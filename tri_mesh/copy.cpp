#include "tri_mesh.h"
#include<assert.h>

void tri_mesh::copy(const tri_mesh& tgt) {
	int i,n;

	if (!initialized) {
		init(tgt);
	}
	else {
		/* CHECK IF BIG ENOUGH */
		if (tgt.nseg > maxpst) {
			*gbl->log << "mesh is too big to copy" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}

	/* COPY VERTEX INFO OVER */
	npnt = tgt.npnt;
	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			pnts(i)(n) = tgt.pnts(i)(n);

	for(i=0;i<npnt;++i) {
		lngth(i) = tgt.lngth(i);
		pnt(i).info = tgt.pnt(i).info;
		pnt(i).nnbor = tgt.pnt(i).nnbor;
		pnt(i).tri = tgt.pnt(i).tri;
	}

	/* COPY VERTEX BOUNDARY INFO */
	for(i=0;i<nvbd;++i)
		vbdry(i)->copy(*tgt.vbdry(i));

/* COPY SIDE INFORMATION */
	nseg = tgt.nseg;
	for(i=0;i<nseg;++i) {
		for(n=0;n<2;++n)
			seg(i).pnt(n) = tgt.seg(i).pnt(n);
		for(n=0;n<2;++n)
			seg(i).tri(n) = tgt.seg(i).tri(n);
		seg(i).info = tgt.seg(i).info;
	}

	/* COPY SIDE BOUNDARY INFO */
	for(i=0;i<nebd;++i)
		ebdry(i)->copy(*tgt.ebdry(i));

	/* COPY ELEMENT DATA */
	ntri = tgt.ntri;
	for(i=0;i<ntri;++i) {
		for(n=0;n<3;++n) {
			tri(i).pnt(n) = tgt.tri(i).pnt(n);
			tri(i).tri(n) = tgt.tri(i).tri(n);
			tri(i).seg(n) = tgt.tri(i).seg(n);
			tri(i).sgn(n) = tgt.tri(i).sgn(n);
		}
		tri(i).info = tgt.tri(i).info;
	}

	qtree.copy(tgt.qtree);
	qtree.change_vptr(pnts);

	return;
}

void tri_mesh::append(const tri_mesh &z) {
	int vrt,sind,flip;
	int nvrtxold, nsideold,ntriold;
	int sind1,tind1,v1b;
	int v1a = 0; // To avoid may be used uninitialized warning
	int sind2,tind2,v2a,v2b;

	for(int i=0;i<z.npnt;++i) {
		for(int n=0;n<ND;++n)
			pnts(i+npnt)(n) = z.pnts(i)(n);
	}

	for(int i=0;i<z.npnt;++i)
		lngth(i+npnt) = z.lngth(i);

	for(int i=0;i<z.npnt;++i) {
		pnt(i+npnt).tri = z.pnt(i).tri +ntri;
	}

	/* MOVE BOUNDARY INFO */
	if (z.nvbd > 0) vbdry.resizeAndPreserve(nvbd+z.nvbd);
	for(int i=0;i<z.nvbd;++i) {
		vbdry(nvbd) = z.vbdry(i)->create(*this);
		vbdry(nvbd)->alloc(4);
		vbdry(nvbd)->pnt = z.vbdry(i)->pnt +npnt;
		++nvbd;
	}

	for(int i=0;i<z.nseg;++i) {
		for(int n=0;n<2;++n) {
			seg(i+nseg).pnt(n) = z.seg(i).pnt(n) +npnt;
			seg(i+nseg).tri(n) = z.seg(i).tri(n) +ntri;
		}
	}

	/* MOVE BOUNDARY INFO */
	ebdry.resizeAndPreserve(nebd+z.nebd);
	for(int i=0;i<z.nebd;++i) {
		ebdry(nebd++) = z.ebdry(i)->create(*this);
		ebdry(nebd-1)->alloc(z.ebdry(i)->maxseg);
		ebdry(nebd-1)->nseg = z.ebdry(i)->nseg;
		for(int j=0;j<z.ebdry(i)->nseg;++j) {
			ebdry(nebd-1)->seg(j) = z.ebdry(i)->seg(j) +nseg;
		}
	}

	/* MOVE TRI INFO */
	for(int i=0;i<z.ntri;++i) {
		for(int n=0;n<3;++n) {
			tri(i+ntri).pnt(n) = z.tri(i).pnt(n) +npnt;
			tri(i+ntri).seg(n) = z.tri(i).seg(n)+nseg;
			tri(i+ntri).sgn(n) = z.tri(i).sgn(n);
			tri(i+ntri).tri(n) = z.tri(i).tri(n)+ntri;
		}
	}

	nvrtxold = npnt;
	nsideold = nseg;
	ntriold = ntri;

	npnt += z.npnt;
	nseg += z.nseg;
	ntri += z.ntri;
	
	treeinit();

	bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT

	setup_for_adapt();

	/* FIND MATCHING COMMUNICATION BOUNDARY */
	for(int i=0;i<nebd -z.nebd;++i) {
		if (!ebdry(i)->is_comm()) continue;

		for(int j=0;j<z.nebd;++j) {
			if (ebdry(i)->idnum == z.ebdry(j)->idnum) {
				int bseg = ebdry(i)->nseg;
				
				assert(bseg == z.ebdry(j)->nseg);

				/* DO FIRST POINT OF EDGE (IF NOT ALREADY DONE) */
				sind1 = ebdry(i)->seg(bseg-1);
				tind1 = seg(sind1).tri(0);
				v1b = seg(sind1).pnt(1);
				sind2 = z.ebdry(j)->seg(0);
				tind2 = z.seg(sind2).tri(0) +ntriold;
				v2a = z.seg(sind2).pnt(0) +nvrtxold;
				if (!(tri(v2a).info&PDLTE)) {
					tri(v2a).info |= PDLTE;
					qtree.dltpt(v2a);
					do {
						for(vrt=0;vrt<3;++vrt)
							if (tri(tind2).pnt(vrt) == v2a) break;
						assert(vrt != 3);

						tri(tind2).pnt(vrt) = v1b;
						sind = tri(tind2).seg((vrt +1)%3);
						flip = (1 +tri(tind2).sgn((vrt +1)%3))/2;
						assert(seg(sind).pnt(flip) == v2a);
						seg(sind).pnt(flip) = v1b;
						tind2 = tri(tind2).tri((vrt +1)%3);
					} while(tind2 > 0);
				}


				/* MERGE SIDES AND VERTICES ALONG MATCHING BOUNDARIES */
				for(int k=0;k<bseg;++k) {
					/* Go backwards on this boundary */
					sind1 = ebdry(i)->seg(bseg-k-1);
					tind1 = seg(sind1).tri(0);
					v1a = seg(sind1).pnt(0);
					/* and forwards on this one */
					sind2 = z.ebdry(j)->seg(k);
					tind2 = z.seg(sind2).tri(0) +ntriold;
					v2b = z.seg(sind2).pnt(1) +nvrtxold;

					/* Get neighboring side and tri information filled in for first mesh */
					for(vrt=0;vrt<3;++vrt)
						if (tri(tind1).seg(vrt) == sind1) break;
					assert(vrt < 3);
					tri(tind1).tri(vrt) = tind2;
					seg(sind1).tri(1) = tind2;
					sind2 += nsideold;

					/* Get neighboring side and tri information filled in for second mesh */
					for(vrt=0;vrt<3;++vrt)
						if (tri(tind2).seg(vrt) == sind2) break;
					assert(vrt < 3);
					tri(tind2).seg(vrt) = sind1;
					tri(tind2).sgn(vrt) = -1;
					tri(tind2).tri(vrt) = tind1;
					tri(sind2).info |= SDLTE;
					
					/* Change all references to v2b to v1a */
					if (!(tri(v2b).info&PDLTE)) {
						vrt = (vrt+2)%3;
						tri(v2b).info |= PDLTE;
						qtree.dltpt(v2b);
						/* Find all segs/tri's referring to v2b and change to v1a */
						for(;;) {
							/* First one */
							tri(tind2).pnt(vrt) = v1a;
							sind = tri(tind2).seg((vrt +2)%3);
							flip = (1 -tri(tind2).sgn((vrt +2)%3))/2;
							assert(seg(sind).pnt(flip) == v2b);
							seg(sind).pnt(flip) = v1a;
							
							/* find next possibility */
							tind2 = tri(tind2).tri((vrt +2)%3);
							if (tind2 < 0) break;

							for(vrt=0;vrt<3;++vrt)
								if (tri(tind2).pnt(vrt) == v2b) break;
							assert(vrt != 3);
						}
					}
				}
				
//				/* Delete 1st boundary */
//				delete ebdry(i);
//				for(int k=i;k<nebd-1;++k)
//					ebdry(k) = ebdry(k+1);
//				nebd -= 1;
//				
//				/* Delete 2nd boundary */
//				delete ebdry(nebd-z.nebd+j);
//				for(int k=nebd-z.nebd+j;k<nebd-1;++k)
//					ebdry(k) = ebdry(k+1);
//				nebd -= 1;
				
				
				
				/* Move deleted ebdry to end then delete */
				/* Avoid deleting z halo mesh before being done with it. */
				edge_bdry *temp = ebdry(i);
				for(int k=i;k<nebd-1;++k)
					ebdry(k) = ebdry(k+1);
				ebdry(nebd-1) = temp;
				nebd -= 1;
				
				/* Move 2nd boundary to end */
				temp = ebdry(nebd-z.nebd+j);
				for(int k=nebd-z.nebd+j;k<nebd;++k)
					ebdry(k) = ebdry(k+1);
				ebdry(nebd) = temp;
				nebd -= 1;
				
				goto DONE;
			}
		}
	}
	
DONE:
	/* DELETE VERTEX BOUNDARY CONDITION THAT REFERS TO DELETED POINT */
	for(int i=nvbd-z.nvbd;i<nvbd;++i) {
		if (tri(vbdry(i)->pnt).info&PDLTE) {
			delete vbdry(i);
			for(int j=i+1;j<nvbd;++j)
				vbdry(j-1) =vbdry(j);
			--nvbd;
		}
	}
	
	delete ebdry(nebd);
	delete ebdry(nebd+1);
	
	// cleanup_after_adapt();

	return;
}


