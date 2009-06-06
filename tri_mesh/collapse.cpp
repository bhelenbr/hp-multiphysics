#include "tri_mesh.h"
#include <utilities.h>
#include <assert.h>
#include <float.h>
#include <blitz/tinyvec-et.h>

void tri_mesh::collapse(int sind, int delt) {
	int ntsrnd, nssrnd, nperim;
	int i,j,vn,pnear,prev,tind,tind1,sind1,stoptri,dir;
	int p0,p1,pt,sd1,sd2,sd3,t1,t2;

	/* gbl->i2wk_lst1 = triangles surrounding delete point */
	/* gbl->i2wk_lst2 = sides surrounding delete point */
	/* -1 indx lists how many */

	/* FIND TRIANGLES / SIDES SURROUNDING ENDPOINT */
	/* EXCLUDES TRIANGLES ADJACENT TO DELETED SIDE */
	/* EXCLUDES SIDES OF TRIANGLES ADJACENT TO DELETED SIDE */
	pnear = seg(sind).pnt(delt);
	tind = pnt(pnear).tri;
	if (tind != seg(sind).tri(0) && tind != seg(sind).tri(1))
		prev = 0;
	else
		prev = 1;
	stoptri = tind;
	dir = 1;
	ntsrnd = 0;
	nssrnd = 0;
	nperim = 0;
	do {
		for(vn=0;vn<3;++vn)
			if (tri(tind).pnt(vn) == pnear) break;

		tind1 = tri(tind).tri((vn +dir)%3);
		if (tind1 < 0) {
			if (dir > 1) break;
			/* REVERSE DIRECTION AND GO BACK TO START */
			++dir;
			tind1 = pnt(pnear).tri;
			prev = 1;
			stoptri = -1;
		}

		if (tind1 != seg(sind).tri(0) && tind1 != seg(sind).tri(1)) {
			gbl->i2wk_lst1(ntsrnd++) = tind1;
			if (!prev) {
				gbl->i2wk_lst2(nssrnd++) = tri(tind).seg((vn +dir)%3);
			}
			prev = 0;
		}
		else {
			prev = 1;
		}

		tind = tind1;

	} while(tind != stoptri);
	gbl->i2wk_lst1(-1) = ntsrnd;
	gbl->i2wk_lst2(-1) = nssrnd;

	/* UPDATE TRI.PNT & SEG.PNT */
	p0 = seg(sind).pnt(1-delt);
	p1 = seg(sind).pnt(delt);

	for(i=0;i<ntsrnd;++i) {
		tind = gbl->i2wk_lst1(i);
		tri(tind).info |= TTOUC;
		for(j=0;j<3;++j) {
			if (tri(tind).pnt(j) == p1) {
				tri(tind).pnt(j) = p0;
				sd3 = tri(tind).seg((j+1)%3);
				tri(sd3).info |= STOUC;
				pt = (1 +tri(tind).sgn((j+1)%3))/2;
				assert(seg(sd3).pnt(pt) == p1 || seg(sd3).pnt(pt) == p0);
				seg(sd3).pnt(pt) = p0;
				sd3 = tri(tind).seg((j+2)%3);
				tri(sd3).info |= STOUC;
				pt = (1 -tri(tind).sgn((j+2)%3))/2;
				assert(seg(sd3).pnt(pt) == p1 || seg(sd3).pnt(pt) == p0);
				seg(sd3).pnt(pt) = p0;
				gbl->i2wk_lst3(nperim++) = tri(tind).seg(j);
				break;
			}
		}
	}

	/* MARK SIDE AS DELETED */
	tri(sind).info |= SDLTE;

	/* CLOSE THE GAP */
	for(i=0;i<2;++i) {
		tind = seg(sind).tri(i);
		if (tind < 0) continue;

		/* MARK TRI AS DELETED TO BE REMOVED LATER */
		tri(tind).info |= TDLTE;

		for(sd1=0;sd1<3;++sd1)
			if(tri(tind).pnt(sd1) == p1) break;

		assert(sd1 != 3);

		if (tri(tind).seg((sd1+1)%3) == sind)
			sd2 = (sd1+2)%3;
		else
			sd2 = (sd1+1)%3;

		sind1 = tri(tind).seg(sd2);
		if (seg(sind1).tri(1) < 0) {
			/* BOUNDARY SIDE SO DELETE OTHER ONE */
			/* ANOMALY OF TRIANGLE WITH 2 EDGES ON BOUNDARY */
			tri(sind1).info |= STOUC;
			if (seg(sind1).pnt(0) == p1)
				seg(sind1).pnt(0) = p0;
			if (seg(sind1).pnt(1) == p1)
				seg(sind1).pnt(1) = p0;
			sind1 = sd2;
			sd2 = sd1;
			sd1 = sind1;
			sind1 = tri(tind).seg(sd2);
		}
		tri(sind1).info |= SDLTE;
		gbl->i2wk_lst2(nssrnd++) = sind1;

		t1 = tri(tind).tri(sd1);
		t2 = tri(tind).tri(sd2);

		/* UPDATE TTRI FOR T1 */
		if (t1 > -1) {
			for(j=0;j<3;++j) {
				if(tri(t1).tri(j) == tind) {
					tri(t1).tri(j) = t2;
					break;
				}
			}
			for(j=0;j<3;++j)
				pnt(tri(tind).pnt(j)).tri = t1;
		}
		/* UPDATE STRI FOR KEPT SIDE */
		sind1 = tri(tind).seg(sd1);
		pt = (1 -tri(tind).sgn(sd1))/2;
		seg(sind1).tri(pt) = t2;
		gbl->i2wk_lst3(nperim++) = sind1;

		/* UPDATE TTRI/TSIDE FOR T2 */
		if (t2 > -1) {
			for(j=0;j<3;++j) {
				if(tri(t2).tri(j) == tind) {
					tri(t2).tri(j) = t1;
					tri(t2).seg(j) = sind1;
					tri(t2).sgn(j) = tri(tind).sgn(sd1);
					break;
				}
			}
			for(j=0;j<3;++j)
				pnt(tri(tind).pnt(j)).tri = t2;
		}
	}

	/* NEED TO REMOVE LEFTOVERS */
	qtree.dltpt(p1);
	tri(p1).info |= PDLTE;

	/* SWAP AFFECTED SIDES */
	swap(gbl->i2wk_lst2(-1),&gbl->i2wk_lst2(0));
	gbl->i2wk_lst2(-1) = nssrnd;
	gbl->i2wk_lst3(-1) = nperim;

	return;
}

/* DELETE UNREFERENCED TRIANGLE */
void tri_mesh::dlttri(int tind) {
	int i,j,p0,t1,sind,flip;


	/* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED TRI */
	while(tri(ntri-1).info&TDLTE)
		--ntri;

	if (ntri <=  tind)  {
		ntri = tind;
		return;
	}

	--ntri;

	for (j=0;j<3;++j) {
		p0 = tri(ntri).pnt(j);
		tri(tind).pnt(j) = p0;
		pnt(p0).tri = tind;
		tri(tind).seg(j) = tri(ntri).seg(j);
		tri(tind).sgn(j) = tri(ntri).sgn(j);
		sind = tri(ntri).seg(j);
		flip = (1 -tri(ntri).sgn(j))/2;
		seg(sind).tri(flip) = tind;

		t1 = tri(ntri).tri(j);
		tri(tind).tri(j) = t1;
		if (t1 > -1) {
			for(i=0;i<3;++i) {
				if(tri(t1).tri(i) == ntri) {
					tri(t1).tri(i) = tind;
					break;
				}
			}
		}
	}

	if (tri(ntri).info&TTOUC) {
		tri(tind).info = -2;
		updatetdata(tind);
	}
	else {
		tri(tind).info = ntri;
		movetdata(ntri,tind);
	}
	/* THIS IS TO PREVENT ROLL BACK NEAR END */
	tri(tind).info &= ~TDLTE;

	tri(ntri).info = tind;

	return;
}

/* DELETE UNREFERENCED SIDE */
void tri_mesh::dltseg(int sind) {
	int j,k,tind;

	/* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED SIDE */
	while(tri(nseg-1).info&SDLTE)
		--nseg;

	if (nseg <= sind) {
		nseg = sind;
		return;
	}

	/* DELETE SIDE */
	--nseg;

	seg(sind).pnt(0) = seg(nseg).pnt(0);
	seg(sind).pnt(1) = seg(nseg).pnt(1);
	seg(sind).tri(0) = seg(nseg).tri(0);
	seg(sind).tri(1) = seg(nseg).tri(1);

	for(j=0;j<2;++j) {
		tind = seg(nseg).tri(j);
		if (tind > -1) {
			for(k=0;k<3;++k)
				if (tri(tind).seg(k) == nseg) break;
			tri(tind).seg(k) = sind;
		}
	}

	if (tri(nseg).info&STOUC) {
		seg(sind).info = -2;
		updatesdata(sind);
	}
	else {
		seg(sind).info = nseg;
		movesdata(nseg,sind);
	}
	seg(nseg).info = sind;

	return;
}

/* DELETE UNREFERENCED VERTEX */
void tri_mesh::dltpnt(int p0) {
	int vn,stoptri,dir;
	int tind, sind, flip;

	/* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED VERTEX */
	while(tri(npnt-1).info&PDLTE) {
		--npnt;
	}

	if (npnt <= p0) {
		npnt = p0;
		return;
	}

	--npnt;

	pnts(p0)(0) = pnts(npnt)(0);
	pnts(p0)(1) = pnts(npnt)(1);
	lngth(p0) = lngth(npnt);
	pnt(p0).nnbor = pnt(npnt).nnbor;
	pnt(p0).tri = pnt(npnt).tri;

	qtree.movept(npnt,p0);

	tind = pnt(npnt).tri;
	stoptri = tind;
	dir = 1;
	do {
		for(vn=0;vn<3;++vn)
			if (tri(tind).pnt(vn) == npnt) break;

		assert(vn != 3);

		tri(tind).pnt(vn) = p0;
		sind = tri(tind).seg((vn +dir)%3);
		flip = (1 +(3-2*dir)*tri(tind).sgn((vn +dir)%3))/2;
		assert(seg(sind).pnt(flip) == npnt);
		seg(sind).pnt(flip) = p0;

		tind = tri(tind).tri((vn +dir)%3);
		if (tind < 0) {
			if (dir > 1) break;
			/* REVERSE DIRECTION AND GO BACK TO START */
			++dir;
			tind = pnt(npnt).tri;
			for(vn=0;vn<3;++vn) {
				if (tri(tind).pnt(vn) == p0) {
					tri(tind).pnt(vn) = npnt;
					break;
				}
			}
			stoptri = -1;
		}

	} while(tind != stoptri);

	if (tri(npnt).info&PTOUC) {
		pnt(p0).info = -2;
		updatepdata(p0);
	}
	else {
		pnt(p0).info = npnt;
		movepdata(npnt,p0);
	}
	pnt(npnt).info = p0;

	return;
}
