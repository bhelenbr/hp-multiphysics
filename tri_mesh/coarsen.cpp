#include "tri_mesh.h"
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <assert.h>
#include <typeinfo>


void tri_mesh::coarsen(FLT factor, const class tri_mesh& inmesh) {
	int i,j,n,sind;
	int p0, p1, odd;
	FLT mindist;

	if (!initialized) {
		/* VERTEX STORAGE ALLOCATION */
		init(inmesh,duplicate,1.9);
	}

    for(i=0;i<maxpst;++i)
		tri(i).info = 0;

	/* PREPARE FOR COARSENING */
	for(i=0;i<inmesh.npnt;++i)
		tri_gbl->fltwk(i) = 1.0e8;

	for(i=0;i<inmesh.nseg;++i) {
		p0 = inmesh.seg(i).pnt(0);
		p1 = inmesh.seg(i).pnt(1);
		tri_gbl->fltwk(p0) = MIN(inmesh.distance(p0,p1),tri_gbl->fltwk(p0));
		tri_gbl->fltwk(p1) = MIN(inmesh.distance(p0,p1),tri_gbl->fltwk(p1));
	}

	for(i=0;i<inmesh.npnt;++i)
		tri_gbl->fltwk(i) *= factor; /* SHOULD BE BETWEEN 1.5 and 2.0 */

	npnt = 0;
	nseg = 0;
	ntri  = 0;

	/* USE tri_gbl->intwk TO KEEP TRACK OF INDICES */

	int inmeshbdrysides = 2*inmesh.nseg -3*inmesh.ntri;
	if (maxpst < inmeshbdrysides/2) {
		*gbl->log << "coarse mesh is not big enough: " << inmeshbdrysides/2 << ' ' << maxpst << std::endl;
	}

	/* COARSEN SIDES    */
	for(i=0;i<nebd;++i) {
		ebdry(i)->nseg = 0;
		if (typeid(ebdry(i)) != typeid(inmesh.ebdry(i))) {
			*gbl->log << "can't coarsen into object with different boundaries" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}

		/* CHECK IF FIRST POINT INSERTED*/
		p0 = inmesh.seg(inmesh.ebdry(i)->seg(0)).pnt(0);
		if (tri_gbl->intwk(p0) < 0) {
			for(n=0;n<ND;++n)
				pnts(npnt)(n) = inmesh.pnts(p0)(n);
			tri_gbl->intwk(p0) = npnt;
			seg(nseg).pnt(0) = npnt;
			++npnt;
		}
		else {
			seg(nseg).pnt(0) = tri_gbl->intwk(p0);
		}

		odd = inmesh.ebdry(i)->nseg%2;
		if (odd) {
			for(j=2;j<inmesh.ebdry(i)->nseg/2;j+=2) {
				p0 = inmesh.seg(inmesh.ebdry(i)->seg(j)).pnt(0);
				for(n=0;n<ND;++n)
					pnts(npnt)(n) = inmesh.pnts(p0)(n);
				tri_gbl->intwk(p0) = npnt;
				seg(nseg).pnt(1) = npnt;
				ebdry(i)->seg(ebdry(i)->nseg) = nseg;
				seg(nseg).tri(1) = -1;
				++nseg;
				++ebdry(i)->nseg;
				seg(nseg).pnt(0) = npnt;
				++npnt;
			}

			/* MIDDLE POINT OF ODD NUMBERED SIDE */
			if (inmesh.ebdry(i)->nseg > 1) {
				j = inmesh.ebdry(i)->nseg/2;
				p0 = inmesh.seg(inmesh.ebdry(i)->seg(j)).pnt(0);
				p1 = inmesh.seg(inmesh.ebdry(i)->seg(j)).pnt(1);
				for(n=0;n<ND;++n)
					pnts(npnt)(n) = 0.5*(inmesh.pnts(p0)(n) +inmesh.pnts(p1)(n));
				tri_gbl->intwk(p0) = npnt;
				tri_gbl->intwk(p1)= npnt;
				seg(nseg).pnt(1) = npnt;
				ebdry(i)->seg(ebdry(i)->nseg) = nseg;
				seg(nseg).tri(1) = -1;
				++nseg;
				++ebdry(i)->nseg;
				seg(nseg).pnt(0) = npnt;
				++npnt;
			}

			for(j = inmesh.ebdry(i)->nseg -((inmesh.ebdry(i)->nseg-2)/4)*2;j<inmesh.ebdry(i)->nseg;j+=2) {
				p0 = inmesh.seg(inmesh.ebdry(i)->seg(j)).pnt(0);
				for(n=0;n<ND;++n)
					pnts(npnt)(n) = inmesh.pnts(p0)(n);
				tri_gbl->intwk(p0) = npnt;
				seg(nseg).pnt(1) = npnt;
				ebdry(i)->seg(ebdry(i)->nseg) = nseg;
				seg(nseg).tri(1) = -1;
				++nseg;
				++ebdry(i)->nseg;
				seg(nseg).pnt(0) = npnt;
				++npnt;
			}
		}
		else {
			for(j=2;j<inmesh.ebdry(i)->nseg;j+=2) {
				p0 = inmesh.seg(inmesh.ebdry(i)->seg(j)).pnt(0);
				for(n=0;n<ND;++n)
					pnts(npnt)(n) = inmesh.pnts(p0)(n);
				tri_gbl->intwk(p0) = npnt;
				seg(nseg).pnt(1) = npnt;
				ebdry(i)->seg(ebdry(i)->nseg) = nseg;
				seg(nseg).tri(1) = -1;
				++nseg;
				++ebdry(i)->nseg;
				seg(nseg).pnt(0) = npnt;
				++npnt;
			}
		}

		/* INSERT LAST POINT */
		p0 = inmesh.seg(inmesh.ebdry(i)->seg(inmesh.ebdry(i)->nseg-1)).pnt(1);
		if (tri_gbl->intwk(p0) < 0) {
			for(n=0;n<ND;++n)
				pnts(npnt)(n) = inmesh.pnts(p0)(n);
			tri_gbl->intwk(p0) = npnt;
			seg(nseg).pnt(1) = npnt;
			ebdry(i)->seg(ebdry(i)->nseg) = nseg;
			seg(nseg).tri(1) = -1;
			++nseg;
			++ebdry(i)->nseg;
			++npnt;
		}
		else {
			seg(nseg).pnt(1) = tri_gbl->intwk(p0);
			ebdry(i)->seg(ebdry(i)->nseg) = nseg;
			seg(nseg).tri(1) = -1;
			++nseg;
			++ebdry(i)->nseg;
		}
	}
    
    for(int i=0;i<npnt;++i)
        pnt(i).info = 0;

	/* MOVE VERTEX BDRY INFORMATION */
	for(i=0;i<inmesh.nvbd;++i) {
		vbdry(i)->copy(*inmesh.vbdry(i));
		vbdry(i)->pnt = tri_gbl->intwk(inmesh.vbdry(i)->pnt);
        pnt(vbdry(i)->pnt).info = vbdry(i)->idnum;
	}

	if (maxpst < nseg/2+nseg) {
		*gbl->log << "coarse mesh is not large enough: " << nseg/2 +nseg << ' ' << maxpst << std::endl;
	}


	treeinit();

	if (nseg > maxpst/3) {
		*gbl->log << gbl->idprefix << " use bigger growth factor, lists aren't long enough " << nseg << ' ' << maxpst/3 << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	for(i=0;i<nseg;++i)
		tri_gbl->i2wk_lst1(i) = i+1;
    
    /* For debugging trianguate problems */
    // output("coarsen.d");

	triangulate(nseg);

/* FUNNY WAY OF MARKING BECAUSE OF -1 AS UNMARKED */
/* PSPEC IS 0x4 */
#if ((-1)&0x4)
#define ISSPEC(A) (!((A)&PSPEC))
#define SETSPEC ((-1)&(~PSPEC))
#else
#define ISSPEC(A) (((A)&PSPEC))
#define SETSPEC ((-1)|(PSPEC))
#endif



	/****************************************************/
	/* Boyer-Watson Algorithm to insert interior points */
	/****************************************************/
	/* MARK BOUNDARY SO DON'T GET INSERTED */
	/* FUNNY WAY OF MARKING SO CAN LEAVE tri_gbl->intwk initialized to -1 */
	for(i=0;i<inmesh.nebd;++i) {
		for(j=0;j<inmesh.ebdry(i)->nseg;++j) {
			sind = inmesh.ebdry(i)->seg(j);
			tri_gbl->intwk(inmesh.seg(sind).pnt(0)) = SETSPEC;
			tri_gbl->intwk(inmesh.seg(sind).pnt(1)) = SETSPEC;
		}
	}

	/* maxsrch must be high for triangulated domain with no interior points */
	tri_gbl->maxsrch = MIN(maxpst,inmesh.ntri);

	for(i=0;i<inmesh.npnt;++i) {
		if (ISSPEC(tri_gbl->intwk(i))) continue;

		mindist = qtree.nearpt(inmesh.pnts(i),j);
		if (sqrt(mindist) < tri_gbl->fltwk(i)) continue;

		insert(inmesh.pnts(i));
	}
	cnt_nbor();

	/* reset maxsrch */
	tri_gbl->maxsrch = 1000;

	/* RESET tri_gbl->intwk */
	for(i=0;i<inmesh.nebd;++i) {
		for(j=0;j<inmesh.ebdry(i)->nseg;++j) {
			sind = inmesh.ebdry(i)->seg(j);
			tri_gbl->intwk(inmesh.seg(sind).pnt(0)) = -1;
			tri_gbl->intwk(inmesh.seg(sind).pnt(1)) = -1;
		}
	}
	bdrylabel();
	initlngth();

	/* PRINT SOME GENERAL DEBUGGING INFO */
	*gbl->log << "#" << std::endl << "#COARSE MESH " << std::endl;
	*gbl->log << "#MAXVST:" << maxpst << " POINTS:" << npnt << " SIDES:" << nseg << " ELEMENTS:" << ntri << std::endl;
	/* PRINT BOUNDARY INFO */
	for(i=0;i<nebd;++i)
		*gbl->log << "#" << ebdry(i)->idprefix << " " << ebdry(i)->mytype << " " << ebdry(i)->nseg << std::endl;

	return;
}


void tri_mesh::coarsen2(FLT factor, const class tri_mesh &inmesh, FLT size_reduce) {
	int i;

	if (!initialized) {
		init(inmesh,duplicate,size_reduce);
	}
	copy(inmesh);
	initlngth();
	for(i=0;i<npnt;++i)
		lngth(i) = factor*lngth(i);
	setup_for_adapt();

	bdry_yaber(1.414);

	for(i=0;i<nebd;++i)
		ebdry(i)->comm_exchange(boundary::all,0,boundary::master_slave);

	/* COARSEN MATCHING BOUNDARIES */
	bdry_yaber1();

	/* INTERIOR SIDE COARSENING */
	yaber(1.414);

	cleanup_after_adapt();
	qtree.reinit();  // REMOVES UNUSED QUADS
	
	/* PRINT SOME GENERAL DEBUGGING INFO */
	*gbl->log << "#" << std::endl << "#COARSE MESH " << std::endl;
	*gbl->log << "#MAXVST:" << maxpst << " POINTS:" << npnt << " SIDES:" << nseg << " ELEMENTS:" << ntri << std::endl;
	/* PRINT BOUNDARY INFO */
	for(i=0;i<nebd;++i)
		*gbl->log << "#" << ebdry(i)->idprefix << " " << ebdry(i)->mytype << " " << ebdry(i)->nseg << std::endl;

	return;
}

void tri_mesh::coarsen3() {
	int i,j,p0,tind,sind,node,cnt=0;

	/* SET-UP ADAPTION TRACKING STUFF */
	/* For vertices 0 = untouched, 1 = touched, 2 = deleted, 3 = special */
	/* For sides 0 = untouced, 1 = touched, 2 = deleted */
	/* For triangles 0 = untouched, 1 = touched, 2 = deleted, 3 = searched */
	for(i=0;i<maxpst;++i)
		tri(i).info = 0;

	/* COARSEN SIDE EDGES FIRST */
	for(i=0;i<nebd;++i) {
		*gbl->log << "coarsening boundary " << i << ": type " << ebdry(i)->mytype << " sides " << ebdry(i)->nseg << std::endl;
		for(j=0;j<ebdry(i)->nseg;++j) {
			sind = ebdry(i)->seg(j);
			p0 = seg(sind).pnt(0);
			if (pnt(p0).info > 0) {
				/* COARSEN SIDE */
				collapse(sind,0);
				pnt(p0).info = -1;
				++cnt;
			}
		}
	}


	for(i=0;i<npnt;++i) {
		if (pnt(i).info > 0) {
			tind = pnt(i).tri;
			for(j=0;j<3;++j) {
				if (tri(tind).pnt(j) == i) {
					break;
				}
			}
			j = (j+1)%3;
			sind = tri(tind).seg(j);
			node = (1+tri(tind).sgn(j))/2;
			collapse(sind,node);
			pnt(i).info = -1;
			++cnt;
		}
	}

	cleanup_after_adapt();

	*gbl->log << "#Coarsen finished: " << cnt << " sides coarsened" << std::endl;
	for(i=0;i<nebd;++i)
		*gbl->log << "boundary " << i << ": type " << ebdry(i)->mytype << " sides " << ebdry(i)->nseg << std::endl;

	return;
}



