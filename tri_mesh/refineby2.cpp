/*
 *  refineby2.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_mesh.h"
#include <assert.h>

void tri_mesh::refineby2(const class tri_mesh& inmesh) {
	int i,j,n,sind,tind,p0,p1,count,pnear,err,initialsidenumber;
	bool found;
	TinyVector<FLT,ND> xpt;

	/* INPUT MESH MUST HAVE GROWTH FACTOR OF 4 */
	/* BECAUSE OF gbl->intwk USAGE */
	if (!initialized) {
		/* VERTEX STORAGE ALLOCATION */
		init(inmesh,duplicate,0.5);
	}

	this->copy(inmesh);

	/* CALCULATE LOCATION OF NEW INTERIOR POINTS */
	count = npnt;
	for(sind=0;sind<nseg;++sind) {

		if (seg(sind).tri(1) < 0) continue;

		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);

		/* MIDPOINT */
		for(n=0;n<ND;++n)
			xpt(n) = 0.5*(pnts(p0)(n) +pnts(p1)(n));

		/* INSERT POINT */
		for(n=0;n<ND;++n)
			pnts(count)(n) = xpt(n);

		++count;
	}

	/*	INSERT INTERIOR POINTS */
	for(i=npnt;i<count;++i) {
		qtree.addpt(npnt);
		qtree.nearpt(npnt,pnear);
		found = findtri(pnts(npnt),pnear,tind);
		assert(found);
		err = insert(npnt,tind);
		++npnt;
	}

	/* INSERT BOUNDARY POINTS */
	for(i=0;i<nebd;++i) {
		initialsidenumber = ebdry(i)->nseg;
		for(j=0;j<initialsidenumber;++j) {
			sind = ebdry(i)->seg(j);
			pnts(npnt) = 0.5*(pnts(seg(sind).pnt(0)) +pnts(seg(sind).pnt(1)));
			ebdry(i)->mvpttobdry(j,0.0,pnts(npnt));
			bdry_insert(npnt,sind);
			++npnt;
		}
	}

	for (i=0;i<nebd;++i)
		ebdry(i)->reorder();

	cnt_nbor();

	return;
}

void tri_mesh::refineby2() {
    /* SET FLAGS ETC... */
    setup_for_adapt();

    /* CALCULATE LOCATION OF NEW INTERIOR POINTS */
    int count = npnt;
    TinyVector<FLT,ND> xpt;
    for(int sind=0;sind<nseg;++sind) {

        if (seg(sind).tri(1) < 0) continue;

        int p0 = seg(sind).pnt(0);
        int p1 = seg(sind).pnt(1);

        /* MIDPOINT */
        for(int n=0;n<ND;++n)
            xpt(n) = 0.5*(pnts(p0)(n) +pnts(p1)(n));

        /* INSERT POINT */
        for(int n=0;n<ND;++n)
            pnts(count)(n) = xpt(n);
        
        ++count;
    }

    /*    INSERT INTERIOR POINTS */
    for(int i=npnt;i<count;++i) {
        int pnear, tind;
        qtree.addpt(npnt);
        qtree.nearpt(npnt,pnear);
        bool found = findtri(pnts(npnt),pnear,tind);
        assert(found);
        bool err = insert(npnt,tind);
        assert(!err);
        tri(npnt).info |= PTOUC;
        ++npnt;
    }
    
    /* REFINE FIRST EDGES */
    bdry_refineby2();

    /* REMOVE DELETED ENTITIES */
    cleanup_after_adapt();
    
    /* Calculate halo's for adaptation interpolation */
    calculate_halo();
    
}

void tri_mesh::bdry_refineby2() {
    TinyVector<FLT,tri_mesh::ND> endpt;

    /* REFINE BOUNDARY SIDES */
    for(int bnum=0;bnum<nebd;++bnum) {
        
        /* INSERT BOUNDARY POINTS */
        int initialsidenumber = ebdry(bnum)->nseg;
        for(int j=0;j<initialsidenumber;++j) {
            int sind = ebdry(bnum)->seg(j);
            pnts(npnt) = 0.5*(pnts(seg(sind).pnt(0)) +pnts(seg(sind).pnt(1)));
            ebdry(bnum)->mvpttobdry(j,0.0,pnts(npnt));
            bdry_insert(npnt,sind);
            ++npnt;

            /* UPDATE NEXT & PREV POINTERS */
            ebdry(bnum)->prev(ebdry(bnum)->next(j)) = ebdry(bnum)->nseg-1;
            ebdry(bnum)->next(ebdry(bnum)->nseg-1) = ebdry(bnum)->next(j);
            ebdry(bnum)->next(j) = ebdry(bnum)->nseg-1;
            ebdry(bnum)->prev(ebdry(bnum)->nseg-1) = j;
        }
    }

    return;
}
