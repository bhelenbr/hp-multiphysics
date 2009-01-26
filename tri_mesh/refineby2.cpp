/*
 *  refineby2.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_mesh.h"
#include <utilities.h>
#include <assert.h>

void tri_mesh::refineby2(const class tri_mesh& inmesh) {
    int i,j,n,sind,tind,p0,p1,count,pnear,err,initialsidenumber,ierr;
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
        ierr = findtri(pnts(npnt),pnear,tind);
        assert(!ierr);
        err = insert(npnt,tind);
        ++npnt;
    }
    
    

    
    /* INSERT BOUNDARY POINTS */
    for(i=0;i<nebd;++i) {
        initialsidenumber = ebdry(i)->nseg;
        for(j=0;j<initialsidenumber;++j) {
            sind = ebdry(i)->seg(j);
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

