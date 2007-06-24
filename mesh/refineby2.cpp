/*
 *  refineby2.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include "mesh.h"
#include <utilities.h>
#include <assert.h>

void mesh::refineby2(const class mesh& inmesh) {
    int i,j,n,sind,tind,v0,v1,count,vnear,err,initialsidenumber;
    TinyVector<FLT,ND> xpt;
    
    /* INPUT MESH MUST HAVE GROWTH FACTOR OF 4 */
    /* BECAUSE OF INTWK USAGE */
     if (!initialized) {
        /* VERTEX STORAGE ALLOCATION */
        allocate_duplicate(0.5,inmesh);
        initialized = 1;
    }
    
    this->copy(inmesh);
    
    /* CALCULATE LOCATION OF NEW INTERIOR POINTS */
    count = nvrtx;
    for(sind=0;sind<nside;++sind) {
    
        if (sd(sind).tri(1) < 0) continue;
        
        v0 = sd(sind).vrtx(0);
        v1 = sd(sind).vrtx(1);
        
        /* MIDPOINT */
        for(n=0;n<ND;++n)
            xpt(n) = 0.5*(vrtx(v0)(n) +vrtx(v1)(n));
                    
        /* INSERT POINT */
        for(n=0;n<ND;++n)
            vrtx(count)(n) = xpt(n);
        
        ++count;
    }
    
    /*	INSERT INTERIOR POINTS */
    for(i=nvrtx;i<count;++i) { 
        qtree.addpt(nvrtx);
        qtree.nearpt(nvrtx,vnear);
        tind = findtri(vrtx(nvrtx),vnear);
        assert(tind > -1);
        err = insert(nvrtx,tind);
        ++nvrtx;
    }
    
    

    
    /* INSERT BOUNDARY POINTS */
    for(i=0;i<nsbd;++i) {
        initialsidenumber = sbdry(i)->nel;
        for(j=0;j<initialsidenumber;++j) {
            sind = sbdry(i)->el(j);
            sbdry(i)->mvpttobdry(j,0.0,vrtx(nvrtx));
            bdry_insert(nvrtx,sind);
            ++nvrtx;
        }
    }
        
    for (i=0;i<nsbd;++i)
        sbdry(i)->reorder();
        
    cnt_nbor();
    
    return;
}

