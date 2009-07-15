/*
 *  intgrt1d.cpp
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_basis.h"

#ifndef BZ_DEBUG
#define f(i) f1[i]
#define rslt(i) rslt1[i]
#endif

void tri_basis::intgrt1d(FLT *rslt1,const FLT *f1) const {
    const int lgpx = gpx, lnmodx = nmodx;
    FLT lcl0;
#ifdef BZ_DEBUG
    const Array<FLT,1> f((FLT *) f1, shape(gpx), neverDeleteData);
    Array<FLT,1> rslt(rslt1, shape(nmodx), neverDeleteData);
#endif
    
    lcl0 = 0.0;
    for(int i=0;i<lgpx;++i)
        lcl0 += gxwtx(1,i)*f(i);
    rslt(0) = lcl0;
    
    lcl0 = 0.0;
    for(int i=0;i<lgpx;++i)
        lcl0 += gxwtx(2,i)*f(i);
    rslt(1) = lcl0;

    for(int n=3;n<lnmodx;++n) {
        lcl0 = 0.0;
        for(int i=0;i<lgpx;++i)
            lcl0 += gxwtx(n,i)*f(i);
        rslt(n-1) = lcl0;
    }
    
    return;
}
    
void tri_basis::intgrtx1d(FLT *rslt1,const FLT *f1) const {
    const int lgpx = gpx, lnmodx = nmodx;
    FLT lcl0;
#ifdef BZ_DEBUG
    const Array<FLT,1> f((FLT *) f1, shape(gpx), neverDeleteData);
    Array<FLT,1> rslt(rslt1, shape(nmodx), neverDeleteData);
#endif
    
    lcl0 = rslt(0);
    for(int i=0;i<lgpx;++i)
        lcl0 -= dgxwtx(1,i)*f(i);
    rslt(0) = lcl0;
    
    lcl0 = rslt(1);
    for(int i=0;i<lgpx;++i)
        lcl0 -= dgxwtx(2,i)*f(i);
    rslt(1) = lcl0;

    for(int n=3;n<lnmodx;++n) {
        lcl0 = rslt(n-1);
        for(int i=0;i<lgpx;++i)
            lcl0 -= dgxwtx(n,i)*f(i);
        rslt(n-1) = lcl0;
    }

    return;
}



