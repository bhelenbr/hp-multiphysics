/*
 *  proj1d.cpp
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "hpbasis.h"
// #include "array2d.h"

#ifndef BZ_DEBUG
#define f(i) f1[(i)]
#define dx(i) dx1[(i)]
#define lin(i) lin1[(i)]
#endif

void hpbasis::proj1d(FLT *lin1, FLT *f1, FLT *dx1) {
    const int lgpx = gpx, lsm2 = sm+2;
    FLT lcl0, lcl1;
#ifdef BZ_DEBUG
    Array<FLT,1> lin(lin1, shape(lsm2), neverDeleteData);
    Array<FLT,1> f(f1, shape(gpx), neverDeleteData);
    Array<FLT,1> dx(dx1, shape(gpx), neverDeleteData);
#endif
    
    for(int i=0;i<lgpx;++i) {
        lcl0 = lin(0)*gx(i,1) +lin(1)*gx(i,2);
        lcl1 = lin(0)*dgx(i,1) +lin(1)*dgx(i,2);
        for(int n=2;n<lsm2;++n) {
            lcl0 += lin(n)*gx(i,n+1);
            lcl1 += lin(n)*dgx(i,n+1);
        }
        f(i) = lcl0;
        dx(i) = lcl1;
    }
            
    return;
}

void hpbasis::proj1d(FLT *lin1, FLT *f1) {
    const int lgpx = gpx, lsm2 = sm+2;
    FLT lcl0;
#ifdef BZ_DEBUG
    Array<FLT,1> lin(lin1, shape(lsm2), neverDeleteData);
    Array<FLT,1> f(f1, shape(gpx), neverDeleteData);
#endif
    
    for(int i=0;i<lgpx;++i) {
        lcl0 = lin(0)*gx(i,1) +lin(1)*gx(i,2);
        for(int n=2;n<lsm2;++n) {
            lcl0 += lin(n)*gx(i,n+1);
        }
        f(i) = lcl0;
    }
            
    return;
}

void hpbasis::proj1d(FLT u1, FLT u2, FLT *f1) {
    const int lgpx = gpx;
#ifdef BZ_DEBUG
    Array<FLT,1> f(f1, shape(gpx), neverDeleteData);
#endif
    
    for(int i=0;i<lgpx;++i)
        f(i) = u1*gx(i,1) +u2*gx(i,2);
        
    return;
}

void hpbasis::derivx1d(FLT *f1, FLT *dx1) {
    TinyVector<FLT,MXGP> wk0;
    const int lgpx = gpx;
    FLT lcl0;
#ifdef BZ_DEBUG
    Array<FLT,1> f(f1, shape(gpx), neverDeleteData);
    Array<FLT,1> dx(dx1, shape(gpx), neverDeleteData);
#endif
    
    for (int i=0;i<lgpx;++i) 
        wk0(i) = f(i)*dltx(i);
        
    for (int i=0;i<lgpx;++i) {
        lcl0 = dx(i);
        for(int n=0;n<i;++n) 
            lcl0 += (wk0(i) + wk0(n))*dltx1(i,n);
        for(int n=i+1;n<lgpx;++n) 
            lcl0 += (wk0(i) + wk0(n))*dltx1(i,n);
        dx(i) = lcl0;
    }

    return;
}

void hpbasis::proj1d_leg(FLT *lin1, FLT *f1) {
    int i,m;
#ifdef BZ_DEBUG
    Array<FLT,1> lin(lin1, shape(sm+2), neverDeleteData);
    Array<FLT,1> f(f1, shape(sm+1), neverDeleteData);
#endif

    for (i=1;i<sm+1;++i)
        f(i) = lin(0)*lgrnge1d(0,i)+lin(1)*lgrnge1d(1,i);
    
    for(m=2;m<sm+2;++m)
        for(i=1;i<sm+1;++i)
            f(i) += lin(m)*lgrnge1d(m,i);
    
    return;
}

void hpbasis::proj1d_leg(FLT u1, FLT u2, FLT *f1) {
    int i;
#ifdef BZ_DEBUG
    Array<FLT,1> f(f1, shape(sm+1), neverDeleteData);
#endif
    
    for (i=1;i<sm+1;++i)
        f(i) = u1*lgrnge1d(0,i) +u2*lgrnge1d(1,i);
        
    return;
}

