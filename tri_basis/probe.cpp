/*
 *  probe.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 18 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"hpbasis.h"
#include<math.h>
#include<utilities.h>

#ifndef BZ_DEBUG
#define lin(n,m) lin1[(n)*stride +m]
#define sin(n,k) sin1[(n)*stride +k]
#define f(i) f1[(i)]
#define dx(i) dx1[(i)]
#define dy(i) dy1[(i)]
#endif

void hpbasis::ptprobe(int nv, FLT *f1, FLT *lin1, int stride) {
    TinyVector<FLT,MXTM> wk0;
    const int bs1 = sm+3, bs2 = 2*sm+3, bint = bm;
    const int lnmodx = nmodx;  
    FLT lcl0;
#ifdef BZ_DEBUG
    Array<FLT,2> lin(lin1, shape(nv,stride), neverDeleteData);
    Array<FLT,1> f(f1, shape(nv), neverDeleteData);
#endif
    
    for(int n=0;n<nv;++n) {
    
        /* SUM ALL S MODE CONTRIBUTIONS */
        /* VERTEX 0            */
        wk0(0) = lin(n,0)*pgn(0);
        
        /* VERTEX 1            */
        lcl0 = lin(n,1)*pgn(1);
        /* SIDE 2        */
        for(int m = bs2; m < bint; ++m )
            lcl0 += lin(n,m)*pgn(m);
        wk0(1) = lcl0;
  
        /* VERTEX 2            */
        lcl0 = lin(n,2)*pgn(2);
        /* SIDE 1        */
        for (int m = bs1; m < bs2; ++m)
            lcl0 += lin(n,m)*pgn(m);
        wk0(2) = lcl0;

        /* LOOP FOR INTERIOR MODES        */
        int ind = bint;
        int bint1 = sm+2-3;
        for(int m = 3; m < bs1; ++m) {
            /* SIDE 0        */
            lcl0 = lin(n,m)*pgn(m);
        
            /* INTERIOR MODES        */
            for(int k = 0; k < bint1; ++k) {
                lcl0 += lin(n,ind+k)*pgn(ind+k);
            }
            ind += bint1--;
            wk0(m) = lcl0;
        }

    /* SUM OVER N X MODES    */      
        lcl0    = 0.0;
        for(int k=0; k < lnmodx; ++k )  
            lcl0 += wk0(k)*pgx(k);
        f(n) = lcl0;
    }
    
    return;
}

void hpbasis::ptprobe_bdry(int nv, FLT *f1, FLT *lin1, int stride) {
    TinyVector<FLT,MXTM> wk0;
    const int bs1 = sm+3, bs2 = 2*sm+3, bint = bm;
    const int lnmodx = nmodx;  
    FLT lcl0;
#ifdef BZ_DEBUG
    Array<FLT,2> lin(lin1, shape(nv,stride), neverDeleteData);
    Array<FLT,1> f(f1, shape(nv), neverDeleteData);
#endif
    
    for(int n=0;n<nv;++n) {
    
        /* SUM ALL S MODE CONTRIBUTIONS */
        /* VERTEX 0            */
        wk0(0) = lin(n,0)*pgn(0);
        
        /* VERTEX 1            */
        lcl0 = lin(n,1)*pgn(1);
        /* SIDE 2        */
        for(int m = bs2; m < bint; ++m )
            lcl0 += lin(n,m)*pgn(m);
        wk0(1) = lcl0;
  
        /* VERTEX 2            */
        lcl0 = lin(n,2)*pgn(2);
        /* SIDE 1        */
        for (int m = bs1; m < bs2; ++m)
            lcl0 += lin(n,m)*pgn(m);
        wk0(2) = lcl0;

        /* SIDE 0        */
        for(int m = 3; m < bs1; ++m) {
            wk0(m) = lin(n,m)*pgn(m);
        }

    /* SUM OVER N X MODES    */      
        lcl0    = 0.0;
        for(int k=0; k < lnmodx; ++k )  
            lcl0 += wk0(k)*pgx(k);
        f(n) = lcl0;
    }
    
    return;
}

void hpbasis::ptprobe_bdry(int nv, FLT *f1, FLT *dx1, FLT *dy1, FLT r, FLT s, FLT *lin1, int stride) {
    TinyVector<FLT,MXTM> wk0,wk1,wk2;
    const int bs1 = sm+3, bs2 = 2*sm+3, bint = bm;
    FLT lcl0, lcl1, lcl2;
    FLT xp1,oeta; 
    FLT x,eta;
#ifdef BZ_DEBUG
    Array<FLT,2> lin(lin1, shape(nv,stride), neverDeleteData);
    Array<FLT,1> f(f1, shape(nv), neverDeleteData);
    Array<FLT,1> dx(dx1, shape(nv), neverDeleteData);
    Array<FLT,1> dy(dy1, shape(nv), neverDeleteData);
#endif
    
    s = MIN(1.0-10.*EPSILON,s);
    x = 2.0*(1+r)/(1-s) -1.0;
    eta = s;
    
    ptvalues_deriv_bdry(x,eta);
    
    for(int n=0;n<nv;++n) {

        /* PART I - sum u*g_mn for each n, s_j    */
        oeta = 2./(1 -eta);
        
        /* VERTEX 0            */
        wk0(0) = lin(n,0)*pgn(0);
        wk1(0) = lin(n,0)*dpgn(0);
        wk2(0) = wk0(0)*oeta;

        /* VERTEX 1            */
        lcl0 = lin(n,1)*pgn(1);
        lcl1 = lin(n,1)*dpgn(1);
        /* SIDE 2        */
        for(int m = bs2; m < bint; ++m ) {
            lcl0 += lin(n,m)*pgn(m);
            lcl1 += lin(n,m)*dpgn(m);
        }  
        wk0(1) = lcl0;
        wk1(1) = lcl1;         
        wk2(1) = lcl0*oeta;
            
  
        /* VERTEX 2            */
        lcl0 = lin(n,2)*pgn(2);
        lcl1 = lin(n,2)*dpgn(2);
        /* SIDE 1        */
        for (int m = bs1; m < bs2; ++m) {
            lcl0 += lin(n,m)*pgn(m);
            lcl1 += lin(n,m)*dpgn(m);
        }            
        wk0(2) = lcl0;
        wk1(2) = lcl1;         
        wk2(2) = lcl0*oeta;

        for(int m = 3; m < bs1; ++m) {
            /* SIDE 0        */
            wk0(m) = lin(n,m)*pgn(m);
            wk1(m) = lin(n,m)*dpgn(m);
            wk2(m) = wk0(m)*oeta;
        }

        /* SUM OVER N AT EACH I,J POINT    */      
        xp1 = 0.5*(1+x);
        lcl0 = 0.0;
        lcl1 = 0.0;
        lcl2 = 0.0;
        for(int k=0; k < nmodx; ++k ) {      
            lcl0 += wk0(k)*pgx(k);
            lcl1 += wk1(k)*pgx(k);
            lcl2 += wk2(k)*dpgx(k);
        }
        f(n) = lcl0;
        dx(n) = lcl2;
        dy(n) = lcl1 +xp1*lcl2;
    }
    return;
}

void hpbasis::ptprobe1d(int nv, FLT *f1, FLT *sin1, int stride) {
    const int pp=p+1;
    FLT lcl0;
#ifdef BZ_DEBUG
    Array<FLT,2> sin(sin1, shape(nv,stride), neverDeleteData);
    Array<FLT,1> f(f1, shape(nv), neverDeleteData);
#endif
    
    for(int n=0;n<nv;++n) {
        lcl0 = 0.0;
        for(int k=0; k < pp; ++k )  
            lcl0 += sin(n,k)*pgx(k);
        f(n) = lcl0;

    }
    
    return;
}

void hpbasis::ptprobe1d(int nv, FLT *f1, FLT *dx1, FLT *sin1, int stride) {
    const int pp = p+1;
    FLT lcl0,lcl1;
#ifdef BZ_DEBUG
    Array<FLT,2> sin(sin1, shape(nv,stride), neverDeleteData);
    Array<FLT,1> f(f1, shape(nv), neverDeleteData);
    Array<FLT,1> dx(dx1, shape(nv), neverDeleteData);
#endif
    
    for(int n=0;n<nv;++n) {
        lcl0 = 0.0;
        lcl1 = 0.0;
        for(int k=0; k < pp; ++k ) {
            lcl0 += sin(n,k)*pgx(k);
            lcl1 += sin(n,k)*dpgx(k);
        }
        f(n) = lcl0;
        dx(n) = lcl1;
    }
    
    return;
}
    
