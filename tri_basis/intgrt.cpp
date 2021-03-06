/*
 *	intgrt.cpp
 *	planar++
 *
 *	Created by helenbrk on Mon Oct 08 2001.
 *	Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "tri_basis.h"
#include <stdio.h>
#include <stdlib.h>

/* 2D to 1D access */
#ifndef BZ_DEBUG
#define f(i,j) f1[(i)*stride +j]
#define dx(i,j) dx1[(i)*stride +j]
#define dy(i,j) dy1[(i)*stride +j]
#define rslt(i) rslt1[i]
#endif

template<int _p, int ep> void tri_basis<_p,ep>::intgrt(FLT *rslt1, const FLT *f1, int stride) const {
	TinyMatrix<FLT,_nmodx,_gpn> wk0;
	const int bs1 = _sm+3, bs2 = 2*_sm+3, bint = _bm, lsm2 = _sm+2;
	const int lgpn=_gpn,lgpx=_gpx,lnmodx=_nmodx;
	FLT lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,2> f((FLT *) f1, shape(_gpx,stride), neverDeleteData);
	Array<FLT,1> rslt(rslt1, shape(_tm), neverDeleteData);
#endif
		
	for (int n = 0; n < lnmodx; ++n) {
		for (int j = 0; j < lgpn; ++j) {
			lcl0 = f(0,j)*_gxwtx(n,0);
			for(int i = 1; i < lgpx; ++i) {
				lcl0 += f(i,j)*_gxwtx(n,i);
			}
			wk0(n,j) = lcl0;
		}
	}
	
	/* SIDE 0 */
	for (int m=0; m < bs1; ++m) {
		lcl0 = 0.0;
		for (int j=0; j < lgpn; ++j ) 
			lcl0 += _gnwtn(m,j)*wk0(m,j);
		rslt(m) = lcl0;
	}
	
	/* SIDE 1 */
	for (int m=bs1;m<bs2;++m) {
		lcl0 = 0.0;
		for (int j=0; j < lgpn; ++j ) 
			lcl0 += _gnwtn(m,j)*wk0(2,j);
		rslt(m) = lcl0;
	}
	
	/* SIDE 2 */
	for (int m=bs2;m<bint;++m) {
		lcl0 = 0.0;
		for (int j=0; j < lgpn; ++j ) 
			lcl0 += _gnwtn(m,j)*wk0(1,j);
		rslt(m) = lcl0;
	}
	
	int indx = bint;
	for(int m = 3; m < lsm2; ++m) {
		for(int n = 0; n < lsm2-m; ++n) {
			lcl0 = 0.0;
			for (int j=0; j < lgpn; ++j )
				lcl0 += _gnwtn(indx,j)*wk0(m,j);
			rslt(indx) = lcl0;
			++indx;
		}
	}
	
	return;
}

/* WARNING THIS ADDS INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
template<int _p, int ep> void tri_basis<_p,ep>::intgrtrs(FLT *rslt1, const FLT *dx1, const FLT *dy1, int stride) const {
	TinyMatrix<FLT,_nmodx,_gpn> wk0,wk1;
	TinyMatrix<FLT,_gpx,_gpn> wk2;
	const int bs1 = _sm+3, bs2 = 2*_sm+3, bint = _bm, lsm2 = _sm+2;
	const int lgpn=_gpn,lgpx=_gpx,lnmodx=_nmodx;
	FLT lcl0, lcl1;
#ifdef BZ_DEBUG
	const Array<FLT,2> dx((FLT *) dx1, shape(_gpx,stride), neverDeleteData);
	const Array<FLT,2> dy((FLT *) dy1, shape(_gpx,stride), neverDeleteData);
	Array<FLT,1> rslt(rslt1, shape(_tm), neverDeleteData);
#endif
		
	for (int j = 0; j < lgpn; ++j) {
		wk2(0,j) = dx(0,j) +dy(0,j)*_x0(0);
		lcl0 = wk2(0,j)*_dgxwtx(0,0);
		lcl1 = dy(0,j)*_gxwtx(0,0);	
		for(int i = 1; i < lgpx; ++i) {
			wk2(i,j) = dx(i,j) +dy(i,j)*_x0(i);
			lcl0 += wk2(i,j)*_dgxwtx(0,i);
			lcl1 += dy(i,j)*_gxwtx(0,i);
		}
		wk0(0,j) = lcl0;
		wk1(0,j) = lcl1;
	}

	for (int n = 1; n < lnmodx; ++n) {
		for (int j = 0; j < lgpn; ++j) {
			lcl0 = wk2(0,j)*_dgxwtx(n,0);
			lcl1 = dy(0,j)*_gxwtx(n,0);
			for(int i = 1; i < lgpx; ++i) {
				lcl0 += wk2(i,j)*_dgxwtx(n,i);
				lcl1 += dy(i,j)*_gxwtx(n,i);
			}
			wk0(n,j) = lcl0;
			wk1(n,j) = lcl1;
		}
	}
	
	/* SIDE 0 */
	for (int m=0; m < bs1; ++m) {
		lcl0 = rslt(m);
		for (int j=0; j < lgpn; ++j ) 
			lcl0 -= _gnwtnn0(m,j)*wk0(m,j) +_dgnwtn(m,j)*wk1(m,j);
		rslt(m) = lcl0;
	}
	
	for (int j=0; j < lgpn; ++j ) {
		/* SIDE 1 */
		lcl0 = wk0(2,j);
		lcl1 = wk1(2,j);
		for (int m=bs1; m < bs2; ++m) 
			rslt(m) -= _gnwtnn0(m,j)*lcl0 +_dgnwtn(m,j)*lcl1;
		
		/* SIDE 2 */
		lcl0 = wk0(1,j);
		lcl1 = wk1(1,j);
		for (int m=bs2; m < bint; ++m)
			rslt(m) -= _gnwtnn0(m,j)*lcl0 +_dgnwtn(m,j)*lcl1;
	}
	
	int indx = bint;
	for(int m = 3; m < lsm2;++m) {
		int nmax = lsm2-m;
		for (int j=0; j < lgpn; ++j ) {
			lcl0 = wk0(m,j);
			lcl1 = wk1(m,j);
			for(int n = 0; n < nmax; ++n) {
				rslt(indx) -= _gnwtnn0(indx,j)*lcl0+_dgnwtn(indx,j)*lcl1;
				++indx;
			}
			indx -= nmax;
		}
		indx += nmax;
	}
	
	
	return;
}

#ifndef BZ_DEBUG
#undef f
#undef dx
#undef dy
#undef rslt
#endif

