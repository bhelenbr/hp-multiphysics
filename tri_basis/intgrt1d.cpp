/*
 *	intgrt1d.cpp
 *	planar++
 *
 *	Created by helenbrk on Fri Oct 12 2001.
 *	Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_basis.h"

#ifndef BZ_DEBUG
#define f(i) f1[i]
#define rslt(i) rslt1[i]
#endif

template<int _p, int ep> void tri_basis<_p,ep>::intgrt1d(FLT *rslt1, const FLT *f1) const {
	const int lgpx = _gpx, lnmodx = _nmodx;
	FLT lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,1> f((FLT *) f1, shape(_gpx), neverDeleteData);
	Array<FLT,1> rslt(rslt1, shape(_nmodx-1), neverDeleteData);
#endif

	lcl0 = 0.0;
	for(int i=0;i<lgpx;++i)
		lcl0 += _gxwtx(1,i)*f(i);
	rslt(0) = lcl0;
	
	lcl0 = 0.0;
	for(int i=0;i<lgpx;++i)
		lcl0 += _gxwtx(2,i)*f(i);
	rslt(1) = lcl0;

	for(int n=3;n<lnmodx;++n) {
		lcl0 = 0.0;
		for(int i=0;i<lgpx;++i)
			lcl0 += _gxwtx(n,i)*f(i);
		rslt(n-1) = lcl0;
	}
	
	return;
}
	
template<int _p, int ep> void tri_basis<_p,ep>::intgrtx1d(FLT *rslt1, const FLT *f1) const {
	const int lgpx = _gpx, lnmodx = _nmodx;
	FLT lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,1> f((FLT *) f1, shape(_gpx), neverDeleteData);
	Array<FLT,1> rslt(rslt1, shape(_nmodx-1), neverDeleteData);
#endif
	
	lcl0 = rslt(0);
	for(int i=0;i<lgpx;++i)
		lcl0 -= _dgxwtx(1,i)*f(i);
	rslt(0) = lcl0;
	
	lcl0 = rslt(1);
	for(int i=0;i<lgpx;++i)
		lcl0 -= _dgxwtx(2,i)*f(i);
	rslt(1) = lcl0;

	for(int n=3;n<lnmodx;++n) {
		lcl0 = rslt(n-1);
		for(int i=0;i<lgpx;++i)
			lcl0 -= _dgxwtx(n,i)*f(i);
		rslt(n-1) = lcl0;
	}

	return;
}

#ifndef BZ_DEBUG
#undef f
#undef rslt
#endif

