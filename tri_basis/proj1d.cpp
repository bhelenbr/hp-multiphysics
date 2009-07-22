/*
 *	proj1d.cpp
 *	planar++
 *
 *	Created by helenbrk on Fri Oct 12 2001.
 *	Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_basis.h"
// #include "array2d.h"

#ifndef BZ_DEBUG
#define f(i) f1[(i)]
#define dx(i) dx1[(i)]
#define lin(i) lin1[(i)]
#endif

template<int _p, int ep> void tri_basis<_p,ep>::proj1d(const FLT *lin1, FLT *f1, FLT *dx1) const {
	const int lgpx = _gpx, lsm2 = _sm+2;
	FLT lcl0, lcl1;
#ifdef BZ_DEBUG
	Array<FLT,1> lin((double *) lin1, shape(lsm2), neverDeleteData);
	Array<FLT,1> f(f1, shape(_gpx), neverDeleteData);
	Array<FLT,1> dx(dx1, shape(_gpx), neverDeleteData);
#endif
	
	for(int i=0;i<lgpx;++i) {
		lcl0 = lin(0)*_gx(i,1) +lin(1)*_gx(i,2);
		lcl1 = lin(0)*_dgx(i,1) +lin(1)*_dgx(i,2);
		for(int n=2;n<lsm2;++n) {
			lcl0 += lin(n)*_gx(i,n+1);
			lcl1 += lin(n)*_dgx(i,n+1);
		}
		f(i) = lcl0;
		dx(i) = lcl1;
	}
			
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj1d(const FLT *lin1, FLT *f1) const {
	const int lgpx = _gpx, lsm2 = _sm+2;
	FLT lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(lsm2), neverDeleteData);
	Array<FLT,1> f(f1, shape(_gpx), neverDeleteData);
#endif
	
	for(int i=0;i<lgpx;++i) {
		lcl0 = lin(0)*_gx(i,1) +lin(1)*_gx(i,2);
		for(int n=2;n<lsm2;++n) {
			lcl0 += lin(n)*_gx(i,n+1);
		}
		f(i) = lcl0;
	}
			
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj1d(FLT u1, FLT u2, FLT *f1) const {
	const int lgpx = _gpx;
#ifdef BZ_DEBUG
	Array<FLT,1> f(f1, shape(_gpx), neverDeleteData);
#endif
	
	for(int i=0;i<lgpx;++i)
		f(i) = u1*_gx(i,1) +u2*_gx(i,2);
		
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::derivx1d(const FLT *f1, FLT *dx1) const {
	TinyVector<FLT,_gpx> wk0;
	const int lgpx = _gpx;
	FLT lcl0;
#ifdef BZ_DEBUG
	Array<FLT,1> f((FLT *) f1, shape(_gpx), neverDeleteData);
	Array<FLT,1> dx(dx1, shape(_gpx), neverDeleteData);
#endif
	
	for (int i=0;i<lgpx;++i) 
		wk0(i) = f(i)*_dltx(i);
		
	for (int i=0;i<lgpx;++i) {
		lcl0 = dx(i);
		for(int n=0;n<i;++n) 
			lcl0 += (wk0(i) + wk0(n))*_dltx1(i,n);
		for(int n=i+1;n<lgpx;++n) 
			lcl0 += (wk0(i) + wk0(n))*_dltx1(i,n);
		dx(i) = lcl0;
	}

	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj1d_leg(const FLT *lin1, FLT *f1) const {
	int i,m;
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(_sm+2), neverDeleteData);
	Array<FLT,1> f(f1, shape(_sm+1), neverDeleteData);
#endif

	for (i=1;i<_sm+1;++i)
		f(i) = lin(0)*_lgrnge1d(0,i)+lin(1)*_lgrnge1d(1,i);
	
	for(m=2;m<_sm+2;++m)
		for(i=1;i<_sm+1;++i)
			f(i) += lin(m)*_lgrnge1d(m,i);
	
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj1d_leg(FLT u1, FLT u2, FLT *f1) const {
	int i;
#ifdef BZ_DEBUG
	Array<FLT,1> f(f1, shape(_sm+1), neverDeleteData);
#endif
	
	for (i=1;i<_sm+1;++i)
		f(i) = u1*_lgrnge1d(0,i) +u2*_lgrnge1d(1,i);
		
	return;
}

#ifndef BZ_DEBUG
#undef f
#undef dx
#undef lin
#endif


