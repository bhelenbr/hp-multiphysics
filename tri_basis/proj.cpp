/*
 *	proj.cpp
 *	planar++
 *
 *	Created by helenbrk on Mon Oct 08 2001.
 *	Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "tri_basis.h"
#include <stdio.h>

#ifndef BZ_DEBUG
#define lin(i) lin1[i]
#define f(i,j) f1[(i)*stride +j]
#define dx(i,j) dx1[(i)*stride +j]
#define dy(i,j) dy1[(i)*stride +j]
#define dnrm(i) dnrm1[i]
#define dtan(i) dtan1[i]
#endif

template<int _p, int ep> void tri_basis<_p,ep>::proj(const FLT *lin1, FLT *f1, FLT *dx1, FLT *dy1, int stride) const {
	TinyMatrix<FLT,_gpn,_nmodx> wk0,wk1,wk2;
	const int bs1 = _sm+3, bs2 = 2*_sm+3, bint = _bm;
	const int lgpx = _gpx, lgpn = _gpn, lnmodx = _nmodx;  
	FLT lcl0, lcl1, lcl2;
	FLT xp1,oeta; 
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(_tm), neverDeleteData);
	Array<FLT,2> f(f1, shape(_gpx,stride), neverDeleteData);
	Array<FLT,2> dx(dx1, shape(_gpx,stride), neverDeleteData);
	Array<FLT,2> dy(dy1, shape(_gpx,stride), neverDeleteData);
#endif
	
	/* DETERMINE U VALUES, GRAD U VALUES
		AT COLLOCATION POINTS
		SUM HAT(U) FOR DU/DX AND DU/DY
	*/
	
	/* GENERAL FORMULA
		dg/dr = 2.0/(1-n) g(s) dg/dx
		dg/ds = g(x)dg/dn +(1+x)/2 dg/dr = g(x)dg/dn +(1+x)/(1-s) g(s) dg/dx 
	*/

	/* PART I - sum u*g_mn for each n, s_j	  */
	for(int j = 0; j < lgpn; ++j ) {
		oeta = _n0(j);
		
		/* VERTEX 0 */
		wk0(j,0) = lin(0)*_gn(j,0);
		wk1(j,0) = lin(0)*_dgn(j,0);
		wk2(j,0) = wk0(j,0)*oeta;
		

		/* VERTEX 1 */
		lcl0 = lin(1)*_gn(j,1);
		lcl1 = lin(1)*_dgn(j,1);
		/* SIDE 2 */
		for(int m = bs2; m < bint; ++m ) {
			lcl0 += lin(m)*_gn(j,m);
			lcl1 += lin(m)*_dgn(j,m);
		}			 
		wk0(j,1) = lcl0;
		wk1(j,1) = lcl1;
		wk2(j,1) = lcl0*oeta;
		

		/* VERTEX 2			   */
		lcl0 = lin(2)*_gn(j,2);
		lcl1 = lin(2)*_dgn(j,2);
		/* SIDE 1		 */
		for (int m = bs1; m < bs2; ++m) {
			lcl0 += lin(m)*_gn(j,m);
			lcl1 += lin(m)*_dgn(j,m);
		}		   
		wk0(j,2) = lcl0;
		wk1(j,2) = lcl1; 
		wk2(j,2) = lcl0*oeta;

		int ind = bint;
		int bint1 = _sm+2-3;
		for(int m = 3; m < bs1; ++m) {
			/* SIDE 0		 */
			lcl0 = lin(m)*_gn(j,m);
			lcl1 = lin(m)*_dgn(j,m);
		
			/* INTERIOR MODES		 */
			for(int n = 0; n < bint1; ++n) {
				lcl0 += lin(ind+n)*_gn(j,ind+n);
				lcl1 += lin(ind+n)*_dgn(j,ind+n);
			}
			ind += bint1--;
			wk0(j,m) = lcl0;
			wk1(j,m) = lcl1;
			wk2(j,m) = lcl0*oeta;
		}
	}

	/* SUM OVER N AT EACH I,J POINT	   */  
	 for (int i=0; i < lgpx; ++i ) {
		xp1 = _x0(i);
		for (int j=0; j < lgpn; ++j) {
			lcl0 = wk0(j,0)*_gx(i,0);
			lcl1 = wk1(j,0)*_gx(i,0);
			lcl2 = wk2(j,0)*_dgx(i,0);
			
			for(int n=1; n < lnmodx; ++n ) {	  
				lcl0 += wk0(j,n)*_gx(i,n);
				lcl1 += wk1(j,n)*_gx(i,n);
				lcl2 += wk2(j,n)*_dgx(i,n);
			}
			f(i,j)	= lcl0;
			dy(i,j) = lcl1 +xp1*lcl2;
			dx(i,j) = lcl2;
		}
	}

	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj(const FLT *lin1, FLT *f1, int stride) const {
	TinyMatrix<FLT,_gpn,_nmodx> wk0;
	const int bs1 = _sm+3, bs2 = 2*_sm+3, bint = _bm;
	const int lgpx = _gpx, lgpn = _gpn, lnmodx = _nmodx;  
	FLT lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(_tm), neverDeleteData);
	Array<FLT,2> f(f1, shape(_gpx,stride), neverDeleteData);
#endif

	/* DETERMINE U VALUES, GRAD U VALUES
		AT COLLOCATION POINTS
		SUM HAT(U) FOR DU/DX AND DU/DY
	*/

	/* PART I - sum u*g_mn for each n, s_j	  */
	for(int j = 0; j < lgpn; ++j ) {
		
		/* VERTEX 0			   */
		wk0(j,0) = lin(0)*_gn(j,0);

		/* VERTEX 1			   */
		lcl0 = lin(1)*_gn(j,1);
		/* SIDE 2		 */
		for(int m = bs2; m < bint; ++m )
			lcl0 += lin(m)*_gn(j,m);
		wk0(j,1) = lcl0;

		/* VERTEX 2			   */
		lcl0 = lin(2)*_gn(j,2);
		/* SIDE 1		 */
		for (int m = bs1; m < bs2; ++m)
			lcl0 += lin(m)*_gn(j,m);
		wk0(j,2) = lcl0;


		/* LOOP FOR INTERIOR MODES		  */
		int ind = bint;
		int bint1 = _sm+2-3;
		for(int m = 3; m < bs1; ++m) {
			/* SIDE 0		 */
			lcl0 = lin(m)*_gn(j,m);
			
		  /* INTERIOR MODES		   */
			for(int n = 0; n < bint1; ++n) {
				lcl0 += lin(ind+n)*_gn(j,ind+n);
			}
			ind += bint1--;
			wk0(j,m) = lcl0;
		}
	}
		 
	/* SUM OVER N AT EACH I,J POINT	   */	   
	 for (int i=0; i < lgpx; ++i ) {
		 for (int j=0; j < lgpn; ++j) {
			 lcl0 = 0.0;

			for(int n=0; n < lnmodx; ++n )	
				lcl0 += wk0(j,n)*_gx(i,n);
				
			f(i,j) = lcl0;
		}
	}

	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj(FLT u1, FLT u2, FLT u3, FLT *f1, int stride) const {
	const int lgpx = _gpx, lgpn = _gpn;
#ifdef BZ_DEBUG
	Array<FLT,2> f(f1, shape(_gpx,stride), neverDeleteData);
#endif

	for (int i=0; i < lgpx; ++i ) 
		for (int j=0; j < lgpn; ++j ) 
			f(i,j) = u1*_gx(i,0)*_gn(j,0) +u2*_gx(i,1)*_gn(j,1) +u3*_gx(i,2)*_gn(j,2);
	
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj_bdry(const FLT *lin1, FLT *f1, FLT *dx1, FLT *dy1, int stride) const {
	TinyMatrix<FLT,_gpn,_nmodx> wk0,wk1,wk2;
	const int bs1 = _sm+3, bs2 = 2*_sm+3, bint = _bm;
	const int lgpx = _gpx, lgpn = _gpn, lnmodx = _nmodx;
	FLT lcl0, lcl1, lcl2;
	FLT xp1,oeta;
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(_tm), neverDeleteData);
	Array<FLT,2> f(f1, shape(_gpx,stride), neverDeleteData);
	Array<FLT,2> dx(dx1, shape(_gpx,stride), neverDeleteData);
	Array<FLT,2> dy(dy1, shape(_gpx,stride), neverDeleteData);
#endif
				 
	/* DETERMINE U VALUES, GRAD U VALUES
		AT COLLOCATION POINTS
		SUM HAT(U) FOR DU/DX AND DU/DY
	*/

	/* PART I - sum u*g_mn for each n, s_j	  */
	for(int j = 0; j < lgpn; ++j ) {
		oeta = _n0(j);
		
		/* VERTEX 0			   */
		wk0(j,0) = lin(0)*_gn(j,0);
		wk1(j,0) = lin(0)*_dgn(j,0);
		wk2(j,0) = wk0(j,0)*oeta;
		
		/* VERTEX 1			   */
		lcl0 = lin(1)*_gn(j,1);
		lcl1 = lin(1)*_dgn(j,1);
		/* SIDE 2		 */
		for(int m = bs2; m < bint; ++m ) {
			lcl0 += lin(m)*_gn(j,m);
			lcl1 += lin(m)*_dgn(j,m);
		}	   
		wk0(j,1) = lcl0;
		wk1(j,1) = lcl1;
		wk2(j,1) = lcl0*oeta;
			

		/* VERTEX 2			   */
		lcl0 = lin(2)*_gn(j,2);
		lcl1 = lin(2)*_dgn(j,2);
		/* SIDE 1		 */
		for (int m = bs1; m < bs2; ++m) {
			lcl0 += lin(m)*_gn(j,m);
			lcl1 += lin(m)*_dgn(j,m);
		}			 
		wk0(j,2) = lcl0;
		wk1(j,2) = lcl1;
		wk2(j,2) = lcl0*oeta;

		for(int m = 3; m < bs1; ++m) {
			/* SIDE 0		 */
			wk0(j,m) = lin(m)*_gn(j,m);
			wk1(j,m) = lin(m)*_dgn(j,m);
			wk2(j,m) = wk0(j,m)*oeta;
		}
	}
		 
	/* SUM OVER N AT EACH I,J POINT	   */	   
	 for (int i=0; i < lgpx; ++i ) {
		xp1 = _x0(i);
		for (int j=0; j < lgpn; ++j) {
			lcl0 = 0.0;
			lcl1 = 0.0;
			lcl2 = 0.0;

			for(int n=0; n < lnmodx; ++n ) {	  
				lcl0 += wk0(j,n)*_gx(i,n);
				lcl1 += wk1(j,n)*_gx(i,n);
				lcl2 += wk2(j,n)*_dgx(i,n);
			}
			f(i,j)	= lcl0;
			dy(i,j) = lcl1 +xp1*lcl2;
			dx(i,j) = lcl2;
		}
	}

	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj_bdry(const FLT *lin1, FLT *f1, int stride) const {
	TinyMatrix<FLT,_gpn,_nmodx> wk0;
	const int bs1 = _sm+3, bs2 = 2*_sm+3, bint = _bm;
	const int lgpx = _gpx, lgpn = _gpn, lnmodx = _nmodx;  
	FLT lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(_tm), neverDeleteData);
	Array<FLT,2> f(f1, shape(_gpx,stride), neverDeleteData);
#endif

	/* DETERMINE U VALUES, GRAD U VALUES
		AT COLLOCATION POINTS
		SUM HAT(U) FOR DU/DX AND DU/DY
	*/

	/* PART I - sum u*g_mn for each n, s_j	  */
	for(int j = 0; j < lgpn; ++j ) {
		
		/* VERTEX 0			   */
		wk0(j,0) = lin(0)*_gn(j,0);

		/* VERTEX 1			   */
		lcl0 = lin(1)*_gn(j,1);
		/* SIDE 2		 */
		for(int m = bs2; m < bint; ++m )
			lcl0 += lin(m)*_gn(j,m);
		wk0(j,1) = lcl0;

		/* VERTEX 2			   */
		lcl0 = lin(2)*_gn(j,2);
		/* SIDE 1		 */
		for (int m = bs1; m < bs2; ++m)
			lcl0 += lin(m)*_gn(j,m);
		wk0(j,2) = lcl0;

		/* SIDE 0		 */
		for(int m = 3; m < bs1; ++m) {
			wk0(j,m) = lin(m)*_gn(j,m);
		}
	}
		 
	/* SUM OVER N AT EACH I,J POINT	   */	   
	 for (int i=0; i < lgpx; ++i ) {
		 for (int j=0; j < lgpn; ++j) {
			 lcl0 = 0.0;

			for(int n=0; n < lnmodx; ++n )	
				lcl0 += wk0(j,n)*_gx(i,n);
				
			f(i,j) = lcl0;
		}
	}

	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::derivr(const FLT *f1, FLT *dx1, int stride) const {
	TinyVector<FLT,_gpx> wk0;
	const int lgpn = _gpn, lgpx = _gpx;
	FLT ln0, lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,2> f((FLT *) f1, shape(_gpx,stride), neverDeleteData);
	Array<FLT,2> dx(dx1, shape(_gpx,stride), neverDeleteData);
#endif

	for (int j=0;j<lgpn;++j) {
		ln0 = _n0(j);
		for (int i=0;i<lgpx;++i)
			wk0(i) = f(i,j)*_dltx(i)*ln0;
		for (int i=0;i<lgpx;++i) {
			lcl0 = dx(i,j);
			for(int n=0;n<i;++n) 
				lcl0 += (wk0(i) + wk0(n))*_dltx1(i,n);
			for(int n=i+1;n<lgpx;++n) 
				lcl0 += (wk0(i) + wk0(n))*_dltx1(i,n);
			dx(i,j) = lcl0;
		}
	}

	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::derivs(const FLT *f1, FLT *dy1, int stride) const {
	TinyVector<FLT,_gpx> wk0;
	const int lgpn = _gpn, lgpx = _gpx;
	FLT ln0, lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,2> f((FLT *) f1, shape(_gpx,stride), neverDeleteData);
	Array<FLT,2> dy(dy1, shape(_gpx,stride), neverDeleteData);
#endif
	
	for (int j=0;j<lgpn;++j) {
		ln0 = _n0(j);
		for (int i=0;i<lgpx;++i)
			wk0(i) = f(i,j)*_dltx(i)*ln0;
		for (int i=0;i<lgpx;++i) {
			lcl0 = dy(i,j);
			for(int n=0;n<i;++n) 
				lcl0 += (wk0(i) + wk0(n))*_dltn1(i,n);
			for(int n=i+1;n<lgpx;++n) 
				lcl0 += (wk0(i) + wk0(n))*_dltn1(i,n);
			dy(i,j) = lcl0;
		}
	}

	for (int i=0;i<lgpx;++i) {
		for (int j=0;j<lgpn;++j)
			wk0(j) = f(i,j)*_dltn(j);
		for (int j=0;j<_gpn;++j) {
			lcl0 = dy(i,j);
			for(int n=0;n<j;++n) 
				lcl0 += (wk0(j) + wk0(n))*_dltn2(j,n);
			for(int n=j+1;n<lgpn;++n) 
				lcl0 += (wk0(j) + wk0(n))*_dltn2(j,n);
			dy(i,j) = lcl0;
		}
	}

	return;
}
	
template<int _p, int ep> void tri_basis<_p,ep>::proj_leg(const FLT *lin1, FLT *f1, int stride) const {
	const int lsm = _sm, ltm = _tm;
	FLT lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(_tm), neverDeleteData);
	Array<FLT,2> f(f1, shape(_gpx,stride), neverDeleteData);
#endif
	
	/* INTERIOR */
	for(int i=1;i<lsm;++i) {
		for(int j=1;j<lsm-(i-1);++j) {
			lcl0 = 0.0;
			for(int k=0;k<ltm;++k)
				lcl0 += lin(k)*_lgrnge(k,i,j);
			f(i,j) = lcl0;
		}
	}
	
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj_leg(FLT u1, FLT u2, FLT u3, FLT *f1, int stride) const {
	int lsm = _sm;
#ifdef BZ_DEBUG
	Array<FLT,2> f(f1, shape(_gpx,stride), neverDeleteData);
#endif
	
	/* INTERIOR */
	for(int i=1;i<lsm;++i)
		for(int j=1;j<lsm-(i-1);++j)
			f(i,j) = u1*_lgrnge(0,i,j) +u2*_lgrnge(1,i,j) +u3*_lgrnge(2,i,j);
	
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj_bdry_leg(const FLT *lin1, FLT *f1, int stride) const {
	int lsm = _sm, lbm = _bm;
	FLT lcl0;
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(_tm), neverDeleteData);
	Array<FLT,2> f(f1, shape(_gpx,stride), neverDeleteData);
#endif
	
	/* INTERIOR */
	for(int i=1;i<lsm;++i) {
		for(int j=1;j<lsm-(i-1);++j) {
			lcl0 = 0.0;
			for(int k=0;k<lbm;++k)
				lcl0 += lin(k)*_lgrnge(k,i,j);
			f(i,j) = lcl0;
		}
	}
	
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::proj_side(int side, const FLT *lin1, FLT *f1, FLT *dtan1, FLT *dnrm1) const {
	TinyVector<FLT,_sm+2> wk0;
	const int lgpx = _gpx, lsm = _sm, ltm = _tm;
#ifdef BZ_DEBUG
	const Array<FLT,1> lin((FLT *) lin1, shape(_tm), neverDeleteData);
	Array<FLT,1> dtan(dtan1, shape(_gpx), neverDeleteData);
	Array<FLT,1> dnrm(dnrm1, shape(_gpx), neverDeleteData);
#endif

	/* NORMAL DERIVATIVE */
	for(int i=0;i<lgpx;++i) {
		dnrm(i) = 0.0;
		for(int m=0;m<ltm;++m)
			dnrm(i) += lin(m)*_dgnorm(side,m,i);
	}

/*	MOVE SIDE & VERTEX MODES TO WK & THEN PROJECT VALUES & TANGENT DERIVS */
	wk0(0) = lin((side+1)%3);
	wk0(1) = lin((side+2)%3);

	int base = side*lsm;
	for(int m=0;m<lsm;++m)
		wk0(m+2) = lin(base +3 +m);
		
	proj1d(&wk0(0),f1,dtan1);
		
	return;
}

#ifndef BZ_DEBUG
#undef lin
#undef f
#undef dx
#undef dy
#undef dnrm
#undef dtan
#endif
