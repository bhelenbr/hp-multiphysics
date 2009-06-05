/*
 *  proj2d.cpp
 *  tet_basis
 *
 *  Created by Michael Brazell on 5/14/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_basis.h"
#include <stdio.h>

#ifndef BZ_DEBUG
#define lin(i) lin1[i]
#define f(i,j) f1[(i)*stride +j]
#define dx(i,j) dx1[(i)*stride +j]
#define dy(i,j) dy1[(i)*stride +j]
#endif

void tet_basis::proj2d(FLT *lin1, FLT *f1, FLT *dx1, FLT *dy1, int stride) {
	Array<FLT,2> wk0(gpy,3+em);
	Array<FLT,2> wk1(gpy,3+em);
	Array<FLT,2> wk2(gpy,3+em);
   const int be2 = em+3, be3 = 2*em+3, bint = 3+3*em;
   const int lgpx = gpx, lgpy = gpy, lnmodx = nmodx;  
   FLT lcl0, lcl1, lcl2;
   FLT xp1,oeta; 
   int sign;
#ifdef BZ_DEBUG
   Array<FLT,1> lin(lin1, shape(3+3*em+fm), neverDeleteData);
   Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
   Array<FLT,2> dx(dx1, shape(gpx,stride), neverDeleteData);
   Array<FLT,2> dy(dy1, shape(gpx,stride), neverDeleteData);
#endif
   
   /* DETERMINE U VALUES, GRAD U VALUES
      AT COLLOCATION POINTS
      SUM HAT(U) FOR DU/DX AND DU/DY
   */
   
   /* GENERAL FORMULA
   	dg/dr = 2.0/(1-n) g(s) dg/dx
   	dg/ds = g(x)dg/dn +(1+x)/2 dg/dr = g(x)dg/dn +(1+x)/(1-s) g(s) dg/dx 
   */

   /* PART I - sum u*g_mn for each n, s_j   */
	for(int j = 0; j < lgpy; ++j ) {
		oeta = y0(j);
      
		/* VERTEX 1 */
		wk0(j,0) = lin(0)*gy(j,1);
		wk1(j,0) = lin(0)*dgy(j,1);
		wk2(j,0) = wk0(j,0)*oeta;

		/* VERTEX 2, EDGE 3 */
		sign = 1;
		lcl0 = lin(1)*gy(j,2);
		lcl1 = lin(1)*dgy(j,2);
		for(int i = 0; i < em; ++i){
			lcl0 += sign*lin(be3+i)*gy(j,3+em+i);	
			lcl1 += sign*lin(be3+i)*dgy(j,3+em+i);
			sign*=-1;
		}
		wk0(j,1) = lcl0;
		wk1(j,1) = lcl1;
		wk2(j,1) = wk0(j,1)*oeta;

		/* VERTEX 3, EDGE 2 */
		lcl0 = lin(2)*gy(j,2);
		lcl1 = lin(2)*dgy(j,2);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be2+i)*gy(j,3+em+i);	
			lcl1 += lin(be2+i)*dgy(j,3+em+i);
		}
		wk0(j,2) = lcl0;
		wk1(j,2) = lcl1;
		wk2(j,2) = wk0(j,2)*oeta;
		
		/* EDGE 1, FACE 0 */		
		int ind2 = 0;		
		for(int p = 3; p < em+3; ++p){
			lcl0=lin(p)*gy(j,p);	
			lcl1=lin(p)*gy(j,p);		
			for(int i = 1; i <= em-p+2; ++i){	
				lcl0 += lin(bint+ind2)*gy(j,3+2*em+ind2);
				lcl1 += lin(bint+ind2)*gy(j,3+2*em+ind2);
				++ind2;
			}			
			wk0(j,p) = lcl0;
			wk1(j,p) = lcl1;
			wk2(j,p) = wk0(j,p)*oeta;
		}
	}
		
		
	/* SUM OVER N AT EACH I,J POINT   */  
	for (int i = 0; i < lgpx; ++i ) {
		xp1 = x0(i);
		for (int j = 0; j < lgpy; ++j) {
			lcl0 = wk0(j,0)*gx(i,0);
			lcl1 = wk1(j,0)*gx(i,0);
			lcl2 = wk2(j,0)*dgx(i,0);
			
			for(int n = 1; n < lnmodx; ++n ) {     
				lcl0 += wk0(j,n)*gx(i,n);
				lcl1 += wk1(j,n)*gx(i,n);
				lcl2 += wk2(j,n)*dgx(i,n);
			}
			f(i,j)  = lcl0;
			dy(i,j) = lcl1 +xp1*lcl2;
			dx(i,j) = lcl2;
		}
	}

   return;
}

/* added routine */
void tet_basis::proj2d_bdry(FLT *lin1, FLT *f1, FLT *dx1, FLT *dy1, int stride) {
	Array<FLT,2> wk0(gpy,3+em);
	Array<FLT,2> wk1(gpy,3+em);
	Array<FLT,2> wk2(gpy,3+em);
	const int be2 = em+3, be3 = 2*em+3, bint = 3+3*em;
	const int lgpx = gpx, lgpy = gpy, lnmodx = nmodx;  
	FLT lcl0, lcl1, lcl2;
	FLT xp1,oeta; 
	int sign;
#ifdef BZ_DEBUG
	Array<FLT,1> lin(lin1, shape(3+3*em+fm), neverDeleteData);
	Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
	Array<FLT,2> dx(dx1, shape(gpx,stride), neverDeleteData);
	Array<FLT,2> dy(dy1, shape(gpx,stride), neverDeleteData);
#endif
	
	/* DETERMINE U VALUES, GRAD U VALUES
	 AT COLLOCATION POINTS
	 SUM HAT(U) FOR DU/DX AND DU/DY
	 */
	
	/* GENERAL FORMULA
	 dg/dr = 2.0/(1-n) g(s) dg/dx
	 dg/ds = g(x)dg/dn +(1+x)/2 dg/dr = g(x)dg/dn +(1+x)/(1-s) g(s) dg/dx 
	 */
	
	/* PART I - sum u*g_mn for each n, s_j   */
	for(int j = 0; j < lgpy; ++j ) {
		oeta = y0(j);
		
		/* VERTEX 1 */
		wk0(j,0) = lin(0)*gy(j,1);
		wk1(j,0) = lin(0)*dgy(j,1);
		wk2(j,0) = wk0(j,0)*oeta;
		
		/* VERTEX 2, EDGE 3 */
		sign = 1;
		lcl0 = lin(1)*gy(j,2);
		lcl1 = lin(1)*dgy(j,2);
		for(int i = 0; i < em; ++i){
			lcl0 += sign*lin(be3+i)*gy(j,3+em+i);	
			lcl1 += sign*lin(be3+i)*dgy(j,3+em+i);
			sign*=-1;
		}
		wk0(j,1) = lcl0;
		wk1(j,1) = lcl1;
		wk2(j,1) = wk0(j,1)*oeta;
		
		/* VERTEX 3, EDGE 2 */
		lcl0 = lin(2)*gy(j,2);
		lcl1 = lin(2)*dgy(j,2);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be2+i)*gy(j,3+em+i);	
			lcl1 += lin(be2+i)*dgy(j,3+em+i);
		}
		wk0(j,2) = lcl0;
		wk1(j,2) = lcl1;
		wk2(j,2) = wk0(j,2)*oeta;
		
		/* EDGE 1, FACE 0 */		
		int ind2 = 0;		
		for(int p = 3; p < em+3; ++p){
			lcl0=lin(p)*gy(j,p);	
			lcl1=lin(p)*gy(j,p);		
			for(int i = 1; i <= em-p+2; ++i){	
				lcl0 += lin(bint+ind2)*gy(j,3+2*em+ind2);
				lcl1 += lin(bint+ind2)*gy(j,3+2*em+ind2);
				++ind2;
			}			
			wk0(j,p) = lcl0;
			wk1(j,p) = lcl1;
			wk2(j,p) = wk0(j,p)*oeta;
		}
	}
	
	
	/* SUM OVER N AT EACH I,J POINT   */  
	for (int i = 0; i < lgpx; ++i ) {
		xp1 = x0(i);
		for (int j = 0; j < lgpy; ++j) {
			lcl0 = wk0(j,0)*gx(i,0);
			lcl1 = wk1(j,0)*gx(i,0);
			lcl2 = wk2(j,0)*dgx(i,0);
			
			for(int n = 1; n < lnmodx; ++n ) {     
				lcl0 += wk0(j,n)*gx(i,n);
				lcl1 += wk1(j,n)*gx(i,n);
				lcl2 += wk2(j,n)*dgx(i,n);
			}
			f(i,j)  = lcl0;
			dy(i,j) = lcl1 +xp1*lcl2;
			dx(i,j) = lcl2;
		}
	}
	
	return;
}

void tet_basis::proj2d(FLT *lin1, FLT *f1, int stride) {
	Array<FLT,2> wk0(gpy,3+em);
   const int be2 = em+3, be3 = 2*em+3, bint = 3+3*em;
   const int lgpx = gpx, lgpy = gpy, lnmodx = nmodx;  
   int sign;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,1> lin(lin1, shape(3+3*em+fm), neverDeleteData);
   Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
#endif

   /* DETERMINE U VALUES, GRAD U VALUES
      AT COLLOCATION POINTS
      SUM HAT(U) FOR DU/DX AND DU/DY
   */

   /* PART I - sum u*g_mn for each n, s_j   */
  	for(int j = 0; j < lgpy; ++j ) {
     
		/* VERTEX 1 */
		wk0(j,0) = lin(0)*gy(j,1);

		/* VERTEX 2, EDGE 3 */
		sign = 1;
		lcl0 = lin(1)*gy(j,2);
		for(int i = 0; i < em; ++i){
			lcl0 += sign*lin(be3+i)*gy(j,3+em+i);	
			sign*=-1;
		}
		wk0(j,1) = lcl0;

		/* VERTEX 3, EDGE 2 */
		lcl0 = lin(2)*gy(j,2);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be2+i)*gy(j,3+em+i);
		}
		wk0(j,2) = lcl0;
		
		/* EDGE 1, FACE 0 */		
		int ind2 = 0;	
		int bint1 = em-3+2;	
		for(int p = 3; p < em+3; ++p){
			lcl0=lin(p)*gy(j,p);	
			for(int i = 1; i <= em-p+2; ++i){	
				lcl0 += lin(bint+ind2)*gy(j,3+2*em+ind2);
				//cout << bint+i+ind2<< endl;
				++ind2;
			}						
			wk0(j,p) = lcl0;
	
		}
	}		
	
		
	/* SUM OVER N AT EACH I,J POINT   */  
	for (int i = 0; i < lgpx; ++i ) {
		for (int j = 0; j < lgpy; ++j) {
			lcl0 = wk0(j,0)*gx(i,0);
			for(int n = 1; n < lnmodx; ++n ) {     
				lcl0 += wk0(j,n)*gx(i,n);
			}
			f(i,j)  = lcl0;
		}
	}
	
   return;
}

void tet_basis::proj2d_bdry(FLT *lin1, FLT *f1, int stride) {
	Array<FLT,2> wk0(gpy,3+em);
   const int be2 = em+3, be3 = 2*em+3;
   const int lgpx = gpx, lgpy = gpy, lnmodx = nmodx;  
   int sign;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,1> lin(lin1, shape(3+3*em+fm), neverDeleteData);
   Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
#endif

   /* DETERMINE U VALUES, GRAD U VALUES
      AT COLLOCATION POINTS
      SUM HAT(U) FOR DU/DX AND DU/DY
   */

   /* PART I - sum u*g_mn for each n, s_j   */
  	for(int j = 0; j < lgpy; ++j ) {
     
		/* VERTEX 1 */
		wk0(j,0) = lin(0)*gy(j,1);

		/* VERTEX 2, EDGE 3 */
		sign = 1;
		lcl0 = lin(1)*gy(j,2);
		for(int i = 0; i < em; ++i){
			lcl0 += sign*lin(be3+i)*gy(j,3+em+i);	
			sign*=-1;
		}
		wk0(j,1) = lcl0;

		/* VERTEX 3, EDGE 2 */
		lcl0 = lin(2)*gy(j,2);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be2+i)*gy(j,3+em+i);
		}
		wk0(j,2) = lcl0;
		
		/* EDGE 1 */		
		for(int p = 3; p < em+3; ++p){
			wk0(j,p)=lin(p)*gy(j,p);	
	
		}
	}		
		
	/* SUM OVER N AT EACH I,J POINT   */  
	for (int i = 0; i < lgpx; ++i ) {
		for (int j = 0; j < lgpy; ++j) {
			lcl0 = wk0(j,0)*gx(i,0);
			for(int n = 1; n < lnmodx; ++n ) {     
				lcl0 += wk0(j,n)*gx(i,n);
			}
			f(i,j)  = lcl0;
		}
	}
	
   return;
}


void tet_basis::proj2d(FLT u1, FLT u2, FLT u3, FLT *f1, int stride) {
   const int lgpx = gpx, lgpy = gpy;
#ifdef BZ_DEBUG
   Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
#endif

	for (int i=0; i < lgpx; ++i ){ 
		for (int j=0; j < lgpy; ++j ){ 
			f(i,j) = u1*gx(i,0)*gy(j,1) +u2*gx(i,1)*gy(j,2) +u3*gx(i,2)*gy(j,2);
		}
	}
	return;
}


void tet_basis::proj2d_leg(FLT *lin1, FLT *f1, int stride) {
   const int ltm = 3+3*em+fm;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,1> lin(lin1, shape(3+3*em+fm), neverDeleteData);
   Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
#endif
   
   /* INTERIOR */
	for(int i = 1; i < em; ++i) {
		for(int j = 1; j < em+1-i; ++j) {	
			lcl0 = 0.0;
			for(int m = 0; m < ltm; ++m)
				lcl0 += lin(m)*lgrnge2d(m,i,j);
			f(i,j) = lcl0;
		}
	}
   
   return;
}

void tet_basis::proj2d_bdry_leg(FLT *lin1, FLT *f1, int stride) {
   const int ltm = 3+3*em;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,1> lin(lin1, shape(3+3*em+fm), neverDeleteData);
   Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
#endif
   
   /* INTERIOR */
	for(int i = 1; i < em; ++i) {
		for(int j = 1; j < em+1-i; ++j) {	
			lcl0 = 0.0;
			for(int m = 0; m < ltm; ++m)
				lcl0 += lin(m)*lgrnge2d(m,i,j);
			f(i,j) = lcl0;
		}
	}
   
   return;
}


void tet_basis::proj2d_leg(FLT u1, FLT u2, FLT u3, FLT *f1, int stride) {
#ifdef BZ_DEBUG
   Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
#endif
   
   /* INTERIOR */
	for(int i = 1; i < em; ++i)
		for(int j = 1; j < em+1-i; ++j)
			f(i,j) = u1*lgrnge2d(0,i,j) +u2*lgrnge2d(1,i,j) +u3*lgrnge2d(2,i,j);
   
   return;
}




