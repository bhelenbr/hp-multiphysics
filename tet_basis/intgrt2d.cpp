/*
 *  intgrt2d.cpp
 *  tet_basis
 *
 *  Created by Michael Brazell on 5/8/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_basis.h"
#include <stdio.h>
#include <stdlib.h>

/* 2D to 1D access */
#ifndef BZ_DEBUG
#define f(i,j) f1[(i)*stride +j]
#define dx(i,j) dx1[(i)*stride +j]
#define dy(i,j) dy1[(i)*stride +j]
#define rslt(i) rslt1[i]
#endif

/* Integrates Face 0, Edge 1,2,3, Vertex 1,2,3*/
void tet_basis::intgrt2d(FLT *rslt1, FLT *f1, int stride) {
   Array<FLT,2> wk0(nmodx,gpy);
   const int lgpy=gpy,lgpx=gpx,lnmodx=nmodx;
   FLT lcl0;
   int sign;
#ifdef BZ_DEBUG
   Array<FLT,2> f(f1, shape(gpx,stride), neverDeleteData);
   Array<FLT,1> rslt(rslt1, shape(3+3*em+fm), neverDeleteData);
#endif
    
	/* INTEGRATE IN X-DIRECTION */
	
	for (int n = 0; n < lnmodx; ++n) {
		for (int j = 0; j < lgpy; ++j) {
			lcl0 = f(0,j)*gxwtx(n,0); 				
			for(int i = 1; i < lgpx; ++i) {
				lcl0 += f(i,j)*gxwtx(n,i);
			}
			wk0(n,j) = lcl0;
		}
	}  	
  
	/* INTEGRATE IN Y-DIRECTION */	
	 
	/* VERTEX 1*/	
	lcl0 = wk0(0,0)*gywty(1,0);
	for (int j = 1; j < lgpy; ++j) {
		lcl0 += wk0(0,j)*gywty(1,j);
	}
	rslt(0)=lcl0;		
	
	/* VERTEX 2*/	
	lcl0 = wk0(1,0)*gywty(2,0);
	for (int j = 1; j < lgpy; ++j) {
		lcl0 += wk0(1,j)*gywty(2,j);
	}
	rslt(1)=lcl0;
	
	/* VERTEX 3*/	
	lcl0 = wk0(2,0)*gywty(2,0);
	for (int j = 1; j < lgpy; ++j) {
		lcl0 += wk0(2,j)*gywty(2,j);
	}
	rslt(2)=lcl0;		
	
	/* EDGE 1*/
	int ind = 2;
	for (int n = 3; n < em+3; ++n){
		lcl0 = wk0(n,0)*gywty(n,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(n,j)*gywty(n,j);
		}
		rslt(++ind)=lcl0;
	}
	
	/* EDGE 2*/
	for (int n = 0; n < em; ++n){
		lcl0 = wk0(2,0)*gywty(n+em+3,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(2,j)*gywty(n+em+3,j);
		}
		rslt(++ind)=lcl0;
	}
	
	/* EDGE 3*/
	sign = 1;
	for (int n = 0; n < em; ++n){
		lcl0 = wk0(1,0)*gywty(n+em+3,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(1,j)*gywty(n+em+3,j);
		}
		rslt(++ind)=sign*lcl0;
		sign*=-1;
	}
	
	/* FACE 0*/   
	int ind2 = 0;
	for (int p = 1; p <= em; ++p){
		for(int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = wk0(p+2,0)*gywty(2+2*em+ind2,0);
			for (int j = 1; j < lgpy; ++j){
				lcl0 += wk0(p+2,j)*gywty(2+2*em+ind2,j);
			}
			rslt(++ind)=lcl0;	
		}
	}			

   return;
}

