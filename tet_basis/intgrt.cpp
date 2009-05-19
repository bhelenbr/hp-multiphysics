/*
 *  intgrt3d.cpp
 *  tet_basis
 *
 *  Created by Michael Brazell on 10/19/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_basis.h"
#include <stdio.h>
#include <stdlib.h>

/* 3D to 1D access */
#ifndef BZ_DEBUG
#define f(i,j,k) f1[(i)*stridex +(j)*stridey +k]
#define dx(i,j,k) dx1[(i)*stridex +(j)*stridey +k]
#define dy(i,j,k) dy1[(i)*stridex +(j)*stridey +k]
#define dz(i,j,k) dz1[(i)*stridex +(j)*stridey +k]
#define rslt(i) rslt1[i]
#endif

void tet_basis::intgrt(FLT *rslt1, FLT *f1, int stridex, int stridey) {
   Array<FLT,3> wk0(nmodx,gpz,gpy);
   Array<FLT,2> wk1(4+3*em+fm,gpz);
   const int lgpz=gpz,lgpy=gpy,lgpx=gpx,lnmodx=nmodx;
   int sign;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
   Array<FLT,1> rslt(rslt1, shape(tm), neverDeleteData);
#endif
    
	/* INTEGRATE IN X-DIRECTION */
	
	for (int n = 0; n < lnmodx; ++n) {
		for (int k = 0; k < lgpz; ++k){	  
			for (int j = 0; j < lgpy; ++j) {
				lcl0 = f(0,j,k)*gxwtx(n,0); 				
				for(int i = 1; i < lgpx; ++i) {
					lcl0 += f(i,j,k)*gxwtx(n,i);
				}
				wk0(n,k,j) = lcl0;
			}
		}
	}  
	
  
	/* INTEGRATE IN Y-DIRECTION */
	
	/* VERTEX 0 */
	for (int k = 0; k < lgpz; ++k){
		lcl0 = wk0(0,k,0)*gywty(0,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(0,k,j)*gywty(0,j);
		}
		wk1(0,k)=lcl0;
	}
	 
	/* VERTEX 1, EDGE 4*/	
	for (int k = 0; k < lgpz; ++k){
		lcl0 = wk0(0,k,0)*gywty(1,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(0,k,j)*gywty(1,j);
		}
		wk1(1,k)=lcl0;
	}
	
	
	/* VERTEX 2, EDGE 5 */	
	for (int k = 0; k < lgpz; ++k){
		lcl0 = wk0(1,k,0)*gywty(2,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(1,k,j)*gywty(2,j);
		}
		wk1(2,k)=lcl0;
	}
	
	/* VERTEX 3, EDGE 6 */	
	for (int k = 0; k < lgpz; ++k){
		lcl0 = wk0(2,k,0)*gywty(2,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(2,k,j)*gywty(2,j);
		}
		wk1(3,k)=lcl0;
	}		
	
	/* EDGE 1, FACE 1 */
	int ind = 4;
	for (int n = 0; n < em; ++n){
		for (int k = 0; k < lgpz; ++k){
			lcl0 = wk0(n+3,k,0)*gywty(n+3,0);
			for (int j = 1; j < lgpy; ++j) {
				lcl0 += wk0(n+3,k,j)*gywty(n+3,j);
			}
			wk1(ind,k)=lcl0;
		}
		++ind;
	}
	
	/* EDGE 2, FACE 2 */
	for (int n = 0; n < em; ++n){
		for (int k = 0; k < lgpz; ++k){
			lcl0 = wk0(2,k,0)*gywty(n+em+3,0);
			for (int j = 1; j < lgpy; ++j) {
				lcl0 += wk0(2,k,j)*gywty(n+em+3,j);
			}
			wk1(ind,k)=lcl0;
		}
		++ind;
	}
	
	/* EDGE 3, FACE 3 */
	for (int n = 0; n < em; ++n){
		for (int k = 0; k < lgpz; ++k){
			lcl0 = wk0(1,k,0)*gywty(n+em+3,0);
			for (int j = 1; j < lgpy; ++j) {
				lcl0 += wk0(1,k,j)*gywty(n+em+3,j);
			}
			wk1(ind,k)=lcl0;
		}
		++ind;
	}
	
	/* FACE 0, INTERIOR */   
	int ind2 = 0;
	for (int p = 1; p <= em; ++p){
		for(int q = 1; q <= em-p; ++q){
			++ind2;
			for (int k = 0; k < lgpz; ++k){
				lcl0 = wk0(p+2,k,0)*gywty(2+2*em+ind2,0);
				for (int j = 1; j < lgpy; ++j) {
					lcl0 += wk0(p+2,k,j)*gywty(2+2*em+ind2,j);
				}
				wk1(ind,k)=lcl0;
			}
			++ind;		
		}
	}	
		
	/* INTEGRATE IN Z-DIRECTION */
	
	/* VERTEX 0 */	
	lcl0 = 0.0;
	for (int k=0; k < lgpz; ++k) 
         lcl0 += gzwtz(0,k)*wk1(0,k);
	rslt(0)=lcl0;
	
	/* VERTEX 1 */	
	lcl0 = 0.0;
	for (int k=0; k < lgpz; ++k) 
         lcl0 += gzwtz(1,k)*wk1(1,k);
	rslt(1)=lcl0;
	
	/* VERTEX 2 */	
	lcl0 = 0.0;
	for (int k=0; k < lgpz; ++k) 
         lcl0 += gzwtz(1,k)*wk1(2,k);
	rslt(2)=lcl0;
	
	/* VERTEX 3 */	
	lcl0 = 0.0;
	for (int k=0; k < lgpz; ++k) 
         lcl0 += gzwtz(1,k)*wk1(3,k);
	rslt(3)=lcl0;
	
	/* EDGE 1 */
	ind = 4;
	for (int n = 2; n < em+2;++n){
		lcl0 = 0.0;
		for (int k=0; k < lgpz; ++k) 
			lcl0 += gzwtz(n,k)*wk1(n+2,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 2 */
	for (int n = 2; n < em+2;++n){
		lcl0 = 0.0;
		for (int k=0; k < lgpz; ++k) 
			lcl0 += gzwtz(n,k)*wk1(n+em+2,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 3 */
	for (int n = 2; n < em+2;++n){
		lcl0 = 0.0;
		for (int k=0; k < lgpz; ++k) 
			lcl0 += gzwtz(n,k)*wk1(n+2*em+2,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 4 */
	for (int n = 0; n < em;++n){
		lcl0 = 0.0;
		for (int k=0; k < lgpz; ++k) 
			lcl0 += gzwtz(n+em+2,k)*wk1(1,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 5 */
	for (int n = 0; n < em;++n){
		lcl0 = 0.0;
		for (int k=0; k < lgpz; ++k) 
			lcl0 += gzwtz(n+em+2,k)*wk1(2,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 6 */
	for (int n = 0; n < em;++n){
		lcl0 = 0.0;
		for (int k=0; k < lgpz; ++k) 
			lcl0 += gzwtz(n+em+2,k)*wk1(3,k);
		rslt(ind++) = lcl0;
	}
	
	/* FACE 0 */  
	ind2 = 0;
	for (int p = 1; p <= em-1; ++p){
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = 0.0;
			for (int k=0; k < lgpz; ++k) 
				lcl0 += gzwtz(p+q+1,k)*wk1(3*em+3+ind2,k);
			rslt(ind++) = lcl0;
		}
	}	

	/* FACE 1 */
	ind2 = 0;
	sign = 1;
	for (int p = 1; p <= em-1;++p){
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = 0.0;
			for (int k=0; k < lgpz; ++k) 
				lcl0 += gzwtz(1+2*em+ind2,k)*wk1(3+p,k);
			rslt(ind++) = sign*lcl0;
		}
		sign*=-1;
	}

	/* FACE 2 */
	ind2 = 0;
	sign = 1;
	for (int p = 1; p <= em-1;++p){
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = 0.0;
			for (int k=0; k < lgpz; ++k) 
				lcl0 += gzwtz(1+2*em+ind2,k)*wk1(em+3+p,k);
			rslt(ind++) = sign*lcl0;
		}
		sign*=-1;
	}
	
	/* FACE 3 */
	ind2 = 0;
	for (int p = 1; p <= em-1;++p){
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = 0.0;
			for (int k=0; k < lgpz; ++k) 
				lcl0 += gzwtz(1+2*em+ind2,k)*wk1(2*em+3+p,k);
			rslt(ind++) = lcl0;
		}
	}

	

	/* INTERIOR */	
	ind2 = 0;
	int ind3 = 0;
	for (int p = 1; p <= em-1; ++p) {
		for (int q = 1; q <= em-p; ++q) {
			++ind2;
			for (int r = 1; r <= em-p-q; ++r) {
				++ind3;
				lcl0 = 0.0;
				for (int k=0; k < lgpz; ++k) 
					lcl0 += gzwtz(1+2*em+fm+ind3,k)*wk1(3*em+3+ind2,k);
				rslt(ind++) = lcl0;
			}
		}
	}
	

   return;
}


///* WARNING THIS ADDS INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
void tet_basis::intgrtrst(FLT *rslt1, FLT *dx1, FLT *dy1, FLT *dz1, int stridex, int stridey) {
	Array<FLT,3> wk0(nmodx,gpz,gpy);
	Array<FLT,3> wk1(nmodx,gpz,gpy);
	Array<FLT,3> wk2(nmodx,gpz,gpy);
	Array<FLT,3> wk3(gpx,gpy,gpz);
	Array<FLT,3> wk4(gpx,gpy,gpz);
	Array<FLT,2> wk5(4+3*em+fm,gpz);
	Array<FLT,2> wk6(4+3*em+fm,gpz);
	const int lgpz=gpz,lgpy=gpy,lgpx=gpx,lnmodx=nmodx;
	int sign;
	FLT lcl0, lcl1, lcl2;
#ifdef BZ_DEBUG
	Array<FLT,3> dx(dx1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,3> dy(dy1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,3> dz(dz1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,1> rslt(rslt1, shape(tm), neverDeleteData);
#endif


	/* INTEGRATE IN X-DIRECTION */

	for (int k = 0; k < lgpz; ++k){	  
		for (int j = 0; j < lgpy; ++j) {
			wk3(0,j,k) = y0(j)*z0(k)*(dx(0,j,k)+x0(0)*(dy(0,j,k)+dz(0,j,k)));
			wk4(0,j,k) = z0(k)*(dy(0,j,k)+y1(j)*dz(0,j,k));
			lcl0 = wk3(0,j,k)*dgxwtx(0,0);
			lcl1 = wk4(0,j,k)*gxwtx(0,0);
			lcl2 = dz(0,j,k)*gxwtx(0,0); 				
			for(int i = 1; i < lgpx; ++i) {
				wk3(i,j,k) = y0(j)*z0(k)*(dx(i,j,k)+x0(i)*(dy(i,j,k)+dz(i,j,k)));
				wk4(i,j,k) = z0(k)*(dy(i,j,k)+y1(j)*dz(i,j,k));
				lcl0 += wk3(i,j,k)*dgxwtx(0,i);
				lcl1 += wk4(i,j,k)*gxwtx(0,i);
				lcl2 += dz(i,j,k)*gxwtx(0,i); 
			}
			wk0(0,k,j) = lcl0;
			wk1(0,k,j) = lcl1;
			wk2(0,k,j) = lcl2;
		}
	}
	
	for (int n = 1; n < lnmodx; ++n) {
		for (int k = 0; k < lgpz; ++k){	  
			for (int j = 0; j < lgpy; ++j) {
			lcl0 = wk3(0,j,k)*dgxwtx(n,0);
			lcl1 = wk4(0,j,k)*gxwtx(n,0);
			lcl2 = dz(0,j,k)*gxwtx(n,0); 				
			for(int i = 1; i < lgpx; ++i) {
				lcl0 += wk3(i,j,k)*dgxwtx(n,i);
				lcl1 += wk4(i,j,k)*gxwtx(n,i);
				lcl2 += dz(i,j,k)*gxwtx(n,i); 
			}
			wk0(n,k,j) = lcl0;
			wk1(n,k,j) = lcl1;
			wk2(n,k,j) = lcl2;
			}
		}
	}  
   
   
   	/* INTEGRATE IN Y-DIRECTION */
	
	/* VERTEX 0 */
	for (int k = 0; k < lgpz; ++k){
		lcl0 = wk0(0,k,0)*gywty(0,0)+wk1(0,k,0)*dgywty(0,0);
		lcl1 = wk2(0,k,0)*gywty(0,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(0,k,j)*gywty(0,j)+wk1(0,k,j)*dgywty(0,j);
			lcl1 += wk2(0,k,j)*gywty(0,j);
		}
		wk5(0,k)=lcl0;
		wk6(0,k)=lcl1;
	}
	 
	/* VERTEX 1, EDGE 4 */	
	for (int k = 0; k < lgpz; ++k){
		lcl0 = wk0(0,k,0)*gywty(1,0) + wk1(0,k,0)*dgywty(1,0);
		lcl1 = wk2(0,k,0)*gywty(1,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(0,k,j)*gywty(1,j)+ wk1(0,k,j)*dgywty(1,j);
			lcl1 += wk2(0,k,j)*gywty(1,j);
		}
		wk5(1,k)=lcl0;
		wk6(1,k)=lcl1;
	}
	
	
	/* VERTEX 2, EDGE 5 */	
	for (int k = 0; k < lgpz; ++k){
		lcl0 = wk0(1,k,0)*gywty(2,0) + wk1(1,k,0)*dgywty(2,0);
		lcl1 = wk2(1,k,0)*gywty(2,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(1,k,j)*gywty(2,j) + wk1(1,k,j)*dgywty(2,j);
			lcl1 += wk2(1,k,j)*gywty(2,j);
		}
		wk5(2,k)=lcl0;
		wk6(2,k)=lcl1;
	}
	
	/* VERTEX 3, EDGE 6 */	
	for (int k = 0; k < lgpz; ++k){
		lcl0 = wk0(2,k,0)*gywty(2,0) + wk1(2,k,0)*dgywty(2,0);
		lcl1 = wk2(2,k,0)*gywty(2,0);
		for (int j = 1; j < lgpy; ++j) {
			lcl0 += wk0(2,k,j)*gywty(2,j) + wk1(2,k,j)*dgywty(2,j);
			lcl1 += wk2(2,k,j)*gywty(2,j);
		}
		wk5(3,k)=lcl0;
		wk6(3,k)=lcl1;
	}	
	
	
	/* EDGE 1, FACE 1 */
	int ind = 4;
	for (int n = 0; n < em; ++n){
		for (int k = 0; k < lgpz; ++k){
			lcl0 = wk0(n+3,k,0)*gywty(n+3,0) + wk1(n+3,k,0)*dgywty(n+3,0);
			lcl1 = wk2(n+3,k,0)*gywty(n+3,0);
			for (int j = 1; j < lgpy; ++j) {
				lcl0 += wk0(n+3,k,j)*gywty(n+3,j) + wk1(n+3,k,j)*dgywty(n+3,j);
				lcl1 += wk2(n+3,k,j)*gywty(n+3,j);
			}
			wk5(ind,k)=lcl0;
			wk6(ind,k)=lcl1;
		}
		++ind;
	}
	
	/* EDGE 2, FACE 2 */
	for (int n = 0; n < em; ++n){
		for (int k = 0; k < lgpz; ++k){
			lcl0 = wk0(2,k,0)*gywty(n+em+3,0) + wk1(2,k,0)*dgywty(n+em+3,0);
			lcl1 = wk2(2,k,0)*gywty(n+em+3,0);
			for (int j = 1; j < lgpy; ++j) {
				lcl0 += wk0(2,k,j)*gywty(n+em+3,j) + wk1(2,k,j)*dgywty(n+em+3,j);
				lcl1 += wk2(2,k,j)*gywty(n+em+3,j);
			}
			wk5(ind,k)=lcl0;
			wk6(ind,k)=lcl1;
		}
		++ind;
	}
	
	/* EDGE 3, FACE 3 */
	for (int n = 0; n < em; ++n){
		for (int k = 0; k < lgpz; ++k){
			lcl0 = wk0(1,k,0)*gywty(n+em+3,0) + wk1(1,k,0)*dgywty(n+em+3,0);
			lcl1 = wk2(1,k,0)*gywty(n+em+3,0);
			for (int j = 1; j < lgpy; ++j) {
				lcl0 += wk0(1,k,j)*gywty(n+em+3,j) + wk1(1,k,j)*dgywty(n+em+3,j);
				lcl1 += wk2(1,k,j)*gywty(n+em+3,j);
			}
			wk5(ind,k)=lcl0;
			wk6(ind,k)=lcl1;
		}
		++ind;
	}

	/* FACE 0, INTERIOR */   
	int ind2 = 0;
	for (int p = 1; p <= em-1; ++p){
		for(int q = 1; q <= em - p; ++q){
			++ind2;
			for (int k = 0; k < lgpz; ++k){
				lcl0 = wk0(p+2,k,0)*gywty(2+2*em+ind2,0)+ wk1(p+2,k,0)*dgywty(2+2*em+ind2,0);
				lcl1 = wk2(p+2,k,0)*gywty(2+2*em+ind2,0);
				for (int j = 1; j < lgpy; ++j) {
					lcl0 += wk0(p+2,k,j)*gywty(2+2*em+ind2,j)+ wk1(p+2,k,j)*dgywty(2+2*em+ind2,j);
					lcl1 += wk2(p+2,k,j)*gywty(2+2*em+ind2,j);
				}
				wk5(ind,k)=lcl0;
				wk6(ind,k)=lcl1;
			}
			++ind;		
		}
	}	
	
	
	
		/* INTEGRATE IN Z-DIRECTION */
	
	/* VERTEX 0 */	
	for (int k=0; k < lgpz; ++k) 
         rslt(0) -= gzwtz(0,k)*wk5(0,k)+dgzwtz(0,k)*wk6(0,k);
	
	/* VERTEX 1 */	
	for (int k=0; k < lgpz; ++k) 
         rslt(1) -= gzwtz(1,k)*wk5(1,k)+dgzwtz(1,k)*wk6(1,k);
	
	/* VERTEX 2 */	
	for (int k=0; k < lgpz; ++k) 
        rslt(2) -= gzwtz(1,k)*wk5(2,k)+dgzwtz(1,k)*wk6(2,k);
	
	/* VERTEX 3 */	
	for (int k=0; k < lgpz; ++k) 
		rslt(3) -= gzwtz(1,k)*wk5(3,k)+dgzwtz(1,k)*wk6(3,k);
	
	/* EDGE 1 */
	ind = 4;
	for (int n = 2; n < em+2;++n){
		lcl0 = rslt(ind);
		for (int k=0; k < lgpz; ++k) 
			lcl0 -= gzwtz(n,k)*wk5(n+2,k)+dgzwtz(n,k)*wk6(n+2,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 2 */
	for (int n = 2; n < em+2;++n){
		lcl0 = rslt(ind);
		for (int k=0; k < lgpz; ++k) 
			lcl0 -= gzwtz(n,k)*wk5(n+em+2,k)+dgzwtz(n,k)*wk6(n+em+2,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 3 */
	for (int n = 2; n < em+2;++n){
		lcl0 = rslt(ind);
		for (int k=0; k < lgpz; ++k) 
			lcl0 -= gzwtz(n,k)*wk5(n+2*em+2,k)+dgzwtz(n,k)*wk6(n+2*em+2,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 4 */
	for (int n = 0; n < em;++n){
		lcl0 = rslt(ind);
		for (int k=0; k < lgpz; ++k) 
			lcl0 -= gzwtz(n+em+2,k)*wk5(1,k)+dgzwtz(n+em+2,k)*wk6(1,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 5 */
	for (int n = 0; n < em;++n){
		lcl0 = rslt(ind);
		for (int k=0; k < lgpz; ++k) 
			lcl0 -= gzwtz(n+em+2,k)*wk5(2,k)+dgzwtz(n+em+2,k)*wk6(2,k);
		rslt(ind++) = lcl0;
	}
	
	/* EDGE 6 */
	for (int n = 0; n < em;++n){
		lcl0 = rslt(ind);
		for (int k=0; k < lgpz; ++k) 
			lcl0 -= gzwtz(n+em+2,k)*wk5(3,k)+dgzwtz(n+em+2,k)*wk6(3,k);
		rslt(ind++) = lcl0;
	}
	
	/* FACE 0 */  
	ind2 = 0;
	for (int p = 1; p <= em-1; ++p){
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = rslt(ind);
			for (int k=0; k < lgpz; ++k) 
				lcl0 -= gzwtz(p+q+1,k)*wk5(3*em+3+ind2,k)+dgzwtz(p+q+1,k)*wk6(3*em+3+ind2,k);
			rslt(ind++) = lcl0;
		}
	}	

	/* FACE 1 */
	ind2 = 0;
	sign = 1;
	for (int p = 1; p <= em-1;++p){
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = rslt(ind);
			for (int k=0; k < lgpz; ++k) 
				lcl0 -= sign*(gzwtz(1+2*em+ind2,k)*wk5(3+p,k)+dgzwtz(1+2*em+ind2,k)*wk6(3+p,k));
			rslt(ind++) = lcl0;
		}
		sign*=-1; 
	}

	/* FACE 2 */
	ind2 = 0;
	sign = 1;
	for (int p = 1; p <= em-1;++p){
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = rslt(ind);
			for (int k=0; k < lgpz; ++k) 
				lcl0 -= sign*(gzwtz(1+2*em+ind2,k)*wk5(em+3+p,k)+dgzwtz(1+2*em+ind2,k)*wk6(em+3+p,k));
			rslt(ind++) = lcl0;
		}
		sign*=-1;
	}
	
	/* FACE 3 */
	ind2 = 0;
	for (int p = 1; p <= em-1;++p){
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			lcl0 = rslt(ind);
			for (int k=0; k < lgpz; ++k) 
				lcl0 -= gzwtz(1+2*em+ind2,k)*wk5(2*em+3+p,k)+dgzwtz(1+2*em+ind2,k)*wk6(2*em+3+p,k);
			rslt(ind++) = lcl0;
		}
	}

	/* INTERIOR */	
	ind2 = 0;
	int ind3 = 0;
	for (int p = 1; p <= em-1; ++p) {
		for (int q = 1; q <= em-p; ++q){
			++ind2;
			for (int r = 1; r <= em-p-q; ++r){
				++ind3;
				lcl0 = rslt(ind);
				for (int k=0; k < lgpz; ++k) 
					lcl0 -= gzwtz(1+2*em+fm+ind3,k)*wk5(3*em+3+ind2,k)+dgzwtz(1+2*em+fm+ind3,k)*wk6(3*em+3+ind2,k);
				rslt(ind++) = lcl0;
			}
		}
	}

  
   return;
}

