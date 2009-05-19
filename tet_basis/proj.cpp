/*
 *  proj.cpp
 *
 *  Created by helenbrk on Mon Oct 08 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "tet_basis.h"
#include <stdio.h>

#ifndef BZ_DEBUG
#define lin(i) lin1[i]
#define f(i,j,k) f1[(i)*stridex +(j)*stridey +k]
#define dx(i,j,k) dx1[(i)*stridex +(j)*stridey +k]
#define dy(i,j,k) dy1[(i)*stridex +(j)*stridey +k]
#define dz(i,j,k) dz1[(i)*stridex +(j)*stridey +k]
//#define dnrm(i,j) dnrm1[(i)*stridex +j]
//#define dtan(i,j) dtan1[(i)*stridex +j]
#endif

void tet_basis::proj(FLT *lin1, FLT *f1, FLT *dx1, FLT *dy1, FLT *dz1, int stridex, int stridey) {
	Array<FLT,2> wk0(4+3*em+fm,gpz);
	Array<FLT,2> wk1(4+3*em+fm,gpz);
	Array<FLT,2> wk2(4+3*em+fm,gpz);
	Array<FLT,3> wk3(3+em,gpy,gpz);
	Array<FLT,3> wk4(3+em,gpy,gpz);
	Array<FLT,3> wk5(3+em,gpy,gpz);
	Array<FLT,3> wk6(3+em,gpy,gpz);
	Array<FLT,3> wk7(3+em,gpy,gpz);
	const int be1 = 4,be2 = be1+em,be3 = be2+em,be4 = be3+em,be5 = be4+em,be6 = be5+em;
	const int bf0 = be6+em, bf1=bf0+fm,bf2=bf1+fm,bf3=bf2+fm,bi=bf3+fm;
	const int lgpx = gpx, lgpy = gpy, lgpz = gpz; 
	int sign; 
	FLT lcl0, lcl1, lcl2, lcl3;
	FLT c0, b0, b1, a0; 
#ifdef BZ_DEBUG
	Array<FLT,1> lin(lin1, shape(tm), neverDeleteData);
	Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,3> dx(dx1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,3> dy(dy1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,3> dz(dz1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif
   
   /* DETERMINE U VALUES, GRAD U VALUES
      AT COLLOCATION POINTS
      SUM HAT(U) FOR DU/DX AND DU/DY
   */  
   
   	/* Z-DIRECTION  */
	for(int n = 0; n < lgpz; ++n){
		c0=z0(n);
	
	/* VERTEX 0 */
		wk0(0,n) = lin(0)*gz(n,0);
		wk1(0,n) = lin(0)*dgz(n,0);
		wk2(0,n) = wk0(0,n)*c0;
			 
	/* VERTEX 1, EDGE 4 */
		lcl0 = lin(1)*gz(n,1);
		lcl1 = lin(1)*dgz(n,1);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be4+i)*gz(n,2+em+i);	
			lcl1 += lin(be4+i)*dgz(n,2+em+i);
		}
		wk0(1,n) = lcl0;
		wk1(1,n) = lcl1;
		wk2(1,n) = lcl0*c0;
					
	/* VERTEX 2, EDGE 5 */
		lcl0 = lin(2)*gz(n,1);
		lcl1 = lin(2)*dgz(n,1);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be5+i)*gz(n,2+em+i);	
			lcl1 += lin(be5+i)*dgz(n,2+em+i);	
		}		
		wk0(2,n) = lcl0;
		wk1(2,n) = lcl1;
		wk2(2,n) = lcl0*c0;
		
	/* VERTEX 3, EDGE 6 */
		lcl0 = lin(3)*gz(n,1);
		lcl1 = lin(3)*dgz(n,1);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be6+i)*gz(n,2+em+i);	
			lcl1 += lin(be6+i)*dgz(n,2+em+i);		
		}
		wk0(3,n) = lcl0;
		wk1(3,n) = lcl1;
		wk2(3,n) = lcl0*c0;
		
	/* EDGE 1, FACE 1 */
		int ind = 3;
		int ind2 = 0;
		sign = 1;
		for(int i = 1; i <= em; ++i){
			lcl0 = lin(i+be1-1)*gz(n,i+1);
			lcl1 = lin(i+be1-1)*dgz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				lcl0 += sign*lin(bf1+ind2)*gz(n,2+2*em+ind2);
				lcl1 += sign*lin(bf1+ind2)*dgz(n,2+2*em+ind2);
				++ind2;
			}
			sign*=-1;
			wk0(++ind,n) = lcl0;
			wk1(ind,n) = lcl1;
			wk2(ind,n) = lcl0*c0;
		}	

	/* EDGE 2, FACE 2 */
	    ind2 = 0;
		sign = 1;
		for(int i = 1; i <= em; ++i){
			lcl0 = lin(i+be2-1)*gz(n,i+1);
			lcl1 = lin(i+be2-1)*dgz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				lcl0 += sign*lin(bf2+ind2)*gz(n,2+2*em+ind2);
				lcl1 += sign*lin(bf2+ind2)*dgz(n,2+2*em+ind2);
				++ind2;
			}
			sign*=-1;
			wk0(++ind,n) = lcl0;
			wk1(ind,n) = lcl1;
			wk2(ind,n) = lcl0*c0;
		}	
	
	/* EDGE 3, FACE 3 */
	    ind2 = 0;
		for(int i = 1; i <= em; ++i){
			lcl0 = lin(i+be3-1)*gz(n,i+1);
			lcl1 = lin(i+be3-1)*dgz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				lcl0 += lin(bf3+ind2)*gz(n,2+2*em+ind2);
				lcl1 += lin(bf3+ind2)*dgz(n,2+2*em+ind2);
				++ind2;
			}
			wk0(++ind,n) = lcl0;
			wk1(ind,n) = lcl1;
			wk2(ind,n) = lcl0*c0;
		}		
	
	/* FACE 0, INTERIOR */  
		ind2 = 0;
		int ind3 = 0;
		for(int i = 1; i <= em-1; ++i){
			for(int j = 1; j <= em-i; ++j){
				lcl0 = lin(bf0+ind2)*gz(n,i+j+1);
				lcl1 = lin(bf0+ind2)*dgz(n,i+j+1);
				++ind2;
				for(int k = 1; k <= em-i-j; ++k){
					lcl0 += lin(bi+ind3)*gz(n,2+2*em+fm+ind3);
					lcl1 += lin(bi+ind3)*dgz(n,2+2*em+fm+ind3);
					++ind3;				
				}
				wk0(++ind,n) = lcl0;
				wk1(ind,n) = lcl1;
				wk2(ind,n) = lcl0*c0;				
			}
		}			
		
	}	
	
	/* Y-DIRECTION  */
	for(int n = 0; n < lgpz; ++n){
		for(int m = 0; m < lgpy; ++m){
			b0=y0(m);
			b1=y1(m);
		
		/* VERTEX 0,1 EDGE 4 */
			lcl0=gy(m,0);
			lcl1=gy(m,1);
			wk3(0,m,n) = wk0(0,n)*lcl0+wk0(1,n)*lcl1;		
			wk4(0,m,n) = b0*(wk2(0,n)*lcl0 + wk2(1,n)*lcl1);
			wk5(0,m,n) = wk1(0,n)*lcl0+wk1(1,n)*lcl1;
			wk6(0,m,n) = wk2(0,n)*dgy(m,0)+wk2(1,n)*dgy(m,1);
			wk7(0,m,n) = b1*wk6(0,m,n);
					 
		/* VERTEX 2, EDGE 3,5 FACE 3 */	
			lcl0 = wk0(2,n)*gy(m,2);
			lcl1 = wk2(2,n)*gy(m,2);
			lcl2 = wk1(2,n)*gy(m,2);
			lcl3 = wk2(2,n)*dgy(m,2);
			for(int i = 0; i < em; ++i){
				lcl0 += wk0(4+2*em+i,n)*gy(m,3+em+i);
				lcl1 += wk2(4+2*em+i,n)*gy(m,3+em+i);
				lcl2 += wk1(4+2*em+i,n)*gy(m,3+em+i);
				lcl3 += wk2(4+2*em+i,n)*dgy(m,3+em+i);		
			}			
			wk3(1,m,n) = lcl0;		
			wk4(1,m,n) = b0*lcl1;
			wk5(1,m,n) = lcl2;
			wk6(1,m,n) = lcl3;
			wk7(1,m,n) = b1*wk6(1,m,n);
			
	
		/* VERTEX 3, EDGE 2,6 FACE 2 */	
			lcl0 = wk0(3,n)*gy(m,2);
			lcl1 = wk2(3,n)*gy(m,2);
			lcl2 = wk1(3,n)*gy(m,2);
			lcl3 = wk2(3,n)*dgy(m,2);
			for(int i = 0; i < em; ++i){
				lcl0 += wk0(4+em+i,n)*gy(m,3+em+i);
				lcl1 += wk2(4+em+i,n)*gy(m,3+em+i);
				lcl2 += wk1(4+em+i,n)*gy(m,3+em+i);
				lcl3 += wk2(4+em+i,n)*dgy(m,3+em+i);
			}
			wk3(2,m,n) = lcl0;		
			wk4(2,m,n) = b0*lcl1;
			wk5(2,m,n) = lcl2;
			wk6(2,m,n) = lcl3;
			wk7(2,m,n) = b1*wk6(2,m,n);
			
		/* EDGE 1, FACE 0,1, INTERIOR */
			int ind2 = 0;			
			for(int p = 3; p < em+3; ++p){
				lcl0=wk0(p+1,n)*gy(m,p);	
				lcl1=wk2(p+1,n)*gy(m,p);
				lcl2=wk1(p+1,n)*gy(m,p);
				lcl3=wk2(p+1,n)*dgy(m,p);			
				for(int i = 1; i <= em-p+2; ++i){	
					lcl0 += wk0(4+3*em+ind2,n)*gy(m,3+2*em+ind2);
					lcl1 += wk2(4+3*em+ind2,n)*gy(m,3+2*em+ind2);
					lcl2 += wk1(4+3*em+ind2,n)*gy(m,3+2*em+ind2);
					lcl3 += wk2(4+3*em+ind2,n)*dgy(m,3+2*em+ind2);
					++ind2;
				}
				wk3(p,m,n) = lcl0;		
				wk4(p,m,n) = b0*lcl1;
				wk5(p,m,n) = lcl2;
				wk6(p,m,n) = lcl3;
				wk7(p,m,n) = b1*wk6(p,m,n);				
			}				
			
		}
	}
	
	
	/* X-DIRECTION  */
	for(int l = 0; l < lgpx; ++l){
		a0=x0(l);
		for(int m = 0; m < lgpy; ++m){
			for(int n = 0; n < lgpz; ++n){				
							/* V0,E4,V1 */      /* V2,E3,E5,F3 */     /* V3,E2,E6,F2 */
				f(l,m,n) = wk3(0,m,n)*gx(l,0) + wk3(1,m,n)*gx(l,1) + wk3(2,m,n)*gx(l,2);
				dx(l,m,n) = wk4(0,m,n)*dgx(l,0) + wk4(1,m,n)*dgx(l,1) + wk4(2,m,n)*dgx(l,2);
				dy(l,m,n) = wk6(0,m,n)*gx(l,0) + wk6(1,m,n)*gx(l,1) + wk6(2,m,n)*gx(l,2);
				dz(l,m,n) = (wk5(0,m,n)+wk7(0,m,n))*gx(l,0) + (wk5(1,m,n)+wk7(1,m,n))*gx(l,1) +(wk5(2,m,n)+wk7(2,m,n))*gx(l,2);					
				for(int i = 0; i < em; ++i){
								/* E1,F0,F1,INTERIOR */
					f(l,m,n) += wk3(3+i,m,n)*gx(l,3+i);
					dx(l,m,n) += wk4(3+i,m,n)*dgx(l,3+i);
					dy(l,m,n) += wk6(3+i,m,n)*gx(l,3+i);
					dz(l,m,n) += (wk5(3+i,m,n)+wk7(3+i,m,n))*gx(l,3+i);
				
				}
				dy(l,m,n) += a0*dx(l,m,n);
				dz(l,m,n) += a0*dx(l,m,n);
			}
		}
	}

   return;
}


void tet_basis::proj(FLT *lin1, FLT *f1, int stridex, int stridey){
	Array<FLT,2> wk0(4+3*em+fm,gpz);
	Array<FLT,3> wk1(3+em,gpy,gpz);
	const int be1 = 4,be2 = be1+em,be3 = be2+em,be4 = be3+em,be5 = be4+em,be6 = be5+em;
	const int bf0 = be6+em, bf1=bf0+fm,bf2=bf1+fm,bf3=bf2+fm,bi=bf3+fm;
	const int lgpx = gpx, lgpy = gpy, lgpz=gpz;  
	int sign;
	
#ifdef BZ_DEBUG
	Array<FLT,1> lin(lin1, shape(tm), neverDeleteData);
	Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif
	
	/* Z-DIRECTION  */
	for(int n = 0; n < lgpz; ++n){
	
	/* VERTEX 0 */
		wk0(0,n) = lin(0)*gz(n,0);
			 
	/* VERTEX 1, EDGE 4 */
		wk0(1,n) = lin(1)*gz(n,1);
		for(int i = 0; i < em; ++i)
			wk0(1,n) += lin(be4+i)*gz(n,2+em+i);	
					
	/* VERTEX 2, EDGE 5 */
		wk0(2,n) = lin(2)*gz(n,1);
		for(int i = 0; i < em; ++i)
			wk0(2,n) += lin(be5+i)*gz(n,2+em+i);	
					
	/* VERTEX 3, EDGE 6 */
		wk0(3,n) = lin(3)*gz(n,1);
		for(int i = 0; i < em; ++i)
			wk0(3,n) += lin(be6+i)*gz(n,2+em+i);			
		
	/* EDGE 1, FACE 1 */
		int ind = 3;
		int ind2 = 0;
		sign = 1;
		for(int i = 1; i <= em; ++i){
			wk0(++ind,n) = lin(i+be1-1)*gz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				wk0(ind,n) += sign*lin(bf1+ind2)*gz(n,2+2*em+ind2);
				++ind2;
			}
			sign*=-1;
		}	

	/* EDGE 2, FACE 2 */
	    ind2 = 0;
		sign = 1;
		for(int i = 1; i <= em; ++i){
			wk0(++ind,n) = lin(i+be2-1)*gz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				wk0(ind,n) += sign*lin(bf2+ind2)*gz(n,2+2*em+ind2);
				++ind2;
			}
			sign*=-1;
		}	
	
	/* EDGE 3, FACE 3 */
	    ind2 = 0;
		for(int i = 1; i <= em; ++i){
			wk0(++ind,n) = lin(i+be3-1)*gz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				wk0(ind,n) += lin(bf3+ind2)*gz(n,2+2*em+ind2);
				++ind2;
			}
		}		
	
	/* FACE 0, INTERIOR */  
		ind2 = 0;
		int ind3 = 0;
		for(int i = 1; i <= em-1; ++i){
			for(int j = 1; j <= em-i; ++j){
				wk0(++ind,n) = lin(bf0+ind2)*gz(n,i+j+1);
				++ind2;
				for(int k = 1; k <= em-i-j; ++k){
					wk0(ind,n) += lin(bi+ind3)*gz(n,2+2*em+fm+ind3);
					++ind3;				
				}
			}
		}			
		
	}	

	/* Y-DIRECTION  */
	for(int n = 0; n < lgpz; ++n){
		for(int m = 0; m < lgpy; ++m){
		
		/* VERTEX 0,1 EDGE 4 */
			wk1(0,m,n) = wk0(0,n)*gy(m,0)+wk0(1,n)*gy(m,1);			
					 
		/* VERTEX 2, EDGE 3,5 FACE 3 */	
			wk1(1,m,n) = wk0(2,n)*gy(m,2);
			for(int i = 0; i < em; ++i)
				wk1(1,m,n) += wk0(4+2*em+i,n)*gy(m,3+em+i);			
	
		/* VERTEX 3, EDGE 2,6 FACE 2 */	
			wk1(2,m,n) = wk0(3,n)*gy(m,2);
			for(int i = 0; i < em; ++i)
				wk1(2,m,n) += wk0(4+em+i,n)*gy(m,3+em+i);
		
		/* EDGE 1, FACE 0,1, INTERIOR */
			int ind2 = 0;			
			for(int p = 3; p < em+3; ++p){
				wk1(p,m,n)=wk0(p+1,n)*gy(m,p);				
				for(int i = 1; i <= em-p+2; ++i){	
					wk1(p,m,n) += wk0(4+3*em+ind2,n)*gy(m,3+2*em+ind2);
					++ind2;
				}
			}				
			
		}
	}
	
		
	/* X-DIRECTION  */
	for(int l = 0; l < lgpx; ++l){
		for(int m = 0; m < lgpy; ++m){
			for(int n = 0; n < lgpz; ++n){
							/* V0,E4,V1 */      /* V2,E3,E5,F3 */     /* V3,E2,E6,F2 */
				f(l,m,n) = wk1(0,m,n)*gx(l,0) + wk1(1,m,n)*gx(l,1) + wk1(2,m,n)*gx(l,2);					
				for(int i = 0; i < em; ++i){
								/* E1,F0,F1,INTERIOR */
					f(l,m,n) += wk1(3+i,m,n)*gx(l,3+i);
				}
			}
		}
	}
		
	
	return;
}


void tet_basis::proj(FLT u1, FLT u2, FLT u3, FLT u4, FLT *f1, int stridex, int stridey) {
   const int lgpx = gpx, lgpy = gpy, lgpz = gpz;
#ifdef BZ_DEBUG
   Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif

	for (int i=0; i < lgpx; ++i ){ 
		for (int j=0; j < lgpy; ++j ){ 
			for (int k=0; k < lgpz; ++k ){
				f(i,j,k) = u1*gz(k,0) +u2*gy(j,1)*gz(k,1) +u3*gx(i,1)*gy(j,2)*gz(k,1)+u4*gx(i,2)*gy(j,2)*gz(k,1);
			}
		}
	}

	return;
}

void tet_basis::proj_bdry(FLT *lin1, FLT *f1, FLT *dx1, FLT *dy1, FLT *dz1, int stridex, int stridey) {
	Array<FLT,2> wk0(4+3*em+fm,gpz);
	Array<FLT,2> wk1(4+3*em+fm,gpz);
	Array<FLT,2> wk2(4+3*em+fm,gpz);
	Array<FLT,3> wk3(3+em,gpy,gpz);
	Array<FLT,3> wk4(3+em,gpy,gpz);
	Array<FLT,3> wk5(3+em,gpy,gpz);
	Array<FLT,3> wk6(3+em,gpy,gpz);
	Array<FLT,3> wk7(3+em,gpy,gpz);
	const int be1 = 4,be2 = be1+em,be3 = be2+em,be4 = be3+em,be5 = be4+em,be6 = be5+em;
	const int bf0 = be6+em, bf1=bf0+fm,bf2=bf1+fm,bf3=bf2+fm;
	const int lgpx = gpx, lgpy = gpy, lgpz = gpz;  
	int sign;
	FLT lcl0, lcl1, lcl2, lcl3;
	FLT c0, b0, b1, a0; 
#ifdef BZ_DEBUG
	Array<FLT,1> lin(lin1, shape(bm), neverDeleteData);
	Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,3> dx(dx1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,3> dy(dy1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
	Array<FLT,3> dz(dz1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif
   
   /* DETERMINE U VALUES, GRAD U VALUES
      AT COLLOCATION POINTS
      SUM HAT(U) FOR DU/DX AND DU/DY
   */  
   
   	/* Z-DIRECTION  */
	for(int n = 0; n < lgpz; ++n){
		c0=z0(n);
	
	/* VERTEX 0 */
		wk0(0,n) = lin(0)*gz(n,0);
		wk1(0,n) = lin(0)*dgz(n,0);
		wk2(0,n) = wk0(0,n)*c0;
			 
	/* VERTEX 1, EDGE 4 */
		lcl0 = lin(1)*gz(n,1);
		lcl1 = lin(1)*dgz(n,1);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be4+i)*gz(n,2+em+i);	
			lcl1 += lin(be4+i)*dgz(n,2+em+i);
		}
		wk0(1,n) = lcl0;
		wk1(1,n) = lcl1;
		wk2(1,n) = lcl0*c0;
					
	/* VERTEX 2, EDGE 5 */
		lcl0 = lin(2)*gz(n,1);
		lcl1 = lin(2)*dgz(n,1);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be5+i)*gz(n,2+em+i);	
			lcl1 += lin(be5+i)*dgz(n,2+em+i);	
		}		
		wk0(2,n) = lcl0;
		wk1(2,n) = lcl1;
		wk2(2,n) = lcl0*c0;
		
	/* VERTEX 3, EDGE 6 */
		lcl0 = lin(3)*gz(n,1);
		lcl1 = lin(3)*dgz(n,1);
		for(int i = 0; i < em; ++i){
			lcl0 += lin(be6+i)*gz(n,2+em+i);	
			lcl1 += lin(be6+i)*dgz(n,2+em+i);		
		}
		wk0(3,n) = lcl0;
		wk1(3,n) = lcl1;
		wk2(3,n) = lcl0*c0;
		
	/* EDGE 1, FACE 1 */
		int ind = 3;
		int ind2 = 0;
		sign = 1;
		for(int i = 1; i <= em; ++i){
			lcl0 = lin(i+be1-1)*gz(n,i+1);
			lcl1 = lin(i+be1-1)*dgz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				lcl0 += sign*lin(bf1+ind2)*gz(n,2+2*em+ind2);
				lcl1 += sign*lin(bf1+ind2)*dgz(n,2+2*em+ind2);
				++ind2;
			}
			sign*=-1;
			wk0(++ind,n) = lcl0;
			wk1(ind,n) = lcl1;
			wk2(ind,n) = lcl0*c0;
		}	

	/* EDGE 2, FACE 2 */
	    ind2 = 0;
		sign = 1;
		for(int i = 1; i <= em; ++i){
			lcl0 = lin(i+be2-1)*gz(n,i+1);
			lcl1 = lin(i+be2-1)*dgz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				lcl0 += sign*lin(bf2+ind2)*gz(n,2+2*em+ind2);
				lcl1 += sign*lin(bf2+ind2)*dgz(n,2+2*em+ind2);
				++ind2;
			}
			sign*=-1;
			wk0(++ind,n) = lcl0;
			wk1(ind,n) = lcl1;
			wk2(ind,n) = lcl0*c0;
		}	
	
	/* EDGE 3, FACE 3 */
	    ind2 = 0;
		for(int i = 1; i <= em; ++i){
			lcl0 = lin(i+be3-1)*gz(n,i+1);
			lcl1 = lin(i+be3-1)*dgz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				lcl0 += lin(bf3+ind2)*gz(n,2+2*em+ind2);
				lcl1 += lin(bf3+ind2)*dgz(n,2+2*em+ind2);
				++ind2;
			}
			wk0(++ind,n) = lcl0;
			wk1(ind,n) = lcl1;
			wk2(ind,n) = lcl0*c0;
		}		
	
	/* FACE 0 */  
		ind2 = 0;
		for(int i = 1; i <= em-1; ++i){
			for(int j = 1; j <= em-i; ++j){ 
				wk0(++ind,n) = lin(bf0+ind2)*gz(n,i+j+1);
				wk1(ind,n) = lin(bf0+ind2)*dgz(n,i+j+1);
				wk2(ind,n) = wk0(ind,n)*c0;
				++ind2;				
			}
		}			
		
	}	
	
	/* Y-DIRECTION  */
	for(int n = 0; n < lgpz; ++n){
		for(int m = 0; m < lgpy; ++m){
			b0=y0(m);
			b1=y1(m);
		
		/* VERTEX 0,1 EDGE 4 */
			lcl0=gy(m,0);
			lcl1=gy(m,1);
			wk3(0,m,n) = wk0(0,n)*lcl0+wk0(1,n)*lcl1;		
			wk4(0,m,n) = b0*(wk2(0,n)*lcl0 + wk2(1,n)*lcl1);
			wk5(0,m,n) = wk1(0,n)*lcl0+wk1(1,n)*lcl1;
			wk6(0,m,n) = wk2(0,n)*dgy(m,0)+wk2(1,n)*dgy(m,1);
			wk7(0,m,n) = b1*wk6(0,m,n);
					 
		/* VERTEX 2, EDGE 3,5 FACE 3 */	
			lcl0 = wk0(2,n)*gy(m,2);
			lcl1 = wk2(2,n)*gy(m,2);
			lcl2 = wk1(2,n)*gy(m,2);
			lcl3 = wk2(2,n)*dgy(m,2);
			for(int i = 0; i < em; ++i){
				lcl0 += wk0(4+2*em+i,n)*gy(m,3+em+i);
				lcl1 += wk2(4+2*em+i,n)*gy(m,3+em+i);
				lcl2 += wk1(4+2*em+i,n)*gy(m,3+em+i);
				lcl3 += wk2(4+2*em+i,n)*dgy(m,3+em+i);		
			}			
			wk3(1,m,n) = lcl0;		
			wk4(1,m,n) = b0*lcl1;
			wk5(1,m,n) = lcl2;
			wk6(1,m,n) = lcl3;
			wk7(1,m,n) = b1*wk6(1,m,n);
			
	
		/* VERTEX 3, EDGE 2,6 FACE 2 */	
			lcl0 = wk0(3,n)*gy(m,2);
			lcl1 = wk2(3,n)*gy(m,2);
			lcl2 = wk1(3,n)*gy(m,2);
			lcl3 = wk2(3,n)*dgy(m,2);
			for(int i = 0; i < em; ++i){
				lcl0 += wk0(4+em+i,n)*gy(m,3+em+i);
				lcl1 += wk2(4+em+i,n)*gy(m,3+em+i);
				lcl2 += wk1(4+em+i,n)*gy(m,3+em+i);
				lcl3 += wk2(4+em+i,n)*dgy(m,3+em+i);
			}
			wk3(2,m,n) = lcl0;		
			wk4(2,m,n) = b0*lcl1;
			wk5(2,m,n) = lcl2;
			wk6(2,m,n) = lcl3;
			wk7(2,m,n) = b1*wk6(2,m,n);
			
		/* EDGE 1, FACE 0,1*/
			int ind2 = 0;			
			for(int p = 3; p < em+3; ++p){
				lcl0=wk0(p+1,n)*gy(m,p);	
				lcl1=wk2(p+1,n)*gy(m,p);
				lcl2=wk1(p+1,n)*gy(m,p);
				lcl3=wk2(p+1,n)*dgy(m,p);			
				for(int i = 1; i <= em-p+2; ++i){	
					lcl0 += wk0(4+3*em+ind2,n)*gy(m,3+2*em+ind2);
					lcl1 += wk2(4+3*em+ind2,n)*gy(m,3+2*em+ind2);
					lcl2 += wk1(4+3*em+ind2,n)*gy(m,3+2*em+ind2);
					lcl3 += wk2(4+3*em+ind2,n)*dgy(m,3+2*em+ind2);
					++ind2;
				}
				wk3(p,m,n) = lcl0;		
				wk4(p,m,n) = b0*lcl1;
				wk5(p,m,n) = lcl2;
				wk6(p,m,n) = lcl3;
				wk7(p,m,n) = b1*wk6(p,m,n);				
			}				
			
		}
	}
	
	
	/* X-DIRECTION  */
	for(int l = 0; l < lgpx; ++l){
		a0=x0(l);
		for(int m = 0; m < lgpy; ++m){
			for(int n = 0; n < lgpz; ++n){				
							/* V0,E4,V1 */      /* V2,E3,E5,F3 */     /* V3,E2,E6,F2 */
				f(l,m,n) = wk3(0,m,n)*gx(l,0) + wk3(1,m,n)*gx(l,1) + wk3(2,m,n)*gx(l,2);
				dx(l,m,n) = wk4(0,m,n)*dgx(l,0) + wk4(1,m,n)*dgx(l,1) + wk4(2,m,n)*dgx(l,2);
				dy(l,m,n) = wk6(0,m,n)*gx(l,0) + wk6(1,m,n)*gx(l,1) + wk6(2,m,n)*gx(l,2);
				dz(l,m,n) = (wk5(0,m,n)+wk7(0,m,n))*gx(l,0) + (wk5(1,m,n)+wk7(1,m,n))*gx(l,1) +(wk5(2,m,n)+wk7(2,m,n))*gx(l,2);					
				for(int i = 0; i < em; ++i){
								/* E1,F0,F1 */
					f(l,m,n) += wk3(3+i,m,n)*gx(l,3+i);
					dx(l,m,n) += wk4(3+i,m,n)*dgx(l,3+i);
					dy(l,m,n) += wk6(3+i,m,n)*gx(l,3+i);
					dz(l,m,n) += (wk5(3+i,m,n)+wk7(3+i,m,n))*gx(l,3+i);
				
				}
				dy(l,m,n) += a0*dx(l,m,n);
				dz(l,m,n) += a0*dx(l,m,n);
			}
		}
	}

   return;
}

void tet_basis::proj_bdry(FLT *lin1, FLT *f1, int stridex, int stridey) {
	Array<FLT,2> wk0(4+3*em+fm,gpz);
	Array<FLT,3> wk1(3+em,gpy,gpz);
	const int be1 = 4,be2 = be1+em,be3 = be2+em,be4 = be3+em,be5 = be4+em,be6 = be5+em;
	const int bf0 = be6+em, bf1=bf0+fm,bf2=bf1+fm,bf3=bf2+fm;
	const int lgpx = gpx, lgpy = gpy, lgpz=gpz;  
	int sign;
#ifdef BZ_DEBUG
	Array<FLT,1> lin(lin1, shape(bm), neverDeleteData);
	Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif
	
	/* Z-DIRECTION  */
	for(int n = 0; n < lgpz; ++n){
	
	/* VERTEX 0 */
		wk0(0,n) = lin(0)*gz(n,0);
			 
	/* VERTEX 1, EDGE 4 */
		wk0(1,n) = lin(1)*gz(n,1);
		for(int i = 0; i < em; ++i)
			wk0(1,n) += lin(be4+i)*gz(n,2+em+i);	
					
	/* VERTEX 2, EDGE 5 */
		wk0(2,n) = lin(2)*gz(n,1);
		for(int i = 0; i < em; ++i)
			wk0(2,n) += lin(be5+i)*gz(n,2+em+i);	
					
	/* VERTEX 3, EDGE 6 */
		wk0(3,n) = lin(3)*gz(n,1);
		for(int i = 0; i < em; ++i)
			wk0(3,n) += lin(be6+i)*gz(n,2+em+i);			
		
	/* EDGE 1, FACE 1 */
		int ind = 3;
		int ind2 = 0;
		sign = 1;
		for(int i = 1; i <= em; ++i){
			wk0(++ind,n) = lin(i+be1-1)*gz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				wk0(ind,n) += sign*lin(bf1+ind2)*gz(n,2+2*em+ind2);
				++ind2;
			}
			sign*=-1;
		}	

	/* EDGE 2, FACE 2 */
	    ind2 = 0;
		sign = 1;
		for(int i = 1; i <= em; ++i){
			wk0(++ind,n) = lin(i+be2-1)*gz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				wk0(ind,n) += sign*lin(bf2+ind2)*gz(n,2+2*em+ind2);
				++ind2;
			}
			sign*=-1;
		}	
	
	/* EDGE 3, FACE 3 */
	    ind2 = 0;
		for(int i = 1; i <= em; ++i){
			wk0(++ind,n) = lin(i+be3-1)*gz(n,i+1);
			for(int j = 1; j <= em-i; ++j){
				wk0(ind,n) += lin(bf3+ind2)*gz(n,2+2*em+ind2);
				++ind2;
			}
		}		
	
	/* FACE 0 */  
		ind2 = 0;
		for(int i = 1; i <= em-1; ++i){
			for(int j = 1; j <= em-i; ++j){
				wk0(++ind,n) = lin(bf0+ind2)*gz(n,i+j+1);
				++ind2;
			}
		}			
		
	}	

	/* Y-DIRECTION  */
	for(int n = 0; n < lgpz; ++n){
		for(int m = 0; m < lgpy; ++m){
		
		/* VERTEX 0,1 EDGE 4 */
			wk1(0,m,n) = wk0(0,n)*gy(m,0)+wk0(1,n)*gy(m,1);			
					 
		/* VERTEX 2, EDGE 3,5 FACE 3 */	
			wk1(1,m,n) = wk0(2,n)*gy(m,2);
			for(int i = 0; i < em; ++i)
				wk1(1,m,n) += wk0(4+2*em+i,n)*gy(m,3+em+i);			
	
		/* VERTEX 3, EDGE 2,6 FACE 2 */	
			wk1(2,m,n) = wk0(3,n)*gy(m,2);
			for(int i = 0; i < em; ++i)
				wk1(2,m,n) += wk0(4+em+i,n)*gy(m,3+em+i);
		
		/* EDGE 1, FACE 0,1 */
			int ind2 = 0;			
			for(int p = 3; p < em+3; ++p){
				wk1(p,m,n)=wk0(p+1,n)*gy(m,p);				
				for(int i = 1; i <= em-p+2; ++i){	
					wk1(p,m,n) += wk0(4+3*em+ind2,n)*gy(m,3+2*em+ind2);
					++ind2;
				}
			}				
			
		}
	}
	
		
	/* X-DIRECTION  */
	for(int p = 0; p < lgpx; ++p){
		for(int m = 0; m < lgpy; ++m){
			for(int n = 0; n < lgpz; ++n){
							/* V0,E4,V1 */      /* V2,E3,E5,F3 */     /* V3,E2,E6,F2 */
				f(p,m,n) = wk1(0,m,n)*gx(p,0) + wk1(1,m,n)*gx(p,1) + wk1(2,m,n)*gx(p,2);					
				for(int i = 0; i < em; ++i){
								/* E1,F0,F1 */
					f(p,m,n) += wk1(3+i,m,n)*gx(p,3+i);
				}
			}
		}
	}
		
	
	return;
}

void tet_basis::derivr(FLT *f1, FLT *dx1, int stridex, int stridey) {
   Array<FLT,1> wk0(gpx);
   const int lgpy = gpy, lgpx = gpx, lgpz = gpz;
   FLT lx0, lcl0;
#ifdef BZ_DEBUG
   Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
   Array<FLT,3> dx(dx1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif

	for (int j = 0; j < lgpy; ++j){
		for (int k = 0; k < lgpz; ++k){		
			lx0=y0(j)*z0(k);
			for (int i = 0; i < lgpx; ++i)
				wk0(i) = f(i,j,k)*dltx(i)*lx0;
			for (int i = 0; i < lgpx; ++i){
				lcl0 = dx(i,j,k);
				for (int n = 0; n < i; ++n)
					lcl0 += (wk0(i) + wk0(n))*dltx2(i,n);
				for(int n=i+1;n<lgpx;++n) 
					lcl0 += (wk0(i) + wk0(n))*dltx2(i,n);
				dx(i,j,k) = lcl0;
			}
		}
	}

   return;
}

void tet_basis::derivs(FLT *f1, FLT *dy1, int stridex, int stridey) {
    Array<FLT,1> wk0(gpx);
   const int lgpy = gpy, lgpx = gpx, lgpz = gpz;
   FLT ly0, lcl0;
#ifdef BZ_DEBUG
   Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
   Array<FLT,3> dy(dy1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif

   for (int j = 0; j < lgpy; ++j){
		for (int k = 0; k < lgpz; ++k){		
			ly0=y0(j)*z0(k);
			for (int i = 0; i < lgpx; ++i) 
				wk0(i) = f(i,j,k)*dltx(i)*ly0;
			for (int i = 0; i < lgpx; ++i){
				lcl0 = dy(i,j,k);
				for (int n = 0; n < i; ++n)
					lcl0 += (wk0(i) + wk0(n))*dltx3(i,n);
				for(int n=i+1;n<lgpx;++n) 
					lcl0 += (wk0(i) + wk0(n))*dltx3(i,n);
				dy(i,j,k) = lcl0;
			}
		}
	}
	
	for (int i = 0; i < lgpx; ++i){	
		for (int k = 0; k < lgpz; ++k){
			ly0=z0(k);				
			for (int j = 0; j < lgpy; ++j)
				wk0(j) = f(i,j,k)*dlty(j)*ly0;
			for (int j = 0; j < lgpy; ++j){
				lcl0 = dy(i,j,k);
				for (int n = 0; n < j; ++n)
					lcl0 += (wk0(j) + wk0(n))*dlty2(j,n);
				for(int n = j+1; n < lgpy; ++n) 
					lcl0 += (wk0(j) + wk0(n))*dlty2(j,n);
				dy(i,j,k) = lcl0;
			}
		}
	}

   return;
}

void tet_basis::derivt(FLT *f1, FLT *dz1, int stridex, int stridey) {
    Array<FLT,1> wk0(gpx);
   const int lgpy = gpy, lgpx = gpx, lgpz = gpz;
   FLT lz0, lcl0;
#ifdef BZ_DEBUG
   Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
   Array<FLT,3> dz(dz1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif

	for (int j = 0; j < lgpy; ++j){
		for (int k = 0; k < lgpz; ++k){		
			lz0=y0(j)*z0(k);
			for (int i = 0; i < lgpx; ++i)
				wk0(i) = f(i,j,k)*dltx(i)*lz0;
			for (int i = 0; i < lgpx; ++i){
				lcl0 = dz(i,j,k);
				for (int n = 0; n < i; ++n)
					lcl0 += (wk0(i) + wk0(n))*dltx3(i,n);
				for(int n = i+1; n < lgpx;++n) 
					lcl0 += (wk0(i) + wk0(n))*dltx3(i,n);
				dz(i,j,k) = lcl0;
			}
		}
	}
	
	for (int i = 0; i < lgpx; ++i){
		for (int k = 0; k < lgpz; ++k){
			lz0=z0(k);		
			for (int j = 0; j < lgpy; ++j)
				wk0(j) = f(i,j,k)*dlty(j)*lz0;
			for (int j = 0; j < lgpy; ++j){
				lcl0 = dz(i,j,k);
				for (int n = 0; n < j; ++n)
					lcl0 += (wk0(j) + wk0(n))*dlty3(j,n);
				for(int n=j+1; n < lgpy; ++n) 
					lcl0 += (wk0(j) + wk0(n))*dlty3(j,n);
				dz(i,j,k) = lcl0;
			}
		}
	}
	
	for (int i = 0; i < lgpx; ++i){
		for (int j = 0; j < lgpy; ++j){	
			for (int k = 0; k < lgpz; ++k)
				wk0(k) = f(i,j,k)*dltz(k);
			for (int k = 0; k < lgpz; ++k){
				lcl0 = dz(i,j,k);
				for (int n = 0; n < k; ++n)
					lcl0 += (wk0(k) + wk0(n))*dltz2(k,n);
				for(int n = k+1; n < lgpx; ++n) 
					lcl0 += (wk0(k) + wk0(n))*dltz2(k,n);
				dz(i,j,k) = lcl0;
			}
		}
	}

   return;
}

   
void tet_basis::proj_leg(FLT *lin1, FLT *f1, int stridex, int stridey) {
   const int ltm = tm;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,1> lin(lin1, shape(tm), neverDeleteData);
   Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif
   
   /* INTERIOR */
	for(int i = 1; i < em; ++i) {
		for(int j = 1; j < em+1-i; ++j) {
			for(int k = 1; k < em+1-i-j; ++k){
				lcl0 = 0.0;
				for(int m = 0; m < ltm; ++m)
					lcl0 += lin(m)*lgrnge3d(m,i,j,k);
				f(i,j,k) = lcl0;
			}
		}
	}
   
   return;
}

void tet_basis::proj_leg(FLT u1, FLT u2, FLT u3, FLT u4, FLT *f1, int stridex, int stridey) {
#ifdef BZ_DEBUG
   Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif
   
   /* INTERIOR */
	for(int i = 1; i < em; ++i)
		for(int j = 1; j < em+1-i; ++j)
			for(int k = 1; k < em+1-i-j; ++k)
				f(i,j,k) = u1*lgrnge3d(0,i,j,k) +u2*lgrnge3d(1,i,j,k) +u3*lgrnge3d(2,i,j,k)+u4*lgrnge3d(3,i,j,k);
   
   return;
}

void tet_basis::proj_bdry_leg(FLT *lin1, FLT *f1, int stridex, int stridey) {
   int lbm = bm;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,1> lin(lin1, shape(lbm), neverDeleteData);
   Array<FLT,3> f(f1, shape(gpx,stridex/stridey,stridey), neverDeleteData);
#endif
   
   /* INTERIOR */
	for(int i = 1; i < em; ++i) {
		for(int j = 1; j < em+1-i; ++j) {
			for(int k = 1; k < em+1-i-j; ++k){
				lcl0 = 0.0;
				for(int m = 0; m < lbm; ++m)
					lcl0 += lin(m)*lgrnge3d(m,i,j,k);
				f(i,j,k) = lcl0;
			}
		}
	}
   
   return;
}

//void tet_basis::proj_face(FLT *lin1, FLT *f1, FLT *dtan1, FLT dtan2, FLT *dnrm1) {
//	const int lgpx = gpx, lsm = em, ltm = 3+3*em+fm, stride = gpx;
//	Array<FLT,1> wk0(ltm);
//#ifdef BZ_DEBUG
//	Array<FLT,1> lin(lin1, shape(ltm), neverDeleteData);
//	Array<FLT,2> dtan(dtan1, shape(gpx,gpy), neverDeleteData);
//	Array<FLT,2> dnrm(dnrm1, shape(gpx,gpy), neverDeleteData);
//#endif
//	
//   /* NORMAL DERIVATIVE */
//   for(int i=0;i<lgpx;++i) {
//      dnrm(i) = 0.0;
//      for(int m=0;m<ltm;++m)
//         dnrm(i) += lin(m)*dgnorm(side,m,i);
//   }
//
///*	MOVE SIDE & VERTEX MODES TO WK & THEN PROJECT VALUES & TANGENT DERIVS */
//   wk0(0) = lin((side+1)%3);
//   wk0(1) = lin((side+2)%3);
//
//   int base = side*lsm;
//   for(int m=0;m<lsm;++m)
//      wk0(m+2) = lin(base +3 +m);
//      
//	proj2d(&wk0(0),f1,dtan,dtan, stride);
//      
//	return;
//}
