/*
 *  probe.cpp
 *  
 *
 *  Created by helenbrk on Thu Oct 18 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"tet_basis.h"
#include<math.h>
#include<utilities.h>

#ifndef BZ_DEBUG
#define lin(n,m) lin1[(n)*stride +m]
#define sin(n,k) sin1[(n)*stride +k]
#define f(i) f1[(i)]
#define dx(i) dx1[(i)]
#define dy(i) dy1[(i)]
#define dz(i) dz1[(i)]

#endif

/*must call ptvalues before using*/
void tet_basis::ptprobe(int nv, FLT *f1, FLT *lin1, int stride) {
	FLT lcl0;
	int ind,ind2,ind3,sign;
#ifdef BZ_DEBUG
	Array<FLT,2> lin(lin1, shape(nv,stride), neverDeleteData);
	Array<FLT,1> f(f1, shape(nv), neverDeleteData);
#endif
   
   for(int n = 0; n < nv; ++n) {
   
		lcl0 = 0.0;
		ind = 0;
		
		/* SUM ALL MODE CONTRIBUTIONS */
		/* VERTEX 0 */
		lcl0 += lin(n,ind)*pgz(0); //note: pgx,pgy=1
      
		/* VERTEX 1 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(1); //pgx=1
		
		/* VERTEX 2 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(2)*pgx(1);
		
		/* VERTEX 3 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(2)*pgx(2);
		
		/* EDGE 1 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(m+1)*pgx(m+1);
		}
		
		/* EDGE 2 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(em+m+1)*pgx(2);
		}
		
		/* EDGE 3 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(em+m+1)*pgx(1);
		}
		
		/* EDGE 4 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(1); //pgx=1
		}
		
		/* EDGE 5 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(2)*pgx(1);
		}
		
		/* EDGE 6 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(2)*pgx(2);
		}
		
		/* FACE 0 */
		ind2 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += lin(n,++ind)*pgz(p+q+1)*pgy(2+2*em+ind2)*pgx(p+2);
			}
		}
		
		/* FACE 1 */
		ind2 = 0;
		sign = 1;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += sign*lin(n,++ind)*pgz(1+2*em+ind2)*pgy(p+2)*pgx(p+2);
			}
			sign*=-1;
		}
		
		/* FACE 2 */
		ind2 = 0;
		sign = 1;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += sign*lin(n,++ind)*pgz(1+2*em+ind2)*pgy(em+2+p)*pgx(2);
			}
			sign*=-1;
		}
		
		/* FACE 3 */
		ind2 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += lin(n,++ind)*pgz(1+2*em+ind2)*pgy(em+2+p)*pgx(1);
			}
		}
		
		/* INTERIOR */
		ind2 = 0;
		ind3 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind3;
				for(int r = 1; r <= em-p-q; ++r){
					++ind2;
					lcl0 += lin(n,++ind)*pgz(1+2*em+fm+ind2)*pgy(2+2*em+ind3)*pgx(p+2);
				}
			}
		}
		
		f(n) = lcl0;
   }
   
   return;
}


/*must call ptvalues before using*/
void tet_basis::ptprobe(int nv, FLT *f1, FLT *dx1, FLT *dy1, FLT *dz1, FLT r, FLT s, FLT t, FLT *lin1, int stride) {
	FLT lcl0,lcl1,lcl2,lcl3,x,y,dxdr,dxds,dxdt,dyds,dydt;
	int ind,ind2,ind3,sign;
#ifdef BZ_DEBUG
	Array<FLT,2> lin(lin1, shape(nv,stride), neverDeleteData);
	Array<FLT,1> f(f1, shape(nv), neverDeleteData);
	Array<FLT,1> dx(dx1, shape(nv), neverDeleteData);
	Array<FLT,1> dy(dy1, shape(nv), neverDeleteData);
	Array<FLT,1> dz(dz1, shape(nv), neverDeleteData);
#endif
   
   x=2.0*(1.0+r)/(-s-t)-1.0;
   y=2.0*(1+s)/(1-t)-1.0;
   dxdr=2.0/(-s-t);
   dxds=(2.0+2.0*r)/((-s-t)*(-s-t));
   dxdt=dxds;
   dyds=2.0/(1-t);
   dydt=(2.0+2.0*s)/((1.0-t)*(1.0-t));
   ptvalues_deriv(x,y,t);
   
   for(int n = 0; n < nv; ++n) {
   
		lcl0 = 0.0;
		lcl1 = 0.0;
		lcl2 = 0.0;
		lcl3 = 0.0;

		ind = 0;
		/* SUM ALL MODE CONTRIBUTIONS */
		/* VERTEX 0 */
		lcl0 += lin(n,ind)*pgz(0); //note: pgx,pgy=1
		lcl3 += lin(n,ind)*dpgz(0); //note: pgx,pgy=1
      
		/* VERTEX 1 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(1); //pgx=1
		lcl2 += lin(n,ind)*pgz(1)*dpgy(1)*dyds; //pgx=1
		lcl3 += lin(n,ind)*(pgz(1)*dpgy(1)*dydt+dpgz(1)*pgy(1)); //pgx=1

		
		/* VERTEX 2 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(2)*pgx(1);
		lcl1 += lin(n,ind)*pgz(1)*pgy(2)*dpgx(1)*dxdr;
		lcl2 += lin(n,ind)*(pgz(1)*pgy(2)*dpgx(1)*dxds+pgz(1)*dpgy(2)*pgx(1)*dyds);
		lcl3 += lin(n,ind)*(pgz(1)*pgy(2)*dpgx(1)*dxdt+pgz(1)*dpgy(2)*pgx(1)*dydt+dpgz(1)*pgy(2)*pgx(1));

		
		/* VERTEX 3 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(2)*pgx(2);
		lcl1 += lin(n,ind)*pgz(1)*pgy(2)*dpgx(2)*dxdr;
		lcl2 += lin(n,ind)*(pgz(1)*pgy(2)*dpgx(2)*dxds+pgz(1)*dpgy(2)*pgx(2)*dyds);
		lcl3 += lin(n,ind)*(pgz(1)*pgy(2)*dpgx(2)*dxdt+pgz(1)*dpgy(2)*pgx(2)*dydt+dpgz(1)*pgy(2)*pgx(2));

		
		/* EDGE 1 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(m+1)*pgx(m+1);
			lcl1 += lin(n,ind)*pgz(m)*pgy(m+1)*dpgx(m+1)*dxdr;
			lcl2 += lin(n,ind)*(pgz(m)*pgy(m+1)*dpgx(m+1)*dxds+pgz(m)*dpgy(m+1)*pgx(m+1)*dyds);
			lcl3 += lin(n,ind)*(pgz(m)*pgy(m+1)*dpgx(m+1)*dxdt+pgz(m)*dpgy(m+1)*pgx(m+1)*dydt+dpgz(m)*pgy(m+1)*pgx(m+1));
		}
		
		/* EDGE 2 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(em+m+1)*pgx(2);		
			lcl1 += lin(n,ind)*pgz(m)*pgy(em+m+1)*dpgx(2)*dxdr;
			lcl2 += lin(n,ind)*(pgz(m)*pgy(em+m+1)*dpgx(2)*dxds+pgz(m)*dpgy(em+m+1)*pgx(2)*dyds);
			lcl3 += lin(n,ind)*(pgz(m)*pgy(em+m+1)*dpgx(2)*dxdt+pgz(m)*dpgy(em+m+1)*pgx(2)*dydt+dpgz(m)*pgy(em+m+1)*pgx(2));
		}
		
		/* EDGE 3 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(em+m+1)*pgx(1);
			lcl1 += lin(n,ind)*pgz(m)*pgy(em+m+1)*dpgx(1)*dxdr;
			lcl2 += lin(n,ind)*(pgz(m)*pgy(em+m+1)*dpgx(1)*dxds+pgz(m)*dpgy(em+m+1)*pgx(1)*dyds);
			lcl3 += lin(n,ind)*(pgz(m)*pgy(em+m+1)*dpgx(1)*dxdt+pgz(m)*dpgy(em+m+1)*pgx(1)*dydt+dpgz(m)*pgy(em+m+1)*pgx(1));
		}
		
		/* EDGE 4 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(1); //pgx=1
			lcl2 += lin(n,ind)*pgz(em+m)*dpgy(1)*dyds; //pgx=1
			lcl3 += lin(n,ind)*(pgz(em+m)*dpgy(1)*dydt+dpgz(em+m)*pgy(1)); //pgx=1

		}
		
		/* EDGE 5 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(2)*pgx(1);
			lcl1 += lin(n,ind)*pgz(em+m)*pgy(2)*dpgx(1)*dxdr;
			lcl2 += lin(n,ind)*(pgz(em+m)*pgy(2)*dpgx(1)*dxds+pgz(em+m)*dpgy(2)*pgx(1)*dyds);
			lcl3 += lin(n,ind)*(pgz(em+m)*pgy(2)*dpgx(1)*dxdt+pgz(em+m)*dpgy(2)*pgx(1)*dydt+dpgz(em+m)*pgy(2)*pgx(1));

		}
		
		/* EDGE 6 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(2)*pgx(2);
			lcl1 += lin(n,ind)*pgz(em+m)*pgy(2)*dpgx(2)*dxdr;
			lcl2 += lin(n,ind)*(pgz(em+m)*pgy(2)*dpgx(2)*dxds+pgz(em+m)*dpgy(2)*pgx(2)*dyds);
			lcl3 += lin(n,ind)*(pgz(em+m)*pgy(2)*dpgx(2)*dxdt+pgz(em+m)*dpgy(2)*pgx(2)*dydt+dpgz(em+m)*pgy(2)*pgx(2));
		}
		
		/* FACE 0 */
		ind2 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += lin(n,++ind)*pgz(p+q+1)*pgy(2+2*em+ind2)*pgx(p+2);
				lcl1 += lin(n,ind)*pgz(p+q+1)*pgy(2+2*em+ind2)*dpgx(p+2)*dxdr;
				lcl2 += lin(n,ind)*(pgz(p+q+1)*pgy(2+2*em+ind2)*dpgx(p+2)*dxds+pgz(p+q+1)*dpgy(2+2*em+ind2)*pgx(p+2)*dyds);
				lcl3 += lin(n,ind)*(pgz(p+q+1)*pgy(2+2*em+ind2)*dpgx(p+2)*dxdt+pgz(p+q+1)*dpgy(2+2*em+ind2)*pgx(p+2)*dydt+dpgz(p+q+1)*pgy(2+2*em+ind2)*pgx(p+2));
			}
		}
		
		/* FACE 1 */
		ind2 = 0;
		sign = 1;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += sign*lin(n,++ind)*pgz(1+2*em+ind2)*pgy(p+2)*pgx(p+2);
				lcl1 += sign*lin(n,ind)*pgz(1+2*em+ind2)*pgy(p+2)*dpgx(p+2)*dxdr;
				lcl2 += sign*lin(n,ind)*(pgz(1+2*em+ind2)*pgy(p+2)*dpgx(p+2)*dxds+pgz(1+2*em+ind2)*dpgy(p+2)*pgx(p+2)*dyds);
				lcl3 += sign*lin(n,ind)*(pgz(1+2*em+ind2)*pgy(p+2)*dpgx(p+2)*dxdt+pgz(1+2*em+ind2)*dpgy(p+2)*pgx(p+2)*dydt+dpgz(1+2*em+ind2)*pgy(p+2)*pgx(p+2));
			}
			sign*=-1;
		}
		
		/* FACE 2 */
		ind2 = 0;
		sign = 1;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += sign*lin(n,++ind)*pgz(1+2*em+ind2)*pgy(em+2+p)*pgx(2);
				lcl1 += sign*lin(n,ind)*pgz(1+2*em+ind2)*pgy(em+2+p)*dpgx(2)*dxdr;
				lcl2 += sign*lin(n,ind)*(pgz(1+2*em+ind2)*pgy(em+2+p)*dpgx(2)*dxds+pgz(1+2*em+ind2)*dpgy(em+2+p)*pgx(2)*dyds);
				lcl3 += sign*lin(n,ind)*(pgz(1+2*em+ind2)*pgy(em+2+p)*dpgx(2)*dxdt+pgz(1+2*em+ind2)*dpgy(em+2+p)*pgx(2)*dydt+dpgz(1+2*em+ind2)*pgy(em+2+p)*pgx(2));
			}
			sign*=-1;
		}
		
		/* FACE 3 */
		ind2 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += lin(n,++ind)*pgz(1+2*em+ind2)*pgy(em+2+p)*pgx(1);
				lcl1 += lin(n,ind)*pgz(1+2*em+ind2)*pgy(em+2+p)*dpgx(1)*dxdr;
				lcl2 += lin(n,ind)*(pgz(1+2*em+ind2)*pgy(em+2+p)*dpgx(1)*dxds+pgz(1+2*em+ind2)*dpgy(em+2+p)*pgx(1)*dyds);
				lcl3 += lin(n,ind)*(pgz(1+2*em+ind2)*pgy(em+2+p)*dpgx(1)*dxdt+pgz(1+2*em+ind2)*dpgy(em+2+p)*pgx(1)*dydt+dpgz(1+2*em+ind2)*pgy(em+2+p)*pgx(1));
			}
		}
		
		/* INTERIOR */
		ind2 = 0;
		ind3 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind3;
				for(int r = 1; r <= em-p-q; ++r){
					++ind2;
					lcl0 += lin(n,++ind)*pgz(1+2*em+fm+ind2)*pgy(2+2*em+ind3)*pgx(p+2);
					lcl1 += lin(n,ind)*pgz(1+2*em+fm+ind2)*pgy(2+2*em+ind3)*dpgx(p+2)*dxdr;
					lcl2 += lin(n,ind)*(pgz(1+2*em+fm+ind2)*pgy(2+2*em+ind3)*dpgx(p+2)*dxds+pgz(1+2*em+fm+ind2)*dpgy(2+2*em+ind3)*pgx(p+2)*dyds);
					lcl3 += lin(n,ind)*(pgz(1+2*em+fm+ind2)*pgy(2+2*em+ind3)*dpgx(p+2)*dxdt+pgz(1+2*em+fm+ind2)*dpgy(2+2*em+ind3)*pgx(p+2)*dydt+dpgz(1+2*em+fm+ind2)*pgy(2+2*em+ind3)*pgx(p+2));
				}
			}
		}
		
		f(n) = lcl0;
		dx(n) = lcl1;
		dy(n) = lcl2;
		dz(n) = lcl3;
   }
   cout << "basis " << f(0) << dx(0) << dy(0) << dz(0) << endl;
   return;
}



/*must call ptvalues before using*/
void tet_basis::ptprobe_bdry(int nv, FLT *f1, FLT *lin1, int stride){
	FLT lcl0;
	int ind,ind2,sign;
#ifdef BZ_DEBUG
	Array<FLT,2> lin(lin1, shape(nv,stride), neverDeleteData);
	Array<FLT,1> f(f1, shape(nv), neverDeleteData);
#endif
   
   
   for(int n = 0; n < nv; ++n) {
   
		lcl0 = 0.0;
		ind = 0;
		
		/* SUM ALL MODE CONTRIBUTIONS */
		/* VERTEX 0 */
		lcl0 += lin(n,ind)*pgz(0); //note: pgx,pgy=1
      
		/* VERTEX 1 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(1); //pgx=1
		
		/* VERTEX 2 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(2)*pgx(1);
		
		/* VERTEX 3 */
		lcl0 += lin(n,++ind)*pgz(1)*pgy(2)*pgx(2);
		
		/* EDGE 1 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(m+1)*pgx(m+1);
		}
		
		/* EDGE 2 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(em+m+1)*pgx(2);
		}
		
		/* EDGE 3 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(m)*pgy(em+m+1)*pgx(1);
		}
		
		/* EDGE 4 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(1); //pgx=1
		}
		
		/* EDGE 5 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(2)*pgx(1);
		}
		
		/* EDGE 6 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgz(em+m)*pgy(2)*pgx(2);
		}
		
		/* FACE 0 */
		ind2 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += lin(n,++ind)*pgz(p+q+1)*pgy(2+2*em+ind2)*pgx(p+2);
			}
		}
		
		/* FACE 1 */
		sign = 1;
		ind2 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += sign*lin(n,++ind)*pgz(1+2*em+ind2)*pgy(p+2)*pgx(p+2);
			}
			sign*=-1;
		}
		
		/* FACE 2 */
		ind2 = 0;
		sign = 1;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += sign*lin(n,++ind)*pgz(1+2*em+ind2)*pgy(em+2+p)*pgx(2);
			}
			sign*=-1;
		}
		
		/* FACE 3 */
		ind2 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += lin(n,++ind)*pgz(1+2*em+ind2)*pgy(em+2+p)*pgx(1);
			}
		}
	}
		
	return;
}

/* PT PROBE ON FACE 0 must call ptvalues2d before using*/
void tet_basis::ptprobe2d(int nv, FLT *f1, FLT *lin1, int stride) {
	FLT lcl0;
	int ind,ind2,sign;
#ifdef BZ_DEBUG
	Array<FLT,2> lin(lin1, shape(nv,stride), neverDeleteData);
	Array<FLT,1> f(f1, shape(nv), neverDeleteData);
#endif
   
   for(int n = 0; n < nv; ++n) {
   
		lcl0 = 0.0;
		ind = 0;
		
		/* SUM ALL MODE CONTRIBUTIONS */
      
		/* VERTEX 1 */
		lcl0 += lin(n,ind)*pgy(0); //pgx=1
		
		/* VERTEX 2 */
		lcl0 += lin(n,++ind)*pgy(1)*pgx(1);
		
		/* VERTEX 3 */
		lcl0 += lin(n,++ind)*pgy(1)*pgx(2);
		
		/* EDGE 1 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgy(m)*pgx(m+1);
		}
		
		/* EDGE 2 */
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += lin(n,++ind)*pgy(em+m)*pgx(2);
		}
		
		/* EDGE 3 */
		sign = 1;
		for(int m = 2; m <= em+1; ++m) {
			lcl0 += sign*lin(n,++ind)*pgy(em+m)*pgx(1);
			sign*=-1;
		}
						
		/* FACE 0 */
		ind2 = 0;
		for(int p = 1; p <= em-1; ++p) {      
			for(int q = 1; q <= em-p; ++q) {
				++ind2;
				lcl0 += lin(n,++ind)*pgy(1+2*em+ind2)*pgx(p+2);
			}
		}
		
		f(n) = lcl0;
   }
   
   return;
}


/* PT PROBE ON EDGE 1 must call ptvalues1d before using*/
void tet_basis::ptprobe1d(int nv, FLT *f1, FLT *sin1, int stride) {
	const int ltm = 2+em;
	FLT lcl0;
#ifdef BZ_DEBUG
	Array<FLT,2> sin(sin1, shape(nv,stride), neverDeleteData);
	Array<FLT,1> f(f1, shape(nv), neverDeleteData);
#endif

	for(int n=0;n<nv;++n) {
		lcl0 = 0.0;
		for(int k=0; k < ltm; ++k )  
			lcl0 += sin(n,k)*pgx(k);
		f(n) = lcl0;
	}
	return;
}




//void tet_basis::ptprobe_bdry(int nv, FLT *f1, FLT *dx1, FLT *dy1, FLT *dz1, FLT r, FLT s, FLT t, FLT *lin1, int stride) {
//   TinyVector<FLT,MXTM> wk0,wk1,wk2;
//   const int bs1 = em+3, bs2 = 2*em+3, bint = vm+em+fm;
//   FLT lcl0, lcl1, lcl2;
//   FLT xp1,oeta; 
//   FLT x,eta;
//#ifdef BZ_DEBUG
//   Array<FLT,2> lin(lin1, shape(nv,stride), neverDeleteData);
//   Array<FLT,1> f(f1, shape(nv), neverDeleteData);
//   Array<FLT,1> dx(dx1, shape(nv), neverDeleteData);
//   Array<FLT,1> dy(dy1, shape(nv), neverDeleteData);
//#endif
//   
//   s = MIN(1.0-10.*EPSILON,s);
//   x = 2.0*(1+r)/(1-s) -1.0;
//   eta = s;
//   
//   ptvalues_deriv_bdry(x,eta);
//   
//   for(int n=0;n<nv;++n) {
//
//      /* PART I - sum u*g_mn for each n, s_j   */
//      oeta = 2./(1 -eta);
//      
//      /* VERTEX 0         */
//      wk0(0) = lin(n,0)*pgn(0);
//      wk1(0) = lin(n,0)*dpgn(0);
//      wk2(0) = wk0(0)*oeta;
//
//      /* VERTEX 1         */
//      lcl0 = lin(n,1)*pgn(1);
//      lcl1 = lin(n,1)*dpgn(1);
//      /* SIDE 2      */
//      for(int m = bs2; m < bint; ++m ) {
//         lcl0 += lin(n,m)*pgn(m);
//         lcl1 += lin(n,m)*dpgn(m);
//      }  
//      wk0(1) = lcl0;
//      wk1(1) = lcl1;       
//      wk2(1) = lcl0*oeta;
//         
//  
//      /* VERTEX 2         */
//      lcl0 = lin(n,2)*pgn(2);
//      lcl1 = lin(n,2)*dpgn(2);
//      /* SIDE 1      */
//      for (int m = bs1; m < bs2; ++m) {
//         lcl0 += lin(n,m)*pgn(m);
//         lcl1 += lin(n,m)*dpgn(m);
//      }         
//      wk0(2) = lcl0;
//      wk1(2) = lcl1;       
//      wk2(2) = lcl0*oeta;
//
//      for(int m = 3; m < bs1; ++m) {
//         /* SIDE 0      */
//         wk0(m) = lin(n,m)*pgn(m);
//         wk1(m) = lin(n,m)*dpgn(m);
//         wk2(m) = wk0(m)*oeta;
//      }
//
//      /* SUM OVER N AT EACH I,J POINT   */     
//      xp1 = 0.5*(1+x);
//      lcl0 = 0.0;
//      lcl1 = 0.0;
//      lcl2 = 0.0;
//      for(int k=0; k < nmodx; ++k ) {     
//         lcl0 += wk0(k)*pgx(k);
//         lcl1 += wk1(k)*pgx(k);
//         lcl2 += wk2(k)*dpgx(k);
//      }
//      f(n) = lcl0;
//      dx(n) = lcl2;
//      dy(n) = lcl1 +xp1*lcl2;
//   }
//   return;
//}





//void tet_basis::ptprobe1d(int nv, FLT *f1, FLT *dx1, FLT *sin1, int stride) {
//   const int pp = p+1;
//   FLT lcl0,lcl1;
//#ifdef BZ_DEBUG
//   Array<FLT,2> sin(sin1, shape(nv,stride), neverDeleteData);
//   Array<FLT,1> f(f1, shape(nv), neverDeleteData);
//   Array<FLT,1> dx(dx1, shape(nv), neverDeleteData);
//#endif
//   
//   for(int n=0;n<nv;++n) {
//      lcl0 = 0.0;
//      lcl1 = 0.0;
//      for(int k=0; k < pp; ++k ) {
//         lcl0 += sin(n,k)*pgx(k);
//         lcl1 += sin(n,k)*dpgx(k);
//      }
//      f(n) = lcl0;
//      dx(n) = lcl1;
//   }
//   
//   return;
//}
//   
