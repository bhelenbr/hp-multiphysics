/*
 *  ptvalues.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on Tue Apr 16 2002.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_basis.h"
#include<math.h>

/* CALCULATES VALUES OF GX POLYNOMIALS & GS POLYNOMIALS AT POINT */
void tet_basis::ptvalues(FLT x, FLT y, FLT z) {
   FLT pkp,pk,pkm,temp;
   int k,m,n,ind,ind2;
   
	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   
	/******************************************/
	/* GENERATE JACOBI POLY FOR X DIRECTION */
	/****************************************/

	ind = 0;
	
	/* VERTEX 0,1 EDGE 4 */
	pgx(ind)= 1.0;

	/* VERTEX 2, EDGE 3,5, FACE 3 */
	pgx(++ind) = .5*(1-x);


	/* VERTEX 3, EDGE 2,6, FACE 2*/
	pgx(++ind) = .5*(1+x);

	/* EDGE 1 FACE 0,1, INTERIOR */
	/* CALCULATE P, P' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;

	for (m = 0; m < em;++m) {
		pgx(++ind) = 0.25*(1.0+x)*(1.0-x)*pk;		
		pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
		pkm = pk;
		pk = pkp;				
	}
   

	/******************************************/
	/* GENERATE JACOBI POLY FOR Y DIRECTION */
	/****************************************/
	ind = 0;

	/* VERTEX 0 */
	pgy(ind) = 1.0;
	
	/* VERTEX 1, EDGE 4   */
	pgy(++ind) = (1+y)*0.5;

	/* VERTEX 2,3, EDGE 5,6   */
	pgy(++ind) = (1-y)*0.5;

	/* EDGE 1, FACE 1  */
	for(m = 2; m <= em+1; ++m) {
		pgy(++ind) = pow(.5*(1-y),m);  
	}
   
	/* EDGE 2,3 FACE 2,3 */
	pk = 1.0;
	pkm = 0.0;
	for(n = 0; n < em; ++n) {
		pgy(++ind) = (1.0-y)*(1.0+y)*0.25*pk;
		pkp = (y-a0(0,n))*pk - b0(0,n)*pkm;
		pkm = pk;
		pk = pkp;
	}
      
   /* FACE 0, INTERIOR */
   for(m = 1; m <= em-1;++m) {      
      pk = 1.0;
      pkm = 0.0;
	  temp=pow(0.5*(1.0-y),m+1)*0.5*(1.0+y);
      for(k = 1; k <= em-m;++k) {
         pgy(++ind) = temp*pk;
         pkp = (y-a0(m,k-1))*pk - b0(m,k-1)*pkm;
         pkm = pk;
         pk = pkp; 
        
      }
   }
   
   /******************************************/
   /* GENERATE JACOBI POLY FOR Z DIRECTION */
   /****************************************/   
   ind = 0;

   /* VERTEX 0  */
   pgz(ind) = (1.0+z)*0.5;

   /* VERTEX 1,2,3  */
   pgz(++ind) = (1.0-z)*0.5;
   
	/* EDGE 1,2,3, FACE 0 */
	for(m = 2; m <= em+1; ++m) {
		pgz(++ind) = pow(0.5*(1.0-z),m);

	}
		
	/* EDGE 4,5,6  */
    pk = 1.0;
    pkm = 0.0;
    for(n = 0; n < em; ++n){
		pgz(++ind) = (1.0-z)*(1.0+z)*.25*pk;
		pkp = (z-a0(0,n))*pk - b0(0,n)*pkm;
		pkm = pk;
		pk = pkp;
	
	}	
	
	/* FACE 1,2,3  */
	for(m = 1; m <= em-1;++m) {      
		pk = 1.0;
		pkm = 0.0;
		temp=pow(0.5*(1.0-z),m+1)*0.5*(1.0+z);
		for(k = 1; k <= em-m;++k) {
			pgz(++ind) = temp*pk;
			pkp = (z-a0(m,k-1))*pk - b0(m,k-1)*pkm;
			pkm = pk;
			pk = pkp;
		}
	}
	
	/* INTERIOR  */
	ind2=0;
	for(n = 1; n <= em-1; ++n){
		for(m = 1; m <= em-n; ++m) {      
			pk = 1.0;
			pkm = 0.0;
			temp=pow(0.5*(1.0-z),m+n+1)*0.5*(1.0+z);
			for(k = 1; k <= em-m-n;++k) {
				pgz(++ind) = temp*pk;
				pkp = (z-a1(ind2,k-1))*pk - b1(ind2,k-1)*pkm;				
				pkm = pk;
				pk = pkp;
			}
			++ind2;
		}
	}   

   return;
}

/* CALCULATE VALUE OF G(X) & DG/DX, G(eta), DG/Deta AT POINT */
void tet_basis::ptvalues_deriv(FLT x, FLT y, FLT z) {
   FLT pkp,pk,pkm,dpk,dpkm,dpkp;
   int k,m,n,ind,ind2;
   
   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   ind = 0;
   
   	/* VERTEX 0,1 EDGE 4 */
	pgx(ind) = 1.0;
	dpgx(ind) = 0.0;
	
   /* VERTEX 2, EDGE 3,5 FACE 3*/
   pgx(++ind) = .5*(1-x);
   dpgx(ind) = -0.5;
   
   /* VERTEX 3, EDGE 2,6 FACE 2*/
   pgx(++ind) = .5*(1+x);
   dpgx(ind) = 0.5;


	/* EDGE 1 FACE 0,1, INTERIOR */
	/* CALCULATE P, P' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for (m = 0; m < em;++m) {
		pgx(++ind) = 0.25*(1.0+x)*(1.0-x)*pk;
		dpgx(ind) = (-x*0.5*pk +(1.0+x)*(1.0-x)*0.25*dpk);
		pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
		dpkp = pk + (x-a0(0,m))*dpk - b0(0,m)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}
   

	/******************************************/
	/* GENERATE JACOBI POLY FOR Y DIRECTION */
	/****************************************/
	ind = 0;
	/* VERTEX 0  */
	pgy(ind) = 1.0;
	dpgy(ind) = 0.0;
	
	/* VERTEX 1, EDGE 4   */
	pgy(++ind) = (1+y)*.5;
	dpgy(ind) = 0.5;

	/* VERTEX 2,3, EDGE 5,6   */
	pgy(++ind) = (1-y)*.5;
	dpgy(ind) = -0.5;

	/* EDGE 1, FACE 1  */
	for(m = 2; m <= em+1 ; ++m) {
		pgy(++ind) = pow(.5*(1-y),m);
		dpgy(ind) = -0.5*m*pow(.5*(1.-y),m-1);
	}
   
	/* EDGE 2,3 FACE 2,3 */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for(n = 0; n < em; ++n) {
		pgy(++ind) = (1.0-y)*(1.0+y)*.25*pk;
		dpgy(ind) = (-0.5*y*pk +(1.0-y)*(1.0+y)*0.25*dpk);
		pkp = (y-a0(0,n))*pk - b0(0,n)*pkm;
		dpkp = pk + (y-a0(0,n))*dpk - b0(0,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
      
	}
      
   /* FACE 0, INTERIOR */
   for(m = 1; m <= em-1;++m) {      
      pk = 1.0;
      pkm = 0.0;
      dpk = 0.0;
      dpkm = 0.0;

      for(k = 1; k <= em-m;++k) {
         pgy(++ind) = (pow(.5*(1.0-y),m+1)*0.5*(1.0+y)*pk);
         dpgy(ind) = (-.25*(m+1.0)*pow(.5*(1.0-y),m)*(1.0+y)*pk +pow(.5*(1.0-y),m+1)*.5*(pk + (1.0+y)*dpk));
         pkp = (y-a0(m,k-1))*pk - b0(m,k-1)*pkm; 
         dpkp = pk + (y-a0(m,k-1))*dpk - b0(m,k-1)*dpkm;
         dpkm = dpk;
         dpk = dpkp;
         pkm = pk;
         pk = pkp;

      }
   }
   
   /******************************************/
   /* GENERATE JACOBI POLY FOR Z DIRECTION */
   /****************************************/   
   ind = 0;

   /* VERTEX 0  */
   pgz(ind) = (1+z)*.5;
   dpgz(ind) = 0.5;

   /* VERTEX 1,2,3  */
   pgz(++ind) = (1-z)*.5;
   dpgz(ind) = -0.5;
   
	/* EDGE 1,2,3, FACE 0 */
	for(m = 2; m <= em+1; ++m) {
		pgz(++ind) = pow(.5*(1-z),m);
		dpgz(ind) = -0.5*m*pow(.5*(1.-z),m-1);
	}
	
	
	/* EDGE 4,5,6  */
    pk = 1.0;
    pkm = 0.0;
    dpk = 0.0;
    dpkm = 0.0;
    for(n = 0; n < em; ++n) {
		pgz(++ind) = (1.0-z)*(1.0+z)*.25*pk;
		dpgz(ind) = (-.5*z*pk +(1.0-z)*(1.0+z)*0.25*dpk);
		pkp = (z-a0(0,n))*pk - b0(0,n)*pkm;
		dpkp = pk + (z-a0(0,n))*dpk - b0(0,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}	
	
	/* FACE 1,2,3  */
	for(m = 1; m <= em-1;++m) {      
		pk = 1.0;
		pkm = 0.0;
		dpk = 0.0;
		dpkm = 0.0;
		for(k = 1; k <= em-m;++k) {
			pgz(++ind) = (pow(0.5*(1.0-z),m+1)*0.5*(1.0+z)*pk);
			dpgz(ind) = (-0.25*(m+1.0)*pow(0.5*(1.0-z),m)*(1.0+z)*pk +pow(.5*(1.0-z),m+1)*.5*(pk + (1.0+z)*dpk));
			pkp = (z-a0(m,k-1))*pk - b0(m,k-1)*pkm;
			dpkp = pk + (z-a0(m,k-1))*dpk - b0(m,k-1)*dpkm;
			dpkm = dpk;
			dpk = dpkp;
			pkm = pk;
			pk = pkp;
		}
	}
	
	
	/* INTERIOR  */
	ind2 = 0;
	for(n = 1; n <= em-1 ; ++n){
		for(m = 1; m <= em-n; ++m) {      
			pk = 1.0;
			pkm = 0.0;
			dpk = 0.0;
			dpkm = 0.0;
			for(k = 1; k <= em-m-n;++k) {
				pgz(++ind) = (pow(0.5*(1.0-z),m+n+1)*0.5*(1.0+z)*pk);
				dpgz(ind) = (-0.25*(m+n+1)*pow(0.5*(1.0-z),m+n)*(1.0+z)*pk +pow(.5*(1.0-z),m+n+1)*.5*(pk + (1.+z)*dpk));
				pkp = (z-a1(ind2,k-1))*pk - b1(ind2,k-1)*pkm;
				dpkp = pk + (z-a1(ind2,k-1))*dpk - b1(ind2,k-1)*dpkm;	
				dpkm = dpk;
				dpk = dpkp;
				pkm = pk;
				pk = pkp;
			}
			++ind2;
		}
	}
	

   

   return;
}


/* CALCULATES VALUES OF GX POLYNOMIALS & GS POLYNOMIALS AT POINT */
/* BOUNDARY MODES ONLY */
void tet_basis::ptvalues_bdry(FLT x, FLT y, FLT z) {
	FLT pkp,pk,pkm,temp;
	int k,m,n,ind;
   
	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   
	/******************************************/
	/* GENERATE JACOBI POLY FOR X DIRECTION */
	/****************************************/

	ind = 0;
	
	/* VERTEX 0,1 EDGE 4 */
	pgx(ind)= 1.0;

	/* VERTEX 2, EDGE 3,5, FACE 3 */
	pgx(++ind) = .5*(1-x);


	/* VERTEX 3, EDGE 2,6, FACE 2*/
	pgx(++ind) = .5*(1+x);

	/* EDGE 1 FACE 0,1 */
	/* CALCULATE P, P' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;

	for (m = 0; m < em;++m) {
		pgx(++ind) = 0.25*(1.0+x)*(1.0-x)*pk;		
		pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
		pkm = pk;
		pk = pkp;				
	}
   

	/******************************************/
	/* GENERATE JACOBI POLY FOR Y DIRECTION */
	/****************************************/
	ind = 0;

	/* VERTEX 0 */
	pgy(ind) = 1.0;
	
	/* VERTEX 1, EDGE 4   */
	pgy(++ind) = (1+y)*0.5;

	/* VERTEX 2,3, EDGE 5,6   */
	pgy(++ind) = (1-y)*0.5;

	/* EDGE 1, FACE 1  */
	for(m = 2; m <= em+1; ++m) {
		pgy(++ind) = pow(.5*(1-y),m);  
	}
   
	/* EDGE 2,3 FACE 2,3 */
	pk = 1.0;
	pkm = 0.0;
	for(n = 0; n < em; ++n) {
		pgy(++ind) = (1.0-y)*(1.0+y)*0.25*pk;
		pkp = (y-a0(0,n))*pk - b0(0,n)*pkm;
		pkm = pk;
		pk = pkp;
	}
      
   /* FACE 0*/
   for(m = 1; m <= em-1;++m) {      
      pk = 1.0;
      pkm = 0.0;
	  temp=pow(0.5*(1.0-y),m+1)*0.5*(1.0+y);
      for(k = 1; k <= em-m;++k) {
         pgy(++ind) = temp*pk;
         pkp = (y-a0(m,k-1))*pk - b0(m,k-1)*pkm;
         pkm = pk;
         pk = pkp; 
        
      }
   }
   
   /******************************************/
   /* GENERATE JACOBI POLY FOR Z DIRECTION */
   /****************************************/   
   ind = 0;

   /* VERTEX 0  */
   pgz(ind) = (1.0+z)*0.5;

   /* VERTEX 1,2,3  */
   pgz(++ind) = (1.0-z)*0.5;
   
	/* EDGE 1,2,3, FACE 0 */
	for(m = 2; m <= em+1; ++m) {
		pgz(++ind) = pow(0.5*(1.0-z),m);

	}
		
	/* EDGE 4,5,6  */
    pk = 1.0;
    pkm = 0.0;
    for(n = 0; n < em; ++n){
		pgz(++ind) = (1.0-z)*(1.0+z)*.25*pk;
		pkp = (z-a0(0,n))*pk - b0(0,n)*pkm;
		pkm = pk;
		pk = pkp;
	
	}	
	
	/* FACE 1,2,3  */
	for(m = 1; m <= em-1;++m) {      
		pk = 1.0;
		pkm = 0.0;
		temp=pow(0.5*(1.0-z),m+1)*0.5*(1.0+z);
		for(k = 1; k <= em-m;++k) {
			pgz(++ind) = temp*pk;
			pkp = (z-a0(m,k-1))*pk - b0(m,k-1)*pkm;
			pkm = pk;
			pk = pkp;
		}
	}   
	
	return;
}



/* CALCULATE VALUE OF G(X) & DG/DX, G(eta), DG/Deta AT POINT FOR ONLY BOUNDARY MODES */
void tet_basis::ptvalues_deriv_bdry(FLT x, FLT y, FLT z) {
 FLT pkp,pk,pkm,dpk,dpkm,dpkp;
   int k,m,n,ind;
   
   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   ind = 0;
   
   	/* VERTEX 0,1 EDGE 4 */
	pgx(ind) = 1.0;
	dpgx(ind) = 0.0;
	
   /* VERTEX 2, EDGE 3,5 FACE 3*/
   pgx(++ind) = .5*(1-x);
   dpgx(ind) = -0.5;
   
   /* VERTEX 3, EDGE 2,6 FACE 2*/
   pgx(++ind) = .5*(1+x);
   dpgx(ind) = 0.5;


	/* EDGE 1 FACE 0,1 */
	/* CALCULATE P, P' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for (m = 0; m < em;++m) {
		pgx(++ind) = 0.25*(1.0+x)*(1.0-x)*pk;
		dpgx(ind) = (-x*0.5*pk +(1.0+x)*(1.0-x)*0.25*dpk);
		pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
		dpkp = pk + (x-a0(0,m))*dpk - b0(0,m)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}   

	/******************************************/
	/* GENERATE JACOBI POLY FOR Y DIRECTION */
	/****************************************/
	ind = 0;
	/* VERTEX 0  */
	pgy(ind) = 1.0;
	dpgy(ind) = 0.0;
	
	/* VERTEX 1, EDGE 4   */
	pgy(++ind) = (1+y)*.5;
	dpgy(ind) = 0.5;

	/* VERTEX 2,3, EDGE 5,6   */
	pgy(++ind) = (1-y)*.5;
	dpgy(ind) = -0.5;

	/* EDGE 1, FACE 1  */
	for(m = 2; m <= em+1 ; ++m) {
		pgy(++ind) = pow(.5*(1-y),m);
		dpgy(ind) = -0.5*m*pow(.5*(1.-y),m-1);
	}
   
	/* EDGE 2,3 FACE 2,3 */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for(n = 0; n < em; ++n) {
		pgy(++ind) = (1.0-y)*(1.0+y)*.25*pk;
		dpgy(ind) = (-0.5*y*pk +(1.0-y)*(1.0+y)*0.25*dpk);
		pkp = (y-a0(0,n))*pk - b0(0,n)*pkm;
		dpkp = pk + (y-a0(0,n))*dpk - b0(0,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
      
	}
      
   /* FACE 0 */
   for(m = 1; m <= em-1;++m) {      
      pk = 1.0;
      pkm = 0.0;
      dpk = 0.0;
      dpkm = 0.0;

      for(k = 1; k <= em-m;++k) {
         pgy(++ind) = (pow(.5*(1.0-y),m+1)*0.5*(1.0+y)*pk);
         dpgy(ind) = (-.25*(m+1.0)*pow(.5*(1.0-y),m)*(1.0+y)*pk +pow(.5*(1.0-y),m+1)*.5*(pk + (1.0+y)*dpk));
         pkp = (y-a0(m,k-1))*pk - b0(m,k-1)*pkm; 
         dpkp = pk + (y-a0(m,k-1))*dpk - b0(m,k-1)*dpkm;
         dpkm = dpk;
         dpk = dpkp;
         pkm = pk;
         pk = pkp;

      }
   }
   
   /******************************************/
   /* GENERATE JACOBI POLY FOR Z DIRECTION */
   /****************************************/   
   ind = 0;

   /* VERTEX 0  */
   pgz(ind) = (1+z)*.5;
   dpgz(ind) = 0.5;

   /* VERTEX 1,2,3  */
   pgz(++ind) = (1-z)*.5;
   dpgz(ind) = -0.5;
   
	/* EDGE 1,2,3, FACE 0 */
	for(m = 2; m <= em+1; ++m) {
		pgz(++ind) = pow(.5*(1-z),m);
		dpgz(ind) = -0.5*m*pow(.5*(1.-z),m-1);
	}
	
	
	/* EDGE 4,5,6  */
    pk = 1.0;
    pkm = 0.0;
    dpk = 0.0;
    dpkm = 0.0;
    for(n = 0; n < em; ++n) {
		pgz(++ind) = (1.0-z)*(1.0+z)*.25*pk;
		dpgz(ind) = (-.5*z*pk +(1.0-z)*(1.0+z)*0.25*dpk);
		pkp = (z-a0(0,n))*pk - b0(0,n)*pkm;
		dpkp = pk + (z-a0(0,n))*dpk - b0(0,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}	
	
	/* FACE 1,2,3  */
	for(m = 1; m <= em-1;++m) {      
		pk = 1.0;
		pkm = 0.0;
		dpk = 0.0;
		dpkm = 0.0;
		for(k = 1; k <= em-m;++k) {
			pgz(++ind) = (pow(0.5*(1.0-z),m+1)*0.5*(1.0+z)*pk);
			dpgz(ind) = (-0.25*(m+1.0)*pow(0.5*(1.0-z),m)*(1.0+z)*pk +pow(.5*(1.0-z),m+1)*.5*(pk + (1.0+z)*dpk));
			pkp = (z-a0(m,k-1))*pk - b0(m,k-1)*pkm;
			dpkp = pk + (z-a0(m,k-1))*dpk - b0(m,k-1)*dpkm;
			dpkm = dpk;
			dpk = dpkp;
			pkm = pk;
			pk = pkp;
		}
	}
	
	return;
}


void tet_basis::ptvalues2d(FLT x, FLT y) {
	FLT pkp,pk,pkm,temp;
	int k,m,n,ind;

	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	/******************************************/
	/* GENERATE JACOBI POLY FOR X DIRECTION */
	/****************************************/

	ind = 0;

	/* VERTEX 1 */
	pgx(ind)= 1.0;

	/* VERTEX 2, EDGE 3 */
	pgx(++ind) = .5*(1-x);


	/* VERTEX 3, EDGE 2*/
	pgx(++ind) = .5*(1+x);

	/* EDGE 1 FACE 0 */
	/* CALCULATE P, P' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;

	for (m = 0; m < em;++m) {
		pgx(++ind) = 0.25*(1.0+x)*(1.0-x)*pk;		
		pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
		pkm = pk;
		pk = pkp;				
	}

	/******************************************/
	/* GENERATE JACOBI POLY FOR Y DIRECTION */
	/****************************************/
	ind = 0;

	/* VERTEX 1   */
	pgy(ind) = (1+y)*0.5;

	/* VERTEX 2,3 */
	pgy(++ind) = (1-y)*0.5;

	/* EDGE 1 */
	for(m = 2; m <= em+1; ++m) {
		pgy(++ind) = pow(.5*(1-y),m);  
	}

	/* EDGE 2,3 */
	pk = 1.0;
	pkm = 0.0;
	for(n = 0; n < em; ++n) {
		pgy(++ind) = (1.0-y)*(1.0+y)*0.25*pk;
		pkp = (y-a0(0,n))*pk - b0(0,n)*pkm;
		pkm = pk;
		pk = pkp;
	}
	  
	/* FACE 0*/
	for(m = 1; m <= em-1;++m) {      
	  pk = 1.0;
	  pkm = 0.0;
	  temp=pow(0.5*(1.0-y),m+1)*0.5*(1.0+y);
	  for(k = 1; k <= em-m;++k) {
		 pgy(++ind) = temp*pk;
		 pkp = (y-a0(m,k-1))*pk - b0(m,k-1)*pkm;
		 pkm = pk;
		 pk = pkp; 
		
	  }
	}

	return;
}

void tet_basis::ptvalues2d_deriv(FLT x, FLT y){
	FLT pkp,pk,pkm,dpk,dpkm,dpkp;
	int k,m,n,ind;

	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	ind = 0;

	/* VERTEX 1 */
	pgx(ind) = 1.0;
	dpgx(ind) = 0.0;

	/* VERTEX 2, EDGE 3*/
	pgx(++ind) = .5*(1-x);
	dpgx(ind) = -0.5;

	/* VERTEX 3, EDGE 2*/
	pgx(++ind) = .5*(1+x);
	dpgx(ind) = 0.5;


	/* EDGE 1 FACE 0 */
	/* CALCULATE P, P' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for (m = 0; m < em;++m) {
		pgx(++ind) = 0.25*(1.0+x)*(1.0-x)*pk;
		dpgx(ind) = (-x*0.5*pk +(1.0+x)*(1.0-x)*0.25*dpk);
		pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
		dpkp = pk + (x-a0(0,m))*dpk - b0(0,m)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}   

	/******************************************/
	/* GENERATE JACOBI POLY FOR Y DIRECTION */
	/****************************************/
	ind = 0;	

	/* VERTEX 1  */
	pgy(ind) = (1+y)*.5;
	dpgy(ind) = 0.5;

	/* VERTEX 2,3 */
	pgy(++ind) = (1-y)*.5;
	dpgy(ind) = -0.5;

	/* EDGE 1  */
	for(m = 2; m <= em+1 ; ++m) {
		pgy(++ind) = pow(.5*(1-y),m);
		dpgy(ind) = -0.5*m*pow(.5*(1.-y),m-1);
	}

	/* EDGE 2,3 */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for(n = 0; n < em; ++n) {
		pgy(++ind) = (1.0-y)*(1.0+y)*.25*pk;
		dpgy(ind) = (-0.5*y*pk +(1.0-y)*(1.0+y)*0.25*dpk);
		pkp = (y-a0(0,n))*pk - b0(0,n)*pkm;
		dpkp = pk + (y-a0(0,n))*dpk - b0(0,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;

	}
	  
	/* FACE 0 */
	for(m = 1; m <= em-1;++m) {      
		pk = 1.0;
		pkm = 0.0;
		dpk = 0.0;
		dpkm = 0.0;

		for(k = 1; k <= em-m;++k) {
			pgy(++ind) = (pow(.5*(1.0-y),m+1)*0.5*(1.0+y)*pk);
			dpgy(ind) = (-.25*(m+1.0)*pow(.5*(1.0-y),m)*(1.0+y)*pk +pow(.5*(1.0-y),m+1)*.5*(pk + (1.0+y)*dpk));
			pkp = (y-a0(m,k-1))*pk - b0(m,k-1)*pkm; 
			dpkp = pk + (y-a0(m,k-1))*dpk - b0(m,k-1)*dpkm;
			dpkm = dpk;
			dpk = dpkp;
			pkm = pk;
			pk = pkp;
		}
	}

	return;
}


void tet_basis::ptvalues1d(FLT x) {
	FLT pkp,pk,pkm;
	int m,ind;

	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	/******************************************/
	/* GENERATE JACOBI POLY FOR X DIRECTION */
	/****************************************/

	ind = 0;
	/* VERTEX 2 */
	pgx(ind) = .5*(1-x);

	/* VERTEX 3*/
	pgx(++ind) = .5*(1+x);

	/* EDGE 1  */
	/* CALCULATE P, P' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;

	for (m = 0; m < em;++m) {
		pgx(++ind) = 0.25*(1.0+x)*(1.0-x)*pk;		
		pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
		pkm = pk;
		pk = pkp;				
	}
   
    return;
}

void tet_basis::ptvalues1d_deriv(FLT x) {
	FLT pkp,pk,pkm,dpk,dpkm,dpkp;
	int m,ind;

	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	ind = 0;
	
	/* VERTEX 2*/
	pgx(ind) = .5*(1-x);
	dpgx(ind) = -0.5;

	/* VERTEX 3*/
	pgx(++ind) = .5*(1+x);
	dpgx(ind) = 0.5;


	/* EDGE 1 */
	/* CALCULATE P, P' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for (m = 0; m < em;++m) {
		pgx(++ind) = 0.25*(1.0+x)*(1.0-x)*pk;
		dpgx(ind) = (-x*0.5*pk +(1.0+x)*(1.0-x)*0.25*dpk);
		pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
		dpkp = pk + (x-a0(0,m))*dpk - b0(0,m)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}   
   
   return;
}
