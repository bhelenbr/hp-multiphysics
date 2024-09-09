/*
 *  initialize.cpp
 * 
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#define NO_DEBUG

#include <math.h>
#include <myblas.h>
#include "tet_basis.h"

const int tet_basis::sbwth;
using namespace blitz;

/*comment this for nodal 3/4*/
Array<tet_basis,1> basis::tet;

void tet_basis::initialize(int pdegree, int gpoints) { 
#ifdef BZ_DEBUG
    std::cerr << "#spline: BZ_DEBUG is set\n";
#endif
#ifdef DEBUG
    std::cerr << "#spline: Running in Xcode's DEBUG Mode\n";
#endif
	if (pdegree < 1) {
		printf("error can't use 0th order basis with vertex based modes\n");
		exit(1);
	}
   
	p = pdegree;
	vm = 4;
	em = p-1;
	fm = MAX((p-1)*(p-2)/2,0);
	tm = (p+1)*(p+2)*(p+3)/6;
	im = MAX(tm-4*fm-6*em-vm,0);
	bm = vm+6*em+4*fm;
	ibwth = (2*(em-1) < im -1? 2*(em-1) : im-1);//fix
	ibwth = (ibwth > 0 ? ibwth : 0);
#ifdef MORTHOGONAL
	ibwth = 1;
#endif
   
	nmodx = em+3;
	nmody = 3+2*em+fm;
	nmodz = 2+2*em+fm+im;
	gpx = gpoints;
	gpy = gpoints;
	gpz = gpoints;

	pgx.resize(nmodx);
	dpgx.resize(nmodx);
	pgy.resize(nmody);
	dpgy.resize(nmody);
	pgz.resize(nmodz);
	dpgz.resize(nmodz);
	
   
   /*****************************/
   /* SETUP VALUES OF FUNCTIONS */
   /*****************************/
	initialize_values(); // SET UP THINGS FOR PROJECT/INTEGRATE/DERIV
//   faceinfoinit(); // SET UP THINGS TO CALCULATE NORMAL DERIVATIVE TO SIDE ALONG SIDE
    lumpinv(); // SET UP THINGS FOR INVERSE OF LUMPED MASS MATRIX
	legpt(); // SET UP PROJECTION TO LEGENDRE POINTS (FOR OUTPUTING)
 
   
   return;
}

void tet_basis::initialize_values(void){

	FLT al,be,x,y,z;
	int i,j,k,m,n,ipoly,ierr,temp,ind;
	Array<FLT,1> e(MXTM),e1(MXTM),e2(MXTM);

	/* ALLOCATE STORAGE FOR RECURSION RELATION COEFFICENTS*/
	temp = MAX(gpx+1,em+1);
	a0.resize(em+2,temp);
	b0.resize(em+2,temp);
	a1.resize(fm+1,temp);
	b1.resize(fm+1,temp);
	
	/* ALLOCATE INTEGRATION, PROJECTION, & DERIVATIVE VARIABLES */
	gx.resize(gpx,nmodx);
	dgx.resize(gpx,nmodx);
	wtx.resize(gpx);
	xp.resize(gpx);
	x0.resize(gpx);
	gxwtx.resize(nmodx,gpx);
	dgxwtx.resize(nmodx,gpx);
	dltx.resize(gpx);
	dltx2.resize(gpx,gpx);
	dltx3.resize(gpx,gpx);

	gy.resize(gpy,nmody);
	dgy.resize(gpy,nmody);
	wty.resize(gpy);
	yp.resize(gpy);
	y0.resize(gpy);
	y1.resize(gpy);
	gywty.resize(tm,gpy);
	gywtyy0.resize(tm,gpy);
	dgywty.resize(tm,gpy);
	dlty.resize(gpy);
	dlty2.resize(gpy,gpy);
	dlty3.resize(gpy,gpy);
		
	gz.resize(gpz,nmodz);
	dgz.resize(gpz,nmodz);
	wtz.resize(gpz);
	zp.resize(gpz);
	z0.resize(gpz);
	gzwtz.resize(tm,gpz);
	gzwtzz0.resize(tm,gpz);
	dgzwtz.resize(tm,gpz);
	dltz.resize(gpz);
	dltz2.resize(gpz,gpz);
	
	norm.resize(tm);
	
	norm=1.0; // delete later
      
	/* GENERATE RECURSION RELATION FOR LEGENDRE
	POLYNOMIALS (USED TO GENERATE GAUSS POINTS)
	IPOLY = 1 SIGNIFIES LEGENDRE POLYNOMIALS
	NEED RECURSION RELATIONS FROM 1 TO N FOR
	N POINT GAUSS FORMULA

	RECURSION RELATION IS OF THE FORM:
	p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
	k=0,1,...,n-1,

	p(-1)(x)=0,  p(0)(x)=1
	*/

	/******************************************/
	/* CALCULATE GAUSS POINTS / COEFFICIENTS  */
	/**************************************   */
	ipoly = 1;
	al = 0.0;
	be = 0.0;
	ierr = recur(gpx+1,ipoly,al,be,&a0(0,0),&b0(0,0));
	if (ierr != 0) {
		printf("recur #1 error %d\n",ierr);
		exit(1);
	}
	
	
	ierr = gauss(gpx,&a0(0,0),&b0(0,0),EPSILON,&x0(0),&wtx(0),&e(0));
	if (ierr != 0) {
		printf("gauss #1 error %d\n",ierr);
		exit(1);
	}

	/***********************************************/
	/* CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
	/***********************************************/
	for (i=0;i<gpx;++i) {
		dltx(i) = 1.0;
		for(k=0;k<i;++k)
			dltx(i) *= (x0(i)-x0(k));
		for(k=i+1;k<gpx;++k)
			dltx(i) *= (x0(i)-x0(k));
	}
   
	for (i=0;i<gpx;++i) {
		for(n=0;n<i;++n){
			dltx2(i,n) = dltx(i)/(x0(i)-x0(n));
			dltx3(i,n) = dltx2(i,n)*(1.+x0(i))/2.;
		}
		for(n=i+1;n<gpx;++n){
			dltx2(i,n) = dltx(i)/(x0(i)-x0(n));
			dltx3(i,n) = dltx2(i,n)*(1.+x0(i))/2.;
		}
	}
      
	for (i=0;i<gpx;++i)
		dltx(i) = 1.0/dltx(i);   

	/*************************************************/ 
	/* NOW CALCULATE VALUES OF G, G' AT GAUSS POINTS */
	/*************************************************/

	/* SIDE 1 IPOLY = 6 FOR JACOBI POLYNOMIALS  */
	ipoly = 6;
	al = 2.0;
	be = 2.0;
	//cout << em+1 << endl;
	ierr = recur(em+1,ipoly,al,be,&a0(0,0),&b0(0,0)); 
	if (ierr != 0) {
		printf("recur #3 error %d\n",ierr);
		exit(1);
	}

	for(i=0;i<tm;++i)
		norm(i) = 1.0;
   
	for(i = 0;i < gpx; ++i) {
		x = x0(i);
		xp(i) = x;
		x0(i) = 0.5*(1+x);
      
		ptvalues_deriv(x,0.0,0.0);

		for (m = 0;m < nmodx;++m) {
			gx(i,m) = pgx(m);
			dgx(i,m) = dpgx(m);
		}
	}

	/*********************************************/
	/* CALCULATE GAUSS POINTS / COEFFICIENTS (Y) */
	/****************************************   */
	ipoly = 6;
	al = 1.0;
	be = 0.0;

   ierr = recur(gpy+1,ipoly,al,be,&a0(1,0),&b0(1,0));
   if (ierr != 0) {
      printf("recur #2 error %d\n",ierr);
      exit(1);
   }

   ierr = gauss(gpy,&a0(1,0),&b0(1,0),EPSILON,&y0(0),&wty(0),&e(0));
   if (ierr != 0) {
      printf("gauss #1 error %d\n",ierr);
      exit(1);
   }
   
   for (i=0;i<gpy;++i)
      wty(i) *= 0.5;

   /***********************************************/
   /* CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
   /***********************************************/
   for (i=0;i<gpy;++i) {
      dlty(i) = 1.0;
      for(k=0;k<i;++k)
         dlty(i) *= (y0(i)-y0(k));
      for(k=i+1;k<gpy;++k)
         dlty(i) *= (y0(i)-y0(k));
   }
   
	for (j=0;j<gpy;++j){
		for(n=0;n<j;++n) {
			dlty2(j,n) = dlty(j)/(y0(j)-y0(n));
			dlty3(j,n) = dlty2(j,n)*(1.+y0(j))/2.;
		}
		for(n=j+1;n<gpy;++n){ 
			dlty2(j,n) = dlty(j)/(y0(j)-y0(n));
			dlty3(j,n) = dlty2(j,n)*(1.+y0(j))/2.;
		}
	}
   
   for (j=0;j<gpy;++j)
      dlty(j) = 1.0/dlty(j);
   

   /******************************************/
   /* GENERATE JACOBI POLY FOR Y DIRECTION */
   /****************************************/

   /*	RECURSION RELATION FOR FACE MODES */
   for(m = 1; m <= em;++m) { //note changed m=2 -> m=1  
		/* CALCULATE RECURSION RELATION FOR P^(2m-3,2)(s) SIDE 1 IPOLY = 6 FOR JACOBI     */
		ipoly = 6;
		al = 2.0*m+3.0;
		be = 2.0;

		ierr = recur(em+1,ipoly,al,be,&a0(m,0),&b0(m,0));
		if (ierr != 0) {
			 printf("recur #4 error %d\n",ierr);
			 exit(1);
		}
   }
      
   for(j=0;j<gpy;++j) {
      y = y0(j);
      yp(j) = y;
      y0(j) = 2.0/(1-y);
	  y1(j) = 0.5*(1+y);

      ptvalues_deriv(0.0,y,0.0);
      
      for(m=0;m<nmody;++m) {
         gy(j,m) = pgy(m);
         dgy(j,m) = dpgy(m);
      }
   }
   
   	/*********************************************/
	/* CALCULATE GAUSS POINTS / COEFFICIENTS (Z) */
	/****************************************   */
	ipoly = 6;
	al = 2.0;
	be = 0.0;

   ierr = recur(gpz+1,ipoly,al,be,&a1(0,0),&b1(0,0));
   if (ierr != 0) {
      printf("recur #2 error %d\n",ierr);
      exit(1);
   }

   ierr = gauss(gpz,&a1(0,0),&b1(0,0),EPSILON,&z0(0),&wtz(0),&e(0));
   if (ierr != 0) {
      printf("gauss #1 error %d\n",ierr);
      exit(1);
   }
   
   for (i=0;i<gpz;++i)
      wtz(i) *= 0.25;

   /***********************************************/
   /* CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
    /***********************************************/  
   
   for (i=0;i<gpz;++i) {
      dltz(i) = 1.0;
      for(k=0;k<i;++k)
         dltz(i) *= (z0(i)-z0(k));
      for(k=i+1;k<gpz;++k)
         dltz(i) *= (z0(i)-z0(k));
   }   
   
   for (j=0;j<gpz;++j) {
      for(n=0;n<j;++n) 
         dltz2(j,n) = dltz(j)/(z0(j)-z0(n));
      for(n=j+1;n<gpz;++n) 
         dltz2(j,n) = dltz(j)/(z0(j)-z0(n));
   }      
   
   for (j=0;j<gpz;++j)
      dltz(j) = 1.0/dltz(j);
   
	/******************************************/
	/* GENERATE JACOBI POLY FOR Z DIRECTION */
	/****************************************/

	/*	RECURSION RELATION FOR INTERIOR MODES */
	ind = 0;
	for(n = 1; n <= em-1; ++n) { 
		for(m = 1; m <= em-n; ++m) {
 				ipoly = 6;
				al = 2.0*n+2.0*m+4.0;
				be = 2.0;
				
				ierr = recur(em+1,ipoly,al,be,&a1(ind,0),&b1(ind,0));
				if (ierr != 0) {
					 printf("recur #4 error %d\n",ierr);
					 exit(1);
				 }
			++ind;
		}
	}   
   
   for(j=0;j<gpz;++j) {
	  z = z0(j);
	  zp(j) = z;
	  z0(j) = 2.0/(1-z);

	  ptvalues_deriv(0.0,0.0,z);
	  
	  for(m = 0; m < nmodz; ++m) {
		 gz(j,m) = pgz(m);
		 dgz(j,m) = dpgz(m);
	  }
   }

//   /******************************/
//   /* CALCULATE NORM             */
//   /* ************************** */
//   /* SIDE & VERTEX MODES */
//   for(n=0;n<3;++n)
//      norm(n) = 1.0;
//   for (n = 3; n < em+3; ++n) {
//      norm(n) = 0.0;
//      for(i = 0; i < gpx; ++i)
//         norm(n) += wtx(i)*gx(i,n)*gx(i,n); //l2 norm
//      norm(n) = 1./sqrt(norm(n));
//      norm(n+em) = norm(n);
//      norm(n+2*em) = norm(n);
//   }
//   
//   /* INTERIOR MODES */
//   for(m = bm; m < tm; ++m) {
//      norm(m) = 0.0;
//      for(j = 0; j < gpy; ++j)
//           norm(m) += wty(j)*gy(j,m)*gy(j,m);
//        norm(m) = 1./sqrt(norm(m));
//   }
//
//   /***************/
//   /* RENORMALIZE */
//   /***************/
//   for (n = 3; n < nmodx; ++n) {         
//      for (i =0;i<gpx;++i) {
//         gx(i,n) *= norm(n);
//         dgx(i,n) *= norm(n);
//#ifdef DEBUG
//         printf("IV1: %d %d %e %e\n",i,n,gx(i,n),dgx(i,n));
//#endif
//      }
//      for (j=0;j<gpy;++j) {
//         /* SIDE 2 & 3 MUST BE RENORMALIZED BY SAME CONSTANT TO MATCH */
//         gy(j,n +em) *= norm(n);
//         dgy(j,n+em) *= norm(n);
//     //    gy(j,n +2*em) *= norm(n);
//     //    dgy(j,n +2*em) *= norm(n);
//#ifdef DEBUG
//         printf("IV2: %d %d %e %e\n",j,n,gy(j,n),dgx(j,n));
//#endif
//      }
//   }
//
//   for(m = bm; m < tm; ++m) {
//      for(j = 0; j < gpy; ++j) {
//         gy(j,m) = gy(j,m)*norm(m);
//         dgy(j,m) = dgy(j,m)*norm(m);
//#ifdef DEBUG
//         printf("IV3: %d %d %e %e\n",j,m,gn(j,m),dgn(j,m));
//#endif
//      }
//   }

	/*****************************************/
	/* PRECALCULATE THINGS FOR FAST INTGRTRS */
	/*****************************************/
	for(m=0;m<nmodx;++m) {
		for(i=0;i<gpx;++i) {
			gxwtx(m,i) = gx(i,m)*wtx(i);
			dgxwtx(m,i) = dgx(i,m)*wtx(i);
#ifdef DEBUG
		printf("%d %d %e %e\n",m,i,gxwtx(m,i),dgxwtx(m,i));
#endif
		}
	}
   
	for(m=0;m<nmody;++m) {
		for(j=0;j<gpy;++j) {
			gywty(m,j) = gy(j,m)*wty(j);
			dgywty(m,j) = dgy(j,m)*wty(j);
#ifdef DEBUG
		printf("%d %d %e %e \n",m,j,gzwty(m,j),dgzwty(m,j));
#endif
		}
	}
   
	for(m=0;m<nmodz;++m) {
		for(j=0;j<gpz;++j) {
			gzwtz(m,j) = gz(j,m)*wtz(j);
			dgzwtz(m,j) = dgz(j,m)*wtz(j);
#ifdef DEBUG
		printf("%d %d %e %e \n",m,j,gzwtz(m,j),dgzwtz(m,j));
#endif
		}
	}
   
	return;
}

      
//void tet_basis::faceinfoinit() { //fix
//	int i,m,n,ind;
//	FLT x,y,z,eta,xp1oeta,xp1,oeta;
//
//	/*	THIS IS TO CALCULATE NORMAL DERIVATIVES TO SIDE */
//	/* SIDES 1 & 2 ARE ROTATED TO SIDE 0 POSITION */
//	dgnorm.resize(tm,gpx,gpy);
//   
//   /*	CALCULATE NORMAL DERIVATIVE VALUES ALONG SIDES */
//   /*	FACE 0 */
//	for(int i = 0; i < gpx; ++i) {
//		for(int j = 0; j < gpy; ++j){
//			z = -1.0;
//			x = xp(i);
//			y = yp(j);
//			ptvalues_deriv(x,y,z);
//			xp1oeta = x0(i)*2.0/(1-eta);
//
//			/* CALCULATE POLYNOMIALS */
//			/* VERTEX 0   */
//			dgnorm(0,i,j) = dpgy(0)*pgx(0) +xp1oeta*pgy(0)*dpgx(0);
//
//			/* VERTEX 1  */
//			dgnorm(1,i,j) = dpgy(1)*pgx(1) +xp1oeta*pgy(1)*dpgx(1);
//
//			/* VERTEX 2    */   
//			dgnorm(2,i,j) = dpgy(2)*pgx(2) +xp1oeta*pgy(2)*dpgx(2);
//
//			for(m = 3; m < em+3; ++m)
//			 dgnorm(m,i,j) = dpgy(m)*pgx(m) +xp1oeta*pgy(m)*dpgx(m);
//			 
//			for(m=em+3;m<2*em+3;++m)
//			 dgnorm(m,i,j) = dpgy(m)*pgx(2) +xp1oeta*pgy(m)*dpgx(2);
//
//			for(m=2*em+3;m<bm;++m)
//			 dgnorm(m,i,j) = dpgy(m)*pgx(1) +xp1oeta*pgy(m)*dpgx(1);
//
//			/*  INTERIOR MODES   */
//			ind = bm;
//			for(m = 3; m< em+2;++m) {      
//				for(n = 1; n < em+3-m;++n) {
//					dgnorm(ind,i,j) = dpgy(ind)*pgx(m) +xp1oeta*pgy(ind)*dpgx(m);
//					++ind;
//				}
//			}
//		}
//	}
//     
//   return;
//}



/************************************************/
/** CALCULATE THINGS FOR LUMPED MASS INVERSION  */
/************************************************/
void tet_basis::lumpinv(void) {

   int info,ind,ind2;
   int ltm=(em+1)*(em+2)*(em+3)/6;
   int vdof = 1+3*em+3*fm+im;
   int lumped = 1;
   Array<int,1> ipiv(MXTM);
   Array<FLT,2> mwk(MXTM,MXTM);//, mm(MXTM,MXTM);
   Array<FLT,1> u(MXTM),l(MXTM),vwk(MXTM),wk1(MXTM*MXGP);

   Array<int,1> lind(ltm);//lower order mode constraints
   Array<int,2> dind(6,vdof);//dofs
   TinyVector<TinyVector<TinyVector<int,2>,3>,6> edgecone;
   char trans[] = "T";
   
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	Array<FLT,3> wk2(MXGP,MXGP,MXGP);
//   /* ALLOCATE MASS MATRIX INVERSION VARIABLES */
    if (em > 0) {
		vfms.resize(14,tm);//vertex from sides
		sfmv.resize(2,em);
		ediag.resize(em);
		diag1d.resize(em);
    }
	if (fm > 0) {
		fdiag.resize(fm);
		diag2d.resize(fm);
		ffmv.resize(3,fm);
		sfms.resize(3);//temp
		ffms.resize(3,fm,2);//temp
	}
    if (im > 0) {
		diag3d.resize(im);
		idiag.resize(im);
		ifmb.resize(bm,im);
    }
	
	mm.resize(MXTM,MXTM);
	
   
//   if (em > 1) sfms.resize(em-1,em,3);
//   if (im > 0) {
//      ifmb.resize(bm,im);
//      bfmi.resize(bm,im);
//      idiag.resize(im,ibwth+1);
//   }
//   msi.resize(bm,bm);
//   
//   /* ALLOCATE 1D MASS MATRIX INVERSION VARIABLES */
//   if (em > 0) {
//      vfms1d.resize(2,em);
//      sfmv1d.resize(2,em);
//      sdiag1d.resize(em,sbwth+1);
//   }
//
//   /********************************************/   
//   /* GENERATE MASS MATRIX                  */
//   /********************************************/
	
	
	mm = 0;
	// 1d mass matrix along edge
	int tm1d = 2+em;
	for(int m=0;m<tm1d;++m) {
		for(int i=0;i<tm1d;++i)
			u(i) = 0.0;
		u(m) = 1.0;

		proj1d(&u(0),&wk2(0,0,0)); 
		intgrt1d(l.data(),wk2.data()); 
         
		for(int i=0;i<tm1d;++i) 
			mm(m,i) = l(i);      
	}
	for(int i = 0; i < em; ++i){
		diag1d(i) = 1.0/mm(i+2,i+2);
	}
	
	
	// 2d mass matrix along face
	int tm2d = 3+3*em+fm;
	for(int m=0;m<tm2d;++m) {
		for(int i=0;i<tm2d;++i)
			u(i) = 0.0;
		u(m) = 1.0;

		proj2d(&u(0),&wk2(0,0,0),stridey); 
		intgrt2d(l.data(),wk2.data(),stridey); 
         
		for(int i=0;i<tm2d;++i) 
			mm(m,i) = l(i);      
	}
	
	int bfm = 3+3*em;
	for(int i = 0; i < fm; ++i){
		diag2d(i) = 1.0/mm(i+bfm,i+bfm);
	}
	
	// full 3d mass matrix
	for(int m=0;m<tm;++m) {
		for(int i=0;i<tm;++i)
			u(i) = 0.0;
		u(m) = 1.0;

		proj(&u(0),&wk2(0,0,0),stridex,stridey); 
		intgrt(l.data(),wk2.data(),stridex,stridey); 
         
		for(int i=0;i<tm;++i) 
			mm(m,i) = l(i);      
	}	

	for(int i = 0; i < im; ++i)
		diag3d(i) = 1.0/mm(i+bm,i+bm);
	
	if(lumped)
		vdiag = 1/(mm(0,0)+mm(0,1)+mm(0,2)+mm(0,3));
	else
		vdiag = 1/mm(0,0);

	
	if(p < 2) return;	
		
	
   /*******************************************************/      
   /*  EQUATIONS TO FIND VERTEX VALUES TO em-1 ACCURACY */
   /*******************************************************/      
	// find lower order modes for constraints
		for(int i = 0; i < 4; ++i) 
			lind(i)=i;

		ind = 3;
	
		for(int m = 0; m < 6; ++m)
			for(int j=4+em*m; j < 3+em*(m+1); ++j)
				lind(++ind)=j;

		for(int m = 0; m < 4; ++m){
			ind2 = 0;
			for(int i = 1; i <= em-1; ++i){
				for(int j = 1; j <= em-i-1; ++j){
					lind(++ind)=4+6*em+fm*m+ind2;
					++ind2;
				}
				++ind2;
			}
		}
		
		ind2=0;
		for(int i = 1; i <= em-1; ++i){
			for(int j = 1; j <= em-i-1; ++j){
				for(int k = 1; k <= em-i-j-1; ++k){
					lind(++ind)=bm+ind2;
					++ind2;
				}
				++ind2;
			}
			++ind2;
		}
		
		
		// find all degrees of freedom for vertices
		// vertex 0
		ind = 0;
		dind(0,ind) = 0;
		for(int i = 4+3*em; i < 4+6*em; ++i)
			dind(0,++ind)=i;// edge 4-6
			
		for(int i = 4+6*em+fm; i < tm; ++i)
			dind(0,++ind)=i; // face 1-interior
			
		// vertex 1
		ind = 0;
		dind(1,ind) = 1;
		for(int i = 4+em; i < 4+4*em; ++i)
			dind(1,++ind)=i; // edge 2-4
			
		for(int i = 4+6*em; i < 4+6*em+fm; ++i)
			dind(1,++ind)=i;// face 0
		
		for(int i = 4+6*em+2*fm; i < tm; ++i)
			dind(1,++ind)=i; // face 2-interior			

		
		// vertex 2
		ind = 0;
		dind(2,ind) = 2;
		for(int i = 4; i < 4+em; ++i)
			dind(2,++ind)=i;// edge 1
			
		for(int i = 4+2*em; i < 4+3*em; ++i)
			dind(2,++ind)=i;// edge 3
			
		for(int i = 4+4*em; i < 4+5*em; ++i)
			dind(2,++ind)=i;// edge 5
			
		for(int i = 4+6*em; i < 4+6*em+2*fm; ++i)
			dind(2,++ind)=i;// face 0-1
		
		for(int i = 4+6*em+3*fm; i < tm; ++i)
			dind(2,++ind)=i;// face 3 - interior		

		
		// vertex 3
		ind = 0;
		dind(3,ind) = 3;
		for(int i = 4; i < 4+2*em; ++i)
			dind(3,++ind)=i;//edge 1-2
			
		for(int i = 4+5*em; i < 4+6*em; ++i)
			dind(3,++ind)=i;//edge 6
			
		for(int i = 4+6*em; i < 4+6*em+3*fm; ++i)
			dind(3,++ind)=i;//face 0-2
		
		for(int i = bm; i < tm; ++i)
			dind(3,++ind)=i;//interior

		for(int i = 0; i < 4; ++i) {
            for(int j = 0; j < ltm; ++j)
				for(int k = 0; k < vdof; ++k)
					mwk(j,k)=mm(lind(j),dind(i,k));		
            vwk = 0;
            vwk(i) = 1;
            for(int j = 0; j < vdof; ++j)
                mwk(i,j) = 0;
            mwk(i,0) = 1;
#ifdef F2CFortran
            GETRF(ltm,vdof,&mwk(0,0),MXTM,&ipiv(0),info);
            GETRS(trans,ltm,1,&mwk(0,0),MXTM,&ipiv(0),&vwk(0),MXTM,info);
#else
            const int one = 1, mx = MXTM;
            dgetrf_(&ltm,&vdof,&mwk(0,0),&mx,&ipiv(0),&info);
            dgetrs_(trans,&ltm,&one,&mwk(0,0),&mx,&ipiv(0),&vwk(0),&mx,&info);
#endif
				
            for(int k=0;k<im;++k)
                ifmb(i,k) = -vwk(k+3*em+3*fm+1);
//				cout << vwk << endl;
//				
//          for(int j = 0; j < ltm; ++j)
//				for(int k = 0; k < vdof; ++k)
//					mwk(j,k)=mm(lind(j),dind(i,k));
//			
//			for(int j = 0; j < vdof; ++j)
//				mwk(i,j) = 0;
//			mwk(i,0) = 1;	
//					
//			for(int k = 0; k < MXTM; ++k){
//				wk = 0;
//				for(int j = 0; j < MXTM; ++j){
//					wk+=mwk(k,j)*vwk(j);
//				}
//				cout << wk << endl;
//			}
//			cout << endl << endl;
        
		}
		
	
		for(int k=0;k<em;++k) {
            sfmv(1,k) = -vwk(k+1);
            sfmv(0,k) = -vwk(k+1+em);      
        }
		
		for(int k=0;k<fm;++k) {
            ffmv(0,k) = 0;
            ffmv(1,k) = -vwk(k+3*em+1+fm); 
			ffmv(2,k) = -vwk(k+3*em+1+2*fm);           
        }
	
	vdiag=mm(0,0);
	for(int k=0;k<em;++k) {
		vdiag -= sfmv(1,k)*mm(4+k+3*em,0);
        vdiag -= sfmv(1,k)*mm(4+k+4*em,0);
        vdiag -= sfmv(1,k)*mm(4+k+5*em,0);
	}
	
	for(int k=0;k<fm;++k) {
		vdiag -= ffmv(0,k)*mm(4+k+6*em+fm,0);
        vdiag -= ffmv(0,k)*mm(4+k+6*em+2*fm,0);
        vdiag -= ffmv(0,k)*mm(4+k+6*em+3*fm,0);
	}
	
	for(int k=0;k<im;++k) {
		vdiag -= ifmb(0,k)*mm(4+k+6*em+4*fm,0);
	}
	vdiag=1/vdiag;
	

				   
   /*******************************************************/      
   /*  EQUATIONS TO FIND SIDE VALUES TO em-2 ACCURACY */
   /*******************************************************/ 
	for(int i = 0; i < 4; ++i){
		for(int j = 0; j < tm; ++j){
			vfms(i,j) = mm(i,j);
		}
	}
			
			
	ind = 4;
	for(int i = 4; i < 6*em+4; i=i+em){
		for(int j=0;j < tm; ++j){
			vfms(ind,j) = mm(i,j);
		}
		++ind;
	}
	
	for(int i = 4+6*em; i < 4+6*em+4*fm; i=i+fm){
		for(int j=0;j < tm; ++j){
			vfms(ind,j) = mm(i,j);
		}
		++ind;
	}
	
	// non lumped pure diagonal terms
	for(int i = 4; i < 4+em; ++i)
		ediag(i-4)=mm(i,i);
	
	
//	if(lumped){
//		//lumped high order terms
//		for(int j = 4; j < 4+em;++j){
//			ediag(j-4) = 0.0;
//			for(int i = 0; i < tm; ++i){
//				ediag(j-4)+=fabs(mm(i,j));
//			}
//		}
//	}

	for(int j = 0; j < em; ++j)
		ediag(j)=1/ediag(j);
	
	if(p == 3){	
		edgecone(0)(0)(0)=2;edgecone(0)(1)(0)=4;edgecone(0)(2)(0)=6;
		edgecone(1)(0)(0)=3;edgecone(1)(1)(0)=4;edgecone(1)(2)(0)=5;
		edgecone(2)(0)(0)=2;edgecone(2)(1)(0)=4;edgecone(2)(2)(0)=6;
		edgecone(3)(0)(0)=1;edgecone(3)(1)(0)=5;edgecone(3)(2)(0)=6;
		edgecone(4)(0)(0)=2;edgecone(4)(1)(0)=4;edgecone(4)(2)(0)=6;
		edgecone(5)(0)(0)=3;edgecone(5)(1)(0)=4;edgecone(5)(2)(0)=5;
		
		edgecone(0)(0)(1)=3;edgecone(0)(1)(1)=4;edgecone(0)(2)(1)=5;
		edgecone(1)(0)(1)=1;edgecone(1)(1)(1)=5;edgecone(1)(2)(1)=6;
		edgecone(2)(0)(1)=1;edgecone(2)(1)(1)=5;edgecone(2)(2)(1)=6;
		edgecone(3)(0)(1)=1;edgecone(3)(1)(1)=2;edgecone(3)(2)(1)=3;
		edgecone(4)(0)(1)=1;edgecone(4)(1)(1)=2;edgecone(4)(2)(1)=3;
		edgecone(5)(0)(1)=1;edgecone(5)(1)(1)=2;edgecone(5)(2)(1)=3;
		
//		ind = 0;
//		for(int m = 0; m < 6; ++m)
//			for(int j=4+em*m; j < 2+em*(m+1); ++j)
//				lind(ind++)=j;
				

//		for(int m = 0; m < 4; ++m){
//			ind2 = 0;
//			for(int i = 1; i <= em-1; ++i){
//				for(int j = 1; j <= em-i-1; ++j){
//					lind(ind++)=4+6*em+fm*m+ind2;
//					++ind2;
//				}
//				++ind2;
//			}
//		}
//		
//		ind2=0;
//		for(int i = 1; i <= em-1; ++i){
//			for(int j = 1; j <= em-i-1; ++j){
//				for(int k = 1; k <= em-i-j-1; ++k){
//					lind(ind++)=bm+ind2;
//					++ind2;
//				}
//				++ind2;
//			}
//			++ind2;
//		}
		//cout << lind(Range(0,9)) << endl;
			
													   
		// find all degrees of freedom for edges
		// edge 1
		ind = 0;
		for(int i = 5; i < 4+em; ++i)
			dind(0,ind++)=i;// edge 1
			
		for(int i = 4+6*em; i < 4+6*em+2*fm; ++i)
			dind(0,ind++)=i; // face 0-1
		
		for(int i = bm; i < tm; ++i)
			dind(0,ind++)=i; // interior
			
		// edge 2
		ind = 0;
		for(int i = 5+em; i < 4+2*em; ++i)
			dind(1,ind++)=i;// edge 2
			
		for(int i = 4+6*em; i < 4+6*em+fm; ++i)
			dind(1,ind++)=i; // face 0
			
		for(int i = 4+6*em+2*fm; i < 4+6*em+3*fm; ++i)
			dind(1,ind++)=i; // face 2
		
		for(int i = bm; i < tm; ++i)
			dind(1,ind++)=i; // interior
			
		// edge 3
		ind = 0;
		for(int i = 5+2*em; i < 4+3*em; ++i)
			dind(2,ind++)=i;// edge 3
			
		for(int i = 4+6*em; i < 4+6*em+fm; ++i)
			dind(2,ind++)=i; // face 0
			
		for(int i = 4+6*em+3*fm; i < tm; ++i)
			dind(2,ind++)=i; // face 3,interior
			
		// edge 4
		ind = 0;
		for(int i = 5+3*em; i < 4+4*em; ++i)
			dind(3,ind++)=i;// edge 4
						
		for(int i = 4+6*em+2*fm; i < tm; ++i)
			dind(3,ind++)=i; // face 2,3, interior
			
		// edge 5
		ind = 0;
		for(int i = 5+4*em; i < 4+5*em; ++i)
			dind(4,ind++)=i;// edge 5
			
		for(int i = 4+6*em+fm; i < 4+6*em+2*fm; ++i)
			dind(4,ind++)=i; // face 1
			
		for(int i = 4+6*em+3*fm; i < tm; ++i)
			dind(4,ind++)=i; // face 3, interior
		
		// edge 6
		ind = 0;
		for(int i = 5+5*em; i < 4+6*em; ++i)
			dind(5,ind++)=i;// edge 6
			
		for(int i = 4+6*em+fm; i < 4+6*em+3*fm; ++i)
			dind(5,ind++)=i; // face 1,2
		
		for(int i = bm; i < tm; ++i)
			dind(5,ind++)=i; // interior
	
			lind = 0;
		//cout << lind << endl << dind(3,Range(0,9)) << endl;
		for(int side = 0; side < 2; ++side){		
			for(int i = 0; i < 6; ++i){
				ind = 0;
				for(int j = 0; j < 3; ++j)
					lind(ind++)=4+(edgecone(i)(j)(side)-1)*em;
				
				for(int j = 0; j < 3; ++j)//p=3 only
					for(int k = 0; k < 3; ++k)
						mwk(j,k)=mm(lind(j),dind(i,k));		
				vwk=0;
					
				for(int j = 0; j < 3; ++j)
					vwk(j) = -mm(lind(j),4+i*em);
#ifdef F2CFortran
				GETRF(3,3,&mwk(0,0),MXTM,&ipiv(0),info);
#else
                const int three = 3, mx = MXTM;
                dgetrf_(&three,&three,&mwk(0,0),&mx,&ipiv(0),&info);
#endif

				if (info != 0) {
					printf("DGETRF FAILED - EDGE MODES info:%d EDGE:%d \n",info,i);
					exit(1);
				}
#ifdef F2CFortran
				GETRS(trans,3,1,&mwk(0,0),MXTM,&ipiv(0),&vwk(0),MXTM,info);
#else
                const int one = 1;
                dgetrs_(trans,&three,&one,&mwk(0,0),&mx,&ipiv(0),&vwk(0),&mx,&info);
#endif

				if (info != 0) {
					printf("DGETRS FAILED - EDGE MODES info:%d EDGE:%d \n",info,i);
					exit(1);
				}
				for(int k=0;k<im;++k)
					ifmb(i+4,k) = -vwk(k+em+2*fm);
				if(i < 3){
					for(int k=0;k<fm;++k)
						ffms(i,k,side) = vwk(k+em-1);  
				}	
								
				//cout << vwk(Range(0,2)) << endl;			
			}
		}
			
		sfms(0) = -vwk(0);	   
		
		int omode = 8;
		mdiag=mm(4,4);
		odiag=mm(4,omode);

		for(int k=1;k<em;++k) {
			mdiag += sfms(k-1)*mm(4+k,4);
			odiag += sfms(k-1)*mm(4+k,omode);
		}
		
		int side = 0;
		
		for(int k=0;k<fm;++k) {
			mdiag += ffms(0,k,side)*mm(4+k+6*em,4);
			mdiag += ffms(0,k,side)*mm(4+k+6*em+fm,4);
			odiag += ffms(0,k,side)*mm(4+k+6*em,omode);
			odiag += ffms(0,k,side)*mm(4+k+6*em+fm,omode);
		}
		
//		cout << ffms(Range(0,2),0,0) << ffms(Range(0,2),0,1) << endl;
//		cout << sfms(0)  << endl;

		/* check all modes to see if odiag and mdiag are correct*/
		//for(int i = 0; i < 6; ++i){
//			for(int j = 0; j < 6; ++j){
//				mwk(i,j)=mm(4+i*em,4+j*em);
//			}
//		}
//		
//		for(int i = 0; i < 6; ++i){
//			for(int j = 0; j < 6; ++j){
//				mwk(i,j)+=sfms(0)*mm(5+i*em,4+j*em);
//			}
//		}
//		
//		side = 0;
//		for(int j = 0; j < 6; ++j){
//			mwk(0,j)+=ffms(0,0,side)*mm(4+6*em,4+j*em);
//			mwk(0,j)+=ffms(0,0,side)*mm(4+6*em+fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(1,j)+=ffms(1,0,side)*mm(4+6*em,4+j*em);
//			mwk(1,j)+=ffms(0,0,side)*mm(4+6*em+2*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(2,j)+=ffms(1,0,side)*mm(4+6*em,4+j*em);
//			mwk(2,j)+=ffms(0,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(3,j)+=ffms(1,0,side)*mm(4+6*em+2*fm,4+j*em);
//			mwk(3,j)+=ffms(1,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(4,j)+=ffms(1,0,side)*mm(4+6*em+fm,4+j*em);
//			mwk(4,j)+=ffms(1,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(5,j)+=ffms(1,0,side)*mm(4+6*em+fm,4+j*em);
//			mwk(5,j)+=ffms(1,0,side)*mm(4+6*em+2*fm,4+j*em);
//		}
//
//		
//		cout << mwk(Range(0,5),Range(0,5))<< endl;
//		
//		for(int i = 0; i < 6; ++i){
//			for(int j = 0; j < 6; ++j){
//				mwk(i,j)=mm(4+i*em,4+j*em);
//			}
//		}
//		
//		for(int i = 0; i < 6; ++i){
//			for(int j = 0; j < 6; ++j){
//				mwk(i,j)-=sfms(0)*mm(5+i*em,4+j*em);
//			}
//		}
//		
//		side = 1;
//		for(int j = 0; j < 6; ++j){
//			mwk(0,j)+=ffms(0,0,side)*mm(4+6*em,4+j*em);
//			mwk(0,j)+=ffms(0,0,side)*mm(4+6*em+fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(1,j)+=ffms(1,0,side)*mm(4+6*em,4+j*em);
//			mwk(1,j)+=ffms(0,0,side)*mm(4+6*em+2*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(2,j)+=ffms(1,0,side)*mm(4+6*em,4+j*em);
//			mwk(2,j)+=ffms(0,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(3,j)+=ffms(1,0,side)*mm(4+6*em+2*fm,4+j*em);
//			mwk(3,j)+=ffms(1,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(4,j)+=ffms(1,0,side)*mm(4+6*em+fm,4+j*em);
//			mwk(4,j)+=ffms(1,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(5,j)+=ffms(1,0,side)*mm(4+6*em+fm,4+j*em);
//			mwk(5,j)+=ffms(1,0,side)*mm(4+6*em+2*fm,4+j*em);
//		}
//
//		
//		cout << mwk(Range(0,5),Range(0,5))<< endl;




//
//		ediag(0)=mm(4,4);
//
//		for(int k=0;k<em;++k) {
//			ediag(0) -= sfms(k)*mm(4+k,4);
//		}
//	
//		for(int k=0;k<fm;++k) {
//			ediag(0) -= ffms(0,k)*mm(4+k+6*em,4);
//			ediag(0) -= ffms(0,k)*mm(4+k+6*em+fm,4);
//		}
//	
//		for(int k=0;k<im;++k) {
//			ediag(0) -= ifmb(4,k)*mm(4+k+6*em+4*fm,4);
//		}
//		ediag(0)=1/ediag(0);
		
//		for(int i = 0; i < 3; ++i)
//			cout << 1/sdiag(i) << endl;

//		cout << sdiag(0) << endl;

		// test edge 2 for sdiag
//		sdiag(0) = mm(7,7);
//		for(int k=0;k<em;++k) {
//			sdiag(0) -= sfms(k)*mm(4+k+em,7);
//		}
//	
//		for(int k=0;k<fm;++k) {
//			sdiag(0) -= ffms(1,k)*mm(4+k+6*em,7);
//			sdiag(0) -= ffms(0,k)*mm(4+k+6*em+2*fm,7);
//		}
//	
//		for(int k=0;k<im;++k) {
//			sdiag(0) -= ifmb(5,k)*mm(4+k+6*em+4*fm,7);
//		}
//		
//		cout << sdiag(0) << endl;
   fdiag(0) = 1/mm(16,16); //for p = 3


}
	/*******************************************************/      
   /*  EQUATIONS TO FIND FACE VALUES TO em-2 ACCURACY */
   /*******************************************************/
   if(p == 4){
   
   		edgecone(0)(0)(0)=2;edgecone(0)(1)(0)=4;edgecone(0)(2)(0)=6;
		edgecone(1)(0)(0)=3;edgecone(1)(1)(0)=4;edgecone(1)(2)(0)=5;
		edgecone(2)(0)(0)=2;edgecone(2)(1)(0)=4;edgecone(2)(2)(0)=6;
		edgecone(3)(0)(0)=1;edgecone(3)(1)(0)=5;edgecone(3)(2)(0)=6;
		edgecone(4)(0)(0)=2;edgecone(4)(1)(0)=4;edgecone(4)(2)(0)=6;
		edgecone(5)(0)(0)=3;edgecone(5)(1)(0)=4;edgecone(5)(2)(0)=5;
		
		edgecone(0)(0)(1)=3;edgecone(0)(1)(1)=4;edgecone(0)(2)(1)=5;
		edgecone(1)(0)(1)=1;edgecone(1)(1)(1)=5;edgecone(1)(2)(1)=6;
		edgecone(2)(0)(1)=1;edgecone(2)(1)(1)=5;edgecone(2)(2)(1)=6;
		edgecone(3)(0)(1)=1;edgecone(3)(1)(1)=2;edgecone(3)(2)(1)=3;
		edgecone(4)(0)(1)=1;edgecone(4)(1)(1)=2;edgecone(4)(2)(1)=3;
       edgecone(5)(0)(1)=1;edgecone(5)(1)(1)=2;edgecone(5)(2)(1)=3;
		
													   
		// find all degrees of freedom for edges
		// edge 1
		ind = 0;
		for(int i = 5; i < 4+em; ++i)
			dind(0,ind++)=i;// edge 1
			
		for(int i = 4+6*em; i < 4+6*em+2*fm; ++i)
			dind(0,ind++)=i; // face 0-1
		
		for(int i = bm; i < tm; ++i)
			dind(0,ind++)=i; // interior
			
		// edge 2
		ind = 0;
		for(int i = 5+em; i < 4+2*em; ++i)
			dind(1,ind++)=i;// edge 2
			
		for(int i = 4+6*em; i < 4+6*em+fm; ++i)
			dind(1,ind++)=i; // face 0
			
		for(int i = 4+6*em+2*fm; i < 4+6*em+3*fm; ++i)
			dind(1,ind++)=i; // face 2
		
		for(int i = bm; i < tm; ++i)
			dind(1,ind++)=i; // interior
			
		// edge 3
		ind = 0;
		for(int i = 5+2*em; i < 4+3*em; ++i)
			dind(2,ind++)=i;// edge 3
			
		for(int i = 4+6*em; i < 4+6*em+fm; ++i)
			dind(2,ind++)=i; // face 0
			
		for(int i = 4+6*em+3*fm; i < tm; ++i)
			dind(2,ind++)=i; // face 3,interior
			
		// edge 4
		ind = 0;
		for(int i = 5+3*em; i < 4+4*em; ++i)
			dind(3,ind++)=i;// edge 4
						
		for(int i = 4+6*em+2*fm; i < tm; ++i)
			dind(3,ind++)=i; // face 2,3, interior
			
		// edge 5
		ind = 0;
		for(int i = 5+4*em; i < 4+5*em; ++i)
			dind(4,ind++)=i;// edge 5
			
		for(int i = 4+6*em+fm; i < 4+6*em+2*fm; ++i)
			dind(4,ind++)=i; // face 1
			
		for(int i = 4+6*em+3*fm; i < tm; ++i)
			dind(4,ind++)=i; // face 3, interior
		
		// edge 6
		ind = 0;
		for(int i = 5+5*em; i < 4+6*em; ++i)
			dind(5,ind++)=i;// edge 6
			
		for(int i = 4+6*em+fm; i < 4+6*em+3*fm; ++i)
			dind(5,ind++)=i; // face 1,2
		
		for(int i = bm; i < tm; ++i)
			dind(5,ind++)=i; // interior
	
		lind = 0;
		//cout << lind << endl << dind(3,Range(0,9)) << endl;
		for(int side = 0; side < 2; ++side){		
			for(int i = 0; i < 6; ++i){
				ind = 0;
				for(int j = 0; j < 3; ++j)
					for(int k = 0; k < em-1; ++k)
						lind(ind++)=4+(edgecone(i)(j)(side)-1)*em+k;
				for(int k = 0; k < 4; ++k)
					lind(ind++)=4+6*em+k*fm;
				
				for(int j = 0; j < 10; ++j)//p=4 only
					for(int k = 0; k < 9; ++k)
						mwk(j,k)=mm(lind(j),dind(i,k));		
				vwk=0;
					
				for(int j = 0; j < 10; ++j)
					vwk(j) = -mm(lind(j),4+i*em);
					
				for(int j = 0; j < 10; ++j)
					vwk(j+12) = -mm(lind(j),4+i*em);
#ifdef F2CFortran
				GETRF(9,10,&mwk(0,0),MXTM,&ipiv(0),info);
				if (info != 0) {
					printf("DGETRF FAILED - EDGE MODES info:%d EDGE:%d \n",info,i);
					exit(1);
				}
				GETRS(trans,10,1,&mwk(0,0),MXTM,&ipiv(0),&vwk(0),MXTM,info);
				if (info != 0) {
					printf("DGETRS FAILED - EDGE MODES info:%d EDGE:%d \n",info,i);
					exit(1);
				}
#else
                const int one = 1, nine = 9, ten = 10, mx = MXTM;
                dgetrf_(&nine,&ten,&mwk(0,0),&mx,&ipiv(0),&info);
                if (info != 0) {
                    printf("DGETRF FAILED - EDGE MODES info:%d EDGE:%d \n",info,i);
                    exit(1);
                }
                dgetrs_(trans,&ten,&one,&mwk(0,0),&mx,&ipiv(0),&vwk(0),&mx,&info);
                if (info != 0) {
                    printf("DGETRS FAILED - EDGE MODES info:%d EDGE:%d \n",info,i);
                    exit(1);
                }
#endif
				
//				cout << vwk << endl;
//				exit(0);
				
				for(int k=0;k<im;++k)
					ifmb(i+4,k) = -vwk(k+em+2*fm);
				if(i < 3){
					for(int k=0;k<fm;++k)
						ffms(i,k,side) = vwk(k+em-1);  
				}	
								
				//cout << vwk(Range(0,2)) << endl;			
			}
		}
			
		sfms(0) = -vwk(0);	   
		
		int omode = 8;
		mdiag=mm(4,4);
		odiag=mm(4,omode);

		for(int k=1;k<em;++k) {
			mdiag += sfms(k-1)*mm(4+k,4);
			odiag += sfms(k-1)*mm(4+k,omode);
		}
		
		int side = 0;
		
		for(int k=0;k<fm;++k) {
			mdiag += ffms(0,k,side)*mm(4+k+6*em,4);
			mdiag += ffms(0,k,side)*mm(4+k+6*em+fm,4);
			odiag += ffms(0,k,side)*mm(4+k+6*em,omode);
			odiag += ffms(0,k,side)*mm(4+k+6*em+fm,omode);
		}
		
//		cout << ffms(Range(0,2),0,0) << ffms(Range(0,2),0,1) << endl;
//		cout << sfms(0)  << endl;

		/* check all modes to see if odiag and mdiag are correct*/
		//for(int i = 0; i < 6; ++i){
//			for(int j = 0; j < 6; ++j){
//				mwk(i,j)=mm(4+i*em,4+j*em);
//			}
//		}
//		
//		for(int i = 0; i < 6; ++i){
//			for(int j = 0; j < 6; ++j){
//				mwk(i,j)+=sfms(0)*mm(5+i*em,4+j*em);
//			}
//		}
//		
//		side = 0;
//		for(int j = 0; j < 6; ++j){
//			mwk(0,j)+=ffms(0,0,side)*mm(4+6*em,4+j*em);
//			mwk(0,j)+=ffms(0,0,side)*mm(4+6*em+fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(1,j)+=ffms(1,0,side)*mm(4+6*em,4+j*em);
//			mwk(1,j)+=ffms(0,0,side)*mm(4+6*em+2*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(2,j)+=ffms(1,0,side)*mm(4+6*em,4+j*em);
//			mwk(2,j)+=ffms(0,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(3,j)+=ffms(1,0,side)*mm(4+6*em+2*fm,4+j*em);
//			mwk(3,j)+=ffms(1,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(4,j)+=ffms(1,0,side)*mm(4+6*em+fm,4+j*em);
//			mwk(4,j)+=ffms(1,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(5,j)+=ffms(1,0,side)*mm(4+6*em+fm,4+j*em);
//			mwk(5,j)+=ffms(1,0,side)*mm(4+6*em+2*fm,4+j*em);
//		}
//
//		
//		cout << mwk(Range(0,5),Range(0,5))<< endl;
//		
//		for(int i = 0; i < 6; ++i){
//			for(int j = 0; j < 6; ++j){
//				mwk(i,j)=mm(4+i*em,4+j*em);
//			}
//		}
//		
//		for(int i = 0; i < 6; ++i){
//			for(int j = 0; j < 6; ++j){
//				mwk(i,j)-=sfms(0)*mm(5+i*em,4+j*em);
//			}
//		}
//		
//		side = 1;
//		for(int j = 0; j < 6; ++j){
//			mwk(0,j)+=ffms(0,0,side)*mm(4+6*em,4+j*em);
//			mwk(0,j)+=ffms(0,0,side)*mm(4+6*em+fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(1,j)+=ffms(1,0,side)*mm(4+6*em,4+j*em);
//			mwk(1,j)+=ffms(0,0,side)*mm(4+6*em+2*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(2,j)+=ffms(1,0,side)*mm(4+6*em,4+j*em);
//			mwk(2,j)+=ffms(0,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(3,j)+=ffms(1,0,side)*mm(4+6*em+2*fm,4+j*em);
//			mwk(3,j)+=ffms(1,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(4,j)+=ffms(1,0,side)*mm(4+6*em+fm,4+j*em);
//			mwk(4,j)+=ffms(1,0,side)*mm(4+6*em+3*fm,4+j*em);
//		}
//		for(int j = 0; j < 6; ++j){
//			mwk(5,j)+=ffms(1,0,side)*mm(4+6*em+fm,4+j*em);
//			mwk(5,j)+=ffms(1,0,side)*mm(4+6*em+2*fm,4+j*em);
//		}
//
//		
//		cout << mwk(Range(0,5),Range(0,5))<< endl;




//
//		ediag(0)=mm(4,4);
//
//		for(int k=0;k<em;++k) {
//			ediag(0) -= sfms(k)*mm(4+k,4);
//		}
//	
//		for(int k=0;k<fm;++k) {
//			ediag(0) -= ffms(0,k)*mm(4+k+6*em,4);
//			ediag(0) -= ffms(0,k)*mm(4+k+6*em+fm,4);
//		}
//	
//		for(int k=0;k<im;++k) {
//			ediag(0) -= ifmb(4,k)*mm(4+k+6*em+4*fm,4);
//		}
//		ediag(0)=1/ediag(0);
		
//		for(int i = 0; i < 3; ++i)
//			cout << 1/sdiag(i) << endl;

//		cout << sdiag(0) << endl;

		// test edge 2 for sdiag
//		sdiag(0) = mm(7,7);
//		for(int k=0;k<em;++k) {
//			sdiag(0) -= sfms(k)*mm(4+k+em,7);
//		}
//	
//		for(int k=0;k<fm;++k) {
//			sdiag(0) -= ffms(1,k)*mm(4+k+6*em,7);
//			sdiag(0) -= ffms(0,k)*mm(4+k+6*em+2*fm,7);
//		}
//	
//		for(int k=0;k<im;++k) {
//			sdiag(0) -= ifmb(5,k)*mm(4+k+6*em+4*fm,7);
//		}
//		
//		cout << sdiag(0) << endl;
   fdiag(0) = 1/mm(16,16); //for p = 3
   
   
		lind=-1;
		ind = 0;
			
		for(int m = 0; m < 4; ++m){
			ind2 = 0;
			for(int i = 1; i <= em-1; ++i){
				for(int j = 1; j <= em-i-1; ++j){
					lind(ind++)=4+6*em+fm*m+ind2;
					++ind2;
				}
				++ind2;
			}
		}
		
		ind2=0;
		for(int i = 1; i <= em-1; ++i){
			for(int j = 1; j <= em-i-1; ++j){
				for(int k = 1; k <= em-i-j-1; ++k){
					lind(ind++)=bm+ind2;
					++ind2;
				}
				++ind2;
			}
			++ind2;
		}
		//cout << lind(Range(0,3)) << endl;
			
													   
		// find all degrees of freedom for faces
		ind = 0;			
		for(int i = 0; i < 4; ++i){
			ind = 0;
			for(int j = 4+6*em+i*fm; j < 4+6*em+(i+1)*fm; ++j){
				dind(i,ind++)=j; 
			}
		}
		
		for(int i = bm; i < tm; ++i)
			dind(Range::all(),ind++)=i; // interior
			
			//cout << dind(Range(0,3),Range(0,3)) << endl;
			
		for(int i = 0; i < 4; ++i){
			for(int j = 0; j < 4; ++j)
				for(int k = 0; k < 4; ++k)
					mwk(j,k)=mm(lind(j),dind(i,k));		
            vwk = 0;
            vwk(i) = 1;
            for(int j = 0; j < 4; ++j)
                mwk(i,j) = 0;
            mwk(i,0) = 1;
            
#ifdef F2CFortran
            GETRF(4,4,&mwk(0,0),MXTM,&ipiv(0),info);
            GETRS(trans,4,1,&mwk(0,0),MXTM,&ipiv(0),&vwk(0),MXTM,info);
#else
            const int one = 1, four = 4, mx = MXTM;
            dgetrf_(&four,&four,&mwk(0,0),&mx,&ipiv(0),&info);
            dgetrs_(trans,&four,&one,&mwk(0,0),&mx,&ipiv(0),&vwk(0),&mx,&info);
#endif
				
            for(int k=0;k<im;++k)
                ifmb(i+4*im+6*im,k) = -vwk(k+fm);
            
            //cout << vwk << endl;
				
//          for(int j = 0; j < ltm; ++j)
//				for(int k = 0; k < vdof; ++k)
//					mwk(j,k)=mm(lind(j),dind(i,k));
//			
//			for(int j = 0; j < vdof; ++j)
//				mwk(i,j) = 0;
//			mwk(i,0) = 1;	
//					
//			for(int k = 0; k < MXTM; ++k){
//				wk = 0;
//				for(int j = 0; j < MXTM; ++j){
//					wk+=mwk(k,j)*vwk(j);
//				}
//				cout << wk << endl;
//			}
//			cout << endl << endl;
				
		}
		//cout << ifmb << endl;
		ind = 0;
		for(int i = 4+6*em; i < 4+6*em+fm; ++i)
			fdiag(ind++)=mm(i,i);			
		idiag(0) = mm(tm-1,tm-1);//interior mode diagonal p4 only
		
		//lumped high order terms
		if(lumped){
			ind = 0;
			for(int i = 4+6*em; i < 4+6*em+fm; ++i){
				fdiag(ind) = 0.0;
				for(int j = 0; j < tm; ++j){
					fdiag(ind)+=fabs(mm(i,j));	
				}
				++ind;
			}
			idiag(ind) = 0.0;
			ind = 0;
			for(int j = 0; j < tm; ++j){
				idiag(ind)+=fabs(mm(tm-1,j));	
			}
		}
		for(int j = 0; j < fm; ++j)
			fdiag(j)=1/fdiag(j);
		for(int j = 0; j < im; ++j)
			idiag(j)=1/idiag(j);
		
//  //    uncomment if including orthogonal face modes
//		tdiag(0) = mm(4+6*em,4+6*em);
//		for(int k=0;k<im;++k) {
//			tdiag(0) -= ifmb(4*im+6*im,k)*mm(4+k+6*em+4*fm,4+6*em);
//		}
//		tdiag(0)=1/tdiag(0);
		//cout << tdiag << endl;
	}

   return;
}         


void tet_basis::legpt(){
   FLT x,y,z,r,s,t;
   int ind,ind2,ind3,sign;
   
   lgrnge1d.resize(em+2,em+1);
   lgrnge2d.resize(3+3*em+fm,em+1,em+1);
   lgrnge3d.resize(tm,em+1,em+1,em+1);

   /* CALCULATE PROJECTION POINTS IN INTERIOR */
	for(int i = 1; i < em; ++i){
		for(int j = 1; j < em-i+1; ++j){
			for(int k = 1; k < em-i-j+2; ++k){	  
				r = -1.0 +2.0*((FLT) i)/(FLT)(em+1);
				s = -1.0 +2.0*((FLT) j)/(FLT)(em+1);
				t = -1.0 +2.0*((FLT) k)/(FLT)(em+1);

				x = 2.0*(1+r)/(-s-t)-1.0;
				y = 2.0*(1+s)/(1-t)-1.0;
				z = t;

				/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
				ptvalues(x,y,z);				
				
				/* VERTEX 0 */
				lgrnge3d(0,i,j,k) =  pgz(0); //note: pgx,pgy=1
			  
				/* VERTEX 1 */
				lgrnge3d(1,i,j,k) =  pgz(1)*pgy(1); //pgx=1
				
				/* VERTEX 2 */
				lgrnge3d(2,i,j,k) =  pgz(1)*pgy(2)*pgx(1);
				
				/* VERTEX 3 */
				lgrnge3d(3,i,j,k) =  pgz(1)*pgy(2)*pgx(2);
				
				ind = 3;
				/* EDGE 1 */
				for(int m = 2; m <= em+1; ++m) {
					lgrnge3d(++ind,i,j,k) =  pgz(m)*pgy(m+1)*pgx(m+1);
				}
				
				/* EDGE 2 */
				for(int m = 2; m <= em+1; ++m) {
					lgrnge3d(++ind,i,j,k) =  pgz(m)*pgy(em+m+1)*pgx(2);
				}
				
				/* EDGE 3 */
				for(int m = 2; m <= em+1; ++m) {
					lgrnge3d(++ind,i,j,k) =  pgz(m)*pgy(em+m+1)*pgx(1);
				}
				
				/* EDGE 4 */
				for(int m = 2; m <= em+1; ++m) {
					lgrnge3d(++ind,i,j,k) =  pgz(em+m)*pgy(1); //pgx=1
				}
				
				/* EDGE 5 */
				for(int m = 2; m <= em+1; ++m) {
					lgrnge3d(++ind,i,j,k) =  pgz(em+m)*pgy(2)*pgx(1);
				}
				
				/* EDGE 6 */
				for(int m = 2; m <= em+1; ++m) {
					lgrnge3d(++ind,i,j,k) =  pgz(em+m)*pgy(2)*pgx(2);
				}
				
				/* FACE 0 */
				ind2 = 0;
				for(int p = 1; p <= em-1; ++p) {      
					for(int q = 1; q <= em-p; ++q) {
						++ind2;
						lgrnge3d(++ind,i,j,k) =  pgz(p+q+1)*pgy(2+2*em+ind2)*pgx(p+2);
					}
				}
				
				/* FACE 1 */
				ind2 = 0;
				sign = 1;
				for(int p = 1; p <= em-1; ++p) {      
					for(int q = 1; q <= em-p; ++q) {
						++ind2;
						lgrnge3d(++ind,i,j,k) =  sign*pgz(1+2*em+ind2)*pgy(p+2)*pgx(p+2);
					}
					sign*=-1;
				}
				
				/* FACE 2 */
				ind2 = 0;
				sign = 1;
				for(int p = 1; p <= em-1; ++p) {      
					for(int q = 1; q <= em-p; ++q) {
						++ind2;
						lgrnge3d(++ind,i,j,k) =  sign*pgz(1+2*em+ind2)*pgy(em+2+p)*pgx(2);
					}
					sign*=-1;
				}
				
				/* FACE 3 */
				ind2 = 0;
				for(int p = 1; p <= em-1; ++p) {      
					for(int q = 1; q <= em-p; ++q) {
						++ind2;
						lgrnge3d(++ind,i,j,k) =  pgz(1+2*em+ind2)*pgy(em+2+p)*pgx(1);
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
							lgrnge3d(++ind,i,j,k) =  pgz(1+2*em+fm+ind2)*pgy(2+2*em+ind3)*pgx(p+2);
						}
					}
				}

				
			}
		}
	}
	
	
	/* CALCULATE PROJECTION POINTS IN FACE */
	
	for(int i = 1; i < em; ++i){
		for(int j = 1; j < em-i+1; ++j){
				  
			r = -1.0 +2.0*((FLT) i)/(FLT)(em+1);
			s = -1.0 +2.0*((FLT) j)/(FLT)(em+1);
	
			x = 2.0*(1+r)/(1.0-s)-1.0;
			y = s;

			/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
			ptvalues2d(x,y);	
					  
			/* VERTEX 1 */
			lgrnge2d(0,i,j) =  pgy(0); //pgx=1
			
			/* VERTEX 2 */
			lgrnge2d(1,i,j) =  pgy(1)*pgx(1);
			
			/* VERTEX 3 */
			lgrnge2d(2,i,j) =  pgy(1)*pgx(2);
			
			ind = 2;
			/* EDGE 1 */
			for(int m = 2; m <= em+1; ++m){
				lgrnge2d(++ind,i,j) =  pgy(m)*pgx(m+1);
			}
			
			/* EDGE 2 */
			for(int m = 2; m <= em+1; ++m){
				lgrnge2d(++ind,i,j) =  pgy(em+m)*pgx(2);
			}
			
			/* EDGE 3 */
			sign = 1;
			for(int m = 2; m <= em+1; ++m){
				lgrnge2d(++ind,i,j) =  sign*pgy(em+m)*pgx(1);
				sign*=-1;
			}
			
			
			/* FACE 0 */
			ind2 = 0;
			for(int p = 1; p <= em-1; ++p){      
				for(int q = 1; q <= em-p; ++q){
					++ind2;
					lgrnge2d(++ind,i,j) =  pgy(1+2*em+ind2)*pgx(p+2);
				}
			}				

		}
	}
	
   
   /* NOW CALCULATE VALUES OF G FOR EDGE PROJECTION */
   for (int i = 1; i < em+1; ++i) {
      x = -1.0+2.0*((FLT) i)/(FLT)(em+1);
      ptvalues1d(x);
      for(int m = 0; m < em+2; ++m)
         lgrnge1d(m,i) = pgx(m);
   }
      
//   cout << lgrnge3d << lgrnge2d << lgrnge1d << endl;
   return;
}

