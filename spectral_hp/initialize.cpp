/*
 *  initialize.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include<math.h>
#include<utilities.h>
#include"hpbasis.h"
#include<myblas.h>

int hpbasis::wkpmax = 0;
const int hpbasis::sbwth;
FLT **hpbasis::wk0, **hpbasis::wk1, **hpbasis::wk2, **hpbasis::wk3;
FLT *hpbasis::pgx, *hpbasis::dpgx, *hpbasis::pgn, *hpbasis::dpgn;


void hpbasis::initialize(int pdegree) {

   if (pdegree < 1) {
      printf("error can't use 0th order basis\n");
      exit(1);
   }
   
   p = pdegree;
   sm = p -1;
   im = (p-2)*(p-1)/2;
   bm = 3*p;
   tm = bm + im;
   ibwth = (2*(sm-1) < im -1? 2*(sm-1) : im-1);
   ibwth = (ibwth > 0 ? ibwth : 0);
   
   nmodx = p+2;
   nmodn = tm;
   gpx = p +1;
   gpn = p +1;
   
/*	WORK VARIABLES */   
   if (p > wkpmax) {
      if (wkpmax != 0) {
         free(wk0);
         free(wk1);
         free(wk2);
         free(wk3);
         free(pgx);
         free(dpgx);
         free(pgn);
         free(dpgn);  
         printf("warning: better to allocate hpbasis from largest to smallest\n");
      }
      mat_alloc(wk0,nmodx,gpn,FLT);
      mat_alloc(wk1,nmodx,gpn,FLT);
      mat_alloc(wk2,nmodx,gpn,FLT);
      mat_alloc(wk3,nmodx,gpn,FLT);
      vect_alloc(pgx,nmodx,FLT);
      vect_alloc(dpgx,nmodx,FLT);
      vect_alloc(pgn,tm,FLT);      
      vect_alloc(dpgn,tm,FLT);
      
      wkpmax = p;
   }
   
/*****************************/
/*	SETUP VALUES OF FuploCTIONS */
/*****************************/
   initialize_values(); // SET UP THINGS FOR PROJECT/INTEGRATE/DERIV
   lumpinv(); // SET UP THINGS FOR INVERSE OF LUMPED MASS MATRIX
   legpt(); // SET UP PROJECTION TO LEGENDRE POINTS (FOR OUTPUTING)
   
   return;
}

void hpbasis::initialize_values(void)
{

	FLT al,be,x,eta;
	FLT e[MXTM],e1[MXTM],e2[MXTM];
   FLT pkp, pk, pkm, dpkp, dpk, dpkm;
	int i,j,k,m,n,ipoly,ierr,ind;
   
/*	ALLOCATE STORAGE FOR RECURSION RELATION COEFFICENTS*/
   ipoly = MAX(gpx+1,sm+1);
   ipoly = MAX(ipoly,gpn+1);
   mat_alloc(a0,sm+2,ipoly,FLT);
   mat_alloc(b0,sm+2,ipoly,FLT);
   
/*	ALLOCATE INTEGRATION, PROJECTION, & DERIVATIVE VARIABLES */
   mat_alloc(gx,nmodx,gpx,FLT);
   mat_alloc(dgx,nmodx,gpx,FLT);
   vect_alloc(wtx,gpx,FLT);
   vect_alloc(x0,gpx,FLT);
   mat_alloc(gxwtx,nmodx,gpx,FLT);
   mat_alloc(dgxwtx,nmodx,gpx,FLT);
   vect_alloc(dltx,gpx,FLT);
   mat_alloc(dltx1,gpx,gpx,FLT);
   
   mat_alloc(gn,nmodn,gpn,FLT);
   mat_alloc(dgn,nmodn,gpn,FLT);	
   vect_alloc(wtn,gpn,FLT);
   vect_alloc(n0,gpn,FLT);
   mat_alloc(gnwtnn0,nmodn,gpn,FLT);
   mat_alloc(dgnwtn,nmodn,gpn,FLT);
   vect_alloc(dltn,gpn,FLT);
   mat_alloc(dltn1,gpn,gpn,FLT);
   mat_alloc(dltn2,gpn,gpn,FLT);
   vect_alloc(norm,tm,FLT);
		
/*	GENERATE RECURSION RELATION FOR LEGENDRE
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
/*	CALCULATE GAUSS POINTS / COEFFICIENTS  */
/*	*************************************	*/
   ipoly = 1;
   al = 0.0;
   be = 0.0;
   ierr = recur(gpx+1,ipoly,al,be,a0[0],b0[0]);
   if (ierr != 0) {
      printf("recur #1 error %d\n",ierr);
      return;
   }

   ierr = gauss(gpx,a0[0],b0[0],EPSILON,x0,wtx,e);
   if (ierr != 0) {
      printf("gauss #1 error %d\n",ierr);
      return;
   }

/***********************************************/
/*	CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
/***********************************************/
   for (i=0;i<gpx;++i) {
      dltx[i] = 1.0;
      for(k=0;k<i;++k)
         dltx[i] *= (x0[i]-x0[k]);
      for(k=i+1;k<gpx;++k)
         dltx[i] *= (x0[i]-x0[k]);
   }
   
   for (i=0;i<gpx;++i) {
      for(n=0;n<i;++n) 
         dltx1[i][n] = dltx[i]/(x0[i]-x0[n]);
      for(n=i+1;n<gpx;++n) 
         dltx1[i][n] = dltx[i]/(x0[i]-x0[n]);
   }

   for (i=0;i<gpx;++i) {
      for(n=0;n<i;++n) 
         dltn1[i][n] = dltx[i]/(x0[i]-x0[n])*.5*(1+x0[i]);
      for(n=i+1;n<gpx;++n) 
         dltn1[i][n] = dltx[i]/(x0[i]-x0[n])*.5*(1+x0[i]);
   }
		
   for (i=0;i<gpx;++i)
      dltx[i] = 1.0/dltx[i];	

/*************************************************/ 
/*	NOW CALCULATE VALUES OF G, G' AT GAUSS POINTS */
/*************************************************/
   for(i = 0;i < gpx; ++i) {
      x = x0[i];
      x0[i] = 0.5*(1+x);

      gx[0][i] = .5*(1-x);
      dgx[0][i] = -.5;
      gx[1][i] = .5*(1+x);
      dgx[1][i] = .5;
      gx[2][i] = 1;
      dgx[2][i] = 0;

      if (sm) {
/*			SIDE 1 IPOLY = 6 FOR JACOBI	*/
         ipoly = 6;
         al = 1.0;
         be = 1.0;
         ierr = recur(sm+1,ipoly,al,be,a0[0],b0[0]);
         if (ierr != 0) {
            printf("recur #3 error %d\n",ierr);
            return;
         }
   
   
/*			CALCULATE P, P' USING RECURSION RELATION */
         pk = 1.0;
         pkm = 0.0;
         dpk = 0.0;
         dpkm = 0.0;
   
         for (m = 1;m < sm+1;++m) {
            gx[m+2][i] = (1.+x)*(1.-x)*.25*pk;
            dgx[m+2][i] = -x*.5*pk +(1.+x)*(1.-x)*.25*dpk;
            pkp = (x-a0[0][m-1])*pk - b0[0][m-1]*pkm;
            dpkp = pk + (x-a0[0][m-1])*dpk - b0[0][m-1]*dpkm;
            dpkm = dpk;
            dpk = dpkp;
            pkm = pk;
            pk = pkp;
         }
      }
   }	


/*********************************************/
/*	CALCULATE GAUSS POINTS / COEFFICIENTS (ETA=N) */
/*	***************************************	*/
   ipoly = 6;
   al = 1.0;
   be = 0.0;
   ierr = recur(gpn,ipoly,al,be,a0[1],b0[1]);
   if (ierr != 0) {
      printf("recur #2 error %d\n",ierr);
      return;
   }	
   
   ierr = radau(gpn-1,a0[1],b0[1],-1.0,n0,wtn,e,e1,e2);		
   if (ierr != 0) {
      printf("gauss #3 error %d\n",ierr);
      return;
   }
   
   for (i=0;i<gpn;++i)
      wtn[i] *= 0.5;

/***********************************************/
/*	CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
/***********************************************/
   for (i=0;i<gpn;++i) {
      dltn[i] = 1.0;
      for(k=0;k<i;++k)
         dltn[i] *= (n0[i]-n0[k]);
      for(k=i+1;k<gpn;++k)
         dltn[i] *= (n0[i]-n0[k]);
   }
   
   for (j=0;j<gpn;++j) {
      for(n=0;n<j;++n) 
         dltn2[j][n] = dltn[j]/(n0[j]-n0[n]);
      for(n=j+1;n<gpn;++n) 
         dltn2[j][n] = dltn[j]/(n0[j]-n0[n]);
   }
   
   for (j=0;j<gpn;++j)
      dltn[j] = 1.0/dltn[j];
   

/******************************************/
/* GENERATE JACOBI POLY FOR S DIRECTION */
/****************************************/
   for(i=0;i<gpn;++i) {
      eta = n0[i];
      n0[i] = 2.0/(1-eta);

/*		VERTEX A	*/
      ind = 0;
      gn[ind][i] = (1-eta)*.5;
      dgn[ind][i] = -.5;

/*	 	VERTEX B  */
      ind = ind+1;
      gn[ind][i] = (1-eta)*.5;
      dgn[ind][i] = -.5;

/*		 VERTEX C	 */	
      ind = ind+1;
      gn[ind][i] = (1+eta)*.5;
      dgn[ind][i] = .5;

      if (sm) {
/*	 	 	SIDE 1 (s)		*/
         for(m = 2; m <= sm+1; ++m) {
            ++ind;
            gn[ind][i] = pow(.5*(1-eta),m);
            dgn[ind][i] = -.5*m*pow(.5*(1.-eta),m-1);
         }
   
   
/*			SIDE 2 IPOLY = 6 FOR JACOBI	*/
         ipoly = 6;
         al = 1.0;
         be = 1.0;
         ierr = recur(sm+1,ipoly,al,be,a0[1],b0[1]);
         if (ierr != 0) {
            printf("recur #3 error %d\n",ierr);
            return;
         }
   
/*			SIDE 2	*/
         pk = 1.0;
         pkm = 0.0;
         dpk = 0.0;
         dpkm = 0.0;
   
         for(n=1;n<=sm;++n) {
            ++ind;
            gn[ind][i] = (1.-eta)*(1.+eta)*.25*pk;		
            dgn[ind][i] = -.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk;
            pkp = (eta-a0[1][n-1])*pk - b0[1][n-1]*pkm;
            dpkp = pk + (eta-a0[1][n-1])*dpk - b0[1][n-1]*dpkm;
            dpkm = dpk;
            dpk = dpkp;
            pkm = pk;
            pk = pkp;
         }
   
/*   	 	SIDE 3	*/
         pk = 1.0;
         pkm = 0.0;
         dpk = 0.0;
         dpkm = 0.0;
   
         for(n=1;n<=sm;++n) {
            ++ind;
            gn[ind][i] = (n % 2 ? 1 : -1)*(1.-eta)*(1.+eta)*.25*pk;
            dgn[ind][i] = (n % 2 ? 1 : -1)*(-.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk);
            pkp = (eta-a0[1][n-1])*pk - b0[1][n-1]*pkm;
            dpkp = pk + (eta-a0[1][n-1])*dpk - b0[1][n-1]*dpkm;
            dpkm = dpk;
            dpk = dpkp;
            pkm = pk;
            pk = pkp;
         }
   
/*	 	 	INTERIOR MODES	*/
         if (im) {
            ind = bm;
            for(m = 2; m< sm+1;++m) {		
/*	 				CALCULATE RECURSION RELATION FOR P^(2m-1,1)(s)
               SIDE 1 IPOLY = 6 FOR JACOBI  	*/
               ipoly = 6;
               al = 2.*m-1;
               be = 1.0;
               ierr = recur(sm+2-m,ipoly,al,be,a0[m],b0[m]);
               if (ierr != 0) {
                  printf("recur #4 error %d\n",ierr);
                  return;
               }
      
               pk = 1.0;
               pkm = 0.0;
               dpk = 0.0;
               dpkm = 0.0;
      
               for(n = 1; n < sm+2-m;++n) {
                  gn[ind][i] = pow(.5*(1.-eta),m)*.5*(1.+eta)*pk;
                  dgn[ind][i] = -.25*m*pow(.5*(1.-eta),m-1)*(1.+eta)*pk 
                           +pow(.5*(1.-eta),m)*.5*(pk + (1.+eta)*dpk);
                  pkp = (eta-a0[m][n-1])*pk - b0[m][n-1]*pkm;
                  dpkp = pk + (eta-a0[m][n-1])*dpk - b0[m][n-1]*dpkm;
                  dpkm = dpk;
                  dpk = dpkp;
                  pkm = pk;
                  pk = pkp;
                  ++ind;
               }
            }
         }
      }
   }

/******************************/
/*	CALCULATE NORM          	*/
/* ************************** */
/*	SIDE & VERTEX MODES */
   for(n=0;n<3;++n)
      norm[n] = 1.0;
	for (n = 3; n < sm+3; ++n) {
		norm[n] = 0.0;
		for(i = 0; i < gpx; ++i)
			norm[n] += wtx[i]*gx[n][i]*gx[n][i];
		norm[n] = 1./sqrt(norm[n]);
      norm[n+sm] = norm[n];
      norm[n+2*sm] = norm[n];
	}
	
/*	INTERIOR MODES */
	for(m = bm; m < tm; ++m) {
		norm[m] = 0.0;
		for(j = 0; j < gpn; ++j)
	  		norm[m] += wtn[j]*gn[m][j]*gn[m][j];
	  	norm[m] = 1./sqrt(norm[m]);
	}

/***************/
/* RENORMALIZE */
/***************/
   for (n = 3; n < nmodx; ++n) {			
      for (i =0;i<gpx;++i) {
         gx[n][i] *= norm[n];
         dgx[n][i] *= norm[n];
      }
      for (j=0;j<gpn;++j) {
/*			SIDE 2 & 3 MUST BE RENORMALIZED BY SAME CONSTANT TO MATCH */
         gn[n +sm][j] *= norm[n];
         dgn[n+sm][j] *= norm[n];
         gn[n +2*sm][j] *= norm[n];
         dgn[n +2*sm][j] *= norm[n];
      }
   }

   for(m = bm; m < tm; ++m) {
      for(i = 0; i < gpn; ++i) {
         gn[m][i] = gn[m][i]*norm[m];
         dgn[m][i] = dgn[m][i]*norm[m];
      }
   }

/*****************************************/
/*	PRECALCULATE THINGS FOR FAST INTGRTRS */
/*****************************************/
   for(m=0;m<nmodx;++m) {
      for(i=0;i<gpx;++i) {
         gxwtx[m][i] = gx[m][i]*wtx[i];
         dgxwtx[m][i] = dgx[m][i]*wtx[i];
      }
   }
   
   for(m=0;m<tm;++m) {
      for(j=0;j<gpn;++j) {
         gnwtnn0[m][j] = gn[m][j]*wtn[j]*n0[j];
         dgnwtn[m][j] = dgn[m][j]*wtn[j];
      }
   }

	return;
}


/************************************************/
/** CALCULATE THINGS FOR LUMPED MASS INVERSION  */
/************************************************/
void hpbasis::lumpinv(void) {
	int i,i1,i2,j,k,m,info,ind,ind1,n,ipiv[2*MXTM];
	FLT mwk[MXTM][MXTM],vwk[MXTM], mm[MXTM][MXTM];
   FLT u[MXTM],l[MXTM];
   FLT rcond=1;
   char trans[] = "T", uplo[] = "U";
   
/*	ALLOCATE MASS MATRIX INVERSION VARIABLES */
   mat_alloc(vfms,3,sm,FLT);
   mat_alloc(sfmv,2,sm,FLT);				
   vect_alloc(sdiag,sm,FLT);
   tens_alloc(sfms,sm,sm-1,3,FLT);
   mat_alloc(ifmb,bm,im,FLT);
   mat_alloc(bfmi,bm,im,FLT);
   mat_alloc(idiag,im,ibwth+1,FLT);
   
/*	ALLOCATE 1D MASS MATRIX INVERSION VARIABLES */
   mat_alloc(vfms1d,2,sm,FLT);
   mat_alloc(sfmv1d,2,sm,FLT);				
   mat_alloc(sdiag1d,sm,sbwth+1,FLT);

/********************************************/	
/*		GENERATE MASS MATRIX 					  */
/********************************************/
   for(m=0;m<tm;++m) {
      for(i=0;i<tm;++i)
         u[i] = 0.0;
      u[m] = 1.0;

      proj(u,wk1);  // PROJECT USES WK0
      intgrt(wk1,l); // INTGRT USES WK0
      for(i=0;i<tm;++i) 
         mm[m][i] = l[i];		
   }

#ifdef SKIP
   printf("\n mm MATRIX SM = %d\n",(sm+2));
   for(i=0;i<tm;++i) {
      printf("%2d:",i);
      for(j=0;j<tm;++j)
         printf("%+.8lf ",mm[i][j]);
      printf("\n");
   }
#endif

/*******************************************************/		
/*	  	EQUATIONS TO FIND VERTEX VALUES TO SM-1 ACCURACY */
/*******************************************************/		
   if (p > 1) {		
      for(i=0;i<3;++i) {
         i1 = (i+1)%3;
         i2 = (i+2)%3;
         
/*			2 CROSS VERTEX CONSTRAINTS */
         for(k=0;k<2;++k) {
            ind = (i+1+k)%3;
            vwk[k] = mm[ind][i];
            for(j=0;j<sm;++j) {
               mwk[k][j] = mm[ind][3+j+i*sm];
               mwk[k][j+sm] = mm[ind][3+j+i2*sm];
            }
            for(j=0;j<im;++j)
               mwk[k][2*sm+j] = mm[ind][j+bm];
         }

/*			3x(SM-3) SIDE CONSTRAINTS */
         for(k=0;k<3;++k) {
            for(m=0;m<sm-1;++m) {	
               vwk[m+2+k*(sm-1)] = mm[3+m+k*sm][i];
               for(j=0;j<sm;++j) {
                  mwk[m+2+k*(sm-1)][j] = mm[3+m+k*sm][3+j+i*sm];
                  mwk[m+2+k*(sm-1)][j+sm] = mm[3+m+k*sm][3+j+i2*sm];
               }
               for(j=0;j<im;++j)
                  mwk[m+2+k*(sm-1)][2*sm+j] = mm[3+m+k*sm][j+bm];
            }
         }
         
         
/*			(SM-3)*(SM-4) INTERNAL MODE CONSTRAINTS  */				
         ind = 3*sm-1;	ind1 = 0;
         for(m=2;m<sm+1;++m) {		
            for(n = 1; n < sm+1-m;++n) {
               vwk[ind] = mm[ind1+bm][i];
               for(j=0;j<sm;++j) {
                  mwk[ind][j] = mm[ind1+bm][3+j+i*sm];
                  mwk[ind][j+sm] = mm[ind1+bm][3+j+i2*sm];
               }
               for(j=0;j<im;++j)
                  mwk[ind][2*sm+j] = mm[ind1+bm][j+bm];
               ++ind;
               ++ind1;
            }
            ++ind1;
         }

         GETRF((sm+1)*(sm+2)/2-1,(sm+1)*(sm+2)/2-1,mwk[0],MXTM,ipiv,info);
         if (info != 0) {
            printf("DGETRF FAILED - VRTX info:%d sm:%d i:%d\n",info,(sm+2),i);
            exit(1);
         }
         GETRS(trans,(sm+1)*(sm+2)/2-1,1,mwk[0],MXTM,ipiv,vwk,MXTM,info);
                                       
/*			STORE INTERIOR VALUES */
         for(k=0;k<im;++k)
            ifmb[i][k] = vwk[k+2*sm];			
      }
      
/*		STORE SIDE VALUES */
      for(k=0;k<sm;++k) {
         sfmv[0][k] = vwk[k];
         sfmv[1][k] = vwk[k+sm];
      }

/*		FIND VERTEX DIAGANOL ELEMENT */
      vdiag = mm[2][2];
      for(k=0;k<sm;++k) {
         vdiag -= sfmv[0][k]*mm[2][3+k+2*sm];
         vdiag -= sfmv[1][k]*mm[2][3+k+sm];
      }

      for(k=0;k<im;++k)
         vdiag -= ifmb[2][k]*mm[i][k+bm];
   }
   else {
      vdiag = mm[0][0] + mm[0][1] + mm[0][2];
   }


/*****************************************************************/
/*		EQUATIONS TO FIND STATIC INVERSION OF INTERIORS FROM SIDES */
/*****************************************************************/
/** WARNING THIS ONLY WORKS UP TO SM = 5!!!! ***/
   if (im > 0) {
      for(i=0;i<3;++i) {
         i1 = (i+1)%3;
         i2 = (i+2)%3;

         for(k=0;k<(sm+2)-3;++k) {
            ind = 0;
/*				CROSS SIDE MODES */
            for(j=k;j<sm-1;++j) {
               vwk[ind] =  mm[i*sm+3+k][i1*sm+3+j];
               for(m=0;m<im;++m)  
                  mwk[ind][m] = mm[bm+m][i1*sm+3+j];
               ++ind;
            }
            
            for(j=k;j<MIN(2*k,sm-1);++j) {
               vwk[ind] =  mm[i*sm+3+k][i2*sm+3+j];
               for(m=0;m<im;++m)  
                  mwk[ind][m] = mm[bm+m][i2*sm+3+j];
               ++ind;
            }							

/*				INTERIOR MODES */	
            ind1 = 0;
            for(m = 2; m< sm+1;++m) {		
               for(n = 1; n < sm+1-m;++n) {
                  vwk[ind] = mm[i*sm+3+k][bm+ind1];
                  for(j=0;j<im;++j) 
                     mwk[ind][j] = mm[bm+j][bm+ind1];
                  ++ind;
                  ++ind1;
               }
               ++ind1;
            }

            GETRF(im,im,mwk[0],MXTM,ipiv,info);
            if (info != 0) {
               printf("DGETRF FAILED SIDE info:%d (sm+2):%d k:%d\n",info,(sm+2),k);
               for(j=0;j<im;++j)
                  vwk[j] = 0.0;
            }
            else
               GETRS(trans,im,1,mwk[0],MXTM,ipiv,vwk,MXTM,info);

            for(j=0;j<im;++j)
               ifmb[i*sm+3+k][j] = vwk[j];
         }

/*			FOR HIGHEST ORDER MODE - STATIC INVERT ALL INTERIOR MODES */
         for(j=0;j<tm;++j) {
            vwk[j] = mm[i*sm +sm +2][j+bm];
            for(k=0;k<im;++k)
               mwk[j][k] = mm[k+bm][j+bm];
         }

         GETRF(im,im,mwk[0],MXTM,ipiv,info);
         if (info != 0) {
            printf("GETRF FAILED info: %d cond: %f\n",info,rcond);
            exit(1);
         }
         GETRS(trans,im,1,mwk[0],MXTM,ipiv,vwk,MXTM,info);

         for(j=0;j<im;++j)
            ifmb[i*sm +sm +2][j] = vwk[j];
      }
   }

/*********************************************/
/*		FIND VERTEX-SIDE & SIDE-SIDE MATRICES  */
/*********************************************/
   for(k=0;k<sm;++k) {
      for(j=0;j<tm;++j)
         u[j] = 0.0;
         
      u[k+3] = 1.0;
      for(j=0;j<im;++j) {
         u[j+bm] = -ifmb[k+3][j];
      }

      proj(u,wk1);
      intgrt(wk1,l);
      
      for(j=0;j<3;++j)
         vfms[j][k] = l[j];
         
      sdiag[k] = 1.0/l[3+k];
      
      for(j=0;j<k;++j) {
         sfms[j][k][0] = l[3+j];
         sfms[j][k][1] = l[3+j+sm];
         sfms[j][k][2] = l[3+j+2*sm];
      }
   }


/************************************************/
/*		FIND MATRICES TO DETERMINE INTERIOR MODES */
/************************************************/
   if (im > 0) {
/*		SETUP DIAGANOL FORM OF MATRIX */
/*		ONLY NECESSARY WHEN USING DPBTRF */	
      for(j=0;j<im;++j) {
         i1 = (0 > j-ibwth ? 0 : j-ibwth);
         for(i=i1;i<=j;++i) {
            k = i - (j-ibwth);
            idiag[j][k] = mm[i+bm][j+bm];
         }
      }
      
      PBTRF(uplo,im,ibwth,idiag[0],ibwth+1,info);
      if (info != 0 || rcond < 10.*EPSILON) {
         printf("PBTRF FAILED (sm+2) : %d info: %d cond: %f\n",(sm+2), info,rcond);
         exit(1);
      }
               
/*		MATRIX TO REMOVE BOUNDARY MODES FROM INTERIOR */
      for(i=0; i<bm; ++i ) {
         for(j=0; j<im; ++j ) {
            bfmi[i][j] = mm[i][j+bm];
         }
         PBTRS(uplo,im,ibwth,1,idiag[0],ibwth+1,bfmi[i],im,info);
      }	
   }


#ifdef SKIP
/*	CHECK TO MAKE SURE PREVIOUS RESULTS ARE RIGHT */
   for(i=0;i<3;++i) {
      i1 = (i+1)%3;
      i2 = (i+2)%3;

      u[i] = 1.0;
      u[i1] = 0.0;
      u[i2] = 0.0;
      
      for(j=0;j<sm;++j) {
         u[3+j+i*sm] = -sfmv[0][j];
         u[3+j+i1*sm] = 0.0;
         u[3+j+i2*sm] = -sfmv[1][j];
      }
      
      for(j=0;j<im;++j)
         u[bm+j] = -ifmb[i][j];
         
      proj(u,wk1);
      intgrt(wk1,l);
      printf("%2d:",i);
      for(j=0;j<tm;++j)
         printf("%+.4le  ",l[j]);
      printf("\n");
   }
   
   for(i=0;i<3;++i) {
      for(k=0;k<sm;++k) {
         for(j=0;j<tm;++j)
            u[j] = 0.0;
         u[i*sm+k+3] = 1.0;
         for(j=0;j<im;++j) {
            u[j+bm] = -ifmb[i*sm+k+3][j];
         }

         proj(u,wk1);
         intgrt(wk1,l);
         
         printf("%2d:",3+i*sm+k);
         for(j=0;j<tm;++j)
            printf("%+.4le  ",l[j]);
         printf("\n");
      }			
   }
#endif	

/********************************************************************/	
/*	NOW SETUP SIMILAR THING FOR 1-D MATRICES 									*/
/********************************************************************/	
/*	CALCULATE 1-D MASS MATRIX */	
   for(m=0;m<sm+2;++m) {
      ind = m + (m > 1 ? 1 : 0);
      for(k=0;k<sm+2;++k) {
         ind1 = k + (k > 1 ? 1 : 0);
         mm[m][k]= 0.0;
         for(i=0;i<gpx;++i)
            mm[m][k] += wtx[i]*gx[ind][i]*gx[ind1][i];
      }
   }

#ifdef SKIP
   printf("\n mm MATRIX (sm+2) = %d\n",(sm+2));
   for(i=0;i<sm+2;++i) {
      printf("%2d:",i);
      for(j=0;j<sm+2;++j)
         printf("%+.4lf ",mm[i][j]);
      printf("\n");
   }
#endif

   if(p > 1) {
/*		STATIC INVERSION SIDE MATRIX  */			
      for(j=0;j<sm;++j) {
         i1 = (0 > j-sbwth ? 0 : j-sbwth);
         for(i=i1;i<=j;++i) {
            k = i - (j-sbwth);
            sdiag1d[j][k] = mm[i+2][j+2];
         }
      }

      PBTRF(uplo,sm,sbwth,sdiag1d[0],sbwth+1,info);
      if (info != 0 || rcond < 100.*EPSILON) {
         printf("PBTRF FAILED - 1D (sm+2) : %d info: %d cond: %f\n",(sm+2), info,rcond);
         exit(1);
      }

/*		MATRIX TO REMOVE VERTEX MODES FROM SIDES */
      for(i=0; i<2; ++i) {
         for(j=0; j<sm; ++j) {
            vfms1d[i][j] = mm[i][j+2];
         }
         PBTRS(uplo,sm,sbwth,1,sdiag1d[0],sbwth+1,bfmi[i],(sm+2)-2,info);
      }

/*		MATRIX TO REMOVE SIDE MODES FROM VERTICES */	
      for(i=0;i<2;++i) {		
/*			VERTEX CONSTRAINT */
         vwk[0] = mm[i][(i+1)%2];
         for(m=0;m<sm;++m)
            mwk[0][m] = mm[(i+1)%2][m+2];

/*			SIDE CONSTRAINTS */
         for(m=0;m<sm-1;++m) {
            vwk[m+1] = mm[i][m+2];
            for(k=0;k<sm;++k) {
               mwk[m+1][k] = mm[m+2][k+2];
            }
         }
         GETRF(sm,sm,mwk[0],MXTM,ipiv,info);
         if (info != 0) {
            printf("DGETRF FAILED - VRTX info:%d (sm+2):%d i:%d\n",info,(sm+2),i);
            exit(1);
         }
         GETRS(trans,sm,1,mwk[0],MXTM,ipiv,vwk,MXTM,info);

/*			STORE MATRIX */
         for(k=0;k<sm;++k)
            sfmv1d[i][k] = vwk[k];
      }

/*		FIND VERTEX DIAGANOL ELEMENT */
      vdiag1d = mm[0][0];
      for(k=0;k<sm;++k)
         vdiag1d -= sfmv1d[0][k]*mm[0][2+k];
   }
   else {
      vdiag1d = mm[0][0]+mm[0][1];
   }

#ifdef SKIP
/*	CHECK TO MAKE SURE PREVIOUS 1D RESULTS ARE RIGHT */
   for(k=0;k<2;++k) {
      u[k] = 1.0;
      u[(k+1)%2] = 0.0;
      for(m=0;m<sm;++m)
         u[m+2] = -sfmv1d[k][m];

      for(i=0;i<sm+2;++i) {
         uht[1][i] = u[0]*gx[0][i];
         uht[1][i] += u[1]*gx[1][i];
         for(m=0;m<sm;++m)
            uht[1][i] += u[m+2]*gx[3+m][i];
      }

      for(m=0;m<sm+2;++m) {
         indm = m + (m > 1 ? 1 : 0);
         l[m] = 0.0;
         for(i=0;i<gpx;++i) {
            l[m] += wtx[i]*gx[indm][i]*uht[1][i];
         }
      }
      printf("%2d:",k);
      for(j=0;j<sm+2;++j)
         printf("%+.4le  ",l[j]);
      printf("%+.4le  ",vdiag1d);
      printf("\n");
   }
#endif

	return;
}			

void hpbasis::legpt()
{
	FLT x,eta,r,s, pkp, pk, pkm;
	int i,j,m,n,ind;
   
   mat_alloc(lgrnge1d,nmodx,sm+1,FLT);
   tens_alloc(lgrnge,tm,sm+1,sm+1,FLT);

/*	CALCULATE PROJECTION POINTS IN INTERIOR */
   for(i=1;i<sm;++i) {
      for(j=1;j<sm-(i-1);++j) {
  			s = -1 +2.0*((FLT) j)/(FLT)(sm+1);
  			r = -1 +2.0*((FLT) i)/(FLT)(sm+1);
  			x = 2.0*(1+r)/(1-s) -1.0;
         eta = s;
  					
/*			CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
  			lgrnge1d[0][i] = .5*(1-x);
			lgrnge1d[1][i] = .5*(1+x);
			lgrnge1d[2][i] = 1.0;

/*			SIDE 1	*/
/*			CALCULATE P, P' USING RECURSION RELATION */
			pk = 1.0;
			pkm = 0.0;
			for (m = 1;m < sm+1;++m) {
				lgrnge1d[m+2][i] = (1.+x)*(1.-x)*.25*pk*norm[m+2];
				pkp = (x-a0[0][m-1])*pk - b0[0][m-1]*pkm;
				pkm = pk;
				pk = pkp;
  			}
  	
		
/*			CALCULATE S POLYNOMIALS */
			ind = 0;
/*			VERTEX A	*/
			lgrnge[ind++][i][j] = (1-eta)*.5*lgrnge1d[0][i];

/*	 		VERTEX B  */
			lgrnge[ind++][i][j] = (1-eta)*.5*lgrnge1d[1][i];

/*			VERTEX C	 */	
			lgrnge[ind++][i][j] = (1+eta)*.5*lgrnge1d[2][i];

/*	  		SIDE 1 (s)		*/
			for(m = 2; m <= sm+1; ++m) {
				lgrnge[ind++][i][j] = pow(.5*(1-eta),m)*lgrnge1d[m+1][i];
			}

/*			SIDE 2	*/
			pk = 1.0;
			pkm = 0.0;
			for(n=1;n<=sm;++n) {
				lgrnge[ind++][i][j] = (1.-eta)*(1.+eta)*.25*pk*lgrnge1d[1][i]*norm[n+2];
				pkp = (eta-a0[1][n-1])*pk - b0[1][n-1]*pkm;
				pkm = pk;
				pk = pkp;
  			}

/*	 		SIDE 3	*/
			pk = 1.0;
			pkm = 0.0;
			for( n = 1; n<= sm;++n) {
				lgrnge[ind++][i][j] = (n % 2 ? 1 : -1)*
				(1.-eta)*(1.+eta)*.25*pk*lgrnge1d[0][i]*norm[n+2];
				pkp = (eta-a0[1][n-1])*pk - b0[1][n-1]*pkm;
				pkm = pk;
				pk = pkp;
  			}

/*	  		INTERIOR MODES	*/
			ind = bm;
			for(m = 2; m< sm+1;++m) {		
				pk = 1.0;
				pkm = 0.0;
				for(n = 1; n < sm+2-m;++n) {
					lgrnge[ind][i][j] = pow(.5*(1.-eta),m)*.5*(1.+eta)*pk*norm[ind]*lgrnge1d[m+1][i];
					pkp = (eta-a0[m][n-1])*pk - b0[m][n-1]*pkm;
					pkm = pk;
					pk = pkp;
					++ind;
				}
			}
		}
	}
	
/*	NOW CALCULATE VALUES OF G FOR SIDE PROJECTION */
	for (i=1;i<sm+1;++i) {
		x = 2.0*(FLT) i/(FLT)(sm+1) -1.0;
		lgrnge1d[0][i] = .5*(1-x);
		lgrnge1d[1][i] = .5*(1+x);

/*		SIDE 1	*/
/*		CALCULATE P, P' USING RECURSION RELATION */
		pk = 1.0;
		pkm = 0.0;
		for (m = 1;m < sm+1;++m) {
			lgrnge1d[m+1][i] = (1.+x)*(1.-x)*.25*pk*norm[m+2];
			pkp = (x-a0[0][m-1])*pk - b0[0][m-1]*pkm;
			pkm = pk;
			pk = pkp;
  		}
  	}
		
	return;
}


