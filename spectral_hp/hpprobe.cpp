/*
 *  probe.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 18 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"hpbasis.h"
#include<math.h>

/* RETURNS VALUES OF GX POLYNOMIALS & GS POLYNOMIALS AT POINT */
void hpbasis::ptvalues(FLT *pgx, FLT *pgn, FLT x, FLT eta) {
	static FLT pkp,pk,pkm;
	static int k,m,ind;
   
/*	CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = .5*(1-x);
   pgx[1] = .5*(1+x);
   pgx[2] = 1.0;

/*	SIDE 1	*/
/*	CALCULATE P, P' USING RECURSION RELATION */
   pk = 1.0;
   pkm = 0.0;
   for (m = 0;m < sm;++m) {
      pgx[m+3] = (1.+x)*(1.-x)*.25*pk*norm[m+3];
      pkp = (x-a0[0][m])*pk - b0[0][m]*pkm;
      pkm = pk;
      pk = pkp;
   }


/*	CALCULATE S POLYNOMIALS */
   ind = 0;
   
/*	VERTEX 0,1,2	*/
   pgn[ind++] = (1-eta)*.5;
   pgn[ind++] = (1-eta)*.5;
   pgn[ind++] = (1+eta)*.5;

/*	SIDE 1 (s)		*/
   for(m = 2; m <= sm+1; ++m)
      pgn[ind++] = pow(.5*(1-eta),m);

/*	SIDE 2	*/
   pk = 1.0;
   pkm = 0.0;
   for(m=0;m<sm;++m) {
      pgn[ind++] = (1.-eta)*(1.+eta)*.25*pk*norm[m+3];
      pkp = (eta-a0[1][m])*pk - b0[1][m]*pkm;
      pkm = pk;
      pk = pkp;
   }

/*	SIDE 3	*/
   pk = 1.0;
   pkm = 0.0;
   for(m = 0;m<sm;++m) {
      pgn[ind++] = (m % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk*norm[m+3];
      pkp = (eta-a0[1][m])*pk - b0[1][m]*pkm;
      pkm = pk;
      pk = pkp;
   }

/*	INTERIOR MODES	*/
   ind = bm;
   for(m = 2; m< sm+1;++m) {		
      pk = 1.0;
      pkm = 0.0;
      for(k = 0; k < sm+1-m;++k) {
         pgn[ind++] = pow(.5*(1.-eta),m)*.5*(1.+eta)*pk*norm[ind];
         pkp = (eta-a0[m][k])*pk - b0[m][k]*pkm;
         pkm = pk;
         pk = pkp;
         ++ind;
      }
   }
   
   return;
}

/*	CALCULATE VALUE OF G(X) & DG/DX, G(eta), DG/Deta AT POINT FOR ONLY BOUNDARY MODES */
void hpbasis::ptvalues_bdry(FLT **pgx, FLT **pgn, FLT x, FLT eta) {
	FLT pkp,pk,pkm,dpk,dpkm,dpkp;
	int m,n,ind;

/*	CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0][0] = .5*(1-x);
   pgx[1][0] = -0.5;
   
   pgx[0][1] = .5*(1+x);
   pgx[1][1] = 0.5;
   
   pgx[0][2] = 1.0;
   pgx[1][2] = 0.0;

/*	SIDE 1	*/
/*	CALCULATE P, P' USING RECURSION RELATION */
   pk = 1.0;
   pkm = 0.0;
   dpk = 0.0;
   dpkm = 0.0;
   for (m = 0;m < sm;++m) {
      pgx[0][m+3] = (1.+x)*(1.-x)*.25*pk*norm[m+3];
      pgx[1][m+3] = (-x*.5*pk +(1.+x)*(1.-x)*.25*dpk)*norm[m+2];
      pkp = (x-a0[0][m])*pk - b0[0][m]*pkm;
      dpkp = pk + (x-a0[0][m])*dpk - b0[0][m]*dpkm;
      dpkm = dpk;
      dpk = dpkp;
      pkm = pk;
      pk = pkp;
   }

/******************************************/
/* GENERATE JACOBI POLY FOR S DIRECTION */
/****************************************/
   ind = 0;

/*	VERTEX A	*/
   pgn[0][ind] = (1-eta)*.5;
   pgn[1][ind] = -.5;
   ++ind;
   
/*	VERTEX B  */
   pgn[0][ind] = (1-eta)*.5;
   pgn[1][ind] = -.5;
   ++ind;

/*	VERTEX C	 */	
   pgn[0][ind] = (1+eta)*.5;
   pgn[1][ind] = .5;
   ++ind;

/*	SIDE 1 (s)		*/
   for(m = 2; m <= sm+1; ++m) {
      pgn[0][ind] = pow(.5*(1-eta),m);
      pgn[1][ind] = -.5*m*pow(.5*(1.-eta),m-1);
      ++ind;
   }

/*	SIDE 2	*/
   pk = 1.0;
   pkm = 0.0;
   dpk = 0.0;
   dpkm = 0.0;
   for(n=0;n<sm;++n) {
      pgn[0][ind] = (1.-eta)*(1.+eta)*.25*pk*norm[n+3];;		
      pgn[1][ind] = -.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk;
      pkp = (eta-a0[1][n])*pk - b0[1][n]*pkm;
      dpkp = pk + (eta-a0[1][n])*dpk - b0[1][n]*dpkm;
      dpkm = dpk;
      dpk = dpkp;
      pkm = pk;
      pk = pkp;
      ++ind;
   }

/*	SIDE 3	*/
   pk = 1.0;
   pkm = 0.0;
   dpk = 0.0;
   dpkm = 0.0;

   for(n=0;n<sm;++n) {
      pgn[0][ind] = (n % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk;
      pgn[1][ind] = (n % 2 ? -1 : 1)*(-.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk);
      pkp = (eta-a0[1][n])*pk - b0[1][n]*pkm;
      dpkp = pk + (eta-a0[1][n])*dpk - b0[1][n]*dpkm;
      dpkm = dpk;
      dpk = dpkp;
      pkm = pk;
      pk = pkp;
      ++ind;
   }

   return;
}

void hpbasis::ptvalues1d(FLT *pgx, FLT x) {
	static FLT pkp,pk,pkm;
	static int m;
   
/*	CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = .5*(1-x);
   pgx[1] = .5*(1+x);

/*	SIDE 1	*/
/*	CALCULATE P, P' USING RECURSION RELATION */
   pk = 1.0;
   pkm = 0.0;
   for (m = 0;m < sm;++m) {
      pgx[m+2] = (1.+x)*(1.-x)*.25*pk*norm[m+3];
      pkp = (x-a0[0][m])*pk - b0[0][m]*pkm;
      pkm = pk;
      pk = pkp;
   }
   
   return;
}


void hpbasis::ptprobe(int nv, FLT **lin, FLT *f, FLT r, FLT s, FLT *pgx, FLT *pgn) {
   static int k,m,n,ind;
   static FLT x,eta;
   
   x = 2.0*(1+r)/(1-s) -1.0;
   eta = s;
   
   ptvalues(pgx,pgn,x,eta);
   
   for(n=0;n<nv;++n) {
   
/* 	SUM ALL S MODE CONTRIBUTIONS */
/*		VERTEX A			*/
		wk0[0][0] = lin[n][0]*pgn[0];

/*		SIDE 3		*/
		for(m = 2*sm+3; m < bm; ++m )
			wk0[0][0] += lin[n][m]*pgn[m];
  
/*		SUM FOR N=2		*/
/*		VERTEX B			*/
		wk0[1][0] = lin[n][1]*pgn[1];

/*		SIDE 2		*/
		for (m = sm+3; m < 2*sm+3; ++m)
			wk0[1][0] += lin[n][m]*pgn[m];

/*		SUM FOR N=3		*/
/*		VERTEX C			*/
		wk0[2][0] = lin[n][2]*pgn[2];

/*		LOOP FOR INTERIOR MODES		*/
		ind = bm;
		for(m = 3; m < sm+3; ++m) {
/*			SIDE 1		*/
			wk0[m][0] = lin[n][m]*pgn[m];
		
/*			INTERIOR MODES		*/
			for(k = 0; k < sm+2-m; ++k) {
				wk0[m][0] += lin[n][ind]*pgn[ind];
				++ind;
			}
		}
	 	
/*	SUM OVER N X MODES	*/  	
      f[n]	= 0.0;

      for(k=0; k < nmodx; ++k )  
         f[n]	+= wk0[k][0]*pgx[k];
   }
   return;
}

void hpbasis::ptprobe_bdry(int nv, FLT **lin, FLT *f, FLT *dx, FLT *dy, FLT r, FLT s, FLT **pgx, FLT **pgn) {
   static int k,m,n;
   static FLT n0,x0,x,eta;
   
   x = 2.0*(1+r)/(1-s) -1.0;
   eta = s;
   
   ptvalues_bdry(pgx,pgn,x,eta);
   
   for(n=0;n<nv;++n) {

/*		PART I - sum u*g_mn for each n, s_j	*/
		n0 = 2./(1-eta);
		
/*		SUM FOR N=1		*/
/*		VERTEX A			*/
		wk0[0][0] = lin[n][0]*pgn[0][0];
		wk1[0][0] = lin[n][0]*pgn[1][0];

/*		SIDE 3		*/
		for(m = 2*sm+3; m < bm; ++m ) {
			wk0[0][0] += lin[n][m]*pgn[0][m];
			wk1[0][0] += lin[n][m]*pgn[1][m];
		}			
		wk2[0][0] = wk0[0][0]*n0;
			
  
/*		SUM FOR N=2		*/
/*		VERTEX B			*/
		wk0[1][0] = lin[n][1]*pgn[0][1];
		wk1[1][0] = lin[n][1]*pgn[1][1];

/*		SIDE 2		*/
		for (m = sm+3; m < 2*sm+3; ++m) {
			wk0[1][0] += lin[n][m]*pgn[0][m];
			wk1[1][0] += lin[n][m]*pgn[1][m];
		}			
 		wk2[1][0] = wk0[1][0]*n0;

/*		SUM FOR N=3		*/
/*		VERTEX C			*/
		wk0[2][0] = lin[n][2]*pgn[0][2];
		wk1[2][0] = lin[n][2]*pgn[1][2];
		wk2[2][0] = wk0[2][0]*n0;

		for(m = 3; m < sm+3; ++m) {
/*			SIDE 1		*/
			wk0[m][0] = lin[n][m]*pgn[0][m];
			wk1[m][0] = lin[n][m]*pgn[1][m];
	 		wk2[m][0] = wk0[m][0]*n0;
		}
	 	
/*		SUM OVER N AT EACH I,J POINT	*/  	
      x0 = 0.5*(1+x);
      f[n]	= 0.0;
      dx[n] = 0.0;
      dy[n] = 0.0;

      for(k=0; k < nmodx; ++k ) {	  
         f[n]	+= wk0[k][0]*pgx[0][k];
         dy[n] += wk1[k][0]*pgx[0][k];
         dx[n] += wk2[k][0]*pgx[1][k];
      }
      dy[n] += x0*dx[n];
	}
   return;
}

void hpbasis::ptprobe1d(int nv, FLT **lin, FLT *f, FLT x, FLT *pgx) {
   static int k,n;
    
   ptvalues1d(pgx,x);
/*	SUM OVER N X MODES	*/  	
   
   for(n=0;n<nv;++n) {
      f[n]	= 0.0;

      for(k=0; k < sm+2; ++k )  
         f[n]	+= lin[n][k]*pgx[k];
   }
   
   return;
}
   

   
