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
#include<utilities.h>

/* RETURNS VALUES OF GX POLYNOMIALS & GS POLYNOMIALS AT POINT */
void hpbasis::ptvalues(FLT x, FLT eta) {
    FLT pkp,pk,pkm;
    int k,m,ind;
   
   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = .5*(1-x);
   pgx[1] = .5*(1+x);
   pgx[2] = 1.0;

   /* SIDE 1   */
   /* CALCULATE P, P' USING RECURSION RELATION */
   pk = 1.0;
   pkm = 0.0;
   for (m = 0;m < sm;++m) {
      pgx[m+3] = (1.+x)*(1.-x)*.25*pk*norm[m+3];
      pkp = (x-a0[0][m])*pk - b0[0][m]*pkm;
      pkm = pk;
      pk = pkp;
   }


   /* CALCULATE S POLYNOMIALS */
   ind = 0;
   
   /* VERTEX 0,1,2   */
   pgn[ind++] = (1-eta)*.5;
   pgn[ind++] = (1-eta)*.5;
   pgn[ind++] = (1+eta)*.5;

   /* SIDE 1 (s)      */
   for(m = 2; m <= sm+1; ++m)
      pgn[ind++] = pow(.5*(1-eta),m);

   /* SIDE 2   */
   pk = 1.0;
   pkm = 0.0;
   for(m=0;m<sm;++m) {
      pgn[ind++] = (1.-eta)*(1.+eta)*.25*pk*norm[m+3];
      pkp = (eta-a0[1][m])*pk - b0[1][m]*pkm;
      pkm = pk;
      pk = pkp;
   }

   /* SIDE 3   */
   pk = 1.0;
   pkm = 0.0;
   for(m = 0;m<sm;++m) {
      pgn[ind++] = (m % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk*norm[m+3];
      pkp = (eta-a0[1][m])*pk - b0[1][m]*pkm;
      pkm = pk;
      pk = pkp;
   }

   /* INTERIOR MODES   */
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

/* CALCULATE VALUE OF G(X) & DG/DX, G(eta), DG/Deta AT POINT FOR ONLY BOUNDARY MODES */
void hpbasis::ptvalues_bdry(FLT x, FLT eta) {
   FLT pkp,pk,pkm,dpk,dpkm,dpkp;
   int m,n,ind;

   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = .5*(1-x);
   dpgx[0] = -0.5;
   
   pgx[1] = .5*(1+x);
   dpgx[1] = 0.5;
   
   pgx[2] = 1.0;
   dpgx[2] = 0.0;

   /* SIDE 1   */
   /* CALCULATE P, P' USING RECURSION RELATION */
   pk = 1.0;
   pkm = 0.0;
   dpk = 0.0;
   dpkm = 0.0;
   for (m = 0;m < sm;++m) {
      pgx[m+3] = (1.+x)*(1.-x)*.25*pk*norm[m+3];
      dpgx[m+3] = (-x*.5*pk +(1.+x)*(1.-x)*.25*dpk)*norm[m+3];
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

   /* VERTEX A   */
   pgn[ind] = (1-eta)*.5;
   dpgn[ind] = -.5;
   ++ind;
   
   /* VERTEX B  */
   pgn[ind] = (1-eta)*.5;
   dpgn[ind] = -.5;
   ++ind;

   /* VERTEX C    */   
   pgn[ind] = (1+eta)*.5;
   dpgn[ind] = .5;
   ++ind;

   /* SIDE 1 (s)      */
   for(m = 2; m <= sm+1; ++m) {
      pgn[ind] = pow(.5*(1-eta),m);
      dpgn[ind] = -.5*m*pow(.5*(1.-eta),m-1);
      ++ind;
   }

   /* SIDE 2   */
   pk = 1.0;
   pkm = 0.0;
   dpk = 0.0;
   dpkm = 0.0;
   for(n=0;n<sm;++n) {
      pgn[ind] = (1.-eta)*(1.+eta)*.25*pk*norm[n+3];
      dpgn[ind] = (-.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk)*norm[n+3];
      pkp = (eta-a0[1][n])*pk - b0[1][n]*pkm;
      dpkp = pk + (eta-a0[1][n])*dpk - b0[1][n]*dpkm;
      dpkm = dpk;
      dpk = dpkp;
      pkm = pk;
      pk = pkp;
      ++ind;
   }

   /* SIDE 3   */
   pk = 1.0;
   pkm = 0.0;
   dpk = 0.0;
   dpkm = 0.0;

   for(n=0;n<sm;++n) {
      pgn[ind] = (n % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk*norm[n+3];
      dpgn[ind] = (n % 2 ? -1 : 1)*(-.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk)*norm[n+3];
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

void hpbasis::ptvalues1d(FLT x) {
    FLT pkp,pk,pkm;
    int m;
   
   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = .5*(1-x);
   pgx[1] = .5*(1+x);

   /* SIDE 1   */
   /* CALCULATE P, P' USING RECURSION RELATION */
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

void hpbasis::ptprobe(int nv, FLT **lin, FLT *f) {
   int k,m,n,ind;
   
   for(n=0;n<nv;++n) {
   
      /* SUM ALL S MODE CONTRIBUTIONS */
      /* VERTEX A         */
      wk0[0][0] = lin[n][0]*pgn[0];

      /* SIDE 3      */
      for(m = 2*sm+3; m < bm; ++m )
         wk0[0][0] += lin[n][m]*pgn[m];
  
      /* SUM FOR N=2      */
      /* VERTEX B         */
      wk0[1][0] = lin[n][1]*pgn[1];

      /* SIDE 2      */
      for (m = sm+3; m < 2*sm+3; ++m)
         wk0[1][0] += lin[n][m]*pgn[m];

      /* SUM FOR N=3      */
      /* VERTEX C         */
      wk0[2][0] = lin[n][2]*pgn[2];

      /* LOOP FOR INTERIOR MODES      */
      ind = bm;
      for(m = 3; m < sm+3; ++m) {
         /* SIDE 1      */
         wk0[m][0] = lin[n][m]*pgn[m];
      
         /* INTERIOR MODES      */
         for(k = 0; k < sm+2-m; ++k) {
            wk0[m][0] += lin[n][ind]*pgn[ind];
            ++ind;
         }
      }
       
   /* SUM OVER N X MODES   */     
      f[n]   = 0.0;

      for(k=0; k < nmodx; ++k )  
         f[n]   += wk0[k][0]*pgx[k];
   }
   return;
}

void hpbasis::ptprobe_bdry(int nv, FLT **lin, FLT *f, FLT *dx, FLT *dy, FLT r, FLT s) {
   int k,m,n;
   FLT n0,x0,x,eta;
   
   s = MIN(1.0-10.*EPSILON,s);
   x = 2.0*(1+r)/(1-s) -1.0;
   eta = s;
   
   ptvalues_bdry(x,eta);
   
   for(n=0;n<nv;++n) {

      /* PART I - sum u*g_mn for each n, s_j   */
      n0 = 2./(1 -eta);
      
      /* SUM FOR N=1      */
      /* VERTEX A         */
      wk0[0][0] = lin[n][0]*pgn[0];
      wk1[0][0] = lin[n][0]*dpgn[0];

      /* SIDE 3      */
      for(m = 2*sm+3; m < bm; ++m ) {
         wk0[0][0] += lin[n][m]*pgn[m];
         wk1[0][0] += lin[n][m]*dpgn[m];
      }         
      wk2[0][0] = wk0[0][0]*n0;
         
  
      /* SUM FOR N=2      */
      /* VERTEX B         */
      wk0[1][0] = lin[n][1]*pgn[1];
      wk1[1][0] = lin[n][1]*dpgn[1];

      /* SIDE 2      */
      for (m = sm+3; m < 2*sm+3; ++m) {
         wk0[1][0] += lin[n][m]*pgn[m];
         wk1[1][0] += lin[n][m]*dpgn[m];
      }         
       wk2[1][0] = wk0[1][0]*n0;

      /* SUM FOR N=3      */
      /* VERTEX C         */
      wk0[2][0] = lin[n][2]*pgn[2];
      wk1[2][0] = lin[n][2]*dpgn[2];
      wk2[2][0] = wk0[2][0]*n0;

      for(m = 3; m < sm+3; ++m) {
         /* SIDE 1      */
         wk0[m][0] = lin[n][m]*pgn[m];
         wk1[m][0] = lin[n][m]*dpgn[m];
          wk2[m][0] = wk0[m][0]*n0;
      }

      /* SUM OVER N AT EACH I,J POINT   */     
      x0 = 0.5*(1+x);
      f[n]   = 0.0;
      dx[n] = 0.0;
      dy[n] = 0.0;

      for(k=0; k < nmodx; ++k ) {     
         f[n]   += wk0[k][0]*pgx[k];
         dy[n] += wk1[k][0]*pgx[k];
         dx[n] += wk2[k][0]*dpgx[k];
      }
      dy[n] += x0*dx[n];
   }
   return;
}

void hpbasis::ptprobe1d(int nv, FLT **lin, FLT *f) {
   int k,n;

   for(n=0;n<nv;++n) {
      f[n]   = 0.0;

      for(k=0; k < sm+2; ++k )  
         f[n]   += lin[n][k]*pgx[k];
   }
   
   return;
}
   
