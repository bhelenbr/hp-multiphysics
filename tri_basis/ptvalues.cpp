/*
 *  ptvalues.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on Tue Apr 16 2002.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "hpbasis.h"
#include<math.h>

/* CALCULATES VALUES OF GX POLYNOMIALS & GS POLYNOMIALS AT POINT */
void hpbasis::ptvalues(FLT x, FLT eta) {
    FLT pkp,pk,pkm;
    int k,m,ind;
   
   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = 1.0;
#ifndef VERTEX
   if (p) {
#endif
      pgx[1] = .5*(1-x);
      pgx[2] = .5*(1+x);
#ifndef VERTEX
   }
#endif

   /* SIDE 0   */
   /* CALCULATE P, P' USING RECURSION RELATION */
   pk = 1.0;
   pkm = 0.0;
   for (m = 0;m < sm;++m) {
      pgx[m+3] = pgx[1]*pgx[2]*pk*norm[m+3];
      pkp = (x-a0[0][m])*pk - b0[0][m]*pkm;
      pkm = pk;
      pk = pkp;
   }


   /* CALCULATE S POLYNOMIALS */
   ind = 0;
   
   /* VERTEX 0,1,2   */
#ifdef VERTEX
   pgn[ind++] = (1+eta)*.5;
   pgn[ind++] = (1-eta)*.5;
   pgn[ind++] = (1-eta)*.5;
#else
   pgn[ind++] = 1.0;
   if (p) {
      pgn[ind++] = (1-eta)*.5;
      pgn[ind++] = (1-eta)*.5;
   }
#endif

   /* SIDE 0 (s)      */
   for(m = 2; m <= sm+1; ++m)
      pgn[ind++] = pow(.5*(1-eta),m);

   /* SIDE 1   */
   pk = 1.0;
   pkm = 0.0;
   for(m=0;m<sm;++m) {
      pgn[ind++] = (1.-eta)*(1.+eta)*.25*pk*norm[m+3];
      pkp = (eta-a0[1][m])*pk - b0[1][m]*pkm;
      pkm = pk;
      pk = pkp;
   }

   /* SIDE 2   */
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
      }
   }
   
   return;
}

/* CALCULATE VALUE OF G(X) & DG/DX, G(eta), DG/Deta AT POINT */
void hpbasis::ptvalues_deriv(FLT x, FLT eta) {
   FLT pkp,pk,pkm,dpk,dpkm,dpkp;
   int k,m,n,ind;
   
   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = 1.0;
   dpgx[0] = 0.0;
#ifndef VERTEX
   if (p) {
#endif
      pgx[1] = .5*(1-x);
      dpgx[1] = -0.5;
      
      pgx[2] = .5*(1+x);
      dpgx[2] = 0.5;
#ifndef VERTEX
   }
#endif

   /* SIDE 0   */
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

   /* VERTEX 0,1,2   */
#ifdef VERTEX
   pgn[ind] = (1+eta)*.5;
   dpgn[ind++] = 0.5;
   
   pgn[ind] = (1-eta)*.5;
   dpgn[ind++] = -0.5;

   pgn[ind] = (1-eta)*.5;
   dpgn[ind++] = -0.5;
#else
   pgn[ind] = 1.0;
   dpgn[ind++] = 0.0;
   if (p) {
      pgn[ind] = (1-eta)*.5;
      dpgn[ind++] = -0.5;
      pgn[ind] = (1-eta)*.5;
      dpgn[ind++] = -0.5;
   }
#endif

   /* SIDE 0 (s)      */
   for(m = 2; m <= sm+1; ++m) {
      pgn[ind] = pow(.5*(1-eta),m);
      dpgn[ind] = -.5*m*pow(.5*(1.-eta),m-1);
      ++ind;
   }

   /* SIDE 1   */
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

   /* SIDE 2   */
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

   
   /*  INTERIOR MODES   */
   for(m = 2; m< sm+1;++m) {      
      pk = 1.0;
      pkm = 0.0;
      dpk = 0.0;
      dpkm = 0.0;

      for(k = 0; k < sm+1-m;++k) {
         pgn[ind] = (pow(.5*(1.-eta),m)*.5*(1.+eta)*pk)*norm[ind];
         dpgn[ind] = (-.25*m*pow(.5*(1.-eta),m-1)*(1.+eta)*pk +pow(.5*(1.-eta),m)*.5*(pk + (1.+eta)*dpk))*norm[ind];
         pkp = (eta-a0[m][k])*pk - b0[m][k]*pkm;
         dpkp = pk + (eta-a0[m][k])*dpk - b0[m][k]*dpkm;
         dpkm = dpk;
         dpk = dpkp;
         pkm = pk;
         pk = pkp;
         ++ind;
      }
   }

   return;
}

/* CALCULATES VALUES OF GX POLYNOMIALS & GS POLYNOMIALS AT POINT */
/* BOUNDARY MODES ONLY */
void hpbasis::ptvalues_bdry(FLT x, FLT eta) {
    FLT pkp,pk,pkm;
    int m,ind;
   
   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = 1.0;
#ifndef VERTEX
   if (p) {
#endif
      pgx[1] = .5*(1-x);
      pgx[2] = .5*(1+x);
#ifndef VERTEX
   }
#endif

   /* SIDE 0   */
   /* CALCULATE P, P' USING RECURSION RELATION */
   pk = 1.0;
   pkm = 0.0;
   for (m = 0;m < sm;++m) {
      pgx[m+3] = pgx[1]*pgx[2]*pk*norm[m+3];
      pkp = (x-a0[0][m])*pk - b0[0][m]*pkm;
      pkm = pk;
      pk = pkp;
   }


   /* CALCULATE S POLYNOMIALS */
   ind = 0;
   
   /* VERTEX 0,1,2   */
#ifdef VERTEX
   pgn[ind++] = (1+eta)*.5;
   pgn[ind++] = (1-eta)*.5;
   pgn[ind++] = (1-eta)*.5;
#else
   pgn[ind++] = 1.0;
   if (p) {
      pgn[ind++] = (1-eta)*.5;
      pgn[ind++] = (1-eta)*.5;
   }
#endif

   /* SIDE 0 (s)      */
   for(m = 2; m <= sm+1; ++m)
      pgn[ind++] = pow(.5*(1-eta),m);

   /* SIDE 1   */
   pk = 1.0;
   pkm = 0.0;
   for(m=0;m<sm;++m) {
      pgn[ind++] = (1.-eta)*(1.+eta)*.25*pk*norm[m+3];
      pkp = (eta-a0[1][m])*pk - b0[1][m]*pkm;
      pkm = pk;
      pk = pkp;
   }

   /* SIDE 2   */
   pk = 1.0;
   pkm = 0.0;
   for(m = 0;m<sm;++m) {
      pgn[ind++] = (m % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk*norm[m+3];
      pkp = (eta-a0[1][m])*pk - b0[1][m]*pkm;
      pkm = pk;
      pk = pkp;
   }
   
   return;
}

/* CALCULATE VALUE OF G(X) & DG/DX, G(eta), DG/Deta AT POINT FOR ONLY BOUNDARY MODES */
void hpbasis::ptvalues_deriv_bdry(FLT x, FLT eta) {
   FLT pkp,pk,pkm,dpk,dpkm,dpkp;
   int m,n,ind;

   /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
   pgx[0] = 1.0;
   dpgx[0] = 0.0;
#ifndef VERTEX
   if (p) {
#endif
      pgx[1] = .5*(1-x);
      dpgx[1] = -0.5;
      
      pgx[2] = .5*(1+x);
      dpgx[2] = 0.5;
#ifndef VERTEX
   }
#endif

   /* SIDE 0   */
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

   /* VERTEX 0,1,2   */
#ifdef VERTEX
   pgn[ind] = (1+eta)*.5;
   dpgn[ind++] = 0.5;
   
   pgn[ind] = (1-eta)*.5;
   dpgn[ind++] = -0.5;

   pgn[ind] = (1-eta)*.5;
   dpgn[ind++] = -0.5;
#else
   pgn[ind] = 1.0;
   dpgn[ind++] = 0.0;
   if (p) {
      pgn[ind] = (1-eta)*.5;
      dpgn[ind++] = -0.5;
      pgn[ind] = (1-eta)*.5;
      dpgn[ind++] = -0.5;
   }
#endif

   /* SIDE 0 (s)      */
   for(m = 2; m <= sm+1; ++m) {
      pgn[ind] = pow(.5*(1-eta),m);
      dpgn[ind] = -.5*m*pow(.5*(1.-eta),m-1);
      ++ind;
   }

   /* SIDE 1   */
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

   /* SIDE 2   */
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
#ifdef VERTEX
   pgx[0] = .5*(1-x);
   pgx[1] = .5*(1+x);
#else
   pgx[0] = 1.0;
   if (!p) return;
   pgx[1] = x;
#endif
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

void hpbasis::ptvalues1d_deriv(FLT x) {
    FLT pkp,pk,pkm,dpkp,dpk,dpkm;
    int m;

#ifdef VERTEX
   pgx[0]  = .5*(1-x);
   dpgx[0] = -.5;
   
   pgx[1]  = .5*(1+x);
   dpgx[1] = .5;
#else
   pgx[0]  = 1.0;
   dpgx[0] = 0.0;
   
   if (!p) return;
   
   pgx[1]  = x;
   dpgx[1] = 1.0;
#endif

   /* SIDE 1 CALCULATE P, P' USING RECURSION RELATION */
   pk = 1.0;
   pkm = 0.0;
   dpk = 0.0;
   dpkm = 0.0;

   for (m = 0;m < sm;++m) {
      pgx[m+2] = (1.+x)*(1.-x)*.25*pk*norm[m+3];
      dpgx[m+2] = (-x*.5*pk +(1.+x)*(1.-x)*.25*dpk)*norm[m+3];
      pkp = (x-a0[0][m])*pk - b0[0][m]*pkm;
      dpkp = pk + (x-a0[0][m])*dpk - b0[0][m]*dpkm;
      dpkm = dpk;
      dpk = dpkp;
      pkm = pk;
      pk = pkp;
   }
   
   return;
}
