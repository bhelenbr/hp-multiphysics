/*
 *  proj.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 08 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "hpbasis.h"
#include<stdio.h>

void hpbasis::proj(FLT *lin, FLT **f, FLT **dx, FLT **dy) {
    int i,j,m,n,ind;
    FLT xp1,oeta;
         
   /* DETERMINE U VALUES, GRAD U VALUES
   AT COLLOCATION POINTS
   SUM HAT(U) FOR DU/DX AND DU/DY
   */

   /* GENERAL FORMULA */
   /* dg/dr = 2.0/(1-n) g(s) dg/dx */
   /* dg/ds = g(x)dg/dn +(1+x)/2 dg/dr = g(x)dg/dn +(1+x)/(1-s) g(s) dg/dx */


   /* PART I - sum u*g_mn for each n, s_j   */
   for(j = 0; j < gpn; ++j ) {
      oeta = n0[j];
      
      /* SUM FOR N=1      */
      /* VERTEX A         */
      wk0[0][j] = lin[0]*gn[0][j];
      wk1[0][j] = lin[0]*dgn[0][j];

      /* SIDE 3      */
      for(m = 2*sm+3; m < bm; ++m ) {
         wk0[0][j] += lin[m]*gn[m][j];
         wk1[0][j] += lin[m]*dgn[m][j];
      }         
      wk2[0][j] = wk0[0][j]*oeta;
         
  
      /* SUM FOR N=2      */
      /* VERTEX B         */
      wk0[1][j] = lin[1]*gn[1][j];
      wk1[1][j] = lin[1]*dgn[1][j];

      /* SIDE 2      */
      for (m = sm+3; m < 2*sm+3; ++m) {
         wk0[1][j] += lin[m]*gn[m][j];
         wk1[1][j] += lin[m]*dgn[m][j];
      }         
       wk2[1][j] = wk0[1][j]*oeta;

      /* SUM FOR N=3      */
      /* VERTEX C         */
      wk0[2][j] = lin[2]*gn[2][j];
      wk1[2][j] = lin[2]*dgn[2][j];
      wk2[2][j] = wk0[2][j]*oeta;

      ind = bm;
      for(m = 3; m < sm+3; ++m) {
         /* SIDE 1      */
         wk0[m][j] = lin[m]*gn[m][j];
         wk1[m][j] = lin[m]*dgn[m][j];
      
         /* INTERIOR MODES      */
         for(n = 0; n < sm+2-m; ++n) {
            wk0[m][j] += lin[ind]*gn[ind][j];
            wk1[m][j] += lin[ind]*dgn[ind][j];
            ++ind;
         }
          wk2[m][j] = wk0[m][j]*oeta;
      }
   }
       
   /* SUM OVER N AT EACH I,J POINT   */     
    for (i=0; i < gpx; ++i ) {
      xp1 = x0[i];
       for (j=0; j < gpn; ++j) {
          f[i][j]   = 0.0;
         dx[i][j] = 0.0;
         dy[i][j] = 0.0;

         for(n=0; n < nmodx; ++n ) {     
            f[i][j]   += wk0[n][j]*gx[n][i];
            dy[i][j] += wk1[n][j]*gx[n][i];
            dx[i][j] += wk2[n][j]*dgx[n][i];
         }
         dy[i][j] += xp1*dx[i][j];
      }
   }

   return;
}


void hpbasis::proj(FLT *lin, FLT **f) {
    int i,j,m,n,ind;
         
   /* DETERMINE U VALUES, GRAD U VALUES
   AT COLLOCATION POINTS
   SUM HAT(U) FOR DU/DX AND DU/DY
   */

   /* PART I - sum u*g_mn for each n, s_j   */
   for(j = 0; j < gpn; ++j ) {
      
      /* SUM FOR N=1      */
      /* VERTEX A         */
      wk0[0][j] = lin[0]*gn[0][j];

      /* SIDE 3      */
      for(m = 2*sm+3; m < bm; ++m )
         wk0[0][j] += lin[m]*gn[m][j];
  
      /* SUM FOR N=2      */
      /* VERTEX B         */
      wk0[1][j] = lin[1]*gn[1][j];

      /* SIDE 2      */
      for (m = sm+3; m < 2*sm+3; ++m)
         wk0[1][j] += lin[m]*gn[m][j];

      /* SUM FOR N=3      */
      /* VERTEX C         */
      wk0[2][j] = lin[2]*gn[2][j];

      /* LOOP FOR INTERIOR MODES      */
      ind = bm;
      for(m = 3; m < sm+3; ++m) {
         /* SIDE 1      */
         wk0[m][j] = lin[m]*gn[m][j];
      
         /* INTERIOR MODES      */
         for(n = 0; n < sm+2-m; ++n) {
            wk0[m][j] += lin[ind]*gn[ind][j];
            ++ind;
         }
      }
   }
       
   /* SUM OVER N AT EACH I,J POINT   */     
    for (i=0; i < gpx; ++i ) {
       for (j=0; j < gpn; ++j) {
          f[i][j]   = 0.0;

         for(n=0; n < nmodx; ++n )  
            f[i][j]   += wk0[n][j]*gx[n][i];
      }
   }

   return;
}

void hpbasis::proj(FLT u1, FLT u2, FLT u3, FLT **f) {
    int i,j;
   
   for (i=0; i < gpx; ++i ) 
      for (j=0; j < gpn; ++j ) 
           f[i][j] = u1*gx[0][i]*gn[0][j] +u2*gx[1][i]*gn[1][j] +u3*gx[2][i]*gn[2][j];
        
   return;
}

void hpbasis::proj_bdry(FLT *lin, FLT **f, FLT **dx, FLT **dy) {
    int i,j,m,n;
    FLT xp1,oeta;
         
   /* DETERMINE U VALUES, GRAD U VALUES
   AT COLLOCATION POINTS
   SUM HAT(U) FOR DU/DX AND DU/DY
   */

   /* PART I - sum u*g_mn for each n, s_j   */
   for(j = 0; j < gpn; ++j ) {
      oeta = n0[j];
      
      /* SUM FOR N=1      */
      /* VERTEX A         */
      wk0[0][j] = lin[0]*gn[0][j];
      wk1[0][j] = lin[0]*dgn[0][j];

      /* SIDE 3      */
      for(m = 2*sm+3; m < bm; ++m ) {
         wk0[0][j] += lin[m]*gn[m][j];
         wk1[0][j] += lin[m]*dgn[m][j];
      }         
      wk2[0][j] = wk0[0][j]*oeta;
         
  
      /* SUM FOR N=2      */
      /* VERTEX B         */
      wk0[1][j] = lin[1]*gn[1][j];
      wk1[1][j] = lin[1]*dgn[1][j];

      /* SIDE 2      */
      for (m = sm+3; m < 2*sm+3; ++m) {
         wk0[1][j] += lin[m]*gn[m][j];
         wk1[1][j] += lin[m]*dgn[m][j];
      }         
       wk2[1][j] = wk0[1][j]*oeta;

      /* SUM FOR N=3      */
      /* VERTEX C         */
      wk0[2][j] = lin[2]*gn[2][j];
      wk1[2][j] = lin[2]*dgn[2][j];
      wk2[2][j] = wk0[2][j]*oeta;
      
      /* LOOP FOR SIDE 1 MODES      */
      for(m = 3; m < sm+3; ++m) {
         wk0[m][j] = lin[m]*gn[m][j];
         wk1[m][j] = lin[m]*dgn[m][j];
         wk2[m][j] = wk0[m][j]*oeta;
      }
   }
       
   /* SUM OVER N AT EACH I,J POINT   */     
    for (i=0; i < gpx; ++i ) {
      xp1 = x0[i];
       for (j=0; j < gpn; ++j) {
          f[i][j]   = 0.0;
         dx[i][j] = 0.0;
         dy[i][j] = 0.0;

         for(n=0; n < nmodx; ++n ) {     
            f[i][j]   += wk0[n][j]*gx[n][i];
            dy[i][j] += wk1[n][j]*gx[n][i];
            dx[i][j] += wk2[n][j]*dgx[n][i];
         }
         dy[i][j] += xp1*dx[i][j];
      }
   }

   return;
}

void hpbasis::proj_bdry(FLT *lin, FLT **f) {
    int i,j,m,n;
         
   /* DETERMINE U VALUES, GRAD U VALUES
   AT COLLOCATION POINTS
   SUM HAT(U) FOR DU/DX AND DU/DY
   */

   /* PART I - sum u*g_mn for each n, s_j   */
   for(j = 0; j < gpn; ++j ) {
      
      /* SUM FOR N=1      */
      /* VERTEX A         */
      wk0[0][j] = lin[0]*gn[0][j];

      /* SIDE 3      */
      for(m = 2*sm+3; m < bm; ++m )
         wk0[0][j] += lin[m]*gn[m][j];
  
      /* SUM FOR N=2      */
      /* VERTEX B         */
      wk0[1][j] = lin[1]*gn[1][j];

      /* SIDE 2      */
      for (m = sm+3; m < 2*sm+3; ++m)
         wk0[1][j] += lin[m]*gn[m][j];

      /* SUM FOR N=3      */
      /* VERTEX C         */
      wk0[2][j] = lin[2]*gn[2][j];

      /* SIDE 1      */
      for(m = 3; m < sm+3; ++m)
         wk0[m][j] = lin[m]*gn[m][j];
   }
       
   /* SUM OVER N AT EACH I,J POINT   */     
    for (i=0; i < gpx; ++i ) {
       for (j=0; j < gpn; ++j) {
          f[i][j]   = 0.0;

         for(n=0; n < nmodx; ++n )  
            f[i][j]   += wk0[n][j]*gx[n][i];
      }
   }

   return;
}



void hpbasis::derivr(FLT **f, FLT **dx) {
    int i,j,n;
    FLT uddlt[MXTM];

   for (j=0;j<gpn;++j) {
      for (i=0;i<gpx;++i)
         uddlt[i] = f[i][j]*dltx[i]*n0[j];
      for (i=0;i<gpx;++i) {
         for(n=0;n<i;++n) 
            dx[i][j] += (uddlt[i] + uddlt[n])*dltx1[i][n];
         for(n=i+1;n<gpx;++n) 
            dx[i][j] += (uddlt[i] + uddlt[n])*dltx1[i][n];
      }
   }

   return;
}

void hpbasis::derivs(FLT **f, FLT **dy) {
    int i,j,n;
    FLT uddlt[MXTM];
   
   for (j=0;j<gpn;++j) {
      for (i=0;i<gpx;++i)
         uddlt[i] = f[i][j]*dltx[i]*n0[j];
      for (i=0;i<gpx;++i) {
         for(n=0;n<i;++n) 
            dy[i][j] += (uddlt[i] + uddlt[n])*dltn1[i][n];
         for(n=i+1;n<gpx;++n) 
            dy[i][j] += (uddlt[i] + uddlt[n])*dltn1[i][n];
      }
   }

   for (i=0;i<gpx;++i) {
      for (j=0;j<gpn;++j)
         uddlt[j] = f[i][j]*dltn[j];
      for (j=0;j<gpn;++j) {
         for(n=0;n<j;++n) 
            dy[i][j] += (uddlt[j] + uddlt[n])*dltn2[j][n];
         for(n=j+1;n<gpn;++n) 
            dy[i][j] += (uddlt[j] + uddlt[n])*dltn2[j][n];
      }
   }

   return;
}
   
void hpbasis::proj_leg(FLT *lin, FLT **f) {
   int i,j,k;
     
   /* INTERIOR */
   for(i=1;i<sm;++i) {
      for(j=1;j<sm-(i-1);++j) {
         f[i][j] = 0.0;
         for(k=0;k<tm;++k)
            f[i][j] += lin[k]*lgrnge[k][i][j];
      }
   }
   
   return;
}

void hpbasis::proj_leg(FLT u1, FLT u2, FLT u3, FLT **f) {
   int i,j;
     
   /* INTERIOR */
   for(i=1;i<sm;++i)
      for(j=1;j<sm-(i-1);++j)
         f[i][j] = u1*lgrnge[0][i][j] +u2*lgrnge[1][i][j] +u3*lgrnge[2][i][j];
   
   return;
}

void hpbasis::proj_bdry_leg(FLT *lin, FLT **f) {
   int i,j,k;

   for(i=1;i<sm;++i) {
      for(j=1;j<sm-(i-1);++j) {
         f[i][j] = 0.0;
         for(k=0;k<bm;++k)
            f[i][j] += lin[k]*lgrnge[k][i][j];
      }
   }
   
   return;
}

void hpbasis::proj_side(int side, FLT *lin, FLT *f, FLT *dtan, FLT *dnrm) {
   int i,m;

   /* NORMAL DERIVATIVE */
   for(i=0;i<gpx;++i) {
      dnrm[i] = 0.0;
      for(m=0;m<tm;++m)
         dnrm[i] += lin[m]*dgnorm[side][m][i];
   }

/*	MOVE SIDE & VERTEX MODES TO WK & THEN PROJECT VALUES & TANGENT DERIVS */   
   wk0[0][0] = lin[side];
   wk0[0][1] = lin[(side+1)%3];
   for(m=0;m<sm;++m)
      wk0[0][m+2] = lin[side*sm +3 +m];
      
   proj1d(wk0[0],f,dtan);
      
   return;
}
