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


void hpbasis::ptprobe(int nv, FLT **lin, FLT *f) {
   int k,m,n,ind;
   
   for(n=0;n<nv;++n) {
   
      /* SUM ALL S MODE CONTRIBUTIONS */
      /* VERTEX 0         */
      wk0[0][0] = lin[n][0]*pgn[0];
      
#ifndef VERTEX
      if (p) {
#endif

      /* VERTEX 1         */
      wk0[1][0] = lin[n][1]*pgn[1];

      /* SIDE 2      */
      for(m = 2*sm+3; m < bm; ++m )
         wk0[1][0] += lin[n][m]*pgn[m];
  
      /* VERTEX 2         */
      wk0[2][0] = lin[n][2]*pgn[2];

      /* SIDE 1      */
      for (m = sm+3; m < 2*sm+3; ++m)
         wk0[2][0] += lin[n][m]*pgn[m];

      /* LOOP FOR INTERIOR MODES      */
      ind = bm;
      for(m = 3; m < sm+3; ++m) {
         /* SIDE 0      */
         wk0[m][0] = lin[n][m]*pgn[m];
      
         /* INTERIOR MODES      */
         for(k = 0; k < sm+2-m; ++k) {
            wk0[m][0] += lin[n][ind]*pgn[ind];
            ++ind;
         }
      }
#ifndef VERTEX
      }
#endif
       
   /* SUM OVER N X MODES   */     
      f[n]   = 0.0;

      for(k=0; k < nmodx; ++k )  
         f[n]   += wk0[k][0]*pgx[k];
   }
   return;
}

void hpbasis::ptprobe_bdry(int nv, FLT **lin, FLT *f) {
   int k,m,n;
   
   for(n=0;n<nv;++n) {
   
      /* SUM ALL S MODE CONTRIBUTIONS */
      /* VERTEX 0         */
      wk0[0][0] = lin[n][0]*pgn[0];
      
#ifndef VERTEX
      if (p) {
#endif

      /* VERTEX 1         */
      wk0[1][0] = lin[n][1]*pgn[1];

      /* SIDE 2      */
      for(m = 2*sm+3; m < bm; ++m )
         wk0[1][0] += lin[n][m]*pgn[m];
  
      /* VERTEX 2         */
      wk0[2][0] = lin[n][2]*pgn[2];

      /* SIDE 1      */
      for (m = sm+3; m < 2*sm+3; ++m)
         wk0[2][0] += lin[n][m]*pgn[m];

      /* SIDE 0      */
      for(m = 3; m < sm+3; ++m) 
         wk0[m][0] = lin[n][m]*pgn[m];
#ifndef VERTEX
      }
#endif       
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
   
   ptvalues_deriv_bdry(x,eta);
   
   for(n=0;n<nv;++n) {

      /* PART I - sum u*g_mn for each n, s_j   */
      n0 = 2./(1 -eta);
      
      /* VERTEX 0         */
      wk0[0][0] = lin[n][0]*pgn[0];
      wk1[0][0] = lin[n][0]*dpgn[0];
      wk2[0][0] = wk0[0][0]*n0;
      
#ifndef VERTEX
      if (p) {
#endif

      /* VERTEX 1         */
      wk0[1][0] = lin[n][1]*pgn[1];
      wk1[1][0] = lin[n][1]*dpgn[1];

      /* SIDE 2      */
      for(m = 2*sm+3; m < bm; ++m ) {
         wk0[1][0] += lin[n][m]*pgn[m];
         wk1[1][0] += lin[n][m]*dpgn[m];
      }         
      wk2[1][0] = wk0[1][0]*n0;
         
  
      /* VERTEX 2         */
      wk0[2][0] = lin[n][2]*pgn[2];
      wk1[2][0] = lin[n][2]*dpgn[2];

      /* SIDE 1      */
      for (m = sm+3; m < 2*sm+3; ++m) {
         wk0[2][0] += lin[n][m]*pgn[m];
         wk1[2][0] += lin[n][m]*dpgn[m];
      }         
       wk2[2][0] = wk0[2][0]*n0;

      for(m = 3; m < sm+3; ++m) {
         /* SIDE 0      */
         wk0[m][0] = lin[n][m]*pgn[m];
         wk1[m][0] = lin[n][m]*dpgn[m];
         wk2[m][0] = wk0[m][0]*n0;
      }
#ifndef VERTEX
      }
#endif

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

      for(k=0; k < p+1; ++k )  
         f[n]   += lin[n][k]*pgx[k];

   }
   
   return;
}

void hpbasis::ptprobe1d(int nv, FLT **lin, FLT *f, FLT *dx) {
   int k,n;

   for(n=0;n<nv;++n) {
      f[n]   = 0.0;
      dx[n]  = 0.0;

      for(k=0; k < p+1; ++k )  {
         f[n]   += lin[n][k]*pgx[k];
         dx[n]  += lin[n][k]*dpgx[k];
      }
   }
   
   return;
}
   
