/*
 *  proj1d.cpp
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include <hpbasis.h>

void hpbasis::proj1d(FLT *lin, FLT *f, FLT *dx) {
   static int i,n;
   
   for(i=0;i<gpx;++i) {
      f[i] = lin[0]*gx[0][i] +lin[1]*gx[1][i];
      dx[i] = lin[0]*dgx[0][i] +lin[1]*dgx[1][i];
   }
  
   for(n=2;n<sm+2;++n) {
      for(i=0;i<gpx;++i) {
         f[i] += lin[n]*gx[n+1][i];
         dx[i] += lin[n]*dgx[n+1][i];
      }
   }
         
   return;
}

void hpbasis::proj1d(FLT *lin, FLT *f) {
   static int i,n;
   
   for(i=0;i<gpx;++i)
      f[i] = lin[0]*gx[0][i] +lin[1]*gx[1][i];
  
   for(n=2;n<sm+2;++n)
      for(i=0;i<gpx;++i)
         f[i] += lin[n]*gx[n+1][i];
         
   return;
}

void hpbasis::proj1d(FLT u1, FLT u2, FLT *f) {
   static int i;
   
   for(i=0;i<gpx;++i)
      f[i] = u1*gx[0][i] +u2*gx[1][i];
      
   return;
}

void hpbasis::derivx1d(FLT *f, FLT *dx) {
	static int i,n;
	static FLT uddlt[MXTM];

   for (i=0;i<gpx;++i)
      uddlt[i] = f[i]*dltx[i];
      
   for (i=0;i<gpx;++i) {
      for(n=0;n<i;++n) 
         dx[i] += (uddlt[i] + uddlt[n])*dltx1[i][n];
      for(n=i+1;n<gpx;++n) 
         dx[i] += (uddlt[i] + uddlt[n])*dltx1[i][n];
   }

	return;
}

void hpbasis::proj1d_leg(FLT *lin, FLT *f) {
   static int i,m;

   for (i=1;i<sm+1;++i)
      f[i] = lin[0]*lgrnge1d[0][i]+lin[1]*lgrnge1d[1][i];
      
   for(m=2;m<sm+2;++m)
      for(i=1;i<sm+1;++i)
         f[i] += lin[m]*lgrnge1d[m][i];
   
   return;
}

void hpbasis::proj1d_leg(FLT u1, FLT u2, FLT *f) {
   static int i;
   
   for (i=1;i<sm+1;++i)
      f[i] = u1*lgrnge1d[0][i] +u2*lgrnge1d[1][i];
      
   return;
}