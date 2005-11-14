/*
 *  edgeswap.cpp
 *  mblock
 *
 *  Created by helenbrk on Mon Sep 17 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include <utilities.h>
#include<assert.h>
#include<float.h>

#ifdef DEBUG_ADAPT
extern int adapt_count;
static std::string adapt_file;
#endif

int mesh::swap(int sind, FLT tol) {
   int j,t1,t2,s1p,s2p,snum1,snum2,tind,sind1,v0,v1,vt1,vt2,dir;
         
   t1 = sd(sind).tri(0);
   t2 = sd(sind).tri(1);
   
   if (t1 < 0 || t2 < 0) return(0);

   for(snum1 = 0; snum1 < 3; ++snum1)
      if (td(t1).side(snum1) == sind) break;

   for(snum2 = 0; snum2 < 3; ++snum2)
      if (td(t2).side(snum2) == sind) break;
   
   v0 = sd(sind).vrtx(0);
   v1 = sd(sind).vrtx(1);
   vt1 = td(t1).vrtx(snum1);
   vt2 = td(t2).vrtx(snum2);
   
   if (MIN(minangle(v0,v1,vt1),minangle(v1,v0,vt2)) >
       MIN(minangle(vt2,vt1,v0),minangle(vt1,vt2,v1)) -tol -10.0*EPSILON) return(0);

   /* MARK TOUCHED */
   td(sind).info |= STOUC; 
   td(t1).info |= TTOUC;
   td(t2).info |= TTOUC;
   
   /* SWAP SIDE */
   sd(sind).vrtx(0) = vt2;
   sd(sind).vrtx(1) = vt1;

   /* KEEP 2/3 VERTICES IN SAME SPOT */
   td(t1).vrtx((snum1 +2)%3) = vt2;
   td(t2).vrtx((snum2 +2)%3) = vt1;

   s1p = (snum1 +1)%3;
   s2p = (snum2 +1)%3;
   
   /* FIX 2 CHANGED EXTERIOR SIDES */
   td(t1).side(snum1) = td(t2).side(s2p);
   td(t1).sign(snum1) = td(t2).sign(s2p);
   td(t1).tri(snum1) = td(t2).tri(s2p);

   td(t2).side(snum2) = td(t1).side(s1p);
   td(t2).sign(snum2) = td(t1).sign(s1p);
   td(t2).tri(snum2) = td(t1).tri(s1p);
   
   /* FIX STRI/TTRI FOR 2 CHANGED EXTERIOR SIDES */
   sind1 = td(t1).side(snum1);
   dir = (1 -td(t1).sign(snum1))/2;
   sd(sind1).tri(dir) = t1;
   tind = td(t1).tri(snum1);
   if (tind > -1) {
      for(j=0;j<3;++j)
         if (td(tind).side(j) == sind1) break;
      td(tind).tri(j) = t1;
   }
   
   sind1 = td(t2).side(snum2);
   dir = (1 -td(t2).sign(snum2))/2;
   sd(sind1).tri(dir) = t2;
   tind = td(t2).tri(snum2);
   if (tind > -1) {
      for(j=0;j<3;++j)
         if (td(tind).side(j) == sind1) break;
      td(tind).tri(j) = t2;
   } 

   vd(v0).tri = t1;
   vd(v1).tri = t2;  
   
   td(t1).tri(s1p) = t2;
   td(t2).tri(s2p) = t1;
   
   td(t1).side(s1p) = sind;
   td(t2).side(s2p) = sind;  
   td(t1).sign(s1p) =  1;
   td(t2).sign(s2p) = -1;
   
#ifdef DEBUG_ADAPT
      std::ostringstream nstr;
      nstr << adapt_count++ << std::flush;
      adapt_file = "adapt" +nstr.str();
      nstr.str("");
      output(adapt_file,grid);
#endif

   return(1);
}
   
   
void mesh::swap(int nswp, int *swp, FLT tol) {
   int i,flag;
   
   do {
      flag = 0;
      for(i=0;i<nswp;++i)
         flag += swap(swp[i],tol);
   } while(flag > 0);

   return;
}

void mesh::swap(FLT tol) {
   int nswap,i;
   
   /* PERFORM EDGE SWAPPING */
   do {
      nswap = 0;
      for(i=0;i<nside;++i) {
         nswap += swap(i,tol);
      }
      *sim::log << "#Swap cycle finished: " << nswap << " sides swapped" << std::endl;
   } while(nswap > 0);
   
   return;
}

