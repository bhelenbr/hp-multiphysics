/*
 *  edgeswap.cpp
 *  mblock
 *
 *  Created by helenbrk on Mon Sep 17 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "utilities.h"
#include<stdio.h>
#include<assert.h>
#include<float.h>

int mesh::swap(int sind, FLT tol = 0.0) {
   /* static */int j,t1,t2,s1p,s2p,snum1,snum2,tind,sd,v0,v1,vt1,vt2,dir;
         
   t1 = stri[sind][0];
   t2 = stri[sind][1];
   
   if (t1 < 0 || t2 < 0) return(0);

   for(snum1 = 0; snum1 < 3; ++snum1)
      if (tside[t1].side[snum1] == sind) break;
   assert(snum1 != 3);

   for(snum2 = 0; snum2 < 3; ++snum2)
      if (tside[t2].side[snum2] == sind) break;
   assert(snum2 != 3); 
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   vt1 = tvrtx[t1][snum1];
   vt2 = tvrtx[t2][snum2];
   
   if (MIN(minangle(v0,v1,vt1),minangle(v1,v0,vt2)) >
       MIN(minangle(vt2,vt1,v0),minangle(vt1,vt2,v1)) -tol) return(0);
       
   sinfo[sind] = -2; /* MARK SIDE AS TOUCHED */   
   
/*	SWAP SIDE */
   svrtx[sind][0] = vt2;
   svrtx[sind][1] = vt1;

/*	KEEP 2/3 VERTICES IN SAME SPOT */
   tvrtx[t1][(snum1 +2)%3] = vt2;
   tvrtx[t2][(snum2 +2)%3] = vt1;

   s1p = (snum1 +1)%3;
   s2p = (snum2 +1)%3;
   
/*	FIX 2 CHANGED EXTERIOR SIDES */
   tside[t1].side[snum1] = tside[t2].side[s2p];
   tside[t1].sign[snum1] = tside[t2].sign[s2p];
   ttri[t1][snum1] = ttri[t2][s2p];

   tside[t2].side[snum2] = tside[t1].side[s1p];
   tside[t2].sign[snum2] = tside[t1].sign[s1p];
   ttri[t2][snum2] = ttri[t1][s1p];
   
/*	FIX STRI/TTRI FOR 2 CHANGED EXTERIOR SIDES */
   sd = tside[t1].side[snum1];
   dir = (1 -tside[t1].sign[snum1])/2;
   stri[sd][dir] = t1;
   tind = ttri[t1][snum1];
   if (tind > -1) {
      for(j=0;j<3;++j)
         if (tside[tind].side[j] == sd) break;
      assert(j != 3);
      ttri[tind][j] = t1;
   }
   
   sd = tside[t2].side[snum2];
   dir = (1 -tside[t2].sign[snum2])/2;
   stri[sd][dir] = t2;
   tind = ttri[t2][snum2];
   if (tind > -1) {
      for(j=0;j<3;++j)
         if (tside[tind].side[j] == sd) break;
      assert(j != 3);
      ttri[tind][j] = t2;
   } 

   vtri[v0] = t1;
   vtri[v1] = t2;  
   
   ttri[t1][s1p] = t2;
   ttri[t2][s2p] = t1;
   
   tside[t1].side[s1p] = sind;
   tside[t2].side[s2p] = sind;  
   tside[t1].sign[s1p] =  1;
   tside[t2].sign[s2p] = -1;

   return(1);
}
   
   
void mesh::swap(int nswp, int *swp, FLT tol = 0.0) {
   /* static */int i,flag;
   
   do {
      flag = 0;
      for(i=0;i<nswp;++i)
         flag += swap(swp[i],tol);
   } while(flag > 0);

   return;
}

void mesh::swap(FLT tol = 0.0) {
   /* static */int nswap,i;
   
/*	PERFORM EDGE SWAPPING */
   do {
      nswap = 0;
      for(i=0;i<nside;++i)
         nswap += swap(i,tol);
      printf("#Swap cycle finished: %d sides swapped\n",nswap);
   } while(nswap > 0);
   
   return;
}

