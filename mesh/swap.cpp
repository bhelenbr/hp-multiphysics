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
   int j,t1,t2,tind,s1,snum1,snum2,v0,v1,vt1,vt2,dir;
   int tempsid[2],tempsgn[2];
      
   t1 = stri[sind][0];
   t2 = stri[sind][1];
   
   if (t1 < 0 || t2 < 0) return(0);

   for(snum1 = 0; snum1 < 3; ++snum1)
      if (tside[t1].side[snum1] == sind) break;
   assert(snum1 != 3);

   for(snum2 = 0; snum2 < 3; ++snum2)
      if (tside[t2].side[snum2] == sind) break;
   assert(snum2 != 3); 
   
   dir = (1 +tside[t1].sign[snum1])/2;
   v0 = svrtx[sind][1-dir];
   v1 = svrtx[sind][dir];
   
   vt1 = tvrtx[t1][snum1];
   vt2 = tvrtx[t2][snum2];
   
   if (MIN(minangle(v0,v1,vt1),minangle(v1,v0,vt2)) >
       MIN(minangle(vt2,vt1,v0),minangle(vt1,vt2,v1)) -tol) return(0);
       
   printf("swapping: %d %d %d\n",sind,t1,t2);
   sinfo[sind] = -2; /* MARK SIDE AS TOUCHED */   
   
/*	SWAP SIDE */
   svrtx[sind][0] = vt2;
   svrtx[sind][1] = vt1;
   stri[sind][0] = t1;
   stri[sind][1] = t2;
   
   tvrtx[t1][0] = vt2;
   tvrtx[t1][1] = vt1;
   tvrtx[t1][2] = v0;
   
   tvrtx[t2][0] = vt1;
   tvrtx[t2][1] = vt2;
   tvrtx[t2][2] = v1;
   
   tempsid[0] = tside[t1].side[(snum1+2)%3];
   tempsgn[0] = tside[t1].sign[(snum1+2)%3];
   
   tempsid[1] = tside[t1].side[(snum1+1)%3];
   tempsgn[1] = tside[t1].sign[(snum1+1)%3];

/*	FIX STRI/TTRI FOR TWO SIDES */   
   s1 = tempsid[1];
   dir = (1 -tempsgn[1])/2;
   stri[s1][dir] = t2;
   tind = stri[s1][1-dir];
   if (tind > -1) {
      for(j=0;j<3;++j)
         if (ttri[tind][j] == t1) break;
      assert(j != 3);
      ttri[tind][j] = t2;
   }
   
   s1 = tside[t2].side[(snum2 +1)%3];
   dir = (1 -tside[t2].sign[(snum2 +1)%3])/2;
   stri[s1][dir] = t1;
   tind = stri[s1][1-dir];
   if (tind > -1) {
      for(j=0;j<3;++j)
         if (ttri[tind][j] == t2) break;
      assert(j != 3);
      ttri[tind][j] = t1;
   }
   
   tside[t1].side[2] = sind;
   tside[t1].sign[2] = 1;
   
   tside[t1].side[1] = tside[t2].side[(snum2 +1)%3];
   tside[t1].sign[1] = tside[t2].sign[(snum2 +1)%3];
   
   tside[t1].side[0] = tempsid[0];
   tside[t1].sign[0] = tempsgn[0];
   
   tempsid[0] = tside[t2].side[(snum2 +2)%3];
   tempsgn[0] = tside[t2].sign[(snum2 +2)%3];
   
   tside[t2].side[2] = sind;
   tside[t2].sign[2] = -1;
   
   tside[t2].side[1] = tempsid[1];
   tside[t2].sign[1] = tempsgn[1];
   
   tside[t2].side[0] = tempsid[0];
   tside[t2].sign[0] = tempsgn[0];
   
   vtri[v0] = t1;
   vtri[v1] = t2;
   
   for(j=0;j<3;++j) {
      s1 = tside[t1].side[j];
      dir = (1 +tside[t1].sign[j])/2;
      ttri[t1][j] = stri[s1][dir];
   }

   for(j=0;j<3;++j) {
      s1 = tside[t2].side[j];
      dir = (1 +tside[t2].sign[j])/2;
      ttri[t2][j] = stri[s1][dir];
   }
   
   return(1);
}
   
   
void mesh::swap(int nswp, int *swp, FLT tol = 0.0) {
   int i,flag;
   
   do {
      flag = 0;
      for(i=0;i<nswp;++i)
         flag += swap(swp[i],tol);
   } while(flag > 0);

   return;
}

void mesh::swap(FLT tol = 0.0) {
   int nswap,i;
   
/*	PERFORM EDGE SWAPPING */
   do {
      nswap = 0;
      for(i=0;i<nside;++i)
         nswap += swap(i,tol);
   } while(nswap > 0);
   
   return;
}