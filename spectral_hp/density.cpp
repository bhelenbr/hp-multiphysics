/*
 *  density.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"spectral_hp.h"
#include<math.h>
#include<utilities.h>

/* THIS FUNCTION WILL SET THE vlngth VALUES BASED ON THE TRUNCATION ERROR */

void spectral_hp::density(FLT trncerr, FLT min, FLT max) {
   int i,m,n,v0,v1,indx;
   FLT sum,denom;
   
//   for(i=0;i<nvrtx;++i)
//      vlngth[i] = 1.05*3.1415/20.0*(1.0  - 0.875*exp(-((vrtx[i][0] - center)*(vrtx[i][0] -center) + vrtx[i][1]*vrtx[i][1]) +0.5*0.5));
//
//   for(i=0;i<nvrtx;++i)
//      vlngth[i] = 0.5;

   
   for(i=0;i<nvrtx;++i)
      fltwk[i] = 0.0;

   switch(b.p) {
      case(1):
         for(i=0;i<nside;++i) {
            v0 = svrtx[i][0];
            v1 = svrtx[i][1];
            sum = 0.0;
            for(n=0;n<NV;++n)
               sum += fabs(vug[v0][n] -vug[v1][n])/(fabs(vug[v0][n]) +fabs(vug[v1][n]) +100.*EPSILON);
            fltwk[v0] += sum;
            fltwk[v1] += sum;
         }
         break;
         
      default:
         indx = 0;
         for(i=0;i<nside;++i) {
            v0 = svrtx[i][0];
            v1 = svrtx[i][1];
            sum = 0.0;
            for(n=0;n<NV;++n) {
               denom = fabs(vug[v0][n]) +fabs(vug[v1][n]) +100.*EPSILON;
               for(m=0;m<b.sm;++m)
                  denom += fabs(sug[indx+m][n]);
               sum += fabs(sug[indx+b.sm -1][n])/denom;
               
            }
            fltwk[v0] += sum;
            fltwk[v1] += sum;
            indx += sm0;
         }
         break;
   }

#ifdef SKIP   
   for(i=0;i<nvrtx;++i) {
      fltwk[i] = pow(fltwk[i]/(nnbor[i]*trncerr),1./b.p);
      if (fltwk[i] < 0.5)
         vlngth[i] = MIN(2.0*vlngth[i],max);
      if (fltwk[i] > 2.0)
         vlngth[i] = MAX(0.5*vlngth[i],min);
   }
#endif
   
   for(i=0;i<nvrtx;++i) {
      fltwk[i] = pow(fltwk[i]/(nnbor[i]*trncerr),1./b.p);
      fltwk[i] = MAX(0.5,fltwk[i]);
      fltwk[i] = MIN(2.0,fltwk[i]);
      vlngth[i] /= fltwk[i];
      vlngth[i] = MIN(vlngth[i],max);
      vlngth[i] = MAX(vlngth[i],min);
   }
   
   for(i=0;i<nside;++i)
      fltwk[i] = (vlngth[svrtx[i][0]] +vlngth[svrtx[i][1]])/
      (2.*distance(svrtx[i][0],svrtx[i][1]));
   
   return;
}


