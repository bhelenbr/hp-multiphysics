/*
 *  intgrt1d.cpp
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"hpbasis.h"

void hpbasis::intgrt1d(FLT *f, FLT *rslt) {
   int i,n;
   
   for(n=0;n<sm+2;++n)
      rslt[n] = 0.0;
      
#ifdef VERTEX
   for(i=0;i<gpx;++i)
      rslt[0] += gxwtx[1][i]*f[i];
   
   for(i=0;i<gpx;++i)
      rslt[1] += gxwtx[2][i]*f[i];
#else
   for(i=0;i<gpx;++i)
      rslt[0] += gxwtx[0][i]*f[i];
      
   if (!p) return;
   
   for(i=0;i<gpx;++i)
      rslt[1] += (gxwtx[2][i] -gxwtx[1][i])*f[i];
#endif

   for(n=3;n<nmodx;++n)
      for(i=0;i<gpx;++i)
         rslt[n-1] += gxwtx[n][i]*f[i];
   
   return;
}
   
void hpbasis::intgrtx1d(FLT *f, FLT *rslt) {
   int i,n;
   
#ifdef VERTEX
   for(i=0;i<gpx;++i)
      rslt[0] -= dgxwtx[1][i]*f[i];
   
   for(i=0;i<gpx;++i)
      rslt[1] -= dgxwtx[2][i]*f[i];
#else
   if (!p) return;
   
   for(i=0;i<gpx;++i)
      rslt[1] -= (dgxwtx[2][i] -dgxwtx[1][i])*f[i];
#endif
   
   for(n=3;n<nmodx;++n)
      for(i=0;i<gpx;++i)
         rslt[n-1] -= dgxwtx[n][i]*f[i];
   
   return;
}



