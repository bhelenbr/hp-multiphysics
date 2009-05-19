/*
 *  intgrt1d.cpp
 *  tet_basis
 *
 *  Created by Michael Brazell on 5/8/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_basis.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef BZ_DEBUG
#define f(i) f1[i]
#define rslt(i) rslt1[i]
#endif

/* Integrates along Edge 1, and Vertices 2,3 */

void tet_basis::intgrt1d(FLT *rslt1, FLT *f1) {
   const int lgpx = gpx, lnmodx = nmodx;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,1> f(f1, shape(gpx), neverDeleteData);
   Array<FLT,1> rslt(rslt1, shape(2+em), neverDeleteData);
#endif
   
   lcl0 = 0.0;
   for(int i = 0; i < lgpx; ++i)
      lcl0 += gxwtx(1,i)*f(i);
   rslt(0) = lcl0;
   
   lcl0 = 0.0;
   for(int i = 0; i < lgpx; ++i)
      lcl0 += gxwtx(2,i)*f(i);
   rslt(1) = lcl0;

   for(int n = 3; n < lnmodx; ++n) {
      lcl0 = 0.0;
      for(int i = 0; i < lgpx; ++i)
         lcl0 += gxwtx(n,i)*f(i);
      rslt(n-1) = lcl0;
   }
   
   return;
}


void tet_basis::intgrtx1d(FLT *rslt1, FLT *f1) {
   const int lgpx = gpx, lnmodx = nmodx;
   FLT lcl0;
#ifdef BZ_DEBUG
   Array<FLT,1> f(f1, shape(gpx), neverDeleteData);
   Array<FLT,1> rslt(rslt1, shape(2+em), neverDeleteData);
#endif
   
   lcl0 = 0.0;
   for(int i = 0; i < lgpx; ++i)
      lcl0 -= dgxwtx(1,i)*f(i);
   rslt(0) = lcl0;
   
   lcl0 = 0.0;
   for(int i = 0; i < lgpx; ++i)
      lcl0 -= dgxwtx(2,i)*f(i);
   rslt(1) = lcl0;

   for(int n = 3; n < lnmodx; ++n) {
      lcl0 = 0.0;
      for(int i = 0; i < lgpx; ++i)
         lcl0 -= dgxwtx(n,i)*f(i);
      rslt(n-1) = lcl0;
   }

   return;
}