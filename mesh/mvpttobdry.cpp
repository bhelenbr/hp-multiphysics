/*
 *  mvpttobdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include<float.h>
#include<math.h>
#include<stdlib.h>

#define R 0.5

extern FLT center;
extern FLT amp;

FLT hgt(int type, FLT x, FLT y) {
   if (type&EULR_MASK) {
      x -= center;
      return(x*x + y*y - R*R);
   }
   if (type&(FSRF_MASK +IFCE_MASK)) {
      return(y - amp*cos(2.*M_PI*(x-0.0)));
   }
   
   return(0.0);
}
FLT dhgtdx(int type, FLT x, FLT y) {
   if (type&EULR_MASK) {
      x -= center;
      return(2*x);
   }
   if (type&(FSRF_MASK +IFCE_MASK)) {
      return(amp*2.*M_PI*sin(2.*M_PI*(x-0.0)));
   }
   return(1.0);
}
FLT dhgtdy(int type, FLT x, FLT y) {
   if (type&EULR_MASK) {
      return(2*y);
   }
   if (type&(FSRF_MASK +IFCE_MASK)) {
      return(1.0);
   }
   return(1.0);
}

void mvpttobdry(int typ, FLT& x, FLT &y) {
   int iter;
   FLT mag, delt_dist;

   iter = 0;
   do {
      mag = sqrt(dhgtdx(typ,x,y)*dhgtdx(typ,x,y) +dhgtdy(typ,x,y)*dhgtdy(typ,x,y));
      delt_dist = -hgt(typ,x,y)/mag;
      x += delt_dist*dhgtdx(typ,x,y)/mag;
      y += delt_dist*dhgtdy(typ,x,y)/mag;
      if (++iter > 100) {
         printf("iterations exceeded curved boundary %d %f %f\n",typ,x,y);
         exit(1);
      }
   } while (fabs(delt_dist) > 10.*EPSILON);
   
   return;
}

