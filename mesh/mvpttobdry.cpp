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

/* CYLINDER OR WAVE */
#define WAVE

#ifdef CYLINDER
#define R 0.5
extern FLT center;
FLT hgt(int type, FLT x, FLT y) {
   x -= center;
	return(x*x + y*y - R*R);
}
FLT dhgtdx(int type, FLT x, FLT y) {
   x -= center;
	return(2*x);
}
FLT dhgtdy(int type, FLT x, FLT y) {
	return(2*y);
}
#endif

#ifdef WAVE
extern FLT amp;

FLT hgt(int type, FLT x, FLT y) {
	return(y - amp*sin(2.*M_PI*(x-0.0)));
}
FLT dhgtdx(int type, FLT x, FLT y) {
	return(-amp*2.*M_PI*cos(2.*M_PI*(x-0.0)));
}
FLT dhgtdy(int type, FLT x, FLT y) {
	return(1.0);
}
#endif

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