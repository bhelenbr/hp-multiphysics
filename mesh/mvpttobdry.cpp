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