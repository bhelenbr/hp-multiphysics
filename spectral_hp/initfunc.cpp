/*
 *  initfunc.cpp
 *  planar++
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

/* THESE ARE FUNCTIONS WHICH MUST BE SET UP FOR INDIVIDUAL PROBLEMS */

#include "spectral_hp.h"
#include "math.h"
#include<assert.h>

/***************************/
/* INITIALIZATION FUNCTION */
/***************************/
FLT f1(int n, FLT x, FLT y) {
   switch(n) {
      case(0):
         return(1 +exp(-10.*x*x));
      case(1):
         return(1 +exp(-10.*(x*x)));
      case(2):
         return(1 +exp(-10.*(x*x)));
   }
   return(0.0);
}

/***************************/
/*	CURVED SIDE DEFINITIONS */
/***************************/
/* FOR A SPHERE */
FLT rad = 0.5;
FLT centerx = 0.0;
FLT centery = 0.0;

/*	FOR A SINE WAVE */
FLT amp = 0.1;

FLT hgt(int type, FLT x, FLT y) {
   if (type&(EULR_MASK +INFL_MASK)) {
      x -= centerx;
      y -= centery;
      return(x*x + y*y - rad*rad);
   }
   
   if (type&(FSRF_MASK +IFCE_MASK)) {
      return(y -amp*sin(2.*M_PI*x));
   }
   
   return(0.0);
      
}
FLT dhgtdx(int type, FLT x, FLT y) {
   if (type&(EULR_MASK +INFL_MASK)) {
      x -= centerx;
      y -= centery;
      return(2*x);
   }
   
   if (type&(FSRF_MASK +IFCE_MASK)) {
      return(-2.*amp*M_PI*cos(2.*M_PI*x));
   }
   
   return(0.0);
}
FLT dhgtdy(int type, FLT x, FLT y) {   
   if (type&(EULR_MASK +INFL_MASK)) {
      x -= centerx;
      y -= centery;
      return(2*y);
   }
   
   if (type&(FSRF_MASK +IFCE_MASK)) {
      return(1.0);
   }
   
   return(0.0);
}

/* TO USE IFCE/FSRF FROM A DIFFERENT MESH */
class spectral_hp *tgt;

void mvpttobdry(int typ, FLT& x, FLT &y) {
   int sind,iter;
   FLT mag, delt_dist;
   
   if (typ&(EULR_MASK +INFL_MASK)) {
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

   if (typ&(FSRF_MASK +IFCE_MASK)) {

      tgt->bdry_locate(typ,x,y,sind);
      
      return;
   }

}
