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


// KOVASNAY TEST CYLINDER FREESTREAM
#define FREESTREAM

/***************************/
/* INITIALIZATION FUNCTION */
/***************************/
#ifdef TEST
FLT f1(int n, FLT x, FLT y) {
   switch(n) {
      case(0):
         return(1.0 +x*(1-x));
      case(1):
         return(x*(1-x)*y*(1-y));
      case(2):
         return(y*(1-y));
   }
   return(0.0);
}
#endif

#ifdef FREESTREAM
extern FLT outertime;

FLT f1(int n, FLT x, FLT y) {
   switch(n) {
      case(0):
         return(1.0 +0.0*outertime);
      case(1):
         return(0.0);
      case(2):
         return(0.0);
   }
   return(0.0);
}
#endif

#ifdef CYLINDER
FLT f1(int n, FLT x, FLT y) {
   FLT r;
   
   r = sqrt(x*x +y*y);
   
   switch(n) {
      case(0):
         if (r < 0.55) 
            return(0.0);
         else
            return(1.0);
      case(1):
         return(0.0);
      case(2):
         return(0.0);
   }
   return(0.0);
}
#endif

#ifdef KOVASZNAY
FLT kovamu = 0.025;

double f1(int n, double x, double y) { 
	double re, lda;

	if (kovamu > 0.0) {
		re = 1/kovamu;
		lda = .5*re - sqrt(re*re*.25 + 4*M_PI*M_PI);
	}
	else
		lda = 0.0;

   switch (n) {
      case(0):
         return(1 - cos(2*M_PI*y)*exp(lda*x));
      case(1):
      	return(lda/(2*M_PI)*sin(2*M_PI*y)*exp(lda*x));
      case(2):
         return(-.5*exp(2.*lda*x)+.5*exp(2.*lda*1.0));
   }
   return(0.0);
}


double df1d(int n, double x, double y, int dir) {
	double re, lda;

	if (kovamu > 0.0) {
		re = 1/kovamu;
		lda = .5*re - sqrt(re*re*.25 + 4*M_PI*M_PI);
	}
	else
		lda = 0.0;

   switch(n) {
      case(0):
         if (dir == 0)
            return( -lda*cos(2*M_PI*y)*exp(lda*x));
         else
            return( 2*M_PI*sin(2*M_PI*y)*exp(lda*x));
      case(1):
      	if (dir == 0) 
            return(lda/(2*M_PI)*lda*sin(2*M_PI*y)*exp(lda*x));
         else
            return(lda*cos(2*M_PI*y)*exp(lda*x));
   }
   return(0.0);
}
#endif

/***************************/
/*	CURVED SIDE DEFINITIONS */
/***************************/
/* FOR A SPHERE */
FLT rad = 0.5;
FLT centerx = 0.0;
FLT centery = 0.0;

/*	FOR A SINE WAVE */
FLT amp = 0.05;

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
int startup = 1;

void mvpttobdry(int typ, FLT& x, FLT &y) {
   int iter;
   FLT mag, delt_dist,psi;
   
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
      if (startup) {
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
      }
      else 
         tgt->findbdrypt(typ,x,y,psi);
      
      return;
   }

}
