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
#include<stdlib.h>

/* INITIAL CONDITIONS */
/* KOVASZNAY TEST CYLINDER FREESTREAM */
#define CONVDIFF

/*	CURVED SURFACES */
/* CIRCLE SIN COS */
#define SIN

/* FOR A CIRCLE */
FLT rad = 0.5;
FLT centerx = 0.0;
FLT centery = 0.0;

/* FOR A SINE/COSINE WAVE */
#ifdef COS
FLT amp = 0.375;
#else
FLT amp = 0.25*0.075;
#endif


/***************************/
/* INITIALIZATION FUNCTION */
/***************************/

/* FOR INITIALIZATION STARTUP WILL BE 1 AFTER THAT IT WILL ALWAYS BE 0 */
int startup = 1;

#ifdef TEST
FLT f1(int n, FLT x, FLT y) {
   switch(n) {
      case(0):
         return(x);
      case(1):
         return(x*x);
      case(2):
         return(x*x*x);
   }
   return(0.0);
}
#endif

#ifdef FREESTREAM
extern FLT outertime;
FLT f1(int n, FLT x, FLT y) {
   switch(n) {
      case(0):
         return(1.0 +0.0*(x+0.5)*(x-3.0)*sin(M_PI*y));
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
         return(1 - cos(2*M_PI*y)*exp(lda*x) +0.01*(x+0.5)*(x-3.0));
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

#ifdef CONVDIFF

// FLT forcing(FLT x,FLT y) { return(-cos(2.*M_PI*x));}
FLT forcing(FLT x,FLT y) {return(0.0);}
FLT blayer = 0.0;
FLT axext = 1.0*cos(M_PI*10.0/180.0), ayext = 1.0*sin(M_PI*10.0/180.0);
FLT nuext = 0.0;

FLT f1(int n, FLT x, FLT y) {
   FLT nux = nuext*4.*M_PI*M_PI;
   FLT axx = axext*2.*M_PI;
   FLT xx = x*2.*M_PI;
   FLT yx = y*2.*M_PI;
   switch(n) {
      case(1):
         return(axx*sin(xx)/(axx*axx +nux*nux)
           +nux*cos(xx)/(axx*axx +nux*nux)
           +(nux > 0.0 ? blayer*exp(axx/nux*xx) : 0.0)
           +startup*(x*(1-x)*y*(1-y)*(sin(xx) +sin(100*xx))*(sin(yx)+sin(100*yx))));
      case(0):
         return(startup*(x*(1-x)*y*(1-y)*(sin(xx) +sin(100*xx))*(sin(yx)+sin(100*yx))));
   }
   return(0.0);
}

FLT df1d(int n, FLT x, FLT y) {
   FLT nux = nuext*4.*M_PI*M_PI;
   FLT axx = axext*2.*M_PI;
   FLT xx = x*2.*M_PI;
   switch(n) {
      case(1):
         return(axx*cos(xx)/(axx*axx +nux*nux) 
           -nuext*sin(x)/(axx*axx +nux*nux)
           +(nux > 0.0 ? blayer*axx/nux*exp(axx/nux*xx) : 0.0));
      case(0):
         return(0.0);
   }
   return(0.0);
}

#endif


/***************************/
/* CURVED SIDE DEFINITIONS */
/***************************/
#ifdef CIRCLE
FLT hgt(int type, FLT x, FLT y) {
   if (type&(CURV_MASK)) {
      x -= centerx;
      y -= centery;
      return(x*x + y*y - rad*rad);
   }
   return(0.0);
}

FLT dhgtdx(int type, FLT x, FLT y) {
   if (type&(CURV_MASK)) {
      x -= centerx;
      y -= centery;
      return(2*x);
   }
   return(0.0);
}

FLT dhgtdy(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      x -= centerx;
      y -= centery;
      return(2*y);
   }
   return(0.0);
}
#endif

#ifdef SIN
FLT hgt(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      return(y -amp*sin(2.*M_PI*x));
   }
   return(0.0);
}

FLT dhgtdx(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      return(-2.*amp*M_PI*cos(2.*M_PI*x));
   }
   return(0.0);
}

FLT dhgtdy(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      return(1.0);
   }
   return(0.0);
}
#endif

#ifdef COS
FLT hgt(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      return(y -amp*cos(2.*M_PI*x));
   }
   return(0.0);
}

FLT dhgtdx(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      return(2.*amp*M_PI*sin(2.*M_PI*x));
   }
   return(0.0);
}

FLT dhgtdy(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      return(1.0);
   }
   return(0.0);
}
#endif

/* TO USE IFCE/FSRF FROM A DIFFERENT MESH */
class spectral_hp *tgt;

void mvpttobdry(int typ, FLT& x, FLT &y) {
   int iter;
   FLT mag, delt_dist,psi;
   
   if (typ&(EULR_MASK +INFL_MASK)) {
      iter = 0;
      do {
         mag = dhgtdx(typ,x,y)*dhgtdx(typ,x,y) +dhgtdy(typ,x,y)*dhgtdy(typ,x,y);
         delt_dist = -hgt(typ,x,y)/mag;
         x += delt_dist*dhgtdx(typ,x,y);
         y += delt_dist*dhgtdy(typ,x,y);
         if (++iter > 100) {
            printf("#Warning: iterations exceeded curved boundary %d %f %f\n",typ,x,y);
            exit(1);
         }
      } while (fabs(delt_dist) > 10.*EPSILON);
      
      return;
   }
   
   if (typ == 1026) {
      if (startup) {
         iter = 0;
         do {
            mag = dhgtdx(typ,x,y)*dhgtdx(typ,x,y) +dhgtdy(typ,x,y)*dhgtdy(typ,x,y);
            delt_dist = -(hgt(typ,x,y)+0.5)/mag;
            x += delt_dist*dhgtdx(typ,x,y);
            y += delt_dist*dhgtdy(typ,x,y);
            if (++iter > 100) {
               printf("#Warning: iterations exceeded curved boundary %d %f %f\n",typ,x,y);
               exit(1);
            }
         } while (fabs(delt_dist) > 10.*EPSILON);
      }
      else 
         tgt->findbdrypt(typ,x,y,psi);
      
      return;
   }
   
   if (typ == 66562) {
       if (startup) {
         iter = 0;
         do {
            mag = dhgtdx(typ,x,y)*dhgtdx(typ,x,y) +dhgtdy(typ,x,y)*dhgtdy(typ,x,y);
            delt_dist = -(hgt(typ,x,y)-0.5)/mag;
            x += delt_dist*dhgtdx(typ,x,y);
            y += delt_dist*dhgtdy(typ,x,y);
            if (++iter > 100) {
               printf("#Warning: iterations exceeded curved boundary %d %f %f\n",typ,x,y);
               exit(1);
            }
         } while (fabs(delt_dist) > 10.*EPSILON);
      }
      else 
         tgt->findbdrypt(typ,x,y,psi);
      
      return;
   }  

   if (typ&(FSRF_MASK +IFCE_MASK)) {
      if (startup) {
         iter = 0;
         do {
            mag = dhgtdx(typ,x,y)*dhgtdx(typ,x,y) +dhgtdy(typ,x,y)*dhgtdy(typ,x,y);
            delt_dist = -hgt(typ,x,y)/mag;
            x += delt_dist*dhgtdx(typ,x,y);
            y += delt_dist*dhgtdy(typ,x,y);
            if (++iter > 100) {
               printf("#Warning: iterations exceeded curved boundary %d %f %f\n",typ,x,y);
               exit(1);
            }
         } while (fabs(delt_dist) > 10.*EPSILON);
      }
      else 
         tgt->findbdrypt(typ,x,y,psi);
      
      return;
   }

}
