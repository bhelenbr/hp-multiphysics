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

/* FOR A CIRCLE */
FLT rad = RADIUS;
FLT centerx = CENTERX;
FLT centery = CENTERY;

FLT amp = AMP;
FLT lam = LAM;

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
         return(1.0 +startup*amp*(x+0.5)*(x-3.0)*sin(M_PI*y));
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

#ifdef SPHERE
FLT f1(int n, FLT x, FLT y) {
   FLT r;
   
   r = sqrt(x*x +y*y);
   
   switch(n) {
      case(0):
         return(0.0);
      case(1):
         if (r < 0.55) 
            return(0.0);
         else
            return(1.0);
      case(2):
         return(0.0);
   }
   return(0.0);
}
#endif

#ifdef UNSTEADY_DROP
extern FLT outertime;

FLT f1(int n, FLT x, FLT y) {
   FLT r;
   
   r = sqrt(x*x +y*y);
   
   switch(n) {
      case(0):
         return(0.0);
      case(1):
         return(0.0+amp*sin(2.*M_PI*outertime/100.0));
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

#if (CASE == 0)
FLT forcing(FLT x,FLT y) { return(-cos(2.*M_PI*x));}
#else
FLT forcing(FLT x,FLT y) {return(0.0);}
#endif
FLT blayer = 0.0;
FLT axext = AMP*cos(M_PI*LAM/180.0), ayext = AMP*sin(M_PI*LAM/180.0);
FLT nuext = (1.0 - AMP);

FLT f1(int n, FLT x, FLT y) {
   FLT nux = nuext*4.*M_PI*M_PI;
   FLT axx = axext*2.*M_PI;
   FLT xx = x*2.*M_PI;
   FLT yx = y*2.*M_PI;
   double eps = 0.05;
   
   switch(n) {
      case((0+CASE)%3):
         return(axx*sin(xx)/(axx*axx +nux*nux)
           +nux*cos(xx)/(axx*axx +nux*nux)
           +(nux > 0.0 ? blayer*exp(axx/nux*xx) : 0.0)
           +startup*((sin(xx) +sin(100*xx))*(sin(yx)+sin(100*yx))));
      case((1+CASE)%3):
         return(startup*((sin(xx) +sin(100*xx))*(sin(yx)+sin(100*yx))));
      case((2+CASE)%3):
         if (x > 2.*eps)
            return(0.0);
         return(1+cos(M_PI*(x-eps)/eps));
   }
   return(0.0);
}

FLT df1d(int n, FLT x, FLT y) {
   FLT nux = nuext*4.*M_PI*M_PI;
   FLT axx = axext*2.*M_PI;
   FLT xx = x*2.*M_PI;
   switch(n) {
      case((0+CASE)%3):
         return(axx*cos(xx)/(axx*axx +nux*nux) 
           -nuext*sin(x)/(axx*axx +nux*nux)
           +(nux > 0.0 ? blayer*axx/nux*exp(axx/nux*xx) : 0.0));
      case(!CASE):
         return(0.0);
   }
   return(0.0);
}

#endif


#ifdef TWOLAYER

FLT body[2];
static FLT h = 2.0;
FLT mux[2];
FLT rhox[2];
FLT theta;


double f1(int n, double x, double y) { 
   FLT bf,re,g1,g2,n1,n2,q1,q2;
   
   /* FOR UIFACE TO BE 1 WITH D = 1, h = h/d */
   /* THETA DEFINED + CLOCKWISE */
   bf = mux[0]/(rhox[0]*(0.5 +(h-1)*rhox[1]/rhox[0])*sin(theta));
   body[0] = bf*sin(theta);
   body[1] = -bf*cos(theta);
   
   re = rhox[0]/mux[0];
   g1 = -bf*sin(theta);
   g2 = -bf*rhox[1]/rhox[0]*sin(theta);
   n1 = 1;
   q1 = 1;
   n2 = mux[1]/mux[0];
   q2 = rhox[1]/rhox[0];
   
   if (y < 1) {
      switch (n) {
         case(0):
            return(0.5*re*g1/n1*y*y +(re*g2*(1-h)-re*g1)*y);
         case(1):
            return(0.0);
         case(2):
            return(-bf*q1*cos(theta)*(y-h));
      }
   }
   else {
      switch (n) {
         case(0):
            return(0.5*re*g1/n2*y*y -re*g2*h/n2*y -0.5*g2*re/n2 +re*g2*h/n2 +re*g2*(1-h)-re*g1/2);
         case(1):
            return(0.0);
         case(2):
            return(-bf*q2*cos(theta)*(y-h));
      }
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
      return(y -amp*sin(2.*M_PI*x/lam));
   }
   return(0.0);
}

FLT dhgtdx(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      return(-2.*amp*M_PI/lam*cos(2.*M_PI*x/lam));
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

#ifdef NACA
FLT hgt(int type, FLT x, FLT y) { 
   double thickness = 0.12;
   if (type&(CURV_MASK)) {
      if (x < 0.0) return(x);
      
      return(thickness*(1.4845*sqrt(x) -0.63*x -1.758*pow(x,2) +1.4215*pow(x,3) -0.5180*pow(x,4)) - fabs(y));
      // return(thickness*(1.4845*sqrt(x) -0.63*x -1.758*pow(x,2) +1.4215*pow(x,3) -0.5075*pow(x,4)) - fabs(y));
   }
   
   return(0.0);
}

FLT dhgtdx(int type, FLT x, FLT y) {
   double thickness = 0.12;
   
   if (type&(CURV_MASK)) {
   
      if (x <= 0.0) return(1.0);
      return(thickness*(0.5*1.4845/sqrt(x) -0.63 -2*1.758*x +3*1.4215*pow(x,2) -4*0.5180*pow(x,3)));
   }
   return(0.0);
}

FLT dhgtdy(int type, FLT x, FLT y) {   
   if (type&(CURV_MASK)) {
      return((y > 0 ? -1 : 1));
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
            printf("#Warning: iterations exceeded curved boundary %d %f %f %f\n",typ,x,y,delt_dist);
            break;
            //exit(1);
         }
      } while (fabs(delt_dist) > 10.*EPSILON);
      
      return;
   }

#ifdef TWOLAYER
   if (typ == 1026) {
      if (startup) {
         iter = 0;
         do {
            mag = dhgtdx(typ,x,y)*dhgtdx(typ,x,y) +dhgtdy(typ,x,y)*dhgtdy(typ,x,y);
            delt_dist = -(hgt(typ,x,y)-1.0)/mag;
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

   if (typ == 1025) {
      if (startup) {
         return;
      }
      else
         tgt->findbdrypt(typ,x,y,psi);

      return;
   }
   
#endif

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
