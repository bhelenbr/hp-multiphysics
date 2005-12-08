/*
 *  initfunc.cpp
 *  planar++
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

/* THESE ARE FUNCTIONS WHICH MUST BE SET UP FOR INDIVIDUAL PROBLEMS */

#include "math.h"
#include "tri_hp.h"

FLT amp,lam,theta;
FLT mux[4];
FLT rhox[4];
FLT sigmax[4];
FLT body[mesh::ND] = {0.0, -9.81};

/***************************/
/* INITIALIZATION FUNCTION */
/***************************/

/* FOR INITIALIZATION STARTUP WILL BE 1 AFTER THAT IT WILL ALWAYS BE 0 */
int startup = 1;

#ifdef TEST
FLT f1(int n, FLT x, FLT y) {
   switch(n) {
      case(0):
         return(pow(x,amp));
      case(1):
         return(pow(y,amp));
      case(2):
         return(pow(x*y,amp));
   }
   return(0.0);
}
#endif

#ifdef FREESTREAM
FLT f1(int n, FLT x, FLT y) {
   FLT xx = x*2.*M_PI;
   FLT yx = y*2.*M_PI;

   switch(n) {
      case(0):
         return(lam*cos(theta) +startup*amp*x*(1.-x)*((sin(xx) +sin(13*xx))*(sin(yx)+sin(5*yx))));
      case(1):
         return(lam*sin(theta) +startup*amp*x*(1.-x)*((sin(xx) +sin(5*xx))*(sin(yx)+sin(12*yx))));
      case(2):
         return(startup*amp*x*(1.-x)*((sin(xx) +sin(8*xx))*(sin(yx)+sin(8*yx))));
   }

   return(0.0);
}
#endif

#ifdef SOURCETERM
FLT f1(int n, FLT x, FLT y) {
   switch(n) {
      case(0):
         return(1.0);
      case(1):
         return(0.0);
      case(2):
         return(pow(x-0.5,amp) +startup*0.01*x*(1.-x));
   }

   return(0.0);
}
#endif

#ifdef ACCELERATING
FLT f1(int n, FLT x, FLT y) {
   switch(n) {
      case(0):
         return(lam +1.0*pow(outertime,amp));
      case(1):
         return(0.0);
      case(2):
         return(-(x-1)*amp*pow(outertime,amp-1.0));
   }

   return(0.0);
}
#endif



#ifdef SHEAR
FLT f1(int n, FLT x, FLT y) {
   FLT xx = x*2.*M_PI;
   FLT yx = y*2.*M_PI;

   switch(n) {
      case(0):
         return(1.0 +lam*y +startup*amp*x*(1.-x)*((sin(xx) +sin(13*xx))*(sin(yx)+sin(5*yx))));
      case(1):
         return(startup*amp*x*(1.-x)*((sin(xx) +sin(5*xx))*(sin(yx)+sin(12*yx))));
      case(2):
         return(startup*amp*x*(1.-x)*((sin(xx) +sin(8*xx))*(sin(yx)+sin(8*yx))));
   }

   return(0.0);
}
#endif

#ifdef TAYLOR
FLT ppipi = -0.5;

FLT f1(int n, FLT x, FLT y) {
   
   x *= 2.*M_PI;
   y *= 2.*M_PI;
   switch(n) {
      case(0):
         return(exp(-8.*M_PI*M_PI*outertime)*sin(y)*cos(x));
      case(1):
         return(-exp(-8.*M_PI*M_PI*outertime)*cos(y)*sin(x));
      case(2):
         return(0.5*exp(-16.*M_PI*M_PI*outertime)*(sin(x)*sin(x)+sin(y)*sin(y)) +ppipi);
   }
}
#endif

#ifdef CYLINDER
FLT f1(int n, FLT x, FLT y) {
   FLT r;
   
   r = sqrt(x*x +y*y);
   
   switch(n) {
      case(0):
         if (r < 0.51) 
            return(0.0);
         else
            return(1.0);
      case(1):
         if (r < 0.51) 
            return(0.0);
         else if (x < 10.0 && x > 0.0 && abs(y) < 10.0)
            return(amp*sin(2.*M_PI*x));
         else
            return(0.0);   
      case(2):
         return(0.0);
   }
   return(0.0);
}
#endif

#ifdef NACA
FLT f1(int n, FLT x, FLT y) {
   FLT r;
   
   r = sqrt(x*x +y*y);
   
   switch(n) {
      case(0):
         if (r < 5.0) 
            return(0.0);
         else if (r < 10.0 && r > 5.0)
            return((r-5.0)/5.0*lam*cos(theta));
         else
            return(lam*cos(theta));
      case(1):
         if (r < 5.0) 
            return(0.0);
         else if (r < 10.0 && r > 5.0)
            return((r-5.0)/5.0*lam*sin(theta));
         else
            return(lam*sin(theta));
      case(2):
         return(0.0);
   }
   return(0.0);
}
#endif



#ifdef NOZZLE
FLT f1(int n, FLT x, FLT y) { 
   switch(n) {
      case(0):
         return(0.0);
      case(1):
         if (y==0.0)
            return(-pow(1.0-x,1./7.));
         else
            return(0.0);
      case(2):
         return(0.0);
   }
   return(0.0);
}
#endif

#ifdef BOUNDARY_LAYER
FLT f1(int n, FLT x, FLT y) {   
   
   switch(n) {
      case(0):
         if (y < 1.0e-4) 
            return(1.0 -(1-y/1.0e-4)*(1-y/1.0e-4));
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

#ifdef DROP
FLT f1(int n, FLT x, FLT y) {
   FLT r,sint,cost,k;
   FLT ur,ut;
   
   k = mux[1]/mux[0];
   r = sqrt(x*x +y*y);

   sint = x/r;
   cost = -y/r;
   ur = -(16.0*r*r*r+16.0*r*r*r*k-8.0*r*r-12.0*r*r*k+k)/(r*r*r)/(1.0+k)*cost/16.0;
   ut = sint*(32.0*r*r*r+32.0*r*r*r*k-8.0*r*r-12.0*r*r*k-k)/(r*r*r)/(1.0+k)/32.0;
   switch(n) {
      case(0):
         if (r < 75.0)
            return(ur*sint+ut*cost);
         else
            return(0.0);
      case(1):
         if (r < 75.0)
            return(-ur*cost+ut*sint);
         else
            return(1.0);
      case(2):
         if (r < 75.0)
            return(mux[0]/2*cost*(2+3*k)/(2*r*r*(1+k))); 
         else
            return(0.0);
   }
   return(0.0);
}

FLT f2(int n, FLT x, FLT y) {
   FLT r,sint,cost,k;
   FLT ur,ut;
   
   k = mux[1]/mux[0];
   r = sqrt(x*x +y*y);

   sint = x/(r+FLT_EPSILON);
   cost = (y > 0.0 ? -1 : 1)*sqrt(1.-sint*sint);
   ur = -(4.0*r*r-1.0)*cost/(1.0+k)/2.0;
   ut = sint*(8.0*r*r-1.0)/(1.0+k)/2.0;
   switch(n) {
      case(0):
         return(ur*sint+ut*cost);
      case(1):
         return(-ur*cost+ut*sint);
      case(2):
         return(4*sigmax[1] -5*mux[1]*r*cost*4/(1+k)); 
   }
   return(0.0);
}
#endif

#ifdef UNSTEADY_DROP
FLT f1(int n, FLT x, FLT y) {
   FLT r;
   
   r = sqrt(x*x +y*y);
   
   switch(n) {
      case(0):
         return(0.0);
      case(1):
         // return(lam+amp*sin(2.*M_PI*outertime/theta));
         if (outertime/lam < 1.0)
            return((1.-cos(M_PI*outertime/lam))*0.5*amp);
         else
            return(amp);
      case(2):
         if (outertime/lam < 1.0)
            return(-y*M_PI/lam*sin(M_PI*outertime/lam)*0.5*amp);
         else
            return(0.0);
   }
   return(0.0);
}
#endif

#ifdef KOVASZNAY
double f1(int n, double x, double y) { 
   double re, lda;
   
   if (kovamu > 0.0) {
      re = 1/mux[0];
      lda = .5*re - sqrt(re*re*.25 + 4*M_PI*M_PI);
   }
   else
      lda = 0.0;

   switch (n) {
      case(0):
         return(1.0 - cos(2*M_PI*y)*exp(lda*x));
      case(1):
         return(lda/(2*M_PI)*sin(2*M_PI*y)*exp(lda*x));
      case(2):
         return(-.5*exp(2.*lda*x));
   }
   return(0.0);
}


double df1d(int n, double x, double y, int dir) {
   double re, lda;
   
   if (kovamu > 0.0) {
      re = 1/mux[0];
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

static FLT h = 2.0;

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

#ifdef ONELAYER
double f1(int n, double x, double y) { 
   FLT bf,re,n1,n2,n3,q1,q2,q3,h1,h2,h3;
   int mid,nonmid;
   
   /* FOR UIFACE TO BE 1 WITH D = 1, h = h/d */
   /* THETA DEFINED + CLOCKWISE */
   bf = 2.*mux[0]/(rhox[0]*sin(theta));
   body[0] = bf*sin(theta);
   body[1] = -bf*cos(theta);
   
   re = rhox[0]/mux[0];
   switch (n) {
      case(0):
         return(-(y+1.0)*(y-1.0) -1.0);
      case(1):
         return(0.0);
      case(2):
         return(-2.*cos(theta)/(sin(theta)*re)*y);
   }

   return(0.0);
}
#endif



#ifdef THREELAYER

double f1(int n, double x, double y) { 
   FLT bf,re,n1,n2,n3,q1,q2,q3,h1,h2,h3;
   int mid,nonmid;
   
   /* FOR UIFACE TO BE 1 WITH D = 1, h = h/d */
   /* THETA DEFINED + CLOCKWISE */
   /* FAIL PROOF TEST */
   if (fabs(mux[2] -mux[1]) < 1.0e-6) mid = 0;
   else if (fabs(mux[2] -mux[0]) < 1.0e-6) mid = 1;
   else mid = 2;
   nonmid = (mid+1)%3;
   
   bf = 2.*mux[nonmid]/(rhox[nonmid]*sin(theta));
   body[0] = bf*sin(theta);
   body[1] = -bf*cos(theta);
   
   re = rhox[nonmid]/mux[nonmid];
   
   h1 = 0.475;
   n1 = 1;
   q1 = 1;
   
   h2 = 0.525;
   n2 = mux[mid]/mux[nonmid];
   q2 = rhox[mid]/rhox[nonmid];
   
   h3 = 1.0;
   n3 = mux[nonmid]/mux[nonmid];
   q3 = rhox[nonmid]/rhox[nonmid];

   switch (n) {
      case(0):
         if (y <= 0.475) {
            double c1 = 2.0*h3/n1;
            double c2 = 0.0;
            return(-1./n1*y*y +c1*y +c2);
         }
         else if (y <= 0.525) {
            double c1 = 2.0*h3/n2;
            double c2 = h1*(h1*n1-h1*n2-2*h3*n1+2*h3*n2)/(n1*n2);
            return(-1./n2*y*y +c1*y +c2);
         }
         else {
            double c1 = 2*h3/n3;
            double c2 = (-h2*h2*n1*n3 +h2*h2*n1*n2 +h1*h1*n1*n3-h1*h1*n2*n3-2*h2*h3*n1*n2+2*h2*h3*n1*n3
                         -2*h1*h3*n1*n3 +2*h1*h3*n2*n3)/(n1*n2*n3);
            return(-1./n3*y*y +c1*y +c2);
         }
      case(1):
         return(0.0);
      case(2):
         return(-2.*cos(theta)/(sin(theta)*re)*y);
   }

   return(0.0);
}
#endif

#ifdef TWOSPHERE
FLT f1(int n, FLT x, FLT y) {
   double r = sqrt(x*x+y*y);
   
   switch(n) {
      case(0):
         if (r < 31.0*12.0*2.54/100.0) 
            return(30.75*12.0*2.54/100.0*amp*2.*M_PI/lam*sin(2.*M_PI/lam*outertime)*x/r);
         break;
      case(1):
         if (r < 31.0*12.0*2.54/100.0) 
            return(30.75*12.0*2.54/100.0*amp*2.*M_PI/lam*sin(2.*M_PI/lam*outertime)*y/r);
         break;
      case(2):
         return(rhox[0]*body[1]*(y-9.9690870000e+00));
         break;
   }

   return(0.0);
}
#endif

#ifdef COLUMN
FLT f1(int n, FLT x, FLT y) {
   FLT wallv,position,rad;
   
   switch(n) {
      case(0):
         rad = 30.75*12.0*2.54/100.0*amp*(1-cos(2.*M_PI*outertime/lam));
         wallv = 30.75*12.0*2.54/100.0*amp*2.*M_PI/lam*sin(2.*M_PI*outertime/lam);
         position = (x-rad)/(1.277112e+00-rad);
         return((1.0-position)*wallv);
      case(1):
         return(0.0);
         // return(1.0  +0.1*(y +1.064971e+01)*(1.277112e+00 -x)*x);  
      case(2):
         return(rhox[0]*body[1]*y);
   }

   return(0.0);
}
#endif

#ifdef IMPINGINGJET
 FLT f1(int n, FLT x, FLT y) {
    switch(n) {
       case(0):
          return(1);
       case(1):
          return(-y/x);
       case(2):
          return(0);
    }
    return(0.0);
}
#endif
   
#include "blocks.h"
#include "block.h"
#include "mgblock.h"
#include "r_mesh.h"
#include "tri_hp_cd.h"

FLT linear_src(FLT x, FLT y) {
   return(x);
}

static FLT apow;
FLT xpower(int n, TinyVector<FLT,2> x) { 
   return(pow(x(0),apow));
}

FLT sinewave(int n, TinyVector<FLT,2> x) { 
   return(sin(M_PI*x(0)));
}

FLT zero_src(FLT x, FLT y) {
   return(0.0);
}

class btype {
   public:
      enum ids {plain=1, cd, ins};
};

block* blocks::getnewblock(int idnum, input_map *blockdata) {
   std::string keyword,val;
   std::istringstream data;
   char idntystring[10];
   int type;        
   
   type = idnum&0xffff;

  if (blockdata) {
      sprintf(idntystring,"b%d",idnum);
      keyword = std::string(idntystring) + ".type";
      if ((*blockdata).get(keyword,val)) {
         data.str(val);
         data >> type;  
         data.clear(); 
      }
      else {
         if (!(*blockdata).get("blocktype",val)) {
            *sim::log << "couldn't find block type" << std::endl;
         }
      }
   }
      
   switch(type) {
      case btype::plain: {
         mgrid<r_mesh> *temp = new mgrid<r_mesh>(idnum);
         return(temp);
      }
      
      case btype::cd: {
         mgrid<tri_hp_cd> *temp = new mgrid<tri_hp_cd>(idnum);
         (*temp).gstorage.src = &zero_src;
         (*temp).gstorage.initfunc = &xpower;
         apow = 3.0;
         return(temp);
      }
      default: {
         std::cout << "unrecognizable block type: " <<  type << std::endl;
         mgrid<r_mesh> *temp = new mgrid<r_mesh>(idnum);
         return(temp);
      }
   } 
      
   return(0);
}

