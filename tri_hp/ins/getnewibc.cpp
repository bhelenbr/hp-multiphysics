/*
 *  getnewibc_ins.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ins.h"

namespace ibc_ins {

   class freestream : public init_bdry_cndtn {
      private:
         FLT alpha, speed,perturb_amp;
         
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            FLT amp = MAX(-sim::tstep,0)*perturb_amp;
            switch(n) {
               case(0):
                  return(speed*cos(alpha) +amp*x(0)*(1.0-x(0)));
               case(1):
                  return(speed*sin(alpha));
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;

            keyword = idnty +".flowspeed";
            if (!blockdata.get(keyword,speed)) 
               blockdata.getwdefault("flowspeed",speed,1.0);
    
            keyword = idnty +".flowangle";
            if (!blockdata.get(keyword,alpha)) 
               blockdata.getwdefault("flowangle",alpha,0.0);  
               
            keyword = idnty +".perturb_amplitude";
            if (!blockdata.get(keyword,alpha)) 
               blockdata.getwdefault("perturb_amplitude",perturb_amp,0.0); 

            alpha *= M_PI/180.0;
         }
   };
   
   class sphere : public init_bdry_cndtn {
      private:
         int orientation;
         FLT speed;
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            FLT r;
   
            r = sqrt(x(0)*x(0) +x(1)*x(1));

            switch(abs(orientation -n)) {
               case(0):
                  return(0.0);
               case(1):
                  if (r < 0.6) 
                     return(0.0);
                  else if (r < 1.1)
                     return(speed*0.5*(1.-cos(M_PI*(r-0.6)/0.5)));
                  else
                     return(speed);
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;

            keyword = idnty +".flowspeed";
            if (!blockdata.get(keyword,speed)) 
               blockdata.getwdefault("flowspeed",speed,1.0);
    
            keyword = idnty +".orientation";
            if (!blockdata.get(keyword,orientation)) 
               blockdata.getwdefault("orientation",orientation,0);   
         }
   };



   
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


   class ibc_type {
      public:
         const static int ntypes = 2;
         enum ids {freestream,sphere};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char ibc_type::names[ntypes][40] = {"freestream","sphere"};

}

 
init_bdry_cndtn *tri_hp_ins::getnewibc(input_map& inmap) {
   std::string keyword,ibcname;
   int type;

   /* FIND INITIAL CONDITION TYPE */
   keyword = std::string(idprefix) + ".ibc";
   if (!inmap.get(keyword,ibcname)) {
      if (!inmap.get("ibc",ibcname)) {
         *sim::log << "couldn't find initial condition type" << std::endl;
      }
   }
   type = ibc_ins::ibc_type::getid(ibcname.c_str());

      
   switch(type) {
      case ibc_ins::ibc_type::freestream: {
         init_bdry_cndtn *temp = new ibc_ins::freestream;
         return(temp);
      }
      case ibc_ins::ibc_type::sphere: {
         init_bdry_cndtn *temp = new ibc_ins::sphere;
         return(temp);
      }
      default: {
         return(tri_hp::getnewibc(inmap));
      }
   }
}

