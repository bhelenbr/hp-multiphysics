/*
 *  getnewibc_ins.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ins.h"
#include "bdry_ins.h"

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
               case(2):
                  return(0.0);
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
            if (!blockdata.get(keyword,perturb_amp)) 
               blockdata.getwdefault("perturb_amplitude",perturb_amp,0.0); 

            alpha *= M_PI/180.0;
         }
   };
   
   class sphere : public init_bdry_cndtn {
      private:
         FLT speed,angle,inner,outer;
         TinyVector<FLT,mesh::ND> vel;
      
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            FLT r;
   
            r = sqrt(x(0)*x(0) +x(1)*x(1));
            switch(n) {
               case(0):case(1): 
                  if (r < inner) 
                     return(0.0);
                  else if (r < outer)
                     return(vel(n)*0.5*(1.-cos(M_PI*(r-inner)/(outer-inner))));
                  else
                     return(vel(n));
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;

            keyword = idnty +".flowspeed";
            if (!blockdata.get(keyword,speed)) 
               blockdata.getwdefault("flowspeed",speed,1.0);
    
            keyword = idnty +".angle";
            if (!blockdata.get(keyword,angle)) 
               blockdata.getwdefault("orientation",angle,0.0);
            angle *= M_PI/180.0;
            
            keyword = idnty +".inner_radius";
            if (!blockdata.get(keyword,inner)) 
               blockdata.getwdefault("inner_radius",inner,1.1);
               
            keyword = idnty +".outer_radius";
            if (!blockdata.get(keyword,outer)) 
               blockdata.getwdefault("outer_radius",outer,2.1);
            
            vel(0) = speed*cos(angle);
            vel(1) = speed*sin(angle);
         }
   };


   class accelerating : public init_bdry_cndtn {
      private:
         FLT speed,c,alpha;
         
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            switch(n) {
               case(0):
                  return(speed +c*pow(sim::time,alpha));
               case(1):
                  return(0.0);
               case(2):
                  return(-(x(0)-1)*c*alpha*pow(sim::time,alpha-1.0));
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;

            keyword = idnty +".speed";
            if (!blockdata.get(keyword,speed)) 
               blockdata.getwdefault("speed",speed,1.0);
    
            keyword = idnty +".coefficient";
            if (!blockdata.get(keyword,c)) 
               blockdata.getwdefault("coefficient",c,0.0);  
               
            keyword = idnty +".power";
            if (!blockdata.get(keyword,alpha)) 
               blockdata.getwdefault("power",alpha,0.0); 
         }
   };   
   

#ifdef TAYLOR
FLT ppipi = -0.5;

FLT f1(int n, FLT x, FLT y) {
   
   x *= 2.*M_PI;
   y *= 2.*M_PI;
   switch(n) {
      case(0):
         return(exp(-8.*M_PI*M_PI*sim::time)*sin(y)*cos(x));
      case(1):
         return(-exp(-8.*M_PI*M_PI*sim::time)*cos(y)*sin(x));
      case(2):
         return(0.5*exp(-16.*M_PI*M_PI*sim::time)*(sin(x)*sin(x)+sin(y)*sin(y)) +ppipi);
   }
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


   class stokes_drop_gas : public init_bdry_cndtn {
      private:
         FLT outer_limit;
         FLT mu_g, kappa;
      
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            FLT r,sint,cost;
            FLT ur,ut;   
            
            r = sqrt(x(0)*x(0) +x(1)*x(1));
            sint = x(0)/r;
            cost = -x(1)/r;
            ur = -(16.0*r*r*r+16.0*r*r*r*kappa-8.0*r*r-12.0*r*r*kappa+kappa)/(r*r*r)/(1.0+kappa)*cost/16.0;
            ut = sint*(32.0*r*r*r+32.0*r*r*r*kappa-8.0*r*r-12.0*r*r*kappa-kappa)/(r*r*r)/(1.0+kappa)/32.0;
            switch(n) {
               case(0):
                  if (r < outer_limit)
                     return(ur*sint+ut*cost);
                  else
                     return(0.0);
               case(1):
                  if (r < outer_limit)
                     return(-ur*cost+ut*sint);
                  else
                     return(1.0);
               case(2):
                  if (r < outer_limit)
                     return(mu_g/2*cost*(2+3*kappa)/(2*r*r*(1+kappa))); 
                  else
                     return(0.0);
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;

            keyword = idnty +".outer_radius";
            if (!blockdata.get(keyword,outer_limit))
               blockdata.getwdefault("outer_radius",outer_limit,75.0);
               
            keyword = idnty +".mu";
            if (!blockdata.get(keyword,mu_g)) {
               *sim::log << "couldn't find mu of gas" << std::endl;
               exit(1);
            }
 
            keyword = idnty +".liquid";
            if (!blockdata.get(keyword,val)) { 
               *sim::log << "couldn't find identity of liquid block" << std::endl;
               exit(1);
            }
            
            FLT mu_l;
            keyword = val +".mu";
            if (!blockdata.get(keyword,mu_l)) {
               *sim::log << "couldn't find mu of liquid" << std::endl;
               exit(1);
            }
            kappa = mu_l/mu_g;
         }
   };

   class stokes_drop_liquid : public init_bdry_cndtn {
      private:
         FLT outer_limit;
         FLT mu_l, kappa, sigma;
      
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            FLT r,sint,cost;
            FLT ur,ut;
            
            r = sqrt(x(0)*x(0) +x(1)*x(1));

            sint = x(0)/(r+FLT_EPSILON);
            cost = (x(1) > 0.0 ? -1 : 1)*sqrt(1.-sint*sint);
            ur = -(4.0*r*r-1.0)*cost/(1.0+kappa)/2.0;
            ut = sint*(8.0*r*r-1.0)/(1.0+kappa)/2.0;
            switch(n) {
               case(0):
                  return(ur*sint+ut*cost);
               case(1):
                  return(-ur*cost+ut*sint);
               case(2):
                  return(4*sigma -5*mu_l*r*cost*4/(1+kappa)); 
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
               
            keyword = idnty +".mu";
            if (!blockdata.get(keyword,mu_l)) {
               *sim::log << "couldn't find mu of liquid" << std::endl;
               exit(1);
            }
            keyword = idnty +".liquid_bdry";
            if (!blockdata.get(keyword,val)) { 
               *sim::log << "couldn't find identity of liquid boundary" << std::endl;
               exit(1);
            }
            
            keyword = val +".sigma";
            if (!blockdata.get(keyword,sigma)) {
               *sim::log << "couldn't find sigma" << std::endl;
               exit(1);
            }            
 
            keyword = idnty +".gas";
            if (!blockdata.get(keyword,val)) { 
               *sim::log << "couldn't find identity of gas block" << std::endl;
               exit(1);
            }
            
            FLT mu_g;
            keyword = val +".mu";
            if (!blockdata.get(keyword,mu_g)) {
               *sim::log << "couldn't find mu of gas" << std::endl;
               exit(1);
            }
            kappa = mu_l/mu_g;
         }
   };
   
   class translating_drop : public mesh_mover {
      private:
         tri_hp_ins &x;
         Array<FLT,1> avg;
         bdry_ins::surface *surf;
         FLT penalty;
         FLT delta_g_factor;
         int delta_interval;

      public:
         translating_drop(tri_hp_ins& xin) : mesh_mover(xin), x(xin) {
            int bnum;
            avg.resize(1+x.ND+x.NV);
            for(bnum=0;bnum<x.nsbd;++bnum) 
               if (surf = dynamic_cast<bdry_ins::surface *>(x.hp_sbdry(bnum))) break;
            assert(bnum < x.nsbd);
         }
         void init(input_map& input, std::string idnty) {
            std::string keyword;
            keyword = idnty + ".penalty_parameter";
            input.getwdefault(keyword,penalty,0.5);
            
            keyword = idnty + ".delta_g_factor";
            input.getwdefault(keyword,delta_g_factor,1.0);
            
            keyword = idnty + ".delta_interval";
            input.getwdefault(keyword,delta_interval,100);
         }
         mesh_mover* create(tri_hp& xin) { return new translating_drop(dynamic_cast<tri_hp_ins&>(xin)); }
         void calculate_stuff() {
            bdry_ins::surface::gbl *surf_gbl = surf->surf_gbl;
            
            /* DETRMINE CORRECTION TO CONSERVE AREA */
            /* IMPORTANT FOR STEADY SOLUTIONS */
            /* SINCE THERE ARE MULTIPLE STEADY-STATES */
            /* TO ENSURE GET CORRECT VOLUME */
            FLT rbar, kc; 
            kc = surf_gbl->sigma/(x.ins_gbl->mu +surf_gbl->mu2);
            x.integrated_averages(avg);
            rbar  = pow(3.*0.5*avg(0),1.0/3.0);
#ifdef DROP
            surf_gbl->vflux =  penalty*kc*(rbar -0.5);
            tri_hp_ins::mesh_ref_vel(1) = penalty*kc*avg(2) +avg(4);   
#endif

            /* C_D TO G CONVERSION REMINDER 
            re = 1.0/surf_gbl->mu2;
            cd = 24./re*(1 +0.1935*pow(re,0.6305));
            cd /= 16.0; // (1/2 rho u^2 * Pi r^2 / 2 pi);
            g = amp*(avg +avg) +12.*cd/(ins_gbl->rho -surf_gbl->rho2);
            */
            return;
         }

         
         block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE) {
            if (ctrl_message == block::begin /* && sim::dti == 0.0 */) calculate_stuff();
            return(block::stop);
         }
         
         block::ctrl setup_preconditioner(block::ctrl ctrl_message) {
            if (ctrl_message == block::begin /* && sim::dti == 0.0 */) calculate_stuff();
            return(block::stop);
         }
         
         block::ctrl tadvance(block::ctrl ctrl_message) {
            if (ctrl_message == block::begin && !x.coarse) {
               calculate_stuff();
#ifdef DROP
//               if (sim::dti > 0.0) surf->surf_gbl->vflux = 0.0;
               if ( (sim::tstep % delta_interval) == 0) {
                  sim::g *= delta_g_factor;
               }
#endif
            }
            return(block::stop);
         }
         
   };


   class mesh_mover_type {
      public:
         const static int ntypes = 1;
         enum ids {translating_drop};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
         }
   };
   const char mesh_mover_type::names[ntypes][40] = {"translating_drop"};

#ifdef UNSTEADY_DROP
FLT f1(int n, FLT x, FLT y) {
   FLT r;
   
   r = sqrt(x*x +y*y);
   
   switch(n) {
      case(0):
         return(0.0);
      case(1):
         // return(lam+amp*sin(2.*M_PI*sim::time/theta));
         if (sim::time/lam < 1.0)
            return((1.-cos(M_PI*sim::time/lam))*0.5*amp);
         else
            return(amp);
      case(2):
         if (sim::time/lam < 1.0)
            return(-y*M_PI/lam*sin(M_PI*sim::time/lam)*0.5*amp);
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
            return(30.75*12.0*2.54/100.0*amp*2.*M_PI/lam*sin(2.*M_PI/lam*sim::time)*x/r);
         break;
      case(1):
         if (r < 31.0*12.0*2.54/100.0) 
            return(30.75*12.0*2.54/100.0*amp*2.*M_PI/lam*sin(2.*M_PI/lam*sim::time)*y/r);
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
         rad = 30.75*12.0*2.54/100.0*amp*(1-cos(2.*M_PI*sim::time/lam));
         wallv = 30.75*12.0*2.54/100.0*amp*2.*M_PI/lam*sin(2.*M_PI*sim::time/lam);
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
         const static int ntypes = 5;
         enum ids {freestream,sphere,accelerating,stokes_drop_gas,stokes_drop_liquid};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char ibc_type::names[ntypes][40] = {"freestream","sphere","accelerating","stokes_drop_gas","stokes_drop_liquid"};

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
      case ibc_ins::ibc_type::accelerating: {
         init_bdry_cndtn *temp = new ibc_ins::accelerating;
         return(temp);
      }
      case ibc_ins::ibc_type::stokes_drop_gas: {
         init_bdry_cndtn *temp = new ibc_ins::stokes_drop_gas;
         return(temp);
      }
      case ibc_ins::ibc_type::stokes_drop_liquid: {
         init_bdry_cndtn *temp = new ibc_ins::stokes_drop_liquid;
         return(temp);
      }
      default: {
         return(tri_hp::getnewibc(inmap));
      }
   }
}

mesh_mover *tri_hp_ins::getnewmesh_mover(input_map& inmap) {
   std::string keyword,movername;
   int type;
   
   /* FIND INITIAL CONDITION TYPE */
   keyword = std::string(idprefix) + ".mesh_mover";
   if (!inmap.get(keyword,movername)) {
      if (!inmap.get("mesh_mover",movername)) {
         type = -1;
      }
   }
   
   type = ibc_ins::mesh_mover_type::getid(movername.c_str());
      
   switch(type) {
      case ibc_ins::mesh_mover_type::translating_drop: {
         mesh_mover *temp = new ibc_ins::translating_drop(*this);
         return(temp);
      }
      default: {
         return(tri_hp::getnewmesh_mover(inmap));
      }
   }
}

