/*
 *  hp_initial_conditions.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include <mathclass.h>

class symbolic_ibc : public init_bdry_cndtn {
   private:
      Array<symbolic_function<2>,1> fcn;
   public:
      FLT f(int n, TinyVector<FLT,mesh::ND> x) {
         return(fcn(n).Eval(x));
      }
      void input(input_map &inmap,std::string idnty) {
         std::string keyword,val;
         std::istringstream data;
         std::ostringstream nstr;
         int nvar;
         
         keyword = idnty +".nvariable";

         if (!inmap.get(keyword,nvar))
            inmap.getwdefault("nvariable",nvar,1);
         
         fcn.resize(nvar);

         for(int n=0;n<nvar;++n) {
            nstr.str("");
            nstr << idnty << ".ibc" << n << std::flush;
            if (inmap.find(nstr.str() +".expression") != inmap.end()) {
               fcn(n).init(inmap,nstr.str());
            }
            else {
               nstr.str("");
               nstr << "ibc" << n << std::flush;
               if (inmap.find(nstr.str() +".expression") != inmap.end()) {
                  fcn(n).init(inmap,nstr.str());
               }
               else {
                  *sim::log << "couldn't find initial condition function\n";
                  exit(1);
               }
            }
         }
      }
};


class polynomial_ibc : public init_bdry_cndtn {
   private:
      int nvar;
      Array<FLT,2> cx;
      Array<FLT,2> cy;  
   
   public:
      FLT f(int n, TinyVector<FLT,mesh::ND> x) {
         FLT fx=0.0,fy=0.0;
         for(int i=0;i<5;++i) {
            fx += cx(n,i)*pow(x(0),i);
            fy += cy(n,i)*pow(x(1),i);
         }
         return(fx*fy);
      }
      
      void input(input_map &inmap,std::string idnty) {
         std::string keyword,val;
         std::istringstream data;
         std::ostringstream nstr;
         
         keyword = idnty +".nvariable";

         if (!inmap.get(keyword,nvar))
            inmap.getwdefault("nvariable",nvar,1);
         
         cx.resize(nvar,5);
         cy.resize(nvar,5);
         cx = 0.0;
         cy = 0.0;
         cx(Range::all(),0) = 1.0;
         cy(Range::all(),0) = 1.0;

         for(int n=0;n<nvar;++n) {
            nstr.str("");
            nstr << idnty << ".coeffx" << n << std::flush;
   
            if (inmap.getline(nstr.str(),val)) {
               data.str(val);
               for(int m=0;m<5;++m)
                  if (!(data >> cx(n,m))) break;
               data.clear(); 
            }
            else {
               nstr.str("");
               nstr << "coeffx" << n << std::flush;

               if (inmap.getline(nstr.str(),val)) {
                  data.str(val);
                  for(int m=0;m<5;++m)
                     if (!(data >> cx(n,m))) break;
                  data.clear(); 
               }
            }
            
            
            
            nstr.str("");
            nstr << idnty << ".coeffy" << n << std::flush;
            
            if (inmap.getline(nstr.str(),val)) {
               data.str(val);
               for(int m=0;m<5;++m)
                  if (!(data >> cy(n,m))) break;
               data.clear(); 
            }
            else {
               nstr.str("");
               nstr << "coeffy" << n << std::flush;
               if (inmap.getline(nstr.str(),val)) {
                  data.str(val);
                  for(int m=0;m<5;++m)
                     if (!(data >> cy(n,m))) break;
                  data.clear(); 
               }
            }
            
         }
         *sim::log << idnty << ".coeffx" << std::endl << cx;
         *sim::log << idnty << ".coeffy" << std::endl << cy;

      return;
      }
};

class power_ibc : public init_bdry_cndtn {
   private:
      Array<FLT,1> c;
      Array<FLT,2> a;
      Array<FLT,2> spd;
   public:
      FLT f(int n, TinyVector<FLT,mesh::ND> x) {
         return(c(n)*pow(x(0)-spd(n,0)*sim::time,a(n,0))*pow(x(1)-spd(n,1)*sim::time,a(n,1)));
      }
      
      void input(input_map &inmap,std::string idnty) {
         std::string keyword,val;
         std::istringstream data;
         std::ostringstream nstr;
         int nvar;
         
         keyword = idnty +".nvariable";

         if (!inmap.get(keyword,nvar))
            inmap.getwdefault("nvariable",nvar,1);
         
         c.resize(nvar);
         a.resize(nvar,2);
         spd.resize(nvar,2);

         for(int n=0;n<nvar;++n) {
            nstr.str("");
            nstr << idnty << ".coeff" << n << std::flush;
            if (!inmap.get(nstr.str(),c(n))) {
               nstr.str("");
               nstr << "coeff" << n << std::flush;
               inmap.getwdefault(nstr.str(),c(n),0.0);
            }

            nstr.str("");
            nstr << idnty << ".powers" << n << std::flush;
            if (!inmap.getline(nstr.str(),val)) {
               nstr.str("");
               nstr << "powers" << n << std::flush;
               inmap.getlinewdefault(nstr.str(),val,"0.0 0.0");
            }            
            data.str(val);
            data >> a(n,0) >> a(n,1);  
            data.clear(); 
            
            nstr.str("");
            nstr << idnty << ".speeds" << n << std::flush;
            if (!inmap.getline(nstr.str(),val)) {
               nstr.str("");
               nstr << "speeds" << n << std::flush;
               inmap.getlinewdefault(nstr.str(),val,"0.0 0.0");
            }            
            data.str(val);
            data >> spd(n,0) >> spd(n,1);  
            data.clear(); 
         }
      }
  
};

class sinusoidal_ibc : public init_bdry_cndtn {
   private:
      Array<FLT,2> lam,amp,offset;
   public:
      FLT f(int n, TinyVector<FLT,mesh::ND> x) {
         return(amp(n,0)*sin((x(0)-offset(n,0))*2*M_PI/lam(n,0)) +amp(n,1)*sin((x(1)-offset(n,1))*2*M_PI/lam(n,1)));
      }
      void input(input_map &inmap,std::string idnty) {
         std::string keyword,val;
         std::istringstream data;
         std::ostringstream nstr;
         int nvar;
         
         keyword = idnty +".nvariable";

         if (!inmap.get(keyword,nvar))
            inmap.getwdefault("nvariable",nvar,1);
         
         lam.resize(nvar,2);
         amp.resize(nvar,2);
         offset.resize(nvar,2);

         for(int n=0;n<nvar;++n) {
            nstr.str("");
            nstr << idnty << ".lambda" << n << std::flush;
            if (!inmap.getline(nstr.str(),val)) {
               nstr.str("");
               nstr << "lambda" << n << std::flush;
               inmap.getlinewdefault(nstr.str(),val,"1.0 1.0");
            }            
            data.str(val);
            data >> lam(n,0) >> lam(n,1);  
            data.clear(); 
            
            nstr.str("");
            nstr << idnty << ".amplitude" << n << std::flush;
            if (!inmap.getline(nstr.str(),val)) {
               nstr.str("");
               nstr << "amplitude" << n << std::flush;
               inmap.getlinewdefault(nstr.str(),val,"1.0 0.0");
            }            
            data.str(val);
            data >> amp(n,0) >> amp(n,1);  
            data.clear(); 
            
            nstr.str("");
            nstr << idnty << ".offset" << n << std::flush;
            if (!inmap.getline(nstr.str(),val)) {
               nstr.str("");
               nstr << "offset" << n << std::flush;
               inmap.getlinewdefault(nstr.str(),val,"0.0 0.0");
            }            
            data.str(val);
            data >> offset(n,0) >> offset(n,1);  
            data.clear(); 
         }
      }
};

class ibc_type {
   public:
      const static int ntypes = 4;
      enum ids {symbolic,polynomial,power,sinusoidal};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         int i;
         for(i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i);
         return(-1);
      }
};
const char ibc_type::names[ntypes][40] = {"symbolic","polynomial","power","sinusoidal"};



init_bdry_cndtn *tri_hp::getnewibc(input_map& inmap) {
   std::string keyword,ibcname;
   int type;

   /* FIND INITIAL CONDITION TYPE */
   keyword = std::string(idprefix) + ".ibc";
   if (!inmap.get(keyword,ibcname)) {
      if (!inmap.get("ibc",ibcname)) {
         *sim::log << "couldn't find initial condition type" << std::endl;
      }
   }
   
   type = ibc_type::getid(ibcname.c_str());
      
   switch(type) {
      case ibc_type::symbolic: {
         init_bdry_cndtn *temp = new symbolic_ibc;
         return(temp);
      }

      case ibc_type::polynomial: {
         init_bdry_cndtn *temp = new polynomial_ibc;
         return(temp);
      }
      case ibc_type::power: {
         init_bdry_cndtn *temp = new power_ibc;
         return(temp);
      }
      case ibc_type::sinusoidal: {
         init_bdry_cndtn *temp = new sinusoidal_ibc;
         return(temp);
      }

      default: {
         *sim::log << "couldn't find initial condition function " << ibcname << std::endl;
         exit(1);
      }
   }
}
      
      
class translating : public mesh_mover {
   public: 
      tri_hp &x;
      TinyVector<FLT,2> velocity;
      translating(tri_hp &xin) :mesh_mover(xin), x(xin) {}
   void move(int nvrtx, Array<TinyVector<FLT,mesh::ND>,1>& vrtx, Array<TinyVector<FLT,mesh::ND>,1>& vrtxbd) {
      int i,n;
      FLT dt;
      
      /* CALCULATE TIME INCREMENT */
      dt = sim::cdirk[sim::substep]/sim::dti;
      for(i=0;i<nvrtx;++i)
         for(n=0;n<mesh::ND;++n)
            //vrtx(i)(n) += dt*velocity(n);
            vrtx(i)(n) = vrtxbd(i)(n) +velocity(n)/(sim::bd[0]);
   
   }
   void input(input_map &inmap, std::string idnty) {
      std::string keyword,val;
      std::istringstream data;
      std::ostringstream nstr;
      
      keyword = idnty +".velocity";
      if (!inmap.getline(keyword,val)) inmap.getlinewdefault("velocity",val,"0.0 0.0");
      data.str(val);
      data >> velocity(0) >> velocity(1);  
      data.clear(); 
   }
};

class mesh_mover_type {
   public:
      const static int ntypes = 1;
      enum ids {translating};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         int i;
         for(i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i);
         return(-1);
      }
};
const char mesh_mover_type::names[ntypes][40] = {"translating"};


mesh_mover *tri_hp::getnewmesh_mover(input_map& inmap) {
   std::string keyword,movername;
   int type;
   
   /* FIND INITIAL CONDITION TYPE */
   keyword = std::string(idprefix) + ".mesh_mover";
   if (!inmap.get(keyword,movername)) {
      if (!inmap.get("mesh_mover",movername)) {
         type = -1;
      }
   }
   
   type = mesh_mover_type::getid(movername.c_str());
      
   switch(type) {
      case mesh_mover_type::translating: {
         mesh_mover *temp = new translating(*this);
         return(temp);
      }
      default: {
         mesh_mover *temp = new mesh_mover(*this);
         return(temp);
      }
   }
}

