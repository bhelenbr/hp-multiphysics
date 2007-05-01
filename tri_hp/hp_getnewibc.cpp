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
         return(fcn(n).Eval(x,sim::time));
      }
      void input(input_map &inmap,std::string idnty) {
         std::string keyword,val;
         std::istringstream data;
         std::ostringstream nstr;
         int nvar;
         
         keyword = idnty +"_nvariable";

         if (!inmap.get(keyword,nvar))
            inmap.getwdefault("nvariable",nvar,1);
         
         fcn.resize(nvar);

         for(int n=0;n<nvar;++n) {
            nstr.str("");
            nstr << idnty << "_ibc" << n << std::flush;
            if (inmap.find(nstr.str() +"_expression") != inmap.end()) {
               fcn(n).init(inmap,nstr.str());
            }
            else {
               nstr.str("");
               nstr << "ibc" << n << std::flush;
               if (inmap.find(nstr.str() +"_expression") != inmap.end()) {
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
         
         keyword = idnty +"_nvariable";

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
            nstr << idnty << "_coeffx" << n << std::flush;
   
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
            nstr << idnty << "_coeffy" << n << std::flush;
            
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
         *sim::log << idnty << "_coeffx" << std::endl << cx;
         *sim::log << idnty << "_coeffy" << std::endl << cy;

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
         
         keyword = idnty +"_nvariable";

         if (!inmap.get(keyword,nvar))
            inmap.getwdefault("nvariable",nvar,1);
         
         c.resize(nvar);
         a.resize(nvar,2);
         spd.resize(nvar,2);

         for(int n=0;n<nvar;++n) {
            nstr.str("");
            nstr << idnty << "_coeff" << n << std::flush;
            if (!inmap.get(nstr.str(),c(n))) {
               nstr.str("");
               nstr << "coeff" << n << std::flush;
               inmap.getwdefault(nstr.str(),c(n),0.0);
            }

            nstr.str("");
            nstr << idnty << "_powers" << n << std::flush;
            if (!inmap.getline(nstr.str(),val)) {
               nstr.str("");
               nstr << "powers" << n << std::flush;
               inmap.getlinewdefault(nstr.str(),val,"0.0 0.0");
            }            
            data.str(val);
            data >> a(n,0) >> a(n,1);  
            data.clear(); 
            
            nstr.str("");
            nstr << idnty << "_speeds" << n << std::flush;
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
         
         keyword = idnty +"_nvariable";

         if (!inmap.get(keyword,nvar))
            inmap.getwdefault("nvariable",nvar,1);
         
         lam.resize(nvar,2);
         amp.resize(nvar,2);
         offset.resize(nvar,2);

         for(int n=0;n<nvar;++n) {
            nstr.str("");
            nstr << idnty << "_lambda" << n << std::flush;
            if (!inmap.getline(nstr.str(),val)) {
               nstr.str("");
               nstr << "lambda" << n << std::flush;
               inmap.getlinewdefault(nstr.str(),val,"1.0 1.0");
            }            
            data.str(val);
            data >> lam(n,0) >> lam(n,1);  
            data.clear(); 
            
            nstr.str("");
            nstr << idnty << "_amplitude" << n << std::flush;
            if (!inmap.getline(nstr.str(),val)) {
               nstr.str("");
               nstr << "amplitude" << n << std::flush;
               inmap.getlinewdefault(nstr.str(),val,"1.0 0.0");
            }            
            data.str(val);
            data >> amp(n,0) >> amp(n,1);  
            data.clear(); 
            
            nstr.str("");
            nstr << idnty << "_offset" << n << std::flush;
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
   keyword = std::string(idprefix) + "_ibc";
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

         
      }
      virtual void init(input_map& input, std::string idnty) {
         std::string keyword,val;
         std::istringstream data;
         std::ostringstream nstr;
         
         FLT vdflt[2] = {0.0, 0.0};
         if (!input.get(idnty +"_velocity",velocity.data(),2)) input.getwdefault("velocity",velocity.data(),2,vdflt);
      }
      
      block::ctrl tadvance(block::ctrl ctrl_message) {
         
         if (x.coarse) return(block::stop);
         
         if (ctrl_message == block::begin) {
            if (sim::substep == 0) x.l2error(x.gbl_ptr->ibc);
            
            TinyVector<FLT,2> dx;
            for (int n=0;n<mesh::ND;++n)
               dx(n) = sim::cdirk[sim::substep]/sim::dti*velocity(n);
            
            for(int i=0;i<x.nvrtx;++i)
               x.vrtx(i) += dx;
         }
         
         return(block::stop);
      }
};

class l2_error : public mesh_mover {
   public: 
      tri_hp &x;
      l2_error(tri_hp &xin) :mesh_mover(xin), x(xin) {}
      block::ctrl tadvance(block::ctrl ctrl_message) {
         
         if (x.coarse) return(block::stop);
         
         if (ctrl_message == block::begin && sim::substep == 0) {
            x.l2error(x.gbl_ptr->ibc);
         }
         return(block::stop);
      }
};

class print_averages : public mesh_mover {
   public: 
      tri_hp &x;
      print_averages(tri_hp &xin) :mesh_mover(xin), x(xin) {}
      block::ctrl tadvance(block::ctrl ctrl_message) {
         
         if (x.coarse) return(block::stop);
         
         if (ctrl_message == block::begin && sim::substep == 0) {
            Array<FLT,1> av(x.NV+3);
            x.integrated_averages(av);
            *sim::log << "# area: " << av(0) << " xbar: " << av(1) << " ybar: " << av(2);
            for(int n=0;n<x.NV;++n) 
               *sim::log << " V" << n << ": " << av(n+3);
            *sim::log << std::endl;
         }
         
         return(block::stop);
      }
};

class mesh_mover_type {
   public:
      const static int ntypes = 3;
      enum ids {translating,print_averages,l2error};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         int i;
         for(i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i);
         return(-1);
      }
};
const char mesh_mover_type::names[ntypes][40] = {"translating","print_averages","l2error"};


mesh_mover *tri_hp::getnewmesh_mover(input_map& inmap) {
   std::string keyword,movername;
   int type;
   
   /* FIND INITIAL CONDITION TYPE */
   keyword = std::string(idprefix) + "_mesh_mover";
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
      case mesh_mover_type::print_averages: {
         mesh_mover *temp = new print_averages(*this);
         return(temp);
      }
      case mesh_mover_type::l2error: {
         mesh_mover *temp = new l2_error(*this);
         return(temp);
      }
      default: {
         mesh_mover *temp = new mesh_mover(*this);
         return(temp);
      }
   }
}

