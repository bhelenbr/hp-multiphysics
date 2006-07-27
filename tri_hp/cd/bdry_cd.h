/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include "hp_boundary.h"
#include "myblas.h"

namespace bdry_cd {
   class dirichlet : public hp_side_bdry {
      tri_hp_cd &x;
      
      public:
         dirichlet(tri_hp_cd &xin, side_bdry &bin) : hp_side_bdry(xin,bin), x(xin) {mytype = "dirichlet";}
         dirichlet(const dirichlet& inbdry, tri_hp_cd &xin, side_bdry &bin) : hp_side_bdry(inbdry,xin,bin), x(xin) {}
         dirichlet* create(tri_hp& xin, side_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
         void vdirichlet() {
            int sind,v0;
                     
            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               v0 = x.sd(sind).vrtx(0);
               x.gbl_ptr->res.v(v0,0) = 0.0;
            }
            v0 = x.sd(sind).vrtx(1);
            x.gbl_ptr->res.v(v0,0) = 0.0;
         }
         
         void sdirichlet(int mode) {
            int sind;

            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               x.gbl_ptr->res.s(sind,mode,0) = 0.0;
            }
         }
            
         block::ctrl tadvance(bool coarse, block::ctrl ctrl_message); 
   };

   class neumann : public hp_side_bdry {
      protected:
         tri_hp_cd &x;
         virtual FLT flux(FLT u, TinyVector<FLT,mesh::ND> x, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm) {return(0.0);}
      
      public:
         neumann(tri_hp_cd &xin, side_bdry &bin) : hp_side_bdry(xin,bin), x(xin) {mytype = "neumann";}
         neumann(const neumann& inbdry, tri_hp_cd &xin, side_bdry &bin) : hp_side_bdry(inbdry,xin,bin), x(xin) {}
         neumann* create(tri_hp& xin, side_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
         block::ctrl rsdl(block::ctrl ctrl_message);
   };


   class characteristic : public neumann {
      public:
         FLT flux(FLT u, TinyVector<FLT,mesh::ND> pt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm) {
            FLT vel;

            vel =  (x.gbl_ptr->ax-mv(0))*norm(0) +(x.gbl_ptr->ay -mv(1))*norm(1);      


            if (vel > 0.0)
               return(vel*u);

            return(x.gbl_ptr->ibc->f(0, pt)*vel);
         }
         characteristic(tri_hp_cd &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
         characteristic(const characteristic &inbdry, tri_hp_cd &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
         characteristic* create(tri_hp& xin, side_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
   };

   class mixed : public neumann {
      public:
         TinyVector<FLT,5> c;
         
         FLT flux(FLT u, TinyVector<FLT,mesh::ND> pt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm) {
            
            FLT fout = 0.0;
            for(int i=0;i<5;++i)
               fout += c(i)*pow(u,i);
            
            return(fout);
         }
         mixed(tri_hp_cd &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "mixed";}
         mixed(const mixed& inbdry, tri_hp_cd &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
         mixed* create(tri_hp& xin, side_bdry &bin) const {return new mixed(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
         
         void init(input_map& inmap, void* &gbl_in) {
            std::string keyword;
            std::istringstream data;
            std::string val;
            
            neumann::init(inmap,gbl_in);

            keyword = base.idprefix + ".cd_mixed_coefficients";
            
            if (inmap.getline(keyword,val)) {
                  data.str(val);
                  data >> c(0) >> c(1) >> c(2) >> c(3) >> c(4);  
                  data.clear(); 
            }
            else {
               *sim::log << "couldn't find coefficients" << std::endl;
            }
            return;
         }
   };
}

