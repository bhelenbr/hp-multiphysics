/*
 *  bdry_ps.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ps.h"
#include "../hp_boundary.h"
#include <myblas.h>

namespace bdry_ps {

   class neumann : public hp_side_bdry {
      protected:
         tri_hp_ps &x;
         virtual void flux(TinyVector<FLT,3> u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> norm, TinyVector<FLT,3>& flx) {
            flx(2) = u(0)*norm(0) +u(1)*norm(1);
            flx(0) = 0.0;
            flx(1) = 0.0;
            
            return;
         }
      
      public:
         neumann(tri_hp_ps &xin, side_bdry &bin) : hp_side_bdry(xin,bin), x(xin) {mytype = "neumann";}
         neumann(const neumann& inbdry, tri_hp_ps &xin, side_bdry &bin) : hp_side_bdry(inbdry,xin,bin), x(xin) {}
         neumann* create(tri_hp& xin, side_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_ps&>(xin),bin);}
         block::ctrl rsdl(block::ctrl ctrl_message);
   };



   class dirichlet : public neumann {      
      void flux(TinyVector<FLT,3> u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> norm, TinyVector<FLT,3>& flx) {
         /* THESE DON'T GET USED */
         flx(0) = 0.0;
         flx(1) = 0.0;
         /* MASS FLUX */
         flx(2) = u(0)*norm(0) +u(1)*norm(1);
         return;
      }
      
      public:
         dirichlet(tri_hp_ps &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "dirichlet";}
         dirichlet(const dirichlet& inbdry, tri_hp_ps &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
         dirichlet* create(tri_hp& xin, side_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_ps&>(xin),bin);}
         void vdirichlet() {
            int sind,v0;
                     
            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               v0 = x.sd(sind).vrtx(0);
               x.gbl_ptr->res.v(v0,Range(0,x.ND-1)) = 0.0;
            }
            v0 = x.sd(sind).vrtx(1);
            x.gbl_ptr->res.v(v0,Range(0,x.ND-1)) = 0.0;
         }
         
         void sdirichlet(int mode) {
            int sind;

            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               x.gbl_ptr->res.s(sind,mode,Range(0,x.ND-1)) = 0.0;
            }
         }
            
         block::ctrl tadvance(bool coarse, block::ctrl ctrl_message);
   };
   
   class friction_wall : public hp_side_bdry {
      protected:
         tri_hp_ps &x;
         int dir;
         FLT muwall;
         
         virtual void flux(TinyVector<FLT,3> u, TinyVector<FLT,2> stress, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> norm, TinyVector<FLT,3>& flx) {
            flx(2) = u(0)*norm(0) +u(1)*norm(1);
            flx(dir) = 0.0;  // not used;
            if (u(1-dir) > 1.0e-4) 
               flx(1-dir) = -muwall*fabs(stress(dir));
            else if (u(1-dir) < -1.0e-4)
               flx(1-dir) =  muwall*fabs(stress(dir));
            else
               flx(1-dir) = 0.0;
               
            flx(1-dir) = -muwall*fabs(stress(dir));

            return;
         }
      
      public:
         friction_wall(tri_hp_ps &xin, side_bdry &bin) : hp_side_bdry(xin,bin), x(xin) {mytype = "friction_wall";}
         friction_wall(const friction_wall& inbdry, tri_hp_ps &xin, side_bdry &bin) : hp_side_bdry(inbdry,xin,bin), x(xin) {}
         friction_wall* create(tri_hp& xin, side_bdry &bin) const {return new friction_wall(*this,dynamic_cast<tri_hp_ps&>(xin),bin);}
         void init(input_map& inmap, void* &gbl_in) {
            std::string keyword;
                        
            hp_side_bdry::init(inmap,gbl_in);
            
            keyword = base.idprefix + "_friction";
            inmap.getwdefault(keyword,muwall,0.2);
            
            keyword = base.idprefix + "_dir";
            inmap.getwdefault(keyword,dir,0);
            
            return;
         }
         
         block::ctrl rsdl(block::ctrl ctrl_message);
         
         void vdirichlet() {
            int sind,v0;
                                 
            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               v0 = x.sd(sind).vrtx(0);
               x.gbl_ptr->res.v(v0,dir) = 0.0;
            }
            v0 = x.sd(sind).vrtx(1);
            x.gbl_ptr->res.v(v0,dir) = 0.0;
         }
         
         void sdirichlet(int mode) {
            int sind;

            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               x.gbl_ptr->res.s(sind,mode,dir) = 0.0;
            }
         }

   };
}
