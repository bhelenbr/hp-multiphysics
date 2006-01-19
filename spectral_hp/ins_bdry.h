/*
 *  ins_bdry.h
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

#include "tri_hp_ins.h"
#include "hp_boundary.h"
#include "myblas.h"

namespace ins_bdry {

   class neumann : public hp_side_bdry {
      protected:
         tri_hp_ins &x;
         virtual void flux(TinyVector<FLT,3> u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, TinyVector<FLT,3>& flx) {
            flx(2) = x.ins_gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
            flx(0) = flx(2)*u(0);
            flx(1) = flx(2)*u(1);
            return;
         }
      
      public:
         neumann(tri_hp_ins &xin, side_bdry &bin) : hp_side_bdry(xin,bin), x(xin) {mytype = "neumann";}
         neumann(const neumann& inbdry, tri_hp_ins &xin, side_bdry &bin) : hp_side_bdry(inbdry,xin,bin), x(xin) {}
         neumann* create(tri_hp& xin, side_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
         void addbflux();
   };



   class inflow : public neumann {      
      void flux(TinyVector<FLT,3> u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, TinyVector<FLT,3>& flx) {
         /* THESE DON'T GET USED */
         flx(0) = 0.0;
         flx(1) = 0.0;
         /* MASS FLUX */
         flx(2) = x.ins_gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
         return;
      }
      
      public:
         inflow(tri_hp_ins &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "inflow";}
         inflow(const inflow& inbdry, tri_hp_ins &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
         inflow* create(tri_hp& xin, side_bdry &bin) const {return new inflow(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
         void vdirichlet() {
            int sind,v0;
                     
            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               v0 = x.sd(sind).vrtx(0);
               x.hp_gbl->res.v(v0,Range(0,x.ND-1)) = 0.0;
            }
            v0 = x.sd(sind).vrtx(1);
            x.hp_gbl->res.v(v0,Range(0,x.ND-1)) = 0.0;
         }
         
         void sdirichlet(int mode) {
            int sind;

            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               x.hp_gbl->res.s(sind,mode,Range(0,x.ND-1)) = 0.0;
            }
         }
            
         block::ctrl tadvance(int excpt);
   };
   
   class euler : public neumann {
      protected:
         void flux(TinyVector<FLT,3> u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, TinyVector<FLT,3>& flx) {
            TinyVector<FLT,3> ub;
            for(int n=0;n<x.NV;++n)
               ub(n) = x.hp_gbl->ibc->f(n,xpt);
            
            flx(2) = x.ins_gbl->rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1));
            flx(0) = flx(2)*ub(0) +u(2)*norm(0);
            flx(1) = flx(2)*ub(1) +u(2)*norm(1);
            
            return;
         }
      public:
         euler(tri_hp_ins &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "euler";}
         euler(const euler& inbdry, tri_hp_ins &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
         euler* create(tri_hp& xin, side_bdry &bin) const {return new euler(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
   };
      
}