/*
 *  lvlset_bdry.h
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
 
#ifndef _bdry_lvlset_h_
#define _bdry_lvlset_h_


#include "tri_hp_lvlset.h"
#include "../hp_boundary.h"
#include "../ins/bdry_ins.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>

namespace bdry_lvlset {

    class neumann : public bdry_ins::neumann {
        protected:
            tri_hp_lvlset &x;

            virtual void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx) {
                
                FLT rho = x.gbl_ptr->rho + (x.gbl_ptr->rho2-x.gbl_ptr->rho)*x.heavyside_if(u(2)/x.gbl_ptr->width);

                /* CONTINUITY */
                flx(x.NV-1) = rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
              
                /* X&Y MOMENTUM */
                for (int n=0;n<mesh::ND;++n)
                    flx(n) = flx(x.NV-1)*u(n) +x.gbl_ptr->ibc->f(x.NV-1, xpt)*norm(n);
              
                /* LEVEL-SET */
                flx(2) = 0.0;
        
                return;
            }
        
        public:
            neumann(tri_hp_lvlset &xin, side_bdry &bin) : bdry_ins::neumann(xin,bin), x(xin) {}
            neumann(const neumann& inbdry, tri_hp_lvlset &xin, side_bdry &bin) : bdry_ins::neumann(inbdry,xin,bin), x(xin) {}
            neumann* create(tri_hp& xin, side_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_lvlset&>(xin),bin);}
    };

    class inflow : public bdry_ins::inflow {  
        protected:
            tri_hp_lvlset &x;
    
        void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm,  Array<FLT,1>& flx) {
            
            /* CONTINUITY */
            FLT rho = x.gbl_ptr->rho + (x.gbl_ptr->rho2-x.gbl_ptr->rho)*x.heavyside_if(u(2)/x.gbl_ptr->width);
            flx(x.NV-1) = rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));

            /* EVERYTHING ELSE */
             for (int n=0;n<x.NV-1;++n)
                flx(n) = 0.0;
                
            return;
        }
        
        public:
            inflow(tri_hp_lvlset &xin, side_bdry &bin) : bdry_ins::inflow(xin,bin), x(xin) {}
            inflow(const inflow& inbdry, tri_hp_lvlset &xin, side_bdry &bin) : bdry_ins::inflow(inbdry,xin,bin), x(xin) {}
            inflow* create(tri_hp& xin, side_bdry &bin) const {return new inflow(*this,dynamic_cast<tri_hp_lvlset&>(xin),bin);}
    };
     
    class flow_inflow : public inflow {                
        public:
            flow_inflow(tri_hp_lvlset &xin, side_bdry &bin) : inflow(xin,bin) {mytype = "flow_inflow";}
            flow_inflow(const flow_inflow& inbdry, tri_hp_lvlset &xin, side_bdry &bin) : inflow(inbdry,xin,bin) {}
            flow_inflow* create(tri_hp& xin, side_bdry &bin) const {return new flow_inflow(*this,dynamic_cast<tri_hp_lvlset&>(xin),bin);}
            block::ctrl tadvance(bool coarse, block::ctrl ctrl_message);

            void vdirichlet() {
                int sind,v0;
                            
                for(int j=0;j<base.nel;++j) {
                    sind = base.el(j);
                    v0 = x.sd(sind).vrtx(0);
                    x.gbl_ptr->res.v(v0,Range(0,x.NV-3)) = 0.0;
                }
                v0 = x.sd(sind).vrtx(1);
                x.gbl_ptr->res.v(v0,Range(0,x.NV-3)) = 0.0;
            }
            
            void sdirichlet(int mode) {
                int sind;

                for(int j=0;j<base.nel;++j) {
                    sind = base.el(j);
                    x.gbl_ptr->res.s(sind,mode,Range(0,x.NV-3)) = 0.0;
                }
            }
    };

    
    class euler : public neumann {
        protected:
            void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm,  Array<FLT,1>& flx) {
                Array<FLT,1> ub(x.NV);
                
                FLT rho = x.gbl_ptr->rho + (x.gbl_ptr->rho2-x.gbl_ptr->rho)*x.heavyside_if(u(2)/x.gbl_ptr->width);

                for(int n=0;n<x.NV;++n)
                    ub(n) = x.gbl_ptr->ibc->f(n,xpt);
                
                flx(x.NV-1) = rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1));
                
                /* X&Y MOMENTUM */
                for (int n=0;n<mesh::ND;++n)
                    flx(n) = flx(x.NV-1)*ub(n) +u(x.NV-1)*norm(n);
              
                /* EVERYTHING ELSE */
                for (int n=mesh::ND;n<x.NV-1;++n)
                    flx(n) = flx(x.NV-1)*ub(n);
                    
                /* LEVEL-SET */
                flx(2) = 0.0;
                
                return;
            }
        public:
            euler(tri_hp_lvlset &xin, side_bdry &bin) : neumann(xin,bin) {mytype="euler";}
            euler(const euler& inbdry, tri_hp_lvlset &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
            euler* create(tri_hp& xin, side_bdry &bin) const {return new euler(*this,dynamic_cast<tri_hp_lvlset&>(xin),bin);}
    };
    


    
    class characteristic : public neumann {
        protected:
            void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx);
        public:
            characteristic(tri_hp_lvlset &xin, side_bdry &bin) : neumann(xin,bin) {mytype="characteristic";}
            characteristic(const characteristic& inbdry, tri_hp_lvlset &xin, side_bdry &bin) : neumann(inbdry,xin,bin)  {}
            characteristic* create(tri_hp& xin, side_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_lvlset&>(xin),bin);}
    };
    
//    class surface_levelset_hybrid : public bdry_ins::surface_slave {        
//        private:
//            int mp_phase,excpt,excpt1,stage;
//                     
//        public:
//            surface_levelset_hybrid(tri_hp_ins &xin, side_bdry &bin) : surface_slave(xin,bin) {mytype = "surface_levelset_hybrid";}
//            surface_levelset_hybrid(const surface_levelset_hybrid& inbdry, tri_hp_ins &xin, side_bdry &bin) : surface_slave(inbdry,xin,bin) {}
//            surface_levelset_hybrid* create(tri_hp& xin, side_bdry &bin) const {return new surface_levelset_hybrid(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//
//            void init(input_map& input,void* &gbl_in); 
//
//            /* FOR COUPLED DYNAMIC BOUNDARIES */
//            block::ctrl rsdl(block::ctrl ctrl_message);
//            block::ctrl update(block::ctrl ctrl_message);
//    };
    
}
#endif