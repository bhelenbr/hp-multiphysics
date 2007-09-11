/*
 *  swe_bdry.h
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

#include "tri_hp_swe.h"
#include "../ins/bdry_ins.h"

namespace bdry_swe {
	
	 class wall : public bdry_ins::neumann {        
        void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm,  Array<FLT,1>& flx) {            
            flx(x.NV-1) = 0.0;
            
            /* X&Y MOMENTUM */
            for (int n=0;n<mesh::ND;++n)
                flx(n) = 0.5*sim::g*u(x.NV-1)*u(x.NV-1)*norm(n);
          
            /* EVERYTHING ELSE */
             for (int n=mesh::ND;n<x.NV-1;++n)
                flx(n) = 0.0;
            
            return;
        }
        
        public:
            wall(tri_hp_swe &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "wall";}
            wall(const wall& inbdry, tri_hp_swe &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
            wall* create(tri_hp& xin, side_bdry &bin) const {return new wall(*this,dynamic_cast<tri_hp_swe&>(xin),bin);}
    };    
        
    
    class characteristic : public bdry_ins::neumann {
        protected:
            tri_hp_swe &x;
            void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx);
        public:
            characteristic(tri_hp_swe &xin, side_bdry &bin) : neumann(xin,bin), x(xin) {mytype = "characteristic";}
            characteristic(const characteristic& inbdry, tri_hp_swe &xin, side_bdry &bin) : neumann(inbdry,xin,bin), x(xin) {}
            characteristic* create(tri_hp& xin, side_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_swe&>(xin),bin);}
    };
}
