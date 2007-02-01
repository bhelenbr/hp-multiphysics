/*
 *  swirl_bdry.h
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
#include <myblas.h>

#include "tri_hp_swirl.h"
#include "../ins/bdry_ins.h"
#include "../hp_boundary.h"

namespace bdry_swirl {
	
	 class symmetry : public bdry_ins::neumann {      
      void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx) {
         /* THESE DON'T GET USED */
         flx(Range(0,x.NV-1)) = 0.0;
         return;
      }
      
      public:
         symmetry(tri_hp_swirl &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "symmetry";}
         symmetry(const symmetry& inbdry, tri_hp_swirl &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
         symmetry* create(tri_hp& xin, side_bdry &bin) const {return new symmetry(*this,dynamic_cast<tri_hp_swirl&>(xin),bin);}
         void vdirichlet() {
            int sind,v0;
                     
            for(int j=0;j<base.nel;++j) {
               sind = base.el(j);
               v0 = x.sd(sind).vrtx(0);
               x.gbl_ptr->res.v(v0,0) = 0.0;
					x.gbl_ptr->res.v(v0,2) = 0.0;
            }
            v0 = x.sd(sind).vrtx(1);
            x.gbl_ptr->res.v(v0,0) = 0.0;
				x.gbl_ptr->res.v(v0,2) = 0.0;
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
}
