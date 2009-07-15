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

#include "tri_hp_swirl.h"
#include "../ins/bdry_ins.h"
#include "../hp_boundary.h"

namespace bdry_swirl {

 class symmetry : public bdry_ins::neumann {        
		void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
			/* THESE DON'T GET USED */
			flx(Range(0,x.NV-1)) = 0.0;
			return;
		}

		public:
			symmetry(tri_hp_swirl &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "symmetry";}
			symmetry(const symmetry& inbdry, tri_hp_swirl &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			symmetry* create(tri_hp& xin, edge_bdry &bin) const {return new symmetry(*this,dynamic_cast<tri_hp_swirl&>(xin),bin);}
			void vdirichlet() {
				int sind,v0;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					x.gbl->res.v(v0,0) = 0.0;
					x.gbl->res.v(v0,2) = 0.0;
				}
				v0 = x.seg(sind).pnt(1);
				x.gbl->res.v(v0,0) = 0.0;
				x.gbl->res.v(v0,2) = 0.0;
			}

			void sdirichlet(int mode) {
				int sind;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					x.gbl->res.s(sind,mode,0) = 0.0;
				}
			}

			void tadvance();
	};        
}
