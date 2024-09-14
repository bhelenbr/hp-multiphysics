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

	class symmetry : public bdry_ins::generic {        
		void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
			/* THESE DON'T GET USED */
			flx(Range(0,x.NV-1)) = 0.0;
			return;
		}

		public:
			/* Fixme: Replace init with a constructor that accepts an input_map? */
			symmetry(tri_hp_swirl &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "symmetry";}
			symmetry(const symmetry& inbdry, tri_hp_swirl &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			symmetry* create(tri_hp& xin, edge_bdry &bin) const {return new symmetry(*this,dynamic_cast<tri_hp_swirl&>(xin),bin);}
			void init(input_map& inmap) {
				generic::init(inmap);
				type[0] = essential;
				essential_indices.push_back(0);
				type[2] = essential;
				essential_indices.push_back(2);
			}
			void tadvance();
	};        
}
