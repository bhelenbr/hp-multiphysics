/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_ps_h_
#define _tri_hp_ps_h_

#include "../tri_hp.h"
#include <blocks.h>

class tri_hp_ps : public tri_hp {
	public:
		/* THINGS SHARED BY ALL tri_hp_ps in same multigrid block */
		struct hp_ps_global {
			/* STABILIZATION */
			Array<FLT,1> tau;

			/* PHYSICAL CONSTANTS */
			FLT mu, lami;
		};
        shared_ptr<hp_ps_global> hp_ps_gbl;

		FLT adis; // DISSIPATION CONSTANT

		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);

    public:
		tri_hp_ps* create() { return new tri_hp_ps(); }

		void init(input_map& inmap,shared_ptr<block_global> gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void length();
		int setup_preconditioner();
        void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();
};
#endif
