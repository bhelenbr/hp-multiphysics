/*
 *  hp_mgrid.h
 *  heat++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_buoyancy_h_
#define _tri_hp_buoyancy_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>
#include <symbolic_function.h>

class tri_hp_buoyancy : public tri_hp_ins {
	public:
		struct hp_buoyancy_global {
			/* PHYSICAL CONSTANTS */
			FLT adapt_energy_scaling;
			FLT kcond,cp;
			symbolic_function<1> rho_vs_T;
		};
        shared_ptr<hp_buoyancy_global> hp_buoyancy_gbl;
		hp_vrtx_bdry* getnewvrtxobject(int bnum, std::string name);
		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);

	public:
		tri_hp_buoyancy* create() { return new tri_hp_buoyancy(); }

		void init(input_map& inmap,shared_ptr<block_global> gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		
		void error_estimator();
		int setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();
};
#endif
