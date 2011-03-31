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
		struct global : public tri_hp_ins::global {

			/* PHYSICAL CONSTANTS */
			FLT kcond,cp;
			symbolic_function<1> rho_vs_T;

		} *gbl;
		hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
		hp_edge_bdry* getnewsideobject(int bnum, input_map &bdrydata);

	public:
		void* create_global_structure() {return new global;}
		tri_hp_buoyancy* create() { return new tri_hp_buoyancy(); }

		void init(input_map& input, void *gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		
		void length();
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();
};
#endif
