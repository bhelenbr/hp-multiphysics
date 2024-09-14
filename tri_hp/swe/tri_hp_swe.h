/*
 *  hp_mgrid.h
 *  swe++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_swe_h_
#define _tri_hp_swe_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>

class tri_hp_swe : public tri_hp_ins {
	public:
		struct hp_swe_global {

			/* PHYSICAL CONSTANTS */
			FLT f0, cbeta, cd, ptest;

			/* BATHYMETRY DATA */
			init_bdry_cndtn *bathy;

		};
        shared_ptr<hp_swe_global> hp_swe_gbl;

		init_bdry_cndtn* getnewibc(std::string name);
		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);

	public:
		tri_hp_swe* create() { return new tri_hp_swe(); }

		void init(input_map& inmap,shared_ptr<block_global> gin);  
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

		int setup_preconditioner();
		void rsdl(int stage);
		void calculate_unsteady_sources() {tri_hp::calculate_unsteady_sources();}
};
#endif
