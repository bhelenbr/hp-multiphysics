/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_ins_h_
#define _tri_hp_ins_h_

#include "../tri_hp.h"
#include <blocks.h>

class tri_hp_ins : public tri_hp {
	public:		
		/* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
		struct global : public tri_hp::global {
			/* STABILIZATION */
			Array<FLT,2> tau;

			/* PHYSICAL CONSTANTS */
			FLT rho, mu;
			Array<FLT,1> D;

		} *gbl;

		FLT adis; // DISSIPATION CONSTANT

		hp_vrtx_bdry* getnewvrtxobject(int bnum, std::string name);
		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);
		init_bdry_cndtn* getnewibc(std::string name);
		tri_hp_helper* getnewhelper(std::string name);

	public:
		void* create_global_structure() {return new global;}
		void delete_global_structure() {tri_hp::delete_global_structure(); delete gbl;}
		tri_hp_ins* create() { return new tri_hp_ins(); }

		void init(input_map& inmap, void *gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

		void error_estimator();
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();

};
#endif
