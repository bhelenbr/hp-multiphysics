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
		enum error_estimator_type {none,energy_norm,scale_independent};
		
		/* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
		struct global : public tri_hp::global {
			/* STABILIZATION */
			Array<FLT,2> tau;

			/* PHYSICAL CONSTANTS */
			FLT rho, mu;
			Array<FLT,1> D;

			/* STORAGE FOR CALCULATION OF ENERGY AND AREA */
			TinyVector<FLT,3> eanda, eanda_recv;
			error_estimator_type error_estimator;

		} *gbl;

#ifdef DROP
		/** Rigid Mesh Motion for Finding Steady Translating Solutions */
		static TinyVector<FLT,ND> mesh_ref_vel;
#endif

		FLT adis; // DISSIPATION CONSTANT

		hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
		hp_edge_bdry* getnewsideobject(int bnum, input_map &bdrydata);
		init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);
		tri_hp_helper* getnewhelper(input_map& inmap);

	public:
		void* create_global_structure() {return new global;}
		tri_hp_ins* create() { return new tri_hp_ins(); }

		void init(input_map& input, void *gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

		void length();
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();

};
#endif
