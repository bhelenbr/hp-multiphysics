/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_cns_h_
#define _tri_hp_cns_h_

#include "../tri_hp.h"
#include <blocks.h>
#include <myblas.h>

class tri_hp_cns : public tri_hp {
	public:
		enum error_estimator_type {none,energy_norm,scale_independent};

		/* THINGS SHARED BY ALL tri_hp_cns in same multigrid block */
		struct global : public tri_hp::global {
			/* STABILIZATION */
			Array<FLT,3> tau;

			/* PHYSICAL CONSTANTS */
			FLT kcond, mu, gamma,R;
			Array<FLT,1> D;
			TinyVector<double,tri_mesh::ND> body;

			/* STORAGE FOR CALCULATION OF ENERGY AND AREA */
			TinyVector<FLT,2> eanda, eanda_recv;
			
			/* SOURCE FUNCTION FOR MMS */
			//init_bdry_cndtn *src;
			
			error_estimator_type error_estimator;
			
			vsi	res_temp;
			
			/* preconditioner could make 2d but keep general for now */
			Array<FLT,3> vpreconditioner,tpreconditioner; 


		} *gbl;

		FLT adis; // DISSIPATION CONSTANT

		hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
		hp_edge_bdry* getnewsideobject(int bnum, input_map &bdrydata);
		init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);
		tri_hp_helper* getnewhelper(input_map& inmap);

    public:
		void* create_global_structure() {return new global;}
		tri_hp_cns* create() { return new tri_hp_cns(); }

		void init(input_map& input, void *gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
	
		//void minvrt();
		void update();
		void length();
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();
		void pennsylvania_peanut_butter(Array<double,1> pvu, FLT h, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep);
		void project_new_variables();
		void switch_variables(Array<double,1> pvu, Array<double,1> &a);

};
#endif
