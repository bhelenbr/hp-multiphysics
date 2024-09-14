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

//#define MMS
#define SUTHERLAND

class tri_hp_cns : public tri_hp {
	public:
		/* THINGS SHARED BY ALL tri_hp_cns in same multigrid block */
		struct hp_cns_global {
			/* STABILIZATION */
			Array<FLT,3> tau;

			/* PHYSICAL CONSTANTS */
			FLT kcond, mu, gamma,R, prandtl;
#ifdef SUTHERLAND
            FLT s1, s2;
#endif
			/* SOURCE FUNCTION FOR MMS */
#ifdef MMS
			init_bdry_cndtn *src;
#endif
			vsi	res_temp;
			
			/* preconditioner could make 2d but keep general for now */
			Array<FLT,3> vpreconditioner, spreconditioner, tpreconditioner;
        };
        shared_ptr<hp_cns_global> hp_cns_gbl;

		FLT adis; // DISSIPATION CONSTANT

		hp_vrtx_bdry* getnewvrtxobject(int bnum, std::string name);
		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);
		init_bdry_cndtn* getnewibc(std::string name);

    public:
		tri_hp_cns* create() { return new tri_hp_cns(); }

		void init(input_map& inmap,shared_ptr<block_global> gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
	
		void update();
		void error_estimator();
		int setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();
		void calculate_tau(Array<double,1> pvu, FLT h, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep);
		void project_new_variables();
		void project_res_vertex();
		void project_res_side(int mode);
		void project_res_interior();
		void switch_variables(Array<double,1> pvu, Array<double,1> &a);
#ifdef SUTHERLAND
        void Sutherland(FLT RT);
#endif

};
#endif
