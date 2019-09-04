/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_cns_explicit_h_
#define _tri_hp_cns_explicit_h_

#include "../tri_hp.h"
#include <blocks.h>
#include <myblas.h>

class tri_hp_cns_explicit : public tri_hp {
	public:
		/* THINGS SHARED BY ALL tri_hp_cns in same multigrid block */
		struct global : public tri_hp::global {
			/* STABILIZATION */
			Array<FLT,3> tau;

			/* PHYSICAL CONSTANTS */
			FLT kcond, mu, gamma,R;
            
			/* STORAGE FOR CALCULATION OF ENERGY AND AREA */
			TinyVector<FLT,2> eanda, eanda_recv;
			
			/* SOURCE FUNCTION FOR MMS */
			//init_bdry_cndtn *src;
			
		} *gbl;

		FLT adis; // DISSIPATION CONSTANT

		hp_vrtx_bdry* getnewvrtxobject(int bnum, std::string name);
		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);
		init_bdry_cndtn* getnewibc(std::string name);

    public:
		void* create_global_structure() {return new global;}
		void delete_global_structure() {tri_hp::delete_global_structure(); delete gbl;}
		tri_hp_cns_explicit* create() { return new tri_hp_cns_explicit(); }

		void init(input_map& inmap, void *gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

		void error_estimator();
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void pennsylvania_peanut_butter(Array<double,1> u, FLT hmax, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep);


};
#endif
