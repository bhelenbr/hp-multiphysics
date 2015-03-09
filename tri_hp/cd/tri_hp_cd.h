/*
 *  tri_hp_cd.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_cd_h_
#define _tri_hp_cd_h_

#include "../tri_hp.h"
#include <blocks.h>

#define CONST_A

class tri_hp_cd : public tri_hp {
	public:
		/* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
		struct global : public tri_hp::global {
			/* STABILIZATION */
			Array<FLT,1> tau;

			/* PHYSICAL CONSTANTS */
			FLT rhocv,kcond,ax,ay;

			/* SOURCE FUNCTION & FLOW FIELDS */
			init_bdry_cndtn *src, *a;

			/* ADAPTATION LIMITS */
			FLT minlngth, maxlngth;
			
			//vsi	stiff_diag; // Stuff for Mike's minvrt

		} *gbl;

		FLT adis; // DISSIPATION CONSTANT

		hp_edge_bdry* getnewsideobject(int bnum, input_map &bdrydata);
		init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);

    public:
		void* create_global_structure() {return new global;}
		tri_hp_cd* create() { return new tri_hp_cd(); }
		void init(input_map& input, void *gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void calculate_unsteady_sources();
		void error_estimator();
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		//		void minvrt(); // Stuff for Mike's minvrt
};
#endif
