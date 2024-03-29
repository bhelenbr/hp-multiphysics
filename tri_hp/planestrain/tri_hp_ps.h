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
		struct global : public tri_hp::global {
			/* STABILIZATION */
			Array<FLT,1> tau;

			/* PHYSICAL CONSTANTS */
			FLT mu, lami;

			/* STORAGE FOR CALCULATION OF ENERGY AND AREA */
			TinyVector<FLT,2> eanda, eanda_recv;

		} *gbl;

		FLT adis; // DISSIPATION CONSTANT

		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);

    public:
		void* create_global_structure() {return new global;}
		void delete_global_structure() {tri_hp::delete_global_structure(); delete gbl;}
		tri_hp_ps* create() { return new tri_hp_ps(); }

		void init(input_map& inmap, void *gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void length();
		int setup_preconditioner();
        void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();
};
#endif
