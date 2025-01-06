/*
 *  tri_hp_cd.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tet_hp_cd_multi_h_
#define _tet_hp_cd_multi_h_

#include "../cd/tet_hp_cd.h"

class tet_hp_cd_multi : public tet_hp_cd {
	public:

		Array<int,1> marks;  // Integer identifying material for each element
	
		/* THINGS SHARED BY ALL tet_hp_cd_multi in same multigrid block */
		struct hp_cd_multi_global {
			int nmaterials;
			Array<FLT,1> kcond, rhocv;
		};
        shared_ptr<hp_cd_multi_global> hp_cd_multi_gbl;
		tet_hp_cd_multi* create() { return new tet_hp_cd_multi(); }
		void init(input_map& inmap, shared_ptr<block_global> gin);
		void init(const multigrid_interface& fine, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void calculate_unsteady_sources();
		int setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void connect(multigrid_interface& in);
};
#endif
