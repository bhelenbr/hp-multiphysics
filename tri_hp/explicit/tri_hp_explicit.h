/*
 *  tri_hp_explicit.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_explicit_h_
#define _tri_hp_explicit_h_

#include "../cd/tri_hp_cd.h"
#include <blocks.h>

class tri_hp_explicit : public tri_hp_cd {
	public:
		/* things shared by all tri_hp_explicit in same multigrid block */
		struct global : public tri_hp_cd::global {
			Array<FLT,3> sprcn2;  // Diagonal side preconditioner that can vary for each mode
			std::map<int,Array<FLT,2> > mass;  // mass matrix for curved elements
		} *gbl;

	public:
		void* create_global_structure() {return new global;}
		tri_hp_explicit* create() { return new tri_hp_explicit(); }
		void init(input_map& inmap, void *gin);
		void init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d);
		void minvrt();
		void setup_preconditioner();
};
#endif
