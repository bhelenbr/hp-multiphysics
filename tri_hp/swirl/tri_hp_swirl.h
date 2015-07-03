/*
 *  hp_mgrid.h
 *  swirl++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_swirl_h_
#define _tri_hp_swirl_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>

class tri_hp_swirl : public tri_hp_ins {
	public:
		FLT dpdz; //< 3D driving pressure gradient
		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);
		init_bdry_cndtn* getnewibc(std::string name);

	public:
		void init(input_map& inmap, void *gin);  
		tri_hp_swirl* create() { return new tri_hp_swirl(); }
		void length();
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
};
#endif
