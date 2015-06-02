/*
 *  tri_hp_nonnewtonian.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/29/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */
 
#ifndef _tri_hp_nonnewtonian_h_
#define _tri_hp_nonnewtonian_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>
#include <symbolic_function.h>


class tri_hp_nonnewtonian : public tri_hp_ins {
	public:
		struct global : public tri_hp_ins::global {
			
			/* PHYSICAL CONSTANTS */
			vector_function mu_of_strain;
			
		} *gbl;
		
	public:
		void* create_global_structure() {return new global;}
		tri_hp_nonnewtonian* create() { return new tri_hp_nonnewtonian(); }
		
		void init(input_map& inmap, void *gin); 
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
};
#endif
