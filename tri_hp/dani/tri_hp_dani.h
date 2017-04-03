/*
 *  tri_hp_dani.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 3/3/17.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_dani_h_
#define _tri_hp_dani_h_

#include "../cd/tri_hp_cd.h"

class tri_hp_dani : public tri_hp_cd {
    public:
		tri_hp_dani* create() { return new tri_hp_dani();}
		void init(input_map& inmap, void *gin);
		void create_mass_matrix(sparse_row_major& mass);
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
};
#endif
