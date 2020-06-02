/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_explicit.h"

void tri_hp_explicit::init(input_map& inmap, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	gbl = static_cast<global *>(gin);

	tri_hp_cd::init(inmap,gin);
    if (!inmap.get(gbl->idprefix + "_sigma",gbl->sigma)) {
        inmap.getwdefault("sigma",gbl->sigma,1.0);
    }
    
	gbl->sprcn2.resize(maxpst,sm0,NV);


	int tm0 = basis::tri(log2pmax)->tm();
	for (int tind=0; tind < ntri;++tind) {
		// tri(tind).info = 0; // FOR TESTING CURVED ALGORITHM
		if (tri(tind).info > -1)
			gbl->mass[tind].resize(tm0,tm0);
	}

	return;
}

void tri_hp_explicit::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	const tri_hp_explicit& inmesh = dynamic_cast<const tri_hp_explicit &>(in);
	gbl = inmesh.gbl;
	tri_hp_cd::init(in,why,sizereduce1d);
	
	return;
}
