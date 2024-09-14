/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swirl.h"
#include "../hp_boundary.h"

 void tri_hp_swirl::init(input_map& inmap, shared_ptr<block_global> gin) { 

	gbl = gin;
	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end())
		inmap[gbl->idprefix + "_nvariable"] = "4";

	tri_hp_ins::init(inmap,gin);

#ifdef AXISYMMETRIC
	dpdz = 0.0;
#else
	if (!inmap.get(gbl->idprefix +"_dpdz",dpdz)) 
		inmap.getwdefault("dpdz",dpdz,0.0);
#endif

	return;
}
