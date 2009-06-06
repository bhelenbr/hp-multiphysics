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

 void tri_hp_swirl::init(input_map& input, void *gin) { 

	gbl = static_cast<global *>(gin);    

    if (input.find(gbl->idprefix + "_nvariable") == input.end())
		input[gbl->idprefix + "_nvariable"] = "4";

    tri_hp_ins::init(input,gin);

#ifdef AXISYMMETRIC
	dpdz = 0.0;
#else
	if (!input.get(gbl->idprefix +"_dpdz",dpdz)) 
		input.getwdefault("dpdz",dpdz,0.0);
#endif

    return;
}
