/*
 *  bdry.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 12/13/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */

#include "bdry_buoyancy.h"
#include <tri_boundary.h>

using namespace bdry_buoyancy;
#include "bdry_buoyancy.h"
#include "melt_buoyancy.h"
#include <myblas.h>

//#define MPDEBUG
//#define DEBUG

using namespace bdry_buoyancy;

void solid_fluid::init(input_map& inmap,void* gbl_in) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;
	
	/* Load in the additional flux for radiation if solid/gas boundary */
	keyword = base.idprefix + "_hp_typelist";
	if (inmap.find(base.idprefix+"_hp_typelist") == inmap.end()) {
		inmap[base.idprefix+"_hp_typelist"] = "0 0 1 1";
	}
	if (inmap.find(base.idprefix+"_flux2") == inmap.end()) {
		inmap[base.idprefix+"_flux2"] = "0.0)";
	}
	if (inmap.find(base.idprefix+"_flux3") == inmap.end()) {
		inmap[base.idprefix+"_flux3"] = "rho*(u0*n0 +u1*n1)";
	}
	
	inmap[base.idprefix+"_c0_indices"] = "2";
	hp_edge_bdry::init(inmap,gbl_in);
	
	return;
}

