/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cns_explicit.h"
#include "../hp_boundary.h"

void tri_hp_cns_explicit::init(input_map& input, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	gbl = static_cast<global *>(gin);

	if (input.find(gbl->idprefix + "_nvariable") == input.end()) {
		input[gbl->idprefix + "_nvariable"] = "4";
	}

	tri_hp::init(input,gin);

	input.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);


	gbl->tau.resize(maxpst,NV,NV);

	if (!input.get(gbl->idprefix + "_gamma",gbl->gamma)) input.getwdefault("gamma",gbl->gamma,1.403);
	if (!input.get(gbl->idprefix + "_mu",gbl->mu)) input.getwdefault("mu",gbl->mu,1.0);
	if (!input.get(gbl->idprefix + "_prandtl",gbl->kcond)) input.getwdefault("prandtl",gbl->kcond,0.75);
	//if (!input.get(gbl->idprefix + "_k",gbl->kcond)) input.getwdefault("k",gbl->kcond,1.0);
	if (!input.get(gbl->idprefix + "_R",gbl->R)) input.getwdefault("R",gbl->R,8.314472);

	gbl->kcond = gbl->R*gbl->mu/gbl->kcond*gbl->gamma/(gbl->gamma-1.);
	
	gbl->body(0) = 0.0;
	gbl->body(1) = 0.001;

	/* source term for MMS */
	//gbl->src = getnewibc("src",input);

	return;
}

void tri_hp_cns_explicit::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	const tri_hp_cns_explicit& inmesh = dynamic_cast<const tri_hp_cns_explicit &>(in);
	gbl = inmesh.gbl;

	tri_hp::init(in,why,sizereduce1d);

	adis = inmesh.adis;

	return;
}



