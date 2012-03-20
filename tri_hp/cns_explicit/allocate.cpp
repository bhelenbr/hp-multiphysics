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

	double prandtl;
	
	double bodydflt[2] = {0.0,0.0};
	if (!input.get(gbl->idprefix +"_body_force",gbl->body.data(),2)) input.getwdefault("body_force",gbl->body.data(),2,bodydflt); 

	if (!input.get(gbl->idprefix + "_gamma",gbl->gamma)) input.getwdefault("gamma",gbl->gamma,1.4);
	if (!input.get(gbl->idprefix + "_mu",gbl->mu)) input.getwdefault("mu",gbl->mu,1.0);
	if (!input.get(gbl->idprefix + "_prandtl",prandtl)) input.getwdefault("prandtl",prandtl,0.713);
	if (!input.get(gbl->idprefix + "_R",gbl->R)) input.getwdefault("R",gbl->R,287.058);
	
	gbl->kcond = gbl->R*gbl->mu/prandtl*gbl->gamma/(gbl->gamma-1.0);

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



