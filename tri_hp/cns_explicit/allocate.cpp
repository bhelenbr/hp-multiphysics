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

void tri_hp_cns_explicit::init(input_map& inmap, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	gbl = static_cast<global *>(gin);

	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end()) {
		inmap[gbl->idprefix + "_nvariable"] = "4";
	}

	tri_hp::init(inmap,gin);

	inmap.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);


	gbl->tau.resize(maxpst,NV,NV);

	double prandtl;
	if (!inmap.get(gbl->idprefix + "_gamma",gbl->gamma)) inmap.getwdefault("gamma",gbl->gamma,1.4);
	if (!inmap.get(gbl->idprefix + "_mu",gbl->mu)) inmap.getwdefault("mu",gbl->mu,1.0);
	if (!inmap.get(gbl->idprefix + "_prandtl",prandtl)) inmap.getwdefault("prandtl",prandtl,0.713);
	if (!inmap.get(gbl->idprefix + "_R",gbl->R)) inmap.getwdefault("R",gbl->R,287.058);
	
	gbl->kcond = gbl->R*gbl->mu/prandtl*gbl->gamma/(gbl->gamma-1.0);

	/* source term for MMS */
//	std::string ibcname;
//	keyword = gbl->idprefix + "_src";
//	if (!inmap.get(keyword,ibcname)) {
//		keyword = "src";
//		if (!inmap.get(keyword,ibcname)) {
//			*gbl->log << "couldn't find mms source" << std::endl;
//		}
//	}
//	gbl->src = getnewibc(ibcname);
//	gbl->src->init(inmap,keyword);
	
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



