/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cns_explicit.h"
#include "../hp_boundary.h"


void tet_hp_cns_explicit::init(input_map& inmap, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	gbl = static_cast<global *>(gin);
	
	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end()) {
		inmap[gbl->idprefix + "_nvariable"] = "5";
	}
	
	tet_hp::init(inmap,gin);
	
	inmap.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);
		
	gbl->tau.resize(maxvst,NV,NV);

	double prandtl;
	
	double bodydflt[3] = {0.0,0.0,0.0};
	if (!inmap.get(gbl->idprefix +"_body_force",gbl->body.data(),3)) inmap.getwdefault("body_force",gbl->body.data(),3,bodydflt); 

	if (!inmap.get(gbl->idprefix + "_gamma",gbl->gamma)) inmap.getwdefault("gamma",gbl->gamma,1.4);
	if (!inmap.get(gbl->idprefix + "_mu",gbl->mu)) inmap.getwdefault("mu",gbl->mu,1.0);
	if (!inmap.get(gbl->idprefix + "_prandtl",prandtl)) inmap.getwdefault("prandtl",prandtl,0.713);
	if (!inmap.get(gbl->idprefix + "_R",gbl->R)) inmap.getwdefault("R",gbl->R,287.058);
	
	gbl->kcond = gbl->R*gbl->mu/prandtl*gbl->gamma/(gbl->gamma-1.0);

	/* source term for MMS */
	//	keyword = gbl->idprefix + "_src";
	//	std::string ibcname;
	//	if (!inmap.get(keyword,ibcname)) {
	//		keyword = "src";
	//		if (!inmap.get(keyword,ibcname)) {
	//			*gbl->log << "couldn't find cd velocity field" << std::endl;
	//		}
	//	}
	//	gbl->src = getnewibc(ibcname);
	//  gbl->src->init(inmap,keyword);
	
	return;
}

void tet_hp_cns_explicit::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	const tet_hp_cns_explicit& inmesh = dynamic_cast<const tet_hp_cns_explicit &>(in);
	gbl = inmesh.gbl;

	tet_hp::init(in,why,sizereduce1d);
	
	adis = inmesh.adis;
		
	return;
}

