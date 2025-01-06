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


void tet_hp_cns_explicit::init(input_map& inmap, shared_ptr<block_global> gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	gbl = gin;
	
	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end()) {
		inmap[gbl->idprefix + "_nvariable"] = "5";
	}
	
	tet_hp::init(inmap,gin);
	
	inmap.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);
		
	hp_cns_explicit_gbl->tau.resize(maxvst,NV,NV);

	double prandtl;
	
	double bodydflt[3] = {0.0,0.0,0.0};
	if (!inmap.get(gbl->idprefix +"_body_force",gbl->body.data(),3)) inmap.getwdefault("body_force",gbl->body.data(),3,bodydflt); 

	if (!inmap.get(gbl->idprefix + "_gamma",hp_cns_explicit_gbl->gamma)) inmap.getwdefault("gamma",hp_cns_explicit_gbl->gamma,1.4);
	if (!inmap.get(gbl->idprefix + "_mu",hp_cns_explicit_gbl->mu)) inmap.getwdefault("mu",hp_cns_explicit_gbl->mu,1.0);
	if (!inmap.get(gbl->idprefix + "_prandtl",prandtl)) inmap.getwdefault("prandtl",prandtl,0.713);
	if (!inmap.get(gbl->idprefix + "_R",hp_cns_explicit_gbl->R)) inmap.getwdefault("R",hp_cns_explicit_gbl->R,287.058);
	
	hp_cns_explicit_gbl->kcond = hp_cns_explicit_gbl->R*hp_cns_explicit_gbl->mu/prandtl*hp_cns_explicit_gbl->gamma/(hp_cns_explicit_gbl->gamma-1.0);

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

