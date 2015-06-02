/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swe.h"
#include "../hp_boundary.h"

void tri_hp_swe::init(input_map& inmap, void *gin) {

	gbl = static_cast<global *>(gin);    
	tri_hp_ins::init(inmap,gin);

	if (!inmap.get(gbl->idprefix + "_f0",gbl->f0)) inmap.getwdefault("f0",gbl->f0,0.0);
	if (!inmap.get(gbl->idprefix + "_beta",gbl->beta)) inmap.getwdefault("beta",gbl->cbeta,0.0);
	if (!inmap.get(gbl->idprefix + "_cd",gbl->cd)) inmap.getwdefault("cd",gbl->cd,0.0);
	if (!inmap.get(gbl->idprefix + "_ptest",gbl->ptest)) inmap.getwdefault("ptest",gbl->ptest,1.0);

	/* FIND INITIAL CONDITION TYPE */
	std::string keyword, ibcname;
	keyword = gbl->idprefix + "_bathy";
	if (!inmap.get(keyword,ibcname)) {
		keyword = "bathy";
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find bathymettry function" << std::endl;
		}
	}
	gbl->bathy = getnewibc(ibcname);
	gbl->bathy->init(inmap,keyword);

	return;
}

void tri_hp_swe::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_hp_swe& inmesh = dynamic_cast<const tri_hp_swe &>(in);
	gbl = inmesh.gbl;
	tri_hp_ins::init(in,why,sizereduce1d);
	return;
}
