/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cd.h"

void tet_hp_cd::init(input_map& input, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	gbl = static_cast<global *>(gin);
	keyword = gbl->idprefix + "_nvariable";
	input[keyword] = "1";
	
	tet_hp::init(input,gin);
	
	keyword = gbl->idprefix + "_dissipation";
	input.getwdefault(keyword,adis,1.0);
	
	keyword = gbl->idprefix + "_ax";
	if (!input.get(keyword,gbl->ax)) input.getwdefault("ax",gbl->ax,1.0);

	keyword = gbl->idprefix + "_ay";
	if (!input.get(keyword,gbl->ay)) input.getwdefault("ay",gbl->ay,0.0);
	
	keyword = gbl->idprefix + "_az";
	if (!input.get(keyword,gbl->ay)) input.getwdefault("az",gbl->az,0.0);

	keyword = gbl->idprefix + "_nu";
	if (!input.get(keyword,gbl->nu)) input.getwdefault("nu",gbl->nu,0.0);

	gbl->tau.resize(maxvst);
	
	gbl->src = getnewibc("src",input);
	
	return;
}

void tet_hp_cd::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	const tet_hp_cd& inmesh = dynamic_cast<const tet_hp_cd &>(in);
	gbl = inmesh.gbl;

	tet_hp::init(in,why,sizereduce1d);
	
	adis = inmesh.adis;
		
	return;
}


