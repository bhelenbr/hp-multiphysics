/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"

void tri_hp_cd::init(input_map& input, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	gbl = static_cast<global *>(gin);
	keyword = gbl->idprefix + "_nvariable";
	input[keyword] = "1";

	tri_hp::init(input,gin);

	keyword = gbl->idprefix + "_dissipation";
	input.getwdefault(keyword,adis,1.0);

#ifdef CONST_A
	keyword = gbl->idprefix + "_ax";
	if (!input.get(keyword,gbl->ax)) input.getwdefault("ax",gbl->ax,1.0);
	
	keyword = gbl->idprefix + "_ay";
	if (!input.get(keyword,gbl->ay)) input.getwdefault("ay",gbl->ay,0.0);
#else
	gbl->a = getnewibc("a",input);
#endif

	keyword = gbl->idprefix + "_nu";
	if (!input.get(keyword,gbl->nu)) input.getwdefault("nu",gbl->nu,0.0);

	keyword = gbl->idprefix + "_minlngth";
	if (!input.get(keyword,gbl->minlngth)) input.getwdefault("minlngth",gbl->minlngth,-1.0);

	keyword = gbl->idprefix + "_maxlngth";
	if (!input.get(keyword,gbl->maxlngth)) input.getwdefault("maxlngth",gbl->maxlngth,1.0e99);    

	gbl->tau.resize(maxpst);

	gbl->src = getnewibc("src",input);
	
	gbl->stiff_diag.v.resize(maxpst,NV);
	gbl->stiff_diag.s.resize(maxpst,sm0,NV);
	gbl->stiff_diag.i.resize(maxpst,im0,NV);
		
	return;
}

void tri_hp_cd::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	const tri_hp_cd& inmesh = dynamic_cast<const tri_hp_cd &>(in);
	gbl = inmesh.gbl;

	tri_hp::init(in,why,sizereduce1d);

	adis = inmesh.adis;

	return;
}

