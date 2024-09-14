/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ps.h"
#include "../hp_boundary.h"

void tri_hp_ps::init(input_map& inmap, shared_ptr<block_global> gin) {
	FLT nu, E;

    gbl = gin;
    hp_ps_gbl = make_shared<hp_ps_global>();
    
	inmap[gbl->idprefix + "_nvariable"] = "3";
	tri_hp::init(inmap,gin);


	hp_ps_gbl->tau.resize(maxpst);
	inmap.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);
	inmap.getwdefault(gbl->idprefix + "_nu",nu,0.3);
	inmap.getwdefault(gbl->idprefix + "_E",E,1.0);

	hp_ps_gbl->mu = E/(2.*(1.+nu));
	hp_ps_gbl->lami = (1.+nu)*(1.-2.*nu)/(E*nu);

	return;
}

void tri_hp_ps::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {

	const tri_hp_ps& inmesh = dynamic_cast<const tri_hp_ps &>(in);
	gbl = inmesh.gbl;

	tri_hp::init(in,why,sizereduce1d);

	adis = inmesh.adis;

	return;
}

/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW TO DO NOTHING */
void tri_hp_ps::calculate_unsteady_sources() {}
