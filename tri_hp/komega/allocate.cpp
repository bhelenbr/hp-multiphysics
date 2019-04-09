/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_komega.h"
#include "../hp_boundary.h"

void tri_hp_komega::init(input_map& inmap, void *gin) {
	gbl = static_cast<global *>(gin);
	inmap[gbl->idprefix + "_nvariable"] = "5";
	tri_hp_ins::init(inmap,gin);
	
	if (!inmap.get(gbl->idprefix + "_linf",gbl->linf)) inmap.getwdefault("linf",gbl->linf,1.0);
	if (!inmap.get(gbl->idprefix + "_uinf",gbl->uinf)) inmap.getwdefault("uinf",gbl->uinf,1.0);
	if (!inmap.get(gbl->idprefix + "_c_mu",gbl->c_mu)) inmap.getwdefault("c_mu",gbl->c_mu,0.09);
	if (!inmap.get(gbl->idprefix + "_epslnk",gbl->epslnk)) inmap.getwdefault("epslnk",gbl->epslnk,1.0);
    
	
	// Derived Constants
	const FLT nu = gbl->mu/gbl->rho; // kinematic viscosity
	const FLT Re = gbl->uinf*gbl->linf/nu; // Reynolds number
	
	// Initial estimate for k and omega based on turbulent length scale and Reynolds number
	const FLT	l_turb = 0.07*gbl->linf/(pow(gbl->c_mu,0.75)); // initial estimate of turbulent length scale
	const FLT I = 0.16*pow(Re,(-1./8.)); // initial estimate of the turbulence intensity
	const FLT k_est = 3./2.*pow(gbl->uinf*I,2); // initial estimate of the mean turbulent kinetic energy
	const FLT omg_est = sqrt(k_est)/gbl->c_mu/l_turb; // initial estimate of mean specific dissipation rate
	
	// Use kinf and omginf based on turbulent lengh scale and Reynolds number
	gbl->kinf = k_est;
	gbl->omginf = omg_est;

	return;
}

void tri_hp_komega::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_hp_komega& inmesh = dynamic_cast<const tri_hp_komega &>(in);
	gbl = inmesh.gbl;
	tri_hp_ins::init(in,why,sizereduce1d);
	return;
}
