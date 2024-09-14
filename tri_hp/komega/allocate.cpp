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

void tri_hp_komega::init(input_map& inmap, shared_ptr<block_global> gin) {
	gbl = gin;
    hp_komega_gbl = make_shared<hp_komega_global>();
	inmap[gbl->idprefix + "_nvariable"] = "5";
	tri_hp_ins::init(inmap,gin);
	
	if (!inmap.get(gbl->idprefix + "_linf",hp_komega_gbl->linf)) inmap.getwdefault("linf",hp_komega_gbl->linf,1.0);
	if (!inmap.get(gbl->idprefix + "_uinf",hp_komega_gbl->uinf)) inmap.getwdefault("uinf",hp_komega_gbl->uinf,1.0);
	if (!inmap.get(gbl->idprefix + "_epslnk",hp_komega_gbl->epslnk)) inmap.getwdefault("epslnk",hp_komega_gbl->epslnk,1.0);
    if (!inmap.get(gbl->idprefix + "_kinf",hp_komega_gbl->kinf)) inmap.getwdefault("kinf",hp_komega_gbl->kinf,1.0);
    if (!inmap.get(gbl->idprefix + "_omginf",hp_komega_gbl->omginf)) inmap.getwdefault("omginf",hp_komega_gbl->omginf,1.0);
    if (!inmap.get(gbl->idprefix + "_Clim",hp_komega_gbl->Clim)) inmap.getwdefault("omginf",hp_komega_gbl->Clim,7.0/8.0);
    if (!inmap.get(gbl->idprefix + "_version",hp_komega_gbl->version)) inmap.getwdefault("version",hp_komega_gbl->version,3);
    if (!inmap.get(gbl->idprefix + "_kmom_on",hp_komega_gbl->kmom_on)) inmap.getwdefault("kmom_on",hp_komega_gbl->kmom_on,1);
    if (!inmap.get(gbl->idprefix + "_sust_on",hp_komega_gbl->sust_on)) inmap.getwdefault("sust_on",hp_komega_gbl->sust_on,1);
    
#ifdef MMS
    /* source term for MMS */
    std::string ibcname, keyword;
    keyword = gbl->idprefix + "_src";
    if (!inmap.get(keyword,ibcname)) {
        keyword = "src";
        if (!inmap.get(keyword,ibcname)) {
            *gbl->log << "couldn't find src" << std::endl;
        }
    }
    hp_komega_gbl->src = getnewibc(ibcname);
    hp_komega_gbl->src->init(inmap,keyword);
#endif
   
    return;
}

void tri_hp_komega::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_hp_komega& inmesh = dynamic_cast<const tri_hp_komega &>(in);
	gbl = inmesh.gbl;
	tri_hp_ins::init(in,why,sizereduce1d);
	return;
}
