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
	if (!inmap.get(gbl->idprefix + "_epslnk",gbl->epslnk)) inmap.getwdefault("epslnk",gbl->epslnk,1.0);
    if (!inmap.get(gbl->idprefix + "_kinf",gbl->kinf)) inmap.getwdefault("kinf",gbl->kinf,1.0);
    if (!inmap.get(gbl->idprefix + "_omginf",gbl->omginf)) inmap.getwdefault("omginf",gbl->omginf,1.0);
    if (!inmap.get(gbl->idprefix + "_Clim",gbl->Clim)) inmap.getwdefault("omginf",gbl->Clim,7.0/8.0);
    if (!inmap.get(gbl->idprefix + "_version",gbl->version)) inmap.getwdefault("version",gbl->version,3);
    if (!inmap.get(gbl->idprefix + "_kmom_on",gbl->kmom_on)) inmap.getwdefault("kmom_on",gbl->kmom_on,1);
    if (!inmap.get(gbl->idprefix + "_sust_on",gbl->sust_on)) inmap.getwdefault("sust_on",gbl->sust_on,1);
    
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
    gbl->src = getnewibc(ibcname);
    gbl->src->init(inmap,keyword);
#endif
   
    return;
}

void tri_hp_komega::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_hp_komega& inmesh = dynamic_cast<const tri_hp_komega &>(in);
	gbl = inmesh.gbl;
	tri_hp_ins::init(in,why,sizereduce1d);
	return;
}
