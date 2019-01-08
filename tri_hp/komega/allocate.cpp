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
    
    return;
}

void tri_hp_komega::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
    const tri_hp_komega& inmesh = dynamic_cast<const tri_hp_komega &>(in);
    gbl = inmesh.gbl;
    tri_hp_ins::init(in,why,sizereduce1d);
    return;
}
