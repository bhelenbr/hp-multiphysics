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

 void tri_hp_swe::init(input_map& input, void *gin) {
    
    gbl = static_cast<global *>(gin);    
    tri_hp_ins::init(input,gin);

    if (!input.get(gbl->idprefix + "_f0",gbl->f0)) input.getwdefault("f0",gbl->f0,0.0);
    if (!input.get(gbl->idprefix + "_beta",gbl->beta)) input.getwdefault("beta",gbl->cbeta,0.0);
    if (!input.get(gbl->idprefix + "_cd",gbl->cd)) input.getwdefault("cd",gbl->cd,0.0);
    if (!input.get(gbl->idprefix + "_ptest",gbl->ptest)) input.getwdefault("ptest",gbl->ptest,1.0);

    gbl->bathy = getnewbathy("bathy",input);
    
    return;
}

void tri_hp_swe::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
    const tri_hp_swe& inmesh = dynamic_cast<const tri_hp_swe &>(in);
    gbl = inmesh.gbl;
    tri_hp_ins::init(in,why,sizereduce1d);
    return;
}
