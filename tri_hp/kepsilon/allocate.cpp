/*
 *  allocate.cpp
 *  spectral_hp
 *
 *  Standard (Fluent) k-epsilon turbulence model.
 *
 */

#include "tri_hp_kepsilon.h"
#include "../hp_boundary.h"

void tri_hp_kepsilon::init(input_map& inmap, shared_ptr<block_global> gin) {
	gbl = gin;
    hp_kepsilon_gbl = make_shared<hp_kepsilon_global>();
	inmap[gbl->idprefix + "_nvariable"] = "5";
	tri_hp_ins::init(inmap,gin);

	if (!inmap.get(gbl->idprefix + "_linf",hp_kepsilon_gbl->linf)) inmap.getwdefault("linf",hp_kepsilon_gbl->linf,1.0);
	if (!inmap.get(gbl->idprefix + "_uinf",hp_kepsilon_gbl->uinf)) inmap.getwdefault("uinf",hp_kepsilon_gbl->uinf,1.0);
    if (!inmap.get(gbl->idprefix + "_kinf",hp_kepsilon_gbl->kinf)) inmap.getwdefault("kinf",hp_kepsilon_gbl->kinf,1.0);
    if (!inmap.get(gbl->idprefix + "_epsinf",hp_kepsilon_gbl->epsinf)) inmap.getwdefault("epsinf",hp_kepsilon_gbl->epsinf,1.0);

    /* Floors to protect divisions by k and epsilon (k and epsilon are solved directly) */
    if (!inmap.get(gbl->idprefix + "_kmin",hp_kepsilon_gbl->kmin)) inmap.getwdefault("kmin",hp_kepsilon_gbl->kmin,1.0e-10);
    if (!inmap.get(gbl->idprefix + "_epsmin",hp_kepsilon_gbl->epsmin)) inmap.getwdefault("epsmin",hp_kepsilon_gbl->epsmin,1.0e-10);

    /* Standard k-epsilon model constants (Fluent defaults) */
    if (!inmap.get(gbl->idprefix + "_Cmu",hp_kepsilon_gbl->Cmu)) inmap.getwdefault("Cmu",hp_kepsilon_gbl->Cmu,0.09);
    if (!inmap.get(gbl->idprefix + "_C1eps",hp_kepsilon_gbl->C1eps)) inmap.getwdefault("C1eps",hp_kepsilon_gbl->C1eps,1.44);
    if (!inmap.get(gbl->idprefix + "_C2eps",hp_kepsilon_gbl->C2eps)) inmap.getwdefault("C2eps",hp_kepsilon_gbl->C2eps,1.92);
    if (!inmap.get(gbl->idprefix + "_sigmak",hp_kepsilon_gbl->sigmak)) inmap.getwdefault("sigmak",hp_kepsilon_gbl->sigmak,1.0);
    if (!inmap.get(gbl->idprefix + "_sigmaeps",hp_kepsilon_gbl->sigmaeps)) inmap.getwdefault("sigmaeps",hp_kepsilon_gbl->sigmaeps,1.3);

    if (!inmap.get(gbl->idprefix + "_kmom_on",hp_kepsilon_gbl->kmom_on)) inmap.getwdefault("kmom_on",hp_kepsilon_gbl->kmom_on,1);

    /* source term for MMS */
    std::string ibcname, keyword;
    keyword = gbl->idprefix + "_src";
    if (!inmap.get(keyword,ibcname)) {
        keyword = "src";
        if (!inmap.get(keyword,ibcname)) {
            ibcname = "zero";
        }
    }
    hp_kepsilon_gbl->src = getnewibc(ibcname);
    hp_kepsilon_gbl->src->init(inmap,keyword);

    return;
}

void tri_hp_kepsilon::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_hp_kepsilon& inmesh = dynamic_cast<const tri_hp_kepsilon &>(in);
	gbl = inmesh.gbl;
    hp_kepsilon_gbl = inmesh.hp_kepsilon_gbl;
	tri_hp_ins::init(in,why,sizereduce1d);
	return;
}
