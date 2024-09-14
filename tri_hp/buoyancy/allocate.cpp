/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_buoyancy.h"
#include "../hp_boundary.h"

void tri_hp_buoyancy::init(input_map& inmap, shared_ptr<block_global> gin) {
    gbl = gin;
    hp_buoyancy_gbl = make_shared<hp_buoyancy_global>();
    
	inmap[gbl->idprefix + "_nvariable"] = "4";
	tri_hp_ins::init(inmap,gin);

	if (!inmap.get(gbl->idprefix + "_conductivity",hp_buoyancy_gbl->kcond)) inmap.getwdefault("conductivity",hp_buoyancy_gbl->kcond,0.7*hp_ins_gbl->mu);
	hp_ins_gbl->D(0) = hp_buoyancy_gbl->kcond;
	if (!inmap.get(gbl->idprefix + "_cp",hp_buoyancy_gbl->cp)) inmap.getwdefault("cp",hp_buoyancy_gbl->cp,1.0);

	if (inmap.find(gbl->idprefix+"_rho_vs_T") != inmap.end()) {
		hp_buoyancy_gbl->rho_vs_T.init(inmap,gbl->idprefix+"_rho_vs_T");
	} 
	else if (inmap.find("rho_vs_T") != inmap.end()){
		hp_buoyancy_gbl->rho_vs_T.init(inmap,"rho_vs_T");
	}
	else {
		*gbl->log << "couldn't find rho_vs_T equation for density" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	if (!inmap.get(gbl->idprefix + "_energy_scaling",hp_buoyancy_gbl->adapt_energy_scaling)) inmap.getwdefault("energy_scaling",hp_buoyancy_gbl->adapt_energy_scaling,1.0);


	return;
}

void tri_hp_buoyancy::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_hp_buoyancy& inmesh = dynamic_cast<const tri_hp_buoyancy &>(in);
	gbl = inmesh.gbl;
	tri_hp_ins::init(in,why,sizereduce1d);
	return;
}



/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tri_hp_buoyancy::calculate_unsteady_sources() {
    int i,j,n,tind;
    FLT lrho;
#ifdef petsc
	int start = log2pmax;
#else
	int start = 0;
#endif
	
	for (log2p=start;log2p<=log2pmax;++log2p) {
		for(tind=0;tind<ntri;++tind) {
			if (tri(tind).info > -1) {
				crdtocht(tind,1);
				for(n=0;n<ND;++n)
					basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
			}
			else {
				for(n=0;n<ND;++n)
					basis::tri(log2p)->proj(vrtxbd(1)(tri(tind).pnt(0))(n),vrtxbd(1)(tri(tind).pnt(1))(n),vrtxbd(1)(tri(tind).pnt(2))(n),&crd(n)(0,0),MXGP);

				for(i=0;i<basis::tri(log2p)->gpx();++i) {
					for(j=0;j<basis::tri(log2p)->gpn();++j) {
						for(n=0;n<ND;++n) {
							dcrd(n,0)(i,j) = 0.5*(vrtxbd(1)(tri(tind).pnt(1))(n) -vrtxbd(1)(tri(tind).pnt(0))(n));
							dcrd(n,1)(i,j) = 0.5*(vrtxbd(1)(tri(tind).pnt(2))(n) -vrtxbd(1)(tri(tind).pnt(0))(n));
						}
					}
				}
			}

			ugtouht(tind,1);
			for(n=0;n<NV;++n)
				basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);

			for(i=0;i<basis::tri(log2p)->gpx();++i) {
				for(j=0;j<basis::tri(log2p)->gpn();++j) {    
					cjcb(i,j) = -gbl->bd(0)*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
					lrho = hp_buoyancy_gbl->rho_vs_T.Eval(u(2)(i,j));                    
					for(n=0;n<NV-2;++n)
						dugdt(log2p)(tind,n,i,j) = lrho*u(n)(i,j)*cjcb(i,j);
					dugdt(log2p)(tind,NV-2,i,j) = lrho*hp_buoyancy_gbl->cp*u(NV-2)(i,j)*cjcb(i,j);
					dugdt(log2p)(tind,NV-1,i,j) = lrho*cjcb(i,j);

					for(n=0;n<ND;++n)
						dxdt(log2p)(tind,n,i,j) = crd(n)(i,j);
				}				
			}
		}
    }
    log2p=log2pmax;

    return;
}
