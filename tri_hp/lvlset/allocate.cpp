/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_lvlset.h"
#include "../hp_boundary.h"

void tri_hp_lvlset::init(input_map& input, void *gin) {
	gbl = static_cast<global *>(gin);    
	input[gbl->idprefix + "_nvariable"] = "4";
	tri_hp_ins::init(input,gin);

	input.getwdefault(gbl->idprefix + "_reinit",reinit_flag,false);
	input.getwdefault(gbl->idprefix + "_reinit_iterations",reinit_iterations,10);
	input.getwdefault(gbl->idprefix + "_rho2",gbl->rho2,1.0);
	input.getwdefault(gbl->idprefix + "_mu2",gbl->mu2,0.0);
	input.getwdefault(gbl->idprefix + "_sigma",gbl->sigma,0.0);
	input.getwdefault(gbl->idprefix + "_width",gbl->width,0.02);

	return;
}

void tri_hp_lvlset::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_hp_lvlset& inmesh = dynamic_cast<const tri_hp_lvlset &>(in);
	gbl = inmesh.gbl;
	reinit_flag = inmesh.reinit_flag;
	reinit_iterations = inmesh.reinit_iterations;
	tri_hp_ins::init(in,why,sizereduce1d);
	return;
}


void tri_hp_lvlset::calculate_unsteady_sources() {
	int i,j,n,tind;
	FLT rho;

	for (log2p=0;log2p<=log2pmax;++log2p) {
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
					rho = gbl->rho +(gbl->rho2 -gbl->rho)*heavyside_if(u(2)(i,j)/gbl->width);
					cjcb(i,j) = -gbl->bd(0)*rho*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
#ifndef CONSERVATIVE
					for(n=0;n<NV-2;++n)
						dugdt(log2p,tind,n)(i,j) = u(n)(i,j)*cjcb(i,j);
					dugdt(log2p,tind,NV-2)(i,j) = -gbl->bd(0)*u(NV-2)(i,j);
#else
					for(n=0;n<NV-1;++n)
						dugdt(log2p,tind,n)(i,j) = u(n)(i,j)*cjcb(i,j);
#endif
					dugdt(log2p,tind,NV-1)(i,j) = cjcb(i,j);

					for(n=0;n<ND;++n)
						dxdt(log2p,tind,n)(i,j) = crd(n)(i,j);
				}
			}
		}
	}
	log2p=log2pmax;

	return;
}
