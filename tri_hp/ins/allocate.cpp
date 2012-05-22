/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ins.h"
#include "../hp_boundary.h"

void tri_hp_ins::init(input_map& input, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	gbl = static_cast<global *>(gin);

	if (input.find(gbl->idprefix + "_nvariable") == input.end()) {
		input[gbl->idprefix + "_nvariable"] = "3";
	}

	tri_hp::init(input,gin);

	input.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);

	gbl->tau.resize(maxpst,NV);

	if (!input.get(gbl->idprefix + "_rho",gbl->rho)) input.getwdefault("rho",gbl->rho,1.0);
	if (!input.get(gbl->idprefix + "_mu",gbl->mu)) input.getwdefault("mu",gbl->mu,0.0);
	
	/* LEAVE UP TO DERIVED CLASSES TO LOAD THESE IF NECESSARY */
	gbl->D.resize(NV);
	if (NV > 3) {
		for (int n=2;n<NV-1;++n) {
			stringstream nstr;
			nstr << n-2;
			if (!input.get(gbl->idprefix + "_D" +nstr.str(),gbl->D(n))) 
				if (!input.get("D" +nstr.str(),gbl->D(n)))
					gbl->D(n) = 0.0;
		}
	}

#ifdef DROP
	*gbl->log << "#DROP is defined" << std::endl;
#endif

	return;
}

void tri_hp_ins::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	const tri_hp_ins& inmesh = dynamic_cast<const tri_hp_ins &>(in);
	gbl = inmesh.gbl;

	tri_hp::init(in,why,sizereduce1d);

	adis = inmesh.adis;

	return;
}



/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tri_hp_ins::calculate_unsteady_sources() {
	int i,j,n,tind;

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
					cjcb(i,j) = -gbl->bd(0)*gbl->rho*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
					for(n=0;n<NV-1;++n)
						dugdt(log2p,tind,n)(i,j) = u(n)(i,j)*cjcb(i,j);
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
