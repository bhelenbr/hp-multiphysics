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

void tri_hp_ins::init(input_map& inmap, shared_ptr<block_global> gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	gbl = gin;
    hp_ins_gbl = make_shared<hp_ins_global>();
    if (!inmap.get(gbl->idprefix + "_rho",hp_ins_gbl->rho)) inmap.getwdefault("rho",hp_ins_gbl->rho,1.0);
    if (!inmap.get(gbl->idprefix + "_mu",hp_ins_gbl->mu)) inmap.getwdefault("mu",hp_ins_gbl->mu,0.0);

	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end()) {
		inmap[gbl->idprefix + "_nvariable"] = "3";
	}
    inmap.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);

	tri_hp::init(inmap,gin);

	hp_ins_gbl->tau.resize(maxpst,NV);
    /* LEAVE UP TO DERIVED CLASSES TO LOAD THESE IF NECESSARY */
    hp_ins_gbl->D.resize(NV);
    if (NV > 3) {
        for (int n=2;n<NV-1;++n) {
            stringstream nstr;
            nstr << n-2;
            if (!inmap.get(gbl->idprefix + "_D" +nstr.str(),hp_ins_gbl->D(n)))
                if (!inmap.get("D" +nstr.str(),hp_ins_gbl->D(n)))
                    hp_ins_gbl->D(n) = 0.0;
        }
    }

	return;
}

void tri_hp_ins::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	const tri_hp_ins& inmesh = dynamic_cast<const tri_hp_ins &>(in);
	gbl = inmesh.gbl;
    hp_ins_gbl = inmesh.hp_ins_gbl;

	tri_hp::init(in,why,sizereduce1d);

	adis = inmesh.adis;

	return;
}



/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tri_hp_ins::calculate_unsteady_sources() {
	int i,j,n,tind;
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
					cjcb(i,j) = -gbl->bd(0)*hp_ins_gbl->rho*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
					for(n=0;n<NV-1;++n)
						dugdt(log2p)(tind,n,i,j) = u(n)(i,j)*cjcb(i,j);
					dugdt(log2p)(tind,NV-1,i,j) = cjcb(i,j);

					for(n=0;n<ND;++n)
						dxdt(log2p)(tind,n,i,j) = crd(n)(i,j);
 				
				}	
			}
		}
	}
	log2p=log2pmax;

	return;
}
