/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cns.h"
#include "../hp_boundary.h"

void tri_hp_cns::init(input_map& inmap, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	gbl = static_cast<global *>(gin);

	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end()) {
		inmap[gbl->idprefix + "_nvariable"] = "4";
	}

	tri_hp::init(inmap,gin);

	inmap.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);

	gbl->tau.resize(maxpst,NV,NV);

	gbl->res_temp.v.resize(maxpst,NV);
	gbl->res_temp.s.resize(maxpst,sm0,NV);
	gbl->res_temp.i.resize(maxpst,im0,NV);
	
	gbl->vpreconditioner.resize(maxpst,NV,NV);
	gbl->spreconditioner.resize(maxpst,NV,NV);
	gbl->tpreconditioner.resize(maxpst,NV,NV);

	double prandtl;
	if (!inmap.get(gbl->idprefix + "_gamma",gbl->gamma)) inmap.getwdefault("gamma",gbl->gamma,1.4);
	if (!inmap.get(gbl->idprefix + "_mu",gbl->mu)) inmap.getwdefault("mu",gbl->mu,1.0);
	if (!inmap.get(gbl->idprefix + "_prandtl",prandtl)) inmap.getwdefault("prandtl",prandtl,0.713);
	if (!inmap.get(gbl->idprefix + "_R",gbl->R)) inmap.getwdefault("R",gbl->R,287.058);

	gbl->kcond = gbl->R*gbl->mu/prandtl*gbl->gamma/(gbl->gamma-1.0);

	/* source term for MMS */
//	std::string ibcname;
//	keyword = gbl->idprefix + "_src";
//	if (!inmap.get(keyword,ibcname)) {
//		keyword = "src";
//		if (!inmap.get(keyword,ibcname)) {
//			*gbl->log << "couldn't find src" << std::endl;
//		}
//	}
//	gbl->ibc = getnewibc(ibcname);
//	gbl->ibc->init(inmap,keyword);

}

void tri_hp_cns::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	const tri_hp_cns& inmesh = dynamic_cast<const tri_hp_cns &>(in);
	gbl = inmesh.gbl;

	tri_hp::init(in,why,sizereduce1d);

	adis = inmesh.adis;

	return;
}



/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tri_hp_cns::calculate_unsteady_sources() {
	int i,j,n,tind;
	FLT	ogm1 = 1.0/(gbl->gamma-1.0);
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
					double rho = u(0)(i,j)/u(NV-1)(i,j);
					cjcb(i,j) = -gbl->bd(0)*rho*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
					dugdt(log2p)(tind,0,i,j) = cjcb(i,j);
					for(n=1;n<NV-1;++n)
						dugdt(log2p)(tind,n,i,j) = u(n)(i,j)*cjcb(i,j);

					double e = ogm1*u(NV-1)(i,j) +0.5*(u(1)(i,j)*u(1)(i,j) +u(2)(i,j)*u(2)(i,j));
					dugdt(log2p)(tind,NV-1,i,j) = e*cjcb(i,j);
					for(n=0;n<ND;++n)
						dxdt(log2p)(tind,n,i,j) = crd(n)(i,j);
				}				
			}
		}
	}
	log2p=log2pmax;

	return;
}
