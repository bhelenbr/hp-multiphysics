/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"

void tri_hp_cd::init(input_map& inmap, shared_ptr<block_global> gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	std::string ibcname;

	gbl = gin;
    hp_cd_gbl = make_shared<hp_cd_global>();
	keyword = gbl->idprefix + "_nvariable";
	inmap[keyword] = "1";

	tri_hp::init(inmap,gin);

	keyword = gbl->idprefix + "_dissipation";
	inmap.getwdefault(keyword,adis,1.0);
	
	FLT rho;
	keyword = gbl->idprefix + "_rho";
	if (!inmap.get(keyword,rho)) inmap.getwdefault("rho",rho,1.0);
	
	FLT cv;
	keyword = gbl->idprefix + "_cv";
	if (!inmap.get(keyword,cv)) inmap.getwdefault("cv",cv,1.0);
	hp_cd_gbl->rhocv = rho*cv;

#ifdef CONST_A
	keyword = gbl->idprefix + "_ax";
	if (!inmap.get(keyword,hp_cd_gbl->ax)) inmap.getwdefault("ax",hp_cd_gbl->ax,1.0);
	
	keyword = gbl->idprefix + "_ay";
	if (!inmap.get(keyword,hp_cd_gbl->ay)) inmap.getwdefault("ay",hp_cd_gbl->ay,0.0);
#else
	/* FIND INITIAL CONDITION TYPE */
	keyword = gbl->idprefix + "_a";
	if (!inmap.get(keyword,ibcname)) {
		keyword = "a";
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find cd velocity field" << std::endl;
		}
	}
	hp_cd_gbl->a = getnewibc(ibcname);
	hp_cd_gbl->a->init(inmap,keyword);
#endif

	if (!inmap.get(gbl->idprefix + "_nu",hp_cd_gbl->kcond)) {
		if (!inmap.get("nu",hp_cd_gbl->kcond)) {
			if (!inmap.get(gbl->idprefix + "_conductivity",hp_cd_gbl->kcond)) {
				inmap.getwdefault("conductivity",hp_cd_gbl->kcond,0.0);
			}
		}
	}
	hp_cd_gbl->tau.resize(maxpst);

	/* FIND INITIAL CONDITION TYPE */
	keyword = gbl->idprefix + "_src";
	if (!inmap.get(keyword,ibcname)) {
		keyword = "src";
		if (!inmap.get(keyword,ibcname)) {
			ibcname = "zero";
		}
	}
	hp_cd_gbl->src = getnewibc(ibcname);
	hp_cd_gbl->src->init(inmap,keyword);
	
	/* Stuff for Mike's minvrt */
//	gbl->stiff_diag.v.resize(maxpst,NV);
//	gbl->stiff_diag.s.resize(maxpst,sm0,NV);
//	gbl->stiff_diag.i.resize(maxpst,im0,NV);
		
	return;
}

void tri_hp_cd::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	const tri_hp_cd& inmesh = dynamic_cast<const tri_hp_cd &>(in);
	gbl = inmesh.gbl;

	tri_hp::init(in,why,sizereduce1d);

	adis = inmesh.adis;

	return;
}

/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tri_hp_cd::calculate_unsteady_sources() {
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
			basis::tri(log2p)->proj(&uht(0)(0),&u(0)(0,0),MXGP);
			
			for(i=0;i<basis::tri(log2p)->gpx();++i) {
				for(j=0;j<basis::tri(log2p)->gpn();++j) {
					cjcb(i,j) = -gbl->bd(0)*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
					dugdt(log2p)(tind,0,i,j) = hp_cd_gbl->rhocv*u(0)(i,j)*cjcb(i,j);
					
					for(n=0;n<ND;++n)
						dxdt(log2p)(tind,n,i,j) = crd(n)(i,j);
				}
			}
		}
	}
	log2p=log2pmax;
	
	return;
}

