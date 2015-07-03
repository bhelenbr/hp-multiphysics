/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cd.h"

void tet_hp_cd::init(input_map& inmap, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	gbl = static_cast<global *>(gin);
	keyword = gbl->idprefix + "_nvariable";
	inmap[keyword] = "1";
	
	tet_hp::init(inmap,gin);
	
	keyword = gbl->idprefix + "_dissipation";
	inmap.getwdefault(keyword,adis,1.0);
	
	FLT rho;
	keyword = gbl->idprefix + "_rho";
	if (!inmap.get(keyword,rho)) inmap.getwdefault("rho",rho,1.0);
	
	FLT cv;
	keyword = gbl->idprefix + "_cv";
	if (!inmap.get(keyword,cv)) inmap.getwdefault("cv",cv,1.0);
	gbl->rhocv = rho*cv;
	
	keyword = gbl->idprefix + "_ax";
	if (!inmap.get(keyword,gbl->ax)) inmap.getwdefault("ax",gbl->ax,1.0);

	keyword = gbl->idprefix + "_ay";
	if (!inmap.get(keyword,gbl->ay)) inmap.getwdefault("ay",gbl->ay,0.0);
	
	keyword = gbl->idprefix + "_az";
	if (!inmap.get(keyword,gbl->ay)) inmap.getwdefault("az",gbl->az,0.0);

	if (!inmap.get(gbl->idprefix + "_nu",gbl->kcond)) {
		if (!inmap.get("nu",gbl->kcond)) {
			if (!inmap.get(gbl->idprefix + "_conductivity",gbl->kcond)) {
				inmap.getwdefault("conductivity",gbl->kcond,0.0);
			}
		}
	}
	
	gbl->tau.resize(maxvst);
	
	keyword = gbl->idprefix + "_src";
	std::string ibcname;
	if (!inmap.get(keyword,ibcname)) {
		keyword = "src";
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find cd src" << std::endl;
		}
	}
	gbl->src = getnewibc(ibcname);
	gbl->src->init(inmap,keyword);
	
	return;
}

void tet_hp_cd::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	const tet_hp_cd& inmesh = dynamic_cast<const tet_hp_cd &>(in);
	gbl = inmesh.gbl;

	tet_hp::init(in,why,sizereduce1d);
	
	adis = inmesh.adis;
		
	return;
}

/* A GENERIC CALCULATION OF SOURCES FOR AN AUTONOMOUS SYSTEM IN STANDARD FORM */
/* WILL NEED TO BE OVERRIDDEN FOR SPECIAL CASES */
void tet_hp_cd::calculate_unsteady_sources() {
	int lgpx = basis::tet(log2p).gpx, lgpy = basis::tet(log2p).gpy, lgpz = basis::tet(log2p).gpz;
	int i,j,k,n,tind;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	TinyVector<int,4> v;
	
	for (log2p=0;log2p<=log2pmax;++log2p) {
		for(tind=0;tind<ntet;++tind) {
			v = tet(tind).pnt;
			
			if (tet(tind).info > -1) {
				crdtocht(tind,1);
				for(n=0;n<ND;++n)
					basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
			}
			else {
				for(n=0;n<ND;++n)
					basis::tet(log2p).proj(vrtxbd(1)(v(0))(n),vrtxbd(1)(v(1))(n),vrtxbd(1)(v(2))(n),vrtxbd(1)(v(3))(n),&crd(n)(0)(0)(0),stridex,stridey);
				
				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpy;++j) {
						for(k=0;k<lgpz;++k) {
							for(n=0;n<ND;++n) {
								dcrd(n)(0)(i)(j)(k) = 0.5*(vrtxbd(1)(v(3))(n) -vrtxbd(1)(v(2))(n));
								dcrd(n)(1)(i)(j)(k) = 0.5*(vrtxbd(1)(v(1))(n) -vrtxbd(1)(v(2))(n));
								dcrd(n)(2)(i)(j)(k) = 0.5*(vrtxbd(1)(v(0))(n) -vrtxbd(1)(v(2))(n));
							}
						}
					}
				}
			}
			
			ugtouht(tind,1);
			basis::tet(log2p).proj(&uht(0)(0),&u(0)(0)(0)(0),stridex, stridey);
			
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						cjcb(i)(j)(k) = -gbl->bd(0)*(dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k)));
						
						dugdt(log2p,tind,0)(i)(j)(k) = gbl->rhocv*u(0)(i)(j)(k)*cjcb(i)(j)(k);
						for(n=0;n<ND;++n)
							dxdt(log2p,tind,n)(i)(j)(k) = crd(n)(i)(j)(k);
					}
				}
			}
		}
		
	}
	log2p = log2pmax;
	
	return;
}


