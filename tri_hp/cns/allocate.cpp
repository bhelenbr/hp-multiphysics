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

void tri_hp_cns::init(input_map& inmap, shared_ptr<block_global> gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	gbl = gin;
    hp_cns_gbl = make_shared<hp_cns_global>();
    
	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end()) {
		inmap[gbl->idprefix + "_nvariable"] = "4";
	}

	tri_hp::init(inmap,gin);

	inmap.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);

	hp_cns_gbl->tau.resize(maxpst,NV,NV);

	hp_cns_gbl->res_temp.v.resize(maxpst,NV);
	hp_cns_gbl->res_temp.s.resize(maxpst,sm0,NV);
	hp_cns_gbl->res_temp.i.resize(maxpst,im0,NV);
	
	hp_cns_gbl->vpreconditioner.resize(maxpst,NV,NV);
	hp_cns_gbl->spreconditioner.resize(maxpst,NV,NV);
	hp_cns_gbl->tpreconditioner.resize(maxpst,NV,NV);

	if (!inmap.get(gbl->idprefix + "_gamma",hp_cns_gbl->gamma)) inmap.getwdefault("gamma",hp_cns_gbl->gamma,1.4);
	if (!inmap.get(gbl->idprefix + "_mu",hp_cns_gbl->mu)) inmap.getwdefault("mu",hp_cns_gbl->mu,1.716e-5);
	if (!inmap.get(gbl->idprefix + "_prandtl",hp_cns_gbl->prandtl)) inmap.getwdefault("prandtl",hp_cns_gbl->prandtl,0.713);
	if (!inmap.get(gbl->idprefix + "_Rgas",hp_cns_gbl->R)) inmap.getwdefault("Rgas",hp_cns_gbl->R,287.058);

    /* Pr = mu/(k/cp) with cp = gamma/(gamma-1)*R */
	hp_cns_gbl->kcond = hp_cns_gbl->mu/hp_cns_gbl->prandtl*hp_cns_gbl->R*hp_cns_gbl->gamma/(hp_cns_gbl->gamma-1.0);

#ifdef SUTHERLAND
    *gbl->log << "#SUTHERLAND is defined" << std::endl;
    /* Inputs for Sutherland Viscsosity */
    /* Formula is mu0*(T/T0)^(3/2)*(T0+C)/(T+C)

    For air from CFD on-line
    https://www.cfd-online.com/Wiki/Sutherland%27s_law
    Also agrees with COMSOL
    https://doc.comsol.com/5.5/doc/com.comsol.help.cfd/cfd_ug_fluidflow_high_mach.08.27.html
    mu0 = 1.716e-5; % [kg/(m s)]
    T0 = 273.15 % [K]
    C = 110.4; % [K]
    R = 8314/28.97;
    */
    
    FLT T0,C;
    if (!inmap.get(gbl->idprefix +"_Sutherland_T0",T0))
        inmap.getwdefault("Sutherland_T0",T0,273.15);
    
    if  (!inmap.get(gbl->idprefix +"_Sutherland_C",C))
        inmap.getwdefault("Sutherland_C",C,110.4);
    
    /* Now convert so they can be used with RT instead of T */
    hp_cns_gbl->s1 = hp_cns_gbl->mu/pow(hp_cns_gbl->R*T0,1.5)*(T0+C)*hp_cns_gbl->R;
    hp_cns_gbl->s2 = hp_cns_gbl->R*C;
#endif
    
    /* source term for MMS */
    std::string ibcname;
    keyword = gbl->idprefix + "_src";
    if (!inmap.get(keyword,ibcname)) {
        keyword = "src";
        if (!inmap.get(keyword,ibcname)) {
            ibcname="zero";
        }
    }
    hp_cns_gbl->src = getnewibc(ibcname);
    hp_cns_gbl->src->init(inmap,keyword);

}

void tri_hp_cns::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	const tri_hp_cns& inmesh = dynamic_cast<const tri_hp_cns &>(in);
	gbl = inmesh.gbl;
    hp_cns_gbl = inmesh.hp_cns_gbl;

	tri_hp::init(in,why,sizereduce1d);

	adis = inmesh.adis;

	return;
}



/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tri_hp_cns::calculate_unsteady_sources() {
	int i,j,n,tind;
	FLT	ogm1 = 1.0/(hp_cns_gbl->gamma-1.0);
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
