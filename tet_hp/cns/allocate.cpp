/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cns.h"
//#include "../hp_boundary.h"


void tet_hp_cns::init(input_map& inmap, shared_ptr<block_global> gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	gbl = gin;
	
	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end()) {
		inmap[gbl->idprefix + "_nvariable"] = "5";
	}
	
	tet_hp::init(inmap,gin);
	
	if (!inmap.get(gbl->idprefix + "_dissipation",adis)) inmap.getwdefault("dissipation",adis,1.0);

	/* no preconditioner = 0, weiss-smith preconditioner = 1, squared preconditioner = 2 */
	inmap.getwdefault("preconditioner",hp_cns_gbl->preconditioner,1);

	hp_cns_gbl->tau.resize(maxvst,NV,NV);

	hp_cns_gbl->vpreconditioner.resize(maxvst,NV,NV);
	hp_cns_gbl->epreconditioner.resize(maxvst,NV,NV);
	hp_cns_gbl->betasquared.resize(ntet);

	double prandtl;
	
	double bodydflt[3] = {0.0,0.0,0.0};
	if (!inmap.get(gbl->idprefix +"_body_force",gbl->body.data(),3)) inmap.getwdefault("body_force",gbl->body.data(),3,bodydflt); 
	
	if (!inmap.get(gbl->idprefix + "_atm_pressure",hp_cns_gbl->atm_pressure)) inmap.getwdefault("atm_pressure",hp_cns_gbl->atm_pressure,0.0);	
	if (!inmap.get(gbl->idprefix + "_density",hp_cns_gbl->density)) inmap.getwdefault("density",hp_cns_gbl->density,0.0);	
	if (!inmap.get(gbl->idprefix + "_gamma",hp_cns_gbl->gamma)) inmap.getwdefault("gamma",hp_cns_gbl->gamma,1.4);
	if (!inmap.get(gbl->idprefix + "_mu",hp_cns_gbl->mu)) inmap.getwdefault("mu",hp_cns_gbl->mu,1.0);
	if (!inmap.get(gbl->idprefix + "_prandtl",prandtl)) inmap.getwdefault("prandtl",prandtl,0.713);
	if (!inmap.get(gbl->idprefix + "_R",hp_cns_gbl->R)) inmap.getwdefault("R",hp_cns_gbl->R,287.058);
	
	hp_cns_gbl->kcond = hp_cns_gbl->R*hp_cns_gbl->mu/prandtl*hp_cns_gbl->gamma/(hp_cns_gbl->gamma-1.0);

	*gbl->log << "#conductivity: " << hp_cns_gbl->kcond << endl;
	
	/* source term for MMS */
//	keyword = gbl->idprefix + "_src";
//	std::string ibcname;
//	if (!inmap.get(keyword,ibcname)) {
//		keyword = "src";
//		if (!inmap.get(keyword,ibcname)) {
//			*gbl->log << "couldn't find cd velocity field" << std::endl;
//		}
//	}
//	gbl->src = getnewibc(ibcname);
//  gbl->src->init(inmap,keyword);
	
	return;
}

void tet_hp_cns::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	const tet_hp_cns& inmesh = dynamic_cast<const tet_hp_cns &>(in);
	gbl = inmesh.gbl;

	tet_hp::init(in,why,sizereduce1d);
	
	adis = inmesh.adis;
		
	return;
}



/* OVERRIDE VIRTUAL FUNCTION FOR COMPRESSIBLE FLOW */
void tet_hp_cns::calculate_unsteady_sources() {
    int i,j,k,n,tind;
	int stridex=MXGP*MXGP;
	int stridey=MXGP;
	FLT	ogm1 = 1.0/(hp_cns_gbl->gamma-1.0);

        
    for (log2p=0;log2p<=log2pmax;++log2p) {
        for(tind=0;tind<ntet;++tind) {
            if (tet(tind).info > -1) {
                crdtocht(tind,1);
                for(n=0;n<ND;++n)
					basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
            }
            else {
                for(n=0;n<ND;++n)
                    basis::tet(log2p).proj(vrtxbd(1)(tet(tind).pnt(0))(n),vrtxbd(1)(tet(tind).pnt(1))(n),vrtxbd(1)(tet(tind).pnt(2))(n),vrtxbd(1)(tet(tind).pnt(3))(n),&crd(n)(0)(0)(0),stridex,stridey);

                for(i=0;i<basis::tet(log2p).gpx;++i) {
                    for(j=0;j<basis::tet(log2p).gpy;++j) {
						for(k=0;k<basis::tet(log2p).gpz;++k) {
							for(n=0;n<ND;++n) {
								dcrd(n)(0)(i)(j)(k) = 0.5*(vrtxbd(1)(tet(tind).pnt(3))(n) -vrtxbd(1)(tet(tind).pnt(2))(n));
								dcrd(n)(1)(i)(j)(k) = 0.5*(vrtxbd(1)(tet(tind).pnt(1))(n) -vrtxbd(1)(tet(tind).pnt(2))(n));
								dcrd(n)(2)(i)(j)(k) = 0.5*(vrtxbd(1)(tet(tind).pnt(0))(n) -vrtxbd(1)(tet(tind).pnt(2))(n));

							}
						}
                    }
                }
            }
            
            ugtouht(tind,1);

			for(n=0;n<NV;++n)
                basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);

            for(i=0;i<basis::tet(log2p).gpx;++i) { 
                for(j=0;j<basis::tet(log2p).gpy;++j) {    
					for(k=0;k<basis::tet(log2p).gpz;++k) {    
						double rho = (u(0)(i)(j)(k)+hp_cns_gbl->atm_pressure)/u(NV-1)(i)(j)(k);

						cjcb(i)(j)(k) = -gbl->bd(0)*rho*(dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k)));
						
						dugdt(log2p,tind,0)(i)(j)(k) = cjcb(i)(j)(k);
						
						for(n=1;n<NV-1;++n)
							dugdt(log2p,tind,n)(i)(j)(k) = u(n)(i)(j)(k)*cjcb(i)(j)(k);

						double e = ogm1*u(NV-1)(i)(j)(k) +0.5*(u(1)(i)(j)(k)*u(1)(i)(j)(k) +u(2)(i)(j)(k)*u(2)(i)(j)(k)+u(3)(i)(j)(k)*u(3)(i)(j)(k));

						dugdt(log2p,tind,NV-1)(i)(j)(k) = e*cjcb(i)(j)(k);

						for(n=0;n<ND;++n)
							dxdt(log2p,tind,n)(i)(j)(k) = crd(n)(i)(j)(k);
					}
                }                
            }
        }
    }
	
    log2p=log2pmax;
    
    return;
}
