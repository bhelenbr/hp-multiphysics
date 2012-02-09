/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cns.h"
#include "../hp_boundary.h"


void tet_hp_cns::init(input_map& input, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	gbl = static_cast<global *>(gin);
	
	if (input.find(gbl->idprefix + "_nvariable") == input.end()) {
		input[gbl->idprefix + "_nvariable"] = "5";
	}
	
	tet_hp::init(input,gin);
	
	input.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);
		
	/* no preconditioner = 0, weiss-smith preconditioner = 1, squared preconditioner = 2 */
	input.getwdefault("preconditioner",gbl->preconditioner,1);

	gbl->tau.resize(maxvst,NV,NV);

	gbl->vpreconditioner.resize(maxvst,NV,NV);
	gbl->epreconditioner.resize(maxvst,NV,NV);
	gbl->betasquared.resize(ntet);

	double prandtl;
	
	double bodydflt[3] = {0.0,0.0,0.0};
	if (!input.get(gbl->idprefix +"_body_force",gbl->body.data(),3)) input.getwdefault("body_force",gbl->body.data(),3,bodydflt); 
	
	if (!input.get(gbl->idprefix + "_atm_pressure",gbl->atm_pressure)) input.getwdefault("atm_pressure",gbl->atm_pressure,0.0);	
	if (!input.get(gbl->idprefix + "_gamma",gbl->gamma)) input.getwdefault("gamma",gbl->gamma,1.4);
	if (!input.get(gbl->idprefix + "_mu",gbl->mu)) input.getwdefault("mu",gbl->mu,1.0);
	if (!input.get(gbl->idprefix + "_prandtl",prandtl)) input.getwdefault("prandtl",prandtl,0.713);
	if (!input.get(gbl->idprefix + "_R",gbl->R)) input.getwdefault("R",gbl->R,287.058);
	
	gbl->kcond = gbl->R*gbl->mu/prandtl*gbl->gamma/(gbl->gamma-1.0);

	/* source term for MMS */
	//gbl->src = getnewibc("src",input);
		
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
	FLT	ogm1 = 1.0/(gbl->gamma-1.0);

        
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
						double rho = (u(0)(i)(j)(k)+gbl->atm_pressure)/u(NV-1)(i)(j)(k);

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
