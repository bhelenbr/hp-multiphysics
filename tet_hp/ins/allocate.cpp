/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_ins.h"
#include "../hp_boundary.h"


void tet_hp_ins::init(input_map& inmap, shared_ptr<block_global> gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	gbl = gin;
	
	if (inmap.find(gbl->idprefix + "_nvariable") == inmap.end()) {
		inmap[gbl->idprefix + "_nvariable"] = "4";
	}
	
	tet_hp::init(inmap,gin);
	
	inmap.getwdefault(gbl->idprefix + "_dissipation",adis,1.0);
		
	hp_ins_gbl->tau.resize(maxvst,NV);

	if (!inmap.get(gbl->idprefix + "_rho",hp_ins_gbl->rho)) inmap.getwdefault("rho",hp_ins_gbl->rho,1.0);
	if (!inmap.get(gbl->idprefix + "_mu",hp_ins_gbl->mu)) inmap.getwdefault("mu",hp_ins_gbl->mu,0.0);
	
	/* LEAVE UP TO DERIVED CLASSES TO LOAD THESE IF NECESSARY */
	hp_ins_gbl->D.resize(NV);
	if (NV > 4) {
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

void tet_hp_ins::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	const tet_hp_ins& inmesh = dynamic_cast<const tet_hp_ins &>(in);
	gbl = inmesh.gbl;

	tet_hp::init(in,why,sizereduce1d);
	
	adis = inmesh.adis;
		
	return;
}



/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tet_hp_ins::calculate_unsteady_sources() {
    int i,j,k,n,tind;
	int stridex=MXGP*MXGP;
	int stridey=MXGP;
        
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
						cjcb(i)(j)(k) = -gbl->bd(0)*hp_ins_gbl->rho*(dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k)));
						for(n=0;n<NV-1;++n)
							dugdt(log2p,tind,n)(i)(j)(k) = u(n)(i)(j)(k)*cjcb(i)(j)(k);
						dugdt(log2p,tind,NV-1)(i)(j)(k) = cjcb(i)(j)(k);

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
