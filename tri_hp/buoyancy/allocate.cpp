/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_buoyancy.h"
#include "../hp_boundary.h"

 void tri_hp_buoyancy::init(input_map& input, gbl *gin) {
    std::string keyword;
    bool adapt_storage;
    bool coarse;
    
    keyword = idprefix + "_nvariable";
    input[keyword] = "4";
    
    tri_hp_ins::init(input,gin);
    
    /* Load pointer to block stuff */
    gbl_ptr = gin;
    
    keyword = idprefix + "_adapt_storage";
    input.getwdefault(keyword,adapt_storage,false);
    if (adapt_storage) return;
    
    keyword = idprefix + "_coarse";
    input.getwdefault(keyword,coarse,false);
    if (coarse) return;
    
    if (!input.get(idprefix + "_conductivity",gbl_ptr->kcond)) input.getwdefault("conductivity",gbl_ptr->kcond,0.7*gbl_ptr->mu);
    gbl_ptr->D(2) = gbl_ptr->kcond;
    if (!input.get(idprefix + "_cp",gbl_ptr->cp)) input.getwdefault("cp",gbl_ptr->cp,1.0);
    
    if (input.find(idprefix+"_rhovsT_expression") != input.end()) {
        gbl_ptr->rhovsT.init(input,idprefix+"_rhovsT");
    } 
    else if (input.find("rhovsT_expression") != input.end()){
        gbl_ptr->rhovsT.init(input,"rhovsT");
    }
    else {
        *sim::log << "couldn't find rhovsT equation for density\n";
        exit(1);
    }


    return;
}

/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tri_hp_buoyancy::calculate_unsteady_sources(bool coarse) {
    int i,j,n,tind;
    FLT lrho;
        
    for (log2p=0;log2p<=log2pmax;++log2p) {
        for(tind=0;tind<ntri;++tind) {
            if (td(tind).info > -1) {
                crdtocht(tind,1);
                for(n=0;n<ND;++n)
                    basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
            }
            else {
                for(n=0;n<ND;++n)
                    basis::tri(log2p).proj(vrtxbd(1)(td(tind).vrtx(0))(n),vrtxbd(1)(td(tind).vrtx(1))(n),vrtxbd(1)(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);

                for(i=0;i<basis::tri(log2p).gpx;++i) {
                    for(j=0;j<basis::tri(log2p).gpn;++j) {
                        for(n=0;n<ND;++n) {
                            dcrd(n,0)(i,j) = 0.5*(vrtxbd(1)(td(tind).vrtx(1))(n) -vrtxbd(1)(td(tind).vrtx(0))(n));
                            dcrd(n,1)(i,j) = 0.5*(vrtxbd(1)(td(tind).vrtx(2))(n) -vrtxbd(1)(td(tind).vrtx(0))(n));
                        }
                    }
                }
            }
            
            ugtouht(tind,1);
            for(n=0;n<NV;++n)
                basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);

            for(i=0;i<basis::tri(log2p).gpx;++i) {
                for(j=0;j<basis::tri(log2p).gpn;++j) {    
                    cjcb(i,j) = -sim::bd[0]*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
                    lrho = gbl_ptr->rhovsT.Eval(u(2)(i,j));                    
                    for(n=0;n<NV-2;++n)
                        dugdt(log2p,tind,n)(i,j) = lrho*u(n)(i,j)*cjcb(i,j);
                    dugdt(log2p,tind,NV-2)(i,j) = lrho*gbl_ptr->cp*u(NV-2)(i,j)*cjcb(i,j);
                    dugdt(log2p,tind,NV-1)(i,j) = lrho*cjcb(i,j);

                    for(n=0;n<ND;++n)
                        dxdt(log2p,tind,n)(i,j) = crd(n)(i,j);
                }                
            }
        }
    }
    log2p=log2pmax;
    
    return;
}
