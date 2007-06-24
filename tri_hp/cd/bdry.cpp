/*
 *  cd_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/13/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "bdry_cd.h"
#include "myblas.h"

using namespace bdry_cd;

block::ctrl dirichlet::tadvance(bool coarse, block::ctrl ctrl_message) {
    int j,k,m,n,v0,v1,sind,indx,info;
    TinyVector<FLT,mesh::ND> pt;
    char uplo[] = "U";
    block::ctrl state;
    
    if (ctrl_message == block::begin) excpt1 = 0;
    
    if (excpt1 == 0) {
        if ((state = hp_side_bdry::tadvance(coarse,ctrl_message)) != block::stop) return(state); 
        ++excpt1;  

        /* UPDATE BOUNDARY CONDITION VALUES */
        for(j=0;j<base.nel;++j) {
            sind = base.el(j);
            v0 = x.sd(sind).vrtx(0);
            x.ug.v(v0,0) = x.gbl_ptr->ibc->f(0,x.vrtx(v0));
        }
        v0 = x.sd(sind).vrtx(1);
        x.ug.v(v0,0) = x.gbl_ptr->ibc->f(0,x.vrtx(v0));
        
        /*******************/    
        /* SET SIDE VALUES */
        /*******************/
        for(j=0;j<base.nel;++j) {
            sind = base.el(j);
            v0 = x.sd(sind).vrtx(0);
            v1 = x.sd(sind).vrtx(1);
            
            if (is_curved()) {
                x.crdtocht1d(sind);
                for(n=0;n<mesh::ND;++n)
                    basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
            }
            else {
                for(n=0;n<mesh::ND;++n) {
                    basis::tri(x.log2p).proj1d(x.vrtx(v0)(n),x.vrtx(v1)(n),&x.crd(n)(0,0));
                    
                    for(k=0;k<basis::tri(x.log2p).gpx;++k)
                        x.dcrd(n,0)(0,k) = 0.5*(x.vrtx(v1)(n)-x.vrtx(v0)(n));
                }
            }

            if (basis::tri(x.log2p).sm) {
                for(n=0;n<x.NV;++n)
                    basis::tri(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
        
                for(k=0;k<basis::tri(x.log2p).gpx; ++k) {
                    pt(0) = x.crd(0)(0,k);
                    pt(1) = x.crd(1)(0,k);
                    for(n=0;n<x.NV;++n)
                        x.res(n)(0,k) -= x.gbl_ptr->ibc->f(n,pt);
                }
                for(n=0;n<x.NV;++n)
                    basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
        
                indx = sind*x.sm0;
                for(n=0;n<x.NV;++n) {
                    PBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.lf(n)(2),basis::tri(x.log2p).sm,info);
                    for(m=0;m<basis::tri(x.log2p).sm;++m) 
                        x.ug.s(sind,m,n) = -x.lf(n)(2+m);
                }
            }
        }
    }
    return(block::stop);
}

block::ctrl neumann::rsdl(block::ctrl ctrl_message) {
    int j,k,n,v0,v1,sind;
    TinyVector<FLT,2> pt,mvel,nrm;

    if (ctrl_message == block::begin) {
        for(j=0;j<base.nel;++j) {
            sind = base.el(j);
            v0 = x.sd(sind).vrtx(0);
            v1 = x.sd(sind).vrtx(1);
            
            x.crdtocht1d(sind);
            for(n=0;n<mesh::ND;++n)
                basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
            
            x.crdtocht1d(sind,1);
            for(n=0;n<mesh::ND;++n)
                basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(1,0));
            
            x.ugtouht1d(sind);
            for(n=0;n<x.NV;++n)
                basis::tri(x.log2p).proj1d(&x.uht(n)(0),&x.u(n)(0,0));

            for(k=0;k<basis::tri(x.log2p).gpx;++k) {
                nrm(0) = x.dcrd(1,0)(0,k);
                nrm(1) = -x.dcrd(0,0)(0,k);                
                for(n=0;n<mesh::ND;++n) {
                    pt(n) = x.crd(n)(0,k);
                    mvel(n) = sim::bd[0]*(x.crd(n)(0,k) -x.crd(n)(1,k));
                }

                x.res(0)(0,k) = RAD(x.crd(0)(0,k))*flux(x.u(0)(0,k),pt,mvel,nrm);
            }

            for(n=0;n<x.NV;++n)
                basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
            
            for(n=0;n<x.NV;++n)
                x.gbl_ptr->res.v(v0,n) += x.lf(n)(0);

            for(n=0;n<x.NV;++n)
                x.gbl_ptr->res.v(v1,n) += x.lf(n)(1);
            
            for(k=0;k<basis::tri(x.log2p).sm;++k) {
                for(n=0;n<x.NV;++n)
                    x.gbl_ptr->res.s(sind,k,n) += x.lf(n)(k+2);
            }
        }
    }
    
    return(block::begin);
}
