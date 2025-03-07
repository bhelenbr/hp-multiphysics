//
//  metric.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 3/2/25.
//

#include <stdio.h>
#include "metric.h"
#include "hp_boundary.h"

void tri_hp::metric::calc_metrics(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND,tri_mesh::ND>& dcrd, int tlvl) const {
    const int log2p = x.log2p;

    /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
    x.crdtocht(tind,tlvl);
    
    /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(log2p)->proj_bdry(&x.cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
}

void tri_hp::metric::calc_metrics1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& dcrd, int tlvl) const {
    x.crdtocht1d(sind,tlvl);
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n)(0),&dcrd(n)(0));
}

void tri_hp::metric::calc_positions(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    x.crdtocht(tind,tlvl);
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(x.log2p)->proj_bdry_leg(&x.cht(n,0),&crd(n)(0,0),MXGP);
}

void tri_hp::metric::calc_positions1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    x.crdtocht1d(sind,tlvl);

    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(x.log2p)->proj1d_leg(&x.cht(n,0),&crd(n)(0));
}




void tri_hp::metric::setinfo() {
    /* SET UP pnts BC INFORMATION FOR OUTPUT */
    for(int i=0;i<x.npnt;++i)
        x.pnt(i).info = -1;

    for(int i=0;i<x.nvbd;++i)
        x.pnt(x.vbdry(i)->pnt).info = 0;

    /* SET UP EDGE BC INFORMATION FOR CURVED SIDES OUTPUT */
    for(int i=0;i<x.nseg;++i)
        x.seg(i).info = -1;

    for(int i=0;i<x.ntri;++i)
        x.tri(i).info = -1;

    if (x.log2p > 0) {
        for(int i=0;i<x.nebd;++i) {
            if (x.hp_ebdry(i)->is_curved()) {
                for(int j=0;j<x.ebdry(i)->nseg;++j) {
                    int sind = x.ebdry(i)->seg(j);
                    x.seg(sind).info = 0;
                    x.tri(x.seg(sind).tri(0)).info = 0;
                }
            }
        }
    }

    return;
}



void mapped_metric::init(input_map& input) {
    std::string mapval;
    if (!input.get(x.gbl->idprefix+"_mapping",mapval)) {
        *x.gbl->log << "Couldn't read mapping " << x.gbl->idprefix +"_mapping" << std::endl;
        sim::abort(__LINE__,__FILE__,x.gbl->log);
    }
    
    if (mapval == "polar") {
        map = make_shared<polar_mapping>();
    }
    else if (mapval == "polar_log") {
        map = make_shared<polar_log_mapping>();
    }
    else {
        *x.gbl->log << "Unrecognized mapping " << mapval << std::endl;
    }
    
    map->init(input,x.gbl->idprefix,x.gbl->log);
}

void mapped_metric::calc_metrics(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND,tri_mesh::ND>& dcrd, int tlvl) const {
    const int log2p = x.log2p;
    const int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
    
    /* LOAD INDICES OF VERTEX POINTS */
    TinyVector<int,3> v;
    v = x.tri(tind).pnt;
    
    /* PROJECT VERTEX COORDINATES TO GAUSS POINTS */
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(log2p)->proj(x.vrtxbd(tlvl)(v(0))(n),x.vrtxbd(tlvl)(v(1))(n),x.vrtxbd(tlvl)(v(2))(n),&crd(n)(0,0),MXGP);
    
    /* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
    for(int i=0;i<lgpx;++i) {
        for(int j=0;j<lgpn;++j) {
            const TinyVector<FLT,tri_mesh::ND> pt(crd(0)(i,j),crd(1)(i,j));
            TinyVector<FLT,tri_mesh::ND> xpt;
            
            map->to_physical_frame(pt, xpt);
            crd(0)(i,j) = xpt(0);
            crd(1)(i,j) = xpt(1);
            
            
            TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn, dtndrs;
            map->calc_metrics(pt, dxdtn);
            
            for(int n=0;n<tri_mesh::ND;++n) {
                dtndrs(n,0) = 0.5*(x.pnts(v(2))(n) -x.pnts(v(1))(n));
                dtndrs(n,1) = 0.5*(x.pnts(v(0))(n) -x.pnts(v(1))(n));
            }
            
            // dx/drs = dx/dtn*dtn/drs
            // dx/drs = [dx/dt, dx/dn]*[dtn/dr, dtn/ds]
            for (int i1 = 0; i1 < tri_mesh::ND; ++i1 ) {
                for (int j1 = 0; j1 < tri_mesh::ND; ++j1 ) {
                    FLT sum = 0.0;
                    for (int k1 = 0; k1 < tri_mesh::ND; ++k1 ) {
                        sum += dxdtn(i1,k1)*dtndrs(k1,j1);
                    }
                    dcrd(i1,j1)(i,j) = sum;
                }
            }
        }
    }
}

void mapped_metric::calc_metrics1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& dcrd, int tlvl) const {
    const int log2p = x.log2p;
    const int lgpx = basis::tri(log2p)->gpx();
    
    x.crdtocht1d(sind,tlvl);
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n)(0),&dcrd(n)(0));
    
    for(int i=0;i<lgpx;++i) {
        const TinyVector<FLT,tri_mesh::ND> pt(crd(0)(i),crd(1)(i));
        
        TinyVector<FLT,tri_mesh::ND> xpt;
        map->to_physical_frame(pt, xpt);
        crd(0)(i) = xpt(0);
        crd(1)(i) = xpt(1);
        
        TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
        map->calc_metrics(pt, dxdtn);
        
        // dx/drs = dx/dtn*dtn/drs
        // dx/drs = [dx/dt, dx/dn]*[dtn/dr, dtn/ds]
        for (int i1 = 0; i1 < tri_mesh::ND; ++i1 ) {
            FLT sum = 0.0;
            for (int k1 = 0; k1 < tri_mesh::ND; ++k1 ) {
                sum += dxdtn(i1,k1)*dcrd(k1)(i);
            }
            dcrd(i1)(i) = sum;
        }
    }
}

void mapped_metric::calc_positions(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;
    
    x.crdtocht(tind,tlvl);
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(x.log2p)->proj_bdry_leg(&x.cht(n,0),&crd(n)(0,0),MXGP);

    for(int i=1;i<basis::tri(log2p)->sm();++i) {
        for(int j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
            const TinyVector<FLT,tri_mesh::ND> pt(crd(0)(i,j),crd(1)(i,j));
            TinyVector<FLT,tri_mesh::ND> xpt;
            map->to_physical_frame(pt, xpt);
            crd(0)(i,j) = xpt(0);
            crd(1)(i,j) = xpt(1);
        }
    }
}

void mapped_metric::calc_positions1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;

    x.crdtocht1d(sind,tlvl);

    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(x.log2p)->proj1d_leg(&x.cht(n,0),&crd(n)(0));
                        
    for(int i=1;i<basis::tri(log2p)->sm()+1;++i) {
        const TinyVector<FLT,tri_mesh::ND> pt(crd(0)(i),crd(1)(i));
        TinyVector<FLT,tri_mesh::ND> xpt;
        map->to_physical_frame(pt, xpt);
        crd(0)(i) = xpt(0);
        crd(1)(i) = xpt(1);
    }
}

void mapped_metric::calc_positions0D(int vind, TinyVector<FLT,tri_mesh::ND>& pt, int tlvl) const {
    map->to_physical_frame(x.vrtxbd(tlvl)(vind), pt);
}

void mapped_metric::setinfo() {
    /* SET UP pnts BC INFORMATION FOR OUTPUT */
    for(int i=0;i<x.npnt;++i)
        x.pnt(i).info = -1;

    for(int i=0;i<x.nvbd;++i)
        x.pnt(x.vbdry(i)->pnt).info = 0;
    
    /* SET UP EDGE BC INFORMATION FOR CURVED SIDES OUTPUT */
    for(int i=0;i<x.nseg;++i)
        x.seg(i).info = -1;
    
    if (x.log2p > 0) {
        for(int i=0;i<x.nebd;++i) {
            if (x.hp_ebdry(i)->is_curved()) {
                for(int j=0;j<x.ebdry(i)->nseg;++j) {
                    int sind = x.ebdry(i)->seg(j);
                    x.seg(sind).info = 0;
                    x.tri(x.seg(sind).tri(0)).info = 0;
                }
            }
        }
    }

    for(int i=0;i<x.ntri;++i)
        x.tri(i).info = 0;
    
    return;
}

void allcurved_metric::setinfo() {
    /* SET UP pnts BC INFORMATION FOR OUTPUT */
    for(int i=0;i<x.npnt;++i)
        x.pnt(i).info = -1;

    for(int i=0;i<x.nvbd;++i)
        x.pnt(x.vbdry(i)->pnt).info = 0;

    /* SET UP EDGE BC INFORMATION FOR CURVED SIDES OUTPUT */
    for(int i=0;i<x.nseg;++i)
        x.seg(i).info = 0;

    for(int i=0;i<x.ntri;++i)
        x.tri(i).info = 0;

    return;
}
