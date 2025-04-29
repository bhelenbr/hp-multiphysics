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
    const int log2p = x.log2p;
    
    /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
    x.crdtocht(tind,tlvl);
    
    /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(log2p)->proj_bdry(&x.cht(n,0), &crd(n)(0,0), MXGP);
}

void tri_hp::metric::calc_positions1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    x.crdtocht1d(sind,tlvl);
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n)(0));
}

void tri_hp::metric::calc_positions_leg(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    x.crdtocht(tind,tlvl);
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(x.log2p)->proj_bdry_leg(&x.cht(n,0),&crd(n)(0,0),MXGP);
}

void tri_hp::metric::calc_positions_leg1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
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
    map = getnewmapping(mapval);
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
    
    metric::calc_metrics1D(sind, crd, dcrd, tlvl);
    
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
        }
    }
}

void mapped_metric::calc_positions1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;
    const int lgpx = basis::tri(log2p)->gpx();
    
    metric::calc_positions1D(sind, crd, tlvl);

    for(int i=0;i<lgpx;++i) {
        const TinyVector<FLT,tri_mesh::ND> pt(crd(0)(i),crd(1)(i));
        
        TinyVector<FLT,tri_mesh::ND> xpt;
        map->to_physical_frame(pt, xpt);
        crd(0)(i) = xpt(0);
        crd(1)(i) = xpt(1);
    }
}



void mapped_metric::calc_positions_leg(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;
    
    metric::calc_positions_leg(tind, crd, tlvl);
    
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

void mapped_metric::calc_positions_leg1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;
    
    metric::calc_positions_leg1D(sind, crd, tlvl);
    
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


void mapped_edge_metric::calc_metrics(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND,tri_mesh::ND>& dcrd, int tlvl) const {
    const int log2p = x.log2p;
    
    /* LOAD INDICES OF VERTEX POINTS */
    TinyVector<int,3> v = x.tri(tind).pnt;
    
    /* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
    /* Linear part only*/
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(log2p)->proj(x.pnts(v(0))(n),x.pnts(v(1))(n),x.pnts(v(2))(n),&crd(n)(0,0),MXGP);
    
    /* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
    for(int n=0;n<tri_mesh::ND;++n) {
        for(int i=1;i<basis::tri(log2p)->sm();++i) {
            for(int j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
                dcrd(n,0)(i,j) = 0.5*(x.pnts(v(2))(n) -x.pnts(v(1))(n));
                dcrd(n,1)(i,j) = 0.5*(x.pnts(v(0))(n) -x.pnts(v(1))(n));
            }
        }
    }
    
    for (int s=0;s<3;++s) {
        const int sind = x.tri(tind).seg(s);
        if (x.seg(sind).tri(1) < 0) {
            const int bnum = x.getbdrynum(x.seg(sind).tri(1));
            const int indx = x.getbdryseg(x.seg(sind).tri(1));
            x.hp_ebdry(bnum)->calc_metrics(indx, s, crd, dcrd, tlvl);
        }
    }
}

void mapped_edge_metric::calc_metrics1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& dcrd, int tlvl) const {
    /* I don't think this function is necessary */
    if (x.seg(sind).tri(1) < 0) {
        const int bnum = x.getbdrynum(x.seg(sind).tri(1));
        const int indx = x.getbdryseg(x.seg(sind).tri(1));
        x.hp_ebdry(bnum)->calc_metrics1D(indx, crd, dcrd, tlvl);
    }
    else {
        /* This is a linear side */
        x.crdtocht1d(sind,tlvl);
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n)(0),&dcrd(n)(0));
    }
}

void mapped_edge_metric::calc_positions(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;
    
    /* LOAD INDICES OF VERTEX POINTS */
    TinyVector<int,3> v = x.tri(tind).pnt;
    
    /* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
    /* Linear part only*/
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(log2p)->proj(x.pnts(v(0))(n),x.pnts(v(1))(n),x.pnts(v(2))(n),&crd(n)(0,0),MXGP);
    
    for (int s=0;s<3;++s) {
        const int sind = x.tri(tind).seg(s);
        if (x.seg(sind).tri(1) < 0) {
            const int bnum = x.getbdrynum(x.seg(sind).tri(1));
            const int indx = x.getbdryseg(x.seg(sind).tri(1));
            x.hp_ebdry(bnum)->calc_positions(indx, s, crd, tlvl);
        }
    }
}

void mapped_edge_metric::calc_positions1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    if (x.seg(sind).tri(1) < 0) {
        const int bnum = x.getbdrynum(x.seg(sind).tri(1));
        const int indx = x.getbdryseg(x.seg(sind).tri(1));
        x.hp_ebdry(bnum)->calc_positions1D(indx, crd, tlvl);
    }
    else {
        /* This is a linear side */
        x.crdtocht1d(sind,tlvl);
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n)(0));
    }
}


void mapped_edge_metric::calc_positions_leg(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;
    
    /* LOAD INDICES OF VERTEX POINTS */
    TinyVector<int,3> v = x.tri(tind).pnt;
    
    /* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
    /* Linear part only*/
    for(int n=0;n<tri_mesh::ND;++n)
        basis::tri(log2p)->proj_leg(x.pnts(v(0))(n),x.pnts(v(1))(n),x.pnts(v(2))(n),&crd(n)(0,0),MXGP);
    
    for (int s=0; s<3;++s) {
        const int sind = x.tri(tind).seg(s);
        if (x.seg(sind).tri(1) < 0) {
            const int bnum = x.getbdrynum(x.seg(sind).tri(1));
            const int indx = x.getbdryseg(x.seg(sind).tri(1));
            x.hp_ebdry(bnum)->calc_positions_leg(indx, s, crd, tlvl);
        }
    }
}

void mapped_edge_metric::calc_positions_leg1D(int sind, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    if (x.seg(sind).tri(1) < 0) {
        const int bnum = x.getbdrynum(x.seg(sind).tri(1));
        const int indx = x.getbdryseg(x.seg(sind).tri(1));
        x.hp_ebdry(bnum)->calc_positions_leg1D(indx, crd, tlvl);
    }
    else {
        /* This is a linear side */
        x.crdtocht1d(sind,tlvl);
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n)(0));
    }
}

void mapped_edge_metric::setinfo() {
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
    
    for(int i=0;i<x.nebd;++i) {
        if (x.hp_ebdry(i)->is_curved()  || x.hp_ebdry(i)->mapped) {
            for(int j=0;j<x.ebdry(i)->nseg;++j) {
                int sind = x.ebdry(i)->seg(j);
                x.seg(sind).info = 0;
                x.tri(x.seg(sind).tri(0)).info = 0;
            }
        }
    }
    
    return;
}


/* These are use to calculate mappings on elements adjacent to the boundary */
void hp_edge_bdry::calc_metrics(int indx, int sd, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND,tri_mesh::ND>& dcrd, int tlvl) const {
    TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND> add_to_crd;
    TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND,tri_mesh::ND> add_to_dcrd;
    const int log2p = x.log2p;
    
    if (curved) {
        TinyMatrix<FLT,tri_mesh::ND,MXTM> cht = 0.0;
        const int sm = basis::tri(x.log2p)->sm();
        
        int cnt = 3 +(sd-1)*sm;
        /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS THIS SIDE ONLY */
        for (int m = 0; m < sm; ++m) {
            for(int n = 0; n < tri_mesh::ND; ++n)
                cht(n,cnt) = crvbd(tlvl)(indx,m)(n);
            ++cnt;
        }
        
        /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(log2p)->proj_bdry(&cht(n,0), &add_to_crd(n)(0,0), &add_to_dcrd(n,0)(0,0), &add_to_dcrd(n,1)(0,0),MXGP);
    }
    if (mapped) {
        /* Determine s location of vertex points */
        /* Reverse map back to parametric coordinates */
        const int sind = base.seg(indx);
        const int v0 = x.seg(sind).pnt(0);
        const int v1 = x.seg(sind).pnt(1);
        TinyVector<FLT,tri_mesh::ND> pt1, pt2, pt, xcrv, xlin;
        map->to_parametric_frame(x.pnts(v0),pt1);
        map->to_parametric_frame(x.pnts(v1),pt2);
        if (abs(pt1(1)-pt2(1)) > 1.0e-7) {
            *x.gbl->log << "#Uh-oh difference in normal positions " << base.idprefix << ' ' << pt1(1) << ' ' << pt2(1) << std::endl;
        }
        
        const int gpx = basis::tri(log2p)->gpx(), gpn = basis::tri(log2p)->gpn();
        switch(sd) {
            case(0): {
                for(int i=0;i<gpx;++i) {
                    pt = pt1*basis::tri(log2p)->gx(i,1) +pt2*basis::tri(log2p)->gx(i,2);
                    xlin = x.pnts(v0)*basis::tri(log2p)->gx(i,1) +x.pnts(v1)*basis::tri(log2p)->gx(i,2);
                    map->to_physical_frame(pt,xcrv);
                    xcrv -= xlin;
                    
                    TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
                    map->calc_metrics(pt, dxdtn);
                    
                    for(int j=0;j<gpn;++j) {
                        for(int n=0;n<tri_mesh::ND;++n) {
                            add_to_crd(n)(i,j) = xcrv(n)*basis::tri(log2p)->gn(j,3);
                            add_to_dcrd(n,0)(i,j) = basis::tri(log2p)->n0(j)*0.5*(dxdtn(n,0)*(pt2(0)-pt1(0)) -(x.pnts(v1)(n)-x.pnts(v0)(n)))*basis::tri(log2p)->gn(j,3);
                            add_to_dcrd(n,1)(i,j) = xcrv(n)*basis::tri(log2p)->dgn(j,3) +basis::tri(log2p)->x0(i)*add_to_dcrd(n,0)(i,j);
                        }
                    }
                }
                break;
            }
            case(1): {
                for(int j=0;j<gpn;++j) {
                    pt = pt1*basis::tri(log2p)->gn(j,1) +pt2*basis::tri(log2p)->gn(j,0);
                    xlin = x.pnts(v0)*basis::tri(log2p)->gn(j,1) +x.pnts(v1)*basis::tri(log2p)->gn(j,0);
                    map->to_physical_frame(pt,xcrv);
                    xcrv -= xlin;
                    
                    TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
                    map->calc_metrics(pt, dxdtn);
                    
                    for(int i=0;i<gpx;++i) {
                        for(int n=0;n<tri_mesh::ND;++n) {
                            add_to_crd(n)(i,j) = xcrv(n)*basis::tri(log2p)->gx(i,2);
                            add_to_dcrd(n,0)(i,j) = 0.5*xcrv(n)*basis::tri(log2p)->n0(j);
                            add_to_dcrd(n,1)(i,j) = (0.5*(dxdtn(n,0)*(pt2(0)-pt1(0)) -(x.pnts(v1)(n)-x.pnts(v0)(n)))*basis::tri(log2p)->gx(i,2) +add_to_dcrd(n,0)(i,j))*basis::tri(log2p)->x0(i);
                        }
                    }
                }
                break;
            }
            case(2): {
                for(int j=0;j<gpn;++j) {
                    pt = pt1*basis::tri(log2p)->gn(j,1) +pt2*basis::tri(log2p)->gn(j,0);
                    xlin = x.pnts(v0)*basis::tri(log2p)->gn(j,1) +x.pnts(v1)*basis::tri(log2p)->gn(j,0);
                    map->to_physical_frame(pt,xcrv);
                    xcrv -= xlin;
                    
                    TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
                    map->calc_metrics(pt, dxdtn);
                    
                    for(int i=0;i<gpx;++i) {
                        for(int n=0;n<tri_mesh::ND;++n) {
                            add_to_crd(n)(i,j) = xcrv(n)*basis::tri(log2p)->gx(i,1);
                            add_to_dcrd(n,0)(i,j) = -0.5*xcrv(n)*basis::tri(log2p)->n0(j);
                            add_to_dcrd(n,1)(i,j) = (0.5*(dxdtn(n,0)*(pt2(0)-pt1(0)) -(x.pnts(v1)(n)-x.pnts(v0)(n)))*basis::tri(log2p)->gx(i,1) +add_to_dcrd(n,0)(i,j))*basis::tri(log2p)->x0(i);
                        }
                    }
                }
                break;
            }
            default:
                sim::abort(__LINE__,__FILE__,x.gbl->log);
        }
        
        for(int i=0;i<gpx;++i) {
            for(int j=0;j<gpn;++j) {
                for(int n=0;n<tri_mesh::ND;++n) {
                    dcrd(n,0)(i,j) += add_to_dcrd(n,0)(i,j);
                    dcrd(n,1)(i,j) += add_to_dcrd(n,1)(i,j);
                    crd(n)(i,j) += add_to_crd(n)(i,j);
                }
            }
        }
    }
}

void hp_edge_bdry::calc_positions(int indx, int sd, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND> add_to_crd;
    const int log2p = x.log2p;
    
    if (curved) {
        TinyMatrix<FLT,tri_mesh::ND,MXTM> cht = 0.0;
        const int sm = basis::tri(x.log2p)->sm();
        
        int cnt = 3 +(sd-1)*sm;
        /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS THIS SIDE ONLY */
        for (int m = 0; m < sm; ++m) {
            for(int n = 0; n < tri_mesh::ND; ++n)
                cht(n,cnt) = crvbd(tlvl)(indx,m)(n);
            ++cnt;
        }
        
        /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(log2p)->proj_bdry(&cht(n,0), &add_to_crd(n)(0,0), MXGP);
    }
    if (mapped) {
        /* Determine s location of vertex points */
        /* Reverse map back to parametric coordinates */
        const int sind = base.seg(indx);
        const int v0 = x.seg(sind).pnt(0);
        const int v1 = x.seg(sind).pnt(1);
        TinyVector<FLT,tri_mesh::ND> pt1, pt2, pt, xcrv, xlin;
        map->to_parametric_frame(x.pnts(v0),pt1);
        map->to_parametric_frame(x.pnts(v1),pt2);
        if (abs(pt1(1)-pt2(1)) > 1.0e-7) {
            *x.gbl->log << "#Uh-oh difference in normal positions " << base.idprefix << ' ' << pt1(1) << ' ' << pt2(1) << std::endl;
        }
        
        const int gpx = basis::tri(log2p)->gpx(), gpn = basis::tri(log2p)->gpn();
        switch(sd) {
            case(0): {
                for(int i=0;i<gpx;++i) {
                    pt = pt1*basis::tri(log2p)->gx(i,1) +pt2*basis::tri(log2p)->gx(i,2);
                    xlin = x.pnts(v0)*basis::tri(log2p)->gx(i,1) +x.pnts(v1)*basis::tri(log2p)->gx(i,2);
                    map->to_physical_frame(pt,xcrv);
                    xcrv -= xlin;
                    
                    TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
                    map->calc_metrics(pt, dxdtn);
                    
                    for(int j=0;j<gpn;++j) {
                        for(int n=0;n<tri_mesh::ND;++n) {
                            add_to_crd(n)(i,j) = xcrv(n)*basis::tri(log2p)->gn(j,3);
                        }
                    }
                }
                break;
            }
            case(1): {
                for(int j=0;j<gpn;++j) {
                    pt = pt1*basis::tri(log2p)->gn(j,1) +pt2*basis::tri(log2p)->gn(j,0);
                    xlin = x.pnts(v0)*basis::tri(log2p)->gn(j,1) +x.pnts(v1)*basis::tri(log2p)->gn(j,0);
                    map->to_physical_frame(pt,xcrv);
                    xcrv -= xlin;
                    
                    TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
                    map->calc_metrics(pt, dxdtn);
                    
                    for(int i=0;i<gpx;++i) {
                        for(int n=0;n<tri_mesh::ND;++n) {
                            add_to_crd(n)(i,j) = xcrv(n)*basis::tri(log2p)->gx(i,2);
                        }
                    }
                }
                break;
            }
            case(2): {
                for(int j=0;j<gpn;++j) {
                    pt = pt1*basis::tri(log2p)->gn(j,1) +pt2*basis::tri(log2p)->gn(j,0);
                    xlin = x.pnts(v0)*basis::tri(log2p)->gn(j,1) +x.pnts(v1)*basis::tri(log2p)->gn(j,0);
                    map->to_physical_frame(pt,xcrv);
                    xcrv -= xlin;
                    
                    TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
                    map->calc_metrics(pt, dxdtn);
                    
                    for(int i=0;i<gpx;++i) {
                        for(int n=0;n<tri_mesh::ND;++n) {
                            add_to_crd(n)(i,j) = xcrv(n)*basis::tri(log2p)->gx(i,1);
                        }
                    }
                }
                break;
            }
            default:
                sim::abort(__LINE__,__FILE__,x.gbl->log);
        }
        
        for(int i=0;i<gpx;++i) {
            for(int j=0;j<gpn;++j) {
                for(int n=0;n<tri_mesh::ND;++n) {
                    crd(n)(i,j) += add_to_crd(n)(i,j);
                }
            }
        }
    }
}
    
void hp_edge_bdry::calc_metrics1D(int indx, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& dcrd, int tlvl) const {
    const int log2p = x.log2p;
    const int gpx = basis::tri(log2p)->gpx();
    const int sind = base.seg(indx);
    
    if (mapped) {
        const int v0 = x.seg(sind).pnt(0);
        const int v1 = x.seg(sind).pnt(1);
        TinyVector<FLT,tri_mesh::ND> pt1, pt2, pt, xcrv, xlin;
        map->to_parametric_frame(x.pnts(v0),pt1);
        map->to_parametric_frame(x.pnts(v1),pt2);
        
        for(int i=0;i<gpx;++i) {
            pt(0) = pt1(0)*basis::tri(log2p)->gx(i,1) +pt2(0)*basis::tri(log2p)->gx(i,2);
            pt(1) = pt1(1)*basis::tri(log2p)->gx(i,1) +pt2(1)*basis::tri(log2p)->gx(i,2);  /* FIXME SIDES MUST BE ALIGNED WITH TANGENT COORDINATE */
            xlin = x.pnts(v0)*basis::tri(log2p)->gx(i,1) +x.pnts(v1)*basis::tri(log2p)->gx(i,2);
            map->to_physical_frame(pt,xcrv);
            xcrv -= xlin;
            
            TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
            map->calc_metrics(pt, dxdtn);
            
            for(int n=0;n<tri_mesh::ND;++n) {
                crd(n)(i) = xcrv(n);
                dcrd(n)(i) = 0.5*dxdtn(n,0)*(pt2(0)-pt1(0));
            }
        }
    }
    else {
        x.crdtocht1d(sind,tlvl);
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(log2p)->proj1d(&x.cht(n,0),&crd(n)(0),&dcrd(n)(0));
    }
}

void hp_edge_bdry::calc_positions1D(int indx, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;
    const int gpx = basis::tri(log2p)->gpx();
    const int sind = base.seg(indx);
    
    if (mapped) {
        const int v0 = x.seg(sind).pnt(0);
        const int v1 = x.seg(sind).pnt(1);
        TinyVector<FLT,tri_mesh::ND> pt1, pt2, pt, xcrv, xlin;
        map->to_parametric_frame(x.pnts(v0),pt1);
        map->to_parametric_frame(x.pnts(v1),pt2);
        
        for(int i=0;i<gpx;++i) {
            pt(0) = pt1(0)*basis::tri(log2p)->gx(i,1) +pt2(0)*basis::tri(log2p)->gx(i,2);
            pt(1) = pt1(1)*basis::tri(log2p)->gx(i,1) +pt2(1)*basis::tri(log2p)->gx(i,2);  /* FIXME SIDES MUST BE ALIGNED WITH TANGENT COORDINATE */
            xlin = x.pnts(v0)*basis::tri(log2p)->gx(i,1) +x.pnts(v1)*basis::tri(log2p)->gx(i,2);
            map->to_physical_frame(pt,xcrv);
            xcrv -= xlin;
            
            TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> dxdtn;
            map->calc_metrics(pt, dxdtn);
            
            for(int n=0;n<tri_mesh::ND;++n) {
                crd(n)(i) = xcrv(n);
            }
        }
    }
    else {
        x.crdtocht1d(sind,tlvl);
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(log2p)->proj1d(&x.cht(n,0),&crd(n)(0));
    }
}

void hp_edge_bdry::calc_positions_leg(int indx, int sd, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_mesh::ND> add_to_crd;
    const int log2p = x.log2p;
    const int sm = basis::tri(log2p)->sm();
    const int sind = base.seg(indx);

    if (mapped) {
        /* Determine s location of vertex points */
        /* Reverse map back to parametric coordinates */
        const int v0 = x.seg(sind).pnt(0);
        const int v1 = x.seg(sind).pnt(1);
        TinyVector<FLT,tri_mesh::ND> pt1, pt2, pt, xcrv, xlin;
        map->to_parametric_frame(x.pnts(v0),pt1);
        map->to_parametric_frame(x.pnts(v1),pt2);
        pt(1) = 0.5*(pt1(1)+pt2(1));
        if (abs(pt1(1)-pt2(1)) > 1.0e-7) {
            *x.gbl->log << "#Uh-oh difference in normal positions " << base.idprefix << ' ' << pt1(1) << ' ' << pt2(1) << std::endl;
        }
        
        /* Calculate xi and eta locations */
        TinyMatrix<FLT,MXGP,MXGP> xi, eta;
        for(int i=1;i<sm;++i) {
            for(int j=1;j<sm-(i-1);++j) {
                eta(i,j) = 2.*basis::tri(log2p)->lgrnge(0,i,j)-1.;
                xi(i,j) = 1. -2.*basis::tri(log2p)->lgrnge(0,i,j)/basis::tri(log2p)->lgrnge(0,i,j);
            }
        }

        const int gpx = basis::tri(log2p)->gpx(), gpn = basis::tri(log2p)->gpn();
        switch(sd) {
            case(0): {
                /* INTERIOR */
                for(int i=1;i<sm;++i) {
                    for(int j=1;j<sm-(i-1);++j) {
                        pt = pt1*(1-xi(i,j))/2. +pt2*(1+xi(i,j))/2.;
                        xlin = x.pnts(v0)*(1-xi(i,j))/2. +x.pnts(v1)*(1+xi(i,j))/2.;
                        map->to_physical_frame(pt,xcrv);
                        xcrv -= xlin;
                        for(int n=0;n<tri_mesh::ND;++n) {
                            crd(n)(i,j) += xcrv(n)*pow((1.-eta(i,j))/2.,2.0);
                        }
                    }
                }
                break;
            }
            case(1): {
                /* INTERIOR */
                for(int i=1;i<sm;++i) {
                    for(int j=1;j<sm-(i-1);++j) {
                        pt = pt1*(1-eta(i,j))/2. +pt2*(1+eta(i,j))/2.;
                        xlin = x.pnts(v0)*(1-eta(i,j))/2. +x.pnts(v1)*(1+xi(i,j))/2.;
                        map->to_physical_frame(pt,xcrv);
                        xcrv -= xlin;
                        
                        for(int n=0;n<tri_mesh::ND;++n) {
                            add_to_crd(n)(i,j) = xcrv(n)*(1+xi(i,j))/2.;
                        }
                    }
                }
                break;
            }
            case(2): {
                for(int i=1;i<sm;++i) {
                    for(int j=1;j<sm-(i-1);++j) {
                        pt = pt1*(1-eta(i,j))/2. +pt2*(1+eta(i,j))/2.;
                        xlin = x.pnts(v0)*(1-eta(i,j))/2. +x.pnts(v1)*(1+xi(i,j))/2.;
                        map->to_physical_frame(pt,xcrv);
                        xcrv -= xlin;
                        
                        for(int n=0;n<tri_mesh::ND;++n) {
                            add_to_crd(n)(i,j) = xcrv(n)*(1-xi(i,j))/2.;
                        }
                    }
                }
                break;
            }
            default:
                sim::abort(__LINE__,__FILE__,x.gbl->log);
        }
        
        for(int i=1;i<sm;++i) {
            for(int j=1;j<sm-(i-1);++j) {
                for(int n=0;n<tri_mesh::ND;++n) {
                    crd(n)(i,j) += add_to_crd(n)(i,j);
                }
            }
        }
    }
    else {
        for(int n = 0; n < tri_mesh::ND; ++n) {
            x.cht(n) = 0.0;
        }
        int ind = 3+sd*basis::tri(log2p)->sm();
        for(int m = 0; m < sm; ++m) {
            for(int n = 0; n < tri_mesh::ND; ++n) {
                x.cht(n,ind) = crv(indx,m)(n);
            }
            ++ind;
        }
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(x.log2p)->proj_bdry_leg(&x.cht(n,0),&crd(n)(0,0),MXGP);
    }
}

void hp_edge_bdry::calc_positions_leg1D(int indx, TinyVector<TinyVector<FLT,MXGP>,tri_mesh::ND>& crd, int tlvl) const {
    const int log2p = x.log2p;
    const int sm = basis::tri(log2p)->sm();
    const int sind = base.seg(indx);
    
    
    if (mapped) {
        const int v0 = x.seg(sind).pnt(0);
        const int v1 = x.seg(sind).pnt(1);
        TinyVector<FLT,tri_mesh::ND> pt1, pt2, pt, xcrv, xlin;
        map->to_parametric_frame(x.pnts(v0),pt1);
        map->to_parametric_frame(x.pnts(v1),pt2);
        
        for (int i=1;i<sm+1;++i) {
            pt = pt1*basis::tri(log2p)->lgrnge1d(0,i) +pt2*basis::tri(log2p)->lgrnge1d(1,i);
            map->to_physical_frame(pt,xcrv);
            
            for(int n=0;n<tri_mesh::ND;++n) {
                crd(n)(i) = xcrv(n);
            }
        }
    }
    else {
        x.crdtocht1d(sind,tlvl);
        for(int n=0;n<tri_mesh::ND;++n)
            basis::tri(log2p)->proj1d_leg(&x.cht(n,0),&crd(n)(0));
    }
}
    
