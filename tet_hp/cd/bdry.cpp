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


/* ---------------------- */
/* vrtx boundary routines */
/* ---------------------- */


void dirichlet_pt::tadvance() {

	hp_vrtx_bdry::tadvance(); 
	
	int v0 = base.pnt;
	x.ug.v(v0,0) = ibc->f(0,x.pnts(v0),x.gbl->time);
	
	return;
}

void neumann_pt::rsdl(int stage){
	
	int v0 = base.pnt;
	
	for(int n=0;n<x.NV;++n)
		x.uht(n)(0) = x.ug.v(v0,n);
	
	element_rsdl(stage);
		
	for(int n=0;n<x.NV;++n)
		x.gbl->res.v(v0,n) += x.lf(n)(0);

	return;
}

void neumann_pt::element_rsdl(int stage) {

	x.lf = 0.0;
	
	return;
}


/* ---------------------- */
/* edge boundary routines */
/* ---------------------- */

void dirichlet_edge::tadvance() {	
	int j,k,m,n,v0,v1,sind;
	TinyVector<FLT,tet_mesh::ND> pt;
	
	hp_edge_bdry::tadvance(); 
	
	/* UPDATE BOUNDARY CONDITION VALUES */
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		x.ug.v(v0,0) = ibc->f(0,x.pnts(v0),x.gbl->time);
	}
	v0 = x.seg(sind).pnt(1);
	x.ug.v(v0,0) = ibc->f(0,x.pnts(v0),x.gbl->time);	
	
    /*******************/    
    /* SET SIDE VALUES */
    /*******************/
    for(j=0;j<base.nseg;++j) {
        sind = base.seg(j).gindx;
        v0 = x.seg(sind).pnt(0);
        v1 = x.seg(sind).pnt(1);
        
        if (is_curved()) {
            x.crdtocht1d(sind);
            for(n=0;n<tet_mesh::ND;++n)
                basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0));
        }
        else {
            for(n=0;n<tet_mesh::ND;++n) {
                basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd1d(n)(0));
            }
        }
		
        if (basis::tet(x.log2p).em) {
            for(n=0;n<x.NV;++n)
                basis::tet(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res1d(n)(0));
			
            for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
                pt(0) = x.crd1d(0)(k);
                pt(1) = x.crd1d(1)(k);
                pt(2) = x.crd1d(2)(k);
                for(n=0;n<x.NV;++n)
                    x.res1d(n)(k) -= x.gbl->ibc->f(n,pt,x.gbl->time);
            }
            for(n=0;n<x.NV;++n)
                basis::tet(x.log2p).intgrt1d(&x.lf(n)(0),&x.res1d(n)(0));
			
            for(n=0;n<x.NV;++n) {
                for(m=0;m<basis::tet(x.log2p).em;++m) 
                    x.ug.e(sind,m,n) = -x.lf(n)(2+m)*basis::tet(x.log2p).diag1d(m);
            }
        }
    }
	
	return;
}

void neumann_edge::rsdl(int stage){
	int sind = -1;
	
	for(int i=0;i<base.nseg;++i){
		sind = base.seg(i).gindx;
		
		x.ugtouht1d(sind);
		
		element_rsdl(sind,stage);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(x.seg(sind).pnt(0),n) += x.lf(n)(0);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(x.seg(sind).pnt(1),n) += x.lf(n)(1);
		
		int indx = 3;

		for(int k=0;k<basis::tet(x.log2p).em;++k) {
			for(int n=0;n<x.NV;++n)
				x.gbl->res.e(sind,k,n) += x.lf(n)(indx);
			++indx;
		}	
	}
	
	return;
}

void neumann_edge::element_rsdl(int eind,int stage) {
	
	x.lf = 0.0;
	
	return;
}

/* ---------------------- */
/* face boundary routines */
/* ---------------------- */

void dirichlet::tadvance() {
	int j,k,m,n,v0,v1,sind;
	TinyVector<FLT,tet_mesh::ND> pt;
	
	hp_face_bdry::tadvance(); 
	
	/* UPDATE BOUNDARY CONDITION VALUES */
	for(j=0;j<base.npnt;++j) {
		v0 = base.pnt(j).gindx;
		x.ug.v(v0,0) = ibc->f(0,x.pnts(v0),x.gbl->time);
	}
	
	
    /*******************/    
    /* SET SIDE VALUES */
    /*******************/
    for(j=0;j<base.nseg;++j) {
        sind = base.seg(j).gindx;
        v0 = x.seg(sind).pnt(0);
        v1 = x.seg(sind).pnt(1);
        
        if (is_curved()) {
            x.crdtocht1d(sind);
            for(n=0;n<tet_mesh::ND;++n)
                basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0));
        }
        else {
            for(n=0;n<tet_mesh::ND;++n) {
                basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd1d(n)(0));
                
				//                for(k=0;k<basis::tet(x.log2p).gpx;++k)
				//                    x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
            }
        }
		
        if (basis::tet(x.log2p).em) {
            for(n=0;n<x.NV;++n)
                basis::tet(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res1d(n)(0));
			
            for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
                pt(0) = x.crd1d(0)(k);
                pt(1) = x.crd1d(1)(k);
                pt(2) = x.crd1d(2)(k);
                for(n=0;n<x.NV;++n)
                    x.res1d(n)(k) -= x.gbl->ibc->f(n,pt,x.gbl->time);
            }
            for(n=0;n<x.NV;++n)
                basis::tet(x.log2p).intgrt1d(&x.lf(n)(0),&x.res1d(n)(0));
			
            for(n=0;n<x.NV;++n) {
                for(m=0;m<basis::tet(x.log2p).em;++m) 
                    x.ug.e(sind,m,n) = -x.lf(n)(2+m)*basis::tet(x.log2p).diag1d(m);
            }
        }
    }
	
	return;
}

void neumann::rsdl(int stage){
	int sind;
	FLT sgn,msgn;
	
	for(int i=0;i<base.ntri;++i){
		
		int find = base.tri(i).gindx;
		x.ugtouht2d(find);
		
		element_rsdl(find,stage);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(x.tri(find).pnt(0),n) += x.lf(n)(0);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(x.tri(find).pnt(1),n) += x.lf(n)(1);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(x.tri(find).pnt(2),n) += x.lf(n)(2);
		
		int indx = 3;
		for(int j=0;j<3;++j) {
			sind=x.tri(find).seg(j);
			sgn = x.tri(find).sgn(j);
			msgn = 1.0;
			for(int k=0;k<basis::tet(x.log2p).em;++k) {
				for(int n=0;n<x.NV;++n)
					x.gbl->res.e(sind,k,n) += msgn*x.lf(n)(indx);
				msgn *= sgn;
				++indx;
			}
		}
		
	    for(int k=0;k<basis::tet(x.log2p).fm;++k) {
		    for(int n=0;n<x.NV;++n)
				x.gbl->res.f(find,k,n) += x.lf(n)(indx);
			++indx;
		}		
	}
	
	return;
}

void neumann::element_rsdl(int find,int stage) {
	
	x.lf = 0.0;
	
	return;
}
