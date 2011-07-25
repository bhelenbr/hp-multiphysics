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

void dirichlet::tadvance() {
	int j,k,m,n,v0,v1,sind,indx,info;
	TinyVector<FLT,tet_mesh::ND> pt;
	char uplo[] = "U";
		
	hp_face_bdry::tadvance(); 
	
	/* UPDATE BOUNDARY CONDITION VALUES */
	for(j=0;j<base.npnt;++j) {
		v0 = base.pnt(j).gindx;
		x.ug.v(v0,0) = ibc->f(0,x.pnts(v0),x.gbl->time);
	}

//
//    /*******************/    
//    /* SET SIDE VALUES */
//    /*******************/
//    for(j=0;j<base.nseg;++j) {
//        sind = base.seg(j);
//        v0 = x.seg(sind).pnt(0);
//        v1 = x.seg(sind).pnt(1);
//        
//        if (is_curved()) {
//            x.crdtocht1d(sind);
//            for(n=0;n<tet_mesh::ND;++n)
//                basis::tet(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
//        }
//        else {
//            for(n=0;n<tet_mesh::ND;++n) {
//                basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
//                
//                for(k=0;k<basis::tet(x.log2p).gpx;++k)
//                    x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
//            }
//        }
//
//        if (basis::tet(x.log2p).em) {
//            for(n=0;n<x.NV;++n)
//                basis::tet(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
//
//            for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
//                pt(0) = x.crd(0)(0,k);
//                pt(1) = x.crd(1)(0,k);
//                for(n=0;n<x.NV;++n)
//                    x.res(n)(0,k) -= x.gbl->ibc->f(n,pt,x.gbl->time);
//            }
//            for(n=0;n<x.NV;++n)
//                basis::tet(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
//
//            indx = sind*x.em0;
//            for(n=0;n<x.NV;++n) {
//                // fix temp PBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.lf(n)(2),basis::tri(x.log2p).sm,info);
//                for(m=0;m<basis::tet(x.log2p).em;++m) 
//                    x.ug.e(sind,m,n) = -x.lf(n)(2+m);
//            }
//        }
//    }
	
	return;
}



void neumann::element_rsdl(int find,int stage) {
	int j,k,n;
	TinyVector<FLT,3> pt,mvel,nrm,vec1,vec2;
	FLT u,flx;
	
	x.lf = 0.0;
	
//	x.crdtocht2d(find);
//	for(n=0;n<tet_mesh::ND;++n)
//		basis::tet(x.log2p).proj2d(&x.cht(n)(0),&x.crd2d(n)(0)(0),&x.dcrd2d(n)(0)(0)(0),&x.dcrd2d(n)(1)(0)(0),MXGP);
//	
//	basis::tet(x.log2p).proj2d(&x.uht(0)(0),&x.u2d(0)(0)(0),MXGP);
//	
//	for(j=0;j<basis::tet(x.log2p).gpx;++j) {
//		for(k=0;k<basis::tet(x.log2p).gpy;++k) {
//			/* co-variant vectors: d vec(x)/dxi crossed with d vec(x)/deta */
//			for(n=0;n<3;++n){
//				vec1(n)=x.dcrd2d(n)(0)(j)(k);
//				vec2(n)=x.dcrd2d(n)(1)(j)(k);
//			}
//			nrm=cross(vec1,vec2);
//			
//			for(n=0;n<tet_mesh::ND;++n) {
//				pt(n) = x.crd2d(n)(j)(k);
//				mvel(n) = 0.0;//x.gbl->bd(0)*(x.crd2d(n)(j)(k) -dxdt(x.log2p,j)(n)(j)(k));
//			}
//			
//			u = x.u2d(0)(j)(k);
//			
//			x.res2d(0)(j)(k) = flux(u,pt,mvel,nrm);
//			
//		}
//	}
	
	
//	basis::tet(x.log2p).intgrt2d(&x.lf(0)(0),&x.res2d(0)(0)(0),MXGP);
	
	x.lf = 0.0;
	
	return;
}
