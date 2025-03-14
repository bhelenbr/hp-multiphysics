/*
 *  l2error.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"

void tri_hp::l2error(init_bdry_cndtn *comparison) {
	int i,j,n,tind;
	FLT err;
	Array<int,1> loc(NV);
	Array<FLT,1> mxr(NV),l2r(NV);
	TinyVector<FLT,2> pt;

	for(n=0;n<NV;++n) {
		mxr(n) = 0.0;
		l2r(n) = 0.0;
	}

	for(tind=0;tind<ntri;++tind) {
        pmetric->calc_metrics(tind,crd,dcrd);
        
        /* Calculate position in parametric coordinates for ibc */
		if (tri(tind).info > -1) {
			crdtocht(tind);
			for(n=0;n<ND;++n)
				basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0) ,MXGP);
		}
		else {
			for(n=0;n<ND;++n)
				basis::tri(log2p)->proj(pnts(tri(tind).pnt(0))(n),pnts(tri(tind).pnt(1))(n),pnts(tri(tind).pnt(2))(n),&crd(n)(0,0),MXGP);
		}

		ugtouht(tind);
		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);

 		for (i=0;i<basis::tri(log2p)->gpx();++i) {	
			for (j=0;j<basis::tri(log2p)->gpn();++j) {
				cjcb(i,j) = (dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
				pt(0) = crd(0)(i,j);
				pt(1) = crd(1)(i,j);
				for(n=0;n<NV;++n) {
					err =  fabs(u(n)(i,j)-comparison->f(n,pt,gbl->time));
					if (err >= mxr(n)) {
						mxr(n) = err;
						loc(n) = tind;
					}
					l2r(n) += err*err*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j)*cjcb(i,j);
				}
			}
		}	
	}

	for(n=0;n<NV;++n) {
		l2r(n) = sqrt(l2r(n)); 
		*gbl->log << "#L_2: " << l2r(n) << " L_inf " << mxr(n) <<  ' ' << loc(n);
	}
	*gbl->log << std::endl;

	return;
}

/* CALCULATE AREA/CIRCUMFERENCE/YBAR */
void tri_hp::integrated_averages(Array<FLT,1> a) {
	int i,j,n,tind;

	/* a(0) = area */
	/* a(1) = xbar */
	/* a(2) = ybar */
	/* a(3-...) variable averages */
	a = 0.0;

	for(tind=0;tind<ntri;++tind) {
        pmetric->calc_metrics(tind, crd, dcrd);

		ugtouht(tind);
		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);

		for(i=0;i<basis::tri(log2p)->gpx();++i) {
			for(j=0;j<basis::tri(log2p)->gpn();++j) {
				cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
				a(0) += RAD(crd(0)(i,j))*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j)*cjcb(i,j);
				a(1) += crd(0)(i,j)*RAD(crd(0)(i,j))*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j)*cjcb(i,j);
				a(2) += crd(1)(i,j)*RAD(crd(0)(i,j))*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j)*cjcb(i,j);
				for(n=0;n<NV;++n) 
					a(3+n) += u(n)(i,j)*RAD(crd(0)(i,j))*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j)*cjcb(i,j);
			}
		}
	}
	a(Range(1,2+NV)) /= a(0);

	return;
}

