/*
 *  tobasis.cpp
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"tri_hp.h"
#include<myblas.h>

void tri_hp::tobasis(init_bdry_cndtn *ibc, int tlvl) {
	int tind,i,j,m,n,indx,v0,v1,sind,info;
	char uplo[] = "U";
	TinyVector<FLT,2> pt;

	/* LOOP THROUGH VERTICES */
	for(i=0;i<npnt;++i)
		for(n=0;n<NV;++n)
			ugbd(tlvl).v(i,n) = ibc->f(n,pnts(i),gbl->time);

    if (basis::tri(log2p)->sm() <= 0) return;

    /* LOOP THROUGH SIDES */    
    for(sind=0;sind<nseg;++sind) {

		v0 = seg(sind).pnt(0);
		v1 = seg(sind).pnt(1);

		if (seg(sind).info < 0) {
			for(n=0;n<ND;++n)
				basis::tri(log2p)->proj1d(pnts(v0)(n),pnts(v1)(n),&crd(n)(0,0));
		}
		else {
			crdtocht1d(sind,tlvl);
			for(n=0;n<ND;++n)
				basis::tri(log2p)->proj1d(&cht(n,0),&crd(n)(0,0));
		}

		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj1d(ugbd(tlvl).v(v0,n),ugbd(tlvl).v(v1,n),&res(n)(0,0));

		for(i=0;i<basis::tri(log2p)->gpx(); ++i) {
			pt(0) = crd(0)(0,i);
			pt(1) = crd(1)(0,i);
			for(n=0;n<NV;++n)
				res(n)(0,i) -= ibc->f(n,pt,gbl->time);
		}

		for(n=0;n<NV;++n)
			basis::tri(log2p)->intgrt1d(&lf(n)(0),&res(n)(0,0));

		indx = sind;
		for(n=0;n<NV;++n) {
			PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&lf(n)(2),basis::tri(log2p)->sm(),info);
			for(m=0;m<basis::tri(log2p)->sm();++m) 
				ugbd(tlvl).s(sind,m,n) = -lf(n)(2+m);
		}
	}

	if (basis::tri(log2p)->im() <= 0) return;

	for(tind = 0; tind < ntri; ++tind) {
		ugtouht_bdry(tind,tlvl);
		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj_bdry(&uht(n)(0),&u(n)(0,0),MXGP);

		if (tri(tind).info < 0) {
			for(n=0;n<ND;++n)
				basis::tri(log2p)->proj(vrtxbd(tlvl)(tri(tind).pnt(0))(n),vrtxbd(tlvl)(tri(tind).pnt(1))(n),vrtxbd(tlvl)(tri(tind).pnt(2))(n),&crd(n)(0,0),MXGP);
		}
		else {
			crdtocht(tind,tlvl);
			for(n=0;n<ND;++n)
				basis::tri(log2p)->proj_bdry(&cht(n,0),&crd(n)(0,0),MXGP);
		}

		for (i=0; i < basis::tri(log2p)->gpx(); ++i ) {
			for (j=0; j < basis::tri(log2p)->gpn(); ++j ) {
				pt(0) = crd(0)(i,j);
				pt(1) = crd(1)(i,j);
				for(n=0;n<NV;++n)
					u(n)(i,j) -= ibc->f(n,pt,gbl->time);
			}
		}

		for(n=0;n<NV;++n) {
			basis::tri(log2p)->intgrt(&lf(n)(0),&u(n)(0,0),MXGP);
			DPBTRS(uplo,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),1,(double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,&lf(n)(basis::tri(log2p)->bm()),basis::tri(log2p)->im(),info);
			for(i=0;i<basis::tri(log2p)->im();++i)
				ugbd(tlvl).i(tind,i,n) = -lf(n)(basis::tri(log2p)->bm()+i);
		}
	}

	return;
}    

