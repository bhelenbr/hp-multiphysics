/*
 *  l2error.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"
#include <utilities.h>
#include "hp_boundary.h"

void tet_hp::l2error(init_bdry_cndtn *comparison) {
	int i,j,k,n,tind;
	FLT err;
	Array<int,1> loc(NV);
	Array<FLT,1> mxr(NV),l2r(NV);
	
	TinyVector<FLT,3> pt;
	int stridey = MXGP;
	int stridex = MXGP*MXGP; 
	
	for(n=0;n<NV;++n) {
		mxr(n) = 0.0;
		l2r(n) = 0.0;
	}
	
	for(tind=0;tind<ntet;++tind) {
		if (tet(tind).info < 0) {
			for(n=0;n<ND;++n)
			basis::tet(log2p).proj(pnts(tet(tind).pnt(0))(n),pnts(tet(tind).pnt(1))(n),pnts(tet(tind).pnt(2))(n),pnts(tet(tind).pnt(3))(n),&crd(n)(0)(0)(0),stridex,stridey);

			for(i=0;i<basis::tet(log2p).gpx;++i) {
			for(j=0;j<basis::tet(log2p).gpy;++j) {
				for(k=0;k<basis::tet(log2p).gpz;++k) {
					for(n=0;n<ND;++n) {
						dcrd(n)(0)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(3))(n) -pnts(tet(tind).pnt(2))(n));
						dcrd(n)(1)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(1))(n) -pnts(tet(tind).pnt(2))(n));
						dcrd(n)(2)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(0))(n) -pnts(tet(tind).pnt(2))(n));
					}
				}
			}
			}
		}
		else {
			crdtocht(tind);
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
		}
		
		ugtouht(tind);
		for(n=0;n<NV;++n)
			basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);

		for(i=0;i<basis::tet(log2p).gpx;++i) {
			for(j=0;j<basis::tet(log2p).gpy;++j) {
				for(k=0;k<basis::tet(log2p).gpz;++k) {
					pt(0) = crd(0)(i)(j)(k);
					pt(1) = crd(1)(i)(j)(k);
					pt(2) = crd(2)(i)(j)(k);
					cjcb(i)(j)(k) = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
					for(n=0;n<NV;++n) {
						err =  fabs(u(n)(i)(j)(k)-comparison->f(n,pt,gbl->time));
						//err = fabs(u(n)(i)(j)(k)-1.0/5.0*exp(-5.0*gbl->time)*(1.0-pt(0))*pt(0));
						if (err > mxr(n)) {
							mxr(n) = err;
							loc(n) = tind;
						}
						l2r(n) += err*err*basis::tet(log2p).wtx(i)*basis::tet(log2p).wty(j)*basis::tet(log2p).wtz(k)*cjcb(i)(j)(k);
					}
				}
			}
		}		
	}
	
	for(n=0;n<NV;++n) {
		l2r(n) = sqrt(l2r(n)); 
		*gbl->log << "#L_2: " << l2r(n) << " L_inf " << mxr(n) <<  ' ' << loc(n) << ' ';
		//*gbl->log  << l2r(n) << ' ' << mxr(n) <<  ' ' ;
	}
	*gbl->log << '\n';

	return;
}

