/*
 *  length.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
#include <utilities.h>

#include "tri_hp_ins.h"
#include "../hp_boundary.h"
#include <blitz/tinyvec-et.h>
#include <vector>

/* THIS FUNCTION WILL SET THE lngth VALUES BASED ON THE TRUNCATION ERROR */

void tri_hp_ins::length() {
	TinyMatrix<FLT,ND,ND> ldcrd;
	Array<TinyMatrix<FLT,MXGP,MXGP>,1> u(NV),ul(NV);
	Array<TinyMatrix<FLT,MXGP,MXGP>,2> du(NV,ND), dul(NV,ND);
		
	int sm = basis::tri(log2p).sm;
	int lgpx = basis::tri(log2p).gpx;
	int lgpn = basis::tri(log2p).gpn;
	std::vector<int> highs;
	highs.push_back(2+sm);
	highs.push_back(2+2*sm);
	highs.push_back(2+3*sm);
	int indx = 3+3*sm;
	for(int m = 1; m < sm; ++m) {
		for(int k = 0; k < sm-m-1; ++k) {
			++indx;
		}
		highs.push_back(indx++);
	}

	/* USING BERNOULLI CONSTANT AS ERROR INDICATOR */
	const FLT alpha = 2.0*(basis::tri(log2p).p-0.5)/static_cast<FLT>(ND);
	FLT denom = 0.0, totalbernoulli = 0.0, totalerror = 0.0;
	for (int tind=0;tind<ntri;++tind) {
			
		/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		TinyVector<int,3> v = tri(tind).pnt;
		/* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
		for(int n=0;n<ND;++n) {
			ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
			ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
		}
			
		ugtouht(tind);
		
		for(int n=0;n<ND;++n)
			basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
		basis::tri(log2p).proj(&uht(NV-1)(0),&u(NV-1)(0,0),MXGP);
		
		for(int n=0;n<NV;++n) 
			for(std::vector<int>::iterator it=highs.begin();it!=highs.end();++it)
				uht(n)(*it) = 0.0;
			
		for(int n=0;n<ND;++n)
			basis::tri(log2p).proj(&uht(n)(0),&ul(n)(0,0),&dul(n,0)(0,0),&dul(n,1)(0,0),MXGP);
		basis::tri(log2p).proj(&uht(NV-1)(0),&ul(NV-1)(0,0),MXGP);

		FLT jcb = 0.25*area(tind); 
		FLT error2 = 0.0; 
		FLT bernoulli, dbernoulli;
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {
				FLT dudx = ldcrd(1,1)*du(0,0)(i,j) -ldcrd(1,0)*du(0,1)(i,j);
				FLT dudy = -ldcrd(0,1)*du(0,0)(i,j) +ldcrd(0,0)*du(0,1)(i,j);
				FLT dvdx = ldcrd(1,1)*du(1,0)(i,j) -ldcrd(1,0)*du(1,1)(i,j);
				FLT dvdy = -ldcrd(0,1)*du(1,0)(i,j) +ldcrd(0,0)*du(1,1)(i,j);	
				
				FLT dudxl = ldcrd(1,1)*dul(0,0)(i,j) -ldcrd(1,0)*dul(0,1)(i,j);
				FLT dudyl = -ldcrd(0,1)*dul(0,0)(i,j) +ldcrd(0,0)*dul(0,1)(i,j);
				FLT dvdxl = ldcrd(1,1)*dul(1,0)(i,j) -ldcrd(1,0)*dul(1,1)(i,j);
				FLT dvdyl = -ldcrd(0,1)*dul(1,0)(i,j) +ldcrd(0,0)*dul(1,1)(i,j);		
				
				/* INVISCID PARTS TO ERROR MEASURE */
				bernoulli = gbl->rho*(u(0)(i,j)*u(0)(i,j) +u(1)(i,j)*u(1)(i,j)) +u(NV-1)(i,j);	
				/* VISCOUS PART TO ERROR MEASURE */
				bernoulli += gbl->mu*(fabs(dudx)+fabs(dudy)+fabs(dvdx)+fabs(dvdy))/jcb;
				
				dbernoulli = gbl->rho*(ul(0)(i,j)*ul(0)(i,j) +ul(1)(i,j)*ul(1)(i,j)) +ul(NV-1)(i,j);
				dbernoulli += gbl->mu*(fabs(dudxl)+fabs(dudyl)+fabs(dvdxl)+fabs(dvdyl))/jcb;
				
				dbernoulli -= bernoulli;
				
				error2 += dbernoulli*dbernoulli*jcb*basis::tri(log2p).wtx(i)*basis::tri(log2p).wtn(j);
				totalbernoulli += bernoulli*bernoulli*jcb*basis::tri(log2p).wtx(i)*basis::tri(log2p).wtn(j);
			}
		}
		totalerror += error2;
		denom += pow(error2,1./(1.+alpha));
		gbl->fltwk(tind) = error2;
	}

	/* Need to all-reduce norm,totalerror,and totalbernoulli */
	gbl->eanda(0) = totalbernoulli;
	gbl->eanda(1) = denom;
	gbl->eanda(2) = totalerror;
	sim::blks.allreduce(gbl->eanda.data(),gbl->eanda_recv.data(),3,blocks::flt_msg,blocks::sum);
	
	/* Determine error target (SEE AEA Paper) */
	FLT etarget2 = gbl->error_target*gbl->error_target*gbl->eanda_recv(0);
	*gbl->log << "# Normalized Error " << sqrt(gbl->eanda_recv(2)/gbl->eanda_recv(0)) << " Target " << gbl->error_target << '\n';
	FLT K = pow(etarget2/denom,1./(ND*alpha));
	gbl->res.v(0,Range(0,npnt-1)) = 1.0;
	gbl->res_r.v(0,Range(0,npnt-1)) = 0.0;
	for(int tind=0;tind<ntri;++tind) {
		FLT error2 = gbl->fltwk(tind)+FLT_EPSILON;
		FLT ri = K*pow(error2, -1./(ND*(1.+alpha)));
		for (int j=0;j<3;++j) {
			int p0 = tri(tind).pnt(j);
			/* Calculate average at vertices */
			gbl->res.v(0,p0) *= ri;
			gbl->res_r.v(0,p0) += 1.0;
		}
	}
	
	/* NOW RESCALE AT VERTICES */
	FLT maxlngth = 50.0;
	FLT minlngth = 0.0;
	for (int pind=0;pind<npnt;++pind) {
		lngth(pind) *= pow(gbl->res.v(0,pind),1.0/gbl->res_r.v(0,pind));
		lngth(pind) = MIN(lngth(pind),maxlngth);
		lngth(pind) = MAX(lngth(pind),minlngth);
	}
	
	/* LIMIT BOUNDARY CURVATURE */
	for(int i=0;i<nebd;++i) {
		if (!(hp_ebdry(i)->is_curved())) continue;

		for(int j=0;j<ebdry(i)->nseg;++j) {
			int sind = ebdry(i)->seg(j);
			int v1 = seg(sind).pnt(0);
			int v2 = seg(sind).pnt(1);

			crdtocht1d(sind);

			/* FIND ANGLE BETWEEN LINEAR SIDES */
			int tind = seg(sind).tri(0);
			int k;
			for(k=0;k<3;++k)
				if (tri(tind).seg(k) == sind) break;

			int v0 = tri(tind).pnt(k);

			TinyVector<FLT,ND> dx0;
			dx0(0) = pnts(v2)(0)-pnts(v1)(0);
			dx0(1) = pnts(v2)(1)-pnts(v1)(1);
			FLT length0 = dx0(0)*dx0(0) +dx0(1)*dx0(1);

			TinyVector<FLT,ND> dx1;
			dx1(0) = pnts(v0)(0)-pnts(v2)(0);
			dx1(1) = pnts(v0)(1)-pnts(v2)(1);
			FLT length1 = dx1(0)*dx1(0) +dx1(1)*dx1(1);

			TinyVector<FLT,ND> dx2;
			dx2(0) = pnts(v1)(0)-pnts(v0)(0);
			dx2(1) = pnts(v1)(1)-pnts(v0)(1);
			FLT length2 = dx2(0)*dx2(0) +dx2(1)*dx2(1);
			
			TinyVector<FLT,2> ep, dedpsi;
			basis::tri(log2p).ptprobe1d(2,&ep(0),&dedpsi(0),-1.0,&cht(0,0),MXTM);
			FLT lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);

			FLT ang1 = acos(-(dx0(0)*dx2(0) +dx0(1)*dx2(1))/sqrt(length0*length2));
			FLT curved1 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));

			basis::tri(log2p).ptprobe1d(2,&ep(0),&dedpsi(0),1.0,&cht(0,0),MXTM);
			lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);

			FLT ang2 = acos(-(dx0(0)*dx1(0) +dx0(1)*dx1(1))/sqrt(length0*length1));
			FLT curved2 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));                            

			// FIXME: end points are wrong for periodic boundary or communication boundary
			FLT sum = gbl->curvature_sensitivity*(fabs(curved1/ang1) +fabs(curved2/ang2));
			lngth(v1) /= 1. +sum;
			lngth(v2) /= 1. +sum;
		}
	}

	/* AVOID HIGH ASPECT RATIOS */
	int nsweep = 0;
	int count;
	do {
		count = 0;
		for(int i=0;i<nseg;++i) {
			int v0 = seg(i).pnt(0);
			int v1 = seg(i).pnt(1);
			FLT ratio = lngth(v1)/lngth(v0);

			if (ratio > 3.0) {
				lngth(v1) = 2.5*lngth(v0);
				++count;
			}
			else if (ratio < 0.333) {
				lngth(v0) = 2.5*lngth(v1);
				++count;
			}
		}
		++nsweep;
		*gbl->log << "#aspect ratio fixes " << nsweep << ' ' << count << std::endl;
	} while(count > 0 && nsweep < 5);

	return;
}
