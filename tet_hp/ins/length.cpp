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

#include "tet_hp_ins.h"
#include "../hp_boundary.h"
//#include <blitz/tinyvec-et.h>
#include <vector>

/* THIS FUNCTION WILL SET THE lngth VALUES BASED ON THE TRUNCATION ERROR */

void tet_hp_ins::length() {
	TinyMatrix<FLT,ND,ND> ldcrd;
	TinyMatrix<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,4,3> du,dul;
	Array<TinyVector<TinyVector<TinyVector<double,MXGP>,MXGP>,MXGP>,1> u,ul;  
	TinyVector<TinyVector<FLT,ND>,ND> d;

	int em = basis::tet(log2p).em;
	int lgpx = basis::tet(log2p).gpx;
	int lgpy = basis::tet(log2p).gpy;
	int lgpz = basis::tet(log2p).gpz;
	int stridey = MXGP;
	int stridex = MXGP*MXGP; 
	
	std::vector<int> highs;
	highs.push_back(2+em);
	highs.push_back(2+2*em);
	highs.push_back(2+3*em);
	int indx = 3+3*em;
	//fix me temp
	for(int m = 1; m < em; ++m) {
		for(int k = 0; k < em-m-1; ++k) {
			++indx;
		}
		highs.push_back(indx++);
	}
	
	/* USING BERNOULLI CONSTANT AS ERROR INDICATOR */
	const FLT alpha = 2.0*(basis::tet(log2p).p-0.5)/static_cast<FLT>(ND);
	FLT denom = 0.0, totalbernoulli = 0.0, totalerror = 0.0;
	for (int tind=0;tind<ntet;++tind) {
		
		/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		TinyVector<int,4> v = tet(tind).pnt;
		/* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
		for(int n=0;n<ND;++n) {
			ldcrd(n,0) = 0.5*(pnts(v(3))(n) -pnts(v(2))(n));
			ldcrd(n,1) = 0.5*(pnts(v(1))(n) -pnts(v(2))(n));
			ldcrd(n,2) = 0.5*(pnts(v(0))(n) -pnts(v(2))(n));

		}
		d(0)(0) =  ldcrd(1,1)*ldcrd(2,2)-ldcrd(1,2)*ldcrd(2,1);//dr/dx
		d(0)(1) = -ldcrd(0,1)*ldcrd(2,2)+ldcrd(0,2)*ldcrd(2,1);//dr/dy
		d(0)(2) =  ldcrd(0,1)*ldcrd(1,2)-ldcrd(0,2)*ldcrd(1,1);//dr/dz
		d(1)(0) = -ldcrd(1,0)*ldcrd(2,2)+ldcrd(1,2)*ldcrd(2,0);//ds/dx
		d(1)(1) =  ldcrd(0,0)*ldcrd(2,2)-ldcrd(0,2)*ldcrd(2,0);//ds/dy
		d(1)(2) = -ldcrd(0,0)*ldcrd(1,2)+ldcrd(0,2)*ldcrd(1,0);//ds/dz
		d(2)(0) =  ldcrd(1,0)*ldcrd(2,1)-ldcrd(1,1)*ldcrd(2,0);//dt/dx
		d(2)(1) = -ldcrd(0,0)*ldcrd(2,1)+ldcrd(0,1)*ldcrd(2,0);//dt/dy
		d(2)(2) =  ldcrd(0,0)*ldcrd(1,1)-ldcrd(0,1)*ldcrd(1,0);//dt/dz
		
		ugtouht(tind);
		
		for(int n=0;n<ND;++n)
			basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),&du(n,0)(0)(0)(0),&du(n,1)(0)(0)(0),&du(n,2)(0)(0)(0),stridex,stridey);
		basis::tet(log2p).proj(&uht(NV-1)(0),&u(NV-1)(0)(0)(0),stridex,stridey);
		
		for(int n=0;n<NV;++n) 
			for(std::vector<int>::iterator it=highs.begin();it!=highs.end();++it)
				uht(n)(*it) = 0.0;
		
		for(int n=0;n<ND;++n)
			basis::tet(log2p).proj(&uht(n)(0),&ul(n)(0)(0)(0),&dul(n,0)(0)(0)(0),&dul(n,1)(0)(0)(0),&dul(n,2)(0)(0)(0),stridex,stridey);
		basis::tet(log2p).proj(&uht(NV-1)(0),&ul(NV-1)(0)(0)(0),stridex,stridey);
		
		FLT jcb = 0.125*volume(tind); 
		FLT error2 = 0.0; 
		FLT bernoulli, dbernoulli;
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpy;++j) {
				for(int k=0;k<lgpz;++k) {

					FLT dudx = d(0)(0)*du(0,0)(i)(j)(k)+d(1)(0)*du(0,1)(i)(j)(k)+d(2)(0)*du(0,2)(i)(j)(k);
					FLT dudy = d(0)(1)*du(0,0)(i)(j)(k)+d(1)(1)*du(0,1)(i)(j)(k)+d(2)(1)*du(0,2)(i)(j)(k);
					FLT dudz = d(0)(2)*du(0,0)(i)(j)(k)+d(1)(2)*du(0,1)(i)(j)(k)+d(2)(2)*du(0,2)(i)(j)(k);

					FLT dvdx = d(0)(0)*du(1,0)(i)(j)(k)+d(1)(0)*du(1,1)(i)(j)(k)+d(2)(0)*du(1,2)(i)(j)(k);
					FLT dvdy = d(0)(1)*du(1,0)(i)(j)(k)+d(1)(1)*du(1,1)(i)(j)(k)+d(2)(1)*du(1,2)(i)(j)(k);
					FLT dvdz = d(0)(2)*du(1,0)(i)(j)(k)+d(1)(2)*du(1,1)(i)(j)(k)+d(2)(2)*du(1,2)(i)(j)(k);

					FLT dwdx = d(0)(0)*du(2,0)(i)(j)(k)+d(1)(0)*du(2,1)(i)(j)(k)+d(2)(0)*du(2,2)(i)(j)(k);
					FLT dwdy = d(0)(1)*du(2,0)(i)(j)(k)+d(1)(1)*du(2,1)(i)(j)(k)+d(2)(1)*du(2,2)(i)(j)(k);
					FLT dwdz = d(0)(2)*du(2,0)(i)(j)(k)+d(1)(2)*du(2,1)(i)(j)(k)+d(2)(2)*du(2,2)(i)(j)(k);
					
					FLT dudxl = d(0)(0)*dul(0,0)(i)(j)(k)+d(1)(0)*dul(0,1)(i)(j)(k)+d(2)(0)*dul(0,2)(i)(j)(k);
					FLT dudyl = d(0)(1)*dul(0,0)(i)(j)(k)+d(1)(1)*dul(0,1)(i)(j)(k)+d(2)(1)*dul(0,2)(i)(j)(k);
					FLT dudzl = d(0)(2)*dul(0,0)(i)(j)(k)+d(1)(2)*dul(0,1)(i)(j)(k)+d(2)(2)*dul(0,2)(i)(j)(k);
					
					FLT dvdxl = d(0)(0)*dul(1,0)(i)(j)(k)+d(1)(0)*dul(1,1)(i)(j)(k)+d(2)(0)*dul(1,2)(i)(j)(k);
					FLT dvdyl = d(0)(1)*dul(1,0)(i)(j)(k)+d(1)(1)*dul(1,1)(i)(j)(k)+d(2)(1)*dul(1,2)(i)(j)(k);
					FLT dvdzl = d(0)(2)*dul(1,0)(i)(j)(k)+d(1)(2)*dul(1,1)(i)(j)(k)+d(2)(2)*dul(1,2)(i)(j)(k);
					
					FLT dwdxl = d(0)(0)*dul(2,0)(i)(j)(k)+d(1)(0)*dul(2,1)(i)(j)(k)+d(2)(0)*dul(2,2)(i)(j)(k);
					FLT dwdyl = d(0)(1)*dul(2,0)(i)(j)(k)+d(1)(1)*dul(2,1)(i)(j)(k)+d(2)(1)*dul(2,2)(i)(j)(k);
					FLT dwdzl = d(0)(2)*dul(2,0)(i)(j)(k)+d(1)(2)*dul(2,1)(i)(j)(k)+d(2)(2)*dul(2,2)(i)(j)(k);	
					
					/* INVISCID PARTS TO ERROR MEASURE */
					bernoulli = gbl->rho*(u(0)(i)(j)(k)*u(0)(i)(j)(k) +u(1)(i)(j)(k)*u(1)(i)(j)(k)+u(2)(i)(j)(k)*u(2)(i)(j)(k)) +u(NV-1)(i)(j)(k);	
					/* VISCOUS PART TO ERROR MEASURE */
					bernoulli += gbl->mu*(fabs(dudx)+fabs(dudy)+fabs(dudz)+fabs(dvdx)+fabs(dvdy)+fabs(dvdz)+fabs(dwdx)+fabs(dwdy)+fabs(dwdz))/jcb;
					
					dbernoulli = gbl->rho*(ul(0)(i)(j)(k)*ul(0)(i)(j)(k) +ul(1)(i)(j)(k)*ul(1)(i)(j)(k)+ul(2)(i)(j)(k)*ul(2)(i)(j)(k)) +ul(NV-1)(i)(j)(k);	
					dbernoulli += gbl->mu*(fabs(dudxl)+fabs(dudyl)+fabs(dudzl)+fabs(dvdxl)+fabs(dvdyl)+fabs(dvdzl)+fabs(dwdxl)+fabs(dwdyl)+fabs(dwdzl))/jcb;
					
					dbernoulli -= bernoulli;
					
					error2 += dbernoulli*dbernoulli*jcb*basis::tet(log2p).wtx(i)*basis::tet(log2p).wty(j)*basis::tet(log2p).wtz(j);
					totalbernoulli += bernoulli*bernoulli*jcb*basis::tet(log2p).wtx(i)*basis::tet(log2p).wty(j)*basis::tet(log2p).wtz(j);
				}
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
	for(int tind=0;tind<ntet;++tind) {
		FLT error2 = gbl->fltwk(tind)+FLT_EPSILON;
		FLT ri = K*pow(error2, -1./(ND*(1.+alpha)));
		for (int j=0;j<4;++j) {
			int p0 = tet(tind).pnt(j);
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
	
//	/* LIMIT BOUNDARY CURVATURE */
//	for(int i=0;i<nebd;++i) {
//		if (!(hp_ebdry(i)->is_curved())) continue;
//		
//		for(int j=0;j<ebdry(i)->nseg;++j) {
//			int sind = ebdry(i)->seg(j);
//			int v1 = seg(sind).pnt(0);
//			int v2 = seg(sind).pnt(1);
//			
//			crdtocht1d(sind);
//			
//			/* FIND ANGLE BETWEEN LINEAR SIDES */
//			int tind = seg(sind).tri(0);
//			int k;
//			for(k=0;k<3;++k)
//				if (tri(tind).seg(k) == sind) break;
//			
//			int v0 = tri(tind).pnt(k);
//			
//			TinyVector<FLT,ND> dx0;
//			dx0(0) = pnts(v2)(0)-pnts(v1)(0);
//			dx0(1) = pnts(v2)(1)-pnts(v1)(1);
//			FLT length0 = dx0(0)*dx0(0) +dx0(1)*dx0(1);
//			
//			TinyVector<FLT,ND> dx1;
//			dx1(0) = pnts(v0)(0)-pnts(v2)(0);
//			dx1(1) = pnts(v0)(1)-pnts(v2)(1);
//			FLT length1 = dx1(0)*dx1(0) +dx1(1)*dx1(1);
//			
//			TinyVector<FLT,ND> dx2;
//			dx2(0) = pnts(v1)(0)-pnts(v0)(0);
//			dx2(1) = pnts(v1)(1)-pnts(v0)(1);
//			FLT length2 = dx2(0)*dx2(0) +dx2(1)*dx2(1);
//			
//			TinyVector<FLT,2> ep, dedpsi;
//			basis::tri(log2p)->ptprobe1d(2,&ep(0),&dedpsi(0),-1.0,&cht(0,0),MXTM);
//			FLT lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);
//			
//			FLT ang1 = acos(-(dx0(0)*dx2(0) +dx0(1)*dx2(1))/sqrt(length0*length2));
//			FLT curved1 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));
//			
//			basis::tri(log2p)->ptprobe1d(2,&ep(0),&dedpsi(0),1.0,&cht(0,0),MXTM);
//			lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);
//			
//			FLT ang2 = acos(-(dx0(0)*dx1(0) +dx0(1)*dx1(1))/sqrt(length0*length1));
//			FLT curved2 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));                            
//			
//			// FIXME: end points are wrong for periodic boundary or communication boundary
//			FLT sum = gbl->curvature_sensitivity*(fabs(curved1/ang1) +fabs(curved2/ang2));
//			lngth(v1) /= 1. +sum;
//			lngth(v2) /= 1. +sum;
//		}
//	}
	
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
