/*
 *  lengthcd.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on Thu May 29 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include "../hp_boundary.h"
#include <math.h>
#include <utilities.h>

void tri_hp_cd::error_estimator() {
	TinyMatrix<FLT,ND,ND> ldcrd;
	Array<TinyMatrix<FLT,MXGP,MXGP>,1> u(NV),ul(NV);
	Array<TinyMatrix<FLT,MXGP,MXGP>,2> du(NV,ND), dul(NV,ND);

	if (gbl->error_estimator == global::none) 
	return;

	int sm = basis::tri(log2p)->sm();
	int lgpx = basis::tri(log2p)->gpx();
	int lgpn = basis::tri(log2p)->gpn();
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

	/* USING energy CONSTANT AS ERROR INDICATOR */
	/* Real convergence rate is p+1/2 (for pressure in L_2) */
	/* Real convergence rate of the error for this norm will be probably p-1/2 because it has derivatives */
	/* This norm is measuring error in p-1 solution not pth order solution */
	/* If solution was optimal converence of derivative in this norm would be p-1 (so this the lower bound) */
	/* alpha includes weighting due to area of element +2 */
	const FLT alpha = 2.0*(basis::tri(log2p)->p()-1.0+ND)/static_cast<FLT>(ND);
	FLT e2to_pow = 0.0, totalenergy2 = 0.0, totalerror2 = 0.0;
	for (int tind=0;tind<ntri;++tind) {
		
		/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		TinyVector<int,3> v = tri(tind).pnt;
		/* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
		for(int n=0;n<ND;++n) {
			ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
			ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
		}
		
		ugtouht(tind);
		
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
			
		for(int n=0;n<NV;++n) 
			for(std::vector<int>::iterator it=highs.begin();it!=highs.end();++it)
				uht(n)(*it) = 0.0;
				
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&ul(n)(0,0),&dul(n,0)(0,0),&dul(n,1)(0,0),MXGP);
		
		FLT jcb = 0.25*area(tind); 
		FLT error2 = 0.0; 
		FLT energy, denergy;
		
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {
	
				//				FLT dtdx = ldcrd(1,1)*du(0,0)(i,j) -ldcrd(1,0)*du(0,1)(i,j);
				//				FLT dtdy = -ldcrd(0,1)*du(0,0)(i,j) +ldcrd(0,0)*du(0,1)(i,j);
				
				//				FLT dtdxl = ldcrd(1,1)*dul(0,0)(i,j) -ldcrd(1,0)*dul(0,1)(i,j);
				//				FLT dtdyl = -ldcrd(0,1)*dul(0,0)(i,j) +ldcrd(0,0)*dul(0,1)(i,j);
				
				/* INVISCID PARTS TO ERROR MEASURE */
				energy = gbl->rhocv*u(0)(i,j);	
				/* VISCOUS PART TO ERROR MEASURE */
				energy += (0 /* +gbl->kcond*(fabs(dtdx)+fabs(dtdy)) */ )/jcb;
				
				
				denergy = gbl->rhocv*ul(0)(i,j);
				denergy += (0 /* +gbl->kcond*(fabs(dtdxl)+fabs(dtdyl)) */)/jcb;
				
				denergy -= energy;
				
				error2 += denergy*denergy*jcb*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j);
				totalenergy2 += energy*energy*jcb*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j);
			}
		}
		totalerror2 += error2;
		e2to_pow += pow(error2,1./(1.+alpha));
		gbl->fltwk(tind) = error2;
	}

	/* Need to all-reduce norm,totalerror2,and totalenergy2 */
	gbl->eanda(0) = totalenergy2;
	gbl->eanda(1) = e2to_pow;
	gbl->eanda(2) = totalerror2;

	return;
}

#ifdef OLDWAY
void tri_hp_cd::length() {
	int i,j,k,v0,v1,v2,sind,tind,count;
	TinyVector<FLT,2> dx0,dx1,dx2,ep,dedpsi;
	FLT sum,ratio;
	FLT length0,length1,length2,lengthept;
	FLT ang1,curved1,ang2,curved2;
	
	return;  // TEMPORARY

	gbl->fltwk(Range(0,npnt-1)) = 0.0;

	switch(basis::tri(log2p)->p()) {
		case(1): {
			for(i=0;i<nseg;++i) {
				v0 = seg(i).pnt(0);
				v1 = seg(i).pnt(1);
				sum = distance2(v0,v1)*fabs(ug.v(v0,0) -ug.v(v1,0));
				gbl->fltwk(v0) += sum;
				gbl->fltwk(v1) += sum;
			}
			break;
		}

		default: {
			for(i=0;i<nseg;++i) {
				v0 = seg(i).pnt(0);
				v1 = seg(i).pnt(1);
				sum = distance2(v0,v1)*fabs(ug.s(i,sm0-1,0));
				gbl->fltwk(v0) += sum;
				gbl->fltwk(v1) += sum;
				/* UNCOMMENT FOR SHOCK DETECTION 
				sum = abs(ug.s(i,sm0-1,0)/(abs(ug.s(i,0,0))+1.0e-5));
				gbl->fltwk(v0) = MAX(sum,gbl->fltwk(v0));
				gbl->fltwk(v1) = MAX(sum,gbl->fltwk(v1)); */
			}

			/* BOUNDARY CURVATURE */
			for(i=0;i<nebd;++i) {
				if (!(hp_ebdry(i)->is_curved())) continue;

				for(j=0;j<ebdry(i)->nseg;++j) {
					sind = ebdry(i)->seg(j);
					v1 = seg(sind).pnt(0);
					v2 = seg(sind).pnt(1);

					crdtocht1d(sind);

					/* FIND ANGLE BETWEEN LINEAR SIDES */
					tind = seg(sind).tri(0);
					for(k=0;k<3;++k)
						if (tri(tind).seg(k) == sind) break;

					v0 = tri(tind).pnt(k);

					dx0(0) = pnts(v2)(0)-pnts(v1)(0);
					dx0(1) = pnts(v2)(1)-pnts(v1)(1);
					length0 = dx0(0)*dx0(0) +dx0(1)*dx0(1);

					dx1(0) = pnts(v0)(0)-pnts(v2)(0);
					dx1(1) = pnts(v0)(1)-pnts(v2)(1);
					length1 = dx1(0)*dx1(0) +dx1(1)*dx1(1);

					dx2(0) = pnts(v1)(0)-pnts(v0)(0);
					dx2(1) = pnts(v1)(1)-pnts(v0)(1);
					length2 = dx2(0)*dx2(0) +dx2(1)*dx2(1);

					basis::tri(log2p)->ptprobe1d(2,&ep(0),&dedpsi(0),-1.0,&cht(0,0),MXTM);
					lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);

					ang1 = acos(-(dx0(0)*dx2(0) +dx0(1)*dx2(1))/sqrt(length0*length2));
					curved1 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));

					basis::tri(log2p)->ptprobe1d(2,&ep(0),&dedpsi(0),1.0,&cht(0,0),MXTM);
					lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);

					ang2 = acos(-(dx0(0)*dx1(0) +dx0(1)*dx1(1))/sqrt(length0*length1));
					curved2 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));                            

					sum = gbl->curvature_sensitivity*(curved1/ang1 +curved2/ang2);
					gbl->fltwk(v0) += sum*gbl->error_target*pnt(v0).nnbor;
					gbl->fltwk(v1) += sum*gbl->error_target*pnt(v1).nnbor;
				}
			}
			break;
		}
	}

	// output_error(); FOR SHOCK DETECTION

	for(i=0;i<npnt;++i) {
		gbl->fltwk(i) = pow(gbl->fltwk(i)/(pnt(i).nnbor*gbl->error_target),1./(basis::tri(log2p)->p()+1+ND));
		lngth(i) /= gbl->fltwk(i);  
		lngth(i) = MAX(lngth(i),gbl->minlngth);
		lngth(i) = MIN(lngth(i),gbl->maxlngth);
    }

    /* AVOID HIGH ASPECT RATIOS */
    int nsweep = 0;
    do {
		count = 0;
		for(i=0;i<nseg;++i) {
			v0 = seg(i).pnt(0);
			v1 = seg(i).pnt(1);
			ratio = lngth(v1)/lngth(v0);

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
#endif
