/*
 *  length.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>
#include "tri_hp_cns_explicit.h"
#include "../hp_boundary.h"
#include <vector>


/* THIS FUNCTION WILL SET THE lngth VALUES BASED ON THE TRUNCATION ERROR */

void tri_hp_cns_explicit::error_estimator() {
	TinyMatrix<FLT,ND,ND> ldcrd;
	Array<TinyMatrix<FLT,MXGP,MXGP>,1> u(NV),ul(NV);
	Array<TinyMatrix<FLT,MXGP,MXGP>,2> du(NV,ND), dul(NV,ND);
	Array<FLT,1> cvu(NV);

	if (hp_gbl->error_estimator == hp_global::none)
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
		
		for(int n = 0; n < NV; ++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);
		
		// switch uht from conservative to primitive variables
		for(int i = 0; i < lgpx; ++i) {
			for(int j = 0; j < lgpn; ++j) {
				
				for(int n = 0; n < NV; ++n)
					cvu(n) = u(n)(i,j);
				
				u(0)(i,j) = (hp_cns_explicit_gbl->gamma-1.0)*(cvu(3)-0.5/cvu(0)*(cvu(1)*cvu(1)+cvu(2)*cvu(2))); // pressure
				u(1)(i,j) = cvu(1)/cvu(0);// x-velocity
				u(2)(i,j) = cvu(2)/cvu(0);// y-velocity
				u(3)(i,j) = u(0)(i,j)/cvu(0);// RT Temperature
			}
		}		

		for(int n = 1; n < NV; ++n) {
			du(n,0) = 0.0;
			du(n,1) = 0.0;
			basis::tri(log2p)->derivr(&u(n)(0,0),&du(n,0)(0,0),MXGP);
			basis::tri(log2p)->derivs(&u(n)(0,0),&du(n,1)(0,0),MXGP);
		}

	
		for(int n=0;n<NV;++n) 
			for(std::vector<int>::iterator it=highs.begin();it!=highs.end();++it)
				uht(n)(*it) = 0.0;
			
		for(int n = 0; n < NV; ++n)
			basis::tri(log2p)->proj(&uht(n)(0),&ul(n)(0,0),MXGP);
		
		// switch uht from conservative to primitive variables
		for(int i = 0; i < lgpx; ++i) {
			for(int j = 0; j < lgpn; ++j) {
				
				for(int n = 0; n < NV; ++n)
					cvu(n) = ul(n)(i,j);
				
				ul(0)(i,j) = (hp_cns_explicit_gbl->gamma-1.0)*(cvu(3)-0.5/cvu(0)*(cvu(1)*cvu(1)+cvu(2)*cvu(2)));// pressure
				ul(1)(i,j) = cvu(1)/cvu(0);// x-velocity
				ul(2)(i,j) = cvu(2)/cvu(0);// y-velocity
				ul(3)(i,j) = ul(0)(i,j)/cvu(0);// RT Temperature
			}
		}		
		
		for(int n = 1; n < NV; ++n) {
			dul(n,0) = 0.0;
			dul(n,1) = 0.0;
			basis::tri(log2p)->derivr(&ul(n)(0,0),&dul(n,0)(0,0),MXGP);
			basis::tri(log2p)->derivs(&ul(n)(0,0),&dul(n,1)(0,0),MXGP);
		}
		
		
		FLT jcb = 0.25*area(tind); 
		FLT error2 = 0.0; 
		FLT energy, denergy;
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {
				FLT dudx = ldcrd(1,1)*du(1,0)(i,j) -ldcrd(1,0)*du(1,1)(i,j);
				FLT dudy = -ldcrd(0,1)*du(1,0)(i,j) +ldcrd(0,0)*du(1,1)(i,j);
				FLT dvdx = ldcrd(1,1)*du(2,0)(i,j) -ldcrd(1,0)*du(2,1)(i,j);
				FLT dvdy = -ldcrd(0,1)*du(2,0)(i,j) +ldcrd(0,0)*du(2,1)(i,j);	
				
				FLT dudxl = ldcrd(1,1)*dul(1,0)(i,j) -ldcrd(1,0)*dul(1,1)(i,j);
				FLT dudyl = -ldcrd(0,1)*dul(1,0)(i,j) +ldcrd(0,0)*dul(1,1)(i,j);
				FLT dvdxl = ldcrd(1,1)*dul(2,0)(i,j) -ldcrd(1,0)*dul(2,1)(i,j);
				FLT dvdyl = -ldcrd(0,1)*dul(2,0)(i,j) +ldcrd(0,0)*dul(2,1)(i,j);		
				
				FLT rho = u(0)(i,j)/u(NV-1)(i,j);
				/* INVISCID PARTS TO ERROR MEASURE */
				energy = rho*(u(1)(i,j)*u(1)(i,j) +u(2)(i,j)*u(2)(i,j)) +u(0)(i,j);	
				/* VISCOUS PART TO ERROR MEASURE */
				energy += hp_cns_explicit_gbl->mu*(fabs(dudx)+fabs(dudy)+fabs(dvdx)+fabs(dvdy))/jcb;
				
				rho = ul(0)(i,j)/ul(NV-1)(i,j);

				denergy = rho*(ul(1)(i,j)*ul(1)(i,j) +ul(2)(i,j)*ul(2)(i,j)) +ul(0)(i,j);
				denergy += hp_cns_explicit_gbl->mu*(fabs(dudxl)+fabs(dudyl)+fabs(dvdxl)+fabs(dvdyl))/jcb;
				
				denergy -= energy;
				
				error2 += denergy*denergy*jcb*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j);
				totalenergy2 += energy*energy*jcb*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j);
			}
		}
		totalerror2 += error2;
		e2to_pow += pow(error2,1./(1.+alpha));
		tri_gbl->fltwk(tind) = error2;
	}

	/* Need to all-reduce norm,totalerror2,and totalenergy2 */
	hp_gbl->eanda(0) = totalenergy2;
	hp_gbl->eanda(1) = e2to_pow;
	hp_gbl->eanda(2) = totalerror2;
		
	return;
}
