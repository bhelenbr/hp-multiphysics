/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cns.h"
#include "../hp_boundary.h"

// #define BODYFORCE
	
void tet_hp_cns::element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im){
	
	int i,j,k,n;
	FLT fluxx,fluxy,fluxz;
	const int NV = 5;
	TinyVector<int,4> v;
	TinyMatrix<FLT,ND,ND> ldcrd;
	TinyMatrix<FLT,NV,NV> A,B,C;
	TinyMatrix<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,NV,ND> du;
	int lgpx = basis::tet(log2p).gpx, lgpy = basis::tet(log2p).gpy, lgpz = basis::tet(log2p).gpz;
	int stridey = MXGP;
	int stridex = MXGP*MXGP; 
	FLT rhobd0 = gbl->rho*gbl->bd(0), lmu = gbl->mu, rhorbd0, cjcb, cjcbi;
	FLT lkcond = gbl->kcond;
	TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-2,NV-2> visc;
	TinyVector<TinyVector<FLT,ND>,ND> d,kcond;// temp fix gbl->kcond and kcond problem?
	TinyMatrix<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,NV,NV> cv, df;
	TinyVector<FLT,NV> tres;
	FLT ogm1 = 1.0/(gbl->gamma-1.0);
	FLT gogm1 = gbl->gamma*ogm1;
	
	/* LOAD INDICES OF VERTEX POINTS */
	v = tet(tind).pnt;
	
	/* IF TINFO > -1 IT IS CURVED ELEMENT */
	if (tet(tind).info > -1 ) {
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		crdtocht(tind);
		
		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(n=0;n<ND;++n)
			basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
	}
	else {
		/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(n=0;n<ND;++n)
			basis::tet(log2p).proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),pnts(v(3))(n),&crd(n)(0)(0)(0),stridex,stridey);
		
		/* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
		for(n=0;n<ND;++n) {
			ldcrd(n,0) = 0.5*(pnts(tet(tind).pnt(3))(n) -pnts(tet(tind).pnt(2))(n));
			ldcrd(n,1) = 0.5*(pnts(tet(tind).pnt(1))(n) -pnts(tet(tind).pnt(2))(n));
			ldcrd(n,2) = 0.5*(pnts(tet(tind).pnt(0))(n) -pnts(tet(tind).pnt(2))(n));
		}
	}
	
	/* CALCULATE MESH VELOCITY */
	for(i=0;i<lgpx;++i) {
		for(j=0;j<lgpy;++j) {
			for(k=0;k<lgpz;++k) {
				mvel(0)(i)(j)(k) = 0.0;//gbl->bd(0)*(crd(0)(i)(j)(k) -dxdt(log2p,tind,0)(i)(j)(k));
				mvel(1)(i)(j)(k) = 0.0;//gbl->bd(0)*(crd(1)(i)(j)(k) -dxdt(log2p,tind,1)(i)(j)(k));
				mvel(2)(i)(j)(k) = 0.0;//gbl->bd(0)*(crd(2)(i)(j)(k) -dxdt(log2p,tind,2)(i)(j)(k));
			}
		}
	}
	

	
	/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */

	if (gbl->beta(stage) > 0.0) {
		basis::tet(log2p).proj(&uht(0)(0),&u(NV-1)(0)(0)(0),stridex,stridey);
		for(n=1;n<NV;++n)
			basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),&du(n,0)(0)(0)(0),&du(n,1)(0)(0)(0),&du(n,2)(0)(0)(0),stridex,stridey);
		
	}
	else {
		for(n=0;n<NV;++n)
			basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);
	}
	
	/* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
	for(n=0;n<NV;++n){
		for(i=0;i<basis::tet(log2p).tm;++i){
			lf_re(n)(i) = 0.0;
			lf_im(n)(i) = 0.0;	
		}
	}
	
	if (tet(tind).info > -1) {
		/* CURVED ELEMENT */
		cout << " curvy element being called in rsdl"<<endl;
		/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpy;++j) {
				for(k=0;k<lgpz;++k) {
					
					d(0)(0) =  dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k);
					d(0)(1) = -dcrd(0)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(0)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k);
					d(0)(2) =  dcrd(0)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)-dcrd(0)(2)(i)(j)(k)*dcrd(1)(1)(i)(j)(k);
					d(1)(0) = -dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
					d(1)(1) =  dcrd(0)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(0)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
					d(1)(2) = -dcrd(0)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)+dcrd(0)(2)(i)(j)(k)*dcrd(1)(0)(i)(j)(k);
					d(2)(0) =  dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
					d(2)(1) = -dcrd(0)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)+dcrd(0)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
					d(2)(2) =  dcrd(0)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k)-dcrd(0)(1)(i)(j)(k)*dcrd(1)(0)(i)(j)(k);
					
					double rho = u(0)(i)(j)(k)/u(NV-1)(i)(j)(k);

					fluxx = rho*(u(1)(i)(j)(k) -mvel(0)(i)(j)(k));
					fluxy = rho*(u(2)(i)(j)(k) -mvel(1)(i)(j)(k));
					fluxz = rho*(u(3)(i)(j)(k) -mvel(2)(i)(j)(k));
					
					/* CONTINUITY EQUATION FLUXES */
					cv(0,0)(i)(j)(k) = d(0)(0)*fluxx+d(0)(1)*fluxy+d(0)(2)*fluxz;
					cv(0,1)(i)(j)(k) = d(1)(0)*fluxx+d(1)(1)*fluxy+d(1)(2)*fluxz;					
					cv(0,2)(i)(j)(k) = d(2)(0)*fluxx+d(2)(1)*fluxy+d(2)(2)*fluxz;
					
					
					/* MOMENTUM & ENERGY FLUXES */
					for(n=1;n<NV-1;++n) {
						cv(n,0)(i)(j)(k) = u(n)(i)(j)(k)*cv(0,0)(i)(j)(k);
						cv(n,1)(i)(j)(k) = u(n)(i)(j)(k)*cv(0,1)(i)(j)(k);
						cv(n,2)(i)(j)(k) = u(n)(i)(j)(k)*cv(0,2)(i)(j)(k);
					}

					
					/* PRESSURE TERMS */
					/* U-MOMENTUM */
					cv(1,0)(i)(j)(k) += d(0)(0)*u(0)(i)(j)(k);
					cv(1,1)(i)(j)(k) += d(1)(0)*u(0)(i)(j)(k);
					cv(1,2)(i)(j)(k) += d(2)(0)*u(0)(i)(j)(k);
					
					/* V-MOMENTUM */
					cv(2,0)(i)(j)(k) += d(0)(1)*u(0)(i)(j)(k);
					cv(2,1)(i)(j)(k) += d(1)(1)*u(0)(i)(j)(k);
					cv(2,2)(i)(j)(k) += d(2)(1)*u(0)(i)(j)(k);
					
					/* W-MOMENTUM */
					cv(3,0)(i)(j)(k) += d(0)(2)*u(0)(i)(j)(k);
					cv(3,1)(i)(j)(k) += d(1)(2)*u(0)(i)(j)(k);
					cv(3,2)(i)(j)(k) += d(2)(2)*u(0)(i)(j)(k);
					
					/* ENERGY FLUX */
					double h = gogm1*u(NV-1)(i)(j)(k)+0.5*(u(1)(i)(j)(k)*u(1)(i)(j)(k)+u(2)(i)(j)(k)*u(2)(i)(j)(k)+u(3)(i)(j)(k)*u(3)(i)(j)(k));
					cv(4,0)(i)(j)(k) = h*cv(0,0)(i)(j)(k);
					cv(4,1)(i)(j)(k) = h*cv(0,1)(i)(j)(k);
					cv(4,2)(i)(j)(k) = h*cv(0,2)(i)(j)(k);

				}
			}
		}
		for(n=0;n<NV;++n)
			basis::tet(log2p).intgrtrst(&lf_im(n)(0),&cv(n,0)(0)(0)(0),&cv(n,1)(0)(0)(0),&cv(n,2)(0)(0)(0),stridex,stridey);
		
		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			/* TIME DERIVATIVE TERMS */ 
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						
						d(0)(0) =  dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k);
						d(0)(1) = -dcrd(0)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(0)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k);
						d(0)(2) =  dcrd(0)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)-dcrd(0)(2)(i)(j)(k)*dcrd(1)(1)(i)(j)(k);
						d(1)(0) = -dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
						d(1)(1) =  dcrd(0)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(0)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
						d(1)(2) = -dcrd(0)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)+dcrd(0)(2)(i)(j)(k)*dcrd(1)(0)(i)(j)(k);
						d(2)(0) =  dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
						d(2)(1) = -dcrd(0)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)+dcrd(0)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
						d(2)(2) =  dcrd(0)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k)-dcrd(0)(1)(i)(j)(k)*dcrd(1)(0)(i)(j)(k);
						
						double rho = u(0)(i)(j)(k)/u(NV-1)(i)(j)(k);
						cjcb = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
						
						rhorbd0 = rho*gbl->bd(0)*cjcb;
						mujcbi = lmu/cjcb;
						kcjcbi = lkcond/cjcb/gbl->R;
						
						/* UNSTEADY TERMS */
						res(0)(i)(j)(k) = rhorbd0+dugdt(log2p,tind,0)(i)(j)(k);
						for(n=1;n<NV-1;++n)
							res(n)(i)(j)(k) = rhorbd0*u(n)(i)(j)(k) +dugdt(log2p,tind,n)(i)(j)(k);
						double e = ogm1*u(NV-1)(i)(j)(k) +0.5*(u(1)(i)(j)(k)*u(1)(i)(j)(k)+u(2)(i)(j)(k)*u(2)(i)(j)(k)+u(3)(i)(j)(k)*u(3)(i)(j)(k));
						res(NV-1)(i)(j)(k) = rhorbd0*e +dugdt(log2p,tind,NV-1)(i)(j)(k);
						
#ifdef BODYFORCE
						for(n=1;n<NV-1;++n)
							res(n)(i)(j)(k) -= rho*cjcb*gbl->body(n-1);			
						
#endif                        
						/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
						/* INDICES ARE 1: EQUATION U V OR W, 2: VARIABLE (U V OR W), 3: EQ. DERIVATIVE (R S OR T) 4: VAR DERIVATIVE (R S OR T)*/
						visc(0,0)(0,0) = -mujcbi*(4./3.*d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1)+d(0)(2)*d(0)(2));
						visc(0,0)(0,1) = -mujcbi*(4./3.*d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1)+d(0)(2)*d(1)(2));
						visc(0,0)(0,2) = -mujcbi*(4./3.*d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1)+d(0)(2)*d(2)(2));
						visc(0,0)(1,0) = visc(0,0)(0,1);
						visc(0,0)(1,1) = -mujcbi*(4./3.*d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1)+d(1)(2)*d(1)(2));
						visc(0,0)(1,2) = -mujcbi*(4./3.*d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1)+d(1)(2)*d(2)(2));
						visc(0,0)(2,0) = visc(0,0)(0,2);
						visc(0,0)(2,1) = visc(0,0)(1,2);
						visc(0,0)(2,2) = -mujcbi*(4./3.*d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1)+d(2)(2)*d(2)(2));
						
						visc(0,1)(0,0) = -mujcbi*1./3.*d(0)(1)*d(0)(0);
						visc(0,1)(0,1) = -mujcbi*(d(0)(1)*d(1)(0)-2./3.*d(0)(0)*d(1)(1));
						visc(0,1)(0,2) = -mujcbi*(d(0)(1)*d(2)(0)-2./3.*d(0)(0)*d(2)(1));
						visc(0,1)(1,0) = -mujcbi*(d(1)(1)*d(0)(0)-2./3.*d(1)(0)*d(0)(1));
						visc(0,1)(1,1) = -mujcbi*1./3.*d(1)(1)*d(1)(0);
						visc(0,1)(1,2) = -mujcbi*(d(1)(1)*d(2)(0)-2./3.*d(1)(0)*d(2)(1));
						visc(0,1)(2,0) = -mujcbi*(d(2)(1)*d(0)(0)-2./3.*d(2)(0)*d(0)(1));
						visc(0,1)(2,1) = -mujcbi*(d(2)(1)*d(1)(0)-2./3.*d(2)(0)*d(1)(1));
						visc(0,1)(2,2) = -mujcbi*1./3.*d(2)(1)*d(2)(0);
						
						visc(0,2)(0,0) = -mujcbi*1./3.*d(0)(2)*d(0)(0);
						visc(0,2)(0,1) = -mujcbi*(d(0)(2)*d(1)(0)-2./3.*d(0)(0)*d(1)(2));
						visc(0,2)(0,2) = -mujcbi*(d(0)(2)*d(2)(0)-2./3.*d(0)(0)*d(2)(2));
						visc(0,2)(1,0) = -mujcbi*(d(1)(2)*d(0)(0)-2./3.*d(1)(0)*d(0)(2));
						visc(0,2)(1,1) = -mujcbi*1./3.*d(1)(2)*d(1)(0);
						visc(0,2)(1,2) = -mujcbi*(d(1)(2)*d(2)(0)-2./3.*d(1)(0)*d(2)(2));
						visc(0,2)(2,0) = -mujcbi*(d(2)(2)*d(0)(0)-2./3.*d(2)(0)*d(0)(2));
						visc(0,2)(2,1) = -mujcbi*(d(2)(2)*d(1)(0)-2./3.*d(2)(0)*d(1)(2));
						visc(0,2)(2,2) = -mujcbi*1./3.*d(2)(2)*d(2)(0);
						
						visc(1,0)(0,0) = visc(0,1)(0,0);
						visc(1,0)(0,1) = visc(0,1)(1,0);
						visc(1,0)(0,2) = visc(0,1)(2,0);
						visc(1,0)(1,0) = visc(0,1)(0,1);
						visc(1,0)(1,1) = visc(0,1)(1,1);
						visc(1,0)(1,2) = visc(0,1)(2,1);
						visc(1,0)(2,0) = visc(0,1)(0,2);
						visc(1,0)(2,1) = visc(0,1)(1,2);
						visc(1,0)(2,2) = visc(0,1)(2,2);
						
						visc(1,1)(0,0) = -mujcbi*(4./3.*d(0)(1)*d(0)(1)+d(0)(0)*d(0)(0)+d(0)(2)*d(0)(2));
						visc(1,1)(0,1) = -mujcbi*(4./3.*d(0)(1)*d(1)(1)+d(0)(0)*d(1)(0)+d(0)(2)*d(1)(2));
						visc(1,1)(0,2) = -mujcbi*(4./3.*d(0)(1)*d(2)(1)+d(0)(0)*d(2)(0)+d(0)(2)*d(2)(2));
						visc(1,1)(1,0) = visc(1,1)(0,1);
						visc(1,1)(1,1) = -mujcbi*(4./3.*d(1)(1)*d(1)(1)+d(1)(0)*d(1)(0)+d(1)(2)*d(1)(2));
						visc(1,1)(1,2) = -mujcbi*(4./3.*d(1)(1)*d(2)(1)+d(1)(0)*d(2)(0)+d(1)(2)*d(2)(2));
						visc(1,1)(2,0) = visc(1,1)(0,2);
						visc(1,1)(2,1) = visc(1,1)(1,2);
						visc(1,1)(2,2) = -mujcbi*(4./3.*d(2)(1)*d(2)(1)+d(2)(0)*d(2)(0)+d(2)(2)*d(2)(2));
						
						visc(1,2)(0,0) = -mujcbi*1./3.*d(0)(2)*d(0)(1);
						visc(1,2)(0,1) = -mujcbi*(d(0)(2)*d(1)(1)-2./3.*d(0)(1)*d(1)(2));
						visc(1,2)(0,2) = -mujcbi*(d(0)(2)*d(2)(1)-2./3.*d(0)(1)*d(2)(2));
						visc(1,2)(1,0) = -mujcbi*(d(1)(2)*d(0)(1)-2./3.*d(1)(1)*d(0)(2));
						visc(1,2)(1,1) = -mujcbi*1./3.*d(1)(2)*d(1)(1);
						visc(1,2)(1,2) = -mujcbi*(d(1)(2)*d(2)(1)-2./3.*d(1)(1)*d(2)(2));
						visc(1,2)(2,0) = -mujcbi*(d(2)(2)*d(0)(1)-2./3.*d(2)(1)*d(0)(2));
						visc(1,2)(2,1) = -mujcbi*(d(2)(2)*d(1)(1)-2./3.*d(2)(1)*d(1)(2));
						visc(1,2)(2,2) = -mujcbi*1./3.*d(2)(2)*d(2)(1);
						
						visc(2,0)(0,0) = visc(0,2)(0,0);
						visc(2,0)(0,1) = visc(0,2)(1,0);
						visc(2,0)(0,2) = visc(0,2)(2,0);
						visc(2,0)(1,0) = visc(0,2)(0,1);
						visc(2,0)(1,1) = visc(0,2)(1,1);
						visc(2,0)(1,2) = visc(0,2)(2,1);
						visc(2,0)(2,0) = visc(0,2)(0,2);
						visc(2,0)(2,1) = visc(0,2)(1,2);
						visc(2,0)(2,2) = visc(0,2)(2,2);
						
						visc(2,1)(0,0) = visc(1,2)(0,0);
						visc(2,1)(0,1) = visc(1,2)(1,0);
						visc(2,1)(0,2) = visc(1,2)(2,0);
						visc(2,1)(1,0) = visc(1,2)(0,1);
						visc(2,1)(1,1) = visc(1,2)(1,1);
						visc(2,1)(1,2) = visc(1,2)(2,1);
						visc(2,1)(2,0) = visc(1,2)(0,2);
						visc(2,1)(2,1) = visc(1,2)(1,2);
						visc(2,1)(2,2) = visc(1,2)(2,2);
						
						visc(2,2)(0,0) = -mujcbi*(4./3.*d(0)(2)*d(0)(2)+d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1));
						visc(2,2)(0,1) = -mujcbi*(4./3.*d(0)(2)*d(1)(2)+d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1));
						visc(2,2)(0,2) = -mujcbi*(4./3.*d(0)(2)*d(2)(2)+d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1));
						visc(2,2)(1,0) = visc(2,2)(0,1);
						visc(2,2)(1,1) = -mujcbi*(4./3.*d(1)(2)*d(1)(2)+d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1));
						visc(2,2)(1,2) = -mujcbi*(4./3.*d(1)(2)*d(2)(2)+d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1));
						visc(2,2)(2,0) = visc(2,2)(0,2);
						visc(2,2)(2,1) = visc(2,2)(1,2);
						visc(2,2)(2,2) = -mujcbi*(4./3.*d(2)(2)*d(2)(2)+d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1));
						
						/* HEAT DIFFUSION TENSOR */
						kcond(0)(0) = -cjcb*(d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1)+d(0)(2)*d(0)(2));
						kcond(1)(1) = -cjcb*(d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1)+d(1)(2)*d(1)(2));
						kcond(2)(2) = -cjcb*(d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1)+d(2)(2)*d(2)(2));
						kcond(0)(1) = -cjcb*(d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1)+d(0)(2)*d(1)(2));
						kcond(0)(2) = -cjcb*(d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1)+d(0)(2)*d(2)(2));
						kcond(1)(2) = -cjcb*(d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1)+d(1)(2)*d(2)(2));
						
#define					kcondI1II0I  kcond(0)(1)
#define					kcondI2II0I  kcond(0)(2)
#define					kcondI2II1I  kcond(1)(2)

						
						/*  MOMENTUM EQUATIONS */
						df(1,0)(i)(j)(k) = +visc(0,0)(0,0)*du(1,0)(i)(j)(k)+visc(0,1)(0,0)*du(2,0)(i)(j)(k)+visc(0,2)(0,0)*du(3,0)(i)(j)(k)
										   +visc(0,0)(0,1)*du(1,1)(i)(j)(k)+visc(0,1)(0,1)*du(2,1)(i)(j)(k)+visc(0,2)(0,1)*du(3,1)(i)(j)(k)
										   +visc(0,0)(0,2)*du(1,2)(i)(j)(k)+visc(0,1)(0,2)*du(2,2)(i)(j)(k)+visc(0,2)(0,2)*du(3,2)(i)(j)(k);
						
						df(1,1)(i)(j)(k) = +visc(0,0)(1,0)*du(1,0)(i)(j)(k)+visc(0,1)(1,0)*du(2,0)(i)(j)(k)+visc(0,2)(1,0)*du(3,0)(i)(j)(k)
										   +visc(0,0)(1,1)*du(1,1)(i)(j)(k)+visc(0,1)(1,1)*du(2,1)(i)(j)(k)+visc(0,2)(1,1)*du(3,1)(i)(j)(k)
										   +visc(0,0)(1,2)*du(1,2)(i)(j)(k)+visc(0,1)(1,2)*du(2,2)(i)(j)(k)+visc(0,2)(1,2)*du(3,2)(i)(j)(k);
						
						df(1,2)(i)(j)(k) = +visc(0,0)(2,0)*du(1,0)(i)(j)(k)+visc(0,1)(2,0)*du(2,0)(i)(j)(k)+visc(0,2)(2,0)*du(3,0)(i)(j)(k)
										   +visc(0,0)(2,1)*du(1,1)(i)(j)(k)+visc(0,1)(2,1)*du(2,1)(i)(j)(k)+visc(0,2)(2,1)*du(3,1)(i)(j)(k)
										   +visc(0,0)(2,2)*du(1,2)(i)(j)(k)+visc(0,1)(2,2)*du(2,2)(i)(j)(k)+visc(0,2)(2,2)*du(3,2)(i)(j)(k);
						
						df(2,0)(i)(j)(k) = +visc(1,0)(0,0)*du(1,0)(i)(j)(k)+visc(1,1)(0,0)*du(2,0)(i)(j)(k)+visc(1,2)(0,0)*du(3,0)(i)(j)(k)
										   +visc(1,0)(0,1)*du(1,1)(i)(j)(k)+visc(1,1)(0,1)*du(2,1)(i)(j)(k)+visc(1,2)(0,1)*du(3,1)(i)(j)(k)
										   +visc(1,0)(0,2)*du(1,2)(i)(j)(k)+visc(1,1)(0,2)*du(2,2)(i)(j)(k)+visc(1,2)(0,2)*du(3,2)(i)(j)(k);
						
						df(2,1)(i)(j)(k) = +visc(1,0)(1,0)*du(1,0)(i)(j)(k)+visc(1,1)(1,0)*du(2,0)(i)(j)(k)+visc(1,2)(1,0)*du(3,0)(i)(j)(k)
										   +visc(1,0)(1,1)*du(1,1)(i)(j)(k)+visc(1,1)(1,1)*du(2,1)(i)(j)(k)+visc(1,2)(1,1)*du(3,1)(i)(j)(k)
										   +visc(1,0)(1,2)*du(1,2)(i)(j)(k)+visc(1,1)(1,2)*du(2,2)(i)(j)(k)+visc(1,2)(1,2)*du(3,2)(i)(j)(k);
						
						df(2,2)(i)(j)(k) = +visc(1,0)(2,0)*du(1,0)(i)(j)(k)+visc(1,1)(2,0)*du(2,0)(i)(j)(k)+visc(1,2)(2,0)*du(3,0)(i)(j)(k)
										   +visc(1,0)(2,1)*du(1,1)(i)(j)(k)+visc(1,1)(2,1)*du(2,1)(i)(j)(k)+visc(1,2)(2,1)*du(3,1)(i)(j)(k)
										   +visc(1,0)(2,2)*du(1,2)(i)(j)(k)+visc(1,1)(2,2)*du(2,2)(i)(j)(k)+visc(1,2)(2,2)*du(3,2)(i)(j)(k);
						
						df(3,0)(i)(j)(k) = +visc(2,0)(0,0)*du(1,0)(i)(j)(k)+visc(2,1)(0,0)*du(2,0)(i)(j)(k)+visc(2,2)(0,0)*du(3,0)(i)(j)(k)
										   +visc(2,0)(0,1)*du(1,1)(i)(j)(k)+visc(2,1)(0,1)*du(2,1)(i)(j)(k)+visc(2,2)(0,1)*du(3,1)(i)(j)(k)
										   +visc(2,0)(0,2)*du(1,2)(i)(j)(k)+visc(2,1)(0,2)*du(2,2)(i)(j)(k)+visc(2,2)(0,2)*du(3,2)(i)(j)(k);
						
						df(3,1)(i)(j)(k) = +visc(2,0)(1,0)*du(1,0)(i)(j)(k)+visc(2,1)(1,0)*du(2,0)(i)(j)(k)+visc(2,2)(1,0)*du(3,0)(i)(j)(k)
										   +visc(2,0)(1,1)*du(1,1)(i)(j)(k)+visc(2,1)(1,1)*du(2,1)(i)(j)(k)+visc(2,2)(1,1)*du(3,1)(i)(j)(k)
										   +visc(2,0)(1,2)*du(1,2)(i)(j)(k)+visc(2,1)(1,2)*du(2,2)(i)(j)(k)+visc(2,2)(1,2)*du(3,2)(i)(j)(k);
						
						df(3,2)(i)(j)(k) = +visc(2,0)(2,0)*du(1,0)(i)(j)(k)+visc(2,1)(2,0)*du(2,0)(i)(j)(k)+visc(2,2)(2,0)*du(3,0)(i)(j)(k)
										   +visc(2,0)(2,1)*du(1,1)(i)(j)(k)+visc(2,1)(2,1)*du(2,1)(i)(j)(k)+visc(2,2)(2,1)*du(3,1)(i)(j)(k)
										   +visc(2,0)(2,2)*du(1,2)(i)(j)(k)+visc(2,1)(2,2)*du(2,2)(i)(j)(k)+visc(2,2)(2,2)*du(3,2)(i)(j)(k);
						
						
						/*  ENGERY EQUATION */
						df(4,0)(i)(j)(k) = kcond(0)(0)*du(NV-1,0)(i)(j)(k)+kcond(0)(1)*du(NV-1,1)(i)(j)(k)+kcond(0)(2)*du(NV-1,2)(i)(j)(k);
						df(4,1)(i)(j)(k) = kcondI1II0I*du(NV-1,0)(i)(j)(k)+kcond(1)(1)*du(NV-1,1)(i)(j)(k)+kcond(1)(2)*du(NV-1,2)(i)(j)(k);
						df(4,2)(i)(j)(k) = kcondI2II0I*du(NV-1,0)(i)(j)(k)+kcondI2II1I*du(NV-1,1)(i)(j)(k)+kcond(2)(2)*du(NV-1,2)(i)(j)(k);
						
						/* VISCOUS DISSIPATION */
						df(4,0)(i)(j)(k) += df(1,0)(i)(j)(k)*u(1)(i)(j)(k) +df(2,0)(i)(j)(k)*u(2)(i)(j)(k)+df(3,0)(i)(j)(k)*u(3)(i)(j)(k);
						df(4,1)(i)(j)(k) += df(1,1)(i)(j)(k)*u(1)(i)(j)(k) +df(2,1)(i)(j)(k)*u(2)(i)(j)(k)+df(3,1)(i)(j)(k)*u(3)(i)(j)(k);
						df(4,2)(i)(j)(k) += df(1,2)(i)(j)(k)*u(1)(i)(j)(k) +df(2,2)(i)(j)(k)*u(2)(i)(j)(k)+df(3,2)(i)(j)(k)*u(3)(i)(j)(k);

						for(n=1;n<NV;++n) {
							cv(n,0)(i)(j)(k) += df(n,0)(i)(j)(k);
							cv(n,1)(i)(j)(k) += df(n,1)(i)(j)(k);
							cv(n,2)(i)(j)(k) += df(n,2)(i)(j)(k);
						}
					}
				}
			}
			for(n=0;n<NV;++n)
				basis::tet(log2p).intgrt(&lf_re(n)(0),&res(n)(0)(0)(0),stridex,stridey);
			
			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(n=0;n<NV;++n) {
				basis::tet(log2p).derivr(&cv(n,0)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
				basis::tet(log2p).derivs(&cv(n,1)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
				basis::tet(log2p).derivt(&cv(n,2)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
			}
			

			
			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						
						d(0)(0) =  dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k);
						d(0)(1) = -dcrd(0)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(0)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k);
						d(0)(2) =  dcrd(0)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)-dcrd(0)(2)(i)(j)(k)*dcrd(1)(1)(i)(j)(k);
						d(1)(0) = -dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
						d(1)(1) =  dcrd(0)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(0)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
						d(1)(2) = -dcrd(0)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)+dcrd(0)(2)(i)(j)(k)*dcrd(1)(0)(i)(j)(k);
						d(2)(0) =  dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
						d(2)(1) = -dcrd(0)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)+dcrd(0)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k);
						d(2)(2) =  dcrd(0)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k)-dcrd(0)(1)(i)(j)(k)*dcrd(1)(0)(i)(j)(k);
						
						df(0,0)(i)(j)(k) = 0.0;
						df(0,1)(i)(j)(k) = 0.0;
						
						tres = 0.0;
						for(int m = 0; m < NV; ++m)
							for(int n = 0; n < NV; ++n)							
								tres(m) += gbl->tau(tind,m,n)*res(n)(i)(j)(k);
												
						FLT pr = u(0)(i)(j)(k);
						FLT uv = u(1)(i)(j)(k);
						FLT vv = u(2)(i)(j)(k);
						FLT wv = u(3)(i)(j)(k);
						FLT rt = u(4)(i)(j)(k);					
						FLT ke = 0.5*(uv*uv+vv*vv+wv*wv);
						
						/* df/dw */
						A = 1.0/rt*uv,               pr/rt,                           0.0,         0.0,         -pr/rt/rt*uv,
						    1.0/rt*uv*uv+1.0,        2.0*pr/rt*uv,                    0.0,         0.0,         -pr/rt/rt*uv*uv,
						    1.0/rt*uv*vv,            pr/rt*vv,                        pr/rt*uv,    0.0,         -pr/rt/rt*uv*vv,
                            1.0/rt*uv*wv,            pr/rt*wv,                        0.0,         pr/rt*uv,    -pr/rt/rt*uv*wv,
						    1.0/rt*uv*(gogm1*rt+ke), pr/rt*(gogm1*rt+ke)+pr/rt*uv*uv, pr/rt*uv*vv, pr/rt*uv*wv, -pr/rt/rt*uv*(gogm1*rt+ke)+pr/rt*uv*gogm1;
						
						
						/* dg/dw */
						B = 1.0/rt*vv,               0.0,         pr/rt,                           0.0,         -pr/rt/rt*vv,
						    1.0/rt*uv*vv,            pr/rt*vv,    pr/rt*uv,                        0.0,         -pr/rt/rt*uv*vv,
						    1.0/rt*vv*vv+1.0,        0.0,         2.0*pr/rt*vv,                    0.0,         -pr/rt/rt*vv*vv,
						    1.0/rt*vv*wv,            0.0,         pr/rt*wv,                        pr/rt*vv,    -pr/rt/rt*vv*wv,
						    1.0/rt*vv*(gogm1*rt+ke), pr/rt*uv*vv, pr/rt*(gogm1*rt+ke)+pr/rt*vv*vv, pr/rt*vv*wv,	-pr/rt/rt*vv*(gogm1*rt+ke)+pr/rt*vv*gogm1;
						
						
						
						/* dh/dw */
						C = 1.0/rt*wv,               0.0,         0.0,         pr/rt,                           -pr/rt/rt*wv,
						    1.0/rt*uv*wv,            pr/rt*wv,    0.0,         pr/rt*uv,                        -pr/rt/rt*uv*wv,
						    1.0/rt*vv*wv,            0.0,         pr/rt*wv,    pr/rt*vv,                        -pr/rt/rt*vv*wv,
						    1.0/rt*wv*wv+1.0,        0.0,         0.0,         2.0*pr/rt*wv,                    -pr/rt/rt*wv*wv,
							1.0/rt*wv*(gogm1*rt+ke), pr/rt*uv*wv, pr/rt*vv*wv, pr/rt*(gogm1*rt+ke)+pr/rt*wv*wv, -pr/rt/rt*wv*(gogm1*rt+ke)+pr/rt*wv*gogm1;				
											
						for(int m = 0; m < NV; ++m) {
							for(int n = 0; n < NV; ++n) {
								df(m,0)(i)(j)(k) -= (d(0)(0)*A(m,n)+d(0)(1)*B(m,n)+d(0)(2)*C(m,n))*tres(n);
								df(m,1)(i)(j)(k) -= (d(1)(0)*A(m,n)+d(1)(1)*B(m,n)+d(1)(2)*C(m,n))*tres(n);
								df(m,2)(i)(j)(k) -= (d(2)(0)*A(m,n)+d(2)(1)*B(m,n)+d(2)(2)*C(m,n))*tres(n);
							}
						}
					}
				}
			}
			
			for(n=0;n<NV;++n)
				basis::tet(log2p).intgrtrst(&lf_re(n)(0),&df(n,0)(0)(0)(0),&df(n,1)(0)(0)(0),&df(n,2)(0)(0)(0),stridex,stridey);
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
					lf_re(n)(i) *= gbl->beta(stage);


		}
	}
	else {
		// drdx dsdx dtdx [  ys*zt-zs*yt, -yr*zt+zr*yt,  yr*zs-zr*ys]
		// drdy dsdy dtdy [ -xs*zt+zs*xt,  xr*zt-zr*xt, -xr*zs+zr*xs]
		// drdz dsdz dtdz [  xs*yt-ys*xt, -xr*yt+yr*xt,  xr*ys-yr*xs]
		
		// cjcb = xr*(ys*zt-zs*yt)-xs*(yr*zt-zr*yt)+xt*(yr*zs-zr*ys)
		
		/* LINEAR ELEMENT */
		d(0)(0) =  ldcrd(1,1)*ldcrd(2,2)-ldcrd(1,2)*ldcrd(2,1);//dr/dx
		d(0)(1) = -ldcrd(0,1)*ldcrd(2,2)+ldcrd(0,2)*ldcrd(2,1);//dr/dy
		d(0)(2) =  ldcrd(0,1)*ldcrd(1,2)-ldcrd(0,2)*ldcrd(1,1);//dr/dz
		d(1)(0) = -ldcrd(1,0)*ldcrd(2,2)+ldcrd(1,2)*ldcrd(2,0);//ds/dx
		d(1)(1) =  ldcrd(0,0)*ldcrd(2,2)-ldcrd(0,2)*ldcrd(2,0);//ds/dy
		d(1)(2) = -ldcrd(0,0)*ldcrd(1,2)+ldcrd(0,2)*ldcrd(1,0);//ds/dz
		d(2)(0) =  ldcrd(1,0)*ldcrd(2,1)-ldcrd(1,1)*ldcrd(2,0);//dt/dx
		d(2)(1) = -ldcrd(0,0)*ldcrd(2,1)+ldcrd(0,1)*ldcrd(2,0);//dt/dy
		d(2)(2) =  ldcrd(0,0)*ldcrd(1,1)-ldcrd(0,1)*ldcrd(1,0);//dt/dz
	
		cjcb = ldcrd(0,0)*(ldcrd(1,1)*ldcrd(2,2)-ldcrd(1,2)*ldcrd(2,1))-ldcrd(0,1)*(ldcrd(1,0)*ldcrd(2,2)-ldcrd(2,0)*ldcrd(1,2))+ldcrd(0,2)*(ldcrd(1,0)*ldcrd(2,1)-ldcrd(2,0)*ldcrd(1,1));
		
		/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpy;++j) {
				for(k=0;k<lgpz;++k) {
					
					double rho = u(0)(i)(j)(k)/u(NV-1)(i)(j)(k);
					
					fluxx = rho*(u(1)(i)(j)(k) -mvel(0)(i)(j)(k));
					fluxy = rho*(u(2)(i)(j)(k) -mvel(1)(i)(j)(k));
					fluxz = rho*(u(3)(i)(j)(k) -mvel(2)(i)(j)(k));
					
					/* CONTINUITY EQUATION FLUXES */
					cv(0,0)(i)(j)(k) = d(0)(0)*fluxx+d(0)(1)*fluxy+d(0)(2)*fluxz;
					cv(0,1)(i)(j)(k) = d(1)(0)*fluxx+d(1)(1)*fluxy+d(1)(2)*fluxz;					
					cv(0,2)(i)(j)(k) = d(2)(0)*fluxx+d(2)(1)*fluxy+d(2)(2)*fluxz;
					
					
					/* MOMENTUM & ENERGY FLUXES */
					for(n=1;n<NV-1;++n) {
						cv(n,0)(i)(j)(k) = u(n)(i)(j)(k)*cv(0,0)(i)(j)(k);
						cv(n,1)(i)(j)(k) = u(n)(i)(j)(k)*cv(0,1)(i)(j)(k);
						cv(n,2)(i)(j)(k) = u(n)(i)(j)(k)*cv(0,2)(i)(j)(k);
					}
					
					
					/* PRESSURE TERMS */
					/* U-MOMENTUM */
					cv(1,0)(i)(j)(k) += d(0)(0)*u(0)(i)(j)(k);
					cv(1,1)(i)(j)(k) += d(1)(0)*u(0)(i)(j)(k);
					cv(1,2)(i)(j)(k) += d(2)(0)*u(0)(i)(j)(k);
					
					/* V-MOMENTUM */
					cv(2,0)(i)(j)(k) += d(0)(1)*u(0)(i)(j)(k);
					cv(2,1)(i)(j)(k) += d(1)(1)*u(0)(i)(j)(k);
					cv(2,2)(i)(j)(k) += d(2)(1)*u(0)(i)(j)(k);
					
					/* W-MOMENTUM */
					cv(3,0)(i)(j)(k) += d(0)(2)*u(0)(i)(j)(k);
					cv(3,1)(i)(j)(k) += d(1)(2)*u(0)(i)(j)(k);
					cv(3,2)(i)(j)(k) += d(2)(2)*u(0)(i)(j)(k);
					
					/* ENERGY FLUX */
					double h = gogm1*u(NV-1)(i)(j)(k)+0.5*(u(1)(i)(j)(k)*u(1)(i)(j)(k)+u(2)(i)(j)(k)*u(2)(i)(j)(k)+u(3)(i)(j)(k)*u(3)(i)(j)(k));
					cv(4,0)(i)(j)(k) = h*cv(0,0)(i)(j)(k);
					cv(4,1)(i)(j)(k) = h*cv(0,1)(i)(j)(k);
					cv(4,2)(i)(j)(k) = h*cv(0,2)(i)(j)(k);
					
				}
			}
		}
		for(n=0;n<NV;++n)
			basis::tet(log2p).intgrtrst(&lf_im(n)(0),&cv(n,0)(0)(0)(0),&cv(n,1)(0)(0)(0),&cv(n,2)(0)(0)(0),stridex,stridey);
		
		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			/* TIME DERIVATIVE TERMS */ 
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						
						double rho = u(0)(i)(j)(k)/u(NV-1)(i)(j)(k);
						cjcb = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
						
						rhorbd0 = rho*gbl->bd(0)*cjcb;
						mujcbi = lmu/cjcb;
						kcjcbi = lkcond/cjcb/gbl->R;
						
						/* UNSTEADY TERMS */
						res(0)(i)(j)(k) = rhorbd0+dugdt(log2p,tind,0)(i)(j)(k);
						for(n=1;n<NV-1;++n)
							res(n)(i)(j)(k) = rhorbd0*u(n)(i)(j)(k) +dugdt(log2p,tind,n)(i)(j)(k);
						double e = ogm1*u(NV-1)(i)(j)(k) +0.5*(u(1)(i)(j)(k)*u(1)(i)(j)(k)+u(2)(i)(j)(k)*u(2)(i)(j)(k)+u(3)(i)(j)(k)*u(3)(i)(j)(k));
						res(NV-1)(i)(j)(k) = rhorbd0*e +dugdt(log2p,tind,NV-1)(i)(j)(k);
						
#ifdef BODYFORCE
						for(n=1;n<NV-1;++n)
							res(n)(i)(j)(k) -= rho*cjcb*gbl->body(n-1);			
						
#endif                        
						/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
						/* INDICES ARE 1: EQUATION U V OR W, 2: VARIABLE (U V OR W), 3: EQ. DERIVATIVE (R S OR T) 4: VAR DERIVATIVE (R S OR T)*/
						visc(0,0)(0,0) = -mujcbi*(4./3.*d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1)+d(0)(2)*d(0)(2));
						visc(0,0)(0,1) = -mujcbi*(4./3.*d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1)+d(0)(2)*d(1)(2));
						visc(0,0)(0,2) = -mujcbi*(4./3.*d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1)+d(0)(2)*d(2)(2));
						visc(0,0)(1,0) = visc(0,0)(0,1);
						visc(0,0)(1,1) = -mujcbi*(4./3.*d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1)+d(1)(2)*d(1)(2));
						visc(0,0)(1,2) = -mujcbi*(4./3.*d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1)+d(1)(2)*d(2)(2));
						visc(0,0)(2,0) = visc(0,0)(0,2);
						visc(0,0)(2,1) = visc(0,0)(1,2);
						visc(0,0)(2,2) = -mujcbi*(4./3.*d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1)+d(2)(2)*d(2)(2));
						
						visc(0,1)(0,0) = -mujcbi*1./3.*d(0)(1)*d(0)(0);
						visc(0,1)(0,1) = -mujcbi*(d(0)(1)*d(1)(0)-2./3.*d(0)(0)*d(1)(1));
						visc(0,1)(0,2) = -mujcbi*(d(0)(1)*d(2)(0)-2./3.*d(0)(0)*d(2)(1));
						visc(0,1)(1,0) = -mujcbi*(d(1)(1)*d(0)(0)-2./3.*d(1)(0)*d(0)(1));
						visc(0,1)(1,1) = -mujcbi*1./3.*d(1)(1)*d(1)(0);
						visc(0,1)(1,2) = -mujcbi*(d(1)(1)*d(2)(0)-2./3.*d(1)(0)*d(2)(1));
						visc(0,1)(2,0) = -mujcbi*(d(2)(1)*d(0)(0)-2./3.*d(2)(0)*d(0)(1));
						visc(0,1)(2,1) = -mujcbi*(d(2)(1)*d(1)(0)-2./3.*d(2)(0)*d(1)(1));
						visc(0,1)(2,2) = -mujcbi*1./3.*d(2)(1)*d(2)(0);
						
						visc(0,2)(0,0) = -mujcbi*1./3.*d(0)(2)*d(0)(0);
						visc(0,2)(0,1) = -mujcbi*(d(0)(2)*d(1)(0)-2./3.*d(0)(0)*d(1)(2));
						visc(0,2)(0,2) = -mujcbi*(d(0)(2)*d(2)(0)-2./3.*d(0)(0)*d(2)(2));
						visc(0,2)(1,0) = -mujcbi*(d(1)(2)*d(0)(0)-2./3.*d(1)(0)*d(0)(2));
						visc(0,2)(1,1) = -mujcbi*1./3.*d(1)(2)*d(1)(0);
						visc(0,2)(1,2) = -mujcbi*(d(1)(2)*d(2)(0)-2./3.*d(1)(0)*d(2)(2));
						visc(0,2)(2,0) = -mujcbi*(d(2)(2)*d(0)(0)-2./3.*d(2)(0)*d(0)(2));
						visc(0,2)(2,1) = -mujcbi*(d(2)(2)*d(1)(0)-2./3.*d(2)(0)*d(1)(2));
						visc(0,2)(2,2) = -mujcbi*1./3.*d(2)(2)*d(2)(0);
						
						visc(1,0)(0,0) = visc(0,1)(0,0);
						visc(1,0)(0,1) = visc(0,1)(1,0);
						visc(1,0)(0,2) = visc(0,1)(2,0);
						visc(1,0)(1,0) = visc(0,1)(0,1);
						visc(1,0)(1,1) = visc(0,1)(1,1);
						visc(1,0)(1,2) = visc(0,1)(2,1);
						visc(1,0)(2,0) = visc(0,1)(0,2);
						visc(1,0)(2,1) = visc(0,1)(1,2);
						visc(1,0)(2,2) = visc(0,1)(2,2);
						
						visc(1,1)(0,0) = -mujcbi*(4./3.*d(0)(1)*d(0)(1)+d(0)(0)*d(0)(0)+d(0)(2)*d(0)(2));
						visc(1,1)(0,1) = -mujcbi*(4./3.*d(0)(1)*d(1)(1)+d(0)(0)*d(1)(0)+d(0)(2)*d(1)(2));
						visc(1,1)(0,2) = -mujcbi*(4./3.*d(0)(1)*d(2)(1)+d(0)(0)*d(2)(0)+d(0)(2)*d(2)(2));
						visc(1,1)(1,0) = visc(1,1)(0,1);
						visc(1,1)(1,1) = -mujcbi*(4./3.*d(1)(1)*d(1)(1)+d(1)(0)*d(1)(0)+d(1)(2)*d(1)(2));
						visc(1,1)(1,2) = -mujcbi*(4./3.*d(1)(1)*d(2)(1)+d(1)(0)*d(2)(0)+d(1)(2)*d(2)(2));
						visc(1,1)(2,0) = visc(1,1)(0,2);
						visc(1,1)(2,1) = visc(1,1)(1,2);
						visc(1,1)(2,2) = -mujcbi*(4./3.*d(2)(1)*d(2)(1)+d(2)(0)*d(2)(0)+d(2)(2)*d(2)(2));
						
						visc(1,2)(0,0) = -mujcbi*1./3.*d(0)(2)*d(0)(1);
						visc(1,2)(0,1) = -mujcbi*(d(0)(2)*d(1)(1)-2./3.*d(0)(1)*d(1)(2));
						visc(1,2)(0,2) = -mujcbi*(d(0)(2)*d(2)(1)-2./3.*d(0)(1)*d(2)(2));
						visc(1,2)(1,0) = -mujcbi*(d(1)(2)*d(0)(1)-2./3.*d(1)(1)*d(0)(2));
						visc(1,2)(1,1) = -mujcbi*1./3.*d(1)(2)*d(1)(1);
						visc(1,2)(1,2) = -mujcbi*(d(1)(2)*d(2)(1)-2./3.*d(1)(1)*d(2)(2));
						visc(1,2)(2,0) = -mujcbi*(d(2)(2)*d(0)(1)-2./3.*d(2)(1)*d(0)(2));
						visc(1,2)(2,1) = -mujcbi*(d(2)(2)*d(1)(1)-2./3.*d(2)(1)*d(1)(2));
						visc(1,2)(2,2) = -mujcbi*1./3.*d(2)(2)*d(2)(1);
						
						visc(2,0)(0,0) = visc(0,2)(0,0);
						visc(2,0)(0,1) = visc(0,2)(1,0);
						visc(2,0)(0,2) = visc(0,2)(2,0);
						visc(2,0)(1,0) = visc(0,2)(0,1);
						visc(2,0)(1,1) = visc(0,2)(1,1);
						visc(2,0)(1,2) = visc(0,2)(2,1);
						visc(2,0)(2,0) = visc(0,2)(0,2);
						visc(2,0)(2,1) = visc(0,2)(1,2);
						visc(2,0)(2,2) = visc(0,2)(2,2);
						
						visc(2,1)(0,0) = visc(1,2)(0,0);
						visc(2,1)(0,1) = visc(1,2)(1,0);
						visc(2,1)(0,2) = visc(1,2)(2,0);
						visc(2,1)(1,0) = visc(1,2)(0,1);
						visc(2,1)(1,1) = visc(1,2)(1,1);
						visc(2,1)(1,2) = visc(1,2)(2,1);
						visc(2,1)(2,0) = visc(1,2)(0,2);
						visc(2,1)(2,1) = visc(1,2)(1,2);
						visc(2,1)(2,2) = visc(1,2)(2,2);
						
						visc(2,2)(0,0) = -mujcbi*(4./3.*d(0)(2)*d(0)(2)+d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1));
						visc(2,2)(0,1) = -mujcbi*(4./3.*d(0)(2)*d(1)(2)+d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1));
						visc(2,2)(0,2) = -mujcbi*(4./3.*d(0)(2)*d(2)(2)+d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1));
						visc(2,2)(1,0) = visc(2,2)(0,1);
						visc(2,2)(1,1) = -mujcbi*(4./3.*d(1)(2)*d(1)(2)+d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1));
						visc(2,2)(1,2) = -mujcbi*(4./3.*d(1)(2)*d(2)(2)+d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1));
						visc(2,2)(2,0) = visc(2,2)(0,2);
						visc(2,2)(2,1) = visc(2,2)(1,2);
						visc(2,2)(2,2) = -mujcbi*(4./3.*d(2)(2)*d(2)(2)+d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1));
						
						/* HEAT DIFFUSION TENSOR */
						kcond(0)(0) = -cjcb*(d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1)+d(0)(2)*d(0)(2));
						kcond(1)(1) = -cjcb*(d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1)+d(1)(2)*d(1)(2));
						kcond(2)(2) = -cjcb*(d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1)+d(2)(2)*d(2)(2));
						kcond(0)(1) = -cjcb*(d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1)+d(0)(2)*d(1)(2));
						kcond(0)(2) = -cjcb*(d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1)+d(0)(2)*d(2)(2));
						kcond(1)(2) = -cjcb*(d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1)+d(1)(2)*d(2)(2));
						
#define					kcondI1II0I  kcond(0)(1)
#define					kcondI2II0I  kcond(0)(2)
#define					kcondI2II1I  kcond(1)(2)
						
						
						/*  MOMENTUM EQUATIONS */
						df(1,0)(i)(j)(k) = +visc(0,0)(0,0)*du(1,0)(i)(j)(k)+visc(0,1)(0,0)*du(2,0)(i)(j)(k)+visc(0,2)(0,0)*du(3,0)(i)(j)(k)
						+visc(0,0)(0,1)*du(1,1)(i)(j)(k)+visc(0,1)(0,1)*du(2,1)(i)(j)(k)+visc(0,2)(0,1)*du(3,1)(i)(j)(k)
						+visc(0,0)(0,2)*du(1,2)(i)(j)(k)+visc(0,1)(0,2)*du(2,2)(i)(j)(k)+visc(0,2)(0,2)*du(3,2)(i)(j)(k);
						
						df(1,1)(i)(j)(k) = +visc(0,0)(1,0)*du(1,0)(i)(j)(k)+visc(0,1)(1,0)*du(2,0)(i)(j)(k)+visc(0,2)(1,0)*du(3,0)(i)(j)(k)
						+visc(0,0)(1,1)*du(1,1)(i)(j)(k)+visc(0,1)(1,1)*du(2,1)(i)(j)(k)+visc(0,2)(1,1)*du(3,1)(i)(j)(k)
						+visc(0,0)(1,2)*du(1,2)(i)(j)(k)+visc(0,1)(1,2)*du(2,2)(i)(j)(k)+visc(0,2)(1,2)*du(3,2)(i)(j)(k);
						
						df(1,2)(i)(j)(k) = +visc(0,0)(2,0)*du(1,0)(i)(j)(k)+visc(0,1)(2,0)*du(2,0)(i)(j)(k)+visc(0,2)(2,0)*du(3,0)(i)(j)(k)
						+visc(0,0)(2,1)*du(1,1)(i)(j)(k)+visc(0,1)(2,1)*du(2,1)(i)(j)(k)+visc(0,2)(2,1)*du(3,1)(i)(j)(k)
						+visc(0,0)(2,2)*du(1,2)(i)(j)(k)+visc(0,1)(2,2)*du(2,2)(i)(j)(k)+visc(0,2)(2,2)*du(3,2)(i)(j)(k);
						
						df(2,0)(i)(j)(k) = +visc(1,0)(0,0)*du(1,0)(i)(j)(k)+visc(1,1)(0,0)*du(2,0)(i)(j)(k)+visc(1,2)(0,0)*du(3,0)(i)(j)(k)
						+visc(1,0)(0,1)*du(1,1)(i)(j)(k)+visc(1,1)(0,1)*du(2,1)(i)(j)(k)+visc(1,2)(0,1)*du(3,1)(i)(j)(k)
						+visc(1,0)(0,2)*du(1,2)(i)(j)(k)+visc(1,1)(0,2)*du(2,2)(i)(j)(k)+visc(1,2)(0,2)*du(3,2)(i)(j)(k);
						
						df(2,1)(i)(j)(k) = +visc(1,0)(1,0)*du(1,0)(i)(j)(k)+visc(1,1)(1,0)*du(2,0)(i)(j)(k)+visc(1,2)(1,0)*du(3,0)(i)(j)(k)
						+visc(1,0)(1,1)*du(1,1)(i)(j)(k)+visc(1,1)(1,1)*du(2,1)(i)(j)(k)+visc(1,2)(1,1)*du(3,1)(i)(j)(k)
						+visc(1,0)(1,2)*du(1,2)(i)(j)(k)+visc(1,1)(1,2)*du(2,2)(i)(j)(k)+visc(1,2)(1,2)*du(3,2)(i)(j)(k);
						
						df(2,2)(i)(j)(k) = +visc(1,0)(2,0)*du(1,0)(i)(j)(k)+visc(1,1)(2,0)*du(2,0)(i)(j)(k)+visc(1,2)(2,0)*du(3,0)(i)(j)(k)
						+visc(1,0)(2,1)*du(1,1)(i)(j)(k)+visc(1,1)(2,1)*du(2,1)(i)(j)(k)+visc(1,2)(2,1)*du(3,1)(i)(j)(k)
						+visc(1,0)(2,2)*du(1,2)(i)(j)(k)+visc(1,1)(2,2)*du(2,2)(i)(j)(k)+visc(1,2)(2,2)*du(3,2)(i)(j)(k);
						
						df(3,0)(i)(j)(k) = +visc(2,0)(0,0)*du(1,0)(i)(j)(k)+visc(2,1)(0,0)*du(2,0)(i)(j)(k)+visc(2,2)(0,0)*du(3,0)(i)(j)(k)
						+visc(2,0)(0,1)*du(1,1)(i)(j)(k)+visc(2,1)(0,1)*du(2,1)(i)(j)(k)+visc(2,2)(0,1)*du(3,1)(i)(j)(k)
						+visc(2,0)(0,2)*du(1,2)(i)(j)(k)+visc(2,1)(0,2)*du(2,2)(i)(j)(k)+visc(2,2)(0,2)*du(3,2)(i)(j)(k);
						
						df(3,1)(i)(j)(k) = +visc(2,0)(1,0)*du(1,0)(i)(j)(k)+visc(2,1)(1,0)*du(2,0)(i)(j)(k)+visc(2,2)(1,0)*du(3,0)(i)(j)(k)
						+visc(2,0)(1,1)*du(1,1)(i)(j)(k)+visc(2,1)(1,1)*du(2,1)(i)(j)(k)+visc(2,2)(1,1)*du(3,1)(i)(j)(k)
						+visc(2,0)(1,2)*du(1,2)(i)(j)(k)+visc(2,1)(1,2)*du(2,2)(i)(j)(k)+visc(2,2)(1,2)*du(3,2)(i)(j)(k);
						
						df(3,2)(i)(j)(k) = +visc(2,0)(2,0)*du(1,0)(i)(j)(k)+visc(2,1)(2,0)*du(2,0)(i)(j)(k)+visc(2,2)(2,0)*du(3,0)(i)(j)(k)
						+visc(2,0)(2,1)*du(1,1)(i)(j)(k)+visc(2,1)(2,1)*du(2,1)(i)(j)(k)+visc(2,2)(2,1)*du(3,1)(i)(j)(k)
						+visc(2,0)(2,2)*du(1,2)(i)(j)(k)+visc(2,1)(2,2)*du(2,2)(i)(j)(k)+visc(2,2)(2,2)*du(3,2)(i)(j)(k);
						
						
						/*  ENGERY EQUATION */
						df(4,0)(i)(j)(k) = kcond(0)(0)*du(NV-1,0)(i)(j)(k)+kcond(0)(1)*du(NV-1,1)(i)(j)(k)+kcond(0)(2)*du(NV-1,2)(i)(j)(k);
						df(4,1)(i)(j)(k) = kcondI1II0I*du(NV-1,0)(i)(j)(k)+kcond(1)(1)*du(NV-1,1)(i)(j)(k)+kcond(1)(2)*du(NV-1,2)(i)(j)(k);
						df(4,2)(i)(j)(k) = kcondI2II0I*du(NV-1,0)(i)(j)(k)+kcondI2II1I*du(NV-1,1)(i)(j)(k)+kcond(2)(2)*du(NV-1,2)(i)(j)(k);
						
						/* VISCOUS DISSIPATION */
						df(4,0)(i)(j)(k) += df(1,0)(i)(j)(k)*u(1)(i)(j)(k) +df(2,0)(i)(j)(k)*u(2)(i)(j)(k)+df(3,0)(i)(j)(k)*u(3)(i)(j)(k);
						df(4,1)(i)(j)(k) += df(1,1)(i)(j)(k)*u(1)(i)(j)(k) +df(2,1)(i)(j)(k)*u(2)(i)(j)(k)+df(3,1)(i)(j)(k)*u(3)(i)(j)(k);
						df(4,2)(i)(j)(k) += df(1,2)(i)(j)(k)*u(1)(i)(j)(k) +df(2,2)(i)(j)(k)*u(2)(i)(j)(k)+df(3,2)(i)(j)(k)*u(3)(i)(j)(k);
						
						for(n=1;n<NV;++n) {
							cv(n,0)(i)(j)(k) += df(n,0)(i)(j)(k);
							cv(n,1)(i)(j)(k) += df(n,1)(i)(j)(k);
							cv(n,2)(i)(j)(k) += df(n,2)(i)(j)(k);
						}
					}
				}
			}
			for(n=0;n<NV;++n)
				basis::tet(log2p).intgrt(&lf_re(n)(0),&res(n)(0)(0)(0),stridex,stridey);
			
			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(n=0;n<NV;++n) {
				basis::tet(log2p).derivr(&cv(n,0)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
				basis::tet(log2p).derivs(&cv(n,1)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
				basis::tet(log2p).derivt(&cv(n,2)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
			}
			
			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {

						df(0,0)(i)(j)(k) = 0.0;
						df(0,1)(i)(j)(k) = 0.0;
						
						tres = 0.0;
						for(int m = 0; m < NV; ++m)
							for(int n = 0; n < NV; ++n)							
								tres(m) += gbl->tau(tind,m,n)*res(n)(i)(j)(k);
						
						FLT pr = u(0)(i)(j)(k);
						FLT uv = u(1)(i)(j)(k);
						FLT vv = u(2)(i)(j)(k);
						FLT wv = u(3)(i)(j)(k);
						FLT rt = u(4)(i)(j)(k);					
						FLT ke = 0.5*(uv*uv+vv*vv+wv*wv);
						
						/* df/dw */
						A = 1.0/rt*uv,               pr/rt,                           0.0,         0.0,         -pr/rt/rt*uv,
						1.0/rt*uv*uv+1.0,        2.0*pr/rt*uv,                    0.0,         0.0,         -pr/rt/rt*uv*uv,
						1.0/rt*uv*vv,            pr/rt*vv,                        pr/rt*uv,    0.0,         -pr/rt/rt*uv*vv,
						1.0/rt*uv*wv,            pr/rt*wv,                        0.0,         pr/rt*uv,    -pr/rt/rt*uv*wv,
						1.0/rt*uv*(gogm1*rt+ke), pr/rt*(gogm1*rt+ke)+pr/rt*uv*uv, pr/rt*uv*vv, pr/rt*uv*wv, -pr/rt/rt*uv*(gogm1*rt+ke)+pr/rt*uv*gogm1;
						
						
						/* dg/dw */
						B = 1.0/rt*vv,               0.0,         pr/rt,                           0.0,         -pr/rt/rt*vv,
						1.0/rt*uv*vv,            pr/rt*vv,    pr/rt*uv,                        0.0,         -pr/rt/rt*uv*vv,
						1.0/rt*vv*vv+1.0,        0.0,         2.0*pr/rt*vv,                    0.0,         -pr/rt/rt*vv*vv,
						1.0/rt*vv*wv,            0.0,         pr/rt*wv,                        pr/rt*vv,    -pr/rt/rt*vv*wv,
						1.0/rt*vv*(gogm1*rt+ke), pr/rt*uv*vv, pr/rt*(gogm1*rt+ke)+pr/rt*vv*vv, pr/rt*vv*wv,	-pr/rt/rt*vv*(gogm1*rt+ke)+pr/rt*vv*gogm1;
						
						
						
						/* dh/dw */
						C = 1.0/rt*wv,               0.0,         0.0,         pr/rt,                           -pr/rt/rt*wv,
						1.0/rt*uv*wv,            pr/rt*wv,    0.0,         pr/rt*uv,                        -pr/rt/rt*uv*wv,
						1.0/rt*vv*wv,            0.0,         pr/rt*wv,    pr/rt*vv,                        -pr/rt/rt*vv*wv,
						1.0/rt*wv*wv+1.0,        0.0,         0.0,         2.0*pr/rt*wv,                    -pr/rt/rt*wv*wv,
						1.0/rt*wv*(gogm1*rt+ke), pr/rt*uv*wv, pr/rt*vv*wv, pr/rt*(gogm1*rt+ke)+pr/rt*wv*wv, -pr/rt/rt*wv*(gogm1*rt+ke)+pr/rt*wv*gogm1;				
						
						for(int m = 0; m < NV; ++m) {
							for(int n = 0; n < NV; ++n) {
								df(m,0)(i)(j)(k) -= (d(0)(0)*A(m,n)+d(0)(1)*B(m,n)+d(0)(2)*C(m,n))*tres(n);
								df(m,1)(i)(j)(k) -= (d(1)(0)*A(m,n)+d(1)(1)*B(m,n)+d(1)(2)*C(m,n))*tres(n);
								df(m,2)(i)(j)(k) -= (d(2)(0)*A(m,n)+d(2)(1)*B(m,n)+d(2)(2)*C(m,n))*tres(n);
							}
						}
					}
				}
			}
			
			for(n=0;n<NV;++n)
				basis::tet(log2p).intgrtrst(&lf_re(n)(0),&df(n,0)(0)(0)(0),&df(n,1)(0)(0)(0),&df(n,2)(0)(0)(0),stridex,stridey);
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
					lf_re(n)(i) *= gbl->beta(stage);
			
			
		}

	}
	
	
	return;
}