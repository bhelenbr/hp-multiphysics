/*
 *  global_matrix.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 9/10/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#include "tet_hp_ins.h"
#include "../hp_boundary.h"

// #define BODYFORCE

void tet_hp_ins::create_sparse_global_stiffness_matrix() {
	int i,j,m,n,indx,gindx,eind,find,sgn,msgn,mode,iind;
	
	int kn = MXTM*NV;
	TinyMatrix<FLT,kn,kn K;
	TinyVector<int,kn> loc_to_glo;
	int ind = 0;

	for(int tind = 0; tind < ntet; ++tind){
		
		create_local_stiffness_matrix(tind,&K(0,0));
		
		for (m = 0; m < 4; ++m) {
			gindx = NV*tet(tind).pnt(m);
			for (n = 0; n < NV; ++n)
				loc_to_glo(ind++) = gindx+n;
		}
		
		
		/* EDGE MODES */
		if (basis::tet(log2p).p > 1) {
			for(i = 0; i < 6; ++i) {
				eind = npnt*NV + tet(tind).seg(i)*basis::tet(log2p).em*NV;
				sgn = tet(tind).sgn(i);
				msgn = 1;
				for (m = 0; m < basis::tet(log2p).em; ++m) {
					for(n = 0; n < NV; ++n) {
						for(j = 0; j < kn; ++j)
							K(ind,j) = msgn*K(ind,j);
						loc_to_glo(ind++) = eind + m*basis::tet(log2p).em + n;
					}
					msgn *= sgn;
					++indx;
				}
			}
		}
		
		/* FACE MODES */
		if (basis::tet(log2p).p > 2) {
			for(i = 0; i < 4; ++i){
				sgn = -tet(tind).rot(i);
				find = npnt*NV+nseg*basis::tet(log2p).em*NV+tet(tind).tri(i)*basis::tet(log2p).em*NV;
				mode = 0;
				msgn = 1;		
				for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
					for(int j = 1; j <= basis::tet(log2p).em-m; ++j){
						for(n = 0; n < NV; ++n){
							for(k = 0; k < kn; ++k)
								K(ind,k) = msgn*K(ind,k);
							loc_to_glo(ind++) = find+mode*basis::tet(log2p).fm+n;
						}
						++mode;
					}
					msgn *= sgn;
				}
			}		
		}		
		
		/* INTERIOR MODES */
		mode = 0;
		iind = npnt*NV+nseg*basis::tet(log2p).em*NV+ntri*basis::tet(log2p).em*NV
		for(m = 0; m < basis::tet(log2p).im; ++m) {
			for(n = 0; n < NV; ++n){
				for(k = 0; k < kn; ++k)
					K(ind,k) = msgn*K(ind,k);
				loc_to_glo(ind++) = iind+mode*basis::tet(log2p).im+n;
			}
			++mode;
		}
		
		for(i = 0; i < kn; ++i)
			for(j = 0; j < kn; ++j)
				insert_sparse(loc_to_glo(i), loc_to_gbl(j), K(i,j))

		
		
	}
	return;
}

void tet_hp_ins::create_local_stiffness_matrix(int tind, FLT *K) {
	
	TinyMatrix<FLT,NV,MXTM> R,Rbar;
	int kcol = 0;
	FLT dw = 0.1;
	
	Rbar = 0;
	tet_hp_ins::residual(tind, &Rbar(0,0));
	
	for(int mode = 0; mode < MXTM; ++mode){
		for(int var = 0; variable < NV; ++var){
			R = 0;
			R(var)(mode) = dw;
			tet_hp_ins::residual(tind, &R(0,0));
			for(i=0;i<basis::tet(log2p).tm;++i)
				for(n=0;n<NV;++n)
					K(++krow,kcol) = (R(n,i)-Rbar(n,i))/dw;
			++kcol;
		}
	}	
	
	return;
}

void tet_hp_ins::residual(int tind, FLT *dw){
	int i,j,k,n,tind;
	FLT fluxx,fluxy,fluxz;
	const int NV = 4;
	TinyVector<int,4> v;
	TinyMatrix<FLT,ND,ND> ldcrd;
	TinyMatrix<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,NV,ND> du;
	int lgpx = basis::tet(log2p).gpx, lgpy = basis::tet(log2p).gpy, lgpz = basis::tet(log2p).gpz;
	int stridey = MXGP;
	int stridex = MXGP*MXGP; 
	FLT rhobd0 = gbl->rho*gbl->bd(0), lmu = gbl->mu, rhorbd0, cjcb, cjcbi, oneminusbeta;
	TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-1,NV-1> visc;
	TinyVector<TinyVector<FLT,ND>,ND> d;
	TinyMatrix<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,NV-1,NV-1> cv, df;
	TinyVector<FLT,NV> tres;
	
	tet_hp::rsdl(stage);
	
	oneminusbeta = 1.0-gbl->beta(stage);
	
	// temporary for now 
	int stage = 0;
	
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
	
	/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
	/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
	ugtouht(tind);
	
	/* pertube by dw */
	for(n=0;n<NV;++n)
		for(i=0;i<basis::tet(log2p).tm;++i)
			uht(n)(i) += dw(n,i);
	
	if (gbl->beta(stage) > 0.0) {
		for(n=0;n<NV-1;++n)
			basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),&du(n,0)(0)(0)(0),&du(n,1)(0)(0)(0),&du(n,2)(0)(0)(0),stridex,stridey);
		basis::tet(log2p).proj(&uht(NV-1)(0),&u(NV-1)(0)(0)(0),stridex,stridey);
		
	}
	else {
		for(n=0;n<NV;++n)
			basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);
	}
	
	/* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
	for(n=0;n<NV;++n)
		for(i=0;i<basis::tet(log2p).tm;++i)
			lf(n)(i) = 0.0;
	
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
					
					fluxx = gbl->rho*(u(0)(i)(j)(k) -mvel(0)(i)(j)(k));
					fluxy = gbl->rho*(u(1)(i)(j)(k) -mvel(1)(i)(j)(k));
					fluxz = gbl->rho*(u(2)(i)(j)(k) -mvel(2)(i)(j)(k));
					
					/* CONTINUITY EQUATION FLUXES */
					du(NV-1,0)(i)(j)(k) = d(0)(0)*fluxx+d(0)(1)*fluxy+d(0)(2)*fluxz;
					du(NV-1,1)(i)(j)(k) = d(1)(0)*fluxx+d(1)(1)*fluxy+d(1)(2)*fluxz;					
					du(NV-1,2)(i)(j)(k) = d(2)(0)*fluxx+d(2)(1)*fluxy+d(2)(2)*fluxz;
					
					
#ifndef INERTIALESS
					/* CONVECTIVE FLUXES */
					for(n=0;n<NV-1;++n) {
						cv(n,0)(i)(j)(k) = u(n)(i)(j)(k)*du(NV-1,0)(i)(j)(k);
						cv(n,1)(i)(j)(k) = u(n)(i)(j)(k)*du(NV-1,1)(i)(j)(k);
						cv(n,2)(i)(j)(k) = u(n)(i)(j)(k)*du(NV-1,2)(i)(j)(k);
					}
#else
					for(n=0;n<NV-1;++n) {
						cv(n,0)(i)(j)(k) = 0.0;
						cv(n,1)(i)(j)(k) = 0.0;
						cv(n,2)(i)(j)(k) = 0.0;
					}
#endif
					
					/* PRESSURE TERMS */
					/* U-MOMENTUM */
					cv(0,0)(i)(j)(k) += d(0)(0)*u(NV-1)(i)(j)(k);
					cv(0,1)(i)(j)(k) += d(1)(0)*u(NV-1)(i)(j)(k);
					cv(0,2)(i)(j)(k) += d(2)(0)*u(NV-1)(i)(j)(k);
					
					/* V-MOMENTUM */
					cv(1,0)(i)(j)(k) += d(0)(1)*u(NV-1)(i)(j)(k);
					cv(1,1)(i)(j)(k) += d(1)(1)*u(NV-1)(i)(j)(k);
					cv(1,2)(i)(j)(k) += d(2)(1)*u(NV-1)(i)(j)(k);
					
					/* W-MOMENTUM */
					cv(2,0)(i)(j)(k) += d(0)(2)*u(NV-1)(i)(j)(k);
					cv(2,1)(i)(j)(k) += d(1)(2)*u(NV-1)(i)(j)(k);
					cv(2,2)(i)(j)(k) += d(2)(2)*u(NV-1)(i)(j)(k);
				}
			}
		}
		for(n=0;n<NV-1;++n)
			basis::tet(log2p).intgrtrst(&lf(n)(0),&cv(n,0)(0)(0)(0),&cv(n,1)(0)(0)(0),&cv(n,2)(0)(0)(0),stridex,stridey);
		basis::tet(log2p).intgrtrst(&lf(NV-1)(0),&du(NV-1,0)(0)(0)(0),&du(NV-1,1)(0)(0)(0),&du(NV-1,2)(0)(0)(0),stridex,stridey);
		
		/* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
		for(n=0;n<NV;++n)
			for(i=0;i<basis::tet(log2p).tm;++i)
				dw(n,i)=lf(n)(i);
		
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
						
						cjcb = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
						
						rhorbd0 = rhobd0*cjcb;
						cjcbi = lmu/cjcb;
						
						/* UNSTEADY TERMS */
						for(n=0;n<NV-1;++n)
							res(n)(i)(j)(k) = rhorbd0*u(n)(i)(j)(k) +dugdt(log2p,tind,n)(i)(j)(k);
						res(NV-1)(i)(j)(k) = rhorbd0 +dugdt(log2p,tind,NV-1)(i)(j)(k);
						
#ifdef INERTIALESS
						for(n=0;n<NV-1;++n)
							res(n)(i)(j)(k) = 0.0;	
#endif
						
#ifdef BODYFORCE
						for(n=0;n<NV-1;++n)
							res(n)(i)(j)(k) -= gbl->rho*cjcb*gbl->body(n);			
						
#endif                        
						/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
						/* INDICES ARE 1: EQUATION U V OR W, 2: VARIABLE (U V OR W), 3: EQ. DERIVATIVE (R S OR T) 4: VAR DERIVATIVE (R S OR T)*/
						visc(0,0)(0,0) = -cjcbi*(2.*d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1)+d(0)(2)*d(0)(2));
						visc(0,0)(0,1) = -cjcbi*(2.*d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1)+d(0)(2)*d(1)(2));
						visc(0,0)(0,2) = -cjcbi*(2.*d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1)+d(0)(2)*d(2)(2));
						visc(0,0)(1,0) = visc(0,0)(0,1);
						visc(0,0)(1,1) = -cjcbi*(2.*d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1)+d(1)(2)*d(1)(2));
						visc(0,0)(1,2) = -cjcbi*(2.*d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1)+d(1)(2)*d(2)(2));
						visc(0,0)(2,0) = visc(0,0)(0,2);
						visc(0,0)(2,1) = visc(0,0)(1,2);
						visc(0,0)(2,2) = -cjcbi*(2.*d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1)+d(2)(2)*d(2)(2));
						
						visc(0,1)(0,0) = -cjcbi*d(0)(1)*d(0)(0);
						visc(0,1)(0,1) = -cjcbi*d(0)(1)*d(1)(0);
						visc(0,1)(0,2) = -cjcbi*d(0)(1)*d(2)(0);
						visc(0,1)(1,0) = -cjcbi*d(1)(1)*d(0)(0);
						visc(0,1)(1,1) = -cjcbi*d(1)(1)*d(1)(0);
						visc(0,1)(1,2) = -cjcbi*d(1)(1)*d(2)(0);
						visc(0,1)(2,0) = -cjcbi*d(2)(1)*d(0)(0);
						visc(0,1)(2,1) = -cjcbi*d(2)(1)*d(1)(0);
						visc(0,1)(2,2) = -cjcbi*d(2)(1)*d(2)(0);
						
						visc(0,2)(0,0) = -cjcbi*d(0)(2)*d(0)(0);
						visc(0,2)(0,1) = -cjcbi*d(0)(2)*d(1)(0);
						visc(0,2)(0,2) = -cjcbi*d(0)(2)*d(2)(0);
						visc(0,2)(1,0) = -cjcbi*d(1)(2)*d(0)(0);
						visc(0,2)(1,1) = -cjcbi*d(1)(2)*d(1)(0);
						visc(0,2)(1,2) = -cjcbi*d(1)(2)*d(2)(0);
						visc(0,2)(2,0) = -cjcbi*d(2)(2)*d(0)(0);
						visc(0,2)(2,1) = -cjcbi*d(2)(2)*d(1)(0);
						visc(0,2)(2,2) = -cjcbi*d(2)(2)*d(2)(0);
						
						visc(1,0)(0,0) = visc(0,1)(0,0);
						visc(1,0)(0,1) = visc(0,1)(1,0);
						visc(1,0)(0,2) = visc(0,1)(2,0);
						visc(1,0)(1,0) = visc(0,1)(0,1);
						visc(1,0)(1,1) = visc(0,1)(1,1);
						visc(1,0)(1,2) = visc(0,1)(2,1);
						visc(1,0)(2,0) = visc(0,1)(0,2);
						visc(1,0)(2,1) = visc(0,1)(1,2);
						visc(1,0)(2,2) = visc(0,1)(2,2);
						
						visc(1,1)(0,0) = -cjcbi*(2.*d(0)(1)*d(0)(1)+d(0)(0)*d(0)(0)+d(0)(2)*d(0)(2));
						visc(1,1)(0,1) = -cjcbi*(2.*d(0)(1)*d(1)(1)+d(0)(0)*d(1)(0)+d(0)(2)*d(1)(2));
						visc(1,1)(0,2) = -cjcbi*(2.*d(0)(1)*d(2)(1)+d(0)(0)*d(2)(0)+d(0)(2)*d(2)(2));
						visc(1,1)(1,0) = visc(1,1)(0,1);
						visc(1,1)(1,1) = -cjcbi*(2.*d(1)(1)*d(1)(1)+d(1)(0)*d(1)(0)+d(1)(2)*d(1)(2));
						visc(1,1)(1,2) = -cjcbi*(2.*d(1)(1)*d(2)(1)+d(1)(0)*d(2)(0)+d(1)(2)*d(2)(2));
						visc(1,1)(2,0) = visc(1,1)(0,2);
						visc(1,1)(2,1) = visc(1,1)(1,2);
						visc(1,1)(2,2) = -cjcbi*(2.*d(2)(1)*d(2)(1)+d(2)(0)*d(2)(0)+d(2)(2)*d(2)(2));
						
						visc(1,2)(0,0) = -cjcbi*d(0)(2)*d(0)(1);
						visc(1,2)(0,1) = -cjcbi*d(0)(2)*d(1)(1);
						visc(1,2)(0,2) = -cjcbi*d(0)(2)*d(2)(1);
						visc(1,2)(1,0) = -cjcbi*d(1)(2)*d(0)(1);
						visc(1,2)(1,1) = -cjcbi*d(1)(2)*d(1)(1);
						visc(1,2)(1,2) = -cjcbi*d(1)(2)*d(2)(1);
						visc(1,2)(2,0) = -cjcbi*d(2)(2)*d(0)(1);
						visc(1,2)(2,1) = -cjcbi*d(2)(2)*d(1)(1);
						visc(1,2)(2,2) = -cjcbi*d(2)(2)*d(2)(1);
						
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
						
						visc(2,2)(0,0) = -cjcbi*(2.*d(0)(2)*d(0)(2)+d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1));
						visc(2,2)(0,1) = -cjcbi*(2.*d(0)(2)*d(1)(2)+d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1));
						visc(2,2)(0,2) = -cjcbi*(2.*d(0)(2)*d(2)(2)+d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1));
						visc(2,2)(1,0) = visc(2,2)(0,1);
						visc(2,2)(1,1) = -cjcbi*(2.*d(1)(2)*d(1)(2)+d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1));
						visc(2,2)(1,2) = -cjcbi*(2.*d(1)(2)*d(2)(2)+d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1));
						visc(2,2)(2,0) = visc(2,2)(0,2);
						visc(2,2)(2,1) = visc(2,2)(1,2);
						visc(2,2)(2,2) = -cjcbi*(2.*d(2)(2)*d(2)(2)+d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1));
						
						
						df(0,0)(i)(j)(k) = +visc(0,0)(0,0)*du(0,0)(i)(j)(k)+visc(0,1)(0,0)*du(1,0)(i)(j)(k)+visc(0,2)(0,0)*du(2,0)(i)(j)(k)
						+visc(0,0)(0,1)*du(0,1)(i)(j)(k)+visc(0,1)(0,1)*du(1,1)(i)(j)(k)+visc(0,2)(0,1)*du(2,1)(i)(j)(k)
						+visc(0,0)(0,2)*du(0,2)(i)(j)(k)+visc(0,1)(0,2)*du(1,2)(i)(j)(k)+visc(0,2)(0,2)*du(2,2)(i)(j)(k);
						
						df(0,1)(i)(j)(k) = +visc(0,0)(1,0)*du(0,0)(i)(j)(k)+visc(0,1)(1,0)*du(1,0)(i)(j)(k)+visc(0,2)(1,0)*du(2,0)(i)(j)(k)
						+visc(0,0)(1,1)*du(0,1)(i)(j)(k)+visc(0,1)(1,1)*du(1,1)(i)(j)(k)+visc(0,2)(1,1)*du(2,1)(i)(j)(k)
						+visc(0,0)(1,2)*du(0,2)(i)(j)(k)+visc(0,1)(1,2)*du(1,2)(i)(j)(k)+visc(0,2)(1,2)*du(2,2)(i)(j)(k);
						
						df(0,2)(i)(j)(k) = +visc(0,0)(2,0)*du(0,0)(i)(j)(k)+visc(0,1)(2,0)*du(1,0)(i)(j)(k)+visc(0,2)(2,0)*du(2,0)(i)(j)(k)
						+visc(0,0)(2,1)*du(0,1)(i)(j)(k)+visc(0,1)(2,1)*du(1,1)(i)(j)(k)+visc(0,2)(2,1)*du(2,1)(i)(j)(k)
						+visc(0,0)(2,2)*du(0,2)(i)(j)(k)+visc(0,1)(2,2)*du(1,2)(i)(j)(k)+visc(0,2)(2,2)*du(2,2)(i)(j)(k);
						
						df(1,0)(i)(j)(k) = +visc(1,0)(0,0)*du(0,0)(i)(j)(k)+visc(1,1)(0,0)*du(1,0)(i)(j)(k)+visc(1,2)(0,0)*du(2,0)(i)(j)(k)
						+visc(1,0)(0,1)*du(0,1)(i)(j)(k)+visc(1,1)(0,1)*du(1,1)(i)(j)(k)+visc(1,2)(0,1)*du(2,1)(i)(j)(k)
						+visc(1,0)(0,2)*du(0,2)(i)(j)(k)+visc(1,1)(0,2)*du(1,2)(i)(j)(k)+visc(1,2)(0,2)*du(2,2)(i)(j)(k);
						
						df(1,1)(i)(j)(k) = +visc(1,0)(1,0)*du(0,0)(i)(j)(k)+visc(1,1)(1,0)*du(1,0)(i)(j)(k)+visc(1,2)(1,0)*du(2,0)(i)(j)(k)
						+visc(1,0)(1,1)*du(0,1)(i)(j)(k)+visc(1,1)(1,1)*du(1,1)(i)(j)(k)+visc(1,2)(1,1)*du(2,1)(i)(j)(k)
						+visc(1,0)(1,2)*du(0,2)(i)(j)(k)+visc(1,1)(1,2)*du(1,2)(i)(j)(k)+visc(1,2)(1,2)*du(2,2)(i)(j)(k);
						
						df(1,2)(i)(j)(k) = +visc(1,0)(2,0)*du(0,0)(i)(j)(k)+visc(1,1)(2,0)*du(1,0)(i)(j)(k)+visc(1,2)(2,0)*du(2,0)(i)(j)(k)
						+visc(1,0)(2,1)*du(0,1)(i)(j)(k)+visc(1,1)(2,1)*du(1,1)(i)(j)(k)+visc(1,2)(2,1)*du(2,1)(i)(j)(k)
						+visc(1,0)(2,2)*du(0,2)(i)(j)(k)+visc(1,1)(2,2)*du(1,2)(i)(j)(k)+visc(1,2)(2,2)*du(2,2)(i)(j)(k);
						
						df(2,0)(i)(j)(k) = +visc(2,0)(0,0)*du(0,0)(i)(j)(k)+visc(2,1)(0,0)*du(1,0)(i)(j)(k)+visc(2,2)(0,0)*du(2,0)(i)(j)(k)
						+visc(2,0)(0,1)*du(0,1)(i)(j)(k)+visc(2,1)(0,1)*du(1,1)(i)(j)(k)+visc(2,2)(0,1)*du(2,1)(i)(j)(k)
						+visc(2,0)(0,2)*du(0,2)(i)(j)(k)+visc(2,1)(0,2)*du(1,2)(i)(j)(k)+visc(2,2)(0,2)*du(2,2)(i)(j)(k);
						
						df(2,1)(i)(j)(k) = +visc(2,0)(1,0)*du(0,0)(i)(j)(k)+visc(2,1)(1,0)*du(1,0)(i)(j)(k)+visc(2,2)(1,0)*du(2,0)(i)(j)(k)
						+visc(2,0)(1,1)*du(0,1)(i)(j)(k)+visc(2,1)(1,1)*du(1,1)(i)(j)(k)+visc(2,2)(1,1)*du(2,1)(i)(j)(k)
						+visc(2,0)(1,2)*du(0,2)(i)(j)(k)+visc(2,1)(1,2)*du(1,2)(i)(j)(k)+visc(2,2)(1,2)*du(2,2)(i)(j)(k);
						
						df(2,2)(i)(j)(k) = +visc(2,0)(2,0)*du(0,0)(i)(j)(k)+visc(2,1)(2,0)*du(1,0)(i)(j)(k)+visc(2,2)(2,0)*du(2,0)(i)(j)(k)
						+visc(2,0)(2,1)*du(0,1)(i)(j)(k)+visc(2,1)(2,1)*du(1,1)(i)(j)(k)+visc(2,2)(2,1)*du(2,1)(i)(j)(k)
						+visc(2,0)(2,2)*du(0,2)(i)(j)(k)+visc(2,1)(2,2)*du(1,2)(i)(j)(k)+visc(2,2)(2,2)*du(2,2)(i)(j)(k);
						
						for(n=0;n<NV-1;++n) {
							cv(n,0)(i)(j)(k) += df(n,0)(i)(j)(k);
							cv(n,1)(i)(j)(k) += df(n,1)(i)(j)(k);
							cv(n,2)(i)(j)(k) += df(n,2)(i)(j)(k);
						}
					}
				}
			}
			for(n=0;n<NV;++n)
				basis::tet(log2p).intgrt(&lf(n)(0),&res(n)(0)(0)(0),stridex,stridey);
			
			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(n=0;n<NV-1;++n) {
				basis::tet(log2p).derivr(&cv(n,0)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
				basis::tet(log2p).derivs(&cv(n,1)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
				basis::tet(log2p).derivt(&cv(n,2)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
			}
			basis::tet(log2p).derivr(&du(NV-1,0)(0)(0)(0),&res(NV-1)(0)(0)(0),stridex,stridey);
			basis::tet(log2p).derivs(&du(NV-1,1)(0)(0)(0),&res(NV-1)(0)(0)(0),stridex,stridey);
			basis::tet(log2p).derivt(&du(NV-1,2)(0)(0)(0),&res(NV-1)(0)(0)(0),stridex,stridey);
			
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
						
						tres(0) = gbl->tau(tind,0)*res(0)(i)(j)(k);    
						tres(1) = gbl->tau(tind,0)*res(1)(i)(j)(k);
						tres(2) = gbl->tau(tind,0)*res(2)(i)(j)(k);
						tres(3) = gbl->tau(tind,NV-1)*res(NV-1)(i)(j)(k);
						
#ifndef INERTIALESS
						df(0,0)(i)(j)(k) -= (d(0)(0)*(2.*u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
											 +d(0)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
											 +d(0)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(0)
						+d(0)(1)*u(0)(i)(j)(k)*tres(1)
						+d(0)(2)*u(0)(i)(j)(k)*tres(2)
						+d(0)(0)*tres(NV-1);
						df(0,1)(i)(j)(k) -= (d(1)(0)*(2.*u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
											 +d(1)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
											 +d(1)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(0)
						+d(1)(1)*u(0)(i)(j)(k)*tres(1)
						+d(1)(2)*u(0)(i)(j)(k)*tres(2)
						+d(1)(0)*tres(NV-1);
						df(0,2)(i)(j)(k) -= (d(2)(0)*(2.*u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
											 +d(2)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
											 +d(2)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(0)
						+d(2)(1)*u(0)(i)(j)(k)*tres(1)
						+d(2)(2)*u(0)(i)(j)(k)*tres(2)
						+d(2)(0)*tres(NV-1);
						df(1,0)(i)(j)(k) -= d(0)(0)*u(1)(i)(j)(k)*tres(0)
						+(d(0)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(0)(1)*(2.*u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(0)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(1)
						+d(0)(2)*u(1)(i)(j)(k)*tres(2)				
						+d(0)(1)*tres(NV-1);		
						df(1,1)(i)(j)(k) -= d(1)(0)*u(1)(i)(j)(k)*tres(0)
						+(d(1)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(1)(1)*(2.*u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(1)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(1)
						+d(1)(2)*u(1)(i)(j)(k)*tres(2)				
						+d(1)(1)*tres(NV-1);		
						df(1,2)(i)(j)(k) -= d(2)(0)*u(1)(i)(j)(k)*tres(0)
						+(d(2)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(2)(1)*(2.*u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(2)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(1)
						+d(2)(2)*u(1)(i)(j)(k)*tres(2)				
						+d(2)(1)*tres(NV-1);
						df(2,0)(i)(j)(k) -= d(0)(0)*u(2)(i)(j)(k)*tres(0)
						+d(0)(1)*u(2)(i)(j)(k)*tres(1)
						+(d(0)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(0)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(0)(2)*(2.*u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(2)
						+d(0)(2)*tres(NV-1);	
						df(2,1)(i)(j)(k) -= d(1)(0)*u(2)(i)(j)(k)*tres(0)
						+d(1)(1)*u(2)(i)(j)(k)*tres(1)
						+(d(1)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(1)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(1)(2)*(2.*u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(2)
						+d(1)(2)*tres(NV-1);	
						df(2,2)(i)(j)(k) -= d(2)(0)*u(2)(i)(j)(k)*tres(0)
						+d(2)(1)*u(2)(i)(j)(k)*tres(1)
						+(d(2)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(2)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(2)(2)*(2.*u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(2)
						+d(2)(2)*tres(NV-1);
#endif
						
						du(NV-1,0)(i)(j)(k) = -d(0)(0)*tres(0)-d(0)(1)*tres(1)-d(0)(2)*tres(2);
						du(NV-1,1)(i)(j)(k) = -d(1)(0)*tres(0)-d(1)(1)*tres(1)-d(1)(2)*tres(2);					
						du(NV-1,2)(i)(j)(k) = -d(2)(0)*tres(0)-d(2)(1)*tres(1)-d(2)(2)*tres(2);
					}
				}
			}
			for(n=0;n<NV-1;++n)
				basis::tet(log2p).intgrtrst(&lf(n)(0),&df(n,0)(0)(0)(0),&df(n,1)(0)(0)(0),&df(n,2)(0)(0)(0),stridex,stridey);
			basis::tet(log2p).intgrtrst(&lf(NV-1)(0),&du(NV-1,0)(0)(0)(0),&du(NV-1,1)(0)(0)(0),&du(NV-1,2)(0)(0)(0),stridex,stridey);
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
					lf(n)(i) *= gbl->beta(stage);
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
					dw(n,i)+=lf(n)(i);
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
					
					fluxx = gbl->rho*(u(0)(i)(j)(k) -mvel(0)(i)(j)(k));
					fluxy = gbl->rho*(u(1)(i)(j)(k) -mvel(1)(i)(j)(k));
					fluxz = gbl->rho*(u(2)(i)(j)(k) -mvel(2)(i)(j)(k));
					
					/* CONTINUITY EQUATION FLUXES */
					du(NV-1,0)(i)(j)(k) = d(0)(0)*fluxx+d(0)(1)*fluxy+d(0)(2)*fluxz;
					du(NV-1,1)(i)(j)(k) = d(1)(0)*fluxx+d(1)(1)*fluxy+d(1)(2)*fluxz;					
					du(NV-1,2)(i)(j)(k) = d(2)(0)*fluxx+d(2)(1)*fluxy+d(2)(2)*fluxz;
					
					
#ifndef INERTIALESS
					/* CONVECTIVE FLUXES */
					for(n=0;n<NV-1;++n) {
						cv(n,0)(i)(j)(k) = u(n)(i)(j)(k)*du(NV-1,0)(i)(j)(k);
						cv(n,1)(i)(j)(k) = u(n)(i)(j)(k)*du(NV-1,1)(i)(j)(k);
						cv(n,2)(i)(j)(k) = u(n)(i)(j)(k)*du(NV-1,2)(i)(j)(k);
					}
#else
					for(n=0;n<NV-1;++n) {
						cv(n,0)(i)(j)(k) = 0.0;
						cv(n,1)(i)(j)(k) = 0.0;
						cv(n,2)(i)(j)(k) = 0.0;
					}
#endif
					
					/* PRESSURE TERMS */
					/* U-MOMENTUM */
					cv(0,0)(i)(j)(k) += d(0)(0)*u(NV-1)(i)(j)(k);
					cv(0,1)(i)(j)(k) += d(1)(0)*u(NV-1)(i)(j)(k);
					cv(0,2)(i)(j)(k) += d(2)(0)*u(NV-1)(i)(j)(k);
					
					/* V-MOMENTUM */
					cv(1,0)(i)(j)(k) += d(0)(1)*u(NV-1)(i)(j)(k);
					cv(1,1)(i)(j)(k) += d(1)(1)*u(NV-1)(i)(j)(k);
					cv(1,2)(i)(j)(k) += d(2)(1)*u(NV-1)(i)(j)(k);
					
					/* W-MOMENTUM */
					cv(2,0)(i)(j)(k) += d(0)(2)*u(NV-1)(i)(j)(k);
					cv(2,1)(i)(j)(k) += d(1)(2)*u(NV-1)(i)(j)(k);
					cv(2,2)(i)(j)(k) += d(2)(2)*u(NV-1)(i)(j)(k);		
				}
			}
		}
		for(n=0;n<NV-1;++n)
			basis::tet(log2p).intgrtrst(&lf(n)(0),&cv(n,0)(0)(0)(0),&cv(n,1)(0)(0)(0),&cv(n,2)(0)(0)(0),stridex,stridey);
		basis::tet(log2p).intgrtrst(&lf(NV-1)(0),&du(NV-1,0)(0)(0)(0),&du(NV-1,1)(0)(0)(0),&du(NV-1,2)(0)(0)(0),stridex,stridey);
		
		/* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
//		lftog(tind,gbl->res);
		for(n=0;n<NV;++n)
			for(i=0;i<basis::tet(log2p).tm;++i)
				dw(n,i)=lf(n)(i);
		
		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			cjcbi = lmu/cjcb;
			rhorbd0 = rhobd0*cjcb;
			
			/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
			/* INDICES ARE 1: EQUATION U V OR W, 2: VARIABLE (U V OR W), 3: EQ. DERIVATIVE (R S OR T) 4: VAR DERIVATIVE (R S OR T)*/
			visc(0,0)(0,0) = -cjcbi*(2.*d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1)+d(0)(2)*d(0)(2));
			visc(0,0)(0,1) = -cjcbi*(2.*d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1)+d(0)(2)*d(1)(2));
			visc(0,0)(0,2) = -cjcbi*(2.*d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1)+d(0)(2)*d(2)(2));
			visc(0,0)(1,0) = visc(0,0)(0,1);
			visc(0,0)(1,1) = -cjcbi*(2.*d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1)+d(1)(2)*d(1)(2));
			visc(0,0)(1,2) = -cjcbi*(2.*d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1)+d(1)(2)*d(2)(2));
			visc(0,0)(2,0) = visc(0,0)(0,2);
			visc(0,0)(2,1) = visc(0,0)(1,2);
			visc(0,0)(2,2) = -cjcbi*(2.*d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1)+d(2)(2)*d(2)(2));
			
			visc(0,1)(0,0) = -cjcbi*d(0)(1)*d(0)(0);
			visc(0,1)(0,1) = -cjcbi*d(0)(1)*d(1)(0);
			visc(0,1)(0,2) = -cjcbi*d(0)(1)*d(2)(0);
			visc(0,1)(1,0) = -cjcbi*d(1)(1)*d(0)(0);
			visc(0,1)(1,1) = -cjcbi*d(1)(1)*d(1)(0);
			visc(0,1)(1,2) = -cjcbi*d(1)(1)*d(2)(0);
			visc(0,1)(2,0) = -cjcbi*d(2)(1)*d(0)(0);
			visc(0,1)(2,1) = -cjcbi*d(2)(1)*d(1)(0);
			visc(0,1)(2,2) = -cjcbi*d(2)(1)*d(2)(0);
			
			visc(0,2)(0,0) = -cjcbi*d(0)(2)*d(0)(0);
			visc(0,2)(0,1) = -cjcbi*d(0)(2)*d(1)(0);
			visc(0,2)(0,2) = -cjcbi*d(0)(2)*d(2)(0);
			visc(0,2)(1,0) = -cjcbi*d(1)(2)*d(0)(0);
			visc(0,2)(1,1) = -cjcbi*d(1)(2)*d(1)(0);
			visc(0,2)(1,2) = -cjcbi*d(1)(2)*d(2)(0);
			visc(0,2)(2,0) = -cjcbi*d(2)(2)*d(0)(0);
			visc(0,2)(2,1) = -cjcbi*d(2)(2)*d(1)(0);
			visc(0,2)(2,2) = -cjcbi*d(2)(2)*d(2)(0);
			
			visc(1,0)(0,0) = visc(0,1)(0,0);
			visc(1,0)(0,1) = visc(0,1)(1,0);
			visc(1,0)(0,2) = visc(0,1)(2,0);
			visc(1,0)(1,0) = visc(0,1)(0,1);
			visc(1,0)(1,1) = visc(0,1)(1,1);
			visc(1,0)(1,2) = visc(0,1)(2,1);
			visc(1,0)(2,0) = visc(0,1)(0,2);
			visc(1,0)(2,1) = visc(0,1)(1,2);
			visc(1,0)(2,2) = visc(0,1)(2,2);
			
			visc(1,1)(0,0) = -cjcbi*(2.*d(0)(1)*d(0)(1)+d(0)(0)*d(0)(0)+d(0)(2)*d(0)(2));
			visc(1,1)(0,1) = -cjcbi*(2.*d(0)(1)*d(1)(1)+d(0)(0)*d(1)(0)+d(0)(2)*d(1)(2));
			visc(1,1)(0,2) = -cjcbi*(2.*d(0)(1)*d(2)(1)+d(0)(0)*d(2)(0)+d(0)(2)*d(2)(2));
			visc(1,1)(1,0) = visc(1,1)(0,1);
			visc(1,1)(1,1) = -cjcbi*(2.*d(1)(1)*d(1)(1)+d(1)(0)*d(1)(0)+d(1)(2)*d(1)(2));
			visc(1,1)(1,2) = -cjcbi*(2.*d(1)(1)*d(2)(1)+d(1)(0)*d(2)(0)+d(1)(2)*d(2)(2));
			visc(1,1)(2,0) = visc(1,1)(0,2);
			visc(1,1)(2,1) = visc(1,1)(1,2);
			visc(1,1)(2,2) = -cjcbi*(2.*d(2)(1)*d(2)(1)+d(2)(0)*d(2)(0)+d(2)(2)*d(2)(2));
			
			visc(1,2)(0,0) = -cjcbi*d(0)(2)*d(0)(1);
			visc(1,2)(0,1) = -cjcbi*d(0)(2)*d(1)(1);
			visc(1,2)(0,2) = -cjcbi*d(0)(2)*d(2)(1);
			visc(1,2)(1,0) = -cjcbi*d(1)(2)*d(0)(1);
			visc(1,2)(1,1) = -cjcbi*d(1)(2)*d(1)(1);
			visc(1,2)(1,2) = -cjcbi*d(1)(2)*d(2)(1);
			visc(1,2)(2,0) = -cjcbi*d(2)(2)*d(0)(1);
			visc(1,2)(2,1) = -cjcbi*d(2)(2)*d(1)(1);
			visc(1,2)(2,2) = -cjcbi*d(2)(2)*d(2)(1);
			
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
			
			visc(2,2)(0,0) = -cjcbi*(2.*d(0)(2)*d(0)(2)+d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1));
			visc(2,2)(0,1) = -cjcbi*(2.*d(0)(2)*d(1)(2)+d(0)(0)*d(1)(0)+d(0)(1)*d(1)(1));
			visc(2,2)(0,2) = -cjcbi*(2.*d(0)(2)*d(2)(2)+d(0)(0)*d(2)(0)+d(0)(1)*d(2)(1));
			visc(2,2)(1,0) = visc(2,2)(0,1);
			visc(2,2)(1,1) = -cjcbi*(2.*d(1)(2)*d(1)(2)+d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1));
			visc(2,2)(1,2) = -cjcbi*(2.*d(1)(2)*d(2)(2)+d(1)(0)*d(2)(0)+d(1)(1)*d(2)(1));
			visc(2,2)(2,0) = visc(2,2)(0,2);
			visc(2,2)(2,1) = visc(2,2)(1,2);
			visc(2,2)(2,2) = -cjcbi*(2.*d(2)(2)*d(2)(2)+d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1));
			
			
			/* TIME DERIVATIVE TERMS */ 
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						
						/* UNSTEADY TERMS */
						for(n=0;n<NV-1;++n)
							res(n)(i)(j)(k) = rhorbd0*u(n)(i)(j)(k) +dugdt(log2p,tind,n)(i)(j)(k);
						res(NV-1)(i)(j)(k) = rhorbd0 +dugdt(log2p,tind,NV-1)(i)(j)(k);
						
#ifdef INERTIALESS
						for(n=0;n<NV-1;++n)
							res(n)(i)(j)(k) = 0.0;	
#endif
						
#ifdef BODYFORCE
						for(n=0;n<NV-1;++n)
							res(n)(i)(j)(k) -= gbl->rho*cjcb*gbl->body(n);			
						
#endif                        
						
						
						df(0,0)(i)(j)(k) = +visc(0,0)(0,0)*du(0,0)(i)(j)(k)+visc(0,1)(0,0)*du(1,0)(i)(j)(k)+visc(0,2)(0,0)*du(2,0)(i)(j)(k)
						+visc(0,0)(0,1)*du(0,1)(i)(j)(k)+visc(0,1)(0,1)*du(1,1)(i)(j)(k)+visc(0,2)(0,1)*du(2,1)(i)(j)(k)
						+visc(0,0)(0,2)*du(0,2)(i)(j)(k)+visc(0,1)(0,2)*du(1,2)(i)(j)(k)+visc(0,2)(0,2)*du(2,2)(i)(j)(k);
						
						df(0,1)(i)(j)(k) = +visc(0,0)(1,0)*du(0,0)(i)(j)(k)+visc(0,1)(1,0)*du(1,0)(i)(j)(k)+visc(0,2)(1,0)*du(2,0)(i)(j)(k)
						+visc(0,0)(1,1)*du(0,1)(i)(j)(k)+visc(0,1)(1,1)*du(1,1)(i)(j)(k)+visc(0,2)(1,1)*du(2,1)(i)(j)(k)
						+visc(0,0)(1,2)*du(0,2)(i)(j)(k)+visc(0,1)(1,2)*du(1,2)(i)(j)(k)+visc(0,2)(1,2)*du(2,2)(i)(j)(k);
						
						df(0,2)(i)(j)(k) = +visc(0,0)(2,0)*du(0,0)(i)(j)(k)+visc(0,1)(2,0)*du(1,0)(i)(j)(k)+visc(0,2)(2,0)*du(2,0)(i)(j)(k)
						+visc(0,0)(2,1)*du(0,1)(i)(j)(k)+visc(0,1)(2,1)*du(1,1)(i)(j)(k)+visc(0,2)(2,1)*du(2,1)(i)(j)(k)
						+visc(0,0)(2,2)*du(0,2)(i)(j)(k)+visc(0,1)(2,2)*du(1,2)(i)(j)(k)+visc(0,2)(2,2)*du(2,2)(i)(j)(k);
						
						df(1,0)(i)(j)(k) = +visc(1,0)(0,0)*du(0,0)(i)(j)(k)+visc(1,1)(0,0)*du(1,0)(i)(j)(k)+visc(1,2)(0,0)*du(2,0)(i)(j)(k)
						+visc(1,0)(0,1)*du(0,1)(i)(j)(k)+visc(1,1)(0,1)*du(1,1)(i)(j)(k)+visc(1,2)(0,1)*du(2,1)(i)(j)(k)
						+visc(1,0)(0,2)*du(0,2)(i)(j)(k)+visc(1,1)(0,2)*du(1,2)(i)(j)(k)+visc(1,2)(0,2)*du(2,2)(i)(j)(k);
						
						df(1,1)(i)(j)(k) = +visc(1,0)(1,0)*du(0,0)(i)(j)(k)+visc(1,1)(1,0)*du(1,0)(i)(j)(k)+visc(1,2)(1,0)*du(2,0)(i)(j)(k)
						+visc(1,0)(1,1)*du(0,1)(i)(j)(k)+visc(1,1)(1,1)*du(1,1)(i)(j)(k)+visc(1,2)(1,1)*du(2,1)(i)(j)(k)
						+visc(1,0)(1,2)*du(0,2)(i)(j)(k)+visc(1,1)(1,2)*du(1,2)(i)(j)(k)+visc(1,2)(1,2)*du(2,2)(i)(j)(k);
						
						df(1,2)(i)(j)(k) = +visc(1,0)(2,0)*du(0,0)(i)(j)(k)+visc(1,1)(2,0)*du(1,0)(i)(j)(k)+visc(1,2)(2,0)*du(2,0)(i)(j)(k)
						+visc(1,0)(2,1)*du(0,1)(i)(j)(k)+visc(1,1)(2,1)*du(1,1)(i)(j)(k)+visc(1,2)(2,1)*du(2,1)(i)(j)(k)
						+visc(1,0)(2,2)*du(0,2)(i)(j)(k)+visc(1,1)(2,2)*du(1,2)(i)(j)(k)+visc(1,2)(2,2)*du(2,2)(i)(j)(k);
						
						df(2,0)(i)(j)(k) = +visc(2,0)(0,0)*du(0,0)(i)(j)(k)+visc(2,1)(0,0)*du(1,0)(i)(j)(k)+visc(2,2)(0,0)*du(2,0)(i)(j)(k)
						+visc(2,0)(0,1)*du(0,1)(i)(j)(k)+visc(2,1)(0,1)*du(1,1)(i)(j)(k)+visc(2,2)(0,1)*du(2,1)(i)(j)(k)
						+visc(2,0)(0,2)*du(0,2)(i)(j)(k)+visc(2,1)(0,2)*du(1,2)(i)(j)(k)+visc(2,2)(0,2)*du(2,2)(i)(j)(k);
						
						df(2,1)(i)(j)(k) = +visc(2,0)(1,0)*du(0,0)(i)(j)(k)+visc(2,1)(1,0)*du(1,0)(i)(j)(k)+visc(2,2)(1,0)*du(2,0)(i)(j)(k)
						+visc(2,0)(1,1)*du(0,1)(i)(j)(k)+visc(2,1)(1,1)*du(1,1)(i)(j)(k)+visc(2,2)(1,1)*du(2,1)(i)(j)(k)
						+visc(2,0)(1,2)*du(0,2)(i)(j)(k)+visc(2,1)(1,2)*du(1,2)(i)(j)(k)+visc(2,2)(1,2)*du(2,2)(i)(j)(k);
						
						df(2,2)(i)(j)(k) = +visc(2,0)(2,0)*du(0,0)(i)(j)(k)+visc(2,1)(2,0)*du(1,0)(i)(j)(k)+visc(2,2)(2,0)*du(2,0)(i)(j)(k)
						+visc(2,0)(2,1)*du(0,1)(i)(j)(k)+visc(2,1)(2,1)*du(1,1)(i)(j)(k)+visc(2,2)(2,1)*du(2,1)(i)(j)(k)
						+visc(2,0)(2,2)*du(0,2)(i)(j)(k)+visc(2,1)(2,2)*du(1,2)(i)(j)(k)+visc(2,2)(2,2)*du(2,2)(i)(j)(k);
						
						for(n=0;n<NV-1;++n) {
							cv(n,0)(i)(j)(k) += df(n,0)(i)(j)(k);
							cv(n,1)(i)(j)(k) += df(n,1)(i)(j)(k);
							cv(n,2)(i)(j)(k) += df(n,2)(i)(j)(k);
						}
					}
				}
			}
			for(n=0;n<NV;++n)
				basis::tet(log2p).intgrt(&lf(n)(0),&res(n)(0)(0)(0),stridex,stridey);
			
			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */				
			for(n=0;n<NV-1;++n) {
				basis::tet(log2p).derivr(&cv(n,0)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
				basis::tet(log2p).derivs(&cv(n,1)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
				basis::tet(log2p).derivt(&cv(n,2)(0)(0)(0),&res(n)(0)(0)(0),stridex,stridey);
			}
			basis::tet(log2p).derivr(&du(NV-1,0)(0)(0)(0),&res(NV-1)(0)(0)(0),stridex,stridey);
			basis::tet(log2p).derivs(&du(NV-1,1)(0)(0)(0),&res(NV-1)(0)(0)(0),stridex,stridey);
			basis::tet(log2p).derivt(&du(NV-1,2)(0)(0)(0),&res(NV-1)(0)(0)(0),stridex,stridey);
			
			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {							
						tres(0) = gbl->tau(tind,0)*res(0)(i)(j)(k);    
						tres(1) = gbl->tau(tind,0)*res(1)(i)(j)(k);
						tres(2) = gbl->tau(tind,0)*res(2)(i)(j)(k);
						tres(3) = gbl->tau(tind,NV-1)*res(NV-1)(i)(j)(k);
						
#ifndef INERTIALESS
						
						/* 1st index is equation, 2nd equation is curvilinear coordinate derivative of test function */
						
						df(0,0)(i)(j)(k) -= (d(0)(0)*(2.*u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
											 +d(0)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
											 +d(0)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(0)
						+d(0)(1)*u(0)(i)(j)(k)*tres(1)
						+d(0)(2)*u(0)(i)(j)(k)*tres(2)
						+d(0)(0)*tres(NV-1);
						df(0,1)(i)(j)(k) -= (d(1)(0)*(2.*u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
											 +d(1)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
											 +d(1)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(0)
						+d(1)(1)*u(0)(i)(j)(k)*tres(1)
						+d(1)(2)*u(0)(i)(j)(k)*tres(2)
						+d(1)(0)*tres(NV-1);
						df(0,2)(i)(j)(k) -= (d(2)(0)*(2.*u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
											 +d(2)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
											 +d(2)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(0)
						+d(2)(1)*u(0)(i)(j)(k)*tres(1)
						+d(2)(2)*u(0)(i)(j)(k)*tres(2)
						+d(2)(0)*tres(NV-1);
						df(1,0)(i)(j)(k) -= d(0)(0)*u(1)(i)(j)(k)*tres(0)
						+(d(0)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(0)(1)*(2.*u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(0)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(1)
						+d(0)(2)*u(1)(i)(j)(k)*tres(2)				
						+d(0)(1)*tres(NV-1);		
						df(1,1)(i)(j)(k) -= d(1)(0)*u(1)(i)(j)(k)*tres(0)
						+(d(1)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(1)(1)*(2.*u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(1)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(1)
						+d(1)(2)*u(1)(i)(j)(k)*tres(2)				
						+d(1)(1)*tres(NV-1);		
						df(1,2)(i)(j)(k) -= d(2)(0)*u(1)(i)(j)(k)*tres(0)
						+(d(2)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(2)(1)*(2.*u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(2)(2)*(u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(1)
						+d(2)(2)*u(1)(i)(j)(k)*tres(2)				
						+d(2)(1)*tres(NV-1);
						df(2,0)(i)(j)(k) -= d(0)(0)*u(2)(i)(j)(k)*tres(0)
						+d(0)(1)*u(2)(i)(j)(k)*tres(1)
						+(d(0)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(0)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(0)(2)*(2.*u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(2)
						+d(0)(2)*tres(NV-1);	
						df(2,1)(i)(j)(k) -= d(1)(0)*u(2)(i)(j)(k)*tres(0)
						+d(1)(1)*u(2)(i)(j)(k)*tres(1)
						+(d(1)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(1)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(1)(2)*(2.*u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(2)
						+d(1)(2)*tres(NV-1);	
						df(2,2)(i)(j)(k) -= d(2)(0)*u(2)(i)(j)(k)*tres(0)
						+d(2)(1)*u(2)(i)(j)(k)*tres(1)
						+(d(2)(0)*(u(0)(i)(j)(k)-mvel(0)(i)(j)(k))
						  +d(2)(1)*(u(1)(i)(j)(k)-mvel(1)(i)(j)(k))
						  +d(2)(2)*(2.*u(2)(i)(j)(k)-mvel(2)(i)(j)(k)))*tres(2)
						+d(2)(2)*tres(NV-1);
#endif
						
						
						du(NV-1,0)(i)(j)(k) = -d(0)(0)*tres(0)-d(0)(1)*tres(1)-d(0)(2)*tres(2);
						du(NV-1,1)(i)(j)(k) = -d(1)(0)*tres(0)-d(1)(1)*tres(1)-d(1)(2)*tres(2);					
						du(NV-1,2)(i)(j)(k) = -d(2)(0)*tres(0)-d(2)(1)*tres(1)-d(2)(2)*tres(2);
					}
				}
			}
			for(n=0;n<NV-1;++n)
				basis::tet(log2p).intgrtrst(&lf(n)(0),&df(n,0)(0)(0)(0),&df(n,1)(0)(0)(0),&df(n,2)(0)(0)(0),stridex,stridey);
			basis::tet(log2p).intgrtrst(&lf(NV-1)(0),&du(NV-1,0)(0)(0)(0),&du(NV-1,1)(0)(0)(0),&du(NV-1,2)(0)(0)(0),stridex,stridey);
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
					lf(n)(i) *= gbl->beta(stage);
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
					dw(n,i)=lf(n)(i);
		}
	}
	
	
	return;
}
