/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cns.h"
#include "../hp_boundary.h"

//#define BODYFORCE

void tri_hp_cns::element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im){

	FLT fluxx,fluxy;
	const int NV = 4;
	TinyVector<int,3> v;
	TinyVector<FLT,ND> pt;
	TinyMatrix<FLT,ND,ND> ldcrd;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
	TinyMatrix<FLT,NV,NV> A,B;
	TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-1,NV-1> visc;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,NV> cv, df;
	TinyMatrix<FLT,ND,ND> kcond;
	TinyVector<FLT,NV> tres;
    TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> mvel; // for local mesh velocity info
    const int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
	const FLT gam = gbl->gamma;
	const FLT gm1 = gam-1.0;
	const FLT ogm1 = 1.0/gm1;
	const FLT gogm1 = gam*ogm1;
#ifndef SUTHERLAND
    const FLT lmu = gbl->mu;
    const FLT lkcond = gbl->kcond;
#endif

	/* LOAD INDICES OF VERTEX POINTS */
	v = tri(tind).pnt;

	/* IF TINFO > -1 IT IS CURVED ELEMENT */
	if (tri(tind).info > -1) {
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		crdtocht(tind);

		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		/* dcrd(i,j) is derivative of physical coordinate i with respect to curvilinear coordinate j */
		/* dxi/dx = dy/deta, dxi/dy = -dx/deta, deta/dx = -dy/dxi, deta/dy = dx/dxi (divided by jacobian) */
		for(int n = 0; n < ND; ++n)
			basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
	}
	else {
		/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(int n = 0; n < ND; ++n)
			basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);

		/* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
		for(int n = 0; n < ND; ++n) {
			ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
			ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
		}
	}

	/* CALCULATE MESH VELOCITY */
	for(int i = 0; i < lgpx; ++i) {
		for(int j = 0; j < lgpn; ++j) {
			mvel(0)(i,j) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
			mvel(1)(i,j) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j));
#ifdef MESH_REF_VEL
			mvel(0)(i,j) += mesh_ref_vel(0);
			mvel(1)(i,j) += mesh_ref_vel(1);
#endif                        
		}
	}
	
	/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
	/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
	if (gbl->beta(stage) > 0.0) {
		basis::tri(log2p)->proj(&uht(0)(0),&u(0)(0,0),MXGP);
		for(int n = 1; n < NV; ++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
	}
	else {
		for(int n = 0; n < NV; ++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);
	}
	
	/* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
	for(int n = 0; n < NV; ++n) {
		for(int m = 0; m < basis::tri(log2p)->tm(); ++m) {
			lf_re(n)(m) = 0.0;
			lf_im(n)(m) = 0.0;	
		}
	}

	if (tri(tind).info > -1) {
		/* CURVED ELEMENT */
		/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
		for(int i = 0; i < lgpx; ++i) {
			for(int j = 0; j < lgpn; ++j) {
				const FLT rho = u(0)(i,j)/u(NV-1)(i,j);
				fluxx = rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(0)(i,j));
				fluxy = rho*RAD(crd(0)(i,j))*(u(2)(i,j) -mvel(1)(i,j));

				/* CONTINUITY EQUATION FLUXES */
				cv(0,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
				cv(0,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;

#ifndef INERTIALESS
				/* MOMENTUM & ENERGY FLUXES */
				for(int n = 1; n < NV-1; ++n) {
					cv(n,0)(i,j) = u(n)(i,j)*cv(0,0)(i,j);
					cv(n,1)(i,j) = u(n)(i,j)*cv(0,1)(i,j);
				}
#else
				for(n=1;n<NV-1;++n) {
					cv(n,0)(i,j) = 0.0;
					cv(n,1)(i,j) = 0.0;
				}
#endif

				/* PRESSURE TERMS */
				/* U-MOMENTUM */
				cv(1,0)(i,j) += dcrd(1,1)(i,j)*RAD(crd(0)(i,j))*u(0)(i,j);
				cv(1,1)(i,j) -= dcrd(1,0)(i,j)*RAD(crd(0)(i,j))*u(0)(i,j);
				/* V-MOMENTUM */
				cv(2,0)(i,j) -=  dcrd(0,1)(i,j)*RAD(crd(0)(i,j))*u(0)(i,j);
				cv(2,1)(i,j) +=  dcrd(0,0)(i,j)*RAD(crd(0)(i,j))*u(0)(i,j);

				/* ENERGY FLUX */
                const FLT E = u(3)(i,j)/gm1+0.5*(u(1)(i,j)*u(1)(i,j)+u(2)(i,j)*u(2)(i,j));
                const FLT new_fluxx = E*fluxx+u(0)(i,j)*RAD(crd(0)(i,j))*u(1)(i,j);
                const FLT new_fluxy = E*fluxy+u(0)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
                
                cv(3,0)(i,j) = +dcrd(1,1)(i,j)*new_fluxx -dcrd(0,1)(i,j)*new_fluxy;
                cv(3,1)(i,j) = -dcrd(1,0)(i,j)*new_fluxx +dcrd(0,0)(i,j)*new_fluxy;
			}
		}
		for(int n = 0; n < NV; ++n)
			basis::tri(log2p)->intgrtrs(&lf_im(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);

		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			/* TIME DERIVATIVE TERMS */ 
			for(int i = 0; i < lgpx; ++i) {
				for(int j = 0; j < lgpn; ++j) {
                    
#ifdef SUTHERLAND
                    Sutherland(u(NV-1)(i,j));
                    const FLT lmu = gbl->mu;
                    const FLT lkcond = gbl->kcond;
#endif
                    
					const FLT rho = u(0)(i,j)/u(NV-1)(i,j);
					const FLT cjcb = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
					const FLT rhorbd0 = rho*gbl->bd(0)*RAD(crd(0)(i,j))*cjcb;
                    const FLT mujcbi = lmu*RAD(crd(0)(i,j))/cjcb;
					const FLT kcjcbi = lkcond*RAD(crd(0)(i,j))/cjcb/gbl->R;

					/* UNSTEADY TERMS */
					res(0)(i,j) = rhorbd0 +dugdt(log2p)(tind,0,i,j);
					for(int n = 1; n < NV-1; ++n)
						res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p)(tind,n,i,j);
					const FLT e = ogm1*u(NV-1)(i,j) +0.5*(u(1)(i,j)*u(1)(i,j) +u(2)(i,j)*u(2)(i,j));
					res(NV-1)(i,j) = rhorbd0*e +dugdt(log2p)(tind,NV-1,i,j);
					
#ifdef BODYFORCE
					res(1)(i,j) -= rho*RAD(crd(0)(i,j))*cjcb*gbl->body(0);
					res(2)(i,j) -= rho*RAD(crd(0)(i,j))*cjcb*gbl->body(1);
					res(3)(i,j) -= rho*cjcb*RAD(crd(0)(i,j))*(u(1)(i,j)*gbl->body(0)+u(2)(i,j)*gbl->body(1));
#endif     
					
#ifdef MMS
					/* source terms for MMS */
					pt(0) = crd(0)(i,j);
					pt(1) = crd(1)(i,j);
					for(int n = 0; n < NV; ++n)
						res(n)(i,j) -= cjcb*gbl->src->f(n,pt,gbl->time);
#endif

					/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
					/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
					visc(0,0)(0,0) = -mujcbi*(4./3.*dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
					visc(0,0)(1,1) = -mujcbi*(4./3.*dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
					visc(0,0)(0,1) =  mujcbi*(4./3.*dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI0II0II1II0I visc(0,0)(0,1)
					
					visc(1,1)(0,0) = -mujcbi*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +4./3.*dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
					visc(1,1)(1,1) = -mujcbi*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +4./3.*dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
					visc(1,1)(0,1) =  mujcbi*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +4./3.*dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI1II1II1II0I visc(1,1)(0,1)
					
					visc(0,1)(0,0) =  mujcbi*1./3.*dcrd(0,1)(i,j)*dcrd(1,1)(i,j);
					visc(0,1)(1,1) =  mujcbi*1./3.*dcrd(0,0)(i,j)*dcrd(1,0)(i,j);
					visc(0,1)(0,1) = -mujcbi*(dcrd(0,1)(i,j)*dcrd(1,0)(i,j)-2./3.*dcrd(0,0)(i,j)*dcrd(1,1)(i,j));
					visc(0,1)(1,0) = -mujcbi*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j)-2./3.*dcrd(0,1)(i,j)*dcrd(1,0)(i,j));
					
					/* OTHER SYMMETRIES     */                
#define				viscI1II0II0II0I visc(0,1)(0,0)
#define				viscI1II0II1II1I visc(0,1)(1,1)
#define				viscI1II0II0II1I visc(0,1)(1,0)
#define				viscI1II0II1II0I visc(0,1)(0,1)                      

					/* HEAT DIFFUSION TENSOR */ 
					kcond(0,0) = -kcjcbi*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
					kcond(1,1) = -kcjcbi*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
					kcond(0,1) =  kcjcbi*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define				kcondI1II0I kcond(0,1)

					/* Momentum equations */
					df(1,0)(i,j) = +visc(0,0)(0,0)*du(1,0)(i,j) +visc(0,1)(0,0)*du(2,0)(i,j)
									+visc(0,0)(0,1)*du(1,1)(i,j) +visc(0,1)(0,1)*du(2,1)(i,j);

					df(1,1)(i,j) = +viscI0II0II1II0I*du(1,0)(i,j) +visc(0,1)(1,0)*du(2,0)(i,j)
									+visc(0,0)(1,1)*du(1,1)(i,j) +visc(0,1)(1,1)*du(2,1)(i,j);

					df(2,0)(i,j) = +viscI1II0II0II0I*du(1,0)(i,j) +visc(1,1)(0,0)*du(2,0)(i,j)
									+viscI1II0II0II1I*du(1,1)(i,j) +visc(1,1)(0,1)*du(2,1)(i,j);

					df(2,1)(i,j) = +viscI1II0II1II0I*du(1,0)(i,j) +viscI1II1II1II0I*du(2,0)(i,j)
									+viscI1II0II1II1I*du(1,1)(i,j) +visc(1,1)(1,1)*du(2,1)(i,j);
				
					/* Energy equation */
					df(3,0)(i,j) = +kcond(0,0)*du(NV-1,0)(i,j) +kcond(0,1)*du(NV-1,1)(i,j);
					df(3,1)(i,j) = +kcondI1II0I*du(NV-1,0)(i,j) +kcond(1,1)*du(NV-1,1)(i,j);

					/* VISCOUS DISSIPATION */
					df(3,0)(i,j) += df(1,0)(i,j)*u(1)(i,j) +df(2,0)(i,j)*u(2)(i,j);
					df(3,1)(i,j) += df(1,1)(i,j)*u(1)(i,j) +df(2,1)(i,j)*u(2)(i,j);

					for(int n = 1; n < NV; ++n) {
						cv(n,0)(i,j) += df(n,0)(i,j);
						cv(n,1)(i,j) += df(n,1)(i,j);
					}
				}
			}
	
			for(int n = 0; n < NV; ++n)
				basis::tri(log2p)->intgrt(&lf_re(n)(0),&res(n)(0,0),MXGP);

			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(int n = 0; n < NV; ++n) {
				basis::tri(log2p)->derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
				basis::tri(log2p)->derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
			}

			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(int i=0;i<lgpx;++i) {
				for(int j=0;j<lgpn;++j) {
					
					df(0,0)(i,j) = 0.0;
					df(0,1)(i,j) = 0.0;
					
					tres = 0.0;
					for(int m = 0; m < NV; ++m)
						for(int n = 0; n < NV; ++n)							
							tres(m) += gbl->tau(tind,m,n)*res(n)(i,j);

					FLT pr = u(0)(i,j);
					FLT uv = u(1)(i,j);
					FLT vv = u(2)(i,j);
					FLT rt = u(3)(i,j);	
		
					FLT rho = pr/rt;
					FLT ke = 0.5*(uv*uv+vv*vv);
					
					/* df/dw derivative of fluxes wrt primitive variables*/
					A = uv/rt,               rho,                         0.0,       -rho*uv/rt,
						uv*uv/rt+1.0,        2.0*rho*uv,                  0.0,       -rho*uv*uv/rt,
						uv*vv/rt,            rho*vv,                      rho*uv,    -rho*uv*vv/rt,
						uv*(gogm1*rt+ke)/rt, rho*(gogm1*rt+ke)+rho*uv*uv, rho*uv*vv, -rho*uv*(gogm1*rt+ke)/rt+rho*uv*gogm1;
										
					/* dg/dw */
					B = vv/rt,               0.0,       rho,                         -rho*vv/rt,
					    uv*vv/rt,            rho*vv,    rho*uv,                      -rho*uv*vv/rt,
						vv*vv/rt+1.0,        0.0,       2.0*rho*vv,                  -rho*vv*vv/rt,
						vv*(gogm1*rt+ke)/rt, rho*uv*vv, rho*(gogm1*rt+ke)+rho*vv*vv, -rho*vv*(gogm1*rt+ke)/rt+rho*vv*gogm1;
					
					for(int m = 0; m < NV; ++m) {
						for(int n = 0; n < NV; ++n) {
							df(m,0)(i,j) -= (dcrd(1,1)(i,j)*A(m,n)-dcrd(0,1)(i,j)*B(m,n))*tres(n);
							df(m,1)(i,j) -= (-dcrd(1,0)(i,j)*A(m,n)+dcrd(0,0)(i,j)*B(m,n))*tres(n);
						}
					}					
				}
			}			
			
			for(int n = 0; n < NV; ++n)
				basis::tri(log2p)->intgrtrs(&lf_re(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
			
			for(int n = 0; n < NV; ++n)
				for(int m = 0; m < basis::tri(log2p)->tm(); ++m)
					lf_re(n)(m) *= gbl->beta(stage);

		}
	}
	else {
		/* STRAIGHT EDGED ELEMENT */
		/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
		for(int i = 0; i < lgpx; ++i) {
			for(int j = 0; j < lgpn; ++j) {
				const FLT rho = u(0)(i,j)/u(NV-1)(i,j);
				fluxx = rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(0)(i,j));
				fluxy = rho*RAD(crd(0)(i,j))*(u(2)(i,j) -mvel(1)(i,j));
				
				/* CONTINUITY EQUATION FLUXES */
				cv(0,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
				cv(0,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;
				
#ifndef INERTIALESS
				/* MOMENTUM & ENERGY FLUXES */
				for(int n = 1; n < NV-1; ++n) {
					cv(n,0)(i,j) = u(n)(i,j)*cv(0,0)(i,j);
					cv(n,1)(i,j) = u(n)(i,j)*cv(0,1)(i,j);
				}
#else
				for(n=1;n<NV-1;++n) {
					cv(n,0)(i,j) = 0.0;
					cv(n,1)(i,j) = 0.0;
				}
#endif
				
				/* PRESSURE TERMS */
				/* U-MOMENTUM */
				cv(1,0)(i,j) += ldcrd(1,1)*RAD(crd(0)(i,j))*u(0)(i,j);
				cv(1,1)(i,j) -= ldcrd(1,0)*RAD(crd(0)(i,j))*u(0)(i,j);
				/* V-MOMENTUM */
				cv(2,0)(i,j) -=  ldcrd(0,1)*RAD(crd(0)(i,j))*u(0)(i,j);
				cv(2,1)(i,j) +=  ldcrd(0,0)*RAD(crd(0)(i,j))*u(0)(i,j);
				
				/* ENERGY FLUX */
                const FLT E = u(3)(i,j)/gm1+0.5*(u(1)(i,j)*u(1)(i,j)+u(2)(i,j)*u(2)(i,j));
                const FLT new_fluxx = E*fluxx+u(0)(i,j)*RAD(crd(0)(i,j))*u(1)(i,j);
                const FLT new_fluxy = E*fluxy+u(0)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
                
                cv(3,0)(i,j) = +ldcrd(1,1)*new_fluxx -ldcrd(0,1)*new_fluxy;
                cv(3,1)(i,j) = -ldcrd(1,0)*new_fluxx +ldcrd(0,0)*new_fluxy;
			}
		}
		for(int n = 0; n < NV; ++n)
			basis::tri(log2p)->intgrtrs(&lf_im(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
		
		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			/* TIME DERIVATIVE TERMS */ 
			for(int i = 0; i < lgpx; ++i) {
				for(int j = 0; j < lgpn; ++j) {
                    
#ifdef SUTHERLAND
                    Sutherland(u(NV-1)(i,j));
                    const FLT lmu = gbl->mu;
                    const FLT lkcond = gbl->kcond;
#endif
                    
					const FLT rho = u(0)(i,j)/u(NV-1)(i,j);
					const FLT cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
					const FLT rhorbd0 = rho*gbl->bd(0)*RAD(crd(0)(i,j))*cjcb;
					const FLT mujcbi = lmu*RAD(crd(0)(i,j))/cjcb;
					const FLT kcjcbi = lkcond*RAD(crd(0)(i,j))/cjcb/gbl->R;
					
					/* UNSTEADY TERMS */
					res(0)(i,j) = rhorbd0 +dugdt(log2p)(tind,0,i,j);
					for(int n = 1; n < NV-1; ++n)
						res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p)(tind,n,i,j);
					const FLT e = ogm1*u(NV-1)(i,j) +0.5*(u(1)(i,j)*u(1)(i,j) +u(2)(i,j)*u(2)(i,j));
					res(NV-1)(i,j) = rhorbd0*e +dugdt(log2p)(tind,NV-1,i,j);
					
#ifdef BODYFORCE
					res(1)(i,j) -= rho*RAD(crd(0)(i,j))*cjcb*gbl->body(0);
					res(2)(i,j) -= rho*RAD(crd(0)(i,j))*cjcb*gbl->body(1);
					res(3)(i,j) -= rho*cjcb*RAD(crd(0)(i,j))*(u(1)(i,j)*gbl->body(0)+u(2)(i,j)*gbl->body(1));

#endif     
					
#ifdef MMS
					/* source terms for MMS */
					pt(0) = crd(0)(i,j);
					pt(1) = crd(1)(i,j);
					for(int n = 0; n < NV; ++n)
						res(n)(i,j) -= cjcb*gbl->src->f(n,pt,gbl->time);
#endif
					
					/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
					/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
					visc(0,0)(0,0) = -mujcbi*(4./3.*ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
					visc(0,0)(1,1) = -mujcbi*(4./3.*ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
					visc(0,0)(0,1) =  mujcbi*(4./3.*ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define                 viscI0II0II1II0I visc(0,0)(0,1)
					
					visc(1,1)(0,0) = -mujcbi*(ldcrd(1,1)*ldcrd(1,1) +4./3.*ldcrd(0,1)*ldcrd(0,1));
					visc(1,1)(1,1) = -mujcbi*(ldcrd(1,0)*ldcrd(1,0) +4./3.*ldcrd(0,0)*ldcrd(0,0));
					visc(1,1)(0,1) =  mujcbi*(ldcrd(1,1)*ldcrd(1,0) +4./3.*ldcrd(0,1)*ldcrd(0,0));
#define                 viscI1II1II1II0I visc(1,1)(0,1)
					
					visc(0,1)(0,0) =  mujcbi*1./3.*ldcrd(0,1)*ldcrd(1,1);
					visc(0,1)(1,1) =  mujcbi*1./3.*ldcrd(0,0)*ldcrd(1,0);
					visc(0,1)(0,1) = -mujcbi*(ldcrd(0,1)*ldcrd(1,0)-2./3.*ldcrd(0,0)*ldcrd(1,1));
					visc(0,1)(1,0) = -mujcbi*(ldcrd(0,0)*ldcrd(1,1)-2./3.*ldcrd(0,1)*ldcrd(1,0));
					
					/* OTHER SYMMETRIES     */                
#define				viscI1II0II0II0I visc(0,1)(0,0)
#define				viscI1II0II1II1I visc(0,1)(1,1)
#define				viscI1II0II0II1I visc(0,1)(1,0)
#define				viscI1II0II1II0I visc(0,1)(0,1)                      
					
					/* HEAT DIFFUSION TENSOR */ 
					kcond(0,0) = -kcjcbi*(ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
					kcond(1,1) = -kcjcbi*(ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
					kcond(0,1) =  kcjcbi*(ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define				kcondI1II0I kcond(0,1)
					
					/* Momentum equations */
					df(1,0)(i,j) = +visc(0,0)(0,0)*du(1,0)(i,j) +visc(0,1)(0,0)*du(2,0)(i,j)
					+visc(0,0)(0,1)*du(1,1)(i,j) +visc(0,1)(0,1)*du(2,1)(i,j);
					
					df(1,1)(i,j) = +viscI0II0II1II0I*du(1,0)(i,j) +visc(0,1)(1,0)*du(2,0)(i,j)
					+visc(0,0)(1,1)*du(1,1)(i,j) +visc(0,1)(1,1)*du(2,1)(i,j);
					
					df(2,0)(i,j) = +viscI1II0II0II0I*du(1,0)(i,j) +visc(1,1)(0,0)*du(2,0)(i,j)
					+viscI1II0II0II1I*du(1,1)(i,j) +visc(1,1)(0,1)*du(2,1)(i,j);
					
					df(2,1)(i,j) = +viscI1II0II1II0I*du(1,0)(i,j) +viscI1II1II1II0I*du(2,0)(i,j)
					+viscI1II0II1II1I*du(1,1)(i,j) +visc(1,1)(1,1)*du(2,1)(i,j);
					
					/* Energy equation */
					df(3,0)(i,j) = +kcond(0,0)*du(NV-1,0)(i,j) +kcond(0,1)*du(NV-1,1)(i,j);
					df(3,1)(i,j) = +kcondI1II0I*du(NV-1,0)(i,j) +kcond(1,1)*du(NV-1,1)(i,j);
					
					/* VISCOUS DISSIPATION */
					df(3,0)(i,j) += df(1,0)(i,j)*u(1)(i,j) +df(2,0)(i,j)*u(2)(i,j);
					df(3,1)(i,j) += df(1,1)(i,j)*u(1)(i,j) +df(2,1)(i,j)*u(2)(i,j);
					
					for(int n = 1; n < NV; ++n) {
						cv(n,0)(i,j) += df(n,0)(i,j);
						cv(n,1)(i,j) += df(n,1)(i,j);
					}
				}
			}
			
			for(int n = 0; n < NV; ++n)
				basis::tri(log2p)->intgrt(&lf_re(n)(0),&res(n)(0,0),MXGP);
			
			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(int n = 0; n < NV; ++n) {
				basis::tri(log2p)->derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
				basis::tri(log2p)->derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
			}
			
			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(int i=0;i<lgpx;++i) {
				for(int j=0;j<lgpn;++j) {
					
					df(0,0)(i,j) = 0.0;
					df(0,1)(i,j) = 0.0;
					
					tres = 0.0;
					for(int m = 0; m < NV; ++m)
						for(int n = 0; n < NV; ++n)							
							tres(m) += gbl->tau(tind,m,n)*res(n)(i,j);
					
					FLT pr = u(0)(i,j);
					FLT uv = u(1)(i,j);
					FLT vv = u(2)(i,j);
					FLT rt = u(3)(i,j);	

					FLT rho = pr/rt;
					FLT ke = 0.5*(uv*uv+vv*vv);
					
					/* df/dw derivative of fluxes wrt primitive variables */
					A = uv/rt,               rho,                         0.0,       -rho*uv/rt,
						uv*uv/rt+1.0,        2.0*rho*uv,                  0.0,       -rho*uv*uv/rt,
						uv*vv/rt,            rho*vv,                      rho*uv,    -rho*uv*vv/rt,
						uv*(gogm1*rt+ke)/rt, rho*(gogm1*rt+ke)+rho*uv*uv, rho*uv*vv, -rho*uv*(gogm1*rt+ke)/rt+rho*uv*gogm1;
					
					/* dg/dw */
					B = vv/rt,               0.0,       rho,                         -rho*vv/rt,
						uv*vv/rt,            rho*vv,    rho*uv,                      -rho*uv*vv/rt,
						vv*vv/rt+1.0,        0.0,       2.0*rho*vv,                  -rho*vv*vv/rt,
						vv*(gogm1*rt+ke)/rt, rho*uv*vv, rho*(gogm1*rt+ke)+rho*vv*vv, -rho*vv*(gogm1*rt+ke)/rt+rho*vv*gogm1;
						
					for(int m = 0; m < NV; ++m) {
						for(int n = 0; n < NV; ++n) {
							df(m,0)(i,j) -= (+ldcrd(1,1)*A(m,n)-ldcrd(0,1)*B(m,n))*tres(n);
							df(m,1)(i,j) -= (-ldcrd(1,0)*A(m,n)+ldcrd(0,0)*B(m,n))*tres(n);
						}
					}					
				}
			}
			
			
			for(int n = 0; n < NV; ++n)
				basis::tri(log2p)->intgrtrs(&lf_re(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
			
			for(int n = 0; n < NV; ++n)
				for(int m = 0; m < basis::tri(log2p)->tm(); ++m)
					lf_re(n)(m) *= gbl->beta(stage);
			
		}
	}



	return;
}


#ifdef SUTHERLAND
void tri_hp_cns::Sutherland(FLT RT) {
    gbl->mu = gbl->s1*pow(RT,1.5)/(RT+gbl->s2);
    gbl->kcond = gbl->R*gbl->mu/gbl->prandtl*gbl->gamma/(gbl->gamma-1.0);
    return;
}
#endif
