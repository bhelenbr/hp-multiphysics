/*
 *  tri_hp_nonnewtonian.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/29/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */

#include "tri_hp_nonnewtonian.h"

void tri_hp_nonnewtonian::init(input_map& inmap, void *gin) {
	gbl = static_cast<global *>(gin);
	tri_hp_ins::init(inmap,gin);
	
	Array<string,1> names(4);
	Array<int,1> dims(4);
	dims = ND;
	names(0) = "u";
	dims(0) = NV;
	names(1) = "x";
	names(2) = "xt";
	dims(3) = 3;
	names(3) = "s";

	gbl->mu_of_strain.set_arguments(4,dims,names);

	if (inmap.find(gbl->idprefix +"_mu_function") != inmap.end()) {
		gbl->mu_of_strain.init(inmap,gbl->idprefix +"_mu_function");
	}
	else if (inmap.find("mu_function") != inmap.end()){
		gbl->mu_of_strain.init(inmap,"mu_function");
	}
	else {
		*gbl->log << "couldn't find nonnewtonian viscosity function " << gbl->idprefix +"_mu_function" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	return;
}

void tri_hp_nonnewtonian::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_hp_nonnewtonian& inmesh = dynamic_cast<const tri_hp_nonnewtonian &>(in);
	gbl = inmesh.gbl;
	tri_hp_ins::init(in,why,sizereduce1d);
	return;
}


void tri_hp_nonnewtonian::element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im) {
	int i,j,n;
	FLT fluxx,fluxy;
	const int NV = 3;
	TinyVector<int,3> v;
	TinyMatrix<FLT,ND,ND> ldcrd;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
	int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
	FLT rhobd0 = gbl->rho*gbl->bd(0), rhorbd0, lrhorbd0, cjcb, cjcbi;
	TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-1,NV-1> visc;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV-1,NV-1> cv, df;
	TinyVector<FLT,NV> tres;
	TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> mvel; // for local mesh velocity info

	
	/* LOAD INDICES OF VERTEX POINTS */
	v = tri(tind).pnt;
	
	/* IF TINFO > -1 IT IS CURVED ELEMENT */
	if (tri(tind).info > -1) {
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		crdtocht(tind);
		
		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		/* dcrd(i,j) is derivative of physical coordinate i with respect to curvilinear coordinate j */
		/* dxi/dx = dy/deta, dxi/dy = -dx/deta, deta/dx = -dy/dxi, deta/dy = dx/dxi (divided by jacobian) */
		for(n=0;n<ND;++n)
			basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
	}
	else {
		/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(n=0;n<ND;++n)
			basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);
		
		/* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
		for(n=0;n<ND;++n) {
			ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
			ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
		}
	}
	
	/* CALCULATE MESH VELOCITY */
	for(i=0;i<lgpx;++i) {
		for(j=0;j<lgpn;++j) {
			mvel(0)(i,j) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
			mvel(1)(i,j) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j));
#ifdef MESH_REF_VEL
			mvel(0)(i,j) += gbl->mesh_ref_vel(0);
			mvel(1)(i,j) += gbl->mesh_ref_vel(1);
#endif                        
		}
	}
	
	/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
	/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
	//ugtouht(tind);
	if (gbl->beta(stage) > 0.0) {
		for(n=0;n<NV-1;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
		basis::tri(log2p)->proj(&uht(NV-1)(0),&u(NV-1)(0,0),MXGP);
	}
	else {
		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);
	}
	
	/* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
	for(n=0;n<NV;++n){
		for(i=0;i<basis::tri(log2p)->tm();++i){
			lf_re(n)(i) = 0.0;
			lf_im(n)(i) = 0.0;	
		}
	}
	
	if (tri(tind).info > -1) {
		/* CURVED ELEMENT */
		/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpn;++j) {
				
				fluxx = gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
				fluxy = gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
				
				/* CONTINUITY EQUATION FLUXES */
				du(NV-1,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
				du(NV-1,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
				
#ifndef INERTIALESS
				/* CONVECTIVE FLUXES */
				for(n=0;n<NV-1;++n) {
					cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
					cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
				}
#else
				for(n=0;n<NV-1;++n) {
					cv(n,0)(i,j) = 0.0;
					cv(n,1)(i,j) = 0.0;
				}
#endif
				/* PRESSURE TERMS */
				/* U-MOMENTUM */
				cv(0,0)(i,j) += dcrd(1,1)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				cv(0,1)(i,j) -= dcrd(1,0)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				/* V-MOMENTUM */
				cv(1,0)(i,j) -=  dcrd(0,1)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				cv(1,1)(i,j) +=  dcrd(0,0)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
			}
		}
		for(n=0;n<NV-1;++n)
			basis::tri(log2p)->intgrtrs(&lf_im(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
		basis::tri(log2p)->intgrtrs(&lf_im(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
		
		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			/* TIME DERIVATIVE TERMS */ 
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					cjcb = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
					
					/* Calculate Non-Newtonian Viscosity */
					Array<FLT,1> axpt(ND), amv(ND), as(3), au(NV);
					for (n=0;n<NV-1;++n)
						au(n) = u(n)(i,j);
					axpt(0) = crd(0)(i,j); axpt(1) = crd(1)(i,j);
					amv(0) = mvel(0)(i,j); amv(1) = mvel(1)(i,j);
					as(0) = dcrd(1,1)(i,j)*du(0,0)(i,j) -dcrd(1,0)(i,j)*du(0,1)(i,j);
					as(1) = dcrd(0,0)(i,j)*du(1,1)(i,j) -dcrd(0,1)(i,j)*du(1,0)(i,j);
					as(2) = dcrd(1,1)(i,j)*du(1,0)(i,j) -dcrd(1,0)(i,j)*du(1,1)(i,j)+ dcrd(0,0)(i,j)*du(0,1)(i,j) -dcrd(0,1)(i,j)*du(0,0)(i,j); 
					as /= cjcb;
					FLT lmu = gbl->mu_of_strain.Eval(au,axpt,amv,as,gbl->time);	
					
					
					rhorbd0 = rhobd0*RAD(crd(0)(i,j))*cjcb;
					cjcbi = lmu*RAD(crd(0)(i,j))/cjcb;
					
					/* UNSTEADY TERMS */
					for(n=0;n<NV-1;++n)
						res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p)(tind,n,i,j);
					res(NV-1)(i,j) = rhorbd0 +dugdt(log2p)(tind,NV-1,i,j);
					
#ifdef INERTIALESS
					res(0)(i,j) = 0.0;
					res(1)(i,j) = 0.0;
#endif
#ifdef AXISYMMETRIC
					res(0)(i,j) -= cjcb*(u(NV-1)(i,j) -2.*lmu*u(0)(i,j)/crd(0)(i,j));
#endif
#ifdef BODYFORCE
					res(0)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(0);
					res(1)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(1);
#endif                  
					
					/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
					/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
					visc(0,0)(0,0) = -cjcbi*(2.*dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
					visc(0,0)(1,1) = -cjcbi*(2.*dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
					visc(0,0)(0,1) =  cjcbi*(2.*dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI0II0II1II0I visc(0,0)(0,1)
					
					visc(1,1)(0,0) = -cjcbi*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +2.*dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
					visc(1,1)(1,1) = -cjcbi*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +2.*dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
					visc(1,1)(0,1) =  cjcbi*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +2.*dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI1II1II1II0I visc(1,1)(0,1)
					
					visc(0,1)(0,0) =  cjcbi*dcrd(0,1)(i,j)*dcrd(1,1)(i,j);
					visc(0,1)(1,1) =  cjcbi*dcrd(0,0)(i,j)*dcrd(1,0)(i,j);
					visc(0,1)(0,1) = -cjcbi*dcrd(0,1)(i,j)*dcrd(1,0)(i,j);
					visc(0,1)(1,0) = -cjcbi*dcrd(0,0)(i,j)*dcrd(1,1)(i,j);
					
					/* OTHER SYMMETRIES     */                
#define                 viscI1II0II0II0I visc(0,1)(0,0)
#define                 viscI1II0II1II1I visc(0,1)(1,1)
#define                 viscI1II0II0II1I visc(0,1)(1,0)
#define                 viscI1II0II1II0I visc(0,1)(0,1) 			
					
					df(0,0)(i,j) = +visc(0,0)(0,0)*du(0,0)(i,j) +visc(0,1)(0,0)*du(1,0)(i,j)
					+visc(0,0)(0,1)*du(0,1)(i,j) +visc(0,1)(0,1)*du(1,1)(i,j);
					
					df(0,1)(i,j) = +viscI0II0II1II0I*du(0,0)(i,j) +visc(0,1)(1,0)*du(1,0)(i,j)
					+visc(0,0)(1,1)*du(0,1)(i,j) +visc(0,1)(1,1)*du(1,1)(i,j);
					
					df(1,0)(i,j) = +viscI1II0II0II0I*du(0,0)(i,j) +visc(1,1)(0,0)*du(1,0)(i,j)
					+viscI1II0II0II1I*du(0,1)(i,j) +visc(1,1)(0,1)*du(1,1)(i,j);
					
					df(1,1)(i,j) = +viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
					+viscI1II0II1II1I*du(0,1)(i,j) +visc(1,1)(1,1)*du(1,1)(i,j);
					
					for(n=0;n<NV-1;++n) {
						cv(n,0)(i,j) += df(n,0)(i,j);
						cv(n,1)(i,j) += df(n,1)(i,j);
					}
				}
			}
			for(n=0;n<NV;++n)
				basis::tri(log2p)->intgrt(&lf_re(n)(0),&res(n)(0,0),MXGP);
			
			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(n=0;n<NV-1;++n) {
				basis::tri(log2p)->derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
				basis::tri(log2p)->derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
			}
			basis::tri(log2p)->derivr(&du(NV-1,0)(0,0),&res(NV-1)(0,0),MXGP);
			basis::tri(log2p)->derivs(&du(NV-1,1)(0,0),&res(NV-1)(0,0),MXGP);
			
			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					
					tres(0) = gbl->tau(tind,0)*res(0)(i,j);
					tres(1) = gbl->tau(tind,0)*res(1)(i,j);
					tres(NV-1) = gbl->tau(tind,NV-1)*res(NV-1)(i,j);
					
#ifndef INERTIALESS
					df(0,0)(i,j) -= (dcrd(1,1)(i,j)*(2*u(0)(i,j)-mvel(0)(i,j))
													 -dcrd(0,1)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(0)
					-dcrd(0,1)(i,j)*u(0)(i,j)*tres(1)
					+dcrd(1,1)(i,j)*tres(NV-1);
					df(0,1)(i,j) -= (-dcrd(1,0)(i,j)*(2*u(0)(i,j)-mvel(0)(i,j))
													 +dcrd(0,0)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(0)
					+dcrd(0,0)(i,j)*u(0)(i,j)*tres(1)
					-dcrd(1,0)(i,j)*tres(NV-1);
					df(1,0)(i,j) -= +dcrd(1,1)(i,j)*u(1)(i,j)*tres(0)
					+(dcrd(1,1)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
						-dcrd(0,1)(i,j)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres(1)
					-dcrd(0,1)(i,j)*tres(NV-1);
					df(1,1)(i,j) -= -dcrd(1,0)(i,j)*u(1)(i,j)*tres(0)
					+(-dcrd(1,0)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
						+dcrd(0,0)(i,j)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres(1)
					+dcrd(0,0)(i,j)*tres(NV-1);
#endif
					
					du(NV-1,0)(i,j) = -(dcrd(1,1)(i,j)*tres(0) -dcrd(0,1)(i,j)*tres(1));
					du(NV-1,1)(i,j) = -(-dcrd(1,0)(i,j)*tres(0) +dcrd(0,0)(i,j)*tres(1));
				}
			}
			for(n=0;n<NV-1;++n)
				basis::tri(log2p)->intgrtrs(&lf_re(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
			basis::tri(log2p)->intgrtrs(&lf_re(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tri(log2p)->tm();++i)
					lf_re(n)(i) *= gbl->beta(stage);
			
		}
	}
	else {
		/* LINEAR ELEMENT */
		/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpn;++j) {
				
				fluxx = gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
				fluxy = gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
				
				/* CONTINUITY EQUATION FLUXES */
				du(NV-1,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
				du(NV-1,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;
				
#ifndef INERTIALESS
				/* CONVECTIVE FLUXES */
				for(n=0;n<NV-1;++n) {
					cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
					cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
				}
#else
				for(n=0;n<NV-1;++n) {
					cv(n,0)(i,j) = 0.0;
					cv(n,1)(i,j) = 0.0;
				}
#endif
				/* PRESSURE TERMS */
				/* U-MOMENTUM */
				cv(0,0)(i,j) += ldcrd(1,1)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				cv(0,1)(i,j) -= ldcrd(1,0)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				/* V-MOMENTUM */
				cv(1,0)(i,j) -=  ldcrd(0,1)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				cv(1,1)(i,j) +=  ldcrd(0,0)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
			}
		}
		for(n=0;n<NV-1;++n)
			basis::tri(log2p)->intgrtrs(&lf_im(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
		basis::tri(log2p)->intgrtrs(&lf_im(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
		
		
		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
			cjcbi = 1./cjcb;
			lrhorbd0 = rhobd0*cjcb;
			
			/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
			/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
			visc(0,0)(0,0) = -cjcbi*(2.*ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
			visc(0,0)(1,1) = -cjcbi*(2.*ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
			visc(0,0)(0,1) =  cjcbi*(2.*ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define         viscI0II0II1II0I visc(0,0)(0,1)
			
			visc(1,1)(0,0) = -cjcbi*(ldcrd(1,1)*ldcrd(1,1) +2.*ldcrd(0,1)*ldcrd(0,1));
			visc(1,1)(1,1) = -cjcbi*(ldcrd(1,0)*ldcrd(1,0) +2.*ldcrd(0,0)*ldcrd(0,0));
			visc(1,1)(0,1) =  cjcbi*(ldcrd(1,1)*ldcrd(1,0) +2.*ldcrd(0,1)*ldcrd(0,0));
#define         viscI1II1II1II0I visc(1,1)(0,1)
			
			visc(0,1)(0,0) =  cjcbi*ldcrd(0,1)*ldcrd(1,1);
			visc(0,1)(1,1) =  cjcbi*ldcrd(0,0)*ldcrd(1,0);
			visc(0,1)(0,1) = -cjcbi*ldcrd(0,1)*ldcrd(1,0);
			visc(0,1)(1,0) = -cjcbi*ldcrd(0,0)*ldcrd(1,1);
			
			/* OTHER SYMMETRIES     */                
#define         viscI1II0II0II0I visc(0,1)(0,0)
#define         viscI1II0II1II1I visc(0,1)(1,1)
#define         viscI1II0II0II1I visc(0,1)(1,0)
#define         viscI1II0II1II0I visc(0,1)(0,1)
			
			/* TIME DERIVATIVE TERMS */ 
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					rhorbd0 = RAD(crd(0)(i,j))*lrhorbd0;
					
					/* Calculate Non-Newtonian Viscosity */
					Array<FLT,1> axpt(ND), amv(ND), as(3), au(NV);
					for (n=0;n<NV-1;++n)
						au(n) = u(n)(i,j);
					axpt(0) = crd(0)(i,j); axpt(1) = crd(1)(i,j);
					amv(0) = mvel(0)(i,j); amv(1) = mvel(1)(i,j);
					as(0) = ldcrd(1,1)*du(0,0)(i,j) -ldcrd(1,0)*du(0,1)(i,j);
					as(1) = ldcrd(0,0)*du(1,1)(i,j) -ldcrd(0,1)*du(1,0)(i,j);
					as(2) = ldcrd(1,1)*du(1,0)(i,j) -ldcrd(1,0)*du(1,1)(i,j)+ ldcrd(0,0)*du(0,1)(i,j) -ldcrd(0,1)*du(0,0)(i,j); 
					as /= cjcb;
					FLT lmu = gbl->mu_of_strain.Eval(au,axpt,amv,as,gbl->time);	
					
					/* UNSTEADY TERMS */
					for(n=0;n<NV-1;++n)
						res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p)(tind,n,i,j);
					res(NV-1)(i,j) = rhorbd0 +dugdt(log2p)(tind,NV-1,i,j);
					
#ifdef INERTIALESS
					res(0)(i,j) = 0.0;
					res(1)(i,j) = 0.0;
#endif
#ifdef AXISYMMETRIC
					res(0)(i,j) -= cjcb*(u(NV-1)(i,j) -2.*lmu*u(0)(i,j)/crd(0)(i,j));
#endif
#ifdef BODYFORCE
					res(0)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(0);
					res(1)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(1);
#endif        
					df(0,0)(i,j) = lmu*RAD(crd(0)(i,j))*(+visc(0,0)(0,0)*du(0,0)(i,j) +visc(0,1)(0,0)*du(1,0)(i,j)
																					 +visc(0,0)(0,1)*du(0,1)(i,j) +visc(0,1)(0,1)*du(1,1)(i,j));
					
					df(0,1)(i,j) = lmu*RAD(crd(0)(i,j))*(+viscI0II0II1II0I*du(0,0)(i,j) +visc(0,1)(1,0)*du(1,0)(i,j)
																					 +visc(0,0)(1,1)*du(0,1)(i,j) +visc(0,1)(1,1)*du(1,1)(i,j));
					
					df(1,0)(i,j) = lmu*RAD(crd(0)(i,j))*(+viscI1II0II0II0I*du(0,0)(i,j) +visc(1,1)(0,0)*du(1,0)(i,j)
																					 +viscI1II0II0II1I*du(0,1)(i,j) +visc(1,1)(0,1)*du(1,1)(i,j));
					
					df(1,1)(i,j) = lmu*RAD(crd(0)(i,j))*(+viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
																					 +viscI1II0II1II1I*du(0,1)(i,j) +visc(1,1)(1,1)*du(1,1)(i,j));
					
					for(n=0;n<NV-1;++n) {
						cv(n,0)(i,j) += df(n,0)(i,j);
						cv(n,1)(i,j) += df(n,1)(i,j);
					} 
				}
			}
			for(n=0;n<NV;++n)
				basis::tri(log2p)->intgrt(&lf_re(n)(0),&res(n)(0,0),MXGP);
			
			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(n=0;n<NV-1;++n) {
				basis::tri(log2p)->derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
				basis::tri(log2p)->derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
			}
			basis::tri(log2p)->derivr(&du(NV-1,0)(0,0),&res(NV-1)(0,0),MXGP);
			basis::tri(log2p)->derivs(&du(NV-1,1)(0,0),&res(NV-1)(0,0),MXGP);
			
			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					tres(0) = gbl->tau(tind,0)*res(0)(i,j);
					tres(1) = gbl->tau(tind,0)*res(1)(i,j);
					tres(NV-1) = gbl->tau(tind,NV-1)*res(NV-1)(i,j);
					
#ifndef INERTIALESS
					df(0,0)(i,j) -= (ldcrd(1,1)*(2*u(0)(i,j)-mvel(0)(i,j))
													 -ldcrd(0,1)*(u(1)(i,j)-mvel(1)(i,j)))*tres(0)
					-ldcrd(0,1)*u(0)(i,j)*tres(1)
					+ldcrd(1,1)*tres(NV-1);
					df(0,1)(i,j) -= (-ldcrd(1,0)*(2*u(0)(i,j)-mvel(0)(i,j))
													 +ldcrd(0,0)*(u(1)(i,j)-mvel(1)(i,j)))*tres(0)
					+ldcrd(0,0)*u(0)(i,j)*tres(1)
					-ldcrd(1,0)*tres(NV-1);
					df(1,0)(i,j) -= +ldcrd(1,1)*u(1)(i,j)*tres(0)
					+(ldcrd(1,1)*(u(0)(i,j)-mvel(0)(i,j))
						-ldcrd(0,1)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres(1)
					-ldcrd(0,1)*tres(NV-1);
					df(1,1)(i,j) -= -ldcrd(1,0)*u(1)(i,j)*tres(0)
					+(-ldcrd(1,0)*(u(0)(i,j)-mvel(0)(i,j))
						+ldcrd(0,0)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres(1)
					+ldcrd(0,0)*tres(NV-1);
#endif
					
					du(NV-1,0)(i,j) = -(ldcrd(1,1)*tres(0) -ldcrd(0,1)*tres(1));
					du(NV-1,1)(i,j) = -(-ldcrd(1,0)*tres(0) +ldcrd(0,0)*tres(1));
				}
			}
			for(n=0;n<NV-1;++n)
				basis::tri(log2p)->intgrtrs(&lf_re(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
			basis::tri(log2p)->intgrtrs(&lf_re(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tri(log2p)->tm();++i)
					lf_re(n)(i) *= gbl->beta(stage);
			
		}
	}
	
	return;
}

void tri_hp_nonnewtonian::setup_preconditioner() {
	/* SET-UP DIAGONAL PRECONDITIONER */
	int tind,i,j,side;
	FLT jcb,h,hmax,q,qmax,lam1,gam;
	TinyVector<int,3> v;
		
	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
	if (basis::tri(log2p)->sm() > 0) {
		gbl->sprcn(Range(0,nseg-1),Range::all()) = 0.0;
	}
	
	for(tind = 0; tind < ntri; ++tind) {
		jcb = 0.25*area(tind);  // area is 2 x triangle area
		v = tri(tind).pnt;
		
		/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
		/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
		ugtouht(tind);
		for(int n=0;n<ND;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
		basis::tri(log2p)->proj(&uht(NV-1)(0),&u(NV-1)(0,0),MXGP);
		
		
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		crdtocht(tind);
		
		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(int n=0;n<ND;++n)
			basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
		
		TinyVector<FLT,ND> mvel;
		qmax = 0.0;
		hmax = 0.0;
		FLT jcbmin = jcb;
		FLT nu = 0.0;
		int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpn;++j) {
				
				mvel(0) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
				mvel(1) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j));
#ifdef MESH_REF_VEL
				mvel(0) += gbl->mesh_ref_vel(0);
				mvel(1) += gbl->mesh_ref_vel(1);
#endif
				FLT cjcb = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);

				jcbmin = MIN(jcbmin,cjcb);
				
				/* Calculate Non-Newtonian Viscosity */
				Array<FLT,1> axpt(ND), amv(ND), as(3), au(NV);
				for (int n=0;n<NV-1;++n)
					au(n) = u(n)(i,j);
				axpt(0) = crd(0)(i,j); axpt(1) = crd(1)(i,j);
				amv(0) = mvel(0); amv(1) = mvel(1);
				as(0) = dcrd(1,1)(i,j)*du(0,0)(i,j) -dcrd(1,0)(i,j)*du(0,1)(i,j);
				as(1) = dcrd(0,0)(i,j)*du(1,1)(i,j) -dcrd(0,1)(i,j)*du(1,0)(i,j);
				as(2) = dcrd(1,1)(i,j)*du(1,0)(i,j) -dcrd(1,0)(i,j)*du(1,1)(i,j)+ dcrd(0,0)(i,j)*du(0,1)(i,j) -dcrd(0,1)(i,j)*du(0,0)(i,j); 
				as /= cjcb;
				FLT lmu = gbl->mu_of_strain.Eval(au,axpt,amv,as,gbl->time);	
				nu = MAX(nu,lmu/gbl->rho);
				
				/* CALCULATE CURVED SIDE LENGTHS */
				h = 0.0;
				for (int n=0;n<ND;++n)
					h += dcrd(n,0)(i,j)*dcrd(n,0)(i,j);
				hmax = MAX(h,hmax);
				
				h = 0.0;
				for (int n=0;n<ND;++n)
					h += dcrd(n,1)(i,j)*dcrd(n,1)(i,j);
				hmax = MAX(h,hmax);
				
				h = 0.0;
				for (int n=0;n<ND;++n)
					h += (dcrd(n,1)(i,j) -dcrd(n,0)(i,j))*(dcrd(n,1)(i,j) -dcrd(n,0)(i,j));
				hmax = MAX(h,hmax);
				
				q = pow(u(0)(i,j)-0.5*mvel(0),2.0)  +pow(u(1)(i,j)-0.5*mvel(1),2.0);
				qmax = MAX(qmax,q);
			}
		}	
		hmax = 2.*sqrt(hmax);
		h = 4.*jcbmin/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
		hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
		
		if (!(h > 0.0)) { 
			*gbl->log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
			*gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1) << std::endl;
			tri_mesh::output("negative",grid);
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		if  (std::isnan(qmax)) { 
			*gbl->log << gbl->idprefix << ' ' << tind << std::endl;
			*gbl->log << "flow solution has nan's " << qmax << std::endl;
			output("nan",tecplot);
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
		if (gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
		q = sqrt(qmax);
		lam1 = q + sqrt(qmax +gam);
		
		/* SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
		gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);
		
		/* SET UP DIAGONAL PRECONDITIONER */
		jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
		jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);
		
		gbl->tprcn(tind,0) = gbl->rho*jcb;   
		gbl->tprcn(tind,1) = gbl->rho*jcb;   
		gbl->tprcn(tind,2) = jcb/gam; 			
		
		for(i=0;i<3;++i) {
			gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
			if (basis::tri(log2p)->sm() > 0) {
				side = tri(tind).seg(i);
				gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
			}
		}
	}
	
	tri_hp::setup_preconditioner();
}
