/*
 *  bdry.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 12/13/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */

#include "bdry_buoyancy.h"
#include <tri_boundary.h>

using namespace bdry_buoyancy;

void surface9::init(input_map& inmap,void* gbl_in) {
	bdry_ins::surface::init(inmap,gbl_in);
	
	if (inmap.find(base.idprefix+"_sigma_vs_T") != inmap.end()) {
		sigma_vs_T.init(inmap,base.idprefix+"_sigma_vs_T");
	} 
	else if (inmap.find("sigma_vs_T") != inmap.end()){
		sigma_vs_T.init(inmap,"sigma_vs_T");
	}
	else {
		*x.gbl->log << "couldn't find sigma_vs_T equation for surface tension" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}

	return;
}



/* Free-Surface that Allows Radiative Heat Flux & Marangoni Effects */
void surface9::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
	int i,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV),flx(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,9,MXGP> res;
	Array<FLT,1> axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND),au(x.NV);

	sind = base.seg(indx);
	v0 = x.seg(sind).pnt(0);
	v1 = x.seg(sind).pnt(1);
	
	x.crdtocht1d(sind);
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
	
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));    
	
	for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
		norm(0) =  dcrd(1,i);
		norm(1) = -dcrd(0,i);
		jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
		
		/* RELATIVE VELOCITY STORED IN MVEL(N)*/
		for(n=0;n<tri_mesh::ND;++n) {
			mvel(n,i) = u(n)(i) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef MESH_REF_VEL
			mvel(n,i) -= x.gbl->mesh_ref_vel(n);
#endif
		}
	
		/* Evaluate Fluxes */
		axpt(0) = crd(0,i); axpt(1) = crd(1,i);
		amv(0) = -(mvel(0,i)-u(0)(i)); amv(1) = -(mvel(1,i)-u(1)(i));
		anorm(0)= norm(0); anorm(1) = norm(1);
		for(int n=0;n<x.NV;++n)
			au(n) = u(n)(i);
		for(int n=0;n<x.NV;++n) {
			flx(n) = fluxes[n].Eval(au,axpt,amv,anorm,x.gbl->time);
		}
		
		/* Evaluate Surface Tension */
		FLT sigma = sigma_vs_T.Eval(au(2));
		/* TANGENTIAL SPACING */                
		res(0,i) = -ksprg(indx)*jcb;
		/* NORMAL FLUX */
		res(1,i) = -RAD(crd(0,i))*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1)) +flx(x.NV-1);     
		/* UPWINDING BASED ON TANGENTIAL VELOCITY */
		res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*gbl->meshc(indx);
		
		/* SURFACE TENSION SOURCE TERM X-DIRECTION */ 
		res(4,i) = +RAD(crd(0,i))*((x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i) +gbl->p_ext)*norm(0) +flx(0)*jcb;
#ifdef AXISYMMETRIC
		res(4,i) += sigma*jcb;
#endif
		/* AND INTEGRATION BY PARTS TERM */
		res(5,i) = +RAD(crd(0,i))*sigma*norm(1)/jcb;
		
		
		/* SURFACE TENSION SOURCE TERM Y-DIRECTION */
		res(6,i) = +RAD(crd(0,i))*((x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i) +gbl->p_ext)*norm(1) +flx(1)*jcb;                
		/* AND INTEGRATION BY PARTS TERM */
		res(7,i) = -RAD(crd(0,i))*sigma*norm(0)/jcb;
		
		/* Auxiliary Fluxes Here */
		res(8,i) = flx(2)*jcb;
	}
	
	lf = 0.0;
	/* INTEGRATE & STORE SURFACE TENSION SOURCE TERM */
	basis::tri(x.log2p)->intgrt1d(&lf(0)(0),&res(4,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(0)(0),&res(5,0));
	basis::tri(x.log2p)->intgrt1d(&lf(1)(0),&res(6,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(1)(0),&res(7,0));
	
	/* Heat Flux Source Term */
	basis::tri(x.log2p)->intgrt1d(&lf(2)(0),&res(8,0));

	
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */                    
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV)(0),&res(0,0));
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+1)(0),&res(1,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV+1)(0),&res(2,0));
	
#ifndef petsc
	/* mass flux preconditioning */
	for(int m=0;m<basis::tri(x.log2p)->sm()+2;++m)
		lf(x.NV-1)(m) = -x.gbl->rho*lf(x.NV+1)(m); 
#ifndef INERTIALESS
	for (n=0;n<x.NV-1;++n) 
		ubar(n) = 0.5*(x.uht(n)(0) +x.uht(n)(1));
	
	for (n=0;n<x.NV-1;++n) {
		lf(n)(0) -= x.uht(n)(0)*(x.gbl->rho -gbl->rho2)*lf(x.NV+1)(0);
		lf(n)(1) -= x.uht(n)(1)*(x.gbl->rho -gbl->rho2)*lf(x.NV+1)(1);
		for(int m=0;m<basis::tri(x.log2p)->sm();++m)
			lf(n)(m+2) -= ubar(n)*(x.gbl->rho -gbl->rho2)*lf(x.NV+1)(m+2);
	}
#endif
#endif
	
	return;
}

