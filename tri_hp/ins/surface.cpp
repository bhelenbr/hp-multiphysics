//
//  surface.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 12/30/14.
//
//

#include "surface.h"

#include "bdry_ins.h"
#include <myblas.h>

//#define MPDEBUG

//#define DEBUG

//#define BODYFORCE

using namespace bdry_ins;

// extern FLT body[ND];



void surface::init(input_map& inmap,void* gin) {
	std::string keyword,matching_block,side_id,master_block,master_id;
	std::istringstream data;
	std::string filename;
	
	find_matching_boundary_name(inmap, matching_block, side_id);
	
	gbl = static_cast<global *>(gin);
	gbl->mu2 = 0.0;
	gbl->rho2 = 0.0;
	master_block = x.gbl->idprefix;
	is_master = true;
	if (base.is_comm()) {
		keyword = matching_block +"_mu";
		if (!inmap.get(keyword,gbl->mu2)) {
			*x.gbl->log << "couldn't find matching blocks viscosity " << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		
		keyword = matching_block +"_rho";
		if (!inmap.get(keyword,gbl->rho2)) {
			*x.gbl->log << "couldn't find matching blocks density" << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		
		// Decide who is the master
		if (fabs(x.gbl->rho - gbl->rho2) < 100.*FLT_EPSILON) {
			if (x.gbl->idprefix < matching_block) {
				is_master = true;
				master_block = x.gbl->idprefix;
			}
			else {
				is_master = false;
				master_block = matching_block;
			}
		}
		else if (x.gbl->rho > gbl->rho2) {
			is_master = true;
			master_block = x.gbl->idprefix;
		}
		else {
			is_master = false;
			master_block = matching_block;
		}
	}
	master_id = master_block +"_" +side_id;
	
	
	/* Let pressure be discontinuous */
	ostringstream vars;
	for(int n=0;n<x.NV-1;++n) {
		vars << n << ' ';
	}
	inmap[base.idprefix +"_c0_indices"] = vars.str();
	hp_coupled_bdry::init(inmap,gin);
	
	inmap.getwdefault(master_id +"_sigma",gbl->sigma,0.0);
	inmap.getwdefault(master_id +"_p_ext",gbl->p_ext,0.0);
	
	return;
}

void surface::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
	
	if (!is_master) return;
	
	int i,n,sind;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,8,MXGP> res;
	
	sind = base.seg(indx);	
	x.crdtocht1d(sind);
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
	
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));
	
//	/* Calculate stabilization constant based on analysis of linear elements and constant tau */
//	int v0 = x.seg(sind).pnt(0);
//	int v1 = x.seg(sind).pnt(1);
//	norm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
//	norm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
//	FLT h = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
//	FLT hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
//
//	mvel(0) = x.uht(0)(0) -(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
//	mvel(1) = x.uht(1)(0)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
//#ifdef MESH_REF_VEL
//	mvel(0) -= x.gbl->mesh_ref_vel(0);
//	mvel(1) -= x.gbl->mesh_ref_vel(1);
//#endif
//	FLT vslp0 = (-mvel(0)*norm(1) +mvel(1)*norm(0))/h;
//
//	mvel(0) = x.uht(0)(1)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
//	mvel(1) = x.uht(1)(1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
//#ifdef MESH_REF_VEL
//	mvel(0) -= x.gbl->mesh_ref_vel(0);
//	mvel(1) -= x.gbl->mesh_ref_vel(1);
//#endif
//	FLT vslp1 = (-mvel(0)*norm(1) +mvel(1)*norm(0))/h;
//	gbl->meshc(indx) = gbl->adis*hsm*(3*(abs(vslp0)+abs(vslp1)) +vslp0-vslp1)/(4*(vslp0*vslp0+vslp0*vslp1+vslp1*vslp1)+10*FLT_EPSILON)*2/h;
	
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
		
		/* TANGENTIAL SPACING */
		res(0,i) = -ksprg(indx)*jcb;
		/* NORMAL FLUX */
        res(1,i) = -RAD(crd(0,i))*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1));        
		/* UPWINDING BASED ON TANGENTIAL VELOCITY */
		res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*gbl->meshc(indx);
		
		/* surface TENSION SOURCE TERM X-DIRECTION */
		res(4,i) = +RAD(crd(0,i))*((x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i) +gbl->p_ext)*norm(0);
#ifdef AXISYMMETRIC
		res(4,i) += gbl->sigma*jcb;
#endif
		/* AND INTEGRATION BY PARTS TERM */
		res(5,i) = +RAD(crd(0,i))*gbl->sigma*norm(1)/jcb;
		
		
		/* surface TENSION SOURCE TERM Y-DIRECTION */
		res(6,i) = +RAD(crd(0,i))*((x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i) +gbl->p_ext)*norm(1);
		/* AND INTEGRATION BY PARTS TERM */
		res(7,i) = -RAD(crd(0,i))*gbl->sigma*norm(0)/jcb;
	}
	
	lf = 0.0;
	/* INTEGRATE & STORE surface TENSION SOURCE TERM */
	basis::tri(x.log2p)->intgrt1d(&lf(0)(0),&res(4,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(0)(0),&res(5,0));
	basis::tri(x.log2p)->intgrt1d(&lf(1)(0),&res(6,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(1)(0),&res(7,0));
	
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV)(0),&res(0,0));
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+1)(0),&res(1,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV+1)(0),&res(2,0));

	return;
}

void surface::setup_preconditioner() {
	int indx,m,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> nrm;
	FLT h, hsm;
	FLT dttang, dtnorm;
	FLT strss;
	FLT drho;
	FLT nu1, nu2;
	FLT gam1, gam2;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd;
	TinyMatrix<FLT,4,MXGP> res;
	TinyMatrix<FLT,4,MXGP> lf;
	TinyVector<FLT,2> mvel;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	int last_phase, mp_phase;
	
	hp_coupled_bdry::setup_preconditioner();

	if (!gbl->symmetric && !is_master) return;
	
	drho = x.gbl->rho -gbl->rho2;
	nu1 = x.gbl->mu/x.gbl->rho;
	if (gbl->rho2 > 0.0)
		nu2 = gbl->mu2/gbl->rho2;
	else
		nu2 = 0.0;
	
	/**************************************************/
	/* DETERMINE surface MOVEMENT TIME STEP              */
	/**************************************************/
	gbl->vdt(0,Range::all(),Range::all()) = 0.0;
	
	for(indx=0; indx < base.nseg; ++indx) {
		sind = base.seg(indx);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		
#ifdef DETAILED_DT
		x.crdtocht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
		
		x.ugtouht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));
		
		dtnorm = 1.0e99;
		dttang = 1.0e99;
		gbl->meshc(indx) = 1.0e99;
		for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
			nrm(0) =  dcrd(1,i)*2;
			nrm(1) = -dcrd(0,i)*2;
			h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
			
			/* RELATIVE VELOCITY STORED IN MVEL(N)*/
			for(n=0;n<tri_mesh::ND;++n) {
				mvel(n) = u(n)(i) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef MESH_REF_VEL
				mvel(n) -= x.gbl->mesh_ref_vel(n);
#endif
			}
			qmax = u(0)(i)*u(0)(i) +u(1)(i)*u(1)(i);
			vslp = fabs(-u(0)(i)*nrm(1)/h +u(1)(i)*nrm(0)/h);
			hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
			
			dttang = MIN(dttang,2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm);
#ifndef BODYFORCE
			strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*gbl->g*nrm(1)/h);
#else
			strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*(-gbl->body(0)*nrm(0) +(gbl->g-gbl->body(1))*nrm(1))/h);
#endif
			
			gam1 = 3.0*qmax +(0.5*hsm*x.gbl->bd(0) + 2.*nu1/hsm)*(0.5*hsm*x.gbl->bd(0) + 2.*nu1/hsm);
			gam2 = 3.0*qmax +(0.5*hsm*x.gbl->bd(0) + 2.*nu2/hsm)*(0.5*hsm*x.gbl->bd(0) + 2.*nu2/hsm);
			
			if (x.gbl->bd(0) + x.gbl->mu == 0.0) gam1 = MAX(gam1,0.1);
			
#ifdef INERTIALESS
			gam1 = (2.*nu1/hsm)*(2.*nu1/hsm);
			gam2 = (2.*nu2/hsm)*(2.*nu2/hsm);
#endif
			dtnorm = MIN(dtnorm,2.*vslp/hsm +x.gbl->bd(0) +1.*strss/(x.gbl->rho*sqrt(qmax +gam1) +gbl->rho2*sqrt(qmax +gam2)));
			
			/* SET UP DISSIPATIVE COEFFICIENT */
			/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
			/* RESIDUAL HAS DX/2 WEIGHTING */
			/* |a| dx/2 dv/dx  dx/2 dpsi */
			/* |a| dx/2 2/dx dv/dpsi  dpsi */
			/* |a| dv/dpsi  dpsi */
			// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5);/* FAILED IN NATES UPSTREAM surface WAVE CASE */
			// gbl->meshc(indx) = MIN(gbl->meshc(indx),gbl->adis/(h*(vslp/hsm +x.gbl->bd(0)))); /* FAILED IN MOVING UP TESTS */
			gbl->meshc(indx) = MIN(gbl->meshc(indx),gbl->adis/(h*(sqrt(qmax)/hsm +x.gbl->bd(0)))); /* SEEMS THE BEST I'VE GOT */
		}
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
#else
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
		h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
		
		mvel(0) = x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
		mvel(1) = x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
#ifdef MESH_REF_VEL
		mvel(0) -= x.gbl->mesh_ref_vel(0);
		mvel(1) -= x.gbl->mesh_ref_vel(1);
#endif
		
		FLT qmax0 = mvel(0)*mvel(0)+mvel(1)*mvel(1);
		FLT vslp0 = (-mvel(0)*nrm(1) +mvel(1)*nrm(0))/h;
		
		mvel(0) = x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
		mvel(1) = x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
#ifdef MESH_REF_VEL
		mvel(0) -= x.gbl->mesh_ref_vel(0);
		mvel(1) -= x.gbl->mesh_ref_vel(1);
#endif
		FLT qmax1 = mvel(0)*mvel(0)+mvel(1)*mvel(1);
		FLT vslp1 = (-mvel(0)*nrm(1) +mvel(1)*nrm(0))/h;
		
		hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
		
		dttang = 2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm;
#ifndef BODYFORCE
		strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*x.gbl->g*nrm(1)/h);
#else
		strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*(-gbl->body(0)*nrm(0) +(gbl->g-gbl->body(1))*nrm(1))/h);
#endif
        FLT qmax = MAX(qmax0,qmax1);
        FLT vslp = MAX(abs(vslp0),abs(vslp1));
		gam1 = 3.0*qmax +(0.5*hsm*x.gbl->bd(0) + 2.*nu1/hsm)*(0.5*hsm*x.gbl->bd(0) + 2.*nu1/hsm);
		gam2 = 3.0*qmax +(0.5*hsm*x.gbl->bd(0) + 2.*nu2/hsm)*(0.5*hsm*x.gbl->bd(0) + 2.*nu2/hsm);
		
		if (x.gbl->bd(0) + x.gbl->mu == 0.0) gam1 = MAX(gam1,0.1);
		
#ifdef INERTIALESS
		gam1 = (2.*nu1/hsm)*(2.*nu1/hsm);
		gam2 = (2.*nu2/hsm)*(2.*nu2/hsm);
#endif
		dtnorm = 2.*vslp/hsm +x.gbl->bd(0) +1.*strss/(x.gbl->rho*sqrt(qmax +gam1) +gbl->rho2*sqrt(qmax +gam2));
		
		/* SET UP DISSIPATIVE COEFFICIENT */
		/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
		/* RESIDUAL HAS DX/2 WEIGHTING */
		/* |a| dx/2 dv/dx  dx/2 dpsi */
		/* |a| dx/2 2/dx dv/dpsi  dpsi */
		/* |a| dv/dpsi  dpsi */
		// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5); /* FAILED IN NATES UPSTREAM surface WAVE CASE */
		// gbl->meshc(indx) = gbl->adis/(h*(vslp/hsm +x.gbl->bd(0))); /* FAILED IN MOVING UP TESTS */
		 gbl->meshc(indx) = gbl->adis/(h*(sqrt(qmax)/hsm +x.gbl->bd(0))); /* SEEMS THE BEST I'VE GOT */
       //gbl->meshc(indx) = gbl->adis*hsm*(3*(abs(vslp0)+abs(vslp1)) +vslp0-vslp1)/(4*(vslp0*vslp0+vslp0*vslp1+vslp1*vslp1)+10*FLT_EPSILON)*2/h;
#endif
		
		dtnorm *= RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));
		
		nrm *= 0.5;
		
		gbl->vdt(indx,0,0) += -dttang*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx,0,1) +=  dttang*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx,1,0) +=  dtnorm*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx,1,1) +=  dtnorm*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,0,0) = -dttang*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,0,1) =  dttang*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,1,0) =  dtnorm*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,1,1) =  dtnorm*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		
		if (basis::tri(x.log2p)->sm()) {
			gbl->sdt(indx,0,0) = -dttang*nrm(1);
			gbl->sdt(indx,0,1) =  dttang*nrm(0);
			gbl->sdt(indx,1,0) =  dtnorm*nrm(0);
			gbl->sdt(indx,1,1) =  dtnorm*nrm(1);
			
#ifdef DETAILED_MINV
			int lsm = basis::tri(x.log2p)->sm();
			x.crdtocht1d(sind);
			for(n=0;n<tri_mesh::ND;++n)
				basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
			
			for(int m = 0; m<lsm; ++m) {
				for(int i=0;i<basis::tri(x.log2p)->gpx();++i) {
					nrm(0) =  dcrd(1,i);
					nrm(1) = -dcrd(0,i);
					res(0,i) = -dttang*nrm(1)*basis::tri(x.log2p)->gx(i,m+3);
					res(1,i) =  dttang*nrm(0)*basis::tri(x.log2p)->gx(i,m+3);
					res(2,i) =  dtnorm*nrm(0)*basis::tri(x.log2p)->gx(i,m+3);
					res(3,i) =  dtnorm*nrm(1)*basis::tri(x.log2p)->gx(i,m+3);
				}
				lf = 0;
				basis::tri(x.log2p)->intgrt1d(&lf(0,0),&res(0,0));
				basis::tri(x.log2p)->intgrt1d(&lf(1,0),&res(1,0));
				basis::tri(x.log2p)->intgrt1d(&lf(2,0),&res(2,0));
				basis::tri(x.log2p)->intgrt1d(&lf(3,0),&res(3,0));
				
				/* CFL = 0 WON'T WORK THIS WAY */
				lf(0) /= gbl->cfl(x.log2p,0);
				lf(1) /= gbl->cfl(x.log2p,0);
				lf(2) /= gbl->cfl(x.log2p,1);
				lf(3) /= gbl->cfl(x.log2p,1);
				
				for (n=0;n<lsm;++n) {
					gbl->ms(indx,2*m,2*n) = lf(0,n+2);
					gbl->ms(indx,2*m,2*n+1) = lf(1,n+2);
					gbl->ms(indx,2*m+1,2*n) = lf(2,n+2);
					gbl->ms(indx,2*m+1,2*n+1) = lf(3,n+2);
				}
				
				/* tang/norm, x/y,  mode,  vert */
				gbl->vms(indx,0,0,m,0) = lf(0,0);
				gbl->vms(indx,0,1,m,0) = lf(1,0);
				gbl->vms(indx,0,0,m,1) = lf(0,1);
				gbl->vms(indx,0,1,m,1) = lf(1,1);
				gbl->vms(indx,1,0,m,0) = lf(2,0);
				gbl->vms(indx,1,1,m,0) = lf(3,0);
				gbl->vms(indx,1,0,m,1) = lf(2,1);
				gbl->vms(indx,1,1,m,1) = lf(3,1);
			}
			
			int info;
			GETRF(2*lsm,2*lsm,&gbl->ms(indx,0,0),2*MAXP,&gbl->ipiv(indx,0),info);
			if (info != 0) {
				*x.gbl->log << "DGETRF FAILED IN SIDE MODE PRECONDITIONER\n";
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			/*
			 \phi_n dx,dy*t = \phi_n Vt
			 \phi_t dx,dy*n = \phi_t Vn
			 */
#endif
		}
	}
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vdt(0,0,0),0,gbl->vdt.length(secondDim)*gbl->vdt.length(secondDim)-1,0);
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vdt(base.nseg,0,0),0,gbl->vdt.length(secondDim)*gbl->vdt.length(secondDim)-1,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		
		last_phase = true;
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(0,0,0),0,gbl->vdt.length(secondDim)*gbl->vdt.length(secondDim)-1,0);
		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(base.nseg,0,0),0,gbl->vdt.length(secondDim)*gbl->vdt.length(secondDim)-1,0);
	}
	
	if (gbl->is_loop) {
		for(m=0;m<tri_mesh::ND;++m)
			for(n=0;n<tri_mesh::ND;++n)
				gbl->vdt(0,m,n) = 0.5*(gbl->vdt(0,m,n) +gbl->vdt(base.nseg+1,m,n));
		gbl->vdt(base.nseg+1,Range::all(),Range::all()) = gbl->vdt(0,Range::all(),Range::all());
	}
	
	FLT jcbi,temp;
	for(indx=0;indx<base.nseg+1;++indx) {
		/* INVERT VERTEX MATRIX */
		jcbi = 1.0/(gbl->vdt(indx,0,0)*gbl->vdt(indx,1,1) -gbl->vdt(indx,0,1)*gbl->vdt(indx,1,0));
		
		temp = gbl->vdt(indx,0,0)*jcbi*gbl->cfl(x.log2p,1);
		gbl->vdt(indx,0,0) = gbl->vdt(indx,1,1)*jcbi*gbl->cfl(x.log2p,0);
		gbl->vdt(indx,1,1) = temp;
		gbl->vdt(indx,0,1) *= -jcbi*gbl->cfl(x.log2p,1);
		gbl->vdt(indx,1,0) *= -jcbi*gbl->cfl(x.log2p,0);
		
		/* DIRECT FORMATION OF vdt^{-1} theta is angle of normal from horizontal */
		//		FLT theta =  100.0*M_PI/180.0;
		//		gbl->vdt(indx,0,0) = -sin(theta);
		//		gbl->vdt(indx,1,1) =  sin(theta);
		//		gbl->vdt(indx,0,1) = cos(theta);
		//		gbl->vdt(indx,1,0) = cos(theta);
	}
	
	/* INVERT SIDE MATRIX */
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx=0;indx<base.nseg;++indx) {
			/* INVERT SIDE MVDT MATRIX */
			jcbi = 1.0/(gbl->sdt(indx,0,0)*gbl->sdt(indx,1,1) -gbl->sdt(indx,0,1)*gbl->sdt(indx,1,0));
			
			temp = gbl->sdt(indx,0,0)*jcbi*gbl->cfl(x.log2p,1);
			gbl->sdt(indx,0,0) = gbl->sdt(indx,1,1)*jcbi*gbl->cfl(x.log2p,0);
			gbl->sdt(indx,1,1) = temp;
			gbl->sdt(indx,0,1) *= -jcbi*gbl->cfl(x.log2p,1);
			gbl->sdt(indx,1,0) *= -jcbi*gbl->cfl(x.log2p,0);
		}
	}

	return;
}

void surface_outflow::init(input_map& inmap,void* gbl_in) {
	hp_deformable_free_pnt::init(inmap,gbl_in);
	
	if (surf->is_master) {
		inmap.getwdefault(base.idprefix +"_contact_angle",contact_angle,90.0);
		contact_angle *= M_PI/180.0;
	}
}

void surface_outflow::element_rsdl(Array<FLT,1> lf) {
	TinyVector<FLT,tri_mesh::ND> tangent, wall_normal;
	surface *surf2 = dynamic_cast<surface *>(surf);

	
	lf = 0.0;
	hp_deformable_free_pnt::element_rsdl(lf);
	
	
	if (surf->is_master) {
		/* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
		/* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
		if (wall_type == curved) {
            /* Calculate normal from geometry definition */
			if (surfbdry == 0) {
				x.ebdry(base.ebdry(1))->bdry_normal(0,-1.0,wall_normal);
			}
			else {
				x.ebdry(base.ebdry(0))->bdry_normal(x.ebdry(base.ebdry(0))->nseg-1,1.0,wall_normal);
			}
		}
        else {
            /* Calculate normal using mesh */
            int bnumwall = base.ebdry(1-surfbdry);
            TinyVector<FLT,tri_mesh::ND> rp;
            if (surfbdry == 0) {
                int sindwall = x.ebdry(bnumwall)->seg(0);
                x.crdtocht1d(sindwall);
                basis::tri(x.log2p)->ptprobe1d(2,&rp(0),&wall_normal(0),-1.0,&x.cht(0,0),MXTM);
            }
            else {
                int sindwall = x.ebdry(bnumwall)->seg(x.ebdry(bnumwall)->nseg-1);
                x.crdtocht1d(sindwall);
                basis::tri(x.log2p)->ptprobe1d(2,&rp(0),&wall_normal(0),1.0,&x.cht(0,0),MXTM);
            }
            FLT length = sqrt(wall_normal(0)*wall_normal(0) +wall_normal(1)*wall_normal(1));
            wall_normal /= length;
            FLT temp = wall_normal(0);
            wall_normal(0) = wall_normal(1);
            wall_normal(1) = -temp;
        }
		
		if (surfbdry == 0) {
			/* Surf-boundary then point then wall (in ccw sense) */
			tangent(0) = wall_normal(0)*sin(contact_angle) +wall_normal(1)*cos(contact_angle);
			tangent(1) = -wall_normal(0)*cos(contact_angle) +wall_normal(1)*sin(contact_angle);
			lf(0) -= RAD(x.pnts(base.pnt)(0))*surf2->gbl->sigma*tangent(0);
			lf(1) -= RAD(x.pnts(base.pnt)(0))*surf2->gbl->sigma*tangent(1);
		}
		else {
			/* Wall then point then Surf-boundary (in ccw sense) */
			tangent(0) = wall_normal(0)*sin(contact_angle) -wall_normal(1)*cos(contact_angle);
			tangent(1) = wall_normal(0)*cos(contact_angle) +wall_normal(1)*sin(contact_angle);
			lf(0) -= RAD(x.pnts(base.pnt)(0))*surf2->gbl->sigma*tangent(0);
			lf(1) -= RAD(x.pnts(base.pnt)(0))*surf2->gbl->sigma*tangent(1);
		}
	}
	
	return;
}
