/*
 *  cd_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/13/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "bdry_cd.h"
#include "myblas.h"
#include "melt_cd.h"

//#define MPDEBUG
//#define DETAILED_DT

using namespace bdry_cd;

void melt_cd::init(input_map& inmap,void* gbl_in) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;
	std::string side_id, matching_block, liquid_block;

	gbl = static_cast<global *>(gbl_in);
	
	if (x.NV == 1) {
		is_master = false; // Liquid boundary is master
		/* Let temperature be continuous */
		inmap[base.idprefix +"_c0_indices"] = "0";
		keyword = base.idprefix + "_hp_typelist";
		inmap[keyword] = "1";
		
		find_matching_boundary_name(inmap, matching_block, side_id);
		liquid_block = matching_block +"_" +side_id;
		
		if (!inmap.get(x.gbl->idprefix + "_rho",gbl->rho_s)) inmap.getwdefault("rho",gbl->rho_s,1.0);
		if (!inmap.get(x.gbl->idprefix + "_cv",gbl->cp_s)) inmap.getwdefault("cv",gbl->cp_s,1.0);
		if (!inmap.get(x.gbl->idprefix + "_conductivity",gbl->k_s)) inmap.getwdefault("conductivity",gbl->k_s,0.0);
		
		if (!inmap.get(matching_block + "_rho",gbl->rho_l)) inmap.getwdefault("rho",gbl->rho_l,1.0);
		if (!inmap.get(matching_block + "_cp",gbl->cp_l)) inmap.getwdefault("cp",gbl->cp_l,1.0);
		FLT mu;
		if (!inmap.get(matching_block +"_mu",mu)) inmap.getwdefault("mu",mu,0.0);
		if (!inmap.get(matching_block + "_conductivity",gbl->k_l)) inmap.getwdefault("conductivity",gbl->k_l,0.7*mu);
	}
	else {
		// This is for the liquid side
		is_master = true;
		/* Let temperature be continuous */
		inmap[base.idprefix +"_c0_indices"] = "2";
		
		/* Make u & v dirichlet B.C.'s */
		keyword = base.idprefix + "_hp_typelist";
		inmap[keyword] = "0 0 1 1";
		liquid_block = base.idprefix;

		if (!inmap.get(matching_block + "_rho",gbl->rho_s)) inmap.getwdefault("rho",gbl->rho_s,1.0);
		if (!inmap.get(matching_block + "_cv",gbl->cp_s)) inmap.getwdefault("cv",gbl->cp_s,1.0);
		if (!inmap.get(matching_block + "_conductivity",gbl->k_s)) inmap.getwdefault("conductivity",gbl->k_s,0.0);
		
		if (!inmap.get(x.gbl->idprefix + "_rho",gbl->rho_l)) inmap.getwdefault("rho",gbl->rho_l,1.0);
		if (!inmap.get(x.gbl->idprefix + "_cp",gbl->cp_l)) inmap.getwdefault("cp",gbl->cp_l,1.0);
		FLT mu;
		if (!inmap.get(x.gbl->idprefix +"_mu",mu)) inmap.getwdefault("mu",mu,0.0);
		if (!inmap.get(x.gbl->idprefix + "_conductivity",gbl->k_l)) inmap.getwdefault("conductivity",gbl->k_l,0.7*mu);
	}
	
	hp_deformable_bdry::init(inmap,gbl_in);
	
	//    gbl->vdt_kinetic.resize(base.maxseg+1);
	//    gbl->sdt_kinetic.resize(base.maxseg);
	
	keyword = liquid_block + "_Lf";
	if (!inmap.get(keyword,gbl->Lf)) {
		*x.gbl->log << "Missing latent heat of fusion " << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	if (!base.is_comm()) {
		keyword = liquid_block + "_cp_s";
		if (!inmap.get(keyword,gbl->cp_s)) {
			*x.gbl->log << "Missing c_p of solid " << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		
		keyword = liquid_block + "_rho_s";
		if (!inmap.get(keyword,gbl->rho_s)) {
			*x.gbl->log << "Missing rho of solid " << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}
	
	inmap.getwdefault(liquid_block + "_Krough",gbl->Krough,0.0);
	inmap.getwdefault(liquid_block + "_Kgt",gbl->Kgt,0.0);
	inmap.getwdefault(liquid_block + "_Ksn",gbl->Ksn,0.0);
	inmap.getwdefault(liquid_block + "_K2Dn",gbl->K2Dn,0.0);
	inmap.getwdefault(liquid_block + "_K2Dn_max",gbl->K2Dn_max,300.0);
	inmap.getwdefault(liquid_block + "_A2Dn",gbl->A2Dn,1.0);
	inmap.getwdefault(liquid_block + "_surge_time",gbl->surge_time,1.0);
	
	FLT angle;
	inmap.getwdefault(liquid_block + "_facet_angle",angle,0.0);
	gbl->facetdir(0) = cos(M_PI*angle/180.0);
	gbl->facetdir(1) = sin(M_PI*angle/180.0);
	
#ifdef SWAP_ROWS
	*x.gbl->log << "SWAP_ROWS is defined\n";
#else
	*x.gbl->log << "SWAP_ROWS is not defined\n";
#endif
	
	return;
}

void melt_cd::setup_preconditioner() {
	int indx,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> nrm;
	FLT h, hsm;
	FLT dttang, dtnorm;
	FLT vslp;
	FLT qmax;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd;
	TinyMatrix<FLT,4,MXGP> res;
	TinyMatrix<FLT,4,MXGP> lf;
	TinyVector<FLT,2> mvel;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	int last_phase, mp_phase;
	TinyVector<FLT,tri_mesh::ND> aPoint;
	aPoint = 0.0;
	const FLT alpha = gbl->k_s/(gbl->rho_s*gbl->cp_s) +gbl->k_l/(gbl->rho_l*gbl->cp_l);
	const FLT Tm = ibc->f(c0_indices.back(), aPoint, 0.0);
	TinyVector<FLT,2> us;
	FLT sign_dir;
	
	if (x.NV > 1) {
		// get solid velocity
		for (int n=0;n<tri_mesh::ND;++n)
			us(n) = ibc->f(n, aPoint, 0.0);
			sign_dir = 1.0;
	}
	else {
		tri_hp_cd& xcd = dynamic_cast<tri_hp_cd&>(x);
#ifdef CONST_A
		us(0) = xcd.gbl->ax;
		us(1) = xcd.gbl->ay;
#else
		for (int n=0;n<tri_mesh::ND;++n)
			us(n) = gbl->a(n, aPoint, 0.0);
#endif
		sign_dir = -1;
	}
	
	/**************************************************/
	/* DETERMINE SURFACE MOVEMENT TIME STEP              */
	/**************************************************/
	gbl->vdt(0) = 0.0;
	
	for(indx=0; indx < base.nseg; ++indx) {
		sind = base.seg(indx);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		
#ifdef DETAILED_DT
		x.crdtocht1d(sind);
		for(int n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
		
		x.ugtouht1d(sind);
		
		dtnorm = 1.0e99;
		dttang = 1.0e99;
		gbl->meshc(indx) = 1.0e99;
		for(int i=0;i<basis::tri(x.log2p)->gpx();++i) {
			nrm(0) =  dcrd(1,i)*2;
			nrm(1) = -dcrd(0,i)*2;
			h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
			
			/* RELATIVE VELOCITY STORED IN MVEL(N)*/
			for(int n=0;n<tri_mesh::ND;++n) {
				mvel(n) = us(n) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef MESH_REF_VEL
				mvel(n) -= x.gbl->mesh_ref_vel(n);
#endif
			}
			
			vslp = fabs(-us(0)*nrm(1)/h +us(1)*nrm(0)/h);
			qmax = mvel(0)*mvel(0)+mvel(1)*mvel(1);
			hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
			
			dttang = MIN(dttang,2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm);
			//dtnorm = MIN(dtnorm,(2.*vslp/hsm +x.gbl->bd(0) +1.5*alpha/(hsm*hsm))*gbl->Lf*(gbl->rho_s +gbl->rho_l));
			
			const FLT alpha = gbl->k_l/(gbl->rho_l*gbl->cp_l);
			dtnorm = MIN(dtnorm,(2.*vslp/hsm +x.gbl->bd(0))*gbl->Lf*gbl->rho_l +1.5*alpha/(hsm*hsm)*gbl->rho_l*gbl->cp_l*Tm);

			/* SET UP DISSIPATIVE COEFFICIENT */
			/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
			/* RESIDUAL HAS DX/2 WEIGHTING */
			/* |a| dx/2 dv/dx  dx/2 dpsi */
			/* |a| dx/2 2/dx dv/dpsi  dpsi */
			/* |a| dv/dpsi  dpsi */
			// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5);/* FAILED IN NATES UPSTREAM SURFACE WAVE CASE */
			//gbl->meshc(indx) = MIN(gbl->meshc(indx),gbl->adis/(h*(vslp/hsm +x.gbl->bd(0)))); /* FAILED IN MOVING UP TESTS */
			gbl->meshc(indx) = MIN(gbl->meshc(indx),gbl->adis/(h*(sqrt(qmax)/hsm +x.gbl->bd(0)))); /* SEEMS THE BEST I'VE GOT */
		}
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
#else
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
		h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
		
		mvel(0) = us(0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
		mvel(1) = us(1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
#ifdef MESH_REF_VEL
		mvel -= x.gbl->mesh_ref_vel;
#endif
		
		nrm *= sign_dir;
		qmax = mvel(0)*mvel(0)+mvel(1)*mvel(1);
		vslp = fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h);
		
		mvel(0) = us(0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
		mvel(1) = us(1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
		
#ifdef MESH_REF_VEL
		mvel -= x.gbl->mesh_ref_vel;
#endif
		
		qmax = MAX(qmax,mvel(0)*mvel(0)+mvel(1)*mvel(1));
		vslp = MAX(vslp,fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h));
		
		hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
		
		dttang = 2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm;
		//dtnorm = (2.*vslp/hsm +x.gbl->bd(0) +1.5*alpha/(hsm*hsm))*gbl->Lf*(gbl->rho_s +gbl->rho_l);

		const FLT alpha = gbl->k_l/(gbl->rho_l*gbl->cp_l);
		dtnorm = (2.*vslp/hsm +x.gbl->bd(0))*gbl->Lf*gbl->rho_l +1.5*alpha/(hsm*hsm)*gbl->rho_l*gbl->cp_l*Tm;
		
		
		/* SET UP DISSIPATIVE COEFFICIENT */
		/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
		/* RESIDUAL HAS DX/2 WEIGHTING */
		/* |a| dx/2 dv/dx  dx/2 dpsi */
		/* |a| dx/2 2/dx dv/dpsi  dpsi */
		/* |a| dv/dpsi  dpsi */
		gbl->meshc(indx) = gbl->adis/(h*(sqrt(qmax)/hsm +x.gbl->bd(0)));
#endif
		dtnorm *= RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));
		nrm *= 0.5;
		
//		*x.gbl->log << dttang << ' ' << dtnorm << std::endl;
		
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
				for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
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
					gbl->ms(indx)(2*m,2*n) = lf(0,n+2);
					gbl->ms(indx)(2*m,2*n+1) = lf(1,n+2);
					gbl->ms(indx)(2*m+1,2*n) = lf(2,n+2);
					gbl->ms(indx)(2*m+1,2*n+1) = lf(3,n+2);
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
			GETRF(2*lsm,2*lsm,&gbl->ms(indx)(0,0),2*MAXP,&gbl->ipiv(indx)(0),info);
			if (info != 0) {
				*x.gbl->log << "DGETRF FAILED IN SIDE MODE PRECONDITIONER\n");
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
		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vdt(0,0,0),0,3,0);
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vdt(base.nseg,0,0),0,3,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		
		last_phase = true;
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(0,0,0),0,3,0);
		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(base.nseg,0,0),0,3,0);
	}
	
	if (gbl->is_loop) {
		for(int m=0;m<tri_mesh::ND;++m)
			for(int n=0;n<tri_mesh::ND;++n)
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


		
		// FOR TESTING
//		*x.gbl->log << indx << ' ' << gbl->vdt(indx,0,0) << ' ' << gbl->vdt(indx,0,1) << ' ' << gbl->vdt(indx,1,0) << ' ' << gbl->vdt(indx,1,1) << ' ' << std::endl;
//		gbl->vdt(indx,Range::all(),Range::all()) = 0.0;
//		gbl->vdt(indx,0,0) = 1.0;
//		gbl->vdt(indx,1,1) = 1.0;
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
			

			
			// FOR TESTING
//			*x.gbl->log << indx << ' ' << gbl->sdt(indx,0,0) << ' ' << gbl->sdt(indx,0,1) << ' ' << gbl->sdt(indx,1,0) << ' ' << gbl->sdt(indx,1,1) << ' ' << std::endl;
//			gbl->sdt(indx,Range::all(),Range::all()) = 0.0;
//			gbl->sdt(indx,0,0) = 1.0;
//			gbl->sdt(indx,1,1) = 1.0;
		}
	}
	
	//	/**************************************************/
	//	/* DETERMINE TEMPRATURE UPDATE TIME STEP          */
	//	/**************************************************/
	//	gbl->vdt_kinetic(0) = 0.0;
	//
	//	for(indx=0; indx < base.nseg; ++indx) {
	//		sind = base.seg(indx);
	//		v0 = x.seg(sind).pnt(0);
	//		v1 = x.seg(sind).pnt(1);
	//
	//		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
	//		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
	//		nrm *= 0.5;
	//		h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
	//		dtnorm = 1.0;
	//		dtnorm *= RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));
	//
	//		gbl->vdt_kinetic(indx) += dtnorm*h*basis::tri(x.log2p)->vdiag1d();
	//		if (basis::tri(x.log2p)->sm()) {
	//			gbl->sdt_kinetic(indx) = dtnorm*h;
	//		}
	//	}
	//
	//	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
	//		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vdt_kinetic(0),0,0,0);
	//		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vdt_kinetic(base.nseg),0,0,0);
	//		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
	//		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
	//
	//		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
	//		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
	//
	//		last_phase = true;
	//		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
	//		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
	//		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt_kinetic(0),0,0,0);
	//		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt_kinetic(base.nseg),0,0,0);
	//	}
	//
	//	if (gbl->is_loop) {
	//		gbl->vdt_kinetic(0) = 0.5*(gbl->vdt_kinetic(0) +gbl->vdt_kinetic(base.nseg+1));
	//		gbl->vdt_kinetic(base.nseg+1) = gbl->vdt_kinetic(0);
	//	}
	//
	//	for(indx=0;indx<base.nseg+1;++indx) {
	//		/* INVERT VERTEX MATRIX */
	//		gbl->vdt_kinetic(indx) = 1./gbl->vdt_kinetic(indx);
	//	}
	//
	//	/* INVERT SIDE MATRIX */
	//	if (basis::tri(x.log2p)->sm() > 0) {
	//		for(indx=0;indx<base.nseg;++indx) {
	//			/* INVERT SIDE MVDT MATRIX */
	//			gbl->sdt_kinetic(indx) = 1./gbl->sdt_kinetic(indx);
	//		}
	//	}
	
	
	return;
}


#if defined(petsc) && defined(SWAP_ROWS)
void melt_cd::non_sparse(Array<int,1> &nnzero) {
	
	hp_deformable_bdry::non_sparse(nnzero);
	
	/* add degrees of freedom to allow swapping of rows */
	const int sm=basis::tri(x.log2p)->sm();
	const int im=basis::tri(x.log2p)->im();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
	// On slave side make more space so vertex equations can be swapped
	for (int i=0;i<base.nseg;++i) {
		int sind = base.seg(i);
		if (!is_master) {
			int pind = x.seg(sind).pnt(0);
			nnzero(Range(pind*vdofs+x.NV,(pind+1)*vdofs-1)) += NV*sm;
			
			pind = x.seg(sind).pnt(1);
			nnzero(Range(pind*vdofs+x.NV,(pind+1)*vdofs-1)) += NV*sm;
		}
	}

	if(x.sm0 > 0) {
		if (is_master) {
			nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) = 3*vdofs +3*x.NV*sm +x.NV*im +tri_mesh::ND*sm;
		}
	}
}

void melt_cd::non_sparse_rcv(Array<int, 1> &nnzero, Array<int, 1> &nnzero_mpi) {
	hp_deformable_bdry::non_sparse_rcv(nnzero,nnzero_mpi);
	
	/* add degrees of freedom to allow swapping of rows */
	const int sm=basis::tri(x.log2p)->sm();
	const int im=basis::tri(x.log2p)->im();
	const int vdofs = 1 +(x.mmovement == tri_hp::coupled_deformable)*x.ND;  // equations in solid have one variable
	
	if(x.sm0 > 0) {
		if (is_master && base.is_comm()) {
			nnzero_mpi(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) = 3*vdofs +3*sm +im +tri_mesh::ND*sm;
		}
	}
}


void melt_cd::petsc_make_1D_rsdl_vector(Array<double,1> res) {
	const int sm = basis::tri(x.log2p)->sm();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
	if (is_master) {
		/* Swap kinetic and energy side residuals */
		for(int j=0;j<base.nseg;++j) {
			int ind1 = x.npnt*vdofs +base.seg(j)*sm*x.NV +x.NV-2;  // Index of heat equation
			for(int m=0;m<sm;++m) {
				FLT temp = res(ind1);
				res(ind1) = gbl->sres(j,m,1);
				gbl->sres(j,m,1) = temp;
				ind1 += x.NV;
			}
		}
	}
	

	/* Swap kinetic and energy vertex residuals */
	int i = 0;
	int sind;
	const int heatindex = c0_indices.back();
	do {
		sind = base.seg(i);
		int row = x.seg(sind).pnt(0)*vdofs;
		double temp = res(row+heatindex);
		res(row+heatindex) = res(row+vdofs-1);
		res(row+vdofs-1) = temp;
	} while (++i < base.nseg);
	int row = x.seg(sind).pnt(1)*vdofs;
	double temp = res(row+heatindex);
	res(row+heatindex) = res(row+vdofs-1);
	res(row+vdofs-1) = temp;

	
	/* now rotate based on normal */
	hp_deformable_bdry::petsc_make_1D_rsdl_vector(res);
}

void melt_cd::petsc_premultiply_jacobian() {
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int heatindex = c0_indices.back();
	
	/* Swap side equations */
	if (is_master) {
		int row = jacobian_start; // Index of motion equations
		const int sm = basis::tri(x.log2p)->sm();
		for(int j=0;j<base.nseg;++j) {
			int row1 = x.npnt*vdofs +base.seg(j)*sm*x.NV +heatindex;  // Index of heat equation
			for(int m=0;m<sm;++m) {
				// Match sparseness of heat equation row by filling in zeros in normal and tangential rows
				x.J.match_patterns(row, row1);
				x.J.match_patterns(row+1, row1);
				x.J.swap_rows(row+1, row1);
				x.J_mpi.match_patterns(row, row1);
				x.J_mpi.match_patterns(row+1, row1);
				x.J_mpi.swap_rows(row+1, row1);
				row += tri_mesh::ND;
				row1 += x.NV;
			}
		}
	}
	
	/* Swap kinetic and energy vertex rows */
#ifdef ONE_SIDED
	if (!is_master) {
		int i = 0;
		int sind;
		do {
			sind = base.seg(i);
			int row = x.seg(sind).pnt(0)*vdofs;
			x.J.match_patterns(row+heatindex, row+vdofs-1);
			x.J.match_patterns(row+heatindex, row+vdofs-2);
			x.J_mpi.match_patterns(row+heatindex, row+vdofs-1);
			x.J_mpi.match_patterns(row+heatindex, row+vdofs-2);
		} while (++i < base.nseg);
		int row = x.seg(sind).pnt(1)*vdofs;
		x.J.match_patterns(row+heatindex, row+vdofs-1);
		x.J.match_patterns(row+heatindex, row+vdofs-2);
		x.J_mpi.match_patterns(row+heatindex, row+vdofs-1);
		x.J_mpi.match_patterns(row+heatindex, row+vdofs-2);
	}
	else {
#endif
		int i = 0;
		int sind;
		do {
			sind = base.seg(i);
			int row = x.seg(sind).pnt(0)*vdofs;
			x.J.match_patterns(row+heatindex, row+vdofs-1);
			x.J.match_patterns(row+heatindex, row+vdofs-2);
			x.J.swap_rows(row+heatindex, row+vdofs-1);
			x.J_mpi.match_patterns(row+heatindex, row+vdofs-1);
			x.J_mpi.match_patterns(row+heatindex, row+vdofs-2);
			x.J_mpi.swap_rows(row+heatindex, row+vdofs-1);
		} while (++i < base.nseg);
		int row = x.seg(sind).pnt(1)*vdofs;
		x.J.match_patterns(row+heatindex, row+vdofs-1);
		x.J.match_patterns(row+heatindex, row+vdofs-2);
		x.J.swap_rows(row+heatindex, row+vdofs-1);
		x.J_mpi.match_patterns(row+heatindex, row+vdofs-1);
		x.J_mpi.match_patterns(row+heatindex, row+vdofs-2);
		x.J_mpi.swap_rows(row+heatindex, row+vdofs-1);
#ifdef ONE_SIDED
	}
#endif
	/* now rotate based on normal */
	hp_deformable_bdry::petsc_premultiply_jacobian();
}
#endif



