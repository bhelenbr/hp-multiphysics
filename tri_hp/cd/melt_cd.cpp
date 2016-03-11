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
	
	hp_coupled_bdry::init(inmap,gbl_in);
#ifdef PRECONDITION
	gbl->field_is_coupled = true;
#endif
	
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
#ifdef OLDKINETICS
	inmap.getwdefault(liquid_block + "_K2Dn_max",gbl->K2Dn_max,300.0);
#else
	inmap.getwdefault(liquid_block + "_K2Dn_DT_min",gbl->K2Dn_DT_min,1.0);
#endif
	inmap.getwdefault(liquid_block + "_A2Dn",gbl->A2Dn,1.0);
	inmap.getwdefault(liquid_block + "_surge_time",gbl->surge_time,1.0);
	
	FLT angle;
	inmap.getwdefault(liquid_block + "_facet_angle",angle,0.0);
	gbl->facetdir(0) = cos(M_PI*angle/180.0);
	gbl->facetdir(1) = sin(M_PI*angle/180.0);
	
	return;
}

#ifdef OLDKINETICS
FLT melt_cd::calculate_kinetic_coefficients(FLT DT,FLT sint) {
	FLT K;
	FLT K2Dn_exp = gbl->K2Dn_max*gbl->K2Dn/(exp(-gbl->A2Dn/(abs(DT) +100.*EPSILON))*gbl->K2Dn_max +gbl->K2Dn);

#ifdef TWOFACETS
	// Hack to get other facet angle
	// Theta is defined as angle between outward liquid normal and facet direction (outward from solid)
	// This inconsistency makes counterclockwise rotations negative
	const int p = 2;
	FLT theta = asin(sint);
	theta -= 70.0*M_PI/180.0;
	K = gbl->Krough*pow(1. + 1./(pow(sint/gbl->Ksn,p) +pow(1./K2Dn_exp,p)) +1./(pow(fabs(sin(theta))/gbl->Ksn,p) +pow(1./K2Dn_exp,p)),1.0/p);
	return(K);
#endif
	
	if (sint == 0.0) {
		K = gbl->Krough*K2Dn_exp;
		
		// Step Oscillation */
		// FLT K1 = gbl->Krough*(1.); // + 2.*gbl->Ksn/((-sint +fabs(sint)) +EPSILON));
		// K = 0.5*(K1-K)*erfc(cos(M_PI*x.gbl->time/gbl->surge_time)*4.0) +K;
		
		//        /* Step Change */
		//        if (M_PI*x.gbl->time/gbl->surge_time < 1.0) {
		//            K = 0.5*(K1-K)*erfc(cos(M_PI*x.gbl->time/gbl->surge_time)*4.0) +K;
		//        }
		//        else {
		//            K = K1;
		//        }
		//
		//		/* For an unsteady perturbation */
		//		const FLT T = 1.5e-3;
		//		if (x.gbl->time < T) {
		//			K = K*(1.0+0.1*sin(2*M_PI*x.gbl->time/T));
		//		}
	}
	else {
		K = gbl->Krough*(1. + 2.*gbl->Ksn/((-sint +fabs(sint)) +EPSILON));
	}
	
	return(K);
}
#else
FLT melt_cd::calculate_kinetic_coefficients(FLT DT,FLT sint) {
	FLT K;
	const int p = 2;
	
	FLT K2Dn_exp = gbl->K2Dn*exp(gbl->A2Dn/(max(abs(DT),gbl->K2Dn_DT_min)));
	
	if (sint == 0.0) {
		K = K2Dn_exp;
	}
	else {
		K = pow(pow(gbl->Krough,p) + pow(gbl->Ksn/(fabs(sint) +EPSILON),p),1.0/p);
	}
	
	return(K);
}
#endif

void melt_cd::output(const std::string& filename, tri_hp::filetype typ,int tlvl) {
	
	hp_coupled_bdry::output(filename,typ,tlvl);
	
	TinyVector<FLT,tri_mesh::ND> norm, aloc;
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	
	switch (typ) {
		case(tri_hp::tecplot): {
			if (!report_flag) break;
			
			std::string fname;
			fname = filename +"kinetics_" +base.idprefix +".dat";
			std::ofstream fout;
			fout.open(fname.c_str());
			
			int indx = 0;
			FLT s = 0.0;
			do {
				int sind = base.seg(indx);
				
				x.crdtocht1d(sind);
				for(int n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
				
				x.ugtouht1d(sind);
				for(int n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));
				
				for(int i=0;i<basis::tri(x.log2p)->gpx();++i) {
					norm(0) =  dcrd(1,i);
					norm(1) = -dcrd(0,i);
					jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
					s += basis::tri(x.log2p)->wtx(i)*jcb;
					FLT sint = (-gbl->facetdir(0)*norm(1) +gbl->facetdir(1)*norm(0))/jcb;
					aloc(0) = crd(0,i); aloc(1) = crd(1,i);
					FLT DT = ibc->f(c0_indices[0], aloc, x.gbl->time) -u(c0_indices[0])(i);
					FLT K = calculate_kinetic_coefficients(DT,sint);
					
					fout << s << ' ' << crd(0,i) << ' ' << crd(1,i) << ' ' << K << ' ' << DT << ' ' << sint << std::endl;
					
				}
			} while (++indx < base.nseg);
			
			/* Calculate exactly what is happening at triple junction */
			const int seg = base.nseg-1;
			const int sind = base.seg(seg);
			const int v0 = x.seg(sind).pnt(1);
			TinyVector<FLT,4> up;
			basis::tri(x.log2p)->ptprobe1d(4,up.data(),1.0,&x.uht(0)(0),MXTM);
			FLT DT = ibc->f(c0_indices[0], x.pnts(v0), x.gbl->time) -up(c0_indices[0]);
			FLT K = calculate_kinetic_coefficients(DT,0.0);
			fout << s << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << ' ' << K << ' ' << DT << ' ' << 0.0 << std::endl;
			fout.close();
			
			break;
		}
		default:
			break;
	}
}

void melt_cd::setup_preconditioner() {
	int indx,last_phase, mp_phase;
	TinyVector<FLT,tri_mesh::ND> nrm,mvel,us,aPoint,dXdXi,dXdEta,dX;
	FLT h,hsm,vslp,vnrm,dttang,dtnorm,sign_dir;
	const int Tvar = c0_indices[0];

	aPoint = 0.0;
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
	gbl->vdt(0,Range::all(),Range::all()) = 0.0;
	for(indx=0; indx < base.nseg; ++indx) {
		const int sind = base.seg(indx);
		const int v0 = x.seg(sind).pnt(0);
		const int v1 = x.seg(sind).pnt(1);
		const int tind = x.seg(sind).tri(0);
		int snum;
		for(snum=0;snum<3;++snum) {
			if (x.tri(tind).seg(snum) == sind) goto found;
		}
		sim::abort(__LINE__, __FILE__, x.gbl->log);
		found: const int v2 = x.tri(tind).pnt(snum);

		const FLT Tm = 0.5*(x.ug.v(v0,Tvar) +x.ug.v(v1,Tvar));
		for(int n=0;n<tri_mesh::ND;++n) {
			dXdXi(n) = x.pnts(v1)(n)-x.pnts(v0)(n);
			dXdEta(n) = x.pnts(v2)(n)-x.pnts(v0)(n);
			aPoint(n) = 0.5*(x.pnts(v1)(n)+x.pnts(v0)(n));
		}
		nrm(0) =  dXdXi(1);
		nrm(1) = -dXdXi(0);
		h = nrm(0)*nrm(0) +nrm(1)*nrm(1);
		const FLT Xi = (dXdEta(0)*dXdXi(0) +dXdEta(1)*dXdXi(1))/h;
		const FLT dT = x.ug.v(v2,Tvar) -(Xi*x.ug.v(v0,Tvar)+(1-Xi)*x.ug.v(v1,Tvar));
		for(int n=0;n<tri_mesh::ND;++n)
			dX(n) = x.pnts(v2)(n) -(Xi*x.pnts(v0)(n)+(1-Xi)*x.pnts(v1)(n));
		const FLT dXnorm =sqrt(dX(0)*dX(0) +dX(1)*dX(1));
		const FLT dTdXnorm = dT/dXnorm;
		h = sqrt(h);
		
		mvel(0) = us(0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
		mvel(1) = us(1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
#ifdef MESH_REF_VEL
		mvel -= x.gbl->mesh_ref_vel;
#endif
		
		nrm *= sign_dir;
		FLT qmax = mvel(0)*mvel(0)+mvel(1)*mvel(1); // TEMPO
		vslp = fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h);
		vnrm = mvel(0)*nrm(0) +mvel(1)*nrm(1);
		mvel(0) = us(0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
		mvel(1) = us(1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
		
#ifdef MESH_REF_VEL
		mvel -= x.gbl->mesh_ref_vel;
#endif
		qmax = MAX(qmax,mvel(0)*mvel(0)+mvel(1)*mvel(1));
		vnrm = MAX(vnrm,mvel(0)*nrm(0)+mvel(1)*nrm(1));
		vslp = MAX(vslp,fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h));
		
		hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
		
		dttang = 2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm;
		// Precondition Interface Rows from 1D Matlab
		// [diagTs+diagTl,dti*rhos*Lf +(rhol*cl-rhos*cs)*dti*T -Kl/dxl*(dT/dx) -Ks/dxs*(dT/dx)]
		// [-dKc*(vg-xv(N/2+1)) +1.0, -Kc*dti]
		dtnorm = (vslp/hsm +x.gbl->bd(0))*gbl->Lf*gbl->rho_s +(gbl->rho_l*gbl->cp_l -gbl->rho_s*gbl->cp_s)*x.gbl->bd(0)*Tm -gbl->k_s/dXnorm*dTdXnorm;
		
		/* This is from Weinsteins expression */
		const FLT sint = (-gbl->facetdir(0)*nrm(1) +gbl->facetdir(1)*nrm(0))/h;
		const FLT DTcool = ibc->f(c0_indices[0], aPoint, x.gbl->time) -Tm;
		const FLT Kc = calculate_kinetic_coefficients(DTcool, sint);
		// Numerically calculate derivatives
		const FLT dKcdT = (calculate_kinetic_coefficients(DTcool +1.0e-6*Tm, sint) -Kc)/(1.0e-6*Tm);
		// const FLT dKcDsint = (calculate_kinetic_coefficients(DTcool, sint +1.0e-6*gbl->Ksn/(gbl->Krough +FLT_EPSILON)) -Kc)/(1.0e-6*gbl->Ksn/(gbl->Krough +FLT_EPSILON));
		
		
		dtnorm *= RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));
		gbl->meshc(indx) = gbl->adis/(h*(vslp/hsm +x.gbl->bd(0) +gbl->k_s/dXnorm*fabs(dTdXnorm)/(gbl->Lf*gbl->rho_s)));
		gbl->meshc(indx) = gbl->adis/(h*(sqrt(qmax)/hsm +x.gbl->bd(0)));
		
		// Indices are T,x,y = 0,1,2
		
		// Energy Equation
		gbl->vdt(indx+1,Range::all(),Range::all()) = 0.0;
		gbl->vdt(indx+1,0,0) = x.gbl->vprcn(v1,c0_indices[0]);  // Energy equation diaganol;
		gbl->vdt(indx+1,0,1) = dtnorm*nrm(0);
		gbl->vdt(indx+1,0,2) = dtnorm*nrm(1);
		
		// Tangent Equation
		gbl->vdt(indx+1,1,0) = 0.0;
		gbl->vdt(indx+1,1,1) = -dttang*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,1,2) =  dttang*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		
		// Kinetic Equation
		gbl->vdt(indx+1,2,0) = (-dKcdT*vnrm +h)*gbl->rho_l*RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));
		gbl->vdt(indx+1,2,1) = (-Kc*x.gbl->bd(0))*nrm(0)*gbl->rho_l*RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));  // Missing term here for dK/dTheta
		gbl->vdt(indx+1,2,2) = (-Kc*x.gbl->bd(0))*nrm(1)*gbl->rho_l*RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0))); // Missing term here for dK/dTheta
		
		gbl->vdt(indx,Range::all(),Range::all()) += gbl->vdt(indx+1,Range::all(),Range::all());
		gbl->vdt(indx,0,0) = x.gbl->vprcn(v0,c0_indices[0]);
		
#ifdef TEST_WITH_DIAG
		gbl->vdt(indx,Range::all(),Range::all()) = 0.0;
		gbl->vdt(indx,0,0) = 1.0;
		gbl->vdt(indx,1,1) = 1.0;
		gbl->vdt(indx,2,2) = 1.0;
		
		gbl->vdt(indx+1,Range::all(),Range::all()) = 0.0;
		gbl->vdt(indx+1,0,0) = 1.0;
		gbl->vdt(indx+1,1,1) = 1.0;
		gbl->vdt(indx+1,2,2) = 1.0;
#endif

		if (basis::tri(x.log2p)->sm()) {
			
			gbl->sdt(indx,0,0) = x.gbl->sprcn(sind,c0_indices[0]); // Energy equation diagonal for side modes;
			gbl->sdt(indx,0,1) = dtnorm*nrm(0);
			gbl->sdt(indx,0,2) = dtnorm*nrm(1);
			
			gbl->sdt(indx,1,0) = 0.0;
			gbl->sdt(indx,1,1) = -dttang*nrm(1);
			gbl->sdt(indx,1,2) =  dttang*nrm(0);
			
			// Kinetic Equation
			gbl->sdt(indx,2,0) = (-dKcdT*vnrm +h)*gbl->rho_l*RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));
			gbl->sdt(indx,2,1) = (-Kc*x.gbl->bd(0))*nrm(0)*gbl->rho_l*RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));  // Missing term here for dK/dTheta
			gbl->sdt(indx,2,2) = (-Kc*x.gbl->bd(0))*nrm(1)*gbl->rho_l*RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0))); // Missing term here for dK/dTheta
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
	
	for(indx=0;indx<base.nseg+1;++indx) {
		/* INVERT 3x3 VERTEX MATRICES */
		int info;
		GETRF(3,3,&gbl->vdt(indx,0,0),x.NV+NV,&gbl->vpiv(indx,0),info);
		if (info != 0) {
			*x.gbl->log << "DGETRF FAILED IN VERTEX MODE PRECONDITIONER\n";
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}
	
	/* INVERT SIDE MATRIX */
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx=0;indx<base.nseg;++indx) {
			/* INVERT 3x3 SIDE MATRIX */
			int info;
			GETRF(3,3,&gbl->sdt(indx,0,0),x.NV+NV,&gbl->spiv(indx,0),info);
			if (info != 0) {
				*x.gbl->log << "DGETRF FAILED IN SIDE MODE PRECONDITIONER\n";
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
		}
	}
	
#ifndef TEST_WITH_DIAG
	// Set Field Preconditioners to 1
	int i = 0;
	int sind;
	do {
		sind = base.seg(i);
		int pnt = x.seg(sind).pnt(0);
		x.r_tri_mesh::gbl->diag(pnt) = 1.0;
		for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n)
			x.gbl->vprcn(pnt,*n) = 1.0;
		x.gbl->vprcn(pnt,c0_indices[0]) = 1.0;
	} while (++i < base.nseg);
	int pnt = x.seg(sind).pnt(1);
	x.r_tri_mesh::gbl->diag(pnt) = 1.0;
	for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n)
		x.gbl->vprcn(pnt,*n) = 1.0;
	x.gbl->vprcn(pnt,c0_indices[0]) = 1.0;
#endif
	
	return;
}

void melt_cd::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
	TinyVector<FLT,tri_mesh::ND> aPoint,us;
	const int Tindx = c0_indices[0];
	int sign_dir;
	
#ifndef SYMMETRIC
	if (!is_master) return;
#endif
	
	aPoint = 0.0;
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
	
	
	int i,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp, aloc, anorm, amv;
	Array<FLT,1> ubar(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,8,MXGP> res;
	
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
		
		Array<FLT,1> au(x.NV), flx(x.NV);
		aloc(0) = crd(0,i) ; aloc(1) = crd(1,i);
		
		amv(0) = (x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i)));
		amv(1) = (x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)));
#ifdef MESH_REF_VEL
		amv(0) += x.gbl->mesh_ref_vel(0);
		amv(1) += x.gbl->mesh_ref_vel(1);
#endif
		anorm(0)= norm(0)/jcb; anorm(1) = norm(1)/jcb;
		
		/* This is from Weinsteins expression */
		FLT sint = sign_dir*(-gbl->facetdir(0)*anorm(1) +gbl->facetdir(1)*anorm(0));
		FLT DT = ibc->f(Tindx, aloc, x.gbl->time) -u(Tindx)(i);
		FLT K = calculate_kinetic_coefficients(DT,sint);
		
		/* Additional heat fluxes for one-sided calculation or transparent w/radiation */
		flux(au,aloc,amv,anorm,jcb,flx);
		
		/* Calculate relative velocity */
		amv = us-amv;
		
		/* TANGENTIAL SPACING */
		res(0,i) = -ksprg(indx)*sign_dir*jcb;
		/* NORMAL FLUX */
		res(1,i) = RAD(crd(0,i))*gbl->rho_s*sign_dir*(amv(0)*norm(0) +amv(1)*norm(1));
		/* Kinetic equation for surface temperature */
		res(2,i) = RAD(crd(0,i))*gbl->rho_s*(-DT)*jcb +K*res(1,i);
		/* Heat source */
		res(3,i) =  -(is_master)*gbl->Lf*res(1,i) +(!base.is_comm())*(RAD(crd(0,i))*flx(Tindx)*jcb +gbl->rho_s*gbl->cp_s*u(Tindx)(i)*res(1,i));
		
		/* UPWINDING BASED ON TANGENTIAL VELOCITY (not used) */
		//		res(4,i) = -res(3,i)*(-norm(1)*amv(0) +norm(0)*amv(1))/jcb*gbl->meshc(indx);
		
		/* UPWINDING KINETIC EQUATION BASED ON TANGENTIAL VELOCITY */
		res(4,i) = -res(2,i)*sign_dir*(-norm(1)*amv(0) +norm(0)*amv(1))/jcb*gbl->meshc(indx);
		
		/* This is a linearization of Gibbs-Thomson for a horizontal surface */
		res(4,i) += RAD(crd(0,i))*gbl->Kgt*anorm(0);
	}
	lf = 0.0;
	
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */
	basis::tri(x.log2p)->intgrt1d(&lf(Tindx)(0),&res(3,0)); // heat flux
	//	basis::tri(x.log2p)->intgrtx1d(&lf(Tindx)(0),&res(4,0)); // surface energy balance upwinded (not used)
	
	if (x.NV > 1) basis::tri(x.log2p)->intgrt1d(&lf(x.NV-1)(0),&res(1,0)); // mass flux
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV)(0),&res(0,0)); // tangent
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+1)(0),&res(2,0)); // kinetic equation
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV+1)(0),&res(4,0)); // kinetic equation upwinded
	//	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV-1)(0),&res(4,0)); // mass flux upwinded
	return;
}

void melt_facet_pt2::element_rsdl(Array<FLT,1> lf) {
	
	const int Tindx = c0_indices[0];
	
	hp_deformable_free_pnt::element_rsdl(lf);
	
#ifndef SYMMETRIC
	if (!surface->is_master) return;
#endif
	
	TinyVector<FLT,tri_mesh::ND> aPoint,us;
	aPoint = 0.0;
	if (x.NV > 1) {
		// get solid velocity
		for (int n=0;n<tri_mesh::ND;++n)
			us(n) = ibc->f(n, aPoint, 0.0);
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
	}

	
	Array<FLT,1> u(x.NV);
	TinyVector<FLT, 2> xp, dxpdpsi, mvel, anorm;
	
	int endpt,seg,sind;
	const int bnum = base.ebdry(surfbdry);
	if (surfbdry == 0) {
		seg = x.ebdry(bnum)->nseg-1;
		sind = x.ebdry(bnum)->seg(seg);
		endpt = x.ebdry(bnum)->nseg;
	}
	else {
		seg = 0;
		sind = x.ebdry(bnum)->seg(seg);
		endpt = 0;
	}
	
	x.crdtocht1d(sind);
	x.ugtouht1d(sind);
	
	/* Calculate Temperature from velocity of point & 2D nucleation coefficient only */
	/* Assuming growth at facet angle */
	int v0 = base.pnt;
	basis::tri(x.log2p)->ptprobe1d(tri_mesh::ND,xp.data(),dxpdpsi.data(),1.0-2.*surfbdry,&x.cht(0,0),MXTM);
	basis::tri(x.log2p)->ptprobe1d(x.NV,u.data(),1.0-2.*surfbdry,&x.uht(0)(0),MXTM);
	
	/* RELATIVE VELOCITY STORED IN MVEL(N)*/
	for(int n=0;n<tri_mesh::ND;++n) {
		mvel(n) = us(n) -x.gbl->bd(0)*(xp(n) -x.vrtxbd(1)(v0)(n));
#ifdef MESH_REF_VEL
		mvel(n) -= x.gbl->mesh_ref_vel(n);
#endif
	}
	
	FLT jcb = sqrt(dxpdpsi(0)*dxpdpsi(0)+dxpdpsi(1)*dxpdpsi(1));
	FLT DT = surface->ibc->f(Tindx, xp, x.gbl->time) -u(Tindx);
	bdry_cd::melt_cd *surf1 = dynamic_cast<bdry_cd::melt_cd *>(surface);
	
	/* This is to allow the general expression at the triple point */
	//     anorm(0)= dxpdpsi(1)/jcb; anorm(1) = -dxpdpsi(0)/jcb;
	//     FLT sint = -surf1->gbl->facetdir(0)*anorm(1) +surf1->gbl->facetdir(1)*anorm(0);
	//     FLT K = surf1->calculate_kinetic_coefficients(DT,sint);
	
	FLT K = surf1->calculate_kinetic_coefficients(DT,0.0);
	FLT res1 = jcb*RAD(xp(0))*surf1->gbl->rho_s*(mvel(0)*surf1->gbl->facetdir(0) +mvel(1)*surf1->gbl->facetdir(1));
	/* Kinetic equation for surface temperature */
	lf(x.NV+1) = RAD(xp(0))*surf1->gbl->rho_s*(-DT)*jcb +res1*K;  // -gbl->K_gt*kappa?;
}



/* Routine to add surface tension stress or endpoint movement residual */
void melt_facet_pt2::rsdl(int stage) {
	
	r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
	r_gbl->res(base.pnt)(0) = 0.0;
	r_gbl->res(base.pnt)(1) = 0.0;
	
#ifndef SYMMETRIC
	if (!surface->is_master) return;
#endif
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	Array<FLT,1> lf(vdofs);
	
	int endpt,seg,sind;
	const int bnum = base.ebdry(surfbdry);
	if (surfbdry == 0) {
		seg = x.ebdry(bnum)->nseg-1;
		sind = x.ebdry(bnum)->seg(seg);
		endpt = x.ebdry(bnum)->nseg;
	}
	else {
		seg = 0;
		sind = x.ebdry(bnum)->seg(seg);
		endpt = 0;
	}
	
	element_rsdl(lf);
	
	x.gbl->res.v(base.pnt,Range::all()) += lf(Range(0,x.NV-1));
	
	// FIXME: These should be deleted eventually
	surface->gbl->vres(endpt,0) = lf(x.NV);
	surface->gbl->vres(endpt,1) = lf(x.NV+1);
	
	r_gbl->res(base.pnt)(0) = lf(x.NV);
	r_gbl->res(base.pnt)(1) = lf(x.NV+1);
	
}

#ifdef petsc

void melt_facet_pt2::petsc_jacobian() {
	
#ifndef SYMMETRIC
	if (!surface->is_master) return;
#endif
	
	const int sm = basis::tri(x.log2p)->sm();
	const int vdofs = x.NV+x.ND;
	int nvars = vdofs*(sm+2);
	Array<FLT,2> K(vdofs,nvars);
	Array<FLT,1> lf(vdofs), Rbar(vdofs);
	Array<int,1> rows(vdofs), cols(nvars);
	
	/* Jacobian for row determine x position */
	/* Variables effecting dT/dx are strictly those on the edge */
	
	int endpt,seg,sind;
	const int bnum = base.ebdry(surfbdry);
	if (surfbdry == 0) {
		seg = x.ebdry(bnum)->nseg-1;
		sind = x.ebdry(bnum)->seg(seg);
		endpt = x.ebdry(bnum)->nseg;
	}
	else {
		seg = 0;
		sind = x.ebdry(bnum)->seg(seg);
		endpt = 0;
	}
	
	element_rsdl(lf);
	Rbar = lf;
	
	Array<FLT,1> dw(x.NV);
	dw = 0.0;
	for(int i=0;i<2;++i)
		for(int n=0;n<x.NV;++n)
			dw(n) += fabs(x.uht(n)(i));
	
	dw *= eps_r;
	dw += eps_a;
	FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;
	
	/* Numerically create Jacobian */
	int kcol = 0;
	for(int mode = 0; mode < 2; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			x.ug.v(x.seg(sind).pnt(mode),var) += dw(var);
			element_rsdl(lf);
			
			int krow = 0;
			for(int n=0;n<vdofs;++n)
				K(krow++,kcol) = (lf(n)-Rbar(n))/dw(var);
			++kcol;
			
			x.ug.v(x.seg(sind).pnt(mode),var) -= dw(var);
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var){
			x.pnts(x.seg(sind).pnt(mode))(var) += dx;
			
			element_rsdl(lf);
			
			int krow = 0;
			for(int n=0;n<vdofs;++n)
				K(krow++,kcol) = (lf(n)-Rbar(n))/dx;
			++kcol;
			
			x.pnts(x.seg(sind).pnt(mode))(var) -= dx;
			
		}
	}
	
	for(int mode = 0; mode < sm; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			x.ug.s(sind,mode,var) += dw(var);
			
			element_rsdl(lf);
			
			int krow = 0;
			for(int n=0;n<vdofs;++n)
				K(krow++,kcol) = (lf(n)-Rbar(n))/dw(var);
			++kcol;
			
			x.ug.s(sind,mode,var) -= dw(var);
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var){
			x.hp_ebdry(bnum)->crv(seg,mode)(var) += dx;
			
			element_rsdl(lf);
			
			int krow = 0;
			for(int n=0;n<vdofs;++n)
				K(krow++,kcol) = (lf(n)-Rbar(n))/dx;
			++kcol;
			
			x.hp_ebdry(bnum)->crv(seg,mode)(var) -= dx;
		}
	}
	
	/* CREATE GLOBAL NUMBERING LIST */
	int ind = 0;
	for(int mode = 0; mode < 2; ++mode) {
		int gindx = vdofs*x.seg(sind).pnt(mode);
		for(int var = 0; var < vdofs; ++var)
			cols(ind++) = gindx++;
	}
	
	int gindxNV = x.npnt*vdofs +x.NV*sind*sm;
	int gindxND = surface->jacobian_start +seg*tri_mesh::ND*sm;
	for(int mode = 0; mode < sm; ++mode) {
		for(int var = 0; var < x.NV; ++var)
			cols(ind++) = gindxNV++;
		
		for(int var = 0; var < tri_mesh::ND; ++var)
			cols(ind++) = gindxND++;
	}
	
#ifdef MY_SPARSE
	// ZERO TANGENT ROW & KINEMATIC T ROW
	for(int n=0;n<tri_mesh::ND;++n)
		rows(n) = base.pnt*vdofs +x.NV +n;
	x.J.zero_rows(tri_mesh::ND,rows);
	x.J_mpi.zero_rows(tri_mesh::ND,rows);
	
	for(int n=0;n<vdofs;++n)
		rows(n) = base.pnt*vdofs +n;
	x.J.add_values(vdofs,rows,nvars,cols,K);
#else
	MatSetValuesLocal(x.petsc_J,vdofs,rows.data(),nvars,cols.data(),K.data(),ADD_VALUES);
#endif
}
#endif

