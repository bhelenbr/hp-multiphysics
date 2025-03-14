#include "melt_buoyancy.h"
#include <myblas.h>

#define NEW_STABILIZATION

using namespace bdry_buoyancy;

void surface_marangoni::init(input_map& inmap) {
	bdry_ins::surface::init(inmap);
	
	if (!is_master) return;
	
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
void surface_marangoni::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
	
	if (!is_master) return;

	int i,n,sind;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV),flx(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,9,MXGP> res;
	Array<FLT,1> axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND),au(x.NV);
	
	sind = base.seg(indx);
	x.crdtocht1d(sind);
    
#ifdef NEW_STABILIZATION
    /* Calculate stabilization constant based on analysis of linear elements and constant tau */
    int v0 = x.seg(sind).pnt(0);
    int v1 = x.seg(sind).pnt(1);
    norm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
    norm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
    FLT h = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
    FLT hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));

    mvel(0) = x.uht(0)(0) -(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
    mvel(1) = x.uht(1)(0)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
#ifdef MESH_REF_VEL
    mvel(0) -= x.hp_gbl->mesh_ref_vel(0);
    mvel(1) -= x.hp_gbl->mesh_ref_vel(1);
#endif
    FLT vslp0 = (-mvel(0)*norm(1) +mvel(1)*norm(0))/h;
    
    mvel(0) = x.uht(0)(1)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
    mvel(1) = x.uht(1)(1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
#ifdef MESH_REF_VEL
    mvel(0) -= x.hp_gbl->mesh_ref_vel(0);
    mvel(1) -= x.hp_gbl->mesh_ref_vel(1);
#endif
    FLT vslp1 = (-mvel(0)*norm(1) +mvel(1)*norm(0))/h;
    hp_bdry_gbl->meshc(indx) = hp_bdry_gbl->adis*hsm*(3*(abs(vslp0)+abs(vslp1)) +vslp0-vslp1)/(4*(vslp0*vslp0+vslp0*vslp1+vslp1*vslp1)+10*FLT_EPSILON)*2/h;
#endif
    
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
			mvel(n,i) -= x.hp_gbl->mesh_ref_vel(n);
#endif
		}
		
		/* Evaluate Fluxes */
		axpt(0) = crd(0,i); axpt(1) = crd(1,i);
		amv(0) = -(mvel(0,i)-u(0)(i)); amv(1) = -(mvel(1,i)-u(1)(i));
		anorm(0)= norm(0)/jcb; anorm(1) = norm(1)/jcb;
		for(int n=0;n<x.NV;++n)
			au(n) = u(n)(i);
		for(int n=0;n<x.NV;++n) {
			flx(n) = fluxes[n]->Eval(au,axpt,amv,anorm,x.gbl->time);
		}
		
		/* Evaluate Surface Tension */
		FLT sigma = sigma_vs_T.Eval(au(2));
		/* TANGENTIAL SPACING */
		res(0,i) = -ksprg(indx)*jcb;
		/* NORMAL FLUX */
		res(1,i) = -RAD(crd(0,i))*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1)) +flx(x.NV-1);
		/* UPWINDING BASED ON TANGENTIAL VELOCITY */
		res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*hp_bdry_gbl->meshc(indx);
		
		/* SURFACE TENSION SOURCE TERM X-DIRECTION */
		res(4,i) = +RAD(crd(0,i))*((x.hp_ins_gbl->rho -surface_gbl->rho2)*x.gbl->g*crd(1,i) +surface_gbl->p_ext)*norm(0) +flx(0)*jcb;
#ifdef AXISYMMETRIC
		res(4,i) += sigma*jcb;
#endif
		/* AND INTEGRATION BY PARTS TERM */
		res(5,i) = +RAD(crd(0,i))*sigma*norm(1)/jcb;
		
		
		/* SURFACE TENSION SOURCE TERM Y-DIRECTION */
		res(6,i) = +RAD(crd(0,i))*((x.hp_ins_gbl->rho -surface_gbl->rho2)*x.gbl->g*crd(1,i) +surface_gbl->p_ext)*norm(1) +flx(1)*jcb;
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
	
	return;
}

void triple_junction::init(input_map& inmap) {
	melt_facet_pt::init(inmap);
	inmap.getwdefault(base.idprefix +"_growth_angle",growth_angle,0.0);
	growth_angle *= M_PI/180.0;
	
	if (!dynamic_cast<melt_buoyancy *>(surf)) {
		*x.gbl->log << "Did not make triple junction point general with respect to orientation of surfaces\n";
		sim::abort(__LINE__, __FILE__, x.gbl->log);
	}
}


void triple_junction::element_rsdl(Array<FLT,1> lf) {
	
	hp_vrtx_bdry::element_rsdl(lf);
	
	if (!surf->is_master) return;

	/* Solidification kinetics at nucleation point */
	/* Calculate temperature from velocity of point & 2D nucleation coefficient only */
	/* assuming growth at facet angle */
	
	Array<FLT,1> u(x.NV);
	TinyVector<FLT, 2> xp, dxpdpsi, mvel, anorm;
	

	int seg,sind;
	int seg_surf,sind_surf;
	const int bnum = base.ebdry(surfbdry);
	const int bnum_surf = base.ebdry(1-surfbdry);
	if (surfbdry == 0) {
		seg = x.ebdry(bnum)->nseg-1;
		sind = x.ebdry(bnum)->seg(seg);
		seg_surf = 0;
		sind_surf = x.ebdry(bnum_surf)->seg(seg_surf);
	}
	else {
		seg = 0;
		sind = x.ebdry(bnum)->seg(seg);
		seg_surf = x.ebdry(bnum_surf)->nseg-1;
		sind_surf = x.ebdry(bnum_surf)->seg(seg_surf);
	}
	
	x.crdtocht1d(sind);
	x.ugtouht1d(sind);
	
	int v0 = base.pnt;
	basis::tri(x.log2p)->ptprobe1d(2,xp.data(),dxpdpsi.data(),1.0-2.*surfbdry,&x.cht(0,0),MXTM);
	basis::tri(x.log2p)->ptprobe1d(4,u.data(),1.0-2.*surfbdry,&x.uht(0)(0),MXTM);
	/* RELATIVE VELOCITY STORED IN MVEL(N)*/
	for(int n=0;n<tri_mesh::ND;++n) {
		mvel(n) = u(n) -x.gbl->bd(0)*(xp(n) -x.vrtxbd(1)(v0)(n));
#ifdef MESH_REF_VEL
		mvel(n) -= x.hp_gbl->mesh_ref_vel(n);
#endif
	}
	
	FLT jcb = sqrt(dxpdpsi(0)*dxpdpsi(0)+dxpdpsi(1)*dxpdpsi(1));
	FLT DT = surf->ibc->f(2, xp, x.gbl->time) -u(2);
	melt_buoyancy *surf1 = dynamic_cast<melt_buoyancy *>(surf);
#ifdef TWOFACETS
	/* This is to allow the general expression at the triple point */
     anorm(0)= dxpdpsi(1)/jcb; anorm(1) = -dxpdpsi(0)/jcb;
     FLT sint = -surf1->bdry_cd_gbl->facetdir(0)*anorm(1) +surf1->bdry_cd_gbl->facetdir(1)*anorm(0);
     FLT K = surf1->calculate_kinetic_coefficients(DT,sint);
#else
	FLT K = surf1->calculate_kinetic_coefficients(DT,0.0);
#endif
	FLT res1 = jcb*RAD(xp(0))*surf1->melt_cd_gbl->rho_s*(mvel(0)*surf1->melt_cd_gbl->facetdir(0) +mvel(1)*surf1->melt_cd_gbl->facetdir(1));
	/* Kinetic equation for surface temperature */
	lf(x.NV+1) = RAD(xp(0))*surf1->melt_cd_gbl->rho_s*(-DT)*jcb +res1*K;  // -gbl->K_gt*kappa?;
	
	/* Motion constraint */
	/* Constraint is that motion of triple junction relative to solid is in direction of free surface + growth angle */
	/* Find tangent to free-surface */
	x.crdtocht1d(sind_surf);
	basis::tri(x.log2p)->ptprobe1d(2,xp.data(),dxpdpsi.data(),1.0-2.*surfbdry,&x.cht(0,0),MXTM);
	jcb = sqrt(dxpdpsi(0)*dxpdpsi(0)+dxpdpsi(1)*dxpdpsi(1));
	
	/* Rotate vector normal to side counterclockwise by growth angle */
	//Reminder
	//anorm(0) =  dxpdpsi(1)/jcb;
	//anorm(1) = -dxpdpsi(0)/jcb;
	// Rotation matrix is [cos(t) -sin(t); sin(t) cos(t)]
	anorm(0) = (cos(growth_angle)*dxpdpsi(1) +sin(growth_angle)*dxpdpsi(0))/jcb;
	anorm(1) = (sin(growth_angle)*dxpdpsi(1) -cos(growth_angle)*dxpdpsi(0))/jcb;
	
	// orthogonality constraint
	lf(x.NV) = mvel(0)*anorm(0) +mvel(1)*anorm(1);
}


#ifdef petsc
void triple_junction::petsc_jacobian() {
	
	if (!surf->is_master) return;
	
	const int sm = basis::tri(x.log2p)->sm();
	const int vdofs = x.NV+x.ND;
	int nvars = vdofs*(2*(sm+1)+1);  // two connecting edges affect residual
	Array<FLT,2> K(vdofs,nvars);
	Array<FLT,1> lf(vdofs), Rbar(vdofs);
	Array<int,1> rows(vdofs), cols(nvars);
	
	/* Jacobian for row determine x position */
	/* Variables effecting dT/dx are strictly those on the edge */
	element_rsdl(lf);
	Rbar = lf;
	
	int endpt,seg,sind;
	int endpt_surf,seg_surf,sind_surf;
	const int bnum = base.ebdry(surfbdry);
	const int bnum_surf = base.ebdry(1-surfbdry);
	if (surfbdry == 0) {
		seg = x.ebdry(bnum)->nseg-1;
		sind = x.ebdry(bnum)->seg(seg);
		endpt = x.ebdry(bnum)->nseg;
		seg_surf = 0;
		sind_surf = x.ebdry(bnum_surf)->seg(seg_surf);
		endpt_surf = 0;
	}
	else {
		seg = 0;
		sind = x.ebdry(bnum)->seg(seg);
		endpt = 0;
		seg_surf = x.ebdry(bnum_surf)->nseg-1;
		sind_surf = x.ebdry(bnum_surf)->seg(seg_surf);
		endpt_surf = x.ebdry(bnum_surf)->nseg;
	}
	
	Array<FLT,1> dw(x.NV);
	dw = 0.0;
	for(int i=0;i<2;++i)
		for(int n=0;n<x.NV;++n)
			dw(n) = dw(n) + fabs(x.uht(n)(i));
	
	dw *= eps_r;
	dw += eps_a;
	FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;
	
	/* Numerically create Jacobian */
	/* These are all the variables on the solidification side */
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
	
	
	/* These are variables on the free-surface side */
	/* Already did the triple junction point so can skip that */
	for(int mode = 1-surfbdry; mode < 2-surfbdry; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			x.ug.v(x.seg(sind_surf).pnt(mode),var) += dw(var);
			element_rsdl(lf);
			
			int krow = 0;
			for(int n=0;n<vdofs;++n)
				K(krow++,kcol) = (lf(n)-Rbar(n))/dw(var);
			++kcol;
			
			x.ug.v(x.seg(sind_surf).pnt(mode),var) -= dw(var);
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var){
			x.pnts(x.seg(sind_surf).pnt(mode))(var) += dx;
			
			element_rsdl(lf);
			
			int krow = 0;
			for(int n=0;n<vdofs;++n)
				K(krow++,kcol) = (lf(n)-Rbar(n))/dx;
			++kcol;
			
			x.pnts(x.seg(sind_surf).pnt(mode))(var) -= dx;
		}
	}
	
	for(int mode = 0; mode < sm; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			x.ug.s(sind_surf,mode,var) += dw(var);
			element_rsdl(lf);
			
			int krow = 0;
			for(int n=0;n<vdofs;++n)
				K(krow++,kcol) = (lf(n)-Rbar(n))/dw(var);
			++kcol;
			
			x.ug.s(sind_surf,mode,var) -= dw(var);
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var){
			x.hp_ebdry(bnum_surf)->crv(seg_surf,mode)(var) += dx;
			element_rsdl(lf);
			
			int krow = 0;
			for(int n=0;n<vdofs;++n)
				K(krow++,kcol) = (lf(n)-Rbar(n))/dx;
			++kcol;
			
			x.hp_ebdry(bnum_surf)->crv(seg_surf,mode)(var) -= dx;
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
	int gindxND = surf->jacobian_start +seg*tri_mesh::ND*sm;
	for(int mode = 0; mode < sm; ++mode) {
		for(int var = 0; var < x.NV; ++var)
			cols(ind++) = gindxNV++;
		
		for(int var = 0; var < tri_mesh::ND; ++var)
			cols(ind++) = gindxND++;
	}
	
	/* Free-surface dependency */
	for(int mode = 1-surfbdry; mode < 2-surfbdry; ++mode) {
		int gindx = vdofs*x.seg(sind_surf).pnt(mode);
		for(int var = 0; var < vdofs; ++var) {
			cols(ind++) = gindx++;
		}
	}
	
	gindxNV = x.npnt*vdofs +x.NV*sind_surf*sm;
	gindxND = x.hp_ebdry(bnum_surf)->jacobian_start +seg_surf*tri_mesh::ND*sm;
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
	This does not work
	MatSetValuesLocal(x.petsc_J,vdofs,rows.data(),nvars,cols.data(),K.data(),ADD_VALUES);
#endif
}
#endif



