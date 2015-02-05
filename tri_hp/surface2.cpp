//
//  surface22.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 12/30/14.
//
//

#include "surface2.h"

#include "bdry_ins.h"
#include <myblas.h>

//#define MPDEBUG

//#define DEBUG

//#define BODYFORCE


using namespace bdry_ins;

// extern FLT body[ND];


void surface2::init(input_map& inmap,void* gin) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;

	if (base.is_comm()) {
		keyword = base.idprefix + "_matching_block";
		if (!inmap.get(keyword,val)) {
			is_master = false;
		}
		else {
			is_master = true;
		}
	}
	else {
		is_master = true;
	}
	
	/* Let pressure be discontinuous */
	ostringstream vars;
	for(int n=0;n<x.NV-1;++n) {
		vars << n << ' ';
	}
	inmap[base.idprefix +"_c0_indices"] = vars.str();
	hp_deformable_bdry::init(inmap,gin);
	gbl = static_cast<global *>(gin);
	
	if (!is_master) return;
	
	keyword = base.idprefix + "_sigma";
	inmap.getwdefault(keyword,gbl->sigma,0.0);
	
	inmap.getwdefault(base.idprefix +"_p_ext",gbl->p_ext,0.0);
	
	gbl->mu2 = 0.0;
	gbl->rho2 = 0.0;
	if (base.is_comm()) {
		keyword = val +"_mu";
		if (!inmap.get(keyword,gbl->mu2)) {
			*x.gbl->log << "couldn't find matching blocks viscosity" << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		
		keyword = val +"_rho";
		if (!inmap.get(keyword,gbl->rho2)) {
			*x.gbl->log << "couldn't find matching blocks density" << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}
	
	return;
}

#ifndef petsc
void surface2::rsdl(int stage) {
	hp_deformable_bdry::rsdl(stage);
	
    if (is_master) {
        /* Communicate mass flux residual for better preconditioning */
        if (base.is_comm()) {
            int count = 0;
            for(int j=0;j<base.nseg+1;++j) {
                base.fsndbuf(count++) = gbl->vres(j,1)*gbl->rho2;
#ifdef MPDEBUG
                *x.gbl->log << gbl->vres(j,1)*gbl->rho2 << '\n';
#endif
            }
            for(int j=0;j<base.nseg;++j) {
                for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
                    base.fsndbuf(count++) = gbl->sres(j,m,1)*gbl->rho2;
                }
            }
            base.sndsize() = count;
            base.sndtype() = boundary::flt_msg;
            base.comm_prepare(boundary::all,0,boundary::master_slave);
            base.comm_exchange(boundary::all,0,boundary::master_slave);
            base.comm_wait(boundary::all,0,boundary::master_slave);
        }
    }
    else {
        int v0,sind;
        base.comm_prepare(boundary::all,0,boundary::master_slave);
        base.comm_exchange(boundary::all,0,boundary::master_slave);
        base.comm_wait(boundary::all,0,boundary::master_slave);
        int count = 0;
        int i = base.nseg-1;
        do {
            sind = base.seg(i);
            v0 = x.seg(sind).pnt(1);
#ifdef MPDEBUG
            *x.gbl->log << x.gbl->res.v(v0,x.NV-1) << ' ' << base.frcvbuf(0,count) << '\n';
#endif
            x.gbl->res.v(v0,x.NV-1) += base.frcvbuf(0,count++);
        } while(--i >= 0);
        v0 = x.seg(sind).pnt(0);
#ifdef MPDEBUG
        *x.gbl->log << x.gbl->res.v(v0,x.NV-1) << ' ' << base.frcvbuf(0,count) << '\n';
#endif
        x.gbl->res.v(v0,x.NV-1) += base.frcvbuf(0,count++);
        
        for(i=base.nseg-1;i>=0;--i) {
            sind = base.seg(i);
            int msgn = 1;
            for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
                x.gbl->res.s(sind,m,x.NV-1) += msgn*base.frcvbuf(0,count++);
                msgn *= -1;
            }
        }
    }
}
#endif

void surface2::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
    
	if (!is_master) return;
    
	int i,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
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
	
	for(n=0;n<tri_mesh::ND;++n)
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
		
		/* TANGENTIAL SPACING */
		res(0,i) = -ksprg(indx)*jcb;
		/* NORMAL FLUX */
		res(1,i) = -RAD(crd(0,i))*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1));
		/* UPWINDING BASED ON TANGENTIAL VELOCITY */
		res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*gbl->meshc(indx);
		
		/* surface2 TENSION SOURCE TERM X-DIRECTION */
		res(4,i) = +RAD(crd(0,i))*((x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i) +gbl->p_ext)*norm(0);
#ifdef AXISYMMETRIC
		res(4,i) += gbl->sigma*jcb;
#endif
		/* AND INTEGRATION BY PARTS TERM */
		res(5,i) = +RAD(crd(0,i))*gbl->sigma*norm(1)/jcb;
		
		
		/* surface2 TENSION SOURCE TERM Y-DIRECTION */
		res(6,i) = +RAD(crd(0,i))*((x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i) +gbl->p_ext)*norm(1);
		/* AND INTEGRATION BY PARTS TERM */
		res(7,i) = -RAD(crd(0,i))*gbl->sigma*norm(0)/jcb;
	}
	
	lf = 0.0;
	/* INTEGRATE & STORE surface2 TENSION SOURCE TERM */
	basis::tri(x.log2p)->intgrt1d(&lf(0)(0),&res(4,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(0)(0),&res(5,0));
	basis::tri(x.log2p)->intgrt1d(&lf(1)(0),&res(6,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(1)(0),&res(7,0));
	
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV)(0),&res(0,0));
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+1)(0),&res(1,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV+1)(0),&res(2,0));
	
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

	return;
}

/* Remember that vertex jacobians need to be written */
void surface2::setup_preconditioner() {
	
	if (!is_master) return;
	
	int indx,m,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> nrm;
	FLT h, hsm;
	FLT dttang, dtnorm;
	FLT vslp, strss;
	FLT drho, srho, smu;
	FLT nu1, nu2;
	FLT qmax, gam1, gam2;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd;
	TinyMatrix<FLT,4,MXGP> res;
	TinyMatrix<FLT,4,MXGP> lf;
	TinyVector<FLT,2> mvel;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	int last_phase, mp_phase;
	
	drho = x.gbl->rho -gbl->rho2;
	srho = x.gbl->rho +gbl->rho2;
	smu = x.gbl->mu +gbl->mu2;
	nu1 = x.gbl->mu/x.gbl->rho;
	if (gbl->rho2 > 0.0)
		nu2 = gbl->mu2/gbl->rho2;
	else
		nu2 = 0.0;
	
	/**************************************************/
	/* DETERMINE surface2 MOVEMENT TIME STEP              */
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
			// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5);/* FAILED IN NATES UPSTREAM surface2 WAVE CASE */
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
		
		qmax = mvel(0)*mvel(0)+mvel(1)*mvel(1);
		vslp = fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h);
		
		mvel(0) = x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
		mvel(1) = x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
#ifdef MESH_REF_VEL
		mvel(0) -= x.gbl->mesh_ref_vel(0);
		mvel(1) -= x.gbl->mesh_ref_vel(1);
#endif
		qmax = MAX(qmax,mvel(0)*mvel(0)+mvel(1)*mvel(1));
		vslp = MAX(vslp,fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h));
		
		hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
		
		dttang = 2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm;
#ifndef BODYFORCE
		strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*x.gbl->g*nrm(1)/h);
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
		dtnorm = 2.*vslp/hsm +x.gbl->bd(0) +1.*strss/(x.gbl->rho*sqrt(qmax +gam1) +gbl->rho2*sqrt(qmax +gam2));
		
		/* SET UP DISSIPATIVE COEFFICIENT */
		/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
		/* RESIDUAL HAS DX/2 WEIGHTING */
		/* |a| dx/2 dv/dx  dx/2 dpsi */
		/* |a| dx/2 2/dx dv/dpsi  dpsi */
		/* |a| dv/dpsi  dpsi */
		// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5); /* FAILED IN NATES UPSTREAM surface2 WAVE CASE */
		// gbl->meshc(indx) = gbl->adis/(h*(vslp/hsm +x.gbl->bd(0))); /* FAILED IN MOVING UP TESTS */
		gbl->meshc(indx) = gbl->adis/(h*(sqrt(qmax)/hsm +x.gbl->bd(0))); /* SEEMS THE BEST I'VE GOT */
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
		for(m=0;m<tri_mesh::ND;++m)
			for(n=0;n<tri_mesh::ND;++n)
				gbl->vdt(0,m,n) = 0.5*(gbl->vdt(0,m,n) +gbl->vdt(base.nseg+1,m,n));
		gbl->vdt(base.nseg+1) = gbl->vdt(0);
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
		
		/* TEMPORARY DIRECT FORMATION OF vdt^{-1} theta is angle of normal from horizontal */
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

void surface_outflow2::rsdl(int stage) {
	int bnum,sind;
	TinyVector<FLT,tri_mesh::ND> ubar, tangent, rp;
	FLT jcb;
	
	/* SET TANGENT RESDIUAL TO 0 */
	surface_fixed_pt2::rsdl(stage);
	
	switch(contact_type) {
		case(free_angle): {
			/* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
			/* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
			bnum = base.ebdry(surfbdry);
			if (surfbdry == 0) {
				sind = x.ebdry(bnum)->seg(x.ebdry(bnum)->nseg-1);
				x.crdtocht1d(sind);
				basis::tri(x.log2p)->ptprobe1d(2,&rp(0),&tangent(0),1.0,&x.cht(0,0),MXTM);
				jcb = sqrt(tangent(0)*tangent(0) +tangent(1)*tangent(1));
				x.gbl->res.v(base.pnt,0) += -RAD(rp(0))*surf->gbl->sigma*tangent(0)/jcb;
				x.gbl->res.v(base.pnt,1) += -RAD(rp(0))*surf->gbl->sigma*tangent(1)/jcb;
			}
			else {
				sind = x.ebdry(bnum)->seg(0);
				x.crdtocht1d(sind);
				basis::tri(x.log2p)->ptprobe1d(2,&rp(0),&tangent(0),-1.0,&x.cht(0,0),MXTM);
				jcb = sqrt(tangent(0)*tangent(0) +tangent(1)*tangent(1));
				x.gbl->res.v(base.pnt,0) -= -RAD(rp(0))*surf->gbl->sigma*tangent(0)/jcb;
				x.gbl->res.v(base.pnt,1) -= -RAD(rp(0))*surf->gbl->sigma*tangent(1)/jcb;
			}
			break;
		}
		case(fixed_angle): {
			/* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
			/* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
			if (wall_type == curved) {
				if (surfbdry == 0) {
					x.ebdry(base.ebdry(1))->bdry_normal(0,-1.0,wall_normal);
				}
				else {
					x.ebdry(base.ebdry(0))->bdry_normal(x.ebdry(base.ebdry(0))->nseg-1,1.0,wall_normal);
				}
			}
			TinyVector<FLT,tri_mesh::ND> tangent;
			if (surfbdry == 0) {
				/* Surf-boundary then point then wall (in ccw sense) */
				tangent(0) = wall_normal(0)*sin(contact_angle) +wall_normal(1)*cos(contact_angle);
				tangent(1) = -wall_normal(0)*cos(contact_angle) +wall_normal(1)*sin(contact_angle);
				x.gbl->res.v(base.pnt,0) -= RAD(x.pnts(base.pnt)(0))*surf->gbl->sigma*tangent(0);
				x.gbl->res.v(base.pnt,1) -= RAD(x.pnts(base.pnt)(0))*surf->gbl->sigma*tangent(1);
			}
			else {
				/* Wall then point then Surf-boundary (in ccw sense) */
				tangent(0) = wall_normal(0)*sin(contact_angle) -wall_normal(1)*cos(contact_angle);
				tangent(1) = wall_normal(0)*cos(contact_angle) +wall_normal(1)*sin(contact_angle);
				x.gbl->res.v(base.pnt,0) -= RAD(x.pnts(base.pnt)(0))*surf->gbl->sigma*tangent(0);
				x.gbl->res.v(base.pnt,1) -= RAD(x.pnts(base.pnt)(0))*surf->gbl->sigma*tangent(1);
			}
		}
		case(prdc): {
			break;
		}
	}
	
	return;
}

#ifdef petsc
//void surface_outflow::petsc_jacobian_dirichlet() {
//				/* Replace tangential residual equation with equation that fixes wall position */
//				int row = (x.NV+tri_mesh::ND)*base.pnt +x.NV;
//				MatZeroRows(x.petsc_J,1,&row,PETSC_NULL);
//				TinyVector<int,2> col(row,row+1);
//				MatSetValuesLocal(x.petsc_J,1,&row,2,col.data(),wall_normal.data(),INSERT_VALUES);
//				MatAssemblyBegin(x.petsc_J,MAT_FINAL_ASSEMBLY);
//				MatAssemblyEnd(x.petsc_J,MAT_FINAL_ASSEMBLY);
//			}


void surface_outflow2::petsc_jacobian() {
	
	if (contact_type == free_angle)
		*x.gbl->log << "Haven't added jacobian of surfacetension free contact add\n";
	
	int indx;
	if (surfbdry == 0) {
		indx = x.ebdry(base.ebdry(0))->nseg;
	}
	else {
		indx = 0;
	}
	
#ifdef MY_SPARSE
	
	/* GET X & Y MESH MOVEMENT ROW */
	/* CONSTRAIN MOTION NORMAL TO BOUNDARY */
	int row = (x.NV+tri_mesh::ND)*base.pnt +x.NV;
	int nnz1 = x.J._cpt(row+1) -x.J._cpt(row);
	int nnz2 = x.J._cpt(row+2) -x.J._cpt(row+1);
	Array<int,1> cols(nnz1);
	Array<FLT,2> vals(2,nnz1);
	
	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	if (nnz1 != nnz2) {
		*x.gbl->log << "zeros problem in deforming mesh on angled boundary\n";
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	int row1 = x.J._cpt(row);
	int row2 = x.J._cpt(row+1);
	for(int col=0;col<nnz1;++col) {
		if (x.J._col(row1++) != x.J._col(row2++)) {
			*x.gbl->log << "zeros indexing problem in deforming mesh on angled boundary\n";
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}
	
	
	FLT J = surf->gbl->vdt(indx,0,0)*surf->gbl->vdt(indx,1,1) -surf->gbl->vdt(indx,1,0)*surf->gbl->vdt(indx,0,1);
	
	vals(0,Range(0,nnz1-1)) = -x.J._val(Range(x.J._cpt(row),x.J._cpt(row+1)-1))*surf->gbl->vdt(indx,1,0)/J;
	vals(0,Range(0,nnz1-1)) += x.J._val(Range(x.J._cpt(row+1),x.J._cpt(row+2)-1))*surf->gbl->vdt(indx,0,0)/J;
	
	/* Replace x equation with tangential position equation */
	/* Replacy y equation with normal displacement equation */
	/* Normal Equation */
	vals(1,Range::all()) = 0.0;
	for(int col=0;col<nnz1;++col) {
		if (x.J._col(row1+col) == row) {
			vals(1,col) = wall_normal(0);
			break;
		}
	}
	for(int col=0;col<nnz1;++col) {
		if (x.J._col(row1+col) == row+1) {
			vals(1,col) = wall_normal(1);
			break;
		}
	}
	
	/* tangent = -sin(theta) i +cos(theta) j */
	/* normal = cos(theta) i + sin(theta) j */
	/* Rotate equations for diagonal dominance to match what is done to residual */
	Array<FLT,2> temp(2,nnz1);
	temp(0,Range::all()) =  vals(0,Range::all())*surf->gbl->vdt(indx,0,1) +vals(1,Range::all())*wall_normal(0);
	temp(1,Range::all()) =  vals(0,Range::all())*surf->gbl->vdt(indx,1,1) +vals(1,Range::all())*wall_normal(1);
	
	TinyVector<int,2> rows(row,row+1);
	x.J._val(Range(x.J._cpt(row),x.J._cpt(row+1)-1)) = temp(0,Range::all());
	x.J._val(Range(x.J._cpt(row+1),x.J._cpt(row+2)-1)) = temp(1,Range::all());
#else
	MatAssemblyBegin(x.petsc_J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(x.petsc_J,MAT_FINAL_ASSEMBLY);
	
	int nnz1, nnz2;
	const PetscInt *cols1, *cols2;
	const PetscScalar *vals1, *vals2;
	int row = (x.NV+tri_mesh::ND)*base.pnt +x.NV;
	
	MatGetRow(x.petsc_J,row,&nnz1,&cols1,&vals1);
	MatGetRow(x.petsc_J,row+1,&nnz2,&cols2,&vals2);
	if (nnz1 != nnz2) {
		*x.gbl->log << "zeros problem in surface_outflow\n";
	}
	
	
	FLT J = surf->gbl->vdt(indx)(0,0)*surf->gbl->vdt(indx,1,1) -surf->gbl->vdt(indx,1,0)*surf->gbl->vdt(indx,0,1);
	Array<int,1> cols(nnz1);
	Array<FLT,2> vals(2,nnz1);
	for (int col=0;col<nnz1;++col) {
		cols(col) = cols1[col];
		vals(0,col) = -vals1[col]*surf->gbl->vdt(indx,1,0)/J;
	}
	for (int col=0;col<nnz1;++col) {
		if (cols(col) != cols2[col]) {
			*x.gbl->log << "zeros problem in deforming mesh on angled boundary\n";
			*x.gbl->log << cols << ' ' << col << ' ' << cols2[col] << std::endl;
		}
		vals(0,col) += vals2[col]*surf->gbl->vdt(indx,0,0)/J;
	}
	MatRestoreRow(x.petsc_J,row,&nnz1,&cols1,&vals1);
	MatRestoreRow(x.petsc_J,row+1,&nnz2,&cols2,&vals2);
	
	/* Replace x equation with tangential position equation */
	/* Replacy y equation with normal displacement equation */
	/* Normal Equation */
	vals(1,Range::all()) = 0.0;
	for(int col=0;col<nnz1;++col) {
		if (cols(col) == row) {
			vals(1,col) = wall_normal(0);
			break;
		}
	}
	for(int col=0;col<nnz1;++col) {
		if (cols(col) == row+1) {
			vals(1,col) = wall_normal(1);
			break;
		}
	}
	
	/* tangent = -sin(theta) i +cos(theta) j */
	/* normal = cos(theta) i + sin(theta) j */
	/* Rotate equations for diagonal dominance to match what is done to residual */
	Array<FLT,2> temp(2,nnz1);
	temp(0,Range::all()) =  vals(0,Range::all())*surf->gbl->vdt(indx,0,1) +vals(1,Range::all())*wall_normal(0);
	temp(1,Range::all()) =  vals(0,Range::all())*surf->gbl->vdt(indx,1,1) +vals(1,Range::all())*wall_normal(1);
	
	TinyVector<int,2> rows(row,row+1);
	MatSetValuesLocal(x.petsc_J,2,rows.data(),nnz1,cols.data(),temp.data(),INSERT_VALUES);
#endif
	
}
#endif
