/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_komega.h"
#include "../hp_boundary.h"

#define CALC_TAU1
//#define CALC_TAU2


#define BODYFORCE

//#define WILCOX1988
//#define WILCOX1988_KL
#define WILCOX2006

double psifunc(double x,double xs,double xe) {
    if (x >= xe) {
        return(1.0);
    }
    else if (x < xs) {
        return(0.0);
    }
    double theta = M_PI/2*(2*x - (xs + xe))/(xe - xs);
    return(0.5*(sin(theta) +1));
}

double dpsifunc(double x,double xs,double xe) {
    if (x >= xe || x < xs) {
        return(0.0);
    }
    double theta = M_PI/2*(2*x - (xs + xe))/(xe - xs);
    return(0.5*M_PI/xe*cos(theta));
}

void tri_hp_komega::element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im) {
	FLT fluxx,fluxy;
	const int NV = 5;
	TinyVector<int,3> v;
    TinyVector<FLT,ND> pt;
	TinyMatrix<FLT,ND,ND> ldcrd;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
	int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
    FLT lmu = gbl->mu, rhorbd0, lrhobd0, cjcb, cjcbi;
    FLT psiktld, psinktld, ktrb, omg, omgMx, tmu, tmuLmtd, mutld, cjcbik, cjcbiomg, strninv;
    TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-1,NV-1> visc;
    TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV-1,NV-1> cv, df;
    TinyVector<FLT,NV> tres;
    TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> mvel; // for local mesh velocity info
    
    /* Turbulence Model Constants*/
    const FLT kinf = gbl->kinf;
    const FLT omginf = gbl->omginf;
    const FLT epslnk = gbl->epslnk;
    const FLT sgmomg = 0.5,  betastr = 0.09;
#ifdef WILCOX1988
    const FLT sgmk = 0.5, betakomg = 0.075, gamma = 5./9.;
#endif
#ifdef WILCOX1988_KL
    FLT vrtctinv;
    const FLT sgmk = 0.5, betakomg = 0.075, gamma = 5./9.;
#endif
#ifdef WILCOX2006
    FLT dktlddx, dktlddy, dkdx, dkdy, CrssD, omgLmtr;
    FLT betakomg = 0.0708;
    const FLT sgmk = 0.6, gamma = 13./25., sgmdo = 1./8., Clim = 7.0/8.0; // The default value from Wilcox(2006) for Clim is 7./8.
#endif
    const FLT k_mom= 1.0; // the term in the momentum equation including k
    const FLT susk = 1.0, susomg = 1.0; // turbulence sustaning terms
    
	/* LOAD INDICES OF VERTEX POINTS */
	v = tri(tind).pnt;
    
	/* IF TINFO > -1 IT IS CURVED ELEMENT */
    if (tri(tind).info > -1) {
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		crdtocht(tind);

		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		/* dcrd(i,j) is derivative of physical coordinate i with respect to curvilinear coordinate j */
		/* dxi/dx = dy/deta, dxi/dy = -dx/deta, deta/dx = -dy/dxi, deta/dy = dx/dxi (divided by jacobian) */
		for(int n=0;n<ND;++n)
			basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
	}
	else {
		/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(int n=0;n<ND;++n)
			basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);

		/* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
		for(int n=0;n<ND;++n) {
			ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
			ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
		}
	}

	/* CALCULATE MESH VELOCITY */
	for(int i=0;i<lgpx;++i) {
		for(int j=0;j<lgpn;++j) {
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
		for(int n=0;n<NV-1;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
		basis::tri(log2p)->proj(&uht(NV-1)(0),&u(NV-1)(0,0),MXGP);
	}
	else {
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);
	}

	/* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
	for(int n=0;n<NV;++n){
		for(int i=0;i<basis::tri(log2p)->tm();++i){
			lf_re(n)(i) = 0.0;
			lf_im(n)(i) = 0.0;	
		}
	}

#ifdef CALC_TAU1
    FLT qmax = 0.0;
    FLT qmax2 = 0.0;
    FLT hmax = 0.0;
    FLT tmuav = 0.0, mutldav = 0.0;
    FLT jcb = 0.25*area(tind);
    FLT jcbmin = jcb;
    FLT h;
#endif

    if (tri(tind).info > -1) {
		/* CURVED ELEMENT */
		/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {

				fluxx = gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
				fluxy = gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));

				/* CONTINUITY EQUATION FLUXES */
				du(NV-1,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
				du(NV-1,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;

				/* CONVECTIVE FLUXES */
				for(int n=0;n<NV-1;++n) {
					cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
					cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
				}

				/* PRESSURE TERMS */
				/* U-MOMENTUM */
				cv(0,0)(i,j) += dcrd(1,1)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				cv(0,1)(i,j) -= dcrd(1,0)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				/* V-MOMENTUM */
				cv(1,0)(i,j) -= dcrd(0,1)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				cv(1,1)(i,j) += dcrd(0,0)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
			}
		}
		for(int n=0;n<NV-1;++n)
            basis::tri(log2p)->intgrtrs(&lf_im(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
        basis::tri(log2p)->intgrtrs(&lf_im(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);

		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			/* TIME DERIVATIVE TERMS */ 
			for(int i=0;i<lgpx;++i) {
				for(int j=0;j<lgpn;++j) {
					cjcb = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                    rhorbd0 = gbl->rho*gbl->bd(0)*RAD(crd(0)(i,j))*cjcb;
                    psiktld = psifunc(u(2)(i,j),0.0,epslnk);
                    psinktld = psifunc(-u(2)(i,j),0.0,epslnk);
                    ktrb = psiktld*u(2)(i,j);
                    omg = exp(u(3)(i,j));
                    tmu = gbl->rho*ktrb/omg;
                    FLT dudx = (dcrd(1,1)(i,j)*du(0,0)(i,j) -dcrd(1,0)(i,j)*du(0,1)(i,j))/cjcb;
                    FLT dvdx = (dcrd(1,1)(i,j)*du(1,0)(i,j) -dcrd(1,0)(i,j)*du(1,1)(i,j))/cjcb;
                    FLT dudy = (dcrd(0,0)(i,j)*du(0,1)(i,j) -dcrd(0,1)(i,j)*du(0,0)(i,j))/cjcb;
                    FLT dvdy = (dcrd(0,0)(i,j)*du(1,1)(i,j) -dcrd(0,1)(i,j)*du(1,0)(i,j))/cjcb;
                    strninv = dudx*dudx +0.5*(dudy +dvdx)*(dudy +dvdx) +dvdy*dvdy;
#ifdef AXISYMMETRIC
                    strninv += u(0)(i,j)*u(0)(i,j)/crd(0)(i,j)/crd(0)(i,j);
#endif
#ifndef WILCOX2006
                    omgMx = omg;
                    tmuLmtd = tmu;
#else
                    /* STRESS-LIMITER MODIFICATION */
                    omgLmtr = Clim*sqrt(2.0*strninv/betastr);
                    omgMx = std::max(omg,omgLmtr);
                    tmuLmtd = gbl->rho*ktrb/omgMx;
#endif
                    mutld = tmu -gbl->rho*psinktld*u(2)(i,j)/omginf;
                    cjcbi = (lmu +tmuLmtd)*RAD(crd(0)(i,j))/cjcb;
                    cjcbik = (lmu +sgmk*mutld)*RAD(crd(0)(i,j))/cjcb;
                    cjcbiomg = (lmu +sgmomg*tmu)*RAD(crd(0)(i,j))/cjcb;

					/* UNSTEADY TERMS */
					for(int n=0;n<NV-1;++n)
						res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p)(tind,n,i,j);
					res(NV-1)(i,j) = rhorbd0 +dugdt(log2p)(tind,NV-1,i,j);

#ifdef AXISYMMETRIC
					res(0)(i,j) -= cjcb*(u(NV-1)(i,j) -2.*(lmu +tmuLmtd)*u(0)(i,j)/crd(0)(i,j) +k_mom*2./3.*gbl->rho*ktrb);
#endif
#ifdef BODYFORCE
					res(0)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(0);
					res(1)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(1);
#endif
                    
#ifdef MMS
                    /* source terms for MMS */
                    pt(0) = crd(0)(i,j);
                    pt(1) = crd(1)(i,j);
                    for(int n = 0; n < NV; ++n)
                        res(n)(i,j) -= cjcb*gbl->src->f(n,pt,gbl->time);
#endif
                            
					/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
					/* INDICES ARE 1: EQUATION U, V, K-TILDE OR ln(OMEGA) 2: VARIABLE (U, V, K-TILDE OR ln(OMEGA)), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
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
                    
                    /* K-TILDE DIFFUSION */
                    visc(2,2)(0,0) = -cjcbik*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
                    visc(2,2)(1,1) = -cjcbik*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
                    visc(2,2)(0,1) =  cjcbik*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI2II2II1II0I visc(2,2)(0,1)
                    
                    /* ln(OMEGA) DIFFUSION */
                    visc(3,3)(0,0) = cjcbiomg/cjcbik*visc(2,2)(0,0);
                    visc(3,3)(1,1) = cjcbiomg/cjcbik*visc(2,2)(1,1);
                    visc(3,3)(0,1) = cjcbiomg/cjcbik*visc(2,2)(0,1);
#define                 viscI3II3II1II0I visc(3,3)(0,1)

					df(0,0)(i,j) = +visc(0,0)(0,0)*du(0,0)(i,j) +visc(0,1)(0,0)*du(1,0)(i,j)
									+visc(0,0)(0,1)*du(0,1)(i,j) +visc(0,1)(0,1)*du(1,1)(i,j) +k_mom*RAD(crd(0)(i,j))*2./3.*gbl->rho*dcrd(1,1)(i,j)*ktrb;

					df(0,1)(i,j) = +viscI0II0II1II0I*du(0,0)(i,j) +visc(0,1)(1,0)*du(1,0)(i,j)
									+visc(0,0)(1,1)*du(0,1)(i,j) +visc(0,1)(1,1)*du(1,1)(i,j) -k_mom*RAD(crd(0)(i,j))*2./3.*gbl->rho*dcrd(1,0)(i,j)*ktrb;

					df(1,0)(i,j) = +viscI1II0II0II0I*du(0,0)(i,j) +visc(1,1)(0,0)*du(1,0)(i,j)
									+viscI1II0II0II1I*du(0,1)(i,j) +visc(1,1)(0,1)*du(1,1)(i,j) -k_mom*RAD(crd(0)(i,j))*2./3.*gbl->rho*dcrd(0,1)(i,j)*ktrb;

					df(1,1)(i,j) = +viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
									+viscI1II0II1II1I*du(0,1)(i,j) +visc(1,1)(1,1)*du(1,1)(i,j) +k_mom*RAD(crd(0)(i,j))*2./3.*gbl->rho*dcrd(0,0)(i,j)*ktrb;
                    
                    df(2,0)(i,j) = visc(2,2)(0,0)*du(2,0)(i,j) +visc(2,2)(0,1)*du(2,1)(i,j);
                    df(2,1)(i,j) = viscI2II2II1II0I*du(2,0)(i,j) +visc(2,2)(1,1)*du(2,1)(i,j);
                    df(3,0)(i,j) = visc(3,3)(0,0)*du(3,0)(i,j) +visc(3,3)(0,1)*du(3,1)(i,j);
                    df(3,1)(i,j) = viscI3II3II1II0I*du(3,0)(i,j) +visc(3,3)(1,1)*du(3,1)(i,j);
                    
					for(int n=0;n<NV-1;++n) {
						cv(n,0)(i,j) += df(n,0)(i,j);
						cv(n,1)(i,j) += df(n,1)(i,j);
                    }

                    /*PRODUCTION TERMS FOR K-TILDE AND ln(omega) */
#ifndef WILCOX1988_KL
                    res(2)(i,j) -= RAD(crd(0)(i,j))*tmuLmtd*2.0*strninv*cjcb;
                    res(3)(i,j) -= RAD(crd(0)(i,j))*gamma*gbl->rho*2.0*strninv/omgMx*cjcb;
#else
                    vrtctinv = 0.5*(dudy -dvdx)*(dudy -dvdx);
                    res(2)(i,j) -= RAD(crd(0)(i,j))*2.0*tmu*sqrt(vrtctinv*strninv)*cjcb;
                    res(3)(i,j) -= RAD(crd(0)(i,j))*gamma*gbl->rho*2.0*sqrt(vrtctinv*strninv)/omg*cjcb;
#endif
                    /* PRODUCTION TERM FOR ln(OMEGA) (DUE TO LOGARITHMIC TRANSFORMATION) */
                    FLT domgtlddx = (dcrd(1,1)(i,j)*du(3,0)(i,j) -dcrd(1,0)(i,j)*du(3,1)(i,j))/cjcb;
                    FLT domgtlddy = (dcrd(0,0)(i,j)*du(3,1)(i,j) -dcrd(0,1)(i,j)*du(3,0)(i,j))/cjcb;
                    res(3)(i,j) -= cjcbiomg*(domgtlddx*domgtlddx + domgtlddy*domgtlddy)*cjcb*cjcb;
                    
                    /* ROUND JET CORRECTION FOR WILCOX2006 */
#ifdef WILCOX2006
#ifdef AXISYMMETRIC
                    FLT chi = fabs(0.25*(dudy -dvdx)*(dudy -dvdx)*u(0)(i,j)/crd(0)(i,j)/pow(betastr*omg,3.0));
                    FLT f_beta = (1.0 +85.0*chi)/(1.0 +100.0*chi);
                    betakomg *= f_beta;
#endif
#endif
                    /* DISSIPATION TERMS FOR K-TILDE and ln(OMEGA)  */
                    res(2)(i,j) += RAD(crd(0)(i,j))*betastr*gbl->rho*(omg*ktrb + omginf*psinktld*u(2)(i,j))*cjcb;
                    res(3)(i,j) += RAD(crd(0)(i,j))*betakomg*gbl->rho*omg*cjcb;
#ifdef WILCOX2006
                    /* CROSS-DIFFUSION TERM FOR ln(OMEGA) */
                    dktlddx = (dcrd(1,1)(i,j)*du(2,0)(i,j) -dcrd(1,0)(i,j)*du(2,1)(i,j))/cjcb;
                    dktlddy = (dcrd(0,0)(i,j)*du(2,1)(i,j) -dcrd(0,1)(i,j)*du(2,0)(i,j))/cjcb;
                    dkdx = dpsifunc(u(2)(i,j),0,epslnk)*dktlddx*u(2)(i,j) +psifunc(u(2)(i,j),0,epslnk)*dktlddx;
                    dkdy = dpsifunc(u(2)(i,j),0,epslnk)*dktlddy*u(2)(i,j) +psifunc(u(2)(i,j),0,epslnk)*dktlddy;
                    CrssD = dkdx*domgtlddx +dkdy*domgtlddy;
                    if ( CrssD > 0.0)
                        res(3)(i,j) -= RAD(crd(0)(i,j))*gbl->rho*CrssD/omg*sgmdo*cjcb;
#endif
                    /* TURBULENCE SUSTAINING TERMS */
                    res(2)(i,j) -= RAD(crd(0)(i,j))*susk*betastr*gbl->rho*kinf*omginf*cjcb;
                    res(3)(i,j) -= RAD(crd(0)(i,j))*susomg*betastr*gbl->rho*omginf*omginf/omg*cjcb;
                      
#ifdef CALC_TAU1
                    jcbmin = MIN(jcbmin,cjcb);
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
                        
                    FLT q = pow(u(0)(i,j)-0.5*mvel(0)(i,j),2.0)  +pow(u(1)(i,j)-0.5*mvel(1)(i,j),2.0);
                    qmax = MAX(qmax,q);
                        
                    FLT q2 = pow(u(0)(i,j)-mvel(0)(i,j),2.0)  +pow(u(1)(i,j)-mvel(1)(i,j),2.0);
                    qmax2 = MAX(qmax2,q2);
                    
                    tmuav = tmuav+tmu;
                    mutldav = mutldav+mutld;
#endif
				}
			}
			for(int n=0;n<NV;++n)
				basis::tri(log2p)->intgrt(&lf_re(n)(0),&res(n)(0,0),MXGP);

			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(int n=0;n<NV-1;++n) {
				basis::tri(log2p)->derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
				basis::tri(log2p)->derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
			}
			basis::tri(log2p)->derivr(&du(NV-1,0)(0,0),&res(NV-1)(0,0),MXGP);
			basis::tri(log2p)->derivs(&du(NV-1,1)(0,0),&res(NV-1)(0,0),MXGP);

#ifdef CALC_TAU2
            FLT h = inscribedradius(tind)/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
#endif
            

            
#ifdef CALC_TAU1
            hmax = 2.*sqrt(hmax);
            h = 4.*jcbmin/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
            hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
            
            tmuav = tmuav/lgpx/lgpn;
            mutldav = mutldav/lgpx/lgpn;
            FLT nu = (lmu + tmuav)/gbl->rho;
            FLT nuk = (lmu + sgmk*mutldav)/gbl->rho;
            FLT nuomg = (lmu + sgmomg*tmuav)/gbl->rho;
            FLT gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
            if (tmuav + gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
            
            FLT q2 = sqrt(qmax2);
            FLT lam2k  = (q2 +1.5*nuk/h +hmax*gbl->bd(0));
            FLT lam2omg  = (q2 +1.5*nuomg/h +hmax*gbl->bd(0));
            
            
            /* SET UP DISSIPATIVE COEFFICIENTS */
            gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
            gbl->tau(tind,2)  = adis*h/(jcb*lam2k);
            gbl->tau(tind,3)  = adis*h/(jcb*lam2omg);
            gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);
#endif
            
			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(int i=0;i<lgpx;++i) {
				for(int j=0;j<lgpn;++j) {

#ifdef CALC_TAU2
                    FLT q = pow(u(0)(i,j)-0.5*mvel(0)(i,j),2.0) +pow(u(1)(i,j)-0.5*mvel(1)(i,j),2.0);
                    FLT q2 = pow(u(0)(i,j)-mvel(0)(i,j),2.0) +pow(u(1)(i,j)-mvel(1)(i,j),2.0);
                    
                    psiktld = psifunc(u(2)(i,j),0.0,epslnk);
                    psinktld = psifunc(-u(2)(i,j),0.0,epslnk);
                    ktrb = psiktld*u(2)(i,j);
                    omg = exp(u(3)(i,j));
                    tmu = gbl->rho*ktrb/omg;
                    mutld = tmu - gbl->rho*psinktld*u(2)(i,j)/omginf;

                    FLT nu = (gbl->mu +tmu)/gbl->rho;
                    FLT nuk = (lmu +sgmk*mutld)/gbl->rho;
                    FLT nuomg = (lmu +sgmomg*tmu)/gbl->rho;
                    
                    FLT gam = 3.0*q +(0*0.5*h*gbl->bd(0) +2.*nu/h)*(0*0.5*h*gbl->bd(0) +2.*nu/h);
                    if (tmu + gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
                    FLT lam2k  = sqrt(q2) +1.5*nuk/h +h*gbl->bd(0);
                    FLT lam2omg  = sqrt(q2) +1.5*nuomg/h +h*gbl->bd(0);
                    
                    /* SET UP DISSIPATIVE COEFFICIENTS */
                    gbl->tau(tind,0) = adis*h/(cjcb*sqrt(gam));
                    gbl->tau(tind,2)  = adis*h/(cjcb*lam2k);
                    gbl->tau(tind,3)  = adis*h/(cjcb*lam2omg);
                    gbl->tau(tind,NV-1) = sqrt(q)*gbl->tau(tind,0);
#endif
                    

					tres(0) = gbl->tau(tind,0)*res(0)(i,j);
					tres(1) = gbl->tau(tind,0)*res(1)(i,j);
                    tres(2) = gbl->tau(tind,2)*res(2)(i,j);
                    tres(3) = gbl->tau(tind,3)*res(3)(i,j);
					tres(NV-1) = gbl->tau(tind,NV-1)*res(NV-1)(i,j);

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
                    
                    df(2,0)(i,j) -= (dcrd(1,1)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                     -dcrd(0,1)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(2);
                    
                    df(2,1)(i,j) -= (-dcrd(1,0)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                     +dcrd(0,0)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(2);
                    
                    df(3,0)(i,j) -= (dcrd(1,1)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                     -dcrd(0,1)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(3);
                    
                    df(3,1)(i,j) -= (-dcrd(1,0)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                     +dcrd(0,0)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(3);
                    
					du(NV-1,0)(i,j) = -(dcrd(1,1)(i,j)*tres(0) -dcrd(0,1)(i,j)*tres(1));
					du(NV-1,1)(i,j) = -(-dcrd(1,0)(i,j)*tres(0) +dcrd(0,0)(i,j)*tres(1));
				}
			}
			for(int n=0;n<NV-1;++n)
				basis::tri(log2p)->intgrtrs(&lf_re(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
			basis::tri(log2p)->intgrtrs(&lf_re(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);

			for(int n=0;n<NV;++n)
				for(int i=0;i<basis::tri(log2p)->tm();++i)
					lf_re(n)(i) *= gbl->beta(stage);

		}
	}
	else {
		/* LINEAR ELEMENT */
		/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {

				fluxx = gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
				fluxy = gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));

				/* CONTINUITY EQUATION FLUXES */
				du(NV-1,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
				du(NV-1,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;

				/* CONVECTIVE FLUXES */
				for(int n=0;n<NV-1;++n) {
					cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
					cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
				}

				/* PRESSURE TERMS */
				/* U-MOMENTUM */
				cv(0,0)(i,j) += ldcrd(1,1)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				cv(0,1)(i,j) -= ldcrd(1,0)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				/* V-MOMENTUM */
				cv(1,0)(i,j) -=  ldcrd(0,1)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
				cv(1,1)(i,j) +=  ldcrd(0,0)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
			}
		}
		for(int n=0;n<NV-1;++n)
			basis::tri(log2p)->intgrtrs(&lf_im(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
		basis::tri(log2p)->intgrtrs(&lf_im(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);


		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {

            cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
			lrhobd0 = gbl->rho*gbl->bd(0)*cjcb;
           

			/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
            /* INDICES ARE 1: EQUATION U, V, K-TILDE OR ln(OMEGA) 2: VARIABLE (U, V, K-TILDE OR ln(OMEGA) ), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
			visc(0,0)(0,0) = -(2.*ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
			visc(0,0)(1,1) = -(2.*ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
			visc(0,0)(0,1) =  (2.*ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define         viscI0II0II1II0I visc(0,0)(0,1)

			visc(1,1)(0,0) = -(ldcrd(1,1)*ldcrd(1,1) +2.*ldcrd(0,1)*ldcrd(0,1));
			visc(1,1)(1,1) = -(ldcrd(1,0)*ldcrd(1,0) +2.*ldcrd(0,0)*ldcrd(0,0));
			visc(1,1)(0,1) =  (ldcrd(1,1)*ldcrd(1,0) +2.*ldcrd(0,1)*ldcrd(0,0));
#define         viscI1II1II1II0I visc(1,1)(0,1)
            
            visc(0,1)(0,0) =  ldcrd(0,1)*ldcrd(1,1);
            visc(0,1)(1,1) =  ldcrd(0,0)*ldcrd(1,0);
            visc(0,1)(0,1) = -ldcrd(0,1)*ldcrd(1,0);
            visc(0,1)(1,0) = -ldcrd(0,0)*ldcrd(1,1);
            
            /* OTHER SYMMETRIES     */
#define         viscI1II0II0II0I visc(0,1)(0,0)
#define         viscI1II0II1II1I visc(0,1)(1,1)
#define         viscI1II0II0II1I visc(0,1)(1,0)
#define         viscI1II0II1II0I visc(0,1)(0,1)
            
            /*K-TILDE DIFFUSION */
            visc(2,2)(0,0) = -(ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
            visc(2,2)(1,1) = -(ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
            visc(2,2)(0,1) =  (ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define         viscI2II2II1II0I visc(2,2)(0,1)
            
            /*ln(OMEGA) DIFFUSION */
            visc(3,3)(0,0) = visc(2,2)(0,0);
            visc(3,3)(1,1) = visc(2,2)(1,1);
            visc(3,3)(0,1) = visc(2,2)(0,1);
#define         viscI3II3II1II0I visc(3,3)(0,1)

            
#ifdef CALC_TAU1
            jcbmin = MIN(jcbmin,cjcb);
            /* CALCULATE CURVED SIDE LENGTHS */
            h = 0.0;
            for (int n=0;n<ND;++n)
                h += ldcrd(n,0)*ldcrd(n,0);
            hmax = MAX(h,hmax);
            
            h = 0.0;
            for (int n=0;n<ND;++n)
                h += ldcrd(n,1)*ldcrd(n,1);
            hmax = MAX(h,hmax);
            
            h = 0.0;
            for (int n=0;n<ND;++n)
                h += (ldcrd(n,1) -ldcrd(n,0))*(ldcrd(n,1) -ldcrd(n,0));
            hmax = MAX(h,hmax);
            
            FLT tmuav = 0.0, mutldav = 0.0;
#endif

			/* TIME DERIVATIVE TERMS */ 
			for(int i=0;i<lgpx;++i) {
				for(int j=0;j<lgpn;++j) {
					rhorbd0 = RAD(crd(0)(i,j))*lrhobd0;
                    psiktld = psifunc(u(2)(i,j),0.0,epslnk);
                    psinktld = psifunc(-u(2)(i,j),0.0,epslnk);
                    ktrb = psiktld*u(2)(i,j);
                    omg = exp(u(3)(i,j));
                    tmu = gbl->rho*ktrb/omg;
                    FLT dudx = (ldcrd(1,1)*du(0,0)(i,j) -ldcrd(1,0)*du(0,1)(i,j))/cjcb;
                    FLT dvdx = (ldcrd(1,1)*du(1,0)(i,j) -ldcrd(1,0)*du(1,1)(i,j))/cjcb;
                    FLT dudy = (ldcrd(0,0)*du(0,1)(i,j) -ldcrd(0,1)*du(0,0)(i,j))/cjcb;
                    FLT dvdy = (ldcrd(0,0)*du(1,1)(i,j) -ldcrd(0,1)*du(1,0)(i,j))/cjcb;
                    strninv = dudx*dudx +0.5*(dudy +dvdx)*(dudy +dvdx) +dvdy*dvdy;
#ifdef AXISYMMETRIC
                    strninv += u(0)(i,j)*u(0)(i,j)/crd(0)(i,j)/crd(0)(i,j);
#endif
#ifndef WILCOX2006
                    omgMx = omg;
                    tmuLmtd = tmu;
#else
                    /* STRESS-LIMITER MODIFICATION */
                    omgLmtr = Clim*sqrt(2.0*strninv/betastr);
                    omgMx = std::max(omg,omgLmtr);
                    tmuLmtd = gbl->rho*ktrb/omgMx;
#endif
                    mutld = tmu -gbl->rho*psinktld*u(2)(i,j)/omginf;
                    cjcbi = (lmu +tmuLmtd)*RAD(crd(0)(i,j))/cjcb;
                    cjcbik = (lmu +sgmk*mutld)*RAD(crd(0)(i,j))/cjcb;
                    cjcbiomg = (lmu +sgmomg*tmu)*RAD(crd(0)(i,j))/cjcb;
                  
					/* UNSTEADY TERMS */
					for(int n=0;n<NV-1;++n)
						res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p)(tind,n,i,j);
					res(NV-1)(i,j) = rhorbd0 +dugdt(log2p)(tind,NV-1,i,j);

#ifdef AXISYMMETRIC
					res(0)(i,j) -= cjcb*(u(NV-1)(i,j) -2.*(lmu +tmuLmtd)*u(0)(i,j)/crd(0)(i,j) +k_mom*2./3.*gbl->rho*ktrb);
#endif
#ifdef BODYFORCE
					res(0)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(0);
					res(1)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(1);
#endif
                    
#ifdef MMS
                    /* source terms for MMS */
                    pt(0) = crd(0)(i,j);
                    pt(1) = crd(1)(i,j);
                    for(int n = 0; n < NV; ++n)
                        res(n)(i,j) -= cjcb*gbl->src->f(n,pt,gbl->time);
#endif
                    df(0,0)(i,j) = cjcbi*(+visc(0,0)(0,0)*du(0,0)(i,j) +visc(0,1)(0,0)*du(1,0)(i,j)
                    +visc(0,0)(0,1)*du(0,1)(i,j) +visc(0,1)(0,1)*du(1,1)(i,j)) +k_mom*RAD(crd(0)(i,j))*2./3.*gbl->rho*ldcrd(1,1)*ktrb;
                    
                    df(0,1)(i,j) = cjcbi*(+viscI0II0II1II0I*du(0,0)(i,j) +visc(0,1)(1,0)*du(1,0)(i,j)
                    +visc(0,0)(1,1)*du(0,1)(i,j) +visc(0,1)(1,1)*du(1,1)(i,j)) -k_mom*RAD(crd(0)(i,j))*2./3.*gbl->rho*ldcrd(1,0)*ktrb;
                    
                    df(1,0)(i,j) = cjcbi*(+viscI1II0II0II0I*du(0,0)(i,j) +visc(1,1)(0,0)*du(1,0)(i,j)
                    +viscI1II0II0II1I*du(0,1)(i,j) +visc(1,1)(0,1)*du(1,1)(i,j)) -k_mom*RAD(crd(0)(i,j))*2./3.*gbl->rho*ldcrd(0,1)*ktrb;
                    
                    df(1,1)(i,j) = cjcbi*(+viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
                    +viscI1II0II1II1I*du(0,1)(i,j) +visc(1,1)(1,1)*du(1,1)(i,j)) +k_mom*RAD(crd(0)(i,j))*2./3.*gbl->rho*ldcrd(0,0)*ktrb;
                    
                    df(2,0)(i,j) = cjcbik*(visc(2,2)(0,0)*du(2,0)(i,j) +visc(2,2)(0,1)*du(2,1)(i,j));
                    df(2,1)(i,j) = cjcbik*(viscI2II2II1II0I*du(2,0)(i,j) +visc(2,2)(1,1)*du(2,1)(i,j));
                    df(3,0)(i,j) = cjcbiomg*(visc(3,3)(0,0)*du(3,0)(i,j) +visc(3,3)(0,1)*du(3,1)(i,j));
                    df(3,1)(i,j) = cjcbiomg*(viscI3II3II1II0I*du(3,0)(i,j) +visc(3,3)(1,1)*du(3,1)(i,j));

					for(int n=0;n<NV-1;++n) {
						cv(n,0)(i,j) += df(n,0)(i,j);
						cv(n,1)(i,j) += df(n,1)(i,j);
                    }
                    
                    /*PRODUCTION TERMS FOR K-TILDE AND ln(omega) */
#ifndef WILCOX1988_KL
                    res(2)(i,j) -= RAD(crd(0)(i,j))*tmuLmtd*2.0*strninv*cjcb;
                    res(3)(i,j) -= RAD(crd(0)(i,j))*gamma*gbl->rho*2.0*strninv/omgMx*cjcb;
#else
                    /* KATO-LAUNDER */
                    vrtctinv = 0.5*(dudy -dvdx)*(dudy -dvdx);
                    res(2)(i,j) -= RAD(crd(0)(i,j))*2.0*tmu*sqrt(vrtctinv*strninv)*cjcb;
                    res(3)(i,j) -= RAD(crd(0)(i,j))*gamma*gbl->rho*2.0*sqrt(vrtctinv*strninv)/omg*cjcb;
#endif
                    /* PRODUCTION TERM FOR ln(OMEGA) (DUE TO LOGARITHMIC TRANSFORMATION) */
                    FLT domgtlddx = (ldcrd(1,1)*du(3,0)(i,j) -ldcrd(1,0)*du(3,1)(i,j))/cjcb;
                    FLT domgtlddy = (ldcrd(0,0)*du(3,1)(i,j) -ldcrd(0,1)*du(3,0)(i,j))/cjcb;
                    res(3)(i,j) -= cjcbiomg*(domgtlddx*domgtlddx + domgtlddy*domgtlddy)*cjcb*cjcb;
                    /* ROUND JET CORRECTION FOR WILCOX2006 */
#ifdef WILCOX2006
#ifdef AXISYMMETRIC
                    FLT chi = fabs(0.25*(dudy -dvdx)*(dudy -dvdx)*u(0)(i,j)/crd(0)(i,j)/pow(betastr*omg,3.0));
                    FLT f_beta = (1.0 +85.0*chi)/(1.0 +100.0*chi);
                    betakomg *= f_beta;
#endif
#endif
                    /* DISSIPATION TERMS FOR K-TILDE and ln(OMEGA)  */
                    res(2)(i,j) += RAD(crd(0)(i,j))*betastr*gbl->rho*(omg*ktrb + omginf*psinktld*u(2)(i,j))*cjcb;
                    res(3)(i,j) += RAD(crd(0)(i,j))*betakomg*gbl->rho*omg*cjcb;
                    
#ifdef WILCOX2006
                    /* CROSS-DIFFUSION TERM FOR ln(OMEGA) */
                    dktlddx = (ldcrd(1,1)*du(2,0)(i,j) -ldcrd(1,0)*du(2,1)(i,j))/cjcb;
                    dktlddy = (ldcrd(0,0)*du(2,1)(i,j) -ldcrd(0,1)*du(2,0)(i,j))/cjcb;
                    dkdx = dpsifunc(u(2)(i,j),0,epslnk)*dktlddx*u(2)(i,j) +psifunc(u(2)(i,j),0,epslnk)*dktlddx;
                    dkdy = dpsifunc(u(2)(i,j),0,epslnk)*dktlddy*u(2)(i,j) +psifunc(u(2)(i,j),0,epslnk)*dktlddy;
                    CrssD = dkdx*domgtlddx +dkdy*domgtlddy;
                    if ( CrssD > 0.0)
                        res(3)(i,j) -= RAD(crd(0)(i,j))*gbl->rho*CrssD/omg*sgmdo*cjcb;
#endif
                    
                    /* TURBULENCE SUSTAINING TERMS */
                    res(2)(i,j) -= RAD(crd(0)(i,j))*susk*betastr*gbl->rho*kinf*omginf*cjcb;
                    res(3)(i,j) -= RAD(crd(0)(i,j))*susomg*betastr*gbl->rho*omginf*omginf/omg*cjcb;

#ifdef CALC_TAU1
                    FLT q = pow(u(0)(i,j)-0.5*mvel(0)(i,j),2.0)  +pow(u(1)(i,j)-0.5*mvel(1)(i,j),2.0);
                    qmax = MAX(qmax,q);
                    
                    FLT q2 = pow(u(0)(i,j)-mvel(0)(i,j),2.0)  +pow(u(1)(i,j)-mvel(1)(i,j),2.0);
                    qmax2 = MAX(qmax2,q2);
                    
                    tmuav = tmuav+tmu;
                    mutldav = mutldav+mutld;
#endif
					 
				}
			}
			for(int n=0;n<NV;++n)
				basis::tri(log2p)->intgrt(&lf_re(n)(0),&res(n)(0,0),MXGP);

			/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
			for(int n=0;n<NV-1;++n) {
				basis::tri(log2p)->derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
				basis::tri(log2p)->derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
			}
			basis::tri(log2p)->derivr(&du(NV-1,0)(0,0),&res(NV-1)(0,0),MXGP);
			basis::tri(log2p)->derivs(&du(NV-1,1)(0,0),&res(NV-1)(0,0),MXGP);
            
#ifdef CALC_TAU2
            FLT h = inscribedradius(tind)/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
#endif
            

            
#ifdef CALC_TAU1
            hmax = 2.*sqrt(hmax);
            h = 4.*jcbmin/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
            hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));

            tmuav = tmuav/lgpx/lgpn;
            mutldav = mutldav/lgpx/lgpn;
            FLT nu = (lmu + tmuav)/gbl->rho;
            FLT nuk = (lmu + sgmk*mutldav)/gbl->rho;
            FLT nuomg = (lmu + sgmomg*tmuav)/gbl->rho;
            FLT gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
            if (tmuav + gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
            
            FLT q2 = sqrt(qmax2);
            FLT lam2k  = (q2 +1.5*nuk/h +hmax*gbl->bd(0));
            FLT lam2omg  = (q2 +1.5*nuomg/h +hmax*gbl->bd(0));
            
            /* SET UP DISSIPATIVE COEFFICIENTS */
            gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
            gbl->tau(tind,2)  = adis*h/(jcb*lam2k);
            gbl->tau(tind,3)  = adis*h/(jcb*lam2omg);
            gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);
#endif

			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(int i=0;i<lgpx;++i) {
				for(int j=0;j<lgpn;++j) {
                    
#ifdef CALC_TAU2
                    FLT q = pow(u(0)(i,j)-0.5*mvel(0)(i,j),2.0) +pow(u(1)(i,j)-0.5*mvel(1)(i,j),2.0);
                    FLT q2 = pow(u(0)(i,j)-mvel(0)(i,j),2.0) +pow(u(1)(i,j)-mvel(1)(i,j),2.0);
                    psiktld = psifunc(u(2)(i,j),0.0,epslnk);
                    psinktld = psifunc(-u(2)(i,j),0.0,epslnk);
                    ktrb = psiktld*u(2)(i,j);
                    omg = exp(u(3)(i,j));
                    tmu = gbl->rho*ktrb/omg;
                    mutld = tmu - gbl->rho*psinktld*u(2)(i,j)/omginf;
                    FLT nu = (gbl->mu +tmu)/gbl->rho;
                    FLT nuk = (lmu +sgmk*mutld)/gbl->rho;
                    FLT nuomg = (lmu +sgmomg*tmu)/gbl->rho;
                    
                    FLT gam = 3.0*q +(0*0.5*h*gbl->bd(0) +2.*nu/h)*(0*0.5*h*gbl->bd(0) +2.*nu/h);
                    if (tmu + gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
                    FLT lam2k  = sqrt(q2) +1.5*nuk/h +h*gbl->bd(0);
                    FLT lam2omg  = sqrt(q2) +1.5*nuomg/h +h*gbl->bd(0);
                    
                    /* SET UP DISSIPATIVE COEFFICIENTS */
                    gbl->tau(tind,0) = adis*h/(cjcb*sqrt(gam));
                    gbl->tau(tind,2)  = adis*h/(cjcb*lam2k);
                    gbl->tau(tind,3)  = adis*h/(cjcb*lam2omg);
                    gbl->tau(tind,NV-1) = sqrt(q)*gbl->tau(tind,0);
#endif

                    tres(0) = gbl->tau(tind,0)*res(0)(i,j);
                    tres(1) = gbl->tau(tind,0)*res(1)(i,j);
                    tres(2) = gbl->tau(tind,2)*res(2)(i,j);
                    tres(3) = gbl->tau(tind,3)*res(3)(i,j);
                    tres(NV-1) = gbl->tau(tind,NV-1)*res(NV-1)(i,j);

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

                    df(2,0)(i,j) -= (ldcrd(1,1)*(u(0)(i,j)-mvel(0)(i,j))
                                     -ldcrd(0,1)*(u(1)(i,j)-mvel(1)(i,j)))*tres(2);
                    
                    df(2,1)(i,j) -= (-ldcrd(1,0)*(u(0)(i,j)-mvel(0)(i,j))
                                     +ldcrd(0,0)*(u(1)(i,j)-mvel(1)(i,j)))*tres(2);
                    
                    df(3,0)(i,j) -= (ldcrd(1,1)*(u(0)(i,j)-mvel(0)(i,j))
                                     -ldcrd(0,1)*(u(1)(i,j)-mvel(1)(i,j)))*tres(3);
                    
                    df(3,1)(i,j) -= (-ldcrd(1,0)*(u(0)(i,j)-mvel(0)(i,j))
                                     +ldcrd(0,0)*(u(1)(i,j)-mvel(1)(i,j)))*tres(3);

					du(NV-1,0)(i,j) = -(ldcrd(1,1)*tres(0) -ldcrd(0,1)*tres(1));
					du(NV-1,1)(i,j) = -(-ldcrd(1,0)*tres(0) +ldcrd(0,0)*tres(1));
				}
			}
			for(int n=0;n<NV-1;++n)
				basis::tri(log2p)->intgrtrs(&lf_re(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
			basis::tri(log2p)->intgrtrs(&lf_re(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);

			for(int n=0;n<NV;++n)
				for(int i=0;i<basis::tri(log2p)->tm();++i)
					lf_re(n)(i) *= gbl->beta(stage);

		}
	}
	
	return;
}
