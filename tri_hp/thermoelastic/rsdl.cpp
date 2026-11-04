/*
 *  rsdl.cpp
 *  tri_hp
 *
 *  Element residual for thermoelastic plane strain.  The elasticity
 *  terms follow planestrain/rsdl.cpp and the temperature terms follow
 *  cd/rsdl.cpp.  Thermal strain enters through the constitutive
 *  (continuity) equation for p:
 *     lami*p +div(u) = alpha*(3 +2*mu*lami)*(T -T_ref)
 *  i.e. p = -lambda*div(u) +(3*lambda +2*mu)*alpha*(T -T_ref)
 *
 */

#include "tri_hp_thermoelastic.h"
#include "../hp_boundary.h"

void tri_hp_thermoelastic::element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im) {
	int i,j,n;
	FLT fluxx,fluxy,rad,cjcbi,kc;
	const int NV = 4, T = 2, P = 3;
	const int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
	const FLT lmu = hp_ps_gbl->mu, lami = hp_ps_gbl->lami;
	/* (3*lambda +2*mu)*alpha/lambda */
	const FLT thermal = hp_te_gbl->alpha*(3.0 +2.0*lmu*lami);
	TinyMatrix<TinyMatrix<FLT,ND,ND>,ND,ND> visc;
	TinyMatrix<FLT,ND,ND> diff;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> cv, df;
	TinyVector<FLT,NV> tres;
	TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> mvel; // for local mesh velocity info

	/* IF TINFO > -1 IT IS CURVED ELEMENT */
	if (tri(tind).info > -1) {
		pmetric->calc_metrics(tind,crd,dcrd);
	}
	else {
		/* LOAD INDICES OF VERTEX POINTS */
		TinyVector<int,3> v;
		v = tri(tind).pnt;

		/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(n=0;n<ND;++n)
			basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);

		/* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpn;++j) {
				for(n=0;n<ND;++n) {
					dcrd(n,0)(i,j) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
					dcrd(n,1)(i,j) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
				}
			}
		}
	}

	/* CALCULATE MESH VELOCITY */
	for(i=0;i<lgpx;++i) {
		for(j=0;j<lgpn;++j) {
			mvel(0)(i,j) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
			mvel(1)(i,j) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j));
#ifdef MESH_REF_VEL
			mvel(0)(i,j) += hp_gbl->mesh_ref_vel(0);
			mvel(1)(i,j) += hp_gbl->mesh_ref_vel(1);
#endif
		}
	}

	/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES FOR DISPLACEMENTS & TEMPERATURE */
	for(n=0;n<NV-1;++n)
		basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
	basis::tri(log2p)->proj(&uht(P)(0),&u(P)(0,0),MXGP);

	/* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
	for(n=0;n<NV;++n){
		for(i=0;i<basis::tri(log2p)->tm();++i){
			lf_re(n)(i) = 0.0;
			lf_im(n)(i) = 0.0;
		}
	}

	/* CONVECTIVE TERMS (IMAGINARY FIRST)*/
	for(i=0;i<lgpx;++i) {
		for(j=0;j<lgpn;++j) {
			rad = RAD(crd(0)(i,j));

			/* CONTINUITY EQUATION FLUXES */
			fluxx = rad*u(0)(i,j);
			fluxy = rad*u(1)(i,j);
			cv(P,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
			cv(P,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;

			/* PRESSURE TERMS */
			/* U-MOMENTUM */
			cv(0,0)(i,j) = +dcrd(1,1)(i,j)*rad*u(P)(i,j);
			cv(0,1)(i,j) = -dcrd(1,0)(i,j)*rad*u(P)(i,j);
			/* V-MOMENTUM */
			cv(1,0)(i,j) = -dcrd(0,1)(i,j)*rad*u(P)(i,j);
			cv(1,1)(i,j) = +dcrd(0,0)(i,j)*rad*u(P)(i,j);

			/* TEMPERATURE CONVECTION */
			fluxx = hp_te_gbl->rhocv*rad*(hp_te_gbl->ax -mvel(0)(i,j))*u(T)(i,j);
			fluxy = hp_te_gbl->rhocv*rad*(hp_te_gbl->ay -mvel(1)(i,j))*u(T)(i,j);
			cv(T,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
			cv(T,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
		}
	}
	for(n=0;n<NV;++n)
		basis::tri(log2p)->intgrtrs(&lf_im(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);

	/* NEGATIVE REAL TERMS */
	if (gbl->beta(stage) > 0.0) {

		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpn;++j) {
				TinyVector<FLT,ND> pt(crd(0)(i,j),crd(1)(i,j));
				rad = RAD(crd(0)(i,j));
				cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);

				/* SOURCE TERMS */
				res(0)(i,j) = 0.0;
				res(1)(i,j) = 0.0;
				/* UNSTEADY & HEAT SOURCE */
				res(T)(i,j) = hp_te_gbl->rhocv*rad*gbl->bd(0)*u(T)(i,j)*cjcb(i,j) +dugdt(log2p)(tind,T,i,j);
				res(T)(i,j) -= rad*cjcb(i,j)*hp_te_gbl->src->f(0,pt,gbl->time);
				/* CONSTITUTIVE EQUATION WITH THERMAL STRAIN */
				res(P)(i,j) = (lami*u(P)(i,j) -thermal*(u(T)(i,j) -hp_te_gbl->T_ref))*cjcb(i,j);

				/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
				/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
				cjcbi = lmu*rad/cjcb(i,j);
				visc(0,0)(0,0) = -cjcbi*(2.*dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
				visc(0,0)(1,1) = -cjcbi*(2.*dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
				visc(0,0)(0,1) =  cjcbi*(2.*dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define         viscI0II0II1II0I visc(0,0)(0,1)

				visc(1,1)(0,0) = -cjcbi*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +2.*dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
				visc(1,1)(1,1) = -cjcbi*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +2.*dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
				visc(1,1)(0,1) =  cjcbi*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +2.*dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define         viscI1II1II1II0I visc(1,1)(0,1)

				visc(0,1)(0,0) =  cjcbi*dcrd(0,1)(i,j)*dcrd(1,1)(i,j);
				visc(0,1)(1,1) =  cjcbi*dcrd(0,0)(i,j)*dcrd(1,0)(i,j);
				visc(0,1)(0,1) = -cjcbi*dcrd(0,1)(i,j)*dcrd(1,0)(i,j);
				visc(0,1)(1,0) = -cjcbi*dcrd(0,0)(i,j)*dcrd(1,1)(i,j);

				/* OTHER SYMMETRIES     */
#define         viscI1II0II0II0I visc(0,1)(0,0)
#define         viscI1II0II1II1I visc(0,1)(1,1)
#define         viscI1II0II0II1I visc(0,1)(1,0)
#define         viscI1II0II1II0I visc(0,1)(0,1)

				df(0,0)(i,j) = +visc(0,0)(0,0)*du(0,0)(i,j) +visc(0,1)(0,0)*du(1,0)(i,j)
								+visc(0,0)(0,1)*du(0,1)(i,j) +visc(0,1)(0,1)*du(1,1)(i,j);

				df(0,1)(i,j) = +viscI0II0II1II0I*du(0,0)(i,j) +visc(0,1)(1,0)*du(1,0)(i,j)
								+visc(0,0)(1,1)*du(0,1)(i,j) +visc(0,1)(1,1)*du(1,1)(i,j);

				df(1,0)(i,j) = +viscI1II0II0II0I*du(0,0)(i,j) +visc(1,1)(0,0)*du(1,0)(i,j)
								+viscI1II0II0II1I*du(0,1)(i,j) +visc(1,1)(0,1)*du(1,1)(i,j);

				df(1,1)(i,j) = +viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
								+viscI1II0II1II1I*du(0,1)(i,j) +visc(1,1)(1,1)*du(1,1)(i,j);

				/* DIFFUSION TENSOR FOR TEMPERATURE */
				kc = hp_te_gbl->kcond*rad/cjcb(i,j);
				diff(0,0) = -kc*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
				diff(1,1) = -kc*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
				diff(0,1) =  kc*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));

				df(T,0)(i,j) = +diff(0,0)*du(T,0)(i,j) +diff(0,1)*du(T,1)(i,j);
				df(T,1)(i,j) = +diff(0,1)*du(T,0)(i,j) +diff(1,1)*du(T,1)(i,j);

				for(n=0;n<NV-1;++n) {
					cv(n,0)(i,j) += df(n,0)(i,j);
					cv(n,1)(i,j) += df(n,1)(i,j);
				}
			}
		}
		for(n=0;n<NV;++n)
			basis::tri(log2p)->intgrt(&lf_re(n)(0),&res(n)(0,0),MXGP);

		/* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
		for(n=0;n<NV;++n) {
			basis::tri(log2p)->derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
			basis::tri(log2p)->derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
		}

		/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpn;++j) {
				/* SET UP DISSIPATIVE COEFFICIENTS */
				/* MOMENTUM RESIDUAL STABILIZES THE CONSTITUTIVE EQUATION */
				tres(0) = hp_ps_gbl->tau(tind)*res(0)(i,j);
				tres(1) = hp_ps_gbl->tau(tind)*res(1)(i,j);

				df(P,0)(i,j) = -(dcrd(1,1)(i,j)*tres(0) -dcrd(0,1)(i,j)*tres(1));
				df(P,1)(i,j) = -(-dcrd(1,0)(i,j)*tres(0) +dcrd(0,0)(i,j)*tres(1));

				/* STREAMWISE STABILIZATION OF TEMPERATURE */
				tres(T) = hp_te_gbl->tau(tind)*res(T)(i,j);
				df(T,0)(i,j) -= (dcrd(1,1)(i,j)*(hp_te_gbl->ax -mvel(0)(i,j))
								-dcrd(0,1)(i,j)*(hp_te_gbl->ay -mvel(1)(i,j)))*tres(T);
				df(T,1)(i,j) -= (-dcrd(1,0)(i,j)*(hp_te_gbl->ax -mvel(0)(i,j))
								+dcrd(0,0)(i,j)*(hp_te_gbl->ay -mvel(1)(i,j)))*tres(T);
			}
		}
		for(n=0;n<NV;++n)
			basis::tri(log2p)->intgrtrs(&lf_re(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);

		for(n=0;n<NV;++n)
			for(i=0;i<basis::tri(log2p)->tm();++i)
				lf_re(n)(i) *= gbl->beta(stage);
	}

	return;
}
