/*
 *  allocate.cpp
 *  tri_hp
 *
 *  Thermoelastic plane strain: tri_hp_ps with an added temperature variable
 *
 */

#include "tri_hp_thermoelastic.h"
#include "../hp_boundary.h"

void tri_hp_thermoelastic::init(input_map& inmap, shared_ptr<block_global> gin) {
	std::string keyword,ibcname;

	gbl = gin;
	hp_te_gbl = make_shared<hp_thermoelastic_global>();

	/* VARIABLES ARE U, V, T, P */
	inmap[gbl->idprefix + "_nvariable"] = "4";
	tri_hp_ps::init(inmap,gin);

	/* ELASTIC CONSTANTS (MU & LAMI) ARE LOADED BY tri_hp_ps */
	inmap.getwdefault(gbl->idprefix + "_E",hp_te_gbl->E,1.0);

	/* THERMAL EXPANSION */
	if (!inmap.get(gbl->idprefix + "_alpha",hp_te_gbl->alpha)) inmap.getwdefault("alpha",hp_te_gbl->alpha,0.0);
	if (!inmap.get(gbl->idprefix + "_T_ref",hp_te_gbl->T_ref)) inmap.getwdefault("T_ref",hp_te_gbl->T_ref,0.0);

	/* HEAT EQUATION CONSTANTS (SAME KEYWORDS AS cd EXCEPT nu, WHICH IS POISSON'S RATIO HERE) */
	FLT rho;
	if (!inmap.get(gbl->idprefix + "_rho",rho)) inmap.getwdefault("rho",rho,1.0);
	FLT cv;
	if (!inmap.get(gbl->idprefix + "_cv",cv)) inmap.getwdefault("cv",cv,1.0);
	hp_te_gbl->rhocv = rho*cv;

	if (!inmap.get(gbl->idprefix + "_conductivity",hp_te_gbl->kcond)) inmap.getwdefault("conductivity",hp_te_gbl->kcond,0.0);

	/* CONVECTION VELOCITY (DEFAULT IS A SOLID AT REST) */
	if (!inmap.get(gbl->idprefix + "_ax",hp_te_gbl->ax)) inmap.getwdefault("ax",hp_te_gbl->ax,0.0);
	if (!inmap.get(gbl->idprefix + "_ay",hp_te_gbl->ay)) inmap.getwdefault("ay",hp_te_gbl->ay,0.0);

	hp_te_gbl->tau.resize(maxpst);

	/* HEAT SOURCE */
	keyword = gbl->idprefix + "_src";
	if (!inmap.get(keyword,ibcname)) {
		keyword = "src";
		if (!inmap.get(keyword,ibcname)) {
			ibcname = "zero";
		}
	}
	hp_te_gbl->src = getnewibc(ibcname);
	hp_te_gbl->src->init(inmap,keyword);

	return;
}

void tri_hp_thermoelastic::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {

	const tri_hp_thermoelastic& inmesh = dynamic_cast<const tri_hp_thermoelastic &>(in);
	gbl = inmesh.gbl;
	hp_te_gbl = inmesh.hp_te_gbl;

	tri_hp_ps::init(in,why,sizereduce1d);

	return;
}

/* UNSTEADY SOURCES FOR THE TEMPERATURE EQUATION ONLY (ELASTICITY IS QUASI-STATIC) */
void tri_hp_thermoelastic::calculate_unsteady_sources() {
	int i,j,n,tind;
	const int T = 2;
#ifdef petsc
	int start = log2pmax;
#else
	int start = 0;
#endif

	for (log2p=start;log2p<=log2pmax;++log2p) {
		for(tind=0;tind<ntri;++tind) {
			if (tri(tind).info > -1) {
				crdtocht(tind,1);
				for(n=0;n<ND;++n)
					basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
			}
			else {
				for(n=0;n<ND;++n)
					basis::tri(log2p)->proj(vrtxbd(1)(tri(tind).pnt(0))(n),vrtxbd(1)(tri(tind).pnt(1))(n),vrtxbd(1)(tri(tind).pnt(2))(n),&crd(n)(0,0),MXGP);

				for(i=0;i<basis::tri(log2p)->gpx();++i) {
					for(j=0;j<basis::tri(log2p)->gpn();++j) {
						for(n=0;n<ND;++n) {
							dcrd(n,0)(i,j) = 0.5*(vrtxbd(1)(tri(tind).pnt(1))(n) -vrtxbd(1)(tri(tind).pnt(0))(n));
							dcrd(n,1)(i,j) = 0.5*(vrtxbd(1)(tri(tind).pnt(2))(n) -vrtxbd(1)(tri(tind).pnt(0))(n));
						}
					}
				}
			}

			ugtouht(tind,1);
			basis::tri(log2p)->proj(&uht(T)(0),&u(T)(0,0),MXGP);

			for(i=0;i<basis::tri(log2p)->gpx();++i) {
				for(j=0;j<basis::tri(log2p)->gpn();++j) {
					cjcb(i,j) = -gbl->bd(0)*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
					for(n=0;n<NV;++n)
						dugdt(log2p)(tind,n,i,j) = 0.0;
					dugdt(log2p)(tind,T,i,j) = hp_te_gbl->rhocv*u(T)(i,j)*cjcb(i,j);

					for(n=0;n<ND;++n)
						dxdt(log2p)(tind,n,i,j) = crd(n)(i,j);
				}
			}
		}
	}
	log2p=log2pmax;

	return;
}
