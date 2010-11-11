#include "tri_hp_cns.h"
#include "../hp_boundary.h"
#include <myblas.h>


void tri_hp_cns::project_new_variables(){
	int info,last_phase, mp_phase;
	char uplo[] = "U";
	Array<double,2> P(NV,NV);
	Array<double,1> lcl(NV);
	Array<TinyVector<double,MXGP>,1> u1d(NV),res1d(NV),temp1d(NV);
	Array<TinyMatrix<double,MXGP,MXGP>,1> u2d(NV),res2d(NV),temp2d(NV);
	Array<TinyVector<double,MXTM>,1> ucoef(NV),rcoef(NV),tcoef(NV);
	
	/* LOOP THROUGH VERTICES */
	for(int i=0;i<npnt;++i){
	
		for(int m = 0; m < NV; ++m)
			lcl(m) = ug.v(i,m);
		
		calculate_preconditioner(lcl, P);

		lcl = 0.0;
		for(int m = 0; m < NV; ++m)
			for(int n = 0; n < NV; ++n)							
				lcl(m) += P(m,n)*gbl->res.v(i,n);
		
		for(int m = 0; m < NV; ++m)
			gbl->res.v(i,m) = lcl(m);
		
// /* maybe forward substitution with Pinv??? */		
//		for(int j=0;j<NV;++j){
//			FLT lcl2 = gbl->res.v(i,j);
//			for(int k=0;k<j;++k){
//				lcl2 -= Pinv(j,k)*gbl->res.v(i,k);
//			}
//			gbl->res.v(i,j) = lcl2/Pinv(j,j);
//		}
		
		
	}
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		vc0load(mp_phase,gbl->res.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= vc0wait_rcv(mp_phase,gbl->res.v.data());
	}
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet();
	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet2d();
	
	if (basis::tri(log2p)->sm() <= 0) return;

	/* LOOP THROUGH SIDES */    
    for(int sind=0;sind<nseg;++sind) {
		
		/* project linears */
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->proj1d(gbl->res.v(seg(sind).pnt(0),n),gbl->res.v(seg(sind).pnt(1),n),&temp1d(n)(0));
		
		/* take global coefficients and put into local vector */
		for(int n=0;n<NV;++n) {
			for (int m=0; m<2; ++m) {
				ucoef(n)(m) = ug.v(seg(sind).pnt(m),n);
				rcoef(n)(m) = gbl->res_temp.v(seg(sind).pnt(m),n);
			}

		}
		
		for(int m=0;m<basis::tri(log2p)->sm();++m) {
			for(int n=0;n<NV;++n) {
				ucoef(n)(m+2) = ug.s(sind,m,n);
				rcoef(n)(m+2) = gbl->res_temp.s(sind,m,n);
			}
		}
		
		/* project sides */
		for(int n=0;n<NV;++n) {
			basis::tri(log2p)->proj1d(&ucoef(n)(0),&u1d(n)(0));
			basis::tri(log2p)->proj1d(&rcoef(n)(0),&res1d(n)(0));
		}
		
		for(int i=0;i<basis::tri(log2p)->gpx(); ++i) {
			
			for(int m = 0; m < NV; ++m)
				lcl(m) = u1d(m)(i);
				
			calculate_preconditioner(lcl, P);
			
			for(int m = 0; m < NV; ++m)
				for(int n = 0; n < NV; ++n)							
					temp1d(m)(i) -= P(m,n)*res1d(n)(i);
	
		}
		
		/* integrate right hand side: phi*P*res */
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->intgrt1d(&rcoef(n)(0),&temp1d(n)(0));
		
		/* invert 1d mass matrix */
		for(int n=0;n<NV;++n) {
			PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&rcoef(n)(2),basis::tri(log2p)->sm(),info);
			for(int m=0;m<basis::tri(log2p)->sm();++m) 
				gbl->res.s(sind,m,n) = -rcoef(n)(2+m);
		}
	}
	
	/* APPLY DIRCHLET B.C.S TO MODE */
	for(int i=0;i<nebd;++i)
		for(int m=0;m<basis::tri(log2p)->sm();++m) 
			hp_ebdry(i)->sdirichlet(m);
	
	if (basis::tri(log2p)->im() <= 0) return;
	
	for(int tind = 0; tind < ntri; ++tind) {
		
		for (int i=0; i<3; ++i) {
			int vrtx = tri(tind).pnt(i);
			for(int n=0; n<NV; ++n){
				ucoef(n)(i) = ug.v(vrtx,n);
				rcoef(n)(i) = gbl->res_temp.v(vrtx,n);
				tcoef(n)(i) = gbl->res.v(vrtx,n);
			}
		}
		
		/* SIDES */
		int cnt = 3;
		for(int i=0;i<3;++i) {
			int sind = tri(tind).seg(i);
			int sign = tri(tind).sgn(i);
			int msgn = 1;
			for (int m = 0; m < basis::tri(log2p)->sm(); ++m) {
				for(int n=0; n<NV; ++n){
					ucoef(n)(cnt) = msgn*ug.s(sind,m,n);
					rcoef(n)(cnt) = msgn*gbl->res_temp.s(sind,m,n);
					tcoef(n)(cnt) = msgn*gbl->res.s(sind,m,n);
				}
				msgn *= sign;
				++cnt;
			}
		}
		
		/* INTERIORS */    
		if (basis::tri(log2p)->im() > 0) {    
			int indx = 0;
			for(int m = 1; m < basis::tri(log2p)->sm(); ++m) {
				for(int k = 0; k < basis::tri(log2p)->sm()-m; ++k) {
					for(int n=0; n<NV; ++n){
						ucoef(n)(cnt) = ug.i(tind,indx,n);
						rcoef(n)(cnt) = gbl->res_temp.i(tind,indx,n);
						tcoef(n)(cnt) = 0.0;
					}
					++cnt; ++indx;
				}
				indx += sm0 -basis::tri(log2p)->sm();
			}
		}
		
		
		for(int n=0;n<NV;++n){			
			basis::tri(log2p)->proj(&ucoef(n)(0),&u2d(n)(0,0),MXGP);
			basis::tri(log2p)->proj(&rcoef(n)(0),&res2d(n)(0,0),MXGP);
			basis::tri(log2p)->proj_bdry(&tcoef(n)(0),&temp2d(n)(0,0),MXGP);
		}
		
		for (int i=0; i < basis::tri(log2p)->gpx(); ++i ) {
			for (int j=0; j < basis::tri(log2p)->gpn(); ++j ) {

				for(int m = 0; m < NV; ++m)
					lcl(m) = u2d(m)(i,j);
				
				calculate_preconditioner(lcl, P);
				
				for(int m = 0; m < NV; ++m)
					for(int n = 0; n < NV; ++n)							
						temp2d(m)(i,j) -= P(m,n)*res2d(n)(i,j);

			}
		}

		for(int n=0;n<NV;++n) {
			basis::tri(log2p)->intgrt(&rcoef(n)(0),&temp2d(n)(0,0),MXGP);
			DPBTRS(uplo,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),1,(double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,&rcoef(n)(basis::tri(log2p)->bm()),basis::tri(log2p)->im(),info);
			for(int i=0;i<basis::tri(log2p)->im();++i)
				gbl->res.i(tind,i,n) = -rcoef(n)(basis::tri(log2p)->bm()+i);
		}
	}
	
	return;
}

void tri_hp_cns::calculate_preconditioner(Array<double,1> pvu, Array<double,2> &P){
	
	double pr = pvu(0),u = pvu(1),v = pvu(2),rt = pvu(3);
	double rho = pr/rt;
	double ke = 0.5*(u*u+v*v);
	double gam = gbl->gamma;
	double gm1 = gam-1.0;
	
//	/* Preconditioner */
//	P = ke*gm1,          -u*gm1,     -v*gm1,      gm1,
//		-u/rho,          1.0/rho,    0.0,         0.0,
//		-v/rho,          0.0,        1.0/rho,     0.0,
//		(gm1*ke-rt)/rho, -u*gm1/rho, -v*gm1/rho, gm1/rho;
	
	FLT umag = sqrt(u*u+v*v);	
	FLT c = sqrt(gam*rt);
	FLT M = MIN(MAX(1.0e-5,umag/c),1.0);
	FLT delta = 1.0;
	FLT omega = gam-gm1*delta;
	FLT k = 1.0; // can use more complicated formula for k see choi and merkle
	FLT beta = k*gam*rt;
	FLT bM2 = beta*M*M;
	FLT E = rt/gm1+ke;
	
	/* Preconditioner */
	P = bM2,                                    0.0,              0.0,              0.0,
		-u/rho,                                 1.0/rho,          0.0,              0.0,
		-v/rho,                                 0.0,              1.0/rho,          0.0,
		gm1*(u*u+v*v-E-rt+delta*bM2)/(rho*gam), -u*gm1/(rho*gam), -v*gm1/(rho*gam), gm1/(rho*gam);	
	
	return;
}

void tri_hp_cns::minvrt() {
	
	if (gbl->diagonal_preconditioner) {

		if (basis::tri(log2p)->sm()) {
			*gbl->log << "MINVRT CNS only works for p=1" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		gbl->res.v(Range(0,npnt-1),Range::all()) *= gbl->vprcn(Range(0,npnt-1),Range::all());

		gbl->res_temp.v = gbl->res.v;
		gbl->res_temp.s = gbl->res.s;
		gbl->res_temp.i = gbl->res.i;
		
		project_new_variables();
	}
	else {
		tri_hp::minvrt();
	}
	
	return;
}


//void tri_hp_cns::minvrt() {
//	int i,j,k,m,n,tind,sind,v0,indx,indx1,indx2,sgn,msgn;
//	TinyVector<int,3> sign,side;
//	Array<FLT,2> tinv(NV,NV),P(NV,NV);
//	Array<FLT,1> temp(NV);
//	int last_phase, mp_phase;
//	FLT gm1 = gbl->gamma-1.0;
//	
//
//	/* only p = 1 for now */
//	/* uses Pinv */
//	for(i=0;i<npnt;++i) {
//		
//		gbl->res.v(i,Range::all()) *= gbl->vprcn(i,Range::all());
//
//		FLT pr = ug.v(i,0);
//		FLT u = ug.v(i,1);
//		FLT v = ug.v(i,2);
//		FLT rt = ug.v(i,3);
//		FLT rho = pr/rt;
//		FLT ke = 0.5*(u*u+v*v);
//		
//		/* Preconditioner */
//		P = ke*gm1,          -u*gm1,     -v*gm1,      gm1,
//		-u/rho,          1.0/rho,    0.0,         0.0,
//		-v/rho,          0.0,        1.0/rho,     0.0,
//		(gm1*ke-rt)/rho, -u*gm1/rho, -v*gm1/rho, gm1/rho;
//		
//		temp = 0.0;
//		for(int m = 0; m < NV; ++m)
//			for(int n = 0; n < NV; ++n)							
//				temp(m) += P(m,n)*gbl->res.v(i,n);
//		
//		for(int m = 0; m < NV; ++m)
//			gbl->res.v(i,m) = temp(m);
//		
//	}
//	
//	
//	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
//		vc0load(mp_phase,gbl->res.v.data());
//		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
//		last_phase = true;
//		last_phase &= vc0wait_rcv(mp_phase,gbl->res.v.data());
//	}
//	
//	/* APPLY VERTEX DIRICHLET B.C.'S */
//	for(i=0;i<nebd;++i)
//		hp_ebdry(i)->vdirichlet();
//	
//	for(i=0;i<nvbd;++i)
//		hp_vbdry(i)->vdirichlet2d();
//	
//	
//	return;
//}
//
