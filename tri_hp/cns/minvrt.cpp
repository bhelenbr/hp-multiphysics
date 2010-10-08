#include "tri_hp_cns.h"
#include "../hp_boundary.h"
#include <myblas.h>


void tri_hp_cns::project_new_variables(){
	int i,j,k,m,n,tind,sind,v0,v1,indx,cnt,sign,msgn,info;
	char uplo[] = "U";
	TinyVector<double,2> pt;
	Array<double,2> P(NV,NV);
	Array<double,1> temp(NV);
	Array<TinyVector<double,MXGP>,1> u1d(NV),res1d(NV),temp1d(NV);
	Array<TinyMatrix<double,MXGP,MXGP>,1> u2d(NV),res2d(NV),temp2d(NV);
	Array<TinyVector<double,MXTM>,1> ucoef(NV),rcoef(NV),tcoef(NV);
	
	
	/* LOOP THROUGH VERTICES */
	for(int i=0;i<npnt;++i){
	
		for(int m = 0; m < NV; ++m)
			temp(m) = gbl->ug0.v(i,m);
		
		calculate_preconditioner(temp, P);

		temp = 0.0;
		for(int m = 0; m < NV; ++m)
			for(int n = 0; n < NV; ++n)							
				temp(m) += P(m,n)*gbl->res.v(i,n);
		
		for(int m = 0; m < NV; ++m)
			gbl->res.v(i,m) = temp(m);
		
	}
	
	if (basis::tri(log2p)->sm() <= 0) return;

	/* LOOP THROUGH SIDES */    
    for(sind=0;sind<nseg;++sind) {
		
		v0 = seg(sind).pnt(0);
		v1 = seg(sind).pnt(1);		
		
		/* project linears*/
		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj1d(gbl->res.v(v0,n),gbl->res.v(v1,n),&temp1d(n)(0));
		
		/* take global coefficients and put into local vector */
		for(n=0;n<NV;++n) {
			ucoef(n)(0) = gbl->ug0.v(v0,n);
			ucoef(n)(1) = gbl->ug0.v(v1,n);
			rcoef(n)(0) = gbl->res.v(v0,n);
			rcoef(n)(1) = gbl->res.v(v1,n);
		}
		
		for(m=0;m<basis::tri(log2p)->sm();++m) {
			for(n=0;n<NV;++n) {
				ucoef(n)(m+2) = gbl->ug0.s(sind,m,n);
				rcoef(n)(m+2) = gbl->res.s(sind,m,n);
			}
		}
		
		/* project sides */
		for(n=0;n<NV;++n) {
			basis::tri(log2p)->proj1d(&ucoef(n)(0),&u1d(n)(0));
			basis::tri(log2p)->proj1d(&rcoef(n)(0),&res1d(n)(0));
		}
		
		for(i=0;i<basis::tri(log2p)->gpx(); ++i) {
			
			for(int m = 0; m < NV; ++m)
				temp(m) = u1d(m)(i);
				
			calculate_preconditioner(temp, P);
			
			temp = 0.0;
			for(int m = 0; m < NV; ++m)
				for(int n = 0; n < NV; ++n)							
					temp(m) += P(m,n)*res1d(n)(i);
			
			for(n=0;n<NV;++n)
				temp1d(n)(i) -= temp(n);
		}
		
		/* integrate right hand side: phi*P*res */
		for(n=0;n<NV;++n)
			basis::tri(log2p)->intgrt1d(&lf(n)(0),&temp1d(n)(0));
		
		/* invert 1d mass matrix */
		for(n=0;n<NV;++n) {
			PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&lf(n)(2),basis::tri(log2p)->sm(),info);
			for(m=0;m<basis::tri(log2p)->sm();++m) 
				gbl->res.s(sind,m,n) = -lf(n)(2+m);
		}
	}
	
	if (basis::tri(log2p)->im() <= 0) return;
	
	for(tind = 0; tind < ntri; ++tind) {
		
		for (i=0; i<3; ++i) {
			indx = tri(tind).pnt(i);
			for(n=0; n<NV; ++n){
				ucoef(n)(i) = gbl->ug0.v(indx,n);
				rcoef(n)(i) = gbl->res.v(indx,n);
				tcoef(n)(i) = gbl->res.v(indx,n);
			}
		}
		
		/* SIDES */
		cnt = 3;
		for(i=0;i<3;++i) {
			sind = tri(tind).seg(i);
			sign = tri(tind).sgn(i);
			msgn = 1;
			for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
				for(n=0; n<NV; ++n){
					ucoef(n)(cnt) = msgn*gbl->ug0.s(sind,m,n);
					rcoef(n)(cnt) = msgn*gbl->res.s(sind,m,n);
					tcoef(n)(cnt) = msgn*gbl->res.s(sind,m,n);
				}
				msgn *= sign;
				++cnt;
			}
		}
		
		/* INTERIORS */    
		if (basis::tri(log2p)->im() > 0) {    
			indx = 0;
			for(m = 1; m < basis::tri(log2p)->sm(); ++m) {
				for(k = 0; k < basis::tri(log2p)->sm()-m; ++k) {
					for(n=0; n<NV; ++n){
						ucoef(n)(cnt) = gbl->ug0.i(tind,indx,n);
						rcoef(n)(cnt) = gbl->res.i(tind,indx,n);
						tcoef(n)(cnt) = 0.0;
					}
					++cnt; ++indx;
				}
				indx += sm0 -basis::tri(log2p)->sm();
			}
		}
		
		
		for(n=0;n<NV;++n){			
			basis::tri(log2p)->proj(&ucoef(n)(0),&u2d(n)(0,0),MXGP);
			basis::tri(log2p)->proj(&rcoef(n)(0),&res2d(n)(0,0),MXGP);
			basis::tri(log2p)->proj_bdry(&tcoef(n)(0),&temp2d(n)(0,0),MXGP);
		}
		
		for (i=0; i < basis::tri(log2p)->gpx(); ++i ) {
			for (j=0; j < basis::tri(log2p)->gpn(); ++j ) {

				for(int m = 0; m < NV; ++m)
					temp(m) = u2d(m)(i,j);
				
				calculate_preconditioner(temp, P);
				
				temp = 0.0;
				for(int m = 0; m < NV; ++m)
					for(int n = 0; n < NV; ++n)							
						temp(m) += P(m,n)*res2d(n)(i,j);
				
				for(n=0;n<NV;++n)
					temp2d(n)(i,j) -= temp(n);
			}
		}
		
		for(n=0;n<NV;++n) {
			basis::tri(log2p)->intgrt(&lf(n)(0),&temp2d(n)(0,0),MXGP);
			DPBTRS(uplo,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),1,(double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,&lf(n)(basis::tri(log2p)->bm()),basis::tri(log2p)->im(),info);
			for(i=0;i<basis::tri(log2p)->im();++i)
				gbl->res.i(tind,i,n) = -lf(n)(basis::tri(log2p)->bm()+i);
		}
	}
	
	return;
}

void tri_hp_cns::calculate_preconditioner(Array<double,1> upv, Array<double,2> &P){
	
	FLT pr = upv(0);
	FLT u = upv(1);
	FLT v = upv(2);
	FLT rt = upv(3);
	FLT rho = pr/rt;
	FLT ke = 0.5*(u*u+v*v);
	FLT gm1 = gbl->gamma-1.0;
	
	/* Preconditioner */
	P = ke*gm1,          -u*gm1,     -v*gm1,      gm1,
		-u/rho,          1.0/rho,    0.0,         0.0,
		-v/rho,          0.0,        1.0/rho,     0.0,
		(gm1*ke-rt)/rho, -u*gm1/rho, -v*gm1/rho, gm1/rho;
	
	return;
}

void tri_hp_cns::minvrt() {
	
	tri_hp::minvrt();
	
	project_new_variables();
	
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
