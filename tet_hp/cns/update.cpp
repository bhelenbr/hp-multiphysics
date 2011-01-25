/*
 *  update.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 1/10/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */

#include "tet_hp_cns.h"
#include "../hp_boundary.h"
#include <myblas.h>

void tet_hp_cns::project_new_variables(){
	int info,last_phase, mp_phase;
	char uplo[] = "U";
	Array<double,1> lcl(NV), lclug(NV),lclres(NV);
	Array<TinyVector<double,MXGP>,2> P(NV,NV);
	Array<TinyMatrix<double,MXGP,MXGP>,2> P2d(NV,NV);
	Array<TinyVector<double,MXGP>,1> u1d(NV),res1d(NV),temp1d(NV);
	Array<TinyMatrix<double,MXGP,MXGP>,1> u2d(NV),res2d(NV),temp2d(NV);
	Array<TinyVector<double,MXTM>,1> ucoef(NV),rcoef(NV),tcoef(NV);
	
	vefi res_temp;
	res_temp.v.resize(maxvst,NV);
	res_temp.e.resize(maxvst,em0,NV);
	res_temp.f.resize(maxvst,fm0,NV);
	res_temp.i.resize(maxvst,im0,NV);
	res_temp.v = gbl->res.v;
	res_temp.e = gbl->res.e;
	res_temp.f = gbl->res.f;
	res_temp.i = gbl->res.i;
	
	/* LOOP THROUGH VERTICES */
	for(int i=0;i<npnt;++i){
		
		for(int n = 0; n < NV; ++n){
			lclug(n) = ug.v(i,n);
			lclres(n) = gbl->res.v(i,n);
		}
		
		switch_variables(lclug,lclres);
		
		for(int j=0;j<NV;++j){
			FLT lcl0 = lclres(j);
			for(int k=0;k<j;++k){
				lcl0 -= gbl->vpreconditioner(i,j,k)*lclres(k);
			}
			lclres(j) = lcl0/gbl->vpreconditioner(i,j,j);
		}
		
		for(int n = 0; n < NV; ++n)
			gbl->res.v(i,n) = lclres(n);
		
	}
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		pc0load(mp_phase,gbl->res.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= pc0wait_rcv(mp_phase,gbl->res.v.data());
	}
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<nfbd;++i)
		hp_fbdry(i)->vdirichlet();
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet3d();	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet3d();
	
	if (!basis::tet(log2p).em) return;
	
//	/* LOOP THROUGH SIDES */    
//    for(int sind=0;sind<nseg;++sind) {
//		
//		/* project linears */
//		for(int n=0;n<NV;++n)
//			basis::tri(log2p)->proj1d(gbl->res.v(seg(sind).pnt(0),n),gbl->res.v(seg(sind).pnt(1),n),&temp1d(n)(0));
//		
//		for(int m=0;m<NV;++m) 
//			for(int n=0;n<NV;++n) 
//				basis::tri(log2p)->proj1d(gbl->vpreconditioner(seg(sind).pnt(0),m,n),gbl->vpreconditioner(seg(sind).pnt(1),m,n),&P(m,n)(0));
//		
//		/* take global coefficients and put into local vector */
//		for(int n=0;n<NV;++n) {
//			for (int m=0; m<2; ++m) {
//				ucoef(n)(m) = ug.v(seg(sind).pnt(m),n);
//				rcoef(n)(m) = gbl->res_temp.v(seg(sind).pnt(m),n);
//			}
//		}
//		
//		for(int n=0;n<NV;++n) {
//			for(int m=0;m<basis::tri(log2p)->sm();++m) {
//				ucoef(n)(m+2) = ug.s(sind,m,n);
//				rcoef(n)(m+2) = gbl->res_temp.s(sind,m,n);
//			}
//		}
//		
//		/* project sides */
//		for(int n=0;n<NV;++n) {
//			basis::tri(log2p)->proj1d(&ucoef(n)(0),&u1d(n)(0));
//			basis::tri(log2p)->proj1d(&rcoef(n)(0),&res1d(n)(0));
//		}
//		
//		for(int i=0;i<basis::tri(log2p)->gpx(); ++i) {
//			
//			for(int m = 0; m < NV; ++m){
//				lclug(m) = u1d(m)(i);
//				lclres(m) = res1d(m)(i);
//			}
//			
//			switch_variables(lclug,lclres);
//			
//			for(int j=0;j<NV;++j){
//				FLT lcl0 = lclres(j);
//				for(int k=0;k<j;++k){
//					lcl0 -= P(j,k)(i)*lclres(k);
//				}
//				lclres(j) = lcl0/P(j,j)(i);
//			}
//			
//			for(int m = 0; m < NV; ++m)
//				temp1d(m)(i) -= lclres(m);	
//			
//		}
//		
//		/* integrate right hand side: phi*P*res */
//		for(int n=0;n<NV;++n)
//			basis::tri(log2p)->intgrt1d(&rcoef(n)(0),&temp1d(n)(0));
//		
//		/* invert 1d mass matrix */
//		for(int n=0;n<NV;++n) {
//			PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&rcoef(n)(2),basis::tri(log2p)->sm(),info);
//			for(int m=0;m<basis::tri(log2p)->sm();++m) 
//				gbl->res.s(sind,m,n) = -rcoef(n)(m+2);
//			
//		}
//		
//	}
//	
//	for(int mode = 0; mode < basis::tri(log2p)->sm(); ++mode) {
//		
//		sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
//		smsgpass(boundary::all,0,boundary::symmetric);
//		sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
//		
//		/* APPLY DIRCHLET B.C.S TO MODE */
//		for(int i=0;i<nebd;++i)
//			hp_ebdry(i)->sdirichlet(mode);
//	}
//	
//	if (!basis::tri(log2p)->im()) return;
//	
//	for(int tind = 0; tind < ntri; ++tind) {
//		
//		for (int i=0; i<3; ++i) {
//			int vrtx = tri(tind).pnt(i);
//			for(int n=0; n<NV; ++n){
//				ucoef(n)(i) = ug.v(vrtx,n);
//				rcoef(n)(i) = gbl->res_temp.v(vrtx,n);
//				tcoef(n)(i) = gbl->res.v(vrtx,n);
//			}
//		}
//		
//		/* SIDES */
//		int cnt = 3;
//		for(int i=0;i<3;++i) {
//			int sind = tri(tind).seg(i);
//			int sign = tri(tind).sgn(i);
//			int msgn = 1;
//			for (int m = 0; m < basis::tri(log2p)->sm(); ++m) {
//				for(int n=0; n<NV; ++n){
//					ucoef(n)(cnt) = msgn*ug.s(sind,m,n);
//					rcoef(n)(cnt) = msgn*gbl->res_temp.s(sind,m,n);
//					tcoef(n)(cnt) = msgn*gbl->res.s(sind,m,n);
//				}
//				msgn *= sign;
//				++cnt;
//			}
//		}
//		
//		/* INTERIORS */    
//		if (basis::tri(log2p)->im() > 0) {    
//			int indx = 0;
//			for(int m = 1; m < basis::tri(log2p)->sm(); ++m) {
//				for(int k = 0; k < basis::tri(log2p)->sm()-m; ++k) {
//					for(int n=0; n<NV; ++n){
//						ucoef(n)(cnt) = ug.i(tind,indx,n);
//						rcoef(n)(cnt) = gbl->res_temp.i(tind,indx,n);
//						tcoef(n)(cnt) = 0.0;
//					}
//					++cnt; ++indx;
//				}
//				indx += sm0 -basis::tri(log2p)->sm();
//			}
//		}
//		
//		
//		for(int n=0;n<NV;++n){			
//			basis::tri(log2p)->proj(&ucoef(n)(0),&u2d(n)(0,0),MXGP);
//			basis::tri(log2p)->proj(&rcoef(n)(0),&res2d(n)(0,0),MXGP);
//			basis::tri(log2p)->proj_bdry(&tcoef(n)(0),&temp2d(n)(0,0),MXGP);
//		}
//		
//		for(int m=0;m<NV;++m)		
//			for(int n=0;n<NV;++n)
//				basis::tri(log2p)->proj(gbl->vpreconditioner(tri(tind).pnt(0),m,n),gbl->vpreconditioner(tri(tind).pnt(1),m,n),gbl->vpreconditioner(tri(tind).pnt(2),m,n),&P2d(m,n)(0,0),MXGP);
//		
//		
//		for (int i=0; i < basis::tri(log2p)->gpx(); ++i ) {
//			for (int j=0; j < basis::tri(log2p)->gpn(); ++j ) {
//				
//				for(int m = 0; m < NV; ++m){
//					lclug(m) = u2d(m)(i,j);
//					lclres(m) = res2d(m)(i,j);
//				}
//				
//				switch_variables(lclug,lclres);
//				
//				for(int m=0;m<NV;++m){
//					FLT lcl0 = lclres(m);
//					for(int n=0;n<m;++n){
//						lcl0 -= P2d(m,n)(i,j)*lclres(n);
//					}
//					lclres(m) = lcl0/P2d(m,m)(i,j);
//				}
//				
//				for(int m = 0; m < NV; ++m)
//					temp2d(m)(i,j) -= lclres(m);	
//				
//				
//			}
//		}
//		
//		for(int n=0;n<NV;++n) {
//			basis::tri(log2p)->intgrt(&rcoef(n)(0),&temp2d(n)(0,0),MXGP);
//			DPBTRS(uplo,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),1,(double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,&rcoef(n)(basis::tri(log2p)->bm()),basis::tri(log2p)->im(),info);
//			for(int i=0;i<basis::tri(log2p)->im();++i)
//				gbl->res.i(tind,i,n) = -rcoef(n)(basis::tri(log2p)->bm()+i);
//		}
//	}
	
	return;
}

void tet_hp_cns::switch_variables(Array<double,1> pvu, Array<double,1> &a){
	
	Array<double,2> dpdc(NV,NV);
	Array<double,1> temp(NV);
	double gm1 = gbl->gamma-1.0;
	
	double pr = pvu(0),u = pvu(1),v = pvu(2),w = pvu(3), rt = pvu(4);
	double rho = pr/rt;
	double ke = 0.5*(u*u+v*v+w*w);
	
	/* jacobian derivative of primitive wrt conservative */
	dpdc = ke*gm1,          -u*gm1,     -v*gm1,     -w*gm1,     gm1,
		   -u/rho,          1.0/rho,    0.0,        0.0,        0.0,
		   -v/rho,          0.0,        1.0/rho,    0.0,        0.0,
		   -w/rho,          0.0,        0.0,        1.0/rho,    0.0,
		   (gm1*ke-rt)/rho, -u*gm1/rho, -v*gm1/rho, -w*gm1/rho, gm1/rho;
	
	temp = 0.0;
	for(int i = 0; i < NV; ++i)
		for(int j = 0; j < NV; ++j)
			temp(i) += dpdc(i,j)*a(j);
	
	a = temp;
	
	return;
}




void tet_hp_cns::update() {
		
	/* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
	gbl->ug0.v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all());
	if (basis::tet(log2p).em) {
		gbl->ug0.e(Range(0,nseg-1),Range(0,em0-1),Range::all()) = ug.e(Range(0,nseg-1),Range::all(),Range::all());
		if (basis::tet(log2p).fm) {
			gbl->ug0.f(Range(0,ntri-1),Range(0,fm0-1),Range::all()) = ug.f(Range(0,ntri-1),Range::all(),Range::all());
			if (basis::tet(log2p).im) {
				gbl->ug0.i(Range(0,ntet-1),Range(0,im0-1),Range::all()) = ug.i(Range(0,ntet-1),Range::all(),Range::all());
			}
		}
	}
	
	for(int i=0;i<nfbd;++i)
		hp_fbdry(i)->update(-1);
	
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->update(-1);
	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->update(-1);
	
	helper->update(-1);
	
	for (int stage = 0; stage < gbl->nstage; ++stage) {
		
		tet_hp::rsdl(stage);	
		
		for(int i=0;i<nfbd;++i)
			hp_fbdry(i)->vdirichlet();
		for(int i=0;i<nfbd;++i)
			hp_fbdry(i)->edirichlet();
		for(int i=0;i<nfbd;++i)
			hp_fbdry(i)->fdirichlet();
		
#ifdef JACOBI
		tet_hp::jacobi_relaxation();
#else
		tet_hp::minvrt();
#endif
		
		for(int i=0;i<nfbd;++i)
			hp_fbdry(i)->modify_boundary_residual();		
		
		project_new_variables();
		
		FLT cflalpha = gbl->alpha(stage)*gbl->cfl(log2p);
		
		ug.v(Range(0,npnt-1),Range::all()) = gbl->ug0.v(Range(0,npnt-1),Range::all()) -cflalpha*gbl->res.v(Range(0,npnt-1),Range::all());
		
		if (basis::tet(log2p).em > 0) {
			ug.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = gbl->ug0.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) -cflalpha*gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
			
			if (basis::tet(log2p).fm > 0) {
				
				for(int i=0;i<ntri;++i) {
					int indx = 0;
					int indx1 = 0;
					for(int m=1;m<=basis::tet(log2p).em;++m) {
						for(int k=1;k<=basis::tet(log2p).em-m;++k) {
							for(int n=0;n<NV;++n) {
								ug.f(i,indx1,n) =  gbl->ug0.f(i,indx1,n) -cflalpha*gbl->res.f(i,indx,n);
							}
							++indx; ++indx1;
						}
						indx1 += em0 -basis::tet(log2p).em;
					}
				}
				if (basis::tet(log2p).im > 0) {
					
					for(int i=0;i<ntet;++i) {
						int indx = 0;
						int indx1 = 0;
						for(int m=1;m<=basis::tet(log2p).em-1;++m) {
							for(int j=1;j<=basis::tet(log2p).em-m;++j) {
								for(int k=1;k<=basis::tet(log2p).em-m-j;++k) {
									for(int n=0;n<NV;++n) {
										ug.i(i,indx1,n) =  gbl->ug0.i(i,indx1,n) -cflalpha*gbl->res.i(i,indx,n);
									}
									++indx; ++indx1;
								}
								indx1 += em0 -basis::tet(log2p).em;
							}
						}
					}
				}
			}
		}
		
		helper->update(stage);
		
		for(int i=0;i<nfbd;++i) {
			hp_fbdry(i)->update(stage);
		}
		
		for(int i=0;i<nebd;++i) {
			hp_ebdry(i)->update(stage);
		}
		
		for(int i=0;i<nvbd;++i) {
			hp_vbdry(i)->update(stage);
		}
	
	}
	return;
}


