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





void tet_hp_cns::update() {
	
	//tet_hp::update();return;
		
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
		
		
//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->vdirichlet();		
//		for(int i=0;i<nebd;++i)
//			hp_ebdry(i)->vdirichlet3d();	
//		for(int i=0;i<nvbd;++i)
//			hp_vbdry(i)->vdirichlet3d();		
//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->edirichlet();		
//		for (int i=0;i<nebd;++i) 
//			hp_ebdry(i)->edirichlet3d();
		

		minvrt();


//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->vdirichlet();		
//		for(int i=0;i<nebd;++i)
//			hp_ebdry(i)->vdirichlet3d();	
//		for(int i=0;i<nvbd;++i)
//			hp_vbdry(i)->vdirichlet3d();		
//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->edirichlet();		
//		for (int i=0;i<nebd;++i) 
//			hp_ebdry(i)->edirichlet3d();
		
		
//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->modify_boundary_residual();		
		
		//project_new_variables();
		
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


void tet_hp_cns::minvrt() {
	int i,j,k,n,tind,msgn,sgn,sind,v0;
	Array<FLT,2> spokemass;
	int last_phase, mp_phase;
	
	Array<double,1> lcl(NV), lclug(NV),lclres(NV),uavg(NV);
	Array<TinyVector<double,MXGP>,2> P(NV,NV);
	Array<TinyVector<double,MXGP>,1> u1d(NV),res1d(NV),temp1d(NV);
	Array<TinyVector<double,MXTM>,1> ucoef(NV),rcoef(NV),tcoef(NV);
	
	if (basis::tet(log2p).p > 2) {
		*gbl->log << "cns minvrt only works for p = 1 and 2" << endl;
		exit(4);
	}
	
	/* LOOP THROUGH EDGES */
	if (basis::tet(log2p).em > 0) {
		for(int eind = 0; eind<nseg;++eind) {
			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */         
			for (k=0; k <basis::tet(log2p).em; ++k) {
				for (i=0; i<2; ++i) {
					v0 = seg(eind).pnt(i);
					for(n=0;n<NV;++n)
						gbl->res.v(v0,n) -= basis::tet(log2p).sfmv(i,k)*gbl->res.e(eind,k,n);
				}
			}
		}		
	}
	
	gbl->res.v(Range(0,npnt-1),Range::all()) *= gbl->vprcn(Range(0,npnt-1),Range::all())*basis::tet(log2p).vdiag;

	/* LOOP THROUGH VERTICES */
	for(int i=0;i<npnt;++i){
		
		for(int n = 0; n < NV; ++n)
			lclres(n) = gbl->res.v(i,n);

		
		if(gbl->preconditioner == 0 || gbl->preconditioner == 1) {
			for(int n = 0; n < NV; ++n)
				lclug(n) = ug.v(i,n);

			switch_variables(lclug,lclres);

			for(int j=0;j<NV;++j){
				FLT lcl0 = lclres(j);
				for(int k=0;k<j;++k){
					lcl0 -= gbl->vpreconditioner(i,j,k)*lclres(k);
				}
				lclres(j) = lcl0/gbl->vpreconditioner(i,j,j);
			}
		}
		else {		
			int info,ipiv[NV];
			Array<double,2> P(NV,NV);
			
			for(int j=0;j<NV;++j)
				for(int k=0;k<NV;++k)
					P(j,k) = gbl->vpreconditioner(i,j,k);
				
			GETRF(NV, NV, P.data(), NV, ipiv, info);

			if (info != 0) {
				*gbl->log << "DGETRF FAILED FOR CNS MINVRT" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			char trans[] = "T";
			GETRS(trans,NV,1,P.data(),NV,ipiv,lclres.data(),NV,info);
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
	for(i=0;i<nfbd;++i)
		hp_fbdry(i)->vdirichlet();
	
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet3d();        
	
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet3d();
	
	if(basis::tet(log2p).em == 0) return;
	
	/* LOOP THROUGH SIDES */    
	for(int sind=0;sind<nseg;++sind) {
		
		for(int n = 0; n < NV; ++n)
			lclres(n) = gbl->res.e(sind,0,n);

		Array<FLT,2> P(NV,NV);
		for(int j=0;j<NV;++j){
			for(int k=0;k<NV;++k){
				P(j,k) = gbl->epreconditioner(sind,j,k);
				//P(j,k) = 0.5*(gbl->vpreconditioner(seg(sind).pnt(0),j,k)+gbl->vpreconditioner(seg(sind).pnt(1),j,k));
			}
		}

		if(gbl->preconditioner == 0 || gbl->preconditioner == 1) {
			for(int n = 0; n < NV; ++n)
				uavg(n) = 0.5*(ug.v(seg(sind).pnt(0),n)+ug.v(seg(sind).pnt(1),n));
				
			switch_variables(uavg,lclres);
			
			for(int j=0;j<NV;++j){
				FLT lcl0 = lclres(j);
				for(int k=0;k<j;++k){
					lcl0 -= P(j,k)*lclres(k);
				}
				lclres(j) = lcl0/P(j,j);
			}
		}
		else {
			int info,ipiv[NV];
			
			GETRF(NV, NV, P.data(), NV, ipiv, info);
			
			if (info != 0) {
				*gbl->log << "DGETRF FAILED FOR CNS MINVRT EDGE" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			char trans[] = "T";
			GETRS(trans,NV,1,P.data(),NV,ipiv,lclres.data(),NV,info);
		}
		
		for(int n = 0; n < NV; ++n)
			gbl->res.e(sind,0,n) = lclres(n);
		
	}
	
	/* REMOVE VERTEX CONTRIBUTION FROM SIDE MODES */
	/* SOLVE FOR SIDE MODES */
	/* PART 1 REMOVE VERTEX CONTRIBUTIONS */
	for(tind=0;tind<ntet;++tind) {         
		for(i=0;i<4;++i) {
			v0 = tet(tind).pnt(i);
			for(n=0;n<NV;++n)
				uht(n)(i) = gbl->res.v(v0,n)*gbl->iprcn(tind,n);
		}
		/* edges */
		for(i=0;i<6;++i) {
			sind = tet(tind).seg(i);
			sgn  = tet(tind).sgn(i);
			for(j=0;j<4;++j) {
				msgn = 1;
				for(k=0;k<basis::tet(log2p).em;++k) {
					for(n=0;n<NV;++n)
						gbl->res.e(sind,k,n) -= msgn*basis::tet(log2p).vfms(j,4+k+i*basis::tet(log2p).em)*uht(n)(j);
					msgn *= sgn;
				}
			}
		}				
	}
	
	
	basis::tet(log2p).ediag(0) = 100.0;//for fast convergence 
	//basis::tet(log2p).ediag(0) = 48.0; //for accuracy mass lumped edge modes
	gbl->res.e(Range(0,nseg-1),0,Range::all()) *= gbl->eprcn(Range(0,nseg-1),Range::all())*basis::tet(log2p).ediag(0);
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		sc0load(mp_phase,gbl->res.e.data(),0,0,gbl->res.e.extent(secondDim));
		smsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= sc0wait_rcv(mp_phase,gbl->res.e.data(),0,0,gbl->res.e.extent(secondDim));
	}
	
	/* APPLY DIRCHLET B.C.S TO MODE */
	for(int i=0;i<nfbd;++i)
		hp_fbdry(i)->edirichlet();
	
	for (int i=0;i<nebd;++i) 
		hp_ebdry(i)->edirichlet3d();	
	
	return;
}

void tet_hp_cns::switch_variables(Array<double,1> pvu, Array<double,1> &a){
	
	Array<double,2> dpdc(NV,NV);
	Array<double,1> temp(NV);
	double gm1 = gbl->gamma-1.0;
	
	double pr = pvu(0);
	double u = pvu(1),v = pvu(2),w = pvu(3), rt = pvu(4);
	double rho = (pr+gbl->atm_pressure)/rt;
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


//void tet_hp_cns::project_new_variables(){
//	int info,last_phase, mp_phase;
//	char uplo[] = "U";
//	Array<double,1> lcl(NV), lclug(NV),lclres(NV);
//	Array<TinyVector<double,MXGP>,2> P(NV,NV);
//	Array<TinyMatrix<double,MXGP,MXGP>,2> P2d(NV,NV);
//	Array<TinyVector<double,MXGP>,1> u1d(NV),res1d(NV),temp1d(NV);
//	Array<TinyMatrix<double,MXGP,MXGP>,1> u2d(NV),res2d(NV),temp2d(NV);
//	Array<TinyVector<double,MXTM>,1> ucoef(NV),rcoef(NV),tcoef(NV);
//	
//	vefi res_temp(gbl->res);
//	
//	/* LOOP THROUGH VERTICES */
//	for(int i=0;i<npnt;++i){
//		
//		for(int n = 0; n < NV; ++n){
//			lclug(n) = ug.v(i,n);
//			lclres(n) = gbl->res.v(i,n);
//		}
//		
//		switch_variables(lclug,lclres);
//		
//		for(int j=0;j<NV;++j){
//			FLT lcl0 = lclres(j);
//			for(int k=0;k<j;++k){
//				lcl0 -= gbl->vpreconditioner(i,j,k)*lclres(k);
//			}
//			lclres(j) = lcl0/gbl->vpreconditioner(i,j,j);
//		}
//		
//		for(int n = 0; n < NV; ++n)
//			gbl->res.v(i,n) = lclres(n);
//		
//	}
//	
//	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
//		pc0load(mp_phase,gbl->res.v.data());
//		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
//		last_phase = true;
//		last_phase &= pc0wait_rcv(mp_phase,gbl->res.v.data());
//	}
//	
//	/* APPLY VERTEX DIRICHLET B.C.'S */
//	for(int i=0;i<nfbd;++i)
//		hp_fbdry(i)->vdirichlet();
//	for(int i=0;i<nebd;++i)
//		hp_ebdry(i)->vdirichlet3d();	
//	for(int i=0;i<nvbd;++i)
//		hp_vbdry(i)->vdirichlet3d();
//	
//	if (!basis::tet(log2p).em) return;
//	
//	/* LOOP THROUGH SIDES */    
//    for(int sind=0;sind<nseg;++sind) {
//		
//		/* project linears */
//		for(int n=0;n<NV;++n)
//			basis::tet(log2p).proj1d(gbl->res.v(seg(sind).pnt(0),n),gbl->res.v(seg(sind).pnt(1),n),&temp1d(n)(0));
//		
//		for(int m=0;m<NV;++m) 
//			for(int n=0;n<NV;++n) 
//				basis::tet(log2p).proj1d(gbl->vpreconditioner(seg(sind).pnt(0),m,n),gbl->vpreconditioner(seg(sind).pnt(1),m,n),&P(m,n)(0));
//		
//		/* take global coefficients and put into local vector */
//		for(int n=0;n<NV;++n) {
//			for (int m=0; m<2; ++m) {
//				ucoef(n)(m) = ug.v(seg(sind).pnt(m),n);
//				rcoef(n)(m) = res_temp.v(seg(sind).pnt(m),n);
//			}
//		}
//		
//		for(int n=0;n<NV;++n) {
//			for(int m=0;m<basis::tet(log2p).em;++m) {
//				ucoef(n)(m+2) = ug.e(sind,m,n);
//				rcoef(n)(m+2) = res_temp.e(sind,m,n);
//			}
//		}
//		
//		/* project sides */
//		for(int n=0;n<NV;++n) {
//			basis::tet(log2p).proj1d(&ucoef(n)(0),&u1d(n)(0));
//			basis::tet(log2p).proj1d(&rcoef(n)(0),&res1d(n)(0));
//		}
//		
//		for(int i=0;i<basis::tet(log2p).gpx; ++i) {
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
//			basis::tet(log2p).intgrt1d(&rcoef(n)(0),&temp1d(n)(0));
//		
//		/* invert 1d mass matrix */
//		for(int n=0;n<NV;++n) {
//			for(int m=0;m<basis::tet(log2p).em;++m) 
//				gbl->res.e(sind,m,n) = -rcoef(n)(m+2)*basis::tet(log2p).diag1d(m);			
//		}
//		
//	}
//	
//	for(int mode = 0; mode < basis::tet(log2p).em; ++mode) {
//		for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
//			sc0load(mp_phase,gbl->res.e.data(),0,0,gbl->res.e.extent(secondDim));
//			smsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
//			last_phase = true;
//			last_phase &= sc0wait_rcv(mp_phase,gbl->res.e.data(),0,0,gbl->res.e.extent(secondDim));
//		}
//		/* APPLY DIRCHLET B.C.S TO MODE */
//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->edirichlet();
//		
//		for (int i=0;i<nebd;++i) 
//			hp_ebdry(i)->edirichlet3d();
//	}
//	
//	
//	return;
//}



