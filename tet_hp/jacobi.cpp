/*
 *  jacobi.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 1/10/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "hp_boundary.h"
#include <myblas.h>

void tet_hp::jacobi_relaxation() {
	int em = basis::tet(log2p).em;
	int fm = basis::tet(log2p).fm;
	int im = basis::tet(log2p).im;

	/* Jacobi's relaxation */
	FLT omega = 0.5;
	hp_gbl->res.v(Range(0,npnt-1),Range::all()) = omega*hp_gbl->res.v(Range(0,npnt-1),Range::all())/hp_gbl->jacob_diag.v(Range(0,npnt-1),Range::all());
	if(em) hp_gbl->res.e(Range(0,nseg-1),Range(0,em-1),Range::all()) = omega*hp_gbl->res.e(Range(0,nseg-1),Range(0,em-1),Range::all())/hp_gbl->jacob_diag.e(Range(0,nseg-1),Range(0,em-1),Range::all());
	if(fm) hp_gbl->res.f(Range(0,ntri-1),Range(0,fm-1),Range::all()) = omega*hp_gbl->res.f(Range(0,ntri-1),Range(0,fm-1),Range::all())/hp_gbl->jacob_diag.f(Range(0,ntri-1),Range(0,fm-1),Range::all());
	if(im) hp_gbl->res.i(Range(0,ntet-1),Range(0,im-1),Range::all()) = omega*hp_gbl->res.i(Range(0,ntet-1),Range(0,im-1),Range::all())/hp_gbl->jacob_diag.i(Range(0,ntet-1),Range(0,im-1),Range::all());
	
	all_dirichlet();
	
	
	return;
}

void tet_hp::jacobian_diagonal() {
	
	Array<TinyVector<FLT,MXTM>,1> Rbar(NV),lf_re(NV),lf_im(NV);
	Array<FLT,1> dw(NV);
#ifdef BZ_DEBUG
	const FLT eps_r = 0.0e-6, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
	const FLT eps_r = 1.0e-6, eps_a = 1.0e-10;  /*<< constants for accurate numerical determination of jacobians */
#endif
	
	int em = basis::tet(log2p).em;
	int fm = basis::tet(log2p).fm;
	int im = basis::tet(log2p).im;
	int tm = basis::tet(log2p).tm;
	
	int dof = tm*NV;
	Array<FLT,1> Kdiag(dof);
	
	hp_gbl->jacob_diag.v = 0.0;
	hp_gbl->jacob_diag.e = 0.0;
	hp_gbl->jacob_diag.f = 0.0;
	hp_gbl->jacob_diag.i = 0.0;
	
	
	/* Calculate diagonal stiffness matrix terms */
	for(int tind = 0; tind < ntet; ++tind) {	
		
		ugtouht(tind);
		
		dw = 0.0;
		for(int i=0;i<4;++i)
			for(int n=0;n<NV;++n)
				dw = dw + fabs(uht(n)(i));
		
		dw = dw*eps_r;
		dw = dw+eps_a;
		
		element_rsdl(tind,0,uht,lf_re,lf_im);
		
		for(int i=0;i<tm;++i) 
			for(int n=0;n<NV;++n) 
				Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);		
		
		int ind = 0;
		for(int i = 0; i < tm; ++i){
			for(int n = 0; n < NV; ++n){
				uht(n)(i) += dw(n);
				
				element_rsdl(tind,0,uht,lf_re,lf_im);
				
				Kdiag(ind++) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw(n);
				
				uht(n)(i) -= dw(n);
			}
		}
		
		ind = 0;	
		for (int i=0; i < 4; ++i) 
			for(int n=0;n<NV;++n)
				hp_gbl->jacob_diag.v(tet(tind).pnt(i),n) += Kdiag(ind++);
		
		for (int i=0; i < 6; ++i) 
			for(int m=0;m<em;++m)
				for(int n=0;n<NV;++n)
					hp_gbl->jacob_diag.e(tet(tind).seg(i),m,n) += Kdiag(ind++);
		
		for (int i=0; i < 4; ++i) 
			for(int m=0;m<fm;++m)
				for(int n=0;n<NV;++n)
					hp_gbl->jacob_diag.f(tet(tind).tri(i),m,n) += Kdiag(ind++);
		
		for(int m=0;m<im;++m)
			for(int n=0;n<NV;++n)
				hp_gbl->jacob_diag.i(tind,m,n) += Kdiag(ind++);			
		
		
	}
	
	return;
}
