/*
 *  minvrt.cpp
 *  tri_hp
 *
 *  Created by michael brazell on 1/4/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */


#include "tri_hp_cd.h"
#include "../hp_boundary.h"
#include <myblas.h>


// Jacobi Relaxation

void tri_hp_cd::minvrt() {
		
	Array<TinyVector<FLT,MXTM>,1> Rbar(NV),lf_re(NV),lf_im(NV);
	Array<FLT,1> dw(NV);
	
	int sm = basis::tri(log2p)->sm();
	int im = basis::tri(log2p)->im();
	int tm = basis::tri(log2p)->tm();
	
	int dof = tm*NV;
	Array<FLT,1> Kdiag(dof);
	
	gbl->stiff_diag.v = 0.0;
	gbl->stiff_diag.s = 0.0;
	gbl->stiff_diag.i = 0.0;
	
	
	/* Calculate diagonal stiffness matrix terms */
	for(int tind = 0; tind < ntri; ++tind) {	
		
		ugtouht(tind);
		
		dw = 0.0;
		for(int i=0;i<3;++i)
			for(int n=0;n<NV;++n)
				dw(n) = dw(n) + fabs(uht(n)(i));

		dw *= eps_r;
		dw = dw+eps_a;
		
		element_rsdl(tind,0,uht,lf_re,lf_im);
		
		for(int i=0;i<basis::tri(log2p)->tm();++i) 
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
		for (int i=0; i < 3; ++i) 
			for(int n=0;n<NV;++n)
				gbl->stiff_diag.v(tri(tind).pnt(i),n) += Kdiag(ind++);
		
		for (int i=0; i < 3; ++i) 
			for(int m=0;m<sm;++m)
				for(int n=0;n<NV;++n)
					gbl->stiff_diag.s(tri(tind).seg(i),m,n) += Kdiag(ind++);
		
		for(int m=0;m<im;++m)
			for(int n=0;n<NV;++n)
				gbl->stiff_diag.i(tind,m,n) += Kdiag(ind++);			
		
		
	}
	
	/* Preconditioner */
	/* vprcn is predivided by vdiag so multiply by vdiag to cancel out */
	hp_gbl->res.v(Range(0,npnt-1),Range::all()) *= hp_gbl->vprcn(Range(0,npnt-1),Range::all())*basis::tri(log2p)->vdiag();
	
	for(int mode=0;mode<sm;++mode)
		hp_gbl->res.s(Range(0,nseg-1),mode,Range::all()) *= hp_gbl->sprcn(Range(0,nseg-1),Range::all());
	
	for(int mode=0;mode<im;++mode)
		hp_gbl->res.i(Range(0,ntri-1),mode,Range::all()) /= hp_gbl->tprcn(Range(0,ntri-1),Range::all());
	
	/* Jacobi's relaxation */
	hp_gbl->res.v(Range(0,npnt-1),Range::all()) /= gbl->stiff_diag.v(Range(0,npnt-1),Range::all());
	if(sm) hp_gbl->res.s(Range(0,nseg-1),Range(0,sm-1),Range::all()) /= gbl->stiff_diag.s(Range(0,nseg-1),Range(0,sm-1),Range::all());
	if(im) hp_gbl->res.i(Range(0,ntri-1),Range(0,im-1),Range::all()) /= gbl->stiff_diag.i(Range(0,ntri-1),Range(0,im-1),Range::all());
	
	/* Apply BC's */
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet();
	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet();
	
	for(int mode=0;mode<sm;++mode)
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->sdirichlet(mode);
	
	return;
}
