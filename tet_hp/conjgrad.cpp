/*
 *  conjgrad.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 1/10/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */


#include "tet_hp.h"
#include "hp_boundary.h"
#include <myblas.h>

void tet_hp::conjugate_gradient() {
	
	rsdl();
	
	FLT omega = 1.0;
	FLT alpha,beta;
	int em = basis::tet(log2p).em;
	int fm = basis::tet(log2p).fm;
	int im = basis::tet(log2p).im;
	
	struct vefi p;
	p.v.resize(maxvst,NV);
	p.e.resize(maxvst,em0,NV);
	p.f.resize(maxvst,fm0,NV);
	p.i.resize(maxvst,im0,NV);
	
	struct vefi z;
	z.v.resize(maxvst,NV);
	z.e.resize(maxvst,em0,NV);
	z.f.resize(maxvst,fm0,NV);
	z.i.resize(maxvst,im0,NV);
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<nfbd;++i)
		hp_fbdry(i)->vdirichlet();
	
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet3d();        
	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet3d();
	
	/* APPLY EDGE DIRICHLET B.C.'S */
    for(int i=0;i<nfbd;++i)
        hp_fbdry(i)->edirichlet();	
	
	for (int i=0;i<nebd;++i) 
		hp_ebdry(i)->edirichlet3d();
	
	/* APPLY FACE DIRICHLET B.C.'S */
	for(int i=0;i<nfbd;++i)
		hp_fbdry(i)->fdirichlet();
	
	z = gbl->res;
	
	jacobi_relaxation();
	//tri_hp::minvrt();
	
	p = gbl->res;
	
	struct vefi Ap;
	Ap.v.resize(maxvst,NV);
	Ap.e.resize(maxvst,em0,NV);
	Ap.f.resize(maxvst,fm0,NV);
	Ap.i.resize(maxvst,im0,NV);
	
	FLT rtrold,rtrnew;
	
	inner_product(rtrold,z,gbl->res);
	
	for(int i = 0; i < 300; ++i) {
		Ap = p;	
		matrix_multiply(Ap);
		FLT ptAp;
		inner_product(ptAp,p,Ap);
		
		alpha = rtrold/ptAp;
		
		ug.v = ug.v - alpha*p.v;
		if(em) ug.e = ug.e - alpha*p.e;
		if(fm) ug.f = ug.f - alpha*p.f;
		if(im) ug.i = ug.i - alpha*p.i;
		
		//		gbl->res.v = gbl->res.v - alpha*Ap.v;
		//		if(sm) gbl->res.s = gbl->res.s - alpha*Ap.s;
		//		if(im) gbl->res.i = gbl->res.i - alpha*Ap.i;
		
		rsdl();		
		
		/* APPLY VERTEX DIRICHLET B.C.'S */
		for(int i=0;i<nfbd;++i)
			hp_fbdry(i)->vdirichlet();
		
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->vdirichlet3d();        
		
		for(int i=0;i<nvbd;++i)
			hp_vbdry(i)->vdirichlet3d();
		
		/* APPLY EDGE DIRICHLET B.C.'S */
		for(int i=0;i<nfbd;++i)
			hp_fbdry(i)->edirichlet();	
		
		for (int i=0;i<nebd;++i) 
			hp_ebdry(i)->edirichlet3d();
		
		/* APPLY FACE DIRICHLET B.C.'S */
		for(int i=0;i<nfbd;++i)
			hp_fbdry(i)->fdirichlet();
		
		z = gbl->res;
		
		jacobi_relaxation();
		//tri_hp::minvrt();
		
		inner_product(rtrnew,z,gbl->res);
		
		beta = rtrnew/rtrold;
		
		p.v = gbl->res.v + beta*p.v;
		if(em) p.e = gbl->res.e + beta*p.e;
		if(fm) p.f = gbl->res.f + beta*p.f;
		if(im) p.i = gbl->res.i + beta*p.i;
		
		rtrold = rtrnew;
		FLT resmax = maxres();
		cout <<  i << ' ' << rtrold << ' ' << resmax << endl;
		if(resmax < 1.0e-10) break;
		
	}
	
	return;
	
}



void tet_hp::matrix_multiply(struct vefi& vec) {
	
#ifdef BZ_DEBUG
	const FLT eps_r = 0.0e-6, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
	const FLT eps_r = 1.0e-6, eps_a = 1.0e-10;  /*<< constants for accurate numerical determination of jacobians */
#endif
	
	int em = basis::tet(log2p).em;
	int fm = basis::tet(log2p).fm;
	int im = basis::tet(log2p).im;
	
	struct vefi temp;
	temp.v.resize(maxvst,NV);
	temp.e.resize(maxvst,em0,NV);
	temp.f.resize(maxvst,fm0,NV);
	temp.i.resize(maxvst,im0,NV);
	
	FLT dw = 0.0;
 	for(int i=0;i<npnt;++i)
		for(int n=0;n<NV;++n)
			dw = MAX(dw,fabs(ug.v(i,n)));
	
	dw = dw*eps_r;
	dw = dw+eps_a;
	
	rsdl();// rsdl may already be called
	
	temp = gbl->res;
	
	ug.v = ug.v + dw*vec.v;
	if(em) ug.e = ug.e + dw*vec.e;
	if(fm) ug.f = ug.f + dw*vec.f;
	if(im) ug.i = ug.i + dw*vec.i;
	
	rsdl();
	
	ug.v = ug.v - dw*vec.v;
	if(em) ug.e = ug.e - dw*vec.e;
	if(fm) ug.f = ug.f - dw*vec.f;
	if(im) ug.i = ug.i - dw*vec.i;
	
	vec.v = (gbl->res.v - temp.v)/dw;
	if(em) vec.e = (gbl->res.e - temp.e)/dw;
	if(fm) vec.f = (gbl->res.f - temp.f)/dw;
	if(im) vec.i = (gbl->res.i - temp.i)/dw;
	
	gbl->res = temp;	
	
	return;
	
}


void tet_hp::inner_product(FLT &alpha, struct vefi vec1,struct vefi vec2 ) {
	
	alpha = 0.0;
	
	for(int i=0; i < npnt; ++i) 
		for(int n=0;n<NV;++n)
			alpha += vec1.v(i,n)*vec2.v(i,n);
	
	for(int i=0; i < nseg; ++i) 
		for(int m=0;m<basis::tet(log2p).em;++m)
			for(int n=0;n<NV;++n)
				alpha += vec1.e(i,m,n)*vec2.e(i,m,n);
		
	for(int i=0; i < ntri; ++i) 
		for(int m=0;m<basis::tet(log2p).fm;++m)
			for(int n=0;n<NV;++n)
				alpha += vec1.f(i,m,n)*vec2.f(i,m,n);
	
	for(int i=0; i < ntet; ++i) 
		for(int m=0;m<basis::tet(log2p).im;++m)
			for(int n=0;n<NV;++n)
				alpha += vec1.i(i,m,n)*vec2.i(i,m,n);
	
	return;	
}