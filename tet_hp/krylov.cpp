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
	
	int max_iter = 50;

	FLT alpha,beta;
	FLT rtrold,rtrnew,ptAp;

	struct vefi p;
	p.v.resize(maxvst,NV);
	p.e.resize(maxvst,em0,NV);
	p.f.resize(maxvst,fm0,NV);
	p.i.resize(maxvst,im0,NV);
	
	struct vefi z(p),Ap(p);
//	z.v.resize(maxvst,NV);
//	z.e.resize(maxvst,em0,NV);
//	z.f.resize(maxvst,fm0,NV);
//	z.i.resize(maxvst,im0,NV);
//	
//	struct vefi Ap;
//	Ap.v.resize(maxvst,NV);
//	Ap.e.resize(maxvst,em0,NV);
//	Ap.f.resize(maxvst,fm0,NV);
//	Ap.i.resize(maxvst,im0,NV);
	
	rsdl();

	all_dirichlet();
	
	z = gbl->res;
	
//#ifdef JACOBI
//	jacobi_relaxation();
//#else
//	tet_hp::minvrt();
//#endif
	
	p = gbl->res;	
	
	inner_product(rtrold,z,gbl->res);

	for(int i = 0; i < max_iter; ++i) {

		matrix_multiply(p,Ap);

		inner_product(ptAp,p,Ap);
		
		alpha = rtrold/ptAp;
		//alpha *= gbl->cfl(log2p);
		
		/* ug = ug - alpha*p */ 
		saxpy(-alpha,1.0,p,ug);
		
		/* res = res-alpha*Ap */
		saxpy(-alpha,1.0,Ap,gbl->res);
		
		//rsdl();		
		
		all_dirichlet();
		
		z = gbl->res;
		
//#ifdef JACOBI
//		jacobi_relaxation();
//#else
//		tet_hp::minvrt();
//#endif
		
		inner_product(rtrnew,z,gbl->res);
		
		beta = rtrnew/rtrold;
		
		/* p = res + beta*p */
		saxpy(1.0,beta,gbl->res,p);
		
		rtrold = rtrnew;
	
		*gbl->log << i << ' ' << rtrold << endl;

		//if(rtrold < 1.0e-15) break;
		
	}
	
	return;
	
}


/* http://math.nist.gov/iml++/bicgstab.h.txt */
/* bi-conjugate gradient stabilized */
void tet_hp::bicgstab() {
	
	int max_iter = 30;
	FLT rho_1, rho_2, alpha, beta, omega;
	
	struct vefi p;
	p.v.resize(maxvst,NV);
	p.e.resize(maxvst,em0,NV);
	p.f.resize(maxvst,fm0,NV);
	p.i.resize(maxvst,im0,NV);
	
	struct vefi phat(p),s(p),shat(p),v(p),vhat(p),t(p),rtilde(p);
//	phat.v.resize(maxvst,NV);
//	phat.e.resize(maxvst,em0,NV);
//	phat.f.resize(maxvst,fm0,NV);
//	phat.i.resize(maxvst,im0,NV);
//	
//	struct vefi s;
//	s.v.resize(maxvst,NV);
//	s.e.resize(maxvst,em0,NV);
//	s.f.resize(maxvst,fm0,NV);
//	s.i.resize(maxvst,im0,NV);
//	
//	struct vefi shat;
//	shat.v.resize(maxvst,NV);
//	shat.e.resize(maxvst,em0,NV);
//	shat.f.resize(maxvst,fm0,NV);
//	shat.i.resize(maxvst,im0,NV);
//	
//	struct vefi v;
//	v.v.resize(maxvst,NV);
//	v.e.resize(maxvst,em0,NV);
//	v.f.resize(maxvst,fm0,NV);
//	v.i.resize(maxvst,im0,NV);
//	
//	struct vefi vhat;
//	vhat.v.resize(maxvst,NV);
//	vhat.e.resize(maxvst,em0,NV);
//	vhat.f.resize(maxvst,fm0,NV);
//	vhat.i.resize(maxvst,im0,NV);
//	
//	struct vefi t;
//	t.v.resize(maxvst,NV);
//	t.e.resize(maxvst,em0,NV);
//	t.f.resize(maxvst,fm0,NV);
//	t.i.resize(maxvst,im0,NV);
//	
//	struct vefi rtilde;
//	rtilde.v.resize(maxvst,NV);
//	rtilde.e.resize(maxvst,em0,NV);
//	rtilde.f.resize(maxvst,fm0,NV);
//	rtilde.i.resize(maxvst,im0,NV);
	
	rsdl();
	
	gbl->res.v = -gbl->res.v;
	gbl->res.e = -gbl->res.e;
	gbl->res.f = -gbl->res.f;
	gbl->res.i = -gbl->res.i;
	rtilde = gbl->res;

	all_dirichlet();
	
	
	p = gbl->res;
	
	inner_product(rho_1,rtilde,gbl->res);

	for (int i = 0; i < max_iter; i++) {
		
		s = gbl->res;//store to use later		

		gbl->res = p;
		
#ifdef JACOBI
		jacobi_relaxation();
#else
		tet_hp::minvrt();
#endif
		
		all_dirichlet();
		
		phat = gbl->res;
		
		matrix_multiply(phat,v);
		
		FLT rtv;
		inner_product(rtv,rtilde,v);
		
		alpha = rho_1/rtv;
		
		saxpy(-alpha, 1.0, v, s);
		
		gbl->res = s;
		
#ifdef JACOBI
		jacobi_relaxation();
#else
		tet_hp::minvrt();
#endif	
		
		all_dirichlet();
		
		shat = gbl->res;
		
		matrix_multiply(shat,t);
		
		FLT ts,tt;
		inner_product(ts,t,s);
		inner_product(tt,t,t);
		
		omega = ts/tt;
		
		ug.v = ug.v + alpha*phat.v + omega*shat.v;
		ug.e = ug.e + alpha*phat.e + omega*shat.e;
		ug.f = ug.f + alpha*phat.f + omega*shat.f;
		ug.i = ug.i + alpha*phat.i + omega*shat.i;
		
//		gbl->res.v = s.v - omega*t.v;
//		gbl->res.e = s.e - omega*t.e;
//		gbl->res.f = s.f - omega*t.f;
//		gbl->res.i = s.i - omega*t.i;
		
		rsdl();		
		gbl->res.v = -gbl->res.v;
		gbl->res.e = -gbl->res.e;
		gbl->res.f = -gbl->res.f;
		gbl->res.i = -gbl->res.i;
		
		all_dirichlet();		
		
		rho_2 = rho_1;	
		
		inner_product(rho_1,rtilde,gbl->res);

		beta = (rho_1/rho_2)*(alpha/omega);
		p.v = gbl->res.v + beta*p.v - beta*omega*v.v;
		p.e = gbl->res.e + beta*p.e - beta*omega*v.e;
		p.f = gbl->res.f + beta*p.f - beta*omega*v.f;
		p.i = gbl->res.i + beta*p.i - beta*omega*v.i;
		
		FLT rtr;
		inner_product(rtr,gbl->res,gbl->res);
		rtr = sqrt(rtr);
		
//		cout << i << ' ' << rtv << ' ' << alpha << ' ' << beta << ' ' << omega << ' ' << rtr << endl;
		cout << i << ' ' <<  rtr << endl;
		
		if (rtr < 1.0e-9) break;
		
	}
	
	return;
	
}


/* http://math.nist.gov/iml++/cgs.h.txt */
/* conjugate gradient squared */
/* still a bug in it */
void tet_hp::cgs() {
	
	int max_iter = 30;
	FLT rho_1, rho_2, alpha, beta;
	
	struct vefi p;
	p.v.resize(maxvst,NV);
	p.e.resize(maxvst,em0,NV);
	p.f.resize(maxvst,fm0,NV);
	p.i.resize(maxvst,im0,NV);
	
	struct vefi phat;
	phat.v.resize(maxvst,NV);
	phat.e.resize(maxvst,em0,NV);
	phat.f.resize(maxvst,fm0,NV);
	phat.i.resize(maxvst,im0,NV);
	
	struct vefi q;
	q.v.resize(maxvst,NV);
	q.e.resize(maxvst,em0,NV);
	q.f.resize(maxvst,fm0,NV);
	q.i.resize(maxvst,im0,NV);
	
	struct vefi qhat;
	qhat.v.resize(maxvst,NV);
	qhat.e.resize(maxvst,em0,NV);
	qhat.f.resize(maxvst,fm0,NV);
	qhat.i.resize(maxvst,im0,NV);
	
	struct vefi vhat;
	vhat.v.resize(maxvst,NV);
	vhat.e.resize(maxvst,em0,NV);
	vhat.f.resize(maxvst,fm0,NV);
	vhat.i.resize(maxvst,im0,NV);
	
	struct vefi u;
	u.v.resize(maxvst,NV);
	u.e.resize(maxvst,em0,NV);
	u.f.resize(maxvst,fm0,NV);
	u.i.resize(maxvst,im0,NV);
	
	struct vefi uhat;
	uhat.v.resize(maxvst,NV);
	uhat.e.resize(maxvst,em0,NV);
	uhat.f.resize(maxvst,fm0,NV);
	uhat.i.resize(maxvst,im0,NV);
	
	struct vefi rtilde;
	rtilde.v.resize(maxvst,NV);
	rtilde.e.resize(maxvst,em0,NV);
	rtilde.f.resize(maxvst,fm0,NV);
	rtilde.i.resize(maxvst,im0,NV);
		
	rsdl();
	
	gbl->res.v = -gbl->res.v;
	gbl->res.e = -gbl->res.e;
	gbl->res.f = -gbl->res.f;
	gbl->res.i = -gbl->res.i;
	
	rtilde = gbl->res;

	all_dirichlet();
	
	
	for (int i = 0; i <= max_iter; i++) {
		
		inner_product(rho_1,rtilde,gbl->res);

		if (i == 0) {
			u = gbl->res;
			p = u;
		} else {
			beta = rho_1/rho_2;
			u.v = gbl->res.v + beta*q.v;
			p.v = u.v + beta*q.v+beta*beta*p.v;
		}
		
		gbl->res = p;
		
#ifdef JACOBI
		jacobi_relaxation();
#else
		tet_hp::minvrt();
#endif
		all_dirichlet();
		
		matrix_multiply(gbl->res,vhat);
		
		gbl->res = vhat;
		all_dirichlet();
		vhat = gbl->res;
		
		FLT rtv;
		inner_product(rtv,rtilde,vhat);
		
		alpha = rho_1/rtv;
		
		q.v = u.v - alpha*vhat.v;
		
		gbl->res.v = u.v + q.v;
		
#ifdef JACOBI
		jacobi_relaxation();
#else
		tet_hp::minvrt();
#endif	
		
		all_dirichlet();
		
		saxpy(alpha,1.0,gbl->res,ug);
		
		rsdl();
		gbl->res.v = -gbl->res.v;
		gbl->res.e = -gbl->res.e;
		gbl->res.f = -gbl->res.f;
		gbl->res.i = -gbl->res.i;
		//gbl->res.v = -gbl->res.v;

		all_dirichlet();
		
//		matrix_multiply(gbl->res, qhat)
//		saxpy(alpha,1.0,qhat,gbl->res_store)
		
		rho_2 = rho_1;	
		
		FLT rtr;
		inner_product(rtr,gbl->res,gbl->res);
		cout << i << ' ' << rtr << endl;

	}
	
	return;
}



//void tet_hp::gmres(){
//	int i,i1,k,k1,its;
//	FLT t,t0,tt,beta,gam,negt;
//	FLT eps1 = 0.0;
//		
//	int im = 20;
//	FLT tol = 1.0e-10;
//	int itmax = 50;
//	int maxits = itmax;
//
//	int em = basis::tet(log2p).em;
//	int fm = basis::tet(log2p).fm;
//	int im = basis::tet(log2p).im;
//	
//	rsdl();
//	
//	/* APPLY VERTEX DIRICHLET B.C.'S */
//	for(int i=0;i<nfbd;++i)
//		hp_fbdry(i)->vdirichlet();
//	
//	for(int i=0;i<nebd;++i)
//		hp_ebdry(i)->vdirichlet3d();        
//	
//	for(int i=0;i<nvbd;++i)
//		hp_vbdry(i)->vdirichlet3d();
//	
//	/* APPLY EDGE DIRICHLET B.C.'S */
//    for(int i=0;i<nfbd;++i)
//        hp_fbdry(i)->edirichlet();	
//	
//	for (int i=0;i<nebd;++i) 
//		hp_ebdry(i)->edirichlet3d();
//	
//	/* APPLY FACE DIRICHLET B.C.'S */
//	for(int i=0;i<nfbd;++i)
//		hp_fbdry(i)->fdirichlet();
//	
//	vefi rhs;
//	rhs.v.resize(npnt,NV);
//	rhs.e.resize(npnt,em,NV);
//	rhs.f.resize(npnt,fm,NV);
//	rhs.i.resize(npnt,im,NV);
//	rhs.v = gbl->res.v;
//	rhs.e = gbl->res.e;
//	rhs.f = gbl->res.f;
//	rhs.i = gbl->res.i;
//	
//	Array<vefi,1> vv(im+1);
//	for(int ind = 0; ind < im+1; ++ind) {	
//		vv(ind).v.resize(npnt,NV);
//		vv(ind).e.resize(npnt,em,NV);
//		vv(ind).f.resize(npnt,fm,NV);
//		vv(ind).i.resize(npnt,im,NV);
//	}
//		
//	Array<Array<double,1>,1> hh(im);
//	for(int ind = 0; ind < im; ++ind) 
//		hh(ind).resize(im+1);
//	
//	Array<vefi,1> z(im);
//	for(int ind = 0; ind < im; ++ind) {	
//		z(ind).v.resize(npnt,NV);
//		z(ind).e.resize(npnt,em,NV);
//		z(ind).f.resize(npnt,fm,NV);
//		z(ind).i.resize(npnt,im,NV);
//	}
//	
//	Array<double,1> c(im);
//	Array<double,1> s(im);
//	Array<double,1> rs(im+1);
//	
//	
//	its = 0;
//	
//	/* outer loop starts here */
//	do{
//		
//		/* compute initial residual */
//		rsdl();
//		
//		/* APPLY VERTEX DIRICHLET B.C.'S */
//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->vdirichlet();
//		
//		for(int i=0;i<nebd;++i)
//			hp_ebdry(i)->vdirichlet3d();        
//		
//		for(int i=0;i<nvbd;++i)
//			hp_vbdry(i)->vdirichlet3d();
//		
//		/* APPLY EDGE DIRICHLET B.C.'S */
//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->edirichlet();	
//		
//		for (int i=0;i<nebd;++i) 
//			hp_ebdry(i)->edirichlet3d();
//		
//		/* APPLY FACE DIRICHLET B.C.'S */
//		for(int i=0;i<nfbd;++i)
//			hp_fbdry(i)->fdirichlet();
//		
//		vv(0).v = gbl->res.v;
//		vv(0).e = gbl->res.e;
//		vv(0).f = gbl->res.f;
//		vv(0).i = gbl->res.i;		
//		
//		vv(0).v = rhs.v - vv(0).v;
//		vv(0).e = rhs.e - vv(0).e;
//		vv(0).f = rhs.f - vv(0).f;
//		vv(0).i = rhs.i - vv(0).i;
//		
//		inner_product(beta, vv(0), vv(0));
//
//		//cout << "initial residual of system Ax-b: " << beta << endl;
//		FLT rtr;
//		inner_product(rtr, gbl->res, gbl->res);
//
//		if ( !(beta > tol * rtr) )
//			break;
//		
//		t = 1.0/beta;
//		
//		/* normalize vv(0) = vv(0)/beta */
//		vv(0).v *= t;
//		vv(0).e *= t;
//		vv(0).f *= t;
//		vv(0).i *= t;
//		
//		if(its == 0)
//			eps1 = tol*beta;
//		
//		/* initialize 1st term of rhs of hessenberg system */
//		rs(0) = beta;
//		for(i = 0; i < im; i++){
//			its++;
//			i1 = i+1;
//			
////			/* Right preconditioning LUz=Mz=v -> z=(LU)^(-1)v*/
////			dpsolve(n, L, U, perm_c, perm_r, z(i), vv(i), stat);
////			/* instead of preconditioning copy z=v */
//			z(i).v = vv(i).v;
//			z(i).e = vv(i).e;
//			z(i).f = vv(i).f;
//			z(i).i = vv(i).i;
//
//			/* matvec operation w = A z_{j} = A M^{-1} v_{j} */
//			sp_dgemv("N",1.0,&A,z(i).data(),1,0.0,vv(i1).data(),1);
//			matrix_multiply(z(i))
//			/*------------------------------------------------------------
//			 |     modified gram - schmidt...
//			 |     h_{i,j} = (w,v_{i})
//			 |     w  = w - h_{i,j} v_{i}
//			 +------------------------------------------------------------*/
//			
//			t0=nrm2(n,vv(i1));
//			
//			for(int j = 0; j <= i; j++){
//				tt=dotprod(n,vv(j),vv(i1));
//				hh(i)(j) = tt;
//				negt = -tt;
//				vv(i1)+=negt*vv(j); //daxpy				
//			}
//			/* h_{j+1,j} = ||w||_{2} */
//			t=nrm2(n,vv(i1));
//			while ( t < 0.5*t0){
//				t0 = t;
//				for(int j = 0; j <=i; j++){
//					tt=dotprod(n,vv(j),vv(i1));
//					hh(i)(j) += tt;
//					negt = -tt;
//					vv(i1)+=negt*vv(j); //daxpy
//				}	
//				t=nrm2(n,vv(i1));
//			}
//			
//			hh(i)(i1) = t;
//			
//			if (t != 0.0){
//				t=1.0/t;
//				vv(i1) *= t;
//			}
//			/*---------------------------------------------------
//			 |     done with modified gram schimdt and arnoldi step
//			 |     now  update factorization of hh
//			 +--------------------------------------------------*/
//			
//			/*--------------------------------------------------------
//			 |   perform previous transformations  on i-th column of h
//			 +-------------------------------------------------------*/
//			
//			for(int k = 1; k <= i; k++){
//				k1 = k-1;
//				tt = hh(i)(k1);
//				hh(i)(k1) = c(k1)*tt+s(k1)*hh(i)(k);
//				hh(i)(k) = -s(k1)*tt+c(k1)*hh(i)(k);
//			}
//			
//			gam = sqrt(hh(i)(i)*hh(i)(i)+hh(i)(i1)*hh(i)(i1));
//			
//			/*---------------------------------------------------
//			 |     if gamma is zero then any small value will do
//			 |     affect only residual estimate
//			 +--------------------------------------------------*/
//			if( gam == 0.0){
//				c(i) = 1.0;
//				s(i) = 0.0;
//			}
//			else {
//				c(i) = hh(i)(i)/gam;
//				s(i) = hh(i)(i1)/gam;
//			}
//			
//			rs(i1) = -s(i)*rs(i);
//			rs(i) = c(i)*rs(i);
//			
//			/*----------------------------------------------------
//			 |   determine residual norm and test for convergence
//			 +---------------------------------------------------*/
//			
//			hh(i)(i) = c(i)*hh(i)(i) + s(i)*hh(i)(i1);
//			beta = fabs(rs(i1));
//			if (beta <= eps1 || its >= itmax){
//				break;
//				cout << "break 0" << endl;
//			}
//		}
//		
//		if (i == im) i--;
//		
//		/* now compute solution. 1st solve upper trianglular system */
//		rs(i) /= hh(i)(i);
//		
//		for(int ii = 1; ii <= i; ii++){
//			k = i-ii;
//			k1 = k+1;
//			tt = rs(k);
//			for(int j = k1;j <= i; j++)
//				tt -= hh(j)(k)*rs(j);
//			rs(k) = tt/hh(k)(k);
//		}
//		/* linear combination of v[i]'s to get sol. */
//		for(int j = 0; j <= i; j++){
//			tt = rs(j);
//			for(int k = 0; k < n; k++)
//				sol(k) += tt*z(j)(k);
//		}
//		
//		/* calculate the residual and output vv=A*sol*/
//		sp_dgemv("N",1.0,&A,sol.data(),1,0.0,vv(0).data(),1);
//		
//		for(int j = 0; j < n; j++)
//			vv(0)(j)=rhs(j)-vv(0)(j);
//		
//		beta=nrm2(n,vv(0));
//		cout << "gmres beta " << beta << endl;
//		
//		/* not sure if I should do this or not */
//		if(beta < tol)
//			break;
//		/*---- restart outer loop if needed ----*/
//		
//		if ( !(beta < eps1 / tol) )	{
//			its = maxits + 10;
//			break;
//			cout << "break 1" << endl;
//		}
//		if (beta <= eps1){
//			break;
//			cout << "break 2" << endl;
//		}
//		
//	} while(its < maxits);
//	
//	itmax = its;	
//	
//	return (its >= maxits);
//}



void tet_hp::matrix_multiply(struct vefi vecx, struct vefi& vecb) {
	
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
	
	saxpy(dw,1.0,vecx,ug);
	
	rsdl();
	
	saxpy(-dw,1.0,vecx,ug);
	
	vecb.v = (gbl->res.v - temp.v)/dw;
	if(em) vecb.e = (gbl->res.e - temp.e)/dw;
	if(fm) vecb.f = (gbl->res.f - temp.f)/dw;
	if(im) vecb.i = (gbl->res.i - temp.i)/dw;
	
	gbl->res = temp;	
	
	return;
	
}

void tet_hp::saxpy(FLT alpha, FLT beta, struct vefi vecx,struct vefi& vecy ) {
	
	vecy.v = alpha*vecx.v + beta*vecy.v;
	if(basis::tet(log2p).em) vecy.e = alpha*vecx.e + beta*vecy.e;
	if(basis::tet(log2p).fm) vecy.f = alpha*vecx.f + beta*vecy.f;
	if(basis::tet(log2p).im) vecy.i = alpha*vecx.i + beta*vecy.i;

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

void tet_hp::all_dirichlet() {

	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int j=0;j<nfbd;++j)
		hp_fbdry(j)->vdirichlet();
	
	for(int j=0;j<nebd;++j)
		hp_ebdry(j)->vdirichlet3d();        
	
	for(int j=0;j<nvbd;++j)
		hp_vbdry(j)->vdirichlet3d();
	
	/* APPLY EDGE DIRICHLET B.C.'S */
	for(int j=0;j<nfbd;++j)
		hp_fbdry(j)->edirichlet();	
	
	for (int j=0;j<nebd;++j) 
		hp_ebdry(j)->edirichlet3d();
	
	/* APPLY FACE DIRICHLET B.C.'S */
	for(int j=0;j<nfbd;++j)
		hp_fbdry(j)->fdirichlet();
}