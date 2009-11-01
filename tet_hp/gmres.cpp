/*
 *  gmres.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 10/31/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "slu_ddefs.h"

void tet_hp::fgmres(int n,SuperMatrix &A, SuperMatrix &L,SuperMatrix &U, int &perm_c, int &perm_r, Array<double,1> &rhs,Array<double,1> &sol,FLT tol,int im,int &itmax,SuperLUStat_t  &stat){
	int i,i1,k1,its;
	FLT t,t0,tt,beta,gam,fltwork,negt;
	FLT eps1 = 0.0;
	int maxits = itmax;
	Array<Array<double,1>,1> vv(im+1);
	for(int ind = 0; ind < im+1; ++ind)	vv(ind).resize(n);
	Array<Array<double,1>,1> hh(im);
	for(int ind = 0; ind < im; ++ind) hh(ind).resize(im+1);

	Array<Array<double,1>,1> z(im);
	for(int ind = 0; ind < im; ++ind) z(ind).resize(n);
	
	Array<double,1> c(im);
	Array<double,1> s(im);
	Array<double,1> rs(im+1);

	its = 0;
	
	do{
//		sp_dgemv (char *, double, SuperMatrix *, double *,
//				  int, double, double *, int);
		/* compute initial residual */
		sp_dgemv("N",1.0, &A, sol.data(),1,0.0,vv(0).data(),1);
		vv(0)=rhs-vv(0);
		beta=0.0;
		for(int ind=0;ind<n;++ind) beta+=vv(0)(ind)*vv(0)(ind);
		beta=sqrt(beta);
		
		fltwork=0.0;
		for(int ind=0;ind<n;++ind) fltwork+=rhs(ind)*rhs(ind);
		fltwork=sqrt(fltwork);
		
		if ( !(beta > tol * fltwork) )
			break;
		
		t = 1.0/beta;
		
		/* normalize vv(0) = vv(0)/beta */
		vv(0) *= t;
		
		if(its == 0)
			eps1 = tol*beta;
		
		/* initialize 1st term of rhs of hessenberg system */
		rs(0) = beta;
		for(i = 0; i < im; i++){
			its++;
			i1 = i+1;
			
			/* Right preconditioning M*z=v */
			dpsolve(n, L, U, perm_c, perm_r, z(i), vv(i), stat);
			
			/* matvec operation w = A z_{j} = A M^{-1} v_{j} */
			sp_dgemv("N",1.0,&A,z(i).data(),1,0.0,vv(i1).data(),1);
			
			/*------------------------------------------------------------
			 |     modified gram - schmidt...
			 |     h_{i,j} = (w,v_{i})
			 |     w  = w - h_{i,j} v_{i}
			 +------------------------------------------------------------*/
			
			t0 = 0.0;
			for(int ind=0;ind<n;++ind) t0+=vv(i1)(ind)*vv(i1)(ind);
			t0=sqrt(t0);
			
			for(int j = 0; j <= i; j++){
				tt = 0.0;
				for(int ind=0;ind<n;++ind) tt+=vv(j)(ind)*vv(i1)(ind);
				hh(i)(j) = tt;
				negt = -tt;
				vv(i1)+=negt*vv(j); //daxpy				
			}
			/* h_{j+1,j} = ||w||_{2} */
			t = 0.0;
			for(int ind=0;ind<n;++ind) t+=vv(i1)(ind)*vv(i1)(ind);
			t = sqrt(t);
			while ( t < 0.5*t0){
				t0 = t;
				for(int j = 0; j <=i; j++){
					tt = 0.0;
					for(int ind=0;ind<n;++ind) t+=vv(j)(ind)*vv(i1)(ind);
					hh(i)(j) += tt;
					negt = -tt;
					vv(i1)+=negt*vv(j);
				}
				t = 0.0;
				for(int ind=0;ind<n;++ind) t+=vv(i1)(ind)*vv(i1)(ind);
				t = sqrt(t);				
			}
			
			hh(i)(i1) = t;
			
			if (t != 0.0){
				t=1.0/t;
				vv(i1) *= t;
			}
			/*---------------------------------------------------
			 |     done with modified gram schimdt and arnoldi step
			 |     now  update factorization of hh
			 +--------------------------------------------------*/
			
			/*--------------------------------------------------------
			 |   perform previous transformations  on i-th column of h
			 +-------------------------------------------------------*/
			
			for(int k = 1; k <= i; k++){
				k1 = k-1;
				tt = hh(i)(k1);
				hh(i)(k1) = c(k1)*tt+s(k1)*hh(i)(k);
				hh(i)(k) = -s(k1)*tt+c(k1)*hh(i)(k);
			}
			
			gam = sqrt(hh(i)(i)*hh(i)(i)+hh(i)(i1)*hh(i)(i1));
			
			/*---------------------------------------------------
			 |     if gamma is zero then any small value will do
			 |     affect only residual estimate
			 +--------------------------------------------------*/
			if( gam == 0.0){
				c(i) = 1.0;
				s(i) = 0.0;
			}
			else {
				c(i) = hh(i)(i)/gam;
				s(i) = hh(i)(i1)/gam;
			}
			
			rs(i1) = -s(i)*rs(i);
			rs(i) = c(i)*rs(i);
			
			/*----------------------------------------------------
			 |   determine residual norm and test for convergence
			 +---------------------------------------------------*/
			
			hh(i)(i) = c(i)*hh(i)(i) + s(i)*hh(i)(i1);
			beta = fabs(rs(i1));
			if (beta <= eps1 || its >= itmax)
				break;
		}
		
		if (i == im) i--;
		
		/* now compute solution. 1st solve upper trianglular system */
		rs(i) /= hh(i)(i);
		
		for(int ii = 1; ii <= i; ii++){
			int k = i-ii;
			k1 = k+1;
			tt = rs(k);
			for(int j = k1;j <= i; j++)
				tt -= hh(j)(k)*rs(j);
			rs(k) = tt/hh(k)(k);
		}
		
		/* linear combination of v[i]'s to get sol. */
		for(int j = 0; j <= i; j++){
			tt = rs(j);
			for(int k = 0; k < n; k++)
				sol(k) += tt*z(j)(k);
		}
		
		/* calculate the residual and output */
		sp_dgemv("N",1.0,&A,sol.data(),1,0.0,vv(0).data(),1);
		for(int j = 0; j < n; j++)
			vv(0)(j)=rhs(j)-vv(0)(j);
		
		beta = 0.0;
		for(int ind=0;ind<n;++ind) beta+=vv(0)(ind)*vv(0)(ind);
		beta = sqrt(beta);
		
		/*---- restart outer loop if needed ----*/

		if ( !(beta < eps1 / tol) )	{
			its = maxits + 10;
			break;
		}
		if (beta <= eps1)
			break;
				
	} while(its < maxits);
	
	itmax = its;	
	
	return;
}

void tet_hp::dpsolve(int n, SuperMatrix &L, SuperMatrix &U, int &perm_c, int &perm_r,Array<double,1> &x, Array<double,1> &y, SuperLUStat_t &stat){
	
    int info;
    static DNformat X;
    static SuperMatrix XX = {SLU_DN, SLU_D, SLU_GE, 1, 1, &X};
	
	x=y;
	
    XX.nrow = n;
    X.lda = n;
    X.nzval = x.data();
    dgstrs(NOTRANS, &L, &U, &perm_c, &perm_r, &XX, &stat, &info);
}

