/*
 *  utilities.cpp
 *  tri_hp
 *
 *  Created by michael brazell on 1/11/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */

#include "tri_hp.h"
#include <math.h>
#include <utilities.h>
#include <myblas.h>

FLT tri_hp::spectral_radius(Array<FLT,2> A, int n){
	Array<FLT,1> x(n),temp(n);
	FLT lambda;
	x=1;
	
	/* Power Iteration to find eigenvector for dominate eigenvalue */
	for(int m=0;m<10;++m){
		temp=0;
		for(int i=0; i<n; ++i){
			for(int j=0; j<n; ++j){
				temp(i)+=A(i,j)*x(j);
			}
		}
		x=temp/l2norm(temp,n);
	}
	
	/* Spectral radius using the Rayleigh quotient */
	lambda = inner_product(x, temp, n)/inner_product(x, x, n);
	
	return(lambda);
}

FLT tri_hp::l2norm(Array<FLT,1> x, int n){
	FLT temp=0.0;
	
	for(int i=0;i<n;++i)
		temp+=x(i)*x(i);	
	
	return(pow(temp, 0.5));
}

FLT tri_hp::inner_product(Array<FLT,1> x, Array<FLT,1> y, int n){
	FLT temp=0.0;
	
	for(int i=0;i<n;++i)
		temp+=x(i)*y(i);	
	
	return(temp);
}

void tri_hp::matrix_absolute_value(Array<FLT,2> &A, int n){
	Array<double,1> lambda_real(n),lambda_imag(n);
	Array<double,2> VR(n,n),temp(n,n);
	int lwork=4*n,info,ipiv[n];
	char vchar[] = "V";
	char nchar[] = "N";
	char trans[] = "T";
	double work[lwork];

	/* find eigenvalues and eigenvectors */
	GEEV(nchar,vchar, n, A.data(), n, lambda_real.data(), lambda_imag.data(), temp.data(), n, VR.data(), n, work, lwork, info);
	
	/* absolute value of eigenvalues */
	for (int i=0; i < n; ++i) 
		lambda_real(i) = sqrt(pow(lambda_real(i),2.0)+pow(lambda_imag(i),2.0));
	
	/* |diag(lambda)|V^T (transpose temp so that GETRF is compatible ) */
	for (int i = 0; i < n; ++i) 
		for (int j = 0; j < n;++j) 
			temp(j,i)=lambda_real(i)*VR(i,j);
	
	/*  LU factorization  */
	GETRF(n, n, VR.data(), n, ipiv, info);

	/* Solve transposed system temp = inv(VR)*temp */
	GETRS(trans,n,n,VR.data(),n,ipiv,temp.data(),n,info);
	
	/* store transpose of temp in A */
	for (int i = 0; i < n; ++i) 
		for (int j = 0; j < n;++j) 
			A(i,j)=temp(j,i);
		
	return;
}

