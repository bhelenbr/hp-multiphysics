/*
 *  lsolver.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 10/8/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#ifndef petsc


#include "tet_hp.h"
#include "gmres.h"

void tet_hp::lsolver(){
	
	int its = -1;
	int max_newton_its = 100;
	int inner_its = 5;
	FLT newton_norm,tol=1.0e-12;
//	Array<double,1> du(size_sparse_matrix);
	double *du = new double[size_sparse_matrix];
	double *b = new double[size_sparse_matrix];

	/* create jacobian struct */
	jacobian_matrix J();
	
	/* send global solution to ug_vec */
	ug_to_vec();
	
	for(int i = 0; i < max_newton_its; ++i) {
		
		/* zero out sparse matrix ija and residual res_vec */
		zero_sparse();
		
		/* insert values into jacobian J (ija) */
		create_jacobian();
		
		/* insert values into residual res_vec */ 
		create_rsdl();
		
		for(int i = 0; i < size_sparse_matrix; ++i)
			b[i]=res_vec(i);
		
		/* solve system with gmres J*du=res_vec */
//		its = gmres(inner_its,size_sparse_matrix,J,res_vec,du,tol);
		//its = gmres(inner_its,size_sparse_matrix,J,b,du,tol);// wont work just compiling
		mult(J,b,du);//testing delete later
		
		/* update solution */
		for(int i = 0; i < size_sparse_matrix; ++i)
			ug_vec(i)-=du[i];

//		ug_vec-=du;
		
		/* send ug_vec to global solution */
		vec_to_ug();
		
		newton_norm = 0.0;			
		for (int j = 0; j < size_sparse_matrix; ++j)
			newton_norm += du[j]*du[j];

//			newton_norm += du(j)*du(j);
		
		if (newton_norm < tol*tol) break;
		
	}

	delete[] du;
	delete[] b;
	
	return;
}

//struct jacobian_matrix{
//	jacobian_matrix();
//};

///* row-indexed matrix multiplication*/
//void tet_hp::mult(const jacobian_matrix &J, const double *vec_in, double *vec_out){
//	for(int i = 0; i < size_sparse_matrix; ++i){
//		vec_out[i] = sa(i)*vec_in[i];
//		for(int j = ija(i); j < ija(i+1); ++j)
//			vec_out[i] += sa(j)*vec_in[ija(j)];
//	}
//}

/* compressed row matrix multiplication*/
void tet_hp::mult(const jacobian_matrix &J, const double *vec_in, double *vec_out){
	for(int i = 0; i < size_sparse_matrix; ++i){
		vec_out[i] = 0.0;
		for(int j = row_ptr(i); j < row_ptr(i+1); ++j)
			vec_out[i] += sval(j)*vec_in[col_ind(j)];
	}
}


#endif