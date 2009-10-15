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
	Array<double,1> du(size_sparse_matrix);

	jacobian_matrix J;
	/* send global solution to ug_vec */
	ug_to_vec();
	
	for(int i = 0; i < max_newton_its; ++i) {
		
		/* zero out sparse matrix ija and residual res_vec */
		zero_sparse();
		
		/* insert values into jacobian J: 0 is compressed row */
		create_jacobian(false);
		
		/* insert values into residual res_vec */ 
		create_rsdl();
		
		/* apply dirichlet boundary conditions to sparse matrix and vector */
		for(int j = 0; j < nfbd; ++j)
			fbdry(j)->apply_sparse_dirichlet(false);
		
		/* solve system with gmres J*du=res_vec */
		its = gmres(inner_its,size_sparse_matrix,J,res_vec.data(),du.data(),tol);
		
		/* update solution */
//		for(int i = 0; i < size_sparse_matrix; ++i)
//			ug_vec(i)-=du[i];
		ug_vec-=du;
		
		/* send ug_vec to global solution */
		vec_to_ug();
		
		newton_norm = 0.0;			
		for (int j = 0; j < size_sparse_matrix; ++j)
			newton_norm += du(j)*du(j);
		
		if (newton_norm < tol*tol) break;
		
	}

	
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
		for(int j = sparse_ptr(i); j < sparse_ptr(i+1); ++j)
			vec_out[i] += sparse_val(j)*vec_in[sparse_ind(j)];
	}
}


#endif