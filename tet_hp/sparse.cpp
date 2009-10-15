/*
 *  sparse.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 9/14/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "ins/tet_hp_ins.h"
#include "gmres.h"


#ifndef petsc
/* compressed row storage */
void tet_hp::insert_sparse(int row, int col, FLT value, bool compressed_col ){

	/* to make compressed column switch row and column */
	if (compressed_col) {
		int temp = row;
		row = col;
		col = temp;
	}

	
	int find_ind = -1;
	
	/* add value to already existing element in sparse */
	for(int i = sparse_ptr(row); i < sparse_ptr(row+1); ++i){
		if(sparse_ind(i) == col){
			sparse_val(i) += value;
			return;
		}
	}
	
	/* if not pre-exisiting element then add to list in correct location */
	for(int i = sparse_ptr(row); i < sparse_ptr(row+1); ++i){
		if(sparse_ind(i) > col){
			find_ind = i;
			break;
		}
	}
	/* no column bigger so insert at end */
	if(find_ind == -1)
		find_ind = sparse_ptr(row+1);
	
	++number_sparse_elements;
	
	/* shift start of column index by one */
	for(int i = row+1; i <= size_sparse_matrix+1; ++i)
		++sparse_ptr(i);
	
	/* slide all non-diagonal entries after insertion column to the right (slow but only done once) */
	for(int i = number_sparse_elements; i > find_ind; --i){
		sparse_val(i)=sparse_val(i-1);
		sparse_ind(i)=sparse_ind(i-1);
	}
	
	/* insert new element into sparse matrix */
	sparse_val(find_ind) = value;
	sparse_ind(find_ind) = col;	
	
	return;	
}

/* zero out sparse matrix but keep same structure */
void tet_hp::zero_sparse(){
	
	for (int i = 0; i < number_sparse_elements; ++i){
		res_vec(i) = 0.0;
		sparse_val(i) = 0.0;
	}
	
	return;
}

void tet_hp::initialize_sparse(){
	
	size_sparse_matrix = (npnt+nseg*em0+ntri*fm0+ntet*im0)*NV;

	/* sparse matrix allocation */
	sparse_ind.resize(MXTM*NV*ntet);//too much storage resize later
	sparse_val.resize(MXTM*NV*ntet);
	sparse_ptr.resize(size_sparse_matrix+1);
	number_sparse_elements = size_sparse_matrix;
	
	/* creates sparse matrix with zeros on diagonal */
	for (int i = 0; i < size_sparse_matrix+1; ++i) {
		sparse_ptr(i) = i;
		sparse_ind(i) = i;
	}
	
	return;
}

/* clears row in sparse matrix, inserts 1.0 on diagonal, and inserts 0.0 in the residual */
void tet_hp::sparse_dirichlet(int ind, bool compressed_col){
	
	res_vec(ind) = 0.0;

	if(compressed_col){
		/* compressed column: only works for symmetric sparsity pattern */
		for(int i = sparse_ptr(ind); i < sparse_ptr(ind+1); ++i){
			for(int j = sparse_ptr(sparse_ind(i)); j < sparse_ptr(sparse_ind(i)+1); ++j){
				if(ind == sparse_ind(j)){
					sparse_val(j) = 0.0;
					break;
				}
			}	
			if(ind == sparse_ind(i))
				sparse_val(i) = 1.0;
		}
	}
	else {
		/* compressed row */
		for(int i = sparse_ptr(ind); i < sparse_ptr(ind+1); ++i){
			sparse_val(i) = 0.0;
			if(ind == sparse_ind(i))
				sparse_val(i) = 1.0;
		}
	}
	return;
}


#endif

#ifdef petsc

//void tet_hp::petsc_dirichlet(int ind){
//	FLT zero = 0.0;
//	FLT one = 1.0;
//	
//	ierr = MatSetValues(petsc_J,1,&ind,1,&loc_to_glo(j),&zero,INSERT_VALUES);
//	return;
//}
#endif

void tet_hp::create_jacobian(bool jac_tran) {
	int indx,gindx,eind,find,iind,sgn,msgn,mode;
	
	int kn = basis::tet(log2p).tm*NV;
	Array<FLT,2> K(kn,kn);
	Array<int,1> loc_to_glo(kn);
	
	for(int tind = 0; tind < ntet; ++tind){	
		
		int ind = 0;

		create_local_jacobian_matrix(tind, K);
		
		for (int m = 0; m < 4; ++m) {
			gindx = NV*tet(tind).pnt(m);
			for (int n = 0; n < NV; ++n)
				loc_to_glo(ind++) = gindx+n;
		}		
		
		/* EDGE MODES */
		if (basis::tet(log2p).p > 1) {
			for(int i = 0; i < 6; ++i) {
				eind = npnt*NV + tet(tind).seg(i)*basis::tet(log2p).em*NV;
				sgn = tet(tind).sgn(i);
				msgn = 1;
				for (int m = 0; m < basis::tet(log2p).em; ++m) {
					for(int n = 0; n < NV; ++n) {
						for(int j = 0; j < kn; ++j) {
							K(ind,j) *= msgn;
							K(j,ind) *= msgn;
						}
						loc_to_glo(ind++) = eind + m*NV + n;
					}
					msgn *= sgn;
				}
			}
		}
		
		/* FACE MODES */
		if (basis::tet(log2p).p > 2) {
			for(int i = 0; i < 4; ++i){
				sgn = -tet(tind).rot(i);
				find = npnt*NV+nseg*basis::tet(log2p).em*NV+tet(tind).tri(i)*basis::tet(log2p).em*NV;
				mode = 0;
				msgn = 1;		
				for(int m = 1; m <= basis::tet(log2p).em-1; ++m) {
					for(int j = 1; j <= basis::tet(log2p).em-m; ++j){
						for(int n = 0; n < NV; ++n){
							for(int k = 0; k < kn; ++k){
								K(ind,k) *= msgn;
								K(k,ind) *= msgn;
							}
							loc_to_glo(ind++) = find+mode*NV+n;
						}
						++mode;
					}
					msgn *= sgn;
				}
			}		
		}		
		
		/* INTERIOR MODES */
		iind = npnt*NV+nseg*basis::tet(log2p).em*NV+ntri*basis::tet(log2p).fm*NV+tind*basis::tet(log2p).im*NV;
		for(int m = 0; m < basis::tet(log2p).im; ++m) {
			for(int n = 0; n < NV; ++n){
				for(int k = 0; k < kn; ++k)
					K(ind,k) = K(ind,k);
				loc_to_glo(ind++) = iind+m*NV+n;
			}
		}

		
#ifdef petsc

		PetscErrorCode ierr;
		for(int i = 0; i < kn; ++i)
			for(int j = 0; j < kn; ++j)
				ierr = MatSetValues(petsc_J,1,&loc_to_glo(i),1,&loc_to_glo(j),&K(i,j),ADD_VALUES);
		// CHKERRQ(ierr);
#endif

#ifndef petsc
		for(int i = 0; i < kn; ++i)
			for(int j = 0; j < kn; ++j)
				insert_sparse(loc_to_glo(j), loc_to_glo(i), K(i,j),jac_tran);

#endif

				
				
	}
	
#ifndef petsc
	/* resize and preserve sparse matrix on first call only */
	if (!sparse_resized) {
		
		//ija.resizeAndPreserve(number_sparse_elements+1);
		//sa.resizeAndPreserve(number_sparse_elements+1);

		sparse_ind.resizeAndPreserve(number_sparse_elements);
		sparse_val.resizeAndPreserve(number_sparse_elements);
		sparse_resized = true;
	}
#endif
	return;
}

void tet_hp::create_local_jacobian_matrix(int tind, Array<FLT,2> &K) {
	Array<TinyVector<FLT,MXTM>,1> R(NV),Rbar(NV),lf_re(NV),lf_im(NV);
	int kcol = 0;
	FLT dw = 0.001;  //dw=sqrt(eps/l2_norm(q))
	
	
	ugtouht(tind);
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tet(log2p).tm;++i)
		for(int n=0;n<NV;++n)
			Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);

	
	for(int mode = 0; mode < basis::tet(log2p).tm; ++mode){
		for(int var = 0; var < NV; ++var){
			uht(var)(mode) += dw;
			element_rsdl(tind,0,uht,lf_re,lf_im);
			for(int i=0;i<basis::tet(log2p).tm;++i)
				for(int n=0;n<NV;++n)
					R(n)(i)=lf_re(n)(i)+lf_im(n)(i);

			int krow = 0;
			for(int i=0;i<basis::tet(log2p).tm;++i)
				for(int n=0;n<NV;++n)
					K(krow++,kcol) = (R(n)(i)-Rbar(n)(i))/dw;
			++kcol;
			
			uht(var)(mode) -= dw;
		}
	}	
	
	return;
}

void tet_hp::create_rsdl() {
	int indx,gindx,eind,find,iind,sgn,msgn,mode;
	
	int kn = NV*basis::tet(log2p).tm;
	Array<FLT,1> lclres(kn);
	Array<int,1> loc_to_glo(kn);
	
	for(int tind = 0; tind < ntet; ++tind){	
		
		int ind = 0;
		
		create_local_rsdl(tind, lclres);
		
		for (int m = 0; m < 4; ++m) {
			gindx = NV*tet(tind).pnt(m);
			for (int n = 0; n < NV; ++n)
				loc_to_glo(ind++) = gindx+n;
		}		
		
		/* EDGE MODES */
		if (basis::tet(log2p).p > 1) {
			for(int i = 0; i < 6; ++i) {
				eind = npnt*NV + tet(tind).seg(i)*basis::tet(log2p).em*NV;
				sgn = tet(tind).sgn(i);
				msgn = 1;
				for (int m = 0; m < basis::tet(log2p).em; ++m) {
					for(int n = 0; n < NV; ++n) {
						lclres(ind) *= msgn;
						loc_to_glo(ind++) = eind + m*NV + n;
					}
					msgn *= sgn;
				}
			}
		}
		
		/* FACE MODES */
		if (basis::tet(log2p).p > 2) {
			for(int i = 0; i < 4; ++i){
				sgn = -tet(tind).rot(i);
				find = npnt*NV+nseg*basis::tet(log2p).em*NV+tet(tind).tri(i)*basis::tet(log2p).em*NV;
				mode = 0;
				msgn = 1;		
				for(int m = 1; m <= basis::tet(log2p).em-1; ++m) {
					for(int j = 1; j <= basis::tet(log2p).em-m; ++j){
						for(int n = 0; n < NV; ++n){
							lclres(ind) *= msgn;							
							loc_to_glo(ind++) = find+mode*NV+n;
						}
						++mode;
					}
					msgn *= sgn;
				}
			}		
		}		
		
		/* INTERIOR MODES */
		iind = npnt*NV+nseg*basis::tet(log2p).em*NV+ntri*basis::tet(log2p).fm*NV+tind*basis::tet(log2p).im*NV;
		for(int m = 0; m < basis::tet(log2p).im; ++m) 
			for(int n = 0; n < NV; ++n)
				loc_to_glo(ind++) = iind+m*NV+n;
			
		
	
#ifdef petsc
		
		PetscErrorCode ierr;
		for(int i = 0; i < kn; ++i)
			ierr = VecSetValues(petsc_f,1,&loc_to_glo(i),&lclres(i),INSERT_VALUES);
		// CHKERRQ(ierr);

#endif
		
#ifndef petsc
		for(int i = 0; i < kn; ++i)
			res_vec(loc_to_glo(i))+=lclres(i);	
#endif
		
		
		
	}
	
	return;
}

void tet_hp::create_local_rsdl(int tind, Array<FLT,1> &lclres) {
	Array<TinyVector<FLT,MXTM>,1> lf_re(NV),lf_im(NV);

	int ind = 0;

	ugtouht(tind);
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tet(log2p).tm;++i)
		for(int n=0;n<NV;++n)
			lclres(ind++)=lf_re(n)(i)+lf_im(n)(i);	
	
	return;
}

#ifndef petsc
void tet_hp::vec_to_ug(){
	
	int ind = 0;
	
	for(int i = 0; i < npnt; ++i)
		for(int n = 0; n < NV; ++n)
			ug.v(i,n) = ug_vec(ind++);
	
	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tet(log2p).em; ++m)
			for(int n = 0; n < NV; ++n)
				ug.e(i,m,n) = ug_vec(ind++);
	
	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tet(log2p).fm; ++m)
			for(int n = 0; n < NV; ++n)
				ug.f(i,m,n) = ug_vec(ind++);		
	
	for(int i = 0; i < ntet; ++i)
		for(int m = 0; m < basis::tet(log2p).im; ++m)
			for(int n = 0; n < NV; ++n)
				ug.i(i,m,n) = ug_vec(ind++);	
	
	return;	
}

void tet_hp::ug_to_vec(){
	
	int ind = 0;
	
	for(int i = 0; i < npnt; ++i)
		for(int n = 0; n < NV; ++n)
			ug_vec(ind++) = ug.v(i,n);
	
	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tet(log2p).em; ++m)
			for(int n = 0; n < NV; ++n)
				ug_vec(ind++) = ug.e(i,m,n);
	
	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tet(log2p).fm; ++m)
			for(int n = 0; n < NV; ++n)
				ug_vec(ind++) = ug.f(i,m,n);;		
	
	for(int i = 0; i < ntet; ++i)
		for(int m = 0; m < basis::tet(log2p).im; ++m)
			for(int n = 0; n < NV; ++n)
				ug_vec(ind++) = ug.i(i,m,n);	
	
	return;	
}
#endif





//#ifdef row_index
///* row-indexed sparse storage */
//void tet_hp::insert_sparse(int row, int col, FLT value){
//	
//	int find_col = -1;
//	
//	/* insert diagonal element */
//	if(row == col){
//		sa(row) += value;
//		return;
//	}
//	
//	/* add value to already existing element in sparse */
//	for(int i = ija(row); i < ija(row+1); ++i){
//		if(ija(i) == col){
//			sa(i) += value;
//			return;
//		}
//	}
//	
//	/* if not diagonal entry or pre-exisiting element then add to list in correct location */
//	for(int i = ija(row); i < ija(row+1); ++i){
//		if(ija(i) > col){
//			find_col = i;
//			break;
//		}
//	}
//	/* no column bigger so insert at end */
//	if(find_col == -1)
//		find_col = ija(row+1);
//		
//	++number_sparse_elements;
//	
//	/* shift start of column index by one */
//	for(int i = row+1; i < size_sparse_matrix+1; ++i)
//		++ija(i);
//		
//	/* slide all non-diagonal entries after insertion column to the right (slow but only done once) */
//		for(int i = number_sparse_elements; i > find_col; --i){
//			sa(i)=sa(i-1);
//			ija(i)=ija(i-1);
//		}
//	
//	/* insert new element into sparse matrix */
//	sa(find_col) = value;
//	ija(find_col) = col;
//	
//	
//	return;	
//}
//
///* zero out sparse matrix but keep same structure */
//void tet_hp::zero_sparse(){
//	
//	for (int i = 0; i < number_sparse_elements; ++i){
//		res_vec(i) = 0.0;
//		sa(i) = 0.0;
//	}
//	sa(number_sparse_elements) = 0.0;
//	
//	return;
//}
//
//void tet_hp::initialize_sparse(){
//	
//	size_sparse_matrix = (npnt+nseg*em0+ntri*fm0+ntet*im0)*NV;
//	
//	/* sparse matrix allocation */
//	ija.resize(MXTM*NV*ntet);//too much storage resize later
//	sa.resize(MXTM*NV*ntet);
//	number_sparse_elements = size_sparse_matrix;
//	sa = 0.0;
//	/* creates sparse matrix with zeros on diagonal */
//	for(int i = 0; i < number_sparse_elements+1; ++i)
//		ija(i) = size_sparse_matrix+1;
//	return;
//}
//#endif


//void tet_hp::sparse_matrix_multiply(FLT *x, FLT *b){
//
//	for(int i = 0; i < size_sparse_matrix; ++i){
//		b(i) = sa(i)*x(i);
//		for(int j = ija(i); j < ija(i+1); ++j)
//			b(i) += sa(j)*x(ija(j));
//	}
//	
//	return;
//}
//
//
//void tet_hp::sparse_matrix_multiply_transpose(FLT *x, FLT *b){
//	
//	for(int i = 0; i < size_sparse_matrix; ++i)
//		b(i) = sa(i)*x(i);
//
//	for(int i = 0; i < size_sparse_matrix; ++i){
//		for(int j = ija(i); j < ija(i+1); ++j){
//			int k = ija(j)
//			b(k) += sa(j)*x(i);
//		}
//	}
//	
//	return;
//}