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

void tet_hp::insert_sparse(int row, int col, FLT value){
	
	int find_col = -1;
	
	/* insert diagonal element */
	if(row == col){
		sa(row) += value;
		return;
	}
	/* add value to already existing element in sparse */
	for(int i = ija(row); i < ija(row+1); ++i){
		if(ija(i) == col){
			sa(i) += value;
			return;
		}
	}
	
	/* if not diagonal entry or pre-exisiting element then add to list in correct location */
	for(int i = ija(row); i < ija(row+1); ++i){
		if(ija(i) > col){
			find_col = i;
			break;
		}
	}
	/* no column bigger so insert at end */
	if(find_col == -1)
		find_col = ija(row+1);
		
	++number_sparse_elements;
	
	/* shift start of column index by one */
	for(int i = row+1; i < size_sparse_matrix+1; ++i)
		++ija(i);
		
	/* slide all non-diagonal entries after insertion column to the right (slow but only done once) */
		for(int i = number_sparse_elements; i > find_col; --i){
			sa(i)=sa(i-1);
			ija(i)=ija(i-1);
		}
	
	/* insert new element into sparse matrix */
	ija(find_col) = col;
	sa(find_col) = value;
	
	
	return;	
}


void tet_hp::create_jacobian() {
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
		
		for(int i = 0; i < kn; ++i)
			for(int j = 0; j < kn; ++j)
				insert_sparse(loc_to_glo(i), loc_to_glo(j), K(i,j));				
				
				
	}
	
	ija.resizeAndPreserve(number_sparse_elements+1);
	sa.resizeAndPreserve(number_sparse_elements+1);

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