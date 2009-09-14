/*
 *  sparse.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 9/14/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#include "tet_hp.h"

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
		
	/* slide all non-diagonal entries after insertion column to the right*/
		for(int i = number_sparse_elements; i > find_col; --i){
			sa(i)=sa(i-1);
			ija(i)=ija(i-1);
		}
	
	/* insert new element into sparse matrix */
	ija(find_col) = col;
	sa(find_col) = value;
	
	
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