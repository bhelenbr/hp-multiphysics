/*
 *  matrix.cpp
 *  utilities
 *
 *  Created by Brian Helenbrook on 1/14/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */

#include "myblas.h"
#define BZ_DEBUG

using namespace blitz;


double spectral_radius(Array<double,2> A) {
	int n = A.extent(firstDim);
	Array<double,1> x(n),temp(n);
	double lambda;
	x=1;
	
	/* Power Iteration to find eigenvector for dominate eigenvalue */
	for(int m=0;m<20;++m){
		temp=0;
		for(int i=0; i<n; ++i){
			for(int j=0; j<n; ++j){
				temp(i)+=A(i,j)*x(j);
			}
		}
		x=temp/l2norm(temp);
	}
	
	/* Spectral radius using the Rayleigh quotient */
	lambda = inner_product(x, temp)/inner_product(x, x);
	
	return(fabs(lambda));
}

double l2norm(Array<double,1> x) {
	int n = x.extent(firstDim);
	double temp=0.0;
	
	for(int i=0;i<n;++i)
		temp+=x(i)*x(i);	
	
	return(pow(temp, 0.5));
}

double inner_product(Array<double,1> x, Array<double,1> y) {
	int n = x.extent(firstDim);
	double temp=0.0;
	
	for(int i=0;i<n;++i)
		temp+=x(i)*y(i);	
	
	return(temp);
}

void matrix_absolute_value(Array<double,2> &A) {
	int n = A.extent(firstDim);
	Array<double,1> lambda_real(n),lambda_imag(n);
	Array<double,2> VR(n,n),temp(n,n);
	int lwork=4*n,info,ipiv[n];
	char vchar[] = "V";
	char nchar[] = "N";
	char trans[] = "T";
	double work[lwork];
	
	/* find eigenvalues and eigenvectors */
	GEEV(nchar,vchar, n, A.data(), n, lambda_real.data(), lambda_imag.data(), temp.data(), n, VR.data(), n, work, lwork, info);
	
	if (info != 0) {
		std::cerr << "DGEEV FAILED FOR MATRIX ABSOLUTE VALUE IN UTILITIES" << std::endl;
		assert(0);
	}
	
	/* absolute value of eigenvalues */
	for (int i=0; i < n; ++i) 
		lambda_real(i) = sqrt(pow(lambda_real(i),2.0)+pow(lambda_imag(i),2.0));
	
	/* |diag(lambda)|V^T (transpose temp so that GETRF is compatible ) */
	for (int i = 0; i < n; ++i) 
		for (int j = 0; j < n;++j) 
			temp(j,i)=lambda_real(i)*VR(i,j);
	
	/*  LU factorization  */
	GETRF(n,  n, VR.data(), n, ipiv, info);
	
	if (info != 0) {
		std::cerr << "DGETRF FAILED FOR MATRIX ABSOLUTE VALUE IN UTILITIES" << std::endl;
		assert(0);
	}
	
	/* Solve transposed system temp = inv(VR)*temp */
	GETRS(trans,n,n,VR.data(),n,ipiv,temp.data(),n,info);
	
	if (info != 0) {
		std::cerr << "DGETRS FAILED FOR MATRIX ABSOLUTE VALUE IN UTILITIES" << std::endl;
		assert(0);
	}
	
	/* store transpose of temp in A */
	for (int i = 0; i < n; ++i) 
		for (int j = 0; j < n;++j) 
			A(i,j)=temp(j,i);
	
	return;
}

sparse_row_major::sparse_row_major(int nrow,const Array<int,1>& nnzero, int offset) {
	resize(nrow,nnzero,offset);
}

void sparse_row_major::resize(int nrow,const Array<int,1>& nnzero, int offset) {
	_nrow = nrow;
	_offset = offset;
	_cpt.resize(nrow+1);
	_cpt(0)=0;
	for (int i=1; i<nrow+1; ++i) 
		_cpt(i) = _cpt(i-1) +nnzero(i-1);
	int number_elements = _cpt(nrow);
	_val.resize(number_elements);
	_col.resize(number_elements);
	_col = INT_MAX-1;
	_cpt.reindexSelf(shape(offset));
}

void sparse_row_major::reset_columns() {
	_col = INT_MAX-1;
}



FLT& sparse_row_major::operator()(int row, int col) {
	int cindx;
	for (cindx=_cpt(row);col>_col(cindx);++cindx);
	
	if (_col(cindx) == col) {
		return(_val(cindx));
	}
	
#ifdef BZ_DEBUG
	if (cindx >= _cpt(row+1) || _col(_cpt(row+1)-1) != INT_MAX-1) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
	}
#endif
	
	/* slide all entries after insertion column to the right*/
	for(int i = _cpt(row+1)-1; i > cindx; --i) {
		_val(i)=_val(i-1);
		_col(i)=_col(i-1);
	}
	
	/* insert new element into sparse matrix */
	_col(cindx) = col;
	return(_val(cindx)=0.0);
}




void sparse_row_major::add_values(int row, int col, double val) {
	/* add value to already existing element in sparse */
	int cindx;
	for (cindx=_cpt(row);col>_col(cindx);++cindx);
	
	if (_col(cindx) == col) {
		_val(cindx) += val;
		return;
	}
	
#ifdef BZ_DEBUG
	if (cindx >= _cpt(row+1) || _col(_cpt(row+1)-1) != INT_MAX-1) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
	}
#endif
	
	/* slide all entries after insertion column to the right*/
	for(int i = _cpt(row+1)-1; i > cindx; --i) {
		_val(i)=_val(i-1);
		_col(i)=_col(i-1);
	}
	
	/* insert new element into sparse matrix */
	_val(cindx) = val;
	_col(cindx) = col;	
	
	return;	
}

void sparse_row_major::add_values(int nrows, const Array<int,1>& rows, int col, const Array<FLT,1>& D) {
	for (int indx = 0; indx < nrows; ++indx) {
		int row = rows(indx);
		add_values(row,col,D(indx));
	}
	return;	
}

void sparse_row_major::add_values(int row,int ncols, const Array<int,1>& cols,const Array<FLT,1>& D) {
	for (int indx = 0; indx < ncols; ++indx) {
		int col = cols(indx);
		add_values(row,col,D(indx));
	}
	return;
}

void sparse_row_major::add_values(int nels,const Array<int,1>& rows, const Array<int,1>& cols,const Array<FLT,1>& D) {
	for (int indx = 0; indx < nels; ++indx) {
		int row = rows(indx);
		int col = cols(indx);
		add_values(row,col,D(indx));
	}
	return;
}

void sparse_row_major::add_values(int nrows,const Array<int,1>& rows, int ncols, const Array<int,1>& cols,const Array<FLT,2>& M) {
	for (int lrow = 0; lrow < nrows; ++lrow) {
		int row = rows(lrow);
		for (int lcol = 0; lcol < ncols; ++lcol) {
			int col = cols(lcol);
			add_values(row,col,M(lrow,lcol));
		}
	}
	return;	
}


void sparse_row_major::set_values(int row, int col, double val) {
	/* add value to already existing element in sparse */
	int cindx;
	for (cindx=_cpt(row);col>_col(cindx);++cindx);
	
	if (_col(cindx) == col) {
		_val(cindx) = val;
		return;
	}
	
#ifdef BZ_DEBUG
	if (cindx >= _cpt(row+1) || _col(_cpt(row+1)-1) != INT_MAX-1) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
	}
#endif
	
	/* slide all entries after insertion column to the right*/
	for(int i = _cpt(row+1)-1; i > cindx; --i) {
		_val(i)=_val(i-1);
		_col(i)=_col(i-1);
	}
	
	/* insert new element into sparse matrix */
	_val(cindx) = val;
	_col(cindx) = col;	
	
	return;	
}

void sparse_row_major::set_values(int nrows, const Array<int,1>& rows, int col, const Array<FLT,1>& D) {
	for (int indx = 0; indx < nrows; ++indx) {
		int row = rows(indx);
		set_values(row,col,D(indx));
	}
	return;	
}

void sparse_row_major::set_values(int row,int ncols, const Array<int,1>& cols,const Array<FLT,1>& D) {
	for (int indx = 0; indx < ncols; ++indx) {
		int col = cols(indx);
		set_values(row,col,D(indx));
	}
	return;	
}


void sparse_row_major::set_values(int nels,const Array<int,1>& rows, const Array<int,1>& cols,const Array<FLT,1>& D) {
	for (int indx = 0; indx < nels; ++indx) {
		int row = rows(indx);
		int col = cols(indx);
		set_values(row,col,D(indx));
	}
	
	return;	
}

void sparse_row_major::set_values(int nrows,const Array<int,1>& rows, int ncols, const Array<int,1>& cols,const Array<FLT,2>& M) {
	for (int lrow = 0; lrow < nrows; ++lrow) {
		int row = rows(lrow);
		for (int lcol = 0; lcol < ncols; ++lcol) {
			int col = cols(lcol);
			set_values(row,col,M(lrow,lcol));
		}
	}
	return;	
}

void sparse_row_major::set_diag(int nels,const Array<int,1>& rows, FLT val, int offset) {	
	for (int indx = 0; indx < nels; ++indx) {
		int row = rows(indx);
		int col = row+offset;
		set_values(row,col,val);
	}
	return;	
}


void sparse_row_major::set_diag(int nels,const Array<int,1>& rows, const Array<FLT,1>& vals, int offset) {
	for (int indx = 0; indx < nels; ++indx) {
		int row = rows(indx);
		int col = row+offset;
		set_values(row,col,vals(indx));
	}
	return;
}

void sparse_row_major::zero_row(int row) {
		for(int j=_cpt(row);j<_cpt(row+1);++j)
			_val(j) = 0.0;
}


void sparse_row_major::zero_rows(int nrows,const Array<int,1>& rows) {
	for(int i=0;i<nrows;++i)
		for(int j=_cpt(rows(i));j<_cpt(rows(i)+1);++j)
			_val(j) = 0.0;
}

void sparse_row_major::check_for_unused_entries() {
	
	for(int i=_offset;i<_offset+_nrow;++i) {
		for (int j=_cpt(i);j<_cpt(i+1);++j) {
			if (_col(j) == INT_MAX-1) {
				std::cerr << "unused entry for row " << i << " allocated " << _cpt(i+1) -_cpt(i) << std::endl;
				Array<int,1> used(_col(Range(_cpt(i),j-1)));
				std::cerr << "used " << used << " unused " << _col(j) << std::endl;
				assert(0);
			}
		}
	}
	return;
}


ostream &operator<<(ostream &stream, sparse_row_major ob) {
	for (int i=ob._offset;i<ob._offset+ob._nrow;++i) {
		stream << i << " [";
		for (int j=ob._cpt(i);j<ob._cpt(i+1);++j)
			stream << '(' << ob._col(j) << ',' << ob._val(j) << ") ";
		stream << "]\n";
	}
  return stream;
}

void sparse_row_major::multiply_row(int row, FLT val) {
	for (int j=_cpt(row);j<_cpt(row+1);++j)
		_val(j) *= val;
	return;
}

void sparse_row_major::swap_rows(int row1, int row2) {
    
    const int nentries1 = _cpt(row1+1)-_cpt(row1);
    const int nentries2 = _cpt(row2+1)-_cpt(row2);
    if (nentries1 != nentries2) {
        std::cerr << "can't swap rows with different numbers of entries\n" << row1 << " allocated " <<   nentries1 << '\n' << row2 << " allocated " <<   nentries2 << std::endl;
        assert(0);
    }
    
    int j1 = _cpt(row1);
    int j2 = _cpt(row2);
    for (int j=0;j<nentries1;++j) {
        int col = _col(j2);
        _col(j2) = _col(j1);
        _col(j1) = col;
        
        FLT val = _val(j2);
        _val(j2) = _val(j1);
        _val(j1) = val;
        ++j1;
        ++j2;
    }
    
    return;
}

void sparse_row_major::combine_rows(int nrows, const Array<int, 1> &rows, int ncols, const Array<int, 1> &cols, const Array<double, 2> &M) {
	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	
	int row0 = cols(0);
	int nnz0 = _cpt(row0+1) -_cpt(row0);
	if (nnz0) {
		Array<FLT,2> _val_store(ncols,nnz0);
		_val_store(0,Range::all()) = _val(Range(_cpt(row0),_cpt(row0+1)-1));
		for(int i=1;i<ncols;++i) {
			int row = cols(i);
			int nnz = _cpt(row+1) -_cpt(row);
			if (nnz != nnz0) {
				std::cerr << "sparseness problem combine rows" << std::endl;
				assert(0);
			}
			
			int col0 = _cpt(row0);
			int col1 = _cpt(row);
			for(int col=0;col<nnz0;++col) {
				if (_col(col0++) != _col(col1++)) {
					std::cerr << "zeros indexing problem in combine rows" << std::endl;
					assert(0);
				}
				_val_store(i,Range::all()) = _val(Range(_cpt(row),_cpt(row+1)-1));
			}
		}
		
		Array<FLT,1> temp(nnz0);
		for(int i=0;i<nrows;++i) {
			temp = 0.0;
			for(int j=0;j<ncols;++j) {
				temp(Range::all()) += M(i,j)*_val_store(j,Range::all());
			}
			int row = rows(i);
			_val(Range(_cpt(row),_cpt(row+1)-1)) = temp(Range::all());
			_col(Range(_cpt(row),_cpt(row+1)-1)) = _col(Range(_cpt(row0),_cpt(row0+1)-1));
		}
	}
}

// Combine rows by multiply by multiplying by inverse of a DGETRF factorized square matrix
void sparse_row_major::combine_rows(int nrows, const Array<int, 1> &rows,const Array<double, 2> &A, int LDA,const Array<int, 1> &ipiv) {
	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	
	int row0 = rows(0);
	int nnz0 = _cpt(row0+1) -_cpt(row0);
	if (nnz0) {
		Array<FLT,2> _val_store(nnz0,nrows);
		_val_store(Range::all(),0) = _val(Range(_cpt(row0),_cpt(row0+1)-1));
		for(int i=1;i<nrows;++i) {
			int row = rows(i);
			int nnz = _cpt(row+1) -_cpt(row);
			if (nnz != nnz0) {
				std::cerr << "sparseness problem combine rows" << std::endl;
				assert(0);
			}
			
			int col0 = _cpt(row0);
			int col1 = _cpt(row);
			for(int col=0;col<nnz0;++col) {
				if (_col(col0++) != _col(col1++)) {
					std::cerr << "zeros indexing problem in combine rows" << std::endl;
					assert(0);
				}
			}
			_val_store(Range::all(),i) = _val(Range(_cpt(row),_cpt(row+1)-1));
		}
		
		int info;
		char trans[] = "T";
		GETRS(trans,nrows,nnz0,const_cast<FLT *>(A.data()),LDA,const_cast<int *>(ipiv.data()),_val_store.data(),nrows,info);
		
		for(int i=0;i<nrows;++i) {
			int row = rows(i);
			_val(Range(_cpt(row),_cpt(row+1)-1)) = _val_store(Range::all(),i);
		}
	}
}

void sparse_row_major::match_patterns(int row1, int row2) {
	
	int nnz1 = _cpt(row1+1) -_cpt(row1);
	int nnz2 = _cpt(row2+1) -_cpt(row2);
	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	if (nnz1 != nnz2) {
		std::cerr << "sparseness problem in match patterns" << nnz1 << ' ' << nnz2 << std::endl;
		assert(0);
	}
	
	int cpt1 = _cpt(row1);
	int cpt2 = _cpt(row2);
	for(int col=0;col<nnz1;++col) {
		/* match rows */
		if (_col(cpt1) < _col(cpt2)) {
			set_values(row2, _col(cpt1), 0.0);
		}
		else if (_col(cpt1) > _col(cpt2)) {
			set_values(row1,_col(cpt2),0.0);
		}
		++cpt1;
		++cpt2;
	}
}

void sparse_row_major::match_patterns(int nrows,const Array<int,1>& rows) {
	int nnz = _cpt(rows(0)+1) -_cpt(rows(0));
	Array<int,1> cols(nnz);
	cols = INT_MAX-1;
	
	if (!nnz) return;
	
	cols = _col(Range(_cpt(rows(0)),_cpt(rows(0)+1)-1));
	for(int row=1;row<nrows;++row) {
		int nnz1 = _cpt(rows(row)+1) -_cpt(rows(row));

		/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
		if (nnz1 != nnz) {
			std::cerr << "sparseness problem in match patterns" << nnz << ' ' << nnz1 << std::endl;
			assert(0);
		}
		
		for(int col_cnt=0;col_cnt<nnz;++col_cnt) {
			int col = _col(_cpt(rows(row))+col_cnt);
			if (col == INT_MAX-1) break;
		
			/* Find insertion spot */
			int cindx;
			for (cindx=0;col>cols(cindx);++cindx);
		
			if (col == cols(cindx)) continue;
			
#ifdef BZ_DEBUG
			if (cindx >= nnz || cols(nnz-1) != INT_MAX-1) {
				std::cerr << "Too many columns in match_patterns " << rows(row) << ' ' << col << std::endl;
				for(int row2=0;row2<nrows;++row2) {
					std::cerr << rows(row2) << ' ' << _col(Range(_cpt(rows(row2)),_cpt(rows(row2)+1)-1)) << std::endl;
				}
				std::cerr << " allocated entries " << nnz << std::endl;
				std::cerr << "Current usage " << cols << std::endl;
				assert(0);
			}
#endif
		
			/* slide all entries after insertion column to the right*/
			for(int i = nnz-1; i > cindx; --i) {
				cols(i)=cols(i-1);
			}

			cols(cindx) = col;
		}
	}
	
	// Insert 0's if entry doesn't exist
	for(int row=0;row<nrows;++row) {
		for(int col=0;col<nnz;++col) {
			add_values(rows(row),cols(col),0.0);
		}
	}

	return;
}


void sparse_row_major::add_rows(int row1, int row2) {
	
    const int nentries1 = _cpt(row1+1)-_cpt(row1);
    const int nentries2 = _cpt(row2+1)-_cpt(row2);
    if (nentries1 != nentries2) {
        std::cerr << "can't swap rows with different numbers of entries\n" << row1 << " allocated " <<   nentries1 << '\n' << row2 << " allocated " <<   nentries2 << std::endl;
        assert(0);
    }
    
    int j1 = _cpt(row1);
    int j2 = _cpt(row2);
    int col_check_sum = 0;
    for (int j=0;j<nentries1;++j) {
        col_check_sum += abs(_col(j1)-_col(j2));
         _val(j1) += _val(j2);
        ++j1;
        ++j2;
    }
    
    if (col_check_sum) {
        std::cerr << "can't add rows with different sparseness patter\n" << row1 << " allocated " <<   nentries1 << '\n' << row2 << " allocated " <<   nentries2 << std::endl;
        assert(0);
    }
    
    return;
}


void sparse_row_major::output_row(ostream &stream,int row) {
	for (int j=_cpt(row);j<_cpt(row+1);++j)
		stream << '(' << _col(j) << ',' << _val(j) << ") ";
	stream << std::endl;
	return;
}

void sparse_row_major::mmult(const Array<FLT,1>& x,Array<FLT,1> ax) {
	
	for (int i=_offset;i<_offset+_nrow;++i) {
		FLT lclax = 0.0;
		for (int j=_cpt(i);j<_cpt(i+1);++j)
			lclax += _val(j)*x(_col(j));
		ax(i-_offset) = lclax;
	}
	return;
}

void sparse_row_major::unpack(Array<FLT,2> tgt) {
    
    tgt = 0.0;
    for (int i=_offset;i<_offset+_nrow;++i) {
        FLT lclax = 0.0;
        for (int j=_cpt(i);j<_cpt(i+1);++j)
            tgt(i,_col(j)) = _val(j);
    }
    return;
}









