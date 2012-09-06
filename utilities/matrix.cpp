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
	GETRF(n, n, VR.data(), n, ipiv, info);
	
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

sparse_row_major::sparse_row_major(int nrow, Array<int,1>& nnzero, int offset) {
	resize(nrow,nnzero,offset);
}

void sparse_row_major::resize(int nrow, Array<int,1>& nnzero, int offset) {
	_nrow = nrow;
	_offset = offset;
	_cpt.resize(nrow+1);
	_cpt(0)=0;
	for (int i=1; i<nrow+1; ++i) 
		_cpt(i) = _cpt(i-1) +nnzero(i-1);
	int number_elements = _cpt(nrow);
	_val.resize(number_elements);
	_col.resize(number_elements);
	_col = INT_MAX;
	_cpt.reindexSelf(shape(offset));
}

void sparse_row_major::reset_columns() {
	_col = INT_MAX;
}



FLT& sparse_row_major::operator()(int row, int col) {
	int cindx;
	for (cindx=_cpt(row);col>_col(cindx);++cindx);
	
	if (_col(cindx) == col) {
		return(_val(cindx));
	}
	
#ifdef BZ_DEBUG
	if (cindx >= _cpt(row+1)) {
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
	if (cindx >= _cpt(row+1)) {
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
		
		/* add value to already existing element in sparse */
		int cindx;
		for (cindx = _cpt(row); col > _col(cindx); ++cindx);
		
		if (_col(cindx) == col) {
			_val(cindx) += D(indx);
			continue;
		}
		
#ifdef BZ_DEBUG
		if (cindx >= _cpt(row+1)) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
		}
#endif
		
		/* slide all entries after insertion column to the right*/
		for(int i = _cpt(row+1)-1; i > cindx; --i){
			_val(i)=_val(i-1);
			_col(i)=_col(i-1);
		}
		
		/* insert new element into sparse matrix */
		_val(cindx) = D(indx);
		_col(cindx) = col;	
	}
	
	return;	
}

void sparse_row_major::add_values(int row,int ncols, const Array<int,1>& cols,const Array<FLT,1>& D) {
	for (int indx = 0; indx < ncols; ++indx) {
		int col = cols(indx);
		
		/* add value to already existing element in sparse */
		int cindx;
		for (cindx = _cpt(row); col > _col(cindx); ++cindx);
		
		if (_col(cindx) == col) {
			_val(cindx) += D(indx);
			continue;
		}
		
#ifdef BZ_DEBUG
		if (cindx >= _cpt(row+1)) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
		}
#endif
		
		/* slide all entries after insertion column to the right*/
		for(int i = _cpt(row+1)-1; i > cindx; --i){
			_val(i)=_val(i-1);
			_col(i)=_col(i-1);
		}
		
		/* insert new element into sparse matrix */
		_val(cindx) = D(indx);
		_col(cindx) = col;	
	}
	
	return;	
}


void sparse_row_major::add_values(int nels,const Array<int,1>& rows, const Array<int,1>& cols,const Array<FLT,1>& D) {
	for (int indx = 0; indx < nels; ++indx) {
		int row = rows(indx);
		int col = cols(indx);
		
		/* add value to already existing element in sparse */
		int cindx;
		for (cindx = _cpt(row); col > _col(cindx); ++cindx);
		
		if (_col(cindx) == col) {
			_val(cindx) += D(indx);
			continue;
		}
		
#ifdef BZ_DEBUG
		if (cindx >= _cpt(row+1)) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
		}
#endif
		
		/* slide all entries after insertion column to the right*/
		for(int i = _cpt(row+1)-1; i > cindx; --i){
			_val(i)=_val(i-1);
			_col(i)=_col(i-1);
		}
		
		/* insert new element into sparse matrix */
		_val(cindx) = D(indx);
		_col(cindx) = col;	
	}
	
	return;	
}

void sparse_row_major::add_values(int nrows,const Array<int,1>& rows, int ncols, const Array<int,1>& cols,const Array<FLT,2>& M) {
	for (int lrow = 0; lrow < nrows; ++lrow) {
		int row = rows(lrow);
		for (int lcol = 0; lcol < ncols; ++lcol) {
			int col = cols(lcol);
			
			/* add value to already existing element in sparse */
			int cindx;
			for (cindx=_cpt(row);col>_col(cindx);++cindx);
			
			if (_col(cindx) == col) {
				_val(cindx) += M(lrow,lcol);
				continue;
			}
			
			
			
#ifdef BZ_DEBUG
			if (cindx >= _cpt(row+1)) {
				std::cerr << "Too many entries for row " << row << " and col " << col ;
				std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
				Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
				std::cerr << "Current usage " << colinds << std::endl;
				assert(0);
			}
#endif
			
			/* slide all entries after insertion column to the right*/
			for(int i = _cpt(row+1)-1; i > cindx; --i){
				_val(i)=_val(i-1);
				_col(i)=_col(i-1);
			}
			
			/* insert new element into sparse matrix */
			_val(cindx) = M(lrow,lcol);
			_col(cindx) = col;	
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
	if (cindx >= _cpt(row+1)) {
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
		
		/* add value to already existing element in sparse */
		int cindx;
		for (cindx = _cpt(row); col > _col(cindx); ++cindx);
		
		if (_col(cindx) == col) {
			_val(cindx) = D(indx);
			continue;
		}
		
#ifdef BZ_DEBUG
		if (cindx >= _cpt(row+1)) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;			assert(0);
		}
#endif
		
		/* slide all entries after insertion column to the right*/
		for(int i = _cpt(row+1)-1; i > cindx; --i){
			_val(i)=_val(i-1);
			_col(i)=_col(i-1);
		}
		
		/* insert new element into sparse matrix */
		_val(cindx) = D(indx);
		_col(cindx) = col;	
	}
	
	return;	
}

void sparse_row_major::set_values(int row,int ncols, const Array<int,1>& cols,const Array<FLT,1>& D) {
	for (int indx = 0; indx < ncols; ++indx) {
		int col = cols(indx);
		
		/* add value to already existing element in sparse */
		int cindx;
		for (cindx = _cpt(row); col > _col(cindx); ++cindx);
		
		if (_col(cindx) == col) {
			_val(cindx) = D(indx);
			continue;
		}
		
#ifdef BZ_DEBUG
		if (cindx >= _cpt(row+1)) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
		}
#endif
		
		/* slide all entries after insertion column to the right*/
		for(int i = _cpt(row+1)-1; i > cindx; --i){
			_val(i)=_val(i-1);
			_col(i)=_col(i-1);
		}
		
		/* insert new element into sparse matrix */
		_val(cindx) = D(indx);
		_col(cindx) = col;	
	}
	
	return;	
}


void sparse_row_major::set_values(int nels,const Array<int,1>& rows, const Array<int,1>& cols,const Array<FLT,1>& D) {
	for (int indx = 0; indx < nels; ++indx) {
		int row = rows(indx);
		int col = cols(indx);
		
		/* add value to already existing element in sparse */
		int cindx;
		for (cindx = _cpt(row); col > _col(cindx); ++cindx);
		
		if (_col(cindx) == col) {
			_val(cindx) = D(indx);
			continue;
		}
		
#ifdef BZ_DEBUG
		if (cindx >= _cpt(row+1)) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
		}
#endif
		
		/* slide all entries after insertion column to the right*/
		for(int i = _cpt(row+1)-1; i > cindx; --i){
			_val(i)=_val(i-1);
			_col(i)=_col(i-1);
		}
		
		/* insert new element into sparse matrix */
		_val(cindx) = D(indx);
		_col(cindx) = col;	
	}
	
	return;	
}

void sparse_row_major::set_values(int nrows,const Array<int,1>& rows, int ncols, const Array<int,1>& cols,const Array<FLT,2>& M) {
	for (int lrow = 0; lrow < nrows; ++lrow) {
		int row = rows(lrow);
		for (int lcol = 0; lcol < ncols; ++lcol) {
			int col = cols(lcol);
			
			/* add value to already existing element in sparse */
			int cindx;
			for (cindx=_cpt(row);col>_col(cindx);++cindx);
			
			if (_col(cindx) == col) {
				_val(cindx) = M(lrow,lcol);
				continue;
			}
			
#ifdef BZ_DEBUG
			if (cindx >= _cpt(row+1)) {
				std::cerr << "Too many entries for row " << row << " and col " << col ;
				std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
				Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
				std::cerr << "Current usage " << colinds << std::endl;
				assert(0);
			}
#endif
			
			/* slide all entries after insertion column to the right*/
			for(int i = _cpt(row+1)-1; i > cindx; --i){
				_val(i)=_val(i-1);
				_col(i)=_col(i-1);
			}
			
			/* insert new element into sparse matrix */
			_val(cindx) = M(lrow,lcol);
			_col(cindx) = col;	
		}
	}
	
	return;	
}

void sparse_row_major::set_diag(int nels,const Array<int,1>& rows, FLT val, int offset) {	
	for (int indx = 0; indx < nels; ++indx) {
		int row = rows(indx);
		int col = row+offset;
		
		/* add value to already existing element in sparse */
		int cindx;
		for (cindx = _cpt(row); col > _col(cindx); ++cindx);
		
		if (_col(cindx) == col) {
			_val(cindx) = val;
			continue;
		}
		
#ifdef BZ_DEBUG
		if (cindx >= _cpt(row+1)) {
			std::cerr << "Too many entries for row " << row << " and col " << col ;
			std::cerr << " allocated entries " << _cpt(row+1) - _cpt(row) << std::endl;
			Array<int,1> colinds(_col(Range(_cpt(row),_cpt(row+1)-1)));
			std::cerr << "Current usage " << colinds << std::endl;
			assert(0);
		}
#endif
		
		/* slide all entries after insertion column to the right*/
		for(int i = _cpt(row+1)-1; i > cindx; --i){
			_val(i)=_val(i-1);
			_col(i)=_col(i-1);
		}
		
		/* insert new element into sparse matrix */
		_val(cindx) = val;
		_col(cindx) = col;	
	}
	
	return;	
}


void sparse_row_major::zero_rows(int nrows,const Array<int,1>& rows) {
	for(int i=0;i<nrows;++i)
		for(int j=_cpt(rows(i));j<_cpt(rows(i)+1);++j)
			_val(j) = 0.0;
}

void sparse_row_major::check_for_unused_entries() {
	
	for(int i=_offset;i<_offset+_nrow;++i) {
		for (int j=_cpt(i);j<_cpt(i+1);++j) {
			if (_col(j) == INT_MAX) {
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

void sparse_row_major::output_row(ostream &stream,int row) {
	for (int j=_cpt(row);j<_cpt(row+1);++j)
		stream << '(' << _col(j) << ',' << _val(j) << ") ";
	stream << std::endl;
	return;
}









