/*
 *  myblas.h
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#ifndef _myblas_h_
#define _myblas_h_

#ifdef F2CFortran
#include<cfortran.h>
#endif

extern "C" {
    void DPBTRSNU2(double *abd, int stride, int ordr, int ofdg, double *b, int rhs);
    void DPBTRSNU1(double *abd, int ordr, int ofdg, double *b, int rhs);
    void DPBTRSNU(double **abd, int ordr, int ofdg, double *b, int rhs);
    void DPBSLN(double **abd, int ordr, int ofdg, double *b, int rhs);
    void DGETLS(double *abd, int stride, int ordr, double *b);
    void DGETUS(double *abd, int stride, int ordr, double *b);
    void BLCKTRI(int nblck,double (*a)[2], double (*b)[2], double (*c)[2], double (*d)[2]);
    void BLCKTRI2(int nblck,double (*a)[2], double (*b)[2], double (*c)[2], double (*d)[2], double (*e)[2]);
    int recur(int n,int ipoly, double al, double be, double *a, double *b);
    int radau(int n,double *alpha,double *beta, double end,double *zero,double *weight,double *e, double *a, double *b);
    int lob(int n,double *alpha,double *beta, double left, double right,double *zero,double *weight,double *e, double *a, double *b);
    int gauss(int n, const double *alpha,const double *beta, double eps, double *zero, double *weight, double *e);
#ifndef F2CFortran
    extern int dgetrf_(const int *,const int *,double *,const int *,int *,int *);
    extern int dgetrs_(const char *,const int *,const int *,const double *,const int *,const int *,double *,const int *,int *);
    extern int dpbtrf_(const char *,const int *,const int *,double *,const int *,int *);
    extern int dpbtrs_(const char *,const int *,const int *,const int *,const double *,const int *,double *,const int *,int *);
    extern int dpbco_(double *,const int *,const int *,const int *,double *,double *,int *);
    extern int dpbsl_(const double *,const int *,const int *,const int *,double *);
    extern int dpoco_(double *,const int *,const int *,double *,double *,int *);
    extern int dposl_(const double *,const int *,const int *,double *);
    // extern int dsics_(int *,int *,int *,int *,double *,int *,int *,int *,int *,double *,double *,double *,int *);
    // extern int sllti2_(int *,double *,double *,int *,int *,int *,double *,double *);
    extern int dgels_(const char *,const int *,const int *,const int *,double *,const int *,double *,const int *,double *,const int *,int *);
    extern int dgeev_(const char *,const char *,const int *,double *,const int *,double *,double *,double *,const int *,double *,const int *,double *,const int *,int *);
    extern int dgecon_(const char *,const int *,const double *,const int *,const double *,double *,double *,int *,int *);
    extern int dpotrf_(const char *,const int *,double *,const int *,int *);
    extern int dpotrs_(const char *,const int *,const int *,const double *,const int *,double *,const int *,int *);
    extern int dspev_(const char *,const char *,const int *,double *,double *,double *,const int *,double *,int *);
#endif
}

#ifdef F2CFortran
PROTOCCALLSFSUB6(DGETRF,dgetrf,INT,INT,DOUBLEV,INT,INTV,PINT)
#define DGETRF(M,N,A,LDA,IPIV,INFO) CCALLSFSUB6(DGETRF,dgetrf,INT,INT,DOUBLEV,INT,INTV,PINT,M,N,A,LDA,IPIV,INFO)

PROTOCCALLSFSUB9(DGETRS,dgetrs,STRING,INT,INT,DOUBLEV,INT,INTV,DOUBLEV,INT,PINT)
#define DGETRS(trans,n,nrhs,a,lda,ipiv,b,ldb,info) CCALLSFSUB9(DGETRS,dgetrs,STRING,INT,INT,DOUBLEV,INT,INTV,DOUBLEV,INT,PINT,trans,n,nrhs,a,lda,ipiv,b,ldb,info)

PROTOCCALLSFSUB6(DPBTRF,dpbtrf,STRING,INT,INT,DOUBLEV,INT,PINT)
#define DPBTRF(UPLO,N,KD,AB,LDAB,INFO) CCALLSFSUB6(DPBTRF,dpbtrf,STRING,INT,INT,DOUBLEV,INT,PINT,UPLO,N,KD,AB,LDAB,INFO)

PROTOCCALLSFSUB9(DPBTRS,dpbtrs,STRING,INT,INT,INT,DOUBLEV,INT,DOUBLEV,INT,PINT)
#define DPBTRS(UPLO,N,KD,NRHS,AB,LDAB,B,LDB,INFO) CCALLSFSUB9(DPBTRS,dpbtrs,STRING,INT,INT,INT,DOUBLEV,INT,DOUBLEV,INT,PINT,UPLO,N,KD,NRHS,AB,LDAB,B,LDB,INFO)

PROTOCCALLSFSUB7(DPBCO,dpbco,DOUBLEV,INT,INT,INT,PDOUBLE,DOUBLEV,PINT)
#define DPBCO(ABD,LDA,N,M,RCOND,Z,INFO) CCALLSFSUB7(DPBCO,dpbco,DOUBLEV,INT,INT,INT,PDOUBLE,DOUBLEV,PINT,ABD,LDA,N,M,RCOND,Z,INFO)

PROTOCCALLSFSUB5(DPBSL,dpbsl,DOUBLEV,INT,INT,INT,DOUBLEV)
#define DPBSL(ABD,LDA,N,M,B) CCALLSFSUB5(DPBSL,dpbsl,DOUBLEV,INT,INT,INT,DOUBLEV,ABD,LDA,N,M,B)

PROTOCCALLSFSUB6(DPOCO,dpoco,DOUBLEV,INT,INT,PDOUBLE,DOUBLEV,PINT)
#define DPOCO(ABD,LDA,N,RCOND,Z,INFO) CCALLSFSUB6(DPOCO,dpoco,DOUBLEV,INT,INT,PDOUBLE,DOUBLEV,PINT,ABD,LDA,N,RCOND,Z,INFO)

PROTOCCALLSFSUB4(DPOSL,dposl,DOUBLEV,INT,INT,DOUBLEV)
#define DPOSL(ABD,LDA,N,B) CCALLSFSUB4(DPOSL,dposl,DOUBLEV,INT,INT,DOUBLEV,ABD,LDA,N,B)

PROTOCCALLSFSUB13(DSICS,dsics,INT,INT,INTV,INTV,DOUBLEV,INT,PINT,INTV,INTV,DOUBLEV,DOUBLEV,DOUBLEV,PINT)
#define DSICS(N,NELT,IA,JA,A,ISYM,NEL,IEL,JEL,EL,D,R,IWARN) CCALLSFSUB13(DSICS,dsics,INT,INT,INTV,INTV,DOUBLEV,INT,PINT,INTV,INTV,DOUBLEV,DOUBLEV,DOUBLEV,PINT,N,NELT,IA,JA,A,ISYM,NEL,IEL,JEL,EL,D,R,IWARN)

PROTOCCALLSFSUB8(SLLTI2,sllti2,INT,DOUBLEV,DOUBLEV,INT,INTV,INTV,DOUBLEV,DOUBLEV)
#define SLLTI2( N, B, X, NEL, IEL, JEL, EL, DINV ) CCALLSFSUB8(SLLTI2,sllti2,INT,DOUBLEV,DOUBLEV,INT,INTV,INTV,DOUBLEV,DOUBLEV,N, B, X, NEL, IEL, JEL, EL, DINV)

PROTOCCALLSFSUB11(DGELS,dgels,STRING,INT,INT,INT,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,INT,PINT)
#define DGELS(trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info) CCALLSFSUB11(DGELS,dgels,STRING,INT,INT,INT,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,INT,PINT,trans,m,n,nrhs,a,lda,b,ldb,work,lwork,info)

PROTOCCALLSFSUB14(DGEEV,dgeev,STRING,STRING,INT,DOUBLEV,INT,DOUBLEV,DOUBLEV,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,INT,PINT)
#define DGEEV(JOBVL,JOBVR,N,A,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,INFO) CCALLSFSUB14(DGEEV,dgeev,STRING,STRING,INT,DOUBLEV,INT,DOUBLEV,DOUBLEV,DOUBLEV,INT,DOUBLEV,INT,DOUBLEV,INT,PINT,JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO)

PROTOCCALLSFSUB9(DGECON,dgecon,STRING,INT,DOUBLEV,INT,DOUBLE,PDOUBLE,DOUBLEV,INTV,PINT)
#define DGECON(norm, n, a, lda, anorm, rcond, work, iwork, info) CCALLSFSUB9(DGECON,dgecon,STRING,INT,DOUBLEV,INT,DOUBLE,PDOUBLE,DOUBLEV,INTV,PINT,norm, n, a, lda, anorm, rcond, work, iwork, info)

PROTOCCALLSFSUB5(DPOTRF,dpotrf,STRING,INT,DOUBLEV,INT,PINT)
#define DPOTRF(UPLO,N,A,LDA,INFO) CCALLSFSUB5(DPOTRF,dpotrf,STRING,INT,DOUBLEV,INT,PINT,UPLO,N,A,LDA,INFO)

PROTOCCALLSFSUB8(DPOTRS,dpotrs,STRING,INT,INT,DOUBLEV,INT,DOUBLEV,INT,PINT)
#define DPOTRS(UPLO,N,NRHS,A,LDA,B,LDB,INFO) CCALLSFSUB8(DPOTRS,dpotrs,STRING,INT,INT,DOUBLEV,INT,DOUBLEV,INT,PINT,UPLO,N,NRHS,A,LDA,B,LDB,INFO)

PROTOCCALLSFSUB9(DSPEV,dspev,STRING,STRING,INT,DOUBLEV,DOUBLEV,DOUBLEV,INT,DOUBLEV,PINT)
#define DSPEV(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO) CCALLSFSUB9(DSPEV,dspev,STRING,STRING,INT,DOUBLEV,DOUBLEV,DOUBLEV,INT,DOUBLEV,PINT,JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, INFO)
#endif


#ifdef SINGLE
#define GETRF SGETRF
#define GETRS SGETRS
#define PBTRF SPBTRF
#define PBTRS SPBTRS
#define GEEV SGEEV
#define EPSILON FLT_EPSILON
#ifndef FLT
#define FLT float
#endif
#else
#define GETRF DGETRF
#define GETRS DGETRS
#define PBTRF DPBTRF
#define PBTRS DPBTRS
#define GEEV DGEEV
#define EPSILON DBL_EPSILON
#ifndef FLT
#define FLT double
#endif
#endif

#include<blitz/array.h>
double spectral_radius(blitz::Array<double,2> A);
double l2norm(blitz::Array<double,1> x);
double inner_product(blitz::Array<double,1> x, blitz::Array<double,1> y);
void matrix_absolute_value(blitz::Array<double,2> &A);

using namespace blitz;

class sparse_row_major {
public:
	int _nrow, _offset;
	Array<int,1> _cpt; //pointer to column indices for each row
	Array<int,1> _col; //sparse list of column indices
	Array<FLT,1> _val; //sparse list of matrix values
	
	sparse_row_major() {}
	sparse_row_major(int nrow,const Array<int,1>& nnzero, int offset = 0);
	void resize(int nrow,const Array<int,1>& nnzero, int offset=0);
	void reset_columns();
	
	void add_values(int row, int col, FLT value);
	void add_values(int nrows, const Array<int,1>& rows, int col, const Array<FLT,1>& D);
	void add_values(int row, int ncols, const Array<int,1>& cols, const Array<FLT,1>& D);
	void add_values(int nels, const Array<int,1>& rows, const Array<int,1>& cols, const Array<FLT,1>& D);
	void add_values(int nrows,const Array<int,1>& rows, int ncols, const Array<int,1>& cols,const Array<FLT,2>& M);
	
	void set_values(int row, int col, FLT value);
	void set_values(int nrows, const Array<int,1>& rows, int col, const Array<FLT,1>& D);
	void set_values(int row, int ncols, const Array<int,1>& cols, const Array<FLT,1>& D);
	void set_values(int nels, const Array<int,1>& rows, const Array<int,1>& cols, const Array<FLT,1>& D);
	void set_values(int nrows,const Array<int,1>& rows, int ncols, const Array<int,1>& cols,const Array<FLT,2>& M);
	
	void set_diag(int nels,const Array<int,1>& rows, FLT val, int offset=0);
	void set_diag(int nels,const Array<int,1>& rows, const Array<FLT,1>& vals, int offset=0);
	void zero_row(int row);
	void zero_rows(int nrows,const Array<int,1>& rows);
	void reset_row(int row);
	void reset_rows(int nrows,const Array<int,1>& rows);
	void multiply_row(int row, FLT val);
	void match_patterns(int row1, int row2);
	void match_patterns(int nrows,const Array<int,1>& rows);
	void swap_rows(int row1, int row2);  // Rows must have same sparseness pattern
	void add_rows(int row1, int row2); // Rows must have same sparseness pattern
	void combine_rows(int nrows,const Array<int,1>& rows, int ncols, const Array<int,1>& cols,const Array<FLT,2>& M); // Rows must have same sparseness pattern
	void combine_rows(int nrows, const Array<int, 1> &rows,const Array<double, 2> &A, int LDA,const Array<int, 1> &ipiv); 	// Combine rows by multiply by multiplying by inverse of a DGETRF factorized square matrix
	void mmult(const Array<FLT,1>& x,Array<FLT,1> rslt);
	void unpack(Array<FLT,2> tgt);
	
	// Array<FLT,1>& operator*(Array<FLT,1>&);
	
	FLT& operator()(int row, int col);
    int nrows() {return(_nrow);}
    int nentries_for_row(int row);
    void get_value_and_col(int row, int entry, FLT &value, int& col);
	void check_for_unused_entries(std::ostream &log);
	friend ostream &operator<<(ostream &stream, sparse_row_major mat);
	void output_row(ostream &stream,int row);
	
};
#endif

