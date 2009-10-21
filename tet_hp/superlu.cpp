/*
 *  superlu.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 10/9/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */


/*
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 *
 */

#include "tet_hp.h"
#include "slu_ddefs.h"
#include "hp_boundary.h"

void tet_hp::superlu(){
	cout << "superlu called" <<endl;
	/*
	 * Purpose
	 * =======
	 *
	 * The driver program DLINSOLX2.
	 * 
	 * This example illustrates how to use DGSSVX to solve systems repeatedly
	 * with the same sparsity pattern of matrix A.
	 * In this case, the column permutation vector perm_c is computed once.
	 * The following data structures will be reused in the subsequent call to
	 * DGSSVX: perm_c, etree
	 * 
	 */
	
//    /* advanced driver dgssvx allocation */ 
//    char           equed[1];
//    SuperMatrix    A, L, U;
//    SuperMatrix    B, X;
//	int            *perm_r; /* row permutations from partial pivoting */
//    int            *perm_c; /* column permutation vector */
//    int            *etree;
//    void           *work;
//    int            info, lwork;
//    double         *R, *C;
//    double         *ferr, *berr;
//    double         rpg, rcond;
//    mem_usage_t    mem_usage;
//    superlu_options_t options;
//    SuperLUStat_t  stat;
//	lwork = 0; // superlu determines memory allocation
	
	
	/* simple driver dgssv allocation */
    SuperMatrix    A, L, U;
    SuperMatrix    B;
	  int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            info;
    superlu_options_t options;
    SuperLUStat_t stat;	
	
	FLT newton_norm,tol=1.0e-14;
	int max_newton_its = 50;
		
    /* Set the default input options:
	 options.Fact = DOFACT;
	 options.Equil = YES;
	 options.ColPerm = COLAMD;
	 options.DiagPivotThresh = 1.0;
	 options.Trans = NOTRANS;
	 options.IterRefine = NOREFINE;
	 options.SymmetricMode = NO;
	 options.PivotGrowth = NO;
	 options.ConditionNumber = NO;
	 options.PrintStat = YES;
     */
    
	set_default_options(&options);	
    options.DiagPivotThresh = 1.0;// .1 for diagonally dominant I think but need to run in symmetric mode
	//options.PrintStat = NO;
	
	/* send global solution to ug_vec */
	ug_to_vec();
	
	bool compressed_column = true;
	if ( !(perm_c = intMalloc(size_sparse_matrix)) ) ABORT("Malloc fails for perm_c[].");

	
	for(int i = 0; i < max_newton_its; ++i) {

		if ( !(perm_r = intMalloc(size_sparse_matrix)) ) ABORT("Malloc fails for perm_r[].");


		/* zero out sparse and residual but keep sparsity pattern */
		zero_sparse();
		
//		dCreate_Dense_Matrix(&X, size_sparse_matrix, 1, res_vec.data(), size_sparse_matrix, SLU_DN, SLU_D, SLU_GE);

		/* create jacobian */
		create_jacobian(compressed_column);

		/* apply neumman bc's */
		apply_neumman(compressed_column);
		
		/* create residual */
		create_rsdl();	
		
		/* apply dirichlet boundary conditions to sparse matrix and vector */
		for(int j = 0; j < nfbd; ++j)
			hp_fbdry(j)->apply_sparse_dirichlet(compressed_column);
		
		FLT max_resid = 0.0;
		for (int j = 0; j < size_sparse_matrix; ++j){
			if(fabs(res_vec(j)) > max_resid)
				max_resid = fabs(res_vec(j));
		}
		cout << "L infinity norm of residual " << max_resid << endl;
		
		if (max_resid < tol) break;		

		/* create super matrix A using compressed column storage */
		dCreate_CompCol_Matrix(&A, size_sparse_matrix, size_sparse_matrix, number_sparse_elements, sparse_val.data(), sparse_ind.data(), sparse_ptr.data(), SLU_NC, SLU_D, SLU_GE);
		//dPrint_CompCol_Matrix("A", &A);

//		if ( !(etree = intMalloc(size_sparse_matrix)) ) ABORT("Malloc fails for etree[].");
//		if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) ) 
//			ABORT("SUPERLU_MALLOC fails for R[].");
//		if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
//			ABORT("SUPERLU_MALLOC fails for C[].");
//		if ( !(ferr = (double *) SUPERLU_MALLOC(1 * sizeof(double))) )
//			ABORT("SUPERLU_MALLOC fails for ferr[].");
//		if ( !(berr = (double *) SUPERLU_MALLOC(1 * sizeof(double))) ) 
//			ABORT("SUPERLU_MALLOC fails for berr[].");
		
		
		dCreate_Dense_Matrix(&B, size_sparse_matrix, 1, res_vec.data(), size_sparse_matrix, SLU_DN, SLU_D, SLU_GE);
		//dPrint_Dense_Matrix("B", &B);
		

		/* Initialize the statistics variables. */
		StatInit(&stat);

//		/* advanced superlu driver */
//		dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,&mem_usage, &stat, &info);		
//		
//		/* get access to solution and store in du */
//		double *du = (double*) ((DNformat*) X.Store)->nzval;
			

		/* solve system Ax=b and stores solution x in b */
		dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
		/* get access to solution and store in du */
		double *du = (double*) ((DNformat*) B.Store)->nzval;

		if (info > 0) {
			cout << "something went wrong with superlu \n";
			exit(1);
		}

		for(int j = 0; j < size_sparse_matrix; ++j)
			ug_vec(j)-=du[j];
		
		/* send ug_vec to global solution */
		vec_to_ug();
		
		if ( options.PrintStat ) StatPrint(&stat);
		StatFree(&stat);
		
//		newton_norm = 0.0;		
//		for (int j = 0; j < size_sparse_matrix; ++j)
//			newton_norm += du[j]*du[j];
//
//		newton_norm=sqrt(newton_norm);
//		cout << "L2 norm of du in newtons method " << newton_norm << endl;
		
//		if (newton_norm < tol) {
//			cout << "total newton iterations " << i+1 << endl;
//			break;
//		}

		/* once factorized can keep same pattern for next iteration */
//		options.Fact = SamePattern;
		SUPERLU_FREE (perm_r);
		//SUPERLU_FREE (xact);
//		SUPERLU_FREE (etree);
//		SUPERLU_FREE (R);
//		SUPERLU_FREE (C);
//		SUPERLU_FREE (ferr);
//		SUPERLU_FREE (berr);


	}

	SUPERLU_FREE (perm_c);

	return;
}

//void tet_hp::ilu_gmres(){
//	cout << "ilu_gmres called" <<endl;
//	/*
//	 * Purpose
//	 * =======
//	 *
//	 * The driver program DLINSOLX2.
//	 * 
//	 * This example illustrates how to use DGSSVX to solve systems repeatedly
//	 * with the same sparsity pattern of matrix A.
//	 * In this case, the column permutation vector perm_c is computed once.
//	 * The following data structures will be reused in the subsequent call to
//	 * DGSSVX: perm_c, etree
//	 * 
//	 */
//	
//    /* advanced driver dgssvx allocation */ 
//    char           equed[1];
//    SuperMatrix    A, L, U;
//    SuperMatrix    B, X;
//	int            *perm_r; /* row permutations from partial pivoting */
//    int            *perm_c; /* column permutation vector */
//    int            *etree;
//    void           *work;
//    int            info, lwork;
//    double         *R, *C;
//    double         *ferr, *berr;
//    double         rpg, rcond;
//    mem_usage_t    mem_usage;
//    superlu_options_t options;
//    SuperLUStat_t  stat;
//	lwork = 0; // superlu determines memory allocation
//	int restrt, iter, maxit;
//    double resid;
//    double *x, *b;
//	
//	//	/* simple driver dgssv allocation */
//	//    SuperMatrix    A, L, U;
//	//    SuperMatrix    B;
//	//	  int            *perm_r; /* row permutations from partial pivoting */
//	//    int            *perm_c; /* column permutation vector */
//	//    int            info;
//	//    superlu_options_t options;
//	//    SuperLUStat_t stat;	
//	
//	FLT newton_norm,tol=1.0e-14;
//	int max_newton_its = 50;
//	
//    /* Set the default input options:
//	 options.Fact = DOFACT;
//	 options.Equil = YES;
//	 options.ColPerm = COLAMD;
//	 options.DiagPivotThresh = 0.1; //different from complete LU
//	 options.Trans = NOTRANS;
//	 options.IterRefine = NOREFINE;
//	 options.SymmetricMode = NO;
//	 options.PivotGrowth = NO;
//	 options.ConditionNumber = NO;
//	 options.PrintStat = YES;
//	 options.RowPerm = LargeDiag;
//	 options.ILU_DropTol = 1e-4;
//	 options.ILU_FillTol = 1e-2;
//	 options.ILU_FillFactor = 10.0;
//	 options.ILU_DropRule = DROP_BASIC | DROP_AREA;
//	 options.ILU_Norm = INF_NORM;
//	 options.ILU_MILU = SMILU_2;
//     */
//    
//	ilu_set_default_options(&options);	
//	//options.PrintStat = NO;
////    options.PivotGrowth = YES;	  /* Compute reciprocal pivot growth */
////    options.ConditionNumber = YES;/* Compute reciprocal condition number */
//
//	/* send global solution to ug_vec */
//	ug_to_vec();
//	
//	bool compressed_column = true;
//	if ( !(perm_c = intMalloc(size_sparse_matrix)) ) ABORT("Malloc fails for perm_c[].");
//	
//	
//	for(int i = 0; i < max_newton_its; ++i) {
//		
//		if ( !(perm_r = intMalloc(size_sparse_matrix)) ) ABORT("Malloc fails for perm_r[].");
//		
//		
//		/* zero out sparse and residual but keep sparsity pattern */
//		zero_sparse();
//		
//		dCreate_Dense_Matrix(&X, size_sparse_matrix, 1, res_vec.data(), size_sparse_matrix, SLU_DN, SLU_D, SLU_GE);
//		
//		/* create jacobian */
//		create_jacobian(compressed_column);
//		
//		/* create residual */
//		create_rsdl();	
//		
//		/* apply dirichlet boundary conditions to sparse matrix and vector */
//		for(int j = 0; j < nfbd; ++j)
//			hp_fbdry(j)->apply_sparse_dirichlet(compressed_column);
//		
//		FLT max_resid = 0.0;
//		for (int j = 0; j < size_sparse_matrix; ++j){
//			if(fabs(res_vec(j)) > max_resid)
//				max_resid = fabs(res_vec(j));
//		}
//		cout << "L infinity norm of residual " << max_resid << endl;
//		
//		if (max_resid < tol) break;		
//		
//		/* create super matrix A using compressed column storage */
//		dCreate_CompCol_Matrix(&A, size_sparse_matrix, size_sparse_matrix, number_sparse_elements, sparse_val.data(), sparse_ind.data(), sparse_ptr.data(), SLU_NC, SLU_D, SLU_GE);
//		
//		if ( !(etree = intMalloc(size_sparse_matrix)) ) ABORT("Malloc fails for etree[].");
//		if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) ) 
//			ABORT("SUPERLU_MALLOC fails for R[].");
//		if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
//			ABORT("SUPERLU_MALLOC fails for C[].");
//		if ( !(ferr = (double *) SUPERLU_MALLOC(1 * sizeof(double))) )
//			ABORT("SUPERLU_MALLOC fails for ferr[].");
//		if ( !(berr = (double *) SUPERLU_MALLOC(1 * sizeof(double))) ) 
//			ABORT("SUPERLU_MALLOC fails for berr[].");
//		
//		
//		dCreate_Dense_Matrix(&B, size_sparse_matrix, 1, res_vec.data(), size_sparse_matrix, SLU_DN, SLU_D, SLU_GE);
//		
//		/* Initialize the statistics variables. */
//		StatInit(&stat);
//		
//		/* Compute the incomplete factorization and compute the condition number and pivot growth using dgsisx. */
//		dgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work,lwork, &B, &X, &rpg, &rcond, &mem_usage, &stat, &info);
//		
//		/* get access to solution and store in du */
//		//double *du = (double*) ((DNformat*) X.Store)->nzval;
//			
//		if (info > 0) {
//			cout << "something went wrong with superlu \n";
//			exit(1);
//		}
//		
//		if (info > 0 || rcond < 1e-8 || rpg > 1e8)
//			printf("WARNING: This preconditioner might be unstable.\n");
//		
//		restrt = SUPERLU_MIN(n / 3 + 1, 50);
//		maxit = 1000;
//		iter = maxit;
//		resid = 1e-8;
//		if (!(b = doubleMalloc(m))) ABORT("Malloc fails for b[].");
//		if (!(x = doubleMalloc(n))) ABORT("Malloc fails for x[].");
//		
//		for (int i = 0; i < size_sparse_matrix; i++)
//			x[i] = 0.0;
//		
//		/* call GMRES */
//		dfgmr(size_sparse_matrix, dmatvec_mult, dpsolve, b, x, resid, restrt, &iter, stdout);
//
//			
//		for(int j = 0; j < size_sparse_matrix; ++j)
//			ug_vec(j)-=x[j];
//		
//		/* send ug_vec to global solution */
//		vec_to_ug();
//		
//		if ( options.PrintStat ) StatPrint(&stat);
//		StatFree(&stat);
//		
//		/* once factorized can keep same pattern for next iteration */
//		options.Fact = SamePattern;
//		SUPERLU_FREE (perm_r);
//		SUPERLU_FREE (etree);
//		SUPERLU_FREE (R);
//		SUPERLU_FREE (C);
//		SUPERLU_FREE (ferr);
//		SUPERLU_FREE (berr);
//		
//		
//	}
//	
//	SUPERLU_FREE (perm_c);
//	
//	return;
//}
//
//
//int *GLOBAL_PERM_C, *GLOBAL_PERM_R;
//SuperMatrix *GLOBAL_A, *GLOBAL_L, *GLOBAL_U;
//SuperLUStat_t *GLOBAL_STAT;
//
//void dmatvec_mult(double alpha, double x[], double beta, double y[])
//{
//    SuperMatrix *A = GLOBAL_A;
//	
//    sp_dgemv("N", alpha, A, x, 1, beta, y, 1);
//}
//
//void dpsolve(int n, double x[], double y[])
//{
//    extern void dcopy_(int *, double [], int *, double [], int *);
//	
//    int i_1 = 1;
//    SuperMatrix *L = GLOBAL_L, *U = GLOBAL_U;
//    SuperLUStat_t *stat = GLOBAL_STAT;
//    int *perm_c = GLOBAL_PERM_C, *perm_r = GLOBAL_PERM_R;
//    int info;
//    static DNformat X;
//    static SuperMatrix XX = {SLU_DN, SLU_D, SLU_GE, 1, 1, &X};
//	
//    dcopy_(&n, y, &i_1, x, &i_1);
//    XX.nrow = n;
//    X.lda = n;
//    X.nzval = x;
//    dgstrs(NOTRANS, L, U, perm_c, perm_r, &XX, stat, &info);
//}
//
//
//
//
///*! @file dfgmr.c
// * \brief flexible GMRES written by Yousef Saad.
// */
//
//#define  epsmac  1.0e-16
//
//extern double ddot_(int *, double [], int *, double [], int *);
//extern double dnrm2_(int *, double [], int *);
//
//
//int dfgmr(int n,
//		  void (*dmatvec) (double, double[], double, double[]),
//		  void (*dpsolve) (int, double[], double[]),
//		  double *rhs, double *sol, double tol, int im, int *itmax, FILE * fits)
//{
//	/*----------------------------------------------------------------------
//	 |                 *** Preconditioned FGMRES ***
//	 +-----------------------------------------------------------------------
//	 | This is a simple version of the ARMS preconditioned FGMRES algorithm.
//	 +-----------------------------------------------------------------------
//	 | Y. S. Dec. 2000. -- Apr. 2008
//	 +-----------------------------------------------------------------------
//	 | on entry:
//	 |----------
//	 |
//	 | rhs     = real vector of length n containing the right hand side.
//	 | sol     = real vector of length n containing an initial guess to the
//	 |           solution on input.
//	 | tol     = tolerance for stopping iteration
//	 | im      = Krylov subspace dimension
//	 | (itmax) = max number of iterations allowed.
//	 | fits    = NULL: no output
//	 |        != NULL: file handle to output " resid vs time and its"
//	 |
//	 | on return:
//	 |----------
//	 | fgmr      int =  0 --> successful return.
//	 |           int =  1 --> convergence not achieved in itmax iterations.
//	 | sol     = contains an approximate solution (upon successful return).
//	 | itmax   = has changed. It now contains the number of steps required
//	 |           to converge --
//	 +-----------------------------------------------------------------------
//	 | internal work arrays:
//	 |----------
//	 | vv      = work array of length [im+1][n] (used to store the Arnoldi
//	 |           basis)
//	 | hh      = work array of length [im][im+1] (Householder matrix)
//	 | z       = work array of length [im][n] to store preconditioned vectors
//	 +-----------------------------------------------------------------------
//	 | subroutines called :
//	 | matvec - matrix-vector multiplication operation
//	 | psolve - (right) preconditionning operation
//	 |	   psolve can be a NULL pointer (GMRES without preconditioner)
//	 +---------------------------------------------------------------------*/
//	
//    int maxits = *itmax;
//    int i, i1, ii, j, k, k1, its, retval, i_1 = 1, i_2 = 2;
//    double beta, eps1 = 0.0, t, t0, gam;
//    double **hh, *c, *s, *rs;
//    double **vv, **z, tt;
//    double zero = 0.0;
//    double one = 1.0;
//	
//    its = 0;
//    vv = (double **)SUPERLU_MALLOC((im + 1) * sizeof(double *));
//    for (i = 0; i <= im; i++) vv[i] = doubleMalloc(n);
//    z = (double **)SUPERLU_MALLOC(im * sizeof(double *));
//    hh = (double **)SUPERLU_MALLOC(im * sizeof(double *));
//    for (i = 0; i < im; i++)
//    {
//		hh[i] = doubleMalloc(i + 2);
//		z[i] = doubleMalloc(n);
//    }
//    c = doubleMalloc(im);
//    s = doubleMalloc(im);
//    rs = doubleMalloc(im + 1);
//	
//    /*---- outer loop starts here ----*/
//    do
//    {
//		/*---- compute initial residual vector ----*/
//		dmatvec(one, sol, zero, vv[0]);
//		for (j = 0; j < n; j++)
//			vv[0][j] = rhs[j] - vv[0][j];	/* vv[0]= initial residual */
//		beta = dnrm2_(&n, vv[0], &i_1);
//		
//		/*---- print info if fits != null ----*/
//		if (fits != NULL && its == 0)
//			fprintf(fits, "%8d   %10.2e\n", its, beta);
//		/*if ( beta <= tol * dnrm2_(&n, rhs, &i_1) )*/
//		if ( !(beta > tol * dnrm2_(&n, rhs, &i_1)) )
//			break;
//		t = 1.0 / beta;
//		
//		/*---- normalize: vv[0] = vv[0] / beta ----*/
//		for (j = 0; j < n; j++)
//			vv[0][j] = vv[0][j] * t;
//		if (its == 0)
//			eps1 = tol * beta;
//		
//		/*---- initialize 1-st term of rhs of hessenberg system ----*/
//		rs[0] = beta;
//		for (i = 0; i < im; i++)
//		{
//			its++;
//			i1 = i + 1;
//			
//			/*------------------------------------------------------------
//			 |  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
//			 +-----------------------------------------------------------*/
//			if (dpsolve)
//				dpsolve(n, z[i], vv[i]);
//			else
//				dcopy_(&n, vv[i], &i_1, z[i], &i_1);
//			
//			/*---- matvec operation w = A z_{j} = A M^{-1} v_{j} ----*/
//			dmatvec(one, z[i], zero, vv[i1]);
//			
//			/*------------------------------------------------------------
//			 |     modified gram - schmidt...
//			 |     h_{i,j} = (w,v_{i})
//			 |     w  = w - h_{i,j} v_{i}
//			 +------------------------------------------------------------*/
//			t0 = dnrm2_(&n, vv[i1], &i_1);
//			for (j = 0; j <= i; j++)
//			{
//				double negt;
//				tt = ddot_(&n, vv[j], &i_1, vv[i1], &i_1);
//				hh[i][j] = tt;
//				negt = -tt;
//				daxpy_(&n, &negt, vv[j], &i_1, vv[i1], &i_1);
//			}
//			
//			/*---- h_{j+1,j} = ||w||_{2} ----*/
//			t = dnrm2_(&n, vv[i1], &i_1);
//			while (t < 0.5 * t0)
//			{
//				t0 = t;
//				for (j = 0; j <= i; j++)
//				{
//					double negt;
//					tt = ddot_(&n, vv[j], &i_1, vv[i1], &i_1);
//					hh[i][j] += tt;
//					negt = -tt;
//					daxpy_(&n, &negt, vv[j], &i_1, vv[i1], &i_1);
//				}
//				t = dnrm2_(&n, vv[i1], &i_1);
//			}
//			
//			hh[i][i1] = t;
//			
//			if (t != 0.0)
//			{
//				/*---- v_{j+1} = w / h_{j+1,j} ----*/
//				t = 1.0 / t;
//				for (k = 0; k < n; k++)
//					vv[i1][k] = vv[i1][k] * t;
//			}
//			/*---------------------------------------------------
//			 |     done with modified gram schimdt and arnoldi step
//			 |     now  update factorization of hh
//			 +--------------------------------------------------*/
//			
//			/*--------------------------------------------------------
//			 |   perform previous transformations  on i-th column of h
//			 +-------------------------------------------------------*/
//			for (k = 1; k <= i; k++)
//			{
//				k1 = k - 1;
//				tt = hh[i][k1];
//				hh[i][k1] = c[k1] * tt + s[k1] * hh[i][k];
//				hh[i][k] = -s[k1] * tt + c[k1] * hh[i][k];
//			}
//			
//			gam = sqrt(pow(hh[i][i], 2) + pow(hh[i][i1], 2));
//			
//			/*---------------------------------------------------
//			 |     if gamma is zero then any small value will do
//			 |     affect only residual estimate
//			 +--------------------------------------------------*/
//			/* if (gam == 0.0) gam = epsmac; */
//			
//			/*---- get next plane rotation ---*/
//			if (gam == 0.0)
//			{
//				c[i] = one;
//				s[i] = zero;
//			}
//            else
//			{
//				c[i] = hh[i][i] / gam;
//				s[i] = hh[i][i1] / gam;
//			}
//			
//			rs[i1] = -s[i] * rs[i];
//			rs[i] = c[i] * rs[i];
//			
//			/*----------------------------------------------------
//			 |   determine residual norm and test for convergence
//			 +---------------------------------------------------*/
//			hh[i][i] = c[i] * hh[i][i] + s[i] * hh[i][i1];
//			beta = fabs(rs[i1]);
//			if (fits != NULL)
//				fprintf(fits, "%8d   %10.2e\n", its, beta);
//			if (beta <= eps1 || its >= maxits)
//				break;
//		}
//		
//		if (i == im) i--;
//		
//		/*---- now compute solution. 1st, solve upper triangular system ----*/
//		rs[i] = rs[i] / hh[i][i];
//		
//		for (ii = 1; ii <= i; ii++)
//		{
//			k = i - ii;
//			k1 = k + 1;
//			tt = rs[k];
//			for (j = k1; j <= i; j++)
//				tt = tt - hh[j][k] * rs[j];
//			rs[k] = tt / hh[k][k];
//		}
//		
//		/*---- linear combination of v[i]'s to get sol. ----*/
//		for (j = 0; j <= i; j++)
//		{
//			tt = rs[j];
//			for (k = 0; k < n; k++)
//				sol[k] += tt * z[j][k];
//		}
//		
//		/* calculate the residual and output */
//		dmatvec(one, sol, zero, vv[0]);
//		for (j = 0; j < n; j++)
//			vv[0][j] = rhs[j] - vv[0][j];	/* vv[0]= initial residual */
//		
//		/*---- print info if fits != null ----*/
//		beta = dnrm2_(&n, vv[0], &i_1);
//		
//		/*---- restart outer loop if needed ----*/
//		/*if (beta >= eps1 / tol)*/
//		if ( !(beta < eps1 / tol) )
//		{
//			its = maxits + 10;
//			break;
//		}
//		if (beta <= eps1)
//			break;
//    } while(its < maxits);
//	
//    retval = (its >= maxits);
//    for (i = 0; i <= im; i++)
//		SUPERLU_FREE(vv[i]);
//    SUPERLU_FREE(vv);
//    for (i = 0; i < im; i++)
//    {
//		SUPERLU_FREE(hh[i]);
//		SUPERLU_FREE(z[i]);
//    }
//    SUPERLU_FREE(hh);
//    SUPERLU_FREE(z);
//    SUPERLU_FREE(c);
//    SUPERLU_FREE(s);
//    SUPERLU_FREE(rs);
//	
//    *itmax = its;
//	
//    return retval;
//} /*----end of fgmr ----*/
