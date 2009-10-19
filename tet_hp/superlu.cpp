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
//	NCformat	*Astore;	
//	double *xact;
//	lwork = 0;// superlu determines memory allocation
	
	
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

	for(int i = 0; i < max_newton_its; ++i) {

		if ( !(perm_c = intMalloc(size_sparse_matrix)) ) ABORT("Malloc fails for perm_c[].");
		if ( !(perm_r = intMalloc(size_sparse_matrix)) ) ABORT("Malloc fails for perm_r[].");
		

		/* zero out sparse and residual but keep sparsity pattern */
		zero_sparse();
		
		/* create jacobian */
		create_jacobian(compressed_column);

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

		//cout << "residual "<<res_vec << endl;
		
		dCreate_Dense_Matrix(&B, size_sparse_matrix, 1, res_vec.data(), size_sparse_matrix, SLU_DN, SLU_D, SLU_GE);
		//dPrint_Dense_Matrix("B", &B);

		/* Initialize the statistics variables. */
		StatInit(&stat);

//		/* advanced superlu driver */
//		dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,&mem_usage, &stat, &info);		
//		/* get access to solution and store in du */
//		double *du = (double*) ((DNformat*) X.Store)->nzval;

			
//		/* not sure what this crap does */
//		xact = doubleMalloc(size_sparse_matrix);
//		dGenXtrue(size_sparse_matrix, 1, xact, size_sparse_matrix);
//		dFillRHS(options.Trans, 1, xact, size_sparse_matrix, &A, &B);

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
		//options.Fact = SamePattern;
		SUPERLU_FREE (perm_r);
		SUPERLU_FREE (perm_c);
//		SUPERLU_FREE (xact);
//		SUPERLU_FREE (etree);
//		SUPERLU_FREE (R);
//		SUPERLU_FREE (C);
//		SUPERLU_FREE (ferr);
//		SUPERLU_FREE (berr);

	}

	//cout << "solution "<<ug_vec << endl;
	
//	ofstream out;
//	out.open("superlu.dat");
//
//	out <<"ZONE F=FEPOINT, ET=TETRAHEDRON, N = "<<npnt <<", E = " << ntet << endl;
//	for(int i = 0; i < npnt; ++i)
//		out << pnts(i)(0) << ' ' << pnts(i)(1) << ' ' << pnts(i)(2) <<' ' << ug_vec(i) << endl;
//	for(int i = 0; i < ntet; ++i)
//		out << tet(i).pnt(0)+1 << ' ' << tet(i).pnt(1)+1 << ' ' <<tet(i).pnt(2)+1 << ' '<<tet(i).pnt(3)+1 << endl;
//	out.close();
	

//	
//	SUPERLU_FREE (xact);
//	SUPERLU_FREE (etree);
//	SUPERLU_FREE (R);
//	SUPERLU_FREE (C);
//	SUPERLU_FREE (ferr);
//	SUPERLU_FREE (berr);

	return;
}



//void superilu_gmres(){
//	
//	
////	/* make all these global */
////	GLOBAL_A = &A;
////    GLOBAL_L = &L;
////    GLOBAL_U = &U;
////    GLOBAL_STAT = &stat;
////    GLOBAL_PERM_C = perm_c;
////    GLOBAL_PERM_R = perm_r;
////	
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
//	
//	lwork = 0;// superlu determines memory allocation
//	
//	
//	/* simple driver dgssv allocation */
//	//    SuperMatrix    A, L, U;
//	//    SuperMatrix    B;
//	//	int            *perm_r; /* row permutations from partial pivoting */
//	//    int            *perm_c; /* column permutation vector */
//	//    int            info;
//	//    superlu_options_t options;
//	//    SuperLUStat_t stat;	
//	
//	FLT newton_norm,tol=1.0e-12;
//	int max_newton_its = 100;
//	
//    ilu_set_default_options(&options);
//	
//    /* Modify the defaults. */
//    options.PivotGrowth = YES;	  /* Compute reciprocal pivot growth */
//    options.ConditionNumber = YES;/* Compute reciprocal condition number */
//	
//	//options.PrintStat = NO;
//	
//	/* send global solution to ug_vec */
//	ug_to_vec();
//	
//	for(int i = 0; i < max_newton_its; ++i) {
//		
//		zero_sparse();
//		create_jacobian(1);
//		
//		/* create super matrix A using compressed column storage */
//		dCreate_CompCol_Matrix(&A, size_sparse_matrix, size_sparse_matrix, number_sparse_elements, sparse_val.data(), sparse_ind.data(), sparse_ptr.data(), SLU_NC, SLU_D, SLU_GE);
//		
//		/* create residual and store in B */
//		create_rsdl();	
//		dCreate_Dense_Matrix(&B, size_sparse_matrix, 1, res_vec.data(), size_sparse_matrix, SLU_DN, SLU_D, SLU_GE);
//		
//		/* Initialize the statistics variables. */
//		StatInit(&stat);
//		
//		/* incomplete lu factorization advanced superlu driver */
//		dgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, &mem_usage, &stat, &info);
//		
//		if (info > 0 || rcond < 1e-8 || rpg > 1e8)
//			printf("WARNING: This preconditioner might be unstable.\n");
//
//		restrt = SUPERLU_MIN(n / 3 + 1, 50);
//		maxit = 1000;
//		iter = maxit;
//		resid = 1e-8;
//
//		sp_dgemv("N", 1.0, &A, xact, 1, 0.0, b, 1);
//		
//		/* call gmres */
//		
//		
//		/* get access to solution and store in du */
//		double *du = (double*) ((DNformat*) X.Store)->nzval;
//		
//		for(int i = 0; i < size_sparse_matrix; ++i)
//			ug_vec(i)-=du[i];
//		
//		/* send ug_vec to global solution */
//		vec_to_ug();
//		
//		if ( options.PrintStat ) StatPrint(&stat);
//		StatFree(&stat);
//		Destroy_CompCol_Matrix(&A);
//		Destroy_Dense_Matrix(&B);
//		Destroy_SuperNode_Matrix(&L);
//		Destroy_CompCol_Matrix(&U);
//		
//		newton_norm = 0.0;			
//		for (int j = 0; j < size_sparse_matrix; ++j)
//			newton_norm += du[j]*du[j];
//		
//		if (newton_norm < tol*tol) break;
//		
//		/* once factorized can keep same pattern for next iteration */
//		options.Fact = SamePattern;
//	}
//	
//    SUPERLU_FREE (xact);
//    SUPERLU_FREE (etree);
//    SUPERLU_FREE (perm_r);
//    SUPERLU_FREE (perm_c);
//    SUPERLU_FREE (R);
//    SUPERLU_FREE (C);
//    SUPERLU_FREE (ferr);
//    SUPERLU_FREE (berr);
//	
//	return;
//}
//
///*! @file dfgmr.c
// * \brief flexible GMRES written by Yousef Saad.
// */
//
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
//
//
//
///*
// * -- SuperLU routine (version 4.0) --
// * Lawrence Berkeley National Laboratory
// * June 30, 2009
// */
//
//int *GLOBAL_PERM_C, *GLOBAL_PERM_R;
//SuperMatrix *GLOBAL_A, *GLOBAL_L, *GLOBAL_U;
//SuperLUStat_t *GLOBAL_STAT;
//
//void dmatvec_mult(double alpha, double x[], double beta, double y[])
//{
//    SuperMatrix *A = GLOBAL_A;
////	y := alpha*A*x + beta*y
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
//int main(int argc, char *argv[])
//{
//    void dmatvec_mult(double alpha, double x[], double beta, double y[]);
//    void dpsolve(int n, double x[], double y[]);
//    extern int dfgmr( int n,
//					 void (*matvec_mult)(double, double [], double, double []),
//					 void (*psolve)(int n, double [], double[]),
//					 double *rhs, double *sol, double tol, int restrt, int *itmax,
//					 FILE *fits);
//    extern int dfill_diag(int n, NCformat *Astore);
//	
//    char     equed[1] = {'B'};
//    yes_no_t equil;
//    trans_t  trans;
//    SuperMatrix A, L, U;
//    SuperMatrix B, X;
//    NCformat *Astore;
//    NCformat *Ustore;
//    SCformat *Lstore;
//    double   *a;
//    int      *asub, *xa;
//    int      *etree;
//    int      *perm_c; /* column permutation vector */
//    int      *perm_r; /* row permutations from partial pivoting */
//    int      nrhs, ldx, lwork, info, m, n, nnz;
//    double   *rhsb, *rhsx, *xact;
//    double   *work = NULL;
//    double   *R, *C;
//    double   u, rpg, rcond;
//    double zero = 0.0;
//    double one = 1.0;
//    mem_usage_t   mem_usage;
//    superlu_options_t options;
//    SuperLUStat_t stat;
//	
//    int restrt, iter, maxit, i;
//    double resid;
//    double *x, *b;
//	
//#ifdef DEBUG
//    extern int num_drop_L, num_drop_U;
//#endif
//	
//#if ( DEBUGlevel>=1 )
//    CHECK_MALLOC("Enter main()");
//#endif
//	
//    /* Defaults */
//    lwork = 0;
//    nrhs  = 1;
//    equil = YES;
//    u	  = 0.1; /* u=1.0 for complete factorization */
//    trans = NOTRANS;
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
//    ilu_set_default_options(&options);
//	
//    /* Modify the defaults. */
//    options.PivotGrowth = YES;	  /* Compute reciprocal pivot growth */
//    options.ConditionNumber = YES;/* Compute reciprocal condition number */
//	
//    if ( lwork > 0 ) {
//		work = SUPERLU_MALLOC(lwork);
//		if ( !work ) ABORT("Malloc fails for work[].");
//    }
//	
//    /* Read matrix A from a file in Harwell-Boeing format.*/
//    if (argc < 2)
//    {
//		printf("Usage:\n%s [OPTION] < [INPUT] > [OUTPUT]\nOPTION:\n"
//			   "-h -hb:\n\t[INPUT] is a Harwell-Boeing format matrix.\n"
//			   "-r -rb:\n\t[INPUT] is a Rutherford-Boeing format matrix.\n"
//			   "-t -triplet:\n\t[INPUT] is a triplet format matrix.\n",
//			   argv[0]);
//		return 0;
//    }
//    else
//    {
//		switch (argv[1][1])
//		{
//			case 'H':
//			case 'h':
//				printf("Input a Harwell-Boeing format matrix:\n");
//				dreadhb(&m, &n, &nnz, &a, &asub, &xa);
//				break;
//			case 'R':
//			case 'r':
//				printf("Input a Rutherford-Boeing format matrix:\n");
//				dreadrb(&m, &n, &nnz, &a, &asub, &xa);
//				break;
//			case 'T':
//			case 't':
//				printf("Input a triplet format matrix:\n");
//				dreadtriple(&m, &n, &nnz, &a, &asub, &xa);
//				break;
//			default:
//				printf("Unrecognized format.\n");
//				return 0;
//		}
//    }
//	
//    dCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
//    Astore = A.Store;
//    dfill_diag(n, Astore);
//    printf("Dimension %dx%d; # nonzeros %d\n", A.nrow, A.ncol, Astore->nnz);
//    fflush(stdout);
//	
//    if ( !(rhsb = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsb[].");
//    if ( !(rhsx = doubleMalloc(m * nrhs)) ) ABORT("Malloc fails for rhsx[].");
//    dCreate_Dense_Matrix(&B, m, nrhs, rhsb, m, SLU_DN, SLU_D, SLU_GE);
//    dCreate_Dense_Matrix(&X, m, nrhs, rhsx, m, SLU_DN, SLU_D, SLU_GE);
//    xact = doubleMalloc(n * nrhs);
//    ldx = n;
//    dGenXtrue(n, nrhs, xact, ldx);
//    dFillRHS(trans, nrhs, xact, ldx, &A, &B);
//	
//    if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
//    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
//    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
//    if ( !(R = (double *) SUPERLU_MALLOC(A.nrow * sizeof(double))) )
//		ABORT("SUPERLU_MALLOC fails for R[].");
//    if ( !(C = (double *) SUPERLU_MALLOC(A.ncol * sizeof(double))) )
//		ABORT("SUPERLU_MALLOC fails for C[].");
//	
//    info = 0;
//#ifdef DEBUG
//    num_drop_L = 0;
//    num_drop_U = 0;
//#endif
//	
//    /* Initialize the statistics variables. */
//    StatInit(&stat);
//	
//    /* Compute the incomplete factorization and compute the condition number
//	 and pivot growth using dgsisx. */
//    dgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work,
//		   lwork, &B, &X, &rpg, &rcond, &mem_usage, &stat, &info);
//	
//    Lstore = (SCformat *) L.Store;
//    Ustore = (NCformat *) U.Store;
//    printf("dgsisx(): info %d\n", info);
//    if (info > 0 || rcond < 1e-8 || rpg > 1e8)
//		printf("WARNING: This preconditioner might be unstable.\n");
//	
//    if ( info == 0 || info == n+1 ) {
//		
//		if ( options.PivotGrowth == YES )
//			printf("Recip. pivot growth = %e\n", rpg);
//		if ( options.ConditionNumber == YES )
//			printf("Recip. condition number = %e\n", rcond);
//		
//    } else if ( info > 0 && lwork == -1 ) {
//		printf("** Estimated memory: %d bytes\n", info - n);
//    }
//    printf("n(A) = %d, nnz(A) = %d\n", n, Astore->nnz);
//    printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
//    printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
//    printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
//    printf("Fill ratio: nnz(F)/nnz(A) = %.3f\n",
//		   ((double)(Lstore->nnz) + (double)(Ustore->nnz) - (double)n)
//		   / (double)Astore->nnz);
//    printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
//		   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
//    fflush(stdout);
//	
//    /* Set the global variables. */
//    GLOBAL_A = &A;
//    GLOBAL_L = &L;
//    GLOBAL_U = &U;
//    GLOBAL_STAT = &stat;
//    GLOBAL_PERM_C = perm_c;
//    GLOBAL_PERM_R = perm_r;
//	
//    /* Set the variables used by GMRES. */
//    restrt = SUPERLU_MIN(n / 3 + 1, 50);
//    maxit = 1000;
//    iter = maxit;
//    resid = 1e-8;
//    if (!(b = doubleMalloc(m))) ABORT("Malloc fails for b[].");
//    if (!(x = doubleMalloc(n))) ABORT("Malloc fails for x[].");
//    sp_dgemv("N", one, &A, xact, 1, zero, b, 1);
//	// b=A*xact
//	
//    if (info <= n + 1)
//    {
//		int i_1 = 1;
//		double maxferr = 0.0, nrmA, nrmB, res, t;
//        double temp;
//		extern double dnrm2_(int *, double [], int *);
//		extern void daxpy_(int *, double *, double [], int *, double [], int *);
//		
//		/* Call GMRES. */
//		/*double *sol = (double*) ((DNformat*) X.Store)->nzval;
//		 for (i = 0; i < n; i++) x[i] = sol[i];*/
//		for (i = 0; i < n; i++) x[i] = zero;
//		
//		t = SuperLU_timer_();
//		
//		dfgmr(n, dmatvec_mult, dpsolve, b, x, resid, restrt, &iter, stdout);
//		
//		t = SuperLU_timer_() - t;
//		
//		/* Output the result. */
//		nrmA = dnrm2_(&(Astore->nnz), (double *)((DNformat *)A.Store)->nzval,
//					  &i_1);
//		nrmB = dnrm2_(&m, b, &i_1);
////		y := alpha*A*x + beta*y
//		sp_dgemv("N", -1.0, &A, x, 1, 1.0, b, 1);
//		//b=-Ax+b;residual
//		res = dnrm2_(&m, b, &i_1);
//		resid = res / nrmB;
//		printf("||A||_F = %.1e, ||B||_2 = %.1e, ||B-A*X||_2 = %.1e, "
//			   "relres = %.1e\n", nrmA, nrmB, res, resid);
//		
//		if (iter >= maxit)
//		{
//			if (resid >= 1.0) iter = -180;
//			else if (resid > 1e-8) iter = -111;
//		}
//		printf("iteration: %d\nresidual: %.1e\nGMRES time: %.2f seconds.\n",
//			   iter, resid, t);
//		
//		for (i = 0; i < m; i++)
//			maxferr = SUPERLU_MAX(maxferr, fabs(x[i] - xact[i]));
//		printf("||X-X_true||_oo = %.1e\n", maxferr);
//    }
//#ifdef DEBUG
//    printf("%d entries in L and %d entries in U dropped.\n",
//		   num_drop_L, num_drop_U);
//#endif
//    fflush(stdout);
//	
//    if ( options.PrintStat ) StatPrint(&stat);
//    StatFree(&stat);
//	
//    SUPERLU_FREE (rhsb);
//    SUPERLU_FREE (rhsx);
//    SUPERLU_FREE (xact);
//    SUPERLU_FREE (etree);
//    SUPERLU_FREE (perm_r);
//    SUPERLU_FREE (perm_c);
//    SUPERLU_FREE (R);
//    SUPERLU_FREE (C);
//    Destroy_CompCol_Matrix(&A);
//    Destroy_SuperMatrix_Store(&B);
//    Destroy_SuperMatrix_Store(&X);
//    if ( lwork >= 0 ) {
//		Destroy_SuperNode_Matrix(&L);
//		Destroy_CompCol_Matrix(&U);
//    }
//    SUPERLU_FREE(b);
//    SUPERLU_FREE(x);
//	
//#if ( DEBUGlevel>=1 )
//    CHECK_MALLOC("Exit main()");
//#endif
//	
//    return 0;
//}
//
//
//
