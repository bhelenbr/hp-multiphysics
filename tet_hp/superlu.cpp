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
	
    /* advanced driver dgssvx allocation */ 
    char           equed[1];
    SuperMatrix    A, L, U;
    SuperMatrix    B, X;
	int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work;
    int            info, lwork;
    double         *R, *C;
    double         *ferr, *berr;
    double         rpg, rcond;
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t  stat;
	lwork = 0; // superlu determines memory allocation
	//double *du;
	//du = doubleMalloc(size_sparse_matrix);
	Array<FLT,1> du(size_sparse_matrix);
	
	du = 0.0;
	
	FLT tol=1.0e-14;
	int max_newton_its = 50;
	
	bool compressed_column = true;
	
	perm_c = intMalloc(size_sparse_matrix);
	perm_r = intMalloc(size_sparse_matrix);
	etree = intMalloc(size_sparse_matrix);

	ferr = (double *) SUPERLU_MALLOC(1 * sizeof(double));
	berr = (double *) SUPERLU_MALLOC(1 * sizeof(double));

	R = (double *) SUPERLU_MALLOC(size_sparse_matrix * sizeof(double));	
	C = (double *) SUPERLU_MALLOC(size_sparse_matrix * sizeof(double));
	
	/* create super matrix A using compressed column storage */
	dCreate_CompCol_Matrix(&A, size_sparse_matrix, size_sparse_matrix, number_sparse_elements, sparse_val.data(), sparse_ind.data(), sparse_ptr.data(), SLU_NC, SLU_D, SLU_GE);
	
	dCreate_Dense_Matrix(&X, size_sparse_matrix, 1, du.data(), size_sparse_matrix, SLU_DN, SLU_D, SLU_GE);

	dCreate_Dense_Matrix(&B, size_sparse_matrix, 1, res_vec.data(), size_sparse_matrix, SLU_DN, SLU_D, SLU_GE);
	   
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
    
	ilu_set_default_options(&options);	
	//set_default_options(&options);	
		
	/* send global solution to ug_vec */
	ug_to_vec();

	for(int i = 0; i < max_newton_its; ++i) {
		
		/* zero out sparse and residual but keep sparsity pattern */
		zero_sparse();		

		/* create jacobian */
		FLT t = SuperLU_timer_();
		create_jacobian(compressed_column);
		t = SuperLU_timer_()-t;
		cout << "time to make jacobian " << t << endl; 
		
		/* create residual */
		create_rsdl();
		
		/* apply neumman bc's */
		apply_neumman(compressed_column);
				
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

		/* Initialize the statistics variables. */
		StatInit(&stat);

		/* advanced superlu driver */
		//dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,&mem_usage, &stat, &info);				
			
		/* ILU driver */
		dgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, &mem_usage, &stat, &info);
		
		if (info > 0) {
			cout << "something went wrong with superlu \n";
			exit(1);
		}
		
		int im = 50;
		int itmax = 1000;
		
		/* Flexible GMRES input approximate solution (du) returns solution du and number of iterations */ 
		fgmres(size_sparse_matrix, A, L, U, perm_c, perm_r, res_vec, du, tol, im, itmax, stat);


		
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
			   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
		
		/* update solution */
		ug_vec -= du;
		
		/* send ug_vec to global solution */
		vec_to_ug();
		
		if ( options.PrintStat ) StatPrint(&stat);
		StatFree(&stat);
		
		/* once factorized can keep same pattern for next iteration */
		options.Fact = SamePattern;
		//options.ColPerm = MY_PERMC;
		
		output("superlu", tecplot);

		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
	}

	SUPERLU_FREE (etree);
	SUPERLU_FREE (R);
	SUPERLU_FREE (C);
	SUPERLU_FREE (ferr);
	SUPERLU_FREE (berr);
	
	SUPERLU_FREE (perm_r);
	SUPERLU_FREE (perm_c);

//	Destroy_CompCol_Matrix(&A);
//	Destroy_Dense_Matrix(&B);
//	Destroy_Dense_Matrix(&X);

	
	return;
}

//void tet_hp::fgmres(int n, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, 
//					double *sol, double *rhs,int *perm_c, int *perm_r,
//					FLT tol, int im, int *itmax){
//	
//	//sol = (double*) ((DNformat*) X.Store)->nzval;
//	//double *rhs = (double*) ((DNformat*) B.Store)->nzval;
//
//	int maxits = *itmax;
//    int i, i1, ii, j, k, k1, its, retval, i_1 = 1, i_2 = 2;
//    double beta, eps1 = 0.0, t, t0, gam;
//    double **hh, *c, *s, *rs;
//    double **vv, **z, tt;
//    double zero = 0.0;
//    double one = 1.0;
//	
//	static DNformat Xdata;
//    static SuperMatrix XX = {SLU_DN, SLU_D, SLU_GE, 1, 1, &Xdata};
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
//		sp_dgemv("N", one, A, sol, 1, zero, vv[0], 1);
//		
//		for (j = 0; j < n; j++)
//			vv[0][j] = rhs[j] - vv[0][j];	/* vv[0]= initial residual */
//		beta = 0.0;
//		for(int ind=0;ind<n;++ind) beta+=vv[0][ind]*vv[0][ind];
//		beta = sqrt(beta);
//		
//		FLT fltwork = 0.0;
//		for(int ind=0;ind<n;++ind) fltwork+=rhs[ind]*rhs[ind];
//		fltwork=sqrt(fltwork);
//
//		if ( !(beta > tol * fltwork) )
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
//			
//			XX.nrow = n;
//			Xdata.lda = n;
//			Xdata.nzval = vv[i];
//			int info;
//			SuperLUStat_t stat;
//			StatInit(&stat);
//			dgstrs(NOTRANS, L, U, perm_c, perm_r, &XX, &stat, &info);
//			StatFree(&stat);
//			z[i] = (double*) ((DNformat*) XX.Store)->nzval;
//
//			/*---- matvec operation w = A z_{j} = A M^{-1} v_{j} ----*/
//			sp_dgemv("N", one, A, z[i], 1, zero, vv[i1], 1);
//
//			/*------------------------------------------------------------
//			 |     modified gram - schmidt...
//			 |     h_{i,j} = (w,v_{i})
//			 |     w  = w - h_{i,j} v_{i}
//			 +------------------------------------------------------------*/
//			t0=0.0;
//			for(int ind=0;ind<n;++ind) t0+=vv[i1][ind]*vv[i1][ind];
//			t0=sqrt(t0);
//			for (j = 0; j <= i; j++)
//			{
//				double negt;
//				tt = 0.0;
//				for(int ind=0;ind<n;++ind) t0+=vv[j][ind]*vv[i1][ind];
//				
//				hh[i][j] = tt;
//				negt = -tt;
//				for(int ind=0;ind<n;++ind) vv[i1][ind]+=vv[i1][ind]+negt*vv[j][ind];
//			}
//			
//			/*---- h_{j+1,j} = ||w||_{2} ----*/
//			t=0.0;
//			for(int ind=0;ind<n;++ind) t+=vv[i1][ind]*vv[i1][ind];
//			t=sqrt(t);
//			
//			while (t < 0.5 * t0)
//			{
//				t0 = t;
//				for (j = 0; j <= i; j++)
//				{
//					double negt;
//					tt = 0.0;
//					for(int ind=0;ind<n;++ind) tt+=vv[j][ind]*vv[i1][ind];
//
//					hh[i][j] += tt;
//					negt = -tt;
//					for(int ind=0;ind<n;++ind) vv[i1][ind]+=vv[i1][ind]+negt*vv[j][ind];
//				}
//				t=0.0;
//				for(int ind=0;ind<n;++ind) t+=vv[i1][ind]*vv[i1][ind];
//				t=sqrt(t);
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
////			if (fits != NULL)
////				fprintf(fits, "%8d   %10.2e\n", its, beta);
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
//		sp_dgemv("N",one,A,sol,1,zero,vv[0],1);
//		for (j = 0; j < n; j++)
//			vv[0][j] = rhs[j] - vv[0][j];	/* vv[0]= initial residual */
//		
//		/*---- print info if fits != null ----*/
//		beta = 0.0;
//		for(int ind=0;ind<n;++ind) beta+=vv[0][ind]*vv[0][ind];
//		beta = sqrt(beta);
//		
//		/*---- restart outer loop if needed ----*/
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
//	return;
//}