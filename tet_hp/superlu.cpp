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

#ifndef petsc

#include "tet_hp.h"
#include "slu_ddefs.h"
#include "hp_boundary.h"

void tet_hp::superilu(){
	cout << "superilu called" <<endl;

	
    /* advanced driver dgssvx allocation */ 
    char           equed[1];
    SuperMatrix    A, L, U;
	//NCformat       *Astore;
    SuperMatrix    B, X;
	trans_t  trans;
	int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work = NULL;
    int            info, lwork;
	int			   n=size_sparse_matrix;
    double         *R, *C;
    double         rpg, rcond;
	//double		   *xact;
	//double *x, *b;

	/* GMRES inputs */
	int im = 50; // krlov subspace
	int itmax = 5000; // total iterations
	FLT gtol = 1.0e-13; // tolerance
	
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t  stat;
	lwork = 0; // superlu determines memory allocation

	Array<double,1> du(n);
	Array<double,1> x(n),b(n);
	du = 0.0;
	trans = NOTRANS;
	FLT tol=1.0e-12;
	int max_newton_its = 50;
	bool compressed_column = true;
	
	perm_c = intMalloc(n);
	perm_r = intMalloc(n);
	etree = intMalloc(n);


	R = (double *) SUPERLU_MALLOC(n * sizeof(double));	
	C = (double *) SUPERLU_MALLOC(n * sizeof(double));
	
	//xact = doubleMalloc(n);
	
	//dGenXtrue(n, 1, xact,n);
    //dFillRHS(trans, 1, xact, n, &A, &B);
	
	
	/* create super matrix A using compressed column storage */
	dCreate_CompCol_Matrix(&A, n, n, number_sparse_elements, sparse_val.data(), sparse_ind.data(), sparse_ptr.data(), SLU_NC, SLU_D, SLU_GE);
	//Astore = (NCformat *) A.Store;
	
	dCreate_Dense_Matrix(&X, n, 1, du.data(), n, SLU_DN, SLU_D, SLU_GE);

	dCreate_Dense_Matrix(&B, n, 1, res_vec.data(), n, SLU_DN, SLU_D, SLU_GE);
	
	/* Set the default input options:
	 options.Fact = DOFACT;
	 options.Equil = YES;
	 options.ColPerm = COLAMD;
	 options.DiagPivotThresh = 0.1; //different from complete LU
	 options.Trans = NOTRANS;
	 options.IterRefine = NOREFINE;
	 options.SymmetricMode = NO;
	 options.PivotGrowth = NO;
	 options.ConditionNumber = NO;
	 options.PrintStat = YES;
	 options.RowPerm = LargeDiag;
	 options.ILU_DropTol = 1e-4;
	 options.ILU_FillTol = 1e-2;
	 options.ILU_FillFactor = 10.0;
	 options.ILU_DropRule = DROP_BASIC | DROP_AREA;
	 options.ILU_Norm = INF_NORM;
	 options.ILU_MILU = SMILU_2;
     */
    
	
//	typedef enum {NO, YES}                                          yes_no_t;
//	typedef enum {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED} fact_t;
//	typedef enum {NOROWPERM, LargeDiag, MY_PERMR}                   rowperm_t;
//	typedef enum {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD, MY_PERMC}colperm_t;
//	typedef enum {NOTRANS, TRANS, CONJ}                             trans_t;
//	typedef enum {NOEQUIL, ROW, COL, BOTH}                          DiagScale_t;
//	typedef enum {NOREFINE, SINGLE=1, DOUBLE, EXTRA}                IterRefine_t;
//	typedef enum {LUSUP, UCOL, LSUB, USUB}                          MemType;
//	typedef enum {HEAD, TAIL}                                       stack_end_t;
//	typedef enum {SYSTEM, USER}                                     LU_space_t;
//	typedef enum {ONE_NORM, TWO_NORM, INF_NORM}			norm_t;
//	typedef enum {SILU, SMILU_1, SMILU_2, SMILU_3}			milu_t;
	
	
	ilu_set_default_options(&options);	
	options.RowPerm = NOROWPERM;
	options.ColPerm = NATURAL;
	options.ColPerm = COLAMD;
	options.Equil = NO;
	
	options.PivotGrowth = YES;	  /* Compute reciprocal pivot growth */
	options.ConditionNumber = YES;/* Compute reciprocal condition number */


	options.DiagPivotThresh = 0.001;
	//options.ILU_DropRule = DROP_AREA;
	//options.ILU_DropTol = 1.0e-4;
	//options.ILU_FillTol = 1.0e-2;
	options.ILU_FillFactor = 3.0;
	options.ILU_MILU = SMILU_3;

	
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
		for (int j = 0; j < n; ++j){
			if(fabs(res_vec(j)) > max_resid)
				max_resid = fabs(res_vec(j));
		}
		cout << "L infinity norm of residual " << max_resid << endl;		

		if (max_resid < tol) break;		

		/* Initialize the statistics variables. */
		StatInit(&stat);
		
		/* ILU driver gives approximate solution X */
		dgsisx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, &mem_usage, &stat, &info);
		cout << " equed     "<<equed << endl;

//		int panel_size = sp_ienv(1);
//		int relax      = sp_ienv(2);		
//		dgsitrf(&options, &A, relax, panel_size,etree, work, lwork, perm_c, perm_r,	&L, &U, &stat, &info);
		
		if (info > 0) {
			cout << "something went wrong with superlu \n";
			//exit(1);
		}

		printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
			   mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
		
		if ( options.PivotGrowth == YES )
			printf("Recip. pivot growth = %e\n", rpg);
		if ( options.ConditionNumber == YES )
			printf("Recip. condition number = %e\n", rcond);
		
		if ( options.PrintStat ) StatPrint(&stat);
		
		b=res_vec;
		sp_dgemv("N", -1.0, &A, du.data(), 1, 1.0, b.data(), 1);//b=-Ax+b
		
		FLT ilu_res=nrm2(n, b);
		cout << endl << endl << "ILU residual "<< ilu_res << endl << endl;
		if (ilu_res > 1.0e-2) du = 0.0;// sometimes ILU answer is too crappy to feed into gmres
		
		/* Flexible GMRES */ 
		t = SuperLU_timer_();
		info = fgmres(n, A, L, U, *perm_c, *perm_r, res_vec, du, gtol, im, itmax, stat);
		t = SuperLU_timer_()-t;
		cout << "GMRES total time " << t << " iterations " << itmax << endl;

		if (info > 0) {
			cout << "something went wrong with gmres \n";
			//exit(1);
		}
		
		b=res_vec;
		sp_dgemv("N", -1.0, &A, du.data(), 1, 1.0, b.data(), 1);//b=-Ax+b
		
		cout << endl << endl << "gmres residual "<< nrm2(n, b) << endl << endl;

		/* update solution */
		ug_vec -= du;

		/* send ug_vec to global solution */
		vec_to_ug();
				
		output("superlu", tecplot);


		cout << endl;
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
		
		StatFree(&stat);

	}
	

	
	SUPERLU_FREE (perm_r);
	SUPERLU_FREE (perm_c);
	SUPERLU_FREE (etree);
	SUPERLU_FREE (R);
	SUPERLU_FREE (C);

	return;
}

void tet_hp::superlu(){
	cout << "superlu called" <<endl;

    char           equed[1];
    SuperMatrix    A, L, U;
    SuperMatrix    B, X;
	trans_t  trans;
	int            *perm_r; /* row permutations from partial pivoting */
    int            *perm_c; /* column permutation vector */
    int            *etree;
    void           *work;
    int            info, lwork=0; /* superlu determines memory allocation */
	int			   n=size_sparse_matrix;
    double         *R, *C;
    double         *ferr, *berr;
    double         rpg, rcond;
	
    mem_usage_t    mem_usage;
    superlu_options_t options;
    SuperLUStat_t  stat;

	Array<double,1> du(n);

	FLT tol=1.0e-14;
	int max_newton_its = 50;
	bool compressed_column = true;
	
	perm_c = intMalloc(n);
	perm_r = intMalloc(n);
	etree = intMalloc(n);
	
	ferr = (double *) SUPERLU_MALLOC(1 * sizeof(double));
	berr = (double *) SUPERLU_MALLOC(1 * sizeof(double));
	
	R = (double *) SUPERLU_MALLOC(n * sizeof(double));	
	C = (double *) SUPERLU_MALLOC(n * sizeof(double));
	
	/* create super matrix A using compressed column storage */
	dCreate_CompCol_Matrix(&A, n, n, number_sparse_elements, sparse_val.data(), sparse_ind.data(), sparse_ptr.data(), SLU_NC, SLU_D, SLU_GE);
	
	dCreate_Dense_Matrix(&X, n, 1, du.data(), n, SLU_DN, SLU_D, SLU_GE);
	
	dCreate_Dense_Matrix(&B, n, 1, res_vec.data(), n, SLU_DN, SLU_D, SLU_GE);
	
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
		for (int j = 0; j < n; ++j){
			if(fabs(res_vec(j)) > max_resid)
				max_resid = fabs(res_vec(j));
		}
		cout << "L infinity norm of residual " << max_resid << endl;
		
		if (max_resid < tol) break;		
		
		/* Initialize the statistics variables. */
		StatInit(&stat);
		
		/* advanced superlu driver */
		dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C, &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,&mem_usage, &stat, &info);				
		
		if (info > 0) {
			cout << "something went wrong with superlu \n";
			exit(1);
		}
		
		printf("L\\U MB %.3f\ttotal MB needed %.3f\n", mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
		
		/* update solution */
		ug_vec -= du;
		
		/* send ug_vec to global solution */
		vec_to_ug();
		
		if ( options.PrintStat ) StatPrint(&stat);
		StatFree(&stat);
		
		output("superlu", tecplot);
		
		Destroy_SuperNode_Matrix(&L);
		Destroy_CompCol_Matrix(&U);
		cout << endl;
	}
	
	SUPERLU_FREE (etree);
	SUPERLU_FREE (R);
	SUPERLU_FREE (C);
	SUPERLU_FREE (ferr);
	SUPERLU_FREE (berr);
	
	SUPERLU_FREE (perm_r);
	SUPERLU_FREE (perm_c);
		
	
	return;
}

#endif