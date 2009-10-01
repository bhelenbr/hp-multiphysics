
/* 
 Include "petscksp.h" so that we can use KSP solvers.  Note that this file
 automatically includes:
 petsc.h       - base PETSc routines   petscvec.h - vectors
 petscsys.h    - system routines       petscmat.h - matrices
 petscis.h     - index sets            petscksp.h - Krylov subspace methods
 petscviewer.h - viewers               petscpc.h  - preconditioners
 */
#include <petscksp.h>
#include "tet_hp.h"

void petsc_initialize(){
	
	PetscInt m=size_sparse_matrix;
	PetscErrorCode ierr;
	PetscTruth mat_nonsymmetric;

	PetscInitialize(&argc,&args,(char *)0,help);

	/* 
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.
	 */
	ierr = MatCreate(PETSC_COMM_WORLD,&petsc_J);CHKERRQ(ierr);
	ierr = MatSetSizes(petsc_J,PETSC_DECIDE,PETSC_DECIDE,m,m);CHKERRQ(ierr);
	ierr = MatSetFromOptions(petsc_J);CHKERRQ(ierr);
	
	/* 
	 Create parallel vectors.
	 - When using VecSetSizes(), we specify only the vector's global
	 dimension; the parallel partitioning is determined at runtime. 
	 - Note: We form 1 vector from scratch and then duplicate as needed.
	 */
	ierr = VecCreate(PETSC_COMM_WORLD,&petsc_u);CHKERRQ(ierr);
	ierr = VecSetSizes(petsc_u,PETSC_DECIDE,m);CHKERRQ(ierr);
	ierr = VecSetFromOptions(petsc_u);CHKERRQ(ierr);
	ierr = VecDuplicate(petsc_u,&petsc_f);CHKERRQ(ierr);
	
	ierr = PetscOptionsHasName(PETSC_NULL,"-mat_nonsym",&mat_nonsymmetric);CHKERRQ(ierr);

	/* choose KSP type */
	KSPSetType(KSP ksp,KSPType KSPGMRES);
	
	/* 
	 Set runtime options (e.g., -ksp_type <type> -pc_type <type>)
	 */
	
	/* move this maybe? */
	ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
	
	/* 
	 Create linear solver context
	 */
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	
	/* 
	 Set operators. Here the matrix that defines the linear system
	 also serves as the preconditioning matrix.
	 */
	
	/* figure out preconditioners */
	ierr = KSPSetOperators(ksp,petsc_J,petsc_J,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
	
	return;
}

void petsc_solve(){
	
	Vec            du;          /* solution, update, function */
	PetscReal      norm=1.0;             /* norm of solution error */
	PetscErrorCode ierr;
	PetscInt       its,newton_its,m=size_sparse_matrix,max_newton_its;
	//PetscMPIInt    size,rank;

	max_newton_its = 1000;
	//ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	//ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	
	/*
     Set flag if we are doing a nonsymmetric problem; the default is symmetric.
	 */
		
	ierr = VecDuplicate(petsc_f,&du);CHKERRQ(ierr);
	
	/* initialize u with ug */
	ug_to_petsc();
	
	
	// -ksp_max_it <its>
	//ierr = KSPSetTolerances(KSP ksp,double rtol,double atol,double dtol,int maxits);CHKERRQ(ierr);
	
	newton_its = 0;
	while(norm > 1e-12){
		
		/* insert values into jacobian matrix J (every processor will do this need to do it a different way) */
		create_jacobian();
		
		/* 
		 Assemble matrix, using the 2-step process:
		 MatAssemblyBegin(), MatAssemblyEnd()
		 Computations can be done while messages are in transition
		 by placing code between these two statements.
		 */
		ierr = MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		

		/* insert values into residual f (every processor will do this need to do it a different way) */
		create_rsdl();		
		
		/* 
		 Assemble vector, using the 2-step process:
		 VecAssemblyBegin(), VecAssemblyEnd()
		 Computations can be done while messages are in transition,
		 by placing code between these two statements.
		 */
		ierr = VecAssemblyBegin(petsc_f);CHKERRQ(ierr);
		ierr = VecAssemblyEnd(petsc_f);CHKERRQ(ierr);
		
	
		/* 
		 Solve linear system.  Here we explicitly call KSPSetUp() for more
		 detailed performance monitoring of certain preconditioners, such
		 as ICC and ILU.  This call is optional, as KSPSetUp() will
		 automatically be called within KSPSolve() if it hasn't been
		 called already.
		 */
		
		/* because nonlinear fix number of iterations? */
		ierr = KSPSetUp(ksp);CHKERRQ(ierr);
		ierr = KSPSolve(ksp,petsc_f,du);CHKERRQ(ierr);
		
		/* u=u-J^-1*f=u-du */
		ierr = VecAXPY(petsc_u,-1.0,du);CHKERRQ(ierr);

		/* send petsc vector u back to ug */
		petsc_to_ug();

		++newton_its;
		
		/* check the error */

		ierr = VecNorm(du,NORM_2,&norm);CHKERRQ(ierr);
		ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Newton iteration %D, Norm of error %A, ksp Iterations %D\n",newton_its,norm,its);CHKERRQ(ierr);
		
		if (newton_its > max_newton_its ) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Max Newton iterations exceeded: break \n");CHKERRQ(ierr);
			break;
		}
	}
	
	ierr = VecDestroy(du);CHKERRQ(ierr);

	return;
}

/* temp fix can I input petsc vectors ? */
void tet_hp::petsc_to_ug(){
	
	int ind = 0;
	
	for(int i = 0; i < npnt; ++i)
		for(int n = 0; n < NV; ++n)
			ug.v(i,n) = petsc_u(ind++);
	
	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tet(log2p).em; ++m)
			for(int n = 0; n < NV; ++n)
				ug.e(i,m,n) = petsc_u(ind++);
	
	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tet(log2p).fm; ++m)
			for(int n = 0; n < NV; ++n)
				ug.f(i,m,n) = petsc_u(ind++);		
	
	for(int i = 0; i < ntet; ++i)
		for(int m = 0; m < basis::tet(log2p).im; ++m)
			for(int n = 0; n < NV; ++n)
				ug.i(i,m,n) = petsc_u(ind++);	
	
	return;	
}

/* temp fix can I input petsc vectors ? */
void tet_hp::ug_to_petsc(){
	
	int ind = 0;
	
	for(int i = 0; i < npnt; ++i)
		for(int n = 0; n < NV; ++n)
			petsc_u(ind++) = ug.v(i,n);
	
	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tet(log2p).em; ++m)
			for(int n = 0; n < NV; ++n)
				petsc_u(ind++) = ug.e(i,m,n);
	
	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tet(log2p).fm; ++m)
			for(int n = 0; n < NV; ++n)
				petsc_u(ind++) = ug.f(i,m,n);;		
	
	for(int i = 0; i < ntet; ++i)
		for(int m = 0; m < basis::tet(log2p).im; ++m)
			for(int n = 0; n < NV; ++n)
				petsc_u(ind++) = ug.i(i,m,n);	
	
	return;	
}

void petsc_finalize(){

	PetscErrorCode ierr;

	/* 
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
	 */
	ierr = KSPDestroy(ksp);CHKERRQ(ierr);
	ierr = VecDestroy(petsc_f);CHKERRQ(ierr);
	ierr = VecDestroy(petsc_u);CHKERRQ(ierr);
	ierr = MatDestroy(petsc_J);CHKERRQ(ierr);
	
	/*
     Indicate to PETSc profiling that we're concluding the second stage 
	 */
	//ierr = PetscLogStagePop();CHKERRQ(ierr);
	
	ierr = PetscFinalize();CHKERRQ(ierr);
	
	return;
}
