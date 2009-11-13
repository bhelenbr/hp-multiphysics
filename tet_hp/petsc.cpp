
/* 
 Include "petscksp.h" so that we can use KSP solvers.  Note that this file
 automatically includes:
 petsc.h       - base PETSc routines   petscvec.h - vectors
 petscsys.h    - system routines       petscmat.h - matrices
 petscis.h     - index sets            petscksp.h - Krylov subspace methods
 petscviewer.h - viewers               petscpc.h  - preconditioners
 */
#ifdef petsc

#include "tet_hp.h"
#include "hp_boundary.h"

void tet_hp::petsc_initialize(){
	
	size_sparse_matrix = (npnt+nseg*em0+ntri*fm0+ntet*im0)*NV;
	
	PetscErrorCode ierr;
	PetscTruth mat_nonsymmetric;

//	PetscInitialize(&argc,&args,(char *)0,help);
//	PetscInitialize(&argc,&args,(char *)0);
//	PetscInitialize();
	PetscInitializeNoArguments();

	/* 
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.
	 */
	//ierr = MatCreate(PETSC_COMM_WORLD,&petsc_J);//CHKERRQ(ierr);
	//ierr = MatSetSizes(petsc_J,PETSC_DECIDE,PETSC_DECIDE,size_sparse_matrix,size_sparse_matrix);//CHKERRQ(ierr);
	
	/* or do this */
	find_sparse_bandwidth();
//	ierr = MatSetFromOptions(petsc_J);//CHKERRQ(ierr);

	/* 
	 Create parallel vectors.
	 - When using VecSetSizes(), we specify only the vector's global
	 dimension; the parallel partitioning is determined at runtime. 
	 - Note: We form 1 vector from scratch and then duplicate as needed.
	 */
	ierr = VecCreate(PETSC_COMM_WORLD,&petsc_u);//CHKERRQ(ierr);
	ierr = VecSetSizes(petsc_u,PETSC_DECIDE,size_sparse_matrix);//CHKERRQ(ierr);
	ierr = VecSetFromOptions(petsc_u);//CHKERRQ(ierr);
	ierr = VecDuplicate(petsc_u,&petsc_f);//CHKERRQ(ierr);

	/*
     Set flag if we are doing a nonsymmetric problem; the default is symmetric.
	 */
	
	ierr = PetscOptionsHasName(PETSC_NULL,"-mat_nonsym",&mat_nonsymmetric);//CHKERRQ(ierr);

	/* 
	 Create linear solver context
	 */
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);//CHKERRQ(ierr);
	
	
	/* choose KSP type */
//	KSPSetType(ksp, KSPGMRES); // GMRES
//	KSPSetType(ksp, KSPCHEBYCHEV); // Chebychev
//	KSPSetType(ksp, KSPRICHARDSON);  // Richardson
//	KSPSetType(ksp, KSPLSQR);  // Least Squares
//	KSPSetType(ksp, KSPBICG);  // BiConjugate Gradient

	
	/*
	Set linear solver defaults for this problem (optional). - By extracting the KSP and PC contexts from the KSP context,
	we can then directly call any KSP and PC routines to set
	various options. - The following four statements are optional; all of these
	parameters could alternatively be specified at runtime via KSPSetFromOptions();
	*/
	KSPGetPC(ksp,&pc);
	
	/* choose preconditioner type */

//	PCSetType(pc, PCLU);     // LU
//	PCSetType(pc, PCILU);    // incomplete LU
//	PCSetType(pc, PCSOR);    // SOR
//	PCSORSetOmega(pc,1.3); // Set Omega in SOR
//	PCSetType(pc, PCJACOBI); // Jacobi
//	PCSetType(pc, PCASM); // Additive Schwarz method
	
	
//	KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

	/*
	Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
	These options will override those specified above as long as KSPSetFromOptions() is called _after_ any other customization routines.
	*/

	ierr = KSPSetFromOptions(ksp);//CHKERRQ(ierr);
	
	return;
}


void tet_hp::petsc_solve(){
	cout << "petsc called "<< endl;
	Vec            resid,du;          /* solution, update, function */
	PetscErrorCode ierr;
	PetscInt       newton_its,its,max_newton_its;
	PetscScalar    petsc_norm;
	//PetscMPIInt    size,rank;
	max_newton_its = 10;
	FLT tol = 1.0e-12;
	//ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);//CHKERRQ(ierr);
	//ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);//CHKERRQ(ierr);
	
	
	VecDuplicate(petsc_f,&du);
	VecDuplicate(petsc_f,&resid);

	/* initialize u with ug */
	ug_to_petsc();
	
	// -ksp_max_it <its>
	//ierr = KSPSetTolerances(KSP ksp,double rtol,double atol,double dtol,int maxits);//CHKERRQ(ierr);

	for(int iter = 0; iter < max_newton_its; ++iter) {
		
		/* creates jacobian and residual with BC's applied */
		create_jacobian_residual();
		
		VecNorm(petsc_f,NORM_2,&petsc_norm);		
		cout << "norm of residual: " << petsc_norm << endl;
		
		if(petsc_norm < tol) break;
		
		/* 
		 Set operators. Here the matrix that defines the linear system
		 also serves as the preconditioning matrix.
		 */
		
		KSPSetOperators(ksp,petsc_J,petsc_J,DIFFERENT_NONZERO_PATTERN);// SAME_NONZERO_PATTERN
		
		/* 
		 Solve linear system.  Here we explicitly call KSPSetUp() for more
		 detailed performance monitoring of certain preconditioners, such
		 as ICC and ILU.  This call is optional, as KSPSetUp() will
		 automatically be called within KSPSolve() if it hasn't been
		 called already.
		 */
		
		KSPSetUp(ksp);
		KSPSolve(ksp,petsc_f,du);
	
		MatMult(petsc_J,du,resid);		
		VecAXPY(resid,-1.0,petsc_f);

		VecNorm(resid,NORM_2,&petsc_norm);
		KSPGetIterationNumber(ksp,&its);

		cout << "linear system residual: " << petsc_norm << endl;
		cout << "iterations: " << its << endl;
		
		KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
		
		/* update: u=u-J^-1*f=u-du */
		VecAXPY(petsc_u,-1.0,du);

		/* send petsc vector u back to ug */
		petsc_to_ug();

		output("petsc", tecplot);
		
		//MatSetOption(petsc_J,MAT_NO_NEW_NONZERO_LOCATIONS,PETSC_TRUE);
	}
	
	
	ierr = VecDestroy(du);//CHKERRQ(ierr);

	return;
}

/* temp fix can I input petsc vectors ? */
void tet_hp::petsc_to_ug(){
	PetscErrorCode ierr;
	PetscScalar *array;
	int ind = 0;
	ierr = VecGetArray(petsc_u,&array);

	for(int i = 0; i < npnt; ++i)
		for(int n = 0; n < NV; ++n)
			ug.v(i,n) = array[ind++];
	
	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tet(log2p).em; ++m)
			for(int n = 0; n < NV; ++n)
				ug.e(i,m,n) = array[ind++];
	
	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tet(log2p).fm; ++m)
			for(int n = 0; n < NV; ++n)
				ug.f(i,m,n) = array[ind++];		
	
	for(int i = 0; i < ntet; ++i)
		for(int m = 0; m < basis::tet(log2p).im; ++m)
			for(int n = 0; n < NV; ++n)
				ug.i(i,m,n) = array[ind++];	
	
	ierr = VecRestoreArray(petsc_u,&array);
	return;	
}

/* temp fix can I input petsc vectors ? */
void tet_hp::ug_to_petsc(){
	PetscErrorCode ierr;
	int ind = 0;
	
	for(int i = 0; i < npnt; ++i){
		for(int n = 0; n < NV; ++n){
			ierr = VecSetValues(petsc_u,1,&ind,&ug.v(i,n),INSERT_VALUES);
			++ind;
		}
	}
	
	for(int i = 0; i < nseg; ++i){
		for(int m = 0; m < basis::tet(log2p).em; ++m){
			for(int n = 0; n < NV; ++n){
				ierr = VecSetValues(petsc_u,1,&ind,&ug.e(i,m,n),INSERT_VALUES);
				++ind;
			}
		}
	}
	
	for(int i = 0; i < ntri; ++i){
		for(int m = 0; m < basis::tet(log2p).fm; ++m){
			for(int n = 0; n < NV; ++n){
				ierr = VecSetValues(petsc_u,1,&ind,&ug.f(i,m,n),INSERT_VALUES);
				++ind;
			}
		}
	}
	
	for(int i = 0; i < ntet; ++i){
		for(int m = 0; m < basis::tet(log2p).im; ++m){
			for(int n = 0; n < NV; ++n){
				ierr = VecSetValues(petsc_u,1,&ind,&ug.i(i,m,n),INSERT_VALUES);
				++ind;	
			}
		}
	}
	
	return;	
}

void tet_hp::petsc_finalize(){

	PetscErrorCode ierr;

	/* 
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
	 */
	ierr = KSPDestroy(ksp);//CHKERRQ(ierr);
	ierr = VecDestroy(petsc_f);//CHKERRQ(ierr);
	ierr = VecDestroy(petsc_u);//CHKERRQ(ierr);
	ierr = MatDestroy(petsc_J);//CHKERRQ(ierr);
	
	/*
     Indicate to PETSc profiling that we're concluding the second stage 
	 */
	//ierr = PetscLogStagePop();//CHKERRQ(ierr);
	
	ierr = PetscFinalize();//CHKERRQ(ierr);
	
	return;
}

#endif