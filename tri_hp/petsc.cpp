
/* 
 Include "petscksp.h" so that we can use KSP solvers.  Note that this file
 automatically includes:
 petsc.h       - base PETSc routines   petscvec.h - vectors
 petscsys.h    - system routines       petscmat.h - matrices
 petscis.h     - index sets            petscksp.h - Krylov subspace methods
 petscviewer.h - viewers               petscpc.h  - preconditioners
 */
#ifdef petsc

#include "tri_hp.h"
#include "hp_boundary.h"

void tri_hp::petsc_initialize(){
	
	size_sparse_matrix = (npnt+nseg*basis::tri(log2p)->sm()+ntri*basis::tri(log2p)->im())*NV;

	PetscTruth mat_nonsymmetric;

	PetscInitializeNoArguments();

	/*
     Set flag if we are doing a nonsymmetric problem; the default is symmetric.
	 */
	
	PetscOptionsHasName(PETSC_NULL,"-mat_nonsym",&mat_nonsymmetric);

	/* 
     Create parallel matrix, specifying only its global dimensions.
     When using MatCreate(), the matrix format can be specified at
     runtime. Also, the parallel partitioning of the matrix is
     determined by PETSc at runtime.
	 */
	MatCreate(PETSC_COMM_WORLD,&petsc_J);
	//MatSetSizes(petsc_J,PETSC_DECIDE,PETSC_DECIDE,size_sparse_matrix,size_sparse_matrix);
	//MatSetFromOptions(petsc_J);

	/* finds number of columns for each row */
	find_sparse_bandwidth(); 

	/* 
	 Create parallel vectors.
	 - When using VecSetSizes(), we specify only the vector's global
	 dimension; the parallel partitioning is determined at runtime. 
	 - Note: We form 1 vector from scratch and then duplicate as needed.
	 */
	VecCreate(PETSC_COMM_WORLD,&petsc_u);
	VecSetSizes(petsc_u,PETSC_DECIDE,size_sparse_matrix);
	VecSetFromOptions(petsc_u);
	VecDuplicate(petsc_u,&petsc_f);

	/* 
	 Create linear solver context
	 */
	KSPCreate(PETSC_COMM_WORLD,&ksp);
	
	
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
//	PCSetType(pc,PCSPAI);
//	PCSORSetOmega(pc,1.3); // Set Omega in SOR
//	PCSetType(pc, PCJACOBI); // Jacobi
//	PCSetType(pc, PCASM); // Additive Schwarz method
	
	
//	KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

	/*
	Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
	These options will override those specified above as long as KSPSetFromOptions() is called _after_ any other customization routines.
	*/

	KSPSetFromOptions(ksp);
	
	/* set preconditioner side */
	//KSPSetPreconditionerSide(ksp,PC_RIGHT);
	
	/* count total degrees of freedom on boundaries */
	int ndofs = 0;
	for(int i=0;i < nebd; ++i)
		ndofs+=(ebdry(i)->nseg*basis::tet(log2p)->sm+ebdry(i)->npnt)*NV;
	
	dirichlet_rows.resize(ndofs);

	return;
}


void tri_hp::petsc_solve(){
	cout << "petsc called "<< endl;
	Vec            resid,du;          /* solution, update, function */
	PetscErrorCode ierr;
	PetscInt       its,max_newton_its;
	PetscScalar    petsc_norm;
	//PetscMPIInt    size,rank;
	max_newton_its = 10;
	FLT tol = 1.0e-12;
	//ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);//CHKERRQ(ierr);
	//ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);//CHKERRQ(ierr);
	
	PetscLogDouble time1,time2;
	
	VecDuplicate(petsc_f,&du);
	VecDuplicate(petsc_f,&resid);

	/* initialize u with ug */
	ug_to_petsc();
	
	for(int iter = 0; iter < max_newton_its; ++iter) {
		

		VecSet(petsc_f,0.0);
		if(iter > 0) MatZeroEntries(petsc_J);
		
		/* insert values into jacobian matrix J */		
		PetscGetTime(&time1);
		create_jacobian();
		PetscGetTime(&time2);
		
		cout << "jacobian made " << time2-time1 << " seconds" << endl;
		
		/* insert values into residual f (every processor will do this need to do it a different way) */
		create_rsdl();		
		
		/* apply neumman bc's */
		apply_neumman();
				
		MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);	
				
		row_counter = 0;

		/* apply dirichlet boundary conditions to sparse matrix and vector */
		for(int j = 0; j < nfbd; ++j)
			hp_fbdry(j)->apply_sparse_dirichlet();		
		
		MatZeroRows(petsc_J,row_counter,dirichlet_rows.data(),1.0);

		
		VecAssemblyBegin(petsc_f);
		VecAssemblyEnd(petsc_f);
		
		VecNorm(petsc_f,NORM_2,&petsc_norm);		
		cout << "norm of residual: " << petsc_norm << endl;
		
		if(petsc_norm < tol) break;
		
		MatSetOption(petsc_J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
		MatSetOption(petsc_J,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
		
		double rtol=1.0e-12; // relative tolerance
		double atol=MAX(petsc_norm*1.0e-3,1.0e-15);// absolute tolerance
		double dtol = 10000; // divergence tolerance
		int maxits = 10000; // maximum iterations
		
		KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
		/* 
		 Set operators. Here the matrix that defines the linear system
		 also serves as the preconditioning matrix.
		 */
		
		KSPSetOperators(ksp,petsc_J,petsc_J,SAME_NONZERO_PATTERN);// SAME_NONZERO_PATTERN
		
		/* 
		 Solve linear system.  Here we explicitly call KSPSetUp() for more
		 detailed performance monitoring of certain preconditioners, such
		 as ICC and ILU.  This call is optional, as KSPSetUp() will
		 automatically be called within KSPSolve() if it hasn't been
		 called already.
		 */
		
		KSPSetUp(ksp);

		PetscGetTime(&time1);
		KSPSolve(ksp,petsc_f,du);
		PetscGetTime(&time2);
	
		MatMult(petsc_J,du,resid);		
		VecAXPY(resid,-1.0,petsc_f);

		VecNorm(resid,NORM_2,&petsc_norm);
		KSPGetIterationNumber(ksp,&its);

		cout << "linear system residual: " << petsc_norm << endl;
		cout << "iterations: " << its << endl;
		cout << "solve time: " << time2-time1 << " seconds" << endl;
		
		//KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
		
		/* update: u=u-J^-1*f=u-du */
		VecAXPY(petsc_u,-1.0,du);

		/* send petsc vector u back to ug */
		petsc_to_ug();

		output("petsc", tecplot);
		
		//MatSetOption(petsc_J,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
		//MatSetOption(petsc_J,MAT_NO_NEW_NONZERO_LOCATIONS,PETSC_TRUE);
		
//		MatSetOption(petsc_J,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);
	}
	
	
	ierr = VecDestroy(du);//CHKERRQ(ierr);

	return;
}

/* temp fix can I input petsc vectors ? */
void tri_hp::petsc_to_ug(){
	PetscScalar *array;
	int ind = 0;
	ierr = VecGetArray(petsc_u,&array);

	for(int i = 0; i < npnt; ++i)
		for(int n = 0; n < NV; ++n)
			ug.v(i,n) = array[ind++];
	
	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tri(log2p)->sm(); ++m)
			for(int n = 0; n < NV; ++n)
				ug.s(i,m,n) = array[ind++];
	
	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tri(log2p)->im(); ++m)
			for(int n = 0; n < NV; ++n)
				ug.i(i,m,n) = array[ind++];		
	
	
	VecRestoreArray(petsc_u,&array);
	return;	
}

/* temp fix can I input petsc vectors ? */
void tet_hp::ug_to_petsc(){
	int ind = 0;
	
	for(int i = 0; i < npnt; ++i){
		for(int n = 0; n < NV; ++n){
			VecSetValues(petsc_u,1,&ind,&ug.v(i,n),INSERT_VALUES);
			++ind;
		}
	}
	
	for(int i = 0; i < nseg; ++i){
		for(int m = 0; m < basis::tri(log2p)->sm(); ++m){
			for(int n = 0; n < NV; ++n){
				VecSetValues(petsc_u,1,&ind,&ug.s(i,m,n),INSERT_VALUES);
				++ind;
			}
		}
	}
	
	for(int i = 0; i < ntri; ++i){
		for(int m = 0; m < basis::tri(log2p)->im(); ++m){
			for(int n = 0; n < NV; ++n){
				VecSetValues(petsc_u,1,&ind,&ug.i(i,m,n),INSERT_VALUES);
				++ind;
			}
		}
	}
	

	
	return;	
}

void tri_hp::petsc_finalize(){


	/* 
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
	 */
	KSPDestroy(ksp);
	VecDestroy(petsc_f);
	VecDestroy(petsc_u);
	MatDestroy(petsc_J);
	
	PetscFinalize();
	
	return;
}

#endif