
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
#include <limits.h>


// #define DEBUG_JAC

void tri_hp::petsc_initialize(){
	int sm=basis::tri(log2p)->sm();
	int im=basis::tri(log2p)->im();
	int tm=basis::tri(log2p)->tm();
	int vdofs;
	
	if (mmovement != coupled_deformable) 
		vdofs = NV;
	else
		vdofs = ND+NV;

	jacobian_size = npnt*vdofs +(nseg*sm +ntri*im)*NV;

	/* count total degrees of freedom on boundaries */
	if (mmovement == coupled_deformable) {
		for(int i=0;i < nebd; ++i) {
			jacobian_size += hp_ebdry(i)->dofs(jacobian_size);
		}
		jacobian_size += helper->dofs(jacobian_size);
	}
	
	
	const int nblock = sim::blks.nblock;
	const int idnum = gbl->idnum;

	Array<int,1> sndsize(nblock), size(nblock);
	sndsize = 0;
	sndsize(idnum) = jacobian_size; 
	sim::blks.allreduce(sndsize.data(),size.data(),nblock,blocks::int_msg,blocks::sum);
	~sndsize;
	
	int total_size = 0;
	for(int i=0;i<nblock;++i) {
		total_size += size(i);
	}

	/* find number of non-zeros for each row */
	Array<int,1> nnzero(jacobian_size);
	int begin_seg = npnt*vdofs;
	int begin_tri = begin_seg+nseg*sm*NV;
	int begin_aux = begin_tri+ntri*im*NV;
		
	/* SELF CONNECTIONS */
	nnzero(Range(0,begin_seg-1)) = vdofs; 
	if (sm) {
		nnzero(Range(begin_seg,begin_tri-1)) = (2*vdofs +sm*NV);
		/* connections of high order isoparametric mappings */
	}
	if (im) nnzero(Range(begin_tri,begin_aux-1)) = 3*vdofs +(tm-3)*NV;
	
	/* edges and vertices connected to avertex */
	for(int i=0; i<npnt; ++i)	
		for(int n=0;n<vdofs;++n)
			nnzero(i*vdofs+n) += pnt(i).nnbor*(vdofs +sm*NV);
	
	for(int i=0; i<ntri; ++i) {	
		/* interior and opposing side mode for each vertex */
		for(int j=0;j<3;++j)
			for(int n=0;n<vdofs;++n)
				nnzero(tri(i).pnt(j)*vdofs+n) += (im+sm)*NV;
				
		
		/* interior modes,opposing side modes, and opposing vertex to each side */
		for(int j=0;j<3;++j)
			for(int m=0;m<sm;++m)
				for(int n=0;n<NV;++n)
					nnzero(begin_seg+tri(i).seg(j)*sm*NV+m*NV+n) += (im +2*sm)*NV +1*vdofs;
	}
			
	/* Connections to extra boundary unknowns & to other blocks */
	Array<int,1> nnzero_mpi(jacobian_size);
	nnzero_mpi = 0;
	for(int i=0;i<nebd;++i) 
		hp_ebdry(i)->non_sparse(nnzero,nnzero_mpi);

	helper->non_sparse(nnzero,nnzero_mpi);
		
	/* CREATE VECTORS & MATRICES */
	PetscErrorCode err;
	
	err = VecCreate(PETSC_COMM_WORLD,&petsc_u);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = VecSetSizes(petsc_u,jacobian_size,total_size);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = VecSetFromOptions(petsc_u);
	
	int high;
	VecGetOwnershipRange(petsc_u,&jacobian_start,&high);
	if (high-jacobian_start != jacobian_size) {
		*gbl->log << "weird size allocation for petsc vector" << std::endl;
		exit(1);
	}

	CHKERRABORT(MPI_COMM_WORLD,err);
	err = VecDuplicate(petsc_u,&petsc_f);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
#ifdef MY_SPARSE
	J.resize(jacobian_size,nnzero); //,jacobian_start);
	J_mpi.resize(jacobian_size,nnzero_mpi); //,jacobian_start);
	err = MatCreate(PETSC_COMM_SELF,&petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
#else
#ifndef MPISRC
	err = MatCreateSeqAIJ(PETSC_COMM_SELF,jacobian_size,jacobian_size,PETSC_NULL,nnzero.data(),&petsc_J);
#else
	err = MatCreateMPIAIJ(PETSC_COMM_WORLD,jacobian_size,jacobian_size,total_size,total_size,PETSC_NULL,nnzero.data(),PETSC_NULL,nnzero_mpi.data(),&petsc_J);
#endif
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = MatSetOption(petsc_J,MAT_SYMMETRIC,PETSC_FALSE); 
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = MatSetOption(petsc_J,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE); 
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	/* Uncomment for debugging */
	//	err = MatSetOption(petsc_J,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE); 
	//	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = MatSetFromOptions(petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
#endif

	/* Create linear solver context */
	err = KSPCreate(PETSC_COMM_WORLD,&ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = KSPSetFromOptions(ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = KSPGetPC(ksp,&pc);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	/* Want default to be LU */
//	err = PCSetType(pc, PCILU);     // LU
//	CHKERRABORT(MPI_COMM_WORLD,err);
	
//	err = PCFactorSetUseInPlace(pc);
//	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = PCSetFromOptions(pc);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	return;
}

void tri_hp::petsc_setup_preconditioner() {
	PetscLogDouble time1,time2;
	PetscErrorCode err;
	
	PetscGetTime(&time1);
	
#ifdef MY_SPARSE
	/* Not sure if I have to delete it each time or not */
	err = MatDestroy(petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	J._val = 0.0;
	J_mpi._val = 0.0;
	J_mpi.reset_columns();  // Fix me: stupid petsc!!!
	petsc_jacobian();
	J.check_for_unused_entries();
	J_mpi.check_for_unused_entries();

#ifndef MPISRC
	err = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,jacobian_size,jacobian_size,J._cpt.data(),J._col.data(),J._val.data(),&petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
#else
	int total_size;
	sim::blks.allreduce(&jacobian_size, &total_size, 1, blocks::int_msg, blocks::sum);
	err =  MatCreateMPIAIJWithSplitArrays(MPI_COMM_WORLD,jacobian_size,jacobian_size,total_size,total_size,
		J._cpt.data(),J._col.data(),J._val.data(),J_mpi._cpt.data(),J_mpi._col.data(),J_mpi._val.data(),&petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);

//	err =  MatCreateMPIAIJWithArrays(MPI_COMM_WORLD,jacobian_size,jacobian_size,total_size,total_size,J._cpt.data(),J._col.data(),J._val.data(),&petsc_J);
//	CHKERRABORT(MPI_COMM_WORLD,err);
#endif

//	err = MatSetOption(petsc_J,MAT_SYMMETRIC,PETSC_FALSE); 
//	CHKERRABORT(MPI_COMM_WORLD,err);
//	
//	err = MatSetOption(petsc_J,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE); 
//	CHKERRABORT(MPI_COMM_WORLD,err);

//	err = MatSetFromOptions(petsc_J);
//	CHKERRABORT(MPI_COMM_WORLD,err);

#else
	/* insert values into jacobian matrix J */		
	MatZeroEntries(petsc_J);
	petsc_jacobian();		
	err = MatSetOption(petsc_J,MAT_NEW_NONZERO_LOCATIONS,PETSC_FALSE);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = MatSetOption(petsc_J,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
//		double rtol=1.0e-12; // relative tolerance
//		double atol=MAX(petsc_norm*1.0e-3,1.0e-15);// absolute tolerance
//		double dtol = 10000; // divergence tolerance
//		int maxits = 10000; // maximum iterations
//		
//		KSPSetTolerances(ksp,rtol,atol,dtol,maxits);
	/* 
	 Set operators. Here the matrix that defines the linear system
	 also serves as the preconditioning matrix.
	 */
#endif
	PetscGetTime(&time2);
	*gbl->log << "jacobian made " << time2-time1 << " seconds" << endl;

#ifdef DEBUG_JAC
	test_jacobian();
#endif

	PetscGetTime(&time1);	 
	err = KSPSetOperators(ksp,petsc_J,petsc_J,SAME_NONZERO_PATTERN);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = KSPSetUp(ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	PetscGetTime(&time2);
	*gbl->log << "matrix inverted " << time2-time1 << " seconds" << endl;

	return;
	
}

void tri_hp::petsc_update() {
	Vec	resid,du;
	PetscInt its;
	PetscErrorCode err;
	
	PetscLogDouble time1,time2;
	err = VecDuplicate(petsc_f,&du);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = VecDuplicate(petsc_f,&resid);
	CHKERRABORT(MPI_COMM_WORLD,err);

	ug_to_petsc();
		
	petsc_rsdl();
			
	VecAssemblyBegin(petsc_f);
	VecAssemblyEnd(petsc_f);
	
	PetscGetTime(&time1);
	err = KSPSolve(ksp,petsc_f,du);
	CHKERRABORT(MPI_COMM_WORLD,err);

	KSPGetIterationNumber(ksp,&its);
	PetscGetTime(&time2);
	*gbl->log << "# iterations " << its << " solve time: " << time2-time1 << " seconds" << endl;
	
	helper->update(-1);
	//KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	
	/* update: u=u-J^-1*f=u-du */
	VecAXPY(petsc_u,-1.0,du);

	/* send petsc vector u back to ug */
	petsc_to_ug();

	err = VecDestroy(du);
	CHKERRABORT(MPI_COMM_WORLD,err);

	return;
}

/* temp fix can I input petsc vectors ? */
void tri_hp::petsc_to_ug(){
	PetscScalar *array;
	PetscErrorCode err;
	int ind = 0;
	
	err = VecGetArray(petsc_u,&array);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	if (mmovement != coupled_deformable) {
		for(int i = 0; i < npnt; ++i)
			for(int n = 0; n < NV; ++n)
				ug.v(i,n) = array[ind++];
	}
	else {
		for(int i = 0; i < npnt; ++i) {
			for(int n = 0; n < NV; ++n)
				ug.v(i,n) = array[ind++];
			for(int n = 0; n < ND; ++n)
				pnts(i)(n) = array[ind++];
		}
	}

	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tri(log2p)->sm(); ++m)
			for(int n = 0; n < NV; ++n)
				ug.s(i,m,n) = array[ind++];
	
	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tri(log2p)->im(); ++m)
			for(int n = 0; n < NV; ++n)
				ug.i(i,m,n) = array[ind++];		
				
	if (mmovement == coupled_deformable) {
		for(int i = 0; i < nebd; ++i) {
			if (!hp_ebdry(i)->curved || !hp_ebdry(i)->coupled) continue;
			
			for(int j = 0; j < ebdry(i)->nseg;++j) 
				for(int m = 0; m < basis::tri(log2p)->sm(); ++m) 
					for(int n = 0; n < ND; ++n) 
						hp_ebdry(i)->crds(j,m,n) = array[ind++];
		}
	}
	
	err = VecRestoreArray(petsc_u,&array);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	return;	
}

/* temp fix can I input petsc vectors ? */
void tri_hp::ug_to_petsc(){
	int ind = jacobian_start;
	
	if (mmovement != coupled_deformable) {
		for(int i = 0; i < npnt; ++i){
			for(int n = 0; n < NV; ++n){
				VecSetValues(petsc_u,1,&ind,&ug.v(i,n),INSERT_VALUES);
				++ind;
			}
		}	
	}
	else {
		for(int i = 0; i < npnt; ++i){
			for(int n = 0; n < NV; ++n){
				VecSetValues(petsc_u,1,&ind,&ug.v(i,n),INSERT_VALUES);
				++ind;
			}
			for(int n = 0; n < ND; ++n){
				VecSetValues(petsc_u,1,&ind,&pnts(i)(n),INSERT_VALUES);
				++ind;
			}			
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
	
	if (mmovement == coupled_deformable) {
		for(int i = 0; i < nebd; ++i) {
			if (!hp_ebdry(i)->curved || !hp_ebdry(i)->coupled) continue;
			
			for(int j = 0; j < ebdry(i)->nseg;++j) {
				for(int m = 0; m < basis::tri(log2p)->sm(); ++m) {
					for(int n = 0; n < ND; ++n) {
						VecSetValues(petsc_u,1,&ind,&hp_ebdry(i)->crds(j,m,n),INSERT_VALUES);
						++ind;					
					}
				}
			}
		}
	}
	
	return;	
}


void tri_hp::petsc_make_1D_rsdl_vector(Array<FLT,1> rv) {
	int ind = 0;	
	if (mmovement != coupled_deformable) {
		for (int i=0;i<npnt;++i)
			for(int n=0;n<NV;++n)
				rv(ind++) = gbl->res.v(i,n);
	}
	else {
		for (int i=0;i<npnt;++i) {
			for(int n=0;n<NV;++n)
				rv(ind++) = gbl->res.v(i,n);
			for(int n=0;n<ND;++n)
				rv(ind++) = dynamic_cast<r_tri_mesh::global *>(gbl)->res(i)(n);
		}
	}
	
	for (int i=0;i<nseg;++i) 
		for(int m=0;m<sm0;++m)
			for(int n=0;n<NV;++n)
				rv(ind++) = gbl->res.s(i,m,n);
	
	for (int i=0;i<ntri;++i) 
		for(int m=0;m<im0;++m)
			for(int n=0;n<NV;++n)
				rv(ind++) = gbl->res.i(i,m,n);
	
	if (sm0) {
		for (int i=0;i<nebd;++i) {
			if (hp_ebdry(i)->curved && hp_ebdry(i)->coupled)
				ind += hp_ebdry(i)->petsc_rsdl(rv(Range(ind,jacobian_size-1)));
		}
	}
}

void tri_hp::petsc_finalize(){

	/* 
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
	 */

	PetscErrorCode err;
	err = KSPDestroy(ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = VecDestroy(petsc_f);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = VecDestroy(petsc_u);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = MatDestroy(petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
		
	return;
}

#endif
