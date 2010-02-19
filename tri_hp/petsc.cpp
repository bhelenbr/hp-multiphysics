
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
	int isodofs = 0;
	if (mmovement == coupled_deformable) {
		for(int i=0;i < nebd; ++i) {
			if (hp_ebdry(i)->curved && hp_ebdry(i)->coupled) {
				isodofs += ebdry(i)->nseg*sm*ND;
			}
		}
		jacobian_size += isodofs;
	}
	
	/* find number of non-zeros for each row */
	Array<int,1> nnzero(jacobian_size);
	int begin_seg = npnt*vdofs;
	int begin_tri = begin_seg+nseg*sm*NV;
	int begin_iso = begin_tri +ntri*im*NV;
	
	/* SELF CONNECTIONS */
	nnzero(Range(0,begin_seg-1)) = vdofs; 
	if (sm) {
		nnzero(Range(begin_seg,begin_tri-1)) = (2 +sm)*NV;
		/* connections of high order isoparametric mappings */
		if (mmovement == coupled_deformable && jacobian_size > begin_iso) nnzero(Range(begin_iso,jacobian_size-1)) = vdofs*(sm+3) +NV*(im+2*sm);
	}
	if (im) nnzero(Range(begin_tri,jacobian_size-1)) = tm*NV;
	
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
	
	/* Connections to isoparametric coordinates */
	if (sm) {
		for(int i=0;i<nebd;++i) {
			if (!hp_ebdry(i)->curved || !hp_ebdry(i)->coupled) continue;
			
			for(int j=0;j<ebdry(i)->nseg;++j) {

				int tind = seg(ebdry(i)->seg(j)).tri(0);
				if (im) nnzero(Range(begin_tri +tind*im*NV,begin_tri +(tind+1)*im*NV-1)) += ND*sm;
				
				for(int j=0;j<3;++j) {
					int sind = tri(tind).seg(j);
					nnzero(Range(begin_seg+sind*NV*sm,begin_seg+(sind+1)*NV*sm-1)) += ND*sm;
					int pind = tri(tind).pnt(j);
					nnzero(Range(pind*vdofs,(pind+1)*vdofs-1))+= ND*sm;
				}
			}
		}
	}
	
	/* FIXME: NEED TO DO PARALLEL COMMUNICATION HERE TO PASS DATA ABOUT NNZ'S ON PARALLEL BOUNDARIES (create nmpizero.data())*/
	
	PetscErrorCode err = MatCreateSeqAIJ(PETSC_COMM_SELF,jacobian_size,jacobian_size,PETSC_NULL,nnzero.data(),&petsc_J);
	//MatCreateMPIAIJ(PETSC_COMM_WORLD,jacobian_size,jacobian_size,PETSC_DECIDE,PETSC_DECIDE,PETSC_NULL,nnzero.data(),PETSC_NULL,nmpizero.data(),&petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = MatSetOption(petsc_J,MAT_SYMMETRIC,PETSC_FALSE); 
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = MatSetFromOptions(petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);

	/* 
	 Create parallel vectors.
	 - When using VecSetSizes(), we specify only the vector's global
	 dimension; the parallel partitioning is determined at runtime. 
	 - Note: We form 1 vector from scratch and then duplicate as needed.
	 */
	err = VecCreate(PETSC_COMM_WORLD,&petsc_u);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = VecSetSizes(petsc_u,jacobian_size,PETSC_DECIDE);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = VecSetFromOptions(petsc_u);
	
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = VecDuplicate(petsc_u,&petsc_f);
	CHKERRABORT(MPI_COMM_WORLD,err);

	/* Create linear solver context */
	err = KSPCreate(PETSC_COMM_WORLD,&ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = KSPSetFromOptions(ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = KSPGetPC(ksp,&pc);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	/* Want default to be LU */
	err = PCSetType(pc, PCLU);     // LU
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = PCFactorSetUseInPlace(pc);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = KSPSetFromOptions(ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	return;
}

void tri_hp::petsc_setup_preconditioner() {
	PetscLogDouble time1,time2;
	PetscErrorCode err;
	
	/* insert values into jacobian matrix J */		
	PetscGetTime(&time1);
	MatZeroEntries(petsc_J);
	petsc_jacobian();
	PetscGetTime(&time2);
	
	*gbl->log << "jacobian made " << time2-time1 << " seconds" << endl;
		
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
	
	err = KSPSetOperators(ksp,petsc_J,petsc_J,SAME_NONZERO_PATTERN);// SAME_NONZERO_PATTERN
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	/* 
	 Solve linear system.  Here we explicitly call KSPSetUp() for more
	 detailed performance monitoring of certain preconditioners, such
	 as ICC and ILU.  This call is optional, as KSPSetUp() will
	 automatically be called within KSPSolve() if it hasn't been
	 called already.
	 */
	PetscGetTime(&time1);
	err = KSPSetUp(ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	PetscGetTime(&time2);
	
	*gbl->log << "matrix inverted " << time2-time1 << " seconds" << endl;

	
	//MatView(petsc_J,0);

	return;
	
}

void tri_hp::petsc_update() {
	Vec            resid,du;          /* solution, update, function */
	PetscErrorCode err;
	PetscInt       its,max_newton_its;
	PetscScalar    petsc_norm;
	//PetscMPIInt    size,rank;
	max_newton_its = 20;
	//err = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);//CHKERRABORT(MPI_COMM_WORLD,err);
	//err = MPI_Comm_size(PETSC_COMM_WORLD,&size);//CHKERRABORT(MPI_COMM_WORLD,err);
	
	PetscLogDouble time1,time2;
	err = VecDuplicate(petsc_f,&du);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = VecDuplicate(petsc_f,&resid);
	CHKERRABORT(MPI_COMM_WORLD,err);

	/* initialize u with ug */
	ug_to_petsc();
		
	/* insert values into residual f (every processor will do this need to do it a different way) */
	petsc_rsdl();
			
	VecAssemblyBegin(petsc_f);
	VecAssemblyEnd(petsc_f);
	
	PetscGetTime(&time1);
	err = KSPSolve(ksp,petsc_f,du);
	CHKERRABORT(MPI_COMM_WORLD,err);
	PetscGetTime(&time2);

	MatMult(petsc_J,du,resid);
	VecAXPY(resid,-1.0,petsc_f);		
	VecNorm(resid,NORM_2,&petsc_norm);
	KSPGetIterationNumber(ksp,&its);

	*gbl->log << "#linear system residual: " << petsc_norm << " iterations " << its << " solve time: " << time2-time1 << " seconds" << endl;
	
	//KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	
	/* update: u=u-J^-1*f=u-du */
	VecAXPY(petsc_u,-1.0,du);

	/* send petsc vector u back to ug */
	petsc_to_ug();
	
	//MatSetOption(petsc_J,MAT_KEEP_ZEROED_ROWS,PETSC_TRUE);
	//MatSetOption(petsc_J,MAT_NO_NEW_NONZERO_LOCATIONS,PETSC_TRUE);
	
//		MatSetOption(petsc_J,MAT_NEW_NONZERO_LOCATION_ERR,PETSC_TRUE);	
	
	err = VecDestroy(du); CHKERRABORT(MPI_COMM_WORLD,err);

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
	int ind = 0;
	
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
