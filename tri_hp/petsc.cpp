
/*
 *  petsc.cpp
 *  tri_hp
 *
 *  Created by michael brazell on 9/14/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

/*
 Include "petscksp.h" so that we can use KSP solvers.  Note that this file
 automatically includes:
 petsc.h       - base PETSc routines   petscvec.h - vectors
 petscsys.h    - system routines       petscmat.h - matrices
 petscis.h     - index sets            petscksp.h - Krylov subspace methods
 petscviewer.h - viewers               petscpc.h  - preconditioners
 */

#ifdef petsc

#include <limits.h>
#include <petsctime.h>
#include <petscviewer.h>

#include "tri_hp.h"
#include "hp_boundary.h"
#include <r_tri_boundary.h>

#define DEBUG_TOL 1.0e-9
#define DEBUG_ABS_TOL 1.0e-9
#define DEBUG_REL_TOL 1.0e-2


void tri_hp::petsc_initialize(){
	int sm=basis::tri(log2p)->sm();
	int im=basis::tri(log2p)->im();
	int tm=basis::tri(log2p)->tm();
	const int vdofs = NV +(mmovement == tri_hp::coupled_deformable)*ND;


	/* Internal degrees of freedom */
	jacobian_size = npnt*vdofs +(nseg*sm +ntri*im)*NV;

	/* Count total degrees of freedom on boundaries */
	if (mmovement == coupled_deformable) {
		for(int i=0;i < nebd; ++i) {
			jacobian_size += hp_ebdry(i)->dofs(jacobian_size);
		}
		jacobian_size += helper->dofs(jacobian_size);
	}
	
	
	/* Assemble list of sizes for each block */
	const int nblock = sim::blks.nblock;
	const int idnum = gbl->idnum;

	Array<int,1> sndsize(nblock), size(nblock);
	sndsize = 0;
	sndsize(idnum) = jacobian_size; 
	sim::blks.allreduce(sndsize.data(),size.data(),nblock,blocks::int_msg,blocks::sum);
	~sndsize;
	
	/* Add up total size */
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
	
	/* Connections to extra boundary unknowns	*/	
	for(int i=0;i<nebd;++i) 
		hp_ebdry(i)->non_sparse(nnzero);
	
	for(int i=0;i<nvbd;++i) 
		hp_vbdry(i)->non_sparse(nnzero);
		
	/* Connections to other blocks */
	/* Jacobian rows for matching degrees of freedom are added together */
	/* so there must be enough space for the entire jacobian row */
	/* Both rows end up being basically the same except the diagonal entries */
	/* one each row are the local degrees of freedom which implicitly enforces */
	/* continuity of the degrees of freedom by: */
	/* a_ii u_l +rest of row = b */
	/* a_ii u_r +rest of row = b */
	/* a_ii (u_r -u_l) = 0 */
	Array<int,1> nnzero_mpi(jacobian_size);
	nnzero_mpi = 0;
	
	
	/* SEND COMMUNICATIONS TO ADJACENT MESHES */
	for(int last_phase = false, phase = 0; !last_phase; ++phase) {
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->non_sparse_snd(nnzero,nnzero_mpi);
		
		for(int i=0;i<nvbd;++i)
			hp_vbdry(i)->non_sparse_snd(nnzero,nnzero_mpi);
		
		for(int i=0;i<nebd;++i)
			ebdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
		
		for(int i=0;i<nvbd;++i)
			vbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
		
		pmsgpass(boundary::all_phased,phase,boundary::symmetric);
		
		last_phase = true;
		for(int i=0;i<nebd;++i) {
			last_phase &= ebdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
		}
		for(int i=0;i<nvbd;++i) {
			last_phase &= vbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
		}
		
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->non_sparse_rcv(phase,nnzero,nnzero_mpi);
		
		for(int i=0;i<nvbd;++i)
			hp_vbdry(i)->non_sparse_rcv(phase,nnzero,nnzero_mpi);
	}

	helper->non_sparse(nnzero,nnzero_mpi);
	
//	*gbl->log << nnzero << std::endl;
//	*gbl->log << nnzero_mpi << std::endl;
	
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
	err = VecDuplicate(petsc_u,&petsc_du);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
#ifdef MY_SPARSE
	J.resize(jacobian_size,nnzero); //,jacobian_start);
	J_mpi.resize(jacobian_size,nnzero_mpi); //,jacobian_start);
	err = MatCreate(PETSC_COMM_WORLD,&petsc_J);
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

void tri_hp::petsc_finalize(){
	
	/*
	 Free work space.  All PETSc objects should be destroyed when they
	 are no longer needed.
	 */
	
	PetscErrorCode err;
	err = KSPDestroy(&ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = VecDestroy(&petsc_f);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = VecDestroy(&petsc_u);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = VecDestroy(&petsc_du);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	err = MatDestroy(&petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
	return;
}


void tri_hp::petsc_setup_preconditioner() {
	PetscLogDouble time1,time2;
	PetscErrorCode err;
	
	PetscTime(&time1);
	
#ifdef MY_SPARSE

	/* Not sure if I have to delete it each time or not */
	err = MatDestroy(&petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
	petsc_jacobian();
	petsc_premultiply_jacobian();

	if (gbl->jac_debug)	{
		streamsize oldprecision = (*gbl->log).precision(2);
//		*gbl->log << "J:\n";
//		*gbl->log << J << std::endl;
//		*gbl->log << "J_mpi:\n";
//		*gbl->log << J_mpi << std::endl;
        std::string desc;
        
        /* A more useful way to output the matrix */
        const int vdofs = NV +(mmovement == tri_hp::coupled_deformable)*ND;
        for (int row = 0; row < jacobian_size; ++row) {
            local_index_to_mesh_descriptor(row,desc);
            *gbl->log << desc << ' ';
            int entries = J.nentries_for_row(row);
            for (int k=0; k<entries; ++k) {
                int col;
                FLT val;
                J.get_value_and_col(row, k, val, col);
                local_index_to_mesh_descriptor(col,desc);
                *gbl->log << '(' << desc << ',' << val << ") ";
            }
            *gbl->log << std::endl;
        }
        
        /* Output J_mpi */
        const int nblock = sim::blks.nblock;
        const int myid = sim::blks.myid;
        assert(sim::blks.myblock == 1);
        
        /* Get counts from other meshes */
        Array<int,1> sndsizes(4*nblock);
        sndsizes = 0;
        sndsizes(myid*4) = npnt;
        sndsizes(myid*4+1) = nseg;
        sndsizes(myid*4+2) = ntri;
        sndsizes(myid*4+3) = jacobian_size;
        
        Array<int,1> sizes(4*nblock);
        sim::blks.allreduce(sndsizes.data(),sizes.data(),4*nblock,blocks::int_msg,blocks::sum);
        ~sndsizes;
        
        for (int row = 0; row < jacobian_size; ++row) {
            local_index_to_mesh_descriptor(row,desc);
            *gbl->log << desc << ' ';
            int entries = J_mpi.nentries_for_row(row);
            for (int k=0; k<entries; ++k) {
                int col;
                FLT val;
                J_mpi.get_value_and_col(row, k, val, col);
                int block;
                for (block = 0; block < nblock; ++block) {
                    if (col < sizes(4*block+3))
                        goto foundv;
                    else
                        col -= sizes(4*block+3);
                }
                *gbl->log << "index out of range" << std::endl;
                sim::abort(__LINE__,__FILE__,gbl->log);
            foundv:
                int const npnt_mpi = sizes(4*block);
                int const nseg_mpi = sizes(4*block+1);
                int const ntri_mpi = sizes(4*block+2);
                int const jacobian_size_mpi = sizes(4*block+3);
                
                
                if (col < npnt_mpi*vdofs) {
                    *gbl->log << "(b" << block << ',' << 'v' << col/vdofs << ',' << col % vdofs << ',' << val << ") ";
                }
                else if ((col -= npnt_mpi*vdofs) < nseg_mpi*sm0*NV) {
                    *gbl->log << "(b" << block << ',' << 's' << col/(NV*sm0) << ',' << (col % (sm0*NV))/NV << ',' << col % NV << ',' << val << ") ";
                }
                else if ((col -= nseg_mpi*NV*sm0) < ntri_mpi*im0*NV) {
                    *gbl->log << "(b" << block << ',' << 'i' << col/(NV*im0) << ',' << (col % (im0*NV))/NV << ',' << col % NV << ',' << val << ") ";
                }
                else if ((col -= ntri_mpi*NV*im0) < jacobian_size_mpi){
                    *gbl->log << "(b" << block << ',' << 'e' << col/(ND*sm0) << ',' << (col % (sm0*ND))/ND << ',' << col % ND << ',' << val << ") ";
                }
                else {
                    *gbl->log << "index out of range" << std::endl;
                    sim::abort(__LINE__,__FILE__,gbl->log);
                }
            }
            *gbl->log << std::endl;
        }
        (*gbl->log).precision(oldprecision);
	}
	J.check_for_unused_entries(*gbl->log);
	J_mpi.check_for_unused_entries(*gbl->log);

#ifndef MPISRC
	err = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,jacobian_size,jacobian_size,J._cpt.data(),J._col.data(),J._val.data(),&petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
#else
	int total_size;
	sim::blks.allreduce(&jacobian_size, &total_size, 1, blocks::int_msg, blocks::sum);
	err =  MatCreateMPIAIJWithSplitArrays(PETSC_COMM_WORLD,jacobian_size,jacobian_size,total_size,total_size,
		J._cpt.data(),J._col.data(),J._val.data(),J_mpi._cpt.data(),J_mpi._col.data(),J_mpi._val.data(),&petsc_J);
	CHKERRABORT(MPI_COMM_WORLD,err);
	
    
//	err =  MatCreateMPIAIJWithArrays(PETSC_COMM_WORLD,jacobian_size,jacobian_size,total_size,total_size,J._cpt.data(),J._col.data(),J._val.data(),&petsc_J);
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
	PetscTime(&time2);
    *gbl->log << "#jacobian made " << time2-time1 << " seconds" << endl;

	if (gbl->jac_debug) {
        
		test_jacobian();
		sim::finalize(__LINE__,__FILE__,gbl->log);
	}

	PetscTime(&time1);	 
	err = KSPSetOperators(ksp,petsc_J,petsc_J);
	//err = KSPSetOperators(ksp,petsc_J,petsc_J,DIFFERENT_NONZERO_PATTERN);
	CHKERRABORT(MPI_COMM_WORLD,err);
	err = KSPSetUp(ksp);
	CHKERRABORT(MPI_COMM_WORLD,err);
	PetscTime(&time2);
    *gbl->log << "#matrix inverted " << time2-time1 << " seconds" << endl;

	return;
	
}


void tri_hp::petsc_jacobian() {
	
#ifdef MY_SPARSE
	J._val = 0.0;
	J.reset_columns();
	J_mpi._val = 0.0;
	J_mpi.reset_columns();
#endif
	
	int gindx;
	
	const int sm = basis::tri(log2p)->sm();
	const int im = basis::tri(log2p)->im();
	const int tm = basis::tri(log2p)->tm();
	
	int vdofs = NV;
	if (mmovement == coupled_deformable) vdofs += ND;
	int kn = 3*vdofs +(tm-3)*NV;
	Array<int,1> loc_to_glo(kn);
	Array<FLT,2> K(kn,kn);
	FLT sgn,msgn;
	Array<FLT,2> Rbar(NV,tm);
	
	/* DO ELEMENTS */
	for(int tind = 0; tind < ntri; ++tind) {
		
		int ind = 0;
		element_jacobian(tind, K);
		
		for (int m = 0; m < 3; ++m) {
			gindx = vdofs*tri(tind).pnt(m);
			for (int n = 0; n < vdofs; ++n)
				loc_to_glo(ind++) = gindx++;
		}
		
		/* EDGE MODES */
		if (sm) {
			for(int i = 0; i < 3; ++i) {
				gindx = npnt*vdofs +tri(tind).seg(i)*sm*NV;
				sgn = tri(tind).sgn(i);
				msgn = 1;
				for (int m = 0; m < sm; ++m) {
					for(int n = 0; n < NV; ++n) {
						for(int j = 0; j < kn; ++j) {
							K(ind,j) *= msgn;
							K(j,ind) *= msgn;
						}
						loc_to_glo(ind++) = gindx++;
					}
					msgn *= sgn;
				}
			}
		}
		
		/* INTERIOR	MODES */
		if (tm) {
			gindx = npnt*vdofs +nseg*sm*NV +tind*im*NV;
			for(int m = 0; m < im; ++m) {
				for(int n = 0; n < NV; ++n){
					loc_to_glo(ind++) = gindx++;
				}
			}
		}
		
#ifdef MY_SPARSE
		J.add_values(kn,loc_to_glo,kn,loc_to_glo,K);
#else
		MatSetValuesLocal(petsc_J,kn,loc_to_glo.data(),kn,loc_to_glo.data(),K.data(),ADD_VALUES);
#endif
	}
	
	/* APPLY ALL MOVING MESH B.C.'s FIRST */
	if (mmovement == coupled_deformable) {
		/* This mostly does nothing, except for angled boundarys */
#ifdef MY_SPARSE
		for(int i=0;i<nebd;++i)
			r_sbdry(i)->jacobian(J,J_mpi,vdofs);
		
		for(int i=0;i<nvbd;++i)
			r_vbdry(i)->jacobian(J,J_mpi,vdofs);
		
		for(int i=0;i<nebd;++i)
			r_sbdry(i)->jacobian_dirichlet(J,J_mpi,vdofs);
		
		for(int i=0;i<nvbd;++i)
			r_vbdry(i)->jacobian_dirichlet(J,J_mpi,vdofs);
#else
		for(int i=0;i<nebd;++i)
			r_sbdry(i)->jacobian(petsc_J,vdofs);
		
		for(int i=0;i<nvbd;++i)
			r_vbdry(i)->jacobian(petsc_J,vdofs);
		
		/* PETSC IS RETARDED */
		FLT zero = 0.0;
		for (int i=jacobian_start +npnt*(NV+ND)+sm*NV*nseg+ntri*im*NV;i<jacobian_start +jacobian_size;++i)
			MatSetValuesLocal(petsc_J,1,&i,1,&i,&zero,ADD_VALUES);
		
		MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);
		
		for(int i=0;i<nebd;++i)
			r_sbdry(i)->jacobian_dirichlet(petsc_J);
		
		for(int i=0;i<nvbd;++i)
			r_vbdry(i)->jacobian_dirichlet(petsc_J);
#endif
	}
	
	/* DO NEUMANN & COUPLED BOUNDARY CONDITIONS */
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->petsc_jacobian();
	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->petsc_jacobian();
	
#ifndef MY_SPARSE
	MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);
#endif
	
	
	/* SEND COMMUNICATIONS TO ADJACENT MESHES */
	for(int last_phase = false, phase = 0; !last_phase; ++phase) {
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->petsc_matchjacobian_snd();
		
		for(int i=0;i<nvbd;++i)
			hp_vbdry(i)->petsc_matchjacobian_snd();
		
		for(int i=0;i<nebd;++i)
			ebdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
		
		for(int i=0;i<nvbd;++i)
			vbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
		
		pmsgpass(boundary::all_phased,phase,boundary::symmetric);
		
		last_phase = true;
		for(int i=0;i<nebd;++i) {
			last_phase &= ebdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
		}
		for(int i=0;i<nvbd;++i) {
			last_phase &= vbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
		}
		
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->petsc_matchjacobian_rcv(phase);
		
		for(int i=0;i<nvbd;++i)
			hp_vbdry(i)->petsc_matchjacobian_rcv(phase);
		
	}
	
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->petsc_jacobian_dirichlet();
	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->petsc_jacobian_dirichlet();
	
	return;
}

void tri_hp::petsc_premultiply_jacobian() {

	for(int i=0;i<nebd;++i) {
		hp_ebdry(i)->petsc_premultiply_jacobian();
	}
	
	for(int i=0;i<nvbd;++i) {
		hp_vbdry(i)->petsc_premultiply_jacobian();
	}
}


void tri_hp::petsc_update() {
	PetscInt its;
	PetscErrorCode err;
	
	PetscLogDouble time1,time2;
	
	ug_to_petsc();
		
	petsc_rsdl();
			
	VecAssemblyBegin(petsc_f);
	VecAssemblyEnd(petsc_f);
	
	VecNorm(petsc_f, NORM_2, &max_residual);
	
	PetscTime(&time1);
	err = KSPSolve(ksp,petsc_f,petsc_du);
	CHKERRABORT(MPI_COMM_WORLD,err);
    
    if (gbl->rsdl_debug == 3) {
        petsc_output_vector(petsc_du);
        sim::finalize(__LINE__,__FILE__,&std::cerr);
    }
	
	double resmax2;
	VecNorm(petsc_du, NORM_2, &resmax2 );
	
	KSPGetIterationNumber(ksp,&its);
	PetscTime(&time2);
	*gbl->log << "# iterations " << its << " residual0 " << max_residual << " du " << resmax2 << " solve time: " << time2-time1 << " seconds" << endl;
	
	helper->update(-1);
	//KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	
	/* update: u=u-J^-1*f=u-du */
	VecAXPY(petsc_u,-under_relax,petsc_du);

	/* send petsc vector u back to ug */
	petsc_to_ug();

	for(int i=0;i<nebd;++i) {
		hp_ebdry(i)->update(-1);
	}
	
	for(int i=0;i<nvbd;++i) {
		hp_vbdry(i)->update(-1);
	}
	
//	std::cerr << "temporary exit" << std::endl;
//	sim::finalize(__LINE__,__FILE__,&std::cerr);

	return;
}

void tri_hp::petsc_rsdl() {
	
	rsdl();
	
	enforce_continuity(gbl->res, r_tri_mesh::gbl->res);
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet();
	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet();
	
	/* APPLY DIRCHLET B.C.S TO MODE */
	for(int m=0;m<basis::tri(log2p)->sm();++m)
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->sdirichlet(m);
	
	if (gbl->rsdl_debug == 1) {
		const int vdofs = NV +(mmovement == tri_hp::coupled_deformable)*ND;
		for(int i=0;i<npnt;++i) {
			*gbl->log << gbl->idprefix << " v: " << i << ' ';
			for(int n=0;n<NV;++n) {
				if (fabs(gbl->res.v(i,n)) > DEBUG_TOL) *gbl->log << gbl->res.v(i,n) << ' ';
				else *gbl->log << "0.0 ";
			}
			
			for(int n=0;n<vdofs-NV;++n) {
				if (fabs(r_tri_mesh::gbl->res(i)(n)) > DEBUG_TOL) *gbl->log << r_tri_mesh::gbl->res(i)(n) << ' ';
				else *gbl->log << "0.0 ";
			}
			*gbl->log << '\n';
		}
		
		for(int i=0;i<nseg;++i) {
			for(int m=0;m<basis::tri(log2p)->sm();++m) {
				*gbl->log << gbl->idprefix << " s: " << i << ' ';
				for(int n=0;n<NV;++n) {
					if (fabs(gbl->res.s(i,m,n)) > DEBUG_TOL) *gbl->log << gbl->res.s(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}
		
		
		for(int i=0;i<ntri;++i) {
			for(int m=0;m<basis::tri(log2p)->im();++m) {
				*gbl->log << gbl->idprefix << " i: " << i << ' ';
				for(int n=0;n<NV;++n) {
					if (fabs(gbl->res.i(i,m,n)) > DEBUG_TOL) *gbl->log << gbl->res.i(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}
		sim::finalize(__LINE__,__FILE__,gbl->log);
	}
	
	
	PetscScalar *array;
	VecGetArray(petsc_f,&array);
	Array<FLT,1> res(array, shape(jacobian_size), neverDeleteData);
	petsc_make_1D_rsdl_vector(res);
	VecRestoreArray(petsc_f, &array);
	
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
				rv(ind++) = r_tri_mesh::gbl->res(i)(n);
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
	
	
	for (int i=0;i<nebd;++i) {
		hp_ebdry(i)->petsc_make_1D_rsdl_vector(rv);
	}
	
	for (int i=0;i<nvbd;++i) {
		hp_vbdry(i)->petsc_make_1D_rsdl_vector(rv);
	}
	
	
	if (gbl->rsdl_debug == 2) {
		const int vdofs = NV +(mmovement == tri_hp::coupled_deformable)*ND;
		ind = 0;
		for(int i=0;i<npnt;++i) {
			*gbl->log << gbl->idprefix << " v: " << i << ' ';
			for(int n=0;n<vdofs;++n) {
				if (fabs(rv(ind)) > DEBUG_TOL) *gbl->log << rv(ind) << ' ';
				else *gbl->log << "0.0 ";
				++ind;
			}
			*gbl->log << '\n';
		}
		
		for(int i=0;i<nseg;++i) {
			for(int m=0;m<basis::tri(log2p)->sm();++m) {
				*gbl->log << gbl->idprefix << " s: " << i << ' ';
				for(int n=0;n<NV;++n) {
					if (fabs(rv(ind)) > DEBUG_TOL) *gbl->log << rv(ind) << ' ';
					else *gbl->log << "0.0 ";
					++ind;
				}
				*gbl->log << '\n';
			}
		}
		
		
		for(int i=0;i<ntri;++i) {
			for(int m=0;m<basis::tri(log2p)->im();++m) {
				*gbl->log << gbl->idprefix << " i: " << i << ' ';
				for(int n=0;n<NV;++n) {
					if (fabs(rv(ind)) > DEBUG_TOL) *gbl->log << rv(ind) << ' ';
					else *gbl->log << "0.0 ";
					++ind;
				}
				*gbl->log << '\n';
			}
		}
		
		for (int i = ind; i< rv.extent(firstDim); ++i) {
			if (fabs(rv(i)) > DEBUG_TOL) *gbl->log << rv(i) << ' ';
			else *gbl->log << "0.0 ";
		}
		
		sim::finalize(__LINE__,__FILE__,gbl->log);
	}
}

void tri_hp::petsc_output_vector(Vec petsc_vec) {
    PetscScalar *array;
    PetscErrorCode err;
    PetscInt local_size;
    int ind = 0;
    err = VecGetArray(petsc_vec,&array);
    CHKERRABORT(MPI_COMM_WORLD,err);
    err = VecGetLocalSize(petsc_vec,&local_size);
    CHKERRABORT(MPI_COMM_WORLD,err);
    
    const int vdofs = NV +(mmovement == tri_hp::coupled_deformable)*ND;
    ind = 0;
    for(int i=0;i<npnt;++i) {
        *gbl->log << gbl->idprefix << " v: " << i << ' ';
        for(int n=0;n<vdofs;++n) {
            if (fabs(array[ind]) > DEBUG_TOL) *gbl->log << array[ind] << ' ';
            else *gbl->log << "0.0 ";
            ++ind;
        }
        *gbl->log << '\n';
    }
    
    for(int i=0;i<nseg;++i) {
        for(int m=0;m<basis::tri(log2p)->sm();++m) {
            *gbl->log << gbl->idprefix << " s: " << i << ' ';
            for(int n=0;n<NV;++n) {
                if (fabs(array[ind]) > DEBUG_TOL) *gbl->log << array[ind] << ' ';
                else *gbl->log << "0.0 ";
                ++ind;
            }
            *gbl->log << '\n';
        }
    }
    
    
    for(int i=0;i<ntri;++i) {
        for(int m=0;m<basis::tri(log2p)->im();++m) {
            *gbl->log << gbl->idprefix << " i: " << i << ' ';
            for(int n=0;n<NV;++n) {
                if (fabs(array[ind]) > DEBUG_TOL) *gbl->log << array[ind] << ' ';
                else *gbl->log << "0.0 ";
                ++ind;
            }
            *gbl->log << '\n';
        }
    }
    
    for (int i = ind; i< local_size; ++i) {
        if (fabs(array[i]) > DEBUG_TOL) *gbl->log << array[i] << '\n';
        else *gbl->log << "0.0\n";
    }
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
	
	for (int i=0;i<nebd;++i) {
		ind += hp_ebdry(i)->petsc_to_ug(&array[ind]);
	}
	
	for (int i=0;i<nvbd;++i) {
		ind += hp_vbdry(i)->petsc_to_ug(&array[ind]);
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
	
	for(int i = 0; i < nebd; ++i)
		hp_ebdry(i)->ug_to_petsc(ind);
	
	for(int i = 0; i < nvbd; ++i)
		hp_vbdry(i)->ug_to_petsc(ind);
	
	return;
}


void tri_hp::enforce_continuity(vsi& ug, Array<TinyVector<FLT,ND>,1>& pnts) {
	int last_phase, mp_phase;
	
	if (mmovement == coupled_deformable) {
		/* Residual for r_mesh vertices */
		for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
			r_tri_mesh::pmsgload(boundary::all_phased,mp_phase, boundary::symmetric,(FLT *) pnts.data(),0,1,2);
			r_tri_mesh::pmsgpass(boundary::all_phased,mp_phase, boundary::symmetric);
			last_phase = true;
			last_phase &= r_tri_mesh::pmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average,(FLT *) pnts.data(),0,1,2);
		}
	}
	
	/* Do flow communication */
	/* Vertices */
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		vc0load(mp_phase,ug.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= vc0wait_rcv(mp_phase,ug.v.data());
	}
	
	/* Sides */
	sc0load(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));
	smsgpass(boundary::all,0,boundary::symmetric);
	sc0wait_rcv(gbl->res.s.data(),0,sm0-1,ug.s.extent(secondDim));
	
	return;
}

//  Rows with dirichlet boundary conditions will not match.  Sparse Jacobian will have a 1 on diagonal only
//  Diagonal entry of communication rows will not match because of equality constraint
//  For test jacobian, unknowns on each side are separately changed but residual is matched
//  For actual jacobian, variables are individual changed as well and rows are added together so this should agree except
//  diagonal entry is switched between equations as follows
//	Shift all entries for this vertex */
//	for(int n_mpi=0;n_mpi<vdofs;++n_mpi) {
//		FLT dval = (*pJ_mpi)(row+n,row_mpi+n_mpi);
//		(*pJ_mpi)(row+n,row_mpi+n_mpi) = 0.0;
//		x.J(row+n,row+n_mpi) += dval;
//	}
//	Thus, the diagonal entries will not match
//  The remote diagonal value will be 0 in the sparse representation
//  and the sparse local diagonal entry will be the sum of the two entries in the test jacobian (remote & local)
//	For sides, all of the continuous mode entries are moved over for all degrees of freedom so there may be a lot of non-matching modes

void tri_hp::test_jacobian() {
    
    /* This outputs a binary file that can be read by matlab */
    /* add to Matlab path: ${HOME}/Packages/petsc/share/petsc/matlab */
     /* J = PetscBinaryRead('Jacobian') */
    /* [U,S,V] = svds(J,1,'smallest') */
    PetscViewer viewer;
    int ierr1 = PetscViewerBinaryOpen(PETSC_COMM_WORLD, "Jacobian", FILE_MODE_WRITE, &viewer);
    CHKERRABORT(MPI_COMM_WORLD,ierr1);
    
    int ierr2 = MatView(petsc_J,viewer);
    CHKERRABORT(MPI_COMM_WORLD,ierr2);
    PetscViewerDestroy(&viewer);
    
    /* Get neighboring block information */
    const int nblock = sim::blks.nblock;
    const int myid = sim::blks.myid;
    const int vdofs = NV +(mmovement == tri_hp::coupled_deformable)*ND;
    assert(sim::blks.myblock == 1);
    
    /* Get counts from other meshes */
    Array<int,1> sndsizes(4*nblock);
    sndsizes = 0;
    sndsizes(myid*4) = npnt;
    sndsizes(myid*4+1) = nseg;
    sndsizes(myid*4+2) = ntri;
    sndsizes(myid*4+3) = jacobian_size;

    
    Array<int,1> sizes(4*nblock);
    sim::blks.allreduce(sndsizes.data(),sizes.data(),4*nblock,blocks::int_msg,blocks::sum);
    ~sndsizes;
    
	/*************** TESTING ROUTINE ***********************/
	/* HARD TEST OF JACOBIAN WITH DIRICHLET B.C.'s APPLIED */
	/*******************************************************/
	Array<FLT,1> dw(NV);
	dw = eps_a;
	FLT dx = eps_a;
	int dof = jacobian_size;
	PetscScalar *array;
    std::string desc;

	petsc_rsdl();
	VecGetArray(petsc_f,&array);
	Array<FLT,1> rbar(array, shape(jacobian_size), duplicateData);
	VecRestoreArray(petsc_f, &array);

	const PetscInt *ranges;
	VecGetOwnershipRanges(petsc_f,&ranges);
	for(int proc=0;proc < sim::blks.nproc;++proc) {
		if (proc == sim::blks.myid) {
			
			*gbl->log << "ON PROCESSOR COMPONENTS FOR PROCESSOR " << proc << std::endl;
			*gbl->log << "COLUMN NUMBERS ARE LOCAL INDICES" << std::endl;
			IS Irows;
			IS Icols;
			ISCreateStride(MPI_COMM_WORLD,jacobian_size,jacobian_start,1,&Irows);
			ISCreateStride(MPI_COMM_WORLD,ranges[proc+1]-ranges[proc],ranges[proc],1,&Icols);
			Mat *submat;
			MatCreateSubMatrices(petsc_J,1,&Irows,&Icols,MAT_INITIAL_MATRIX,&submat);
			ISDestroy(&Irows);
			ISDestroy(&Icols);
			MatDestroyMatrices(1,&submat);
			
			Array<FLT,2> testJ(dof,dof);
			testJ = 0.0;
			
			int ind = 0;
			if (mmovement != coupled_deformable) {
				for (int pind=0;pind<npnt;++pind) {
					for(int n=0;n<NV;++n) {
						FLT stored_value = ug.v(pind,n);
						ug.v(pind,n) += dw(n);
						
						petsc_rsdl();
						VecGetArray(petsc_f,&array);
						Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
						testJ(Range::all(),ind) = (rtemp-rbar)/(ug.v(pind,n)-stored_value);
						VecRestoreArray(petsc_f, &array);

						
						
						++ind;
						ug.v(pind,n) = stored_value;
					}
				}
			}
			else {
				for (int pind=0;pind<npnt;++pind) {
					for(int n=0;n<NV;++n) {
						FLT stored_value = ug.v(pind,n);
						ug.v(pind,n) += dw(n);
						petsc_rsdl();
						VecGetArray(petsc_f,&array);
						Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
						testJ(Range::all(),ind) = (rtemp-rbar)/(ug.v(pind,n)-stored_value);
						VecRestoreArray(petsc_f, &array);
						++ind;
						ug.v(pind,n) = stored_value;
					}
					
					for(int n=0;n<ND;++n) {
						FLT stored_value = pnts(pind)(n);
						pnts(pind)(n) += dx;
						petsc_rsdl();
						VecGetArray(petsc_f,&array);
						Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
						testJ(Range::all(),ind) = (rtemp-rbar)/(pnts(pind)(n)-stored_value);
						VecRestoreArray(petsc_f, &array);
						++ind;
						pnts(pind)(n) = stored_value;
					}
				}
			}
			
			for (int sind=0;sind<nseg;++sind) {
				for (int m=0;m<sm0;++m) {
					for(int n=0;n<NV;++n) {
						FLT stored_value = ug.s(sind,m,n);
						ug.s(sind,m,n) += dw(n);
						petsc_rsdl();
						VecGetArray(petsc_f,&array);
						Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
						testJ(Range::all(),ind) = (rtemp-rbar)/(ug.s(sind,m,n) -stored_value);
						VecRestoreArray(petsc_f, &array);
						++ind;
						ug.s(sind,m,n) = stored_value;
					}
				}
			}
			
			for (int tind=0;tind<ntri;++tind) {
				for (int m=0;m<im0;++m) {
					for(int n=0;n<NV;++n) {
						FLT stored_value = ug.i(tind,m,n);
						ug.i(tind,m,n) += dw(n);
						petsc_rsdl();
						VecGetArray(petsc_f,&array);
						Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
						testJ(Range::all(),ind) = (rtemp-rbar)/(ug.i(tind,m,n)-stored_value);
						VecRestoreArray(petsc_f, &array);
						++ind;
						ug.i(tind,m,n) = stored_value;
					}
				}
			}
			
			/* How to perturb extra unknowns? */
			for(int i=0;i<nebd;++i) {
				if (!(hp_ebdry(i)->curved) || !(hp_ebdry(i)->coupled)) continue;
				for(int j=0;j<ebdry(i)->nseg;++j) {
					for(int m=0;m<sm0;++m) {
						for(int n=0;n<ND;++n) {
							FLT stored_value = hp_ebdry(i)->crds(j,m,n);
							hp_ebdry(i)->crds(j,m,n) += dx;
							petsc_rsdl();
							VecGetArray(petsc_f,&array);
							Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
							testJ(Range::all(),ind) = (rtemp-rbar)/(hp_ebdry(i)->crds(j,m,n)-stored_value);
							VecRestoreArray(petsc_f, &array);
							++ind;
							hp_ebdry(i)->crds(j,m,n) = stored_value;
						}
					}
				}
			}
			assert(ind == jacobian_size);
			
			const PetscScalar *vals;
			const PetscInt *cols;
			int nnz;
			for(int i=0;i<jacobian_size;++i) {
				MatGetRow(petsc_J,i+jacobian_start,&nnz,&cols,&vals);
                local_index_to_mesh_descriptor(i,desc);
				*gbl->log << desc << ' ';
				int cnt = 0;
				/* Skip things on previous processors */
				while (cols[cnt] < jacobian_start) {
					++cnt;
				}
				
				for(int j=0;j<jacobian_size;++j) {
					if (!(fabs(testJ(i,j)) < DEBUG_ABS_TOL)) {
						if (cnt >= nnz) {
                            local_index_to_mesh_descriptor(j,desc);
							*gbl->log << " (F " <<  desc << ' ' << testJ(i,j) << ") ";
							continue;
						}
						if (cols[cnt] == j+jacobian_start) {
                            if (!(fabs(testJ(i,j) -vals[cnt])/(fabs(testJ(i,j)) +fabs(vals[cnt])) < DEBUG_REL_TOL)) {
                                local_index_to_mesh_descriptor(j,desc);
                                *gbl->log << " (FS " << desc << ", "<< testJ(i,j) << ' ' << vals[cnt] << ") ";
                            }
							++cnt;
						}
						else if (cols[cnt] < j+jacobian_start) {
							do {
                                if (!(fabs(vals[cnt]) < DEBUG_ABS_TOL) && cols[cnt] >= jacobian_start) {
                                    local_index_to_mesh_descriptor(cols[cnt]-jacobian_start,desc);
                                    *gbl->log << " (S " << desc << ' ' << vals[cnt] << ") ";

                                }
								++cnt;
							} while (cols[cnt] < j+jacobian_start);
							--j;
						}
						else {
                            local_index_to_mesh_descriptor(j,desc);
							*gbl->log << " (F " <<  desc  << ' ' << testJ(i,j) << ") ";
						}
					}
				}
				if (cnt < nnz) {
					do {
                        if (!(fabs(vals[cnt]) < DEBUG_ABS_TOL) && cols[cnt] < jacobian_start+jacobian_size) {
                            local_index_to_mesh_descriptor(cols[cnt] -jacobian_start,desc);
                            *gbl->log << " (S " << desc << ' ' << vals[cnt] << ") ";
                        }
					} while (++cnt < nnz);
				}
				
				*gbl->log << std::endl;
				MatRestoreRow(petsc_J,i,&nnz,&cols,&vals);
			}
		}
		else {
            int const npnt_mpi = sizes(4*proc);
            int const nseg_mpi = sizes(4*proc+1);
            int const ntri_mpi = sizes(4*proc+2);
            int const jacobian_size_mpi = sizes(4*proc+3);
            
			IS Irows;
			IS Icols;
			ISCreateStride(MPI_COMM_WORLD,jacobian_size,jacobian_start,1,&Irows);
			ISCreateStride(MPI_COMM_WORLD,ranges[proc+1]-ranges[proc],ranges[proc],1,&Icols);
			Mat *submat;
			MatCreateSubMatrices(petsc_J,1,&Irows,&Icols,MAT_INITIAL_MATRIX,&submat);
			ISDestroy(&Irows);
			ISDestroy(&Icols);
			
			Array<FLT,2> testJ(jacobian_size,ranges[proc+1]-ranges[proc]);
			testJ = 0.0;
			vsi ugtemp;
			ugtemp.v.resize(ug.v.extent(firstDim),ug.v.extent(secondDim));
			ugtemp.s.resize(ug.s.extent(firstDim),ug.s.extent(secondDim),ug.s.extent(thirdDim));
			ugtemp.i.resize(ug.i.extent(firstDim),ug.i.extent(secondDim),ug.i.extent(thirdDim));
			ugtemp.v = ug.v;
			ugtemp.s = ug.s;
			ugtemp.i = ug.i;
			
			
			Array<TinyVector<FLT,ND>,1> pntstemp;
			if (coupled_deformable) {
				pntstemp.resize(pnts.extent(firstDim));
				pntstemp = pnts;
			}
			
			
			/* calculate residual enough times */
			for(int i=0;i<ranges[proc+1]-ranges[proc];++i) {
				// enforce_continuity(ug, pnts);
				petsc_rsdl();
				VecGetArray(petsc_f,&array);
				Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
				testJ(Range::all(),i) = (rtemp-rbar)/eps_a;
				VecRestoreArray(petsc_f, &array);
				ug = ugtemp;
				if (coupled_deformable) pnts = pntstemp;
			}
			
			const PetscScalar *vals;
			const PetscInt *cols;
			int nnz;

			*gbl->log << "OFF PROCESSOR COMPONENTS FROM PROCESSOR " << sim::blks.myid << " TO PROCESSOR " << proc << std::endl;
			*gbl->log << "COLUMN NUMBERS ARE GLOBAL INDICES" << std::endl;

			for(int i=0;i<jacobian_size;++i) {
				MatGetRow(*submat,i,&nnz,&cols,&vals);
                local_index_to_mesh_descriptor(i,desc);
                *gbl->log << desc << ' ';
				
				int cnt = 0;
				for(int j=0;j<ranges[proc+1]-ranges[proc];++j){
                    ostringstream nstr;
                    int col = j;
                    if (col < npnt_mpi*vdofs) {
                        nstr << "b" << proc << ',' << 'v' << col/vdofs << ',' << col % vdofs << ',';
                    }
                    else if ((col -= npnt_mpi*vdofs) < nseg_mpi*sm0*NV) {
                        nstr << "b" << proc << ',' << 's' << col/(NV*sm0) << ',' << (col % (sm0*NV))/NV << ',' << col % NV << ',';
                    }
                    else if ((col -= nseg_mpi*NV*sm0) < ntri_mpi*im0*NV) {
                        nstr << "b" << proc << ',' << 'i' << col/(NV*im0) << ',' << (col % (im0*NV))/NV << ',' << col % NV << ',';
                    }
                    else if ((col -= ntri_mpi*NV*im0) < jacobian_size_mpi){
                        nstr << "b" << proc << ',' << 'e' << col/(ND*sm0) << ',' << (col % (sm0*ND))/ND << ',' << col % ND << ',';
                    }
                    desc = nstr.str();
                    
                    if (fabs(testJ(i,j)) > DEBUG_ABS_TOL) {
						if (cnt >= nnz) {
							*gbl->log << " (F " <<  desc << ' ' << testJ(i,j) << ") ";
							continue;
						}
						if (cols[cnt] == j) {
							if (fabs(testJ(i,j) -vals[cnt])/(fabs(testJ(i,j)) +fabs(vals[cnt])) > DEBUG_REL_TOL)
								*gbl->log << " (FS " << desc << ' ' << testJ(i,j) << ' ' << vals[cnt] << ") ";
							++cnt;
						}
						else if (cols[cnt] < j) {
							do {
                                if (fabs(vals[cnt]) > DEBUG_ABS_TOL && cols[cnt] >= ranges[proc]) {
                                    int col = cols[cnt];
                                    if (col < npnt_mpi*vdofs) {
                                        nstr << "b" << proc << ',' << 'v' << col/vdofs << ',' << col % vdofs << ',';
                                    }
                                    else if ((col -= npnt_mpi*vdofs) < nseg_mpi*sm0*NV) {
                                        nstr << "b" << proc << ',' << 's' << col/(NV*sm0) << ',' << (col % (sm0*NV))/NV << ',' << col % NV << ',';
                                    }
                                    else if ((col -= nseg_mpi*NV*sm0) < ntri_mpi*im0*NV) {
                                        nstr << "b" << proc << ',' << 'i' << col/(NV*im0) << ',' << (col % (im0*NV))/NV << ',' << col % NV << ',';
                                    }
                                    else if ((col -= ntri_mpi*NV*im0) < jacobian_size_mpi){
                                        nstr << "b" << proc << ',' << 'e' << col/(ND*sm0) << ',' << (col % (sm0*ND))/ND << ',' << col % ND << ',';
                                    }
                                    *gbl->log << " (S " << nstr.str() << ' ' << vals[cnt] << ") ";

                                }
								++cnt;
							} while (cols[cnt] < j);
							--j;
						}
						else {
							*gbl->log << " (F " <<  desc  << ' ' << testJ(i,j) << ") ";
						}
					}
				}
				if (cnt < nnz) {
					do {
						if (fabs(vals[cnt]) > DEBUG_ABS_TOL && cols[cnt] < ranges[proc+1])
							*gbl->log << " (S " << cols[cnt] +ranges[proc] << ' ' << vals[cnt] << ") ";
					} while (++cnt < nnz);
				}
				
				*gbl->log << std::endl;
				MatRestoreRow(*submat,i,&nnz,&cols,&vals);
			}
			
			
			
			MatDestroyMatrices(1,&submat);
			
		}
	}
}

void tri_hp::local_index_to_mesh_descriptor(int col, std::string& desc) {
    const int block = gbl->idnum;
    const int vdofs = NV +(mmovement == tri_hp::coupled_deformable)*ND;

    ostringstream nstr;
    
    nstr << "b" << block << ',';
    if (col < npnt*vdofs) {
        nstr << 'v' << col/vdofs << ',' << col % vdofs;
    }
    else if ((col -= npnt*vdofs) < nseg*sm0*NV) {
        nstr << 's' << col/(NV*sm0) << ',' << (col % (sm0*NV))/NV << ',' << col % NV;
    }
    else if ((col -= nseg*NV*sm0) < ntri*im0*NV) {
        nstr << 'i' << col/(NV*im0) << ',' << (col % (im0*NV))/NV << ',' << col % NV;
    }
    else {
        col -= ntri*NV*im0;
        nstr << 'e' << col/(ND*sm0) << ',' << (col % (sm0*ND))/ND << ',' << col % ND;
    }
    desc = nstr.str();
}

#endif
