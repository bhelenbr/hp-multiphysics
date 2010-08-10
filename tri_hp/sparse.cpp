/*
 *  sparse.cpp
 *  tri_hp
 *
 *  Created by michael brazell on 9/14/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "ins/tri_hp_ins.h"
#include "hp_boundary.h"
#include <r_tri_boundary.h>

#ifdef petsc

void tri_hp::petsc_jacobian() {
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
		//loc_to_glo += jacobian_start;  /* shift indices for this block */
		
#ifdef MY_SPARSE
		J.add_values(kn,loc_to_glo,kn,loc_to_glo,K);
#else
		MatSetValuesLocal(petsc_J,kn,loc_to_glo.data(),kn,loc_to_glo.data(),K.data(),ADD_VALUES);
#endif
	}
	/* APPLY ALL MOVING MESH B.C.'s FIRST */
	if (mmovement == coupled_deformable) {
		/* This mostly does nothing, except for angled boundarys */
		for(int i=0;i<nebd;++i) 
			r_sbdry(i)->jacobian();  // Fixme: parallel jacobian stuff be done here

#ifndef MY_SPARSE	
		/* PETSC IS RETARDED */
		FLT zero = 0.0;
		for (int i=jacobian_start +npnt*(NV+ND)+sm*NV*nseg+ntri*im*NV;i<jacobian_start +jacobian_size;++i) 
			MatSetValuesLocal(petsc_J,1,&i,1,&i,&zero,ADD_VALUES);
			
		MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);	
#endif

		for(int i=0;i<nebd;++i)
			r_sbdry(i)->jacobian_dirichlet();
	}
					
	/* DO NEUMANN & COUPLED BOUNDARY CONDITIONS */
	for(int i=0;i<nebd;++i) 
		hp_ebdry(i)->petsc_jacobian();
	
	/* This one does nothing */
	for(int i=0;i<nvbd;++i) 
		hp_vbdry(i)->petsc_jacobian();
		
#ifndef MY_SPARSE
	MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);	
#endif
		
	for(int i=0;i<nebd;++i) 
		hp_ebdry(i)->petsc_jacobian_dirichlet();

	for(int i=0;i<nvbd;++i) 
		hp_vbdry(i)->petsc_jacobian_dirichlet();
	
	return;
}

void tri_hp::petsc_rsdl() {


	rsdl();
	
	/* Do communication */
	/* Vertices */
	int last_phase, mp_phase;
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		vc0load(mp_phase,gbl->res.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= vc0wait_rcv(mp_phase,gbl->res.v.data());
	}
	
	/* Sides */
	sc0load(gbl->res.s.data(),0,sm0-1,gbl->res.s.extent(secondDim));
	smsgpass(boundary::all,0,boundary::symmetric);
	sc0wait_rcv(gbl->res.s.data(),0,sm0-1,gbl->res.s.extent(secondDim));  
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet();
		
	for(int i=0;i<nvbd;++i) 
		hp_vbdry(i)->vdirichlet();

	for(int i=0;i<nvbd;++i) 
		hp_vbdry(i)->vdirichlet2d();

		
	/* APPLY DIRCHLET B.C.S TO MODE */
	for(int m=0;m<basis::tri(log2p)->sm();++m)
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->sdirichlet(m);
				
	PetscScalar *array;
	VecGetArray(petsc_f,&array);
	Array<FLT,1> res(array, shape(jacobian_size), neverDeleteData);
	petsc_make_1D_rsdl_vector(res);
	VecRestoreArray(petsc_f, &array);
		
	return;
}


#define DEBUG_TOL 1.0e-7
#define WBC

void tri_hp::test_jacobian() {

	/*************** TESTING ROUTINE ***********************/
	/* HARD TEST OF JACOBIAN WITH DIRICHLET B.C.'s APPLIED */
	/*******************************************************/
	const FLT eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
	Array<FLT,1> dw(NV);
	dw = eps_a;
	FLT dx = eps_a;
	int dof = jacobian_size;
	Array<FLT,2> testJ(dof,dof);
	testJ = 0.0;
	PetscScalar *array;

#ifdef WBC
	petsc_rsdl();
	VecGetArray(petsc_f,&array);
	Array<FLT,1> rbar(array, shape(jacobian_size), duplicateData);
	VecRestoreArray(petsc_f, &array);
#else
	rsdl();
	petsc_make_1D_rsdl_vector(rbar);
#endif

	int ind = 0;
	if (mmovement != coupled_deformable) {
		for (int pind=0;pind<npnt;++pind) {
			for(int n=0;n<NV;++n) {
				FLT stored_value = ug.v(pind,n);
				ug.v(pind,n) += dw(n);
#ifdef WBC
				petsc_rsdl();
				VecGetArray(petsc_f,&array);
				Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
				testJ(Range::all(),ind) = (rtemp-rbar)/dw(n);
				VecRestoreArray(petsc_f, &array);
#else
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-rbar)/dw(n);
#endif
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
#ifdef WBC
				petsc_rsdl();
				VecGetArray(petsc_f,&array);
				Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
				testJ(Range::all(),ind) = (rtemp-rbar)/dw(n);
				VecRestoreArray(petsc_f, &array);
#else
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-rbar)/dw(n);
#endif
				++ind;
				ug.v(pind,n) = stored_value;
			}
			
			for(int n=0;n<ND;++n) {
				FLT stored_value = pnts(pind)(n);
				pnts(pind)(n) += dx;
#ifdef WBC
				petsc_rsdl();
				VecGetArray(petsc_f,&array);
				Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
				testJ(Range::all(),ind) = (rtemp-rbar)/dw(n);
				VecRestoreArray(petsc_f, &array);
#else
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-rbar)/dw(n);
#endif
				++ind;
				pnts(pind)(n) = stored_value;
			}
		}
	}
			
	for (int sind=0;sind<nseg;++sind) {
		for (int m=0;m<sm0;++m) {
			for(int n=0;n<NV;++n) {
				ug.s(sind,m,n) += dw(n);
#ifdef WBC
				petsc_rsdl();
				VecGetArray(petsc_f,&array);
				Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
				testJ(Range::all(),ind) = (rtemp-rbar)/dw(n);
				VecRestoreArray(petsc_f, &array);
#else
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-rbar)/dw(n);
#endif
				++ind;
				ug.s(sind,m,n) -= dw(n);
			}
		}			
	}

	for (int tind=0;tind<ntri;++tind) {
		for (int m=0;m<im0;++m) {
			for(int n=0;n<NV;++n) {
				ug.i(tind,m,n) += dw(n);
#ifdef WBC
				petsc_rsdl();
				VecGetArray(petsc_f,&array);
				Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
				testJ(Range::all(),ind) = (rtemp-rbar)/dw(n);
				VecRestoreArray(petsc_f, &array);
#else
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-rbar)/dw(n);
#endif
				++ind;
				ug.i(tind,m,n) -= dw(n);
			}
		}
	}
	
	/* How to perturb extra unknowns? */
	for(int i=0;i<nebd;++i) {
		if (!(hp_ebdry(i)->curved) || !(hp_ebdry(i)->coupled)) continue;
		for(int j=0;j<ebdry(i)->nseg;++j) {
			for(int m=0;m<sm0;++m) {
				for(int n=0;n<ND;++n) {
					hp_ebdry(i)->crds(j,m,n) += dx;
#ifdef WBC
					petsc_rsdl();
					VecGetArray(petsc_f,&array);
					Array<FLT,1> rtemp(array, shape(jacobian_size), neverDeleteData);
					testJ(Range::all(),ind) = (rtemp-rbar)/dw(n);
					VecRestoreArray(petsc_f, &array);
#else
					rsdl();
					petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
					testJ(Range::all(),ind) = (testJ(Range::all(),ind)-rbar)/dw(n);
#endif
					++ind;
					hp_ebdry(i)->crds(j,m,n) -= dx;
				}
			}
		}
	}
		
	const PetscScalar *vals;
	const PetscInt *cols;
	int nnz;
	
	for(int i=0;i<ind;++i) {
		MatGetRow(petsc_J,i+jacobian_start,&nnz,&cols,&vals);
		*gbl->log << "row " << i+jacobian_start << ": ";
		int cnt = 0;
		for(int j=0;j<ind;++j) {
			if (fabs(testJ(i,j)) > DEBUG_TOL) {
				if (cnt >= nnz) {
					*gbl->log << " (Extra entry in full matrix " <<  j+jacobian_start << ' ' << testJ(i,j) << ") ";
					continue;
				}
				if (cols[cnt] == j+jacobian_start) {
					if (fabs(testJ(i,j) -vals[cnt]) > DEBUG_TOL) 
						*gbl->log << " (Jacobian " << j+jacobian_start << ", "<< testJ(i,j) << ' ' << vals[cnt] << ") ";
					++cnt;
				}
				else if (cols[cnt] < j+jacobian_start) {
					do {
						if (fabs(vals[cnt]) > DEBUG_TOL)
							*gbl->log << " (Extra entry in sparse matrix " << cols[cnt] << ' ' << vals[cnt] << ") ";
						++cnt;
					} while (cols[cnt] < j+jacobian_start);
					--j;
				}
				else {
					*gbl->log << " (Extra entry in full matrix " <<  j+jacobian_start  << ' ' << testJ(i,j) << ") ";
				}
			}
		}
		if (cnt < nnz) {
			do {
				if (fabs(vals[cnt]) > DEBUG_TOL)
					*gbl->log << " (Extra entry in sparse matrix " << cols[cnt] << ' ' << vals[cnt] << ") ";
			} while (++cnt < nnz);
		}
			
		*gbl->log << std::endl;
		MatRestoreRow(petsc_J,i,&nnz,&cols,&vals);
	}
	
	// MatView(petsc_J,0);
}

#endif
