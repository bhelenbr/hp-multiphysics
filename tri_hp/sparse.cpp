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


#ifdef petsc

// #define DEBUG
#define DEBUG_TOL 1.0e-9

/* finds number of nonzeros per row */ 
void tri_hp::find_sparse_bandwidth(){
	
	Array<int,1> bw(size_sparse_matrix);
	
	int sm=basis::tri(log2p)->sm();
	int im=basis::tri(log2p)->im();
	int tm=basis::tri(log2p)->tm();

	int vdofs;
	if (mmovement != coupled_deformable) 
		vdofs = NV;
	else
		vdofs = ND+NV;
		

	int begin_seg = npnt*vdofs;
	int begin_tri = begin_seg+nseg*sm*NV;
	int begin_iso = begin_tri +ntri*im*NV;
	
	/* SELF CONNECTIONS */
	bw(Range(0,begin_seg-1)) = vdofs; 
	if (sm) {
		bw(Range(begin_seg,begin_tri-1)) = (2 +sm)*NV;
		/* connections of high order isoparametric mappings */
		if (mmovement == coupled_deformable && size_sparse_matrix > begin_iso) bw(Range(begin_iso,size_sparse_matrix-1)) = vdofs*(sm+3) +NV*(im+2*sm);
	}
	if (im) bw(Range(begin_tri,size_sparse_matrix-1)) = tm*NV;
	
	/* edges and vertices connected to avertex */
	for(int i=0; i<npnt; ++i)	
		for(int n=0;n<vdofs;++n)
			bw(i*vdofs+n) += pnt(i).nnbor*(vdofs +sm*NV);
	
	for(int i=0; i<ntri; ++i) {	
		/* interior and opposing side mode for each vertex */
		for(int j=0;j<3;++j)
			for(int n=0;n<vdofs;++n)
				bw(tri(i).pnt(j)*vdofs+n) += (im+sm)*NV;
				
		
		/* interior modes,opposing side modes, and opposing vertex to each side */
		for(int j=0;j<3;++j)
			for(int m=0;m<sm;++m)
				for(int n=0;n<NV;++n)
					bw(begin_seg+tri(i).seg(j)*sm*NV+m*NV+n) += (im +2*sm)*NV +1*vdofs;
	}
	
	/* Connections to isoparametric coordinates */
	if (sm) {
		for(int i=0;i<nebd;++i) {
			if (!hp_ebdry(i)->curved || !hp_ebdry(i)->coupled) continue;
			
			for(int j=0;j<ebdry(i)->nseg;++j) {

				int tind = seg(ebdry(i)->seg(j)).tri(0);
				if (im) bw(Range(begin_tri +tind*im*NV,begin_tri +(tind+1)*im*NV-1)) += ND*sm;
				
				for(int j=0;j<3;++j) {
					int sind = tri(tind).seg(j);
					bw(Range(begin_seg+sind*NV*sm,begin_seg+(sind+1)*NV*sm-1)) += ND*sm;
					int pind = tri(tind).pnt(j);
					bw(Range(pind*vdofs,(pind+1)*vdofs-1))+= ND*sm;
				}
			}
		}
	}
	
	MatCreateSeqAIJ(PETSC_COMM_SELF,size_sparse_matrix,size_sparse_matrix,PETSC_NULL,bw.data(),&petsc_J);
	//MatCreateMPIAIJ(comm,size_sparse_matrix,size_sparse_matrix,global_size,global_size,int d nz,int *d nnz, int o nz,int *o nnz,Mat *A);

	return;
}

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
		MatSetValues(petsc_J,kn,loc_to_glo.data(),kn,loc_to_glo.data(),K.data(),ADD_VALUES);
	}
				
	/* DO NEUMANN & COUPLED BOUNDARY CONDITIONS */
	for(int i=0;i<nebd;++i) 
		hp_ebdry(i)->petsc_jacobian();
	
	/* This one does nothing */
	for(int i=0;i<nvbd;++i) 
		hp_vbdry(i)->petsc_jacobian();
		
	MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);	
	
	
	
#ifdef DEBUG
	/*************** TESTING ROUTINE ***********/
	/* HARD TEST OF JACOBIAN ROUTINE */
	FLT dw = 1.0e-4;
	int dof = size_sparse_matrix;
	Array<FLT,2> testJ(dof,dof);
	testJ = 0.0;
	
	rsdl();
	petsc_make_1D_rsdl_vector(testJ(Range::all(),dof-1));
		
	int ind = 0;
	if (mmovement != coupled_deformable) {
		for (int pind=0;pind<npnt;++pind) {
			for(int n=0;n<NV;++n) {
				FLT stored_value = ug.v(pind,n);
				ug.v(pind,n) += dw;
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-testJ(Range::all(),dof-1))/dw;
				++ind;
				ug.v(pind,n) = stored_value;
			}
		}
	}
	else {
		for (int pind=0;pind<npnt;++pind) {
			for(int n=0;n<NV;++n) {
				FLT stored_value = ug.v(pind,n);
				ug.v(pind,n) += dw;
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-testJ(Range::all(),dof-1))/dw;
				++ind;
				ug.v(pind,n) = stored_value;
			}
			
			for(int n=0;n<ND;++n) {
				FLT stored_value = pnts(pind)(n);
				pnts(pind)(n) += dw;
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-testJ(Range::all(),dof-1))/dw;
				++ind;
				pnts(pind)(n) = stored_value;
			}
		}
	}
			
	for (int sind=0;sind<nseg;++sind) {
		for (int m=0;m<sm0;++m) {
			for(int n=0;n<NV;++n) {
				ug.s(sind,m,n) += dw;
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-testJ(Range::all(),dof-1))/dw;
				++ind;
				ug.s(sind,m,n) -= dw;
			}
		}			
	}

	for (int tind=0;tind<ntri;++tind) {
		for (int m=0;m<im0;++m) {
			for(int n=0;n<NV;++n) {
				ug.i(tind,m,n) += dw;
				rsdl();
				petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
				testJ(Range::all(),ind) = (testJ(Range::all(),ind)-testJ(Range::all(),dof-1))/dw;
				++ind;
				ug.i(tind,m,n) -= dw;
			}
		}
	}
	
	/* How to perturb extra unknowns? */
	for(int i=0;i<nebd;++i) {
		if (!(hp_ebdry(i)->curved) || !(hp_ebdry(i)->coupled)) continue;
		for(int j=0;j<ebdry(i)->nseg;++j) {
			for(int m=0;m<sm0;++m) {
				for(int n=0;n<ND;++n) {
					hp_ebdry(i)->crds(j,m,n) += dw;
					rsdl();
					petsc_make_1D_rsdl_vector(testJ(Range::all(),ind));
					testJ(Range::all(),ind) = (testJ(Range::all(),ind)-testJ(Range::all(),dof-1))/dw;
					++ind;
					hp_ebdry(i)->crds(j,m,n) -= dw;
				}
			}
		}
	}
		
	const PetscScalar *vals;
	const PetscInt *cols;
	int nnz;
	
	for(int i=0;i<ind;++i) {
		MatGetRow(petsc_J,i,&nnz,&cols,&vals);
		*gbl->log << "row " << i << ": ";
		int cnt = 0;
		for(int j=0;j<ind;++j) {
			if (fabs(testJ(i,j)) > DEBUG_TOL) {
				if (cnt >= nnz) {
					*gbl->log << " (Extra entry in full matrix " <<  j << ' ' << testJ(i,j) << ") ";
					continue;
				}
				if (cols[cnt] == j) {
					if (fabs(testJ(i,j) -vals[cnt]) > DEBUG_TOL) 
						*gbl->log << " (Jacobian " << j << ", "<< testJ(i,j) << ' ' << vals[cnt] << ") ";
					++cnt;
				}
				else if (cols[cnt] < j) {
					do {
						if (fabs(vals[cnt]) > DEBUG_TOL)
							*gbl->log << " (Extra entry in sparse matrix " << cols[cnt] << ' ' << vals[cnt] << ") ";
						++cnt;
					} while (cols[cnt] < j);
					--j;
				}
				else {
					*gbl->log << " (Extra entry in full matrix " <<  j  << ' ' << testJ(i,j) << ") ";
				}
			}
		}
		*gbl->log << std::endl;
		MatRestoreRow(petsc_J,i,&nnz,&cols,&vals);
	}
	
	// MatView(petsc_J,0);
	
	exit(1);
#endif
		
	for(int i=0;i<nebd;++i) 
		hp_ebdry(i)->petsc_jacobian_dirichlet();
		
	for(int i=0;i<nvbd;++i) 
		hp_vbdry(i)->petsc_jacobian_dirichlet();	
	
	return;
}

void tri_hp::petsc_rsdl() {
	rsdl();
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet();

	for(int i=0;i<nvbd;++i) {
		hp_vbdry(i)->vdirichlet2d();
		
	/* APPLY DIRCHLET B.C.S TO MODE */
	for(int m=0;m<basis::tri(log2p)->sm();++m)
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->sdirichlet(m);
				
	PetscScalar *array;
	VecGetArray(petsc_f,&array);
	Array<FLT,1> res(array, shape(size_sparse_matrix), neverDeleteData);
	petsc_make_1D_rsdl_vector(res);
	VecRestoreArray(petsc_f, &array);
	
	
	return;
}
#endif
