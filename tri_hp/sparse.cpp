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
		if (mmovement == coupled_deformable) bw(Range(begin_iso,size_sparse_matrix-1)) = vdofs*(sm+3) +NV*(im+2*sm);
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
	for(int i=0;i<nebd;++i) {
		if (!hp_ebdry(i)->curved || !hp_ebdry(i)->coupled) continue;
		
		for(int j=0;j<ebdry(i)->nseg;++j) {

			int tind = seg(ebdry(i)->seg(j)).tri(0);
			bw(Range(begin_tri +tind*im*NV,begin_tri +(tind+1)*im*NV-1)) += ND*sm;
			
			for(int j=0;j<3;++j) {
				int sind = tri(tind).seg(j);
				bw(Range(begin_seg+sind*NV*sm,begin_seg+(sind+1)*NV*sm-1)) += ND*sm;
				int pind = tri(tind).pnt(j);
				bw(Range(pind*vdofs,(pind+1)*vdofs-1))+= ND*sm;
			}
		}
	}
	
	MatCreateSeqAIJ(PETSC_COMM_SELF,size_sparse_matrix,size_sparse_matrix,PETSC_NULL,bw.data(),&petsc_J);
	//MatCreateMPIAIJ(comm,size_sparse_matrix,size_sparse_matrix,global_size,global_size,int d nz,int *d nnz, int o nz,int *o nnz,Mat *A);

	return;
}


void tri_hp::petsc_jacobian() {
	int gindx,eind,find;
	int sm = basis::tri(log2p)->sm();
	int tm = basis::tri(log2p)->tm();

	int vdofs = NV;
	if (mmovement == coupled_deformable) vdofs += ND;
	Array<FLT,2> K(NV*(sm+2),NV*(sm+2));
	Array<FLT,2> Kiso(vdofs*(sm+2),vdofs*(sm+2));
	Array<int,1> loc_to_glo(NV*(sm+2));
	Array<int,1> loc_to_glo_iso(vdofs*(sm+2));

	int kcol,krow;
	FLT dw = 1.0e-4;// fix me temp make a global value?
	FLT sgn,msgn;
	Array<FLT,2> Rbar(NV,tm);
	
	/* DO NEUMANN & COUPLED BOUNDARY CONDITIONS */
	for(int i=0;i<nebd;++i){
		for(int j=0;j<ebdry(i)->nseg;++j){
			
			/* Calculate and store initial residual */
			eind = ebdry(i)->seg(j);
			ugtouht1d(eind);
			lf = 0.0;
			hp_ebdry(i)->element_rsdl(j,0);

			int ind = 0;
			for(int k=0;k<2;++k) {
				int gindx = NV*seg(eind).pnt(k);
				for(int n=0;n<NV;++n) {
					Rbar(n,k) = lf(n)(k);
					loc_to_glo(ind++) = gindx+n;
				}
			}
			
			/* EDGE MODES */
			if (sm > 0) {
				int gbl_eind = npnt*NV + eind*sm*NV;
				for (int m = 0; m < sm; ++m) {
					for(int n = 0; n < NV; ++n) {
						Rbar(n,m+2) = lf(n)(m+2);
						loc_to_glo(ind++) = gbl_eind + m*NV + n;
					}
				}
			}
			
			/* Numerically create Jacobian */
			kcol = 0;
			for(int mode = 0; mode < sm+2; ++mode){
				for(int var = 0; var < NV; ++var){
					uht(var)(mode) += dw;
					
					lf = 0.0;
					hp_ebdry(i)->element_rsdl(j,0);

					krow = 0;
					for(int k=0;k<sm+2;++k)
						for(int n=0;n<NV;++n)
							K(krow++,kcol) = (lf(n)(k)-Rbar(n,k))/dw;
					
					++kcol;
					
					uht(var)(mode) -= dw;
				}
			}
			

			MatSetValues(petsc_J,NV*(sm+2),loc_to_glo.data(),NV*(sm+2),loc_to_glo.data(),K.data(),ADD_VALUES);
		}
	}
		
	int kn = tm*NV;
	K.resize(kn,kn);
	loc_to_glo.resize(kn);
		
	/* NOW DO ELEMENTS */
	for(int tind = 0; tind < ntri; ++tind){	
		
		int ind = 0;

		element_jacobian(tind, K);
		
		for (int m = 0; m < 3; ++m) {
			gindx = NV*tri(tind).pnt(m);
			for (int n = 0; n < NV; ++n)
				loc_to_glo(ind++) = gindx+n;
		}		
		
		/* EDGE MODES */
		if (basis::tri(log2p)->p() > 1) {
			for(int i = 0; i < 3; ++i) {
				eind = npnt*NV + tri(tind).seg(i)*basis::tri(log2p)->sm()*NV;
				sgn = tri(tind).sgn(i);
				msgn = 1;
				for (int m = 0; m < basis::tri(log2p)->sm(); ++m) {
					for(int n = 0; n < NV; ++n) {
						for(int j = 0; j < kn; ++j) {
							K(ind,j) *= msgn;
							K(j,ind) *= msgn;
						}
						loc_to_glo(ind++) = eind + m*NV + n;
					}
					msgn *= sgn;
				}
			}
		}
		
		/* INTERIOR	MODES */
		if (basis::tri(log2p)->p() > 2) {
			find = npnt*NV+nseg*basis::tri(log2p)->sm()*NV+tind*basis::tri(log2p)->im()*NV;
			for(int m = 0; m < basis::tri(log2p)->im(); ++m) {
				for(int n = 0; n < NV; ++n){
					loc_to_glo(ind++) = find+m*NV+n;
				}
			}
					
		}		
		
		
		MatSetValues(petsc_J,kn,loc_to_glo.data(),kn,loc_to_glo.data(),K.data(),ADD_VALUES);
			
				
	}	

	return;
}

void tri_hp::petsc_rsdl() {
	rsdl();

	PetscScalar *array;
	VecGetArray(petsc_f,&array);
	Array<FLT,1> res(array, shape(size_sparse_matrix), neverDeleteData);
	petsc_make_1D_rsdl_vector(res);
	VecRestoreArray(petsc_f, &array);
	
	return;
}
#endif
