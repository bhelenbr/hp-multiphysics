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
	
	int begin_seg = npnt*NV;
	int begin_tri = begin_seg+nseg*sm*NV;
	
	/* SELF CONNECTIONS */
	bw(Range(0,begin_seg-1)) = NV;  
	if (sm) bw(Range(begin_seg,begin_tri-1)) = (2 +sm)*NV;
	if (im) bw(Range(begin_tri,size_sparse_matrix-1)) = tm*NV;
	
	/* connected edges and vertices to vertex */
	for(int i=0; i<npnt; ++i)	
		for(int n=0;n<NV;++n)
			bw(i*NV+n) += NV*(pnt(i).nnbor*(sm+1));		
	
	for(int i=0; i<ntri; ++i) {	
		/* add interior modes to each vertex */
		for(int j=0;j<3;++j)
			for(int n=0;n<NV;++n)
				bw(tri(i).pnt(j)*NV+n) += (im+sm)*NV;
				
		
		/* add interior modes to each edge */
		for(int j=0;j<3;++j)
			for(int m=0;m<sm;++m)
				for(int n=0;n<NV;++n)
					bw(begin_seg+tri(i).seg(j)*sm*NV+m*NV+n) += (im +2*sm+1)*NV;
	}
	

	MatCreateSeqAIJ(PETSC_COMM_SELF,size_sparse_matrix,size_sparse_matrix,PETSC_NULL,bw.data(),&petsc_J);
	//MatCreateMPIAIJ(comm,size_sparse_matrix,size_sparse_matrix,global_size,global_size,int d nz,int *d nnz, int o nz,int *o nnz,Mat *A);


	
	
	return;
}


/* clears row in sparse matrix, inserts 1.0 on diagonal, and inserts 0.0 in the residual */
void tri_hp::sparse_dirichlet(int ind){
	
	const PetscInt row = ind;
	PetscScalar zero = 0.0;
	PetscScalar one = 1.0;

	/* apply dirichlet by inserting zero in f */
	VecSetValues(petsc_f,1,&row,&zero,INSERT_VALUES);
	dirichlet_rows(row_counter++)=ind;

	
	return;

// too slow, but easy to implement
//	MatZeroRows(petsc_J,1,&row,1.0);
//
//	return;

}


void tri_hp::create_jacobian() {
	int gindx,eind,find;
	int sm = basis::tri(log2p)->sm();
	int tm = basis::tri(log2p)->tm();

	Array<FLT,2> K(NV*(sm+2),NV*(sm+2));
	Array<int,1> loc_to_glo(NV*(sm+2));

	int kcol,krow;
	FLT dw = 1.0e-4;// fix me temp make a global value?
	FLT sgn,msgn;
	Array<FLT,2> Rbar(NV,tm);
	
	/* DO NEUMANN BOUNDARY CONDITIONS */
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

		create_local_jacobian_matrix(tind, K);
		
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

void tri_hp::create_local_jacobian_matrix(int tind, Array<FLT,2> &K) {
	Array<TinyVector<FLT,MXTM>,1> R(NV),Rbar(NV),lf_re(NV),lf_im(NV);
	FLT dw = 1.0e-4;  //dw=sqrt(eps/l2_norm(q))
	
	
	ugtouht(tind);
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tri(log2p)->tm();++i)
		for(int n=0;n<NV;++n)
			Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);
	
	int kcol = 0;
	for(int mode = 0; mode < basis::tri(log2p)->tm(); ++mode){
		for(int var = 0; var < NV; ++var){
			uht(var)(mode) += dw;
			
			element_rsdl(tind,0,uht,lf_re,lf_im);

			int krow = 0;
			for(int i=0;i<basis::tri(log2p)->tm();++i)
				for(int n=0;n<NV;++n)
					K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw;
			
			++kcol;
			uht(var)(mode) -= dw;
		}
	}	

	return;
}

void tri_hp::create_rsdl() {
	int gindx,eind,find,iind,sgn,msgn,mode;
	
	int kn = NV*basis::tri(log2p)->tm();
	Array<FLT,1> lclres(kn);
	Array<int,1> loc_to_glo(kn);
	
	
	/* DO BOUNDARY CONDITION RESIDUALS */
	int sm = basis::tri(log2p)->sm();
	for(int i=0;i<nebd;++i) {
		for(int j=0;j<ebdry(i)->nseg;++j) {
			
			eind = ebdry(i)->seg(j);
			ugtouht1d(eind);
			lf = 0.0;
			hp_ebdry(i)->element_rsdl(j,0);

			/* Vertex Modes */
			int ind = 0;
			for (int m = 0; m < 2; ++m) {
				int gindx = NV*seg(eind).pnt(m);
				for (int n = 0; n < NV; ++n) {
					lclres(ind) = lf(n)(m);
					loc_to_glo(ind++) = gindx+n;
				}
			}		
			
			/* Edge Modes */
			if (sm > 0) {
				int gbl_eind = npnt*NV + eind*sm*NV;
				for (int m = 0; m < sm; ++m) {
					for(int n = 0; n < NV; ++n) {
						lclres(ind) = lf(n)(m+2);
						loc_to_glo(ind++) = gbl_eind + m*NV + n;
					}
				}
			}
			VecSetValues(petsc_f,(sm+2)*NV,loc_to_glo.data(),lclres.data(),ADD_VALUES);
		}
	}
	
	
	/* NOW DO ELEMENT RESIDUAL */
	for(int tind = 0; tind < ntri; ++tind){	
		
		int ind = 0;
		
		create_local_rsdl(tind, lclres);
		
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
						lclres(ind) *= msgn;
						loc_to_glo(ind++) = eind + m*NV + n;
					}
					msgn *= sgn;
				}
			}
		}
		
		/* FACE MODES */
		if (basis::tri(log2p)->p() > 2) {
			find = npnt*NV+nseg*basis::tri(log2p)->sm()*NV+tind*basis::tri(log2p)->im()*NV;
			for(int m = 0; m < basis::tri(log2p)->im(); ++m) {
				for(int n = 0; n < NV; ++n){
					loc_to_glo(ind++) = find+m*NV+n;
				}					
			}					
		}		
		VecSetValues(petsc_f,kn,loc_to_glo.data(),lclres.data(),ADD_VALUES);	
	}
	
	return;
}

void tri_hp::create_local_rsdl(int tind, Array<FLT,1> &lclres) {
	Array<TinyVector<FLT,MXTM>,1> lf_re(NV),lf_im(NV);

	int ind = 0;

	ugtouht(tind);
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tri(log2p)->tm();++i)
		for(int n=0;n<NV;++n)
			lclres(ind++)=lf_re(n)(i)+lf_im(n)(i);	
	
	return;
}
#endif
