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
	
	bw = 0;
	
	/* opposite edge and vertices from vertex */
	for(int i=0; i<npnt; ++i)	
		for(int n=0;n<NV;++n)
			bw(i*NV+n) += NV*(pnt(i).nnbor*(sm+2));		
	
	
	// temp fix me will this work for parallel BC's
	// may not matter because matrix storage can correct size of bandwidth
	for(int i=0; i<nseg; ++i){		
		/* add side modes to each vertex */
		for(int j=0;j<2;++j)
			for(int n=0;n<NV;++n)
				bw(seg(i).pnt(j)*NV+n) += sm*NV;
		
		/* add side and vertex modes to edges */
		if(seg(i).tri(1) > 0){
			for(int m=0;m<sm;++m)
				for(int n=0;n<NV;++n)
					bw(begin_seg+i*sm*NV+m*NV+n) += NV*(5*sm+4);	
		}
		else{
			for(int m=0;m<sm;++m)
				for(int n=0;n<NV;++n)
					bw(begin_seg+i*sm*NV+m*NV+n) += NV*(3*sm+3);
		}
	}
	
	
	for(int i=0; i<ntri; ++i){	
		/* add interior modes to each vertex */
		for(int j=0;j<3;++j)
			for(int n=0;n<NV;++n)
				bw(tri(i).pnt(j)*NV+n) += im*NV;
		/* add interior modes to each edge */
		for(int j=0;j<3;++j)
			for(int m=0;m<sm;++m)
				for(int n=0;n<NV;++n)
					bw(begin_seg+tri(i).seg(j)*sm*NV+m*NV+n) += im*NV;
		/* add total modes to each interior mode */
		for(int m=0;m<im;++m)
			for(int n=0;n<NV;++n)
				bw(begin_tri+i*im*NV+m*NV+n)+= tm*NV;
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
	int gindx,eind,find,iind,sgn,msgn,mode;
	int kn = basis::tri(log2p)->tm()*NV;
	Array<FLT,2> K(kn,kn);
	Array<int,1> loc_to_glo(kn);
	
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
			for(int m = ; m < basis::tri(log2p)->im(); ++m) {
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
	int kcol = 0;
	FLT dw = 1.0e-6;  //dw=sqrt(eps/l2_norm(q))
	
	
	ugtouht(tind);
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tri(log2p)->tm();++i)
		for(int n=0;n<NV;++n)
			Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);

	for(int mode = 0; mode < basis::tri(log2p)->tm(); ++mode){
		for(int var = 0; var < NV; ++var){
			uht(var)(mode) += dw;
			
			element_rsdl(tind,0,uht,lf_re,lf_im);
			for(int i=0;i<basis::tri(log2p)->tm();++i)
				for(int n=0;n<NV;++n)
					R(n)(i)=lf_re(n)(i)+lf_im(n)(i);

			int krow = 0;
			for(int i=0;i<basis::tri(log2p)->tm();++i)
				for(int n=0;n<NV;++n)
					K(krow++,kcol) = (R(n)(i)-Rbar(n)(i))/dw;
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
	
	for (int m = 0; m < basis::tri(log2p)->tm(); ++m) {
		for (int n = 0; n < NV; ++n) {
			lf_im(n)(m)=0.0;
			lf_re(n)(m)=0.0;
		}
	}
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tri(log2p)->tm();++i)
		for(int n=0;n<NV;++n)
			lclres(ind++)=lf_re(n)(i)+lf_im(n)(i);	
	
	return;
}

void tri_hp::apply_neumman() {
	int tm = 2+basis::tri(log2p)->sm();
	int kn = tm*NV;
	int kcol,krow,eind,ind;
	FLT dw = 1.0e-6;// fix me temp make a global value?
	FLT sgn,msgn;
	Array<FLT,2> R(NV,tm),Rbar(NV,tm);
	Array<FLT,2> K(kn,kn);
	Array<FLT,1> lclres(kn);
	Array<int,1> loc_to_glo(kn);
	
	for(int i=0;i<nebd;++i){
		for(int j=0;j<ebdry(i)->nseg;++j){
			eind = ebdry(i)->seg(j).gindx;
			ugtouht1d(eind);
			hp_ebdry(i)->element_rsdl(eind,0);
			kcol = 0;
			ind = 0;
			for(int k=0;k<tm;++k){
				for(int n=0;n<NV;++n){
					Rbar(n,k) = lf(n)(k);
					lclres(ind++) = lf(n)(k);
				}
			}
			for(int mode = 0; mode < tm; ++mode){
				for(int var = 0; var < NV; ++var){
					uht(var)(mode) += dw;
					
					hp_ebdry(i)->element_rsdl(eind,0);
					for(int k=0;k<tm;++k)
						for(int n=0;n<NV;++n)
							R(n,k)=lf(n)(k);
					
					krow = 0;
					for(int k=0;k<tm;++k)
						for(int n=0;n<NV;++n)
							K(krow++,kcol) = (R(n,k)-Rbar(n,k))/dw;
					++kcol;
					
					uht(var)(mode) -= dw;
				}
			}
			ind = 0;
			for (int m = 0; m < 2; ++m) {
				int gindx = NV*seg(eind).pnt(m);
				for (int n = 0; n < NV; ++n)
					loc_to_glo(ind++) = gindx+n;
			}		
			
			/* EDGE MODES */
			if (basis::tri(log2p)->p() > 1) {
				int gbl_eind = npnt*NV + eind*basis::tri(log2p)->sm()*NV;
				sgn = seg(eind).sgn(k);
				msgn = 1.0;
				for (int m = 0; m < basis::tri(log2p)->sm(); ++m) {
					for(int n = 0; n < NV; ++n) {
						for(int j = 0; j < kn; ++j) {
							K(ind,j) *= msgn;
							K(j,ind) *= msgn;
						}
						lclres(ind) *= msgn;

						loc_to_glo(ind++) = eind + m*NV + n;
					}
					msgn *= sgn;
				}
			}
			

			MatSetValues(petsc_J,kn,loc_to_glo.data(),kn,loc_to_glo.data(),K.data(),ADD_VALUES);
			VecSetValues(petsc_f,kn,loc_to_glo.data(),lclres.data(),ADD_VALUES);

		}
	}
	
	return;
}
#endif