/*
 *  sparse.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 9/14/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "ins/tet_hp_ins.h"
//#include "gmres.h"
#include "hp_boundary.h"

#ifndef petsc
/* compressed row storage */
void tet_hp::insert_sparse(int row, int col, FLT value, bool compressed_col ){
	
//	if (fabs(value)<10e-15) {
//		return;
//	}

	/* to make compressed column switch row and column */
	if (compressed_col) {
		int temp = row;
		row = col;
		col = temp;
	}
		
	/* add value to already existing element in sparse */
	for(int i = sparse_ptr(row); i < sparse_ptr(row+1); ++i){
		if(sparse_ind(i) == col){
			sparse_val(i) += value;
			return;
		}
	}
	
	int find_ind = -1;

	/* if not pre-exisiting element then add to list in correct location */
	for(int i = sparse_ptr(row); i < sparse_ptr(row+1); ++i){
		if(sparse_ind(i) > col){
			find_ind = i;
			break;
		}
	}
	
//	/* bandwidth known so dont do this stuff*/
//	/* no column bigger so insert at end */
//	if(find_ind == -1)
//		find_ind = sparse_ptr(row+1);
//	
//	++number_sparse_elements;
//	
//	/* shift start of column index by one */
//	for(int i = row+1; i < size_sparse_matrix+1; ++i)
//		++sparse_ptr(i);
//	
//	/* slide all non-diagonal entries after insertion column to the right (slow but only done once) */
//	for(int i = number_sparse_elements; i > find_ind; --i){
//		sparse_val(i)=sparse_val(i-1);
//		sparse_ind(i)=sparse_ind(i-1);
//	}
	
	/* slide all entries after insertion column to the right*/
	for(int i = sparse_ptr(row+1)-1; i > find_ind; --i){
		sparse_val(i)=sparse_val(i-1);
		sparse_ind(i)=sparse_ind(i-1);
	}
	
	/* insert new element into sparse matrix */
	sparse_val(find_ind) = value;
	sparse_ind(find_ind) = col;	
	
	return;	
}




/* zero out sparse matrix but keep same structure */
void tet_hp::zero_sparse(){
	
	for (int i = 0; i < number_sparse_elements; ++i)
		sparse_val(i) = 0.0;
	
	for (int i = 0; i < size_sparse_matrix; ++i)
		res_vec(i) = 0.0;

	
	return;
}

void tet_hp::initialize_sparse(){
	
	size_sparse_matrix = (npnt+nseg*basis::tet(log2p).em+ntri*basis::tet(log2p).fm+ntet*basis::tet(log2p).im)*NV;

	sparse_ptr.resize(size_sparse_matrix+1);
	find_sparse_bandwidth();
	number_sparse_elements = sparse_ptr(size_sparse_matrix);
	
	cout << "number of sparse elements "<< number_sparse_elements << endl;
	cout << "number of unknowns "<< size_sparse_matrix << endl;

	sparse_val.resize(number_sparse_elements);
	sparse_ind.resize(number_sparse_elements);
	
	sparse_ind = 10*size_sparse_matrix; // some number bigger than size of matrix

	res_vec.resize(size_sparse_matrix);
	ug_vec.resize(size_sparse_matrix);
	
	return;
}

#endif

/* finds number of nonzeros per row tested up to p=4 */ 
void tet_hp::find_sparse_bandwidth(){
	
	Array<int,1> bw(size_sparse_matrix);
	
	int em=basis::tet(log2p).em;
	int fm=basis::tet(log2p).fm;
	int im=basis::tet(log2p).im;
	int tm=basis::tet(log2p).tm;
	
	int begin_seg = npnt*NV;
	int begin_tri = begin_seg+nseg*em*NV;
	int begin_tet = begin_tri+ntri*fm*NV;
	
	bw = 0;
	
	for(int i=0; i<npnt; ++i){		
		for(int n=0;n<NV;++n)
			bw(i*NV+n) += NV*(pnt(i).nnbor*fm+pnt(i).nspk+pnt(i).ntri*em+1);		
	}
	
	
	for(int i=0; i<nseg; ++i){		
		for(int j=0;j<2;++j)
			for(int n=0;n<NV;++n)
				bw(seg(i).pnt(j)*NV+n) += em*NV;
		
		for(int m=0;m<em;++m)
			for(int n=0;n<NV;++n)
				bw(begin_seg+i*em*NV+m*NV+n) += NV*((2*fm+em)*seg(i).nnbor+em+(2*em+1)*seg(i).nspk+2);		
	}
	
	
	for(int i=0; i<ntri; ++i){		
		for(int j=0;j<3;++j)
			for(int n=0;n<NV;++n)
				bw(tri(i).pnt(j)*NV+n) += fm*NV;
		
		for(int j=0;j<3;++j)
			for(int m=0;m<em;++m)
				for(int n=0;n<NV;++n)
					bw(begin_seg+tri(i).seg(j)*em*NV+m*NV+n) += fm*NV;
		
		if(tri(i).tet(1) > 0){
			for(int m=0;m<fm;++m)
				for(int n=0;n<NV;++n)
					bw(begin_tri+i*fm*NV+m*NV+n) += NV*(7*fm+9*em+5);
		}
		else {
			for(int m=0;m<fm;++m)
				for(int n=0;n<NV;++n)
					bw(begin_tri+i*fm*NV+m*NV+n) += NV*(4*fm+6*em+4);
		}
	}
	
	for(int i=0; i<ntet; ++i){		
		for(int j=0;j<4;++j)
			for(int n=0;n<NV;++n)
				bw(tet(i).pnt(j)*NV+n) += im*NV;
		
		for(int j=0;j<6;++j)
			for(int m=0;m<em;++m)
				for(int n=0;n<NV;++n)
					bw(begin_seg+tet(i).seg(j)*em*NV+m*NV+n) += im*NV;
		
		for(int j=0;j<4;++j)
			for(int m=0;m<fm;++m)
				for(int n=0;n<NV;++n)
					bw(begin_tri+tet(i).tri(j)*fm*NV+m*NV+n) += im*NV;
		
		for(int m=0;m<im;++m)
			for(int n=0;n<NV;++n)
				bw(begin_tet+i*im*NV+m*NV+n)+= tm*NV;
	}
	

#ifdef petsc
	MatCreateSeqAIJ(PETSC_COMM_SELF,size_sparse_matrix,size_sparse_matrix,PETSC_NULL,bw.data(),&petsc_J);
	//MatCreateMPIAIJ(comm,size_sparse_matrix,size_sparse_matrix,global_size,global_size,int d nz,int *d nnz, int o nz,int *o nnz,Mat *A);
#else
	sparse_ptr(0)=0;
	for (int i=1; i<size_sparse_matrix+1; ++i) 
		sparse_ptr(i) = bw(i-1)+sparse_ptr(i-1);
#endif

	
	
	return;
}


/* clears row in sparse matrix, inserts 1.0 on diagonal, and inserts 0.0 in the residual */
void tet_hp::sparse_dirichlet(int ind, bool compressed_col){
	
#ifdef petsc
	const PetscInt row = ind;
	PetscScalar zero = 0.0;
	PetscScalar one = 1.0;

	/* apply dirichlet by inserting zero in f */
	VecSetValues(petsc_f,1,&row,&zero,INSERT_VALUES);
	dirichlet_rows(row_counter++)=ind;
	
	return;
	
//	PetscInt          ncols;
//	const PetscInt    *cols;
//	const PetscScalar *vals;
//	PetscTruth        flg1=PETSC_FALSE;
//	PetscScalar       *zeros;
//	
//	MatGetRow(petsc_J,row,&ncols,&cols,&vals); 
//	
//	MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);
//	MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);	
//	
//	PetscMalloc(sizeof(PetscScalar)*(ncols+1),&zeros);
//	PetscMemzero(zeros,(ncols+1)*sizeof(PetscScalar));
//	PetscOptionsHasName(PETSC_NULL, "-set_row_zero", &flg1);
//
//	MatSetValues(petsc_J,1,&row,ncols,cols,zeros,INSERT_VALUES);
//
//	
//	
//	return;
	
	
	MatZeroRows(petsc_J,1,&row,1.0);

	
	return;
	
//	PetscInt ncols;
//	const PetscInt *cols;
//	const PetscScalar *vals;
//   
//	MatGetRow(petsc_J,row,&ncols,&cols,&vals);
//
//	Array<double,1> values(ncols);
//	values = 0.0;
//	
//	MatSetValues(petsc_J,1,&row,ncols,cols,values.data(),INSERT_VALUES);
//	MatSetValues(petsc_J,1,&row,1,&row,&one,INSERT_VALUES);
	
//	MatRestoreRow(petsc_J,row,&ncols,&cols,&vals);

//	return;
//	MatRestoreRow(petsc_J,row,&ncols,&cols,&values.data());
	
	/* zero out row in Jacobian and insert 1.0 on diagonal */
	//MatZeroRowsLocal(petsc_J,1,&pind,one);
	//MatZeroRowsIS(petsc_J,IS rows,one);
	
//	/* hacked */
//	PetscInt *cols;
//	PetscScalar *vals;
//	PetscInt ncol = size_sparse_matrix;
//
//	for(int i = 0; i < size_sparse_matrix; ++i){
//		cols[i]=i;
//		vals[i]=0.0;
//	}
//	
//	MatSetValues(petsc_J,1,&row,ncol,cols,vals,INSERT_VALUES);
//	MatSetValues(petsc_J,1,&row,1,&row,&one,INSERT_VALUES);

#else
	res_vec(ind) = 0.0;
	
	if(compressed_col){
		/* compressed column: only works for symmetric sparsity pattern */
		for(int i = sparse_ptr(ind); i < sparse_ptr(ind+1); ++i){
			for(int j = sparse_ptr(sparse_ind(i)); j < sparse_ptr(sparse_ind(i)+1); ++j){
				if(ind == sparse_ind(j)){
					sparse_val(j) = 0.0;
					break;
				}
			}	
			if(ind == sparse_ind(i))
				sparse_val(i) = 1.0;
		}
	}
	else {
		/* compressed row */
		for(int i = sparse_ptr(ind); i < sparse_ptr(ind+1); ++i){
			sparse_val(i) = 0.0;
			if(ind == sparse_ind(i))
				sparse_val(i) = 1.0;
		}
	}
#endif
	
	return;
}

#ifdef petsc
void tet_hp::create_jacobian_residual(){
	
	bool compressed_column = false;
	PetscLogDouble time1,time2;
	
	//find_sparse_bandwidth();//shouldn't have to do this everytime
	
	VecSet(petsc_f,0.0);
//	MatZeroEntries(petsc_J);
	cout << "zeroed enries"<< endl;
	/* insert values into jacobian matrix J */
	
	PetscGetTime(&time1);
	create_jacobian(compressed_column);
	PetscGetTime(&time2);

	cout << "jacobian made " << time2-time1 << " seconds" << endl;

	/* insert values into residual f (every processor will do this need to do it a different way) */
	create_rsdl();		
	cout << "residual made"<< endl;

	/* apply neumman bc's */
	apply_neumman(compressed_column);
	
	cout << "neumman BC's applied"<< endl;

	MatAssemblyBegin(petsc_J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(petsc_J,MAT_FINAL_ASSEMBLY);	

	cout << "jacobian assembled"<< endl;

	
	PetscGetTime(&time1);
	/* apply dirichlet boundary conditions to sparse matrix and vector */
	for(int j = 0; j < nfbd; ++j)
		hp_fbdry(j)->apply_sparse_dirichlet(compressed_column);		
	PetscGetTime(&time2);

	cout << "dirichlet BC's applied "<< time2-time1 <<" seconds"<< endl;

	VecAssemblyBegin(petsc_f);
	VecAssemblyEnd(petsc_f);
	cout << "petsc_f assembled "<< endl;

	return;
	
}

#else

void tet_hp::create_jacobian_residual(){
	
	bool compressed_column = true;
	
	zero_sparse();

	/* insert values into jacobian matrix J (every processor will do this need to do it a different way) */
	create_jacobian(compressed_column);
	
	/* insert values into residual f (every processor will do this need to do it a different way) */
	create_rsdl();		
	
	/* apply neumman bc's */
	apply_neumman(compressed_column);

	
	/* apply dirichlet boundary conditions to sparse matrix and vector */
	for(int j = 0; j < nfbd; ++j)
		hp_fbdry(j)->apply_sparse_dirichlet(compressed_column);		

	
	return;
	
}

#endif


void tet_hp::create_jacobian(bool jac_tran) {
	int gindx,eind,find,iind,sgn,msgn,mode;
	int kn = basis::tet(log2p).tm*NV;
	Array<FLT,2> K(kn,kn);
	Array<int,1> loc_to_glo(kn);
	
	for(int tind = 0; tind < ntet; ++tind){	
		
		int ind = 0;

		create_local_jacobian_matrix(tind, K);
		
		for (int m = 0; m < 4; ++m) {
			gindx = NV*tet(tind).pnt(m);
			for (int n = 0; n < NV; ++n)
				loc_to_glo(ind++) = gindx+n;
		}		
		
		/* EDGE MODES */
		if (basis::tet(log2p).p > 1) {
			for(int i = 0; i < 6; ++i) {
				eind = npnt*NV + tet(tind).seg(i)*basis::tet(log2p).em*NV;
				sgn = tet(tind).sgn(i);
				msgn = 1;
				for (int m = 0; m < basis::tet(log2p).em; ++m) {
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
		
		/* FACE MODES */
		if (basis::tet(log2p).p > 2) {
			for(int i = 0; i < 4; ++i){
				sgn = -tet(tind).rot(i);
				find = npnt*NV+nseg*basis::tet(log2p).em*NV+tet(tind).tri(i)*basis::tet(log2p).fm*NV;
				mode = 0;
				msgn = 1;		
				for(int m = 1; m <= basis::tet(log2p).em-1; ++m) {
					for(int j = 1; j <= basis::tet(log2p).em-m; ++j){
						for(int n = 0; n < NV; ++n){
							for(int k = 0; k < kn; ++k){
								K(ind,k) *= msgn;
								K(k,ind) *= msgn;
							}
							loc_to_glo(ind++) = find+mode*NV+n;
						}
						++mode;
					}
					msgn *= sgn;
				}
			}		
		}		
		
		/* INTERIOR MODES */
		iind = npnt*NV+nseg*basis::tet(log2p).em*NV+ntri*basis::tet(log2p).fm*NV+tind*basis::tet(log2p).im*NV;
		for(int m = 0; m < basis::tet(log2p).im; ++m) {
			for(int n = 0; n < NV; ++n){
				for(int k = 0; k < kn; ++k)
					K(ind,k) = K(ind,k);
				loc_to_glo(ind++) = iind+m*NV+n;
			}
		}
		
#ifdef petsc

		PetscErrorCode ierr;
		ierr = MatSetValues(petsc_J,kn,loc_to_glo.data(),kn,loc_to_glo.data(),K.data(),ADD_VALUES);
#else
		for(int i = 0; i < kn; ++i)
			for(int j = 0; j < kn; ++j)
				insert_sparse(loc_to_glo(i), loc_to_glo(j), K(i,j),jac_tran);

#endif

				
				
	}	

	return;
}

void tet_hp::create_local_jacobian_matrix(int tind, Array<FLT,2> &K) {
	Array<TinyVector<FLT,MXTM>,1> R(NV),Rbar(NV),lf_re(NV),lf_im(NV);
	int kcol = 0;
	FLT dw = 1.0e-6;  //dw=sqrt(eps/l2_norm(q))
	
	
	ugtouht(tind);
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tet(log2p).tm;++i)
		for(int n=0;n<NV;++n)
			Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);

	for(int mode = 0; mode < basis::tet(log2p).tm; ++mode){
		for(int var = 0; var < NV; ++var){
			uht(var)(mode) += dw;
			
			element_rsdl(tind,0,uht,lf_re,lf_im);
			for(int i=0;i<basis::tet(log2p).tm;++i)
				for(int n=0;n<NV;++n)
					R(n)(i)=lf_re(n)(i)+lf_im(n)(i);

			int krow = 0;
			for(int i=0;i<basis::tet(log2p).tm;++i)
				for(int n=0;n<NV;++n)
					K(krow++,kcol) = (R(n)(i)-Rbar(n)(i))/dw;
			++kcol;
			
			uht(var)(mode) -= dw;
		}
	}	

	return;
}

void tet_hp::create_rsdl() {
	int gindx,eind,find,iind,sgn,msgn,mode;
	
	int kn = NV*basis::tet(log2p).tm;
	Array<FLT,1> lclres(kn);
	Array<int,1> loc_to_glo(kn);
	
	for(int tind = 0; tind < ntet; ++tind){	
		
		int ind = 0;
		
		create_local_rsdl(tind, lclres);
		
		for (int m = 0; m < 4; ++m) {
			gindx = NV*tet(tind).pnt(m);
			for (int n = 0; n < NV; ++n)
				loc_to_glo(ind++) = gindx+n;
		}		
		
		/* EDGE MODES */
		if (basis::tet(log2p).p > 1) {
			for(int i = 0; i < 6; ++i) {
				eind = npnt*NV + tet(tind).seg(i)*basis::tet(log2p).em*NV;
				sgn = tet(tind).sgn(i);
				msgn = 1;
				for (int m = 0; m < basis::tet(log2p).em; ++m) {
					for(int n = 0; n < NV; ++n) {
						lclres(ind) *= msgn;
						loc_to_glo(ind++) = eind + m*NV + n;
					}
					msgn *= sgn;
				}
			}
		}
		
		/* FACE MODES */
		if (basis::tet(log2p).p > 2) {
			for(int i = 0; i < 4; ++i){
				sgn = -tet(tind).rot(i);
				find = npnt*NV+nseg*basis::tet(log2p).em*NV+tet(tind).tri(i)*basis::tet(log2p).fm*NV;
				mode = 0;
				msgn = 1;		
				for(int m = 1; m <= basis::tet(log2p).em-1; ++m) {
					for(int j = 1; j <= basis::tet(log2p).em-m; ++j){
						for(int n = 0; n < NV; ++n){
							lclres(ind) *= msgn;							
							loc_to_glo(ind++) = find+mode*NV+n;
						}
						++mode;
					}
					msgn *= sgn;
				}
			}		
		}		
		
		/* INTERIOR MODES */
		iind = npnt*NV+nseg*basis::tet(log2p).em*NV+ntri*basis::tet(log2p).fm*NV+tind*basis::tet(log2p).im*NV;
		for(int m = 0; m < basis::tet(log2p).im; ++m) 
			for(int n = 0; n < NV; ++n)
				loc_to_glo(ind++) = iind+m*NV+n;
			
		
	
#ifdef petsc		
		PetscErrorCode ierr;
		ierr = VecSetValues(petsc_f,kn,loc_to_glo.data(),lclres.data(),ADD_VALUES);

#else
		for(int i = 0; i < kn; ++i)
			res_vec(loc_to_glo(i))+=lclres(i);	
#endif
		
		
		
	}
	
	return;
}

void tet_hp::create_local_rsdl(int tind, Array<FLT,1> &lclres) {
	Array<TinyVector<FLT,MXTM>,1> lf_re(NV),lf_im(NV);

	int ind = 0;

	ugtouht(tind);
	
	for (int m = 0; m < basis::tet(log2p).tm; ++m) {
		for (int n = 0; n < NV; ++n) {
			lf_im(n)(m)=0.0;
			lf_re(n)(m)=0.0;
		}
	}
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tet(log2p).tm;++i)
		for(int n=0;n<NV;++n)
			lclres(ind++)=lf_re(n)(i)+lf_im(n)(i);	
	
	return;
}

void tet_hp::apply_neumman(bool jac_tran) {
	int fm = 3+3*basis::tet(log2p).em+basis::tet(log2p).fm;
	int kn = fm*NV;
	int kcol,krow,find,ind;
	FLT dw = 1.0e-6;
	FLT sgn,msgn;
	Array<FLT,2> R(NV,fm),Rbar(NV,fm);
	Array<FLT,2> K(kn,kn);
	Array<FLT,1> lclres(kn);
	Array<int,1> loc_to_glo(kn);
	
	for(int i=0;i<nfbd;++i){
		for(int j=0;j<fbdry(i)->ntri;++j){
			find = fbdry(i)->tri(j).gindx;
			ugtouht2d(find);
			hp_fbdry(i)->element_rsdl(find,0);
			kcol = 0;
			ind = 0;
			for(int k=0;k<fm;++k){
				for(int n=0;n<NV;++n){
					Rbar(n,k) = lf(n)(k);
					lclres(ind++) = lf(n)(k);
				}
			}
			for(int mode = 0; mode < fm; ++mode){
				for(int var = 0; var < NV; ++var){
					uht(var)(mode) += dw;
					
					hp_fbdry(i)->element_rsdl(find,0);
					for(int k=0;k<fm;++k)
						for(int n=0;n<NV;++n)
							R(n,k)=lf(n)(k);
					
					krow = 0;
					for(int k=0;k<fm;++k)
						for(int n=0;n<NV;++n)
							K(krow++,kcol) = (R(n,k)-Rbar(n,k))/dw;
					++kcol;
					
					uht(var)(mode) -= dw;
				}
			}
			ind = 0;
			for (int m = 0; m < 3; ++m) {
				int gindx = NV*tri(find).pnt(m);
				for (int n = 0; n < NV; ++n)
					loc_to_glo(ind++) = gindx+n;
			}		
			
			/* EDGE MODES */
			if (basis::tet(log2p).p > 1) {
				for(int k = 0; k < 3; ++k) {
					int eind = npnt*NV + tri(find).seg(k)*basis::tet(log2p).em*NV;
					sgn = tri(find).sgn(k);
					msgn = 1.0;
					for (int m = 0; m < basis::tet(log2p).em; ++m) {
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
			}
			
			/* FACE MODE */
			int gbl_find = npnt*NV+nseg*basis::tet(log2p).em*NV+find*basis::tet(log2p).fm*NV;
			for(int m = 0; m < basis::tet(log2p).fm; ++m) 
				for(int n = 0; n < NV; ++n)
					loc_to_glo(ind++) = gbl_find+m*NV+n;				
#ifdef petsc				

			MatSetValues(petsc_J,kn,loc_to_glo.data(),kn,loc_to_glo.data(),K.data(),ADD_VALUES);
			VecSetValues(petsc_f,kn,loc_to_glo.data(),lclres.data(),ADD_VALUES);
				
#else
			for(int k = 0; k < kn; ++k){
				for(int m = 0; m < kn; ++m)
					insert_sparse(loc_to_glo(k), loc_to_glo(m), K(k,m),jac_tran);
				
				res_vec(loc_to_glo(k)) += lclres(k);
			}
#endif
		}
	}
	
	return;
}

#ifndef petsc
void tet_hp::vec_to_ug(){
	int ind = 0;
	
	for(int i = 0; i < npnt; ++i)
		for(int n = 0; n < NV; ++n)
			ug.v(i,n) = ug_vec(ind++);
	
	if(!basis::tet(log2p).em) return;
	
	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tet(log2p).em; ++m)
			for(int n = 0; n < NV; ++n)
				ug.e(i,m,n) = ug_vec(ind++);
	
	if(!basis::tet(log2p).fm) return;

	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tet(log2p).fm; ++m)
			for(int n = 0; n < NV; ++n)
				ug.f(i,m,n) = ug_vec(ind++);		

	if(!basis::tet(log2p).im) return;
	
	for(int i = 0; i < ntet; ++i)
		for(int m = 0; m < basis::tet(log2p).im; ++m)
			for(int n = 0; n < NV; ++n)
				ug.i(i,m,n) = ug_vec(ind++);	
	
	return;	
}

void tet_hp::ug_to_vec(){
	
	int ind = 0;
	
	for(int i = 0; i < npnt; ++i)
		for(int n = 0; n < NV; ++n)
			ug_vec(ind++) = ug.v(i,n);
	
	if(!basis::tet(log2p).em) return;

	for(int i = 0; i < nseg; ++i)
		for(int m = 0; m < basis::tet(log2p).em; ++m)
			for(int n = 0; n < NV; ++n)
				ug_vec(ind++) = ug.e(i,m,n);
	
	if(!basis::tet(log2p).fm) return;

	for(int i = 0; i < ntri; ++i)
		for(int m = 0; m < basis::tet(log2p).fm; ++m)
			for(int n = 0; n < NV; ++n)
				ug_vec(ind++) = ug.f(i,m,n);;		
	
	if(!basis::tet(log2p).im) return;

	for(int i = 0; i < ntet; ++i)
		for(int m = 0; m < basis::tet(log2p).im; ++m)
			for(int n = 0; n < NV; ++n)
				ug_vec(ind++) = ug.i(i,m,n);	
	
	return;	
}
#endif





//#ifdef row_index
///* row-indexed sparse storage */
//void tet_hp::insert_sparse(int row, int col, FLT value){
//	
//	int find_col = -1;
//	
//	/* insert diagonal element */
//	if(row == col){
//		sa(row) += value;
//		return;
//	}
//	
//	/* add value to already existing element in sparse */
//	for(int i = ija(row); i < ija(row+1); ++i){
//		if(ija(i) == col){
//			sa(i) += value;
//			return;
//		}
//	}
//	
//	/* if not diagonal entry or pre-exisiting element then add to list in correct location */
//	for(int i = ija(row); i < ija(row+1); ++i){
//		if(ija(i) > col){
//			find_col = i;
//			break;
//		}
//	}
//	/* no column bigger so insert at end */
//	if(find_col == -1)
//		find_col = ija(row+1);
//		
//	++number_sparse_elements;
//	
//	/* shift start of column index by one */
//	for(int i = row+1; i < size_sparse_matrix+1; ++i)
//		++ija(i);
//		
//	/* slide all non-diagonal entries after insertion column to the right (slow but only done once) */
//		for(int i = number_sparse_elements; i > find_col; --i){
//			sa(i)=sa(i-1);
//			ija(i)=ija(i-1);
//		}
//	
//	/* insert new element into sparse matrix */
//	sa(find_col) = value;
//	ija(find_col) = col;
//	
//	
//	return;	
//}
//
///* zero out sparse matrix but keep same structure */
//void tet_hp::zero_sparse(){
//	
//	for (int i = 0; i < number_sparse_elements; ++i){
//		res_vec(i) = 0.0;
//		sa(i) = 0.0;
//	}
//	sa(number_sparse_elements) = 0.0;
//	
//	return;
//}
//
//void tet_hp::initialize_sparse(){
//	
//	size_sparse_matrix = (npnt+nseg*em0+ntri*fm0+ntet*im0)*NV;
//	
//	/* sparse matrix allocation */
//	ija.resize(MXTM*NV*ntet);//too much storage resize later
//	sa.resize(MXTM*NV*ntet);
//	number_sparse_elements = size_sparse_matrix;
//	sa = 0.0;
//	/* creates sparse matrix with zeros on diagonal */
//	for(int i = 0; i < number_sparse_elements+1; ++i)
//		ija(i) = size_sparse_matrix+1;
//	return;
//}
//#endif


//void tet_hp::sparse_matrix_multiply(FLT *x, FLT *b){
//
//	for(int i = 0; i < size_sparse_matrix; ++i){
//		b(i) = sa(i)*x(i);
//		for(int j = ija(i); j < ija(i+1); ++j)
//			b(i) += sa(j)*x(ija(j));
//	}
//	
//	return;
//}
//
//
//void tet_hp::sparse_matrix_multiply_transpose(FLT *x, FLT *b){
//	
//	for(int i = 0; i < size_sparse_matrix; ++i)
//		b(i) = sa(i)*x(i);
//
//	for(int i = 0; i < size_sparse_matrix; ++i){
//		for(int j = ija(i); j < ija(i+1); ++j){
//			int k = ija(j)
//			b(k) += sa(j)*x(i);
//		}
//	}
//	
//	return;
//}