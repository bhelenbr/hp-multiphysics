
/*
 *  pod_generate.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 1/18/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#include "pod_generate.h"
#include <myblas.h>
#ifdef libbinio
#include <libbinio/binfile.h>
#endif
#include <netcdf.h>

#define ERR(e,logp) {*logp << "netCDF error " <<  nc_strerror(e); sim::abort(__LINE__,__FILE__,logp);}

#define LOW_NOISE_DOT


extern "C" {
	double dlamch_(const char *cmach);
	/* Subroutine */ int dspevx_(char *jobz, char *range, char *uplo, int *n,
															 double *ap, double *vl, double *vu, int *il, int *
															 iu, double *abstol, int *m, double *w, double *z__,
															 int *ldz, double *work, int *iwork, int *ifail,
															 int *info);
	
	/* Subroutine */ int dspgvx_(int *ITYPE, char *jobz, char *range, char *uplo, int *n,
															 double *ap, double *ab, double *vl, double *vu, int *il,
															 int *iu, double *abstol, int *m, double *w, double *z__,
															 int *ldz, double *work, int *iwork, int *ifail,
															 int *info);
	
	/* Subroutine */ int dspgv_(int *ITYPE, char *jobz, char *uplo, int *n,
															double *ap, double *ab, double *w, double *z__,
															int *ldz, double *work, int *info);

	/* Subroutine */ int dgesvd_(char *JOBU, char *JOBVT, int *m, int *n, double *A, int *LDA, double *s, double *U, int *LDU,
															double *VT, int *LDVT, double *work, int *iwork, int *info);	

	/* Subroutine */ //int dspev_(char *jobz, char *uplo, int *n, double *ap, double *w, double *z__, int *ldz, double *work, 
								//							int *info);	
}


template<class BASE> void pod_generate<BASE>::init(input_map& inmap, shared_ptr<block_global> gin) {
	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	std::istringstream instr;
	int i;
	
	/* Initialize base class */
	BASE::init(inmap,gin);
	
#ifdef DIRECT_METHOD
	*BASE::gbl->log << "DIRECT_METHOD is defined\n";
#else
	*BASE::gbl->log << "DIRECT_METHOD is not defined\n";
#endif

#ifdef USING_MASS_MATRIX // This uses the mass matrix to form inner products rather than numerically evaluate them
	*BASE::gbl->log << "USING_MASS_MATRIX is defined\n";
#else
	*BASE::gbl->log << "USING_MASS_MATRIX is not defined\n";
#endif
	
#ifdef LOW_NOISE_DOT // This uses the mass matrix to form inner products rather than numerically evaluate them
	*BASE::gbl->log << "LOW_NOISE_DOT is defined\n";
#else
	*BASE::gbl->log << "LOW_NOISE_DOT is not defined\n";
#endif
	
	std::string mystring;
	inmap.getlinewdefault(BASE::gbl->idprefix +"_groups",mystring,std::string("0"));
	istringstream data;
	data.str(mystring);
	data >> pod_id;   // First group is pod_id.  Other groups are for pod boundary conditions
	
	if (!inmap.get("snapshots",nsnapshots)) {
		*BASE::gbl->log << "# Couldn't find number of snapshots" << std::endl;
		sim::abort(__LINE__,__FILE__,BASE::gbl->log);
	}
	inmap.getwdefault("restart_interval",restart_interval,1);
	inmap.getwdefault("restart",restartfile,1);
	
	inmap.getwdefault("coefficient_start",coefficient_start,restartfile);
	inmap.getwdefault("coefficient_interval",coefficient_interval,restart_interval);
	inmap.getwdefault("coefficient_end",coefficient_end,restart_interval*nsnapshots +restartfile);
	
	std::string job_name;
	inmap.getlinewdefault("job",job_name,std::string("snapshot"));  //job = R-POD
	if (job_name == "direct") {
		job_id = direct;
	}
	else if (job_name == "snapshot") {
		job_id = snapshot;
	}
	else if (job_name == "mPOD") {
		job_id = mPOD;
	}
	else {
		*BASE::gbl->log << "Undefined job." << std::endl;
		sim::abort(__LINE__,__FILE__,BASE::gbl->log);
	}
	
	inmap.getwdefault("podmodes",nmodes,nsnapshots);
	if ( nmodes > nsnapshots ) {
		*BASE::gbl->log << "#Number of modes is more than number of snapshots" << std::endl;
		sim::abort(__LINE__,__FILE__,BASE::gbl->log);
	}
	nmodes = MAX(nmodes,2);
	
	inmap.getwdefault("ndeflation",ndeflation,1);  // This is if you want to use the recursive method (1 means no, was LOW_NOISE)
	if ( nmodes % ndeflation != 0 ) {
		*BASE::gbl->log << "# Number of modes and number of deflations are not compatible." << std::endl;
		sim::abort(__LINE__,__FILE__,BASE::gbl->log);
	}
	nmodes_per_deflation = nmodes/ndeflation;
	
	inmap.getwdefault("MinEigen",MinEigen,1e-8);

	inmap.getwdefault("M",M,nsnapshots);
	if ( M > nsnapshots ) {
		*BASE::gbl->log << "#MaxSize can't be more than number of snapshots. MaxSize replaced with snapshots." << std::endl;
		M = nsnapshots;
	}

	
	/* THIS IS TO CHANGE THE WAY SNAPSHOT MATRIX ENTRIES ARE FORMED */
	scaling.resize(BASE::NV);
	scaling = 1;
	if (inmap.getline(BASE::gbl->idprefix + "_scale_vector",linebuff) || inmap.getline("scale_vector",linebuff)) {
		instr.str(linebuff);
		for(i=0;i<BASE::NV;++i)
			instr >> scaling(i);
	}
	
	modes.resize(nmodes);
	for(i=0;i<nmodes;++i) {
		modes(i).v.resize(BASE::maxpst,BASE::NV);
		modes(i).s.resize(BASE::maxpst,BASE::sm0,BASE::NV);
		modes(i).i.resize(BASE::maxpst,BASE::im0,BASE::NV);
	}
	
#ifdef USING_MASS_MATRIX
	const int sm = basis::tri(BASE::log2p)->sm();
	const int im = basis::tri(BASE::log2p)->im();
	const int ndofs = BASE::npnt +BASE::nseg*sm +BASE::ntri*im;
	int ntotal = ndofs*BASE::NV;
	mass_times_snapshot.resize(ntotal);
#endif
	
#ifdef POD_BDRY
	pod_ebdry.resize(BASE::nebd);
	for (int i=0;i<BASE::nebd;++i) {
		pod_ebdry(i) = new pod_gen_edge_bdry<BASE>(*this,*BASE::ebdry(i));
		pod_ebdry(i)->init(inmap);
	}
#endif
	
	return;
}

template<class BASE> void pod_generate<BASE>::tadvance() {
	std::string filename;
	std::ostringstream nstr;
	
	/* MAKE FILES FOR VOLUME MODE SNAPSHOT-PROJECTION */
	/* THIS DOES NOTHING EXCEPT CREATE DUPLICATE FILES IF POD_BDRY IS OFF */
	for(int k=0;k<nsnapshots;++k) {
		nstr.str("");
		nstr << k*restart_interval +restartfile << std::flush;
		filename = "rstrt" +nstr.str() +"_d0";
		BASE::input(filename, BASE::reload_type);
		
#ifdef POD_BDRY
		/* ZERO SNAPSHOTS ON POD BOUNDARY'S */
		for (int i=0;i<BASE::nebd;++i)
			pod_ebdry(i)->zero_bdry(BASE::ug);
#endif
		
		nstr.str("");
		nstr << k << std::flush;
		filename = "temp" +nstr.str();
		BASE::output(filename, BASE::output_type(1));
	}
	
	BASE::tadvance();
	
	int psi1dcounter;
	Array<FLT,1> psimatrix(nsnapshots*nsnapshots),psimatrix_recv(nsnapshots*nsnapshots);
	Array<FLT,1> eigenvalues(nmodes);

#ifdef USING_MASS_MATRIX
	create_mass_matrix(mass);
#endif
	
	
	switch (job_id) {
		/////////////////////
		// DIRECT METHOD   //
		/////////////////////
		case(direct): {
			const int sm = basis::tri(BASE::log2p)->sm();
			const int im = basis::tri(BASE::log2p)->im();
			const int ndofs = BASE::npnt +BASE::nseg*sm +BASE::ntri*im;
			int ntotal = ndofs*BASE::NV;
			Array<FLT,2> mass_times_timeaverage(ntotal,ntotal);
			Array<FLT,1> packed_mass_times_timeaverage((ntotal+1)*ntotal/2),packed_mass((ntotal+1)*ntotal/2);
			
			for (int deflation_count=0; deflation_count<ndeflation;++deflation_count) {
				int start = nmodes_per_deflation*deflation_count;

				
				/* Make packed storage version of block mass matrix for Lapack routine */
				for(int n=0;n<BASE::NV;++n) {
					for(int m=0;m<BASE::NV;++m) {
						Array<FLT,2> mass_block = mass_times_timeaverage(Range(n*ndofs,(n+1)*ndofs-1),Range(m*ndofs,(m+1)*ndofs-1));  // Reference to data, not a copy
						mass.unpack(mass_block);
					}
				}
				
				// Store block mass matrix in packed format
				psi1dcounter = 0;
				for (int k=0;k<ntotal;++k) {
					for(int l=k;l<ntotal;++l) {
						packed_mass(psi1dcounter++) = mass_times_timeaverage(k,l);
					}
				}
				
				mass_times_timeaverage = 0;
				for (int k=0;k<nsnapshots;++k) {
					nstr.str("");
					nstr << k << std::flush;
					filename = "temp" +nstr.str();
					BASE::input(filename, BASE::output_type(1)); // Loads into ug
					time_average(BASE::ug,mass_times_timeaverage);
				}
				mass_times_timeaverage /= static_cast<FLT>(nsnapshots);
				
				/* Left multiply time average blocks by mass matrix */
				Array<FLT,1> subarray(ndofs);
				for(int n=0;n<BASE::NV;++n){
					for(int i=0;i<ndofs*BASE::NV;++i) {
						// Copies data for this section of time average so it can be repeatedly multiplied
						subarray = mass_times_timeaverage(Range(n*ndofs,(n+1)*ndofs-1),i);
						// Makes reference to data so can be stored in correct location
						Array<FLT,1> mass_times_timeaverage_subarray = mass_times_timeaverage(Range(n*ndofs,(n+1)*ndofs-1),i);
						mass.mmult(subarray,mass_times_timeaverage_subarray);
					}
				}
				
				/* Right multiply time average blocks by mass matrix (done by transposing as M is symmetric) */
				mass_times_timeaverage.transposeSelf(secondDim, firstDim);
				for(int n=0;n<BASE::NV;++n){
					for(int i=0;i<ndofs*BASE::NV;++i) {
						// Copies data for this section of time average so it can be repeatedly multiplied
						subarray = mass_times_timeaverage(Range(n*ndofs,(n+1)*ndofs-1),i);
						// Makes reference to data so can be stored in correct location
						Array<FLT,1> mass_times_timeaverage_subarray = mass_times_timeaverage(Range(n*ndofs,(n+1)*ndofs-1),i);
						mass.mmult(subarray,mass_times_timeaverage_subarray);
					}
				}
				
				// Fixme:  not sure how this is going to work in parallel
				// *BASE::gbl->log << "mass_times_timeaverage" << std::endl;
				// sim::blks.allreduce(psimatrix2.data(),psimatrix_recv2.data(),nsnapshots*(nsnapshots+1)/2,blocks::flt_msg,blocks::sum,pod_id);
				
				// Store in packed format
				psi1dcounter = 0;
				for (int k=0;k<ndofs*BASE::NV;++k) {
					for(int l=k;l<ndofs*BASE::NV;++l) {
						packed_mass_times_timeaverage(psi1dcounter++) = mass_times_timeaverage(k,l);
					}
				}
				
				// generalized eigenvalue problem
				int info;
				int ITYPE = 1;
				char jobz[2] = "V", uplo[2] = "L";
				Array<FLT,1> eigenvalues_subarray(nmodes_per_deflation);
				Array<FLT,2> eigenvectors(ntotal,nmodes_per_deflation,ColumnMajorArray<2>());
				// Compute only some eigenvectors
				char range[2] = "I";
				int il=ntotal-nmodes_per_deflation+1, iu=ntotal;
				double vl=0, vu=1e18;
				int neig;
				Array<int,1> iwork(5*ntotal);
				Array<double,1> work(8*ntotal);
				Array<int,1> ifail(ntotal);
				double abstol = 2.*dlamch_("Safe minimum");
				
				dspgvx_(&ITYPE, jobz, range, uplo, &ntotal, packed_mass_times_timeaverage.data(), packed_mass.data(), &vl, &vu, &il, &iu, &abstol,
								&neig, eigenvalues_subarray.data(), eigenvectors.data(), &ntotal, work.data(), iwork.data(), ifail.data(), &info);
				
				if (info != 0) {
					*BASE::gbl->log << "Failed to find eigenmodes " << info << std::endl;
					sim::abort(__LINE__,__FILE__,BASE::gbl->log);
				}
				
				/* Reorder largest to smallest */
				if (eigenvalues_subarray(nmodes_per_deflation-1) > eigenvalues_subarray(0)) {
					eigenvalues_subarray.reverseSelf(firstDim);
					eigenvectors.reverseSelf(secondDim);
				}
				
				// Put eigenvector back into modes */
				for(int k=0;k<nmodes_per_deflation;++k)	{
					eigenvalues(start+k) = eigenvalues_subarray(k);
					int ind = 0;
					for(int n=0;n<BASE::NV;++n) {
						for (int i=0;i<BASE::npnt;++i)
							modes(start+k).v(i,n) = eigenvectors(ind++,k);
						
						for (int i=0;i<BASE::nseg;++i)
							for(int m=0;m<sm;++m)
								modes(start+k).s(i,m,n) = eigenvectors(ind++,k);
						
						for (int i=0;i<BASE::ntri;++i)
							for(int m=0;m<im;++m)
								modes(start+k).i(i,m,n) = eigenvectors(ind++,k);
					}
				}
				
				/* CALCULATE INNER PRODUCT OF MODES FOR RENORMALIZATION */
				for(int k=0;k<nmodes_per_deflation;++k) {
					psimatrix(k) = norm2(modes(start+k));
				}
				sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nmodes_per_deflation,blocks::flt_msg,blocks::sum,pod_id);
				//*BASE::gbl->log << psimatrix_recv(Range(0,nmodes_per_deflation-1)) << std::endl;
				
				/* RENORMALIZE MODES & OUTPUT */
				for (int k=0;k<nmodes_per_deflation;++k) {
					FLT norm = sqrt(psimatrix_recv(k));
					modes(start+k).v(Range(0,BASE::npnt-1)) /= norm;
					modes(start+k).s(Range(0,BASE::nseg-1)) /= norm;
					modes(start+k).i(Range(0,BASE::ntri-1)) /= norm;
					nstr.str("");
					nstr << start+k << std::flush;
					filename = "mode" +nstr.str();
					output(modes(start+k),filename);
				}
				
				/* To test normalization */
				//		psimatrix = 0.0;
				//		psimatrix_recv = 0.0;
				//		/* CALCULATE INNER PRODUCT OF MODES FOR RENORMALIZATION */
				//		for(int k=0;k<nmodes_per_deflation;++k) {
				//			psimatrix(k) = norm2(modes(start+k));
				//		}
				//		sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nmodes_per_deflation,blocks::flt_msg,blocks::sum,pod_id);
				//		*BASE::gbl->log << psimatrix_recv(Range(0,nmodes_per_deflation-1)) << std::endl;
				
				/***************************************/
				/* SUBSTRACT PROJECTION FROM SNAPSHOTS */
				/***************************************/
				if (deflation_count+1 < ndeflation) {
					psi1dcounter = 0;
					for (int k=0;k<nsnapshots;++k) {
						/* LOAD SNAPSHOT */
						nstr.str("");
						nstr << k << std::flush;
						filename = "temp" +nstr.str();
						BASE::input(filename, BASE::output_type(1));
						project_to_gauss(BASE::ug);
						
						for(int l=start;l<nmodes_per_deflation+start;++l) {
							psimatrix(psi1dcounter++) = inner_product_with_projection(modes(l));
						}
					}
					sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nsnapshots*nmodes_per_deflation,blocks::flt_msg,blocks::sum,pod_id);
					
					psi1dcounter = 0;
					for (int k=0;k<nsnapshots;++k) {
						/* LOAD SNAPSHOT */
						nstr.str("");
						nstr << k << std::flush;
						filename = "temp" +nstr.str();
						BASE::input(filename, BASE::output_type(1));
						
						for(int l=start;l<nmodes_per_deflation+start;++l) {
							BASE::ug.v(Range(0,BASE::npnt-1)) -= psimatrix_recv(psi1dcounter)*modes(l).v(Range(0,BASE::npnt-1));
							BASE::ug.s(Range(0,BASE::nseg-1)) -= psimatrix_recv(psi1dcounter)*modes(l).s(Range(0,BASE::nseg-1));
							BASE::ug.i(Range(0,BASE::ntri-1)) -= psimatrix_recv(psi1dcounter)*modes(l).i(Range(0,BASE::ntri-1));
							++psi1dcounter;
						}
						BASE::output(filename, BASE::output_type(1));
						// BASE::output(filename, BASE::output_type(0));
					}
				}
			}
			break;
		}
			
		/////////////////////
		// SNAPSHOT METHOD //
		/////////////////////
		case(snapshot): {
			for (int deflation_count=0; deflation_count<ndeflation;++deflation_count) {
				int start = nmodes_per_deflation*deflation_count;

				psi1dcounter = 0;
				for (int k=0;k<nsnapshots;++k) {
					nstr.str("");
					nstr << k << std::flush;
					filename = "temp" +nstr.str();
					BASE::input(filename, BASE::output_type(1)); // Loads into ug
					project_to_gauss(BASE::ug);
				
					for(int l=k;l<nsnapshots;++l) {
						nstr.str("");
						nstr << l << std::flush;
						filename = "temp" +nstr.str();
						BASE::input(filename, BASE::output_type(1));
						psimatrix(psi1dcounter) = inner_product_with_projection(BASE::ug);
						++psi1dcounter;
					}
				}
				// *BASE::gbl->log << "psimatrix" << psimatrix << std::endl;
				sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nsnapshots*(nsnapshots+1)/2,blocks::flt_msg,blocks::sum,pod_id);
				psimatrix_recv = psimatrix_recv/nsnapshots;
			
				int info;
				char jobz[2] = "V", uplo[2] = "L";
				Array<FLT,1> eigenvalues_subarray(nmodes_per_deflation);
				Array<FLT,2> eigenvectors(nsnapshots,nmodes_per_deflation,ColumnMajorArray<2>());
			
				if(nmodes_per_deflation == nsnapshots) {
					// Compute all eigenvectors
					Array<FLT,1> work(3*nsnapshots);
					// normal eigenvalue problem
#ifdef F2CFortran
					DSPEV(jobz,uplo,nsnapshots,psimatrix_recv.data(),eigenvalues_subarray.data(),eigenvectors.data(),nsnapshots,work.data(),info);
#else
                    dspev_(jobz,uplo,&nsnapshots,psimatrix_recv.data(),eigenvalues_subarray.data(),eigenvectors.data(),&nsnapshots,work.data(),&info);
#endif
				}
				else {
					// Compute only some eigenvectors
					char range[2] = "I";
					int il=nsnapshots-nmodes_per_deflation+1, iu=nsnapshots;
					double vl, vu;
					int neig;
					Array<int,1> iwork(5*nsnapshots);
					Array<double,1> work(8*nsnapshots);
					Array<int,1> ifail(nsnapshots);
					double abstol = 2.*dlamch_("Safe minimum");
				
					/* To compute needed sizes
					 Array<int,1> isuppz(2*nsnapshots);
					 int lwork = 30*nsnapshots;
					 Array<FLT,1> work(lwork);
					 Array<int,1> iwork(lwork);
					 dsyevr_(jobz, range, uplo, &nsnapshots, psimatrix_recv.data(), &nsnapshots, &vl, &vu, &il, &iu, &abstol,
					 &neig, eigenvalues.data(),eigenvectors.data(), &nsnapshots, isuppz.data(), work.data(), &lwork,&iwork, iwork.data(), &info);
					 std::cout << "optimal sizes:" <<  work(0) << ' ' << iwork(0) << std::endl;
					 */
					 
					// normal eigenvalue problem
					dspevx_(jobz, range, uplo, &nsnapshots, psimatrix_recv.data(),&vl, &vu, &il, &iu, &abstol,
									&neig, eigenvalues_subarray.data(), eigenvectors.data(), &nsnapshots, work.data(), iwork.data(), ifail.data(), &info);
				}
			
				if (info != 0) {
					*BASE::gbl->log << "Failed to find eigenmodes " << info << std::endl;
					sim::abort(__LINE__,__FILE__,BASE::gbl->log);
				}
				/* Reorder largest to smallest */
				if (eigenvalues_subarray(nmodes_per_deflation-1) > eigenvalues_subarray(0)) {
					eigenvalues_subarray.reverseSelf(firstDim);
					eigenvectors.reverseSelf(secondDim);
				}
			
				/* construct POD MODES */
				for(int k=0;k<nmodes_per_deflation;++k)	{
					eigenvalues(start+k) = eigenvalues_subarray(k);
					modes(start+k).v(Range(0,BASE::npnt-1)) = 0.0;
					modes(start+k).s(Range(0,BASE::nseg-1)) = 0.0;
					modes(start+k).i(Range(0,BASE::ntri-1)) = 0.0;
				
					/* LOAD SNAPSHOTS AND CALCULATE MODE */
					for(int l=0;l<nsnapshots;++l)	{
						nstr.str("");
						nstr << l << std::flush;
						filename = "temp" +nstr.str();
						BASE::input(filename, BASE::output_type(1));
						modes(start+k).v(Range(0,BASE::npnt-1)) += eigenvectors(l,k)*BASE::ug.v(Range(0,BASE::npnt-1));
						modes(start+k).s(Range(0,BASE::nseg-1)) += eigenvectors(l,k)*BASE::ug.s(Range(0,BASE::nseg-1));
						modes(start+k).i(Range(0,BASE::ntri-1)) += eigenvectors(l,k)*BASE::ug.i(Range(0,BASE::ntri-1));
					}
				}
			
				/* CALCULATE INNER PRODUCT OF MODES FOR RENORMALIZATION */
				for(int k=0;k<nmodes_per_deflation;++k) {
					psimatrix(k) = norm2(modes(start+k));
				}
				sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nmodes_per_deflation,blocks::flt_msg,blocks::sum,pod_id);
				//*BASE::gbl->log << psimatrix_recv(Range(0,nmodes_per_deflation-1)) << std::endl;
			
				/* RENORMALIZE MODES & OUTPUT */
				for (int k=0;k<nmodes_per_deflation;++k) {
					FLT norm = sqrt(psimatrix_recv(k));
					modes(start+k).v(Range(0,BASE::npnt-1)) /= norm;
					modes(start+k).s(Range(0,BASE::nseg-1)) /= norm;
					modes(start+k).i(Range(0,BASE::ntri-1)) /= norm;
					nstr.str("");
					nstr << start+k << std::flush;
					filename = "mode" +nstr.str();
					output(modes(start+k),filename);
				}
			
				/* To test normalization */
				//		psimatrix = 0.0;
				//		psimatrix_recv = 0.0;
				//		/* CALCULATE INNER PRODUCT OF MODES FOR RENORMALIZATION */
				//		for(int k=0;k<nmodes_per_deflation;++k) {
				//			psimatrix(k) = norm2(modes(start+k));
				//		}
				//		sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nmodes_per_deflation,blocks::flt_msg,blocks::sum,pod_id);
				//		*BASE::gbl->log << psimatrix_recv(Range(0,nmodes_per_deflation-1)) << std::endl;
			
				/***************************************/
				/* SUBSTRACT PROJECTION FROM SNAPSHOTS */
				/***************************************/
				if (deflation_count+1 < ndeflation) {
					psi1dcounter = 0;
					for (int k=0;k<nsnapshots;++k) {
						/* LOAD SNAPSHOT */
						nstr.str("");
						nstr << k << std::flush;
						filename = "temp" +nstr.str();
						BASE::input(filename, BASE::output_type(1));
						project_to_gauss(BASE::ug);
					
						for(int l=start;l<nmodes_per_deflation+start;++l) {
							psimatrix(psi1dcounter++) = inner_product_with_projection(modes(l));
						}
					}
					sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nsnapshots*nmodes_per_deflation,blocks::flt_msg,blocks::sum,pod_id);
				
					psi1dcounter = 0;
					for (int k=0;k<nsnapshots;++k) {
						/* LOAD SNAPSHOT */
						nstr.str("");
						nstr << k << std::flush;
						filename = "temp" +nstr.str();
						BASE::input(filename, BASE::output_type(1));
					
						for(int l=start;l<nmodes_per_deflation+start;++l) {
							BASE::ug.v(Range(0,BASE::npnt-1)) -= psimatrix_recv(psi1dcounter)*modes(l).v(Range(0,BASE::npnt-1));
							BASE::ug.s(Range(0,BASE::nseg-1)) -= psimatrix_recv(psi1dcounter)*modes(l).s(Range(0,BASE::nseg-1));
							BASE::ug.i(Range(0,BASE::ntri-1)) -= psimatrix_recv(psi1dcounter)*modes(l).i(Range(0,BASE::ntri-1));
							++psi1dcounter;
						}
						BASE::output(filename, BASE::output_type(1));
						// BASE::output(filename, BASE::output_type(0));
					}
				}
			}
			break;
		}
		case(mPOD): {
			/******************************************************************************/
			/**************************  mPOD         *************************************/
			/******************************************************************************/
			std::string write_name = "temp";
			
			/* For first level eigenvalues are just 1 */
			eigenvalues.resize(nsnapshots);
			eigenvalues = 1./nsnapshots;
			
			int npartitions = ceil(nsnapshots/M);  //  Calculate initial number of partitions
			int nlevel = 0;
			
			while (true) { //  do while loop implemented with a break statement at end
				*BASE::gbl->log << "level = " << nlevel << ", npartitions = " << npartitions << ", M = " << M << ", Lam_min = " << MinEigen << std::endl;

				if (npartitions == 1) write_name = "mode";
				
				 // counter of number of modes for next level
				int next_nsnapshots = 0;
						
				//  Loop over each partition (this loop could be parallelized...)
				for (int np=0;np<npartitions;++np) { 		

					// index of beginning of set by eq. 15
					int I0 = np*M;
					// index of end of set avoiding integer errors
					int I1 = min((np+1)*M,nsnapshots);
					// M for this set, again avoiding integer errors
					int Mp = I1-I0;					
						
					psimatrix.resize(Mp*(Mp+1)/2);  //  this is the snapshot matrix (see eq. 17)
					psimatrix_recv.resize(Mp*(Mp+1)/2);  // % this is the snapshot matrix (see eq. 17) for mpi communications
					psimatrix = 0.0;
					psimatrix_recv = 0.0;
					psi1dcounter = 0;
					for (int k=I0;k<I1;++k) {
						nstr.str("");
						nstr << k << std::flush;
						filename = "temp" +nstr.str();
						BASE::input(filename, BASE::reload_type); 
						project_to_gauss(BASE::ug);
						
						for(int l=k;l<I1;++l) {
							nstr.str("");
							nstr << l << std::flush;
							filename = "temp"+nstr.str();
							BASE::input(filename, BASE::reload_type);
							psimatrix(psi1dcounter) = inner_product_with_projection(BASE::ug);
							++psi1dcounter;
						}
					}
					sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),Mp*(Mp+1)/2,blocks::flt_msg,blocks::sum,pod_id);
					
					psimatrix=0.;
					for (int i = 0;i<Mp;++i) {
						psimatrix(i + i*(2*Mp-i-1)/2) = eigenvalues(I0+i);
					}
					
					// Compute all eigenvectors
					int info;
					char jobz[2] = "V", uplo[2] = "L";
					int ITYPE=3;
					Array<FLT,1> work;
					work.resize(3*Mp);
					Array<FLT,1> eigenvalues_p(Mp);
					Array<FLT,2> eigenvectors_p(Mp,Mp,ColumnMajorArray<2>());
					dspgv_(&ITYPE, jobz, uplo, &Mp, psimatrix_recv.data(), psimatrix.data(), eigenvalues_p.data(), eigenvectors_p.data(), &Mp, work.data(), &info);
					if (info != 0) {
						*BASE::gbl->log << "Failed to find eigenmodes " << info << std::endl;
						sim::abort(__LINE__,__FILE__,BASE::gbl->log);
					}

					/* Reorder largest to smallest */
					if (eigenvalues_p(Mp-1) > eigenvalues_p(0)) {
						eigenvalues_p.reverseSelf(firstDim);
						eigenvectors_p.reverseSelf(secondDim);
					}
					*BASE::gbl->log << "eigenvalues of partition "<<  np << " in decomposition level " << nlevel <<  " are "<<  eigenvalues_p << std::endl;
					int modes_p = first(eigenvalues_p < MinEigen);
					if (modes_p < 0) modes_p = Mp;
					
					*BASE::gbl->log << "the index of first eigenvalue below "<<  MinEigen << " in the partition " <<  np << " is "<<  modes_p << std::endl;

					if (modes_p > 0) {
						eigenvalues(Range(next_nsnapshots,next_nsnapshots+modes_p-1)) = eigenvalues_p(Range(0,modes_p-1));
				
						/* construct POD MODES */
						for(int k=0;k<modes_p;++k) {
							modes(k).v(Range(0,BASE::npnt-1)) = 0.0;
							modes(k).s(Range(0,BASE::nseg-1)) = 0.0;
							modes(k).i(Range(0,BASE::ntri-1)) = 0.0;
						
							/* LOAD SNAPSHOTS AND CALCULATE MODE */
							for(int l=0;l<Mp;++l) {
								nstr.str("");
								nstr << I0 +l << std::flush;
								filename = "temp" +nstr.str();
								BASE::input(filename, BASE::reload_type);
								modes(k).v(Range(0,BASE::npnt-1)) += eigenvectors_p(l,k)*BASE::ug.v(Range(0,BASE::npnt-1));
								modes(k).s(Range(0,BASE::nseg-1)) += eigenvectors_p(l,k)*BASE::ug.s(Range(0,BASE::nseg-1));
								modes(k).i(Range(0,BASE::ntri-1)) += eigenvectors_p(l,k)*BASE::ug.i(Range(0,BASE::ntri-1));
							}
						}
		
						/* CALCULATE INNER PRODUCT OF MODES FOR RENORMALIZATION */
						for(int k=0;k<modes_p;++k) {
							psimatrix(k) = norm2(modes(k));
						}
						sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),modes_p,blocks::flt_msg,blocks::sum,pod_id);
				
						/* RENORMALIZE MODES & OUTPUT */
						for (int k=0;k<modes_p;++k) {
							FLT norm = sqrt(psimatrix_recv(k));
							modes(k).v(Range(0,BASE::npnt-1)) /= norm;
							modes(k).s(Range(0,BASE::nseg-1)) /= norm;
							modes(k).i(Range(0,BASE::ntri-1)) /= norm;
							nstr.str("");
							nstr << k +next_nsnapshots << std::flush;
							filename = write_name +nstr.str();
							output(modes(k),filename);
						}
						next_nsnapshots += modes_p;
					}
				}
				nsnapshots = next_nsnapshots;
				nmodes = nsnapshots;

				if (npartitions == 1)  {
					break;
				}
					
				++nlevel;
				npartitions = ceil(min(npartitions/2,nsnapshots/M));
				M = ceil(nsnapshots/npartitions);
			}
		}
		break;
	}

	/* OUTPUT EIGENVALUES VECTOR */
	*BASE::gbl->log << "eigenvalues: "<<  eigenvalues(Range(0,nmodes-1)) << std::endl;
	filename = "eigenvalues_" +BASE::gbl->idprefix;
	BASE::output(nmodes,eigenvalues,filename,BASE::output_type(1));
	BASE::output(nmodes,eigenvalues,filename,BASE::output_type(0));

	/* CALCULATE POD COEFFICIENTS FOR EXPANSION OF SNAPSHOTS */
	const int ncoefficients = (coefficient_end-coefficient_start)/coefficient_interval*nmodes;
	psimatrix.resize(ncoefficients);
	psimatrix_recv.resize(ncoefficients);
	psimatrix = 0.0;
	psimatrix_recv = 0.0;
	psi1dcounter=0;
	for (int k=coefficient_start;k<coefficient_end;k+=coefficient_interval) {
		/* LOAD SNAPSHOT */
		nstr.str("");
		nstr << k << std::flush;
		filename = "rstrt" +nstr.str() +"_d0";
		BASE::input(filename, BASE::reload_type);
		project_to_gauss(BASE::ug);
		for(int l=0;l<nmodes;++l) {
			psimatrix(psi1dcounter) = inner_product_with_projection(modes(l));
			++psi1dcounter;
		}
	}
	sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),ncoefficients,blocks::flt_msg,blocks::sum,pod_id);

	psi1dcounter = 0;
	for (int k=coefficient_start;k<coefficient_end;k+=coefficient_interval) {
		/* OUTPUT COEFFICIENT VECTOR */
		nstr.str("");
		nstr << k << std::flush;
		filename = "coeff" +nstr.str() + "_" +BASE::gbl->idprefix;
		BASE::output(nmodes,psimatrix_recv(Range(psi1dcounter,psi1dcounter+nmodes-1)),filename,BASE::output_type(1));
		BASE::output(nmodes,psimatrix_recv(Range(psi1dcounter,psi1dcounter+nmodes-1)),filename,BASE::output_type(0));
		psi1dcounter += nmodes;
	}

#ifdef POD_BDRY
	/* NOW GENERATE BDRY POD MODES */
	for (int i=0;i<BASE::nebd;++i)
		pod_ebdry(i)->calculate_modes();
#endif
	
	sim::finalize(__LINE__,__FILE__,BASE::gbl->log);
	
	return;
}


template<class BASE> void pod_generate<BASE>::create_mass_matrix(sparse_row_major& mass) {
	
	/* Calculate Storage for Sparse Mass Matrix */
	const int sm=basis::tri(BASE::log2p)->sm();
	const int im=basis::tri(BASE::log2p)->im();
	const int tm=basis::tri(BASE::log2p)->tm();
	const int lgpx = basis::tri(BASE::log2p)->gpx(), lgpn = basis::tri(BASE::log2p)->gpn();
	
	/* rank of mass matrix */
	const int size = BASE::npnt +BASE::nseg*sm +BASE::ntri*im;
	
	/* find number of non-zeros for each row */
	Array<int,1> nnzero(size);
	const int begin_seg = BASE::npnt;
	const int begin_tri = begin_seg+BASE::nseg*sm;
	const int ndofs = begin_tri +BASE::ntri*im;
	
	/* SELF CONNECTIONS */
	nnzero(Range(0,begin_seg-1)) = 1;
	if (sm) {
		nnzero(Range(begin_seg,begin_tri-1)) = (2 +sm);
	}
	if (im) nnzero(Range(begin_tri,ndofs-1)) = tm;
	
	/* edges and vertices connected to a vertex */
	for(int i=0; i<BASE::npnt; ++i)
		nnzero(i) += BASE::pnt(i).nnbor*(sm+1);
	
	for(int i=0; i<BASE::ntri; ++i) {
		/* interior and opposing side mode for each vertex */
		for(int j=0;j<3;++j)
			nnzero(BASE::tri(i).pnt(j)) += (im+sm);
		
		
		/* interior modes,opposing side modes, and opposing vertex to each side */
		for(int j=0;j<3;++j)
			for(int m=0;m<sm;++m)
				nnzero(begin_seg+BASE::tri(i).seg(j)*sm+m) += im +2*sm +1;
	}
	mass.resize(size,nnzero);
	
	/* Calculate Sparse Mass Matrix */
	BASE::uht(0) = 0.0;
	
	for(int tind=0;tind<BASE::ntri;++tind) {
		
		/* Create vector of global indices */
		Array<int,1> gindx(tm),gsign(tm);
		gsign=1;
		
		/* VERTEX MODES */
		int indx = 0;
		for (int m=0; m<3; ++m) {
			gindx(indx++) = BASE::tri(tind).pnt(m);
		}
		
		if (sm > 0) {
			/* SIDE MODES */
			for(int i=0;i<3;++i) {
				int sind = BASE::npnt +BASE::tri(tind).seg(i)*sm;
				int sgn = BASE::tri(tind).sgn(i);
				int msgn = 1;
				for (int m = 0; m < sm; ++m) {
					gindx(indx) = sind;
					gsign(indx++) = msgn;
					msgn *= sgn;
					++sind;
				}
			}
			
			int iind = BASE::npnt +BASE::nseg*sm +tind*im;
			for(int m=0;m<im;++m) {
				gindx(indx++) = iind++;
			}
		}
		
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		BASE::crdtocht(tind);
		
		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(int n=0;n<BASE::ND;++n)
			basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);
		
		for(int m=0;m<tm;++m) {
			BASE::uht(0)(m) = 1.0*gsign(m);
			basis::tri(BASE::log2p)->proj(&BASE::uht(0)(0),&BASE::u(0)(0,0),MXGP);
			
			
			for(int i=0;i<lgpx;++i) {
				for(int j=0;j<lgpn;++j) {
					BASE::u(0)(i,j) *= RAD(BASE::crd(0)(i,j))*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
				}
			}
			basis::tri(BASE::log2p)->intgrt(&BASE::lf(0)(0),&BASE::u(0)(0,0),MXGP);
			BASE::uht(0)(m) = 0.0;
			
			/* store in mass matrix */
			for(int k=0;k<tm;++k)
				mass.add_values(gindx(m),gindx(k),gsign(k)*BASE::lf(0)(k));
		}
	}
	mass.check_for_unused_entries(*BASE::gbl->log);
	
	return;
}


template<class BASE> void pod_generate<BASE>::test_orthogonality() {
	/* CALCULATE POD MODES ORTHOGONALITY */
	Array<FLT,1> psimatrix(nmodes*nmodes);
	Array<FLT,1> psimatrix_recv(nmodes*nmodes);
	
	psimatrix = 0.0;
	psimatrix_recv = 0.0;
	int psi1dcounter = 0 ;
	
	std::string filename;
	std::ostringstream nstr;
	
	for (int k=0;k<nmodes;++k) {
		/* LOAD MODE k */
		nstr.str("");
		nstr << k << std::flush;
		filename = "mode" +nstr.str();
		BASE::input(filename, BASE::output_type(1));
		project_to_gauss(BASE::ug);
		
		for(int l=0;l<nmodes;++l) {
			/* LOAD MODE l */
			nstr.str("");
			nstr << l << std::flush;
			filename = "mode" +nstr.str();
			BASE::input(filename, BASE::output_type(1));
			psimatrix(psi1dcounter) = inner_product_with_projection(BASE::ug);
			++psi1dcounter;
		}
	}
	sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nmodes*nmodes,blocks::flt_msg,blocks::sum,pod_id);
	
	/* Output Orthogonality Matrix */
	filename = "orthogonality_" +BASE::gbl->idprefix +".txt";
	ofstream out1;
	out1.open(filename.c_str());
	out1.precision(16);
	
	for (int k=0;k<nmodes;++k) {
		for (int l=0;l<nmodes;++l) {
			out1 << k << "   " << l << "   " << psimatrix_recv(k*nmodes +l) << std::endl;
		}
	}
	out1.close();
	
	return;
}

#ifndef USING_MASS_MATRIX
template<class BASE> void pod_generate<BASE>::project_to_gauss(vsi& target) {
	const int lgpx = basis::tri(BASE::log2p)->gpx(), lgpn = basis::tri(BASE::log2p)->gpn();
	
	
	BASE::ug.v.reference(target.v);
	BASE::ug.s.reference(target.s);
	BASE::ug.i.reference(target.i);
	
	for(int tind=0;tind<BASE::ntri;++tind) {
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		BASE::crdtocht(tind);
		
		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(int n=0;n<BASE::ND;++n)
			basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);
		
		/* PROJECT MODE k TO GAUSS POINTS */
		BASE::ugtouht(tind);
		for(int n=0;n<BASE::NV;++n)
			basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);
		
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {
				FLT cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p)->wtx(i)*basis::tri(BASE::log2p)->wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
				for(int n=0;n<BASE::NV;++n) {
					BASE::dugdt(BASE::log2p)(tind,n,i,j) = BASE::u(n)(i,j)*scaling(n)*cjcb;
				}
			}
		}
	}
	/* The default is for ug to refer to ugbd(0) */
	BASE::ug.v.reference(BASE::ugbd(0).v);
	BASE::ug.s.reference(BASE::ugbd(0).s);
	BASE::ug.i.reference(BASE::ugbd(0).i);
	
	return;
}

template<class BASE> FLT pod_generate<BASE>::inner_product_with_projection(vsi& target) {
	const int lgpx = basis::tri(BASE::log2p)->gpx(), lgpn = basis::tri(BASE::log2p)->gpn();
	FLT returnval = 0.0;
	
#ifdef LOW_NOISE_DOT
	Array<FLT,1> low_noise_dot;
	low_noise_dot.resize(BASE::ntri);
#endif
	
	BASE::ug.v.reference(target.v);
	BASE::ug.s.reference(target.s);
	BASE::ug.i.reference(target.i);
	
	
	for(int tind=0;tind<BASE::ntri;++tind) {
		/* PROJECT MODE l TO GAUSS POINTS */
		BASE::ugtouht(tind);
		for(int n=0;n<BASE::NV;++n)
			basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);
		
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {
				for(int n=0;n<BASE::NV;++n) {
					returnval += BASE::dugdt(BASE::log2p)(tind,n,i,j)*BASE::u(n)(i,j);
				}
			}
		}
#ifdef LOW_NOISE_DOT
		low_noise_dot(tind) = returnval;
		returnval = 0.0;
#endif
	}
	
#ifdef LOW_NOISE_DOT
	/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
	int halfcount,remainder;
	for (remainder=BASE::ntri % 2, halfcount = BASE::ntri/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
		for (int tind=0;tind<halfcount;++tind)
			low_noise_dot(tind) += low_noise_dot(tind+halfcount);
		if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
	}
	returnval = low_noise_dot(0);
#endif
	
	/* The default is for ug to refer to ugbd(0) */
	BASE::ug.v.reference(BASE::ugbd(0).v);
	BASE::ug.s.reference(BASE::ugbd(0).s);
	BASE::ug.i.reference(BASE::ugbd(0).i);
	
	return(returnval);
}
#else
template<class BASE> void pod_generate<BASE>::project_to_gauss(vsi& target) {
	const int sm=basis::tri(BASE::log2p)->sm();
	const int im=basis::tri(BASE::log2p)->im();
	const int ndofs = BASE::npnt +BASE::nseg*sm +BASE::ntri*im;
	
	
	/* Make 1D Vector */
	Array<FLT,1> snapshot(ndofs);
	for(int n=0;n<BASE::NV;++n) {
		int ind = 0;
		for (int i=0;i<BASE::npnt;++i)
			snapshot(ind++) = target.v(i,n);
		
		for (int i=0;i<BASE::nseg;++i)
			for(int m=0;m<sm;++m)
				snapshot(ind++) = target.s(i,m,n);
		
		for (int i=0;i<BASE::ntri;++i)
			for(int m=0;m<im;++m)
				snapshot(ind++) = target.i(i,m,n);
		
		Array<FLT,1> subarray = mass_times_snapshot(Range(n*ndofs,(n+1)*ndofs-1));
		mass.mmult(snapshot,subarray);
		subarray *= scaling(n);
	}
	return;
}

template<class BASE> void pod_generate<BASE>::time_average(vsi& target, Array<FLT,2> average_matrix) {
	const int sm=basis::tri(BASE::log2p)->sm();
	const int im=basis::tri(BASE::log2p)->im();
	const int ndofs = BASE::npnt +BASE::nseg*sm +BASE::ntri*im;
	
	/* Make 1D Vector */
	Array<FLT,1> snapshot(ndofs*BASE::NV);
	int ind = 0;
	for(int n=0;n<BASE::NV;++n) {
		for (int i=0;i<BASE::npnt;++i)
			snapshot(ind++) = target.v(i,n);
		
		for (int i=0;i<BASE::nseg;++i)
			for(int m=0;m<sm;++m)
				snapshot(ind++) = target.s(i,m,n);
		
		for (int i=0;i<BASE::ntri;++i)
			for(int m=0;m<im;++m)
				snapshot(ind++) = target.i(i,m,n);
	}
	
	for (int i=0;i<ndofs*BASE::NV;++i){
		for (int j=0;j<ndofs*BASE::NV;++j){
			average_matrix(i,j) += snapshot(i)*snapshot(j);
		}
		
	}
	return;
}

template<class BASE> FLT pod_generate<BASE>::inner_product_with_projection(vsi& target) {
	const int sm=basis::tri(BASE::log2p)->sm();
	const int im=basis::tri(BASE::log2p)->im();
	const int ndofs = BASE::npnt +BASE::nseg*sm +BASE::ntri*im;
	
	
	Array<FLT,1> snapshot(ndofs);
	
	FLT returnval = 0.0;
	
	for(int n=0;n<BASE::NV;++n) {
		int ind = 0;
		for (int i=0;i<BASE::npnt;++i)
			snapshot(ind++) = target.v(i,n);
		
		for (int i=0;i<BASE::nseg;++i)
			for(int m=0;m<sm;++m)
				snapshot(ind++) = target.s(i,m,n);
		
		for (int i=0;i<BASE::ntri;++i)
			for(int m=0;m<im;++m)
				snapshot(ind++) = target.i(i,m,n);
		
		returnval += dot(snapshot,mass_times_snapshot(Range(n*ndofs,(n+1)*ndofs-1)));
	}
	return(returnval);
}
#endif

template<class BASE> FLT pod_generate<BASE>::norm2(vsi& target) {
	const int lgpx = basis::tri(BASE::log2p)->gpx(), lgpn = basis::tri(BASE::log2p)->gpn();
	FLT returnval = 0.0;
	
	BASE::ug.v.reference(target.v);
	BASE::ug.s.reference(target.s);
	BASE::ug.i.reference(target.i);
	
	for(int tind=0;tind<BASE::ntri;++tind) {
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		BASE::crdtocht(tind);
		
		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(int n=0;n<BASE::ND;++n)
			basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);
		
		/* PROJECT MODE TO GAUSS POINTS */
		BASE::ugtouht(tind);
		for(int n=0;n<BASE::NV;++n)
			basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);
		
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {
				FLT cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p)->wtx(i)*basis::tri(BASE::log2p)->wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
				for(int n=0;n<BASE::NV;++n) {
					returnval += BASE::u(n)(i,j)*BASE::u(n)(i,j)*scaling(n)*cjcb;
				}
			}
		}
	}
	
	/* The default is for ug to refer to ugbd(0) */
	BASE::ug.v.reference(BASE::ugbd(0).v);
	BASE::ug.s.reference(BASE::ugbd(0).s);
	BASE::ug.i.reference(BASE::ugbd(0).i);
	
	return(returnval);
	
}

template<class BASE> void pod_generate<BASE>::output(vsi& target,std::string filename) {
	
	/* Output refers directly to ugbd rather than using ug */
	vsi ugstore;
	ugstore.v.reference(BASE::ugbd(0).v);
	ugstore.s.reference(BASE::ugbd(0).s);
	ugstore.i.reference(BASE::ugbd(0).i);
	
	BASE::ugbd(0).v.reference(target.v);
	BASE::ugbd(0).s.reference(target.s);
	BASE::ugbd(0).i.reference(target.i);
	
	BASE::output(filename, BASE::output_type(0));
	BASE::output(filename, BASE::output_type(1));
	
	BASE::ugbd(0).v.reference(ugstore.v);
	BASE::ugbd(0).s.reference(ugstore.s);
	BASE::ugbd(0).i.reference(ugstore.i);
	
	return;
}

#ifdef POD_BDRY
/***********************/
/* POD_BDRY FUNCTION   */
/***********************/
template<class BASE> void pod_gen_edge_bdry<BASE>::init(input_map& inmap) {
	std::string keyword;
	
	keyword = base.idprefix + "_pod";
	inmap.getwdefault(keyword,active,false);
	
	if (!active) return;
	
	keyword = base.idprefix + "_podmodes";
	inmap.getwdefault(keyword,nmodes,x.nsnapshots);
	
	keyword = base.idprefix + "_pod_id";
	if (!inmap.get(keyword,pod_id)) {
		*x.gbl->log << "Must provide a pod id for pod boundary" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	/* ERROR CHECK THAT NUMBER IS LISTED IN GROUP LIST */
	if (!inmap.getline(x.gbl->idprefix +"_groups",keyword)) {
		*x.gbl->log << "group list must contain pod_id for " << base.idprefix << "" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	istringstream mystr(keyword);
	int bnum;
	do {
		if (!(mystr >> bnum)) {
			*x.gbl->log << "group list must contain pod_id for " << base.idprefix << "" << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	} while (bnum != pod_id);
}

template<class BASE> void pod_gen_edge_bdry<BASE>::zero_bdry(tri_hp::vsi& ug) {
	
	if (!active) return;
	
	int sind,v0;
	
	int j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		ug.v(v0,Range(0,x.NV-1)) = 0.0;
		ug.s(sind,Range(0,x.sm0-1),Range(0,x.NV-1)) = 0.0;
	} while(++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	ug.v(v0,Range(0,x.NV-1)) = 0.0;
}

/* 1D VERSION TO GENERATE MODES */
template<class BASE> void pod_gen_edge_bdry<BASE>::calculate_modes() {
	
	if (!active) return;
	
	Array<FLT,1> psimatrix(x.nsnapshots*x.nsnapshots);
	Array<FLT,1> psimatrix_recv(x.nsnapshots*x.nsnapshots);
	Array<FLT,1> low_noise_dot(x.BASE::ntri);
	Array<FLT,1> eigenvalues(x.nsnapshots);
	Array<FLT,1> eigenvector(x.nsnapshots);
	Array<FLT,2> coeff(x.nsnapshots,x.nsnapshots);
	char jobz[2] = "V", range[2] = "I", uplo[2] = "L";
	double vl, vu;
	int il=x.nsnapshots, iu=x.nsnapshots;
	double abstol = 2.*dlamch_("Safe minimum");
	int sind,v0,neig,info;
	Array<int,1> iwork(5*x.nsnapshots);
	Array<double,1> work(8*x.nsnapshots);
	Array<int,1> ifail(x.nsnapshots);
	int i,k,l,n;
	int lgpx = basis::tri(x.log2p)->gpx();
	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	FLT cjcb;
	
	/* MAKE FILES FOR EDGE MODE SNAPSHOT-PROJECTION */
	for(k=0;k<x.nsnapshots;++k) {
		nstr.str("");
		nstr << k*x.restart_interval +x.restartfile << std::flush;
		filename = "rstrt" +nstr.str() +"_d0";
		x.input(filename, x.BASE::reload_type);
		nstr.str("");
		nstr << k << std::flush;
		filename = "temp" +nstr.str();
		x.BASE::output(filename, x.BASE::output_type(1));
	}
	
	
	for (int eig_ct=0; eig_ct<x.nsnapshots;++eig_ct) {
		
		/* ******************************************/
		/* GENERATE POD MODES SNAPSHOT MATRIX       */
		/********************************************/
		int psi1dcounter = 0;
		psimatrix = 0.0;
		for (k=0;k<x.nsnapshots;++k) {
			nstr.str("");
			nstr << k << std::flush;
			filename = "temp" +nstr.str();
			x.input(filename, x.BASE::output_type(1));
			
			for(l=k;l<x.nsnapshots;++l) {
				nstr.str("");
				nstr << l << std::flush;
				filename = "temp" +nstr.str();
				x.input(filename, x.BASE::output_type(1), 1);
				
				/* PERFORM 1D INTEGRATION */
				int bsind = 0;
				do {
					sind = base.seg(bsind);
					
					/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
					x.crdtocht1d(sind);
					
					/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
					for(n=0;n<x.ND;++n)
						basis::tri(x.log2p)->proj1d(&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0));
					
					x.ugtouht1d(sind,0);
					for(n=0;n<x.NV;++n)
						basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.u(n)(0,0));
					
					x.ugtouht1d(sind,1);
					for(n=0;n<x.NV;++n)
						basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.res(n)(0,0));
					
					FLT tmp_store = 0.0;
					for(i=0;i<lgpx;++i) {
						cjcb =  basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
						for(n=0;n<x.NV;++n) {
							tmp_store += x.u(n)(0,i)*x.res(n)(0,i)*x.scaling(n)*cjcb;
						}
					}
					low_noise_dot(bsind) = tmp_store;
				} while(++bsind < base.nseg);
				
				/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
				int halfcount,remainder;
				for (remainder=base.nseg % 2, halfcount = base.nseg/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
					for (int bsind=0;bsind<halfcount;++bsind)
						low_noise_dot(bsind) += low_noise_dot(bsind+halfcount);
					if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
				}
				psimatrix(psi1dcounter) = low_noise_dot(0);
				
				++psi1dcounter;
			}
		}
		sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),x.nsnapshots*(x.nsnapshots+1)/2,blocks::flt_msg,blocks::sum,pod_id);
		
		/*********************************/
		/* TO COMPUTE JUST 1 EIGENVECTOR */
		/*********************************/
		dspevx_(jobz, range, uplo, &x.nsnapshots, psimatrix_recv.data(), &vl, &vu, &il, &iu, &abstol,
						&neig, eigenvalues.data(),eigenvector.data(), &x.nsnapshots, work.data(), iwork.data(), ifail.data(), &info);
		if (info != 0) {
			*x.gbl->log << "Failed to find eigenmodes " << info << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		
		*x.gbl->log << base.idprefix << " eigenvalue "<<  eig_ct << ' ' << eigenvalues(0) << std::endl;
		
		
		/*************************************/
		/* LOAD SNAPSHOTS AND CALCULATE MODE */
		/*************************************/
		/* FOR NOW I AM GOING TO CALCULATE 2D MODE */
		/* EVEN THOUGH THIS IS UNNECESSARY */
		x.ug.v = 0.0;
		x.ug.s = 0.0;
		x.ug.i = 0.0;
		for(l=0;l<x.nsnapshots;++l)	{
			nstr.str("");
			nstr << l << std::flush;
			filename = "temp" +nstr.str();
			x.input(filename, x.output_type(1), 1);
			
			int bsind = 0;
			do {
				sind = base.seg(bsind);
				v0 = x.seg(sind).pnt(0);
				x.ug.v(v0) += eigenvector(l)*x.ugbd(1).v(v0);
				x.ug.s(sind) += eigenvector(l)*x.ugbd(1).s(sind);
			} while(++bsind < base.nseg);
			v0 = x.seg(sind).pnt(1);
			x.ug.v(v0) += eigenvector(l)*x.ugbd(1).v(v0);
		}
		
		
		/********************/
		/* RENORMALIZE MODE */
		/********************/
		for(int bsind=0;bsind<base.nseg;++bsind) {
			sind = base.seg(bsind);
			
			/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
			x.crdtocht1d(sind);
			
			/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(n=0;n<x.ND;++n)
				basis::tri(x.log2p)->proj1d(&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0));
			
			x.ugtouht1d(sind);
			for(n=0;n<x.NV;++n)
				basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.res(n)(0,0));
			
			FLT tmp_store = 0.0;
			for(i=0;i<lgpx;++i) {
				cjcb =  basis::tri(x.log2p)->wtx(i)*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
				for(n=0;n<x.NV;++n) {
					tmp_store += x.res(n)(0,i)*x.res(n)(0,i)*x.scaling(n)*cjcb;
				}
			}
			low_noise_dot(bsind) = tmp_store;
		}
		/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
		int halfcount,remainder;
		for (remainder=base.nseg % 2, halfcount = base.nseg/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
			for (int bsind=0;bsind<halfcount;++bsind)
				low_noise_dot(bsind) += low_noise_dot(bsind+halfcount);
			if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
		}
		double norm = low_noise_dot(0);
		double norm_recv;
		
		sim::blks.allreduce(&norm,&norm_recv,1,blocks::flt_msg,blocks::sum,pod_id);
		
		/* DIVIDE BY NORM */
		norm = sqrt(norm_recv);
		x.ug.v(Range(0,x.npnt-1)) /= norm;
		x.ug.s(Range(0,x.nseg-1)) /= norm;
		
		/* 2D OUTPUT RENORMALIZED MODE */
		nstr.str("");
		nstr << eig_ct << std::flush;
		filename = "mode" +nstr.str();
		x.BASE::output(filename, x.BASE::output_type(0));
		
		/* 1D OUTPUT */
		output(x.ug,filename,x.BASE::output_type(0));
		output(x.ug,filename,x.BASE::output_type(1));
		
		/***************************************/
		/* SUBSTRACT PROJECTION FROM SNAPSHOTS */
		/***************************************/
		double dotp_recv, dotp;
		for (k=0;k<x.nsnapshots;++k) {
			/* LOAD SNAPSHOT */
			nstr.str("");
			nstr << k << std::flush;
			filename = "temp" +nstr.str();
			x.input(filename, x.output_type(1), 1);
			
			/* PERFORM 1D INTEGRATION */
			for(int bsind=0;bsind<base.nseg;++bsind) {
				sind = base.seg(bsind);
				
				/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
				x.crdtocht1d(sind);
				
				/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
				for(n=0;n<x.ND;++n)
					basis::tri(x.log2p)->proj1d(&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0));
				
				x.ugtouht1d(sind,0);
				for(n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.u(n)(0,0));
				
				x.ugtouht1d(sind,1);
				for(n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.res(n)(0,0));
				
				FLT tmp_store = 0.0;
				for(i=0;i<lgpx;++i) {
					cjcb =  basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
					for(n=0;n<x.NV;++n) {
						tmp_store += x.u(n)(0,i)*x.res(n)(0,i)*x.scaling(n)*cjcb;
					}
				}
				low_noise_dot(bsind) = tmp_store;
			}
			
			/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
			int halfcount,remainder;
			for (remainder=base.nseg % 2, halfcount = base.nseg/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
				for (int bsind=0;bsind<halfcount;++bsind)
					low_noise_dot(bsind) += low_noise_dot(bsind+halfcount);
				if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
			}
			dotp = low_noise_dot(0);
			
			sim::blks.allreduce(&dotp,&dotp_recv,1,blocks::flt_msg,blocks::sum,pod_id);
			x.ugbd(1).v(Range(0,x.npnt-1)) -= dotp_recv*x.ug.v(Range(0,x.npnt-1));
			x.ugbd(1).s(Range(0,x.nseg-1)) -= dotp_recv*x.ug.s(Range(0,x.nseg-1));
			x.BASE::output(filename, x.BASE::output_type(1), 1);
			coeff(k,eig_ct) = dotp_recv;
		}
	}
	
	/* Output coefficient vectors */
	for (k=0;k<x.nsnapshots;++k) {
		nstr.str("");
		nstr << k*x.restart_interval +x.restartfile << std::flush;
		filename = "coeff" +nstr.str() + "_" +base.idprefix;
		x.BASE::output(nmodes,coeff(k,Range::all()),filename,x.output_type(1));
		x.BASE::output(nmodes,coeff(k,Range::all()),filename,x.output_type(0));
	}
}

template<class BASE> void pod_gen_edge_bdry<BASE>::output(vsi& target,std::string filename, typename BASE::filetype typ) {
	
	switch(typ) {
		case(tri_hp::text): {
			std::string fname;
			fname = filename +"_" +base.idprefix +".txt";
			ofstream fout;
			fout.open(fname.c_str(),std::ofstream::out | std::ofstream::app);
			
			int bsind = 0;
			int sind,v0;
			do {
				sind = base.seg(bsind);
				v0 = x.seg(sind).pnt(0);
				for (int n=0;n<x.NV;++n)
					fout << x.ug.v(v0,n) << ' ';
				fout << std::endl;
			} while(++bsind < base.nseg);
			v0 = x.seg(sind).pnt(1);
			for (int n=0;n<x.NV;++n)
				fout << x.ug.v(v0,n) << ' ';
			fout << std::endl;
			
			for (int bsind=0;bsind<base.nseg;++bsind) {
				sind = base.seg(bsind);
				for (int m=0;m<x.sm0;++m) {
					for (int n=0;n<x.NV;++n)
						fout << x.ug.s(sind,m,n) << ' ';
					fout << std::endl;
				}
			}
			fout.close();
			
			break;
		}
			
#ifdef libbinio
		case(tri_hp::binary): {
			/* 1D OUTPUT RENORMALIZED MODE */
			std::string fname = filename + "_" + base.idprefix +".bin";
			binofstream bout;
			bout.open(fname.c_str());
			if (bout.error()) {
				*x.gbl->log << "couldn't open coefficient output file " << filename;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
			bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);
			
			int bsind = 0;
			int sind,v0;
			do {
				sind = base.seg(bsind);
				v0 = x.seg(sind).pnt(0);
				for (int n=0;n<x.NV;++n)
					bout.writeFloat(x.ug.v(v0,n),binio::Double);
			} while(++bsind < base.nseg);
			v0 = x.seg(sind).pnt(1);
			for (int n=0;n<x.NV;++n)
				bout.writeFloat(x.ug.v(v0,n),binio::Double);
			
			for (int bsind=0;bsind<base.nseg;++bsind) {
				sind = base.seg(bsind);
				for (int m=0;m<x.sm0;++m)
					for (int n=0;n<x.NV;++n)
						bout.writeFloat(x.ug.s(sind,m,n),binio::Double);
			}
			bout.close();
			break;
		}
#endif
			
		case(tri_hp::netcdf): {
			std::string fname = filename + "_" + base.idprefix +".nc";
			
			/* Create the file. The NC_CLOBBER parameter tells netCDF to
			 * overwrite this file, if it already exists.*/
			int retval, ncid, dims[3];
			int one, two, three;
			if ((retval = nc_create(fname.c_str(), NC_CLOBBER|NC_NETCDF4, &ncid))) ERR(retval,x.gbl->log);
			
			/* some fixed dimensions */
			if ((retval = nc_def_dim(ncid,"1",1,&one))) ERR(retval,x.gbl->log);
			if ((retval = nc_def_dim(ncid,"2",2,&two))) ERR(retval,x.gbl->log);
			if ((retval = nc_def_dim(ncid,"3",3,&three))) ERR(retval,x.gbl->log);
			
			/* Define the dimensions. NetCDF will hand back an ID for each. */
			int npnt_id,nseg_id,sm0_id,NV_id;
			if ((retval = nc_def_dim(ncid, "npnt", base.nseg+1, &npnt_id))) ERR(retval,x.gbl->log);
			if ((retval = nc_def_dim(ncid, "nseg", base.nseg, &nseg_id))) ERR(retval,x.gbl->log);
			if ((retval = nc_def_dim(ncid, "sm0", x.sm0, &sm0_id))) ERR(retval,x.gbl->log);
			if ((retval = nc_def_dim(ncid, "NV", x.NV, &NV_id))) ERR(retval,x.gbl->log);
			
			int ugv_id;
			dims[0] = npnt_id;
			dims[1] = NV_id;
			if ((retval = nc_def_var(ncid, "ugv", NC_DOUBLE, 2, dims, &ugv_id))) ERR(retval,x.gbl->log);
			
			int ugs_id;
			dims[0] = nseg_id;
			dims[1] = sm0_id;
			dims[2] = NV_id;
			if ((retval = nc_def_var(ncid, "ugs", NC_DOUBLE, 3, dims, &ugs_id))) ERR(retval,x.gbl->log);
			
			if ((retval = nc_enddef(ncid))) ERR(retval,x.gbl->log);

			size_t index[3];
			int bsind = 0;
			int sind,v0;
			do {
				sind = base.seg(bsind);
				index[0] = bsind;
				v0 = x.seg(sind).pnt(0);
				for (int n=0;n<x.NV;++n) {
					index[1] = n;
					nc_put_var1_double(ncid,ugv_id,index,&x.ug.v(v0,n));
				}
			} while(++bsind < base.nseg);
			index[0] = bsind;
			v0 = x.seg(sind).pnt(1);
			for (int n=0;n<x.NV;++n) {
				index[1] = n;
				nc_put_var1_double(ncid,ugv_id,index,&x.ug.v(v0,n));
			}
	
			for(int i=0;i<base.nseg;++i) {
				index[0] = i;
				for(int m=0;m<x.sm0;++m) {
					index[1] = m;
					for(int n=0;n<x.NV;++n) {
						index[2] = n;
						nc_put_var1_double(ncid,ugs_id,index,&x.ug.s(base.seg(i),m,n));
					}
				}
			}
			if ((retval = nc_close(ncid))) ERR(retval,x.gbl->log);
			
			break;
		}
			
		case(tri_hp::tecplot): {
			std::string fname = filename + "_" + base.idprefix +".dat";
			std::ofstream fout;
			fout.open(fname.c_str());
			
			fout << "VARIABLES=\"S\",\"X\",\"Y\",";
			for(int n=0;n<x.NV;++n)
				fout << "\"V" << n << "\",";
			fout << "\nTITLE = " << base.idprefix << '\n'<< "ZONE\n";
			
			FLT circumference = 0.0;
			int ind = 0;
			do {
				int sind = base.seg(ind);
				int tind = x.seg(sind).tri(0);
				
				int seg;
				for(seg=0;seg<3;++seg)
					if (x.tri(tind).seg(seg) == sind) break;
				assert(seg != 3);
				
				x.crdtocht(tind);
				for(int m=basis::tri(x.log2p)->bm();m<basis::tri(x.log2p)->tm();++m)
					for(int n=0;n<tri_mesh::ND;++n)
						x.cht(n,m) = 0.0;
				
				for(int n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj_side(seg,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));
				
				x.ugtouht(tind);
				for(int n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj_side(seg,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));
				
				for (int i=0;i<basis::tri(x.log2p)->gpx();++i) {
					FLT arclength = sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
					FLT jcb =  basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*arclength;
					circumference += jcb;
					
					fout << circumference << ' ' << x.crd(0)(0,i) << ' ' << x.crd(1)(0,i) << ' ';
					
					for (int n=0;n<x.NV;++n) {
						/* Output value, tangent, and normal derivatives */
						fout << x.u(n)(0,i) << ' ';
						
					}
					fout << std::endl;
				}
			} while (++ind < base.nseg);
			fout.close();
			
			break;
		}
			
		default:
			break;
	}
	return;
}

#endif

