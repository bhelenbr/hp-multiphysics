/*
 *  pod_generate.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 1/18/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#include "pod_generate.h"

#define LOW_NOISE_DOT

#include <myblas.h>
#include <libbinio/binfile.h>
// #include <veclib/clapack.h>

extern "C" {
	double dlamch_(const char *cmach);
	/* Subroutine */ int dspevx_(char *jobz, char *range, char *uplo, int *n, 
							 double *ap, double *vl, double *vu, int *il, int *
							 iu, double *abstol, int *m, double *w, double *z__, 
							 int *ldz, double *work, int *iwork, int *ifail, 
							 int *info);
}	


template<class BASE> void pod_generate<BASE>::init(input_map& input, void *gin) {
	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	std::istringstream instr;
	int i;

	/* Initialize base class */
	BASE::init(input,gin);

	if (!input.get(BASE::gbl->idprefix + "_snapshots",nsnapshots)) {
		input.getwdefault("snapshots",nsnapshots,10);
	}

	if (!input.get(BASE::gbl->idprefix + "_nmodes",nmodes)) input.getwdefault("nmodes",nmodes,nsnapshots);    

	/* THIS IS TO CHANGE THE WAY SNAPSHOT MATRIX ENTRIES ARE FORMED */
	scaling.resize(BASE::NV);
	scaling = 1;
	if (input.getline(BASE::gbl->idprefix + "_scale_vector",linebuff) || input.getline("scale_vector",linebuff)) {
		instr.str(linebuff);
		for(i=0;i<BASE::NV;++i)
			instr >> scaling(i);
	}

	nmodes = MAX(nmodes,2);


#ifndef LOWNOISE
	modes.resize(nmodes);
	for(i=0;i<nmodes;++i) {
		modes(i).v.resize(BASE::maxpst,BASE::NV);
		modes(i).s.resize(BASE::maxpst,BASE::sm0,BASE::NV);
		modes(i).i.resize(BASE::maxpst,BASE::im0,BASE::NV);
	}
#else
	input.getwdefault(BASE::gbl->idprefix + "_groups",pod_id,0);

#ifdef POD_BDRY
	pod_ebdry.resize(BASE::nebd);
	for (int i=0;i<BASE::nebd;++i) {
		pod_ebdry(i) = new pod_gen_edge_bdry<BASE>(*this,*BASE::ebdry(i));
		pod_ebdry(i)->init(input);
	}
#endif
#endif

	input.getwdefault("restart",restartfile,1);

	return;
}

#ifndef LOWNOISE

template<class BASE> void pod_generate<BASE>::tadvance() {
	int i,j,k,l,n,tind,info;
	int lgpx = basis::tri(BASE::log2p)->gpx(), lgpn = basis::tri(BASE::log2p)->gpn();
	Array<FLT,1> low_noise_dot;
	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	std::istringstream instr;
	FLT cjcb;

	if (BASE::gbl->substep != 0) return;

	BASE::tadvance(); 

	int psi1dcounter = 0;
	vsi ugstore;
	ugstore.v.reference(BASE::ugbd(0).v);
	ugstore.s.reference(BASE::ugbd(0).s);
	ugstore.i.reference(BASE::ugbd(0).i);

	psimatrix.resize(nsnapshots*nsnapshots);
	psimatrix_recv.resize(nsnapshots*nsnapshots);
	low_noise_dot.resize(BASE::ntri);

	/* GENERATE POD MODES SNAPSHOT COEFFICIENTS */
	psimatrix = 0.0;
	for (k=0;k<nsnapshots;++k) {
		nstr.str("");
		nstr << k+restartfile << std::flush;
		filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
		BASE::ugbd(0).v.reference(modes(0).v);
		BASE::ugbd(0).s.reference(modes(0).s);
		BASE::ugbd(0).i.reference(modes(0).i);
		BASE::input(filename, BASE::binary);

		for(l=k;l<nsnapshots;++l) {
			nstr.str("");
			nstr << l+restartfile << std::flush;
			filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
			BASE::ugbd(0).v.reference(modes(1).v);
			BASE::ugbd(0).s.reference(modes(1).s);
			BASE::ugbd(0).i.reference(modes(1).i);
			BASE::input(filename, BASE::binary);

			for(tind=0;tind<BASE::ntri;++tind) {          
				/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
				BASE::crdtocht(tind);

				/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
				for(n=0;n<BASE::ND;++n)
					basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);

				BASE::ugbd(0).v.reference(modes(0).v);
				BASE::ugbd(0).s.reference(modes(0).s);
				BASE::ugbd(0).i.reference(modes(0).i);
				BASE::ugtouht(tind,0);
				for(n=0;n<BASE::NV;++n)
					basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);

				BASE::ugbd(0).v.reference(modes(1).v);
				BASE::ugbd(0).s.reference(modes(1).s);
				BASE::ugbd(0).i.reference(modes(1).i);
				BASE::ugtouht(tind,0);
				for(n=0;n<BASE::NV;++n)
					basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);

				FLT tmp_store = 0.0;
				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpn;++j) {
						cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p)->wtx(i)*basis::tri(BASE::log2p)->wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
						for(n=0;n<BASE::NV;++n) {
						tmp_store += BASE::u(n)(i,j)*BASE::res(n)(i,j)*scaling(n)*cjcb;
						}
					}
				}
				psimatrix(psi1dcounter) += tmp_store;
			}
			++psi1dcounter;
		}	
		BASE::ugbd(0).v.reference(ugstore.v);
		BASE::ugbd(0).s.reference(ugstore.s);
		BASE::ugbd(0).i.reference(ugstore.i);
	}
	sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nsnapshots*(nsnapshots+1)/2,blocks::flt_msg,blocks::sum);

	Array<FLT,1> eigenvalues(nsnapshots);
	Array<FLT,2> eigenvectors(nsnapshots,nsnapshots);

	/* To compute needed sizes 
	 char jobz[2] = "V", range[2] = "A", uplo[2] = "L";
	 double vl, vu;
	 int il, iu;
	 double abstol = dlamch_("Safe minimum");
	 int neig;
	 Array<int,1> isuppz(2*nsnapshots);
	 int lwork = 30*nsnapshots;
	 lwork = -1;
	 Array<FLT,1> work(lwork);
	 Array<int,1> iwork(lwork);
	 int info1, info;

	 info1 = dsyevr_(jobz, range, uplo, &nsnapshots, psimatrix_recv.data(), &nsnapshots, &vl, &vu, &il, &iu, &abstol,
	 &neig, eigenvalues.data(),eigenvectors.data(), &nsnapshots, isuppz.data(), work.data(), &lwork,&iwork, iwork.data(), &info); 
	 std::cout << "optimal sizes:" <<  work(0) << ' ' << iwork(0) << std::endl; 
	 */

	/* TO COMPUTE JUST 1 EIGENVECTOR
	 char jobz[2] = "V", range[2] = "I"; uplo[2] = "L";
	 double vl, vu;
	 int il=nsnapshots, iu=nsnapshots;
	 double abstol = 2.*dlamch_("Safe minimum");
	 int neig,info;
	 Array<int,1> iwork(5*nsnapshots);
	 Array<double,1> work(8*nsnapshots);
	 Array<int,1> ifail(nsnapshots);

	 info = dspevx_(jobz, range, uplo, &nsnapshots, psimatrix_recv.data(), &nsnapshots, &vl, &vu, &il, &iu, &abstol,
	 &neig, eigenvalues.data(),eigenvectors.data(), &nsnapshots, isuppz.data(), work.data(), &lwork,&iwork, iwork.data(), &info); 
	 */

	char jobz[2] = "V", uplo[2] = "L";
	Array<FLT,1> work(3*nsnapshots);
	DSPEV(jobz,uplo,nsnapshots,psimatrix_recv.data(),eigenvalues.data(),eigenvectors.data(),nsnapshots,work.data(),info);

	if (info != 0) {
		*BASE::gbl->log << "Failed to find eigenmodes " << info << std::endl;
		sim::abort(__LINE__,__FILE__,BASE::gbl->log);
	}

	*BASE::gbl->log << "eigenvalues "<<  eigenvalues << std::endl;

	eigenvectors.transposeSelf(secondDim,firstDim);  // FORTRAN ROUTINE RETURNS TRANSPOSE

	// construct POD MODES
	for(k=0;k<nmodes;++k)	{
		modes(k).v(Range(0,BASE::npnt-1)) = 0.0;
		modes(k).s(Range(0,BASE::nseg-1)) = 0.0;
		modes(k).i(Range(0,BASE::ntri-1)) = 0.0;

		/* LOAD SNAPSHOTS AND CALCULATE MODE */
		for(l=0;l<nsnapshots;++l)	{
			nstr.str("");
			nstr << l +restartfile << std::flush;
			filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
			BASE::input(filename, BASE::binary);

			modes(k).v(Range(0,BASE::npnt-1)) += eigenvectors(l,nmodes -k -1)*BASE::ug.v(Range(0,BASE::npnt-1));
			modes(k).s(Range(0,BASE::nseg-1)) += eigenvectors(l,nmodes -k -1)*BASE::ug.s(Range(0,BASE::nseg-1));
			modes(k).i(Range(0,BASE::ntri-1)) += eigenvectors(l,nmodes -k -1)*BASE::ug.i(Range(0,BASE::ntri-1));
		}
	}

	ugstore.v.reference(BASE::ugbd(1).v);
	ugstore.s.reference(BASE::ugbd(1).s);
	ugstore.i.reference(BASE::ugbd(1).i);

	psimatrix = 0.0;
	psimatrix_recv = 0.0;

	/* CALCULATE INNER PRODUCT OF MODES FOR RENORMALIZATION */
	for(l=0;l<nmodes;++l) {
		BASE::ugbd(1).v.reference(modes(l).v);
		BASE::ugbd(1).s.reference(modes(l).s);
		BASE::ugbd(1).i.reference(modes(l).i);

		for(tind=0;tind<BASE::ntri;++tind) {          
			/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
			BASE::crdtocht(tind);

			/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(n=0;n<BASE::ND;++n)
				basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);

			/* PROJECT MODE TO GAUSS POINTS */
			BASE::ugtouht(tind,1);
			for(n=0;n<BASE::NV;++n)
				basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);

			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p)->wtx(i)*basis::tri(BASE::log2p)->wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
					for(n=0;n<BASE::NV;++n) {
						psimatrix(l) += BASE::res(n)(i,j)*BASE::res(n)(i,j)*scaling(n)*cjcb;
					}
				}
			}
		}
	}	
	BASE::ugbd(1).v.reference(ugstore.v);
	BASE::ugbd(1).s.reference(ugstore.s);
	BASE::ugbd(1).i.reference(ugstore.i);
	sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nmodes,blocks::flt_msg,blocks::sum);

	ugstore.v.reference(BASE::ugbd(0).v);
	ugstore.s.reference(BASE::ugbd(0).s);
	ugstore.i.reference(BASE::ugbd(0).i);

	FLT norm;
	/* RENORMALIZE MODES AND OUTPUT COEFFICIENTS */
	for (k=0;k<nmodes;++k) {
		norm = sqrt(psimatrix_recv(k));
		modes(k).v(Range(0,BASE::npnt-1)) /= norm;
		modes(k).s(Range(0,BASE::nseg-1)) /= norm;
		modes(k).i(Range(0,BASE::ntri-1)) /= norm;

		nstr.str("");
		nstr << k << std::flush;
		filename = "mode" +nstr.str() + "_" + BASE::gbl->idprefix;
		BASE::ugbd(0).v.reference(modes(k).v);
		BASE::ugbd(0).s.reference(modes(k).s);
		BASE::ugbd(0).i.reference(modes(k).i);
		BASE::output(filename, BASE::binary);
		BASE::output(filename, BASE::tecplot);
	}
	BASE::ugbd(0).v.reference(ugstore.v);
	BASE::ugbd(0).s.reference(ugstore.s);
	BASE::ugbd(0).i.reference(ugstore.i);


	ugstore.v.reference(BASE::ugbd(1).v);
	ugstore.s.reference(BASE::ugbd(1).s);
	ugstore.i.reference(BASE::ugbd(1).i);

	/* CALCULATE POD COEFFICIENTS FOR EXPANSION OF SNAPSHOTS */
	psimatrix = 0.0;
	psimatrix_recv = 0.0;
	psi1dcounter=0;
	for (k=0;k<nsnapshots;++k) {
		/* LOAD SNAPSHOT */
		nstr.str("");
		nstr << k+restartfile << std::flush;
		filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
		BASE::input(filename, BASE::binary);

		for(l=0;l<nmodes;++l) {
			BASE::ugbd(1).v.reference(modes(l).v);
			BASE::ugbd(1).s.reference(modes(l).s);
			BASE::ugbd(1).i.reference(modes(l).i);

			for(tind=0;tind<BASE::ntri;++tind) {          
				/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
				BASE::crdtocht(tind);

				/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
				for(n=0;n<BASE::ND;++n)
					basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);

				/* PROJECT SNAPSHOT TO GAUSS POINTS */
				BASE::ugtouht(tind);
				for(n=0;n<BASE::NV;++n)
					basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);

				/* PROJECT MODE TO GAUSS POINTS */
				BASE::ugtouht(tind,1);
				for(n=0;n<BASE::NV;++n)
					basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);

				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpn;++j) {
						cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p)->wtx(i)*basis::tri(BASE::log2p)->wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
						for(n=0;n<BASE::NV;++n) {
							psimatrix(psi1dcounter) += BASE::u(n)(i,j)*BASE::res(n)(i,j)*scaling(n)*cjcb;
						}
					}
				}
			}
			++psi1dcounter;
		}
	}
	BASE::ugbd(1).v.reference(ugstore.v);
	BASE::ugbd(1).s.reference(ugstore.s);
	BASE::ugbd(1).i.reference(ugstore.i);

	sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nsnapshots*nmodes,blocks::flt_msg,blocks::sum);

	for (k=0;k<nsnapshots;++k) {
		/* OUTPUT COEFFICIENT VECTOR */
		nstr.str("");
		nstr << k+restartfile << std::flush;
		filename = "coeff" +nstr.str() + "_" +BASE::gbl->idprefix +".bin";
		binofstream bout;
		bout.open(filename.c_str());
		if (bout.error()) {
			*BASE::gbl->log << "couldn't open coefficient output file " << filename;
			sim::abort(__LINE__,__FILE__,BASE::gbl->log);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (l=0;l<nmodes;++l) 
			bout.writeFloat(psimatrix_recv(k*nmodes +l),binio::Double);

		bout.close();
	}
	
	sim::finalize(__LINE__,__FILE__,BASE::gbl->log);
		
	return;
}

#else	
template<class BASE> void pod_generate<BASE>::tadvance() {

	if (BASE::gbl->substep != 0) return;

	BASE::tadvance(); 

	Array<FLT,1> psimatrix(nsnapshots*nsnapshots);
	Array<FLT,1> psimatrix_recv(nsnapshots*nsnapshots);
	Array<FLT,1> low_noise_dot(BASE::ntri);
	Array<FLT,1> eigenvalues(nsnapshots);
	Array<FLT,1> eigenvector(nsnapshots);
	Array<FLT,2> coeff(nsnapshots,nsnapshots);
	char jobz[2] = "V", range[2] = "I", uplo[2] = "L";
	double vl, vu;
	int il=nsnapshots, iu=nsnapshots;
	double abstol = 2.*dlamch_("Safe minimum");
	int neig,info;
	Array<int,1> iwork(5*nsnapshots);
	Array<double,1> work(8*nsnapshots);
	Array<int,1> ifail(nsnapshots);
	int i,j,k,l,n,tind;
	int lgpx = basis::tri(BASE::log2p)->gpx(), lgpn = basis::tri(BASE::log2p)->gpn();
	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	std::istringstream instr;
	FLT cjcb;

	/* MAKE FILES FOR VOLUME MODE SNAPSHOT-PROJECTION */
	for(k=restartfile;k<nsnapshots+restartfile;++k) {
		nstr.str("");
		nstr << k << std::flush;
		filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
		BASE::input(filename, BASE::binary);
		nstr.str("");
		nstr << k-restartfile << std::flush;
		filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;

#ifdef POD_BDRY
		/* ZERO SNAPSHOTS ON POD BOUNDARY'S */
		for (i=0;i<BASE::nebd;++i)
			pod_ebdry(i)->zero_bdry(BASE::ug);
#endif

		BASE::output(filename, BASE::binary);
	}


	for (int eig_ct=0; eig_ct<nmodes;++eig_ct) {

		/* ******************************************/
		/* GENERATE POD MODES SNAPSHOT MATRIX       */
		/********************************************/
		int psi1dcounter = 0;
		psimatrix = 0.0;
		for (k=0;k<nsnapshots;++k) {
			nstr.str("");
			nstr << k << std::flush;
			filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;
			BASE::input(filename, BASE::binary); // Loads into ug/ugbd(0)

			for(l=k;l<nsnapshots;++l) {
				nstr.str("");
				nstr << l << std::flush;
				filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;
				BASE::input(filename, BASE::binary, 1); // Loads into ug/ugbd(1)

				for(tind=0;tind<BASE::ntri;++tind) {          
					/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
					BASE::crdtocht(tind);

					/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
					for(n=0;n<BASE::ND;++n)
						basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);

					BASE::ugtouht(tind,0);
					for(n=0;n<BASE::NV;++n)
						basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);

					BASE::ugtouht(tind,1);
					for(n=0;n<BASE::NV;++n)
						basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);

					FLT tmp_store = 0.0;
					for(i=0;i<lgpx;++i) {
						for(j=0;j<lgpn;++j) {
							cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p)->wtx(i)*basis::tri(BASE::log2p)->wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
							for(n=0;n<BASE::NV;++n) {
								tmp_store += BASE::u(n)(i,j)*BASE::res(n)(i,j)*scaling(n)*cjcb;
							}
						}
					}
					low_noise_dot(tind) = tmp_store;
#ifndef LOW_NOISE_DOT
					psimatrix(psi1dcounter) += tmp_store;
#endif
				}
#ifdef LOW_NOISE_DOT
				/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
				int halfcount,remainder;
				for (remainder=BASE::ntri % 2, halfcount = BASE::ntri/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
					for (tind=0;tind<halfcount;++tind) 
						low_noise_dot(tind) += low_noise_dot(tind+halfcount);
					if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
				}
				psimatrix(psi1dcounter) = low_noise_dot(0);
#endif
				++psi1dcounter;
			}	
		}
		sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nsnapshots*(nsnapshots+1)/2,blocks::flt_msg,blocks::sum,pod_id);

		/*********************************/
		/* TO COMPUTE JUST 1 EIGENVECTOR */
		/*********************************/
		dspevx_(jobz, range, uplo, &nsnapshots, psimatrix_recv.data(), &vl, &vu, &il, &iu, &abstol,
				&neig, eigenvalues.data(),eigenvector.data(), &nsnapshots, work.data(), iwork.data(), ifail.data(), &info); 
		if (info != 0) {
			*BASE::gbl->log << "Failed to find eigenmodes " << info << std::endl;
			sim::abort(__LINE__,__FILE__,BASE::gbl->log);
		}

		*BASE::gbl->log << "eigenvalue "<<  eig_ct << ' ' << eigenvalues(0) << std::endl;


		/*************************************/
		/* LOAD SNAPSHOTS AND CALCULATE MODE */
		/*************************************/
		BASE::ug.v = 0.0;
		BASE::ug.s = 0.0;
		BASE::ug.i = 0.0;
		for(l=0;l<nsnapshots;++l)	{
			nstr.str("");
			nstr << l << std::flush;
			filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;
			BASE::input(filename, BASE::binary, 1);

			BASE::ug.v(Range(0,BASE::npnt-1)) += eigenvector(l)*BASE::ugbd(1).v(Range(0,BASE::npnt-1));
			BASE::ug.s(Range(0,BASE::nseg-1)) += eigenvector(l)*BASE::ugbd(1).s(Range(0,BASE::nseg-1));
			BASE::ug.i(Range(0,BASE::ntri-1)) += eigenvector(l)*BASE::ugbd(1).i(Range(0,BASE::ntri-1));
		}


		/********************/
		/* RENORMALIZE MODE */
		/********************/
		double norm = 0.0;
		for(tind=0;tind<BASE::ntri;++tind) {          
			/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
			BASE::crdtocht(tind);

			/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(n=0;n<BASE::ND;++n)
				basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);

			/* PROJECT MODE TO GAUSS POINTS */
			BASE::ugtouht(tind);
			for(n=0;n<BASE::NV;++n)
				basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);

			FLT tmp_store = 0.0;
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p)->wtx(i)*basis::tri(BASE::log2p)->wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
					for(n=0;n<BASE::NV;++n) {
						tmp_store += BASE::res(n)(i,j)*BASE::res(n)(i,j)*scaling(n)*cjcb;
					}
				}
			}
			low_noise_dot(tind) = tmp_store;
#ifndef LOW_NOISE_DOT
			norm += tmp_store;
#endif
		}
#ifdef LOW_NOISE_DOT
		/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
		int halfcount,remainder;
		for (remainder=BASE::ntri % 2, halfcount = BASE::ntri/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
			for (tind=0;tind<halfcount;++tind) 
				low_noise_dot(tind) += low_noise_dot(tind+halfcount);
			if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
		}
		norm = low_noise_dot(0);
#endif
		double norm_recv;

		sim::blks.allreduce(&norm,&norm_recv,1,blocks::flt_msg,blocks::sum,pod_id);

		/* DIVIDE BY NORM */
		norm = sqrt(norm_recv);
		BASE::ug.v(Range(0,BASE::npnt-1)) /= norm;
		BASE::ug.s(Range(0,BASE::nseg-1)) /= norm;
		BASE::ug.i(Range(0,BASE::ntri-1)) /= norm;

		/* OUTPUT RENORMALIZED MODE */
		nstr.str("");
		nstr << eig_ct << std::flush;
		filename = "mode" +nstr.str() + "_" + BASE::gbl->idprefix;
		BASE::output(filename, BASE::binary);
		BASE::output(filename, BASE::tecplot);

		/***************************************/
		/* SUBSTRACT PROJECTION FROM SNAPSHOTS */
		/***************************************/
		double dotp_recv, dotp;
		for (k=0;k<nsnapshots;++k) {
			/* LOAD SNAPSHOT */
			nstr.str("");
			nstr << k << std::flush;
			filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;
			BASE::input(filename, BASE::binary, 1);
			dotp = 0.0;

			for(tind=0;tind<BASE::ntri;++tind) {          
				/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
				BASE::crdtocht(tind);

				/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
				for(n=0;n<BASE::ND;++n)
					basis::tri(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);

				/* PROJECT SNAPSHOT TO GAUSS POINTS */
				BASE::ugtouht(tind);
				for(n=0;n<BASE::NV;++n)
					basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);

				/* PROJECT MODE TO GAUSS POINTS */
				BASE::ugtouht(tind,1);
				for(n=0;n<BASE::NV;++n)
					basis::tri(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);


				FLT tmp_store = 0.0;
				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpn;++j) {
						cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p)->wtx(i)*basis::tri(BASE::log2p)->wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
						for(n=0;n<BASE::NV;++n) {
						tmp_store += BASE::u(n)(i,j)*BASE::res(n)(i,j)*scaling(n)*cjcb;
						}
					}
				}
				low_noise_dot(tind) = tmp_store;
#ifndef LOW_NOISE_DOT
				dotp += tmp_store;
#endif
			}
#ifdef LOW_NOISE_DOT
			/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
			int halfcount,remainder;
			for (remainder=BASE::ntri % 2, halfcount = BASE::ntri/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
				for (tind=0;tind<halfcount;++tind) 
					low_noise_dot(tind) += low_noise_dot(tind+halfcount);
				if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
			}
			dotp = low_noise_dot(0);
#endif
			sim::blks.allreduce(&dotp,&dotp_recv,1,blocks::flt_msg,blocks::sum,pod_id);

			BASE::ugbd(1).v(Range(0,BASE::npnt-1)) -= dotp_recv*BASE::ug.v(Range(0,BASE::npnt-1));
			BASE::ugbd(1).s(Range(0,BASE::nseg-1)) -= dotp_recv*BASE::ug.s(Range(0,BASE::nseg-1));
			BASE::ugbd(1).i(Range(0,BASE::ntri-1)) -= dotp_recv*BASE::ug.i(Range(0,BASE::ntri-1));
			BASE::output(filename, BASE::binary, 1);
			coeff(k,eig_ct) = dotp_recv;
		}
	}

	/*****************************/
	/* OUTPUT COEFFICIENT VECTOR */
	/*****************************/
	for (k=restartfile;k<restartfile+nsnapshots;++k) {
		nstr.str("");
		nstr << k << std::flush;
		filename = "coeff" +nstr.str() + "_" +BASE::gbl->idprefix +".bin";
		binofstream bout;
		bout.open(filename.c_str());
		if (bout.error()) {
			*BASE::gbl->log << "couldn't open coefficient output file " << filename;
			sim::abort(__LINE__,__FILE__,BASE::gbl->log);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (l=0;l<nmodes;++l) 
			bout.writeFloat(coeff(k-restartfile,l),binio::Double);

		bout.close();
	}

#ifdef POD_BDRY
	/* NOW GENERATE BDRY POD MODES */
	for (int i=0;i<BASE::nebd;++i)
		pod_ebdry(i)->calculate_modes();
#endif

#ifdef petsc
	PetscFinalize();
#endif
#ifdef PTH
	pth_exit(NULL);
#endif    
#ifdef MPISRC
	MPI_Finalize();
#endif
	sim::abort(__LINE__,__FILE__,BASE::gbl->log);
}
#endif

#ifdef POD_BDRY
template<class BASE> void pod_gen_edge_bdry<BASE>::init(input_map& input) {
	std::string keyword;

	keyword = base.idprefix + "_pod";
	input.getwdefault(keyword,active,false);

	if (!active) return;

	keyword = base.idprefix + "_podmodes";
	input.getwdefault(keyword,nmodes,x.nsnapshots);

	keyword = base.idprefix + "_pod_id";
	if (!input.get(keyword,pod_id)) {
		*x.gbl->log << "Must provide a pod id for pod boundary" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}

	/* ERROR CHECK THAT NUMBER IS LISTED IN GROUP LIST */
	if (!input.getline(x.gbl->idprefix +"_groups",keyword)) {
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

template<class BASE> void pod_gen_edge_bdry<BASE>::zero_bdry(tri_hp::vsi ug) {

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
	std::istringstream instr;
	FLT cjcb;

	/* MAKE FILES FOR EDGE MODE SNAPSHOT-PROJECTION */
	for(k=0;k<x.nsnapshots;++k) {
		nstr.str("");
		nstr << k+x.restartfile << std::flush;
		filename = "rstrt" +nstr.str() + "_" + x.gbl->idprefix +".d0";
		x.input(filename, BASE::binary);
		nstr.str("");
		nstr << k << std::flush;
		filename = "temp" +nstr.str() + "_" + base.idprefix;
		x.output(filename, BASE::binary);
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
			filename = "temp" +nstr.str() + "_" + base.idprefix;
			x.input(filename, x.binary);

			for(l=k;l<x.nsnapshots;++l) {
				nstr.str("");
				nstr << l << std::flush;
				filename = "temp" +nstr.str() + "_" + base.idprefix;
				x.input(filename, x.binary, 1);

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

		*x.gbl->log << "eigenvalue "<<  eig_ct << ' ' << eigenvalues(0) << std::endl;


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
			filename = "temp" +nstr.str() + "_" + base.idprefix;
			x.input(filename, x.binary, 1);

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
		filename = "mode" +nstr.str() + "_" + base.idprefix;
		x.output(filename, x.tecplot);

		/* 1D OUTPUT RENORMALIZED MODE */
		filename = "mode" +nstr.str() + "_" + base.idprefix +".bin";
		binofstream bout;
		bout.open(filename.c_str());
		if (bout.error()) {
			*x.gbl->log << "couldn't open coefficient output file " << filename;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		int bsind = 0;
		do {
			sind = base.seg(bsind);
			v0 = x.seg(sind).pnt(0);
			for (n=0;n<x.NV;++n)
				bout.writeFloat(x.ug.v(v0,n),binio::Double);
		} while(++bsind < base.nseg);
		v0 = x.seg(sind).pnt(1);
		for (n=0;n<x.NV;++n)
			bout.writeFloat(x.ug.v(v0,n),binio::Double);		

		for (int bsind=0;bsind<base.nseg;++bsind) {
			sind = base.seg(bsind);
			for (int m=0;m<x.sm0;++m)
				for (n=0;n<x.NV;++n)
					bout.writeFloat(x.ug.s(sind,m,n),binio::Double);
		}
		bout.close();


		/***************************************/
		/* SUBSTRACT PROJECTION FROM SNAPSHOTS */
		/***************************************/
		double dotp_recv, dotp;
		for (k=0;k<x.nsnapshots;++k) {
			/* LOAD SNAPSHOT */
			nstr.str("");
			nstr << k << std::flush;
			filename = "temp" +nstr.str() + "_" + base.idprefix;
			x.input(filename, x.binary, 1);

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
			x.output(filename, x.binary, 1);
			coeff(k,eig_ct) = dotp_recv;
		}
	}

	/*****************************/
	/* OUTPUT COEFFICIENT VECTOR */
	/*****************************/
	for (k=0;k<x.nsnapshots;++k) {
		nstr.str("");
		nstr << k+x.restartfile << std::flush;
		filename = "coeff" +nstr.str() + "_" +base.idprefix +".bin";
		binofstream bout;
		bout.open(filename.c_str());
		if (bout.error()) {
			*x.gbl->log << "couldn't open coefficient output file " << filename << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (l=0;l<nmodes;++l) 
			bout.writeFloat(coeff(k,l),binio::Double);

		bout.close();
	}
}
#endif
