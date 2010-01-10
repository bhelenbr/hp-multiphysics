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
#include <libbinio/binfile.h>
// #include <veclib/clapack.h>

extern "C" {
	double dlamch_(char *cmach);
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

	if (!input.get(BASE::gbl->idprefix + "_podmodes",nmodes)) input.getwdefault("podmodes",nmodes,nsnapshots);    

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
		modes(i).v.resize(BASE::maxvst,BASE::NV);
		modes(i).e.resize(BASE::maxvst,BASE::em0,BASE::NV);
		modes(i).f.resize(BASE::maxvst,BASE::fm0,BASE::NV);
		modes(i).i.resize(BASE::maxvst,BASE::im0,BASE::NV);
	}
	coeffs.resize(nmodes);
#else
	input.getwdefault(BASE::gbl->idprefix + "_groups",pod_id,0);

	pod_fbdry.resize(BASE::nfbd);
	for (int i=0;i<BASE::nfbd;++i) {
		pod_fbdry(i) = new pod_gen_face_bdry<BASE>(*this,*BASE::fbdry(i));
		pod_fbdry(i)->init(input);
	}
	pod_ebdry.resize(BASE::nebd);
	for (int i=0;i<BASE::nebd;++i) {
		pod_ebdry(i) = new pod_gen_edge_bdry<BASE>(*this,*BASE::ebdry(i));
		pod_ebdry(i)->init(input);
	}
#endif

	return;
}

#ifndef LOWNOISE

template<class BASE> void pod_generate<BASE>::tadvance() {
	int i,j,k,s,l,n,tind,info;
	int lgpx = basis::tet(BASE::log2p)->gpx(), lgpy = basis::tet(BASE::log2p)->gpy(),lgpz = basis::tet(BASE::log2p)->gpz();
	int stridex = MXGP*MXGP, stridey = MXGP;
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
	ugstore.e.reference(BASE::ugbd(0).e);
	ugstore.f.reference(BASE::ugbd(0).f);
	ugstore.i.reference(BASE::ugbd(0).i);

	psimatrix.resize(nsnapshots*nsnapshots);
	psimatrix_recv.resize(nsnapshots*nsnapshots);
	low_noise_dot.resize(BASE::ntet);

	/* GENERATE POD MODES SNAPSHOT COEFFICIENTS */
	psimatrix = 0.0;
	for (s=0;s<nsnapshots;++s) {
		nstr.str("");
		nstr << s+1 << std::flush;
		filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
		BASE::ugbd(0).v.reference(modes(0).v);
		BASE::ugbd(0).e.reference(modes(0).e);
		BASE::ugbd(0).f.reference(modes(0).f);
		BASE::ugbd(0).i.reference(modes(0).i);
		BASE::input(filename, BASE::binary);

		for(l=s;l<nsnapshots;++l) {
			nstr.str("");
			nstr << l+1 << std::flush;
			filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
			BASE::ugbd(0).v.reference(modes(1).v);
			BASE::ugbd(0).e.reference(modes(1).e);
			BASE::ugbd(0).f.reference(modes(1).f);
			BASE::ugbd(0).i.reference(modes(1).i);
			BASE::input(filename, BASE::binary);

			for(tind=0;tind<BASE::ntet;++tind) {          
				/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
				BASE::crdtocht(tind);

				/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
				for(n=0;n<BASE::ND;++n)
					basis::tet(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0)(0)(0), &BASE::dcrd(n,0)(0)(0)(0), &BASE::dcrd(n,1)(0)(0)(0), &BASE::dcrd(n,2)(0)(0)(0),stridex,stridey);

				BASE::ugbd(0).v.reference(modes(0).v);
				BASE::ugbd(0).e.reference(modes(0).e);
				BASE::ugbd(0).f.reference(modes(0).f);
				BASE::ugbd(0).i.reference(modes(0).i);
				BASE::ugtouht(tind,0);
				for(n=0;n<BASE::NV;++n)
					basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0)(0)(0),stridex,stridey);

				BASE::ugbd(0).v.reference(modes(1).v);
				BASE::ugbd(0).e.reference(modes(1).e);
				BASE::ugbd(0).f.reference(modes(1).f);
				BASE::ugbd(0).i.reference(modes(1).i);
				BASE::ugtouht(tind,0);
				for(n=0;n<BASE::NV;++n)
					basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0)(0)(0),stridex, stridey);

				FLT tmp_store = 0.0;
				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpy;++j) {
						for(k=0;k<lgpz;++k) {
							cjcb = basis::tet(BASE::log2p)->wtx(i)*basis::tet(BASE::log2p)->wty(j)*basis::tet(BASE::log2p)->wtz(k)*(BASE::dcrd(0)(0)(i)(j)(k)*(BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k))-BASE::dcrd(0)(1)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k))+BASE::dcrd(0)(2)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k)-BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k)));
							for(n=0;n<BASE::NV;++n) {
								tmp_store += BASE::u(n)(i)(j)(k)*BASE::res(n)(i)(j)(k)*scaling(n)*cjcb;
							}
						}
					}
				}
				low_noise_dot(tind) = tmp_store;
				
				psimatrix(psi1dcounter) += tmp_store; // TO NOT USE LOW NOISE DOT
			}
			/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
			int halfcount,remainder;
			for (remainder=BASE::ntet % 2, halfcount = BASE::ntet/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
				for (tind=0;tind<halfcount;++tind) 
					low_noise_dot(tind) += low_noise_dot(tind+halfcount);
				if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
			}
			// psimatrix(psi1dcounter) = low_noise_dot(0);  // TO USE LOW NOISE DOT!! TEMPO!!
			
			++psi1dcounter;
		}	
		BASE::ugbd(0).v.reference(ugstore.v);
		BASE::ugbd(0).e.reference(ugstore.e);
		BASE::ugbd(0).f.reference(ugstore.f);
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
	Array<FLT,1> work(3*nsnapshots);//3->4 fix me temp?
	DSPEV(jobz,uplo,nsnapshots,psimatrix_recv.data(),eigenvalues.data(),eigenvectors.data(),nsnapshots,work.data(),info);

	if (info != 0) {
		*BASE::gbl->log << "Failed to find eigenmodes " << info << std::endl;
		exit(1);
	}

	*BASE::gbl->log << "eigenvalues "<<  eigenvalues << std::endl;

	eigenvectors.transposeSelf(secondDim,firstDim);  // FORTRAN ROUTINE RETURNS TRANSPOSE

	//reconstruct POD MODES
	for(k=0;k<nmodes;++k)	{
		modes(k).v(Range(0,BASE::npnt-1)) = 0.0;
		modes(k).e(Range(0,BASE::nseg-1)) = 0.0;
		modes(k).f(Range(0,BASE::ntri-1)) = 0.0;
		modes(k).i(Range(0,BASE::ntet-1)) = 0.0;

		/* LOAD SNAPSHOTS AND CALCULATE MODE */
		for(l=0;l<nsnapshots;++l)	{
			nstr.str("");
			nstr << l +1 << std::flush;
			filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
			BASE::input(filename, BASE::binary);

			modes(k).v(Range(0,BASE::npnt-1)) += eigenvectors(l,nmodes -k -1)*BASE::ug.v(Range(0,BASE::npnt-1));
			modes(k).e(Range(0,BASE::nseg-1)) += eigenvectors(l,nmodes -k -1)*BASE::ug.e(Range(0,BASE::nseg-1));
			modes(k).f(Range(0,BASE::ntri-1)) += eigenvectors(l,nmodes -k -1)*BASE::ug.f(Range(0,BASE::ntri-1));
			modes(k).i(Range(0,BASE::ntet-1)) += eigenvectors(l,nmodes -k -1)*BASE::ug.i(Range(0,BASE::ntet-1));
		}
	}

	ugstore.v.reference(BASE::ugbd(1).v);
	ugstore.e.reference(BASE::ugbd(1).e);
	ugstore.f.reference(BASE::ugbd(1).f);
	ugstore.i.reference(BASE::ugbd(1).i);

	psimatrix = 0.0;
	psimatrix_recv = 0.0;

	for(l=0;l<nmodes;++l) {
		BASE::ugbd(1).v.reference(modes(l).v);
		BASE::ugbd(1).e.reference(modes(l).e);
		BASE::ugbd(1).f.reference(modes(l).f);
		BASE::ugbd(1).i.reference(modes(l).i);

		for(tind=0;tind<BASE::ntet;++tind) {          
			/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
			BASE::crdtocht(tind);

			/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(n=0;n<BASE::ND;++n)
				basis::tet(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0)(0)(0), &BASE::dcrd(n,0)(0)(0)(0), &BASE::dcrd(n,1)(0)(0)(0),&BASE::dcrd(n,2)(0)(0)(0),stridex, stridey);

			/* PROJECT MODE TO GAUSS POINTS */
			BASE::ugtouht(tind,1);
			for(n=0;n<BASE::NV;++n)
				basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0)(0)(0),stridex, stridey);

			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						cjcb = basis::tet(BASE::log2p)->wtx(i)*basis::tet(BASE::log2p)->wty(j)*basis::tet(BASE::log2p)->wtz(k)*(BASE::dcrd(0)(0)(i)(j)(k)*(BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k))-BASE::dcrd(0)(1)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k))+BASE::dcrd(0)(2)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k)-BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k)));
						for(n=0;n<BASE::NV;++n) {
							psimatrix(l) += BASE::res(n)(i)(j)(k)*BASE::res(n)(i)(j)(k)*scaling(n)*cjcb;
						}
					}
				}
			}
		}
	}	
	BASE::ugbd(1).v.reference(ugstore.v);
	BASE::ugbd(1).e.reference(ugstore.e);
	BASE::ugbd(1).f.reference(ugstore.f);
	BASE::ugbd(1).i.reference(ugstore.i);
	sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nmodes,blocks::flt_msg,blocks::sum);

	ugstore.v.reference(BASE::ugbd(0).v);
	ugstore.e.reference(BASE::ugbd(0).e);
	ugstore.f.reference(BASE::ugbd(0).f);
	ugstore.i.reference(BASE::ugbd(0).i);

	FLT norm;
	/* RENORMALIZE MODES AND OUTPUT COEFFICIENTS */
	for (k=0;k<nmodes;++k) {
		norm = sqrt(psimatrix_recv(k));
		modes(k).v(Range(0,BASE::npnt-1)) /= norm;
		modes(k).e(Range(0,BASE::nseg-1)) /= norm;
		modes(k).f(Range(0,BASE::ntri-1)) /= norm;
		modes(k).i(Range(0,BASE::ntet-1)) /= norm;

		nstr.str("");
		nstr << k << std::flush;
		filename = "mode" +nstr.str() + "_" + BASE::gbl->idprefix;
		BASE::ugbd(0).v.reference(modes(k).v);
		BASE::ugbd(0).e.reference(modes(k).e);
		BASE::ugbd(0).f.reference(modes(k).f);
		BASE::ugbd(0).i.reference(modes(k).i);
		BASE::output(filename, BASE::binary);
		BASE::output(filename, BASE::tecplot);
	}
	BASE::ugbd(0).v.reference(ugstore.v);
	BASE::ugbd(0).e.reference(ugstore.e);
	BASE::ugbd(0).f.reference(ugstore.f);
	BASE::ugbd(0).i.reference(ugstore.i);


	ugstore.v.reference(BASE::ugbd(1).v);
	ugstore.e.reference(BASE::ugbd(1).e);
	ugstore.f.reference(BASE::ugbd(1).f);
	ugstore.i.reference(BASE::ugbd(1).i);

	psimatrix = 0.0;
	psimatrix_recv = 0.0;
	psi1dcounter=0;
	for (s=0;s<nsnapshots;++s) {
		/* LOAD SNAPSHOT */
		nstr.str("");
		nstr << s+1 << std::flush;
		filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
		BASE::input(filename, BASE::binary);

		for(l=0;l<nmodes;++l) {
			BASE::ugbd(1).v.reference(modes(l).v);
			BASE::ugbd(1).e.reference(modes(l).e);
			BASE::ugbd(1).f.reference(modes(l).f);
			BASE::ugbd(1).i.reference(modes(l).i);

			for(tind=0;tind<BASE::ntet;++tind) {          
				/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
				BASE::crdtocht(tind);

				/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
				for(n=0;n<BASE::ND;++n)
					basis::tet(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0)(0)(0), &BASE::dcrd(n,0)(0)(0)(0), &BASE::dcrd(n,1)(0)(0)(0), &BASE::dcrd(n,2)(0)(0)(0), stridex, stridey);

				/* PROJECT SNAPSHOT TO GAUSS POINTS */
				BASE::ugtouht(tind);
				for(n=0;n<BASE::NV;++n)
					basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0)(0)(0),stridex, stridey);

				/* PROJECT MODE TO GAUSS POINTS */
				BASE::ugtouht(tind,1);
				for(n=0;n<BASE::NV;++n)
					basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0)(0)(0),stridex, stridey);

				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpy;++j) {
						for(k=0;k<lgpz;++k) {
							cjcb = basis::tet(BASE::log2p)->wtx(i)*basis::tet(BASE::log2p)->wty(j)*basis::tet(BASE::log2p)->wtz(k)*(BASE::dcrd(0)(0)(i)(j)(k)*(BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k))-BASE::dcrd(0)(1)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k))+BASE::dcrd(0)(2)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k)-BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k)));
							for(n=0;n<BASE::NV;++n) {
								psimatrix(psi1dcounter) += BASE::u(n)(i)(j)(k)*BASE::res(n)(i)(j)(k)*scaling(n)*cjcb;
							}
						}
					}
				}
			}
			++psi1dcounter;
		}
	}
	BASE::ugbd(1).v.reference(ugstore.v);
	BASE::ugbd(1).e.reference(ugstore.e);
	BASE::ugbd(1).f.reference(ugstore.f);
	BASE::ugbd(1).i.reference(ugstore.i);

	sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nsnapshots*nmodes,blocks::flt_msg,blocks::sum);

	for (k=0;k<nsnapshots;++k) {
		/* OUTPUT COEFFICIENT VECTOR */
		nstr.str("");
		nstr << k+1 << std::flush;
		filename = "coeff" +nstr.str() + "_" +BASE::gbl->idprefix +".bin";
		binofstream bout;
		bout.open(filename.c_str());
		if (bout.error()) {
			*BASE::gbl->log << "couldn't open coefficient output file " << filename;
			exit(1);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (l=0;l<nmodes;++l) 
			bout.writeFloat(psimatrix_recv(k*nmodes +l),binio::Double);

		bout.close();
	}

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
	int i,j,k,s,l,n,tind;
	int lgpx = basis::tet(BASE::log2p)->gpx(), lgpy = basis::tet(BASE::log2p)->gpy(),lgpz = basis::tet(BASE::log2p)->gpz();
	int stridex = MXGP*MXGP, stridey = MXGP;	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	std::istringstream instr;
	FLT cjcb;

	/* MAKE FILES FOR VOLUME MODE SNAPSHOT-PROJECTION */
	for(k=0;k<nsnapshots;++k) {
		nstr.str("");
		nstr << k+1 << std::flush;
		filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
		BASE::input(filename, BASE::binary);
		nstr.str("");
		nstr << k << std::flush;
		filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;

		/* ZERO SNAPSHOTS ON POD BOUNDARY'S */
		for (i=0;i<BASE::nebd;++i)
			pod_ebdry(i)->zero_bdry(BASE::ug);

		BASE::output(filename, BASE::binary);
	}


	for (int eig_ct=0; eig_ct<nsnapshots;++eig_ct) {

		/* ******************************************/
		/* GENERATE POD MODES SNAPSHOT MATRIX       */
		/********************************************/
		int psi1dcounter = 0;
		psimatrix = 0.0;
		for (s=0;s<nsnapshots;++s) {
			nstr.str("");
			nstr << s << std::flush;
			filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;
			BASE::input(filename, BASE::binary);

			for(l=s;l<nsnapshots;++l) {
				nstr.str("");
				nstr << l << std::flush;
				filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;
				BASE::input(filename, BASE::binary, 1);

				for(tind=0;tind<BASE::ntet;++tind) {          
					/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
					BASE::crdtocht(tind);

					/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
					for(n=0;n<BASE::ND;++n)
						basis::tet(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0)(0)(0), &BASE::dcrd(n,0)(0)(0)(0), &BASE::dcrd(n,1)(0)(0)(0), &BASE::dcrd(n,2)(0)(0)(0),stridex,stridey);

					BASE::ugtouht(tind,0);
					for(n=0;n<BASE::NV;++n)
						basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0)(0)(0),stridex,stridey);

					BASE::ugtouht(tind,1);
					for(n=0;n<BASE::NV;++n)
						basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0)(0)(0),stridex,stridey);

					FLT tmp_store = 0.0;
					for(i=0;i<lgpx;++i) {
						for(j=0;j<lgpy;++j) {
							for(k=0;k<lgpz;++k) {
								cjcb = basis::tet(BASE::log2p)->wtx(i)*basis::tet(BASE::log2p)->wty(j)*basis::tet(BASE::log2p)->wtz(k)*(BASE::dcrd(0)(0)(i)(j)(k)*(BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k))-BASE::dcrd(0)(1)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k))+BASE::dcrd(0)(2)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k)-BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k)));
								for(n=0;n<BASE::NV;++n) {
									tmp_store += BASE::u(n)(i)(j)(k)*BASE::res(n)(i)(j)(k)*scaling(n)*cjcb;
								}
							}
						}
					}
					low_noise_dot(tind) = tmp_store;
				}

				/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
				int halfcount,remainder;
				for (remainder=BASE::ntet % 2, halfcount = BASE::ntet/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
					for (tind=0;tind<halfcount;++tind) 
						low_noise_dot(tind) += low_noise_dot(tind+halfcount);
					if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
				}
				psimatrix(psi1dcounter) = low_noise_dot(0);

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
			exit(1);
		}

		*BASE::gbl->log << "eigenvalue "<<  eig_ct << ' ' << eigenvalues(0) << std::endl;


		/*************************************/
		/* LOAD SNAPSHOTS AND CALCULATE MODE */
		/*************************************/
		BASE::ug.v = 0.0;
		BASE::ug.e = 0.0;
		BASE::ug.f = 0.0;
		BASE::ug.i = 0.0;
		for(l=0;l<nsnapshots;++l)	{
			nstr.str("");
			nstr << l << std::flush;
			filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;
			BASE::input(filename, BASE::binary, 1);

			BASE::ug.v(Range(0,BASE::npnt-1)) += eigenvector(l)*BASE::ugbd(1).v(Range(0,BASE::npnt-1));
			BASE::ug.e(Range(0,BASE::nseg-1)) += eigenvector(l)*BASE::ugbd(1).e(Range(0,BASE::nseg-1));
			BASE::ug.f(Range(0,BASE::ntri-1)) += eigenvector(l)*BASE::ugbd(1).f(Range(0,BASE::ntri-1));
			BASE::ug.i(Range(0,BASE::ntet-1)) += eigenvector(l)*BASE::ugbd(1).i(Range(0,BASE::ntet-1));
		}


		/********************/
		/* RENORMALIZE MODE */
		/********************/
		for(tind=0;tind<BASE::ntet;++tind) {          
			/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
			BASE::crdtocht(tind);

			/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(n=0;n<BASE::ND;++n)
				basis::tet(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0)(0)(0), &BASE::dcrd(n,0)(0)(0)(0), &BASE::dcrd(n,1)(0)(0)(0),&BASE::dcrd(n,2)(0)(0)(0),stridex,stridey);

			/* PROJECT MODE TO GAUSS POINTS */
			BASE::ugtouht(tind);
			for(n=0;n<BASE::NV;++n)
				basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0)(0)(0),stridex,stridey);

			FLT tmp_store = 0.0;
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						cjcb = basis::tet(BASE::log2p)->wtx(i)*basis::tet(BASE::log2p)->wty(j)*basis::tet(BASE::log2p)->wtz(k)*(BASE::dcrd(0)(0)(i)(j)(k)*(BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k))-BASE::dcrd(0)(1)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k))+BASE::dcrd(0)(2)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k)-BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k)));
						for(n=0;n<BASE::NV;++n) {
							tmp_store += BASE::res(n)(i)(j)(k)*BASE::res(n)(i)(j)(k)*scaling(n)*cjcb;
						}
					}
				}
			}
			low_noise_dot(tind) = tmp_store;
		}
		/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
		int halfcount,remainder;
		for (remainder=BASE::ntet % 2, halfcount = BASE::ntet/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
			for (tind=0;tind<halfcount;++tind) 
				low_noise_dot(tind) += low_noise_dot(tind+halfcount);
			if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
		}
		double norm = low_noise_dot(0);
		double norm_recv;

		sim::blks.allreduce(&norm,&norm_recv,1,blocks::flt_msg,blocks::sum,pod_id);

		/* DIVIDE BY NORM */
		norm = sqrt(norm_recv);
		BASE::ug.v(Range(0,BASE::npnt-1)) /= norm;
		BASE::ug.e(Range(0,BASE::nseg-1)) /= norm;
		BASE::ug.f(Range(0,BASE::ntri-1)) /= norm;
		BASE::ug.i(Range(0,BASE::ntet-1)) /= norm;

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
		for (s=0;s<nsnapshots;++s) {
			/* LOAD SNAPSHOT */
			nstr.str("");
			nstr << s << std::flush;
			filename = "temp" +nstr.str() + "_" + BASE::gbl->idprefix;
			BASE::input(filename, BASE::binary, 1);

			for(tind=0;tind<BASE::ntet;++tind) {          
				/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
				BASE::crdtocht(tind);

				/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
				for(n=0;n<BASE::ND;++n)
					basis::tet(BASE::log2p)->proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0)(0)(0), &BASE::dcrd(n,0)(0)(0)(0), &BASE::dcrd(n,1)(0)(0)(0),&BASE::dcrd(n,2)(0)(0)(0),stridex,stridey);

				/* PROJECT SNAPSHOT TO GAUSS POINTS */
				BASE::ugtouht(tind);
				for(n=0;n<BASE::NV;++n)
					basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::u(n)(0)(0)(0),stridex,stridey);

				/* PROJECT MODE TO GAUSS POINTS */
				BASE::ugtouht(tind,1);
				for(n=0;n<BASE::NV;++n)
					basis::tet(BASE::log2p)->proj(&BASE::uht(n)(0),&BASE::res(n)(0)(0)(0),stridex,stridey);


				FLT tmp_store = 0.0;
				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpy;++j) {
						for(k=0;k<lgpz;++k) {
							cjcb = basis::tet(BASE::log2p)->wtx(i)*basis::tet(BASE::log2p)->wty(j)*basis::tet(BASE::log2p)->wtz(k)*(BASE::dcrd(0)(0)(i)(j)(k)*(BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k))-BASE::dcrd(0)(1)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(2)(i)(j)(k)-BASE::dcrd(1)(2)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k))+BASE::dcrd(0)(2)(i)(j)(k)*(BASE::dcrd(1)(0)(i)(j)(k)*BASE::dcrd(2)(1)(i)(j)(k)-BASE::dcrd(1)(1)(i)(j)(k)*BASE::dcrd(2)(0)(i)(j)(k)));
							for(n=0;n<BASE::NV;++n) {
								tmp_store += BASE::u(n)(i)(j)(k)*BASE::res(n)(i)(j)(k)*scaling(n)*cjcb;
							}
						}
					}
				}
				low_noise_dot(tind) = tmp_store;
			}
			/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
			int halfcount,remainder;
			for (remainder=BASE::ntet % 2, halfcount = BASE::ntet/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
				for (tind=0;tind<halfcount;++tind) 
					low_noise_dot(tind) += low_noise_dot(tind+halfcount);
				if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
			}
			dotp = low_noise_dot(0);

			sim::blks.allreduce(&dotp,&dotp_recv,1,blocks::flt_msg,blocks::sum,pod_id);

			BASE::ugbd(1).v(Range(0,BASE::npnt-1)) -= dotp_recv*BASE::ug.v(Range(0,BASE::npnt-1));
			BASE::ugbd(1).e(Range(0,BASE::nseg-1)) -= dotp_recv*BASE::ug.e(Range(0,BASE::nseg-1));
			BASE::ugbd(1).f(Range(0,BASE::ntri-1)) -= dotp_recv*BASE::ug.f(Range(0,BASE::ntri-1));
			BASE::ugbd(1).i(Range(0,BASE::ntet-1)) -= dotp_recv*BASE::ug.i(Range(0,BASE::ntet-1));
			BASE::output(filename, BASE::binary, 1);
			coeff(k,eig_ct) = dotp_recv;
		}
	}

	/*****************************/
	/* OUTPUT COEFFICIENT VECTOR */
	/*****************************/
	for (k=0;k<nsnapshots;++k) {
		nstr.str("");
		nstr << k+1 << std::flush;
		filename = "coeff" +nstr.str() + "_" +BASE::gbl->idprefix +".bin";
		binofstream bout;
		bout.open(filename.c_str());
		if (bout.error()) {
			*BASE::gbl->log << "couldn't open coefficient output file " << filename;
			exit(1);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (l=0;l<nmodes;++l) 
			bout.writeFloat(coeff(k,l),binio::Double);

		bout.close();
	}

	/* NOW GENERATE BDRY POD MODES */
	for (int i=0;i<BASE::nfbd;++i)
		pod_fbdry(i)->calculate_modes();
	for (int i=0;i<BASE::nebd;++i)
		pod_ebdry(i)->calculate_modes();

	return;
}
#endif

template<class BASE> void pod_gen_face_bdry<BASE>::init(input_map& input) {
	std::string keyword;

	keyword = base.idprefix + "_pod";
	input.getwdefault(keyword,active,false);

	if (!active) return;

	keyword = base.idprefix + "_podmodes";
	input.getwdefault(keyword,nmodes,x.nsnapshots);

	keyword = base.idprefix + "_pod_id";
	if (!input.get(keyword,pod_id)) {
		*x.gbl->log << "Must provide a pod id for pod boundary\n";
		exit(1);
	}

	/* ERROR CHECK THAT NUMBER IS LISTED IN GROUP LIST */
	if (!input.getline(x.gbl->idprefix +"_groups",keyword)) {
		*x.gbl->log << "group list must contain pod_id for " << base.idprefix << "\n";
		exit(1);
	}

	istringstream mystr(keyword);
	int bnum;
	do {
		if (!(mystr >> bnum)) {
			*x.gbl->log << "group list must contain pod_id for " << base.idprefix << "\n";
			exit(1);
		}
	} while (bnum != pod_id);
}

template<class BASE> void pod_gen_edge_bdry<BASE>::init(input_map& input) {
	std::string keyword;
	
	keyword = base.idprefix + "_pod";
	input.getwdefault(keyword,active,false);
	
	if (!active) return;
	
	keyword = base.idprefix + "_podmodes";
	input.getwdefault(keyword,nmodes,x.nsnapshots);
	
	keyword = base.idprefix + "_pod_id";
	if (!input.get(keyword,pod_id)) {
		*x.gbl->log << "Must provide a pod id for pod boundary\n";
		exit(1);
	}
	
	/* ERROR CHECK THAT NUMBER IS LISTED IN GROUP LIST */
	if (!input.getline(x.gbl->idprefix +"_groups",keyword)) {
		*x.gbl->log << "group list must contain pod_id for " << base.idprefix << "\n";
		exit(1);
	}
	
	istringstream mystr(keyword);
	int bnum;
	do {
		if (!(mystr >> bnum)) {
			*x.gbl->log << "group list must contain pod_id for " << base.idprefix << "\n";
			exit(1);
		}
	} while (bnum != pod_id);
}

template<class BASE> void pod_gen_face_bdry<BASE>::zero_bdry(tet_hp::vefi ug) {

	if (!active) return;

	for(int i=0;i<base.npnt;++i)
		ug.v(i,Range(0,x.NV-1)) = 0.0;	
	for(int i=0;i<base.nseg;++i)
		ug.e(i,Range(0,x.em0-1),Range(0,x.NV-1)) = 0.0;
	for(int i=0;i<base.ntri;++i)
		ug.f(i,Range(0,x.fm0-1),Range(0,x.NV-1)) = 0.0;

}

template<class BASE> void pod_gen_edge_bdry<BASE>::zero_bdry(tet_hp::vefi ug) {
	if (!active) return;
	
	int sind,v0;
	
	for(int j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		ug.v(v0,Range(0,x.NV-1)) = 0.0;
		ug.e(sind,Range(0,x.em0-1),Range(0,x.NV-1)) = 0.0;
	}
	v0 = x.seg(sind).pnt(1);
	ug.v(v0,Range(0,x.NV-1)) = 0.0;
}

/* 2D VERSION TO GENERATE MODES */
template<class BASE> void pod_gen_face_bdry<BASE>::calculate_modes() {

	if (!active) return;

	Array<FLT,1> psimatrix(x.nsnapshots*x.nsnapshots);
	Array<FLT,1> psimatrix_recv(x.nsnapshots*x.nsnapshots);
	Array<FLT,1> low_noise_dot(x.BASE::ntet);
	Array<FLT,1> eigenvalues(x.nsnapshots);
	Array<FLT,1> eigenvector(x.nsnapshots);
	Array<FLT,2> coeff(x.nsnapshots,x.nsnapshots);
	char jobz[2] = "V", range[2] = "I", uplo[2] = "L";
	double vl, vu;
	int il=x.nsnapshots, iu=x.nsnapshots;
	double abstol = 2.*dlamch_("Safe minimum");
	int find,v0,neig,info;
	Array<int,1> iwork(5*x.nsnapshots);
	Array<double,1> work(8*x.nsnapshots);
	Array<int,1> ifail(x.nsnapshots);
	int i,j,k,l,n;
	int lgpx = basis::tet(x.log2p)->gpx(),lgpy = basis::tet(x.log2p)->gpy();
	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	std::istringstream instr;
	FLT cjcb;

	/* MAKE FILES FOR EDGE MODE SNAPSHOT-PROJECTION */
	for(k=0;k<x.nsnapshots;++k) {
		nstr.str("");
		nstr << k+1 << std::flush;
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
				for(int bfind=0;bfind<base.ntri;++bfind) {
					find = base.tri(bfind).gindx;

					/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
					x.crdtocht2d(find);

					/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
					for(n=0;n<x.ND;++n)
						basis::tet(x.log2p)->proj2d(&x.cht(n,0), &x.crd2d(n)(0)(0), &x.dcrd2d(n)(0)(0)(0),&x.dcrd2d(n)(1)(0)(0),MXGP);

					x.ugtouht2d(find,0);
					for(n=0;n<x.NV;++n)
						basis::tet(x.log2p)->proj2d(&x.uht(n)(0),&x.u2d(n)(0)(0),MXGP);

					x.ugtouht2d(find,1);
					for(n=0;n<x.NV;++n)
						basis::tet(x.log2p)->proj2d(&x.uht(n)(0),&x.res2d(n)(0)(0),MXGP);

					FLT tmp_store = 0.0;
					for(i=0;i<lgpx;++i) {
						for(j=0;j<lgpy;++j) {
							cjcb =  basis::tet(x.log2p)->wtx(i)*basis::tet(x.log2p)->wty(j)*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));//fix me temp
							for(n=0;n<x.NV;++n) {
								tmp_store += x.u2d(n)(i)(j)*x.res2d(n)(i)(j)*x.scaling(n)*cjcb;
							}
						}
					}
					low_noise_dot(bsind) = tmp_store;
				}

				/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
				int halfcount,remainder;
				for (remainder=base.nseg % 2, halfcount = base.ntri/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
					for (int bfind=0;bfind<halfcount;++bfind) 
						low_noise_dot(bfind) += low_noise_dot(bfind+halfcount);
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
			exit(1);
		}

		*x.gbl->log << "eigenvalue "<<  eig_ct << ' ' << eigenvalues(0) << std::endl;


		/*************************************/
		/* LOAD SNAPSHOTS AND CALCULATE MODE */
		/*************************************/
		/* FOR NOW I AM GOING TO CALCULATE 2D MODE */
		/* EVEN THOUGH THIS IS UNNECESSARY */
		x.ug.v = 0.0;
		x.ug.e = 0.0;
		x.ug.f = 0.0;
		x.ug.i = 0.0;
		for(l=0;l<x.nsnapshots;++l)	{
			nstr.str("");
			nstr << l << std::flush;
			filename = "temp" +nstr.str() + "_" + base.idprefix;
			x.input(filename, x.binary, 1);
			
			x.ug.v(Range(0,x.npnt-1)) += eigenvector(l)*x.ugbd(1).v(Range(0,x.npnt-1));
			x.ug.e(Range(0,x.nseg-1)) += eigenvector(l)*x.ugbd(1).e(Range(0,x.nseg-1));
			x.ug.f(Range(0,x.ntri-1)) += eigenvector(l)*x.ugbd(1).f(Range(0,x.ntri-1));
		}


		/********************/
		/* RENORMALIZE MODE */
		/********************/
		for(int bfind=0;bfind<base.ntri;++bfind) {
			find = base.tri(bfind).gindx;

			/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
			x.crdtocht2d(find);

			/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(n=0;n<x.ND;++n)
				basis::tet(x.log2p)->proj2d(&x.cht(n,0), &x.crd2d(n)(0)(0), &x.dcrd2d(n)(0)(0)(0),&x.dcrd2d(n)(1)(0)(0),MXGP);

			x.ugtouht2d(find);
			for(n=0;n<x.NV;++n)
				basis::tet(x.log2p)->proj2d(&x.uht(n)(0),&x.res2d(n)(0)(0),MXGP);

			FLT tmp_store = 0.0;
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					cjcb =  basis::tet(x.log2p)->wtx(i)* basis::tet(x.log2p)->wty(j)*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));//temp fix me
					for(n=0;n<x.NV;++n) {
						tmp_store += x.res2d(n)(i,j)*x.res2d(n)(i,j)*x.scaling(n)*cjcb;
					}
				}
			}
			low_noise_dot(bsind) = tmp_store;
		}
		/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
		int halfcount,remainder;
		for (remainder=base.ntri % 2, halfcount = base.ntri/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
			for (int bfind=0;bfind<halfcount;++bfind) 
				low_noise_dot(bfind) += low_noise_dot(bfind+halfcount);
			if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
		}
		double norm = low_noise_dot(0);
		double norm_recv;

		sim::blks.allreduce(&norm,&norm_recv,1,blocks::flt_msg,blocks::sum,pod_id);

		/* DIVIDE BY NORM */
		norm = sqrt(norm_recv);
		x.ug.v(Range(0,x.npnt-1)) /= norm;
		x.ug.e(Range(0,x.nseg-1)) /= norm;
		x.ug.f(Range(0,x.ntri-1)) /= norm;

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
			exit(1);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (int bvind=0;bvind<base.npnt;++bvind) {
			v0 = x.pnt(bvind).gindx;
			for (n=0;n<x.NV;++n)
				bout.writeFloat(x.ug.v(v0,n),binio::Double);
		}		

		for (int bsind=0;bsind<base.nseg;++bsind) {
			sind = base.seg(bsind).gindx;
			for (int m=0;m<x.em0;++m)
				for (n=0;n<x.NV;++n)
					bout.writeFloat(x.ug.e(sind,m,n),binio::Double);
		}
		for (int bfind=0;bfind<base.ntri;++bfind) {
			find = base.tri(bfind).gindx;
			for (int m=0;m<x.fm0;++m)
				for (n=0;n<x.NV;++n)
					bout.writeFloat(x.ug.f(find,m,n),binio::Double);
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
			for(int bfind=0;bfind<base.ntri;++bfind) {
				find = base.tri(bfind).gindx;

				/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
				x.crdtocht2d(find);

				/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
				for(n=0;n<x.ND;++n)
					basis::tet(x.log2p)->proj2d(&x.cht(n,0), &x.crd2d(n)(0)(0), &x.dcrd2d(n)(0)(0)(0),&x.dcrd2d(n)(1)(0)(0),MXGP);

				x.ugtouht2d(find,0);
				for(n=0;n<x.NV;++n)
					basis::tet(x.log2p)->proj2d(&x.uht(n)(0),&x.u2d(n)(0)(0),MXGP);

				x.ugtouht2d(find,1);
				for(n=0;n<x.NV;++n)
					basis::tet(x.log2p)->proj2d(&x.uht(n)(0),&x.res2d(n)(0)(0),MXGP);

				FLT tmp_store = 0.0;
				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpy;++j) {
						cjcb =  basis::tet(x.log2p)->wtx(i)*basis::tet(x.log2p)->wty(j)*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));//temp fix me
						for(n=0;n<x.NV;++n) {
							tmp_store += x.u2d(n)(i)(j)*x.res2d(n)(i)(j)*x.scaling(n)*cjcb;
						}
					}
				}
				low_noise_dot(bsind) = tmp_store;
			}

			/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
			int halfcount,remainder;
			for (remainder=base.ntri % 2, halfcount = base.ntri/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
				for (int bfind=0;bfind<halfcount;++bfind) 
					low_noise_dot(bfind) += low_noise_dot(bfind+halfcount);
				if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
			}
			dotp = low_noise_dot(0);

			sim::blks.allreduce(&dotp,&dotp_recv,1,blocks::flt_msg,blocks::sum,pod_id);
			x.ugbd(1).v(Range(0,x.npnt-1)) -= dotp_recv*x.ug.v(Range(0,x.npnt-1));
			x.ugbd(1).e(Range(0,x.nseg-1)) -= dotp_recv*x.ug.e(Range(0,x.nseg-1));
			x.ugbd(1).f(Range(0,x.ntri-1)) -= dotp_recv*x.ug.f(Range(0,x.ntri-1));
			x.output(filename, x.binary, 1);
			coeff(k,eig_ct) = dotp_recv;
		}
	}

	/*****************************/
	/* OUTPUT COEFFICIENT VECTOR */
	/*****************************/
	for (k=0;k<x.nsnapshots;++k) {
		nstr.str("");
		nstr << k+1 << std::flush;
		filename = "coeff" +nstr.str() + "_" +base.idprefix +".bin";
		binofstream bout;
		bout.open(filename.c_str());
		if (bout.error()) {
			*x.gbl->log << "couldn't open coefficient output file " << filename;
			exit(1);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (l=0;l<nmodes;++l) 
			bout.writeFloat(coeff(k,l),binio::Double);

		bout.close();
	}
}

/* 1D VERSION TO GENERATE MODES */
template<class BASE> void pod_gen_edge_bdry<BASE>::calculate_modes() {
	
	if (!active) return;

	Array<FLT,1> psimatrix(x.nsnapshots*x.nsnapshots);
	Array<FLT,1> psimatrix_recv(x.nsnapshots*x.nsnapshots);
	Array<FLT,1> low_noise_dot(x.BASE::ntet);
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
	int lgpx = basis::tet(x.log2p)->gpx();
	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	std::istringstream instr;
	FLT cjcb;

	/* MAKE FILES FOR EDGE MODE SNAPSHOT-PROJECTION */
	for(k=0;k<x.nsnapshots;++k) {
		nstr.str("");
		nstr << k+1 << std::flush;
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
				for(int bsind=0;bsind<base.nseg;++bsind) {
					sind = base.seg(bsind);

					/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
					x.crdtocht1d(sind);

					/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
					for(n=0;n<x.ND;++n)
						basis::tet(x.log2p)->proj1d(&x.cht(n,0), &x.crd1d(n)(0), &x.dcrd1d(n)(0));

					x.ugtouht1d(sind,0);
					for(n=0;n<x.NV;++n)
						basis::tet(x.log2p)->proj1d(&x.uht(n)(0),&x.u1d(n)(0));

					x.ugtouht1d(sind,1);
					for(n=0;n<x.NV;++n)
						basis::tet(x.log2p)->proj1d(&x.uht(n)(0),&x.res1d(n)(0));

					FLT tmp_store = 0.0;
					for(i=0;i<lgpx;++i) {
						cjcb =  basis::tet(x.log2p)->wtx(i)*sqrt(x.dcrd1d(0)(i)*x.dcrd1d(0)(i) +x.dcrd1d(1)(i)*x.dcrd1d(1)(i)+x.dcrd1d(2)(i)*x.dcrd1d(2)(i));//fix me temp
						for(n=0;n<x.NV;++n) {
							tmp_store += x.u1d(n)(i)*x.res1d(n)(i)*x.scaling(n)*cjcb;
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
			exit(1);
		}

		*x.gbl->log << "eigenvalue "<<  eig_ct << ' ' << eigenvalues(0) << std::endl;


		/*************************************/
		/* LOAD SNAPSHOTS AND CALCULATE MODE */
		/*************************************/
		/* FOR NOW I AM GOING TO CALCULATE 2D MODE */
		/* EVEN THOUGH THIS IS UNNECESSARY */
		x.ug.v = 0.0;
		x.ug.e = 0.0;
		x.ug.f = 0.0;
		x.ug.i = 0.0;
		for(l=0;l<x.nsnapshots;++l)	{
			nstr.str("");
			nstr << l << std::flush;
			filename = "temp" +nstr.str() + "_" + base.idprefix;
			x.input(filename, x.binary, 1);

			for (int bsind=0;bsind<base.nseg;++bsind) {
				sind = base.seg(bsind);
				v0 = x.seg(sind).pnt(0);
				x.ug.v(v0) += eigenvector(l)*x.ugbd(1).v(v0);
				x.ug.e(sind) += eigenvector(l)*x.ugbd(1).e(sind);
			}
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
				basis::tet(x.log2p)->proj1d(&x.cht(n,0), &x.crd1d(n)(0), &x.dcrd1d(n)(0));

			x.ugtouht1d(sind);
			for(n=0;n<x.NV;++n)
				basis::tet(x.log2p)->proj1d(&x.uht(n)(0),&x.res1d(n)(0));

			FLT tmp_store = 0.0;
			for(i=0;i<lgpx;++i) {
				cjcb =  basis::tet(x.log2p)->wtx(i)*sqrt(x.dcrd1d(0)(i)*x.dcrd1d(0)(i) +x.dcrd1d(1)(i)*x.dcrd1d(1)(i)+x.dcrd1d(2)(i)*x.dcrd1d(2)(i));//fix me temp
				for(n=0;n<x.NV;++n) {
					tmp_store += x.res1d(n)(i)*x.res1d(n)(i)*x.scaling(n)*cjcb;
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
		x.ug.e(Range(0,x.nseg-1)) /= norm;

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
			exit(1);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (int bsind=0;bsind<base.nseg;++bsind) {
			sind = base.seg(bsind);
			v0 = x.seg(sind).pnt(0);
			for (n=0;n<x.NV;++n)
				bout.writeFloat(x.ug.v(v0,n),binio::Double);
		}
		v0 = x.seg(sind).pnt(1);
		for (n=0;n<x.NV;++n)
			bout.writeFloat(x.ug.v(v0,n),binio::Double);		

		for (int bsind=0;bsind<base.nseg;++bsind) {
			sind = base.seg(bsind);
			for (int m=0;m<x.em0;++m)
				for (n=0;n<x.NV;++n)
					bout.writeFloat(x.ug.e(sind,m,n),binio::Double);
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
					basis::tet(x.log2p)->proj1d(&x.cht(n,0), &x.crd1d(n)(0), &x.dcrd1d(n)(0));

				x.ugtouht1d(sind,0);
				for(n=0;n<x.NV;++n)
					basis::tet(x.log2p)->proj1d(&x.uht(n)(0),&x.u1d(n)(0));

				x.ugtouht1d(sind,1);
				for(n=0;n<x.NV;++n)
					basis::tet(x.log2p)->proj1d(&x.uht(n)(0),&x.res1d(n)(0));

				FLT tmp_store = 0.0;
				for(i=0;i<lgpx;++i) {
					cjcb =  basis::tet(x.log2p)->wtx(i)*sqrt(x.dcrd1d(0)(i)*x.dcrd1d(0)(i) +x.dcrd1d(1)(i)*x.dcrd1d(1)(i)+x.dcrd1d(2)(i)*x.dcrd1d(2)(i));//fix me temp
					for(n=0;n<x.NV;++n) {
						tmp_store += x.u1d(n)(i)*x.res1d(n)(i)*x.scaling(n)*cjcb;
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
			x.ugbd(1).e(Range(0,x.nseg-1)) -= dotp_recv*x.ug.e(Range(0,x.nseg-1));
			x.output(filename, x.binary, 1);
			coeff(k,eig_ct) = dotp_recv;
		}
	}

	/*****************************/
	/* OUTPUT COEFFICIENT VECTOR */
	/*****************************/
	for (k=0;k<x.nsnapshots;++k) {
		nstr.str("");
		nstr << k+1 << std::flush;
		filename = "coeff" +nstr.str() + "_" +base.idprefix +".bin";
		binofstream bout;
		bout.open(filename.c_str());
		if (bout.error()) {
			*x.gbl->log << "couldn't open coefficient output file " << filename;
			exit(1);
		}
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
		bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

		for (l=0;l<nmodes;++l) 
			bout.writeFloat(coeff(k,l),binio::Double);

		bout.close();
	}
}
