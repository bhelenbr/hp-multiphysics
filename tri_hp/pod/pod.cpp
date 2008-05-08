/*
 *  pod.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/26/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include <myblas.h>
// #include <veclib/clapack.h>

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
    
    if (!input.get(BASE::gbl->idprefix + "_podmodes",nsnapshots)) input.getwdefault("podmodes",nmodes,nsnapshots);    
    
    /* THIS IS TO CHANGE THE WAY SNAPSHOT MATRIX ENTRIES ARE FORMED */
    scaling.resize(BASE::NV);
    scaling = 1;
    if (input.getline(BASE::gbl->idprefix + "_scale_vector",linebuff) || input.getline("scale_vector",linebuff)) {
        instr.str(linebuff);
        for(i=0;i<BASE::NV;++i)
            instr >> scaling(i);
    }
    
    nmodes = MAX(nmodes,2);
    modes.resize(nmodes);
    for(i=0;i<nmodes;++i) {
        modes(i).v.resize(BASE::maxpst,BASE::NV);
        modes(i).s.resize(BASE::maxpst,BASE::sm0,BASE::NV);
        modes(i).i.resize(BASE::maxpst,BASE::im0,BASE::NV);
    }
    coeffs.resize(nmodes);
    
    return;
}
    
template<class BASE> void pod_generate<BASE>::tadvance() {
    int i,j,k,l,n,tind,info;
    int lgpx = basis::tri(BASE::log2p).gpx, lgpn = basis::tri(BASE::log2p).gpn;
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
	low_noise_dot.resize(ntri);
    
    /* GENERATE POD MODES SNAPSHOT COEFFICIENTS */
    psimatrix = 0.0;
    for (k=0;k<nsnapshots;++k) {
        nstr.str("");
        nstr << k+1 << std::flush;
        filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
        BASE::ugbd(0).v.reference(modes(0).v);
        BASE::ugbd(0).s.reference(modes(0).s);
        BASE::ugbd(0).i.reference(modes(0).i);
        BASE::input(filename, BASE::text);
        
        for(l=k;l<nsnapshots;++l) {
            nstr.str("");
            nstr << l+1 << std::flush;
            filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
            BASE::ugbd(0).v.reference(modes(1).v);
            BASE::ugbd(0).s.reference(modes(1).s);
            BASE::ugbd(0).i.reference(modes(1).i);
            BASE::input(filename, BASE::text);
            
            for(tind=0;tind<BASE::ntri;++tind) {          
                /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
                BASE::crdtocht(tind);

                /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
                for(n=0;n<BASE::ND;++n)
                    basis::tri(BASE::log2p).proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);
                                 
                BASE::ugbd(0).v.reference(modes(0).v);
                BASE::ugbd(0).s.reference(modes(0).s);
                BASE::ugbd(0).i.reference(modes(0).i);
                BASE::ugtouht(tind,0);
                for(n=0;n<BASE::NV;++n)
                    basis::tri(BASE::log2p).proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);
                    
                BASE::ugbd(0).v.reference(modes(1).v);
                BASE::ugbd(0).s.reference(modes(1).s);
                BASE::ugbd(0).i.reference(modes(1).i);
                BASE::ugtouht(tind,0);
                for(n=0;n<BASE::NV;++n)
                    basis::tri(BASE::log2p).proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);
                
				FLT tmp_store = 0.0;
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                        cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p).wtx(i)*basis::tri(BASE::log2p).wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
                        for(n=0;n<BASE::NV;++n) {
							tmp_store += BASE::u(n)(i,j)*BASE::res(n)(i,j)*scaling(n)*cjcb;
                        }
                    }
                }
				low_noise_dot(tind) = tmp_store;
            }
			/* BALANCED ADDITION FOR MINIMAL ROUNDOFF */
			int halfcount,remainder;
			while (oldhalfcount > 1) {
			for (remainder=ntri % 2, halfcount = ntri/2; halfcount>0; remainder = halfcount % 2, halfcount /= 2) {
				for (tind=0;tind<halfcount;++tind) 
					low_noise_dot(tind) += low_noise_dot(tind+halfcount);
				if (remainder) low_noise_dot(halfcount-1) += low_noise_dot(2*halfcount);
			}
			psimatrix(psi1dcounter) = low_noise_dot(0);
            ++psi1dcounter;
        }	
        BASE::ugbd(0).v.reference(ugstore.v);
        BASE::ugbd(0).s.reference(ugstore.s);
        BASE::ugbd(0).i.reference(ugstore.i);
    }
    sim::blks.allreduce(psimatrix.data(),psimatrix_recv.data(),nsnapshots*(nsnapshots+1)/2,blocks::flt_msg,blocks::sum);
     
    Array<FLT,1> eigenvalues(nsnapshots);
    Array<FLT,2> eigenvectors(nsnapshots,nsnapshots);
    
/*         
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
    exit(1);
*/
    char jobz[2] = "V", uplo[2] = "L";
    Array<FLT,1> work(3*nsnapshots);
    DSPEV(jobz,uplo,nsnapshots,psimatrix_recv.data(),eigenvalues.data(),eigenvectors.data(),nsnapshots,work.data(),info);
    
    if (info != 0) {
        *BASE::gbl->log << "Failed to find eigenmodes " << info << std::endl;
        exit(1);
    }
    
    *BASE::gbl->log << "eigenvalues "<<  eigenvalues << std::endl;
    
    eigenvectors.transposeSelf(secondDim,firstDim);  // FORTRAN ROUTINE RETURNS TRANSPOSE

    //reconstruct POD MODES
    for(k=0;k<nmodes;++k)	{
        modes(k).v = 0.0;
        modes(k).s = 0.0;
        modes(k).i = 0.0;
        
        /* LOAD SNAPSHOTS AND CALCULATE MODE */
        for(l=0;l<nsnapshots;++l)	{
            nstr.str("");
            nstr << l +1 << std::flush;
            filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
            BASE::input(filename, BASE::text);

            modes(k).v += eigenvectors(l,nmodes -k -1)*BASE::ug.v;
            modes(k).s += eigenvectors(l,nmodes -k -1)*BASE::ug.s;
            modes(k).i += eigenvectors(l,nmodes -k -1)*BASE::ug.i;
        }
    }

    ugstore.v.reference(BASE::ugbd(1).v);
    ugstore.s.reference(BASE::ugbd(1).s);
    ugstore.i.reference(BASE::ugbd(1).i);
    
    psimatrix = 0.0;
    psimatrix_recv = 0.0;

    for(l=0;l<nmodes;++l) {
        BASE::ugbd(1).v.reference(modes(l).v);
        BASE::ugbd(1).s.reference(modes(l).s);
        BASE::ugbd(1).i.reference(modes(l).i);
        
        for(tind=0;tind<BASE::ntri;++tind) {          
            /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
            BASE::crdtocht(tind);

            /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
            for(n=0;n<BASE::ND;++n)
                basis::tri(BASE::log2p).proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);
                             
            /* PROJECT MODE TO GAUSS POINTS */
            BASE::ugtouht(tind,1);
            for(n=0;n<BASE::NV;++n)
                basis::tri(BASE::log2p).proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);
            
            for(i=0;i<lgpx;++i) {
                for(j=0;j<lgpn;++j) {
                    cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p).wtx(i)*basis::tri(BASE::log2p).wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
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
        modes(k).v /= norm;
        modes(k).s /= norm;
        modes(k).i /= norm;
        
        nstr.str("");
        nstr << k << std::flush;
        filename = "mode" +nstr.str() + "_" + BASE::gbl->idprefix;
        BASE::ugbd(0).v.reference(modes(k).v);
        BASE::ugbd(0).s.reference(modes(k).s);
        BASE::ugbd(0).i.reference(modes(k).i);
        BASE::output(filename, BASE::text);
        BASE::output(filename, BASE::tecplot);
    }
    BASE::ugbd(0).v.reference(ugstore.v);
    BASE::ugbd(0).s.reference(ugstore.s);
    BASE::ugbd(0).i.reference(ugstore.i);
    

    ugstore.v.reference(BASE::ugbd(1).v);
    ugstore.s.reference(BASE::ugbd(1).s);
    ugstore.i.reference(BASE::ugbd(1).i);

    psimatrix = 0.0;
    psimatrix_recv = 0.0;
	psi1dcounter=0.0;
    for (k=0;k<nsnapshots;++k) {
        /* LOAD SNAPSHOT */
        nstr.str("");
        nstr << k+1 << std::flush;
        filename = "rstrt" +nstr.str() + "_" + BASE::gbl->idprefix +".d0";
        BASE::input(filename, BASE::text);
        
        for(l=0;l<nmodes;++l) {
            BASE::ugbd(1).v.reference(modes(l).v);
            BASE::ugbd(1).s.reference(modes(l).s);
            BASE::ugbd(1).i.reference(modes(l).i);
            
            for(tind=0;tind<BASE::ntri;++tind) {          
                /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
                BASE::crdtocht(tind);

                /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
                for(n=0;n<BASE::ND;++n)
                    basis::tri(BASE::log2p).proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);
                                 
                /* PROJECT SNAPSHOT TO GAUSS POINTS */
                BASE::ugtouht(tind);
                for(n=0;n<BASE::NV;++n)
                    basis::tri(BASE::log2p).proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);
                    
                /* PROJECT MODE TO GAUSS POINTS */
                BASE::ugtouht(tind,1);
                for(n=0;n<BASE::NV;++n)
                    basis::tri(BASE::log2p).proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);
                
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                        cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p).wtx(i)*basis::tri(BASE::log2p).wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
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
        nstr << k+1 << std::flush;
        filename = "coeff" +nstr.str() + "_" +BASE::gbl->idprefix +".txt";
        std::ofstream out;
        out.open(filename.c_str());
        if (!out) {
            *BASE::gbl->log << "couldn't open coefficient output file " << filename;
            exit(1);
        }
        
        for (l=0;l<nmodes;++l) 
            out << psimatrix_recv(k*nmodes +l) << '\n';
            
        out.close();
    }

    return;
}
    
    
template<class BASE> void pod_simulate<BASE>::init(input_map& input, void *gin) {
    std::string filename,keyword,linebuff;
    std::ostringstream nstr;
    std::istringstream instr;
    int i;

    /* Initialize base class */
    BASE::init(input,gin);
        
    if (!input.get(BASE::gbl->idprefix + "_podmodes",nmodes)) input.getwdefault("podmodes",nmodes,5);    
    
    vsi ugstore;
    ugstore.v.reference(BASE::ugbd(0).v);
    ugstore.s.reference(BASE::ugbd(0).s);
    ugstore.i.reference(BASE::ugbd(0).i);

    modes.resize(nmodes);
    for(i=0;i<nmodes;++i) {
        nstr.str("");
        nstr << i << std::flush;
        filename = "mode" +nstr.str() +"_" +BASE::gbl->idprefix;
        nstr.clear();
        modes(i).v.resize(BASE::maxpst,BASE::NV);
        modes(i).s.resize(BASE::maxpst,BASE::sm0,BASE::NV);
        modes(i).i.resize(BASE::maxpst,BASE::im0,BASE::NV);
        BASE::ugbd(0).v.reference(modes(i).v);
        BASE::ugbd(0).s.reference(modes(i).s);
        BASE::ugbd(0).i.reference(modes(i).i);
        BASE::input(filename, BASE::text);
    }
    BASE::ugbd(0).v.reference(ugstore.v);
    BASE::ugbd(0).s.reference(ugstore.s);
    BASE::ugbd(0).i.reference(ugstore.i);
            
    coeffs.resize(nmodes);
    rsdls.resize(nmodes);
    rsdls_recv.resize(nmodes);
    jacobian.resize(nmodes,nmodes);
    ipiv.resize(nmodes);
    
    int initfile;
    input.getwdefault("initfile",initfile,1);
    nstr.str("");
    nstr << initfile << std::flush;
    filename = "coeff" +nstr.str() +"_" +BASE::gbl->idprefix +".txt";
    std::ifstream in;
    in.open(filename.c_str());
    if (!in) {
        *BASE::gbl->log << "couldn't open coefficient input file " << filename;
        exit(1);
    }

    /* CONSTRUCT INITIAL SOLUTION DESCRIPTION */
    BASE::ug.v = 0;
    BASE::ug.s = 0;
    BASE::ug.i = 0;
                
    for (int l=0;l<nmodes;++l) {
        in >> coeffs(l);
        BASE::ug.v(Range(0,BASE::npnt-1)) += coeffs(l)*modes(l).v(Range(0,BASE::npnt-1));
        BASE::ug.s(Range(0,BASE::nseg-1)) += coeffs(l)*modes(l).s(Range(0,BASE::nseg-1));
        BASE::ug.i(Range(0,BASE::ntri-1)) += coeffs(l)*modes(l).i(Range(0,BASE::ntri-1));
    }
    
    return;
}

template<class BASE> void pod_simulate<BASE>::rsdl(int stage) {

    BASE::rsdl(stage);

    rsdls = 0.0;
    
    /* APPLY VERTEX DIRICHLET B.C.'S */
    for(int i=0;i<BASE::nebd;++i)
        BASE::hp_ebdry(i)->vdirichlet();
    
    for(int i=0;i<BASE::nvbd;++i)
        BASE::hp_vbdry(i)->vdirichlet2d();
        
    /* APPLY DIRCHLET B.C.S TO MODE */
    for(int i=0;i<BASE::nebd;++i)
        for(int sm=0;sm<basis::tri(BASE::log2p).sm;++sm)
            BASE::hp_ebdry(i)->sdirichlet(sm);
    
    for (int k = 0; k < nmodes; ++k) {
        rsdls(k) = 0.0;
        for(int i=0; i<BASE::npnt;++i)
            for(int n=0;n<BASE::NV;++n)
                rsdls(k) += modes(k).v(i,n)*BASE::gbl->res.v(i,n);
                
        for(int i=0; i<BASE::nseg;++i)
            for(int sm=0;sm<basis::tri(BASE::log2p).sm;++sm)
                for(int n=0;n<BASE::NV;++n)
                    rsdls(k) += modes(k).s(i,sm,n)*BASE::gbl->res.s(i,sm,n);
                    
        for(int i=0; i<BASE::ntri;++i)
            for(int im=0;im<basis::tri(BASE::log2p).im;++im)
                for(int n=0;n<BASE::NV;++n)
                    rsdls(k) += modes(k).i(i,im,n)*BASE::gbl->res.i(i,im,n);
                    
    }
    
    sim::blks.allreduce(rsdls.data(),rsdls_recv.data(),nmodes,blocks::flt_msg,blocks::sum);

    return;
}


template<class BASE> void pod_simulate<BASE>::setup_preconditioner() {
    BASE::setup_preconditioner();
    rsdl(sim::NSTAGE);
    
    /* STORE BASELINE IN LAST COLUMN */
    jacobian(Range(0,nmodes-1),nmodes-1) = rsdls_recv;
    BASE::gbl->ug0.v = BASE::ug.v;
    BASE::gbl->ug0.s = BASE::ug.s;
    BASE::gbl->ug0.i = BASE::ug.i;

    for (int modeloop = 0; modeloop < nmodes; ++modeloop) {
        /* PERTURB EACH COEFFICIENT */
        BASE::gbl->res.v(Range(0,BASE::npnt-1),Range::all()) = 1.0e-4*modes(modeloop).v(Range(0,BASE::npnt-1),Range::all());
        BASE::gbl->res.s(Range(0,BASE::nseg-1),Range::all(),Range::all()) = 1.0e-4*modes(modeloop).s(Range(0,BASE::nseg-1),Range::all(),Range::all());
        BASE::gbl->res.i(Range(0,BASE::ntri-1),Range::all(),Range::all()) = 1.0e-4*modes(modeloop).i(Range(0,BASE::ntri-1),Range::all(),Range::all());
        
        /* APPLY VERTEX DIRICHLET B.C.'S */
        for(int i=0;i<BASE::nebd;++i)
            BASE::hp_ebdry(i)->vdirichlet();
        
        for(int i=0;i<BASE::nvbd;++i)
            BASE::hp_vbdry(i)->vdirichlet2d();
            
        /* APPLY DIRCHLET B.C.S TO MODE */
        for(int i=0;i<BASE::nebd;++i)
            for(int sm=0;sm<basis::tri(BASE::log2p).sm;++sm)
                BASE::hp_ebdry(i)->sdirichlet(sm);

        BASE::ug.v(Range(0,BASE::npnt-1),Range::all()) = BASE::gbl->ug0.v(Range(0,BASE::npnt-1),Range::all()) +BASE::gbl->res.v(Range(0,BASE::npnt-1),Range::all());
        BASE::ug.s(Range(0,BASE::nseg-1),Range::all(),Range::all()) = BASE::gbl->ug0.s(Range(0,BASE::nseg-1),Range::all()) +BASE::gbl->res.s(Range(0,BASE::nseg-1),Range::all(),Range::all());
        BASE::ug.i(Range(0,BASE::ntri-1),Range::all(),Range::all()) = BASE::gbl->ug0.i(Range(0,BASE::ntri-1),Range::all(),Range::all()) +BASE::gbl->res.i(Range(0,BASE::ntri-1),Range::all(),Range::all());
        
        rsdl(sim::NSTAGE);
        
         /* STORE IN ROW */
        jacobian(Range(0,nmodes-1),modeloop) = (rsdls_recv -jacobian(Range(0,nmodes-1),nmodes-1))/1.0e-4;
    }
            
    /* RESTORE UG */
    BASE::ug.v = BASE::gbl->ug0.v;
    BASE::ug.s = BASE::gbl->ug0.s;
    BASE::ug.i = BASE::gbl->ug0.i;
                
    /* FACTORIZE PRECONDITIONER */
    int info;
    GETRF(nmodes,nmodes,jacobian.data(),nmodes,ipiv.data(),info);
    if (info != 0) {
        printf("DGETRF FAILED FOR POD JACOBIAN\n");
        exit(1);
    }
    return;
}

template<class BASE> void pod_simulate<BASE>::update() {
    char trans[] = "T";
    int info;
    
    rsdl(sim::NSTAGE);
  
    GETRS(trans,nmodes,1,jacobian.data(),nmodes,ipiv.data(),rsdls_recv.data(),nmodes,info);
    if (info != 0) {
        printf("DGETRS FAILED FOR POD UPDATE\n");
        exit(1);
    }
    coeffs -= rsdls_recv;
    
    BASE::gbl->res.v = 0.0;
    BASE::gbl->res.s = 0.0;
    BASE::gbl->res.i = 0.0;
    
    for (int m=0;m<nmodes;++m) {
        BASE::gbl->res.v(Range(0,BASE::npnt-1),Range::all()) += rsdls_recv(m)*modes(m).v(Range(0,BASE::npnt-1),Range::all());
        BASE::gbl->res.s(Range(0,BASE::nseg-1),Range::all(),Range::all()) += rsdls_recv(m)*modes(m).s(Range(0,BASE::nseg-1),Range::all(),Range::all());
        BASE::gbl->res.i(Range(0,BASE::ntri-1),Range::all(),Range::all()) += rsdls_recv(m)*modes(m).i(Range(0,BASE::ntri-1),Range::all(),Range::all());
    }
    
    /* APPLY VERTEX DIRICHLET B.C.'S */
    for(int i=0;i<BASE::nebd;++i)
        BASE::hp_ebdry(i)->vdirichlet();
    
    for(int i=0;i<BASE::nvbd;++i)
        BASE::hp_vbdry(i)->vdirichlet2d();
        
    /* APPLY DIRCHLET B.C.S TO MODE */
    for(int i=0;i<BASE::nebd;++i)
        for(int sm=0;sm<basis::tri(BASE::log2p).sm;++sm)
            BASE::hp_ebdry(i)->sdirichlet(sm);
            

    BASE::ug.v(Range(0,BASE::npnt-1),Range::all()) -= BASE::gbl->res.v(Range(0,BASE::npnt-1),Range::all());
    BASE::ug.s(Range(0,BASE::nseg-1),Range::all(),Range::all()) -= BASE::gbl->res.s(Range(0,BASE::nseg-1),Range::all(),Range::all());
    BASE::ug.i(Range(0,BASE::ntri-1),Range::all(),Range::all()) -= BASE::gbl->res.i(Range(0,BASE::ntri-1),Range::all(),Range::all());
    
    return;
}


template<class BASE> FLT pod_simulate<BASE>::maxres() {
    int i;
    FLT mxr;
        
    mxr = 0.0;
        
    for(i=0;i<nmodes;++i)
        mxr = MAX(fabs(rsdls_recv(i)),mxr);
        
    *BASE::gbl->log << ' ' << mxr << ' ';

    return(mxr);
}

