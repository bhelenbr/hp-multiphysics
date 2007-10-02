#include "mesh.h"
#include <cfloat>
#include <iostream>
#include <stdlib.h>
#include <utilities.h>
#include <string.h>
#include "block.h"
// #include <boost/bind.hpp>
// #include <boost/thread.hpp>

void mesh::connect(multigrid_interface& in) {
    mesh &tgt = dynamic_cast<mesh&>(in);
    
    /* SET UP MULTIGRID STUFF */
    fine = &in;
    tgt.coarse = this;    
    if (fcnnct.ubound(firstDim) < maxvst-1) fcnnct.resize(maxvst);
    if (tgt.ccnnct.ubound(firstDim) < tgt.maxvst-1) tgt.ccnnct.resize(tgt.maxvst);
    
#define OLDRECONNECT
#ifdef OLDRECONNECT
    coarsen(1.6,tgt);
#else
    coarsen2(1.5,tgt);
#endif
    mgconnect(tgt,fcnnct);
    tgt.mgconnect(*this,tgt.ccnnct);
    
    /* THIS IS FOR DIAGNOSIS OF MULTI-GRID FAILURES */
    checkintegrity();
              
    if (gbl->adapt_output) {
        std::string name, fname;
        std::string adapt_file;
        std::ostringstream nstr;
        nstr.str("");
        nstr << gbl->tstep << std::flush;
        name = "coarsen" +nstr.str() +"_" +gbl->idprefix +".";
        nstr.str("");
        nstr << coarse_level << flush;
        fname = name +nstr.str();
        output(fname,mesh::grid);
        fname = name +nstr.str() + "_ft_to_cv";
        tgt.testconnect(fname,tgt.ccnnct,this);
        fname = name +nstr.str() + "_cv_to_ft";
        testconnect(fname,fcnnct,&tgt);
    }
    setinfo();
}

void mesh::mgconnect(mesh &tgt, Array<transfer,1> &cnnct) {   
    int i,bnum,v0;

    /* LOOP THROUGH VERTICES AND FIND SURROUNDING TRIANGLE */
    for(i=0;i<nvrtx;++i) {
        tgt.qtree.nearpt(vrtx(i).data(),v0);
        cnnct(i).tri = abs(tgt.findtri(vrtx(i),v0));
        tgt.getwgts(cnnct(i).wt);
    }
    
    //Array<boost::function0<void>,1> thread_func(nsbd);
    //boost::thread_group threads;
    
    /* REDO BOUNDARY SIDES TO DEAL WITH CURVATURE */
    for(bnum=0;bnum<nsbd;++bnum) {
        /* CHECK TO MAKE SURE THESE ARE THE SAME SIDES */
        if(sbdry(bnum)->idnum != tgt.sbdry(bnum)->idnum) {
            *gbl->log << "error: sides are not numbered the same" << std::endl;
            exit(1);
        }
        sbdry(bnum)->mgconnect(cnnct,tgt,bnum);
        //thread_func(bnum) = boost::bind(&side_bdry::mgconnect,sbdry(bnum),boost::ref(cnnct),boost::ref(tgt),bnum);
        //threads.create_thread(thread_func(bnum));
    }
    //threads.join_all();
    
}


/*	 THIS ROUTINE DETERMINES THE POSITION OF COARSE VERTICES  */
/*  TO TEST USING MULTI-GRID CONNECTION */
/* THE MULTIGRID CONNECTIONS */
 void mesh::testconnect(const std::string &fname,Array<transfer,1> &cnnct, mesh *cmesh) {
    int i,j,n,tind;
    Array<TinyVector<FLT,2>,1> work;

    work.resize(maxvst);
        
    if (cmesh != NULL) {
                
        /* LOOP THROUGH VERTICES TO TO CALCULATE POSITION OF COARSE VERTICES  */
        for(i=0;i<nvrtx;++i) {
            tind = cnnct(i).tri;
            
            for(n=0;n<ND;++n)
                work(i)(n) = 0.0;
            
            for(j=0;j<3;++j) {
                for(n=0;n<ND;++n)
                    work(i)(n) += cnnct(i).wt(j)*cmesh->vrtx(cmesh->td(tind).vrtx(j))(n);
            }
        }
        Array<TinyVector<FLT,2>,1> storevrtx(vrtx);
        vrtx.reference(work);
        output(fname, grid);
        vrtx.reference(storevrtx);
    }
        
    return;
}
    
