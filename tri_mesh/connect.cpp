#include "tri_mesh.h"
#include <cfloat>
#include <iostream>
#include <stdlib.h>
#include <utilities.h>
#include <string.h>
#include "block.h"
// #include <boost/bind.hpp>
// #include <boost/thread.hpp>

void tri_mesh::connect(multigrid_interface& in) {
    tri_mesh &tgt = dynamic_cast<tri_mesh&>(in);
    
    /* SET UP MULTIGRID STUFF */
    fine = &in;
    tgt.coarse = this;    
    if (fcnnct.ubound(firstDim) < maxpst-1) fcnnct.resize(maxpst);
    if (tgt.ccnnct.ubound(firstDim) < tgt.maxpst-1) tgt.ccnnct.resize(tgt.maxpst);
    
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
        output(fname,tri_mesh::grid);
        fname = name +nstr.str() + "_ft_to_cv";
        tgt.testconnect(fname,tgt.ccnnct,this);
        fname = name +nstr.str() + "_cv_to_ft";
        testconnect(fname,fcnnct,&tgt);
    }
    setinfo();
}

void tri_mesh::mgconnect(tri_mesh &tgt, Array<transfer,1> &cnnct) {   
    int i,bnum,p0;

    /* LOOP THROUGH POINTS AND FIND SURROUNDING TRIANGLE */
    for(i=0;i<npnt;++i) {
        tgt.qtree.nearpt(pnts(i).data(),p0);
        cnnct(i).tri = abs(tgt.findtri(pnts(i),p0));
        tgt.getwgts(cnnct(i).wt);
    }
    
    //Array<boost::function0<void>,1> thread_func(nebd);
    //boost::thread_group threads;
    
    /* REDO BOUNDARY SIDES TO DEAL WITH CURVATURE */
    for(bnum=0;bnum<nebd;++bnum) {
        /* CHECK TO MAKE SURE THESE ARE THE SAME SIDES */
        if(ebdry(bnum)->idnum != tgt.ebdry(bnum)->idnum) {
            *gbl->log << "error: sides are not numbered the same" << std::endl;
            exit(1);
        }
        ebdry(bnum)->mgconnect(cnnct,tgt,bnum);
        //thread_func(bnum) = boost::bind(&edge_bdry::mgconnect,ebdry(bnum),boost::ref(cnnct),boost::ref(tgt),bnum);
        //threads.create_thread(thread_func(bnum));
    }
    //threads.join_all();
    
}


/*	 THIS ROUTINE DETERMINES THE POSITION OF COARSE POINTS  */
/*  TO TEST USING MULTI-GRID CONNECTION */
/* THE MULTIGRID CONNECTIONS */
 void tri_mesh::testconnect(const std::string &fname,Array<transfer,1> &cnnct, tri_mesh *cmesh) {
    int i,j,n,tind;
    Array<TinyVector<FLT,2>,1> work;

    work.resize(maxpst);
        
    if (cmesh != NULL) {
                
        /* LOOP THROUGH POINTS TO TO CALCULATE POSITION OF COARSE POINTS  */
        for(i=0;i<npnt;++i) {
            tind = cnnct(i).tri;
            
            for(n=0;n<ND;++n)
                work(i)(n) = 0.0;
            
            for(j=0;j<3;++j) {
                for(n=0;n<ND;++n)
                    work(i)(n) += cnnct(i).wt(j)*cmesh->pnts(cmesh->tri(tind).pnt(j))(n);
            }
        }
        Array<TinyVector<FLT,2>,1> storevrtx(pnts);
        pnts.reference(work);
        output(fname, grid);
        pnts.reference(storevrtx);
    }
        
    return;
}
    
