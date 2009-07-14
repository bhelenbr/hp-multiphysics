#include "tet_mesh.h"
#include <cfloat>
#include <iostream>
#include <stdlib.h>
#include <utilities.h>
#include <string.h>
#include "block.h"
#ifdef USING_MADLIB
#include "MAdLibInterface.h"
#endif

// #include <boost/bind.hpp>
// #include <boost/thread.hpp>

void tet_mesh::coarsen(FLT factor, const tet_mesh& tgt) {
	copy(tgt);
#ifdef USING_MADLIB
	MAdLibInterface::coarsenMesh(factor,this);
#endif
	return;
}

void tet_mesh::connect(multigrid_interface& in) {
	tet_mesh &tgt = dynamic_cast<tet_mesh&>(in);
	
	/* SET UP MULTIGRID STUFF */
	fine = &in;
	tgt.coarse = this;    
	if (fcnnct.ubound(firstDim) < maxvst-1) fcnnct.resize(maxvst);
	if (tgt.ccnnct.ubound(firstDim) < tgt.maxvst-1) tgt.ccnnct.resize(tgt.maxvst);
	
	// copy(tgt);
	coarsen(2.0,tgt);
	mgconnect(tgt,fcnnct);
	tgt.mgconnect(*this,tgt.ccnnct);

	/* THIS IS FOR DIAGNOSIS OF MULTI-GRID FAILURES */
	if (gbl->adapt_output) {
		checkintegrity();
		tgt.checkintegrity();
		std::string name, fname;
		std::string adapt_file;
		std::ostringstream nstr;
		nstr.str("");
		nstr << gbl->tstep << std::flush;
		name = "coarsen" +nstr.str() +"_" +gbl->idprefix +".";
		nstr.str("");
		nstr << coarse_level << flush;
		fname = name +nstr.str();
		output(fname,tet_mesh::grid);
		fname = name +nstr.str() + "_ft_to_cv";
		tgt.testconnect(fname,tgt.ccnnct,this);
		fname = name +nstr.str() + "_cv_to_ft";
		testconnect(fname,fcnnct,&tgt);
	}
	setinfo();
}

void tet_mesh::mgconnect(tet_mesh &tgt, Array<transfer,1> &cnnct) {   
	int i,bnum,p0;

	/* LOOP THROUGH POINTS AND FIND SURROUNDING TRIANGLE */
	for(i=0;i<npnt;++i) {
		tgt.otree.nearpt(pnts(i).data(),p0);
		tgt.findtet(pnts(i),p0,cnnct(i).tet);
		tgt.getwgts(cnnct(i).wt);
	}
	
	//Array<boost::function0<void>,1> thread_func(nebd);
	//boost::thread_group threads;
	
	/* REDO BOUNDARY SIDES TO DEAL WITH CURVATURE */
	for(bnum=0;bnum<nfbd;++bnum) {
		/* CHECK TO MAKE SURE THESE ARE THE SAME SIDES */
		if(fbdry(bnum)->idnum != tgt.fbdry(bnum)->idnum) {
			*gbl->log << "error: sides are not numbered the same" << std::endl;
			exit(1);
		}
		fbdry(bnum)->mgconnect(cnnct,tgt,bnum);
		//thread_func(bnum) = boost::bind(&edge_bdry::mgconnect,ebdry(bnum),boost::ref(cnnct),boost::ref(tgt),bnum);
		//threads.create_thread(thread_func(bnum));
	}
	//threads.join_all();
	
	
	for(bnum=0;bnum<nebd;++bnum) {
		/* CHECK TO MAKE SURE THESE ARE THE SAME SIDES */
		if(ebdry(bnum)->idnum != tgt.ebdry(bnum)->idnum) {
			*gbl->log << "error: sides are not numbered the same" << std::endl;
			exit(1);
		}
		ebdry(bnum)->mgconnect(cnnct,tgt,bnum);
	}	
}


/*	 THIS ROUTINE DETERMINES THE POSITION OF COARSE POINTS  */
/*  TO TEST USING MULTI-GRID CONNECTION */
/* THE MULTIGRID CONNECTIONS */
 void tet_mesh::testconnect(const std::string &fname,Array<transfer,1> &cnnct, tet_mesh *cmesh) {
	int i,j,n,ttind;
	Array<TinyVector<FLT,ND>,1> work;

	work.resize(maxvst);
		
	if (cmesh != NULL) {
				
		/* LOOP THROUGH POINTS TO TO CALCULATE POSITION OF COARSE POINTS  */
		for(i=0;i<npnt;++i) {
			ttind = cnnct(i).tet;
			
			for(n=0;n<ND;++n)
				work(i)(n) = 0.0;
			
			for(j=0;j<4;++j) {
				for(n=0;n<ND;++n)
					work(i)(n) += cnnct(i).wt(j)*cmesh->pnts(cmesh->tet(ttind).pnt(j))(n);
			}
		}
		Array<TinyVector<FLT,ND>,1> storevrtx(pnts);
		pnts.reference(work);
		output(fname, grid);
		pnts.reference(storevrtx);
	}
		
	return;
}
	
