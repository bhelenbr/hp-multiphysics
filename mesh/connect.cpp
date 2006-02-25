#include "mesh.h"
#include "boundary.h"
#include <cfloat>
#include <iostream>
#include <stdlib.h>
#include <utilities.h>
#include <string.h>
#include "block.h"

block::ctrl mesh::mgconnect(block::ctrl ctrl_message, Array<transfer,1> &cnnct, const class mesh& tgt) {
   int i,bnum,v0;
   block::ctrl state = block::stop;
   
   if (ctrl_message == block::begin) excpt = 0;
   else ++excpt;
   
   switch(excpt) {
      case(0):
         /* LOOP THROUGH VERTICES AND FIND SURROUNDING TRIANGLE */
         for(i=0;i<nvrtx;++i) {
            tgt.qtree.nearpt(vrtx(i).data(),v0);
            cnnct(i).tri = abs(tgt.findtri(vrtx(i),v0));
            tgt.getwgts(cnnct(i).wt);
         }
      default:
         /* REDO BOUNDARY SIDES TO DEAL WITH CURVATURE */
         for(bnum=0;bnum<nsbd;++bnum) {
            /* CHECK TO MAKE SURE THESE ARE THE SAME SIDES */
            if(sbdry(bnum)->idnum != tgt.sbdry(bnum)->idnum) {
               *sim::log << "error: sides are not numbered the same" << std::endl;
               exit(1);
            }
            state &= sbdry(bnum)->mgconnect(ctrl_message,cnnct,tgt,bnum);
         }
   }
   return(state);
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
   
