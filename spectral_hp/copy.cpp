/*
 *  copy.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "utilities.h"
#include "hp_boundary.h"

 void tri_hp::copy_data(const tri_hp& tgt) {
   int i,n,t,state;
   
   /* COPY MESH INFORMATION */
   state = mesh::initialized;
   mesh::copy(tgt);

   if (!state) {
      map<std::string,std::string> input;
      std::string prefix = "blk";
      std::string keyword = prefix + ".adapt_storage";
      input[keyword] = "1";
      
      ostringstream myStream;
      myStream << log2p << flush;
      input[prefix+".log2p"] = myStream.str();
      init(input, prefix, tgt.hp_gbl);
   }
      
   for(t=0;t<sim::nhist+1;++t) {
      ugbd(t).v(Range(0,nvrtx),Range::all()) = tgt.ugbd(t).v(Range(0,nvrtx),Range::all());
      ugbd(t).s(Range(0,nside),Range::all(),Range::all()) = tgt.ugbd(t).s(Range(0,nside),Range::all(),Range::all());
      ugbd(t).i(Range(0,ntri),Range::all(),Range::all()) = tgt.ugbd(t).i(Range(0,ntri),Range::all(),Range::all());
      
      for(i=0;i<nvrtx;++i)
         for(n=0;n<ND;++n)
            vrtxbd(t)(i)(n) = tgt.vrtxbd(t)(i)(n);
   }

   for(i=0;i<nsbd;++i)
      hp_sbdry(i)->copy_data(tgt.hp_sbdry(i));
      
   return;
}
         
   


