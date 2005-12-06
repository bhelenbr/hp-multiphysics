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
   int i,n,t;
   
   /* COPY MESH INFORMATION */
   mesh::copy(tgt);

   for(t=0;t<sim::nadapt+1;++t) {
      ugbd(t).v(Range(0,nvrtx-1),Range::all()) = tgt.ugbd(t).v(Range(0,nvrtx-1),Range::all());
      ugbd(t).s(Range(0,nside-1),Range::all(),Range::all()) = tgt.ugbd(t).s(Range(0,nside-1),Range::all(),Range::all());
      ugbd(t).i(Range(0,ntri-1),Range::all(),Range::all()) = tgt.ugbd(t).i(Range(0,ntri-1),Range::all(),Range::all());
      
      for(i=0;i<nvrtx;++i)
         for(n=0;n<ND;++n)
            vrtxbd(t)(i)(n) = tgt.vrtxbd(t)(i)(n);
   }

   for(i=0;i<nsbd;++i)
      hp_sbdry(i)->copy_data(*tgt.hp_sbdry(i));
      
   return;
}
         
   


