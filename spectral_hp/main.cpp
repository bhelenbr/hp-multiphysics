/*
 *  main.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "blocks.h"
#include<stdio.h>
#include<utilities.h>

extern FLT f1(int n, FLT x, FLT y);

int main() {
   
   int i,j;
   class spectral_hp x,y;
   class hpbasis base;
   class blocks myblock;
   class mesh m1,m2;

#ifdef SKIP
//   x.in_mesh("/private/Network/Servers/shelob.camp.clarkson.edu/home/helenbrk/codes/grids/KOVA/WKOVA/kova3",easymesh,10);
   m1.in_mesh("../../grids/TIM/tim",easymesh,10);
   m2.coarsen(m1);
   m2.bcinfo();
   m2.out_mesh("tim2",easymesh);
   exit(1);
   
   x.allocate(base);
   x.tobasis(&f1);
   exit(0);
   
   for(i=0;i<2;++i) {
      x.output("hope",tecplot);
      x.density(0.001,0.001,1.0);
      x.adapt(y,0.66);
      y.output("hope1",tecplot);
      x.output("new",tecplot);
   }

   return(0);
#endif
   
//   myblock.init(1, 1, 0, "../../grids/TIM/tim2",easymesh,5);

//   myblock.init(1, 1, 1, "../../grids/KOVA/kova16",easymesh,5);
//   9  7.569496e-05  3.783905e-05  1.604564e-04

   myblock.init(1,1,0,"../../grids/WAVE/wave3",easymesh,5);
   myblock.tadvance();
   for(i=0;i<1;++i) {
      myblock.cycle(1);
      printf("%d ",i);
      myblock.print_maxres();
      printf("\n");
   }
   
   myblock.output("hope1",tecplot);

   return(0);
}



