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

#ifdef SKIP
   FLT u[3];
 
   base.initialize(2);

//   x.in_mesh("/private/Network/Servers/shelob.camp.clarkson.edu/home/helenbrk/codes/grids/KOVA/WKOVA/kova3",easymesh,10);
   x.in_mesh("../../grids/TIM/tim",easymesh,10);
   x.allocate(base);
   x.tobasis(&f1);
   x.output("hope",tecplot);
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
   
   myblock.init(1, 3, 1, "../../grids/TIM/tim",easymesh,10);
   myblock.output("hope0",tecplot);

   myblock.tadvance();
   for(i=0;i<1;++i) {
      myblock.cycle(1);
      printf("%d ",i);
      myblock.maxres();
      printf("\n");
   }
   
   myblock.output("hope1",tecplot);

   return(0);
}



