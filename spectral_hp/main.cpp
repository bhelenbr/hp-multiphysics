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
extern FLT rhotemporary;

int main(int argc, char **argv) {
   class blocks myblock;

#ifdef SKIP
   class mesh test1, test2;
   
   test1.in_mesh("error1",easymesh);
   test2.coarsen(test1);
   test2.setbcinfo();
   test2.out_mesh("error2");
   
   exit(1);
#endif

   myblock.init(argv[1]);
   myblock.go();
   
   return(0);

   myblock.init(argv[1]);
   myblock.tadvance();
   myblock.output(1,tecplot);
   myblock.adaptation();
   myblock.tadvance();
   myblock.output(2,tecplot);

   
   
   return(0);
}



