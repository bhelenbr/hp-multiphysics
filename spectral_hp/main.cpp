/*
 *  main.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */


#include "blocks.h"
#include"pV3.h"
#include<stdio.h>
#include<utilities.h>
#include<time.h>

#ifdef PV3
extern "C" int MAINPROG(int argc, char **argv);
int MAINPROG(int argc, char**argv) {
#else
int main(int argc, char **argv) {
#endif
   class blocks myblock;
   clock_t cpu_time;
   class hp_mgrid test;
   class hpbasis w;
   struct hp_mgrid_glbls gbl;
      
   myblock.init(argv[1]);

   /* START CLOCK TIMER */
   clock();
   myblock.go();
   cpu_time = clock();
   printf("that took %ld cpu time\n",cpu_time);
   
   return(0);


   /* FOR TEST ADAPTATION */
   myblock.init(argv[1]);
   myblock.output(10,tecplot);
   
   return(0);
   
   /*	TO OUTPUT TRUNCATION ERROR */
   /* 
   test.mesh::in_mesh("mesh530.0",grid);
   w.initialize(4);
   test.mesh::setbcinfo();
   test.spectral_hp::allocate(w);
   gbl.rho = 1.0;
   test.allocate(0,&gbl);
   test.spectral_hp::input("data530.0",text);
   test.outlength("terror",tecplot);
   return(0);
   
   
   

}



