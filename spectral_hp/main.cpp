/*
 *  main.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#include "defines.h"
#include "blocks.h"
#include "pV3.h"
#include <stdio.h>
#include <utilities.h>
#include <string.h>
#include <time.h>

#ifdef PV3
extern "C" int MAINPROG(int argc, char **argv);
int MAINPROG(int argc, char**argv) {
#else
int main(int argc, char **argv) {
#endif
   clock_t cpu_time;
   class blocks myblock;

#ifdef SIMULATION      
   /* NORMAL SIMULATION */
   myblock.init(argv[1]);
   myblock.go();
   cpu_time = clock();
   printf("#that took %ld cpu time\n",cpu_time);
   return(0);
#endif

#ifdef DEBUG
   extern FLT f1(int n, FLT x, FLT y);

   myblock.init(argv[1]);
   myblock.minvrt_test(10,f1);
   myblock.output(-2,text);
   myblock.go();
#endif

#ifdef PV3VIEWER   
   /* PV3 STATIC VIEWER */
   myblock.init(argv[1],0);
   return(0);
#endif

#ifdef FINDMAX
   hpbasis b;
   spectral_hp hp;
   b.initialize(4,5);
   hp.in_mesh(argv[1],grid,1.0);
   hp.allocate(b);
   hp.input(argv[1],tecplot);
   hp.output("test",tecplot);
   hp.output("test",text);
   hp.findmaxy(IFCE_MASK);
#endif
   
   
}



