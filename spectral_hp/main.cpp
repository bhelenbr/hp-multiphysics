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
#include <stdio.h>
#include <utilities.h>
#include <string.h>
#include <time.h>
#include <signal.h>
#ifdef PV3
#include "pV3.h"
#endif

void ctrlc(int signal);
static class blocks myblock;

#ifdef PV3
extern "C" int MAINPROG(int argc, char **argv);
int MAINPROG(int argc, char**argv) {
#else
int main(int argc, char **argv) {
#endif
   clock_t cpu_time;
   struct sigaction action;
	struct sigaction o_action;

/*	INTERRUPT HANDLER FOR GRACEFUL EXIT ON CTRL-C   	*/	
	action.sa_handler = ctrlc;
	sigemptyset(&action.sa_mask);
	action.sa_flags = 0;
	if (sigaction(SIGTERM,&action,&o_action)) printf("interrupt handler failed\n");

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
   hp_mgrid hp;
   b.initialize(4,5);
   hp.in_mesh(argv[1],grid,1.0);
   hp.spectral_hp::allocate(b);
   hp.input(argv[1],tecplot);
//   hp.findmaxy(CURV_MASK);
   hp.findmaxx(CURV_MASK);
   FLT avg[5];
   hp.integrated_averages(avg);
   printf("%e %e %e %e %e\n",avg[0],avg[1],avg[2],avg[3],avg[4]);
#endif

}

void ctrlc(int signal)
{  
   /* THIS ALLOWS FILES TO CLOSE */
   /* AND OUTPUTS SOLUTION AT TIME OF INTERRUPT */
   printf("# exiting gracefully\n");
   myblock.loadbasis();
   myblock.output(-1,text,1);
   myblock.output(-1,tecplot);

   exit(1);
}

