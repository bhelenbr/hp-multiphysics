/*
 *  main.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#include <blocks.h>
#include <signal.h>
#ifdef PV3
#include <pV3.h>
#endif
#ifdef MPISRC
#include <mpi.h>
#endif

void ctrlc(int signal);

int main(int argc, char **argv) {
   struct sigaction action;
	struct sigaction o_action;

/*	INTERRUPT HANDLER FOR GRACEFUL EXIT ON CTRL-C   	*/	
	action.sa_handler = ctrlc;
	sigemptyset(&action.sa_mask);
	action.sa_flags = 0;
	if (sigaction(SIGTERM,&action,&o_action)) printf("interrupt handler failed\n");
   

#ifdef MPISRC
   int myid;
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
   
   /* NORMAL SIMULATION */
   sim::blks.init(argv[1]);
   sim::blks.go();

#ifdef MPISRC
   MPI_Finalize();
#endif

//
//#ifdef DEBUG
//   extern FLT f1(int n, FLT x, FLT y);
//   
//   myblock.init(argv[1]);
//   myblock.minvrt_test(10,f1);
//   myblock.output(-2,text);
//   myblock.go();
//#endif
//
//
//#ifdef FINDMAX
//   hpbasis b;
//   hp_mgrid hp;
//   b.initialize(4,5);
//   hp.in_mesh(argv[1],grid,1.0);
//   hp.tri_hp::allocate(&b);
//   hp.input(argv[1],tecplot);
//   hp.findintercept(CURV_MASK,&ydist);
//   //hp.findmaxy(CURV_MASK);
//   // hp.findmaxx(CURV_MASK);
//   //FLT avg[5];
//   //hp.integrated_averages(avg);
//   //printf("%e %e %e %e %e\n",avg[0],avg[1],avg[2],avg[3],avg[4]);
//#endif

}

void ctrlc(int signal)
{  
   /* THIS ALLOWS FILES TO CLOSE */
   /* AND OUTPUTS SOLUTION AT TIME OF INTERRUPT */
   *sim::log << "# exiting gracefully" << std::endl;
   sim::blks.output("interrupt",block::restart);
   sim::blks.output("interrupt",block::display);
   exit(1);
}

