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
#ifdef MPISRC
#include <mpi.h>
#endif

void ctrlc(int signal);

int main(int argc, char **argv) {
    struct sigaction action;
	struct sigaction o_action;

#ifdef MPISRC
	int myid;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
#ifdef PTH
	// For debugging put interrupt here
	// On interrupt type this into gdb console: handle SIGUSR1 nostop print pass
	// Then continue
	int rc = pth_init();
	if (!rc) {
		std::cerr << "couldn't start pth environment" << std::endl;
	}
#endif

/*	INTERRUPT HANDLER FOR GRACEFUL EXIT ON CTRL-C    	*/	
	action.sa_handler = ctrlc;
	sigemptyset(&action.sa_mask);
	action.sa_flags = 0;
	if (sigaction(SIGTERM,&action,&o_action)) printf("interrupt handler failed\n");

	/* NORMAL SIMULATION */
	if (argc < 2) {
		std::cerr << "# Need to specify input file" << std::endl;
		exit(1);
	}
	sim::blks.go(argv[1]);

#ifdef PTH
	pth_exit(NULL);
#endif    
#ifdef MPISRC
	MPI_Finalize();
#endif

//
//#ifdef DEBUG
//    extern FLT f1(int n, FLT x, FLT y);
//    
//    myblock.init(argv[1]);
//    myblock.minvrt_test(10,f1);
//    myblock.output(-2,text);
//    myblock.go();
//#endif
//
//
//#ifdef FINDMAX
//    tri_basis b;
//    hp_mgrid hp;
//    b.initialize(4,5);
//    hp.in_mesh(argv[1],grid,1.0);
//    hp.tri_hp::allocate(&b);
//    hp.input(argv[1],tecplot);
//    hp.findintercept(CURV_MASK,&ydist);
//    //hp.findmaxy(CURV_MASK);
//    // hp.findmaxx(CURV_MASK);
//    //FLT avg[5];
//    //hp.integrated_averages(avg);
//    //printf("%e %e %e %e %e\n",avg[0],avg[1],avg[2],avg[3],avg[4]);
//#endif

}

void ctrlc(int signal)
{  
    /* THIS ALLOWS FILES TO CLOSE */
    /* AND OUTPUTS SOLUTION AT TIME OF INTERRUPT */
//    sim::blks.output("interrupt",block::restart);
//    sim::blks.output("interrupt",block::display);
#ifdef PTH
	pth_exit(NULL);
#endif    
#ifdef MPISRC
	MPI_Finalize();
#endif
	exit(1);
}

