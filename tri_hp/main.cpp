/*
 *  main.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#ifdef MPISRC
#include <mpi.h>
#endif

#include <blocks.h>
#include <signal.h>
#ifdef petsc
#include <petscksp.h>
static char help[] = "How am I supposed to know???\n\n";
#endif

#include <parseargs.h>
static GBool Debugger = gFalse;
GBool printHelp = gFalse;

static ArgDesc argDesc[] = {
	{"-h",        argFlag,      &printHelp,      0,
		"print usage information"},
	{"-help",    argFlag,      &printHelp,      0,
		"print usage information"},
	{"-stop_for_debugger"  ,argFlag,     &Debugger,            0,
		"Stop for Debugger"},
	{NULL}
};

#include <thread>
#include <unistd.h>


void ctrlc(int signal);

int main(int argc, char **argv) {
	struct sigaction action;
	struct sigaction o_action;

#ifdef MPISRC
	int myid;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	
	std::cout << "Running processes " << getpid() << std::endl;
	// system("sleep 20");
#endif
#ifdef PTH
	// For debugging put interrupt here
	// Then on interrupt in console window type:
    // for gdb console: handle SIGUSR1 nostop print pass
    // for lldb: pro hand -p true -s false SIGUSR1
	// Then hit continue
	int rc = pth_init();
	if (!rc) {
		std::cerr << "couldn't start pth environment" << std::endl;
	}
#endif
#ifdef petsc
	PetscErrorCode err = PetscInitialize(&argc,&argv,(char *)0,help);
	CHKERRABORT(MPI_COMM_WORLD,err);
#endif
	
	// parse args
	GBool ok = parseArgs(argDesc, &argc, argv);
	if (!ok || printHelp) {
		fprintf(stderr, "tri_hp ");
		printUsage("tri_hp", "<inputfile>]", argDesc);
		sim::abort(__LINE__,__FILE__,&std::cerr);
	}
	
	
#if (defined(MPISRC) && !defined(petsc))
	if (Debugger) {
		int size;
		/*
		 we have to make sure that all processors have opened
		 connections to all other processors, otherwise once the
		 debugger has stated it is likely to receive a SIGUSR1
		 and kill the program.
		 */
		int ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);
		if (ierr != MPI_SUCCESS) {
			exit(1);
		}
		if (size > 2) {
			int dummy = 0;
			MPI_Status status;
			for (int i=0; i<size; i++) {
				if (myid != i) {
					ierr = MPI_Send(&dummy,1,MPI_INT,i,109,MPI_COMM_WORLD);
				}
			}
			for (int i=0; i<size; i++) {
				if (myid != i) {
					ierr = MPI_Recv(&dummy,1,MPI_INT,i,109,MPI_COMM_WORLD,&status);				}
			}
		}
		std::cout << "Waiting for debugger.  Process id is " << getpid() << std::endl;
#if __cplusplus >= 199711L
        std::this_thread::sleep_for(std::chrono::milliseconds(10000));
#endif
	}
#endif


/*	INTERRUPT HANDLER FOR GRACEFUL EXIT ON CTRL-C    	*/	
	action.sa_handler = ctrlc;
	sigemptyset(&action.sa_mask);
	action.sa_flags = 0;
	if (sigaction(SIGTERM,&action,&o_action))
		std::cerr << "interrupt handler failed" << std::endl;

	/* NORMAL SIMULATION */
	if (argc < 2) {
		std::cerr << "# Need to specify input file" << std::endl;
		sim::abort(__LINE__,__FILE__,&std::cerr);
	}
	sim::blks.go(argv[1]);

#ifdef petsc
	PetscFinalize();
#endif
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
//#endif
	
	return 0;
}

void ctrlc(int signal)
{  
	sim::abort(__LINE__,__FILE__,&std::cerr);

}

