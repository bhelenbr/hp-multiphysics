/* THIS IS A MULTIBLOCK MESH */
/* CAN DRIVE MULTIGRID OR STANDARD ITERATION */

#ifndef _blocks_h_
#define _blocks_h_

#include "block.h"
#include <fstream>
#include <blitz/array.h>
#include <input_map.h>
#include <utilities.h>

#define PTH
// #define BOOST
#ifdef PTH
#include <pth.h>
#endif

#ifdef MPISRC
#include <mpi.h>
#endif

#ifdef BOOST
#include <boost/thread.hpp>
#include <boost/bind.hpp>
#endif

#define DIRK 4
#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

/* THE PURPOSE OF THE BLOCKS OBJECT IS TO INITIATE THREADS AND COORDINATE THREADED MESSAGE PASSING */
class blocks {
    public:
        int nproc; /**< Number of processors in simulation */
        int myid;  /**< Processor number of myself for MPI */
        int nblock; /**< Number of blocks in simulation */
        int myblock; /**< Number of blocks on this processor */
        Array<block *,1> blk; /**< Array containing pointers to blocks */
    private:
        block* getnewblock(int idnum, input_map& input) {
            return new block(idnum);  /* Only 1 kind of block for now */
        }
        int max_mycommbdry; /**< Maximum number of communication boundaries for any block */
        int blkdig, bdrydig; /**< Number of digits required to make communication tags */

        /** @name allreduce buffer pointer storage
         *  used in allreduce for block communication
         */
        //@{
		class all_reduce_data {
			public: 
				int nentries, buf_cnt;
				blitz::Array<void *,1> sndbufs, rcvbufs;
#ifdef MPISRC
				MPI_Comm comm;
#endif
				all_reduce_data() : nentries(1), buf_cnt(0)
#ifdef MPISRC 
			, comm(MPI_COMM_WORLD) 
#endif
			{}
		};
		std::map<int,all_reduce_data> group_data; /**< Associaciates group numbers to their data */ 
        //@}
        
         /** @name variables needed for thread message passing
         */
        //@{
        std::map<int,bool> message_list; /**< keeps track of whether communication messages have been received */
#if defined(PTH)
        blitz::Array<pth_t,1> threads;
        pth_mutex_t list_mutex, allreduce_mutex;
        pth_cond_t list_change, allreduce_change;
        pth_attr_t attr;
#elif defined(BOOST)
        blitz::Array<boost::function0<void>,1> thread_func;
        boost::thread_group threads;
        boost::mutex list_mutex, allreduce_mutex;
        boost::condition list_change, allreduce_change;
#endif
        //@}

    public:
        /** Initialize multiblock/mgrid mesh */
        blocks() : nproc(1), myid(0) {}
        
        /** Initializes blocks using data from map */
        void go(input_map& input);
        
        /** Loads map from file then initialize */
        void go(const std::string &infile, const std::string &outfile = std::string());
                
        /** Function to allow blocks to perform global sums, averages, etc.. */
        enum operations {sum,max}; /** Only summing and max for now */
        enum msg_type {flt_msg, int_msg};
        void allreduce(void *sendbuf, void *recvbuf, int count, msg_type datatype, operations op, int group = 0);
        
        /** Functions for thread communication */
#if defined(PTH)
        pth_mutex_t data_mutex;
#elif defined(BOOST)
        boost::mutex data_mutex;
#endif

        void setdigits(int ncomm_bdry) {
#if defined(PTH)
            pth_mutex_acquire(&data_mutex,false,NULL);
#elif defined(BOOST)
            boost::mutex::scoped_lock lock(data_mutex);
#endif
            max_mycommbdry = ncomm_bdry;
            blkdig = 1;
            while ((nblock>>blkdig) > 0) ++blkdig;
            bdrydig = 1;
            while ((max_mycommbdry>>bdrydig) > 0) ++bdrydig;
            if (2*(blkdig+bdrydig)+2 > 32) std::cerr << "can't guarantee unique tags\n"; 
#if defined(PTH)
            pth_mutex_release(&data_mutex);
#endif
            return;
        }

        int tagid(int vsf,int b1,int b2,int n1,int n2) {
            int tag,lefttag,righttag;        
            lefttag = (b1<<(bdrydig)) +n1;
            righttag = (b2<<(bdrydig)) +n2;
            tag = (vsf<<(2*(blkdig+bdrydig))) +(lefttag<<(blkdig+bdrydig)) +righttag;  
            return(tag);
        }

        void waitforslot(int msgid, bool set) {
            std::map<int,bool>::iterator mi;
#if defined(PTH)
            pth_mutex_acquire(&list_mutex,false,NULL);
#elif defined(BOOST)
            boost::mutex::scoped_lock lock(list_mutex);
#endif
            mi = message_list.find(msgid);
            if (mi == message_list.end()) message_list[msgid] = false;
            
            while(message_list[msgid] != set) {
#if defined(PTH)
                pth_cond_await(&list_change, &list_mutex, NULL);
#elif defined(BOOST)
                list_change.wait(lock);
#endif
            }
#if defined(PTH)
            pth_mutex_release(&list_mutex);
#endif
        }
        
        void notify_change(int msgid, bool set) {
#if defined(PTH)
            pth_mutex_acquire(&list_mutex,false,NULL);
#elif defined(BOOST)
            boost::mutex::scoped_lock lock(list_mutex);
#endif
            message_list[msgid] = set;
#if defined(PTH)
            pth_cond_notify(&list_change, true);
            pth_mutex_release(&list_mutex);
#elif defined(BOOST)
            list_change.notify_all();
#endif
        }
};


/** \brief Global variables for simulation
 *
 * This namespace contains global variables needed by all blocks of the simulation 
 */
namespace sim {
    extern blocks blks; /**< Contains all blocks for this processor */
    extern int nproc;
    extern int myid;
    
    /** Time stepping for simulation */
#ifdef BACKDIFF
    /** @name backdiff Backwards difference constants
    *  These are constant for backwards difference timestepping
    */
    //@{
    const int nhist = BACKDIFF; /**< number of backwards difference steps */
    const int nadapt = BACKDIFF; /**< number of backwards difference steps that require adaptation */
    const int stepsolves = 1;
    //@}
#endif
#ifdef DIRK
    /** @name DIRK3 scheme constants
     *  Constants used to define the DIRK3 scheme
    */
    //@{
    const FLT GRK3 = 0.43586652150845899941601945;
    const FLT C2RK3 = (2.-9*GRK3+6.*GRK3*GRK3)/(3*(1-4*GRK3+2*GRK3*GRK3));
    const FLT B2RK3 = -3*(1-4*GRK3+2*GRK3*GRK3)*(1-4*GRK3+2*GRK3*GRK3)/(4*(-1+6*GRK3-9*GRK3*GRK3+3*GRK3*GRK3*GRK3));
    //@}

    /** @name DIRK4 scheme constants
     *  Constants used to define the DIRK4 scheme
    */
    //@{
    // const FLT GRK4 0.5  // FOR ERROR PREDICTION
    const FLT GRK4 = 0.43586652150845899941601945; // FOR L-STABILITY
    const FLT C3RK4 = 2*GRK4*(GRK4-1./4.)*(GRK4-1.)/((GRK4-0.5)*(GRK4-0.5)-1./12.);
    const FLT A32RK4 = (C3RK4*C3RK4-2*C3RK4*GRK4)/(4.*GRK4);
    const FLT B1RK4 = (1./6. +GRK4*GRK4-GRK4*GRK4*C3RK4+3./2.*GRK4*C3RK4-GRK4-1./4.*C3RK4)/(GRK4*C3RK4);
    const FLT B2RK4 = (1./3.-GRK4-1./2.*C3RK4+GRK4*C3RK4)/(2.*GRK4*(2.*GRK4-C3RK4));
    const FLT B3RK4 =( 1./3.-2*GRK4+2*GRK4*GRK4)/(C3RK4*(C3RK4-2*GRK4));
    const int nadapt = 2; /**< number of time steps that require adaptation */
    const int nhist = 4; /**< number of time steps to be stored during a time step */
    const int stepsolves = 3; /**< Number of DIRK implicit solutions required */
    //@}
#endif

    /** @name Multistage Explicit Scheme
     *  Constants for multistage iteration scheme 
     */
    //@{
    const int NSTAGE = 5; /**< Number of stages */
    const FLT alpha[NSTAGE+1] = {0.25, 1./6., .375, .5, 1.0, 1.0}; /**< Multistage time step constants (imaginary) */
    const FLT beta[NSTAGE+1] = {1.0, 0.0, 5./9., 0.0, 4./9., 1.0}; /**< Multistage time step constants (real) */
    //@}
    
}

#endif




