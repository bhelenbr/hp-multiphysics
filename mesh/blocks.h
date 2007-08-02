/* THIS IS A MULTIBLOCK MESH */
/* CAN DRIVE MULTIGRID OR STANDARD ITERATION */

#ifndef _blocks_h_
#define _blocks_h_

#include "block.h"
#include <fstream>
#include <blitz/array.h>
#include <input_map.h>
#include <utilities.h>

// #define PTH
// #define BOOST


#ifdef PTH
#include <pth.h>
/* COMPLICATED WAY TO GET CORRECT FORM FOR FUNCTION */
template<class T, void(T::*mem_fn)()> void* thread_proxy(void* p) {
   (static_cast<T*>(p)->*mem_fn)();
   return 0;
}

extern int pth_int1, pth_int2;
template<class T, void(T::*mem_fn)(int)> void* thread_proxy(void* p) {
   (static_cast<T*>(p)->*mem_fn)(pth_int1);
   return 0;
}

template<class T, void(T::*mem_fn)(int,int)> void* thread_proxy(void* p) {
   (static_cast<T*>(p)->*mem_fn)(pth_int1,pth_int2);
   return 0;
}

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

class blocks {
    protected:
        int nproc; /**< Number of processors in simulation */
        int myid;  /**< Processor number of myself for MPI */
        int nblock; /**< Number of blocks on this processor */

        /** @name Multigrid parameters 
         *  constants defining multigrid iteratoin 
        */
        //@{
        int mglvls; /**< Total number of levels of multigrid */
        int ngrid; /**< Number of grids (could be more or less than mglvls) */
        int extra_finest_levels; /**< Number of extra levels to included on finest grid */
        int extra_coarsest_levels; /**< Number of extra levels to be include on coarsest grid */
        int ncycle; /**< Number of iterations per timestep */
        FLT absolute_tolerance; /**< Absolute error tolerance for iterative loop (negative means don't use) */
        FLT relative_tolerance; /**< Relative error tolerance for iterative loop (negative means don't use) */
        /** This for running a two-level iteration on performing extra iterations on coarsest mesh **/
        int error_control_level; /**< Level of multigrid to repeat iterations until convergence (negative means don't use) */
        FLT error_control_tolerance; /**< Relative error tolerance for error_control_level */
        /** Number of cycles between re-evaluation of preconditioner.
         *  negative mean reevaluate for both refinement and coarsening sweeps 
         */
        int prcndtn_intrvl; 
        int itercrsn; /**< Number of iterations to performing on coarsening */
        int iterrfne; /**< Number of iteration to perform on moving to finer mesh */
        int vw; /**< V-cycle = 1, W-cycle = 2 */
        //@}
        int ntstep;  /**< Number of time steps to perform */
        int nstart; /**< Starting step (for restart from file */
        block **blk; /**< Array containing pointers to blocks */
        
        /** @name Output parameters
         *  constants determining when and what to output 
         */
        //@{
        int out_intrvl; /**< Number of time-steps between data outputs */
        int rstrt_intrvl; /**< Number of output intervals between restart files */
        bool debug_output; /**< Output file every iteration */
        //@}
        
        /** @name Adaptation parameters 
         *  constants controlling adaptation 
         */
        //@{
        int adapt_flag; /**< Completely turns adaptation off (can also shut off blocks individiually)  */
        //@}
        
        /** This is a data structure for storing communication info
         *  only used temporarily to sort things out in findmatch
         */
        class commblocks {
            public:
                int nblock;
                struct commblock {
                    int nvcomm;
                    struct vid {
                        int nvbd, idnum;
                    } *vcomm;
                    
                    int nscomm;
                    struct sid {
                        int nsbd, idnum; // , vid0, vid1;
                    } *scomm;
                    
                    int nfcomm;
                    struct fid {
                        int nfbd, idnum; // , nsds;
                        // int *sids;
                    } *fcomm;
                } *blkinfo;
                
                commblocks() : nblock(0) {}
                void output();                 
                int unpack(blitz::Array<int,1> entitylist); 
                ~commblocks() {
                    for(int i=0;i<nblock;++i) {
                        delete []blkinfo[i].vcomm;
                        delete []blkinfo[i].scomm;
                        // for(k=0;k<blkinfo[i].nfcomm;++k)
                        //     delete []blkinfo[i].fcomm[k].sids;
                        delete []blkinfo[i].fcomm;
                    }
                    delete []blkinfo;
                }
        };
        
        /** @name allreduce buffer pointer storage
         *  used in allreduce for block communication
         */
        //@{
        blitz::Array<void *,1> sndbufs, rcvbufs;
        //@}
        
         /** @name variables needed for thread message passing
         */
        //@{
        std::map<int,bool> message_list;  
#if defined(PTH)
        blitz::Array<pth_t,1> threads;
        pth_mutex_t list_mutex;
        pth_cond_t list_change;
        pth_attr_t attr;
#define THREAD_RUN0(f) {for (int b=0;b<nblock;++b) {threads(b) = pth_spawn(attr, thread_proxy<block,&block::f>, blk[b]);} }
#define THREAD_RUN1(f,arg1) {pth_int1 = arg1; for (int b=0;b<nblock;++b) {threads(b) = pth_spawn(attr, thread_proxy<block,&block::f>, blk[b]);}}
#define THREAD_RUN2(f,arg1,arg2) {pth_int1 = arg1; pth_int2 = arg2; for (int b=0;b<nblock;++b) {threads(b) = pth_spawn(attr, thread_proxy<block,&block::f>, blk[b]);}}
#define THREADS_JOIN() {for (int i = 0; i<nblock; ++i) pth_join(threads(i),NULL);}
#elif defined(BOOST)
        blitz::Array<boost::function0<void>,1> thread_func;
        boost::thread_group threads;
        boost::mutex list_mutex;
        boost::condition list_change;
#define THREAD_RUN0(f) {for (int b=0;b<nblock;++b) {thread_func(b) = boost::bind(&block::f,blk[b]); threads.create_thread(thread_func(b));}}
#define THREAD_RUN1(f,arg1) {for (int b=0;b<nblock;++b) {thread_func(b) = boost::bind(&block::f,blk[b],arg1); threads.create_thread(thread_func(b));}}
#define THREAD_RUN2(f,arg1,arg2) {for (int b=0;b<nblock;++b) {thread_func(b) = boost::bind(&block::f,blk[b],arg1,arg2); threads.create_thread(thread_func(b));}}
#define THREADS_JOIN() threads.join_all()
#else
#define THREAD_RUN0(f) {blk[0]->f();}
#define THREAD_RUN1(f,arg1) {blk[0]->f(arg1);}
#define THREAD_RUN2(f,arg1,arg2) {blk[0]->f(arg1,arg2);}
#define THREADS_JOIN() 
#endif
        //@}

    public:
        /** Initialize multiblock/mgrid mesh */
        blocks() : nproc(1), myid(0) {}
        
        /** Initializes blocks using data from map */
        void init(input_map input);
        
        /** Loads map from file then initialize */
        void init(const std::string &infile, const std::string &outfile = std::string());
        
        /** Makes sure vertex positions on boundaries coinside */
        void matchboundaries(int lvl);
        
        /** Outputs solution in various filetypes */
        void output(const std::string &filename, block::output_purpose = block::display, int level = 0);
        
        /** Inputs solution */
        void input(std::string &filename) {}
        
        /** Routine that starts simulation */
        void go();
        
        /** Adapt mesh */
        void restructure(); 
        
        /** Function to allow blocks to perform global sums, averages, etc.. */
        enum operations {sum,max}; /** Only summing and max for now */
        enum msg_type {flt_msg, int_msg};
        void allreduce1(void *sendbuf, void *recvbuf);
        void allreduce2(int count, msg_type datatype, operations op);
        
        /** Functions for thread communication */
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
        
    protected:
        /** Allocates blocks, called by init to generate blocks from initialization file */
        block* getnewblock(int idnum, input_map& blockdata);
                
        /** Sets-up parallel communications, called by init */
        void findmatch();
 
        /** Calculate residuals */
        void rsdl(int);
        
        /** Setup preconditioner */
        void setup_preconditioner(int mglvl);

        /** Iterate on all blocks */  
        void iterate(int mglvl, int niter);

        /** Multigrid cycle */
        void cycle(int vw, int lvl = 0);
        
        /** Shift to next implicit time step */
        void tadvance(); 
                
        /** Print errors */
        FLT maxres(int lvl = 0); 
};


/** \brief Global variables for simulation
 *
 * This namespace contains global variables needed by all blocks of the simulation 
 */
namespace sim {
    extern blocks blks; /**< Contains all blocks for this processor */
    extern FLT dti; /**< Inverse time step */
    extern FLT time; /**< Simulation time */
    extern int tstep; /**< Simulation time step */
    extern int substep; /**< For schemes requiring multiple solves per step */
    extern FLT g;  /**< gravity */
    extern blitz::TinyVector<FLT,2> body; /**< General way for body forces */
    extern std::ostream *log; /**< log file stream */
    
    /** Time stepping for simulation */
#ifdef BACKDIFF
    /** @name backdiff Backwards difference constants
    *  These are constant for backwards difference timestepping
    */
    //@{
    extern FLT bd[BACKDIFF+1];  /**< backwards difference constants */
    const int nhist = BACKDIFF; /**< number of backwards difference steps */
    const int nadapt = BACKDIFF; /**< number of backwards difference steps that require adaptation */
    const int stepsolves = 1;
    //@}
#endif
#ifdef DIRK
    /** @name DIRK variables
     *  These are arrays for diagonally implicit RK Timestepping 
     */
    //@{
    extern FLT bd[1]; /**< Diagonal coefficient */
    extern FLT adirk[DIRK][DIRK]; /**< ``a'' coefficient matrix */
    extern FLT cdirk[DIRK]; /**< ``c'' coefficient matrix */
    //@}
    
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
    extern bool esdirk; /**< Flag to be set when using an explicit 1'st stage */
    //@}
#endif

    /** @name Multistage Explicit Scheme
     *  Constants for multistage iteration scheme 
     */
    //@{
    const int NSTAGE = 5; /**< Number of stages */
    extern blitz::Array<FLT,1> cfl;  /**< Global array of cfl numbers (can use how you see fit) */
    const FLT alpha[NSTAGE+1] = {0.25, 1./6., .375, .5, 1.0, 1.0}; /**< Multistage time step constants (imaginary) */
    const FLT beta[NSTAGE+1] = {1.0, 0.0, 5./9., 0.0, 4./9., 1.0}; /**< Multistage time step constants (real) */
    //@}
    
    extern bool adapt_output; /**< Flag to tell whether to give detailed adaptation data */
}

#endif




