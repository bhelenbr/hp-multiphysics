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
				int myblock; /**< number of blocks on this processor that are in group */
				int nblock; /**< total number of blocks that are in group */
				int buf_cnt;  /**< counter for checking when everyone has called in */
				std::map<int,int> mylocalid;  /**< Map from block number to global sequential list for group */
				blitz::Array<void *,1> sndbufs, rcvbufs;  /**< arrays to keep track of location of send and receive buffers */
#ifdef MPISRC
				MPI_Comm comm;	/**< MPI version of all_reduce_data */
#endif
#if defined(PTH)
				pth_cond_t allreduce_change; /**< conditional await for mpiallreduce */
#elif defined(BOOST)
				boost::condition allreduce_change; /**< conditional await for mpiallreduce */
#endif
				/** allreduce data constructor */
				all_reduce_data() : myblock(1), buf_cnt(0)
#ifdef MPISRC 
					, comm(MPI_COMM_WORLD) 
#endif
					{}
		};
		std::map<int,all_reduce_data *> group_data; /**< Associaiates group numbers to their data */
#if defined(PTH)
		pth_mutex_t allreduce_mutex;
#elif defined(BOOST)
		boost::mutex allreduce_mutex;
#endif
        //@}
        
         /** @name variables needed for threads and thread message passing
         */
        //@{
        std::map<int,bool> message_list; /**< keeps track of whether communication messages have been received */
#if defined(PTH)
        blitz::Array<pth_t,1> threads;
		pth_mutex_t list_mutex;
		pth_cond_t list_change;
        pth_attr_t attr;
#elif defined(BOOST)
        blitz::Array<boost::function0<void>,1> thread_func;
        boost::thread_group threads;
        boost::mutex list_mutex;
        boost::condition list_change;
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
		int allreduce_local_id(int group, int bid) { 
#if defined(PTH)
			pth_mutex_acquire(&allreduce_mutex,false,NULL);
#elif defined(BOOST)
			boost::mutex::scoped_lock lock(allreduce_mutex);
#endif
			int localid = group_data[group]->mylocalid[bid];
#if defined(PTH)
			pth_mutex_release(&allreduce_mutex);
#endif
			return(localid);
		}
	
		int allreduce_nmember(int group) { 
#if defined(PTH)
			pth_mutex_acquire(&allreduce_mutex,false,NULL);
#elif defined(BOOST)
			boost::mutex::scoped_lock lock(allreduce_mutex);
#endif
			int nmember = group_data[group]->nblock;
#if defined(PTH)
			pth_mutex_release(&allreduce_mutex);
#endif
			return(nmember);
		}
	
			        
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
 * This namespace contains global variables and constants needed by all blocks of the simulation 
 */
namespace sim {
    extern blocks blks; /**< Contains all blocks for this processor */
}

#endif




