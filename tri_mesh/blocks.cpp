/*
 *  blocks.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Thu Sep 26 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */


#include "blocks.h"
#include "block.h"
#include <time.h>
#include <input_map.h>
#include <iostream>
#include <sys/param.h>
#include <string>
#include <sstream>
#ifdef MPISRC
#include <mpi.h>
#include <set>
#endif
#include <blitz/array.h>

using namespace std;
using namespace blitz;

#ifdef CAPRI
#include <capri.h>
#endif

#ifdef petsc
#include <petsc.h>
#endif

blocks sim::blks;

#ifdef PTH
struct gostruct {
	block *blk;
	input_map input;
};

static void* thread_go(void* ptr) {
	gostruct *p = static_cast<gostruct *>(ptr);
	p->blk->go(p->input);
	return 0;
}
#endif

static void my_new_handler()
{
	std::cerr << "Out of memory" << endl;
	std::cerr.flush();
	sim::abort(__LINE__,__FILE__,&std::cerr);
}

void blocks::go(const std::string &infile, const std::string &outfile) {
	input_map maptemp;
	std::string name;

	name = infile + ".inpt";
	if (!maptemp.input(name)) {
		if (!maptemp.input(infile)) {
			std::cerr << "Couldn't open input file " << infile << std::endl;
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
	}

	if (!outfile.empty()) {
		maptemp["logfile"] = outfile;
	}

	go(maptemp);

	return;
}

void blocks::go(input_map& input) {
	int i,nb,groupid;
	std::string blkstring,mystring;
	std::ostringstream nstr;
	std::istringstream data;
	std::map<int,all_reduce_data *>::iterator mi;

#ifdef MPISRC
	MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

	input.echo = true;
	input.echoprefix = "#";

	std::set_new_handler(my_new_handler);

	/* EXTRACT NBLOCKS FOR MYID */
	/* SPACE DELIMITED ARRAY OF NBLOCKS FOR EACH PROCESSOR */
	int bstart = 0;
	input.getlinewdefault("nblock",blkstring,std::string("1"));
	data.str(blkstring);
	for (i=0;i<myid;++i) {
		if (!(data >> nb)) {
			std::cerr << "error reading blocks\n";
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
		bstart += nb;
	}
	if (!(data >> myblock)) {
		std::cerr << "error reading blocks\n";
		sim::abort(__LINE__,__FILE__,&std::cerr);
	}
	data.clear();

	nblock = bstart +myblock;
	for (i=myid+1;i<nproc;++i) {
		if (!(data >> nb)) {
			std::cerr << "error reading blocks\n";
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
		nblock += nb;
	}
	blk.resize(myblock);

	/* SET-UP MPI_COMM_WORLD STRUCTURE */
	group_data[0] = new all_reduce_data();
	group_data[0]->myblock = myblock;

	for(i=0;i<myblock;++i) {
		blk(i) = getnewblock(i+bstart,input);
		input.getlinewdefault(blk(i)->idprefix +"_groups",mystring,std::string(""));
		data.str(mystring);
		while (data >> groupid) {
			if (groupid == 0) {
				std::cerr << "Can't use 0 as a group number\n";
				sim::abort(__LINE__,__FILE__,&std::cerr);
			}
			if ((mi = group_data.find(groupid)) != group_data.end()) {
				++mi->second->myblock;
				mi->second->mylocalid[i+bstart] = 1;  // Temporarily use mylocalid's key values to mark blocks in group
			}
			else {
				/* New entry */
				group_data[groupid] = new all_reduce_data();
				group_data[groupid]->mylocalid[i+bstart] = 1; // Temporarily use mylocalid's key values to mark blocks in group
			}
		}
		data.clear();
	}

	/* RESIZE ALLREDUCE BUFFER POINTERS ARRAYS */
	for (mi=group_data.begin();mi != group_data.end(); ++mi) {
		mi->second->sndbufs.resize(mi->second->myblock);
		mi->second->rcvbufs.resize(mi->second->myblock);
	}

#ifdef MPISRC
	/* SET-UP COMMUNICATORS FOR ALL_REDUCE COMMUNICATION */
	/* FOR EACH GROUP NEED LIST OF PROCESSORS */
	std::map<int,std::set<int> > grp_rnk;

	/* Go through list of # of blocks per processor again */
	bstart = 0;
	istringstream nblkstr;
	nblkstr.str(blkstring);
	for (i=0;i<nproc;++i) {
		nblkstr >> nb;
		for (int bn=bstart;bn<bstart+nb;++bn) {
			nstr.str("");
			nstr << "b" << bn << "_groups";
			input.getlinewdefault(nstr.str(),mystring,std::string(""));
			data.str(mystring);
			while (data >> groupid) {
				grp_rnk[groupid].insert(i);
				if ((mi = group_data.find(groupid)) != group_data.end()) {
					/* This group is on this processor */
					mi->second->mylocalid[bn] = 1;  // Temporarily use mylocalid's key values to mark blocks in group
				}
			}
		}
		bstart += nb;
	}

	int ierr,n;
	MPI_Group base_grp, temp_grp;
	MPI_Comm comm_out;
	blitz::Array<int,1> ranks(nproc);

	ierr = MPI_Comm_group (MPI_COMM_WORLD,&base_grp);

	for (std::map<int,std::set<int> >::iterator gri = grp_rnk.begin(); gri != grp_rnk.end(); ++gri) {
		/* Create array of ranks */
		n = gri->second.size();
		for (std::set<int>::iterator si = gri->second.begin(); si != gri->second.end();++si) {
			ranks(i) = *si;
		}

		ierr += MPI_Group_incl(base_grp, n, ranks.data(), &temp_grp);
		ierr += MPI_Comm_create (MPI_COMM_WORLD, temp_grp, &comm_out);
		if ((mi = group_data.find(gri->first)) != group_data.end()) {
			mi->second->comm = comm_out;
		}
	}
#endif
	
	int partition;
	input.getwdefault("partition",partition,0);
	/* Allocate memory for partitioning */
	if (partition) allocate_shared_memory(partition*nblock*100, sizeof(int));

	/* Give blocks in localid for each group sequential numbering */
	for (mi=group_data.begin();mi != group_data.end(); ++mi) {
		int count = 0;
		for (std::map<int,int>::iterator lid = mi->second->mylocalid.begin(); lid != mi->second->mylocalid.end(); ++lid) {
			lid->second = count++;
		}
		mi->second->nblock = count;
	}


#ifdef PTH
	threads.resize(myblock);
	pth_mutex_init(&list_mutex);
	pth_cond_init(&list_change);
	pth_mutex_init(&allreduce_mutex);
	for (mi=group_data.begin();mi != group_data.end(); ++mi)
		pth_cond_init(&(mi->second->allreduce_change));
	pth_mutex_init(&data_mutex);
	attr = PTH_ATTR_DEFAULT;
	pth_attr_set(attr, PTH_ATTR_JOINABLE, true);

	Array<gostruct,1> myGo(myblock);
	for (int b=0;b<myblock;++b) {
		myGo(b).input = input;
		myGo(b).blk = blk(b);
		threads(b) = pth_spawn(attr, thread_go,static_cast<void *>(&myGo(b)));
	}
	
	for (int b=0;b<myblock;++b)
		pth_join(threads(b),NULL);

#elif defined(BOOST)
	thread_func.resize(myblock);
	for (int b=0;b<myblock;++b) {
		thread_func(b) = boost::bind(&block::go,blk(b),input);
		threads.create_thread(thread_func(b));
	}

	threads.join_all();
#else
	if (myblock != 1) {
		std::cerr << "Need pth or boost::threads to run mulitple blocks on single processor\n";
		sim::abort(__LINE__,__FILE__,&std::cerr);
	}
	blk(0)->go(input);
#endif
	
	if (partition) free_shared_memory();


}


/* This routine waits for everyone to exit nicely */
void sim::finalize(int line,const char *file, std::ostream *log) {
	*log << "Exiting at line " << line << " of file " << file << std::endl;
#ifdef PTH
	pth_exit(NULL);
#endif
#ifdef BOOST
	throw boost::thread_interrupted();
#endif
#ifdef PTH
	pth_kill();
#endif
#ifdef petsc
	PetscFinalize();
#endif
#ifdef MPISRC
	MPI_Finalize();
#endif
    
    std::exit(0);
}

/* This routine forces everyone to die */
void sim::abort(int line,const char *file, std::ostream *log) {
	*log << "Exiting at line " << line << " of file " << file << std::endl;
	for (int b=0;b<blks.myblock;++b) {
		sim::blks.blk(b)->output("aborted_solution", block::display);
		sim::blks.blk(b)->output("aborted_solution", block::restart);
	}
#ifdef petsc
	PetscFinalize();
#endif
#ifdef MPI
	MPI_Abort(MPI_COMM_WORLD,1);
#endif
    
	/* Terminates all threads */
	std::exit(1);
}

/* each block has a list of group #'s that it belongs to: integer array of size "n_comm_purposes" */
/* typicaally there will be 2 comm_purposes: everyone to everyone, and only 1 with user defined groups in it */
/* Each entry in array corresponds to different communication purpose entry 1 always has value 1 and that */
/* corresponds to MPI_COMM_WORLD, -1 means block is not active for that communication purpose */
/* blocks has to know how many different groups there are for each purpose: array<int, 1> n_comm_entries_for_purpose(n_comm_purposes) */
/* and must be able to associaciate those group numbers to a sequential integer count: array<map<int,int>,1> buffer_number(n_comm_purposes) */
/* Then blocks has to have buffer pointer arrays for the total # of groups in each purpose: array<1, array<1, void *> > sndbuf(n_comm_purposes), rcvbuf(n_comm_purposes)  */

void blocks::allreduce(void *sendbuf, void *recvbuf, int count, msg_type datatype, operations op, int group) {
	int i,j;
	int *ircvbuf;
	FLT *frcvbuf;

#if defined(PTH)
	pth_mutex_acquire(&allreduce_mutex,false,NULL);
#elif defined(BOOST)
	boost::mutex::scoped_lock lock(allreduce_mutex);
#endif

	std::map<int,all_reduce_data *>::iterator gi = group_data.find(group);
	if (gi == group_data.end()) {
		std::cerr << "couldn't find group " << group << " in all_reduce\n";
		sim::abort(__LINE__,__FILE__,&std::cerr);
	}
	all_reduce_data& gd(*gi->second);
	gd.sndbufs(gd.buf_cnt) = sendbuf;
	gd.rcvbufs(gd.buf_cnt) = recvbuf;
	++gd.buf_cnt;

	if (gd.buf_cnt == gd.myblock) {
		switch(datatype) {
			case(int_msg): {

				Array<int,1> isendbuf(count);

				switch(op) {
					case(sum): {
						isendbuf = 0;
						for(j=0;j<gd.buf_cnt;++j) {
							for(i=0;i<count;++i) {
								isendbuf(i) += static_cast<int *>(gd.sndbufs(j))[i];
							}
						}
#ifdef MPISRC
						MPI_Allreduce(isendbuf.data(),gd.rcvbufs(0),count,MPI_INT,MPI_SUM,gd.comm);

						for(i=1;i<gd.buf_cnt;++i) {
							ircvbuf = static_cast<int *>(gd.rcvbufs(i));
							for(j=0;j<count;++j)
								ircvbuf[j] =  static_cast<int *>(gd.rcvbufs(0))[j];
						}
#else
						for(i=0;i<gd.buf_cnt;++i) {
							ircvbuf = static_cast<int *>(gd.rcvbufs(i));
							for(j=0;j<count;++j)
								ircvbuf[j] = isendbuf(j);
						}
#endif
						break;
					}

					case(max): {
						for(j=0;j<count;++j) {
							isendbuf(j) = static_cast<int *>(gd.sndbufs(0))[j];
						}

						for(i=1;i<gd.buf_cnt;++i) {
							for(j=0;j<count;++j) {
								isendbuf(j) = MAX(static_cast<int *>(gd.sndbufs(i))[j],isendbuf(j));
							}
						}
#ifdef MPISRC
						MPI_Allreduce(isendbuf.data(),gd.rcvbufs(0),count,MPI_INT,MPI_MAX,gd.comm);

						for(i=1;i<gd.buf_cnt;++i) {
							ircvbuf = static_cast<int *>(gd.rcvbufs(i));
							for(j=0;j<count;++j)
								ircvbuf[j] =  static_cast<int *>(gd.rcvbufs(0))[j];
						}
#else
						for(i=0;i<gd.buf_cnt;++i) {
							ircvbuf = static_cast<int *>(gd.rcvbufs(i));
							for(j=0;j<count;++j)
								ircvbuf[j] = isendbuf(j);
						}
#endif
						break;
					}
				}
				break;
			}


			case(flt_msg): {

				Array<FLT,1> fsendbuf(count);

				switch(op) {
					case(sum): {
						fsendbuf = 0;
						for(i=0;i<gd.buf_cnt;++i) {
							for(j=0;j<count;++j) {
								fsendbuf(j) += static_cast<FLT *>(gd.sndbufs(i))[j];
							}
						}
#ifdef MPISRC
						MPI_Allreduce(fsendbuf.data(),gd.rcvbufs(0),count,MPI_DOUBLE,MPI_SUM,gd.comm);

						for(i=1;i<gd.buf_cnt;++i) {
							frcvbuf = static_cast<FLT *>(gd.rcvbufs(i));
							for(j=0;j<count;++j)
								frcvbuf[j] =  static_cast<FLT *>(gd.rcvbufs(0))[j];
						}
#else
						for(i=0;i<gd.buf_cnt;++i) {
							frcvbuf = static_cast<FLT *>(gd.rcvbufs(i));
							for(j=0;j<count;++j)
								frcvbuf[j] = fsendbuf(j);
						}
#endif
						break;
					}

					case(max): {
						for(j=0;j<count;++j) {
							fsendbuf(j) = static_cast<FLT *>(gd.sndbufs(0))[j];
						}

						for(i=1;i<gd.buf_cnt;++i) {
							for(j=0;j<count;++j) {
								fsendbuf(j) = MAX(static_cast<FLT *>(gd.sndbufs(i))[j],fsendbuf(j));
							}
						}
#ifdef MPISRC
						MPI_Allreduce(fsendbuf.data(),gd.rcvbufs(0),count,MPI_DOUBLE,MPI_MAX,gd.comm);

						for(i=1;i<gd.buf_cnt;++i) {
							frcvbuf = static_cast<FLT *>(gd.rcvbufs(i));
							for(j=0;j<count;++j)
								frcvbuf[j] =  static_cast<FLT *>(gd.rcvbufs(0))[j];
						}
#else
						for(i=0;i<gd.buf_cnt;++i) {
							frcvbuf = static_cast<FLT *>(gd.rcvbufs(i));
							for(j=0;j<count;++j)
								frcvbuf[j] = fsendbuf(j);
						}
#endif
						break;
					}
				}
				break;
			}
		}
		gd.buf_cnt = 0;

#if defined(PTH)
		pth_cond_notify(&gd.allreduce_change, true);
#elif defined(BOOST)
		gd.allreduce_change.notify_all();
#endif
	}
	else {
#if defined(PTH)
		pth_cond_await(&gd.allreduce_change, &allreduce_mutex, NULL);
#elif defined(BOOST)
		gd.allreduce_change.wait(lock);
#endif
	}
#if defined(PTH)
	pth_mutex_release(&allreduce_mutex);
#endif
	return;
}

void blocks::allocate_shared_memory(int nentry, size_t size) {
#if defined(PTH)
	pth_mutex_acquire(&shared_mem_mutex,false,NULL);
#elif defined(BOOST)
	boost::mutex::scoped_lock lock(shared_mem_mutex);
#endif
//	if (shared_mem_call_count == 0) {
		shared_mem_size = size*nentry;
#ifndef MPISRC
		shared_mem = malloc(shared_mem_size);
        if (!shared_mem) {
            std::cerr << "Could not allocate memory!" << std::endl;
            sim::abort(__LINE__,__FILE__,&std::cerr);
        }
#else
		int ierr = MPI_Alloc_mem(shared_mem_size, MPI_INFO_NULL, &shared_mem);
		if (ierr != MPI_SUCCESS) {
			std::cerr << "couldn't allocate shared memory for partitioning " << ierr << std::endl;
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
#endif
		Array<FLT,1> temp(static_cast<FLT *>(shared_mem), (shared_mem_size)/sizeof(FLT), neverDeleteData);
		fshared_mem.reference(temp);
		Array<int,1> temp1(static_cast<int *>(shared_mem), (shared_mem_size)/sizeof(int), neverDeleteData);
		ishared_mem.reference(temp1);
		ishared_mem = 0;
		
#ifdef MPISRC
		if (myid == 0) {
			/* Rank 0 will hold the master of the data */
			int ierr = MPI_Win_create(shared_mem,shared_mem_size,1,MPI_INFO_NULL,MPI_COMM_WORLD,&shared_mem_win);
			if (ierr != MPI_SUCCESS) {
				std::cerr << "couldn't create rank 0 window " << ierr << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
			}
		}
		else {
			/* Others do not receive asynchronous messages */
			int ierr = MPI_Win_create(NULL,0,1,MPI_INFO_NULL, MPI_COMM_WORLD, &shared_mem_win);
			if (ierr != MPI_SUCCESS) {
				std::cerr << "couldn't create non rank 0 window " << ierr << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
			}
		}
#endif
//	}
//	++shared_mem_call_count;
//	
//	if(shared_mem_call_count == myblock)
//		shared_mem_call_count = 0;

#if defined(PTH)
	pth_mutex_release(&shared_mem_mutex);
#endif
	return;
}

void blocks::begin_use_shared_memory() {
#if defined(PTH)
	pth_mutex_acquire(&shared_mem_mutex,false,NULL);
#elif defined(BOOST)
	shared_mem_mutex.lock();
#endif
#ifdef MPISRC
	if (myid == 0) {
		/* Request lock of process 1 */
		int ierr = MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,shared_mem_win);
		if (ierr != MPI_SUCCESS) {
			std::cerr << "couldn't get lock rank 0" << ierr << std::endl;
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
	}
	else {
		int ierr = MPI_Win_lock(MPI_LOCK_EXCLUSIVE,0,0,shared_mem_win);
		if (ierr != MPI_SUCCESS) {
			std::cerr << "couldn't get lock rank != 0" << ierr << std::endl;
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
		ierr = MPI_Get(shared_mem,shared_mem_size,MPI_CHAR,0,0,shared_mem_size,MPI_CHAR,shared_mem_win);
		if (ierr != MPI_SUCCESS) {
			std::cerr << "couldn't get shared_mem" << ierr << std::endl;
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
		
		ierr =  MPI_Win_flush(0,shared_mem_win);
		if (ierr != MPI_SUCCESS) {
			std::cerr << "couldn't get shared_mem" << ierr << std::endl;
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
	}
	std::cout << "myid " << myid << "entering shared memory" << std::endl;
	std::cout << "nentry " << ishared_mem(0) << std::endl;
#endif
}

void blocks::end_use_shared_memory() {
	
	std::cout << "myid " << myid << "leaving shared memory" << std::endl;
	std::cout << "nentry " << ishared_mem(0) << std::endl;

#ifdef MPISRC
	if (myid != 0) {
		/* Request lock of process 1 */
		int ierr = MPI_Put(shared_mem,shared_mem_size,MPI_CHAR,0,0,shared_mem_size,MPI_CHAR,shared_mem_win);
		if (ierr != MPI_SUCCESS) {
			std::cerr << "couldn't put shared_mem" << ierr << std::endl;
			sim::abort(__LINE__,__FILE__,&std::cerr);
		}
	}
	/* Block until put succeeds */
	int ierr = MPI_Win_unlock(0,shared_mem_win);
	if (ierr != MPI_SUCCESS) {
		std::cerr << "couldn't unlock shared_mem" << ierr << std::endl;
		sim::abort(__LINE__,__FILE__,&std::cerr);
	}
#endif
	
#if defined(PTH)
	pth_mutex_release(&shared_mem_mutex);
#elif defined(BOOST)
	shared_mem_mutex.unlock();
#endif

}


void blocks::free_shared_memory() {
#if defined(PTH)
	pth_mutex_acquire(&shared_mem_mutex,false,NULL);
#elif defined(BOOST)
	boost::mutex::scoped_lock lock(shared_mem_mutex);
#endif
	
//	if (shared_mem_call_count == 0) {
		
#ifndef MPISRC
		free(shared_mem);
#else
		MPI_Free_mem(shared_mem);
		MPI_Win_free(&shared_mem_win);
#endif
//	}
//	++shared_mem_call_count;
//	
//	if(shared_mem_call_count == myblock)
//		shared_mem_call_count = 0;

#if defined(PTH)
	pth_mutex_release(&shared_mem_mutex);
#endif
	return;
}





/** This is a data structure for storing communication info
 *  only used temporarily to sort things out in findmatch
 */
class multigrid_interface::comm_info {
	public:
		int proc;
		int blk;

		int nvcomm;
		struct vid {
			int nvbd, idnum;
		} *vcomm;

		int nscomm;
		struct sid {
			int nebd, idnum;
		} *ecomm;

		int nfcomm;
		struct fid {
			int nfbd, idnum;
		} *fcomm;

		comm_info() {}
		int unpack(blitz::Array<int,1> entitylist);
		~comm_info() {
			delete []vcomm;
			delete []ecomm;
			delete []fcomm;
		}
};

int multigrid_interface::comm_info::unpack(blitz::Array<int,1> entitylist) {
	int count = 0;
	proc = entitylist(count++);
	blk = entitylist(count++);
	nvcomm = entitylist(count++);
	vcomm = new comm_info::vid[nvcomm];
	for(int i=0;i<nvcomm;++i) {
		vcomm[i].nvbd = entitylist(count++);
		vcomm[i].idnum = entitylist(count++);
	}

	nscomm = entitylist(count++);
	ecomm = new comm_info::sid[nscomm];
	for(int i=0;i<nscomm;++i) {
		ecomm[i].nebd = entitylist(count++);
		ecomm[i].idnum = entitylist(count++);
	}

	nfcomm = entitylist(count++);
	fcomm = new comm_info::fid[nfcomm];
	for(int i=0;i<nfcomm;++i) {
		fcomm[i].nfbd = entitylist(count++);
		fcomm[i].idnum = entitylist(count++);
	}
	return(count);
}

void multigrid_interface::findmatch(shared_ptr<block_global> gbl, int grdlvl) {

	/* GET DATA FROM BLOCKS IN A THREAD SAFE WAY */
#if defined(PTH)
	pth_mutex_acquire(&sim::blks.data_mutex,false,NULL);
#elif defined(BOOST)
	boost::mutex::scoped_lock datalock(sim::blks.data_mutex);
#endif
	const int nblock = sim::blks.nblock;
	const int myid = sim::blks.myid;
	const int myblock = sim::blks.myblock;
	const int idnum = gbl->idnum;

	/* FIGURE OUT MY LOCAL BLOCK NUMBER & GRID LEVEL? */
	int b1;
	for (b1=0;b1<myblock;++b1)
		if (sim::blks.blk(b1)->idprefix == gbl->idprefix) break;
	if (b1 >= myblock) {
		*gbl->log << "Didn't find myself in block list?\n";
		sim::abort(__LINE__,__FILE__,gbl->log);
	}




#if defined(PTH)
	pth_mutex_release(&sim::blks.data_mutex);
#elif defined(BOOST)
	datalock.unlock();
#endif

	/* FIRST DETERMINE TOTAL SIZE OF LIST */
	Array<int,1> sndsize(nblock), size(nblock);
	sndsize = 0;
	sndsize(idnum) = 2 +comm_entity_size(); // 2 is for myid & local block number
	sim::blks.allreduce(sndsize.data(),size.data(),nblock,blocks::int_msg,blocks::sum);
	~sndsize;

	int tsize = 0;
	for(int i=0;i<nblock;++i)
		tsize += size(i);

	Array<int,1> sndentitylist(tsize);
	sndentitylist = 0;

	/* CALCULATE BEGINNING LOCATION FOR THIS BLOCK*/
	int count = 0;
	for(int i=0;i<idnum;++i)
		count += size(i);

	/* START ASSEMBLING LIST */
	sndentitylist(count++) = myid;
	sndentitylist(count++) = b1;

	Array<int,1> sublist;
	sublist.reference(sndentitylist(Range(count,toEnd)));
	comm_entity_list(sublist);
	~sublist;

	Array<int,1> entitylist(tsize);
	entitylist = 0;
	sim::blks.allreduce(sndentitylist.data(),entitylist.data(),tsize,blocks::int_msg,blocks::sum);
	~sndentitylist;

	/* UNPACK ENTITY LIST INTO MORE USABLE FORM */
	int max_comm_bdry_num = 0; // SO I CAN GENERATE UNIQUE TAGS
	count = 0;
	multigrid_interface::comm_info *binfo = new multigrid_interface::comm_info[nblock];
	for(int i=0;i<nblock;++i) {
		count += binfo[i].unpack(entitylist(Range(count,toEnd)));
		max_comm_bdry_num = MAX(max_comm_bdry_num,binfo[i].nvcomm);
		max_comm_bdry_num = MAX(max_comm_bdry_num,binfo[i].nscomm);
		max_comm_bdry_num = MAX(max_comm_bdry_num,binfo[i].nfcomm);
	}
	/* CALCULATE NUMBER OF BINARY DIGITS NEEDED TO MAKE TAGS */
	sim::blks.setdigits(max_comm_bdry_num);

	/* OUTPUT LIST FOR DEBUGGING */
	*gbl->log << "# block " << gbl->idprefix << " block data" << std::endl;
	*gbl->log << "# myid: " << binfo[idnum].proc << "\t local block number: " << binfo[idnum].blk << std::endl;
	*gbl->log << "#\t\tnvcomm: " << binfo[idnum].nvcomm << std::endl;
	for (int i=0;i<binfo[idnum].nvcomm;++i)
		*gbl->log << "#\t\t\tnvbd: " << binfo[idnum].vcomm[i].nvbd << " idnum: " << binfo[idnum].vcomm[i].idnum << std::endl;

	*gbl->log << "#\t\tnscomm: " << binfo[idnum].nscomm << std::endl;
	for (int i=0;i<binfo[idnum].nscomm;++i)
		*gbl->log << "#\t\t\tnsbd: " << binfo[idnum].ecomm[i].nebd << " idnum: " << binfo[idnum].ecomm[i].idnum << std::endl;

	*gbl->log << "#\t\tnfcomm: " << binfo[idnum].nfcomm << std::endl;
	for (int i=0;i<binfo[idnum].nfcomm;++i)
		*gbl->log << "#\t\t\tnfbd: " << binfo[idnum].fcomm[i].nfbd << " idnum: " << binfo[idnum].fcomm[i].idnum << std::endl;

	~entitylist;


	/* CAN NOW START TO LOOK FOR MATCHES */
	/* GOING TO DO THIS THREAD SAFE I HOPE BECAUSE OF ACCESS TO SIM */
#if defined(PTH)
	pth_mutex_acquire(&sim::blks.data_mutex,false,NULL);
#elif defined(BOOST)
	datalock.lock();
#endif
	b1 = idnum;

	*gbl->log << "# finding matches for block " << b1 << std::endl;
	/* LOOK FOR VERTEX MATCHES */
	for(int i=0;i<binfo[b1].nvcomm;++i) {
		bool first_found = false;
		*gbl->log << "#\tvertex " << binfo[b1].vcomm[i].nvbd << " idnum: " << binfo[b1].vcomm[i].idnum << std::endl;
		for(int b2=0;b2<nblock;++b2) {
			for(int j=0;j<binfo[b2].nvcomm;++j) {
				if (binfo[b1].vcomm[i].idnum == binfo[b2].vcomm[j].idnum) {
					if (b1 == b2 && i == j) {
						if (!first_found) first_found = true;  // Leave first flag alone
						continue;  // CAN"T MATCH TO MYSELF
					}

					boundary *bp1 = getvbdry(binfo[b1].vcomm[i].nvbd);
					if (binfo[b1].proc == binfo[b2].proc) {
						/* local match */
						boundary *bp2 = sim::blks.blk(binfo[b2].blk)->grd(grdlvl)->getvbdry(binfo[b2].vcomm[j].nvbd);
						bp1->local_cnnct(bp2,sim::blks.tagid(1,b1,b2,i,j),sim::blks.tagid(1,b2,b1,j,i));
						if (!first_found) {
							bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
							first_found = true;
						}
						*gbl->log <<  "#\t\tlocal match to block: " << b2 << " pnt: " << binfo[b2].vcomm[j].nvbd << " tag: " << sim::blks.tagid(1,b1,b2,i,j) << ' ' << sim::blks.tagid(1,b2,b1,j,i) << " idnum: " << binfo[b2].vcomm[j].idnum << " first: " << bp1->is_frst() << std::endl;
					}
#ifdef MPISRC
					else {
						ostringstream nstr;
						nstr << 'b' << b2 << "_v" << binfo[b2].vcomm[j].idnum;
						bp1->mpi_cnnct(binfo[b2].proc,sim::blks.tagid(1,b1,b2,i,j),sim::blks.tagid(1,b2,b1,j,i),nstr.str());
						if (!first_found) {
							bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
							first_found = true;
						}
						*gbl->log <<  "#\t\t  mpi match to  block: " << b2 << " pnt: " << binfo[b2].vcomm[j].nvbd << " tag: " << sim::blks.tagid(1,b1,b2,i,j) << ' ' << sim::blks.tagid(1,b2,b1,j,i) << " idnum: " << binfo[b2].vcomm[j].idnum << " first: " << bp1->is_frst() << std::endl;
					}
#endif

				}
			}
		}
	}

	/* LOOK FOR SIDE MATCHES */
	for(int i=0;i<binfo[b1].nscomm;++i) {
		bool first_found = false;
		*gbl->log << "#\tside " << binfo[b1].ecomm[i].nebd << " idnum: " << binfo[b1].ecomm[i].idnum << std::endl;
		for(int b2=0;b2<nblock;++b2) {
			for(int j=0;j<binfo[b2].nscomm;++j) {
				if (binfo[b1].ecomm[i].idnum == binfo[b2].ecomm[j].idnum) {
					if (b1 == b2 && i == j) {
						if (!first_found) first_found = true;  // Leave first flag alone
						continue;  // CAN"T MATCH TO MYSELF
					}

					boundary *bp1 = getebdry(binfo[b1].ecomm[i].nebd);
					if (binfo[b1].proc == binfo[b2].proc) {
						/* local match */
						boundary *bp2 = sim::blks.blk(binfo[b2].blk)->grd(grdlvl)->getebdry(binfo[b2].ecomm[j].nebd);
						bp1->local_cnnct(bp2,sim::blks.tagid(2,b1,b2,i,j),sim::blks.tagid(2,b2,b1,j,i));
						if (!first_found) {
							bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
							first_found = true;
						}
						*gbl->log <<  "#\t\tlocal match to block: " << b2 << " edge: " << binfo[b2].ecomm[j].nebd << " tag: " << sim::blks.tagid(2,b1,b2,i,j) << ' ' << sim::blks.tagid(2,b2,b1,j,i) << " idnum: " << binfo[b2].ecomm[j].idnum << " first: " << bp1->is_frst() << std::endl;
					}
#ifdef MPISRC
					else {
						ostringstream nstr;
						nstr << 'b' << b2 << "_s" << binfo[b2].ecomm[j].idnum;
						bp1->mpi_cnnct(binfo[b2].proc,sim::blks.tagid(2,b1,b2,i,j),sim::blks.tagid(2,b2,b1,j,i),nstr.str());
						if (!first_found) {
							bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
							first_found = true;
						}
						*gbl->log <<  "#\t\t  mpi match to  block: " << b2 << " edge: " << binfo[b2].ecomm[j].nebd << " tag: " << sim::blks.tagid(2,b1,b2,i,j) << ' ' << sim::blks.tagid(2,b2,b1,j,i) << " idnum: " << binfo[b2].ecomm[j].idnum << " first: " << bp1->is_frst() << std::endl;
					}
#endif

				}
			}
		}
	}



	/* LOOK FOR FACE MATCHES */
	for(int i=0;i<binfo[b1].nfcomm;++i) {
		bool first_found = false;
		*gbl->log << "#\tface " << binfo[b1].fcomm[i].nfbd << " idnum: " << binfo[b1].fcomm[i].idnum << std::endl;
		for(int b2=0;b2<nblock;++b2) {
			for(int j=0;j<binfo[b2].nfcomm;++j) {
				if (binfo[b1].fcomm[i].idnum == binfo[b2].fcomm[j].idnum) {
					if (b1 == b2 && i == j) {
						if (!first_found) first_found = true;  // Leave first flag alone
						continue;  // CAN"T MATCH TO MYSELF
					}

					boundary *bp1 = getfbdry(binfo[b1].fcomm[i].nfbd);
					if (binfo[b1].proc == binfo[b2].proc) {
						/* local match */
						boundary *bp2 = sim::blks.blk(binfo[b2].blk)->grd(grdlvl)->getfbdry(binfo[b2].fcomm[j].nfbd);
						bp1->local_cnnct(bp2,sim::blks.tagid(3,b1,b2,i,j),sim::blks.tagid(3,b2,b1,j,i));
						if (!first_found) {
							bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
							first_found = true;
						}
						*gbl->log <<  "#\t\tlocal match to block: " << b2 << " face: " << binfo[b2].fcomm[j].nfbd << " tag: " << sim::blks.tagid(3,b1,b2,i,j) << ' ' << sim::blks.tagid(3,b2,b1,j,i) << " idnum: " << binfo[b2].fcomm[j].idnum << " first: " << bp1->is_frst() << std::endl;
					}
#ifdef MPISRC
					else {
						ostringstream nstr;
						nstr << 'b' << b2 << "_f" << binfo[b2].fcomm[j].idnum;
						bp1->mpi_cnnct(binfo[b2].proc,sim::blks.tagid(3,b1,b2,i,j),sim::blks.tagid(3,b2,b1,j,i),nstr.str());
						if (!first_found) {
							bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
							first_found = true;
						}
						*gbl->log <<  "#\t\t  mpi match to  block: " << b2 << " face: " << binfo[b2].fcomm[j].nfbd << " tag: " << sim::blks.tagid(3,b1,b2,i,j) << ' ' << sim::blks.tagid(3,b2,b1,j,i) << " idnum: " << binfo[b2].fcomm[j].idnum << " first: " << bp1->is_frst() << std::endl;
					}
#endif

				}
			}
		}
	}

#if defined(PTH)
	pth_mutex_release(&sim::blks.data_mutex);
#elif defined(BOOST)
	datalock.unlock();
#endif

	/* DELETE DATA STRUCTURE */
	delete []binfo;


	/* Now match internal numbering system of boundaries (if necessary) */
	match_bdry_numbering();
	
	/* Update halo information */
	calculate_halo();

	return;
}

void block::init(input_map &input) {
	std::string mystring,stem;

	/* NEED TO BOOTSTRAP UNTIL I CAN GET LOGFILE OPENED */
	input.getwdefault("ngrid",ngrid,1);
	ngrid = MAX(ngrid,1);
	
	grd.resize(ngrid);
	for (int i=0;i<ngrid;++i)
		grd(i) = getnewlevel(input);
	gbl = make_shared<block_global>();
	gbl->idnum = idnum;
	gbl->idprefix = idprefix;

	/* OPEN LOGFILES FOR EACH BLOCK */
	if (input.get("logfile",stem)) {
		mystring = stem +"_" +idprefix +".log";
		std::ofstream *filelog = new std::ofstream;
		filelog->setf(std::ios::scientific, std::ios::floatfield);
		filelog->precision(3);	
		
		/* Check for pre-existing log file and rotate old log out of the way */
		
		FILE *f = fopen(mystring.c_str(),"r");
		if (f) {
			/* File exists, lets move it to mystring_idprefix.#.log */
			ostringstream nstr;
			int num = 0;
			do {
				++num;
				fclose(f);
				nstr.str("");
				nstr << stem  << "_"  << idprefix << "." << num << ".log";
				f = fopen(nstr.str().c_str(),"r");
			} while (f);
			rename(mystring.c_str(),nstr.str().c_str());
		}
		
		filelog->open(mystring.c_str());
		gbl->log = filelog;
		
		/* GOING TO MAKE STD::COUT POINT TO FILE AS WELL */
//		streambuf *psbuf, *backup;		
//		backup = cout.rdbuf();     // back up cout's streambuf
//		psbuf = filelog->rdbuf();   // get file's streambuf
//		cout.rdbuf(psbuf);         // assign streambuf to cout
//		cerr.rdbuf(psbuf);				// assing streambuf to cerr as well
	}
	else {
		std::cout.setf(std::ios::scientific, std::ios::floatfield);
		std::cout.precision(3);
		gbl->log = &std::cout;
	}
	input.echo = true;
	input.log = gbl->log;
	input.echoprefix = "#";
    
#ifdef BZ_DEBUG
    *gbl->log << "#BZ_DEBUG is set\n";
#endif
#ifdef DEBUG
    *gbl->log << "#Running in Xcode's DEBUG Mode\n";
#endif

#ifdef CAPRI
	int status = gi_uStart();
	*gbl->log << "gi_uStart status = ", status, "\n";
	if (status != CAPRI_SUCCESS) sim::abort(__LINE__,__FILE__,gbl->log);

	if (input.get("BRep",mystring)) {
		status = gi_uLoadPart(mystring.c_str());
		*gbl->log << mystring << ": gi_uLoadPart status = " << status << std::endl;
		if (status != CAPRI_SUCCESS) sim::abort(__LINE__,__FILE__,gbl->log);
	}
#endif
    
    int time_scheme;
    input.getwdefault("time_scheme",time_scheme,static_cast<int>(block_global::time_schemes::DIRK4));
    gbl->time_scheme = static_cast<block_global::time_schemes>(time_scheme);
	
    /* SET TIME STEPPING INFO */
    switch (gbl->time_scheme) {
        case(block_global::time_schemes::DIRK4):
            /* ESDIRK SCHEME */
            gbl->nadapt = 2;
            gbl->nhist = 4;
            gbl->stepsolves = 3;
            gbl->bd.resize(1);
            gbl->bd = 0.0;
            gbl->adirk.resize(4,4);
            gbl->cdirk.resize(4);
            break;
        case(block_global::time_schemes::DIRK3):
            /* DIRK SCHEME */
            gbl->nadapt = 2;
            gbl->nhist = 4;
            gbl->stepsolves = 3;
            gbl->bd.resize(1);
            gbl->bd = 0.0;
            gbl->adirk.resize(3,3);
            gbl->cdirk.resize(3);
            break;
        case(block_global::time_schemes::DIRK2):
            /* BD1 WITH EXTRAPOLATION */
            gbl->nadapt = 2;
            gbl->nhist = 2;
            gbl->stepsolves = 1;
            gbl->bd.resize(1);
            gbl->bd = 0.0;
            gbl->adirk.resize(2,2);
            gbl->cdirk.resize(2);
            break;
        case(block_global::time_schemes::DIRK1):
            /* BD1 WITH NO EXTRAPOLATION */
            gbl->nadapt = 1;
            gbl->nhist = 1;
            gbl->stepsolves = 1;
            gbl->bd.resize(1);
            gbl->bd = 0.0;
            gbl->adirk.resize(1,1);
            gbl->cdirk.resize(1);
            break;
        case(block_global::time_schemes::AM1):
            /* This is the AM1 with extrapolation scheme */
            gbl->nadapt = 2;
            gbl->nhist = 2;
            gbl->stepsolves = 1;
            gbl->bd.resize(1);
            gbl->bd = 0.0;
            gbl->adirk.resize(2,2);
            gbl->cdirk.resize(2);
            break;
        case(block_global::time_schemes::BD1):
            /* This is the BD1 scheme */
            /* This is the 2-step backwards difference scheme */
            gbl->nhist = 1;
            gbl->nadapt = 1;
            gbl->bd.resize(1);
            gbl->stepsolves = 1;
            break;
        case(block_global::time_schemes::BD2):
            /* This is the 2-step backwards difference scheme */
            gbl->nhist = 2;
            gbl->nadapt = 2;
            gbl->bd.resize(2);
            gbl->stepsolves = 1;
            break;
        case(block_global::time_schemes::BD3):
            /* This is the 2-step backwards difference scheme */
            gbl->nhist = 3;
            gbl->nadapt = 3;
            gbl->bd.resize(3);
            gbl->stepsolves = 1;
            break;
    }
    
	gbl->tstep = -1; // Simulation starts at t = 0, This is set negative until first tadvance to alllow change between initialization & B.C.'s
	gbl->substep = -1;
	input.getwdefault("dtinv",gbl->dti,0.0);
	if (gbl->dti < 0.0) {
		*gbl->log << "Inverse time step must be positive number" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	input.getwdefault("dtinv_prev",gbl->dti_prev,gbl->dti);
	if (gbl->dti_prev < 0.0) {
		*gbl->log << "Inverse time step must be positive number" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	input.getwdefault("extrapolate",gbl->extrapolate,0.0);
	if (gbl->extrapolate < 0.0 || gbl->extrapolate > 1.0) {
		*gbl->log << "Guess extrapolation constant should between 0 and 1" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
    
    input.getwdefault("auto_timestep_tries",gbl->auto_timestep_tries,0);
    gbl->recursive_timestep_levels = 0;
    if (gbl->auto_timestep_tries) {
        if (gbl->dti <= 0.0) {
            *gbl->log << "Auto timestepping requires an initial time step"  << std::endl;
            sim::abort(__LINE__,__FILE__,gbl->log);
        }
        if (gbl->time_scheme > 4) {
            *gbl->log << "Auto timestepping can only be used with DIRK schemes"  << std::endl;
            sim::abort(__LINE__,__FILE__,gbl->log);
        }
        input.getwdefault("auto_timestep_ratio",gbl->auto_timestep_ratio,2.0);
        input.getwdefault("auto_dti_min",gbl->auto_dti_min,0.0);
        input.getwdefault("auto_dti_max",gbl->auto_dti_max,gbl->dti*128.0);
        input.getwdefault("auto_timestep_maxtime",gbl->auto_timestep_maxtime,-1.0);
    }
    else {
        input.getwdefault("recursive_timestep_levels",gbl->recursive_timestep_levels,0);
        input.getwdefault("recursive_start_level",gbl->recursive_level,0);
        input.getwdefault("recursive_nsuccesses",gbl->recursive_nsuccesses,1);
        gbl->recursive_successes = 0;
        gbl->recursive_fraction = 0;
    }
	
	input.getwdefault("ntstep",ntstep,1);
	if (ntstep < 0) {
		*gbl->log << "Number of time steps must be positive" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	input.getwdefault("restart",nstart,0);
	if (nstart < 0) {
		*gbl->log << "Restart must be positive" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	
	ntstep += nstart +1;
	gbl->time = 0.0;
	if (gbl->dti > 0.0) gbl->time = nstart/gbl->dti;

	/* OTHER UNIVERSAL CONSTANTS */
	input.getwdefault("gravity",gbl->g,0.0);
    double bodydflt[2] = {0.0,0.0};
    input.getwdefault("body_force",gbl->body.data(),2,bodydflt);

#ifndef petsc
	/* LOAD CONSTANTS FOR MULTISTAGE ITERATIVE SCHEME */
	input.getwdefault("nstage",gbl->nstage,5);
	if (gbl->nstage < 0) {
		*gbl->log << "Number of stages must be positive" << std::endl;
	}
    gbl->alpha.resize(gbl->nstage+1);
    gbl->beta.resize(gbl->nstage+1);
    istringstream datastream;

    /* LOAD COEFFICIENTS FOR IMAGINARY TERMS */
    const FLT alpha_dflt[5] = {0.25, 1./6., .375, .5, 1.0};
    input.getwdefault("alpha",gbl->alpha.data(),gbl->nstage,alpha_dflt);
    gbl->alpha(gbl->nstage) = 1.0;

    const FLT beta_dflt[5] = {1.0, 0.0, 5./9., 0.0, 4./9.};
    input.getwdefault("beta",gbl->beta.data(),gbl->nstage,beta_dflt);
    gbl->beta(gbl->nstage) = 1.0;
#else
    gbl->nstage = 1;
    gbl->alpha.resize(gbl->nstage+1);
    gbl->beta.resize(gbl->nstage+1);
    gbl->alpha = 1.0;
    gbl->beta = 1.0;
#endif


	/* LOAD BASIC CONSTANTS FOR MULTIGRID */
	input.getwdefault("itercrsn",itercrsn,1);
	if (itercrsn < 0) {
		*gbl->log << "Number of iterations must be positive" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}

	input.getwdefault("iterrfne",iterrfne,0);
	if (iterrfne < 0) {
		*gbl->log << "Number of iterations must be positive" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}

	input.getwdefault("ncycle",ncycle,0);
	if (ncycle < 0) {
		*gbl->log << "Number of iterations must be positive" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}

	input.getwdefault("preconditioner_interval",prcndtn_intrvl,1);

	input.getwdefault("vwcycle",vw,1);
	if (vw < 1) {
		*gbl->log << "vwcycle must be positive" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}

	input.getwdefault("absolute_tolerance",absolute_tolerance,1.0e-12);

	input.getwdefault("relative_tolerance",relative_tolerance,-1.0);

	input.getwdefault("ngrid",ngrid,1);  // JUST TO HAVE IN LOG FILE
	ngrid = MAX(ngrid,1);
	
	input.getwdefault("extra_coarsest_levels",extra_coarsest_levels,0);

	input.getwdefault("extra_finest_levels",extra_finest_levels,0);
	mglvls = ngrid+extra_coarsest_levels+extra_finest_levels;
	
	input.getwdefault("error_control_level",error_control_level,-1);
	if (error_control_level > -1 && error_control_level != mglvls-1) {
		*gbl->log << "Error being controlled on other than coarsest mesh" << std::endl;
	}
	
	input.getwdefault("error_control_tolerance",error_control_tolerance,0.33);
	
	input.getwdefault("output_interval", out_intrvl,1);
	if (out_intrvl < 1) {
		*gbl->log << "Output interval must be positive" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}

	input.getwdefault("restart_interval",rstrt_intrvl,1);
	if (rstrt_intrvl < 1) {
		*gbl->log << "Restart interval must be positive" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}

	input.getwdefault("debug_output",debug_output,0);

	input.getwdefault("adapt_output",gbl->adapt_output,false);

	if (!input.get(idprefix+"_adapt",gbl->adapt_interval)) {
		input.getwdefault("adapt",gbl->adapt_interval,0);
	}
    input.getwdefault(idprefix+"_adaptable",gbl->adaptable,true);

	if (!input.get(idprefix + "_tolerance",gbl->tolerance)) {
		input.getwdefault("tolerance",gbl->tolerance,sqrt(2.0));
	}

	if (!input.get(idprefix + "_error_target",gbl->error_target)) {
		input.getwdefault("error_target",gbl->error_target,1.0e-4);
	}
	
	if (!input.get(idprefix+"_length_smoothing_steps",gbl->length_smoothing_steps)) {
		input.getwdefault("length_smoothing_steps",gbl->length_smoothing_steps,0);
	}
	
	input.getwdefault("rsdl_debug",gbl->rsdl_debug,0);
	
	input.getwdefault("jac_debug",gbl->jac_debug,0);

	/* LOAD FINE MESH INFORMATION */
	grd(0)->init(input,gbl);
	if (gbl->adapt_output) {
		output("before_matched",block::display,0);
	}
	grd(0)->matchboundaries();
	if (gbl->adapt_output) {
		output("matched0",block::display,0);
	}

	ostringstream nstr;
	for(int lvl=1;lvl<ngrid;++lvl) {
//#define OLDRECONNECT
#ifdef OLDRECONNECT
			grd(lvl)->init(*grd(lvl-1),multigrid_interface::multigrid,2.0);
#else
		FLT size_reduce = 1.0;
		if (lvl > 1) size_reduce = 2.0;
		grd(lvl)->init(*grd(lvl-1),multigrid_interface::multigrid,size_reduce);
#endif
		grd(lvl)->connect(*grd(lvl-1));
		if (gbl->adapt_output) {
			nstr.str("");
			nstr << lvl << flush;
			mystring = "matched" +nstr.str();
			output(mystring,block::display,lvl);
		}
	}

	return;
}



FLT block::cycle(int vw, int lvl, bool evaluate_preconditioner) {
	int vcount,extra_count = 0;
	int gridlevel,gridlevelp;
	FLT error = 0.0,maxerror = 0.0;
    
	/* THIS ALLOWS FOR EXTRA LEVELS FOR BOTTOM AND TOP GRID */
	/* SO I CAN DO P-MULTIGRID & ALGEBRAIC STUFF */
	gridlevel = MIN(MAX(lvl-extra_finest_levels,0),ngrid-1);
	gridlevelp = MIN(MAX(lvl-extra_finest_levels+1,0),ngrid-1);

	for (vcount=0;vcount<vw;++vcount) {

        if (evaluate_preconditioner) {
            int err = grd(gridlevel)->setup_preconditioner();
            int errorsum = 0;
            sim::blks.allreduce(&err,&errorsum, 1, blocks::int_msg, blocks::sum);
            if (errorsum)
                throw 1;
        }
        
		for(int iter=0;iter<itercrsn;++iter)
			grd(gridlevel)->update();
		
		if (lvl+vcount == 0) {
			*gbl->log << gbl->iteration << ' ';
			error = maxres(0);
			*gbl->log << std::endl;
            /* Check for nan's */
            if (!(error > 0.0))
                throw 1;
		}
    

		/* THIS IS TO RUN A TWO-LEVEL ITERATION */
		/* OR TO CONVERGE THE SOLUTION ON THE COARSEST MESH */
		if (error_control_level == lvl) {
			*gbl->log << "#error_control " << extra_count << ' ';
			error = maxres(gridlevel);
            /* Check for nan's */
            if (!(error > 0.0))
                throw 1;
			maxerror = MAX(error,maxerror);
			*gbl->log << ' ' << error/maxerror << std::endl;
			if (error/maxerror > error_control_tolerance && error > absolute_tolerance && extra_count < ncycle) vcount = vw-2;
			if (debug_output) {
				if (extra_count % debug_output == 0) {
					std::string outname;
					std::ostringstream nstr("");
					nstr.str("");
					nstr << gbl->tstep << '_' << gbl->substep << '_' << extra_count << std::flush;
					outname = "coarse_debug" +nstr.str();
					output(outname,block::debug,gridlevel);
				}
			}
			++extra_count;
			if (lvl == mglvls-1) continue;
		}

		if (lvl == mglvls-1) return(error);

		grd(gridlevel)->rsdl();

		grd(gridlevelp)->mg_restrict();

		cycle(vw, lvl+1, evaluate_preconditioner);

		grd(gridlevelp)->mg_prolongate();

        if (evaluate_preconditioner && iterrfne) {
            int err = grd(gridlevel)->setup_preconditioner();
            int errorsum = 0;
            sim::blks.allreduce(&err,&errorsum, 1, blocks::int_msg, blocks::sum);
            if (errorsum)
                throw 1;
        }

		for(int iter=0;iter<iterrfne;++iter)
			grd(gridlevel)->update();
	}

	return(error);
}

void block::go(input_map input) {
	std::string outname;
	std::ostringstream nstr;
	clock_t begin_time, end_time;

	init(input);

    begin_time = clock();

#ifdef METIS
	/* partition utility */
	int nparts;
	if (input.get("partition",nparts)) {
		const int nblock = sim::blks.nblock;
		Array<int,1> part_list(nparts*nblock+2);
		grd(0)->setpartition(nparts,part_list);
		
		Array<int,1> part_count(nblock+1);
		part_count(0) = 0;
		for(int i=0;i<nblock;++i) {
			part_count(i+1) = part_count(i);
			for(int j=0;j<nparts;++j) {
				part_count(i+1) += part_list(i*nparts+j);
			}
		}
		
		int myparts = part_count(gbl->idnum);
		for(int i=0;i<nparts;++i) {
			if (part_list(nparts*gbl->idnum +i)) {
				multigrid_interface *apart = getnewlevel(input);
				apart->init(*grd(0),multigrid_interface::adapt_storage,1.0);
				apart->partition(*grd(0),myparts,part_list(nparts*nblock),part_list(nparts*nblock+1));
				ostringstream nstr;
				nstr << "b" << myparts++ << std::flush;
				*gbl->log << nstr.str() << "_matches: " << gbl->idprefix << std::endl;
				std::string old_idprefix = gbl->idprefix;
				gbl->idprefix = nstr.str();
				nstr.clear();
				apart->output("partition",block::display);
				apart->output("partition",block::restart);
				gbl->idprefix = old_idprefix;
			}
		}
		return;
	}
#endif
    
	/* OUTPUT INITIAL CONDITION */
	if (nstart == 0) {
		// TEMPORARY FIX FOR INCONSITENCY BETWEEN BOUNDARY OUTPUT (BY NUMBER) & BLOCK OUTPUT (BY STRING)
		gbl->tstep = 0;
		output("data0",block::display);
	}
    
    int rb2;
    if (input.get("refineby2",rb2)) {
        gbl->tstep=nstart+1;
        grd(0)->refineby2();
        nstr.str("");
        nstr << nstart+1 << std::flush;
        outname = "data" +nstr.str();
        output(outname,block::display);
        
        outname = "rstrt" +nstr.str();
        output(outname,block::restart);
        return;
    }
    
    int offset;
    if (input.get("offset",offset)) {
        gbl->tstep=nstart+1;
        grd(0)->offset_geometry(input);
        return;
    }
    
	for(gbl->tstep=nstart+1;gbl->tstep<ntstep;++gbl->tstep) {
        if (gbl->recursive_timestep_levels) {
            gbl->recursive_fraction = 0;
            recursive_time_step();
        }
        else {
            time_step();
        }
        
		/* OUTPUT DISPLAY FILES */
		if (!((gbl->tstep)%out_intrvl)) {
			nstr.str("");
			nstr << gbl->tstep << std::flush;
			outname = "data" +nstr.str();
			output(outname,block::display);
		}

		/* ADAPT MESH */
		if (gbl->adapt_interval && !(gbl->tstep%gbl->adapt_interval)) {
			adapt();
		}

		/* OUTPUT RESTART FILES */
		if (!((gbl->tstep)%(rstrt_intrvl))) {
			outname = "rstrt" +nstr.str();
			output(outname,block::restart);
            if (gbl->auto_timestep_tries) {
                gbl->dti = (gbl->dti/gbl->auto_timestep_ratio>gbl->auto_dti_min) ? gbl->dti/gbl->auto_timestep_ratio : gbl->auto_dti_min;
                *gbl->log << "Setting time step to " << gbl->dti << '\n';
            }
        }
        
        if (gbl->auto_timestep_tries) {
            if (gbl->auto_timestep_maxtime > 0.0 && gbl->time > gbl->auto_timestep_maxtime) {
                *gbl->log << "Auto timestepping reached the maximum time specified\n";
                break;
            }
        }
	}
	end_time = clock();
	*gbl->log << "that took " << static_cast<double>((end_time - begin_time))/ CLOCKS_PER_SEC << " cpu time" << std::endl;

	return;
}

void block::recursive_time_step(int this_level) {
    std::string outname;
    std::ostringstream nstr;
    FLT maxerror,error;
    
    if (gbl->recursive_level != this_level) {
        /* Jump directly to lower level */
        for(int halfstep_count=0;halfstep_count<2;++halfstep_count) {
            recursive_time_step(this_level+1);
        }
        return;
    }
    
    try {
        int denom = 1<<this_level;
        for(gbl->substep=0;gbl->substep<gbl->stepsolves;++gbl->substep) {
            *gbl->log << "#TIMESTEP: " << gbl->tstep << " " << (gbl->recursive_fraction>>(gbl->recursive_timestep_levels-this_level))+1 << '/' << denom << " SUBSTEP: " << gbl->substep << std::endl;
            tadvance();
            
            maxerror = 0.0;
            for(gbl->iteration=0;gbl->iteration<ncycle;++gbl->iteration) {
                error = cycle(vw,0,!(gbl->iteration%prcndtn_intrvl));
                maxerror = MAX(error,maxerror);
                
                if (debug_output) {
                    if (gbl->iteration % debug_output == 0) {
                        nstr.str("");
                        nstr << gbl->tstep << '_' << gbl->substep << '_' << gbl->iteration << std::flush;
                        outname = "debug" +nstr.str();
                        output(outname,block::debug);
                    }
                }
                if (error/maxerror < relative_tolerance || error < absolute_tolerance) break;
            }
            if (gbl->iteration == ncycle) {
                throw 1;
            }
        }
        ++gbl->recursive_successes;
        gbl->recursive_fraction += 1<<(gbl->recursive_timestep_levels -this_level);
    }
    catch(int e) {
        *gbl->log << "#Convergence error for time step " << gbl->tstep << " at level " << this_level << '\n';
        
        if (this_level < gbl->recursive_timestep_levels) {
            gbl->recursive_successes = 0;
            ++gbl->recursive_level;
            reset_timestep();
            gbl->dti *= 2.0;
            *gbl->log << "#Setting recursive time step to " << gbl->dti << '\n';
            
            /* Keep track if the sublevel actually successfully converged */
            for(int halfstep_count=0;halfstep_count<2;++halfstep_count) {
                recursive_time_step(this_level+1);
            }
            if (gbl->recursive_level == this_level+1 && gbl->recursive_successes == 2*gbl->recursive_nsuccesses) {
                gbl->recursive_successes = 0;
                --gbl->recursive_level;
                gbl->dti /= 2.0;
            }
        }
        else {
            *gbl->log << "#convergence error at finest time step\n";
            sim::abort(__LINE__,__FILE__,&std::cerr);
        }
    }

}

void block::time_step() {
    std::string outname;
    std::ostringstream nstr;
    FLT maxerror,error;

    for (int auto_timestep_failures = 0; auto_timestep_failures < gbl->auto_timestep_tries+1; ++auto_timestep_failures) {
        try {
            for(gbl->substep=0;gbl->substep<gbl->stepsolves;++gbl->substep) {
                *gbl->log << "#TIMESTEP: " << gbl->tstep << " SUBSTEP: " << gbl->substep << std::endl;
                tadvance();

                maxerror = 0.0;
                for(gbl->iteration=0;gbl->iteration<ncycle;++gbl->iteration) {
                    error = cycle(vw,0,!(gbl->iteration%prcndtn_intrvl));
                    maxerror = MAX(error,maxerror);

                    if (debug_output) {
                        if (gbl->iteration % debug_output == 0) {
                            nstr.str("");
                            nstr << gbl->tstep << '_' << gbl->substep << '_' << gbl->iteration << std::flush;
                            outname = "debug" +nstr.str();
                            output(outname,block::debug);
                        }
                    }
                    if (error/maxerror < relative_tolerance || error < absolute_tolerance) break;
                }
                if (gbl->auto_timestep_tries && gbl->iteration == ncycle) {
                    throw 1;
                }
            }
        }
        catch(int e) {
            *gbl->log << "#Convergence error for time step " << gbl->tstep << '\n';
            if (gbl->auto_timestep_tries) {
                reset_timestep();
                gbl->dti = gbl->auto_timestep_ratio*gbl->dti;
                *gbl->log << "#Setting time step to " << gbl->dti << '\n';
                if (gbl->dti > gbl->auto_dti_max) {
                    *gbl->log << "dti is too large. Aborting" << std::endl;
                    sim::abort(__LINE__,__FILE__,&std::cerr);
                }
                continue;
            }
            else {
                sim::abort(__LINE__,__FILE__,&std::cerr);
            }
        }
        return;
    }
    sim::abort(__LINE__,__FILE__,&std::cerr);
}


FLT block::maxres(int lvl) {
	FLT error = 0.0,globalmax;

	error = grd(lvl)->maxres();
	sim::blks.allreduce(&error,&globalmax, 1, blocks::flt_msg, blocks::max);

	return(globalmax);
}


void block::tadvance() {
	int lvl;
    
    if (gbl->substep == 0) {
        /* create backup of dti_prev in case we need to reset the time step */
        *gbl->log << "#time " << gbl->time << ", dti " << gbl->dti << ", dti_prev " << gbl->dti_prev << std::endl;
        gbl->dti_prev_store = gbl->dti_prev;
    }
    
    switch(gbl->time_scheme) {
        case(block_global::time_schemes::DIRK4): {
            switch(gbl->tstep) {     /* STARTUP SEQUENCE */
                case(1): {
                    /* THIS IS THE STANDARD FORM */
                    // FLT gbl->adirk(DIRK,DIRK) = {{GRK3,0.0,0.0},{C2RK3-GRK3,GRK3,0.0},{1-B2RK3-GRK3,B2RK3,GRK3}}
                    // FLT gbl->bdirk[DIRK] = {1-B2RK3-GRK3,B2RK3,GRK3};
                    // FLT gbl->cdirk(DIRK) = {GRK3,C2RK3,1.0};
                    /* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
                    gbl->adirk(0,0) = 1./sim::GRK3;                 gbl->adirk(0,1) = 0.0;             gbl->adirk(0,2) = 0.0;
                    gbl->adirk(1,0) = sim::C2RK3-sim::GRK3;      gbl->adirk(1,1) = 1./sim::GRK3; gbl->adirk(1,2) = 0.0;
                    gbl->adirk(2,0) = 1.-sim::B2RK3-sim::C2RK3; gbl->adirk(2,1) = sim::B2RK3;    gbl->adirk(2,2) = 1./sim::GRK3;
                    gbl->cdirk(0) = sim::GRK3; gbl->cdirk(1) = sim::C2RK3-sim::GRK3; gbl->cdirk(2) = 1.0-sim::C2RK3;
                    gbl->esdirk = false;
                    break;
                }
                case(2): {
                    gbl->adirk(0,0) = 1./sim::GRK3;                         gbl->adirk(0,1) = 0.0;             gbl->adirk(0,2) = 0.0;      gbl->adirk(0,3) = 0.0;
                    gbl->adirk(1,0) = sim::GRK4;                             gbl->adirk(1,1) = 1./sim::GRK4;        gbl->adirk(1,2) = 0.0;      gbl->adirk(1,3) = 0.0;
                    gbl->adirk(2,0) = sim::C3RK4-sim::A32RK4-2.*sim::GRK4;        gbl->adirk(2,1) = sim::A32RK4;         gbl->adirk(2,2) = 1./sim::GRK4; gbl->adirk(2,3) = 0.0;
                    gbl->adirk(3,0) = sim::B1RK4-(sim::C3RK4-sim::A32RK4-sim::GRK4); gbl->adirk(3,1) = sim::B2RK4-sim::A32RK4; gbl->adirk(3,2) = sim::B3RK4;    gbl->adirk(3,3) = 1./sim::GRK4;
                    gbl->cdirk(0) = 2.*sim::GRK4; gbl->cdirk(1) = sim::C3RK4-2.*sim::GRK4; gbl->cdirk(2) = 1.0-sim::C3RK4;
                    gbl->esdirk = true;
                    break;
                }
                default : {
                    /* THIS IS THE STANDARD FORM */
                    // FLT gbl->adirk(DIRK,DIRK) = {{0.0,0.0,0.0,0.0},{GRK4,GRK4,0.0,0.0},{C3RK4-A32RK4-GRK4,A32RK4,GRK4,0.0},{B1RK4,B2RK4,B3RK4,GRK4}}
                    // FLT gbl->bdirk[DIRK] = {B1RK4,B2RK4,B3RK4,GRK4};
                    // FLT gbl->cdirk(DIRK) = {0.0,2.*GRK4,C3RK4,1.0};
                    /* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
                    gbl->adirk(0,0) = 1./sim::GRK4;                         gbl->adirk(0,1) = 0.0;             gbl->adirk(0,2) = 0.0;      gbl->adirk(0,3) = 0.0;
                    gbl->adirk(1,0) = sim::GRK4;                             gbl->adirk(1,1) = 1./sim::GRK4;        gbl->adirk(1,2) = 0.0;      gbl->adirk(1,3) = 0.0;
                    gbl->adirk(2,0) = sim::C3RK4-sim::A32RK4-2.*sim::GRK4;        gbl->adirk(2,1) = sim::A32RK4;         gbl->adirk(2,2) = 1./sim::GRK4; gbl->adirk(2,3) = 0.0;
                    gbl->adirk(3,0) = sim::B1RK4-(sim::C3RK4-sim::A32RK4-sim::GRK4); gbl->adirk(3,1) = sim::B2RK4-sim::A32RK4; gbl->adirk(3,2) = sim::B3RK4;    gbl->adirk(3,3) = 1./sim::GRK4;
                    gbl->cdirk(0) = 2.*sim::GRK4; gbl->cdirk(1) = sim::C3RK4-2.*sim::GRK4; gbl->cdirk(2) = 1.0-sim::C3RK4;
                    gbl->esdirk = true;
                }
            }
            gbl->bd(0) = gbl->dti*gbl->adirk(gbl->substep +gbl->esdirk,gbl->substep +gbl->esdirk);
            if (gbl->dti > 0.0) {
                gbl->time += gbl->cdirk(gbl->substep)/gbl->dti;
                if (gbl->esdirk) {
                    /* ALLOWS CHANGES OF TIME STEP BETWEEN RESTARTS */
                    gbl->adirk(0,0) *= gbl->dti_prev/gbl->dti;
                }
                gbl->dti_prev = gbl->dti;
            }
            else {
                gbl->time = gbl->tstep;
            }
            break;
        }
            
        case(block_global::time_schemes::DIRK3): {
            /* THIS IS THE STANDARD FORM */
            // FLT gbl->adirk(DIRK,DIRK) = {{GRK3,0.0,0.0},{C2RK3-GRK3,GRK3,0.0},{1-B2RK3-GRK3,B2RK3,GRK3}}
            // FLT gbl->bdirk[DIRK] = {1-B2RK3-GRK3,B2RK3,GRK3};
            // FLT gbl->cdirk(DIRK) = {GRK3,C2RK3,1.0};
            /* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
            gbl->adirk(0,0) = 1./sim::GRK3;                 gbl->adirk(0,1) = 0.0;             gbl->adirk(0,2) = 0.0;
            gbl->adirk(1,0) = sim::C2RK3-sim::GRK3;      gbl->adirk(1,1) = 1./sim::GRK3; gbl->adirk(1,2) = 0.0;
            gbl->adirk(2,0) = 1.-sim::B2RK3-sim::C2RK3; gbl->adirk(2,1) = sim::B2RK3;    gbl->adirk(2,2) = 1./sim::GRK3;
            gbl->cdirk(0) = sim::GRK3; gbl->cdirk(1) = sim::C2RK3-sim::GRK3; gbl->cdirk(2) = 1.0-sim::C2RK3;
            gbl->esdirk = false;
            
            gbl->bd(0) = gbl->dti*gbl->adirk(gbl->substep,gbl->substep);
            if (gbl->dti > 0.0) {
                gbl->time += gbl->cdirk(gbl->substep)/gbl->dti;
                gbl->dti_prev = gbl->dti;
            }
            else {
                gbl->time = gbl->tstep;
            }
            break;
        }
            
        case(block_global::time_schemes::DIRK2): {
            /* STARTUP SEQUENCE */
            switch(gbl->tstep) {
                case(1): {
                    /* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
                    gbl->adirk(0,0) = 1.;                 gbl->adirk(0,1) = 0.0;
                    gbl->adirk(1,0) = 0.;                 gbl->adirk(1,1) = 0.;
                    gbl->cdirk(0) = 1;
                    gbl->esdirk = false;
                    break;
                }
                default: {
                    gbl->adirk(0,0) = 1.0;                        gbl->adirk(0,1) = 0.0;
                    gbl->adirk(1,0) = 0.0;                        gbl->adirk(1,1) = 1.0;
                    gbl->cdirk(0) = 1;
                    gbl->esdirk = true;
                    break;
                }
            }
            gbl->bd(0) = gbl->dti*gbl->adirk(gbl->substep +gbl->esdirk,gbl->substep +gbl->esdirk);
            if (gbl->dti > 0.0) {
                gbl->time += gbl->cdirk(gbl->substep)/gbl->dti;
                if (gbl->esdirk) {
                    /* ALLOWS CHANGES OF TIME STEP BETWEEN RESTARTS */
                    gbl->adirk(0,0) *= gbl->dti_prev/gbl->dti;
                }
                gbl->dti_prev = gbl->dti;
            }
            else {
                gbl->time = gbl->tstep;
            }
            break;
        }
        
        case(block_global::time_schemes::DIRK1): {
            gbl->adirk(0,0) = 1.0;
            gbl->cdirk(0) = 1.0;
            gbl->esdirk = false;
            gbl->bd(0) = gbl->dti*gbl->adirk(gbl->substep,gbl->substep);
            if (gbl->dti > 0.0) {
                gbl->time += gbl->cdirk(gbl->substep)/gbl->dti;
                gbl->dti_prev = gbl->dti;
            }
            else {
                gbl->time = gbl->tstep;
            }
            break;
        }
            
        case(block_global::time_schemes::AM1): {
            /* STARTUP SEQUENCE */
            switch(gbl->tstep) {
                case(1): {
                    /* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
                    gbl->adirk(0,0) = 1.;                 gbl->adirk(0,1) = 0.0;
                    gbl->adirk(1,0) = 0.;                 gbl->adirk(1,1) = 0.;
                    gbl->cdirk(0) = 1;
                    gbl->esdirk = false;
                    break;
                }
                case(2): {
                    gbl->adirk(0,0) = 1.0;                        gbl->adirk(0,1) = 0.0;
                    gbl->adirk(1,0) = 0.5;                        gbl->adirk(1,1) = 2.0;
                    gbl->cdirk(0) = 1;
                    gbl->esdirk = true;
                    break;
                }
                default : {
                    /* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
                    gbl->adirk(0,0) = 2.0;                        gbl->adirk(0,1) = 0.0;
                    gbl->adirk(1,0) = 0.5;                        gbl->adirk(1,1) = 2.0;
                    gbl->cdirk(0) = 1;
                    gbl->esdirk = true;
                    break;
                }
            }
            gbl->bd(0) = gbl->dti*gbl->adirk(gbl->substep +gbl->esdirk,gbl->substep +gbl->esdirk);
            if (gbl->dti > 0.0) {
                gbl->time += gbl->cdirk(gbl->substep)/gbl->dti;
                if (gbl->esdirk) {
                    /* ALLOWS CHANGES OF TIME STEP BETWEEN RESTARTS */
                    gbl->adirk(0,0) *= gbl->dti_prev/gbl->dti;
                }
                gbl->dti_prev = gbl->dti;
            }
            else {
                gbl->time = gbl->tstep;
            }
            break;
        }
            
        case(block_global::time_schemes::BD3): {
            switch(gbl->tstep) {
                case(1): {
                    gbl->bd(0) =  gbl->dti;
                    gbl->bd(1) = -gbl->dti;
                    break;
                }
                case(2): {
                    gbl->bd(0) =  1.5*gbl->dti;
                    gbl->bd(1) = -2.0*gbl->dti;
                    gbl->bd(2) =  0.5*gbl->dti;
                    break;
                }
                default: {
                    gbl->bd(0) = 11./6*gbl->dti;
                    gbl->bd(1) = -3.*gbl->dti;
                    gbl->bd(2) = 1.5*gbl->dti;
                    gbl->bd(3) = -1./3.*gbl->dti;
                    break;
                }
            }
            if (gbl->dti > 0.0)
                gbl->time += 1./gbl->dti;
            else
                gbl->time = gbl->tstep;
            break;
        }
            
        case(block_global::time_schemes::BD2): {
            switch(gbl->tstep) {
                case(1): {
                    gbl->bd(0) =  gbl->dti;
                    gbl->bd(1) = -gbl->dti;
                    break;
                }
                default: {
                    gbl->bd(0) =  1.5*gbl->dti;
                    gbl->bd(1) = -2.0*gbl->dti;
                    gbl->bd(2) =  0.5*gbl->dti;
                    break;
                }
            }
            if (gbl->dti > 0.0)
                gbl->time += 1./gbl->dti;
            else
                gbl->time = gbl->tstep;
            break;
        }
            
        case(block_global::time_schemes::BD1): {
            switch(gbl->tstep) {
                case(1): {
                    gbl->bd(0) =  gbl->dti;
                    gbl->bd(1) = -gbl->dti;
                    break;
                }
                default: {
                    gbl->bd(0) =  1.5*gbl->dti;
                    gbl->bd(1) = -2.0*gbl->dti;
                    gbl->bd(2) =  0.5*gbl->dti;
                    break;
                }
            }
            if (gbl->dti > 0.0)
                gbl->time += 1./gbl->dti;
            else
                gbl->time = gbl->tstep;
            break;
        }
	}

	for (lvl=0;lvl<ngrid;++lvl) {
		grd(lvl)->tadvance();
	}

	grd(0)->matchboundaries();

	return;
}

void block::reset_timestep() {

    for(int substep = gbl->substep; substep >= 0; --substep)
        gbl->time -= gbl->cdirk(substep)/gbl->dti;
    
    gbl->dti_prev = gbl->dti_prev_store;

    grd(0)->reset_timestep();
    
    return;
}

void block::output(const std::string &filename, output_purpose why, int level) {
	grd(level)->output(filename,why);
}

void block::adapt() {
	int lvl;

	grd(0)->matchboundaries();
	grd(0)->adapt();

	for(lvl=1;lvl<ngrid;++lvl) {
		grd(lvl)->connect(*grd(lvl-1));
	}

	return;
}

block::~block() {
	for(int i=ngrid-1;i>=0;--i) {
		delete grd(i);
	}
    delete gbl->log;
}
