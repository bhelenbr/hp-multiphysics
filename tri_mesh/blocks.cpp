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
#include "boundary.h"
#include <time.h>
#include <input_map.h>
#include <utilities.h>
#include <iostream>
#include <string>
#include <sstream>
#include <utilities.h>
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

blocks sim::blks;

#ifdef PTH
struct gostruct {
    block *blk;
    input_map input;
};

void* thread_go(void* ptr) {
    gostruct *p = static_cast<gostruct *>(ptr);
    p->blk->go(p->input);
    return 0;
}
#endif

void my_new_handler() 
{ 
	std::cerr << "Out of memory" << endl; 
	std::cerr.flush(); 
	std::abort(); 
} 

void blocks::go(const std::string &infile, const std::string &outfile) {
    input_map maptemp;
    std::string name;
    
    name = infile + ".inpt";
    maptemp.input(name);
    
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
			exit(1);
		}
		bstart += nb;
	}
	if (!(data >> myblock)) {
		std::cerr << "error reading blocks\n"; 
		exit(1);
	}
	data.clear();
	
	nblock = bstart +myblock;
	for (i=myid+1;i<nproc;++i) {
		if (!(data >> nb)) {
			std::cerr << "error reading blocks\n"; 
			exit(1);
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
				exit(1);
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
    if (nblock != 1) {
        std::cerr << "Need pth or boost::threads to run mulitple blocks on single processor\n";
        exit(1);
    }
    blk(0)->go(input);
#endif

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
		std::cerr << "couldn't find group in all_reduce\n";
		exit(1);
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
   
/** This is a data structure for storing communication info
 *  only used temporarily to sort things out in findmatch
 */
class block::comm_info {
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

int block::comm_info::unpack(blitz::Array<int,1> entitylist) {
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

void block::findmatch(int grdlvl) {

    *gbl->log << "# finding matches at multigrid level " << grdlvl << " for block " << idnum << std::endl;
    
    /* GET DATA FROM BLOCKS IN A THREAD SAFE WAY */
#if defined(PTH)
    pth_mutex_acquire(&sim::blks.data_mutex,false,NULL);
#elif defined(BOOST)
    boost::mutex::scoped_lock datalock(sim::blks.data_mutex);
#endif       
    int nblock = sim::blks.nblock;
    int myid = sim::blks.myid;
    int myblock = sim::blks.myblock;
    /* FIGURE OUT MY LOCAL BLOCK NUMBER */
    int b1;
    for (b1=0;b1<myblock;++b1)
        if (sim::blks.blk(b1) == this) break;
    if (b1 >= myblock) {
        *gbl->log << "Didn't find myself in block list?\n";
        exit(1);
    }
#if defined(PTH)
    pth_mutex_release(&sim::blks.data_mutex);
#elif defined(BOOST)
    datalock.unlock();
#endif   
    
    /* FIRST DETERMINE TOTAL SIZE OF LIST */
    Array<int,1> sndsize(nblock), size(nblock);
    sndsize = 0;  
    sndsize(idnum) = 2 +grd(grdlvl)->comm_entity_size(); // 2 is for myid & local block number
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
    grd(grdlvl)->comm_entity_list(sublist);
    ~sublist;
        
    Array<int,1> entitylist(tsize);
    entitylist = 0;
    sim::blks.allreduce(sndentitylist.data(),entitylist.data(),tsize,blocks::int_msg,blocks::sum);
    ~sndentitylist;
        
    /* UNPACK ENTITY LIST INTO MORE USABLE FORM */
    int max_comm_bdry_num = 0; // SO I CAN GENERATE UNIQUE TAGS
    count = 0;
    block::comm_info *binfo = new block::comm_info[nblock];
    for(int i=0;i<nblock;++i) {
        count += binfo[i].unpack(entitylist(Range(count,toEnd)));
        max_comm_bdry_num = MAX(max_comm_bdry_num,binfo[i].nvcomm);
        max_comm_bdry_num = MAX(max_comm_bdry_num,binfo[i].nscomm);
        max_comm_bdry_num = MAX(max_comm_bdry_num,binfo[i].nfcomm);
    }
    /* CALCULATE NUMBER OF BINARY DIGITS NEEDED TO MAKE TAGS */      
    sim::blks.setdigits(max_comm_bdry_num);
    
    /* OUTPUT LIST FOR DEBUGGING */
    *gbl->log << "# block " << idprefix << " block data" << std::endl;
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
        *gbl->log << "#\tvertex " << binfo[b1].vcomm[i].nvbd << " b1: " << binfo[b1].vcomm[i].idnum << std::endl;
        for(int b2=0;b2<nblock;++b2) {
            for(int j=0;j<binfo[b2].nvcomm;++j) {
                if (binfo[b1].vcomm[i].idnum == binfo[b2].vcomm[j].idnum) {
                    if (b1 == b2 && i == j) {
                        if (!first_found) first_found = true;  // Leave first flag alone
                        continue;  // CAN"T MATCH TO MYSELF
                    }
                    
                    boundary *bp1 = grd(grdlvl)->getvbdry(binfo[b1].vcomm[i].nvbd);
                    if (binfo[b1].proc == binfo[b2].proc) {
                        /* local match */
                        boundary *bp2 = sim::blks.blk(binfo[b2].blk)->grd(grdlvl)->getvbdry(binfo[b2].vcomm[j].nvbd);
                        bp1->local_cnnct(bp2,sim::blks.tagid(1,b1,b2,i,j),sim::blks.tagid(1,b2,b1,j,i));
                        if (!first_found) {
                            bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
                            first_found = true;
                        }
                        *gbl->log <<  "#\t\tlocal match to block: " << b2 << " pnt: " << binfo[b2].vcomm[j].nvbd << " tag: " << sim::blks.tagid(1,b1,b2,i,j) << ' ' << sim::blks.tagid(1,b2,b1,j,i) << " idnum: " << binfo[b2].vcomm[j].idnum << std::endl;
                    }
#ifdef MPISRC
                    else {
                        bp1->mpi_cnnct(binfo[b2].proc,sim::blks.tagid(1,b1,b2,i,j),sim::blks.tagid(1,b2,b1,j,i));
                        if (!first_found) {
                            bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
                            first_found = true;
                        }
                        *gbl->log <<  "#\t\t  mpi match to  block: " << b2 << " pnt: " << binfo[b2].vcomm[j].nvbd << " tag: " << sim::blks.tagid(1,b1,b2,i,j) << ' ' << sim::blks.tagid(1,b2,b1,j,i) << " idnum: " << binfo[b2].vcomm[j].idnum << std::endl;
                    }
#endif
                    
                }
            }
        }
    }
    
    /* LOOK FOR SIDE MATCHES */
    for(int i=0;i<binfo[b1].nscomm;++i) {
        bool first_found = false;
        *gbl->log << "#\tside " << binfo[b1].ecomm[i].nebd << " b1: " << binfo[b1].ecomm[i].idnum << std::endl;
        for(int b2=0;b2<nblock;++b2) {
            for(int j=0;j<binfo[b2].nscomm;++j) {
                if (binfo[b1].ecomm[i].idnum == binfo[b2].ecomm[j].idnum) {
                    if (b1 == b2 && i == j) {
                        if (!first_found) first_found = true;  // Leave first flag alone
                        continue;  // CAN"T MATCH TO MYSELF
                    }
                    
                    boundary *bp1 = grd(grdlvl)->getebdry(binfo[b1].ecomm[i].nebd);
                    if (binfo[b1].proc == binfo[b2].proc) {
                        /* local match */
                        boundary *bp2 = sim::blks.blk(binfo[b2].blk)->grd(grdlvl)->getebdry(binfo[b2].ecomm[j].nebd);
                        bp1->local_cnnct(bp2,sim::blks.tagid(2,b1,b2,i,j),sim::blks.tagid(2,b2,b1,j,i));
                        if (!first_found) {
                            bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
                            first_found = true;
                        }
                        *gbl->log <<  "#\t\tlocal match to block: " << b2 << " edge: " << binfo[b2].ecomm[j].nebd << " tag: " << sim::blks.tagid(2,b1,b2,i,j) << ' ' << sim::blks.tagid(2,b2,b1,j,i) << " idnum: " << binfo[b2].ecomm[j].idnum << std::endl;
                    }
#ifdef MPISRC
                    else {
                        bp1->mpi_cnnct(binfo[b2].proc,sim::blks.tagid(2,b1,b2,i,j),sim::blks.tagid(2,b2,b1,j,i));
                        if (!first_found) {
                            bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
                            first_found = true;
                        }
                        *gbl->log <<  "#\t\t  mpi match to  block: " << b2 << " edge: " << binfo[b2].ecomm[j].nebd << " tag: " << sim::blks.tagid(2,b1,b2,i,j) << ' ' << sim::blks.tagid(2,b2,b1,j,i) << " idnum: " << binfo[b2].ecomm[j].idnum << std::endl;
                    }
#endif
                    
                }
            }
        }
    }
    
            
    
    /* LOOK FOR FACE MATCHES */
    for(int i=0;i<binfo[b1].nfcomm;++i) {
        bool first_found = false;
        *gbl->log << "#\tface " << binfo[b1].fcomm[i].nfbd << " b1: " << binfo[b1].fcomm[i].idnum << std::endl;
        for(int b2=0;b2<nblock;++b2) {
            for(int j=0;j<binfo[b2].nfcomm;++j) {
                if (binfo[b1].fcomm[i].idnum == binfo[b2].fcomm[j].idnum) {
                    if (b1 == b2 && i == j) {
                        if (!first_found) first_found = true;  // Leave first flag alone
                        continue;  // CAN"T MATCH TO MYSELF
                    }
                    
                    boundary *bp1 = grd(grdlvl)->getfbdry(binfo[b1].fcomm[i].nfbd);
                    if (binfo[b1].proc == binfo[b2].proc) {
                        /* local match */
                        boundary *bp2 = sim::blks.blk(binfo[b2].blk)->grd(grdlvl)->getfbdry(binfo[b2].fcomm[j].nfbd);
                        bp1->local_cnnct(bp2,sim::blks.tagid(3,b1,b2,i,j),sim::blks.tagid(3,b2,b1,j,i));
                        if (!first_found) {
                            bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
                            first_found = true;
                        }
                        *gbl->log <<  "#\t\tlocal match to block: " << b2 << " face: " << binfo[b2].fcomm[j].nfbd << " tag: " << sim::blks.tagid(3,b1,b2,i,j) << ' ' << sim::blks.tagid(3,b2,b1,j,i) << " idnum: " << binfo[b2].fcomm[j].idnum << std::endl;
                    }
#ifdef MPISRC
                    else {
                        bp1->mpi_cnnct(binfo[b2].proc,sim::blks.tagid(3,b1,b2,i,j),sim::blks.tagid(3,b2,b1,j,i));
                        if (!first_found) {
                            bp1->is_frst() = !bp1->is_frst(); // Switches true to false by default
                            first_found = true;
                        }
                        *gbl->log <<  "#\t\t  mpi match to  block: " << b2 << " face: " << binfo[b2].fcomm[j].nfbd << " tag: " << sim::blks.tagid(3,b1,b2,i,j) << ' ' << sim::blks.tagid(3,b2,b1,j,i) << " idnum: " << binfo[b2].fcomm[j].idnum << std::endl;
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

    return;
}

void block::init(input_map &input) {
    std::string mystring;
    
    /* NEED TO BOOTSTRAP UNTIL I CAN GET LOGFILE OPENED */
    input.getwdefault("ngrid",ngrid,1);
    grd.resize(ngrid);
    for (int i=0;i<ngrid;++i) 
        grd(i) = getnewlevel(input);
    gbl = static_cast<block_global *>(grd(0)->create_global_structure());
	gbl->idnum = idnum;
    gbl->idprefix = idprefix;
    
    /* OPEN LOGFILES FOR EACH BLOCK */
    if (input.get("logfile",mystring)) {
        mystring += "." +idprefix +".log";
        std::ofstream *filelog = new std::ofstream;
        filelog->setf(std::ios::scientific, std::ios::floatfield);
        filelog->precision(3);
        filelog->open(mystring.c_str());
        gbl->log = filelog;
    }
    else {
        std::cout.setf(std::ios::scientific, std::ios::floatfield);
        std::cout.precision(3);
        gbl->log = &std::cout;
    }
    input.echo = true;
    input.log = gbl->log;
    input.echoprefix = "#";
    
#ifdef CAPRI
    int status = gi_uStart();
    *gbl->log << "gi_uStart status = ", status, "\n";
    if (status != CAPRI_SUCCESS) exit(1);
    
    if (input.get("BRep",mystring)) {
        status = gi_uLoadPart(mystring.c_str());
        *gbl->log << mystring << ": gi_uLoadPart status = " << status << std::endl;
        if (status != CAPRI_SUCCESS) exit(1);
    }
#endif 

    /* SET TIME STEPPING INFO */
    gbl->tstep = -1; // Simulation starts at t = 0, This is set negative until first tadvance to alllow change between initialization & B.C.'s
    gbl->substep = -1;
    input.getwdefault("dtinv",gbl->dti,0.0);
    input.getwdefault("dtinv_prev",gbl->dti_prev,gbl->dti); 
    input.getwdefault("ntstep",ntstep,1);
    input.getwdefault("restart",nstart,0);
    ntstep += nstart +1;
	gbl->time = 0.0;
    if (gbl->dti > 0.0) gbl->time = nstart/gbl->dti;
        
    /* OTHER UNIVERSAL CONSTANTS */
    input.getwdefault("gravity",gbl->g,0.0);
    
    /* LOAD BASIC CONSTANTS FOR MULTIGRID */
    input.getwdefault("itercrsn",itercrsn,1);
    
    input.getwdefault("iterrfne",iterrfne,0);
        
    input.getwdefault("ncycle",ncycle,0);
    
    input.getwdefault("preconditioner_interval",prcndtn_intrvl,-1);
    
    input.getwdefault("vwcycle",vw,2);
    
    input.getwdefault("absolute_tolerance",absolute_tolerance,1.0e-12);
    
    input.getwdefault("relative_tolerance",relative_tolerance,-1.0);
    
    input.getwdefault("error_control_level",error_control_level,-1);
    
    input.getwdefault("error_control_tolerance",error_control_tolerance,0.33);
    
    input.getwdefault("ngrid",ngrid,1);  // JUST TO HAVE IN LOG FILE

    input.getwdefault("extra_coarsest_levels",extra_coarsest_levels,0);
    
    input.getwdefault("extra_finest_levels",extra_finest_levels,0);
    mglvls = ngrid+extra_coarsest_levels+extra_finest_levels;
    
    input.getwdefault("output_interval", out_intrvl,1);
    
    input.getwdefault("restart_interval",rstrt_intrvl,1);
    rstrt_intrvl = MAX(1,rstrt_intrvl);
    
    input.getwdefault("debug_output",debug_output,false);
    
    input.getwdefault("adapt_output",gbl->adapt_output,false);
    
    if (!input.get(idprefix+"_adapt",gbl->adapt_flag)) {
        input.getwdefault("adapt",gbl->adapt_flag,false);
    }
    
    if (!input.get(idprefix + "_tolerance",gbl->tolerance)) {
        input.getwdefault("tolerance",gbl->tolerance,1.25);
    }   
    
    if (!input.get(idprefix + "_error_target",gbl->error_target)) {
        input.getwdefault("error_target",gbl->error_target,1.0e-4);
    }  
         
    /* LOAD FINE MESH INFORMATION */
    grd(0)->init(input,gbl);
    findmatch(0);
    grd(0)->matchboundaries();
    if (gbl->adapt_output) {
        output("matched0",block::display,0);
    }
    
    ostringstream nstr;
#define OLDRECONNECT
    for(int lvl=1;lvl<ngrid;++lvl) {
#ifdef OLDRECONNECT
        grd(lvl)->init(*grd(lvl-1),multigrid_interface::multigrid,2.0);
#else
        FLT size_reduce = 1.0;
        if (lvl > 1) size_reduce = 2.0;
        grd(lvl)->init(*grd(lvl-1),multigrid_interface::multigrid,size_reduce);
#endif
        findmatch(lvl);             // Because this is done second
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



void block::cycle(int vw, int lvl) {
    int vcount,extra_count = 0; 
    int gridlevel,gridlevelp;
    FLT error,maxerror = 0.0;
    
    /* THIS ALLOWS FOR EXTRA LEVELS FOR BOTTOM AND TOP GRID */
    /* SO I CAN DO P-MULTIGRID & ALGEBRAIC STUFF */
    gridlevel = MIN(MAX(lvl-extra_finest_levels,0),ngrid-1);
    gridlevelp = MIN(MAX(lvl-extra_finest_levels+1,0),ngrid-1);
        
    for (vcount=0;vcount<vw;++vcount) {
    
        if (!(vcount%abs(prcndtn_intrvl)) && itercrsn) grd(gridlevel)->setup_preconditioner();
                        
        for(int iter=0;iter<itercrsn;++iter)
            grd(gridlevel)->update();
        
        /* THIS IS TO RUN A TWO-LEVEL ITERATION */
        /* OR TO CONVERGE THE SOLUTION ON THE COARSEST MESH */
        if (error_control_level == lvl) {
            *gbl->log << "#error_control ";
            error = maxres(gridlevel);
            maxerror = MAX(error,maxerror);
            *gbl->log << ' ' << error/maxerror << std::endl;
            if (error/maxerror > error_control_tolerance) vcount = vw-2;
			if (debug_output) {
			    std::string outname;
				std::ostringstream nstr("");
				nstr.str("");
				nstr << gbl->tstep << '_' << gbl->substep << '_' << extra_count++ << std::flush;
				outname = "coarse_debug" +nstr.str();
				output(outname,block::debug,gridlevel);
			}
            if (lvl == mglvls-1) continue;
        }
        
        if (lvl == mglvls-1) return;
        
        grd(gridlevel)->rsdl();
        
        grd(gridlevelp)->mg_restrict();
        
        cycle(vw, lvl+1);

        grd(gridlevelp)->mg_prolongate();
        
        if (!(vcount%abs(prcndtn_intrvl)) && prcndtn_intrvl < 0 && iterrfne) grd(gridlevel)->setup_preconditioner();
        
        for(int iter=0;iter<iterrfne;++iter)
            grd(gridlevel)->update();
    }

    return;
}

void block::go(input_map input) {
    int i;
    std::string outname;
    std::ostringstream nstr;
    clock_t begin_time, end_time;
    FLT maxerror,error;

    
    init(input);
    
    begin_time = clock();

    /* OUTPUT INITIAL CONDITION */
    if (nstart == 0) output("data0",block::display);
    
    for(gbl->tstep=nstart+1;gbl->tstep<ntstep;++gbl->tstep) {
        for(gbl->substep=0;gbl->substep<sim::stepsolves;++gbl->substep) {
            *gbl->log << "#TIMESTEP: " << gbl->tstep << " SUBSTEP: " << gbl->substep << std::endl;
            tadvance();

            maxerror = 0.0;
            for(i=0;i<ncycle;++i) {
                cycle(vw);
                *gbl->log << i << ' ';
                error = maxres();
                maxerror = MAX(error,maxerror);
                *gbl->log << std::endl << std::flush;
                if (debug_output) {
                    nstr.str("");
                    nstr << gbl->tstep << '_' << gbl->substep << '_' << i << std::flush;
                    outname = "debug" +nstr.str();
                    output(outname,block::debug);
                }
                if (error/maxerror < relative_tolerance || error < absolute_tolerance) break;
            }
        }
        
        /* OUTPUT DISPLAY FILES */
        if (!((gbl->tstep)%out_intrvl)) {
            nstr.str("");
            nstr << gbl->tstep << std::flush;
            outname = "data" +nstr.str();
            output(outname,block::display);
        }
    
        /* ADAPT MESH */
        if (gbl->adapt_flag) adapt();
    
        /* OUTPUT RESTART FILES */
        if (!((gbl->tstep)%(rstrt_intrvl*out_intrvl))) {
            outname = "rstrt" +nstr.str();
            output(outname,block::restart);
        }
    }
    end_time = clock();
    *gbl->log << "that took " << static_cast<double>((end_time - begin_time)/ CLOCKS_PER_SEC) << " cpu time" << std::endl;
    
    return;
}

FLT block::maxres(int lvl) {
    FLT error = 0.0,globalmax;
    
    error = grd(lvl)->maxres();
    sim::blks.allreduce(&error,&globalmax, 1, blocks::flt_msg, blocks::max);
    
    return(globalmax);
}


void block::tadvance() {
    int lvl;    

#ifdef BACKDIFF
    if (gbl->dti > 0.0) gbl->time += 1./gbl->dti;

    for(i=0;i<BACKDIFF+1;++i)
        gbl->bd[i] = 0.0;
    
    switch(gbl->tstep) {
        case(1):
            gbl->bd[0] =  gbl->dti;
            gbl->bd[1] = -gbl->dti;
            break;
#if (BACKDIFF > 1)
        case(2):
            gbl->bd[0] =  1.5*gbl->dti;
            gbl->bd[1] = -2.0*gbl->dti;
            gbl->bd[2] =  0.5*gbl->dti;
            break;
#endif
#if (BACKDIFF > 2)
        case(3):
            gbl->bd[0] = 11./6*gbl->dti;
            gbl->bd[1] = -3.*gbl->dti;
            gbl->bd[2] = 1.5*gbl->dti;
            gbl->bd[3] = -1./3.*gbl->dti;
            break;
#endif
    }
#endif

#ifdef DIRK
#if (DIRK == 4)
    /* STARTUP SEQUENCE */
    switch(gbl->tstep) {
        case(1): {
            /* THIS IS THE STANDARD FORM */
            // FLT gbl->adirk[DIRK][DIRK] = {{GRK3,0.0,0.0},{C2RK3-GRK3,GRK3,0.0},{1-B2RK3-GRK3,B2RK3,GRK3}} 
            // FLT gbl->bdirk[DIRK] = {1-B2RK3-GRK3,B2RK3,GRK3};
            // FLT gbl->cdirk[DIRK] = {GRK3,C2RK3,1.0};
            /* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
            gbl->adirk[0][0] = 1./sim::GRK3;                 gbl->adirk[0][1] = 0.0;             gbl->adirk[0][2] = 0.0;
            gbl->adirk[1][0] = sim::C2RK3-sim::GRK3;      gbl->adirk[1][1] = 1./sim::GRK3; gbl->adirk[1][2] = 0.0;
            gbl->adirk[2][0] = 1.-sim::B2RK3-sim::C2RK3; gbl->adirk[2][1] = sim::B2RK3;    gbl->adirk[2][2] = 1./sim::GRK3;
            gbl->cdirk[0] = sim::GRK3; gbl->cdirk[1] = sim::C2RK3-sim::GRK3; gbl->cdirk[2] = 1.0-sim::C2RK3;
            gbl->esdirk = false;
            break;
        }
        case(2): {
            gbl->adirk[0][0] = 1./sim::GRK3;                         gbl->adirk[0][1] = 0.0;             gbl->adirk[0][2] = 0.0;      gbl->adirk[0][3] = 0.0;
            gbl->adirk[1][0] = sim::GRK4;                             gbl->adirk[1][1] = 1./sim::GRK4;        gbl->adirk[1][2] = 0.0;      gbl->adirk[1][3] = 0.0;
            gbl->adirk[2][0] = sim::C3RK4-sim::A32RK4-2.*sim::GRK4;        gbl->adirk[2][1] = sim::A32RK4;         gbl->adirk[2][2] = 1./sim::GRK4; gbl->adirk[2][3] = 0.0;
            gbl->adirk[3][0] = sim::B1RK4-(sim::C3RK4-sim::A32RK4-sim::GRK4); gbl->adirk[3][1] = sim::B2RK4-sim::A32RK4; gbl->adirk[3][2] = sim::B3RK4;    gbl->adirk[3][3] = 1./sim::GRK4; 
            gbl->cdirk[0] = 2.*sim::GRK4; gbl->cdirk[1] = sim::C3RK4-2.*sim::GRK4; gbl->cdirk[2] = 1.0-sim::C3RK4;
            gbl->esdirk = true;
            break;
        }
        default : {
            /* THIS IS THE STANDARD FORM */
            // FLT gbl->adirk[DIRK][DIRK] = {{0.0,0.0,0.0,0.0},{GRK4,GRK4,0.0,0.0},{C3RK4-A32RK4-GRK4,A32RK4,GRK4,0.0},{B1RK4,B2RK4,B3RK4,GRK4}} 
            // FLT gbl->bdirk[DIRK] = {B1RK4,B2RK4,B3RK4,GRK4};
            // FLT gbl->cdirk[DIRK] = {0.0,2.*GRK4,C3RK4,1.0};
            /* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
            gbl->adirk[0][0] = 1./sim::GRK4;                         gbl->adirk[0][1] = 0.0;             gbl->adirk[0][2] = 0.0;      gbl->adirk[0][3] = 0.0;
            gbl->adirk[1][0] = sim::GRK4;                             gbl->adirk[1][1] = 1./sim::GRK4;        gbl->adirk[1][2] = 0.0;      gbl->adirk[1][3] = 0.0;
            gbl->adirk[2][0] = sim::C3RK4-sim::A32RK4-2.*sim::GRK4;        gbl->adirk[2][1] = sim::A32RK4;         gbl->adirk[2][2] = 1./sim::GRK4; gbl->adirk[2][3] = 0.0;
            gbl->adirk[3][0] = sim::B1RK4-(sim::C3RK4-sim::A32RK4-sim::GRK4); gbl->adirk[3][1] = sim::B2RK4-sim::A32RK4; gbl->adirk[3][2] = sim::B3RK4;    gbl->adirk[3][3] = 1./sim::GRK4; 
            gbl->cdirk[0] = 2.*sim::GRK4; gbl->cdirk[1] = sim::C3RK4-2.*sim::GRK4; gbl->cdirk[2] = 1.0-sim::C3RK4;
            gbl->esdirk = true;
        }
    }
    gbl->bd[0] = gbl->dti*gbl->adirk[gbl->substep +gbl->esdirk][gbl->substep +gbl->esdirk];
    if (gbl->dti > 0.0) {
        gbl->time += gbl->cdirk[gbl->substep]/gbl->dti;
        if (gbl->esdirk) {
            /* ALLOWS CHANGES OF TIME STEP BETWEEN RESTARTS */
            gbl->   adirk[0][0] *= gbl->dti_prev/gbl->dti;
        }
        gbl->dti_prev = gbl->dti;
    }
    else {
        gbl->time = gbl->tstep;
    }
#endif
#endif


    for (lvl=0;lvl<ngrid;++lvl) {  
        grd(lvl)->tadvance();
    }
    
    grd(0)->matchboundaries();

    return;
}

void block::output(const std::string &filename, output_purpose why, int level) {
    std::string fapp;
    fapp = filename +"_" +idprefix; 
    grd(level)->output(fapp,why);
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
