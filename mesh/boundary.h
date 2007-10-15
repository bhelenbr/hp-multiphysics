/*
 *  boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Fri Jun 07 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include <utilities.h>
#include <stdio.h>
#include <input_map.h>
#include <mathclass.h>
#include <float.h>

#ifdef MPISRC
#include <mpi.h>
#endif

#ifndef _boundary_h_
#define _boundary_h_

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

/** \defgroup boundary Mesh Boundary Objects 
 *  This group contains object for dealing with mesh boundaries
 */
 
/** \brief Generic interface for a b.c. 
 *
 *  \ingroup boundary
 *  Contains all the basic functions for parallel communications
 */
class boundary {
    public:
        int idnum;
        std::string idprefix; 
        std::string mytype;

        boundary(int idin) : idnum(idin) {
            char buffer[100];
            std::string keyname;
            sprintf(buffer,"%d",idnum);
            idprefix = std::string(buffer);
            mytype = "boundary";
        }
        virtual void alloc(int n) {}
        virtual void output(std::ostream& fout) {
            fout << idprefix << "_type: " << mytype << std::endl;            
        }
        virtual void input(input_map& bdrydata) {}
        
        /* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
        enum msg_type {flt_msg, int_msg};
        enum groups {all,all_phased,partitions,manifolds};
        enum comm_type {symmetric,master_slave,slave_master};
        enum operation {average,sum,maximum,replace};
        union  {
            bool bdum;
            int idum;
            FLT fdum;
            msg_type mdum;
        } dummy;
        virtual bool is_comm() {return(false);}
        virtual bool& is_frst() {return(dummy.bdum=true);}
        virtual int& group() {return(dummy.idum=1);}
        virtual int matches() {return(0);}
        virtual int local_cnnct(boundary *bin, int snd_tag, int rcv_tag) {return 1;}
#ifdef MPISRC
        virtual int mpi_cnnct(int proc_tgt, int snd_tag, int rcv_tag) {return 1;}
#endif
        virtual int& matchphase(boundary::groups group, int matchnum) {return(dummy.idum=0);}
        virtual void resize_buffers(int size) {}
        virtual void *psndbuf() {return(&dummy);}
        virtual int& isndbuf(int indx) {return(dummy.idum);}
        virtual FLT& fsndbuf(int indx) {return(dummy.fdum);}
        virtual int& ircvbuf(int m,int indx) {return(dummy.idum);}
        virtual FLT& frcvbuf(int m,int indx) {return(dummy.fdum);}
        virtual int& sndsize() {return(dummy.idum=0);}
        virtual boundary::msg_type& sndtype() {return(dummy.mdum);}
        virtual void comm_prepare(boundary::groups group, int phase, comm_type type) {}
        virtual void comm_exchange(boundary::groups group, int phase, comm_type type) {}
        virtual int comm_wait(boundary::groups group, int phase, comm_type type) {return 1;}
        virtual int comm_nowait(boundary::groups group, int phase, comm_type type) {return 1;}
        virtual ~boundary() {}
};

/** \brief Template class to make a communciation boundary 
 *
 * \ingroup boundary
 * Contains variables and routines for communciation boundaries.
 * Can be applied to vrtx_bdry or edge_bdry to make specific type
 */
template<class BASE,class MESH> class comm_bdry : public BASE {
    protected:
        static const int maxmatch = 8;
        bool first; //!< For master-slave communication. Only one boundary of matching boundaries is first
        int maxgroup;
        int groupmask;  //!< To make groups that only communicate in restricted situations group 0 all, group 1 all phased, group 2 partitions, group 3 manifolds 
        Array<int,1> maxphase; //!<  For phased symmetric message passing for each group
        Array<TinyVector<int,maxmatch>,1> phase;  //!< To set-up staggered sequence of symmetric passes for each group (-1 means skip)

        int buffsize; //!< Size in bytes of buffer
        void *sndbuf; //!< Raw memory for outgoing message buffer
        Array<FLT,1> fsndbufarray; //!< Access to outgoing message buffer for floats
        Array<int,1> isndbufarray; //!< Access to outgoing message buffer for ints
        int msgsize; //!< Outgoing size
        boundary::msg_type msgtype; //!< Outgoing type
        
        /** Different types of matching boundaries,
         * local is same processor same thread
         * mpi is different processor
         * someday have threads but not yet
         */
#ifdef MPISRC
        enum matchtype {local, mpi};  
#else
        enum matchtype {local};
#endif
        int nmatch; //!< Number of matching boundaries
        TinyVector<matchtype,maxmatch> mtype; //!< Local or mpi or ?
        TinyVector<boundary *,maxmatch> local_match; //!< Pointers to local matches
        TinyVector<int,maxmatch> snd_tags; //!< Identifies each connection uniquely
        TinyVector<int,maxmatch> rcv_tags; //!< Identifies each connection uniquely
        TinyVector<void *,maxmatch> rcvbuf; //!< Raw memory to store incoming messages
        TinyVector<Array<FLT,1>,maxmatch> frcvbufarray; //!< Access to incoming message buffer for floats
        TinyVector<Array<int,1>,maxmatch> ircvbufarray; //!< Access to incoming message buffer for ints
        
#ifdef MPISRC
        TinyVector<int,maxmatch> mpi_match; //!< Processor numbers for mpi
        TinyVector<MPI_Request,maxmatch> mpi_rcvrqst; //!< Identifier returned from mpi to monitor success of recv
        TinyVector<MPI_Request,maxmatch> mpi_sndrqst; //!< Identifier returned from mpi to monitor success of send
#endif
                
    public:
        comm_bdry(int inid, MESH &xin) : BASE(inid,xin), first(1), maxgroup(1), groupmask(0x3), buffsize(0), nmatch(0) {
            maxphase.resize(maxgroup+1);
            phase.resize(maxgroup+1);
            for(int m=0;m<maxmatch;++m) phase(0)(m) = 0;
        }
        comm_bdry(const comm_bdry<BASE,MESH> &inbdry, MESH& xin) : BASE(inbdry,xin), first(inbdry.first), maxgroup(inbdry.maxgroup), groupmask(inbdry.groupmask), buffsize(0), nmatch(0) {    
            maxphase.resize(maxgroup+1);
            phase.resize(maxgroup+1);
            maxphase = inbdry.maxphase;
            for(int k=0;k<maxgroup+1;++k)
                phase(k) = inbdry.phase(k);
            
            /* COPY THESE, BUT WILL HAVE TO BE RESET TO NEW MATCHING SIDE */
            mtype = inbdry.mtype;
            local_match = local_match;
            snd_tags = inbdry.snd_tags;
            rcv_tags = inbdry.rcv_tags;

#ifdef MPISRC
            mpi_match = inbdry.mpi_match;
#endif
            return;
        }
              
        comm_bdry<BASE,MESH>* create(MESH &xin) const {return(new comm_bdry<BASE,MESH>(*this,xin));}
        bool is_comm() {return(true);}
        bool& is_frst() {return(first);}
        int& group() {return(groupmask);}
        int& sndsize() {return(msgsize);}
        boundary::msg_type& sndtype() {return(msgtype);}
        int matches() {return(nmatch);}
        int& matchphase(boundary::groups group, int matchnum) {return(phase(group)(matchnum));}  
        int& isndbuf(int indx) {return(isndbufarray(indx));}
        FLT& fsndbuf(int indx) {return(fsndbufarray(indx));}
        int& ircvbuf(int m,int indx) {return(ircvbufarray(m)(indx));}
        FLT& frcvbuf(int m,int indx) {return(frcvbufarray(m)(indx));}
        
        void resize_buffers(int nfloats) {
            if (buffsize) free(sndbuf);
            buffsize = nfloats*sizeof(FLT);
            sndbuf = xmalloc(buffsize); 
            Array<FLT,1> temp(static_cast<FLT *>(sndbuf), buffsize/sizeof(FLT), neverDeleteData);
            fsndbufarray.reference(temp);
            Array<int,1> temp1(static_cast<int *>(sndbuf), buffsize/sizeof(int), neverDeleteData);
            isndbufarray.reference(temp1);
        }
        
        void alloc(int nels) {
            BASE::alloc(nels);
            resize_buffers(nels*3);
        }

        void input(input_map& inmap) {
            int j,k,m;
            std::string keyword,val;
            std::map<std::string,std::string>::const_iterator mi;
            std::istringstream data;
            std::ostringstream nstr;
            
            BASE::input(inmap);
            
            keyword = BASE::idprefix +"_first";
            inmap.getwdefault(keyword,first,true);
            
            /* SET GROUP MEMBERSHIP FLAGS */
            maxgroup = 0;
            groupmask = 0;
            inmap.getlinewdefault(BASE::idprefix + "_group",val,"0 1"); // DEFAULT IS FIRST 2 GROUPS
            data.str(val);
            while(data >> m) {
                groupmask = groupmask|(1<<m);
                maxgroup = MAX(maxgroup,m);
            }
            data.clear();
                        
            /* LOAD PHASES */
            maxphase.resize(maxgroup+1);
            maxphase = 0;
            phase.resize(maxgroup+1);
            phase = 0;
            /* SKIP GROUP 0 BECAUSE THAT GROUP IS NOT PHASED */
            for(k=1;k<maxgroup+1;++k) {
                if (!(groupmask&(1<<k))) continue;
                
                nstr.str("");
                nstr << BASE::idprefix << "_phase" << k << std::flush;
                if (inmap.getline(nstr.str(),val)) {
                    data.str(val);
                    m = 0;
                    while(data >> j) {
                        phase(k)(m) = j;
                        maxphase(k) = MAX(maxphase(k),phase(k)(m));
                        ++m;
                    }
                    data.clear();
                }
            }
        }
        
        void output(std::ostream& fout) {
            BASE::output(fout);
            
            fout << BASE::idprefix << "_group" << ": ";
            for(int k=0;k<maxgroup+1;++k)
                if (groupmask&(1<<k)) fout << k << ' ';
            fout << std::endl;  
            
            for(int k=0;k<maxgroup+1;++k) {
                if (groupmask&(1<<k)) {
                    fout << BASE::idprefix << "_phase (not set yet so this is dumb)" << k << ": ";
                    for (int m=0;m<nmatch;++m)
                        fout << phase(k)(m) << " ";
                    fout << std::endl;
                }
            }
        }

        int local_cnnct(boundary *bin, int snd_tag, int rcv_tag) {
            if (bin->idnum == BASE::idnum) {
                mtype(nmatch) = local;
                local_match(nmatch) = bin;
                snd_tags(nmatch) = snd_tag;
                rcv_tags(nmatch) = rcv_tag;
                rcvbuf(nmatch) = xmalloc(buffsize);
                Array<FLT,1> temp(static_cast<FLT *>(rcvbuf(nmatch)), buffsize/sizeof(FLT), neverDeleteData);
                frcvbufarray(nmatch).reference(temp);
                Array<int,1> temp1(static_cast<int *>(rcvbuf(nmatch)), buffsize/sizeof(int), neverDeleteData);
                ircvbufarray(nmatch).reference(temp1);
                ++nmatch;
                return(0);
            }
            *BASE::x.gbl->log << "error: not local match" << BASE::idnum << bin->idnum << std::endl;
            return(1);
        }
        
#ifdef MPISRC
        int mpi_cnnct(int nproc, int snd_tag, int rcv_tag) {
            mtype(nmatch) = mpi;
            mpi_match(nmatch) = nproc;
            snd_tags(nmatch) = snd_tag;
            rcv_tags(nmatch) = rcv_tag;
            rcvbuf(nmatch) = xmalloc(buffsize);
            Array<FLT,1> temp(static_cast<FLT *>(rcvbuf(nmatch)), buffsize/sizeof(FLT), neverDeleteData);
            frcvbufarray(nmatch).reference(temp);
            Array<int,1> temp1(static_cast<int *>(rcvbuf(nmatch)), buffsize/sizeof(int), neverDeleteData);
            ircvbufarray(nmatch).reference(temp1);
            ++nmatch;
            return(0);
        }
#endif

        /* MECHANISM FOR SYMMETRIC SENDING */
        void comm_prepare(boundary::groups grp, int phi, boundary::comm_type type) {            
            int m;
            int nrecvs_to_post;
            int nsends_to_post;
            
#ifdef MPISRC
            int err;
#endif
            
            if (!((1<<grp)&groupmask)) return;
            
            /* SWITCHES FOR MASTER_SLAVE */
            switch(type) {
                case(boundary::master_slave): {
                    if (first) {
                        nsends_to_post = nmatch;
                        nrecvs_to_post = 0;
                    }
                    else {
                        nrecvs_to_post = 1;
                        nsends_to_post = 0;
                    }
                    break;
                }
                case(boundary::slave_master): {
                    if (!first) {
                        nsends_to_post = 1;
                        nrecvs_to_post = 0;
                    }
                    else {
                        nrecvs_to_post = nmatch;
                        nsends_to_post = 0;
                    }
                    break;
                }
                case(boundary::symmetric): {
                    nrecvs_to_post = nmatch;
                    nsends_to_post = nmatch;
                    break;
                }
            }
            
            /* MPI POST RECEIVES FIRST */
            for(m=0;m<nrecvs_to_post;++m) {
                if (phi != phase(grp)(m)) continue;
                
#ifdef MPDEBUG
                *(gbl->) << "preparing for receipt of message: " << BASE::idprefix << " with Group, Phase, Type " << grp << ',' << phi << ',' <<  type << " tag " << rcv_tags(m) << " first:" <<  is_frst()  << " with type: " ;
#endif      
                
                switch(mtype(m)) {
                    case(local):
                         /* NOTHING TO DO FOR LOCAL RECEIVES */
#ifdef MPDEBUG
                        *gbl->log << "local" << std::endl << std::flush;
#endif
                        break;  
#ifdef MPISRC
                    case(mpi):
#ifdef MPDEBUG
                        *gbl->log << "mpi to process " << mpi_match(m) << std::endl << std::flush;
#endif
                        switch(msgtype) {
                            case(boundary::flt_msg): {
#ifdef SINGLE
                                err = MPI_Irecv(&frcvbuf(m,0), buffsize/sizeof(FLT), MPI_FLOAT, 
                                    mpi_match(m), rcv_tags(m), MPI_COMM_WORLD, &mpi_rcvrqst(m));
#else
                                err = MPI_Irecv(&frcvbuf(m,0), buffsize/sizeof(FLT), MPI_DOUBLE, 
                                    mpi_match(m), rcv_tags(m), MPI_COMM_WORLD, &mpi_rcvrqst(m));
#endif    
                                break;
                            }
                            case(boundary::int_msg): {
                                err = MPI_Irecv(&ircvbuf(m,0), buffsize/sizeof(int), MPI_INT, 
                                    mpi_match(m), rcv_tags(m), MPI_COMM_WORLD, &mpi_rcvrqst(m));
                                break;
                            }
                        }
#endif        
                }
            }
            
            /* LOCAL POST SENDS FIRST */
            for(m=0;m<nsends_to_post;++m) {
                if (phi != phase(grp)(m)) continue;
                
#ifdef MPDEBUG
                *(gbl->) << "preparing for send of message: " << BASE::idprefix << " with Group, Phase, Type " << grp << ',' << phi << ',' <<  type << " tag " << snd_tags(m) << " first:" <<  is_frst()  << " with type: " ;
#endif      
                
                switch(mtype(m)) {
                    case(local):
#ifdef MPDEBUG
                        *gbl->log << "local" << std::endl << std::flush;
#endif
                        sim::blks.notify_change(snd_tags(m),true);
                        break;  
#ifdef MPISRC
                    case(mpi):
                        /* NOTHING TO DO FOR MPI SENDS */
#ifdef MPDEBUG
                        *gbl->log << "mpi to process " << mpi_match(m) << std::endl << std::flush;
#endif
                        break;
#endif        
                }
            }
        }
        
        void comm_exchange(boundary::groups grp, int phi, boundary::comm_type type) {
            int i,m;
#ifdef MPISRC
            int err;
#endif
            int nrecvs_to_post = nmatch, nsends_to_post = nmatch;
            
            if (!((1<<grp)&groupmask)) return;
            
            switch(type) {
                case(boundary::master_slave): {
                    if (!first) {
                        nrecvs_to_post = 1;
                        nsends_to_post = 0;
                    }
                    else {
                        nrecvs_to_post = 0;
                        nsends_to_post = nmatch;
                    }
                    break;
                }
                case(boundary::slave_master): {
                    if (first) {
                        nrecvs_to_post = nmatch;
                        nsends_to_post = 0;
                    }
                    else {
                        nrecvs_to_post = 0;
                        nsends_to_post = 1;
                    }
                    break;
                }
                case(boundary::symmetric): {
                    nrecvs_to_post = nmatch;
                    nsends_to_post = nmatch;
                    break;
                }
            }     
            

            /* LOCAL PASSES */
            for(m=0;m<nrecvs_to_post;++m) {
                if (phi != phase(grp)(m) || mtype(m) != local) continue;
#ifdef MPDEBUG
                *(gbl->) << "exchanging local float message: " << BASE::idprefix << " with Group, Phase, Type " << grp << ',' << phi << ',' <<  type << " tag " << rcv_tags(m) << " first:" <<  is_frst()  << " with type: " << mtype(m) << std::endl;
//                        for(i=0;i<local_match(m)->sndsize();++i) 
//                            *gbl->log << "\t" << local_match(m)->fsndbuf(i) << std::endl;
#endif      
                sim::blks.waitforslot(rcv_tags(m),true);
                switch(msgtype) {
                    case(boundary::flt_msg): {
                        for(i=0;i<local_match(m)->sndsize();++i) 
                            frcvbuf(m,i) = local_match(m)->fsndbuf(i);
                        break;
                    }
                    case(boundary::int_msg): {
                        for(i=0;i<local_match(m)->sndsize();++i) 
                            ircvbuf(m,i) = local_match(m)->isndbuf(i);
                        break;   
                    }
                }
                sim::blks.notify_change(rcv_tags(m),false);
            }
  
#ifdef MPISRC
            /* MPI PASSES */
            for(m=0;m<nsends_to_post;++m) {
                if (phi != phase(grp)(m) || mtype(m) != mpi) continue;
                
#ifdef MPDEBUG
                *(gbl->) << "exchanging mpi float message with process: " << mpi_match(m) << ' ' << BASE::idprefix << " with Group, Phase, Type " << grp << ',' << phi << ',' <<  type << " tag " << snd_tags(m) << " first:" <<  is_frst()  << " with type: " << mtype(m) << std::endl;
#endif  
                switch(msgtype) {
                    case(boundary::flt_msg): {
#ifdef SINGLE
                        err = MPI_Isend(&fsndbuf(0), msgsize, MPI_FLOAT, 
                            mpi_match(m), snd_tags(m), MPI_COMM_WORLD, &mpi_sndrqst(m));
#else
                        err = MPI_Isend(&fsndbuf(0), msgsize, MPI_DOUBLE, 
                            mpi_match(m), snd_tags(m), MPI_COMM_WORLD, &mpi_sndrqst(m));                 
#endif
                        break;
                    }
                    case(boundary::int_msg): {
                        err = MPI_Isend(&isndbuf(0), msgsize, MPI_INT, 
                            mpi_match(m), snd_tags(m), MPI_COMM_WORLD, &mpi_sndrqst(m));
                        break;
                    }
                }
            }
#endif
            return;
        }
        
        int comm_wait(boundary::groups grp, int phi, boundary::comm_type type) {
            int nrecvs_to_post;
            int nsends_to_post;
            
            if (!((1<<grp)&groupmask)) return(1);
            
            switch(type) {
                case(boundary::master_slave): {
                    if (first) {
                        nrecvs_to_post = 0;
                        nsends_to_post = nmatch;
                    }
                    else {
                        nrecvs_to_post = 1;
                        nsends_to_post = 0;
                    }
                    break;
                }
                
                case(boundary::slave_master): {
                    if (!first) {
                        nrecvs_to_post = 0;
                        nsends_to_post = 1;
                    }
                    else {
                        nrecvs_to_post = nmatch;
                        nsends_to_post = 0;
                    }
                    break;
                }
                
                case(boundary::symmetric): {
                    nrecvs_to_post = nmatch;
                    nsends_to_post = nmatch;
                    break;
                }
            }
            
            for(int m=0;m<nsends_to_post;++m) {
                if (phi != phase(grp)(m)) continue;
                
                switch(mtype(m)) {
                    case(local): {
                        sim::blks.waitforslot(snd_tags(m),false);
                        break;
                    }
#ifdef MPISRC
                    case(mpi): {
                        MPI_Status status;
#ifdef MPDEBUG
                        *gbl->log << "waiting for mpi send to complete: " << mpi_match(m) << ' ' << BASE::idprefix << " phase " << phi << " tag " << snd_tags(m) << " with type: " << mtype(m) << std::endl;
#endif
                        MPI_Wait(&mpi_sndrqst(m), &status); 
                        break;
                    }
#endif
                }
            }

         
            for(int m=0;m<nrecvs_to_post;++m) {
                if (phi != phase(grp)(m)) continue;
                
                switch(mtype(m)) {
                    case(local): {
                        break;
                    }
#ifdef MPISRC
                    case(mpi): {
                        MPI_Status status;
#ifdef MPDEBUG
                        *gbl->log << "waiting for mpi message from process: " << mpi_match(m) << ' ' << BASE::idprefix << " phase " << phi << " tag " << rcv_tags(m) << " with type: " << mtype(m)  << std::endl;
#endif
                        MPI_Wait(&mpi_rcvrqst(m), &status); 
                        break;
                    }
#endif
                }
                
#ifdef MPDEBUG
                if (msgtype == boundary::flt_msg) {
                    *(gbl->) << "received float message: " << BASE::idprefix << " with Group, Phase, Type " << grp << ',' << phi << ',' <<  type << " tag " << rcv_tags(m) << " with type: " << mtype(m) << std::endl;
                }
                else {
                    *(gbl->) << "received int message: " << BASE::idprefix << " with Group, Phase, Type " << grp << ',' << phi << ',' <<  type << " tag " << rcv_tags(m) << " with type: " << mtype(m) << std::endl;
                }
#endif  
            }
                        
            /* ONE MEANS FINISHED 0 MEANS MORE TO DO */
            return((phi-maxphase(grp) >= 0 ? 1 : 0));
        }
};


/* SOME GENERIC SHAPES */
template<int ND> class geometry {
    protected:
        virtual FLT hgt(TinyVector<FLT,ND> pt, FLT time = 0.0) {return(0.0);}
        virtual FLT dhgt(int dir, TinyVector<FLT,ND> pt, FLT time = 0.0) {return(1.0);}

    public:        
        geometry* create() const {return(new geometry);}
        virtual void mvpttobdry(TinyVector<FLT,2> &pt) {
            int iter,n;
            FLT mag, delt_dist;
                
            /* FOR AN ANALYTIC SURFACE */
            iter = 0;
            do {
                mag = 0.0;
                for(n=0;n<ND;++n)
                    mag += pow(dhgt(n,pt.data()),2);
                mag = sqrt(mag);
                delt_dist = -hgt(pt.data())/mag;
                for(n=0;n<ND;++n)
                    pt(n) += delt_dist*dhgt(n,pt.data())/mag;
                if (++iter > 100) {
                    std::cout << "curved iterations exceeded curved boundary " << pt(0) << ' ' << pt(1) << '\n';  // TEMPORARY NEED TO FIX
                    exit(1);
                }
            } while (fabs(delt_dist) > 10.*EPSILON);
            
            return;
        }
        virtual void input(input_map& inmap, std::string idprefix, std::ostream& log) {}
        virtual void output(std::ostream& fout, std::string idprefix) {}
        virtual ~geometry() {}
};

template<int ND> class vgeometry_interface {
    public:
        virtual void mvpttobdry(TinyVector<FLT,ND> &pt) {}
        virtual ~vgeometry_interface() {}
};

template<int ND> class egeometry_interface {
    public:
        virtual void mvpttobdry(int nel, FLT psi, TinyVector<FLT,ND> &pt) {}
        virtual ~egeometry_interface() {}
};

template<int ND> class fgeometry_interface {
    public:
        virtual void mvpttobdry(int nel, FLT psi, FLT eta, TinyVector<FLT,ND> &pt) {}
        virtual ~fgeometry_interface() {}
};


template<int ND> class symbolic_shape : public geometry<ND> {
    protected:
        symbolic_function<ND> h, dhdx0, dhdx1;
        FLT hgt(TinyVector<FLT,ND> pt, FLT time = 0.0) {
            return(h.Eval(pt,time));
        } 
        FLT dhgt(int dir, TinyVector<FLT,ND> pt, FLT time = 0.0) {
            if (dir) return(dhdx1.Eval(pt,time));
            return(dhdx0.Eval(pt,time));
        }
    public:
        void input(input_map& inmap, std::string idprefix, std::ostream& log) {    
            geometry<ND>::input(inmap,idprefix,log);
            if (inmap.find(idprefix +"_h_expression") != inmap.end()) {
                h.init(inmap,idprefix+"_h");
            }
            else {
                log << "couldn't find shape function for " << idprefix << std::endl;
                exit(1);
            }
            
            if (inmap.find(idprefix +"_dhdx0_expression") != inmap.end()) {
                dhdx0.copy_consts(h);
                dhdx0.init(inmap,idprefix+"_dhdx0");
            }
            else {
                log << "couldn't find shape function for " << idprefix << std::endl;
                exit(1);
            }
            
            if (inmap.find(idprefix +"_dhdx1_expression") != inmap.end()) {
                dhdx1.copy_consts(h);
                dhdx1.init(inmap,idprefix+"_dhdx1");
            }
            else {
                log << "couldn't find shape function for " << idprefix << std::endl;
                exit(1);
            }
        }
};

class circle : public geometry<2> {
    public:
        FLT center[2];
        FLT radius;
        FLT hgt(TinyVector<FLT,2> pt,FLT time = 0.0) {
            return(radius*radius -pow(pt[0]-center[0],2) -pow(pt[1]-center[1],2));
        }
        FLT dhgt(int dir, TinyVector<FLT,2> pt,FLT time = 0.0) {
            return(-2.*(pt[dir]-center[dir]));
        }
        
        circle() : geometry<2>(), radius(0.5) {center[0] = 0.0; center[1] = 0.0;}
        circle(const circle &inbdry) : geometry<2>(inbdry), radius(inbdry.radius) {center[0] = inbdry.center[0]; center[1] = inbdry.center[1];}
        circle* create() const {return(new circle(*this));}

        void output(std::ostream& fout,std::string idprefix) {
            geometry<2>::output(fout,idprefix);
            fout << idprefix << "_center: " << center[0] << '\t' << center[1] << std::endl;
            fout << idprefix << "_radius: " << radius << std::endl;
        }
      
         void input(input_map& inmap,std::string idprefix, std::ostream& log) {
            geometry<2>::input(inmap,idprefix,log);
            
            FLT dflt[2] = {0.0, 0.0};
            inmap.getwdefault(idprefix+"_center",center,2,dflt);
            inmap.getwdefault(idprefix+"_radius",radius,0.5);
        }
};  

class ellipse : public geometry<2> {
    public:
        TinyVector<FLT,2> axes;
        FLT hgt(TinyVector<FLT,2> pt, FLT time = 0.0) {
            return(1 -pow(pt[0]/axes(0),2) -pow(pt[1]/axes(1),2));
        }
        FLT dhgt(int dir, TinyVector<FLT,2> pt, FLT time = 0.0) {
              return(-2.*pt[dir]/pow(axes(dir),2));            
        }
        
            public:
        ellipse() {}
        ellipse(const ellipse &inbdry) : geometry<2>(inbdry), axes(inbdry.axes) {}
        ellipse* create() const {return(new ellipse(*this));}

        void output(std::ostream& fout,std::string idprefix) {
            geometry<2>::output(fout,idprefix);
            fout << idprefix << "_a" << axes(0) << std::endl;
            fout << idprefix << "_b" << axes(1) << std::endl;
        }
        void input(input_map& inmap, std::string idprefix,std::ostream& log) {
            geometry<2>::input(inmap,idprefix,log);
            inmap.getwdefault(idprefix+"_a",axes(0),1.0);
            inmap.getwdefault(idprefix+"_b",axes(1),1.0);
        }
};


class naca : public geometry<2> {
    public:
        FLT sign;
        TinyVector<FLT,5> coeff;
        FLT scale;
        FLT theta;
        TinyVector<FLT,2> pos;
        
        FLT hgt(TinyVector<FLT,2> x, FLT time = 0.0) {
            TinyVector<FLT,2> pt;
            for(int n=0;n<2;++n)
                pt[n] = x[n] -pos(n);
            
            FLT temp = pt[0]*cos(theta) -pt[1]*sin(theta);
            pt[1] = pt[0]*sin(theta) +pt[1]*cos(theta);
            pt[0] = temp;
            pt *= scale;
            FLT poly = coeff[1]*pt[0] +coeff[2]*pow(pt[0],2) +coeff[3]*pow(pt[0],3) +coeff[4]*pow(pt[0],4) - sign*pt[1];            
            return(coeff[0]*pt[0] -poly*poly/coeff[0]);
        }
        FLT dhgt(int dir, TinyVector<FLT,2> x, FLT time = 0.0) {
            TinyVector<FLT,2> pt;
            for(int n=0;n<2;++n)
                pt[n] = x[n] -pos(n);
            
            FLT temp = pt[0]*cos(theta) -pt[1]*sin(theta);
            pt[1] = pt[0]*sin(theta) +pt[1]*cos(theta);
            pt[0] = temp;
            pt *= scale;
            
            TinyVector<FLT,2> ddx; 
            FLT poly = coeff[1]*pt[0] +coeff[2]*pow(pt[0],2) +coeff[3]*pow(pt[0],3) +coeff[4]*pow(pt[0],4) - sign*pt[1];            
            FLT dpolydx = coeff[1] +2*coeff[2]*pt[0] +3*coeff[3]*pow(pt[0],2) +4*coeff[4]*pow(pt[0],3);
            FLT dpolydy = -sign;
            ddx(0) = coeff[0] -2*poly*dpolydx/coeff[0];
            ddx(1) = -2*poly*dpolydy/coeff[0];
            ddx *= scale;
            
            if (dir == 0) return(ddx(0)*cos(theta) +ddx(1)*sin(theta));
            return(ddx(0)*(-sin(theta)) +ddx(1)*cos(theta));
        }
        
        naca() : geometry<2>(), sign(1.0), scale(1.0), theta(0.0) {
            /* NACA 0012 is the default */
            sign = 1;
            coeff[0] = 1.4845; coeff[1] = -0.63; coeff[2] = -1.758; coeff[3] = 1.4215; coeff[4] = -0.5180;
            coeff *= 0.12;
            pos = 0.0;
        }
        naca(const naca &inbdry) : geometry<2>(inbdry), sign(inbdry.sign), scale(inbdry.scale), theta(inbdry.theta) {
            for(int i=0;i<5;++i) 
                coeff[i] = inbdry.coeff[i];
                
            pos = inbdry.pos;
        }
        naca* create() const {return(new naca(*this));}

        void output(std::ostream& fout,std::string idprefix) {
            geometry<2>::output(fout,idprefix);
            fout << idprefix << "_sign: " << sign << std::endl;
            fout << idprefix << "_coeff: ";
            for(int i=0;i<5;++i) 
                fout << coeff[i] << " ";
            fout << std::endl;
            fout << idprefix << "_scale: " << scale << std::endl;
            fout << idprefix << "_theta: " << theta << std::endl;
            fout << idprefix << "_center: " << pos << std::endl;
        }
      
         void input(input_map& inmap,std::string idprefix,std::ostream& log) {
            geometry<2>::input(inmap,idprefix,log);

            FLT thickness;
            inmap.getwdefault(idprefix+"_sign",sign,1.0);
            inmap.getwdefault(idprefix+"_thickness",thickness,0.12);
            inmap.getwdefault(idprefix+"_scale",scale,1.0);
            inmap.getwdefault(idprefix+"_theta",theta,0.0);
            scale = 1./scale;
            theta = theta*M_PI/180.0;
            
            std::string linebuf;
            istringstream datastream;
            inmap.getwdefault(idprefix+"_coeff",linebuf,std::string("1.4845 -0.63 -1.758 1.4215 -0.5180"));
            datastream.str(linebuf);
            for(int i=0;i<5;++i)
                datastream >> coeff[i];
            datastream.clear();
            coeff *= thickness;
            
            inmap.getwdefault(idprefix+"_center",linebuf,std::string("0.0 0.0"));
            datastream.str(linebuf);
            for(int i=0;i<2;++i)
                datastream >> pos(i);
            datastream.clear();            
        }
};         






#endif


