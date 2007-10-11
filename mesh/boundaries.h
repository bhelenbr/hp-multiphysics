#include "tri_mesh.h"

#include <iostream>
#include <input_map.h>
#include <string>
#include <sstream>
#include <fstream>
#include <mathclass.h>

//#define MPDEBUG

/** \brief Template class to make a communciation boundary 
 *
 * \ingroup boundary
 * Contains variables and routines for communciation boundaries.
 * Can be applied to vrtx_bdry or side_bdry to make specific type
 */
template<class BASE> class comm_bdry : public BASE {
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
        comm_bdry(int inid, tri_mesh &xin) : BASE(inid,xin), first(1), maxgroup(1), groupmask(0x3), buffsize(0), nmatch(0) {
            maxphase.resize(maxgroup+1);
            phase.resize(maxgroup+1);
            for(int m=0;m<maxmatch;++m) phase(0)(m) = 0;
        }
        comm_bdry(const comm_bdry<BASE> &inbdry, tri_mesh&xin) : BASE(inbdry,xin), first(inbdry.first), maxgroup(inbdry.maxgroup), groupmask(inbdry.groupmask), buffsize(0), nmatch(0) {    
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
              
        comm_bdry<BASE>* create(tri_mesh &xin) const {return(new comm_bdry<BASE>(*this,xin));}
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
            resize_buffers(nels*3);    // should be tri_mesh::ND;
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

/** \brief Specialization for a communiation vertex 
 *
 * \ingroup boundary
 */
class vcomm : public comm_bdry<vrtx_bdry> {
    public:
        vcomm(int inid, tri_mesh& xin) : comm_bdry<vrtx_bdry>(inid,xin) {mytype="comm";}
        vcomm(const vcomm &inbdry, tri_mesh& xin) : comm_bdry<vrtx_bdry>(inbdry,xin) {}
        
        vcomm* create(tri_mesh& xin) const {return new vcomm(*this,xin);}

        /** Generic routine to load buffers from array */
        void vloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride);
        /** Generic routine to receive into array */
        void vfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride);
};

/** \brief Specialization for a communiation side 
 *
 * \ingroup boundary
 */
class scomm : public comm_bdry<side_bdry> {
    public:                
        /* CONSTRUCTOR */
        scomm(int inid, tri_mesh& xin) : comm_bdry<side_bdry>(inid,xin) {mytype="comm";}
        scomm(const scomm &inbdry, tri_mesh& xin) : comm_bdry<side_bdry>(inbdry,xin) {}
        
        scomm* create(tri_mesh& xin) const {return new scomm(*this,xin);}
        
        /* GENERIC COMMUNICATIONS */
        void vloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride);
        void vfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride);
        void sloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride);
        void sfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride);
};

/** \brief Specialization for a parition side 
 *
 * \ingroup boundary
 * This is specifically for partitions that are not smooth boundaries
 * as typically created by metis.  Because they are not smooth it 
 * changes what happens between different levels of multigrid.
 */
class spartition : public scomm {
    public:
        /* CONSTRUCTOR */
        spartition(int inid, tri_mesh& xin) : scomm(inid,xin) {groupmask = 1;mytype="partition";}
        spartition(const spartition &inbdry, tri_mesh& xin) : scomm(inbdry,xin) {}

        spartition* create(tri_mesh& xin) const {return new spartition(*this,xin);}
        void mgconnect(Array<tri_mesh::transfer,1> &cnnct, tri_mesh& tgt, int bnum);
};


/* CAN DO PERIODIC IN X OR Z (dir = 0/1) IN 3D */
template<class BASE> class prdc_template : public BASE {
    protected:
        int dir;
    public:        
        /* CONSTRUCTOR */
        prdc_template(int idin, tri_mesh &xin) : BASE(idin,xin), dir(0) {BASE::mytype="prdc";}
        prdc_template(const prdc_template<BASE> &inbdry, tri_mesh &xin) : BASE(inbdry,xin), dir(inbdry.dir) {}
        
        prdc_template<BASE>* create(tri_mesh& xin) const {return(new prdc_template<BASE>(*this,xin));}

        int& setdir() {return(dir);}
        void output(std::ostream& fout) {
            BASE::output(fout);
            fout << BASE::idprefix << "_dir" << ": " << dir << std::endl;  
        }
        void input(input_map& inmap) {
            std::string keyword;
            std::map<std::string,std::string>::const_iterator mi;
            std::istringstream data;
            
            BASE::input(inmap);

            keyword = BASE::idprefix + "_dir";
            mi = inmap.find(keyword);
            if (mi != inmap.end()) {
                data.str(mi->second);
                data >> dir;
                data.clear();
            }
        }
        
        /* SEND/RCV VRTX POSITION */
        void loadpositions() { BASE::vloadbuff(BASE::all,&(BASE::x.vrtx(0)(0)),1-dir,1-dir +tri_mesh::ND-2,tri_mesh::ND); }
        void rcvpositions(int phase) { BASE::vfinalrcv(BASE::all_phased,phase,BASE::master_slave,boundary::replace,&(BASE::x.vrtx(0)(0)),1-dir,1-dir +tri_mesh::ND-2,tri_mesh::ND); }
};

class curved_analytic_interface {
    protected:
        virtual FLT hgt(TinyVector<FLT,tri_mesh::ND> pt, FLT time = 0.0) {return(0.0);}
        virtual FLT dhgt(int dir, TinyVector<FLT,tri_mesh::ND> pt, FLT time = 0.0) {return(1.0);}

    public:        
        curved_analytic_interface* create() const {return(new curved_analytic_interface);}
        void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt);
        virtual void input(input_map& inmap, std::string idprefix, std::ostream& log) {}
        virtual void output(std::ostream& fout, std::string idprefix) {}
        virtual ~curved_analytic_interface() {}
};


class symbolic_shape : public curved_analytic_interface {
    protected:
        symbolic_function<2> h, dhdx0, dhdx1;
        FLT hgt(TinyVector<FLT,tri_mesh::ND> pt, FLT time = 0.0) {
            return(h.Eval(pt,time));
        } 
        FLT dhgt(int dir, TinyVector<FLT,tri_mesh::ND> pt, FLT time = 0.0) {
            if (dir) return(dhdx1.Eval(pt,time));
            return(dhdx0.Eval(pt,time));
        }
    public:
        void input(input_map& inmap, std::string idprefix, std::ostream& log) {    
            curved_analytic_interface::input(inmap,idprefix,log);
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

typedef prdc_template<vcomm> vprdc;
typedef prdc_template<scomm> sprdc;


template<class BASE,class GEOM> class analytic_geometry : public BASE {
    public:
        GEOM geometry_object;
        
    public: 
        analytic_geometry(int inid, tri_mesh &xin) : BASE(inid,xin) {BASE::mytype=BASE::mytype+"analytic";}
        analytic_geometry(const analytic_geometry<BASE,GEOM> &inbdry, tri_mesh &xin) : BASE(inbdry,xin), geometry_object(inbdry.geometry_object) {}
        analytic_geometry* create(tri_mesh& xin) const {return(new analytic_geometry<BASE,GEOM>(*this,xin));}

        void output(std::ostream& fout) {
            BASE::output(fout);
            geometry_object.output(fout,BASE::idprefix);
        }
        void input(input_map& inmap) {
            BASE::input(inmap);
            geometry_object.input(inmap,BASE::idprefix,*BASE::x.gbl->log);
        }
        
        void mvpttobdry(int nel,FLT psi, TinyVector<FLT,tri_mesh::ND> &pt) {
            /* GET LINEAR APPROXIMATION FIRST */
            BASE::mvpttobdry(nel,psi,pt);
            geometry_object.mvpttobdry(pt);
            return;
        }
};


template<class BASE> class ssolution_geometry : public sgeometry_pointer, public BASE {
  public: 
        ssolution_geometry(int inid, tri_mesh &xin) : BASE(inid,xin) {BASE::mytype=BASE::mytype+"coupled";}
        ssolution_geometry(const ssolution_geometry<BASE> &inbdry, tri_mesh &xin) : BASE(inbdry,xin) {}
        ssolution_geometry* create(tri_mesh& xin) const {return(new ssolution_geometry<BASE>(*this,xin));}

        virtual void mvpttobdry(int nel,FLT psi, TinyVector<FLT,tri_mesh::ND> &pt) {
            if (BASE::x.gbl->tstep < 0) BASE::mvpttobdry(nel,psi,pt);
            else solution_data->mvpttobdry(nel,psi,pt);
            return;
        }
};


class circle : public curved_analytic_interface {
    public:
        FLT center[tri_mesh::ND];
        FLT radius;
        FLT hgt(TinyVector<FLT,tri_mesh::ND> pt,FLT time = 0.0) {
            return(radius*radius -pow(pt[0]-center[0],2) -pow(pt[1]-center[1],2));
        }
        FLT dhgt(int dir, TinyVector<FLT,tri_mesh::ND> pt,FLT time = 0.0) {
            return(-2.*(pt[dir]-center[dir]));
        }
        
        circle() : curved_analytic_interface(), radius(0.5) {center[0] = 0.0; center[1] = 0.0;}
        circle(const circle &inbdry) : curved_analytic_interface(inbdry), radius(inbdry.radius) {center[0] = inbdry.center[0]; center[1] = inbdry.center[1];}
        circle* create() const {return(new circle(*this));}

        void output(std::ostream& fout,std::string idprefix) {
            curved_analytic_interface::output(fout,idprefix);
            fout << idprefix << "_center: " << center[0] << '\t' << center[1] << std::endl;
            fout << idprefix << "_radius: " << radius << std::endl;
        }
      
         void input(input_map& inmap,std::string idprefix, std::ostream& log) {
            curved_analytic_interface::input(inmap,idprefix,log);
            
            FLT dflt[2] = {0.0, 0.0};
            inmap.getwdefault(idprefix+"_center",center,tri_mesh::ND,dflt);
            inmap.getwdefault(idprefix+"_radius",radius,0.5);
        }
};  

class ellipse : public curved_analytic_interface {
    public:
        TinyVector<FLT,2> axes;
        FLT hgt(TinyVector<FLT,tri_mesh::ND> pt, FLT time = 0.0) {
            return(1 -pow(pt[0]/axes(0),2) -pow(pt[1]/axes(1),2));
        }
        FLT dhgt(int dir, TinyVector<FLT,tri_mesh::ND> pt, FLT time = 0.0) {
              return(-2.*pt[dir]/pow(axes(dir),2));            
        }
        
            public:
        ellipse() {}
        ellipse(const ellipse &inbdry) : curved_analytic_interface(inbdry), axes(inbdry.axes) {}
        ellipse* create() const {return(new ellipse(*this));}

        void output(std::ostream& fout,std::string idprefix) {
            curved_analytic_interface::output(fout,idprefix);
            fout << idprefix << "_a" << axes(0) << std::endl;
            fout << idprefix << "_b" << axes(1) << std::endl;
        }
        void input(input_map& inmap, std::string idprefix,std::ostream& log) {
            curved_analytic_interface::input(inmap,idprefix,log);
            inmap.getwdefault(idprefix+"_a",axes(0),1.0);
            inmap.getwdefault(idprefix+"_b",axes(1),1.0);
        }
};


class naca : public curved_analytic_interface {
    public:
        FLT sign;
        TinyVector<FLT,5> coeff;
        FLT scale;
        FLT theta;
        TinyVector<FLT,2> pos;
        
        FLT hgt(TinyVector<FLT,tri_mesh::ND> x, FLT time = 0.0) {
            TinyVector<FLT,tri_mesh::ND> pt;
            for(int n=0;n<tri_mesh::ND;++n)
                pt[n] = x[n] -pos(n);
            
            FLT temp = pt[0]*cos(theta) -pt[1]*sin(theta);
            pt[1] = pt[0]*sin(theta) +pt[1]*cos(theta);
            pt[0] = temp;
            pt *= scale;
            FLT poly = coeff[1]*pt[0] +coeff[2]*pow(pt[0],2) +coeff[3]*pow(pt[0],3) +coeff[4]*pow(pt[0],4) - sign*pt[1];            
            return(coeff[0]*pt[0] -poly*poly/coeff[0]);
        }
        FLT dhgt(int dir, TinyVector<FLT,tri_mesh::ND> x, FLT time = 0.0) {
            TinyVector<FLT,tri_mesh::ND> pt;
            for(int n=0;n<tri_mesh::ND;++n)
                pt[n] = x[n] -pos(n);
            
            FLT temp = pt[0]*cos(theta) -pt[1]*sin(theta);
            pt[1] = pt[0]*sin(theta) +pt[1]*cos(theta);
            pt[0] = temp;
            pt *= scale;
            
            TinyVector<FLT,tri_mesh::ND> ddx; 
            FLT poly = coeff[1]*pt[0] +coeff[2]*pow(pt[0],2) +coeff[3]*pow(pt[0],3) +coeff[4]*pow(pt[0],4) - sign*pt[1];            
            FLT dpolydx = coeff[1] +2*coeff[2]*pt[0] +3*coeff[3]*pow(pt[0],2) +4*coeff[4]*pow(pt[0],3);
            FLT dpolydy = -sign;
            ddx(0) = coeff[0] -2*poly*dpolydx/coeff[0];
            ddx(1) = -2*poly*dpolydy/coeff[0];
            ddx *= scale;
            
            if (dir == 0) return(ddx(0)*cos(theta) +ddx(1)*sin(theta));
            return(ddx(0)*(-sin(theta)) +ddx(1)*cos(theta));
        }
        
        naca() : curved_analytic_interface(), sign(1.0), scale(1.0), theta(0.0) {
            /* NACA 0012 is the default */
            sign = 1;
            coeff[0] = 1.4845; coeff[1] = -0.63; coeff[2] = -1.758; coeff[3] = 1.4215; coeff[4] = -0.5180;
            coeff *= 0.12;
            pos = 0.0;
        }
        naca(const naca &inbdry) : curved_analytic_interface(inbdry), sign(inbdry.sign), scale(inbdry.scale), theta(inbdry.theta) {
            for(int i=0;i<5;++i) 
                coeff[i] = inbdry.coeff[i];
                
            pos = inbdry.pos;
        }
        naca* create() const {return(new naca(*this));}

        void output(std::ostream& fout,std::string idprefix) {
            curved_analytic_interface::output(fout,idprefix);
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
            curved_analytic_interface::input(inmap,idprefix,log);

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
            for(int i=0;i<tri_mesh::ND;++i)
                datastream >> pos(i);
            datastream.clear();            
        }
};         










