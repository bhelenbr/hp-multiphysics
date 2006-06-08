#include "mesh.h"

#include <iostream>
#include <input_map.h>
#include <string>
#include <sstream>
#include <fstream>

#define NO_MPDEBUG

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
      TinyVector<int,maxmatch> tags; //!< Identifies each connection uniquely
      TinyVector<void *,maxmatch> rcvbuf; //!< Raw memory to store incoming messages
      TinyVector<Array<FLT,1>,maxmatch> frcvbufarray; //!< Access to incoming message buffer for floats
      TinyVector<Array<int,1>,maxmatch> ircvbufarray; //!< Access to incoming message buffer for ints
      
#ifdef MPISRC
      TinyVector<int,maxmatch> mpi_match; //!< Processor numbers for mpi
      TinyVector<MPI_Request,maxmatch> mpi_rcvrqst; //!< Identifier returned from mpi to monitor success of recv
      TinyVector<MPI_Request,maxmatch> mpi_sndrqst; //!< Identifier returned from mpi to monitor success of send
#endif
            
   public:
      comm_bdry(int inid, mesh &xin) : BASE(inid,xin), first(1), maxgroup(1), groupmask(0x3), buffsize(0), nmatch(0) {
         maxphase.resize(maxgroup+1);
         phase.resize(maxgroup+1);
         for(int m=0;m<maxmatch;++m) phase(0)(m) = 0;
      }
      comm_bdry(const comm_bdry<BASE> &inbdry, mesh&xin) : BASE(inbdry,xin), first(inbdry.first), maxgroup(inbdry.maxgroup), groupmask(inbdry.groupmask), buffsize(0), nmatch(0) {         
         maxphase.resize(maxgroup+1);
         phase.resize(maxgroup+1);
         maxphase = inbdry.maxphase;
         for(int k=0;k<maxgroup+1;++k)
            phase(k) = inbdry.phase(k);
         
         /* COPY THESE, BUT WILL HAVE TO BE RESET TO NEW MATCHING SIDE */
         mtype = inbdry.mtype;
         local_match = local_match;
         tags = inbdry.tags;
#ifdef MPISRC
         mpi_match = inbdry.mpi_match;
#endif
         return;
      }
           
      comm_bdry<BASE>* create(mesh &xin) const {return(new comm_bdry<BASE>(*this,xin));}
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
         resize_buffers(nels*3);   // should be mesh::ND;
      }

      void input(input_map& inmap) {
         int j,k,m,maxgroup;
         std::string keyword,val;
         std::map<std::string,std::string>::const_iterator mi;
         std::istringstream data;
         std::ostringstream nstr;
         
         BASE::input(inmap);
         
         keyword = BASE::idprefix +".first";
         inmap.getwdefault(keyword,first,true);
         
         /* SET GROUP MEMBERSHIP FLAGS */
         /* ALWAYS BELONG TO FIRST 2 GROUPS */
         maxgroup = 1;
         groupmask = 0x3;
         keyword = BASE::idprefix + ".group";
         mi = inmap.find(keyword);
         if (inmap.getline(keyword,val)) {
            data.str(val);
            while(data >> m) {
               groupmask = groupmask|(1<<m);
               maxgroup = MAX(maxgroup,m);
            }
            data.clear();
         }
                  
         /* LOAD PHASES */
         maxphase.resize(maxgroup+1);
         maxphase = 0;
         phase.resize(maxgroup+1);
         phase = 0;
         /* SKIP GROUP 0 BECAUSE THAT GROUP IS NOT PHASED */
         for(k=1;k<maxgroup+1;++k) {
            if (!(groupmask&(1<<k))) continue;
            
            nstr.str("");
            nstr << BASE::idprefix << ".phase" << k << std::flush;
            mi = inmap.find(keyword);
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
         
         fout << BASE::idprefix << ".group" << ": ";
         for(int k=0;k<maxgroup+1;++k)
            if (groupmask&(1<<k)) fout << k << ' ';
         fout << std::endl;  
         
         for(int k=0;k<maxgroup+1;++k) {
            if (groupmask&(1<<k)) {
               fout << BASE::idprefix << ".phase (not set yet so this is dumb)" << k << ": ";
               for (int m=0;m<nmatch;++m)
                  fout << phase(k)(m) << " ";
               fout << std::endl;
            }
         }
      }

      int local_cnnct(boundary *bin, int msg_tag) {
         if (bin->idnum == BASE::idnum) {
            mtype(nmatch) = local;
            local_match(nmatch) = bin;
            tags(nmatch) = msg_tag;
            rcvbuf(nmatch) = xmalloc(buffsize);
            Array<FLT,1> temp(static_cast<FLT *>(rcvbuf(nmatch)), buffsize/sizeof(FLT), neverDeleteData);
            frcvbufarray(nmatch).reference(temp);
            Array<int,1> temp1(static_cast<int *>(rcvbuf(nmatch)), buffsize/sizeof(int), neverDeleteData);
            ircvbufarray(nmatch).reference(temp1);
            ++nmatch;
            return(0);
         }
         *sim::log << "error: not local match" << BASE::idnum << bin->idnum << std::endl;
         return(1);
      }
      
#ifdef MPISRC
      int mpi_cnnct(int nproc, int msg_tag) {
         mtype(nmatch) = mpi;
         mpi_match(nmatch) = nproc;
         tags(nmatch) = msg_tag;
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
         int nrecvs_to_post = nmatch;
         
         if (!((1<<grp)&groupmask)) return;
         
         /* SWITCHES FOR MASTER_SLAVE */
         switch(type) {
            case(boundary::master_slave): {
               if (first) 
                  return;
               else
                  nrecvs_to_post = 1;
               break;
            }
            case(boundary::slave_master): {
               if (!first) return;
               break;
            }
            case(boundary::symmetric): {
               break;
            }
         }
         
         switch(msgtype) {
            case(boundary::flt_msg): {
               /* MPI POST RECEIVES FIRST */
               for(m=0;m<nrecvs_to_post;++m) {
                  if (phi != phase(grp)(m)) continue;
                  
#ifdef MPDEBUG
                  *(sim::log) << "preparing for message: " << BASE::idnum << " " << tags(m) << " first:" <<  is_frst()  << " with type: " << mtype(m) << std::endl;
#endif     
                  
                  switch(mtype(m)) {
                     case(local):
                        /* NOTHING TO DO FOR LOCAL PASSES (BUFFER ALREADY LOADED) */
                        break;  
#ifdef MPISRC
                     case(mpi):
#ifdef SINGLE
                        MPI_Irecv(&frcvbuf(m,0), buffsize/sizeof(FLT), MPI_FLOAT, 
                           mpi_match(m), tags(m), MPI_COMM_WORLD, &mpi_rcvrqst(m));
#else
                        MPI_Irecv(&frcvbuf(m,0), buffsize/sizeof(FLT), MPI_DOUBLE, 
                           mpi_match(m), tags(m), MPI_COMM_WORLD, &mpi_rcvrqst(m)); 
                        break;
#endif   
#endif      
                  }
               }
               break;
            }
               
            case(boundary::int_msg): {
               /* MPI POST RECEIVES FIRST */
               for(m=0;m<nrecvs_to_post;++m) {
                  if (phi != phase(grp)(m)) continue;
                  
                  switch(mtype(m)) {
                     case(local):
                        /* NOTHING TO DO FOR LOCAL PASSES (BUFFER ALREADY LOADED) */
                        break;  
#ifdef MPISRC
                     case(mpi):
                        MPI_Irecv(&ircvbuf(m,0), buffsize/sizeof(int), MPI_INT, 
                           mpi_match(m), tags(m),MPI_COMM_WORLD, &mpi_rcvrqst(m));
                        break;
#endif         
                  }

               }
               break; 
            }
         }
      }
      
      void comm_exchange(boundary::groups grp, int phi, boundary::comm_type type) {
         int i,m;
         int nlocalmessages = nmatch, nmpimessages = nmatch;
         
         if (!((1<<grp)&groupmask)) return;
         
         switch(type) {
            case(boundary::master_slave): {
               if (!first) {
                  nlocalmessages = 1;
                  nmpimessages = 0;
               }
               else {
                  nlocalmessages = 0;
                  nmpimessages = nmatch;
               }
               break;
            }
            case(boundary::slave_master): {
               if (first) {
                  nlocalmessages = nmatch;
                  nmpimessages = 0;
               }
               else {
                  nlocalmessages = 0;
                  nmpimessages = 1;
               }
               break;
            }
            case(boundary::symmetric): {
               break;
            }
         }    
         
         switch(msgtype) {
            case(boundary::flt_msg): {
               /* LOCAL PASSES */
               for(m=0;m<nlocalmessages;++m) {
                  if (phi != phase(grp)(m) || mtype(m) != local) continue;
#ifdef MPDEBUG
                  *(sim::log) << "exchanging message: " << BASE::idnum << " " << tags(m) << " first:" <<  is_frst()  << " with type: " << mtype(m) << std::endl;
                  for(i=0;i<local_match(m)->sndsize();++i) 
                     *sim::log << "\t" << local_match(m)->fsndbuf(i) << std::endl;
#endif     
                  for(i=0;i<local_match(m)->sndsize();++i) 
                     frcvbuf(m,i) = local_match(m)->fsndbuf(i);
               }
  
#ifdef MPISRC
               /* MPI PASSES */
               for(m=0;m<nmpimessages;++m) {
                  if (phi != phase(grp)(m) || mtype(m) != mpi) continue;
                  
#ifdef MPDEBUG
                  *(sim::log) << "exchanging message: " << BASE::idnum << " " << tags(m) << " first:" <<  is_frst()  << " with type: " << mtype(m) << std::endl;
                  for(i=0;i<msgsize;++i) 
                     *(sim::log) << "\t" << fsndbuf(i) << std::endl;
#endif  

#ifdef SINGLE
                  MPI_Isend(&fsndbuf(0), msgsize, MPI_FLOAT, 
                     mpi_match(m), tags(m),MPI_COMM_WORLD, &mpi_sndrqst(m));
#else
                  MPI_Isend(&fsndbuf(0), msgsize, MPI_DOUBLE, 
                     mpi_match(m), tags(m),MPI_COMM_WORLD, &mpi_sndrqst(m));             
#endif
               }
#endif
               break;
            }
               
            case(boundary::int_msg): {
               /* LOCAL PASSES */
               for(m=0;m<nlocalmessages;++m) {
                  if (phi != phase(grp)(m) || mtype(m) != local) continue;
  
                  for(i=0;i<local_match(m)->sndsize();++i) 
                     ircvbuf(m,i) = local_match(m)->isndbuf(i);
               }
  
#ifdef MPISRC
               /* MPI PASSES */
               for(m=0;m<nmpimessages;++m) {
                  if (phi != phase(grp)(m) || mtype(m) != mpi) continue;

                  MPI_Isend(&isndbuf(0), msgsize, MPI_INT, 
                     mpi_match(m), tags(m),MPI_COMM_WORLD, &mpi_sndrqst(m));
               }
#endif
               break;     
            }
         }
         
         return;
      }
      
      int comm_wait(boundary::groups grp, int phi, boundary::comm_type type) {
         int nwait = nmatch;
         
         if (!((1<<grp)&groupmask)) return(1);
         
         switch(type) {
            case(boundary::master_slave): {
               if (first) 
                  return(1);
               else
                  nwait = 1;
               break;
            }
            
            case(boundary::slave_master): {
               if (!first) return(1);
               break;
            }
            
            case(boundary::symmetric): {
               break;
            }
         }
                  
         for(int m=0;m<nwait;++m) {
            if (phi != phase(grp)(m)) continue;
            
            switch(mtype(m)) {
               case(local): {
                  break;
               }
#ifdef MPISRC
               case(mpi): {
                  MPI_Status status;
                  MPI_Wait(&mpi_rcvrqst(m), &status); 
                  break;
               }
#endif
            }
            
#ifdef MPDEBUG
            *(sim::log) << "received message: " << BASE::idnum << " " << tags(m) << " with type: " << mtype(m) << std::endl;
            for(int i=0;i<msgsize;++i) 
               *(sim::log) << "\t" << frcvbuf(m,i) << std::endl;
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
      vcomm(int inid, mesh& xin) : comm_bdry<vrtx_bdry>(inid,xin) {mytype="comm";}
      vcomm(const vcomm &inbdry, mesh& xin) : comm_bdry<vrtx_bdry>(inbdry,xin) {}
      
      vcomm* create(mesh& xin) const {return new vcomm(*this,xin);}

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
      scomm(int inid, mesh& xin) : comm_bdry<side_bdry>(inid,xin) {mytype="comm";}
      scomm(const scomm &inbdry, mesh& xin) : comm_bdry<side_bdry>(inbdry,xin) {}
      
      scomm* create(mesh& xin) const {return new scomm(*this,xin);}
      
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
      spartition(int inid, mesh& xin) : scomm(inid,xin) {groupmask = 1;mytype="partition";}
      spartition(const spartition &inbdry, mesh& xin) : scomm(inbdry,xin) {}

      spartition* create(mesh& xin) const {return new spartition(*this,xin);}
      block::ctrl mgconnect(block::ctrl ctrl_message, Array<mesh::transfer,1> &cnnct, const class mesh& tgt, int bnum);
};


/* CAN DO PERIODIC IN X OR Z (dir = 0/1) IN 3D */
template<class BASE> class prdc_template : public BASE {
   protected:
      int dir;
   public:      
      /* CONSTRUCTOR */
      prdc_template(int idin, mesh &xin) : BASE(idin,xin), dir(0) {BASE::mytype="prdc";}
      prdc_template(const prdc_template<BASE> &inbdry, mesh &xin) : BASE(inbdry,xin), dir(inbdry.dir) {}
      
      prdc_template<BASE>* create(mesh& xin) const {return(new prdc_template<BASE>(*this,xin));}

      int& setdir() {return(dir);}
      void output(std::ostream& fout) {
         BASE::output(fout);
         fout << BASE::idprefix << ".dir" << ": " << dir << std::endl;  
      }
      void input(input_map& inmap) {
         std::string keyword;
         std::map<std::string,std::string>::const_iterator mi;
         std::istringstream data;
         
         BASE::input(inmap);

         keyword = BASE::idprefix + ".dir";
         mi = inmap.find(keyword);
         if (mi != inmap.end()) {
            data.str(mi->second);
            data >> dir;
            data.clear();
         }
      }
      
      /* SEND/RCV VRTX POSITION */
      void loadpositions() { BASE::vloadbuff(BASE::all,&(BASE::x.vrtx(0)(0)),1-dir,1-dir +mesh::ND-2,mesh::ND); }
      void rcvpositions(int phase) { BASE::vfinalrcv(BASE::all_phased,phase,BASE::master_slave,boundary::replace,&(BASE::x.vrtx(0)(0)),1-dir,1-dir +mesh::ND-2,mesh::ND); }
};

class curved_analytic_interface {
   protected:
      virtual FLT hgt(FLT x[mesh::ND]) {return(0.0);}
      virtual FLT dhgt(int dir, FLT x[mesh::ND]) {return(1.0);}

   public:      
      curved_analytic_interface* create() const {return(new curved_analytic_interface);}
      void mvpttobdry(TinyVector<FLT,mesh::ND> &pt);
      virtual void input(input_map& inmap, std::string idprefix) {}
      virtual void output(std::ostream& fout, std::string idprefix) {}
      virtual ~curved_analytic_interface() {}
};

class sinewave : public curved_analytic_interface {
   protected:
      int h_or_v;
      FLT amp, lam, phase, offset;
      FLT hgt(FLT pt[mesh::ND]) {
         return(pt[1-h_or_v] -offset -amp*sin(2.*M_PI*pt[h_or_v]/lam +phase));
      }
      FLT dhgt(int dir, FLT pt[mesh::ND]) {
         if (dir == h_or_v) 
            return(-amp*2.*M_PI/lam*cos(2.*M_PI*pt[h_or_v]/lam+phase));
         
         return(1.0);
      }
   public:
      sinewave() : h_or_v(0), amp(0.0), lam(1.0), phase(0.0), offset(0.0) {}
      sinewave(const sinewave &inbdry) : curved_analytic_interface(inbdry), h_or_v(inbdry.h_or_v), amp(inbdry.amp), lam(inbdry.lam), phase(inbdry.phase), offset(inbdry.offset) {}
      sinewave* create() const {return(new sinewave(*this));}

      void output(std::ostream& fout,std::string idprefix) {         
         curved_analytic_interface::output(fout,idprefix);
         fout << idprefix << ".h_or_v" << h_or_v << std::endl;
         fout << idprefix << ".amp: " << amp << std::endl;
         fout << idprefix << ".lam: " << lam << std::endl;
         fout << idprefix << ".phase: " << phase << std::endl;
         fout << idprefix << ".offset: " << offset << std::endl;
      }
         
      void input(input_map& inmap, std::string idprefix) {   
         curved_analytic_interface::input(inmap,idprefix);
         
         std::istringstream data(inmap[idprefix+".h_or_v"]);
         if (!(data >> h_or_v)) h_or_v = 0;
         data.clear();

         data.str(inmap[idprefix+".amp"]);
         if (!(data >> amp)) amp = 0.0;
         data.clear();
         
         data.str(inmap[idprefix +".lam"]);
         if (!(data >> lam)) lam = 1.0;
         data.clear();
         
         data.str(inmap[idprefix +".phase"]);
         if (!(data >> phase)) phase = 0.0;
         phase = phase/360.0*2.0*M_PI;
         data.clear();
         
         data.str(inmap[idprefix +".offset"]);
         if (!(data >> offset)) offset = 0.0;
         data.clear();
      }
};

class circle : public curved_analytic_interface {
   public:
      FLT center[mesh::ND];
      FLT radius;
      FLT hgt(FLT pt[mesh::ND]) {
         return(radius*radius -pow(pt[0]-center[0],2) -pow(pt[1]-center[1],2));
      }
      FLT dhgt(int dir, FLT pt[mesh::ND]) {
         return(-2.*(pt[dir]-center[dir]));
      }
      
      circle() : curved_analytic_interface(), radius(0.5) {center[0] = 0.0; center[1] = 0.0;}
      circle(const circle &inbdry) : curved_analytic_interface(inbdry), radius(inbdry.radius) {center[0] = inbdry.center[0]; center[1] = inbdry.center[1];}
      circle* create() const {return(new circle(*this));}

      void output(std::ostream& fout,std::string idprefix) {
         curved_analytic_interface::output(fout,idprefix);
         fout << idprefix << ".center: " << center[0] << '\t' << center[1] << std::endl;
         fout << idprefix << ".radius: " << radius << std::endl;
      }
     
       void input(input_map& inmap,std::string idprefix) {
         curved_analytic_interface::input(inmap,idprefix);
         
         std::istringstream data(inmap[idprefix+".center"]);
         if (!(data >> center[0] >> center[1])) {center[0] = 0.0; center[1] = 0.0;}
         data.clear();
         
         data.str(inmap[idprefix+".radius"]);
         if (!(data >> radius)) radius = 0.5;
         data.clear();
      }
};  

class naca : public curved_analytic_interface {
   public:
      FLT sign;
      FLT thickness;
      FLT coeff[5];
      
      FLT hgt(FLT pt[mesh::ND]) {
         return(thickness*(coeff[0]*sqrt(pt[0]) +coeff[1]*pt[0] +coeff[2]*pow(pt[0],2) +coeff[3]*pow(pt[0],3) +coeff[4]*pow(pt[0],4)) - sign*pt[1]);
      }
      FLT dhgt(int dir, FLT pt[mesh::ND]) {
         if (dir == 0) {
            if (pt[0] <= 0.0) return(1.0);
            return(thickness*(0.5*coeff[0]/sqrt(pt[0]) +coeff[1] +2*coeff[2]*pt[0] +3*coeff[3]*pow(pt[0],2) +4*coeff[4]*pow(pt[0],3)));
         }
         else {
            return(-sign);
         }
         return(0.0);
      }
      
      naca() : curved_analytic_interface(), sign(1.0), thickness(0.12) {
         /* NACA 0012 is the default */
         sign = 1;
         coeff[0] = 1.4845; coeff[1] = -0.63; coeff[2] = -1.758; coeff[3] = 1.4215; coeff[4] = -0.5180;
      }
      naca(const naca &inbdry) : curved_analytic_interface(inbdry), sign(inbdry.sign), thickness(inbdry.thickness) {
         for(int i=0;i<5;++i) 
            coeff[i] = inbdry.coeff[i];
      }
      naca* create() const {return(new naca(*this));}

      void output(std::ostream& fout,std::string idprefix) {
         curved_analytic_interface::output(fout,idprefix);
         fout << idprefix << ".sign: " << sign << std::endl;
         fout << idprefix << ".thickness: " << thickness << std::endl;
         fout << idprefix << ".coeff: ";
         for(int i=0;i<5;++i) 
            fout << coeff[i] << " ";
         fout << std::endl;
      }
     
       void input(input_map& inmap,std::string idprefix) {
         curved_analytic_interface::input(inmap,idprefix);
         
         std::istringstream data(inmap[idprefix+".sign"]);
         if (!(data >> sign)) sign = 1.0;
         data.clear();
         
         data.str(inmap[idprefix+".thickness"]);
         if (!(data >> thickness)) thickness = 0.12;
         data.clear();
         
         std::map<std::string,std::string>::const_iterator mi;
         std::string keyword;
         keyword = idprefix + ".coeff";
         mi = inmap.find(keyword);
         if (mi != inmap.end()) {
            data.str(mi->second);
            for(int i=0;i<5;++i)
               data >> coeff[i];
            data.clear();
         }
      }
};       

class gaussian : public curved_analytic_interface {
   public:
      FLT width,amp,power;
      TinyVector<FLT,2> intercept, normal;
      
      FLT hgt(FLT pt[mesh::ND]) {
         FLT vert = pt[0]*normal(0) +pt[1]*normal(1);
         FLT horz = (pt[0]*normal(1) -pt[1]*normal(0))/width;
         return(vert -amp*exp(-horz*horz));
      }
      FLT dhgt(int dir, FLT pt[mesh::ND]) {
         FLT horz = (pt[0]*normal(1) -pt[1]*normal(0))/width;
         if (dir == 0) {
            return(normal(0) -amp*exp(-horz*horz)*2*horz*normal(1)/width);
         }
         else {
            return(normal(1) +amp*exp(-horz*horz)*2*horz*normal(0)/width);
         }
         return(0.0);
      }
      
      gaussian() : curved_analytic_interface(), width(1.0), amp(0.1), intercept(0.0,0.0), normal(0.0,1.0) {}
      gaussian(const gaussian &inbdry) : curved_analytic_interface(inbdry), width(inbdry.width), amp(inbdry.amp), intercept(inbdry.intercept), normal(inbdry.normal) {}
      gaussian* create() const {return(new gaussian(*this));}

      void output(std::ostream& fout,std::string idprefix) {
         curved_analytic_interface::output(fout,idprefix);
         fout << idprefix << ".amp: " << amp << std::endl;
         fout << idprefix << ".width: " << width << std::endl;
         fout << idprefix << ".intercept: " << intercept << std::endl;
         fout << idprefix << ".normal: " << normal << std::endl;
      }
     
       void input(input_map& inmap,std::string idprefix) {
         curved_analytic_interface::input(inmap,idprefix);
         
         std::istringstream data(inmap[idprefix+".amp"]);
         if (!(data >> amp)) amp = 0.1;
         data.clear();
         
         data.str(inmap[idprefix+".width"]);
         if (!(data >> width)) width = 1.0;
         data.clear();
    
         data.str(inmap[idprefix+".intercept"]);
         if (!(data >> intercept)) intercept = 0.0;
         data.clear();     
         
         normal = 0.0;
         data.str(inmap[idprefix+".normal"]);
         if (!(data >> intercept)) normal(1) = 1.0;
         data.clear();     
         
         FLT length = sqrt(normal(0)*normal(0) +normal(1)*normal(1));
         normal /= length;
      }
};       



typedef prdc_template<vcomm> vprdc;
typedef prdc_template<scomm> sprdc;


template<class BASE,class GEOM> class analytic_geometry : public BASE {
   public:
      GEOM geometry_object;
      
   public: 
      analytic_geometry(int inid, mesh &xin) : BASE(inid,xin) {BASE::mytype=BASE::mytype+"analytic";}
      analytic_geometry(const analytic_geometry<BASE,GEOM> &inbdry, mesh &xin) : BASE(inbdry,xin), geometry_object(inbdry.geometry_object) {}
      analytic_geometry* create(mesh& xin) const {return(new analytic_geometry<BASE,GEOM>(*this,xin));}

      void output(std::ostream& fout) {
         BASE::output(fout);
         geometry_object.output(fout,BASE::idprefix);
      }
      void input(input_map& inmap) {
         BASE::input(inmap);
         geometry_object.input(inmap,BASE::idprefix);
      }
      
      void mvpttobdry(int nel,FLT psi, TinyVector<FLT,mesh::ND> &pt) {
         /* GET LINEAR APPROXIMATION FIRST */
         BASE::mvpttobdry(nel,psi,pt);
         geometry_object.mvpttobdry(pt);
         return;
      }
};


template<class BASE> class ssolution_geometry : public sgeometry_pointer, public BASE {
  public: 
      ssolution_geometry(int inid, mesh &xin) : BASE(inid,xin) {BASE::mytype=BASE::mytype+"coupled";}
      ssolution_geometry(const ssolution_geometry<BASE> &inbdry, mesh &xin) : BASE(inbdry,xin) {}
      ssolution_geometry* create(mesh& xin) const {return(new ssolution_geometry<BASE>(*this,xin));}

      virtual void mvpttobdry(int nel,FLT psi, TinyVector<FLT,mesh::ND> &pt) {
         if (sim::tstep <= 0) BASE::mvpttobdry(nel,psi,pt);
         else solution_data->mvpttobdry(nel,psi,pt);
         return;
      }
};









