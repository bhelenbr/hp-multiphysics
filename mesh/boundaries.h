#include "boundary.h"

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

#define NO_MPDEBUG

/* TEMPLATE CLASS TO MAKE A COMMUNCIATION BOUNDARY */
template<class BASE> class comm_bdry : public BASE {
   protected:
      static const int maxmatch = 8;
      bool first; // FOR MASTER-SLAVE COMMUNICATION
      int grouping;  // TO MAKE GROUPS THAT ONLY COMMUNICATE IN RESTRICTED SITUATIONS GROUP 1 ALL, GROUP 2 PARTITIONS 
      int maxphase; //  FOR SYMMETRIC MESSAGE PASSING
      int buffsize; // SIZE OF BUFFER (TIMES SIZEOF(FLT))
      void *sndbuf; // OUTGOING MESSAGE BUFFER
      int msgsize; // OUTGOING SIZE
      boundary::msg_type msgtype; // OUTGOING TYPE
      
      /* MATCHES */
      int nmatch; // NUMBER OF LOCAL MATCHING BOUNDARIES
#ifdef MPISRC
      enum matchtype {local, mpi};
#else
      enum matchtype {local};
#endif
      matchtype mtype[maxmatch]; // LOCAL OR MPI OR ?
      void *local_match[maxmatch]; // POINTERS TO REMOTE BUFFERS
      int tags[maxmatch]; //IDENTIFIES EACH CONNECTION UNIQUELY
      int phase[maxmatch];  // TO SET-UP STAGGERED SEQUENCE OF SYMMETRIC PASSES
      void *rcvbuf[maxmatch]; //LOCAL BUFFERS TO STORE INCOMING MESSAGES
#ifdef MPISRC
      int mpi_match[maxmatch]; //PROCESSOR NUMBERS FOR MPI
      MPI_Request mpi_rcvrqst[maxmatch]; // IDENTIFIER RETURNED FROM MPI TO MONITOR SUCCESS OF RECV
      MPI_Request mpi_sndrqst[maxmatch]; // IDENTIFIER RETURNED FROM MPI TO MONITOR SUCCESS OF SEND
#endif
            
   public:
      comm_bdry(int inid, mesh &xin) : BASE(inid,xin), first(1), grouping(0), maxphase(0), buffsize(0), nmatch(0) {for(int m=0;m<maxmatch;++m) phase[m] = 0;}
      comm_bdry(const comm_bdry<BASE> &inbdry, mesh&xin) : BASE(inbdry,xin), first(1), grouping(inbdry.grouping), maxphase(inbdry.maxphase), buffsize(0), nmatch(0) {for(int m=0;m<maxmatch;++m) phase[m] = inbdry.phase[m];} 
     
      comm_bdry<BASE>* create(mesh &xin) const {return(new comm_bdry<BASE>(*this,xin));}

      void output(std::ostream& fout) {
         BASE::output(fout);
         fout << BASE::idprefix << ".group" << ": " << grouping << std::endl;  
         
         fout << BASE::idprefix << ".phase" << ": ";
         for (int m=0;m<nmatch;++m)
            fout << phase[m] << " ";
         fout << std::endl;
      }
      void input(std::map <std::string,std::string>& inmap) {
         int m;
         std::string keyword;
         std::map<std::string,std::string>::const_iterator mi;
         std::istringstream data;
         
         BASE::input(inmap);
         grouping = 0;
         keyword = BASE::idprefix + ".group";
         mi = inmap.find(keyword);
         if (mi != inmap.end()) {
            data.str(mi->second);
            data >> grouping;
            data.clear();
         }
         
         for (m=0;m<maxmatch;++m)
            phase[m] = 0;
            
         keyword = BASE::idprefix + ".phase";
         mi = inmap.find(keyword);
         if (mi != inmap.end()) {
            data.str(mi->second);
            m = 0;
            maxphase = 0;
            while(data >> phase[m]) {
               maxphase = MAX(maxphase,phase[m]);
               ++m;
            }
            data.clear();
         }
      }
      bool is_comm() {return(true);}
      bool& is_frst() {return(first);}
      int& group() {return(grouping);}
      int& sndsize() {return(msgsize);}
      boundary::msg_type& sndtype() {return(msgtype);}
      void *psndbuf() {return(sndbuf);}
      int& isndbuf(int n) {return(static_cast<int *>(sndbuf)[n]);}
      FLT& fsndbuf(int n) {return(static_cast<FLT *>(sndbuf)[n]);}
      void *prcvbuf(int match) {return(rcvbuf);}
      int& ircvbuf(int match,int n) {return(static_cast<int *>(rcvbuf[match])[n]);}
      FLT& frcvbuf(int match,int n) {return(static_cast<FLT *>(rcvbuf[match])[n]);}
      
      void alloc(int size) {
         BASE::alloc(size);
         buffsize = size*3;   // should be mesh::ND;
         sndbuf = xmalloc(buffsize*sizeof(FLT)); 
      }
      void resize_buffers(int size) {
         if (sndbuf) free(sndbuf);
         buffsize = size;
         sndbuf = xmalloc(buffsize*sizeof(FLT)); 
      }
      
      void setphase(int newphase, int msg_tag) {
         int i;
         for(i=0;i<nmatch;++i) {
            if (tags[i] == msg_tag) {
               phase[i] = newphase;
               maxphase = MAX(maxphase,newphase);
            }
         }
         return;
      }

      int local_cnnct(boundary *bin, int msg_tag) {
         if (bin->idnum == BASE::idnum) {
            mtype[nmatch] = local;
            local_match[nmatch] = bin->psndbuf();
            tags[nmatch] = msg_tag;
            rcvbuf[nmatch] = xmalloc(buffsize*sizeof(FLT));
            
            ++nmatch;
            return(0);
         }
         *BASE::x.log << "error: not local match" << BASE::idnum << bin->idnum << std::endl;
         return(1);
      }
      
#ifdef MPISRC
      int mpi_cnnct(int nproc, int msg_tag) {
         mtype[nmatch] = mpi;
         mpi_match[nmatch] = nproc;
         tags[nmatch] = msg_tag;
         rcvbuf[nmatch] = xmalloc(buffsize*sizeof(FLT));
         ++nmatch;
         return(0);
      }
#endif

      /* MECHANISM FOR SYMMETRIC SENDING */
      void comm_prepare(int phi) {         
         int m;
         switch(msgtype) {
            case(boundary::flt_msg):
               /* MPI POST RECEIVES FIRST */
               for(m=0;m<nmatch;++m) {
                  if (phi != phase[m]) continue;
                  
                  switch(mtype[m]) {
                     case(local):
                        /* NOTHING TO DO FOR LOCAL PASSES (BUFFER ALREADY LOADED) */
                        break;  
#ifdef MPISRC
                     case(mpi):
#ifdef SINGLE
                        MPI_Irecv(&frcvbuf(m,0), msgsize, MPI_FLOAT, 
                           mpi_match[m], tags[m], MPI_COMM_WORLD, &mpi_rcvrqst[m]);
#else
                        MPI_Irecv(&frcvbuf(m,0), msgsize, MPI_DOUBLE, 
                           mpi_match[m], tags[m], MPI_COMM_WORLD, &mpi_rcvrqst[m]); 
                        break;
#endif   
#endif      
                  }
               }
               break;
               
            case(boundary::int_msg):
               /* MPI POST RECEIVES FIRST */
               for(m=0;m<nmatch;++m) {
                  if (phi != phase[m]) continue;
                  
                  switch(mtype[m]) {
                     case(local):
                        /* NOTHING TO DO FOR LOCAL PASSES (BUFFER ALREADY LOADED) */
                        break;  
#ifdef MPISRC
                     case(mpi):
                        MPI_Irecv(&ircvbuf(m,0), msgsize, MPI_INT, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);
                        break;
#endif         
                  }

               }
               break;              
         }
      }
      
      void comm_transmit(int phi) {
         int i,m;
         
         switch(msgtype) {
            case(boundary::flt_msg):
               /* LOCAL PASSES */
               for(m=0;m<nmatch;++m) {
                  if (phi != phase[m]) continue;
#ifdef MPDEBUG
                  *(BASE::x.log) << "sending message: " << BASE::idnum << " " << tags[m] << " first:" <<  is_frst()  << " with type: " << mtype[m] << std::endl;
                  for(i=0;i<msgsize;++i) 
                     *(BASE::x.log) << "\t" << fsndbuf(i) << std::endl;
#endif     
                  switch(mtype[m]) {
                  
                     case(local):
                        for(i=0;i<msgsize;++i) 
                           frcvbuf(m,i) = static_cast<FLT *>(local_match[m])[i];
                        break;
#ifdef MPISRC
                     case(mpi):
#ifdef SINGLE
                        MPI_Isend(&fsndbuf(0), msgsize, MPI_FLOAT, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
#else
                        MPI_Isend(&fsndbuf(0), msgsize, MPI_DOUBLE, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);             
#endif
                        break;
#endif
                  }
               }
               break;
               
            case(boundary::int_msg):
               /* LOCAL PASSES */
               for(m=0;m<nmatch;++m) {
                  if (phi != phase[m]) continue;
                  
                  switch(mtype[m]) {
                     case(local):
                        for(i=0;i<msgsize;++i)
                           ircvbuf(m,i) = static_cast<int *>(local_match[m])[i];
                        break;
#ifdef MPISRC
                     case(mpi):
                     /* MPI PASSES */                  
                        MPI_Isend(&isndbuf(0), msgsize, MPI_INT, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
                        break;
#endif
                  }
               }
               break;              
         }
         
         return;
      }
      
      int comm_wait(int phi) {

         for(int m=0;m<nmatch;++m) {
            if (phi != phase[m]) continue;
            
            switch(mtype[m]) {
               case(local):
                  break;
#ifdef MPISRC
               case(mpi):
                  MPI_Status status;
                  MPI_Wait(&mpi_rcvrqst[m], &status); 
                  break;
#endif
            }
            
#ifdef MPDEBUG
            *(BASE::x.log) << "received message: " << BASE::idnum << " " << tags[m] << " with type: " << mtype[m] << std::endl;
            for(int i=0;i<msgsize;++i) 
               *(BASE::x.log) << "\t" << frcvbuf(m,i) << std::endl;
#endif  
         }
         /* ONE MEANS FINISHED 0 MEANS MORE TO DO */
         return((phi-maxphase >= 0 ? 1 : 0));
      }
      
      /* MECHANISM FOR MASTER SLAVE COMMUNICATIONS */
      void master_slave_prepare() {
         if (first) return;

         switch(msgtype) {
            case(boundary::flt_msg):
                switch(mtype[0]) {
                  case(local):
                     /* NOTHING TO DO FOR LOCAL PASSES (BUFFER ALREADY LOADED) */
                     return;
#ifdef MPISRC
                  case(mpi):
                     /* MPI POST RECEIVES FIRST */
#ifdef SINGLE
                     MPI_Irecv(&frcvbuf(0,0), msgsize, MPI_FLOAT, 
                        mpi_match[0], tags[0],MPI_COMM_WORLD, &mpi_rcvrqst[0]);
#else
                     MPI_Irecv(&frcvbuf(0,0), msgsize, MPI_DOUBLE, 
                        mpi_match[0], tags[0],MPI_COMM_WORLD, &mpi_rcvrqst[0]);             
#endif
                     break;
#endif
               }

                     
            case(boundary::int_msg):
               switch(mtype[0]) {
                  case(local):
                     /* NOTHING TO DO FOR LOCAL PASSES (BUFFER ALREADY LOADED) */
                     return;
#ifdef MPISRC
                  case(mpi):
                     /* MPI POST RECEIVES FIRST */
                     MPI_Irecv(&ircvbuf(0,0), msgsize, MPI_INT, 
                        mpi_match[0], tags[0],MPI_COMM_WORLD, &mpi_rcvrqst[0]);
                     break;    
#endif
               }
         }
      }
      
      void master_slave_transmit() {
         if (!first) return;
         
         int i,m;
         switch(msgtype) {
            case(boundary::flt_msg):
               /* LOCAL PASSES */
               for(m=0;m<nmatch;++m) {   
                  switch(mtype[m]) {
                     case(local):
                        for(i=0;i<msgsize;++i)
                           frcvbuf(m,i) = static_cast<FLT *>(local_match[m])[i];
                        break;
#ifdef MPISRC
                     case(mpi):
               /* MPI PASSES */
#ifdef SINGLE
                        MPI_Isend(&fsndbuf(0), msgsize, MPI_FLOAT, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
#else
                        MPI_Isend(&fsndbuf(0), msgsize, MPI_DOUBLE, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
#endif
                        break;
#endif         
                  }
               }
               break;
               
            case(boundary::int_msg):
               /* LOCAL PASSES */
               for(m=0;m<nmatch;++m) {   
                  switch(mtype[m]) {
                     case(local):
                        for(i=0;i<msgsize;++i)
                           ircvbuf(m,i) = static_cast<int *>(local_match[m])[i];
                        break;
#ifdef MPISRC
                     case(mpi):
                        /* MPI PASSES */
                        MPI_Isend(&isndbuf(0), msgsize, MPI_INT, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
                        break;
#endif
                  }
               }
               break;              
         }
         return;
      }

      void master_slave_wait() {
         if (first) return;

         switch(mtype[0]) {
            case(local):
               return;
#ifdef MPISRC
            case(mpi):
               MPI_Status status;
               MPI_Wait(&mpi_rcvrqst[0], &status); 
               return;
#endif
         }
         return;
      }

      /* MECHANISM FOR MASTER SLAVE COMMUNICATIONS */
      void slave_master_prepare() {
         
         if (!first) return;
 
         int m;
         switch(msgtype) {
            case(boundary::flt_msg):
               /* MPI POST RECEIVES FIRST */
               for(m=0;m<nmatch;++m) {
                  switch(mtype[m]) {
                     case(local):
                        /* NOTHING TO DO */
                        break;
#ifdef MPISRC
                     case(mpi):
#ifdef SINGLE
                        MPI_Irecv(&frcvbuf(m,0), msgsize, MPI_FLOAT, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);
#else
                        MPI_Irecv(&frcvbuf(m,0), msgsize, MPI_DOUBLE, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);  
#endif
                        break;
#endif
                  }
               }
               break;
               
            case(boundary::int_msg):
               /* MPI POST RECEIVES FIRST */
               for(m=0;m<nmatch;++m) {
                  switch(mtype[m]) {
                     case(local):
                        /* NOTHING TO DO */
                        break;
#ifdef MPISRC
                     case(mpi):
                        MPI_Irecv(&frcvbuf(m,0), msgsize, MPI_INT, 
                           mpi_match[m], tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);
                        break;
#endif
                  }
               }
               break;              
         }
      }
            
      void slave_master_transmit() {
         if (first) return;
         
         int i;
         
         switch(msgtype) {
            case(boundary::flt_msg):
               switch(mtype[0]) {
                  case(local):
                     /* LOCAL PASSES */
                     for(i=0;i<msgsize;++i)
                        frcvbuf(0,i) = static_cast<FLT *>(local_match[0])[i];
                     break;
#ifdef MPISRC
                  case(mpi):
               /* MPI PASSES */
#ifdef SINGLE
                     MPI_Isend(&fsndbuf(0), msgsize, MPI_FLOAT, 
                        mpi_match[0], tags[0],MPI_COMM_WORLD, &mpi_sndrqst[0]);
#else
                     MPI_Isend(&fsndbuf(0), msgsize, MPI_DOUBLE, 
                        mpi_match[0], tags[0],MPI_COMM_WORLD, &mpi_sndrqst[0]);             
#endif
                     break;
#endif
               }
               break;
               
            case(boundary::int_msg):
               switch(mtype[0]) {
                  case(local):
                     /* LOCAL PASSES */
                     for(i=0;i<msgsize;++i)
                        ircvbuf(0,0) = static_cast<int *>(local_match[0])[i];
                     break;
#ifdef MPISRC
                  case(mpi):
                     /* MPI PASSES */
                     MPI_Isend(&isndbuf(0), msgsize, MPI_INT, mpi_match[0], tags[0],MPI_COMM_WORLD, &mpi_sndrqst[0]);
                     break;
#endif
               }
               break;              
         }
         return;
      }

      void slave_master_wait() {
         if (!first) return;

         for(int m=0;m<nmatch;++m) {
            switch(mtype[m]) {
               case(local):
                  break;
#ifdef MPISRC
               case(mpi):
                  MPI_Status status;
                  MPI_Wait(&mpi_rcvrqst[m], &status); 
                  break;
#endif
            }
         }
         return;
      }
};

class vcomm : public comm_bdry<vrtx_bdry> {
   public:
      vcomm(int inid, mesh& xin) : comm_bdry<vrtx_bdry>(inid,xin) {mytype="comm";}
      vcomm(const vcomm &inbdry, mesh& xin) : comm_bdry<vrtx_bdry>(inbdry,xin) {}
      
      vcomm* create(mesh& xin) const {return new vcomm(*this,xin);}

      /* GENERIC VERTEX COMMUNICATIONS */
      void loadbuff(FLT *base,int bgn,int end, int stride);
      void finalrcv(int phase, FLT *base,int bgn,int end, int stride);
};

class scomm : public comm_bdry<class side_bdry> {
   public:            
      /* CONSTRUCTOR */
      scomm(int inid, mesh& xin) : comm_bdry<side_bdry>(inid,xin) {mytype="comm";}
      scomm(const scomm &inbdry, mesh& xin) : comm_bdry<side_bdry>(inbdry,xin) {}
      
      scomm* create(mesh& xin) const {return new scomm(*this,xin);}
      
      /* GENERIC COMMUNICATIONS */
      virtual void loadbuff(FLT *base,int bgn,int end, int stride);
      virtual void finalrcv(int phase, FLT *base,int bgn,int end, int stride);
};

class spartition : public scomm {
   public:
      /* CONSTRUCTOR */
      spartition(int inid, mesh& xin) : scomm(inid,xin) {grouping = 1;mytype="partition";}
      spartition(const spartition &inbdry, mesh& xin) : scomm(inbdry,xin) {}

      spartition* create(mesh& xin) const {return new spartition(*this,xin);}
      block::ctrl mgconnect(int excpt, Array<mesh::transfer,1> &cnnct, const class mesh& tgt, int bnum);
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
      void input(std::map <std::string,std::string>& inmap) {
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
      void copy(const boundary& bin) {
         BASE::copy(bin);
         dir = dynamic_cast<const prdc_template<BASE>& >(bin).dir;
      } 

      /* SEND/RCV VRTX POSITION */
      void loadpositions() { loadbuff(&(BASE::x.vrtx(0)(0)),1-dir,1-dir +mesh::ND-2,mesh::ND); }
      void rcvpositions(int phase) { finalrcv(phase,&(BASE::x.vrtx(0)(0)),1-dir,1-dir +mesh::ND-2,mesh::ND); }
};

class curved_analytic : public side_bdry {
   protected:
      virtual FLT hgt(FLT x[mesh::ND]) {return(0.0);}
      virtual FLT dhgt(int dir, FLT x[mesh::ND]) {return(1.0);}

   public:      
      /* CONSTRUCTOR */
      curved_analytic(int idin, mesh &xin) : side_bdry(idin,xin) {mytype="curved_analytic";}
      curved_analytic(const curved_analytic &inbdry, mesh &xin) : side_bdry(inbdry,xin) {}

      curved_analytic* create(mesh& xin) const {return new curved_analytic(*this,xin);}

      void mvpttobdry(int nel,FLT psi, TinyVector<FLT,mesh::ND> &pt);
};

class sinewave : public curved_analytic {
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
      sinewave(int inid, mesh &xin) : curved_analytic(inid,xin), h_or_v(0), amp(0.0), lam(1.0), phase(0.0), offset(0.0) {mytype="sinewave";}
      sinewave(const sinewave &inbdry, mesh &xin) : curved_analytic(inbdry,xin), h_or_v(inbdry.h_or_v), amp(inbdry.amp), lam(inbdry.lam), phase(inbdry.phase), offset(inbdry.offset) {}

      sinewave* create(mesh& xin) const {return(new sinewave(*this,xin));}

      void output(std::ostream& fout) {         
         curved_analytic::output(fout);
         fout << idprefix << ".h_or_v" << h_or_v << std::endl;
         fout << idprefix << ".amp: " << amp << std::endl;
         fout << idprefix << ".lam: " << lam << std::endl;
         fout << idprefix << ".phase: " << phase << std::endl;
         fout << idprefix << ".offset: " << offset << std::endl;
      }
         
      void input(std::map <std::string,std::string>& inmap) {   
         curved_analytic::input(inmap);
         
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
         
         // *(x.log) << amp << ' ' << lam << ' ' << phase << ' ' << offset/M_PI*180.0 << std::endl;
      }
};

class circle : public curved_analytic {
   public:
      FLT center[mesh::ND];
      FLT radius;
      FLT hgt(FLT pt[mesh::ND]) {
         return(radius*radius -pow(pt[0]-center[0],2) -pow(pt[1]-center[1],2));
      }
      FLT dhgt(int dir, FLT pt[mesh::ND]) {
         return(-2.*(pt[dir]-center[dir]));
      }
      
      circle(int inid, mesh &xin) : curved_analytic(inid,xin), radius(0.5) {center[0] = 0.0; center[1] = 0.0; mytype="circle";}
      circle(const circle &inbdry, mesh &xin) : curved_analytic(inbdry,xin), radius(inbdry.radius) {center[0] = inbdry.center[0]; center[1] = inbdry.center[1];}
 
      circle* create(mesh& xin) const {return(new circle(idnum,xin));}

      void output(std::ostream& fout) {
         curved_analytic::output(fout);
         fout << idprefix << ".center: " << center[0] << '\t' << center[1] << std::endl;
         fout << idprefix << ".radius: " << radius << std::endl;
      }
     
       void input(std::map <std::string,std::string>& inmap) {
         curved_analytic::input(inmap);
         
         std::istringstream data(inmap[idprefix+".center"]);
         if (!(data >> center[0] >> center[1])) {center[0] = 0.0; center[1] = 0.0;}
         data.clear();
         
         data.str(inmap[idprefix+".radius"]);
         if (!(data >> radius)) radius = 0.5;
         data.clear();
      }
};         

typedef prdc_template<vcomm> vprdc;
typedef prdc_template<scomm> sprdc;

