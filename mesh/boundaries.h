#include "boundary.h"

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <fstream>

/* TEMPLATE CLASS TO MAKE A COMMUNCIATION BOUNDARY */
template<class BASE> class comm_bdry : public BASE {
   protected:
      static const int maxmatch = 8;
      int first;
      int maxphase;
      int buffsize;
      void *sndbuff;
      int *isndalias;
      FLT *fsndalias;
      int msgsize;
      boundary::msg_type msgtype;
      
      /* LOCAL MATCHES */
      int nlocal_match;
      void *local_match[maxmatch];
      int local_tags[maxmatch];
      int local_phase[maxmatch];
      void *local_rcv_buf[maxmatch];
      
      int nmpi_match;
#ifdef MPISRC
      /* MPI MATCHES */
      int mpi_match[maxmatch];
      int mpi_tags[maxmatch];
      int mpi_phase[maxmatch];
      MPI_Request mpi_rcvrqst[maxmatch];
      MPI_Request mpi_sndrqst[maxmatch];
      void *mpi_rcv_buf[maxmatch];
#endif
            
   public:
      comm_bdry(int inid, mesh &xin) : BASE(inid,xin), first(1), maxphase(0), buffsize(0), nlocal_match(0), nmpi_match(0)  {}
      comm_bdry<BASE>* create(mesh &xin) const {return(new comm_bdry<BASE>(idnum,xin));}
      bool is_comm() {return(true);}
      bool is_frst() {return(first);}
      void setfrst(bool tf) {first = tf;}; 
      int& sndsize() {return(msgsize);}
      boundary::msg_type& sndtype() {return(msgtype);}
      void *mysndbuff() {return(sndbuff);}
      int& isndbuf(int n) {return(isndalias[n]);}
      FLT& fsndbuf(int n) {return(fsndalias[n]);}
      
      void alloc(int size) {
         BASE::alloc(size);
         buffsize = size*3;   // should be mesh::DIM;
         sndbuff = xmalloc(buffsize*sizeof(FLT)); 
         isndalias = static_cast<int *>(sndbuff);
         fsndalias = static_cast<FLT *>(sndbuff);        
      }
      void resize_buffers(int size) {
         if (sndbuff) free(sndbuff);
         buffsize = size;
         sndbuff = xmalloc(buffsize*sizeof(FLT)); 
         isndalias = static_cast<int *>(sndbuff);
         fsndalias = static_cast<FLT *>(sndbuff); 
      }
      
      void setphase(int phase, int msg_tag) {
         int i;
         for(i=0;i<nlocal_match;++i) {
            if (local_tags[i] == msg_tag) {
               local_phase[i] = phase;
               maxphase = MAX(maxphase,phase);
            }
         }

#ifdef MPISRC
         for(i=0;i<nmpi_match;++i) {
            if (mpi_tags[i] == msg_tag) {
               mpi_phase[i] = phase;
               maxphase = MAX(maxphase,phase);
            }
         }
#endif
         return;
      }

      int local_cnnct(boundary *bin, int msg_tag) {
         if (bin->idnum == idnum) {
            local_match[nlocal_match] = bin->mysndbuff();
            local_tags[nlocal_match] = msg_tag;
            local_phase[nlocal_match] = 0; // DEFAULT ALL MESSAGE PASSING IN ONE PHASE
            local_rcv_buf[nlocal_match] = xmalloc(buffsize*sizeof(FLT));
            ++nlocal_match;
            return(0);
         }
         *x.log << "error: not local match" << idnum << bin->idnum << std::endl;
         return(1);
      }
      
#ifdef MPISRC
      int mpi_cnnct(int nproc, int msg_tag) {
         mpi_match[nmpi_match] = nproc;
         mpi_tags[nmpi_match] = msg_tag;
         mpi_phase[nmpi_match] = 0; // DEFAULT ALL MESSAGE PASSING IN ONE PHASE
         mpi_rcv_buf[nmpi_match] = xmalloc(buffsize*sizeof(FLT));
         ++nmpi_match;
         return(0);
      }
#endif

      /* GENERIC MECHANISM FOR SENDING */
      int rcv(int phase) {
         int i,m;
         
         switch(msgtype) {
            case(boundary::flt_msg):
               /* LOCAL PASSES */
               for(m=0;m<nlocal_match;++m) {
                  if (phase != local_phase[m]) continue;
             
                  for(i=0;i<msgsize;++i)
                     static_cast<FLT *>(local_rcv_buf[m])[i] = static_cast<FLT *>(local_match[m])[i];
               }

#ifdef MPISRC
               /* MPI PASSES */
               for(m=0;m<nmpi_match;++m) {
                  if (phase != mpi_phase[m]) continue;
#ifdef SINGLE
                  MPI_Isend(fsndalias, msgsize, MPI_FLOAT, 
                     mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
#else
                  MPI_Isend(fsndalias, msgsize, MPI_DOUBLE, 
                     mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);             
#endif
               }
#endif
               break;
               
            case(boundary::int_msg):
               /* LOCAL PASSES */
               for(m=0;m<nlocal_match;++m) {
                  if (phase != local_phase[m]) continue;
                  
                  for(i=0;i<msgsize;++i)
                     static_cast<int *>(local_rcv_buf[m])[i] = static_cast<int *>(local_match[m])[i];
               }

#ifdef MPISRC
               /* MPI PASSES */
               for(m=0;m<nmpi_match;++m) {
                  if (phase != mpi_phase[m]) continue;
                  
                  MPI_Isend(isndalias, msgsize, MPI_INT, 
                     mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
               }
#endif

               break;              
         }
         
         /* ONE MEANS FINISHED 0 MEANS MORE TO DO */
         return((phase-maxphase >= 0 ? 1 : 0));
      }

      void snd(int phase) {
#ifdef MPISRC
         int m;
#endif
         /* NOTHING TO DO FOR LOCAL PASSES (BUFFER ALREADY LOADED) */

#ifdef MPISRC
         switch(msgtype) {
            case(flt_msg):
               /* MPI POST RECEIVES FIRST */
               for(m=0;m<nmpi_match;++m) {
                  if (phase != mpi_phase[m]) continue;
#ifdef SINGLE
                  MPI_Irecv(static_cast<FLT *>(mpi_rcv_buf[m]), msgsize, MPI_FLOAT, 
                     mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);
#else
                  MPI_Irecv(static_cast<FLT *>(mpi_rcv_buf[m]), msgsize, MPI_DOUBLE, 
                     mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);             
#endif
               }
               break;
               
            case(int_msg):
               /* MPI POST RECEIVES FIRST */
               for(m=0;m<nmpi_match;++m) {
                  if (phase != mpi_phase[m]) continue;
                  
                  MPI_Irecv(static_cast<int *>(mpi_rcv_buf[m]), msgsize, MPI_INT, 
                     mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);
               }
               break;              
         }
#endif
      }

};

class vcomm : public comm_bdry<vrtx_bdry> {
   public:
      vcomm(int inid, mesh& xin) : comm_bdry<vrtx_bdry>(inid,xin) {}
      vcomm* create(mesh& xin) const {return new vcomm(idnum,xin);}

      /* GENERIC VERTEX COMMUNICATIONS */
      void loadbuff(FLT *base,int bgn,int end, int stride);
      void finalrcv(FLT *base,int bgn,int end, int stride);
};

class scomm : public comm_bdry<class side_bdry> {
   public:            
      /* CONSTRUCTOR */
      scomm(int inid, mesh& xin) : comm_bdry<side_bdry>(inid,xin) {}
      scomm* create(mesh& xin) const {return new scomm(idnum,xin);}
      
      /* GENERIC COMMUNICATIONS */
      virtual void loadbuff(FLT *base,int bgn,int end, int stride);
      virtual void finalrcv(FLT *base,int bgn,int end, int stride);
};

/* CAN DO PERIODIC IN X OR Z (dir = 0/1) IN 3D */
template<class BASE> class prdc_template : public BASE {
   protected:
      int dir;
   public:      
      /* CONSTRUCTOR */
      prdc_template(int idin, mesh &xin) : BASE(idin,xin), dir(0) {}
      prdc_template<BASE>* create(mesh& xin) const {return new prdc_template<BASE>(idnum,xin);}

      int& setdir() {return(dir);}
      void output(std::ostream& fout) {
         BASE::output(fout);
         fout << BASE::idprefix << ".dir" << ": " << dir << std::endl;  
      }
      void input(std::map <std::string,std::string>& inmap) {
         char idntystring[100];
         std::string keyword;
         std::map<std::string,std::string>::const_iterator mi;
         std::istringstream data;
         
         BASE::input(inmap);

         sprintf(idntystring,"%d",idnum);
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
      void sndpositions() { loadbuff(&(x.vrtx[0][0]),1-dir,1-dir +mesh::DIM-2,mesh::DIM); }
      void rcvpositions() { finalrcv(&(x.vrtx[0][0]),1-dir,1-dir +mesh::DIM-2,mesh::DIM); }
};

class curved_analytic : public side_bdry {
   protected:
      virtual FLT hgt(FLT x[mesh::DIM]) {return(0.0);}
      virtual FLT dhgt(int dir, FLT x[mesh::DIM]) {return(1.0);}

   public:      
      /* CONSTRUCTOR */
      curved_analytic(int idin, mesh &xin) : side_bdry(idin,xin) {}
      curved_analytic* create(mesh& xin) const {return new curved_analytic(idnum,xin);}

      void mvpttobdry(int nel,FLT psi, FLT pt[mesh::DIM]);
};

class sinewave : public curved_analytic {
   protected:
      FLT hgt(FLT pt[mesh::DIM]) {
         return(pt[1] -offset -amp*sin(2.*M_PI*pt[0]/lam +phase));
      }
      FLT dhgt(int dir, FLT pt[mesh::DIM]) {
         switch(dir) {
            case(0):
               return(-amp*2.*M_PI/lam*cos(2.*M_PI*pt[0]/lam+phase));
            case(1):
               return(1.0);
         }
         return(1.0);
      }
   public:
      FLT amp, lam, phase, offset;
      sinewave(int inid, mesh &xin) : curved_analytic(inid,xin), amp(0.0), lam(1.0), phase(0.0), offset(0.0) {}
      sinewave* create(mesh& xin) const {return new sinewave(idnum,xin);}

      void output(std::ostream& fout) {         
         curved_analytic::output(fout);
         fout << idprefix << ".amp: " << amp << std::endl;
         fout << idprefix << ".lam: " << lam << std::endl;
         fout << idprefix << ".phase: " << phase << std::endl;
         fout << idprefix << ".offset: " << offset << std::endl;
      }
         
      void input(std::map <std::string,std::string>& inmap) {   
         curved_analytic::input(inmap);

         std::istringstream data(inmap[idprefix+".amp"]);
         if (!(data >> amp)) amp = 0.0;
         data.clear();
         
         data.str(inmap[idprefix +".lam"]);
         if (!(data >> lam)) lam = 1.0;
         data.clear();
         
         data.str(inmap[idprefix +".phase"]);
         if (!(data >> phase)) phase = 0.0;
         data.clear();
         
         data.str(inmap[idprefix +".offset"]);
         if (!(data >> offset)) offset = 0.0;
         data.clear();
      }
      
};

class circle : public curved_analytic {
   public:
      FLT center[mesh::DIM];
      FLT radius;
      FLT hgt(FLT pt[mesh::DIM]) {
         return(radius*radius -pow(pt[0]-center[0],2) -pow(pt[1]-center[1],2));
      }
      FLT dhgt(int dir, FLT pt[mesh::DIM]) {
         return(-2.*(pt[dir]-center[dir]));
      }
      
      circle(int inid, mesh &xin) : curved_analytic(inid,xin), radius(0.5) {center[0] = 0.0; center[1] = 0.0;}
      circle* create(mesh& xin) const {return new circle(idnum,xin);}

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

