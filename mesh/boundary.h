/*
 *  boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Fri Jun 07 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include "mesh.h"
#include <utilities.h>
#include <stdio.h>
#include <map>
#include <typeinfo>

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

/* GENERIC INTERFACE FOR A B.C. */
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
      }
      virtual boundary* create(mesh &xin) const = 0;
      virtual void alloc(int n) {}
      virtual void copy(const boundary& bin) {}
      virtual void output(std::ostream& fout) {
         fout << idprefix << ".type: " << mytype << std::endl;         
      }
      virtual void input(std::map<std::string,std::string>& bdrydata) {}
      virtual void setupcoordinates() {}
      
      /* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
      enum msg_type {flt_msg, int_msg};
      union  {
         bool bdum;
         int idum;
         FLT fdum;
         msg_type mdum;
      } dummy;
      virtual bool is_comm() {return(false);}
      virtual bool& is_frst() {return(dummy.bdum=true);}
      virtual int& group() {return(dummy.idum=0);}
      virtual int matches() {return(0);}
      virtual int local_cnnct(boundary *bin, int msg_tag) {return 1;}
#ifdef MPISRC
      virtual int mpi_cnnct(int proc_tgt, int msg_tag) {return 1;}
#endif
      virtual void setphase(int phase, int msg_tag) {}
      virtual void resize_buffers(int size) {}
      virtual void *psndbuf() {return(&dummy);}
      virtual int& isndbuf(int indx) {return(dummy.idum);}
      virtual FLT& fsndbuf(int indx) {return(dummy.fdum);}
      virtual int& sndsize() {return(dummy.idum=0);}
      virtual boundary::msg_type& sndtype() {return(dummy.mdum);}
      virtual void comm_prepare(int phase) {}
      virtual int comm_transmit(int phase) {return 1;}
      virtual void comm_wait(int phase) {}
      virtual void master_slave_prepare() {}
      virtual void master_slave_transmit() {}
      virtual void master_slave_wait() {}
      virtual void slave_master_prepare() {}
      virtual void slave_master_transmit() {}
      virtual void slave_master_wait() {}
      virtual void sndpositions() {}
      virtual void rcvpositions() {}
      virtual ~boundary() {}
};

/* SPECIALIZATION FOR A VERTEX */
class vrtx_bdry : public boundary {
   
   public:
      mesh &x;
      int v0;
      
      /* CONSTRUCTOR */
      vrtx_bdry(int intype, mesh &xin) : boundary(intype), x(xin) {idprefix = "v" +idprefix; mytype="plain";}
      
      /* OTHER USEFUL STUFF */
      vrtx_bdry* create(mesh &xin) const { return(new vrtx_bdry(idnum,xin));}
      void copy(const boundary& bin) {
         boundary::copy(bin);
         v0 = dynamic_cast<const vrtx_bdry&>(bin).v0;
      }
      
      /* MORE SPECIFIC SENDING FOR VERTICES */
      virtual void loadbuff(FLT *base,int bgn,int end, int stride) {}
      virtual void finalrcv(FLT *base,int bgn,int end, int stride) {}
      void sndpositions() {loadbuff(&(x.vrtx[0][0]),0,1,mesh::DIM);}
      void rcvpositions() {finalrcv(&(x.vrtx[0][0]),0,1,mesh::DIM);}
};


/* SPECIALIAZATION FOR A SIDE */
class side_bdry : public boundary {
   public:
      mesh &x;
      int maxel;
      int nel;
      int *el;
      
      /* CONSTRUCTOR */
      side_bdry(int inid, mesh &xin) : boundary(inid), x(xin), maxel(0)  {idprefix = "s" +idprefix; mytype="plain";}
      
      /* BASIC B.C. STUFF */
      void alloc(int n);
      virtual side_bdry* create(mesh &xin) const {return(new side_bdry(idnum,xin));}
      void copy(const boundary& bin);
      
      /* ADDITIONAL STUFF FOR SIDES */
      virtual void swap(int s1, int s2);
      virtual void reorder();
      virtual void mvpttobdry(int nel,FLT psi, FLT pt[mesh::DIM]);
      virtual block::ctrl mgconnect(int excpt, mesh::transfer *cnnct, const class mesh& tgt, int bnum);
      virtual void findbdryside(FLT *xpt, int &sidloc, FLT &psiloc) const;
      
      /* DEFAULT SENDING FOR SIDE VERTICES */
      virtual void loadbuff(FLT *base,int bgn,int end, int stride) {}
      virtual void finalrcv(FLT *base,int bgn,int end, int stride) {}
      void sndpositions() {loadbuff(&(x.vrtx[0][0]),0,mesh::DIM-1,mesh::DIM);}
      void rcvpositions() {finalrcv(&(x.vrtx[0][0]),0,mesh::DIM-1,mesh::DIM); }
};

class vtype {
   public:
      static const int ntypes = 3;
      enum ids {plain=1,comm,prdc};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

class stype {
   public:
      static const int ntypes = 7;
      enum ids {plain=1, comm, prdc, sinewave, circle, spline, partition};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

#endif

