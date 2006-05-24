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
#include <input_map.h>

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
      virtual boundary* create(mesh &xin) const = 0;
      virtual void alloc(int n) {}
      virtual void output(std::ostream& fout) {
         fout << idprefix << ".type: " << mytype << std::endl;         
      }
      virtual void input(input_map& bdrydata) {}
      
      /* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
      enum msg_type {flt_msg, int_msg};
      enum groups {all,partitions,manifolds};
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
      virtual int local_cnnct(boundary *bin, int msg_tag) {return 1;}
#ifdef MPISRC
      virtual int mpi_cnnct(int proc_tgt, int msg_tag) {return 1;}
#endif
      virtual int& matchphase(int group, int matchnum) {return(dummy.idum=0);}
      virtual void resize_buffers(int size) {}
      virtual void *psndbuf() {return(&dummy);}
      virtual int& isndbuf(int indx) {return(dummy.idum);}
      virtual FLT& fsndbuf(int indx) {return(dummy.fdum);}
      virtual int& ircvbuf(int m,int indx) {return(dummy.idum);}
      virtual FLT& frcvbuf(int m,int indx) {return(dummy.fdum);}
      virtual int& sndsize() {return(dummy.idum=0);}
      virtual boundary::msg_type& sndtype() {return(dummy.mdum);}
      virtual void comm_prepare(int group, int phase) {}
      virtual void comm_transmit(int group, int phase) {}
      virtual int comm_wait(int group, int phase) {return 1;}
      virtual int comm_nowait(int group, int phase) {return 1;}
      virtual void master_slave_prepare(int group=0) {}
      virtual void master_slave_transmit(int group=0) {}
      virtual void master_slave_wait(int group=0) {}
      virtual void master_slave_nowait(int group=0) {}
      virtual void slave_master_prepare(int group=0) {}
      virtual void slave_master_transmit(int group=0) {}
      virtual void slave_master_wait(int group=0) {}
      virtual void slave_master_nowait(int group=0) {}
      virtual ~boundary() {}
   protected:
      int excpt;
};

class vgeometry_interface {
   public:
      virtual void mvpttobdry(TinyVector<FLT,mesh::ND> &pt) {}
      virtual ~vgeometry_interface() {}
};

class vgeometry_pointer {
   public:
      vgeometry_interface *solution_data;
};


/** \brief Specialization for a vertex 
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for a vertex boundary 
 */
class vrtx_bdry : public boundary, public vgeometry_interface {
   public:
      mesh &x;
      TinyVector<int,2> sbdry;
      int v0;
      
      /* CONSTRUCTOR */
      vrtx_bdry(int intype, mesh &xin) : boundary(intype), x(xin) {idprefix = "v" +idprefix; mytype="plain";}
      vrtx_bdry(const vrtx_bdry &inbdry, mesh &xin) : boundary(inbdry.idnum), x(xin)  {idprefix = inbdry.idprefix; mytype = inbdry.mytype; sbdry = inbdry.sbdry;}

      /* OTHER USEFUL STUFF */
      vrtx_bdry* create(mesh &xin) const { return(new vrtx_bdry(*this,xin));}
      virtual void copy(const vrtx_bdry& bin) {
         v0 = bin.v0;
         sbdry = bin.sbdry;
      }
      virtual void vloadbuff(int group,FLT *base,int bgn,int end, int stride) {}
      virtual void vfinalrcv(int group,int phase,FLT *base,int bgn,int end, int stride) {}
      virtual void mvpttobdry(TinyVector<FLT,mesh::ND> &pt) {}
      virtual void loadpositions() {vloadbuff(all,&(x.vrtx(0)(0)),0,mesh::ND-1,mesh::ND);}
      virtual void rcvpositions(int phase) {vfinalrcv(all,phase,&(x.vrtx(0)(0)),0,mesh::ND-1,mesh::ND);}
};

class sgeometry_interface {
   public:
      virtual void mvpttobdry(int nel, FLT psi, TinyVector<FLT,mesh::ND> &pt) {}
      virtual ~sgeometry_interface() {}
};

class sgeometry_pointer {
   public:
      sgeometry_interface *solution_data;
};


/** \brief Specialization for a side 
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for a side boundary 
 */
class side_bdry : public boundary, public sgeometry_interface {
   public:
      mesh &x;
      TinyVector<int,2> vbdry;
      int maxel;
      int nel;
      Array<int,1> el;
      
      /* CONSTRUCTOR */
      side_bdry(int inid, mesh &xin) : boundary(inid), x(xin), maxel(0)  {idprefix = "s" +idprefix; mytype="plain";}
      side_bdry(const side_bdry &inbdry, mesh &xin) : boundary(inbdry.idnum), x(xin), maxel(0)  {idprefix = inbdry.idprefix; mytype = inbdry.mytype; vbdry = inbdry.vbdry;}
      
      /* BASIC B.C. STUFF */
      void alloc(int n);
      side_bdry* create(mesh &xin) const {
         return(new side_bdry(*this,xin));
      }
      virtual void copy(const side_bdry& bin);
      
      /* ADDITIONAL STUFF FOR SIDES */
      virtual void swap(int s1, int s2);
      virtual void reorder();
      virtual block::ctrl mgconnect(block::ctrl ctrl_message, Array<mesh::transfer,1> &cnnct, const class mesh& tgt, int bnum);
      virtual void mvpttobdry(int nel, FLT psi, TinyVector<FLT,mesh::ND> &pt);
      virtual void findbdrypt(const TinyVector<FLT,2> xpt, int &sidloc, FLT &psiloc) const;
      
      /* DEFAULT SENDING FOR SIDE VERTICES */
      virtual void vloadbuff(int group,FLT *base,int bgn,int end, int stride) {}
      virtual void vfinalrcv(int group,int phase,FLT *base,int bgn,int end, int stride) {}
      virtual void sloadbuff(int group,FLT *base,int bgn,int end, int stride) {}
      virtual void sfinalrcv(int group,int phase,FLT *base,int bgn,int end, int stride) {}
      virtual void loadpositions() {vloadbuff(all,&(x.vrtx(0)(0)),0,mesh::ND-1,mesh::ND);}
      virtual void rcvpositions(int phase) {vfinalrcv(all,phase,&(x.vrtx(0)(0)),0,mesh::ND-1,mesh::ND);}
};


#endif


