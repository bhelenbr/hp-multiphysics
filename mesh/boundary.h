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

#ifdef MPISRC
#include <mpi.h>
#endif

#ifndef _boundary_h_
#define _boundary_h_

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#define FLT double
#define EPSILON DBL_EPSILON
#endif

/* GENERIC INTERFACE FOR A B.C. */
class boundary {
   private:
      int idnum;
   
   public:
      enum msg_type {flt_msg, int_msg};
      boundary(int idin) : idnum(idin) {}
      inline int idnty() const {return(idnum);}
      virtual void alloc(int size) {}
      virtual void copy(const boundary& bin) {} 
      virtual void output(std::ostream& fout) {fout << "idnum: " << idnty() << '\n';}
      virtual void input(FILE *fin, FLT grwfac = 1.0) {}
      virtual void getgeometryfrommesh() {}
      virtual void summarize(std::ostream& fout) {fout << "idnum: " << idnty() << '\n';}
      virtual void tadvance() {}
      
      /* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
      virtual bool is_comm() {return(false);}
      virtual bool is_frst() {return(true);}
      virtual void setfrst(bool tf) {}
      virtual int local_cnnct(boundary *bin, int msg_tag) {return 1;}
      virtual void *mysndbuff() {return 0;}
#ifdef MPISRC
      virtual int mpi_cnnct(int proc_tgt, int msg_tag) {return 1;}
#endif
      virtual void setphase(int phase, int msg_tag) {}
      
      /* GENERIC MECHANISM FOR SENDING? */
      virtual void snd(int phase) {}
      virtual int rcv(int phase) {return 1;}
      virtual void sndpositions() {}
      virtual void rcvpositions() {}
};


/* TEMPLATE CLASS TO MAKE A COMMUNCIATION BOUNDARY */
template<class BASE, class MESH> class comm_boundary : public BASE {
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
      comm_boundary(int inid, MESH &xin) : BASE(inid,xin), first(1), maxphase(0), buffsize(0), nlocal_match(0), nmpi_match(0)  {}
      bool is_comm() {return(true);}
      bool is_frst() {return(first);}
      void setfrst(bool tf) {first = tf;}; 
      int& sndsize() {return(msgsize);}
      boundary::msg_type& sndtype() {return(msgtype);}
      void *mysndbuff() {return(sndbuff);}
      int& isndbuf(int n) {return(isndalias[n]);}
      FLT& fsndbuf(int n) {return(fsndalias[n]);}
      
      void alloc(int size);
      void setphase(int phase, int msg_tag);
      int local_cnnct(boundary *bin, int msg_tag); 
#ifdef MPISRC
      int mpi_cnnct(int nproc, int msg_tag);
#endif
      void snd(int phase);
      int rcv(int phase);
};


/* SPECIALIZATION FOR A VERTEX */
template <class MESH> class vrtx_template : public boundary {
   protected:
      MESH &x;
      int v0;
      FLT pt[MESH::DIM];
      
   public:
      /* CONSTRUCTOR */
      vrtx_template(int inid, MESH &xin) : boundary(inid), x(xin) {}
      void alloc(int n) { boundary::alloc(MESH::DIM*n);}
      
      /* ACCESS FUNCTIONS */
      inline MESH& b() {return(x);} // return base mesh
      inline int& v() {return(v0);} // return vertex number in mesh
      
      /* OTHER USEFUL STUFF */
      void copy(const boundary& bin) {
         boundary::copy(bin);
         const vrtx_template<MESH>& temp = dynamic_cast<const vrtx_template<MESH>& >(bin);
         v0 = temp.v0;
         for(int n=0;n<MESH::DIM;++n)
            pt[n] = temp.pt[n];
      }
      void output(std::ostream& fout) {
         boundary::output(fout);
         fout << "point: " << v0 << '\n';
         for(int n=0;n<MESH::DIM;++n)
            fout << pt[n] << ' ';
         fout << "\n";
      }
      void input(FILE *fin, FLT grwfac=1.0) {
         boundary::input(fin);
         std::fscanf(fin,"point: %d\n",&v0);
         for(int n=0;n<MESH::DIM;++n)
            std::fscanf(fin,"%lf ",&pt[n]);
         fscanf(fin,"\n");
      } 
      void summarize(std::ostream& fout) {
         fout << "#VBDRY "<< idnty() << " VRTX " << v() << '\n';
      }
      
      void getgeometryfrommesh() {
         for(int n=0;n<MESH::DIM;++n) {
            pt[n] = b().vrtx[v0][n];
         }
      }
      
      /* MORE SPECIFIC SENDING FOR VERTICES */
      virtual void loadbuff(FLT *base,int bgn,int end, int stride) {}
      virtual void finalrcv(FLT *base,int bgn,int end, int stride) {}
      void sndpositions() {loadbuff(&(b().vrtx[0][0]),0,1,MESH::DIM);}
      void rcvpositions() {finalrcv(&(b().vrtx[0][0]),0,1,MESH::DIM);}
};

template<class MESH> class vcomm : public comm_boundary<vrtx_template<MESH>,MESH> {
   public:
      vcomm(int inid, MESH& xin) : comm_boundary<vrtx_template<MESH>,MESH>(inid,xin) {}
      
      void alloc(int n) {
         vrtx_template<MESH>::alloc(n);
         comm_boundary<vrtx_template<MESH>,MESH>::alloc(MESH::DIM*n);
      }
         
      /* GENERIC VERTEX COMMUNICATIONS */
      void loadbuff(FLT *base,int bgn,int end, int stride) {
         int i,offset;
         
         sndsize()=end-bgn+1;
         sndtype()=flt_msg;
         
         /* LOAD SEND BUFFER */   
         offset = v()*stride +bgn;
         for (i=0;i<end-bgn+1;++i) 
            fsndalias[i] = base[offset+i];
      }
      
      void finalrcv(FLT *base,int bgn,int end, int stride) {
         int i,m,offset;

         /* REMINDER DON'T OVERWRITE SEND BUFFER */
#ifdef MPISRC
         MPI_Status status;
         /* MPI PASSES */
         for(m=0;m<nmpi_match;++m)
            MPI_Wait(&mpi_rcvrqst[m], &status);         
#endif
         
         offset = v()*stride +bgn;
         for(i=0;i<end-bgn+1;++i) {
            base[offset] = fsndalias[i];  // ELIMINATES UNDESIRABLE SIDE/VRTX INTERACTIONS
            for(m=0;m<nlocal_match;++m)
               base[offset] += static_cast<FLT *>(local_rcv_buf[m])[i];

#ifdef MPISRC
            for(m=0;m<nmpi_match;++m)
               base[offset] += static_cast<FLT *>(mpi_rcv_buf[m])[i];
#endif
            base[offset++] /= (1 +nlocal_match +nmpi_match);
         }
      }
};

/* SPECIALIAZATION FOR A SIDE */
template <class MESH> class side_template : public boundary {
   private:
      MESH &x;
      int maxel;
      int nel;
      int *el;
      FLT *s;
      
   public:
      /* CONSTRUCTOR */
      side_template(int inid, MESH &xin) : boundary(inid), x(xin), maxel(0)  {};
      
      /* ACCESS FUNCTIONS */
      inline int mxsz() const {return(maxel);}
      inline int& nsd() {return(nel);}
      inline int& sd(int ind) {return(el[ind]);}
      inline MESH& b() {return(x);} // return base mesh
      inline FLT& spt(int ind) {return(s[ind]);}
      
      /* BASIC B.C. STUFF */
      void alloc(int n);
      void copy(const boundary& bin);
      void output(std::ostream& fout);
      void input(FILE *fin, FLT grwfac=1.0);
      void getgeometryfrommesh();
      
      /* ADDITIONAL STUFF FOR SIDES */
      virtual void swap(int s1, int s2);
      virtual void reorder();
      virtual void mvpttobdry(int nel,FLT psi, FLT pt[MESH::DIM]);
      virtual void findbdrypt(const boundary *tgt,int ntgt,FLT psitgt,int *nout, FLT *psiout);
      virtual void summarize(std::ostream& fout) {
         fout << "#BDRY " << idnty() << " MAX " << mxsz() << " SIDES " << nsd() << '\n';
      }
      
      /* DEFAULT SENDING FOR SIDE VERTICES */
      virtual void loadbuff(FLT *base,int bgn,int end, int stride) {}
      virtual void finalrcv(FLT *base,int bgn,int end, int stride) {}
      void sndpositions() {loadbuff(&(b().vrtx[0][0]),0,MESH::DIM-1,MESH::DIM);}
      void rcvpositions() {finalrcv(&(b().vrtx[0][0]),0,MESH::DIM-1,MESH::DIM); }
};

template<class MESH> class scomm : public comm_boundary<side_template<MESH>,MESH> {
   public:            
      /* CONSTRUCTOR */
      scomm(int inid, MESH& xin) : comm_boundary<side_template<MESH>,MESH>(inid,xin) {}
      
      /* INITIALIZE STORAGE */
      void alloc(int n) {
         side_template<MESH>::alloc(n);
         comm_boundary<side_template<MESH>,MESH>::alloc(MESH::DIM*n);
      }

      
      /* GENERIC VERTEX COMMUNICATIONS */
      virtual void loadbuff(FLT *base,int bgn,int end, int stride) {
         int j,k,count,sind,offset;

         count = 0;
         for(j=0;j<nsd();++j) {
            sind = sd(j);
            offset = b().svrtx[sind][0]*stride;
            for (k=bgn;k<=end;++k) 
               fsndalias[count++] = base[offset+k];
         }
         offset = b().svrtx[sind][1]*stride;
         for (k=bgn;k<=end;++k) 
            fsndalias[count++] = base[offset+k]; 
            
         sndsize() = count;
         sndtype() = boundary::flt_msg;
      }
      
      
      
      
      virtual void finalrcv(FLT *base,int bgn,int end, int stride) {
         int j,k,m,count,offset,sind;
#ifdef MPISRC
         MPI_Status status;
#endif
         /* ASSUMES REVERSE ORDERING OF SIDES */
         /* WON'T WORK IN 3D */
         
         /* RELOAD FROM BUFFER */
         /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
         /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */
         count = 0;
         for(j=0;j<nsd();++j) {
            sind = sd(j);
            offset = b().svrtx[sind][0]*stride;
            for (k=bgn;k<=end;++k)
               base[offset+k] = fsndalias[count++];
         }
         offset = b().svrtx[sind][1]*stride;
         for (k=bgn;k<=end;++k)
            base[offset+k] = fsndalias[count++];
         
         for(m=0;m<nlocal_match;++m) {            
            count = 0;
            for(j=nsd()-1;j>=0;--j) {
               sind = sd(j);
               offset = b().svrtx[sind][1]*stride;
               for (k=bgn;k<=end;++k) 
                  base[offset+k] += static_cast<FLT *>(local_rcv_buf[m])[count++];
            }
            offset = b().svrtx[sind][0]*stride;
            for (k=bgn;k<=end;++k) 
                  base[offset+k] += static_cast<FLT *>(local_rcv_buf[m])[count++];            
         }
         
      
#ifdef MPISRC
         /* MPI PASSES */
         for(m=0;m<nmpi_match;++m) {
            MPI_Wait(&mpi_rcvrqst[m], &status);
            count = 0;
            for(j=nsd()-1;j>=0;--j) {
               sind = sd(j);
               offset = b().svrtx[sind][1]*stride;
               for (k=bgn;k<=end;++k) 
                  base[offset+k] += static_cast<FLT *>(mpi_rcv_buf[m])[count++];
            }
            offset = b().svrtx[sind][0]*stride;
            for (k=bgn;k<=end;++k) 
                  base[offset+k] += static_cast<FLT *>(mpi_rcv_buf[m])[count++];
         }
#endif

         count = 0;
         for(j=0;j<nsd();++j) {
            sind = sd(j);
            offset = b().svrtx[sind][0]*stride;
            for (k=bgn;k<=end;++k)
               base[offset+k] /= (1 +nlocal_match +nmpi_match);
         }
         offset = b().svrtx[sind][1]*stride;
         for (k=bgn;k<=end;++k)
            base[offset+k] /= (1 +nlocal_match +nmpi_match);
      }
};

/* CAN DO PERIODIC IN X OR Y (dir = 0/1) IN 2D */
/* CAN DO PERIODIC IN X OR Z (dir = 0/1) IN 3D */
template<class BASE, class MESH> class prdc_template : public BASE {
   protected:
      int dir;
   public:      
      /* CONSTRUCTOR */
      prdc_template(int idin, MESH &xin) : BASE(idin,xin), dir(0) {}
      int& setdir() {return(dir);}
      
      /* SEND/RCV VRTX POSITION */
      void sndpositions() { loadbuff(&(b().vrtx[0][0]),1-dir,1-dir +MESH::DIM-2,MESH::DIM); }
      void rcvpositions() { finalrcv(&(b().vrtx[0][0]),1-dir,1-dir +MESH::DIM-2,MESH::DIM); }
};

template<class BASE, class MESH> class curv_template : public BASE {
   protected:
      virtual FLT hgt(FLT x[MESH::DIM]) = 0;
      virtual FLT dhgt(int dir, FLT x[MESH::DIM]) = 0;
   public:      
      /* CONSTRUCTOR */
      curv_template(int idin, MESH &xin) : BASE(idin,xin) {}
      void mvpttobdry(FLT pt[MESH::DIM]) {
         int iter,n;
         FLT mag, delt_dist;
            
         /* FOR AN ANALYTIC SURFACE */
         iter = 0;
         do {
            mag = 0.0;
            for(n=0;n<MESH::DIM;++n)
               mag += pow(dhgt(n,pt),2);
            mag = sqrt(mag);
            delt_dist = -hgt(pt)/mag;
            for(n=0;n<MESH::DIM;++n)
               pt[n] += delt_dist*dhgt(n,pt)/mag;
            if (++iter > 100) {
               *b().log << "iterations exceeded curved boundary " << idnty() << ' ' << pt[0] << ' ' << pt[1] << '\n';
               exit(1);
            }
         } while (fabs(delt_dist) > 10.*EPSILON);
         
         return;
      }
};

#ifdef CAPRI
#include <capri.h>
#endif

#include <rbd.h>

/* 3D to 2D BOUNDARY */
template<class MESH> class threetotwo : public scomm<MESH> {
   int face;
   public:
      threetotwo(int inid, MESH& xin) : scomm<MESH>(inid,xin) {
         if ((inid&0xFFFF)==16) face = 7;
         else face = 8;
      }
      void sndpositions() {loadbuff(&(b().vrtx[0][0]),0,1,3);}
      void rcvpositions() {
         int i,sind,v0;
         for(i=0;i<nsd();++i) {
            sind = sd(i);
            v0 = b().svrtx[sind][0];
            b().vrtx[v0][2] = 0.0;
         }
         v0 = b().svrtx[sind][1];
         b().vrtx[v0][2] = 0.0;
      }
      
      void tadvance() {
         int i,sind,v0,n,iter;
         double x[3],x1[3];
#ifdef CAPRI
         int icode;
         double t = -100.0;
         double uv[2] = {-100.0,100.0};
#endif
                  
         for(i=0;i<nsd();++i) {
            sind = sd(i);
            v0 = b().svrtx[sind][0];
            for(n=0;n<3;++n)
               x[n] = b().vrtx[v0][n];
            rigidbodyrmv(x);
            for (iter = 0; iter < 5; ++iter) {
               for(n=0;n<3;++n)
                  x1[n] = x[n];
               rigidbodyadd(x);
               x1[2] -= x[2];
#ifdef CAPRI
               icode = gi_qNearestOnFace(1, face, x1, uv, x);
#endif
            }
            rigidbodyadd(x);
            for(n=0;n<3;++n)
               b().vrtx[v0][n] = x[n];
         }
         v0 = b().svrtx[sind][1];
         for(n=0;n<3;++n)
            x[n] = b().vrtx[v0][n];
         rigidbodyrmv(x);
         for (iter = 0; iter < 5; ++iter) {
            for(n=0;n<3;++n)
               x1[n] = x[n];
            rigidbodyadd(x);
            x1[2] -= x[2];
#ifdef CAPRI
            icode = gi_qNearestOnFace(1, face, x1, uv, x);
#endif
         }
         rigidbodyadd(x);
         for(n=0;n<3;++n)
            b().vrtx[v0][n] = x[n];
      }
};

template<class MESH> class twotothree : public scomm<MESH> {
   public:
      twotothree(int inid, MESH& xin) : scomm<MESH>(inid,xin) {}
      void finalrcv(FLT *base,int bgn,int end, int stride) {}
      void rcvpositions() {
         FLT *base = &(b().vrtx[0][0]);
         int bgn = 0, end = 1, stride = 2;
         int j,k,m,count,offset,sind;
#ifdef MPISRC
         MPI_Status status;
#endif
         /* NOT TOTALLY SURE HOW 1 TO MANY WILL WORK */
         /* THIS IS ONLY FOR ONE TO ONE MATCHES */
         
         /* ASSUMES REVERSE ORDERING OF SIDES */
         /* WON'T WORK IN 3D */
         for(m=0;m<nlocal_match;++m) {   
            count = 0;
            for(j=nsd()-1;j>=0;--j) {
               sind = sd(j);
               offset = b().svrtx[sind][1]*stride;
               for (k=bgn;k<=end;++k)
                  base[offset+k] = static_cast<FLT *>(local_rcv_buf[m])[count++];
            }
            offset = b().svrtx[sind][0]*stride;
            for (k=bgn;k<=end;++k) 
                  base[offset+k] = static_cast<FLT *>(local_rcv_buf[m])[count++];            
         }
         
      
#ifdef MPISRC
         /* MPI PASSES */
         for(m=0;m<nmpi_match;++m) {
            MPI_Wait(&mpi_rcvrqst[m], &status);
            count = 0;
            for(j=nsd()-1;j>=0;--j) {
               sind = sd(j);
               offset = b().svrtx[sind][1]*stride;
               for (k=bgn;k<=end;++k) 
                  base[offset+k] = static_cast<FLT *>(mpi_rcv_buf[m])[count++];
            }
            offset = b().svrtx[sind][0]*stride;
            for (k=bgn;k<=end;++k) 
                  base[offset+k] = static_cast<FLT *>(mpi_rcv_buf[m])[count++];
         }
#endif
      }
};

#ifdef CAPRI
template<class MESH> class capri_edge : public side_template<MESH> {
   private:
      int vol, edge;
   public:
      capri_edge(int inid, MESH& xin) : side_template<MESH>(inid,xin), vol(1), edge(inid&0xFFFF) {}

      void tadvance() {
         int i,n,icode,sind,v0;
         double x1[3],x2[3];
         double t = -10.0;
         
         for(i=0;i<nsd();++i) {
            sind = sd(i);
            v0 = b().svrtx[sind][0];
            icode = gi_qPointOnEdge(1, edge, spt(i), b().vrtx[v0], 0, x1, x1);
            rigidbodyadd(b().vrtx[v0]);
         }
         v0 = b().svrtx[sind][1];
         icode = gi_qPointOnEdge(1, edge, spt(nsd()), b().vrtx[v0], 0, x1, x1);
         rigidbodyadd(b().vrtx[v0]);
         
         return;
      }
      
      void getgeometryfrommesh() {
         int i,n,icode,sind,v0;
         double x[3],x1[3];
         int cpri_npt;
         double *cpri_pts, *t;
         double dist[2];
         
         icode = gi_dTesselEdge(1, edge, &cpri_npt, &cpri_pts, &t);
         
         if (cpri_npt != nsd()+1) {
            *b().log << "error in cpri_edge: " << edge << ' ' << cpri_npt << ' ' << nsd()+1 << std::endl;
            exit(1);
         }
         
         /* FIGURE OUT DIRECTION */
         icode = gi_qPointOnEdge(1, edge, t[0], x, 0, x1, x1);
         dist[0] = 0.0;
         for(n=0;n<3;++n)
            dist[0] += fabs(x[n] - b().vrtx[b().svrtx[sd(0)][0]][n]);
         
         dist[1] = 0.0;
         for(n=0;n<3;++n)
            dist[1] += fabs(x[n] - b().vrtx[b().svrtx[sd(nsd()-1)][1]][n]); 
            
         if (dist[0] < dist[1]) { 
            for(i=0;i<nsd()+1;++i)
               spt(i) = t[i];
         }
         else {
            for(i=0;i<nsd()+1;++i)
               spt(i) = t[cpri_npt-1-i];
         }
      }
};
#endif
   

#include "boundary.cpp"

#endif

