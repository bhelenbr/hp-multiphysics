/*
 *  boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Fri Jun 07 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef _mesh_h_
ONLY INCLUDE FROM WITH MESH
#endif

#include <stdio.h>
#include <utilities.h>

class vrtx_boundary {
   private:
      mesh &x;
      int idnum;
      FLT pt[ND];
      int v0;
      
   public:
      /* CONSTRUCTOR */
      vrtx_boundary(mesh &xin, int inid) : x(xin), idnum(inid) { alloc(); }
      
      /* ACCESS FUNCTIONS */
      inline int idnty() const {return(idnum);}
      inline int& v() {return(v0);}
      inline mesh& b() {return(x);} // return base mesh
      
      virtual void alloc() {}
      /* OTHER USEFUL STUFF */
      virtual void copy(const vrtx_boundary &bin) {
         v0 = bin.v0;
         for(int n=0;n<ND;++n)
            pt[n] = bin.pt[n];
      }
      virtual void output(FILE *out) {
         fprintf(out,"vbdry idnum: %d point: %d\n",idnum,v0);
         for(int n=0;n<ND;++n)
            fprintf(out,"%f ",pt[n]);
         fprintf(out,"\n");
      }
      virtual void input(FILE *in){
         fscanf(in,"vbdry idnum: %*d point: %d\n",&v0);
         for(int n=0;n<ND;++n)
            fscanf(in,"%lf ",&pt[n]);
         fscanf(in,"\n");
      } 
      virtual void getgeometryfrommesh() {
         for(int n=0;n<ND;++n) {
            pt[n] = b().vrtx[v0][n];
         }
      }
      virtual void summarize() {
         std::printf("#VBDRY %d VRTX %d\n",idnty(),v());
      }
      
      /* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
      virtual int match(vrtx_boundary *in) {return 1;}
      virtual void setphase(int phase, int match) {}
      virtual int send(int phase, FLT *base,int bgn,int end, int stride) {return(1);}
      virtual void rcv(int phase, FLT *base,int bgn,int end, int stride) {}
      virtual int sendposition(int phase) {return(1);}
      virtual void rcvposition(int phase) {return;}
};

class vcom_boundary : public vrtx_boundary {
   public:
      static const int maxmatch = 8;
      int nmatch;
      class vcom_boundary *vmatch[maxmatch];
      int myphase[maxmatch],maxphase;
      FLT *vbuff;
      int msgsize;
      
      vcom_boundary(mesh &xin, int inid) : vrtx_boundary(xin, inid), nmatch(0), maxphase(0) {}
      void alloc() {vbuff = new FLT[4];}
      
      /* ZERO FIRST VALUE */
      void setphase(int phase, int match) {myphase[match] = phase; maxphase = MAX(maxphase,phase);}

      /* MATCH BOUNDARIES */
      int match(vrtx_boundary *in) {
         if (in->idnty() == idnty()) {
            vmatch[nmatch] = dynamic_cast<vcom_boundary *>(in);
            setphase(0,nmatch); // NOT SURE HOW TO DO THIS YET
            ++nmatch;
            return(1);
         }
         return(0);
      }
      
      /* SEND/RCV VRTX POSITION */
      int send(int phase, FLT *base,int bgn,int end, int stride);
      void rcv(int phase, FLT *base,int bgn,int end, int stride);
      int sendpositions(int phase) { return(send(phase,&(b().vrtx[0][0]),0,1,2));}
      void rcvpositions(int phase) {rcv(phase,&(b().vrtx[0][0]),0,1,2); }
   
};

/* INTERFACE FOR A BOUNDARY CONDITION */   
class side_boundary {
   private:
      mesh &x;
      int idnum;
      int maxel;
      int nel;
      int *el;
      FLT (*s)[2];
      
   public:
      /* CONSTRUCTOR */
      side_boundary(mesh &xin, int inid) : x(xin), idnum(inid), maxel(0)  {};
      
      /* ACCESS FUNCTIONS */
      inline int idnty() const {return(idnum);}
      inline int mxsz() const {return(maxel);}
      inline int& nsd() {return(nel);}
      inline int& sd(int ind) {return(el[ind]);}
      inline mesh& b() {return(x);} // return base mesh
      
      /* OTHER USEFUL STUFF */
      virtual void alloc(int n);
      virtual void copy(const side_boundary &bin);
      virtual void output(FILE *out);
      virtual void input(FILE *in, FLT grwfac);
      virtual void swap(int s1, int s2);
      virtual void reorder();
      virtual void mvpttobdry(int nel,FLT psi, FLT pt[2]);
      virtual void getgeometryfrommesh();
      virtual void findbdrypt(const class side_boundary *tgt,int ntgt,FLT psitgt,int *nout, FLT *psiout);
      virtual void summarize() {
         std::printf("#BDRY %d MAX %d SIDES %d\n",idnty(),mxsz(),nsd());
      }
      
      /* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
      virtual int isfrst() {return(1);}
      virtual int match(side_boundary *in) {return 1;}
      virtual void setphase(int phase) {}
      virtual int send(int phase, FLT *base,int bgn,int end, int stride) {return(1);}
      virtual void rcv(int phase, FLT *base,int bgn,int end, int stride) {}
      virtual int sendpositions(int phase) {return(1);}
      virtual void rcvpositions(int phase) {return;}
};

class comm_boundary : public side_boundary {
   public:
      /* ADDITIONAL STUFF FOR COMMUNICATION BOUNDARIES */
      class comm_boundary *bdrymatch;
      int frst;
      int myphase;
      FLT *sbuff;
      int msgsize;
      
      /* CONSTRUCTOR */
      comm_boundary(mesh &xin, int inid) : side_boundary(xin,inid) , frst(0), myphase(0) {}
      
      /* INITIALIZE STORAGE */
      void alloc(int nside) {
         side_boundary::alloc(nside);
         sbuff = new FLT[2*ND*mxsz()];
      }
      
      /* ZERO FIRST VALUE */
      void setphase(int i) {myphase = i;}

      /* MATCH BOUNDARIES */
      int match(side_boundary *in) {
         if (in->idnty() == idnty()) {
            bdrymatch = dynamic_cast<comm_boundary *>(in);
            frst = 1 -bdrymatch->isfrst();
            return(1);
         }
         return(0);
      }
      
      /* TEST FIRST */
      int isfrst() {return(frst);}
      
      /* SEND/RCV VRTX POSITION */
      int send(int phase, FLT *base,int bgn,int end, int stride);
      void rcv(int phase, FLT *base,int bgn,int end, int stride);
      int sendpositions(int phase) { return(send(phase,&(b().vrtx[0][0]),0,1,2));}
      void rcvpositions(int phase) {rcv(phase,&(b().vrtx[0][0]),0,1,2); }
      
};

class prdx_boundary : public comm_boundary {
   public:      
      /* CONSTRUCTOR */
      prdx_boundary(mesh &xin, int idin) : comm_boundary(xin,idin) {setphase(0);}

      /* SEND/RCV Y VRTX POSITION */
      int sendpositions(int phase) { return(send(phase,&(b().vrtx[0][0]),1,1,2)); }
      void rcvpositions(int phase) { rcv(phase,&(b().vrtx[0][0]),1,1,2); }
};

class prdy_boundary : public comm_boundary {
   public:      
      /* CONSTRUCTOR */
      prdy_boundary(mesh &xin, int idin) : comm_boundary(xin,idin) {setphase(1);}

      /* SEND/RCV X VRTX POSITION */
      int sendpositions(int phase){ return(send(phase,&(b().vrtx[0][0]),0,0,2)); }
      void rcvpositions(int phase) { rcv(phase,&(b().vrtx[0][0]),0,0,2); }
};

template<class BASE> class curv_template : public BASE {
   protected:
      virtual FLT hgt(FLT x[ND]) = 0;
      virtual FLT dhgt(int dir, FLT x[ND]) = 0;
   public:      
      /* CONSTRUCTOR */
      curv_template(mesh &xin, int idin) : BASE(xin,idin) {}
      void mvpttobdry(int nel,FLT psi, FLT pt[ND]);
};

typedef curv_template<side_boundary> curv_boundary;
typedef curv_template<comm_boundary> ifce_boundary;

