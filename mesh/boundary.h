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

/* FAT INTERFACE FOR A BOUNDARY CONDITION */

class side_boundary {
   private:
      mesh &x;
      const int idnum;
      int maxel;
      int nel;
      int *el;
      FLT (*s)[2];
      
   public:
      /* PUBLIC VIRTUAL FUNCTIONS */
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
      virtual side_boundary* reorder();
      virtual void mvpttobdry(int nel,FLT psi, FLT pt[2]);
      virtual void getgeometryfrommesh();
      virtual void findbdrypt(const class side_boundary *tgt,int ntgt,FLT psitgt,int *nout, FLT *psiout);
      virtual void length() {}
      virtual void summarize() {
         std::printf("#BDRY %d MAX %d SIDES %d\n",idnty(),mxsz(),nsd());
      }
      
      /* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
      virtual void init_comm_buf(int factor) {}
      virtual void zerofrst() {return;}
      virtual int isfrst() {return(1);}
      virtual int match(side_boundary *in) {return 1;}
      virtual void sendpositions() {return;}
      virtual void rcvpositions() {return;}
      
      /* VIRTUAL FUNCTIONS FOR BOUNDARY DEFORMATION */
      virtual void dirichlet(FLT (*)[ND]) {}
      virtual void fixdx2(FLT (*)[ND]) {}
      virtual void sendx(FLT *base,int bgn,int end, int stride) {}
      virtual void sendy(FLT *base,int bgn,int end, int stride) {}
      virtual void rcvx(FLT *base,int bgn,int end, int stride) {}
      virtual void rcvy(FLT *base,int bgn,int end, int stride) {}
      virtual void tadvance() {}
      
      /* HP BOUNDARY FUNCTIONS */
      virtual void setcrv() {}
      virtual void hpinput(FILE *) {}
      virtual void hpoutput(FILE *) {}
      virtual void setbcinfo() {}
      virtual void loadcrv(int ind, FLT *u[ND], int offset, int ibd = 0) {}
      
      /* HP MGRID BOUNDARY FUNCTIONS */
      virtual void addbflux() {}
      virtual void vsndx() {}
      virtual void vsndy() {}
      virtual void vrcvx() {}
      virtual void vrcvy() {}
      virtual void vzero() {}
      virtual void ssnd(int m) {}
      virtual void srcv(int m) {}
      virtual void szero(int m) {}
      virtual void tstep_sndx() {}
      virtual void tstep_rcvx() {}
      virtual void tstep_sndy() {}
      virtual void tstep_rcvy() {}
      virtual void setinflow() {}
      
      /* FOR DYNAMIC BOUNDARIES */
      virtual void tstep1() {}
      virtual void tstep2() {}
      virtual void minvrt1() {}
      virtual void minvrt2() {}
      virtual void getfres() {}
      virtual void getcchng() {}
      virtual void nstage1() {}
      virtual void nstage2() {}
      virtual void vrttoug() {}
      virtual void ugtovrt() {}
};

class comm_boundary : public side_boundary {
   public:
      /* ADDITIONAL STUFF FOR COMMUNICATION BOUNDARIES */
      class comm_boundary *bdrymatch;
      int frst;
      FLT *sbuff;
      int msgsize;
      
      /* CONSTRUCTOR */
      comm_boundary(mesh &xin, int inid) : side_boundary(xin,inid) , frst(0) {}
      
      /* INITIALIZE STORAGE */
      void init_comm_buf(int factor) {sbuff = new FLT[factor*mxsz()];}
      
      /* ZERO FIRST VALUE */
      void zerofrst() {frst = 0;}
      
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
      
      /* SOME GENERIC VERTEX COMMUNICATION FUNCTIONS */
      /* SHOULD NOT BE CALLED DIRECTLY FROM MESH */
      void send (FLT *base,int bgn,int end, int stride);
      void rcv (FLT *base,int bgn,int end, int stride);
      
      /* SEND/RCV VRTX POSITION */
      void sendpositions() { send(&(b().vrtx[0][0]),0,1,2);}
      void rcvpositions() {rcv(&(b().vrtx[0][0]),0,1,2); }
      
};

class prdx_boundary : public comm_boundary {
   public:      
      /* CONSTRUCTOR */
      prdx_boundary(mesh &xin, int idin) : comm_boundary(xin,idin) {}

      /* SEND/RCV Y VRTX POSITION */
      void sendpositions() { send(&(b().vrtx[0][0]),1,1,2); }
      void rcvpositions() { rcv(&(b().vrtx[0][0]),1,1,2); }
};

class prdy_boundary : public comm_boundary {
   public:      
      /* CONSTRUCTOR */
      prdy_boundary(mesh &xin, int idin) : comm_boundary(xin,idin) {}

      /* SEND/RCV X VRTX POSITION */
      void sendvrtxposition(){ send(&(b().vrtx[0][0]),0,0,2); }
      void rcvvrtxposition() { rcv(&(b().vrtx[0][0]),0,0,2); }
};

class curv_boundary : public side_boundary {
   protected:
      virtual FLT hgt(FLT x[ND]) = 0;
      virtual FLT dhgt(int dir, FLT x[ND]) = 0;
   public:      
      /* CONSTRUCTOR */
      curv_boundary(mesh &xin, int idin) : side_boundary(xin,idin) {}
      void mvpttobdry(int nel,FLT psi, FLT pt[ND]);
};


      