/*
 *  rboundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Mon Jun 10 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include"mesh.h"

class r_mesh;

/* A GENERIC PROTOCOL FOR SENDING */
class rcomm {
   public:
      static void send (class comm_boundary *bin,FLT *base,int bgn,int end, int stride) {
         bin->send(base,bgn,end,stride);
      }
      
      static void rcv (class comm_boundary *bin,FLT *base,int bgn,int end, int stride) {
         bin->rcv(base,bgn,end,stride);
      }
};

class no_rcomm {
   public:
      static void send (class comm_boundary *bin,FLT *base,int bgn,int end, int stride) {}
      static void rcv (class comm_boundary *bin,FLT *base,int bgn,int end, int stride) {}
};

/* GENERIC PROTOCALS FOR MOVING BOUNDARIES */
class no_rfix {
   public:
      static void dirichlet(class side_boundary *bin, FLT (*res)[ND]) {}
};

class rfixx {
   public:
      static void dirichlet(class side_boundary *bin, FLT (*res)[ND]) {
      int j,sind;
      
      for(j=0;j<bin->nsd();++j) {
         sind = bin->sd(j);
         res[bin->b().svrtx[sind][0]][0] = 0.0;
         res[bin->b().svrtx[sind][1]][0] = 0.0;
      }
      return;
   }
};

class rfixy {
   public:
      static void dirichlet(class side_boundary *bin, FLT (*res)[ND]) {
      int j,sind;
      
      for(j=0;j<bin->nsd();++j) {
         sind = bin->sd(j);
         res[bin->b().svrtx[sind][0]][1] = 0.0;
         res[bin->b().svrtx[sind][1]][1] = 0.0;
      }
      return;
   }
};

template<class BASE, class SENDX, class SENDY, class FIXX, class FIXY> class rcomm_generic 
   : public BASE {
   
   private:
      r_mesh &x;
      SENDX xsend;
      SENDY ysend;
      FIXX xfix;
      FIXY yfix;
      
   public:
      rcomm_generic(class r_mesh& xin, int type) : BASE(xin,type) , x(x) {}
      inline r_mesh& b() {return(x);}
      void dirichlet(FLT (*res)[ND]) {
         xfix.dirichlet(this,res);
         yfix.dirichlet(this,res);
      }
         
      void sendx(FLT *base,int bgn,int end, int stride) {
         xsend.send(this,base,bgn,end,stride);
      }
      void sendy(FLT *base,int bgn,int end, int stride) {
         ysend.send(this,base,bgn,end,stride);
      }
      void rcvx(FLT *base,int bgn,int end, int stride) {
         xsend.rcv(this,base,bgn,end,stride);
      }
      void rcvy(FLT *base,int bgn,int end, int stride) {
         ysend.rcv(this,base,bgn,end,stride);
      }
      void tadvance() {}
};

template<class BASE, class FIXX, class FIXY> class rgeneric 
   : public BASE {
   
   private:
   	r_mesh &x;
      FIXX xfix;
      FIXY yfix;
      
   public:
      rgeneric(class r_mesh& xin, int type) : BASE(xin,type) , x(xin) {}
      inline r_mesh& b() {return(x);}
      void dirichlet(FLT (*res)[ND]) {
         xfix.dirichlet(this,res);
         yfix.dirichlet(this,res);
      }
      void sendx(FLT *base,int bgn,int end, int stride) {}
      void sendy(FLT *base,int bgn,int end, int stride) {}
      void rcvx(FLT *base,int bgn,int end, int stride) {}
      void rcvy(FLT *base,int bgn,int end, int stride) {}
      void tadvance() {}
};

typedef rcomm_generic<prdx_boundary,rcomm,no_rcomm,rfixx,no_rfix> rprdx;
typedef rcomm_generic<prdy_boundary,no_rcomm,rcomm,no_rfix,rfixy> rprdy;
typedef rcomm_generic<comm_boundary,rcomm,no_rcomm,no_rfix,no_rfix> rcomx;
typedef rcomm_generic<comm_boundary,no_rcomm,rcomm,no_rfix,no_rfix> rcomy;
typedef rgeneric<side_boundary,rfixx,rfixy> rfixd;
typedef rgeneric<curv_boundary,rfixx,rfixy> rfixd_curved;
typedef rcomm_generic<curv_boundary,no_rcomm,rcomx,rfixx,rfixy> rifce;
typedef rgeneric<side_boundary,rfixx,no_rfix> rsymm;

