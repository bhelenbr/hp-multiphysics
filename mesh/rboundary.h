/*
 *  rboundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Mon Jun 10 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
 #ifndef _r_mesh_h_
ONLY INCLUDE FROM WITH R_MESH.H
#endif

/* VIRTUAL FUNCTIONS */
class rbdry_interface {
   public:
      /* VIRTUAL FUNCTIONS FOR BOUNDARY DEFORMATION */
      virtual void dirichlet(FLT (*)[mesh::ND]) {}
#ifdef FOURTH
      virtual void fixdx2(FLT (*)[mesh::ND]) {}
#endif
      virtual void tadvance() {}
};

/* GENERIC PROTOCALS FOR MOVING BOUNDARIES */
class no_rfix {
   public:
      static void dirichlet(class side_boundary *bin, FLT (*res)[mesh::ND]) {}
};

class rfixx {
   public:
      static void dirichlet(class side_boundary *bin, FLT (*res)[mesh::ND]) {
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
      static void dirichlet(class side_boundary *bin, FLT (*res)[mesh::ND]) {
      int j,sind;
      
      for(j=0;j<bin->nsd();++j) {
         sind = bin->sd(j);
         res[bin->b().svrtx[sind][0]][1] = 0.0;
         res[bin->b().svrtx[sind][1]][1] = 0.0;
      }
      return;
   }
};

template<class BASE, class FIXX, class FIXY> class rgeneric 
   : public BASE, public rbdry_interface {
   
   private:
      r_mesh &x;
      FIXX xfix;
      FIXY yfix;
      
   public:
      rgeneric(class r_mesh& xin, int type) : BASE(xin,type), x(xin) {}
      inline r_mesh& b() {return(x);}
      void dirichlet(FLT (*res)[mesh::ND]) {
         xfix.dirichlet(this,res);
         yfix.dirichlet(this,res);
      }
      void tadvance() {}
};

typedef rgeneric<prdx_boundary,rfixx,no_rfix> rprdx;
typedef rgeneric<prdy_boundary,no_rfix,rfixy> rprdy;
typedef rgeneric<comm_boundary,no_rfix,no_rfix> rcomm;
typedef rgeneric<side_boundary,rfixx,rfixy> rfixd;
typedef rgeneric<curv_boundary,rfixx,rfixy> rfixd_curved;
typedef rgeneric<ifce_boundary,rfixx,rfixy> rifce;
typedef rgeneric<side_boundary,rfixx,no_rfix> rsymm;

