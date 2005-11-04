/*
 *  tri_hp.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include <r_mesh.h>
#include <float.h>
#include <hpbasis.h>
#include <blocks.h>

class hp_vrtx_bdry;
class hp_side_bdry;

#ifdef AXISYMMETRIC
#define RAD(I,J) crd(0)(I,J)
#define RAD1D(I) crd(0)(0,I)
#else
#define RAD(I,J) 1
#define RAD1D(I) 1
#endif

/** This class is just the data storage and nothing for multigrid */
class tri_hp : public r_mesh  {
   protected:
      int NV;
      int p0, sm0, im0;  /**> Initialization values */
      int log2p; /**> index of basis to use in global basis::tri array */
      
      /* STATIC WORK ARRAYS */
      static Array<TinyMatrix<FLT,MXGP,MXGP>,1> u,res;
      static Array<TinyMatrix<FLT,MXGP,MXGP>,2> du;
      static TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> crd;
      static TinyMatrix<FLT,MXGP,MXGP> cjcb;
      static TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,ND,ND> dcrd;
      static Array<TinyVector<FLT,MXTM>,1> uht,lf;
      static TinyMatrix<FLT,ND,MXTM> cht, cf;
      static TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> mvel; // for local mesh velocity info
      static Array<TinyMatrix<FLT,MXGP,MXGP>,2> bdwk;

      /** Stores vertex, side and interior coefficients of solution */
      struct vsi {
         Array<FLT,2> v;
         Array<FLT,3> s;
         Array<FLT,3> i;
      } ug;
      
      /** vertex boundary information */
      Array<hp_vrtx_bdry *,1> hp_vbdry;
      hp_vrtx_bdry* getnewvrtxobject(int bnum, std::map<std::string,std::string> *bdrydata);

      /** side boundary information */
      Array<hp_side_bdry *,1> hp_sbdry;
      hp_side_bdry* getnewsideobject(int bnum, std::map<std::string,std::string> *bdrydata);
      
      /** Array for time history information */
      TinyVector<vsi,sim::nhist+1> ugbd;
      TinyVector<Array<TinyVector<FLT,ND>,1>,sim::nhist+1> vrtxbd;
      Array<TinyVector<FLT,r_mesh::ND>,1> dvrtdt; //!< Precalculated backwards difference mesh info at vertices
      Array<TinyMatrix<FLT,MXGP,MXGP>,2> dugdt; //!< Precalculated unsteady sources at Gauss points
      
#ifdef PV3
      /** Variables to understand iterative convergence using pV3 */
      struct vsi ugpv3; // STORAGE FOR pV3 to see mode change
      Array<TinyVector<FLT,r_mesh::ND>,1> vrtxpv3; // STORAGE FOR pV3 to see mode change
#endif
      
      /* Multigrid stuff needed on each mesh */
      bool isfrst; // FLAG TO SET ON FIRST ENTRY TO COARSE MESH
      Array<vsi,1> dres; //!< Driving term for multigrid
      Array<FLT,2> vug_frst; //!< Solution on first entry to coarse mesh
      
      /** Adaptation constants */
      FLT trncerr, bdrysensitivity, vlngth_tol;  //!<   Adaptation constants  
      
      /* THESE THINGS ARE SHARED BY MESHES OF THE SAME BLOCK */
      struct gbl : public r_mesh::gbl {
         
         /**< Pointer to adaptation solution storage 
          * Also used for backwards difference storage in tadvance 
          * could be used for ug0 res and res_r as well? 
          */
         tri_hp *pstr;  

         /* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
         vsi ug0;

         /** Residual storage for equations */
         vsi res; 

         /* REAL PART FOR RESIDUAL STORAGE */
         vsi res_r;  
         
         /* RESIDUAL STORAGE FOR ENTRY TO MULTIGRID */
         vsi res0;

         /* MATRIX PRECONDITIONER  */
#ifndef MATRIX_PRECONDITIONER
         Array<FLT,2> vprcn, sprcn, tprcn;
#else
         Array<FLT,3> vprcn, sprcn, tprcn;
#endif
         
         /* INITIALIZATION AND BOUNDARY CONDITION FUNCTION */
         FLT (*func)(int n, FLT x, FLT y);
      } *hp_gbl;

      /* FUNCTIONS FOR MOVING GLOBAL TO LOCAL */
      void ugtouht(int tind);
      void ugtouht(int tind,int nhist);
      void ugtouht_bdry(int tind);
      void ugtouht_bdry(int tind, int nhist);
      void ugtouht1d(int sind);
      void ugtouht1d(int sind, int nhist);
      void crdtocht(int tind);
      void crdtocht(int tind, int nhist);
      void crdtocht1d(int sind);
      void crdtocht1d(int sind, int nhist);
      void restouht_bdry(int tind); // USED IN MINVRT
      
      /* THIS FUNCTION ADDS LF TO GLOBAL VECTORS */
      void lftog(int tind, vsi gvect);

      /* SETUP V/S/T INFO */
      void setbcinfo();
      
   public:
      tri_hp() : r_mesh() {}
      virtual ~tri_hp();
      void init(std::map <std::string,std::string>& input, std::string prefix, gbl *rgin);
      void copy_data(const tri_hp &tgt);
      
      /* Initialization functions */
      inline void tobasis(FLT (*func)(int var, TinyVector<FLT,ND> &x), int tlvl = 0);
      void curvinit();
      
      /* Input / Output functions */
      enum filetype {tecplot, text, binary, adapt_diagnostic};
      void input(char *name, filetype type, int tlvl = 0);
      void output(char *name, filetype type = tecplot, int tlvl = 0);

      /* Some other handy utilities */
      void l2error(FLT (*func)(int, TinyVector<FLT,ND> &x));
      block::ctrl findmax(int excpt, int bnum, FLT (*fxy)(TinyVector<FLT,ND> &x));
      void findintercept(int bnum, FLT (*fxy)(TinyVector<FLT,ND> &x));
      
      /* Routines for point probe */
      void ptprobe(TinyVector<FLT,ND> xp, Array<FLT,1> uout, int tlvl);
      void ptprobe_bdry(int bnum, TinyVector<FLT,ND> xp, Array<FLT,1> uout, int tlvl);
      
      /* Adaptation Routines */
      void findinteriorpt(TinyVector<FLT,2> pt, int &tind, FLT &r, FLT &s);
      void findandmvptincurved(TinyVector<FLT,2> pt,int &tind, FLT &r, FLT &s);
      
      /* FUNCTIONS FOR ADAPTION */ 
      virtual block::ctrl length(int excpt) {return(block::stop);}
      block::ctrl adapt(int excpt,FLT tol);
      void movevdata(int frm, int to);
      void movevdata_bdry(int bnum,int bel,int endpt);
      void updatevdata(int v);
      void updatevdata_bdry(int bnum,int bel,int endpt);
      void movesdata(int frm, int to);
      void movesdata_bdry(int bnum,int bel);
      void updatesdata(int s);
      void updatesdata_bdry(int bnum,int bel);
      void movetdata(int frm, int to);
      void updatetdata(int t);

      /* Routines to do explicit update of solution */
      virtual block::ctrl setup_preconditioner(int excpt);
      virtual block::ctrl rsdl(int excpt) {return(block::stop);}
      block::ctrl update(int excpt);
      block::ctrl minvrt(int excpt);
      block::ctrl minvrt_test(int excpt);
      
      /* MGRID TRANSFER */
      inline void setlog2p(int value) { tri_hp::log2p = value; } /* To switch polynomial degree */
      block::ctrl mg_getfres(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh);
      block::ctrl mg_getcchng(int excpt,Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *cmesh);

      /* ADVANCE TIME SOLUTION */
      virtual block::ctrl tadvance(bool coarse,int execpoint,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh) {return(block::stop);}

      /* MESSAGE PASSING ROUTINES SPECIALIZED FOR SOLUTION CONTINUITY */
      void vmsgload(int phase, FLT *vdata);
      int vmsgwait_rcv(int phase,FLT *vdata);
      int vmsgrcv(int phase,FLT *vdata);
      void smsgload(int phase,FLT *sdata, int bgnmode, int endmode, int modestride);
      int smsgwait_rcv(int phase,FLT *sdata, int bgnmode, int endmode, int modestride);
      int smsgrcv(int phase,FLT *sdata, int bgnmode, int endmode, int modestride);

#ifdef BACKDIFF
      void unsteady_sources(int mgrid);
      void shift();
#else
      void unsteady_sources(int stage, int mgrid);
#endif             
      FLT maxres(FLT *err);
      
#ifdef PV3
      void pvstruc(int& knode, int& kequiv, int& kcel1, int& kcel2, int& kcel3, int& kcel4, int& knptet, int &kptet,int& knblock,int &blocks,int &kphedra, int& ksurf,int& knsurf,int& hint);
      void pvcell(int &kn, int &kpoffset, int cel1[][4], int cel2[][5], int cel3[][6], int cel4[][8], int nptet[][8], int ptet[]);
      void pvgrid(int &kn, float (*xyz)[3]);
      void pvsurface(int snum, int &offset, int nsurf[][3], int scon[], int scel[][4], char tsurf[][20]);
      void pvvect(int &offset, float v[][3]);
      void flotov(int &kn, struct vsi flo,int nvar, float *v);
      void meshtov(int &kn, FLT (*vin)[r_mesh::ND], struct bistruct **bin, int nvar, float *v);
      void pv3freeze();
      void pv3subtract(int frozen);
#endif

};


