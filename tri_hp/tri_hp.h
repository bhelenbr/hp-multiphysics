/*
 *  tri_hp.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_h_
#define _tri_hp_h_

#include <r_mesh.h>
#include <float.h>
#include <hpbasis.h>
#include <blocks.h>

#ifdef AXISYMMETRIC
#define RAD(r) (r)
#else
#define RAD(r) 1
#endif

class hp_vrtx_bdry;
class hp_side_bdry;

class init_bdry_cndtn {
   public:
      virtual FLT f(int n, TinyVector<FLT,mesh::ND> x) = 0;
      virtual void input(input_map &blkdata, std::string idnty) {};
};

class mesh_mover;

/** This class is just the data storage and nothing for multigrid */
class tri_hp : public r_mesh  {
   public:
      int NV;
      int p0, sm0, im0;  /**> Initialization values */
      int log2p; /**> index of basis to use in global basis::tri array */
      int log2pmax; /**> Initialization value of log2p */
      bool coarse; /**> tells whether we are coarse level of multigrid or not */
      enum movementtype {fixed,uncoupled_rigid,coupled_rigid,uncoupled_deformable,coupled_deformable} mmovement;
      bool adapt_flag;
      
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
      virtual hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
      /** side boundary information */
      Array<hp_side_bdry *,1> hp_sbdry;
      virtual hp_side_bdry* getnewsideobject(int bnum, input_map &bdrydata); 
      /** object to perform rigid mesh movement */
      mesh_mover *mover;
      
      /** Array for time history information */
      TinyVector<vsi,sim::nhist+1> ugbd;
      TinyVector<Array<TinyVector<FLT,ND>,1>,sim::nhist+1> vrtxbd; //!< Highest level contains pre-summed unsteady mesh velocity source
      Array<TinyMatrix<FLT,MXGP,MXGP>,3> dugdt; //!< Precalculated unsteady sources at Gauss points
      Array<TinyMatrix<FLT,MXGP,MXGP>,3> dxdt; //!< Precalculated mesh velocity sources at Gauss points
      
      /* Multigrid stuff needed on each mesh */
      bool isfrst; // FLAG TO SET ON FIRST ENTRY TO COARSE MESH
      Array<vsi,1> dres; //!< Driving term for multigrid
      Array<FLT,2> vug_frst; //!< Solution on first entry to coarse mesh
      FLT fadd; //!< Controls addition of residuals on coarse mesh
      
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

         /* PRECONDITIONER  */
         bool diagonal_preconditioner;
         Array<FLT,2> vprcn, sprcn, tprcn;  // Diagonal preconditioner
         Array<FLT,3> vprcn_ut, sprcn_ut, tprcn_ut; // Lower triangle preconditioner
         
         /* INITIALIZATION AND BOUNDARY CONDITION FUNCTION */
         init_bdry_cndtn *ibc;
                  
         /* Pointers to block storage objects for side boundary conditions */
         Array<void *,1> sbdry_gbls;
         
         /* Pointers to block storage objects for vrtx boundary conditions */
         Array<void *,1> vbdry_gbls;
         
         /* Time step factor for different polynomial degree */
         TinyVector<FLT,MXGP> cfl;
         
      } *gbl_ptr;
      virtual init_bdry_cndtn* getnewibc(input_map& inmap);
      virtual mesh_mover* getnewmesh_mover(input_map& inmap);

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
      void setinfo();
      
   public:
      tri_hp() : r_mesh() {}
      virtual ~tri_hp();
      void init(input_map& input, gbl *rgin);
      void copy_data(const tri_hp &tgt);
      virtual tri_hp* create() = 0;
      
      /* Initialization functions */
      void tobasis(init_bdry_cndtn *ibc, int tlvl = 0);
      void curvinit();
      
      /* Input / Output functions */
      enum filetype {tecplot, text, binary, adapt_diagnostic};
      TinyVector<filetype,3> output_type;
      void input(const std::string &name);
      void input(const std::string &name, filetype type, int tlvl = 0);
      void output(const std::string &name, block::output_purpose why);
      void output(const std::string &name, filetype type = tecplot, int tlvl = 0);

      /* Some other handy utilities */
      void l2error(init_bdry_cndtn *toCompare);
      block::ctrl findmax(block::ctrl ctrl_message, int bnum, FLT (*fxy)(TinyVector<FLT,ND> &x));
      void findintercept(int bnum, FLT (*fxy)(TinyVector<FLT,ND> &x));
      void integrated_averages(Array<FLT,1> a);
      
      /* Routines for point probe */
      void ptprobe(TinyVector<FLT,ND> xp, Array<FLT,1> uout, int tlvl);
      void ptprobe_bdry(int bnum, TinyVector<FLT,ND> xp, Array<FLT,1> uout, int tlvl);
      
      /* Adaptation Routines */
      void findinteriorpt(TinyVector<FLT,2> pt, int &tind, FLT &r, FLT &s);
      void findandmvptincurved(TinyVector<FLT,2>& pt,int &tind, FLT &r, FLT &s);
      
      /* FUNCTIONS FOR ADAPTION */ 
      virtual block::ctrl length(block::ctrl ctrl_message) {*sim::log << "here\n"; return(block::stop);}
      block::ctrl adapt(block::ctrl ctrl_message,FLT tol);
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
      virtual block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      virtual block::ctrl rsdl(block::ctrl ctrl_message, int stage = sim::NSTAGE) {   
         /* ONLY NEED TO CALL FOR MOVEMENT BETWEEN MESHES */
         if (mmovement == coupled_deformable && stage == sim::NSTAGE && log2p == 0) {
            return(r_mesh::rsdl(ctrl_message));            
         }
         return(block::stop);
      }
      FLT maxres();      
      block::ctrl update(block::ctrl ctrl_message);
      block::ctrl minvrt(block::ctrl ctrl_message);
      block::ctrl minvrt_test(block::ctrl ctrl_message);
      
      /* MGRID TRANSFER */
      inline void setlog2p(int value) { tri_hp::log2p = value; } /* To switch polynomial degree */
      block::ctrl mg_getfres(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh);
      block::ctrl mg_getcchng(block::ctrl ctrl_message,Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *cmesh);

      /* ADVANCE TIME SOLUTION */
      block::ctrl tadvance(bool coarse,block::ctrl ctrl_message,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh);
      virtual void calculate_unsteady_sources(bool coarse);
      
      /* MESSAGE PASSING ROUTINES SPECIALIZED FOR SOLUTION CONTINUITY */
      void vc0load(int phase, FLT *vdata, int vrtstride=1);
      int vc0wait_rcv(int phase,FLT *vdata, int vrtsride=1);
      int vc0rcv(int phase,FLT *vdata, int vrtstride=1);
      void sc0load(FLT *sdata, int bgnmode, int endmode, int modestride);
      int sc0wait_rcv(FLT *sdata, int bgnmode, int endmode, int modestride);
      int sc0rcv(FLT *sdata, int bgnmode, int endmode, int modestride);
      block::ctrl matchboundaries(block::ctrl ctrl_message);
      
   private:
      int excpt,excpt1,stage,mp_phase,mode;
};

class mesh_mover {
   public:
      mesh_mover(tri_hp& xin) {}
      virtual void init(input_map& input, std::string idnty) {}
      virtual mesh_mover* create(tri_hp& xin) { return new mesh_mover(xin); }
      virtual block::ctrl tadvance(block::ctrl ctrl_message) {return(block::stop);}
      virtual block::ctrl setup_preconditioner(block::ctrl ctrl_message) {return(block::stop);}
      virtual block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE) {return(block::stop);}
      virtual block::ctrl update(block::ctrl ctrl_message) {return(block::stop);}
};


#endif
