#ifndef _mesh_h_
#define _mesh_h_

#include <math.h>
#include <quadtree.h>
#include <ftype.h>
#include <iostream>
#include <float.h>
#include <utilities.h>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <blitz/array.h>
#include "block.h"

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

#define NO_DEBUG_ADAPT

using namespace blitz;

class side_bdry;
class vrtx_bdry;

class mesh {

   /***************/
   /* DATA        */
   /***************/
   public:
      int maxvst;
      static const int ND = 2;
            
      /* VERTEX DATA */
      int nvrtx;
      Array<TinyVector<FLT,ND>,1> vrtx;
      Array<FLT,1> vlngth;
      struct vstruct {
         int tri;
         int nnbor;
         int info;
      };
      Array<vstruct,1> vd;
      class quadtree<ND> qtree;
   
      /* VERTEX BOUNDARY INFO */
      int nvbd;
      Array<vrtx_bdry *,1> vbdry;
      vrtx_bdry* getnewvrtxobject(int idnum, std::map<std::string,std::string> *bdrydata);

      /* SIDE DATA */      
      int nside;
      struct sstruct {
         TinyVector<int,2> vrtx;
         TinyVector<int,2> tri;
         int info;
      };
      Array<sstruct,1> sd;
      
      /* SIDE BOUNDARY INFO */
      int nsbd;
      Array<side_bdry *,1> sbdry;
      side_bdry* getnewsideobject(int idnum, std::map<std::string,std::string> *bdrydata);

      /* TRIANGLE DATA */      
      int ntri;
      struct tstruct {
         TinyVector<int,3> vrtx;
         TinyVector<int,3> side;
         TinyVector<int,3> sign;
         TinyVector<int,3> tri;
         int info;
      };
      Array<tstruct,1> td;
      
      /* SOME SCRATCH VARIABLES */
      /* ANY ROUTINE THAT USES i1wk SHOULD RESET IT TO -1 */
      static Array<int,1> i1wk,i2wk, i2wk_lst1, i2wk_lst2, i2wk_lst3;

      /* THIS IS SHARED WORK SPACE */
      /* ONLY SHARE WITH BLOCKS THAT WORK TOTALLY INDEPENDENTLY */
      /* (DIFFERENT LEVELS OF MULTIGRID) */
      sharedmem scratch;
      Array<FLT,1> fscr1;
      
      int initialized;
      std::ostream *log;
      static int maxsrch;
      
   public:
      /**************/
      /*  INTERFACE */
      /**************/
      /* INITIALIZATION & ALLOCATION */
      mesh() : nvbd(0), nsbd(0), initialized(0), log(&std::cout)  {}
      sharedmem* allocate(int mxsize, const sharedmem *wkin = 0);
      void allocate_duplicate(FLT sizereduce1d,const class mesh& xmesh);
      void get_scratch_pointers() {
         Array<FLT,1> tmp(static_cast<FLT *>(scratch.data()), maxvst, neverDeleteData);
         fscr1.reference(tmp);
      }
      void copy(const mesh& tgt);
      ~mesh();
      
      /* INPUT/OUTPUT MESH (MAY MODIFY VINFO/SINFO/TINFO) */
      sharedmem* input(const char *filename, ftype::name filetype = ftype::easymesh,  FLT grwfac = 1, const char *bdrymap = 0,sharedmem *win = 0);
      int output(const char *filename, ftype::name filetype = ftype::easymesh) const;
      void bdry_output(const char *filename) const;
      void setbcinfo();  // FOR EASYMESH OUTPUT (NOT USED)

      /* MESH MODIFICATION UTILTIES */
      void coarsen_substructured(const class mesh &tgt,int p);
	   void symmetrize();
      void append(const mesh &z);
      void shift(TinyVector<FLT,ND>& s);
      void scale(TinyVector<FLT,ND>& s);
      int smooth_cofa(int niter);
      void refineby2(const class mesh& xmesh);
      void settrim();
      void initvlngth();
      block::ctrl adapt(int excpt, FLT tolsize);
      int coarsen(FLT factor, const class mesh& xmesh);
      void coarsen2(FLT factor, const class mesh& inmesh, FLT size_reduce = 1.0);
      void coarsen3();

      /* UTILITIES FOR PARALLEL COMPUTATIONS */
#ifdef METIS
      void setpartition(int nparts); 
#endif 
      void partition(class mesh& xmesh, int npart);
      int comm_entity_size();
      int comm_entity_list(Array<int,1>& list);
      void msgload(int phase,FLT *base,int bgn, int end, int stride);
      void msgpass(int phase);
      int msgwait_rcv(int phase,FLT *base,int bgn, int end, int stride);
      int msgrcv(int phase,FLT *base,int bgn, int end, int stride);
      void matchboundaries1(int phase);
      int matchboundaries2(int phase);
      
      
      /* UTILITIES FOR INTERPOLATION BETWEEN MESHES */
      struct transfer {
         int tri;
         TinyVector<FLT,3> wt;
      };
      block::ctrl mgconnect(int excpt, Array<transfer,1> &cnnct, const class mesh& tgt);
      void testconnect(char *fname,Array<transfer,1> &cnnct, mesh *cmesh);
      
      /* SOME DEGUGGING FUNCTIONS */
      void checkintegrity();
      void checkintwk() const;
      
   protected:
   
      /*******************/
      /* INTERNAL FUNCTIONS */
      /*******************/  
      void cnt_nbor(void);
      void bdrylabel(void);
      void createsideinfo(void);
      void createtdstri(void);
      void createttri(void);
      void createvtri(void);
      void treeinit();
      void treeinit(FLT x1[ND], FLT x2[ND]);

      /* MESH MODIFICATION FUNCTIONS */
      /* 4 binary digits for vertex, side, tri */
      const static int VSPEC = 0x4, VDLTE = 0x2, VTOUC=0x1;
      const static int SDLTE = 0x10*0x2, STOUC=0x10*0x1;
      const static int TSRCH = 0x100*0x4, TDLTE = 0x100*0x2, TTOUC=0x100*0x1;
      void setup_for_adapt();
      void swap(FLT swaptol = EPSILON);
      void bdry_yaber(FLT tolsize);
      void bdry_yaber1();
      void yaber(FLT tolsize);
      void bdry_rebay(FLT tolsize);
      void bdry_rebay1();
      void rebay(FLT tolsize);
      void cleanup_after_adapt();

      FLT stgt_to_actual(int i) {   
         /* CALCULATE SIDE TO TARGET LENGTH RATIO */
         /* ERROR ON THE CONSERVATIVE SIDE (OVER-REFINED) BY TAKING MAX VLNGTH */
         return(MAX(vlngth(sd(i).vrtx(0)),vlngth(sd(i).vrtx(1)))/distance(sd(i).vrtx(0),sd(i).vrtx(1)));
      }      
      /* TO CREATE AN INITIAL TRIANGUlATION */
      void triangulate(int nside);
      void addtri(int v0,int v1,int v2,int sind,int dir);

      /* TO INSERT A POINT */
      int insert(const TinyVector<FLT,2> &x);
      int insert(int vnum, int tnum);
      void bdry_insert(int vnum, int sind, int endpt = 0);
      int findtri(TinyVector<FLT,ND> x, int vnear) const;

      /* FOR COARSENING A SIDE */
      void collapse(int sind, int endpt);

      void dltvrtx(int vind);
      void dltsd(int sind);
      void dlttri(int tind);
      
      /* ORDERED LIST FUNCTIONS */
      void putinlst(int sind);
      void tkoutlst(int sind);
      
      /* TO EDGE SWAP FOR IMPROVED QUALITY */
      int swap(int sind, FLT tol = 0.0);  //SWAPS A SINGLE SIDE 
      void swap(int nswp, int *swp, FLT tol = 0.0); //SWAPS A LIST OF SIDES

      /* PRIMITIVE FUNCTIONS */
   public:
      inline FLT distance(int v0, int v1) const {
         FLT d = 0.0;
         for(int n = 0; n<ND;++n)
            d += pow(vrtx(v0)(n) -vrtx(v1)(n),2);
         return(sqrt(d));
      }
      inline FLT distance2(int v0, int v1) const {
         FLT d = 0.0;
         for(int n = 0; n<ND;++n)
            d += pow(vrtx(v0)(n) -vrtx(v1)(n),2);
         return(d);
      }  
      FLT incircle(int tind, const TinyVector<FLT,2> &x) const;
      FLT insidecircle(int sind, const TinyVector<FLT,2> &x) const;
      FLT area(int sind, int vind) const;
      FLT area(int v0, int v1, int v2) const;
      FLT area(int tind) const;
      FLT minangle(int v0, int v1, int v2) const;
      FLT angle(int v0, int v1, int v2) const;
      FLT circumradius(int tind) const;
      void circumcenter(int tind, TinyVector<FLT,2> &x) const;
      FLT inscribedradius(int tind) const;
      FLT aspect(int tind) const;
      FLT intri(int tind, const TinyVector<FLT,2> &x) const;
      void getwgts(TinyVector<FLT,3> &wgt) const;
      /* tri numbers at boundary point to el and group */
      int getbdrynum(int trinum) const { return((-trinum>>16) -1);}
      int getbdryel(int trinum) const { return(-trinum&0xFFFF);}
      int trinumatbdry(int bnum, int bel) const { return(-(((bnum+1)<<16) +bel));}
};
#endif
