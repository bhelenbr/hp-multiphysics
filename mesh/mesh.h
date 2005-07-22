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

using namespace blitz;

class side_bdry;
class vrtx_bdry;

class mesh {

   /***************/
   /* DATA        */
   /***************/
   public:
      int initialized, maxvst;
      static const int ND = 2;
      std::ostream *log;
      
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
      
      /* THESE ARE NOT SHARED BECAUSE THEY HAVE TO BE KEPT INITIALIZED TO -1 */
      static int maxsrch, maxlst;
      static Array<int,1> i1wk,i2wk,i3wk;
      /* THIS IS SHARED WORK SPACE */
      /* ONLY SHARE WITH BLOCKS THAT WORK TOTALLY INDEPENDENTLY */
      /* (DIFFERENT LEVELS OF MULTIGRID) */
      sharedmem scratch;
      Array<FLT,1> fscr1;
      void get_scratch_pointers() {
         Array<FLT,1> tmp(static_cast<FLT *>(scratch.data()), maxvst, neverDeleteData);
         fscr1.reference(tmp);
      }

      /**************/
      /*  INTERFACE */
      /**************/
      /* DEFAULT INITIALIZATION */
      mesh() : initialized(0), log(&std::cout), nvbd(0), nsbd(0) {}
      sharedmem* allocate(int mxsize, const sharedmem *wkin = 0);
      void allocate_duplicate(FLT sizereduce1d,const class mesh& xmesh);
      void load_scratch_pointers();
      void copy(const mesh& tgt);
      ~mesh();
      void coarsen_substructured(const class mesh &tgt,int p);
	   void symmetrize();
      void append(const mesh &z);
      void shift(TinyVector<FLT,ND>& s);
      void scale(TinyVector<FLT,ND>& s);

      /* INPUT/OUTPUT MESH (MAY MODIFY VINFO/SINFO/TINFO) */
      sharedmem* input(const char *filename, ftype::name filetype = ftype::easymesh,  FLT grwfac = 1, const char *bdrymap = 0,sharedmem *win = 0);
      int output(const char *filename, ftype::name filetype = ftype::easymesh) const;
      void bdry_output(const char *filename) const;

      void setbcinfo();

      /* MESH MODIFICATION */  
      /* TO SET UP ADAPTATION VLENGTH */
      void initvlngth();
      void length1(int phase) {msgload(phase,vlngth.data(),0,0,1);}
      int length2(int phase) {return(msgwait_rcv(phase,vlngth.data(),0,0,1));} 
      int coarsen(FLT factor, const class mesh& xmesh);
      void coarsen2(FLT factor, const class mesh& inmesh, FLT size_reduce = 1.0);
      void coarsen3();
      void refineby2(const class mesh& xmesh);
      void swap(FLT swaptol = EPSILON);
      void yaber(FLT tolsize, int yes_swap, FLT swaptol = EPSILON);
      inline void treeupdate() { qtree.update(0,nvrtx);}
      void rebay(FLT tolsize);
      void trebay(FLT tolsize);
      inline void adapt(FLT tolerance) {
         treeinit(); // TEMPORARY
         yaber(1.0/tolerance,1,EPSILON);
         treeupdate();
         rebay(tolerance);
         return;
      }
      
      /* CENTER OF AREA MESH SMOOTHING */
      int smooth_cofa(int niter);

      /* FIND MATCHING MESH BOUNDARIES */
      void settrim();
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
      
      /* THIS IS AN INTERPOLATION DATA STRUCTURE */
      struct transfer {
         int tri;
         TinyVector<FLT,3> wt;
      };
      block::ctrl mgconnect(int excpt, Array<transfer,1> &cnnct, const class mesh& tgt);
      void testconnect(char *fname,Array<transfer,1> &cnnct, mesh *cmesh);
      
      /* SOME DEGUGGING FUNCTIONS */
      void checkintegrity();
      void checkintwk() const;
      
      /* NEW ADAPTATION ROUTINES */
      void setup_for_adapt();
      void coarsen_bdry(FLT tolsize);
      FLT side_lngth_ratio(int i) {   
         /* CALCULATE SIDE TO TARGET LENGTH RATIO */
         /* ERROR ON THE CONSERVATIVE SIDE (OVER-REFINED) BY TAKING MAX VLNGTH */
         return(MAX(vlngth(sd(i).vrtx(0)),vlngth(sd(i).vrtx(1)))/distance(sd(i).vrtx(0),sd(i).vrtx(1)));
      }
      
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
      /* TO CREATE AN INITIAL TRIANGUlATION */
      void triangulate(int nside);
      void addtri(int v0,int v1,int v2,int sind,int dir);

      /* TO INSERT A POINT */
      int insert(const TinyVector<FLT,2> &x);
      int insert(int tind, int vnum, FLT theta, int &ntdel, int *tdel, int &nsdel, int *sdel);
      void bdry_insert(int tind, int snum, int vnum, int &ntdel, int *tdel, int &nsdel, int *sdel);
      int findtri(TinyVector<FLT,ND> x, int vnear) const;
      void fltwkreb(int sind);  //FUNCTIONS FOR SETTUPING UP FLTWK FOR REBAY
      void fltwkreb();
      
      /* FOR COARSENING A SIDE */
      int collapse(int sind,int& ntdel,int *tdel,int& nsdel,int *sdel);
      int collapse1(int sind,int bgnend, int& ntdel,int *tdel,int& nsdel,int *sdel);

      void dltvrtx(int vind);
      void dltsd(int sind);
      void dlttri(int tind);
      void fltwkyab(int sind);  //FUNCTIONS FOR SETTUPING UP FLTWK FOR YABER
      void fltwkyab();
      
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
      FLT tradius(int tind) const;
      void tcenter(int tind, TinyVector<FLT,2> &x) const;
      FLT aspect(int tind) const;
      FLT intri(int tind, const TinyVector<FLT,2> &x) const;
      void getwgts(TinyVector<FLT,3> &wgt) const;
      /* tri numbers at boundary point to el and group */
      int getbdrynum(int trinum) const { return((-trinum>>16) -1);}
      int getbdryel(int trinum) const { return(-trinum&0xFFFF);}
      int trinumatbdry(int bnum, int bel) { return(-(((bnum+1)<<16) +bel));}
};
#endif
