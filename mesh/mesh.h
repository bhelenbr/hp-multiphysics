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


class side_bdry;
class vrtx_bdry;

class mesh {

   /***************/
   /* DATA        */
   /***************/
   public:
      int initialized, maxvst;
      static const int ND = 2;
      static const int DIM = ND;
      std::ostream *log;
      
      /* VERTEX DATA */
      int nvrtx;
      FLT (*vrtx)[ND];
      FLT *vlngth;
      struct vstruct {
         int tri;
         int nnbor;
         int info;
      } *vd;
      class quadtree<ND> qtree;
   
      /* VERTEX BOUNDARY INFO */
      static const int MAXVB = 8;
      int nvbd;
      vrtx_bdry *vbdry[MAXVB];
      vrtx_bdry* getnewvrtxobject(int idnum, std::map<std::string,std::string> *bdrydata);

      /* SIDE DATA */      
      int nside;
      struct sstruct {
         int vrtx[2];
         int tri[2];
         int info;
      } *sd;
      
      /* SIDE BOUNDARY INFO */
      static const int MAXSB = 8;
      int nsbd;
      side_bdry *sbdry[MAXSB]; 
      side_bdry* getnewsideobject(int idnum, std::map<std::string,std::string> *bdrydata);

      /* TRIANGLE DATA */      
      int ntri;
      struct tstruct {
         int vrtx[3];
         int side[3];
         int sign[3];
         int tri[3];
         int info;
      } *td;
      
      /* THESE ARE NOT SHARED BECAUSE THEY HAVE TO BE KEPT INITIALIZED TO -1 */
      static int maxsrch, maxlst;
      static int *i1wk, *i2wk, *i3wk;
      /* THIS SHARED FLOATING WORK */
      sharedmem *scratch;
      FLT *fwk;
      
      /**************/
      /*  INTERFACE */
      /**************/
      /* DEFAULT INITIALIZATION */
      mesh() : initialized(0), log(&std::cout) {}
      sharedmem* allocate(int mxsize, sharedmem *wkin = 0);
      void allocate_duplicate(FLT sizereduce1d,const class mesh& xmesh);
      void load_scratch_pointers();
      void copy(const mesh& tgt);
      void coarsen_substructured(const class mesh &tgt,int p);
      void shift(FLT s[ND]);
      void scale(FLT s[ND]);

      /* INPUT/OUTPUT MESH (MAY MODIFY VINFO/SINFO/TINFO) */
      sharedmem* input(const char *filename, ftype::name filetype = ftype::easymesh,  FLT grwfac = 1, const char *bdrymap = 0,sharedmem *win = 0);
      int output(const char *filename, ftype::name filetype = ftype::easymesh) const;
      void setbcinfo();

      /* MESH MODIFICATION */  
      /* TO SET UP ADAPTATION VLENGTH */
      void initvlngth();
      void length1();
      void length2(); 
      int coarsen(FLT factor, const class mesh& xmesh);
      void coarsen2(FLT factor, const class mesh& inmesh, FLT size_reduce = 1.0);
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
      void partition(const class mesh& xmesh, int npart);
      int comm_entity_size();
      int comm_entity_list(int *list);
      int msgpass(int phase);
      void matchboundaries1();
      void matchboundaries2();
      
      /* THIS IS AN INTERPOLATION DATA STRUCTURE */
      struct transfer {
         int tri;
         FLT wt[3];
      };
      block::ctrl mgconnect(int excpt, transfer *cnnct, const class mesh& tgt);
      void testconnect(char *fname,transfer *cnnct, mesh *cmesh);
      
      /* SOME DEGUGGING FUNCTIONS */
      void checkintegrity() const;
      void checkintwk() const;
      
      /* NEW ADAPTATION ROUTINES */
      void setup_for_adapt();
      void coarsen_bdry(FLT tolsize);
      FLT side_lngth_ratio(int i) {   
         /* CALCULATE SIDE TO TARGET LENGTH RATIO */
         /* ERROR ON THE CONSERVATIVE SIDE (OVER-REFINED) BY TAKING MAX VLNGTH */
         return(MAX(vlngth[sd[i].vrtx[0]],vlngth[sd[i].vrtx[1]])/distance(sd[i].vrtx[0],sd[i].vrtx[1]));
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

      /* MESH MODIFICATION FUNCTIONS */
      /* TO CREATE AN INITIAL TRIANGUlATION */
      void triangulate(int nside);
      void addtri(int v0,int v1,int v2,int sind,int dir);

      /* TO INSERT A POINT */
      int insert(FLT x[ND]);
      int insert(int tind, int vnum, FLT theta, int &ntdel, int *tdel, int &nsdel, int *sdel);
      void bdry_insert(int tind, int snum, int vnum, int &ntdel, int *tdel, int &nsdel, int *sdel);
      int findtri(FLT x[ND], int vnear) const;
      void fltwkreb(int sind);  //FUNCTIONS FOR SETTUPING UP FLTWK FOR REBAY
      void fltwkreb();
      
      /* FOR COARSENING A SIDE */
      int collapse(int sind,int& ntdel,int *tdel,int& nsdel,int *nsdel);
      void dltvrtx(int vind);
      void dltd(int sind);
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
            d += pow(vrtx[v0][n] -vrtx[v1][n],2);
         return(sqrt(d));
      }
      inline FLT distance2(int v0, int v1) const {
         FLT d = 0.0;
         for(int n = 0; n<ND;++n)
            d += pow(vrtx[v0][n] -vrtx[v1][n],2);
         return(d);
      }  
      FLT incircle(int tind, FLT *a) const;
      FLT insidecircle(int sind, FLT *a) const;
      FLT area(int sind, int vind) const;
      FLT area(int v0, int v1, int v2) const;
      FLT area(int tind) const;
      FLT minangle(int v0, int v1, int v2) const;
      FLT angle(int v0, int v1, int v2) const;
      FLT tradius(int tind) const;
      void tcenter(int tind, FLT x[ND]) const;
      FLT aspect(int tind) const;
      FLT intri(int tind, FLT x[ND]) const;
      void getwgts(FLT *) const;
};
#endif
