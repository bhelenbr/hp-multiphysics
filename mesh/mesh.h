#ifndef _mesh_h_
#define _mesh_h_

#include <math.h>
#include <quadtree.h>
#include <ftype.h>
#include <iostream>

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#define FLT double
#define EPSILON DBL_EPSILON
#endif

template<class T> class side_template;
template<class T> class vrtx_template;

template<int ND> class mesh {

   /***************/
   /* DATA        */
   /***************/
   protected:
      int initialized, maxvst;

   public:
      std::ostream *log;
      static const int DIM = ND;
#ifdef CAPRI
      int cpri_vol, cpri_face;
#endif
      
      /* STORES THINGS NEEDED FOR A BLOCK OF MESHES */
      /* AT THIS POINT DON'T NEED ANYTHING */
      struct gbl {
      } *rg;
      
      /* VERTEX DATA */
      int nvrtx;
      FLT (*vrtx)[ND];
      FLT *vlngth;
      int *vinfo;
      int *nnbor;
      int *vtri;
      class quadtree<ND> qtree;
   
      /* VERTEX BOUNDARY INFO */
      static const int MAXVB = 8;
      int nvbd;
      vrtx_template<mesh<ND> > *vbdry[MAXVB];
      virtual void getnewvrtxobject(int bnum, int type);
     
      /* SIDE DATA */      
      int nside;
      int (*svrtx)[2];
      int (*stri)[2];
      int *sinfo;
      
      /* SIDE BOUNDARY INFO */
      static const int MAXSB = 8;
      int nsbd;
      side_template<mesh<ND> > *sbdry[MAXSB]; 
      virtual void getnewsideobject(int bnum, int type);

      /* TRIANGLE DATA */      
      int ntri;
      int (*tvrtx)[3];
      int (*ttri)[3];
      struct tsidedata {
         int side[3];
         int sign[3];
      } *tside;
      int *tinfo;

      /* THIS IS AN INTERPOLATION STRUCTURE IF USING MULTIGRID */
      /* UNALLOCATED UNLESS CONNECT IS CALLED */      
      struct mg_trans {
         int tri;
         double wt[3];
      };
      void *fmpt;  // CAN BE USED TO POINT TO ANY OBJECT INHERITED FROM MESH
      void *cmpt;  // DITTO
      struct mg_trans *fine;
      struct mg_trans *coarse;
      void mgconnect(struct mg_trans *cnnct, const class mesh& tgt);
      
      /* THESE ARE STATIC WORK ARRAYS FOR MESH REFINEMENT ROUTINES */
      static FLT *fltwk;
      static int gblmaxvst, *intwk1, *intwk2, *intwk3;
      static int maxsrch;  // MAX TRIS TO SEARCH IN FINDTRI ROUTINE
      
      /**************/
      /*  INTERFACE */
      /**************/
      /* DEFAULT INITIALIZATION */
      mesh() : initialized(0), log(&std::cout), fmpt(NULL), cmpt(NULL), fine(NULL), coarse(NULL) {}
      void allocate(bool coarse, mesh::gbl *ginit) {}
      void copy(const mesh& tgt);

      /* INPUT/OUTPUT MESH (MAY MODIFY VINFO/SINFO/TINFO) */
      int in_mesh(FLT (*vin)[ND], const char *filename, FTYPE filetype = easymesh, FLT grwfac = 1);
      inline int in_mesh(const char *filename, FTYPE filetype = easymesh, FLT grwfac = 1) {return(in_mesh(vrtx,filename,filetype,grwfac));}
      int out_mesh(FLT (*vin)[ND], const char *filename, FTYPE filetype = easymesh) const;
      inline int out_mesh(const char *filename, FTYPE filetype = easymesh) const {return(out_mesh(vrtx,filename,filetype));}
      void setbcinfo();
      
      /* ACCESS TO SIMPLE MESH DATA */
      int max() {return maxvst;}

      /* MESH MODIFICATION */  
      /* TO SET UP ADAPTATION VLENGTH */
      void initvlngth();
      virtual void setlength() {
        // for(int i=0;i<nvrtx;++i)
        //    vlngth[i] = 0.0982;
         *log << "in generic setlength" << std::endl;
      }
      void length1();
      void length2(); 
      int coarsen(FLT factor, const class mesh& xmesh);
      void coarsen2(FLT factor, const class mesh& inmesh, class mesh& work);
      void refineby2(const class mesh& xmesh);
      void swap(FLT swaptol = 0.0);
      void yaber(FLT tolsize, int yes_swap, FLT swaptol = 0.0);
      inline void treeupdate() { qtree.update(0,nvrtx);}
      void rebay(FLT tolsize);
      inline void adapt(FLT tolerance) {
         yaber(1.0/tolerance,1,0.0);
         out_mesh("coarse",tecplot);
         treeupdate();
         rebay(tolerance);
         out_mesh("refine",tecplot);
         return;
      }
      
      /* CENTER OF AREA MESH SMOOTHING */
      int smooth_cofa(int niter);

      /* FIND MATCHING MESH BOUNDARIES */
      int comm_entity_size();
      int comm_entity_list(int *list);
      int msgpass(int phase);
      void matchboundaries1();
      void matchboundaries2();

      /* FUNCTION TO ALLOCATE & SET INTERPOLATION WEIGHTS BETWEEN TWO MESHES */
      void setfine(class mesh& tgt);
      void setcoarse(class mesh& tgt);
      void testconnect(char *fname);

      /* SOME DEGUGGING FUNCTIONS */
      void checkintegrity() const;
      void checkintwk() const;
      
   protected:
      /*******************/
      /* INTERNAL FUNCTIONS */
      /*******************/  
      void allocate(int mxsize);
       
      void cnt_nbor(void);
      void bdrylabel(void);
      void createsideinfo(void);
      void createtsidestri(void);
      void createttri(void);
      void createvtri(void);
      void treeinit();

      /* MESH MODIFICATION FUNCTIONS */
      /* TO CREATE AN INITIAL TRIANGUlATION */
      void triangulate(int **sidelst, int *nside, int nloop, int cknbor = 1);
      void findpt(int *vrtxlst,int nv,int *v,int chkadj,int good[], int &ngood);
      void addtri(int v0,int v1,int v2,int sind,int dir);

      /* TO INSERT A POINT */
      int insert(FLT x[ND]);
      int insert(int tind, int vnum, FLT theta = 0.0);
      void bdry_insert(int tind, int snum, int vnum);
      int findtri(FLT x[ND], int vnear) const;
      void fltwkreb(int sind);  //FUNCTIONS FOR SETTUPING UP FLTWK FOR REBAY
      void fltwkreb();
      
      /* FOR COARSENING A SIDE */
      int collapse(int sind);
      void dltvrtx(int vind);
      void dltside(int sind);
      void dlttri(int tind);
      void fltwkyab(int sind);  //FUNCTIONS FOR SETTUPING UP FLTWK FOR YABER
      void fltwkyab();
      
      /* ORDERED LIST FUNCTIONS */
      static const int MAXLST = 1000;
      static int nslst;
      static int ntdel, tdel[1000+1];
      static int nsdel, sdel[1000+1];
      void putinlst(int sind);
      void tkoutlst(int sind);
      
      /* TO EDGE SWAP FOR IMPROVED QUALITY */
      int swap(int sind, FLT tol = 0.0);  //SWAPS A SINGLE SIDE 
      void swap(int nswp, int *swp, FLT tol = 0.0); //SWAPS A LIST OF SIDES

      /* PRIMITIVE FUNCTIONS */
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

#include "boundary.h"

/* REMINDER FIND AND REPLACE TO INSERT FUNCTION INSTANTIATIONS 
template\<int ND\> ([a-z]*) mesh\<ND\>::([^?]*) \{

template \1 mesh\<2\>::\2;
template \1 mesh\<3\>::\2;

template\<int ND\> \1 mesh\<ND\>::\2 \{
*/
#endif
