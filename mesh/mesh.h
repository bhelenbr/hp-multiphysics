#ifndef _mesh_h_
#define _mesh_h_

#include <math.h>
#include <quadtree.h>
#include <cstdio>


#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#define FLT double
#define EPSILON DBL_EPSILON
#endif


enum FILETYPE {easymesh, gambit, tecplot, grid, text, binary, mavriplis};

class side_boundary;
class vrtx_boundary;

class mesh {

   
   /***************/
   /* DATA        */
   /***************/
   protected:
      int initialized, maxvst;

   public:
      static const int ND = 2;
      
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
      vrtx_boundary *vbdry[MAXVB];
      virtual void getnewvrtxobject(int bnum, int type);
     
      /* SIDE DATA */      
      int nside;
      int (*svrtx)[2];
      int (*stri)[2];
      int *sinfo;
      
      /* SIDE BOUNDARY INFO */
      static const int MAXSB = 8;
      int nsbd;
      side_boundary *sbdry[MAXSB]; 
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
      mesh() : initialized(0), fmpt(NULL), cmpt(NULL), fine(NULL), coarse(NULL) {}
      void copy(const mesh& tgt);

      /* INPUT/OUTPUT MESH (MAY MODIFY VINFO/SINFO/TINFO) */
      int in_mesh(FLT (*vin)[ND], const char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1);
      inline int in_mesh(const char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1) {return(in_mesh(vrtx,filename,filetype,grwfac));}
      int out_mesh(FLT (*vin)[ND], const char *filename, FILETYPE filetype = easymesh) const;
      inline int out_mesh(const char *filename, FILETYPE filetype = easymesh) const {return(out_mesh(vrtx,filename,filetype));}
      void setbcinfo();
      
      /* ACCESS TO SIMPLE MESH DATA */
      int max() {return maxvst;}

      /* MESH MODIFICATION */  
      /* TO SET UP ADAPTATION VLENGTH */
      void initvlngth();
      virtual void setlength() {printf("here\n");}
      void length1();
      int length_mp(int phase);
      void length2(int phase); 
      int coarsen(FLT factor, const class mesh& xmesh);
      void coarsen2(FLT factor, const class mesh& inmesh, class mesh& work);
      void refineby2(const class mesh& xmesh);
      void swap(FLT swaptol = 0.0);
      void yaber(FLT tolsize, int yes_swap, FLT swaptol = 0.0);
      inline void treeupdate() { qtree.update(0,nvrtx);}
      void rebay(FLT tolsize);
      inline void adapt(FLT tolerance) {
         yaber(1.0/tolerance,1,0.0);
//         setbcinfo();
//         out_mesh("coarse",grid);
         treeupdate();
         rebay(tolerance);
//        setbcinfo();
//        out_mesh("refine",grid);
         return;
      }
      
      /* CENTER OF AREA MESH SMOOTHING */
      int smooth_cofa(int niter);

      /* FIND MATCHING MESH BOUNDARIES */
      void findmatch(class mesh& tgt);
      FLT *vbuff[MAXSB];
      void matchboundaries1();
      int matchboundaries_mp(int phase);
      void matchboundaries2(int phase);

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
      static int ntdel, tdel[MAXLST+1];
      static int nsdel, sdel[MAXLST+1];
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

/*	MAKE SURE MATCHING BOUNDARIES ARE AT EXACTLY THE SAME POSITIONS */
void inline mesh::matchboundaries1() {
   for(int i=0;i<nsbd;++i) 
      sbdry[i]->sendpositions(0);
}

int inline mesh::matchboundaries_mp(int phase) {
   int stop;
   
   for(int i=0;i<nsbd;++i) 
      sbdry[i]->rcvpositions(phase);
   
   stop = 1;
   for(int i=0;i<nsbd;++i) 
      stop &= sbdry[i]->sendpositions(phase+1);
   
   return(stop);
}

void inline mesh::matchboundaries2(int phase) {
   for(int i=0;i<nsbd;++i)
      sbdry[i]->rcvpositions(phase);
}

/* SOME CONVENIENT BOUNDARY TYPES */
#define FSRF_MASK (1<<0)
#define IFCE_MASK (1<<1)
#define INFL_MASK (1<<2)
#define OUTF_MASK (1<<3)
#define SYMM_MASK (1<<4)
#define EULR_MASK (1<<5)
#define PRDX_MASK (1<<6)   
#define PRDY_MASK (1<<7)
#define COMX_MASK (1<<8)
#define COMY_MASK (1<<9)
#define CURV_MASK (1<<10)
#define PRMT_MASK (1<<11)
#define XDIR_MP COMX_MASK
#define YDIR_MP (COMY_MASK +IFCE_MASK)
#define ALLD_MP (XDIR_MP +YDIR_MP)

#endif
