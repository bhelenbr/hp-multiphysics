#define FLT double

#define ND 2
#define MAXSB 16
#define MAXLST 1000

/*	FIRST 16 BITS DETERMINE TYPE */
/* LAST 16 BITS GIVE UNIQUE IDENTITY NUMBER */
/*	IDENTITY NUMBER SHOULD BE THE SAME FOR SIDES WHICH MATCH */
/*	IDENTITIY NUMBER IS UNIQUE ONLY FOR EACH TYPE */
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

/* MASK FOR ALL MESSAGE PASSING BOUNDARIES */
/* X/Y SPLIT ONLY WORKS FOR FACE TO FACE BOUNDARIES & 2X2 BOUNDARY INTERSECTIONS */
#define XDIR_MP COMX_MASK
#define YDIR_MP (COMY_MASK +IFCE_MASK)
#define ALLD_MP (XDIR_MP +YDIR_MP)
                             
#if (FLT == double)
#define EPSILON DBL_EPSILON
#else
#define EPSILON FLT_EPSILON
#endif

#include<quadtree.h>


enum FILETYPE {easymesh, gambit, tecplot, text, binary};

/*	THIS ROUTINE MUST BE SUPPLIED IF USING CURVED BOUNDARIES */
/*	GIVEN AN INPUT X,Y SHOULD MOVE POINT ONTO CURVED SURFACE */
extern void mvpttobdry(int typ, FLT& x, FLT &y);

class mesh {

   protected:
      int initialized, maxvst;

/*    VERTEX DATA */
      int nvrtx;
      FLT (*vrtx)[ND];
      FLT *vlngth;
      int *vinfo;
      int *nnbor;
      int *vtri;
      class quadtree qtree;
   
/*    VERTEX BOUNDARY INFO */
      int nvbd;
      int maxvbel;  
      struct boundary {
         int type;
         int num;
         int *el;
/*			STUFF FOR COMMUNICATION BOUNDARIES */
         class mesh *adjmesh;
         int adjbnum;
         void *misc; // VOID POINTER TO HANG ANY OTHER WEIRD STUFF
      } vbdry[MAXSB];

/*    SIDE DATA */      
      int nside;
      int (*svrtx)[ND];
      int (*stri)[ND];
      int *sinfo;
      
/*    SIDE BOUNDARY INFO */
      int nsbd;
      int maxsbel;
      struct boundary sbdry[MAXSB];

/*    TRIANGLE DATA */      
      int ntri;
      int (*tvrtx)[3];
      int (*ttri)[3];
      struct tsidedata {
         int side[3];
         int sign[3];
      } *tside;
      int *tinfo;

/*		THIS IS AN INTERPOLATION STRUCTURE IF USING MULTIGRID */
/*		UNALLOCATED UNLESS CONNECT IS CALLED */      
      struct mg_trans {
         int tri;
         double wt[3];
      };
      struct mg_trans *fine;
      struct mg_trans *coarse;
      int mgconnect(struct mg_trans *cnnct, const class mesh& tgt);
      
/*		THESE ARE STATIC WORK ARRAYS FOR MESH REFINEMENT ROUTINES */
      static FLT *fltwk;
      static int gblmaxvst, *intwk1, *intwk2, *intwk3;
      static int maxsrch;  // MAX TRIS TO SEARCH IN FINDTRI ROUTINE

/*    SETUP FUNCTIONS */   
      void cnt_nbor(void);
      void bdrysidereorder(int sidenum);
      void bdrygroupreorder(void);
      void bdrylabel(void);
      void createsideinfo(void);
      void createttri(void);
      void createvtri(void);
      void treeinit();
      void initvlngth();

/*		MESH MODIFICATION FUNCTIONS */
/*		TO CREATE AN INITIAL TRIANGUlATION */
      void triangulate(int **sidelst, int *nside, int nloop, int cknbor = 1);
      void findpt(int *vrtxlst,int nv,int *v,int chkadj,int good[], int &ngood);
      void addtri(int v0,int v1,int v2,int sind,int dir);

/*		TO INSERT A POINT */
      void insert(FLT x, FLT y);
      void insert(int tind, int vnum);
      void bdry_insert(int tind, int snum, int vnum);
      int findtri(FLT x, FLT y, int vnear) const;
      
/*		TO DELETE A SIDE */
      int collapse(int sind);
      void dltvrtx(int vind);
      void dltside(int sind);
      void dlttri(int tind);
      
/*		ORDERED LIST FUNCTIONS */
      void putinlst(int sind);
      void tkoutlst(int sind);
      
/*		TO EDGE SWAP FOR IMPROVED QUALITY */
      int swap(int sind, FLT tol = 0.0);  //SWAPS A SINGLE SIDE 
      void swap(int nswp, int *swp, FLT tol = 0.0); //SWAPS A LIST OF SIDES

/*		PRIMITIVE FUNCTIONS */
      FLT distance(const int v0, const int v1) const;
      FLT incircle(int tind, double *a) const;
      FLT area(int sind, int vind) const;
      FLT area(int v0, int v1, int v2) const;
      FLT area(int tind) const;
      FLT minangle(int v0, int v1, int v2) const;
      FLT angle(int v0, int v1, int v2) const;
      FLT tradius(int tind) const;
      void tcenter(int tind, FLT &xpt, FLT &ypt) const;
      FLT aspect(int tind) const;
      FLT intri(int tind, FLT x, FLT y) const;
      void getwgts(FLT *) const;

   public:
/*		DEFAULT INITIALIZATION */
      mesh() : initialized(0), fine(NULL), coarse(NULL) {};
/*		COPY CONSTRUCTOR TO DUPLICATE MESHES */
      mesh(const mesh& copy) : initialized(0), fine(NULL), coarse(NULL) { printf("in mesh copy\n"); *this = copy;}
      mesh& operator=(const mesh& tgt) {this->copy(tgt); return(*this);}
      void copy(const mesh& tgt);
      
/*		ALLOCATE COMMUNICATION BUFFERS */
/*		THESE HAVE TO BE PUBLIC TO RECEIVE MESSAGES */
      FLT *vbuff[MAXSB];
      FLT *sbuff[MAXSB];
      void inline init_comm_buf(int factor) {
         for(int i=0;i<nvbd;++i)
            if (vbdry[i].type & ALLD_MP)
               vbuff[i] = new FLT[factor];

         for(int i=0;i<nsbd;++i)
            if (sbdry[i].type & ALLD_MP)
               sbuff[i] = new FLT[factor*maxsbel];
      }

/*    INPUT/OUTPUT MESH (MAY MODIFY VINFO/SINFO/TINFO) */
      int in_mesh(char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1);
      int out_mesh(const char *filename, FILETYPE filetype = easymesh) const;
      void bcinfo();
      
/*		ACCESS TO SIMPLE MESH DATA */
      int max() {return maxvst;}

/*    MESH MODIFICATION */   
      int coarsen(const class mesh& xmesh);
      void density();
      void swap(FLT tolsize = 0.0);
      void yaber(FLT tolsize);
      void rebay(FLT tolsize);

      
/*    CENTER OF AREA MESH SMOOTHING */
      int smooth_cofa(int niter);

/*		FUNCTION TO FIND MATCHING MESH BOUNDARIES */
      int findmatch(class mesh& tgt);
/*		TELLS HOW MANY COMMUNICATION BOUNDARIES THERE ARE */
      int alld_mp();
      
/*		FUNCTION TO ALLOCATE & SET INTERPOLATION WEIGHTS BETWEEN TWO MESHES */
      int setfine(const class mesh& tgt);
      int setcoarse(const class mesh& tgt);
};
   
