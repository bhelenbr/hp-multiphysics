/*	AN R-DEFORMABLE MULTI-GRID MESH OBJECT */
/*	GOOD COMBINATIONS ARE: FIXX & FIX2X, FIXX & !FIX2X, AND !FIXX and !FIX2X */
/* DON'T DO !FIXX and FIX2X */
#include"mesh.h"

#define FOURTH

#define FIXX_MASK (0xFF -PRDY_MASK)
#define FIXY_MASK (0xFF -PRDX_MASK -SYMM_MASK) 

#define FIX2X_MASK (0xFF -FSRF_MASK -IFCE_MASK -PRDX_MASK -PRDY_MASK -EULR_MASK)
#define FIX2Y_MASK (0xFF -FSRF_MASK -IFCE_MASK -PRDX_MASK -PRDY_MASK -EULR_MASK -SYMM_MASK)

/*	MESH INDEPENDENT VARIABLES FOR MGRID SEQUENCE */
struct r_mesh_glbls {
   FLT (*work)[ND];
   FLT (*res)[ND];
   FLT *diag;
   FLT vnn;
   FLT fadd;
};

class r_mesh :public mesh {
   /* MESH DEFORMATION VARIABLES */
      protected:
/*			STRUCTURE CONTAINING ALL MESH INDEPENDENT VARIABLES */
         struct r_mesh_glbls *rg;
         
/*    	MESH VARIABLES */
         FLT *ksprg;
         FLT (*src)[ND];
         FLT *kvol;
         FLT (*vrtx_frst)[ND];
         bool isfrst;
         
/*			MGRID MESH POINTERS */
         class r_mesh *cmesh;
         class r_mesh *fmesh;

      public:
         void init(bool coarse, struct r_mesh_glbls *rginit);

/*    	SETUP SPRING CONSTANTS */
/*			LAPLACE CONSTANTS */
         void rklaplace();

/*			NEEDED FOR BIHARMONIC METHOD */
         void calc_kvol();
         void kvol_mp();
         void kvoli();         

/*			SPRING METHOD */
         void rksprg();
         
/*			CALCULATE COARSE SPRING CONSTANTS */
/*			USING INTERPOLATION OPERATORS */
/*			SHOULD WORK??? FOR BIHARMONIC, LAPLACIAN, & SPRING */
         void rkmgrid();
         void rkmgrid_mp();
         void rkmgridi();
   
/*    	CALCULATE RESIDUAL */
         void rsdl();
         void rsdl_mp();
         void update(int lvl);

/*			CALCLATE SOURCE TERM */
         void source();
         void sumsrc();

/*			CALCULATE VOL*DT */
         void vddt();
         void vddt_mp();
         void vddti();

#ifdef FOURTH
/*			TWO STEP PROCEDURE FOR 4TH ORDER */
         void rsdl1();
         void rsdl1_mp();
         void vddt1();
         void vddt1_mp();
#endif
         
/*    	MGRID TRANSFER */
         void mg_getfres();
         void mg_getcchng();
         int setfine(class r_mesh& tgt);
         int setcoarse(class r_mesh& tgt);  
         
/*			COMMUNICATION BOUNDARIES */
         void send(int MASK, FLT *base,int bgn,int end, int stride);
         void rcv(int MASK, FLT *base,int bgn,int end, int stride);
           
/*   		TESTS */
         void perturb(int step);
};
