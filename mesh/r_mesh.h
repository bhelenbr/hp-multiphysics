#ifndef _r_mesh_h_
#define _r_mesh_h_

#include"mesh.h"

/* AN R-DEFORMABLE MULTI-GRID MESH OBJECT */
/* GOOD COMBINATIONS ARE: FIXX & FIX2X, FIXX & !FIX2X, AND !FIXX and !FIX2X */
/* DON'T DO !FIXX and FIX2X */

#define NO_FOURTH

/* MESH INDEPENDENT VARIABLES FOR MGRID SEQUENCE */
struct r_mesh_glbls {
   FLT (*work)[ND];
   FLT (*res)[ND];
   FLT *diag;
};

class r_mesh :public mesh {
   /* MESH DEFORMATION VARIABLES */
      private:
         /* THINGS SHARED BY ALL BLOCKS */
         static FLT vnn, fadd;
         static int fixx_mask, fixy_mask, fix2x_mask, fix2y_mask;
         
         /* STRUCTURE CONTAINING ALL MESH INDEPENDENT VARIABLES */
         struct r_mesh_glbls *rg;
         
         /* MESH VARIABLES */
         FLT *ksprg;
         FLT (*src)[ND];
         FLT *kvol;
         FLT (*vrtx_frst)[ND];
         bool isfrst;
         
         /* SETUP FUNCTION */
         void gbl_alloc(struct r_mesh_glbls *store);

      public:
         void allocate(bool coarse, struct r_mesh_glbls *rginit);

         /* SETUP SPRING CONSTANTS */
         /* LAPLACE CONSTANTS */
         void rklaplace();

         /* NEEDED FOR BIHARMONIC METHOD */
         void calc_kvol();
         void kvol_mp();
         void kvoli();         

         /* SPRING METHOD */
         void rksprg();
         
         /* CALCULATE COARSE SPRING CONSTANTS */
         /* USING INTERPOLATION OPERATORS */
         /* SHOULD WORK??? FOR BIHARMONIC, LAPLACIAN, & SPRING */
         void rkmgrid();
         void rkmgrid_mp();
         void rkmgridi();
   
         /* CALCULATE RESIDUAL */
         void rsdl();
         void rsdl_mp();
         void update();
         void maxres();

         /* CALCLATE SOURCE TERM */
         void source();
         void sumsrc();

         /* CALCULATE VOL*DT */
         void vddt();
         void vddt_mp();
         void vddti();

#ifdef FOURTH
         /* TWO STEP PROCEDURE FOR 4TH ORDER */
         void rsdl1();
         void rsdl1_mp();
         void vddt1();
         void vddt1_mp();
#endif
         
         /* MGRID TRANSFER */
         void mg_getfres();
         void mg_getcchng();
         void setfine(class r_mesh& tgt);
         void setcoarse(class r_mesh& tgt);  
         
         /* COMMUNICATION BOUNDARIES */
         void send(int MASK, FLT *base,int bgn,int end, int stride);
         void rcv(int MASK, FLT *base,int bgn,int end, int stride);
           
         /* TESTS */
         void perturb();
         
         /* TO SET UP ADAPTATION VLENGTH */
         void length1();
         void length_mp();
         void length2();
         
         friend class block;
         friend class blocks;
};
#endif
