#ifndef _r_mesh_h_
#define _r_mesh_h_

#include "mesh.h"
#include "rboundary.h"
/* AN R-DEFORMABLE MULTI-GRID MESH OBJECT */
/* GOOD COMBINATIONS ARE: FIXX & FIX2X, FIXX & !FIX2X, AND !FIXX and !FIX2X */
/* DON'T DO !FIXX and FIX2X */

#define NO_FOURTH
#define GEOMETRIC

/* MESH DEFORMATION VARIABLES */

class r_mesh :public mesh {
      public:
         /* MESH INDEPENDENT VARIABLES FOR MGRID SEQUENCE */
         struct r_gbl {
            FLT (*work)[ND];
            FLT (*res)[ND];
            FLT *diag;
         } *rg;
         
      private:
         /* THINGS SHARED BY ALL BLOCKS */
         static FLT vnn, fadd;
         
         /* MESH VARIABLES */
         FLT *ksprg;
         FLT (*src)[ND];
         FLT *kvol;
         FLT (*vrtx_frst)[ND];
         bool isfrst;
         
         side_boundary* getnewsideobject(int type);
                  
      public:
         /* SETUP FUNCTION */
         void gbl_alloc(r_mesh::r_gbl *store);
         void allocate(bool coarse, r_mesh::r_gbl *rginit);

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
         
         /* TESTS */
         void tadvance() {
            for(int i=0;i<nsbd;++i) 
               sbdry[i]->tadvance();
         }
         
         /* TO SET UP ADAPTATION VLENGTH */
         void length1();
         void length_mp();
         void length2();
         
         friend class block;
         friend class blocks;
};
#endif
