#ifndef _r_mesh_h_
#define _r_mesh_h_

#include "mesh.h"
/* AN R-DEFORMABLE MULTI-GRID MESH OBJECT */
/* GOOD COMBINATIONS ARE: FIXX & FIX2X, FIXX & !FIX2X, AND !FIXX and !FIX2X */
/* DON'T DO !FIXX and FIX2X */

#define NO_FOURTH
#define GEOMETRIC

/* MESH DEFORMATION VARIABLES */
class rbdry_interface;

class r_mesh :public mesh {
      public:
         /* THINGS SHARED BY ALL BLOCKS */
         static FLT vnn, fadd;
         
         /* MESH INDEPENDENT VARIABLES FOR MGRID SEQUENCE */
         struct gbl {
            FLT (*work)[ND];
            FLT (*res)[ND];
            FLT *diag;
         } *rg;
         
      private:         
         /* MESH VARIABLES */
         FLT *ksprg;
         FLT (*src)[ND];
         FLT *kvol;
         FLT (*vrtx_frst)[ND];
         bool isfrst;
         
         rbdry_interface *rbdry[MAXSB];
         void getnewsideobject(int bnum, int type);
         
         /* SETUP FUNCTION */
         /* CALLED BY ALLOCATE WHEN COARSE IS FALSE */
         void gbl_alloc(r_mesh::gbl *store);
                  
      public:
         void allocate(bool coarse, r_mesh::gbl *rginit);

         /* SETUP SPRING CONSTANTS */
         /* LAPLACE CONSTANTS */
         void rklaplace();

         /* NEEDED FOR BIHARMONIC METHOD */
         void calc_kvol();
         int kvol_mp(int phase);
         void kvoli(int phase);         

         /* SPRING METHOD */
         void rksprg();
         
         /* CALCULATE COARSE SPRING CONSTANTS */
         /* USING INTERPOLATION OPERATORS */
         /* SHOULD WORK??? FOR BIHARMONIC, LAPLACIAN, & SPRING */
         void rkmgrid();
         int rkmgrid_mp(int phase);
         void rkmgridi(int phase);
   
         /* CALCULATE RESIDUAL */
         void rsdl();
         int rsdl_snd(int phase);
         void rsdl_rcv(int phase);
         void update();
         void maxres();

         /* CALCLATE SOURCE TERM */
         void zero_source();
         void sumsrc();

         /* CALCULATE VOL*DT */
         void vddt();
         int vddt_mp(int phase);
         void vddti(int phase);

#ifdef FOURTH
         /* TWO STEP PROCEDURE FOR 4TH ORDER */
         void rsdl1();
         int rsdl1_mp(int phase);
         void vddt1();
         int vddt1_mp(int phase);
#endif

         /* MGRID TRANSFER */
         void mg_getfres();
         void mg_getcchng();
         
         /* TESTS */
         void tadvance();
};

#include "rboundary.h"

inline void r_mesh::tadvance() {
   for(int i=0;i<nsbd;++i) 
      rbdry[i]->tadvance();
}
#endif
