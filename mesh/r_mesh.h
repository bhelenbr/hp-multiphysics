#ifndef _r_mesh_h_
#define _r_mesh_h_

#include "mesh.h"
#include "block.h"
#include <map>
#include <string>
#include <sstream>
#include <iostream>

/* AN R-DEFORMABLE MULTI-GRID MESH OBJECT */
/* GOOD COMBINATIONS ARE: FIXX & FIX2X, FIXX & !FIX2X, AND !FIXX and !FIX2X */
/* DON'T DO !FIXX and FIX2X */

#define NO_FOURTH
#define GEOMETRIC

class r_side_bdry;

class r_mesh : public mesh {
      protected:         
         FLT vnn;
         FLT fadd;
         
         /* MESH VARIABLES */
         FLT *ksprg;
         FLT (*src)[ND];
         FLT *kvol;
         FLT (*vrtx_frst)[ND];
         bool isfrst;
         int mp_phase;
         
         r_side_bdry *r_sbdry[MAXSB];
         r_side_bdry* getnewsideobject(int bnum, std::map<std::string,std::string> *bdrydata);
         
          /* SETUP SPRING CONSTANTS */
         /* LAPLACE CONSTANTS */
         void rklaplace();

         /* NEEDED FOR BIHARMONIC METHOD */
         void calc_kvol();
         void kvoli();         
         /* SPRING METHOD */
         void rksprg();
         
         /* CALCULATE COARSE SPRING CONSTANTS */
         /* USING INTERPOLATION OPERATORS */
         /* SHOULD WORK??? FOR BIHARMONIC, LAPLACIAN, & SPRING */
         void rkmgrid(mesh::transfer *fv_to_ct, r_mesh *fmesh);
   
         /* CALCULATE RESIDUAL */
         void rsdl();

         /* CALCLATE SOURCE TERM */
         void zero_source();
         void sumsrc();

         /* CALCULATE VOL*DT */
         void vddt();
         void vddti();

#ifdef FOURTH
         /* TWO STEP PROCEDURE FOR 4TH ORDER */
         void rsdl1();
         void vddt1();
#endif
         void moveboundaries();

      public:
         /* POINTER TO STABLE THINGS SHARED IN BLOCK CONTAINER */
         /* EMPTY BECAUSE R_MESH DOESN'T NEED TO SHARE ANY STABLE INFORMATION */
         struct gbl {} *rg;  
         /* SCRATCH VARIABLES FOR MGRID SEQUENCE */
         /* CAN BE SHARED BETWEEN MGLEVELS BUT NOT DIFFERENT BLOCKS */
         /* ALLOCATED FROM scratch */
         /* NEED TO BE PUBLIC SO THEY CAN BE MANIPULATED BY B.C.'s */
         FLT (*res)[ND];
         FLT (*work)[ND];
         FLT *diag;
         
         /* ACCESSOR FUNCTIONS FOR COMPATIBILITY WITH MGBLOCK */
         sharedmem* init(bool coarse, std::map <std::string,std::string>& input, std::string prefix, gbl *rgin, sharedmem *wkin = 0);
         void load_scratch_pointers();
         void bdry_output(const char *filename) const;
         block::ctrl mg_getfres(int excpt,mesh::transfer *fv_to_ct, mesh::transfer *cv_to_ft, r_mesh *fmesh);
         block::ctrl mg_getcchng(int excpt,mesh::transfer *fv_to_ct, mesh::transfer *cv_to_ft, r_mesh *cmesh);
         block::ctrl tadvance(bool coarse,int execpoint,mesh::transfer *fv_to_ct,mesh::transfer *cv_to_ft, r_mesh *fmesh);
         block::ctrl rsdl(int excpt);
         block::ctrl update(int excpt);
         block::ctrl setup_preconditioner(int excpt);
         void maxres();
};
#endif


