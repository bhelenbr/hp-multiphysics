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
         mesh::filetype output_type;
         
         /* MESH VARIABLES */
         Array<FLT,1> ksprg;
         Array<FLT,1> kvol;
         Array<TinyVector<FLT,ND>,1> src;
         Array<TinyVector<FLT,ND>,1> vrtx_frst;
         bool isfrst;
         int mp_phase;
         
         Array<r_side_bdry *,1> r_sbdry;
         r_side_bdry* getnewsideobject(int bnum, input_map *bdrydata);
         
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
         void rkmgrid(Array<mesh::transfer,1> &fv_to_ct, r_mesh *fmesh);
   
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
         /* POINTER TO THINGS SHARED IN BLOCK CONTAINER */
         /* SCRATCH VARIABLES FOR MGRID SEQUENCE */
         /* CAN BE SHARED BETWEEN MGLEVELS BUT NOT DIFFERENT BLOCKS */
         /* NEED TO BE PUBLIC SO THEY CAN BE MANIPULATED BY B.C.'s */
         struct gbl {
            Array<FLT,1> diag;
            Array<TinyVector<FLT,2>,1> res;
            Array<TinyVector<FLT,2>,1> res1;
         } *rg;  
         ~r_mesh();
         
         /* ACCESSOR FUNCTIONS FOR COMPATIBILITY WITH MGBLOCK */
         void init(input_map& input, gbl *rgin);
         void input(const std::string &fname) {mesh::input(fname,mesh::grid); }
         void output(const std::string &outname,block::output_purpose why) {mesh::output(outname,output_type); }
         void bdry_output(const std::string &filename) const;
         block::ctrl mg_getfres(int excpt,Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, r_mesh *fmesh);
         block::ctrl mg_getcchng(int excpt,Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, r_mesh *cmesh);
         block::ctrl tadvance(bool coarse,int execpoint,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, r_mesh *fmesh);
         block::ctrl rsdl(int excpt);
         block::ctrl update(int excpt);
         block::ctrl setup_preconditioner(int excpt);
         block::ctrl length(int excpt) {return(block::stop);}
         void maxres();
};
#endif


