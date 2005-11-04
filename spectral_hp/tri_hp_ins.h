/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#include"tri_hp.h"
// #include"surface.h"


#define MXLG2P 5
#define NSTAGE 5

#ifdef AXISYMMETRIC
#define RAD(I,J) crd(0)(I,J)
#define RAD1D(I) crd(0)(0,I)
#else
#define RAD(I,J) 1
#define RAD1D(I) 1
#endif

class tri_hp_ins : public tri_hp {
   private:
      /* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
      struct gbl : public tri_hp::gbl {
         /* STABILIZATION */
         Array<FLT,1> tau,delt;
            
         /* PHYSICAL CONSTANTS */
         FLT rho, mu, nu;
         
         /* STORAGE FOR CALCULATION OF ENERGY AND AREA */
         TinyVector<FLT,2> eanda, eanda_recv;
         
      } *ins_gbl;

      static FLT cv00[MXGP][MXGP],cv01[MXGP][MXGP],cv10[MXGP][MXGP],cv11[MXGP][MXGP]; // LOCAL WORK ARRAYS
      static FLT e00[MXGP][MXGP],e01[MXGP][MXGP],e10[MXGP][MXGP],e11[MXGP][MXGP]; // LOCAL WORK ARRAYS
      
      FLT adis; // DISSIPATION CONSTANT
      
   public:
      void init(std::map <std::string,std::string>& input, std::string prefix, gbl *gin);
      block::ctrl length(int excpt);
      block::ctrl setup_preconditioner(int excpt);
      block::ctrl rsdl(int excpt);
      block::ctrl tadvance(bool coarse,int execpoint,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp_ins *fmesh);
      
      
      
      /* BOUNDARY CONDITION ROUTINES */
      // void drag(int bdry_id);
      // void integrated_averages(FLT a[]);
};

