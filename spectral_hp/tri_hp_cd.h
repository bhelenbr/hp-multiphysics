/*
 *  tri_hp_cd.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp.h"
#include <blocks.h>

#ifdef AXISYMMETRIC
#define RAD(I,J) crd(0)(I,J)
#define RAD1D(I) crd(0)(0,I)
#else
#define RAD(I,J) 1
#define RAD1D(I) 1
#endif

class tri_hp_cd : public tri_hp {
   private:
      /* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
      struct gbl : public tri_hp::gbl {
         /* STABILIZATION */
         Array<FLT,1> tau;
            
         /* PHYSICAL CONSTANTS */
         FLT ax, ay, nu;
         
         /* SOURCE FUNCTION */
         FLT (*src)(FLT x, FLT y);

      } *cd_gbl;
      
      FLT adis; // DISSIPATION CONSTANT
      
   public:
      void init(std::map <std::string,std::string>& input, std::string prefix, gbl *gin); 
      block::ctrl length(int excpt);
      block::ctrl setup_preconditioner(int excpt);
      block::ctrl rsdl(int excpt, int stage);
      void addbflux();
      block::ctrl tadvance(bool coarse,int execpoint,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp_cd *fmesh);
};
