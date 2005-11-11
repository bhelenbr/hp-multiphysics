/*
 *  tri_hp_cd.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_cd_h_
#define _tri_hp_cd_h_

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
   public:
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
      
      hp_vrtx_bdry* tri_hp_cd::getnewvrtxobject(int bnum, std::map<std::string,std::string> *bdrydata);
      hp_side_bdry* tri_hp_cd::getnewsideobject(int bnum, std::map<std::string,std::string> *bdrydata);
      
   public:
      void init(std::map <std::string,std::string>& input, std::string prefix, gbl *gin); 
      tri_hp_cd* create() { return new tri_hp_cd(); }
      block::ctrl length(int excpt);
      block::ctrl setup_preconditioner(int excpt);
      block::ctrl rsdl(int excpt, int stage=sim::NSTAGE);
};

#endif
