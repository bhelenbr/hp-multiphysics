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

class tri_hp_cd : public tri_hp {
   public:
      /* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
      struct gbl : public tri_hp::gbl {
         /* STABILIZATION */
         Array<FLT,1> tau;
            
         /* PHYSICAL CONSTANTS */
         FLT ax, ay, nu;
         
         /* SOURCE FUNCTION */
         init_bdry_cndtn *src;

      } *cd_gbl;
      
      FLT adis; // DISSIPATION CONSTANT
      
      hp_side_bdry* getnewsideobject(int bnum, input_map &bdrydata);
      init_bdry_cndtn* getnewibc(input_map& inmap);
      init_bdry_cndtn *getnewsrc(input_map &inmap);
   
   private:
      int excpt;
      
   public:
      void init(input_map& input, gbl *gin); 
      tri_hp_cd* create() { return new tri_hp_cd(); }
      block::ctrl length(block::ctrl ctrl_message);
      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
};
#endif
