/*
 *  hp_mgrid.h
 *  lvlset++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_lvlset_h_
#define _tri_hp_lvlset_h_

#include "tri_hp_ins.h"
#include <blocks.h>

class tri_hp_lvlset : public tri_hp_ins {
   public:
      struct gbl : public tri_hp_ins::gbl {
         /* STABILIZATION */
         Array<FLT,1> tau_lvlset;
            
         /* PHYSICAL CONSTANTS */
         FLT sigma, width;
         FLT rho2, mu2;

      } *gbl_ptr;

		init_bdry_cndtn* getnewibc(input_map& inmap);
   
   private:
      int excpt;
      
   public:
      void init(input_map& input, gbl *gin); 
      tri_hp_lvlset* create() { return new tri_hp_lvlset(); }
      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
      void calculate_unsteady_sources(bool coarse);

};
#endif
