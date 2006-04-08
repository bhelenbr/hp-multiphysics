/*
 *  hp_mgrid.h
 *  swirl++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_swirl_h_
#define _tri_hp_swirl_h_

#include "tri_hp.h"
#include <blocks.h>

class tri_hp_swirl : public tri_hp {
   public:
      /* THINGS SHARED BY ALL tri_hp_swirl in same multigrid block */
      struct gbl : public tri_hp::gbl {
         /* STABILIZATION */
         Array<FLT,1> tau,delt;
            
         /* PHYSICAL CONSTANTS */
         FLT rho, mu, nu;
         
         /* STORAGE FOR CALCULATION OF ENERGY AND AREA */
         TinyVector<FLT,2> eanda, eanda_recv;

      } *swirl_gbl;
      
      FLT adis; // DISSIPATION CONSTANT
		hp_side_bdry* getnewsideobject(int bnum, input_map &bdrydata);
		init_bdry_cndtn* getnewibc(input_map& inmap);
   
   private:
      int excpt;
      
   public:
      void init(input_map& input, gbl *gin); 
      tri_hp_swirl* create() { return new tri_hp_swirl(); }
      block::ctrl length(block::ctrl ctrl_message);
      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
      void calculate_unsteady_sources(bool coarse);
};
#endif