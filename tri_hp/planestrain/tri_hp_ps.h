/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_ps_h_
#define _tri_hp_ps_h_

#include "tri_hp.h"
#include <blocks.h>

class tri_hp_ps : public tri_hp {
   public:
      /* THINGS SHARED BY ALL tri_hp_ps in same multigrid block */
      struct gbl : public tri_hp::gbl {
         /* STABILIZATION */
         Array<FLT,1> tau;
            
         /* PHYSICAL CONSTANTS */
         FLT mu, lami;
         
         /* STORAGE FOR CALCULATION OF ENERGY AND AREA */
         TinyVector<FLT,2> eanda, eanda_recv;

      } *gbl_ptr;
      
      FLT adis; // DISSIPATION CONSTANT
      
      hp_side_bdry* getnewsideobject(int bnum, input_map &bdrydata);
   
   private:
      int excpt;
      
   public:
      void init(input_map& input, gbl *gin); 
      tri_hp_ps* create() { return new tri_hp_ps(); }
      block::ctrl length(block::ctrl ctrl_message);
      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
      void calculate_unsteady_sources(bool coarse);
};
#endif