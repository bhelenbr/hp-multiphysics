/*
 *  hp_mgrid.h
 *  swe++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_swe_h_
#define _tri_hp_swe_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>

class tri_hp_swe : public tri_hp_ins {
   public:
      struct gbl : public tri_hp_ins::gbl {
      
         /* PHYSICAL CONSTANTS */
         FLT f0, beta;
         
         /* BATHYMETRY DATA */
         init_bdry_cndtn *bathy;
         
      } *gbl_ptr;

      init_bdry_cndtn* getnewbathy(input_map& inmap);
		hp_side_bdry* getnewsideobject(int bnum, input_map &bdrydata);
		init_bdry_cndtn* getnewibc(input_map& inmap);
   
   private:
      int excpt;
      
   public:
      void init(input_map& input, gbl *gin); 
      tri_hp_swe* create() { return new tri_hp_swe(); }
      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
      void calculate_unsteady_sources(bool coarse) {tri_hp::calculate_unsteady_sources(coarse);}
};
#endif
