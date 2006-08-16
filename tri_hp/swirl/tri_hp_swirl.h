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

#include "tri_hp_ins.h"
#include <blocks.h>

class tri_hp_swirl : public tri_hp_ins {
   public:
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
};
#endif