/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_ins_h_
#define _tri_hp_ins_h_

#include "../tri_hp.h"
#include <blocks.h>

class tri_hp_ins : public tri_hp {
   public:
      /* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
      struct gbl : public tri_hp::gbl {
         /* STABILIZATION */
         Array<FLT,2> tau;
            
         /* PHYSICAL CONSTANTS */
         FLT rho, mu;
         Array<FLT,1> D;
         
         /* STORAGE FOR CALCULATION OF ENERGY AND AREA */
         TinyVector<FLT,2> eanda, eanda_recv;

      } *gbl_ptr;

#ifdef DROP
      /** Rigid Mesh Motion for Finding Steady Translating Solutions */
      static TinyVector<FLT,ND> mesh_ref_vel;
#endif
      
      FLT adis; // DISSIPATION CONSTANT
      
      hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
      hp_side_bdry* getnewsideobject(int bnum, input_map &bdrydata);
      init_bdry_cndtn* getnewibc(input_map& inmap);
      mesh_mover* getnewmesh_mover(input_map& inmap);

   
   private:
      int excpt;
      
   public:
      void init(input_map& input, gbl *gin); 
      tri_hp_ins* create() { return new tri_hp_ins(); }
      block::ctrl length(block::ctrl ctrl_message);
      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
      void calculate_unsteady_sources(bool coarse);
};
#endif
