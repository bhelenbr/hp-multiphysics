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

#include "../ins/tri_hp_ins.h"
#include <blocks.h>

class tri_hp_lvlset : public tri_hp_ins {
   public:
      struct gbl : public tri_hp_ins::gbl {
         /* PHYSICAL CONSTANTS */
         FLT sigma, width;
         FLT rho2, mu2;

      } *gbl_ptr;
      hp_side_bdry* getnewsideobject(int bnum, input_map &bdrydata);

   
   private:
      int excpt;
      
   public:
      void init(input_map& input, gbl *gin); 
      tri_hp_lvlset* create() { return new tri_hp_lvlset(); }
      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
      void calculate_unsteady_sources(bool coarse);
      
      FLT heavyside(FLT phidw) {
         return(0.5*(phidw +sin(M_PI*phidw)/M_PI));
      }
      
      FLT delta(FLT phidw) {
         return(0.5/gbl_ptr->width*(1.+cos(M_PI*phidw)));
      }
      
      FLT heavyside_if(FLT phidw) {
         if (phidw < -1.0) return(0.0);
         if (phidw >  1.0) return(1.0);
         return(heavyside(phidw));
      }
      
      void heavyside_and_delta_if(FLT phidw, FLT& h, FLT& d) {
         if (phidw < -1.0) {
            h = 0.0;
            d = 0.0;
         }
         else if (phidw >  1.0) {
            h = 1.0;
            d = 0.0;
         }
         else {
            h = heavyside(phidw);
            d = delta(phidw);
         }
         return;
      }         

};
#endif
