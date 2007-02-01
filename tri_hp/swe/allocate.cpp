/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swe.h"
#include "../hp_boundary.h"

 void tri_hp_swe::init(input_map& input, gbl *gin) {
   std::string keyword;
   
   tri_hp_ins::init(input,gin);
   
   /* Load pointer to block stuff */
   gbl_ptr = gin;

   if (!input.get(idprefix + ".f0",gbl_ptr->f0)) input.getwdefault("rho",gbl_ptr->f0,0.0);
   if (!input.get(idprefix + ".beta",gbl_ptr->beta)) input.getwdefault("mu",gbl_ptr->beta,1.0);
   
   gbl_ptr->bathy = getnewbathy(input);
   gbl_ptr->bathy->input(input,idprefix);
   
   return;
}
