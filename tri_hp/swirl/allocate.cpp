/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swirl.h"
#include "../hp_boundary.h"

 void tri_hp_swirl::init(input_map& input, gbl *gin) {
   std::string keyword;
   
   keyword = idprefix + "_nvariable";
   input[keyword] = "4";
   
   tri_hp_ins::init(input,gin);
   
   return;
}
