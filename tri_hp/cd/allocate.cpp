/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"

void tri_hp_cd::init(input_map& input, gbl *gin) {
   int coarse;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   keyword = idprefix + ".nvariable";
   input[keyword] = "1";
   
   tri_hp::init(input,gin);
   
   /* Load pointer to block stuff */
   cd_gbl = gin;
   
   keyword = idprefix + ".coarse";
   input.getwdefault(keyword,coarse,0);
   
   keyword = idprefix + ".dissipation";
   input.getwdefault(keyword,adis,1.0);
   *sim::log << "#" << keyword << ": " << adis << std::endl;
   
   if (coarse) return;
  
   keyword = idprefix + ".ax";
   input.getwdefault(keyword,cd_gbl->ax,1.0);
   *sim::log << "#" << keyword << ": " << cd_gbl->ax << std::endl;

   keyword = idprefix + ".ay";
   input.getwdefault(keyword,cd_gbl->ay,0.0);
   *sim::log << "#" << keyword << ": " << cd_gbl->ay << std::endl;

   keyword = idprefix + ".nu";
   input.getwdefault(keyword,cd_gbl->nu,1.0);
   *sim::log << "#" << keyword << ": " << cd_gbl->nu << std::endl;

   cd_gbl->tau.resize(maxvst);
   
   cd_gbl->src = getnewsrc(input);
   cd_gbl->src->input(input,idprefix);
   
   return;
}