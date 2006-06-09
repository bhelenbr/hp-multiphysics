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
   bool coarse, adapt_storage;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   keyword = idprefix + ".nvariable";
   input[keyword] = "1";
   
   tri_hp::init(input,gin);
   
   /* Load pointer to block stuff */
   cd_gbl = gin;
     
   keyword = idprefix + ".adapt_storage";
   input.getwdefault(keyword,adapt_storage,false);
   if (adapt_storage) return;
   
   keyword = idprefix + ".coarse";
   input.getwdefault(keyword,coarse,false);
   
   keyword = idprefix + ".dissipation";
   input.getwdefault(keyword,adis,1.0);
   
   if (coarse) return;
  
   keyword = idprefix + ".ax";
   input.getwdefault(keyword,cd_gbl->ax,1.0);

   keyword = idprefix + ".ay";
   input.getwdefault(keyword,cd_gbl->ay,0.0);

   keyword = idprefix + ".nu";
   input.getwdefault(keyword,cd_gbl->nu,1.0);

   cd_gbl->tau.resize(maxvst);
   
   cd_gbl->src = getnewsrc(input);
   cd_gbl->src->input(input,idprefix);
   
   return;
}