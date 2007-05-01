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
   
   keyword = idprefix + "_nvariable";
   input[keyword] = "1";
   
   tri_hp::init(input,gin);
   
   /* Load pointer to block stuff */
   gbl_ptr = gin;
     
   keyword = idprefix + "_adapt_storage";
   input.getwdefault(keyword,adapt_storage,false);
   if (adapt_storage) return;
   
   keyword = idprefix + "_coarse";
   input.getwdefault(keyword,coarse,false);
   
   keyword = idprefix + "_dissipation";
   input.getwdefault(keyword,adis,1.0);
   
   if (coarse) return;
  
   keyword = idprefix + "_ax";
   if (!input.get(keyword,gbl_ptr->ax)) input.getwdefault("ax",gbl_ptr->ax,1.0);

   keyword = idprefix + "_ay";
   if (!input.get(keyword,gbl_ptr->ay)) input.getwdefault("ay",gbl_ptr->ay,0.0);

   keyword = idprefix + "_nu";
   if (!input.get(keyword,gbl_ptr->nu)) input.getwdefault("nu",gbl_ptr->nu,0.0);

   gbl_ptr->tau.resize(maxvst);
   
   gbl_ptr->src = getnewsrc(input);
   gbl_ptr->src->input(input,idprefix);
   
   return;
}
