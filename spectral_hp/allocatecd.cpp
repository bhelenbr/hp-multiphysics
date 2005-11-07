/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"

 void tri_hp_cd::init(std::map <std::string,std::string>& input, std::string prefix, gbl *gin) {
   int coarse;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   
   tri_hp::init(input,prefix,gin);
   
   keyword = prefix + ".coarse";
   data.str(input[keyword]);
   if (!(data >> coarse)) {
      coarse = 0;
   }
   data.clear();
   
   if (coarse) return;
   
   /* Load pointer to block stuff */
   cd_gbl = gin;
  
   keyword = prefix + ".ax";
   data.str(input[keyword]);
   if (!(data >> cd_gbl->ax)) {
      cd_gbl->ax = 1.0;
   }
   data.clear();
   
   keyword = prefix + ".ay";
   data.str(input[keyword]);
   if (!(data >> cd_gbl->ay)) {
      cd_gbl->ax = 0.0;
   }
   data.clear();
   
   keyword = prefix + ".nu";
   data.str(input[keyword]);
   if (!(data >> cd_gbl->nu)) {
      cd_gbl->nu = 1.0;
   }
   data.clear();
   
   cd_gbl->tau.resize(maxvst);
   
   return;
}
   
   
         