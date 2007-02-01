/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ps.h"
#include "../hp_boundary.h"

 void tri_hp_ps::init(input_map& input, gbl *gin) {
   FLT nu, E;
   bool coarse;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   keyword = idprefix + ".nvariable";
   input[keyword] = "3";
   
   tri_hp::init(input,gin);
   
   /* Load pointer to block stuff */
   gbl_ptr = gin;
   
   keyword = idprefix + ".coarse";
   input.getwdefault(keyword,coarse,false);
   
   keyword = idprefix + ".dissipation";
   input.getwdefault(keyword,adis,1.0);
   
   if (coarse) return;
   
   gbl_ptr->tau.resize(maxvst);

   keyword = idprefix + ".nu";
   input.getwdefault(keyword,nu,0.0);
   
   keyword = idprefix + ".E";
   input.getwdefault(keyword,E,0.0);
   
   gbl_ptr->mu = E/(2.*(1.+nu));
   gbl_ptr->lami = (1.+nu)*(1.-2.*nu)/(E*nu);
  
   return;
}

/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW TO DO NOTHING */
void tri_hp_ps::calculate_unsteady_sources(bool coarse) {}
