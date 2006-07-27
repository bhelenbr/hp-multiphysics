/*
 *  getnewswirl_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_swirl.h"
#include "bdry_swirl.h"
#include <input_map.h>
#include "bdry_ins.h"

using namespace bdry_swirl;

/** \brief Helper object for side_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_swirl_stype {
   public:
      static const int ntypes = 5;
      enum ids {unknown=-1,plain,inflow,outflow,euler,symmetry};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i);
         return(unknown);
      }
};

const char tri_hp_swirl_stype::names[ntypes][40] = {"plain","inflow","outflow","euler","symmetry"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_side_bdry* tri_hp_swirl::getnewsideobject(int bnum, input_map& bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   int type;        
   hp_side_bdry *temp;  
   
   keyword =  sbdry(bnum)->idprefix + ".swirl_type";   
   if (bdrydata.get(keyword,val)) {
      type = tri_hp_swirl_stype::getid(val.c_str());
      if (type == tri_hp_swirl_stype::unknown)  {
         *sim::log << "unknown side type:" << val << std::endl;
         exit(1);
      }
   }
   else {
      type = tri_hp_swirl_stype::unknown;
   }

   switch(type) {
      case tri_hp_swirl_stype::plain: {
         temp = new hp_side_bdry(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_swirl_stype::inflow: {
         temp = new inflow(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_swirl_stype::outflow: {
         temp = new neumann(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_swirl_stype::euler: {
         temp = new euler(*this,*sbdry(bnum));
         break;
      }
		case tri_hp_swirl_stype::symmetry: {
         temp = new symmetry(*this,*sbdry(bnum));
         break;
      }	
      default: {
         temp = tri_hp_ins::getnewsideobject(bnum,bdrydata);
         break;
      }
   }
   
   return(temp);
}


