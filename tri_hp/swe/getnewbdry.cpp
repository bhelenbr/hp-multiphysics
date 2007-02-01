/*
 *  getnewswe_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_swe.h"
#include "bdry_swe.h"
#include <input_map.h>

using namespace bdry_swe;

/** \brief Helper object for side_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_swe_stype {
   public:
      static const int ntypes = 1;
      enum ids {unknown=-1,wall};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i);
         return(unknown);
      }
};

const char tri_hp_swe_stype::names[ntypes][40] = {"wall"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_side_bdry* tri_hp_swe::getnewsideobject(int bnum, input_map& bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   int type;        
   hp_side_bdry *temp;  
   
   keyword =  sbdry(bnum)->idprefix + ".swe_type";   
   if (bdrydata.get(keyword,val)) {
      type = tri_hp_swe_stype::getid(val.c_str());
      if (type == tri_hp_swe_stype::unknown)  {
         *sim::log << "unknown side type:" << val << std::endl;
         exit(1);
      }
   }
   else {
      type = tri_hp_swe_stype::unknown;
   }

   switch(type) {
		case tri_hp_swe_stype::wall: {
         temp = new wall(*this,*sbdry(bnum));
         break;
      }	
      default: {
         temp = tri_hp_ins::getnewsideobject(bnum,bdrydata);
         break;
      }
   }
   
   return(temp);
}


