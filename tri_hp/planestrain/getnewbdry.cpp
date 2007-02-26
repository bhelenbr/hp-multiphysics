/*
 *  getnewbdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_ps.h"
#include "bdry_ps.h"
#include <input_map.h>

using namespace bdry_ps;

/** \brief Helper object for side_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_ps_stype {
   public:
      static const int ntypes = 4;
      enum ids {unknown=-1,plain,dirichlet,neumann,friction_wall};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i);
         return(unknown);
      }
};

const char tri_hp_ps_stype::names[ntypes][40] = {"plain","dirichlet","neumann","friction_wall"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_side_bdry* tri_hp_ps::getnewsideobject(int bnum, input_map &bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   int type;        
   hp_side_bdry *temp;  
   
   keyword =  sbdry(bnum)->idprefix + "_ps_type";
   if (bdrydata.get(keyword,val)) {
      type = tri_hp_ps_stype::getid(val.c_str());
      if (type == tri_hp_ps_stype::unknown)  {
         *sim::log << "unknown side type:" << val << std::endl;
         exit(1);
      }
   }
   else {
      type = tri_hp_ps_stype::unknown;
   }
   
   switch(type) {
      case tri_hp_ps_stype::plain: {
         temp = new hp_side_bdry(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ps_stype::dirichlet: {
         temp = new dirichlet(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ps_stype::neumann: {
         temp = new neumann(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ps_stype::friction_wall: {
         temp = new friction_wall(*this,*sbdry(bnum));  // TEMPORARY NOT WORKING YET
         break;
      }
      default: {
         temp = tri_hp::getnewsideobject(bnum,bdrydata);
         break;
      }
   }
   
   return(temp);
}


