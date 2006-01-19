/*
 *  getnewins_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_ins.h"
#include "ins_bdry.h"
#include <input_map.h>

using namespace ins_bdry;

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_ins_vtype {
   public:
      static const int ntypes = 1;
      enum ids {plain=1};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char tri_hp_ins_vtype::names[ntypes][40] = {"plain"};

hp_vrtx_bdry* tri_hp_ins::getnewvrtxobject(int bnum, input_map *bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   char idntystring[10];
   int type;        
   hp_vrtx_bdry *temp;  
   int idnum = vbdry(bnum)->idnum;
   
   type = idnum&0xffff;

   if (bdrydata) {
      sprintf(idntystring,"v%d",idnum);
      keyword = std::string(idntystring) + ".ins_type";
      
      if ((*bdrydata).get(keyword,val)) {
         type = tri_hp_ins_vtype::getid(val.c_str());
         if (type < 0)  {
            *sim::log << "unknown vertex type:" << val << std::endl;
            exit(1);
         }
      }
   }
   
   switch(type) {
      case tri_hp_ins_vtype::plain: {
         temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
         break;
      }
      default: {
         std::cout << "unrecognizable vrtx type: " <<  type << " idnum: " << idnum << std::endl;
         temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
         break;
      }
   } 
   
   if (bdrydata) temp->init(*bdrydata);
   
   return(temp);
}


/** \brief Helper object for side_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_ins_stype {
   public:
      static const int ntypes = 5;
      enum ids {plain=1,inflow,outflow,characteristic,euler};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char tri_hp_ins_stype::names[ntypes][40] = {"plain","inflow","outflow","characteristic","euler"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_side_bdry* tri_hp_ins::getnewsideobject(int bnum, input_map *bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   char idntystring[10];
   int type;        
   hp_side_bdry *temp;  
   int idnum = sbdry(bnum)->idnum;

   type = idnum&0xff;
   
   if (bdrydata) {
      sprintf(idntystring,"s%d",idnum);
      keyword = std::string(idntystring) + ".ins_type";
      
      if ((*bdrydata).get(keyword,val)) {
         type = tri_hp_ins_stype::getid(val.c_str());
         if (type < 0)  {
            *sim::log << "unknown side type:" << val << std::endl;
            exit(1);
         }
      }
      else {
         type = tri_hp_ins_stype::plain;
      }
   }

   switch(type) {
      case tri_hp_ins_stype::plain: {
         temp = new hp_side_bdry(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::inflow: {
         temp = new inflow(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::outflow: {
         temp = new neumann(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::characteristic: {
         temp = new neumann(*this,*sbdry(bnum));  // TEMPORARY NOT WORKING YET
         break;
      }
      case tri_hp_ins_stype::euler: {
         temp = new euler(*this,*sbdry(bnum));
         break;
      }
   }

   if (bdrydata) temp->init(*bdrydata);
   
   return(temp);
}


