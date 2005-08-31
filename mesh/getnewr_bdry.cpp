/*
 *  getnewr_bdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2004.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "r_mesh.h"
#include "r_boundary.h"

/* OLD BOUNDARY TYPES */
#define FSRF_MASK (1<<0)
#define IFCE_MASK (1<<1)
#define INFL_MASK (1<<2)
#define OUTF_MASK (1<<3)
#define SYMM_MASK (1<<4)
#define EULR_MASK (1<<5)
#define PRDX_MASK (1<<6)   
#define PRDY_MASK (1<<7)
#define COMX_MASK (1<<8)
#define COMY_MASK (1<<9)
#define CURV_MASK (1<<10)

#define MAPPING

const char r_stype::names[ntypes][40] = {"plain", "fixed", "translating", "oscillating"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
r_side_bdry* r_mesh::getnewsideobject(int bnum, std::map<std::string,std::string> *bdrydata) {
   std::istringstream data;
   std::map<std::string,std::string>::const_iterator mi;
   r_side_bdry *temp;  
   int type = 2;

   if (bdrydata) {
      mi = (*bdrydata).find(sbdry(bnum)->idprefix + ".r_type");
      if (mi != (*bdrydata).end()) {
         type = r_stype::getid((*mi).second.c_str());
         if (type < 0)  {
            *sim::log << "unknown type:" << (*mi).second << std::endl;
            exit(1);
         }
      }
      else {
         *sim::log << "couldn't find type for r_side: " << sbdry(bnum)->idnum << std::endl;
      }
   }

   // *sim::log << "making side " << idnum << std::endl;

   switch(type) {
      case r_stype::plain: {
         temp = new r_side_bdry(*this,*sbdry(bnum));
         break;
      }
      case r_stype::fixed: {
         temp = new r_fixed(*this,*sbdry(bnum)); 
         break;
      }
      case r_stype::translating: {
         temp = new r_translating(*this,*sbdry(bnum)); 
         break;
      }
      case r_stype::oscillating: {
         temp = new r_oscillating(*this,*sbdry(bnum)); 
         break;
      }      
      default: {
         temp = new r_fixed(*this,*sbdry(bnum));
         std::cout << "Don't know this r_side_bdry type\n";
      }
   }
   
   if (bdrydata) temp->input(*bdrydata);
   
   temp->output(*sim::log);

   
   return(temp);
}
