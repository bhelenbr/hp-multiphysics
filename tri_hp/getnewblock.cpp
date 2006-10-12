/*
 *  initfunc.cpp
 *  planar++
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#include "blocks.h"
#include "block.h"
#include "mgblock.h"

#include "r_mesh.h"
#include "tri_hp_cd.h"
#include "tri_hp_ins.h"
#include "tri_hp_ps.h"
#include "tri_hp_swirl.h"
#include "pod/pod.h"

class btype {
   public:
      const static int ntypes = 7;
      enum ids {r_mesh,cd,ins,ps,swirl,pod_ins,pod_cd};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         int i;
         for(i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i);
         return(-1);
      }
};
const char btype::names[ntypes][40] = {"r_mesh","cd","ins","ps","swirl","pod_ins","pod_cd"};


block* blocks::getnewblock(int idnum, input_map& blockdata) {
   std::string keyword,val,ibcname,srcname;
   std::istringstream data;
   char idntystring[10];
   int type;        
   
   /* FIND BLOCK TYPE */
   sprintf(idntystring,"b%d",idnum);
   keyword = std::string(idntystring) + ".type";

   if (blockdata.get(keyword,val)) {
      type = btype::getid(val.c_str());
   }
   else {
      if (!blockdata.get("blocktype",val)) {
         *sim::log << "couldn't find block type" << std::endl;
         exit(1);
      }
      type = btype::getid(val.c_str());
   }
         
   switch(type) {
      case btype::r_mesh: {
         mgrid<r_mesh> *temp = new mgrid<r_mesh>(idnum);
         return(temp);
      }
      
      case btype::cd: {
         mgrid<tri_hp_cd> *temp = new mgrid<tri_hp_cd>(idnum);
         return(temp);
      }
      
      case btype::ins: {
         mgrid<tri_hp_ins> *temp = new mgrid<tri_hp_ins>(idnum);
         return(temp);
      }
      
      case btype::ps: {
         mgrid<tri_hp_ps> *temp = new mgrid<tri_hp_ps>(idnum);
         return(temp);
      }
		
		case btype::swirl: {
         mgrid<tri_hp_swirl> *temp = new mgrid<tri_hp_swirl>(idnum);
         return(temp);
      }
      case btype::pod_ins: {
         mgrid<pod<tri_hp_ins> > *temp = new mgrid<pod<tri_hp_ins> >(idnum);
         return(temp);
      }
      
      case btype::pod_cd: {
         mgrid<pod<tri_hp_cd> > *temp = new mgrid<pod<tri_hp_cd> >(idnum);
         return(temp);
      }

      default: {
         std::cout << "unrecognizable block type: " <<  type << std::endl;
         mgrid<r_mesh> *temp = new mgrid<r_mesh>(idnum);
         return(temp);
      }
   } 
      
   return(0);
}

