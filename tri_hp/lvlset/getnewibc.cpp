/*
 *  getnewibc.cpp
 *  tri_hp
 *
 *  Created by Erik Durkish on 2/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_lvlset.h"

namespace ibc_lvlset {
	
	class stationary : public init_bdry_cndtn {
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            switch(n) {
               case(0):
                  return(0.0);
               case(1):
                  return(0.0);
					case(2):
						return(0.0);
					case(3):
						return(0.0);
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
         }
   };	
	
	
   class ibc_type {
      public:
         const static int ntypes = 1;
         enum ids {stationary};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char ibc_lvlset::ibc_type::names[ntypes][40] = {"stationary"};

}


init_bdry_cndtn *tri_hp_lvlset::getnewibc(input_map& inmap) {
	std::string keyword,ibcname;
	int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = std::string(idprefix) + ".ibc";
	if (!inmap.get(keyword,ibcname)) {
		if (!inmap.get("ibc",ibcname)) {
			*sim::log << "couldn't find initial condition type" << std::endl;
		}
	}
	
	type = ibc_lvlset::ibc_type::getid(ibcname.c_str());
		
	switch(type) {
		case ibc_lvlset::ibc_type::stationary: {
         init_bdry_cndtn *temp = new ibc_lvlset::stationary;
         return(temp);
		}
		
		default: {
         return(tri_hp_ins::getnewibc(inmap));
      }
	}
}
