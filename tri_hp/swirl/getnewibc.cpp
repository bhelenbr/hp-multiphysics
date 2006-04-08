/*
 *  getnewibc.cpp
 *  tri_hp
 *
 *  Created by Erik Durkish on 2/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swirl.h"

namespace ibc_swirl {

   class spinning : public init_bdry_cndtn {
      private:
         FLT omega;
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            switch(n) {
               case(0):
                  return(0);
               case(1):
                  return(0);
					case(2):
						return(x(0)*omega);
					case(3):
						return(0.99*x(0)*x(0)*omega*omega/2);
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
				
				keyword = idnty +".rotationalspeed";
            if (!blockdata.get(keyword,omega)) 
               blockdata.getwdefault("rotationalspeed",omega,1.0); 
         }
   };

   class ibc_type {
      public:
         const static int ntypes = 1;
         enum ids {spinning};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char ibc_swirl::ibc_type::names[ntypes][40] = {"spinning"};

}


init_bdry_cndtn *tri_hp_swirl::getnewibc(input_map& inmap) {
	std::string keyword,ibcname;
	int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = std::string(idprefix) + ".ibc";
	if (!inmap.get(keyword,ibcname)) {
		if (!inmap.get("ibc",ibcname)) {
			*sim::log << "couldn't find initial condition type" << std::endl;
		}
	}
	
	type = ibc_swirl::ibc_type::getid(ibcname.c_str());
		
	switch(type) {
		case ibc_swirl::ibc_type::spinning: {
			init_bdry_cndtn *temp = new ibc_swirl::spinning;
			return(temp);
		}

		default: {
			return(0);
		}
	}
}
