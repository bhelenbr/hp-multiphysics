/*
 *  getnewibc.cpp
 *  tri_hp
 *
 *  Created by Erik Durkish on 2/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swe.h"

namespace ibc_swe {

   class rossby : public init_bdry_cndtn {
      private:
         FLT amplitude;
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            FLT A = 0.771*amplitude*amplitude;
            FLT phi = A*pow(cosh(amplitude*x(0)),-2);
            FLT dphidx = -2*amplitude*tanh(amplitude*x(0))*phi;
            FLT h = 1.0 +phi*(3.+6.*x(1)*x(1))/4.*exp(-x(1)*x(1)/2);
            switch(n) {
               case(0):
                  return(h*phi*(-9.+6.*x(1)*x(1))/4.*exp(-x(1)*x(1)/2.));
               case(1):
                  return(h*dphidx*2*x(1)*exp(-x(1)*x(1)/2.));
					case(2):
                  return(h +0.00*exp(-x(1)*x(1)));

            }
            return(0.0);
         }
         
         void input(input_map &in_map,std::string idnty) {
            
            if (!in_map.get(idnty + "_rossby_amp",amplitude)) in_map.getwdefault("rossby_amp",amplitude,0.395);
         }
   };
	
   class ibc_type {
      public:
         const static int ntypes = 1;
         enum ids {rossby};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char ibc_swe::ibc_type::names[ntypes][40] = {"rossby"};
   
   class flat : public init_bdry_cndtn {
      private:
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            return(0.0);
         }
   };
   
   
   
   class bathy_type {
      public:
         const static int ntypes = 1;
         enum ids {flat};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char ibc_swe::bathy_type::names[ntypes][40] = {"flat"};

};


init_bdry_cndtn *tri_hp_swe::getnewibc(input_map& inmap) {
	std::string keyword,ibcname;
	int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = std::string(idprefix) + "_ibc";
	if (!inmap.get(keyword,ibcname)) {
		if (!inmap.get("ibc",ibcname)) {
			*sim::log << "couldn't find initial condition type" << std::endl;
		}
	}	
	type = ibc_swe::ibc_type::getid(ibcname.c_str());   
		
	switch(type) {
		case ibc_swe::ibc_type::rossby: {
			init_bdry_cndtn *temp = new ibc_swe::rossby;
			return(temp);
		}
		
		default: {
         return(tri_hp_ins::getnewibc(inmap));
      }
	}
}

init_bdry_cndtn *tri_hp_swe::getnewbathy(input_map& inmap) {
	std::string keyword,ibcname;
	int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = std::string(idprefix) + "_bathy";
	if (!inmap.get(keyword,ibcname)) {
		if (!inmap.get("bathy",ibcname)) {
			*sim::log << "couldn't find bathy  type" << std::endl;
		}
	}
	
	type = ibc_swe::bathy_type::getid(ibcname.c_str());
		
	switch(type) {
		case ibc_swe::bathy_type::flat: {
			init_bdry_cndtn *temp = new ibc_swe::flat;
			return(temp);
		}
      default: {
         return(tri_hp::getnewibc(inmap));
      }
	}
   
}


