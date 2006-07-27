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
                  return(0.0);
               case(1):
                  return(0.0);
					case(2):
						return(x(0)*omega);
					//case(3):
						//return(0.995*x(0)*x(0)*omega*omega/2);
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
	
	class jet : public init_bdry_cndtn {
      private:
         FLT omega,speed;
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            switch(n) {
               case(0):
                  return(x(0)*speed/2);
               case(1):
                  return(-speed*x(1));
					case(2):
						return(x(0)*omega);
					case(3):
						return(1*(x(0)*x(0)/2*(omega*omega-speed*speed/4)-x(1)*x(1)*speed*speed/2));
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
				
				keyword = idnty +".rotationalspeed";
            if (!blockdata.get(keyword,omega)) 
               blockdata.getwdefault("rotationalspeed",omega,1.0); 
				
				keyword = idnty +".jetspeed";
            if (!blockdata.get(keyword,omega)) 
               blockdata.getwdefault("jetspeed",speed,1.0); 
         }
   };
	
	class spinninglid : public init_bdry_cndtn {
      private:
         FLT omega,epsilon;
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
				double R,H,A,B,C,D;
				R=1.0;
				H=1.75;
				A=-omega*R*(-R*R*R+3*R*R*epsilon-3*R*epsilon*epsilon+epsilon*epsilon*epsilon)/(epsilon*epsilon*epsilon);
				B=omega*(-3*R*R*R+6*R*R*epsilon-3*R*epsilon*epsilon+epsilon*epsilon*epsilon)/(epsilon*epsilon*epsilon);
				C=-3*(-R+epsilon)*R*omega/(epsilon*epsilon*epsilon);
				D=-R*omega/(epsilon*epsilon*epsilon);
            switch(n) {
               case(0):
                  return(0.0);
               case(1):
                  return(0.0);
					case(2):
						if(x(1)<H) {
							return(0.0);
						}
						else {
							if(x(0)<R-epsilon) {
								return(x(0)*omega);
							}
							else {
								return(A+B*x(0)+C*x(0)*x(0)+D*x(0)*x(0)*x(0));
							}
						}
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
				
				keyword = idnty +".rotationalspeed";
            if (!blockdata.get(keyword,omega)) 
               blockdata.getwdefault("rotationalspeed",omega,1.0); 
				
				keyword = idnty +".offset";
            if (!blockdata.get(keyword,epsilon)) 
               blockdata.getwdefault("offset",epsilon,0.01); 
         }
   };
	
	class spinninglid2 : public init_bdry_cndtn {
      private:
         FLT w,H;
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            switch(n) {
               case(0):
                  return(0.0);
               case(1):
                  return(0.0);
					case(2):
						if(x(1)<H) {
							return(0.0);
						}
						else {
							return(x(0)*w);
						}
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
				
				keyword = idnty +".rotationalspeed";
            if (!blockdata.get(keyword,w)) 
               blockdata.getwdefault("rotationalspeed",w,1.0); 
				
				keyword = idnty +".height";
            if (!blockdata.get(keyword,H)) 
               blockdata.getwdefault("height",H,2.0); 
         }
   };

	class testFunc : public init_bdry_cndtn {
			public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            switch(n) {
               case(0):
                  return(x(0)*x(0));
               case(1):
                  return(-3*x(0)*x(1));
					case(2):
						return(x(0)*x(0));
					case(3):
						return(x(0)*x(0));
			}
            return(0.0);
         }
			
		void input(input_map &blockdata,std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
		}
   };
	
	class freespin : public init_bdry_cndtn {
      private:
         FLT w,R;
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            switch(n) {
               case(0):
                  return(0.0);
               case(1):
                  return(0.0);
					case(2):
						//if(x(1)==-2.0 || x(0)==R) {
							return(x(0)*w);
						//}
						//else {
							//return(0.0);
						//}
					case(3):
						return(x(0)*x(0)*w*w/2.0);
						
            }
            return(0.0);
         }
         
         void input(input_map &blockdata,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
				
				keyword = idnty +".rotationalspeed";
            if (!blockdata.get(keyword,w)) 
               blockdata.getwdefault("rotationalspeed",w,1.0); 
				
				keyword = idnty +".radius";
            if (!blockdata.get(keyword,R)) 
               blockdata.getwdefault("radius",R,1.0); 
         }
   };
	
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
         const static int ntypes = 7;
         enum ids {spinning,jet,spinninglid,spinninglid2,testFunc,freespin,stationary};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char ibc_swirl::ibc_type::names[ntypes][40] = {"spinning","jet","spinninglid","spinninglid2","testFunc","freespin","stationary"};

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
		case ibc_swirl::ibc_type::jet: {
			init_bdry_cndtn *temp = new ibc_swirl::jet;
			return(temp);
			
		}
		case ibc_swirl::ibc_type::spinninglid: {
			init_bdry_cndtn *temp = new ibc_swirl::spinninglid;
			return(temp);
			
		}
		
		case ibc_swirl::ibc_type::spinninglid2: {
		init_bdry_cndtn *temp = new ibc_swirl::spinninglid2;
		return(temp);
			
		}
		
		case ibc_swirl::ibc_type::testFunc: {
		init_bdry_cndtn *temp = new ibc_swirl::testFunc;
		return(temp);
			
		}
		
		case ibc_swirl::ibc_type::freespin: {
		init_bdry_cndtn *temp = new ibc_swirl::freespin;
		return(temp);
			
		}
		
		case ibc_swirl::ibc_type::stationary: {
		init_bdry_cndtn *temp = new ibc_swirl::stationary;
		return(temp);
			
		}
		
		default: {
         return(tri_hp::getnewibc(inmap));
      }
	}
}
