/*
 *  getnewibc.cpp
 *  tri_hp
 *
 *  Created by Erik Durkish on 2/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swirl.h"
#include <cstring>

namespace ibc_swirl {

	class spinning : public init_bdry_cndtn {
		private:
			FLT omega;
		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
				switch(n) {
					case(0):
						return(0.0);
					case(1):
						return(0.0);
					case(2):
						return(x(0)*omega);
					case(3):
						return(x(0)*x(0)*omega*omega/2);
				}
				return(0.0);
			}

			void init(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_rotationalspeed";
				if (!blockdata.get(keyword,omega)) 
					blockdata.getwdefault("rotationalspeed",omega,1.0); 
			}
	};

	class jet : public init_bdry_cndtn {
		private:
			FLT omega,speed;
		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
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

			void init(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_rotationalspeed";
				if (!blockdata.get(keyword,omega)) 
					blockdata.getwdefault("rotationalspeed",omega,1.0); 

				keyword = idnty +"_jetspeed";
				if (!blockdata.get(keyword,omega)) 
					blockdata.getwdefault("jetspeed",speed,1.0); 
			}
	};

	class spinninglid : public init_bdry_cndtn {
		private:
			FLT omega,epsilon;
		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
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

			void init(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_rotationalspeed";
				if (!blockdata.get(keyword,omega)) 
					blockdata.getwdefault("rotationalspeed",omega,1.0); 

				keyword = idnty +"_offset";
				if (!blockdata.get(keyword,epsilon)) 
					blockdata.getwdefault("offset",epsilon,0.01); 
			}
	};

	class spinninglid2 : public init_bdry_cndtn {
		private:
			FLT w,H;
		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
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

			void init(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_rotationalspeed";
				if (!blockdata.get(keyword,w)) 
					blockdata.getwdefault("rotationalspeed",w,1.0); 

				keyword = idnty +"_height";
				if (!blockdata.get(keyword,H)) 
					blockdata.getwdefault("height",H,2.0); 
			}
	};

	class ibc_type {
		public:
			const static int ntypes = 4;
			enum ids {spinning,jet,spinninglid,spinninglid2};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
					if (!strcmp(nin,names[i])) return(i);
				return(-1);
		}
	};
	const char ibc_swirl::ibc_type::names[ntypes][40] = {"spinning","jet","spinninglid","spinninglid2"};

}


init_bdry_cndtn *tri_hp_swirl::getnewibc(std::string name) {
	std::string keyword,ibcname;
	init_bdry_cndtn *temp;
	int type;

	type = ibc_swirl::ibc_type::getid(name.c_str());
	switch(type) {
		case ibc_swirl::ibc_type::spinning: {
			temp = new ibc_swirl::spinning;
			break;

		}
		case ibc_swirl::ibc_type::jet: {
			temp = new ibc_swirl::jet;
			break;

		}
		case ibc_swirl::ibc_type::spinninglid: {
			temp = new ibc_swirl::spinninglid;
			break;
		}

		case ibc_swirl::ibc_type::spinninglid2: {
			temp = new ibc_swirl::spinninglid2;
			break;
		}

		default: {
			return(tri_hp::getnewibc(name));
		}
	}
	return(temp);
}
