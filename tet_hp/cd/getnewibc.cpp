		
/*
 *  getnewibc_ins.cpp
 *  tet_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cd.h"
#include <spline.h>

namespace ibc_cd {
	
	class zero_src : public init_bdry_cndtn {
		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x,FLT time) {
				return(0.0);
			}
			void init(input_map &inmap, std::string idnty) {}
	};
	
	class unsteady_spline : public init_bdry_cndtn {
		spline<1> time_history;
		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x,FLT time) {
				TinyVector<FLT,1> source;
				time_history.interpolate(time, source);
				return(source(0));
			}
			void init(input_map &inmap, std::string idnty) {
				std::string filename;
				if (!inmap.get(idnty+"_filename",filename)) {
					std::cerr << "couldn't find unsteady_spline filename" << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);
				}
				time_history.read(filename);
			}
	};
	
	class unsteady_spline_regions : public init_bdry_cndtn {
		int nregion;
		Array<FLT,3>  bbox;
		Array<spline<1>,1> time_history;
		
	public:
		FLT f(int n, TinyVector<FLT,tet_mesh::ND> x,FLT time) {
			TinyVector<FLT,1> source;
			
			source(0) = 0.0;
			for(int r=0;r<nregion;++r) {
				int sign_sum = 0;
				for (int n=0;n<tet_mesh::ND;++n) {
					sign_sum += ((x(n)-bbox(r,0,n))*(x(n)-bbox(r,1,n)) < 0.0 ? -1 : 1);
				}
				if (sign_sum == -3) {
					time_history(r).interpolate(time, source);
				}
			}
			return(source(0));
		}
		
		void init(input_map &inmap, std::string idnty) {
			
			if (!inmap.get(idnty+"_nregions",nregion)){
				std::cerr << "couldn't find " << idnty+"nregions" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
			}
			bbox.resize(nregion,2,tet_mesh::ND);
			time_history.resize(nregion);
			
			for(int r=0;r<nregion;++r) {
				std::string filename;
				stringstream nstr;
				nstr << r << std::flush;
				if (!inmap.get(idnty+"_region" +nstr.str() +"_filename",filename)) {
					if (!inmap.get(idnty +"_filename",filename)) {
						std::cerr << "couldn't find unsteady_spline filename" << std::endl;
						sim::abort(__LINE__,__FILE__,&std::cerr);
					}
				}
				time_history(r).read(filename);
				
				string bbox_string;
				inmap.getline(idnty+"_region" +nstr.str() +"_bbox0",bbox_string);
				istringstream bbox_stream(bbox_string);
				for (int n=0;n<tet_mesh::ND;++n)
					bbox_stream >> bbox(r,0,n);
	
				inmap.getline(idnty+"_region" +nstr.str() +"_bbox1",bbox_string);
				bbox_stream.clear();
				bbox_stream.str(bbox_string);
				for (int n=0;n<tet_mesh::ND;++n) {
					bbox_stream >> bbox(r,1,n);
				}
			}
		}
	};
	
	

	class ibc_type {
		public:
			const static int ntypes = 3;
			enum ids {zero,unsteady_spline,unsteady_spline_regions};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
			}
	};
	const char ibc_type::names[ntypes][40] = {"zero","unsteady_spline","unsteady_spline_regions"};
}


init_bdry_cndtn *tet_hp_cd::getnewibc(std::string ibcname) {
	init_bdry_cndtn *temp;
	int type;

	type = ibc_cd::ibc_type::getid(ibcname.c_str());
	switch(type) {
		case(ibc_cd::ibc_type::zero): {
			temp = new ibc_cd::zero_src;
			break;
		}
		case(ibc_cd::ibc_type::unsteady_spline): {
			temp = new ibc_cd::unsteady_spline;
			break;
		}
		case(ibc_cd::ibc_type::unsteady_spline_regions): {
			temp = new ibc_cd::unsteady_spline_regions;
			break;
		}
		default: {
			return(tet_hp::getnewibc(ibcname));
		}
	}
	return(temp);
}
