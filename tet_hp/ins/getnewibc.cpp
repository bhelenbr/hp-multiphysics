/*
 *  getnewibc_ins.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_ins.h"
#include "bdry_ins.h"
#include <tet_boundary.h>

namespace ibc_ins {

	class freestream : public init_bdry_cndtn {
		private:
			FLT alpha, speed,perturb_amp;

		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
				FLT amp = (time > 0.0 ? 0.0 : perturb_amp); 
				switch(n) {
					case(0):
						return(0*speed*cos(alpha) +amp*x(0)*(1.0-x(0)));
					case(1):
						return(speed*sin(alpha));
					case(2):
						return(16/0.08333333333/0.08333333333/0.1875/0.1875*(0.1875/2.0+x(0)-0.375)*(0.1875/2.0-x(0)+0.375)*(0.08333333333/2.0+x(1)-0.104166666665)*(0.08333333333/2.0-x(1)+0.104166666665));
					case(3):
						return(0.0);
				}
				return(0.0);
			}

			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_flowspeed";
				if (!blockdata.get(keyword,speed)) 
						blockdata.getwdefault("flowspeed",speed,1.0);

				keyword = idnty +"_flowangle";
				if (!blockdata.get(keyword,alpha)) 
						blockdata.getwdefault("flowangle",alpha,0.0);  

				keyword = idnty +"_perturb_amplitude";
				if (!blockdata.get(keyword,perturb_amp)) 
						blockdata.getwdefault("perturb_amplitude",perturb_amp,0.0); 

				alpha *= M_PI/180.0;
			}
	};
	


	class ibc_type {
		public:
			const static int ntypes = 1;
			enum ids {freestream};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
						if (!strcmp(nin,names[i])) return(i);
				return(-1);
		}
	};
	const char ibc_type::names[ntypes][40] = {"freestream"};
	
	
	class parameter_changer : public tet_hp_helper {
	protected:
		tet_hp_ins &x;
		//bdry_ins::surface *surf;
		FLT delta_rho, rho_factor;
		FLT delta_mu, mu_factor;
		FLT delta_g, delta_g_factor;
		FLT delta_rho2, rho2_factor;
		FLT delta_mu2, mu2_factor;
		FLT delta_sigma, sigma_factor;
		int interval;
		
	public:
		parameter_changer(tet_hp_ins& xin) : tet_hp_helper(xin), x(xin) {
//			int bnum;
//			
//			for(bnum=0;bnum<x.nebd;++bnum) 
//				if (surf = dynamic_cast<bdry_ins::surface *>(x.hp_ebdry(bnum))) break;
//			
//			if (bnum > x.nebd -1) surf = 0;
		}
		void init(input_map& input, std::string idnty) {
			std::string keyword, val;
			
			keyword = idnty + "_delta_rho";
			if (!input.get(keyword,delta_rho)) {
				input.getwdefault("delta_rho",delta_rho,0.0);
			}
			
			keyword = idnty + "_rho_factor";
			if (!input.get(keyword,rho_factor)) {
				input.getwdefault("rho_factor",rho_factor,1.0);
			}
			
			keyword = idnty + "_delta_mu";
			if (!input.get(keyword,delta_mu)) {
				input.getwdefault("delta_mu",delta_mu,0.0);
			}
			
			keyword = idnty + "_mu_factor";
			if (!input.get(keyword,mu_factor)) {
				input.getwdefault("mu_factor",mu_factor,1.0);
			}
			
			keyword = "delta_g";
			input.getwdefault(keyword,delta_g,0.0);
			input[keyword] = "0.0"; // SO ONLY ONE BLOCK PER PROCESSOR CHANGES THIS
			
			keyword = "delta_g_factor";
			input.getwdefault(keyword,delta_g_factor,1.0);
			input[keyword] = "1.0"; // SO ONLY ONE BLOCK PER PROCESSOR CHANGES THIS
			
			input.getwdefault("parameter_interval",interval,1);
			
//			if (surf) {
//				std::string surfidnty = surf->base.idprefix;
//				
//				keyword = surfidnty + "_delta_sigma";
//				input.getwdefault(keyword,delta_sigma,0.0);
//				
//				keyword = surfidnty + "_sigma_factor";
//				input.getwdefault(keyword,sigma_factor,1.0);
//				
//				keyword = surfidnty + "_matching_block";
//				if (!input.get(keyword,val)) {
//					delta_rho2 = 0.0;
//					rho2_factor = 1.0;
//					delta_mu2 = 0.0;
//					mu2_factor = 1.0;
//				}
//				else {                         
//					keyword = val + "_delta_rho";                    
//					if (!input.get(keyword,delta_rho2)) {
//						input.getwdefault("delta_rho",delta_rho2,0.0);
//					}
//					
//					keyword = val + "_rho_factor";
//					if (!input.get(keyword,rho2_factor)) {
//						input.getwdefault("rho_factor",rho2_factor,1.0);
//					}
//					
//					keyword = val + "_delta_mu";
//					if (!input.get(keyword,delta_mu2)) {
//						input.getwdefault("delta_mu",delta_mu2,0.0);
//					}
//					
//					keyword = val + "_mu_factor";
//					if (!input.get(keyword,mu2_factor)) {
//						input.getwdefault("mu_factor",mu2_factor,1.0);
//					}
//				}
//			}
		}
		tet_hp_helper* create(tet_hp& xin) { return new parameter_changer(dynamic_cast<tet_hp_ins&>(xin)); }
		
		
		
		void tadvance() {
			if (!x.coarse_level) {
				if ( (x.gbl->tstep % interval) +x.gbl->substep == 0) {
					
					x.gbl->rho += delta_rho;
					x.gbl->rho *= rho_factor;
					
					x.gbl->mu  += delta_mu;
					x.gbl->mu  *= mu_factor;
					
					x.gbl->g += delta_g;
					x.gbl->g *= delta_g_factor;
					
					*x.gbl->log << "new density, viscosity, and gravity are " << x.gbl->rho << ' ' << x.gbl->mu << ' ' << x.gbl->g << std::endl;
					
//					
//					if (surf) {
//						surf->gbl->rho2 += delta_rho2;
//						surf->gbl->rho2 *= rho2_factor;
//						
//						surf->gbl->mu2  += delta_mu2;
//						surf->gbl->mu2  *= mu2_factor;
//						
//						surf->gbl->sigma  += delta_sigma;
//						surf->gbl->sigma  *= sigma_factor;
//						
//						*x.gbl->log << "matching block density, viscosity, and surface tension are " << surf->gbl->rho2 << ' ' << surf->gbl->mu2 << ' ' << surf->gbl->sigma << std::endl;
//					}
				}
			}
			return;
		}
	};
	
	class helper_type {
	public:
		const static int ntypes = 1;
		enum ids {parameter_changer};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
	};
	const char helper_type::names[ntypes][40] = {"parameter_changer"};
	

}

tet_hp_helper *tet_hp_ins::getnewhelper(input_map& inmap) {
	std::string keyword,movername;
	int type;
	
	/* FIND INITIAL CONDITION TYPE */
	keyword = std::string(gbl->idprefix) + "_helper";
	if (!inmap.get(keyword,movername)) {
		if (!inmap.get("helper",movername)) {
			type = -1;
		}
	}
	
	type = ibc_ins::helper_type::getid(movername.c_str());
	
	switch(type) {
		case ibc_ins::helper_type::parameter_changer: {
			tet_hp_helper *temp = new ibc_ins::parameter_changer(*this);
			return(temp);
		}

		default: {
			return(tet_hp::getnewhelper(inmap));
		}
	}
}


init_bdry_cndtn *tet_hp_ins::getnewibc(std::string suffix, input_map& inmap) {
	std::string keyword,ibcname;
	init_bdry_cndtn *temp;
	int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = gbl->idprefix + "_" +suffix;
	if (!inmap.get(keyword,ibcname)) {
		keyword = suffix;
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find initial condition type" << std::endl;
		}
	}
	type = ibc_ins::ibc_type::getid(ibcname.c_str());

		
	switch(type) {
		case ibc_ins::ibc_type::freestream: {
			temp = new ibc_ins::freestream;
			break;
		}
		default: {
			return(tet_hp::getnewibc(suffix,inmap));
		}
	}
	temp->input(inmap,keyword);
	return(temp);
}






