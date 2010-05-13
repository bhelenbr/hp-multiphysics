/*
 *  getnewibc_cns.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cns_explicit.h"
#include "bdry_cns_explicit.h"
#include <tri_boundary.h>

namespace ibc_cns_explicit {

	class freestream : public init_bdry_cndtn {
		private:
			FLT alpha, speed,perturb_amp,gamma,pr,RT;
			TinyVector<FLT,tri_mesh::ND> vel;

		public:

			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				FLT amp = (time > 0.0 ? 0.0 : perturb_amp); 
			
				switch(n) {
					case(0):
						return(pr/RT);
					case(1):
						return(pr/RT*vel(0) +amp*x(0)*(1.0-x(0)));
					case(2):
						return(pr/RT*vel(1));
					case(3):
						return(pr/(gamma-1.0)+0.5*pr/RT*(vel(0)*vel(0)+vel(1)*vel(1)));
				}
				return(0.0);
			}

			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_flowspeed";
				if (!blockdata.get(keyword,speed)) 
					blockdata.getwdefault("flowspeed",speed,0.1);

				keyword = idnty +"_flowangle";
				if (!blockdata.get(keyword,alpha)) 
					blockdata.getwdefault("flowangle",alpha,0.0);  

				keyword = idnty +"_perturb_amplitude";
				if (!blockdata.get(keyword,perturb_amp)) 
					blockdata.getwdefault("perturb_amplitude",perturb_amp,0.0); 

				keyword = idnty +"_gamma";
				if (!blockdata.get(keyword,gamma))
					blockdata.getwdefault("gamma",gamma,1.403);
				
				keyword = idnty +"_RT";
				if (!blockdata.get(keyword,RT)) 
					blockdata.getwdefault("RT",RT,1.0/gamma);
				
				keyword = idnty +"_pressure";
				if (!blockdata.get(keyword,pr)) 
					blockdata.getwdefault("pressure",pr,1.0/gamma);
				
				alpha *= M_PI/180.0;
				vel(0) = speed*cos(alpha);
				vel(1) = speed*sin(alpha);
			}
    };

	class sphere : public init_bdry_cndtn {
		private:
			FLT speed,angle,inner,outer,RT,pr,gamma;
			TinyVector<FLT,tri_mesh::ND> vel;

		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				FLT r = sqrt(x(0)*x(0) +x(1)*x(1));
								
				TinyVector<FLT,2> v;
				if (r < inner) v = 0.0;
				else if (r < outer) v = vel*0.5*(1.-cos(M_PI*(r-inner)/(outer-inner)));
				else v = vel;
				
				switch(n) {
					case(0):
						return(pr/RT);
					case(1):case(2): 
						return(pr/RT*v(n-1));
					case(3):
						return(pr/(gamma-1.0)+0.5*pr/RT*(v(0)*v(0)+v(1)*v(1)));
				}
				return(0.0);
			}

			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_flowspeed";
				if (!blockdata.get(keyword,speed)) 
					blockdata.getwdefault("flowspeed",speed,0.1);

				keyword = idnty +"_angle";
				if (!blockdata.get(keyword,angle)) 
					blockdata.getwdefault("angle",angle,0.0);
				angle *= M_PI/180.0;

				keyword = idnty +"_inner_radius";
				if (!blockdata.get(keyword,inner)) 
					blockdata.getwdefault("inner_radius",inner,1.1);

				keyword = idnty +"_outer_radius";
				if (!blockdata.get(keyword,outer)) 
					blockdata.getwdefault("outer_radius",outer,2.1);
					
				keyword = idnty +"_gamma";
				if (!blockdata.get(keyword,gamma))
					blockdata.getwdefault("gamma",gamma,1.403);

				keyword = idnty +"_RT";
				if (!blockdata.get(keyword,RT)) 
					blockdata.getwdefault("RT",RT,1.0/gamma);
				
				keyword = idnty +"_pressure";
				if (!blockdata.get(keyword,pr)) 
					blockdata.getwdefault("pressure",pr,1.0/gamma);
				
				vel(0) = speed*cos(angle);
				vel(1) = speed*sin(angle);
			}
	};

	
	class ringleb : public init_bdry_cndtn {
	
		private:
			FLT angle,xshift,yshift;
			FLT gam; // temp figure out how to load gbl->gamma
			double q0;  // Initial guess for q in iteration //
			double theta0; // Initial guess for theta in iteration //
			double scale; // Scale point before calculating solution //
			double shift[2]; // Spatial shift of point before calculating solution //
		
		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				Array<double,1> val(4);
				bool error;
			
				error = eval(x,time,val);			
				
				return(val(n));
			}
		
			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;
				
				keyword = idnty +"_q0";
				if (!blockdata.get(keyword,q0)) 
					blockdata.getwdefault("q0",q0,0.02);
				
				keyword = idnty +"_theta0";
				if (!blockdata.get(keyword,theta0)) 
					blockdata.getwdefault("theta0",theta0,M_PI/2.);
				angle *= M_PI/180.0;
				
				keyword = idnty +"_scale";
				if (!blockdata.get(keyword,scale)) 
					blockdata.getwdefault("scale",scale,20.0);
				
				keyword = idnty +"_xshift";
				if (!blockdata.get(keyword,xshift)) 
					blockdata.getwdefault("xshift",xshift,-40.0);
				
				keyword = idnty +"_yshift";
				if (!blockdata.get(keyword,yshift)) 
					blockdata.getwdefault("yshift",yshift,0.0);
				
				keyword = idnty +"_gamma";
				if (!blockdata.get(keyword,gam)) 
					blockdata.getwdefault("gamma",gam,1.403);
				
				shift[0] = xshift;
				shift[1] = yshift;

				
			}
		
			bool eval(TinyVector<FLT,tri_mesh::ND> x, double t, Array<double,1> &val) const {
				
				double q = q0;
				double theta = theta0;
				double dtheta = -1.0e-4;
				double dq = 1.0e-4;
				double jac[2][2], deti, delta[2]; 
				bool conservative = true;
				
				double xyz[3];
				xyz[0] = scale*x(0) +shift[0];
				xyz[1] = scale*x(1) +shift[1];
				
				if (fabs(xyz[1]) < 1.0e-8)
					theta = M_PI/2.0;
				
				double error[2] = {1.0, 1.0};
				int niter = 0;
				while (fabs(error[0])+fabs(error[1]) > 1.0e-15 && niter++ < 30) {
					calculateError(xyz, error, q, theta);
					calculateError(xyz, delta, q+dq, theta);
					jac[0][0] = (delta[0]-error[0])/dq;
					jac[1][0] = (delta[1]-error[1])/dq;
					calculateError(xyz, delta, q, theta+dtheta);
					jac[0][1] = (delta[0]-error[0])/dtheta;
					jac[1][1] = (delta[1]-error[1])/dtheta;
					deti = 1./(jac[0][0]*jac[1][1] -jac[0][1]*jac[1][0]);
					
					if (fabs(xyz[1]) > 1.0e-6) {
						delta[0] = deti*(error[0]*jac[1][1] -error[1]*jac[0][1]);
						delta[1] = deti*(error[1]*jac[0][0] -error[0]*jac[1][0]);
						q -= delta[0];
						theta -= delta[1];
					}
					else {
						q -= error[0]/jac[0][0];
						theta = M_PI/2.;
					}
					
					if (fabs(theta) > M_PI/2.) {
						theta = M_PI/2. -(theta-M_PI/2.);
					}
				}
				if (niter > 9999 || !(q >= 0.0)) {  
					std::cout << "RINGLEB NOT_CONVERGED: " << xyz[0] << ' ' << xyz[1] << niter << ' ' << q << ' ' << theta << std::endl;
					return false;
				}
				
				double c = sqrt(1. - (gam-1.)/2.*q*q);
				double r = pow(c,(2./(gam-1.)));
				
				if (conservative) {
					val(0) = r;
					val(1) = r*q*cos(theta);
					val(2) = r*q*sin(theta);
					val(3) = r*(c*c/(gam*(gam-1.)) + 0.5*q*q);
				}
				else {
					val(0) = r*c*c/gam; 
					val(1) = q*cos(theta);
					val(2) = q*sin(theta);
					val(3) = val(0)/r; 
				}
				
				return true;
				
			}
		
			void calculateError(const double xyz[3], double error[2], double q, double theta) const {
				double psi,c,J,r,xp,yp;
				
				psi = 1./q*sin(theta);
				c = sqrt(1. - (gam-1.)/2.*q*q);
				J = 1./c +1./(3*c*c*c) +1./(5.*c*c*c*c*c) -1./2.*log((1.+c)/(1.-c));
				r = pow(c,(2./(gam-1.)));
				xp = 1./(2.*r)*(1./(q*q)-2.*psi*psi) +J/2.;
				
				error[0] = xyz[0]-xp;
				
				if (fabs(xyz[1]) > 1.0e-6) {
					yp = xyz[1]/fabs(xyz[1])*psi/(q*r)*sqrt(1.-psi*psi*q*q);
					error[1] = xyz[1]-yp;
				}
				else {
					error[1] = 0.0;
				}
				
				return;
			}
		
			
	};
	
	
	class helper_type {
		public:
			const static int ntypes = 4;
			enum ids {translating_drop,parameter_changer,unsteady_body_force,force_coupling};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
					if (!strcmp(nin,names[i])) return(i);
				return(-1);
			}
	};
	const char helper_type::names[ntypes][40] = {"translating_drop","parameter_changer","unsteady_body_force","force_coupling"};

	class ibc_type {
		public:
			const static int ntypes = 7;
			enum ids {freestream,sphere,accelerating,impinge,stokes_drop_gas,stokes_drop_liquid,ringleb};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
					if (!strcmp(nin,names[i])) return(i);
				return(-1);
			}
	};
	const char ibc_type::names[ntypes][40] = {"freestream","sphere","accelerating","impinge","stokes_drop_gas","stokes_drop_liquid","ringleb"};

}


init_bdry_cndtn *tri_hp_cns_explicit::getnewibc(std::string suffix, input_map& inmap) {
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
	type = ibc_cns_explicit::ibc_type::getid(ibcname.c_str());


	switch(type) {
		case ibc_cns_explicit::ibc_type::freestream: {
			temp = new ibc_cns_explicit::freestream;
			break;
		}
		case ibc_cns_explicit::ibc_type::sphere: {
			temp = new ibc_cns_explicit::sphere;
			break;
		}
		case ibc_cns_explicit::ibc_type::ringleb: {
			temp = new ibc_cns_explicit::ringleb;
			break;
		}
//		case ibc_cns_explicit::ibc_type::accelerating: {
//			temp = new ibc_cns_explicit::accelerating;
//			break;
//		}
//		case ibc_cns_explicit::ibc_type::impinge: {
//			temp = new ibc_cns_explicit::impinge;
//			break;
//		}
//		case ibc_cns_explicit::ibc_type::stokes_drop_gas: {
//			temp = new ibc_cns_explicit::stokes_drop_gas;
//			break;
//		}
//		case ibc_cns_explicit::ibc_type::stokes_drop_liquid: {
//			temp = new ibc_cns_explicit::stokes_drop_liquid;
//			break;
//		}
		default: {
			return(tri_hp::getnewibc(suffix,inmap));
		}
	}
	temp->input(inmap,keyword);
	return(temp);
}

tri_hp_helper *tri_hp_cns_explicit::getnewhelper(input_map& inmap) {
	std::string keyword,movername;
	int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = std::string(gbl->idprefix) + "_helper";
	if (!inmap.get(keyword,movername)) {
		if (!inmap.get("helper",movername)) {
			type = -1;
		}
	}

	type = ibc_cns_explicit::helper_type::getid(movername.c_str());

	switch(type) {
//		case ibc_cns_explicit::helper_type::translating_drop: {
//			tri_hp_helper *temp = new ibc_cns_explicit::translating_drop(*this);
//			return(temp);
//		}
//		case ibc_cns_explicit::helper_type::parameter_changer: {
//			tri_hp_helper *temp = new ibc_cns_explicit::parameter_changer(*this);
//			return(temp);
//		}
//		case ibc_cns_explicit::helper_type::unsteady_body_force: {
//			tri_hp_helper *temp = new ibc_cns_explicit::unsteady_body_force(*this);
//			return(temp);
//		}
//		case ibc_cns_explicit::helper_type::force_coupling: {
//			tri_hp_helper *temp = new ibc_cns_explicit::force_coupling(*this);
//			return(temp);
//		}
		default: {
			return(tri_hp::getnewhelper(inmap));
		}
	}
}
