/*
 *  getnewibc.cpp
 *  tet_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cns.h"
#include "bdry_cns.h"
#include <tet_boundary.h>

namespace ibc_cns {

	class freestream : public init_bdry_cndtn {
		private:
		FLT angle1,angle2, speed,perturb_amp,gamma;
		TinyVector<FLT,tet_mesh::ND> vel;

		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
				FLT amp = (time > 0.0 ? 0.0 : perturb_amp); 
				amp = perturb_amp;
				
				for (int i=0; i<3; ++i) 
					amp *= 4.0*x(i)*(1.0-x(i))*(sin(2.0*M_PI*x(i))+sin(16*M_PI*x(i)));

				switch(n) {
					case(0):
						return(1.0/gamma+amp);
					case(1):
						return(vel(0)+amp);
					case(2):
						return(vel(1)+amp);
					case(3):
						return(vel(2)+amp);
					case(4):
						return(1.0/gamma+amp);
				}
				return(0.0);
			}

			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_flowspeed";
				if (!blockdata.get(keyword,speed)) 
					blockdata.getwdefault("flowspeed",speed,0.1);

				keyword = idnty +"_flowangle1";
				if (!blockdata.get(keyword,angle1)) 
					blockdata.getwdefault("flowangle1",angle1,0.0);  

				keyword = idnty +"_flowangle2";
				if (!blockdata.get(keyword,angle2)) 
					blockdata.getwdefault("flowangle2",angle2,90.0); 
				
				keyword = idnty +"_perturb_amplitude";
				if (!blockdata.get(keyword,perturb_amp)) 
					blockdata.getwdefault("perturb_amplitude",perturb_amp,0.0); 

				keyword = idnty +"_gamma";
				if (!blockdata.get(keyword,gamma))
					blockdata.getwdefault("gamma",gamma,1.4);
				
				angle1 *= M_PI/180.0;
				angle2 *= M_PI/180.0;
				
				/* spherical coordinates */
				vel(0) = speed*cos(angle1)*sin(angle2);
				vel(1) = speed*sin(angle1)*sin(angle2);
				vel(2) = speed*cos(angle2);
			}
	};

	class sphere : public init_bdry_cndtn {
	private:
		FLT speed,angle1,angle2,inner,outer,RT,pr,gamma;
		TinyVector<FLT,tet_mesh::ND> vel;
		
	public:
		FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
			FLT r;
			
			r = sqrt(x(0)*x(0) +x(1)*x(1)+x(2)*x(2));
			switch(n) {
				case(0):
					return(pr);
				case(1):case(2):case(3): 
					if (r < inner) 
						return(0.0);
					else if (r < outer)
						return(vel(n-1)*0.5*(1.-cos(M_PI*(r-inner)/(outer-inner))));
					else
						return(vel(n-1));		
				case(4):
					return(RT);
			}
			return(0.0);
		}
		
		void input(input_map &blockdata,std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
			
			keyword = idnty +"_flowspeed";
			if (!blockdata.get(keyword,speed)) 
				blockdata.getwdefault("flowspeed",speed,0.1);
			
			keyword = idnty +"_angle1";
			if (!blockdata.get(keyword,angle1)) 
				blockdata.getwdefault("angle1",angle1,0.0);
			angle1 *= M_PI/180.0;
			
			keyword = idnty +"_angle2";
			if (!blockdata.get(keyword,angle2)) 
				blockdata.getwdefault("angle2",angle2,90.0);
			angle2 *= M_PI/180.0;
						
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
			
			/* spherical coordinates */
			vel(0) = speed*cos(angle1)*sin(angle2);
			vel(1) = speed*sin(angle1)*sin(angle2);
			vel(2) = speed*cos(angle2);
		}
	};
	
	
	class ringleb : public init_bdry_cndtn {
		
	private:
		FLT angle,xshift,yshift;
		FLT gam; 
		double q0;  // Initial guess for q in iteration //
		double theta0,theta1; // Initial guess for theta in iteration //
		double scale; // Scale point before calculating solution //
		double shift[2]; // Spatial shift of point before calculating solution //
		
		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
				Array<double,1> val(5);
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
					blockdata.getwdefault("theta0",theta0,90.0);
				theta0 *= M_PI/180.0;
				
				keyword = idnty +"_theta1";
				if (!blockdata.get(keyword,theta1)) 
					blockdata.getwdefault("theta1",theta1,90.0);
				theta1 *= M_PI/180.0;
				
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
		
		bool eval(TinyVector<FLT,tet_mesh::ND> x, double t, Array<double,1> &val) const {
			
			double q = q0;
			double theta = theta0;
			double dtheta = -1.0e-4;
			double dq = 1.0e-4;
			double jac[2][2], deti, delta[2]; 
			bool conservative = false;
			
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
				val(1) = r*q*cos(theta)*sin(theta1);
				val(2) = r*q*sin(theta)*sin(theta1);
				val(3) = r*q*cos(theta1);
				val(4) = r*(c*c/(gam*(gam-1.)) + 0.5*q*q);
			}
			else {
				val(0) = r*c*c/gam; 
				val(1) = q*cos(theta)*sin(theta1);
				val(2) = q*sin(theta)*sin(theta1);
				val(3) = q*cos(theta1);
				val(4) = val(0)/r; 
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
	
	class hydrostatic : public init_bdry_cndtn {
	private:
		FLT gamma,pressure,temperature,shiftz,gravity;
		
	public:
		FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {

			switch(n) {
				case(0):
					return(pressure*exp(gravity*(x(2)-shiftz)/temperature));
				case(1):case(2):case(3):
					return(0.0);
				case(4):
					return(temperature);
			}
			return(0.0);
		}
		
		void input(input_map &blockdata,std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
			
			keyword = idnty +"_gamma";
			if (!blockdata.get(keyword,gamma))
				blockdata.getwdefault("gamma",gamma,1.4);
			
			keyword = idnty +"_pressure";
			if (!blockdata.get(keyword,pressure)) 
				blockdata.getwdefault("pressure",pressure,1.0/gamma);
			
			keyword = idnty +"_temperature";
			if (!blockdata.get(keyword,temperature)) 
				blockdata.getwdefault("temperature",temperature,1.0/gamma);  
			
			keyword = idnty +"_shiftz";
			if (!blockdata.get(keyword,shiftz)) 
				blockdata.getwdefault("shiftz",shiftz,0.0); 
			
			keyword = idnty +"_gravity";
			if (!blockdata.get(keyword,gravity)) 
				blockdata.getwdefault("gravity",gravity,-9.81); 
			
		}
	};
	
	
	class spinninglid : public init_bdry_cndtn {

	private:
		FLT omega,epsilon,radius,height,gamma;
	public:
		FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
				
			FLT A = -omega*radius*(-radius*radius*radius+3.0*radius*radius*epsilon-3.0*radius*epsilon*epsilon+epsilon*epsilon*epsilon)/(epsilon*epsilon*epsilon);
			FLT B = omega*(-3.0*radius*radius*radius+6.0*radius*radius*epsilon-3.0*radius*epsilon*epsilon+epsilon*epsilon*epsilon)/(epsilon*epsilon*epsilon);
			FLT C = -3.0*(-radius+epsilon)*radius*omega/(epsilon*epsilon*epsilon);
			FLT D = -radius*omega/(epsilon*epsilon*epsilon);
			
			FLT r = sqrt(x(1)*x(1)+x(2)*x(2));		
			FLT vtheta = A+B*r+C*r*r+D*r*r*r;
			
			switch(n) {
				case(0):
					return(1.0/gamma);
				case(1):
					return(0.0);
				case(2):
					if(x(0)<height) {
						return(0.0);
					}
					else {
						if(r<radius-epsilon) {
							return(x(2)*omega);
						}
						else {
							return(vtheta*x(2)/r);
						}
					}
				case(3):
					if(x(0)<height) {
						return(0.0);
					}
					else {
						if(r<radius-epsilon) {
							return(-x(1)*omega);
						}
						else {
							return(-vtheta*x(1)/r);
						}
					}				
				case(4):
					return(1.0/gamma);
			}
			return(0.0);
		}
		
		void input(input_map &blockdata,std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
			
			keyword = idnty +"_rotationalspeed";
			if (!blockdata.get(keyword,omega)) 
				blockdata.getwdefault("rotationalspeed",omega,0.1); 
			
			keyword = idnty +"_offset";
			if (!blockdata.get(keyword,epsilon)) 
				blockdata.getwdefault("offset",epsilon,0.05); 
			
			keyword = idnty +"_height";
			if (!blockdata.get(keyword,height)) 
				blockdata.getwdefault("height",height,2.5); 
			
			keyword = idnty +"_radius";
			if (!blockdata.get(keyword,radius)) 
				blockdata.getwdefault("radius",radius,1.0); 
			
			keyword = idnty +"_gamma";
			if (!blockdata.get(keyword,gamma)) 
				blockdata.getwdefault("gamma",gamma,1.4);
			
		}
	};
	
	
	class ibc_type {
		public:
			const static int ntypes = 5;
			enum ids {freestream,sphere,ringleb,hydrostatic,spinninglid};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
						if (!strcmp(nin,names[i])) return(i);
				return(-1);
		}
	};
	const char ibc_type::names[ntypes][40] = {"freestream","sphere","ringleb","hydrostatic","spinninglid"};
	
	
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
	

}

tet_hp_helper *tet_hp_cns::getnewhelper(input_map& inmap) {
	std::string keyword,movername;
	int type;
	
	/* FIND INITIAL CONDITION TYPE */
	keyword = std::string(gbl->idprefix) + "_helper";
	if (!inmap.get(keyword,movername)) {
		if (!inmap.get("helper",movername)) {
			type = -1;
		}
	}
	
	type = ibc_cns::helper_type::getid(movername.c_str());
	
	switch(type) {
//		case ibc_cns::helper_type::parameter_changer: {
//			tet_hp_helper *temp = new ibc_cns::parameter_changer(*this);
//			return(temp);
//		}

		default: {
			return(tet_hp::getnewhelper(inmap));
		}
	}
}


init_bdry_cndtn *tet_hp_cns::getnewibc(std::string suffix, input_map& inmap) {
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
	type = ibc_cns::ibc_type::getid(ibcname.c_str());

		
	switch(type) {
		case ibc_cns::ibc_type::freestream: {
			temp = new ibc_cns::freestream;
			break;
		}
		case ibc_cns::ibc_type::sphere: {
			temp = new ibc_cns::sphere;
			break;
		}
		case ibc_cns::ibc_type::ringleb: {
			temp = new ibc_cns::ringleb;
			break;
		}
		case ibc_cns::ibc_type::hydrostatic: {
			temp = new ibc_cns::hydrostatic;
			break;
		}	
		case ibc_cns::ibc_type::spinninglid: {
			temp = new ibc_cns::spinninglid;
			break;
		}	
		default: {
			return(tet_hp::getnewibc(suffix,inmap));
		}
	}
	temp->input(inmap,keyword);
	return(temp);
}






