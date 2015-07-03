/*
 *  getnewibc.cpp
 *  tet_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cns_explicit.h"
#include "bdry_cns_explicit.h"
#include <tet_boundary.h>

namespace ibc_cns_explicit {

	class freestream : public init_bdry_cndtn {
		private:
		FLT angle1,angle2, speed,perturb_amp,gamma;
		TinyVector<FLT,tet_mesh::ND> vel;

		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
				FLT amp = (time > 0.0 ? 0.0 : perturb_amp); 
				switch(n) {
					case(0):
						return(1.0);
					case(1):
						return(vel(0)+amp*x(0)*(1.0-x(0)));
					case(2):
						return(vel(1));
					case(3):
						return(vel(2));
					case(4):
						return(1.0/gamma/(gamma-1.0)+.5*speed*speed);
				}
				return(0.0);
			}

			void init(input_map &inmap, std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_flowspeed";
				if (!inmap.get(keyword,speed)) 
					inmap.getwdefault("flowspeed",speed,0.1);

				keyword = idnty +"_flowangle1";
				if (!inmap.get(keyword,angle1)) 
					inmap.getwdefault("flowangle1",angle1,0.0);  

				keyword = idnty +"_flowangle2";
				if (!inmap.get(keyword,angle2)) 
					inmap.getwdefault("flowangle2",angle2,90.0); 
				
				keyword = idnty +"_perturb_amplitude";
				if (!inmap.get(keyword,perturb_amp)) 
					inmap.getwdefault("perturb_amplitude",perturb_amp,0.0); 

				keyword = idnty +"_gamma";
				if (!inmap.get(keyword,gamma))
					inmap.getwdefault("gamma",gamma,1.4);
				
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
		FLT speed,angle1,angle2,inner,outer,rho,rhoE,gamma;
		TinyVector<FLT,tet_mesh::ND> vel;
		
	public:
		FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
			FLT r;
			
			r = sqrt(x(0)*x(0)+x(1)*x(1)+x(2)*x(2));
			switch(n) {
				case(0):
					return(rho);
				case(1):case(2):case(3): 
					if (r < inner) 
						return(0.0);
					else if (r < outer)
						return(vel(n-1)*0.5*(1.-cos(M_PI*(r-inner)/(outer-inner))));
					else
						return(vel(n-1));		
				case(4):
					return(rhoE);
			}
			return(0.0);
		}
		
		void init(input_map &inmap, std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
			
			keyword = idnty +"_flowspeed";
			if (!inmap.get(keyword,speed)) 
				inmap.getwdefault("flowspeed",speed,0.1);
			
			keyword = idnty +"_angle1";
			if (!inmap.get(keyword,angle1)) 
				inmap.getwdefault("angle1",angle1,0.0);
			angle1 *= M_PI/180.0;
			
			keyword = idnty +"_angle2";
			if (!inmap.get(keyword,angle2)) 
				inmap.getwdefault("angle2",angle2,90.0);
			angle2 *= M_PI/180.0;
						
			keyword = idnty +"_inner_radius";
			if (!inmap.get(keyword,inner)) 
				inmap.getwdefault("inner_radius",inner,1.1);
			
			keyword = idnty +"_outer_radius";
			if (!inmap.get(keyword,outer)) 
				inmap.getwdefault("outer_radius",outer,2.1);
			
			keyword = idnty +"_gamma";
			if (!inmap.get(keyword,gamma))
				inmap.getwdefault("gamma",gamma,1.403);
			
			keyword = idnty +"_rho";
			if (!inmap.get(keyword,rho)) 
				inmap.getwdefault("rho",rho,1.0);
			
			keyword = idnty +"_rhoE";
			if (!inmap.get(keyword,rhoE)) 
				inmap.getwdefault("rhoE",rhoE,1.0/gamma/(gamma-1.0)+0.5*speed*speed);
			
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
			
			void init(input_map &inmap, std::string idnty) {
				std::string keyword,val;
				std::istringstream data;
				
				keyword = idnty +"_q0";
				if (!inmap.get(keyword,q0)) 
					inmap.getwdefault("q0",q0,0.02);
				
				keyword = idnty +"_theta0";
				if (!inmap.get(keyword,theta0)) 
					inmap.getwdefault("theta0",theta0,90.0);
				theta0 *= M_PI/180.0;
				
				keyword = idnty +"_theta1";
				if (!inmap.get(keyword,theta1)) 
					inmap.getwdefault("theta1",theta1,90.0);
				theta1 *= M_PI/180.0;
				
				keyword = idnty +"_scale";
				if (!inmap.get(keyword,scale)) 
					inmap.getwdefault("scale",scale,20.0);
				
				keyword = idnty +"_xshift";
				if (!inmap.get(keyword,xshift)) 
					inmap.getwdefault("xshift",xshift,-40.0);
				
				keyword = idnty +"_yshift";
				if (!inmap.get(keyword,yshift)) 
					inmap.getwdefault("yshift",yshift,0.0);
				
				keyword = idnty +"_gamma";
				if (!inmap.get(keyword,gam))
					inmap.getwdefault("gamma",gam,1.403);
				
				shift[0] = xshift;
				shift[1] = yshift;
				
				
			}
		
		bool eval(TinyVector<FLT,tet_mesh::ND> x, double t, Array<double,1> &val) const {
			
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
	
	
	class ibc_type {
		public:
			const static int ntypes = 3;
			enum ids {freestream,sphere,ringleb};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
						if (!strcmp(nin,names[i])) return(i);
				return(-1);
		}
	};
	const char ibc_type::names[ntypes][40] = {"freestream","sphere","ringleb"};

}

init_bdry_cndtn *tet_hp_cns_explicit::getnewibc(std::string ibcname) {
	init_bdry_cndtn *temp;
	int type;

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
		default: {
			return(tet_hp::getnewibc(ibcname));
		}
	}
	return(temp);
}






