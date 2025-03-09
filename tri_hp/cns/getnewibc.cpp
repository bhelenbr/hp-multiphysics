/*
 *  getnewibc_cns.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cns.h"
#include "bdry_cns.h"
#include <tri_boundary.h>

#include <fstream>
namespace ibc_cns {

	class freestream : public init_bdry_cndtn {
		private:
		FLT alpha, speed,perturb_amp, gamma;

		public:

			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				FLT amp = (time > 0.0 ? 0.0 : perturb_amp); 

				switch(n) {
					case(0):
						return(1.0/gamma);
					case(1):
						return(speed*cos(alpha) +amp*x(0)*(1.0-x(0)));
					case(2):
						return(speed*sin(alpha));
					case(3):
						return(1.0/gamma);
				}
				return(0.0);
			}

			void init(input_map &blockdata,std::string idnty) {
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
					blockdata.getwdefault("gamma",gamma,1.4);

				alpha *= M_PI/180.0;
			}
    };

	class sphere : public init_bdry_cndtn {
		private:
			FLT speed,angle,inner,outer,RT,pr,gamma;
			TinyVector<FLT,tri_mesh::ND> vel;

		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				FLT r;

				r = sqrt(x(0)*x(0) +x(1)*x(1));
				switch(n) {
					case(0):
						return(pr);
					case(1):case(2): 
						if (r < inner) 
							return(0.0);
						else if (r < outer)
							return(vel(n-1)*0.5*(1.-cos(M_PI*(r-inner)/(outer-inner))));
						else
							return(vel(n-1));
					case(3):
						return(RT);
				}
				return(0.0);
			}

			void init(input_map &blockdata,std::string idnty) {
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
					blockdata.getwdefault("gamma",gamma,1.4);
				
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
		protected:
			tri_hp_cns &x;
		
		private:
			FLT angle,xshift,yshift;
			FLT gam; // temp figure out how to load hp_cns_gbl->gamma
			double q0;  // Initial guess for q in iteration //
			double theta0; // Initial guess for theta in iteration //
			double scale; // Scale point before calculating solution //
			double shift[2]; // Spatial shift of point before calculating solution //
		
		public:
            ringleb(tri_hp_cns& cns) : x(cns) {}
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				Array<double,1> val(4);
				bool error;
			
				error = eval(x,time,val);			

				return(val(n));
			}
		
			void init(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;
				
				keyword = idnty +"_q0";
				if (!blockdata.get(keyword,q0)) 
					blockdata.getwdefault("q0",q0,0.02);
				
				keyword = idnty +"_theta0";
				if (!blockdata.get(keyword,theta0)) 
					blockdata.getwdefault("theta0",theta0,0.99*M_PI/2.);
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
				
				shift[0] = xshift;
				shift[1] = yshift;
                
                if (!blockdata.get(idnty+ "_gamma",gam)) blockdata.getwdefault("gamma",gam,1.4);

			}
		
			bool eval(TinyVector<FLT,tri_mesh::ND> xp, double t, Array<double,1> &val) const {
				
				double q = q0;
				double theta = theta0;
				double dtheta = -1.0e-4;
				double dq = 1.0e-4;
				double jac[2][2], deti, error[2];
				bool conservative = false;
				const int maxiter = 30;
				
				double xyz[3];
				xyz[0] = scale*xp(0) +shift[0];
				xyz[1] = scale*xp(1) +shift[1];
                
				
				if (fabs(xyz[1]) < 1.0e-8)
					theta = M_PI/2.0;
				
				double delta[2] = {1.0, 1.0};
				int niter = 0;
				while (fabs(delta[0])+fabs(delta[1]) > 1.0e-10 && niter++ < maxiter) {
					if (fabs(xyz[1]) > 1.0e-10) {
						calculateError(xyz, error, q, theta);
						calculateError(xyz, delta, q+dq, theta);
						jac[0][0] = (delta[0]-error[0])/dq;
						jac[1][0] = (delta[1]-error[1])/dq;
						calculateError(xyz, delta, q, theta+dtheta);
						jac[0][1] = (delta[0]-error[0])/dtheta;
						jac[1][1] = (delta[1]-error[1])/dtheta;
						deti = 1./(jac[0][0]*jac[1][1] -jac[0][1]*jac[1][0]);
                    
						delta[0] = deti*(error[0]*jac[1][1] -error[1]*jac[0][1]);
						delta[1] = deti*(error[1]*jac[0][0] -error[0]*jac[1][0]);
						q -= delta[0];
						theta -= delta[1];
					}
					else {
						theta = M_PI/2.;
						calculateError(xyz, error, q, theta);
						calculateError(xyz, delta, q+dq, theta);
						jac[0][0] = (delta[0]-error[0])/dq;
                        delta[0] = error[0]/jac[0][0];
						q -= delta[0];
					}
                    
					
					if (fabs(theta) > M_PI/2.) {
						theta = theta +(theta > 0 ? -M_PI/2. : M_PI/2.);
					}
                    // *x.gbl->log << niter <<  ' ' << xyz[0] << ' ' << xyz[1] << ' ' << q << ' ' << theta << ' ' << fab*/s(delta[0]) +fabs(delta[1]) << std::endl;
				}
                
				if (niter >= maxiter || !(q >= 0.0)) {
					*x.gbl->log << "RINGLEB NOT_CONVERGED: " << xyz[0] << ' ' << xyz[1] << ' ' << niter << ' ' << q << ' ' << theta << ' ' << fabs(delta[0])+fabs(delta[1]) << std::endl;
                    exit(0);
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
					yp = xyz[1]/fabs(xyz[1])*psi/(q*r)*sqrt(max(1.-psi*psi*q*q,0.0));
					error[1] = xyz[1]-yp;
				}
				else {
					error[1] = xyz[1]-0;
				}
                
                
				
				return;
			}
		
			
	};


    class nozzle : public init_bdry_cndtn {
//        protected:
//            tri_hp_cns &x;
        
        private:
            FLT RTo, Po, xthroat;
            FLT gam; // temp figure out how to load hp_cns_gbl->gamma
            symbolic_function<2> Aratio;
            symbolic_function<2> dArdx;
        
        public:
            FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
                Array<double,1> val(4);
                bool error;
            
                error = eval(x,time,val);

                return(val(n));
            }
        
            void init(input_map &blockdata,std::string idnty) {
                std::string keyword,val;
                std::istringstream data;
                
                keyword = idnty +"_Aratio";
                Aratio.init(blockdata,keyword);
                
                keyword = idnty +"_dArdx";
                dArdx.init(blockdata,keyword);

                
                keyword = idnty +"_RTo";
                if (!blockdata.get(keyword,RTo))
                    blockdata.getwdefault("RTo",RTo,287.058);
                
                keyword = idnty +"_Po";
                if (!blockdata.get(keyword,Po))
                    blockdata.getwdefault("Po",Po,101325.0);
                
                keyword = idnty +"_xthroat";
                if (!blockdata.get(keyword,xthroat))
                    blockdata.getwdefault("xthroat",xthroat,5.0);
                
                if (!blockdata.get(idnty+ "_gamma",gam)) blockdata.getwdefault("gamma",gam,1.4);

            }
        
            bool eval(TinyVector<FLT,tri_mesh::ND> x, double t, Array<double,1> &val) const {
                
                int it = 0;
                int maxit = 1000;
                FLT tol = 1.0e-10;
//                FLT lam = 1.0;
//                FLT fold;
//                int armijo = 0;
//                FLT M_armijo, f_armijo, fp_armijo;
//                FLT alpha = 1.0e-4;
                
                FLT Ar = Aratio.Eval(x,t);
                FLT M;
                
                if (x(0)<xthroat){
                    M = 0.5;
                }
                else{
                    M = 1.5;
                }

                FLT f = M*Ar-pow(((2.0+(gam-1.0)*M*M)/(gam+1.0)),((gam+1.0)/(2.0*(gam-1.0))));
                FLT fp = Ar - (2.0*M*pow((((gam - 1.0)*M*M + 2.0)/(gam + 1.0)),((gam + 1.0)/(2.0*gam - 2.0) - 1.0))*(gam - 1.0))/(2.0*gam - 2.0);
                

//                fold = f;
                while(fabs(f)>tol && it<maxit){
                    
                  /*  Armijo Line Search (no longer necessary) */
//                    lam = 1.0;
//
//                    M_armijo = M-lam*(f/fp);
//
//                    f_armijo = M_armijo*Ar-(pow(((2.0+(gam-1.0)*M_armijo*M_armijo)/(gam+1.0)),((gam+1.0)/(2.0*(gam-1.0)))));
//                    fp_armijo = Ar - (2.0*M_armijo*pow((((gam - 1.0)*M_armijo*M_armijo + 2.0)/(gam + 1.0)),((gam + 1.0)/(2.0*gam - 2.0) - 1.0))*(gam - 1.0))/(2.0*gam - 2.0);
//
//                    while(fabs(f_armijo)>(1.0-alpha*lam)*fabs(fold) && armijo<100){
//
//                        lam = 0.5*lam;
//                        M_armijo = M-lam*(f_armijo/fp_armijo);
//                        f_armijo = M_armijo*Ar-(pow(((2.0+(gam-1.0)*M_armijo*M_armijo)/(gam+1.0)),((gam+1.0)/(2.0*(gam-1.0)))));
//                        fp_armijo = Ar - (2.0*M_armijo*pow((((gam - 1.0)*M_armijo*M_armijo + 2.0)/(gam + 1.0)),((gam + 1.0)/(2.0*gam - 2.0) - 1.0))*(gam - 1.0))/(2.0*gam - 2.0);
//
//                        armijo++;
//                    }
//                    armijo = 0;
//
//                    M = M_armijo;
//
//                    f = f_armijo;
//                    fp = fp_armijo;
                    
                    M = M-(f/fp);
                    f = M*Ar-pow(((2.0+(gam-1.0)*M*M)/(gam+1.0)),((gam+1.0)/(2.0*(gam-1.0))));
                    fp = Ar - (2.0*M*pow((((gam - 1.0)*M*M + 2.0)/(gam + 1.0)),((gam + 1.0)/(2.0*gam - 2.0) - 1.0))*(gam - 1.0))/(2.0*gam - 2.0);

                    it++;

                }
                
//                std::cout << setprecision(14) << x << " " << M << " " << f << std::endl;

                if (it == maxit){
                    printf("Newton Solver in ibc nozzle exceeded iteration limit \n");
                    if (fabs(f)>1.0e-10){
                        std::cerr << "Nozzle not converged" << std::endl;
                        sim::abort(__LINE__,__FILE__,&std::cerr);
                    }
                }
                
                
                FLT RT = RTo/(1.0+0.5*(gam-1.0)*M*M);
                FLT p = Po/(pow((1.0+0.5*(gam-1.0)*M*M),((gam/(gam-1.0)))));
                
                
                FLT c = sqrt(gam*RT);
                FLT u = M*c;
                FLT v;
                FLT rho = p/RT;
                
                FLT dArdM = (2.0*pow((((gam - 1.0)*M*M + 2.0)/(gam + 1.0)),((gam + 1.0)/(2.0*gam - 2.0) - 1.0))*(gam - 1.0))/(2.0*gam - 2.0) - pow((((gam - 1.0)*M*M + 2.0)/(gam + 1.0)),((gam + 1)/(2*gam - 2)))/(M*M);
                FLT dAdx = dArdx.Eval(x,t);
                FLT dMdx = pow(dArdM,-1.0)*dAdx;
                
                FLT dcdM = -sqrt(gam)*(M*RTo*(gam/2.0 - 1.0/2.0))/(pow(((gam/2.0 - 1.0/2.0)*M*M + 1.0),2)*sqrt(RTo/((gam/2.0 - 1.0/2.0)*M*M + 1.0)));
                FLT dudx = dMdx*(sqrt(gam*RT)+M*dcdM);
                
                FLT dpdM = -(2.0*M*Po*gam*(gam/2.0 - 1.0/2.0))/(pow(((gam/2.0 - 1.0/2.0)*M*M + 1.0),(gam/(gam - 1.0) + 1.0))*(gam - 1.0));
                FLT dRTdM = -(2.0*M*RTo*(gam/2.0 - 1.0/2.0))/pow(((gam/2.0 - 1.0/2.0)*M*M + 1.0),2);
                FLT drhodx = (1/(RT*RT))*(RT*dpdM*dMdx-p*dRTdM*dMdx);
                
                v = -(x(1)/rho)*(rho*dudx+u*drhodx);
                
                    
                val(0) = p;
                val(1) = u;
                val(2) = v;
                val(3) = RT;
                
                
                return true;
                
            }
        
            
    };
	
	class helper_type {
		public:
			const static int ntypes = 3;
			enum ids {parameter_changer,unsteady_body_force,force_coupling};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
					if (!strcmp(nin,names[i])) return(i);
				return(-1);
			}
	};
	const char helper_type::names[ntypes][40] = {"parameter_changer","unsteady_body_force","force_coupling"};

	class ibc_type {
		public:
			const static int ntypes = 8;
			enum ids {freestream,sphere,accelerating,impinge,stokes_drop_gas,stokes_drop_liquid,ringleb,nozzle};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
					if (!strcmp(nin,names[i])) return(i);
				return(-1);
			}
	};
	const char ibc_type::names[ntypes][40] = {"freestream","sphere","accelerating","impinge","stokes_drop_gas","stokes_drop_liquid","ringleb","nozzle"};

}


init_bdry_cndtn *tri_hp_cns::getnewibc(std::string name) {
	std::string keyword,ibcname;
	init_bdry_cndtn *temp;
	int type;

	type = ibc_cns::ibc_type::getid(name.c_str());
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
			temp = new ibc_cns::ringleb(*this);
			break;
		}
        case ibc_cns::ibc_type::nozzle: {
            temp = new ibc_cns::nozzle;
            break;
        }
		default: {
			return(tri_hp::getnewibc(name));
		}
	}
	return(temp);
}
