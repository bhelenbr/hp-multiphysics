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
					blockdata.getwdefault("flowspeed",speed,1.0);

				keyword = idnty +"_flowangle";
				if (!blockdata.get(keyword,alpha)) 
					blockdata.getwdefault("flowangle",alpha,0.0);  

				keyword = idnty +"_perturb_amplitude";
				if (!blockdata.get(keyword,perturb_amp)) 
					blockdata.getwdefault("perturb_amplitude",perturb_amp,0.0); 

				keyword = idnty +"_RT";
				if (!blockdata.get(keyword,RT)) 
					blockdata.getwdefault("RT",RT,1.0);
				
				keyword = idnty +"_pressure";
				if (!blockdata.get(keyword,pr)) 
					blockdata.getwdefault("pressure",pr,1.0);
				
				keyword = idnty +"_gamma";
				if (!blockdata.get(keyword,gamma))
					blockdata.getwdefault("gamma",gamma,1.403);
					
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
								
				switch(n) {
					case(1):case(2): 
						if (r < inner) 
							return(0.0);
						else if (r < outer)
							return(pr/RT*vel(n-1)*0.5*(1.-cos(M_PI*(r-inner)/(outer-inner))));
						else
							return(pr/RT*vel(n-1));
					case(0):
						return(pr/RT);
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
					blockdata.getwdefault("flowspeed",speed,1.0);

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

				keyword = idnty +"_RT";
				if (!blockdata.get(keyword,RT)) 
					blockdata.getwdefault("RT",RT,1.0);
				
				keyword = idnty +"_pressure";
				if (!blockdata.get(keyword,pr)) 
					blockdata.getwdefault("pressure",pr,1.0);
				
				keyword = idnty +"_gamma";
				if (!blockdata.get(keyword,gamma))
					blockdata.getwdefault("gamma",gamma,1.403);
				
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
	
	
	

	

//	class accelerating : public init_bdry_cndtn {
//		private:
//			FLT speed,c,alpha;
//
//		public:
//			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
//				switch(n) {
//					case(0):
//						return(speed +c*pow(time,alpha));
//					case(1):
//						return(0.0);
//					case(2):
//						return(-(x(0)-1)*c*alpha*pow(time,alpha-1.0));
//				}
//				return(0.0);
//			}
//
//			void input(input_map &blockdata,std::string idnty) {
//				std::string keyword,val;
//				std::istringstream data;
//
//				keyword = idnty +"_speed";
//				if (!blockdata.get(keyword,speed)) 
//					blockdata.getwdefault("speed",speed,1.0);
//
//				keyword = idnty +"_coefficient";
//				if (!blockdata.get(keyword,c)) 
//					blockdata.getwdefault("coefficient",c,0.0);  
//
//				keyword = idnty +"_power";
//				if (!blockdata.get(keyword,alpha)) 
//					blockdata.getwdefault("power",alpha,0.0); 
//			}
//	};    

//	class impinge : public init_bdry_cndtn {
//		public:
//			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
//				switch(n) {
//					case(0):
//						return(1.0);
//					case(1):
//						return(-x(1)/x(0));
//				case(2):
//					return(0.0);
//
//				}
//				return(0.0);
//			}
//	};

//	class parameter_changer : public tri_hp_helper {
//		protected:
//			tri_hp_cns &x;
//			bdry_cns::surface *surf;
//			FLT delta_rho, rho_factor;
//			FLT delta_mu, mu_factor;
//			FLT delta_g, delta_g_factor;
//			FLT delta_rho2, rho2_factor;
//			FLT delta_mu2, mu2_factor;
//			FLT delta_sigma, sigma_factor;
//			int interval;
//
//		public:
//			parameter_changer(tri_hp_cns& xin) : tri_hp_helper(xin), x(xin) {
//				int bnum;
//
//				for(bnum=0;bnum<x.nebd;++bnum) 
//					if (surf = dynamic_cast<bdry_cns::surface *>(x.hp_ebdry(bnum))) break;
//
//				if (bnum > x.nebd -1) surf = 0;
//			}
//			void init(input_map& input, std::string idnty) {
//				std::string keyword, val;
//
//				keyword = idnty + "_delta_rho";
//				if (!input.get(keyword,delta_rho)) {
//					input.getwdefault("delta_rho",delta_rho,0.0);
//				}
//
//				keyword = idnty + "_rho_factor";
//				if (!input.get(keyword,rho_factor)) {
//					input.getwdefault("rho_factor",rho_factor,1.0);
//				}
//
//				keyword = idnty + "_delta_mu";
//				if (!input.get(keyword,delta_mu)) {
//					input.getwdefault("delta_mu",delta_mu,0.0);
//				}
//
//				keyword = idnty + "_mu_factor";
//				if (!input.get(keyword,mu_factor)) {
//					input.getwdefault("mu_factor",mu_factor,1.0);
//				}
//
//				keyword = "delta_g";
//				input.getwdefault(keyword,delta_g,0.0);
//				input[keyword] = "0.0"; // SO ONLY ONE BLOCK PER PROCESSOR CHANGES THIS
//
//				keyword = "delta_g_factor";
//				input.getwdefault(keyword,delta_g_factor,1.0);
//				input[keyword] = "1.0"; // SO ONLY ONE BLOCK PER PROCESSOR CHANGES THIS
//
//				input.getwdefault("parameter_interval",interval,1);
//
//				if (surf) {
//					std::string surfidnty = surf->base.idprefix;
//
//					keyword = surfidnty + "_delta_sigma";
//					input.getwdefault(keyword,delta_sigma,0.0);
//
//					keyword = surfidnty + "_sigma_factor";
//					input.getwdefault(keyword,sigma_factor,1.0);
//
//					keyword = surfidnty + "_matching_block";
//					if (!input.get(keyword,val)) {
//						delta_rho2 = 0.0;
//						rho2_factor = 1.0;
//						delta_mu2 = 0.0;
//						mu2_factor = 1.0;
//					}
//					else {                         
//						keyword = val + "_delta_rho";                    
//						if (!input.get(keyword,delta_rho2)) {
//							input.getwdefault("delta_rho",delta_rho2,0.0);
//						}
//
//						keyword = val + "_rho_factor";
//						if (!input.get(keyword,rho2_factor)) {
//							input.getwdefault("rho_factor",rho2_factor,1.0);
//						}
//
//						keyword = val + "_delta_mu";
//						if (!input.get(keyword,delta_mu2)) {
//							input.getwdefault("delta_mu",delta_mu2,0.0);
//						}
//
//						keyword = val + "_mu_factor";
//						if (!input.get(keyword,mu2_factor)) {
//							input.getwdefault("mu_factor",mu2_factor,1.0);
//						}
//					}
//				}
//			}
//			tri_hp_helper* create(tri_hp& xin) { return new parameter_changer(dynamic_cast<tri_hp_cns&>(xin)); }
//
//
//
//			void tadvance() {
//				if (!x.coarse_level) {
//					if ( (x.gbl->tstep % interval) +x.gbl->substep == 0) {
//
//						x.gbl->rho += delta_rho;
//						x.gbl->rho *= rho_factor;
//
//						x.gbl->mu  += delta_mu;
//						x.gbl->mu  *= mu_factor;
//
//						x.gbl->g += delta_g;
//						x.gbl->g *= delta_g_factor;
//
//						*x.gbl->log << "new density, viscosity, and gravity are " << x.gbl->rho << ' ' << x.gbl->mu << ' ' << x.gbl->g << std::endl;
//
//
//						if (surf) {
//							surf->gbl->rho2 += delta_rho2;
//							surf->gbl->rho2 *= rho2_factor;
//
//							surf->gbl->mu2  += delta_mu2;
//							surf->gbl->mu2  *= mu2_factor;
//
//							surf->gbl->sigma  += delta_sigma;
//							surf->gbl->sigma  *= sigma_factor;
//
//							*x.gbl->log << "matching block density, viscosity, and surface tension are " << surf->gbl->rho2 << ' ' << surf->gbl->mu2 << ' ' << surf->gbl->sigma << std::endl;
//						}
//					}
//				}
//				return;
//			}
//	};


//	class unsteady_body_force : public tri_hp_helper {
//		protected:
//			TinyVector<symbolic_function<1>,2> fcn;
//
//		public:
//			unsteady_body_force(tri_hp_cns& xin) : tri_hp_helper(xin) {}
//
//			void init(input_map& input, std::string idnty) {                
//				std::string keyword,val;
//				std::ostringstream nstr;
//
//				for(int n=0;n<2;++n) {
//					nstr.str("");
//					nstr << idnty << "_forcing" << n << std::flush;
//					if (input.find(nstr.str()) != input.end()) {
//						fcn(n).init(input,nstr.str());
//					}
//					else {
//						nstr.str("");
//						nstr << "forcing" << n << std::flush;
//						if (input.find(nstr.str()) != input.end()) {
//							fcn(n).init(input,nstr.str());
//						}
//						else {
//							std::cerr << "couldn't find forcing function" << std::endl;
//							exit(1);
//						}
//					}
//				}
//			}
//			tri_hp_helper* create(tri_hp& xin) { return new unsteady_body_force(dynamic_cast<tri_hp_cns&>(xin));}            
//
//			void tadvance() {
//				x.gbl->body(0) = fcn(0).Eval(0,x.gbl->time);
//				x.gbl->body(1) = fcn(1).Eval(0,x.gbl->time);
//				return;
//			}
//	};

	FLT xmax(TinyVector<FLT,2> &pt) {return(pt(0));}

//	class force_coupling : public tri_hp_helper {
//		protected:
//			FLT k_linear, k_torsion;
//			FLT mass, I;
//			int nboundary;
//			Array<bdry_cns::force_coupling *,1> hp_ebdry;
//			Array<eboundary_with_geometry<edge_bdry,naca> *,1> ebdry;
//			FLT w_a, w_b, w_c, w_d;
//			FLT w_a_tilda, w_b_tilda, w_c_tilda, w_d_tilda;
//			FLT k_a[gbl->nhist], k_b[gbl->nhist], k_c[gbl->nhist], k_d[gbl->nhist];
//			/* This is to dump output from the generic boundaries */
//			struct nullstream : std::ostream {
//				nullstream(): std::ios(0), std::ostream(0) {}
//			} ns;
//
//
//		public:
//			force_coupling(tri_hp_cns& xin) : tri_hp_helper(xin), w_a(0.0), w_b(0.0), w_c(0.0), w_d(0.0) {}
//			force_coupling(const force_coupling &in_fc, tri_hp_cns& xin) : k_linear(in_fc.k_linear), k_torsion(in_fc.k_torsion),
//				mass(in_fc.mass), I(in_fc.I), nboundary(in_fc.nboundary), tri_hp_helper(xin) {
//				int i,j;
//
//				hp_ebdry.resize(nboundary);
//				ebdry.resize(nboundary);
//				for (i = 0; i < nboundary; ++i) {
//					for (j=0;j<x.nebd;++j) {
//						if (x.ebdry(j)->idnum == in_fc.ebdry(i)->idnum) {
//						goto found;
//						}
//					}
//					std::cerr << "List of force boundaries is wrong" << std::endl;
//					exit(1);
//
//					found:
//					if (!(hp_ebdry(i) = dynamic_cast<bdry_cns::force_coupling *>(x.hp_ebdry(j)))) {
//						std:cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
//						exit(1);
//					}
//					if (!(ebdry(i) = dynamic_cast<eboundary_with_geometry<edge_bdry,naca> *>(x.ebdry(j)))) {
//						std::cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
//						exit(1);
//					}
//
//				}
//			}
//
//
//			tri_hp_helper* create(tri_hp& xin) {return new force_coupling(*this, dynamic_cast<tri_hp_cns&>(xin));}
//
//			void init(input_map& input, std::string idnty) {     
//				std::string bdrys;
//				std::istringstream bdryin;
//
//				if (!input.get(x.gbl->idprefix +"_y0",w_a)) 
//					input.getwdefault("y0",w_a,0.0);
//
//				if (!input.get(x.gbl->idprefix +"_dydt0",w_b)) 
//					input.getwdefault("dydt0",w_b,0.0);
//
//				if (!input.get(x.gbl->idprefix +"_theta0",w_c)) 
//					input.getwdefault("theta0",w_c,0.0);
//
//				if (!input.get(x.gbl->idprefix +"_dthetadt0",w_d)) 
//					input.getwdefault("dthetadt0",w_d,0.0);
//
//				if (!input.get(x.gbl->idprefix +"_k_linear",k_linear)) 
//					input.getwdefault("k_linear",k_linear,1.0);
//
//				if (!input.get(x.gbl->idprefix +"_k_torsion",k_torsion)) 
//					input.getwdefault("k_torsion",k_torsion,1.0);
//
//				if (!input.get(x.gbl->idprefix +"_mass",mass)) 
//					input.getwdefault("mass",mass,1.0);
//
//				if (!input.get(x.gbl->idprefix +"_I",I)) 
//					input.getwdefault("I",I,1.0);
//
//				if (!input.get(x.gbl->idprefix +"_nboundary",nboundary))
//					input.getwdefault("nboundary",nboundary,1);
//
//				if (!input.getline(x.gbl->idprefix +"_force_boundaries",bdrys)) {
//					if (!input.getline("force_boundaries",bdrys)) {
//						std::cerr << "No boundary number list" << std::endl;
//						exit(1);
//					}
//				}
//
//				hp_ebdry.resize(nboundary);
//				ebdry.resize(nboundary);
//				bdryin.str(bdrys);
//				int bnum;
//				int j;
//				for (int i = 0; i < nboundary; ++i) {
//					bdryin >> bnum;
//					for (j=0;j<x.nebd;++j) {
//						if (x.ebdry(j)->idnum == bnum) {
//							goto found;
//						}
//					}
//					std::cerr << "List of force boundaries is wrong" << std::endl;
//					exit(1);
//
//					found:
//						if (!(hp_ebdry(i) = dynamic_cast<bdry_cns::force_coupling *>(x.hp_ebdry(j)))) {
//							std:cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
//							exit(1);
//						}
//						if (!(ebdry(i) = dynamic_cast<eboundary_with_geometry<edge_bdry,naca> *>(x.ebdry(j)))) {
//							std::cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
//							exit(1);
//						}
//
//				}
//
//			}
//
//			void tadvance() {
//				if (x.coarse_flag) return;
//
//				int stage = x.gbl->substep +x.gbl->esdirk;
//				if (stage > 0) {
//					k_a[stage-1]=(w_a-w_a_tilda)/(1/x.gbl->adirk(stage-1,stage-1));
//					k_b[stage-1]=(w_b-w_b_tilda)/(1/x.gbl->adirk(stage-1,stage-1)); 
//					k_c[stage-1]=(w_c-w_c_tilda)/(1/x.gbl->adirk(stage-1,stage-1));
//					k_d[stage-1]=(w_d-w_d_tilda)/(1/x.gbl->adirk(stage-1,stage-1));
//				}
//
//				if (x.gbl->substep == 0) {
//					w_a_tilda=w_a;
//					w_b_tilda=w_b;   
//					w_c_tilda=w_c;
//					w_d_tilda=w_d;  
//				}
//
//				for (int s=0;s<stage;++s) {
//					w_a_tilda=w_a_tilda+(x.gbl->adirk(stage,s))*k_a[s];
//					w_b_tilda=w_b_tilda+(x.gbl->adirk(stage,s))*k_b[s];
//					w_c_tilda=w_c_tilda+(x.gbl->adirk(stage,s))*k_c[s];
//					w_d_tilda=w_d_tilda+(x.gbl->adirk(stage,s))*k_d[s];
//				}
//
//				/* EXTRAPOLATE */
//				if (stage  && x.gbl->dti > 0.0) {
//					FLT constant =  x.gbl->cdirk(x.gbl->substep);
//					w_a += constant*k_a[stage-1]; 
//					w_b += constant*k_b[stage-1];
//					w_c += constant*k_c[stage-1];
//					w_d += constant*k_d[stage-1];
//
//					/* FIX POSITIONS */
//					FLT dy = constant*k_a[stage-1];
//					FLT dtheta = constant*k_c[stage-1];
//					FLT w_a_first = w_a -dy;
//					TinyVector<FLT,2> ctr(0.0,w_a_first);
//					TinyVector<FLT,2> vel(0.0,w_b);
//					TinyVector<FLT,2> disp(0.0,dy);
//					rigid_body(dtheta,w_d,ctr,disp,vel);	
//				}
//
//
//				return;
//			}
//
//			void mg_restrict() {
//				if(x.coarse_level) {
//					tri_hp *fmesh = dynamic_cast<tri_hp *>(x.fine);
//					force_coupling *fine_helper = dynamic_cast<force_coupling *>(fmesh->helper);
//
//					/* LOOP THROUGH POINTS TO TO CALCULATE POSITION OF COARSE POINTS  */
//					int i,j,n,tind;
//					for(i=0;i<x.npnt;++i) {
//						tind = x.fcnnct(i).tri;
//
//						for(n=0;n<x.ND;++n)
//						x.pnts(i)(n) = 0.0;
//
//						for(j=0;j<3;++j) {
//						for(n=0;n<x.ND;++n)
//							x.pnts(i)(n) += x.fcnnct(i).wt(j)*fmesh->pnts(fmesh->tri(tind).pnt(j))(n);
//						}
//					}
//
//
//					for (int i=0; i<nboundary; ++i) {
//						ebdry(i)->geometry_object.theta = fine_helper->ebdry(i)->geometry_object.theta; 
//						ebdry(i)->geometry_object.pos = fine_helper->ebdry(i)->geometry_object.pos; 
//					}
//				}
//
//				return;
//			}
//
//
//			void rigid_body(FLT dtheta, FLT omega, TinyVector<FLT,2> ctr, TinyVector<FLT,2> disp, TinyVector<FLT,2> vel) {
//				TinyVector<FLT,2> dx;
//
//				/* UPDATE MESH POSITION */
//				FLT r,cost,sint; 
//				FLT cosdt = cos(dtheta);    
//				FLT sindt = sin(dtheta);    
//				for (int i=0;i<x.npnt;++i) {
//					dx = x.pnts(i) -ctr;
//						r = sqrt(dx(0)*dx(0) +dx(1)*dx(1));
//					cost = dx(0)/r;
//					sint = dx(1)/r;
//					x.pnts(i)(0) += -(r-r*cosdt)*cost -r*sindt*sint +disp(0);
//					x.pnts(i)(1) += -(r-r*cosdt)*sint +r*sindt*cost +disp(1);						
//				}
//
//				for (int i=0; i<nboundary; ++i) {
//					hp_ebdry(i)->set_ctr_rot(ctr);
//					hp_ebdry(i)->set_vel(vel);
//					hp_ebdry(i)->set_omega(w_d);
//					ebdry(i)->geometry_object.theta += -dtheta; 
//					ebdry(i)->geometry_object.pos(0) += disp(0) -0.25*(cos(w_c) -cos(w_c -dtheta)); 
//					ebdry(i)->geometry_object.pos(1) += disp(1) -0.25*(sin(w_c) -sin(w_c -dtheta)); 
//
//					/* CAN FIX ENDPOINTS TOO (NOT NECESSARY) */
////						v0 = x.seg(ebdry(i)->seg(0)).pnt(0);
////						FLT distance2 = (1-i)*0.75 -i*0.25;
////						FLT distance1 = (1-i)*(-0.25) +i*0.75;
////						x.pnts(v0)(0) = distance1*cos(w_c);
////						x.pnts(v0)(1) = w_a +distance1*sin(w_c);
////						v0 = x.seg(ebdry(i)->seg(ebdry(i)->nseg-1)).pnt(1);
////						x.pnts(v0)(0) = distance2*cos(w_c);
////						x.pnts(v0)(1) = w_a +distance2*sin(w_c);
//
//					hp_ebdry(i)->curv_init();
//
//				}
//			}
//
//
//
//			void update(int stage) {
//				if (!x.coarse_flag) {
//
//					/* GET FORCE & TORQUE */
//					FLT force_y = 0.0;
//					// FLT force_x = 0.0;
//					FLT moment = 0.0;
//					for (int i=0; i<nboundary; ++i) {                                  
//					     /* FORCE BOUNDARY TO CALCULATE ALL FLUXES */
//						hp_ebdry(i)->output(ns,tri_hp::tecplot);
//						force_y += hp_ebdry(i)->diff_flux(1);
//						moment += hp_ebdry(i)->moment;
//					}   
//
//					FLT dt = 1.0/x.gbl->dti;
//					int stage = x.gbl->substep +x.gbl->esdirk;
//
//					/* CALCULATE NEW POSITION */
//					double w_a_first = w_a;
//					double w_c_first = w_c;
//					w_a=(w_a_tilda+((1/x.gbl->adirk(stage,stage)))*dt*(w_b_tilda+((1/x.gbl->adirk(stage,stage)))*dt*force_y))/(1+k_linear*(dt*dt*(1/x.gbl->adirk(stage,stage))*(1/x.gbl->adirk(stage,stage))));         // translational displacement
//					w_b=(w_b_tilda-((1/x.gbl->adirk(stage,stage)))*dt*(w_a_tilda*k_linear-force_y))/(1+k_linear*(dt*dt*(1/x.gbl->adirk(stage,stage))*(1/x.gbl->adirk(stage,stage))));                                    // translational velocity
//					w_c=(w_c_tilda+((1/x.gbl->adirk(stage,stage)))*dt*(w_d_tilda+((1/x.gbl->adirk(stage,stage)))*dt*moment))/(1+k_torsion*(dt*dt*(1/x.gbl->adirk(stage,stage))*(1/x.gbl->adirk(stage,stage))));         // rotational displacement
//					w_d=(w_d_tilda-((1/x.gbl->adirk(stage,stage)))*dt*(w_c_tilda*k_torsion-moment))/(1+k_torsion*(dt*dt*(1/x.gbl->adirk(stage,stage))*(1/x.gbl->adirk(stage,stage))));                                   // rotational velocity
//
//
//					/* FIX POSITIONS */
//					FLT dy = w_a -w_a_first;
//					FLT dtheta = w_c-w_c_first;
//					TinyVector<FLT,2> ctr(0.0,w_a_first);
//					TinyVector<FLT,2> vel(0.0,w_b);
//					TinyVector<FLT,2> disp(0.0,dy);
//					rigid_body(dtheta,w_d,ctr,disp,vel);	
//				}
//				return;
//			}
//	};

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