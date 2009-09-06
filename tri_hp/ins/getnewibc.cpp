/*
 *  getnewibc_ins.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ins.h"
#include "bdry_ins.h"
#include <tri_boundary.h>

namespace ibc_ins {

	class freestream : public init_bdry_cndtn {
		private:
			FLT alpha, speed,perturb_amp;

		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				FLT amp = (time > 0.0 ? 0.0 : perturb_amp); 
				switch(n) {
					case(0):
						return(speed*cos(alpha) +amp*x(0)*(1.0-x(0)));
					case(1):
						return(speed*sin(alpha));
					case(2):
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

	class sphere : public init_bdry_cndtn {
		private:
			FLT speed,angle,inner,outer;
			TinyVector<FLT,tri_mesh::ND> vel;

		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				FLT r;

				r = sqrt(x(0)*x(0) +x(1)*x(1));
				switch(n) {
					case(0):case(1): 
						if (r < inner) 
							return(0.0);
						else if (r < outer)
							return(vel(n)*0.5*(1.-cos(M_PI*(r-inner)/(outer-inner))));
						else
							return(vel(n));
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

				vel(0) = speed*cos(angle);
				vel(1) = speed*sin(angle);
			}
	};


	class accelerating : public init_bdry_cndtn {
		private:
			FLT speed,c,alpha;

		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
				switch(n) {
					case(0):
						return(speed +c*pow(time,alpha));
					case(1):
						return(0.0);
					case(2):
						return(-(x(0)-1)*c*alpha*pow(time,alpha-1.0));
				}
				return(0.0);
			}

			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_speed";
				if (!blockdata.get(keyword,speed)) 
					blockdata.getwdefault("speed",speed,1.0);

				keyword = idnty +"_coefficient";
				if (!blockdata.get(keyword,c)) 
					blockdata.getwdefault("coefficient",c,0.0);  

				keyword = idnty +"_power";
				if (!blockdata.get(keyword,alpha)) 
					blockdata.getwdefault("power",alpha,0.0); 
			}
	};    

	class impinge : public init_bdry_cndtn {
		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
				switch(n) {
					case(0):
						return(1.0);
					case(1):
						return(-x(1)/x(0));
				case(2):
					return(0.0);

				}
				return(0.0);
			}
	};

	class stokes_drop_gas : public init_bdry_cndtn {
		private:
			FLT outer_limit;
			FLT mu_g, kappa;
			TinyVector<FLT,2> center;
			FLT frame_vel;

		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
				FLT r,sint,cost;
				FLT ur,ut;    

				x -= center;

				r = sqrt(x(0)*x(0) +x(1)*x(1));
				sint = x(0)/r;
				cost = -x(1)/r;
				ur = -(16.0*r*r*r+16.0*r*r*r*kappa-8.0*r*r-12.0*r*r*kappa+kappa)/(r*r*r)/(1.0+kappa)*cost/16.0;
				ut = sint*(32.0*r*r*r+32.0*r*r*r*kappa-8.0*r*r-12.0*r*r*kappa-kappa)/(r*r*r)/(1.0+kappa)/32.0;
				switch(n) {
					case(0):
						if (r < outer_limit)
							return(ur*sint+ut*cost);
						else
							return(0.0);
					case(1):
						if (r < outer_limit)
							return(-ur*cost+ut*sint -frame_vel);
						else
							return(1.0 -frame_vel);
					case(2):
						if (r < outer_limit)
							return(mu_g/2*cost*(2+3*kappa)/(2*r*r*(1+kappa))); 
						else
							return(0.0);
				}
				return(0.0);
			}

			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_outer_radius";
				if (!blockdata.get(keyword,outer_limit))
					blockdata.getwdefault("outer_radius",outer_limit,75.0);

				keyword = idnty +"_frame_velocity";
				if (!blockdata.get(keyword,frame_vel))
					blockdata.getwdefault("frame_velocity",frame_vel,0.0);

				keyword = idnty +"_center";
				center = 0.0;
				if (!blockdata.get(keyword,center.data(),2))
					blockdata.getwdefault("center",center.data(),2,center.data());

				string blkname = idnty.substr(0,idnty.find('_'));
				keyword = blkname +"_mu";
				if (!blockdata.get(keyword,mu_g)) {
					std::cerr << "couldn't find mu of gas " << keyword << std::endl;
					exit(1);
				}

				keyword = blkname +"_liquid";
				if (!blockdata.get(keyword,val)) { 
					std::cerr << "couldn't find identity of liquid block" << std::endl;
					exit(1);
				}

				FLT mu_l;
				keyword = val +"_mu";
				if (!blockdata.get(keyword,mu_l)) {
					std::cerr << "couldn't find mu of liquid" << std::endl;
					exit(1);
				}
				kappa = mu_l/mu_g;
			}
	};

	class stokes_drop_liquid : public init_bdry_cndtn {
		private:
			FLT outer_limit;
			FLT rho_l, mu_l, kappa, sigma;
			FLT frame_vel;
			TinyVector<FLT,2> center;

		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
				FLT r,sint,cost;
				FLT ur,ut;

				x -= center;

				r = sqrt(x(0)*x(0) +x(1)*x(1));

				sint = x(0)/(r+FLT_EPSILON);
				cost = (x(1) > 0.0 ? -1 : 1)*sqrt(1.-sint*sint);
				ur = -(4.0*r*r-1.0)*cost/(1.0+kappa)/2.0;
				ut = sint*(8.0*r*r-1.0)/(1.0+kappa)/2.0;
				switch(n) {
					case(0):
						return(ur*sint+ut*cost);
					case(1):
						return(-ur*cost+ut*sint -frame_vel);
					case(2):
						return(4*sigma -5*mu_l*r*cost*4/(1+kappa)); 
				}
				return(0.0);
			}

			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_frame_velocity";
				if (!blockdata.get(keyword,frame_vel))
					blockdata.getwdefault("frame_velocity",frame_vel,0.0);

				keyword = idnty +"_center";
				center = 0.0;
				if (!blockdata.get(keyword,center.data(),2))
					blockdata.getwdefault("center",center.data(),2,center.data());

				string blkname = idnty.substr(0,idnty.find('_'));
				keyword = blkname +"_rho";
				if (!blockdata.get(keyword,rho_l)) {
					std::cerr << "couldn't find rho of liquid" << std::endl;
					exit(1);
				}

				keyword = blkname +"_mu";
				if (!blockdata.get(keyword,mu_l)) {
					std::cerr << "couldn't find mu of liquid" << std::endl;
					exit(1);
				}
				keyword = blkname +"_liquid_bdry";
				if (!blockdata.get(keyword,val)) { 
					std::cerr << "couldn't find identity of liquid boundary" << std::endl;
					exit(1);
				}

				keyword = val +"_sigma";
				if (!blockdata.get(keyword,sigma)) {
					std::cerr << "couldn't find sigma" << std::endl;
					exit(1);
				}				

				keyword = blkname +"_gas";
				if (!blockdata.get(keyword,val)) { 
					std::cerr << "couldn't find identity of gas block" << std::endl;
					exit(1);
				}

				FLT mu_g;
				keyword = val +"_mu";
				if (!blockdata.get(keyword,mu_g)) {
					std::cerr << "couldn't find mu of gas " << keyword << std::endl;
					exit(1);
				}
				kappa = mu_l/mu_g;
			}
	};

	class parameter_changer : public tri_hp_helper {
		protected:
			tri_hp_ins &x;
			bdry_ins::surface *surf;
			FLT delta_rho, rho_factor;
			FLT delta_mu, mu_factor;
			FLT delta_g, delta_g_factor;
			FLT delta_rho2, rho2_factor;
			FLT delta_mu2, mu2_factor;
			FLT delta_sigma, sigma_factor;
			int interval;

		public:
			parameter_changer(tri_hp_ins& xin) : tri_hp_helper(xin), x(xin) {
				int bnum;

				for(bnum=0;bnum<x.nebd;++bnum) 
					if (surf = dynamic_cast<bdry_ins::surface *>(x.hp_ebdry(bnum))) break;

				if (bnum > x.nebd -1) surf = 0;
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

				if (surf) {
					std::string surfidnty = surf->base.idprefix;

					keyword = surfidnty + "_delta_sigma";
					input.getwdefault(keyword,delta_sigma,0.0);

					keyword = surfidnty + "_sigma_factor";
					input.getwdefault(keyword,sigma_factor,1.0);

					keyword = surfidnty + "_matching_block";
					if (!input.get(keyword,val)) {
						delta_rho2 = 0.0;
						rho2_factor = 1.0;
						delta_mu2 = 0.0;
						mu2_factor = 1.0;
					}
					else {                         
						keyword = val + "_delta_rho";                    
						if (!input.get(keyword,delta_rho2)) {
							input.getwdefault("delta_rho",delta_rho2,0.0);
						}

						keyword = val + "_rho_factor";
						if (!input.get(keyword,rho2_factor)) {
							input.getwdefault("rho_factor",rho2_factor,1.0);
						}

						keyword = val + "_delta_mu";
						if (!input.get(keyword,delta_mu2)) {
							input.getwdefault("delta_mu",delta_mu2,0.0);
						}

						keyword = val + "_mu_factor";
						if (!input.get(keyword,mu2_factor)) {
							input.getwdefault("mu_factor",mu2_factor,1.0);
						}
					}
				}
			}
			tri_hp_helper* create(tri_hp& xin) { return new parameter_changer(dynamic_cast<tri_hp_ins&>(xin)); }



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


						if (surf) {
							surf->gbl->rho2 += delta_rho2;
							surf->gbl->rho2 *= rho2_factor;

							surf->gbl->mu2  += delta_mu2;
							surf->gbl->mu2  *= mu2_factor;

							surf->gbl->sigma  += delta_sigma;
							surf->gbl->sigma  *= sigma_factor;

							*x.gbl->log << "matching block density, viscosity, and surface tension are " << surf->gbl->rho2 << ' ' << surf->gbl->mu2 << ' ' << surf->gbl->sigma << std::endl;
						}
					}
				}
				return;
			}
	};


	class unsteady_body_force : public tri_hp_helper {
		protected:
			TinyVector<symbolic_function<1>,2> fcn;

		public:
			unsteady_body_force(tri_hp_ins& xin) : tri_hp_helper(xin) {}

			void init(input_map& input, std::string idnty) {                
				std::string keyword,val;
				std::ostringstream nstr;

				for(int n=0;n<2;++n) {
					nstr.str("");
					nstr << idnty << "_forcing" << n << std::flush;
					if (input.find(nstr.str()) != input.end()) {
						fcn(n).init(input,nstr.str());
					}
					else {
						nstr.str("");
						nstr << "forcing" << n << std::flush;
						if (input.find(nstr.str()) != input.end()) {
							fcn(n).init(input,nstr.str());
						}
						else {
							std::cerr << "couldn't find forcing function" << std::endl;
							exit(1);
						}
					}
				}
			}
			tri_hp_helper* create(tri_hp& xin) { return new unsteady_body_force(dynamic_cast<tri_hp_ins&>(xin));}            

			void tadvance() {
				x.gbl->body(0) = fcn(0).Eval(0,x.gbl->time);
				x.gbl->body(1) = fcn(1).Eval(0,x.gbl->time);
				return;
			}
    };

    FLT xmax(TinyVector<FLT,2> &pt) {return(pt(0));}

    class translating_drop : public parameter_changer {
			private:
				Array<FLT,1> avg;
				FLT penalty1,penalty2;

			public:
				translating_drop(tri_hp_ins& xin) : parameter_changer(xin) {
					avg.resize(1+x.ND+x.NV);
				}

			void init(input_map& input, std::string idnty) {
				parameter_changer::init(input,idnty);

				std::string keyword;
				keyword = idnty + "_penalty1_parameter";
				input.getwdefault(keyword,penalty1,0.5);

				keyword = idnty + "_penalty2_parameter";
				input.getwdefault(keyword,penalty2,0.5);
			}

			tri_hp_helper* create(tri_hp& xin) { return new translating_drop(dynamic_cast<tri_hp_ins&>(xin)); }
			void calculate_stuff() {
				bdry_ins::surface::global *gbl = surf->gbl;

				/* DETRMINE CORRECTION TO CONSERVE AREA */
				/* IMPORTANT FOR STEADY SOLUTIONS */
				/* SINCE THERE ARE MULTIPLE STEADY-STATES */
				/* TO ENSURE GET CORRECT VOLUME */
				FLT rbar, kc; 
				kc = gbl->sigma/(x.gbl->mu +gbl->mu2);
				x.integrated_averages(avg);
				rbar  = pow(3.*0.5*avg(0),1.0/3.0);
#ifdef DROP
				if (!x.coarse_flag) gbl->vflux =  penalty1*kc*(rbar -0.5);
				if (!x.coarse_flag) tri_hp_ins::mesh_ref_vel(1) = penalty2*kc*avg(2) +avg(4);    
#endif

				/* C_D TO G CONVERSION REMINDER 
				re = 1.0/gbl->mu2;
				cd = 24./re*(1 +0.1935*pow(re,0.6305));
				cd /= 16.0; // (1/2 rho u^2 * Pi r^2 / 2 pi);
				g = amp*(avg +avg) +12.*cd/(gbl->rho -gbl->rho2);
				*/
				return;
			}


			void rsdl(int stage) {
				/* if (x.gbl->dti == 0) */ calculate_stuff();
				return;
			}

			void setup_preconditioner() {
				/* if (gbl->dti == 0.0) */ calculate_stuff();
				return;
			}

			void tadvance() {

				if (x.coarse_flag) return;

				calculate_stuff();
#ifdef DROP
				if ( (x.gbl->tstep % interval) +x.gbl->substep == 0) {
					*x.gbl->log << "#gravity, velocity, height: " << x.gbl->g << ' ' << tri_hp_ins::mesh_ref_vel(1) << ' ';                        
					int v0 = x.seg(surf->base.seg(0)).pnt(0);
					int v1 = x.seg(surf->base.seg(surf->base.nseg-1)).pnt(1);
					FLT height = x.pnts(v0)(1)-x.pnts(v1)(1);
					*x.gbl->log << height << std::endl;
				}

				parameter_changer::tadvance();
#endif
				surf->findmax(xmax);
				return;
			}
	};
	
	



	class force_coupling : public tri_hp_helper {
		protected:
			FLT k_linear, k_torsion;
			FLT mass, I;
			int nboundary;
			Array<bdry_ins::force_coupling *,1> hp_ebdry;
			Array<eboundary_with_geometry<edge_bdry,naca> *,1> ebdry;
			FLT w_a, w_b, w_c, w_d;
			FLT w_a_tilda, w_b_tilda, w_c_tilda, w_d_tilda;
			TinyVector<FLT,5> k_a, k_b, k_c, k_d;
			/* This is to dump output from the generic boundaries */
			struct nullstream : std::ostream {
				nullstream(): std::ios(0), std::ostream(0) {}
			} ns;


		public:
			force_coupling(tri_hp_ins& xin) : tri_hp_helper(xin), w_a(0.0), w_b(0.0), w_c(0.0), w_d(0.0) {}
			force_coupling(const force_coupling &in_fc, tri_hp_ins& xin) : k_linear(in_fc.k_linear), k_torsion(in_fc.k_torsion),
				mass(in_fc.mass), I(in_fc.I), nboundary(in_fc.nboundary), tri_hp_helper(xin) {
				int i,j;

				hp_ebdry.resize(nboundary);
				ebdry.resize(nboundary);
				for (i = 0; i < nboundary; ++i) {
					for (j=0;j<x.nebd;++j) {
						if (x.ebdry(j)->idnum == in_fc.ebdry(i)->idnum) {
						goto found;
						}
					}
					std::cerr << "List of force boundaries is wrong" << std::endl;
					exit(1);

					found:
					if (!(hp_ebdry(i) = dynamic_cast<bdry_ins::force_coupling *>(x.hp_ebdry(j)))) {
						std:cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
						exit(1);
					}
					if (!(ebdry(i) = dynamic_cast<eboundary_with_geometry<edge_bdry,naca> *>(x.ebdry(j)))) {
						std::cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
						exit(1);
					}

				}
			}


			tri_hp_helper* create(tri_hp& xin) {return new force_coupling(*this, dynamic_cast<tri_hp_ins&>(xin));}

			void init(input_map& input, std::string idnty) {     
				std::string bdrys;
				std::istringstream bdryin;

				if (!input.get(x.gbl->idprefix +"_y0",w_a)) 
					input.getwdefault("y0",w_a,0.0);

				if (!input.get(x.gbl->idprefix +"_dydt0",w_b)) 
					input.getwdefault("dydt0",w_b,0.0);

				if (!input.get(x.gbl->idprefix +"_theta0",w_c)) 
					input.getwdefault("theta0",w_c,0.0);

				if (!input.get(x.gbl->idprefix +"_dthetadt0",w_d)) 
					input.getwdefault("dthetadt0",w_d,0.0);

				if (!input.get(x.gbl->idprefix +"_k_linear",k_linear)) 
					input.getwdefault("k_linear",k_linear,1.0);

				if (!input.get(x.gbl->idprefix +"_k_torsion",k_torsion)) 
					input.getwdefault("k_torsion",k_torsion,1.0);

				if (!input.get(x.gbl->idprefix +"_mass",mass)) 
					input.getwdefault("mass",mass,1.0);

				if (!input.get(x.gbl->idprefix +"_I",I)) 
					input.getwdefault("I",I,1.0);

				if (!input.get(x.gbl->idprefix +"_nboundary",nboundary))
					input.getwdefault("nboundary",nboundary,1);

				if (!input.getline(x.gbl->idprefix +"_force_boundaries",bdrys)) {
					if (!input.getline("force_boundaries",bdrys)) {
						std::cerr << "No boundary number list" << std::endl;
						exit(1);
					}
				}

				hp_ebdry.resize(nboundary);
				ebdry.resize(nboundary);
				bdryin.str(bdrys);
				int bnum;
				int j;
				for (int i = 0; i < nboundary; ++i) {
					bdryin >> bnum;
					for (j=0;j<x.nebd;++j) {
						if (x.ebdry(j)->idnum == bnum) {
							goto found;
						}
					}
					std::cerr << "List of force boundaries is wrong" << std::endl;
					exit(1);

					found:
						if (!(hp_ebdry(i) = dynamic_cast<bdry_ins::force_coupling *>(x.hp_ebdry(j)))) {
							std:cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
							exit(1);
						}
						if (!(ebdry(i) = dynamic_cast<eboundary_with_geometry<edge_bdry,naca> *>(x.ebdry(j)))) {
							std::cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
							exit(1);
						}

				}

			}

			void tadvance() {
				if (x.coarse_flag) return;

				int stage = x.gbl->substep +x.gbl->esdirk;
				if (stage > 0) {
					k_a[stage-1]=(w_a-w_a_tilda)/(1/x.gbl->adirk(stage-1,stage-1));
					k_b[stage-1]=(w_b-w_b_tilda)/(1/x.gbl->adirk(stage-1,stage-1)); 
					k_c[stage-1]=(w_c-w_c_tilda)/(1/x.gbl->adirk(stage-1,stage-1));
					k_d[stage-1]=(w_d-w_d_tilda)/(1/x.gbl->adirk(stage-1,stage-1));
				}

				if (x.gbl->substep == 0) {
					w_a_tilda=w_a;
					w_b_tilda=w_b;   
					w_c_tilda=w_c;
					w_d_tilda=w_d;  
				}

				for (int s=0;s<stage;++s) {
					w_a_tilda=w_a_tilda+(x.gbl->adirk(stage,s))*k_a[s];
					w_b_tilda=w_b_tilda+(x.gbl->adirk(stage,s))*k_b[s];
					w_c_tilda=w_c_tilda+(x.gbl->adirk(stage,s))*k_c[s];
					w_d_tilda=w_d_tilda+(x.gbl->adirk(stage,s))*k_d[s];
				}

				/* EXTRAPOLATE */
				if (stage  && x.gbl->dti > 0.0) {
					FLT constant =  x.gbl->cdirk(x.gbl->substep);
					w_a += constant*k_a[stage-1]; 
					w_b += constant*k_b[stage-1];
					w_c += constant*k_c[stage-1];
					w_d += constant*k_d[stage-1];

					/* FIX POSITIONS */
					FLT dy = constant*k_a[stage-1];
					FLT dtheta = constant*k_c[stage-1];
					FLT w_a_first = w_a -dy;
					TinyVector<FLT,2> ctr(0.0,w_a_first);
					TinyVector<FLT,2> vel(0.0,w_b);
					TinyVector<FLT,2> disp(0.0,dy);
					rigid_body(dtheta,w_d,ctr,disp,vel);	
				}


				return;
			}

			void mg_restrict() {
				if(x.coarse_level) {
					tri_hp *fmesh = dynamic_cast<tri_hp *>(x.fine);
					force_coupling *fine_helper = dynamic_cast<force_coupling *>(fmesh->helper);

					/* LOOP THROUGH POINTS TO TO CALCULATE POSITION OF COARSE POINTS  */
					int i,j,n,tind;
					for(i=0;i<x.npnt;++i) {
						tind = x.fcnnct(i).tri;

						for(n=0;n<x.ND;++n)
						x.pnts(i)(n) = 0.0;

						for(j=0;j<3;++j) {
						for(n=0;n<x.ND;++n)
							x.pnts(i)(n) += x.fcnnct(i).wt(j)*fmesh->pnts(fmesh->tri(tind).pnt(j))(n);
						}
					}


					for (int i=0; i<nboundary; ++i) {
						ebdry(i)->geometry_object.theta = fine_helper->ebdry(i)->geometry_object.theta; 
						ebdry(i)->geometry_object.pos = fine_helper->ebdry(i)->geometry_object.pos; 
					}
				}

				return;
			}


			void rigid_body(FLT dtheta, FLT omega, TinyVector<FLT,2> ctr, TinyVector<FLT,2> disp, TinyVector<FLT,2> vel) {
				TinyVector<FLT,2> dx;

				/* UPDATE MESH POSITION */
				FLT r,cost,sint; 
				FLT cosdt = cos(dtheta);    
				FLT sindt = sin(dtheta);    
				for (int i=0;i<x.npnt;++i) {
					dx = x.pnts(i) -ctr;
						r = sqrt(dx(0)*dx(0) +dx(1)*dx(1));
					cost = dx(0)/r;
					sint = dx(1)/r;
					x.pnts(i)(0) += -(r-r*cosdt)*cost -r*sindt*sint +disp(0);
					x.pnts(i)(1) += -(r-r*cosdt)*sint +r*sindt*cost +disp(1);						
				}

				for (int i=0; i<nboundary; ++i) {
					hp_ebdry(i)->set_ctr_rot(ctr);
					hp_ebdry(i)->set_vel(vel);
					hp_ebdry(i)->set_omega(w_d);
					ebdry(i)->geometry_object.theta += -dtheta; 
					ebdry(i)->geometry_object.pos(0) += disp(0) -0.25*(cos(w_c) -cos(w_c -dtheta)); 
					ebdry(i)->geometry_object.pos(1) += disp(1) -0.25*(sin(w_c) -sin(w_c -dtheta)); 

						/* CAN FIX ENDPOINTS TOO (NOT NECESSARY) */
//						v0 = x.seg(ebdry(i)->seg(0)).pnt(0);
//						FLT distance2 = (1-i)*0.75 -i*0.25;
//						FLT distance1 = (1-i)*(-0.25) +i*0.75;
//						x.pnts(v0)(0) = distance1*cos(w_c);
//						x.pnts(v0)(1) = w_a +distance1*sin(w_c);
//						v0 = x.seg(ebdry(i)->seg(ebdry(i)->nseg-1)).pnt(1);
//						x.pnts(v0)(0) = distance2*cos(w_c);
//						x.pnts(v0)(1) = w_a +distance2*sin(w_c);

					hp_ebdry(i)->curv_init();

				}
			}



			void update(int stage) {
				if (!x.coarse_flag) {

					/* GET FORCE & TORQUE */
					FLT force_y = 0.0;
					// FLT force_x = 0.0;
					FLT moment = 0.0;
					for (int i=0; i<nboundary; ++i) {                                  
					     /* FORCE BOUNDARY TO CALCULATE ALL FLUXES */
						hp_ebdry(i)->output(ns,tri_hp::tecplot);
						force_y += hp_ebdry(i)->diff_flux(1);
						moment += hp_ebdry(i)->moment;
					}   

					FLT dt = 1.0/x.gbl->dti;
					int stage = x.gbl->substep +x.gbl->esdirk;

					/* CALCULATE NEW POSITION */
					double w_a_first = w_a;
					double w_c_first = w_c;
					w_a=(w_a_tilda+((1/x.gbl->adirk(stage,stage)))*dt*(w_b_tilda+((1/x.gbl->adirk(stage,stage)))*dt*force_y))/(1+k_linear*(dt*dt*(1/x.gbl->adirk(stage,stage))*(1/x.gbl->adirk(stage,stage))));         // translational displacement
					w_b=(w_b_tilda-((1/x.gbl->adirk(stage,stage)))*dt*(w_a_tilda*k_linear-force_y))/(1+k_linear*(dt*dt*(1/x.gbl->adirk(stage,stage))*(1/x.gbl->adirk(stage,stage))));                                    // translational velocity
					w_c=(w_c_tilda+((1/x.gbl->adirk(stage,stage)))*dt*(w_d_tilda+((1/x.gbl->adirk(stage,stage)))*dt*moment))/(1+k_torsion*(dt*dt*(1/x.gbl->adirk(stage,stage))*(1/x.gbl->adirk(stage,stage))));         // rotational displacement
					w_d=(w_d_tilda-((1/x.gbl->adirk(stage,stage)))*dt*(w_c_tilda*k_torsion-moment))/(1+k_torsion*(dt*dt*(1/x.gbl->adirk(stage,stage))*(1/x.gbl->adirk(stage,stage))));                                   // rotational velocity


					/* FIX POSITIONS */
					FLT dy = w_a -w_a_first;
					FLT dtheta = w_c-w_c_first;
					TinyVector<FLT,2> ctr(0.0,w_a_first);
					TinyVector<FLT,2> vel(0.0,w_b);
					TinyVector<FLT,2> disp(0.0,dy);
					rigid_body(dtheta,w_d,ctr,disp,vel);	
				}
				return;
			}
	};
	
	
	
	class streamlines : public tri_hp_helper {
		protected:
			tri_hp_ins &x;
			int ndiv;
			TinyMatrix<FLT,2,tri_mesh::ND> rake_pts;
			int nterm_lines;
			Array<TinyMatrix<FLT,2,tri_mesh::ND>,1> term_lines;
			int maxtsteps;
			bool forward, backward;


		public:
			streamlines(tri_hp_ins& xin) : tri_hp_helper(xin), x(xin) {};
			streamlines(streamlines& tgt, tri_hp_ins& xin) : tri_hp_helper(xin), x(xin), ndiv(tgt.ndiv), nterm_lines(tgt.nterm_lines), maxtsteps(tgt.maxtsteps), forward(tgt.forward), backward(tgt.backward) {
				term_lines.resize(nterm_lines);
				term_lines = tgt.term_lines;
			}
			tri_hp_helper* create(tri_hp& xin) { return new streamlines(*this,dynamic_cast<tri_hp_ins&>(xin)); }
			void init(input_map& input, std::string idnty) {
			
				if (!input.get(x.gbl->idprefix +"_divisions",ndiv)) 
					input.getwdefault("divisions",ndiv,0);
					
				if (!input.get(x.gbl->idprefix +"_maxtsteps",maxtsteps)) 
					input.getwdefault("maxtsteps",maxtsteps,100000);

				if (!input.get(x.gbl->idprefix +"_pt0",&rake_pts(0,0),tri_mesh::ND)) {
					if (!input.get("pt0",&rake_pts(0,0),tri_mesh::ND)) {
						*x.gbl->log << "Couldn't find endpoint 0 of streamline rake" << std::endl;
						exit(1);
					}
				}
				
				if (!input.get(x.gbl->idprefix +"_pt1",&rake_pts(1,0),tri_mesh::ND)) {
					if (!input.get("pt1",&rake_pts(1,0),tri_mesh::ND)) {
						*x.gbl->log << "Couldn't find endpoint 1 of streamline rake" << std::endl;
						exit(1);
					}
				}
				
				/* Termination lines */
				if (!input.get(x.gbl->idprefix +"_nterm_line",nterm_lines)) 
					input.getwdefault("nterm_line",nterm_lines,0);

				term_lines.resize(nterm_lines);

				for (int n=0;n<nterm_lines;++n) {
					ostringstream nstr;
					nstr << "term_line" << n << std::endl;
					if (!input.get(x.gbl->idprefix +nstr.str() +"_pt0",&term_lines(n)(0,0),tri_mesh::ND)) {
						*x.gbl->log << "Couldn't find endpoint 0 of streamline rake" << std::endl;
						exit(1);
					}
					
					if (!input.get(x.gbl->idprefix +nstr.str() +"_pt1",&term_lines(n)(1,0),tri_mesh::ND)) {
						*x.gbl->log << "Couldn't find endpoint 1 of streamline rake" << std::endl;
						exit(1);
					}
				}
				
				if (!input.get(x.gbl->idprefix +"_forward",forward)) 
					input.getwdefault("forward",forward,true);
					
				if (!input.get(x.gbl->idprefix +"_backward",backward)) 
					input.getwdefault("backward",backward,false);
			}
			
			void output() {
				std::ofstream out;
				std:ostringstream filename;
				filename << "streamlines" << x.gbl->tstep << ".dat";
				out.open(filename.str().c_str());
			
				TinyVector<FLT,tri_mesh::ND> dx;
				dx(0) = (rake_pts(1,0)-rake_pts(0,0))/ndiv;
				dx(1) = (rake_pts(1,1)-rake_pts(0,1))/ndiv;
				
				TinyVector<FLT,tri_mesh::ND> ppt0, ppt1;
				for (int i=0;i<ndiv+1;++i) {
					ppt0(0) = rake_pts(0,0) +dx(0)*i;
					ppt0(1) = rake_pts(0,1) +dx(1)*i;
					
					int tind = -1;
					int counter;
					timeloop: for (counter = 0; counter < maxtsteps; ++counter) {
						
						const int nstage = 4;			
						Array<FLT,2> vels(nstage,x.NV);
						
						Array<FLT,1> vel = vels(0,Range::all());
						if (!x.ptprobe(ppt0,vel,tind)) break;
						out << ppt0(0) << ' ' << ppt0(1) << '\n';
						FLT dt = x.inscribedradius(tind)/(fabs(vel(0))+fabs(vel(1)));
						
						double a[nstage] = {0.0, 0.5, 0.5, 1.0};
						for (int s = 1; s < nstage;++s) {
							ppt1(0) = ppt0(0) +dt*a[s]*vels(s-1,0);
							ppt1(1) = ppt0(1) +dt*a[s]*vels(s-1,1);
							
							Array<FLT,1> vel = vels(s,Range::all());
							if (!x.ptprobe(ppt0,vel,tind)) goto exit;
						}
						ppt0(0) += dt/6.*(vels(0,0) +2.*vels(1,0) +2.*vels(2,0) +vels(3,0));
						ppt0(1) += dt/6.*(vels(0,1) +2.*vels(1,1) +2.*vels(2,1) +vels(3,1));
					}
					exit: out << counter << '\n';
					continue;
				}
			}

	};
	
	static bool tempflip = true;
	class static_particles : public streamlines {
		FLT rho, d;

		public:
			static_particles(tri_hp_ins& xin) : streamlines(xin) {};
			static_particles(static_particles& tgt, tri_hp_ins& xin) : streamlines(tgt), rho(tgt.rho), d(tgt.d) {}
			tri_hp_helper* create(tri_hp& xin) { return new static_particles(*this,dynamic_cast<tri_hp_ins&>(xin)); }
			
			void init(input_map& input, std::string idnty) {
				
				streamlines::init(input,idnty);
				
				if (!input.get(x.gbl->idprefix +"_particle_density",rho)) 
					input.getwdefault("particle_density",rho,1.0);
				
				if (!input.get(x.gbl->idprefix +"_particle_diameter",d)) 
					input.getwdefault("particle_diameter",d,1.0);
			}
			
			void tadvance() {
				d *= 2.0;
			}
				
			void output() {
				std::ofstream out;
				std:ostringstream filename;
				filename << "static_particles" << x.gbl->tstep << ".dat";
				out.precision(10);
				out.open(filename.str().c_str());
				
				TinyVector<FLT,tri_mesh::ND> dx;
				dx(0) = (rake_pts(1,0)-rake_pts(0,0))/ndiv;
				dx(1) = (rake_pts(1,1)-rake_pts(0,1))/ndiv;
				FLT re_over_v = x.gbl->rho*d/x.gbl->mu;
				FLT vol = M_PI*d*d*d/6.;
				FLT m = rho*vol;
				FLT buoy = (rho-x.gbl->rho)*vol*x.gbl->g/m;
				FLT rhoAo2m = 0.5*M_PI*d*d/4.*x.gbl->rho/m;

									
				TinyVector<FLT,tri_mesh::ND> ppt0, ppt1, vel0, vel1;
				for (int i=0;i<ndiv+1;++i) {
					ppt0(0) = rake_pts(0,0) +dx(0)*i;
					ppt0(1) = rake_pts(0,1) +dx(1)*i;
					vel0 = 0.0;

					int tind = -1;
					out << "Curve Begin\n";
					for (int counter = 0; counter < maxtsteps; ++counter) {
						const int nstage = 4;			
						Array<FLT,2> vels(nstage,x.NV);
						Array<FLT,2> fs(nstage,x.NV);
						Array<FLT,1> vel(x.NV);
						
						if (!x.ptprobe(ppt0,vel,tind)) break;
						out << ppt0(0) << ", " << ppt0(1) << '\n';
						
						TinyVector<FLT,2> vslip;
						vslip(0) = vel0(0) -vel(0);
						vslip(1) = vel0(1) -vel(1);
						
						FLT vmag = sqrt(vslip(0)*vslip(0) + vslip(1)*vslip(1));
						FLT re =vmag*re_over_v;
						FLT nonlinear = 0.125*pow(re,0.72);
						FLT cd = 24./re_over_v*(1.+nonlinear);  // Avoid dividing by v
						FLT dcdre = 24./re_over_v*nonlinear*.72/re;


						FLT dt = 1./((fabs(vel0(0))+fabs(vel0(1)))/x.inscribedradius(tind) +rhoAo2m*(cd +dcdre*re));
							
						vels(0,0) = vel0(0);
						vels(0,1) = vel0(1);
						fs(0,0) = -rhoAo2m*cd*vslip(0);
						fs(0,1) = -rhoAo2m*cd*vslip(1) -buoy;
												
						double a[nstage] = {0.0, 0.5, 0.5, 1.0};
						for (int s = 1; s < nstage;++s) {
							ppt1(0) = ppt0(0) +dt*a[s]*vels(s-1,0);
							ppt1(1) = ppt0(1) +dt*a[s]*vels(s-1,1);
							vel1(0) = vel0(0) +dt*a[s]*fs(s-1,0);
							vel1(1) = vel0(1) +dt*a[s]*fs(s-1,1);
							
							if (!x.ptprobe(ppt1,vel,tind)) goto nextpt;
							
							vslip(0) = vel1(0) -vel(0);
							vslip(1) = vel1(1) -vel(1);
							
							vmag = sqrt(vslip(0)*vslip(0) + vslip(1)*vslip(1));
							re =vmag*re_over_v;
							nonlinear = 0.125*pow(re,0.72);
							cd = 24./re_over_v*(1.+nonlinear);  // Avoid dividing by v		
							
							vels(s,0) = vel1(0);
							vels(s,1) = vel1(1);					
							fs(s,0) = -rhoAo2m*cd*vslip(0);
							fs(s,1) = -rhoAo2m*cd*vslip(1) -buoy;
							
						}
						ppt0(0) += dt/6.*(vels(0,0) +2.*vels(1,0) +2.*vels(2,0) +vels(3,0));
						ppt0(1) += dt/6.*(vels(0,1) +2.*vels(1,1) +2.*vels(2,1) +vels(3,1));
						vel0(0) += dt/6.*(fs(0,0) +2.*fs(1,0) +2.*fs(2,0) +fs(3,0));
						vel0(1) += dt/6.*(fs(0,1) +2.*fs(1,1) +2.*fs(2,1) +fs(3,1));
					}
					nextpt: continue;
				}
				
				out.close();
			}
		
	};	
	


	class helper_type {
		public:
			const static int ntypes = 6;
			enum ids {translating_drop,parameter_changer,unsteady_body_force,force_coupling,streamlines,static_particles};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
					if (!strcmp(nin,names[i])) return(i);
				return(-1);
			}
	};
	const char helper_type::names[ntypes][40] = {"translating_drop","parameter_changer","unsteady_body_force","force_coupling","streamlines","static_particles"};

	class ibc_type {
		public:
			const static int ntypes = 6;
			enum ids {freestream,sphere,accelerating,impinge,stokes_drop_gas,stokes_drop_liquid};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
					if (!strcmp(nin,names[i])) return(i);
				return(-1);
		}
	};
	const char ibc_type::names[ntypes][40] = {"freestream","sphere","accelerating","impinge","stokes_drop_gas","stokes_drop_liquid"};

}


init_bdry_cndtn *tri_hp_ins::getnewibc(std::string suffix, input_map& inmap) {
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
		case ibc_ins::ibc_type::sphere: {
			temp = new ibc_ins::sphere;
			break;
		}
		case ibc_ins::ibc_type::accelerating: {
			temp = new ibc_ins::accelerating;
			break;
		}
		case ibc_ins::ibc_type::impinge: {
			temp = new ibc_ins::impinge;
			break;
		}
		case ibc_ins::ibc_type::stokes_drop_gas: {
			temp = new ibc_ins::stokes_drop_gas;
			break;
		}
		case ibc_ins::ibc_type::stokes_drop_liquid: {
			temp = new ibc_ins::stokes_drop_liquid;
			break;
		}
		default: {
			return(tri_hp::getnewibc(suffix,inmap));
		}
	}
	temp->input(inmap,keyword);
	return(temp);
}

tri_hp_helper *tri_hp_ins::getnewhelper(input_map& inmap) {
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
		case ibc_ins::helper_type::translating_drop: {
			tri_hp_helper *temp = new ibc_ins::translating_drop(*this);
			return(temp);
		}
		case ibc_ins::helper_type::parameter_changer: {
			tri_hp_helper *temp = new ibc_ins::parameter_changer(*this);
			return(temp);
		}
		case ibc_ins::helper_type::unsteady_body_force: {
			tri_hp_helper *temp = new ibc_ins::unsteady_body_force(*this);
			return(temp);
		}
		case ibc_ins::helper_type::force_coupling: {
			tri_hp_helper *temp = new ibc_ins::force_coupling(*this);
			return(temp);
		}
		case ibc_ins::helper_type::streamlines: {
			tri_hp_helper *temp = new ibc_ins::streamlines(*this);
			return(temp);
		}
		case ibc_ins::helper_type::static_particles: {
			tri_hp_helper *temp = new ibc_ins::static_particles(*this);
			return(temp);
		}
		default: {
			return(tri_hp::getnewhelper(inmap));
		}
	}
}


#ifdef TWOLAYER

static FLT h = 2.0;

double f1(int n, double x, double y) { 
	FLT bf,re,g1,g2,n1,n2,q1,q2;

	/* FOR UIFACE TO BE 1 WITH D = 1, h = h/d */
	/* THETA DEFINED + CLOCKWISE */
	bf = mux[0]/(rhox[0]*(0.5 +(h-1)*rhox[1]/rhox[0])*sin(theta));
	body[0] = bf*sin(theta);
	body[1] = -bf*cos(theta);

	re = rhox[0]/mux[0];
	g1 = -bf*sin(theta);
	g2 = -bf*rhox[1]/rhox[0]*sin(theta);
	n1 = 1;
	q1 = 1;
	n2 = mux[1]/mux[0];
	q2 = rhox[1]/rhox[0];

	if (y < 1) {
		switch (n) {
			case(0):
				return(0.5*re*g1/n1*y*y +(re*g2*(1-h)-re*g1)*y);
			case(1):
				return(0.0);
			case(2):
				return(-bf*q1*cos(theta)*(y-h));
		}
	}
	else {
		switch (n) {
			case(0):
				return(0.5*re*g1/n2*y*y -re*g2*h/n2*y -0.5*g2*re/n2 +re*g2*h/n2 +re*g2*(1-h)-re*g1/2);
			case(1):
				return(0.0);
			case(2):
				return(-bf*q2*cos(theta)*(y-h));
		}
	}    

	return(0.0);
}
#endif

#ifdef ONELAYER
double f1(int n, double x, double y) { 
	FLT bf,re,n1,n2,n3,q1,q2,q3,h1,h2,h3;
	int mid,nonmid;

	/* FOR UIFACE TO BE 1 WITH D = 1, h = h/d */
	/* THETA DEFINED + CLOCKWISE */
	bf = 2.*mux[0]/(rhox[0]*sin(theta));
	body[0] = bf*sin(theta);
	body[1] = -bf*cos(theta);

	re = rhox[0]/mux[0];
	switch (n) {
		case(0):
			return(-(y+1.0)*(y-1.0) -1.0);
		case(1):
			return(0.0);
		case(2):
			return(-2.*cos(theta)/(sin(theta)*re)*y);
	}

	return(0.0);
}
#endif



#ifdef THREELAYER

double f1(int n, double x, double y) { 
	FLT bf,re,n1,n2,n3,q1,q2,q3,h1,h2,h3;
	int mid,nonmid;

	/* FOR UIFACE TO BE 1 WITH D = 1, h = h/d */
	/* THETA DEFINED + CLOCKWISE */
	/* FAIL PROOF TEST */
	if (fabs(mux[2] -mux[1]) < 1.0e-6) mid = 0;
	else if (fabs(mux[2] -mux[0]) < 1.0e-6) mid = 1;
	else mid = 2;
	nonmid = (mid+1)%3;

	bf = 2.*mux[nonmid]/(rhox[nonmid]*sin(theta));
	body[0] = bf*sin(theta);
	body[1] = -bf*cos(theta);

	re = rhox[nonmid]/mux[nonmid];

	h1 = 0.475;
	n1 = 1;
	q1 = 1;

	h2 = 0.525;
	n2 = mux[mid]/mux[nonmid];
	q2 = rhox[mid]/rhox[nonmid];

	h3 = 1.0;
	n3 = mux[nonmid]/mux[nonmid];
	q3 = rhox[nonmid]/rhox[nonmid];

	switch (n) {
		case(0):
			if (y <= 0.475) {
				double c1 = 2.0*h3/n1;
				double c2 = 0.0;
				return(-1./n1*y*y +c1*y +c2);
			}
			else if (y <= 0.525) {
				double c1 = 2.0*h3/n2;
				double c2 = h1*(h1*n1-h1*n2-2*h3*n1+2*h3*n2)/(n1*n2);
				return(-1./n2*y*y +c1*y +c2);
			}
			else {
				double c1 = 2*h3/n3;
				double c2 = (-h2*h2*n1*n3 +h2*h2*n1*n2 +h1*h1*n1*n3-h1*h1*n2*n3-2*h2*h3*n1*n2+2*h2*h3*n1*n3
					-2*h1*h3*n1*n3 +2*h1*h3*n2*n3)/(n1*n2*n3);
				return(-1./n3*y*y +c1*y +c2);
			}
		case(1):
			return(0.0);
		case(2):
			return(-2.*cos(theta)/(sin(theta)*re)*y);
	}

	return(0.0);
}
#endif




