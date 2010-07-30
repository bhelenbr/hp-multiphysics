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
					sim::abort(__LINE__,__FILE__,&std::cerr);
				}

				keyword = blkname +"_liquid";
				if (!blockdata.get(keyword,val)) { 
					std::cerr << "couldn't find identity of liquid block" << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);
				}

				FLT mu_l;
				keyword = val +"_mu";
				if (!blockdata.get(keyword,mu_l)) {
					std::cerr << "couldn't find mu of liquid" << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);
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
					sim::abort(__LINE__,__FILE__,&std::cerr);
				}

				keyword = blkname +"_mu";
				if (!blockdata.get(keyword,mu_l)) {
					std::cerr << "couldn't find mu of liquid" << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);
				}
				keyword = blkname +"_liquid_bdry";
				if (!blockdata.get(keyword,val)) { 
					std::cerr << "couldn't find identity of liquid boundary" << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);
				}

				keyword = val +"_sigma";
				if (!blockdata.get(keyword,sigma)) {
					std::cerr << "couldn't find sigma" << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);
				}				

				keyword = blkname +"_gas";
				if (!blockdata.get(keyword,val)) { 
					std::cerr << "couldn't find identity of gas block" << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);
				}

				FLT mu_g;
				keyword = val +"_mu";
				if (!blockdata.get(keyword,mu_g)) {
					std::cerr << "couldn't find mu of gas " << keyword << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);
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
					if ((surf = dynamic_cast<bdry_ins::surface *>(x.hp_ebdry(bnum)))) break;

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
							sim::abort(__LINE__,__FILE__,&std::cerr);
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
			tri_hp_ins& x;
			
			/* Physical data for solid body */
			bool horizontal, vertical, rotational;
			TinyVector<FLT,tri_mesh::ND> k_linear;
			FLT k_torsion;
			FLT mass, I;
			
			/* Pointers to boundary information */
			int nboundary;
			Array<bdry_ins::generic *,1> hp_ebdry;
			Array<bdry_ins::rigid *, 1> hp_ebdry_rigid;
			Array<rigid_movement_interface2D *,1> ebdry_rigid;
			Array<edge_bdry *,1> ebdry;
			
			/* Data for DIRK scheme */
			TinyVector<FLT,6> w;  // x,xdot,y,ydot,theta,thetadot 
			TinyVector<FLT,6> w_tilda;  //  x,xdot,y,ydot,theta,thetadot 
			TinyVector<TinyVector<FLT,6>,5> ki;
			
			/* Residual */
			TinyVector<FLT,6> res;
			
			/* For multigrid way */
			bool isfrst;
			TinyMatrix<FLT,6,6> J;
			TinyVector<int,12> ipiv;
			TinyVector<FLT,6> mg_res0;
			TinyVector<FLT,6> mg_w0;
			TinyVector<TinyVector<FLT,6>,3> dres;
			
			/* For full jacobian way (not working) */
			int jacobian_start;

			
			/* This is to dump output from the generic boundaries */
			struct nullstream : std::ostream {
				nullstream(): std::ios(0), std::ostream(0) {}
			} ns;

		public:
			force_coupling(tri_hp_ins& xin) : tri_hp_helper(xin), x(xin), horizontal(false), vertical(false), rotational(false), w(0.0) {}
			force_coupling(const force_coupling &in_fc, tri_hp_ins& xin) : x(xin), horizontal(in_fc.horizontal), vertical(in_fc.vertical),  
				rotational(in_fc.rotational), k_linear(in_fc.k_linear), k_torsion(in_fc.k_torsion),
				mass(in_fc.mass), I(in_fc.I), nboundary(in_fc.nboundary), tri_hp_helper(xin) {
				int i,j;

				hp_ebdry.resize(nboundary);
				hp_ebdry_rigid.resize(nboundary);
				ebdry.resize(nboundary);
				ebdry_rigid.resize(nboundary);
				for (i = 0; i < nboundary; ++i) {
					for (j=0;j<x.nebd;++j) {
						if (x.ebdry(j)->idnum == in_fc.ebdry(i)->idnum) {
						goto found;
						}
					}
					std::cerr << "List of force boundaries is wrong" << std::endl;
					sim::abort(__LINE__,__FILE__,&std::cerr);

					found:
					hp_ebdry(i) = dynamic_cast<bdry_ins::generic *>(x.hp_ebdry(j));
					ebdry(i) = x.ebdry(j);
					if (!(hp_ebdry_rigid(i) = dynamic_cast<bdry_ins::rigid *>(x.hp_ebdry(j)))) {
						std::cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
						sim::abort(__LINE__,__FILE__,&std::cerr);
					}
					if (!(ebdry_rigid(i) = dynamic_cast<rigid_movement_interface2D *>(x.ebdry(j)))) {
						std::cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
						sim::abort(__LINE__,__FILE__,&std::cerr);
					}
				}
			}


			tri_hp_helper* create(tri_hp& xin) {return new force_coupling(*this, dynamic_cast<tri_hp_ins&>(xin));}
			
			/* DECOUPLING THE FLOW & BOUNDARY MOVEMENT EQUATIONS */
			int dofs(int dofs) {
				return(0);  // Nothing in Jacobian
				/* For full jacobian (not working yet) */
				jacobian_start = dofs;
				return(2*(vertical+horizontal+rotational));
			}
			void non_sparse(Array<int,1> &nnzero) {		
				return;
				/* For full jacobian (not working yet) */
				if (x.mmovement == tri_hp::coupled_deformable) {
					const int vdofs = x.NV +x.ND;
					const int sm=basis::tri(x.log2p)->sm();
					const int begin_seg = x.npnt*vdofs;
					
					for (int i=0; i<nboundary; ++i) {
						if (hp_ebdry(i)->is_curved())
							nnzero(Range(hp_ebdry(i)->jacobian_start,hp_ebdry(i)->jacobian_start+hp_ebdry(i)->base.nseg*sm*tri_mesh::ND-1)) = 2*(vertical+horizontal+rotational);

						for(int j=0;j<ebdry(i)->nseg;++j) {
							int sind = ebdry(i)->seg(j);
							int v0 = x.seg(sind).pnt(0);

							nnzero(Range(v0*vdofs,(v0+1)*vdofs-1)) += 2*(vertical+horizontal+rotational);
							nnzero(Range(begin_seg+sind*x.NV*sm,begin_seg+(sind+1)*x.NV*sm-1)) += 2*(vertical+horizontal+rotational);
						}
						int v0 = x.seg(ebdry(i)->seg(ebdry(i)->nseg-1)).pnt(1);
						for(int n=0;n<vdofs;++n) {
							nnzero(v0*vdofs+n) += 2*(vertical+horizontal+rotational);
						}
					}
				}
				
				if (x.mmovement == tri_hp::coupled_rigid) {
					nnzero(Range(0,x.npnt*x.NV-1)) += 2*(vertical+horizontal+rotational);
					/* FIXME: NOT SURE WHAT TO DO WITH HIGH-ORDER MODES ON CURVED BOUNDARY YET */
				}
			}
				
			/* FUll jacobian is not working yet */
			void jacobian() {}
			void jacobian_dirichlet() {}

			void init(input_map& input, std::string idnty) {     
				std::string bdrys;
				std::istringstream bdryin;

				if (!input.get(x.gbl->idprefix +"_horizontal",horizontal)) 
					input.getwdefault("horizontal",horizontal,false);
					
				if (!input.get(x.gbl->idprefix +"_vertical",vertical)) 
					input.getwdefault("vertical",vertical,false);

				if (!input.get(x.gbl->idprefix +"_rotational",rotational)) 
					input.getwdefault("rotational",rotational,false);
								
				if (horizontal) {
					if (!input.get(x.gbl->idprefix +"_x0",w(0))) 
					input.getwdefault("x0",w(0),0.0);

					if (!input.get(x.gbl->idprefix +"_dx0dt",w(1)))
						input.getwdefault("dx0dt",w(1),0.0);
						
					if (!input.get(x.gbl->idprefix +"_k_linear0",k_linear(0))) 
						input.getwdefault("k_linear0",k_linear(0),1.0);
				}
				
				if (vertical) {
					if (!input.get(x.gbl->idprefix +"_x1",w(2))) 
					input.getwdefault("x1",w(2),0.0);

					if (!input.get(x.gbl->idprefix +"_dx1dt",w(3)))
						input.getwdefault("dx1dt",w(3),0.0);
						
					if (!input.get(x.gbl->idprefix +"_k_linear1",k_linear(0))) 
						input.getwdefault("k_linear1",k_linear(1),1.0);
				}
				
				if (horizontal || vertical) {
					if (!input.get(x.gbl->idprefix +"_mass",mass)) 
						input.getwdefault("mass",mass,1.0);
				}
				
				if (rotational) {
					if (!input.get(x.gbl->idprefix +"_theta",w(4))) 
					input.getwdefault("theta",w(4),0.0);

					if (!input.get(x.gbl->idprefix +"_dthetadt",w(5)))
						input.getwdefault("dthetadt",w(5),0.0);
						
					if (!input.get(x.gbl->idprefix +"_I",I)) 
						input.getwdefault("I",I,1.0);
					
					if (!input.get(x.gbl->idprefix +"_k_torsion",k_torsion)) 
						input.getwdefault("k_torsion",k_torsion,1.0);
				}
				w_tilda = w;		

				if (!input.get(x.gbl->idprefix +"_nboundary",nboundary))
					input.getwdefault("nboundary",nboundary,1);

				if (!input.getline(x.gbl->idprefix +"_force_boundaries",bdrys)) {
					if (!input.getline("force_boundaries",bdrys)) {
						std::cerr << "No boundary number list" << std::endl;
						sim::abort(__LINE__,__FILE__,&std::cerr);
					}
				}

				hp_ebdry.resize(nboundary);
				ebdry.resize(nboundary);
				hp_ebdry_rigid.resize(nboundary);
				ebdry_rigid.resize(nboundary);
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
					sim::abort(__LINE__,__FILE__,&std::cerr);

					found:
					hp_ebdry(i) = dynamic_cast<bdry_ins::generic *>(x.hp_ebdry(j));
						ebdry(i) = x.ebdry(j);
						if (!(hp_ebdry_rigid(i) = dynamic_cast<bdry_ins::rigid *>(x.hp_ebdry(j)))) {
							std::cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
							sim::abort(__LINE__,__FILE__,&std::cerr);
						}
						if (!(ebdry_rigid(i) = dynamic_cast<rigid_movement_interface2D *>(x.ebdry(j)))) {
							std::cerr << "Boundary in list of force boundaries is of wrong type" << std::endl;
							sim::abort(__LINE__,__FILE__,&std::cerr);
						}

				}

			}

			void tadvance() {
				if (x.coarse_flag) return;

				int stage = x.gbl->substep +x.gbl->esdirk;
				if (stage > 0)
					ki(stage-1) = (w-w_tilda)/(1/x.gbl->adirk(stage-1,stage-1));

				if (x.gbl->substep == 0)
					w_tilda = w;

				for (int s=0;s<stage;++s) {
					w_tilda = w_tilda +x.gbl->adirk(stage,s)*ki(s);
				}

//				/* EXTRAPOLATE */
//				if (stage  && x.gbl->dti > 0.0) {
//					FLT constant =  x.gbl->cdirk(x.gbl->substep);
//					w += constant*ki(stage-1); 
//					
//					/* FIX POSITIONS */
//					FLT dx = constant*ki(stage-1)(0);
//					FLT dy = constant*ki(stage-1)(2);
//					FLT dtheta = constant*ki(stage-1)(4);
//					FLT x_first = w(0) -dx;
//					FLT y_first = w(2) -dy;
//					TinyVector<FLT,2> ctr(x_first,y_first);
//					TinyVector<FLT,2> vel(w(1),w(3));
//					TinyVector<FLT,2> disp(dx,dy);
//					rigid_body(dtheta,w(5),ctr,disp,vel);	
//				}


				return;
			}

			void mg_restrict() {
				isfrst = true;
				
				if(x.coarse_level) {
					isfrst = true;
					tri_hp *fmesh = dynamic_cast<tri_hp *>(x.fine);
					force_coupling *fine_helper = dynamic_cast<force_coupling *>(fmesh->helper);

					if (x.mmovement != tri_hp::coupled_deformable) {
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
					}
				
					w = fine_helper->w;
					w_tilda = fine_helper->w_tilda;  //  x,xdot,y,ydot,theta,thetadot 
					ki = fine_helper->ki;
					mg_w0 = w;
					mg_res0 = fine_helper->res;
					for (int i=0; i<nboundary; ++i) {
						ebdry_rigid(i)->theta = fine_helper->ebdry_rigid(i)->theta; 
						ebdry_rigid(i)->pos = fine_helper->ebdry_rigid(i)->pos; 
						hp_ebdry_rigid(i)->ctr = 	fine_helper->hp_ebdry_rigid(i)->ctr;
						hp_ebdry_rigid(i)->vel = 	fine_helper->hp_ebdry_rigid(i)->vel;
						hp_ebdry_rigid(i)->omega = fine_helper->hp_ebdry_rigid(i)->omega;
					}					
				}
				else {
					mg_res0 = res;
				}

				return;
			}
			
			void mg_prolongate() {
				if(!x.coarse_level) {
					return;
				}
				/* CALCULATE CORRECTIONS */
				mg_w0 -= w;
				
				/* LOOP THROUGH FINE VERTICES    */
				/* TO DETERMINE CHANGE IN SOLUTION */
				tri_hp *fmesh = dynamic_cast<tri_hp *>(x.fine);
				force_coupling *fine_helper = dynamic_cast<force_coupling *>(fmesh->helper);
				fine_helper->w -= mg_w0;
			}

			void rigid_body(TinyVector<FLT,6> w0, TinyVector<FLT,6> dw) {
				FLT dtheta = dw(4);
				FLT omega = w0(5)+dw(5);
				TinyVector<FLT,2> ctr(w0(0),w0(2));
				TinyVector<FLT,2> disp(dw(0),dw(2));
				TinyVector<FLT,2> vel(w0(1)+dw(1),w0(3)+dw(3));				
				TinyVector<FLT,2> dx;

				/* UPDATE MESH POSITION */
				FLT cosdt = cos(dtheta);    
				FLT sindt = sin(dtheta);
				if (x.mmovement != tri_hp::coupled_deformable) {
					FLT r,cost,sint; 
    
					for (int i=0;i<x.npnt;++i) {
						dx = x.pnts(i) -ctr;
						r = sqrt(dx(0)*dx(0) +dx(1)*dx(1));
						cost = dx(0)/r;
						sint = dx(1)/r;
						x.pnts(i)(0) += -(r-r*cosdt)*cost -r*sindt*sint +disp(0);
						x.pnts(i)(1) += -(r-r*cosdt)*sint +r*sindt*cost +disp(1);						
					}
				}

				for (int i=0; i<nboundary; ++i) {
					hp_ebdry_rigid(i)->ctr = ctr+disp;
					hp_ebdry_rigid(i)->vel = vel;
					hp_ebdry_rigid(i)->omega = omega;
					
					FLT r = 0.25;
					FLT cost = -cos(-ebdry_rigid(i)->theta);
					FLT sint = -sin(-ebdry_rigid(i)->theta);
					ebdry_rigid(i)->theta -= dtheta; 
					/* FIX ME HACK FOR 1/4 CHORD */
					ebdry_rigid(i)->pos(0) += -(r-r*cosdt)*cost -r*sindt*sint +disp(0);
					ebdry_rigid(i)->pos(1) += -(r-r*cosdt)*sint +r*sindt*cost +disp(1);	

						/* CAN FIX ENDPOINTS TOO (NOT NECESSARY) */
//						v0 = x.seg(ebdry(i)->seg(0)).pnt(0);
//						FLT distance2 = (1-i)*0.75 -i*0.25;
//						FLT distance1 = (1-i)*(-0.25) +i*0.75;
//						x.pnts(v0)(0) = distance1*cos(theta;
//						x.pnts(v0)(1) = w(0) +distance1*sin(theta);
//						v0 = x.seg(ebdry(i)->seg(ebdry(i)->nseg-1)).pnt(1);
//						x.pnts(v0)(0) = distance2*cos(theta);
//						x.pnts(v0)(1) = w(0) +distance2*sin(theta);

					hp_ebdry(i)->curv_init();

				}
			}
			
			void rsdl(int s = 0) {
				int stage = x.gbl->substep +x.gbl->esdirk;
				FLT bd0 = x.gbl->dti*x.gbl->adirk(stage,stage);
					
				/* GET FORCE & TORQUE */
				TinyVector<FLT,tri_mesh::ND> force = 0.0;
				FLT moment = 0.0;
				for (int i=0; i<nboundary; ++i) { 
					/* FORCE BOUNDARY TO CALCULATE ALL FLUXES */
					hp_ebdry(i)->output(*x.gbl->log,tri_hp::tecplot);
					force(0) -= hp_ebdry(i)->diff_flux(0);
					force(1) -= hp_ebdry(i)->diff_flux(1);					
					moment += hp_ebdry(i)->moment;
				}
//				std::cout << w << ' ' << w_tilda << std::endl;
//				std::cout << force << std::endl;

				double L = 0.05;  // HACK HACK HACK FIX ME FOR GETTING HYDROSTATIC PRESSURE CORRECT
				force(1) += -x.gbl->rho*x.gbl->g*w(2)*L; // -mass*x.gbl->g;				
				
				if (horizontal) {
					// translational displacement
					res(0) = bd0*(w(0) -w_tilda(0)) -w(1);
					// translational velocity
					res(1) = mass*bd0*(w(1) -w_tilda(1)) -force(0) +k_linear(0)*w(0);
				} 
				else {
					res(0) = 0.0;
					res(1) = 0.0;
				}
				
				if (vertical) {
					// translational displacement
					res(2) = bd0*(w(2) -w_tilda(2)) -w(3);
					// translational velocity
					res(3) = mass*bd0*(w(3) -w_tilda(3)) -force(1) +k_linear(1)*w(2);
				} 
				else {
					res(2) = 0.0;
					res(3) = 0.0;
				}		
				
				if (rotational) {
					// rotational displacement
					res(4) = bd0*(w(4) -w_tilda(4)) -w(5);
					// translational velocity
					res(5) = I*bd0*(w(5) -w_tilda(5)) -moment +k_torsion*w(4);
				} 
				else {
					res(4) = 0.0;
					res(5) = 0.0;
				}
				if(x.coarse_flag) {
					/* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
					if(isfrst) {
						isfrst = false;
						dres(x.log2p) = x.fadd*mg_res0 -res;
					}
					res += dres(x.log2p); 
				}
					
			}
			
			void setup_preconditioner() {
				FLT delta = 1.0e-6;
				
				rsdl(0);
				
				TinyVector<FLT,6> res0, dw;
				res0 = res;
				dw = 0.0;
				J = 0.0;
				 				
				/* PERTURB VARIABLES AND SEE HOW FORCES & MOMENT CHANGES */
				if (horizontal) {
					/* Displacement */
					dw(0) = delta;
					rigid_body(w,dw);
					w(0) += delta;	
					rsdl(0);
					for(int i=0;i<6;++i)
						J(i,0) = (res(i)-res0(i))/delta;
					/* Restore */
					dw(0) = -delta;
					rigid_body(w,dw);
					w(0) -= delta;
					dw(0) = 0.0;
					
					/* Velocity */
					dw(1) = delta;
					rigid_body(w,dw);	
					w(1) += delta;
					rsdl(0);
					for(int i=0;i<6;++i)
						J(i,1) = (res(i)-res0(i))/delta;
					/* Restore */
					dw(1) -= delta;
					rigid_body(w,dw);	
					w(1) -= delta;
					dw(1) = 0.0;
				}
				
				if (vertical) {
					/* Displacement */
					dw(2) = delta;
					rigid_body(w,dw);
					w(2) += delta;	
					rsdl(0);
					for(int i=0;i<6;++i)
						J(i,2) = (res(i)-res0(i))/delta;
					/* Restore */
					dw(2) = -delta;
					rigid_body(w,dw);
					w(2) -= delta;
					dw(2) = 0.0;

					/* Velocity */
					dw(3) = delta;
					rigid_body(w,dw);	
					w(3) += delta;
					rsdl(0);
					for(int i=0;i<6;++i)
						J(i,3) = (res(i)-res0(i))/delta;
					/* Restore */
					dw(3) = -delta;
					rigid_body(w,dw);
					w(3) -= delta;
					dw(3) = 0.0;
				}

				if (rotational) {
					/* Rotation */
					dw(4) = delta;
					rigid_body(w,dw);	
					w(4) += delta;
					rsdl(0);
					for(int i=0;i<6;++i)
						J(i,4) = (res(i)-res0(i))/delta;
					/* Restore */
					dw(4) = -delta;
					rigid_body(w,dw);
					w(4) -= delta;
					dw(4) = 0.0;
					
					/* Angular velocity */
					dw(5) = delta;
					rigid_body(w,dw);	
					w(5) += delta;
					rsdl(0);
					for(int i=0;i<6;++i)
						J(i,5) = (res(i)-res0(i))/delta;
					/* Restore */
					dw(5) = -delta;
					rigid_body(w,dw);
					w(5) -= delta;
					dw(5) = 0.0;
				}
				
				/* Constrain motions */
				if (!horizontal) {
					for(int i=0;i<6;++i) {
						J(0,i) = 0.0;
						J(1,i) = 0.0;
					}
					J(0,0) = 1.0;
					J(1,1) = 1.0;
				}
				if (!vertical) {
					for(int i=0;i<6;++i) {
						J(2,i) = 0.0;
						J(3,i) = 0.0;
					}	
					J(2,2) = 1.0;
					J(3,3) = 1.0;
				}
				if (!rotational) {
					for(int i=0;i<6;++i) {
						J(4,i) = 0.0;
						J(5,i) = 0.0;
					}
					J(4,4) = 1.0;
					J(5,5) = 1.0;
				}
				
				/* Invert J */
				int info;
				int size = 6;				
				GETRF(size,size,J.data(),size,ipiv.data(),info);
				if (info) {
					*x.gbl->log << "Error inverting Jacobian for solid body coupling " << info << std::endl;
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
				
			}

			void update(int stage) {
				if (!x.coarse_flag) {
				
					/* Calculate residual */
					rsdl(0);
					
					int info;
					char trans[] = "T";	
					int size = 6;
					GETRS(trans,size,1,J.data(),size,ipiv.data(),res.data(),size,info);
					if (info) {
						*x.gbl->log << "Error with Jacobian product for solid body coupling" << std::endl;
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}
					
					res *= -1.0;
					
					/* Add correction */
					rigid_body(w, res);
					w += res;
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
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}
				}
				
				if (!input.get(x.gbl->idprefix +"_pt1",&rake_pts(1,0),tri_mesh::ND)) {
					if (!input.get("pt1",&rake_pts(1,0),tri_mesh::ND)) {
						*x.gbl->log << "Couldn't find endpoint 1 of streamline rake" << std::endl;
						sim::abort(__LINE__,__FILE__,x.gbl->log);
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
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}
					
					if (!input.get(x.gbl->idprefix +nstr.str() +"_pt1",&term_lines(n)(1,0),tri_mesh::ND)) {
						*x.gbl->log << "Couldn't find endpoint 1 of streamline rake" << std::endl;
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}
				}
				
				if (!input.get(x.gbl->idprefix +"_forward",forward)) 
					input.getwdefault("forward",forward,true);
					
				if (!input.get(x.gbl->idprefix +"_backward",backward)) 
					input.getwdefault("backward",backward,false);
			}
			
			void output() {
				std::ofstream out;
				std::ostringstream filename;
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
					for (counter = 0; counter < maxtsteps; ++counter) {
						
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
				std::ostringstream filename;
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




