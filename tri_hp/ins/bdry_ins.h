/*
 *  ins_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _bdry_ins_h_
#define _bdry_ins_h_


#include "tri_hp_ins.h"
#include "../hp_boundary.h"
#include <myblas.h>
#include <blitz/array.h>
#include <symbolic_function.h>

// #define DROP

using namespace blitz;

//#define DETAILED_DT
//#define DETAILED_MINV

//#define L2_ERROR

namespace bdry_ins {

	class generic : public hp_edge_bdry {
		protected:
			tri_hp_ins &x;
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
				/* CONTINUITY */
				flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
				/* X&Y MOMENTUM */
#ifdef INERTIALESS
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#else
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = flx(x.NV-1)*u(n) +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#endif			
				/* EVERYTHING ELSE */
				for (int n=tri_mesh::ND;n<x.NV-1;++n)
					flx(n) = flx(x.NV-1)*u(n);
				return;
			}

		public:
			Array<FLT,1> total_flux,diff_flux,conv_flux;
			FLT circumference,moment,convect,circulation;

		public:
			generic(tri_hp_ins &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "generic";}
			generic(const generic& inbdry, tri_hp_ins &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin) {
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);
				total_flux = 0.0;
				diff_flux = 0.0;
				conv_flux = 0.0;
			}
			generic* create(tri_hp& xin, edge_bdry &bin) const {return new generic(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in) {
				hp_edge_bdry::init(inmap,gbl_in);
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);
				total_flux = 0.0;
				diff_flux = 0.0;
				conv_flux = 0.0;
			}
			void output(const std::string& fname, tri_hp::filetype typ,int tlvl = 0);
	};


	class inflow : public generic {  
		public:
			inflow(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "inflow";}
			inflow(const inflow& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			inflow* create(tri_hp& xin, edge_bdry &bin) const {return new inflow(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in) {
				generic::init(inmap,gbl_in);
				for (int n=0;n<x.NV-1;++n) {
					essential_indices.push_back(n);
					type[n] = essential;
				}
			}
	};


	class euler : public generic {
		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
				Array<FLT,1> ub(x.NV);

				for(int n=0;n<x.NV;++n)
					ub(n) = ibc->f(n,xpt,x.gbl->time);

				flx(x.NV-1) = x.gbl->rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1));

				/* X&Y MOMENTUM */
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = flx(x.NV-1)*ub(n) +u(x.NV-1)*norm(n);

				/* EVERYTHING ELSE */
				for (int n=tri_mesh::ND;n<x.NV-1;++n)
					flx(n) = flx(x.NV-1)*ub(n);

				return;
			}
		public:
			euler(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "euler";}
			euler(const euler& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			euler* create(tri_hp& xin, edge_bdry &bin) const {return new euler(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
	};
	
	class rigid : public init_bdry_cndtn {
		public:
			FLT omega;
			TinyVector<FLT,2> vel, ctr;
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				FLT r,cost,sint,v; 
				TinyVector<FLT,2> dx;
				
				dx = x-ctr;
				r = sqrt(dx(0)*dx(0) +dx(1)*dx(1));
				cost = dx(0)/r;
				sint = dx(1)/r;
				if(n == 0)
					v = -r*omega*sint +vel(0);
				else if (n == 1)
					v = r*omega*cost +vel(1);
				else
					v = 0.0;
				
				return(v);
			}
			rigid() : omega(0.0), vel(0.0), ctr(0.0) {}
			void init(input_map &inmap,std::string idnty) {
				inmap.getwdefault(idnty+"_omega",omega,0.0);
				FLT dflt[2] = {0.0, 0.0};
				inmap.getwdefault(idnty+"_center",ctr.data(),2,dflt);
				inmap.getwdefault(idnty+"_velocity",vel.data(),2,dflt);
			}
	};		
	
	class force_coupling : public inflow, public rigid {
		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
				for (int n=0;n<x.NV;++n)
					flx(n) = 0.0;
				return;
			}
			
		public:
			force_coupling(tri_hp_ins &xin, edge_bdry &bin) : inflow(xin,bin), rigid() {mytype = "force_coupling"; /* ibc = this; */}
			force_coupling(const force_coupling& inbdry, tri_hp_ins &xin, edge_bdry &bin) : inflow(inbdry,xin,bin), rigid(inbdry) {/*ibc = this;*/}
			force_coupling* create(tri_hp& xin, edge_bdry &bin) const {return new force_coupling(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in) {
				inflow::init(inmap,gbl_in);
				rigid::init(inmap,base.idprefix);
				report_flag = true;
			}
			void tadvance() {hp_edge_bdry::tadvance();}
			void update(int stage) {
				if (!x.coarse_flag) setvalues(this,essential_indices);
			}
			void output(const std::string& fname, tri_hp::filetype typ,int tlvl = 0) {
				inflow::output(fname,typ,tlvl);
				if (typ == tri_hp::text) {
					*x.gbl->log << omega << ' ' << vel << ' ' << ctr << '\n';
				}
		}
	};
	
	class friction_slip : public generic, public rigid {
		protected:
			FLT slip_length;
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
				Array<FLT,1> urb(x.NV), ub(x.NV);

				/* RIGID COMPONENT */
				for(int n=0;n<x.NV;++n)
					urb(n) = f(n,xpt,x.gbl->time);
	
				/* Additional Components (for belts, sliding objects etc..) */
			  /* Calculate position in object's frame */
				/* i,j are world, hat's are in objects frame */
				/* i = cos(theta) hat{i} -sin(theta) hat{j} */
				/* j = sin(theta) hat{i} +cos(theta) hat{j} */
				
//				xpt -= my_ibc.ctr;
//				FLT temp = xpt(0)*cos(my_ibc.theta) +xpt(1)*sin(my_ibc.theta);
//				xpt(1) = -xpt(0)*sin(my_ibc.theta) +xpt(1)*cos(my_ibc.theta);
//				xpt(0) = temp;
				
				for(int n=0;n<x.NV;++n)
					ub(n) = ibc->f(n,xpt,x.gbl->time);
					
				/* Rotate back to cartesian */
				/* \hat{i} = cos(theta) i +sin(theta) j */
				/* \hat{j} = -sin(theta) i +cos(theta) j */
//				temp = ub(0)*cos(my_ibc.theta) -ub(1)*sin(my_ibc.theta);
//				ub(1) = ub(0)*sin(my_ibc.theta) +ub(1)*cos(my_ibc.theta);
//				ub(0) = temp;
				
				ub += urb;

				// FLT un = (ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1);
				FLT un = 0.0;			// No flux through the surface	
				TinyVector<FLT,tri_mesh::ND> tang(-norm(1),norm(0));
				FLT slip_stress = ((ub(0) -u(0))*tang(0) +(ub(1) -u(1))*tang(1))*x.gbl->mu/(slip_length);
				
				flx(x.NV-1) = x.gbl->rho*un;

				/* X&Y MOMENTUM */
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = flx(x.NV-1)*ub(n) +u(x.NV-1)*norm(n) -slip_stress*tang(n);

				/* EVERYTHING ELSE */
				for (int n=tri_mesh::ND;n<x.NV-1;++n)
					flx(n) = flx(x.NV-1)*ub(n);

				return;
			}
		public:
			friction_slip(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin), rigid(), slip_length(0.0) {mytype = "friction_slip";}
			friction_slip(const friction_slip& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin), rigid(inbdry), slip_length(inbdry.slip_length) {}
			friction_slip* create(tri_hp& xin, edge_bdry &bin) const {return new friction_slip(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in) {
				generic::init(inmap,gbl_in);
				rigid::init(inmap,base.idprefix);
				std::string keyword = base.idprefix +"_sliplength";
				if (!inmap.get(keyword,slip_length)) {
					*x.gbl->log << "Couldn't find slip length for " << keyword << std::endl;
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
			}
//			void input(ifstream& fin,tri_hp::filetype typ,int tlvl) {
//				generic::input(fin,typ,tlvl);
//				if (typ == tri_hp::text) {
//					fin >> omega >> vel >> ctr;
//				}
//			}
//			void output(const std::string& filename, tri_hp::filetype typ,int tlvl = 0) {
//				generic::output(filename,typ,tlvl);
//				if (typ == tri_hp::text) {
//					*x.gbl->log << omega << ' ' << vel << ' ' << ctr << '\n';
//				}
//				if (typ == tri_hp::tecplot) {
//					/* Calculate total flux, (this is getting weird) */
//					diff_flux = 0.0;
//					for(int j=0;j<base.nseg;++j) {
//						int sind = base.seg(j);
//						
//						x.ugtouht1d(sind);
//						element_rsdl(j,x.lf);
//						
//						for(int n=0;n<x.NV;++n)
//							diff_flux(n) += x.lf(n)(0);
//						
//						for(int n=0;n<x.NV;++n)
//							diff_flux(n) += x.lf(n)(1);
//					}
//				}
//			}
	};


	class characteristic : public generic {
		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx);
		public:
			characteristic(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
    };

    class symmetry : public generic {
		int dir;

		public:
			symmetry(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "symmetry";}
			symmetry(const symmetry& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin), dir(inbdry.dir) {}
			symmetry* create(tri_hp& xin, edge_bdry &bin) const {return new symmetry(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in) {
				generic::init(inmap,gbl_in);
				std::string keyword = base.idprefix +"_dir";
				inmap.getwdefault(keyword,dir,0);
				essential_indices.push_back(dir);
				type[dir] = essential;
			}
			void tadvance();
	};

	class applied_stress : public generic {
		Array<symbolic_function<2>,1> stress;

		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {

				/* CONTINUITY */
				flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
				/* X&Y MOMENTUM */
#ifdef INERTIALESS
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = -stress(n).Eval(xpt,x.gbl->time) +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#else
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = flx(x.NV-1)*u(n) -stress(n).Eval(xpt,x.gbl->time) +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#endif

				/* EVERYTHING ELSE */
				for (int n=tri_mesh::ND;n<x.NV-1;++n)
					flx(n) = flx(x.NV-1)*u(n);

				return;
			}
		public:
			applied_stress(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "applied_stress";}
			applied_stress(const applied_stress& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin), stress(inbdry.stress) {}
			applied_stress* create(tri_hp& xin, edge_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in);
	};
	
	class actuator_disc : public generic {
		bool start_pt_open, end_pt_open;
		symbolic_function<3> dp;
		
		void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
			/* CONTINUITY */
			flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
			
			if (base.is_frst()) {
				FLT norm_vel = flx(x.NV-1)/x.gbl->rho;
				TinyVector<FLT,3> inpt(xpt(0),xpt(1),norm_vel);
				FLT delta_p = dp.Eval(inpt,x.gbl->time);			
				
				/* X&Y MOMENTUM */
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = norm(n)*delta_p; 
			}
			else {
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = 0.0; 
			}
				
			/* EVERYTHING ELSE */
			for (int n=tri_mesh::ND;n<x.NV-1;++n)
				flx(n) = flx(x.NV-1)*u(n);

			return;
		}
		public:
			actuator_disc(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "actuator_disc";}
			actuator_disc(const actuator_disc& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin), dp(inbdry.dp) {}
			actuator_disc* create(tri_hp& xin, edge_bdry &bin) const {return new actuator_disc(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			
			/* To read in data */
			void init(input_map& inmap,void* gbl_in) {
				generic::init(inmap,gbl_in);
				
				/* LOAD PRESSURE JUMP FUNCTION */
				dp.init(inmap,base.idprefix+"_jump");
				
				/* Default is for axisymmetric case with one point on centerline */
				inmap.getwdefault(base.idprefix + "_start_open",start_pt_open,false);
				inmap.getwdefault(base.idprefix + "_end_open",start_pt_open,true);
				
				if (base.is_frst()) {
					bool temp = start_pt_open;
					start_pt_open = end_pt_open;
					end_pt_open = temp;
				}
			}
			void output(const std::string& fname, tri_hp::filetype typ,int tlvl = 0);

#ifdef petsc
			void petsc_matchjacobian_snd();
			int petsc_matchjacobian_rcv(int phase);
			void non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
			int non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
#endif
			void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride) {base.vloadbuff(boundary::all,pdata,0,x.NV-2,vrtstride*x.NV);}
			void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,0,x.NV-2, x.NV*vrtstride);}
			void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride); 
			int smatchsolution_rcv(FLT *sdata, int bgnmode, int endmode, int modestride);
	};
}
#endif
