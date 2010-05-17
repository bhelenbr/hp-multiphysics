/*
 *  ins_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _bdry_ins_h_
#define _bdry_ins_h_


#include "tri_hp_ins.h"
#include "../hp_boundary.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array.h>
#include <symbolic_function.h>


using namespace blitz;

//#define DETAILED_DT
//#define DETAILED_MINV

//#define L2_ERROR

namespace bdry_ins {

	class generic : public hp_edge_bdry {
		protected:
			tri_hp_ins &x;
			bool report_flag;
#ifdef L2_ERROR
			symbolic_function<2> l2norm;
#endif
			enum bctypes {ess, nat, mix};

		public:
			Array<FLT,1> total_flux,diff_flux,conv_flux;
			FLT circumference,moment,convect,circulation;

		public:
			generic(tri_hp_ins &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "generic";}
			generic(const generic& inbdry, tri_hp_ins &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {
				if (report_flag) {
#ifdef L2_ERROR
					l2norm = inbdry.l2norm;
#endif
					total_flux.resize(x.NV);
					diff_flux.resize(x.NV);
					conv_flux.resize(x.NV);  
				}
			}
			generic* create(tri_hp& xin, edge_bdry &bin) const {return new generic(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& input,void* gbl_in) {
				hp_edge_bdry::init(input,gbl_in);
				std::string keyword = base.idprefix +"_report";
				input.getwdefault(keyword,report_flag,false);

#ifdef L2_ERROR
				l2norm.init(input,base.idprefix+"_norm");
#endif
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);            
			}
			void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);
	};


	class neumann : public generic {
		protected:
			virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {

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
			neumann(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "neumann";}
			neumann(const neumann& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			neumann* create(tri_hp& xin, edge_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void element_rsdl(int eind,int stage);

	};



	class inflow : public neumann {  
		protected:
			Array<int,1> dirichlets;
			int ndirichlets;
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {

				/* CONTINUITY */
				flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));

				/* EVERYTHING ELSE DOESN'T MATTER */
				for (int n=0;n<x.NV-1;++n)
					flx(n) = 0.0;

				return;
			}

		public:
			inflow(tri_hp_ins &xin, edge_bdry &bin) : neumann(xin,bin) {
				mytype = "inflow";
				ndirichlets = x.NV-1;
				dirichlets.resize(x.NV-1);
				for (int n=0;n<x.NV-1;++n)
					dirichlets(n) = n;
			}
			inflow(const inflow& inbdry, tri_hp_ins &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
			inflow* create(tri_hp& xin, edge_bdry &bin) const {return new inflow(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}

			void vdirichlet() {
				int sind,j,v0;

				j = 0;
				do {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					x.gbl->res.v(v0,Range(0,x.NV-2)) = 0.0;
				} while (++j < base.nseg);
				v0 = x.seg(sind).pnt(1);
				x.gbl->res.v(v0,Range(0,x.NV-2)) = 0.0;
			}

			void sdirichlet(int mode) {
				int sind;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					x.gbl->res.s(sind,mode,Range(0,x.NV-2)) = 0.0;
				}
			}
		// temp fix me not sure how to make compatible
#ifdef petsc			
			void petsc_jacobian_dirichlet() {
				hp_edge_bdry::petsc_jacobian_dirichlet();  // Apply deforming mesh stuff
				
				int sm=basis::tri(x.log2p)->sm();
				Array<int,1> indices((base.nseg+1)*(x.NV-1) +base.nseg*sm*(x.NV-1));
				
				int vdofs;
				if (x.mmovement == x.coupled_deformable)
					vdofs = x.NV +tri_mesh::ND;
				else
					vdofs = x.NV;
					
				/* only works if pressure is 4th variable */
				int gind,v0,sind;
				int counter = 0;
				
				int j = 0;
				do {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					gind = v0*vdofs;
					for(int n=0;n<x.NV-1;++n) {						
						indices(counter++)=gind+n;
					}
				} while (++j < base.nseg);
				v0 = x.seg(sind).pnt(1);
				gind = v0*vdofs;
				for(int n=0;n<x.NV-1;++n) {
					indices(counter++)=gind+n;
				}

				for(int i=0;i<base.nseg;++i) {
					gind = x.npnt*vdofs+base.seg(i)*sm*x.NV;
					for(int m=0; m<sm; ++m) {
						for(int n=0;n<x.NV-1;++n) {
							indices(counter++)=gind+m*x.NV+n;
						}
					}
				}	
#ifdef MY_SPARSE
				x.my_zero_rows(counter,indices);
				x.my_set_diag(counter,indices,1.0);
#else
				MatZeroRows(x.petsc_J,counter,indices.data(),1.0);
#endif
			}
#endif

			void tadvance() {
				hp_edge_bdry::tadvance();
				setvalues(ibc,dirichlets,ndirichlets);
			}
	};

	class flexible : public inflow { 
		protected:
			Array<bctypes,1> type;
			init_bdry_cndtn *ibc;

			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {
				double flux = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));

				for (int n=0;n<x.NV-1;++n) {
					switch(type(n)) {
						case(ess): {
							flx(n) = 0.0;
						}
						case (nat): {
							FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
							flx(n) = ibc->f(n, xpt, x.gbl->time)*length;
						}
						case (mix): {
							FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
							for (int n=0;n<tri_mesh::ND;++n)
								flx(n) = flux*u(n) +ibc->f(n, xpt, x.gbl->time)*length;     
						}
					}
				}

				flx(x.NV-1) = flux;

				return;
			}

		public:
			flexible(tri_hp_ins &xin, edge_bdry &bin) : inflow(xin,bin) {mytype = "flexible"; type.resize(x.NV-1);}
			flexible(const flexible& inbdry, tri_hp_ins &xin, edge_bdry &bin) : inflow(inbdry,xin,bin) {type.resize(x.NV-1); type = inbdry.type; ibc = inbdry.ibc;}
			flexible* create(tri_hp& xin, edge_bdry &bin) const {return new flexible(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& input,void* gbl_in);
	};

	class euler : public neumann {
		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {
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
			euler(tri_hp_ins &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "euler";}
			euler(const euler& inbdry, tri_hp_ins &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
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
	};		

	class force_coupling : public inflow {
		rigid my_ibc;

		void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {
			for (int n=0;n<x.NV;++n)
				flx(n) = 0.0;
			return;
		}

		public:
			force_coupling(tri_hp_ins &xin, edge_bdry &bin) : inflow(xin,bin) {mytype = "force_coupling";}
			force_coupling(const force_coupling& inbdry, tri_hp_ins &xin, edge_bdry &bin) : inflow(inbdry,xin,bin) {}
			force_coupling* create(tri_hp& xin, edge_bdry &bin) const {return new force_coupling(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void tadvance() {hp_edge_bdry::tadvance();}
			void update(int stage) {
				if (!x.coarse_flag) setvalues(&my_ibc,dirichlets,ndirichlets);
			}
			void set_theta(FLT thetain) {};
			void set_omega(FLT dwdt) {my_ibc.omega = dwdt;}
			void set_ctr_rot(TinyVector<FLT,2> c) {my_ibc.ctr = c;}
			void set_vel(TinyVector<FLT,2> v) {my_ibc.vel = v;}
			void input(ifstream& fin,tri_hp::filetype typ,int tlvl) {
				inflow::input(fin,typ,tlvl);
				if (typ == tri_hp::text) {
					fin >> my_ibc.omega >> my_ibc.vel >> my_ibc.ctr;
				}
			}
			void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
				inflow::output(fout,typ,tlvl);
				if (typ == tri_hp::text) {
					fout << my_ibc.omega << ' ' << my_ibc.vel << ' ' << my_ibc.ctr << '\n';
				}
			}

	};
	
	class friction_slip : public neumann {
		rigid my_ibc;

		protected:
			FLT slip_length;
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {
				Array<FLT,1> urb(x.NV), ub(x.NV);

				/* RIGID COMPONENT */
				for(int n=0;n<x.NV;++n)
					urb(n) = my_ibc.f(n,xpt,x.gbl->time);
	
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

				FLT un = (ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1);				
				TinyVector<FLT,tri_mesh::ND> tang(-norm(1),norm(0));
				FLT slip_stress = ((ub(0) -u(0))*tang(1) +(ub(1) -u(1))*tang(0))*x.gbl->mu/(slip_length*sqrt(tang(0)*tang(0)+tang(1)*tang(1)));

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
			friction_slip(tri_hp_ins &xin, edge_bdry &bin) : neumann(xin,bin), slip_length(0.0) {mytype = "friction_slip";}
			friction_slip(const friction_slip& inbdry, tri_hp_ins &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), slip_length(inbdry.slip_length) {}
			friction_slip* create(tri_hp& xin, edge_bdry &bin) const {return new friction_slip(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& input,void* gbl_in) {
				neumann::init(input,gbl_in);
				std::string keyword = base.idprefix +"_sliplength";
				if (!input.get(keyword,slip_length)) {
					*x.gbl->log << "Couldn't find slip length for " << keyword << std::endl;
					exit(1);
				}
			}
			void set_omega(FLT dwdt) {my_ibc.omega = dwdt;}
			void set_ctr_rot(TinyVector<FLT,2> c) {my_ibc.ctr = c;}
			void set_vel(TinyVector<FLT,2> v) {my_ibc.vel = v;}
			void input(ifstream& fin,tri_hp::filetype typ,int tlvl) {
				neumann::input(fin,typ,tlvl);
				if (typ == tri_hp::text) {
					fin >> my_ibc.omega >> my_ibc.vel >> my_ibc.ctr;
				}
			}
			void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
				// neumann::output(fout,typ,tlvl);
				if (typ == tri_hp::text) {
					fout << my_ibc.omega << ' ' << my_ibc.vel << ' ' << my_ibc.ctr << '\n';
				}
				if (typ == tri_hp::tecplot) {
					/* Calculate total flux, (this is getting weird) */
					total_flux = 0.0;
					for(int j=0;j<base.nseg;++j) {
						int sind = base.seg(j);
						
						x.ugtouht1d(sind);
						element_rsdl(j,0);
						
						for(int n=0;n<x.NV;++n)
							total_flux(n) += x.lf(n)(0);
						
						for(int n=0;n<x.NV;++n)
							total_flux(n) += x.lf(n)(1);
					}
				}
			}
	};


	class characteristic : public neumann {
		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx);
		public:
			characteristic(tri_hp_ins &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic& inbdry, tri_hp_ins &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
    };

    class symmetry : public generic {
		int dir;

		public:
			symmetry(tri_hp_ins &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "symmetry";}
			symmetry(const symmetry& inbdry, tri_hp_ins &xin, edge_bdry &bin) : generic(inbdry,xin,bin), dir(inbdry.dir) {}
			symmetry* create(tri_hp& xin, edge_bdry &bin) const {return new symmetry(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& input,void* gbl_in) {
				generic::init(input,gbl_in);
				std::string keyword = base.idprefix +"_dir";
				input.getwdefault(keyword,dir,0);
			}

			void vdirichlet() {
				int sind,j,v0;
				j = 0;
				do {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					x.gbl->res.v(v0,dir) = 0.0;
				} while (++j < base.nseg);
				v0 = x.seg(sind).pnt(1);
				x.gbl->res.v(v0,dir) = 0.0;
			}

			void sdirichlet(int mode) {
				int sind;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					x.gbl->res.s(sind,mode,dir) = 0.0;
				}
			}
			
#ifdef petsc			
			void petsc_jacobian_dirichlet() {
				hp_edge_bdry::petsc_jacobian_dirichlet();  // Apply deforming mesh stuff
				
				int sm=basis::tri(x.log2p)->sm();
				Array<int,1> indices((base.nseg+1)*(x.NV-1) +base.nseg*sm*(x.NV-1));
				
				int vdofs;
				if (x.mmovement == x.coupled_deformable)
					vdofs = x.NV +tri_mesh::ND;
				else
					vdofs = x.NV;
					
				/* only works if pressure is 4th variable */
				int gind,v0,sind;
				int counter = 0;
				
				int j = 0;
				do {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					gind = v0*vdofs;
					indices(counter++)=gind+dir;
				} while (++j < base.nseg);
				v0 = x.seg(sind).pnt(1);
				gind = v0*vdofs;
				indices(counter++)=gind+dir;

				for(int i=0;i<base.nseg;++i) {
					gind = x.npnt*vdofs+base.seg(i)*sm*x.NV;
					for(int m=0; m<sm; ++m) {
							indices(counter++)=gind+m*x.NV+dir;
					}
				}	
#ifdef MY_SPARSE
				x.my_zero_rows(counter,indices);
				x.my_set_diag(counter,indices,1.0);
#else
				MatZeroRows(x.petsc_J,counter,indices.data(),1.0);
#endif
			}
#endif

			void tadvance();
	};

	class applied_stress : public neumann {
		Array<symbolic_function<2>,1> stress;

		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {

				/* CONTINUITY */
				flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));

				FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
				/* X&Y MOMENTUM */
#ifdef INERTIALESS
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = -stress(n).Eval(xpt,x.gbl->time)*length +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#else
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = flx(x.NV-1)*u(n) -stress(n).Eval(xpt,x.gbl->time)*length +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#endif

				/* EVERYTHING ELSE */
				for (int n=tri_mesh::ND;n<x.NV-1;++n)
					flx(n) = flx(x.NV-1)*u(n);

				return;
			}
		public:
			applied_stress(tri_hp_ins &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "applied_stress";}
			applied_stress(const applied_stress& inbdry, tri_hp_ins &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), stress(inbdry.stress) {}
			applied_stress* create(tri_hp& xin, edge_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in);
	};

	class surface_slave : public neumann {
		public:
			surface_slave(tri_hp_ins &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "surface_slave";}
			surface_slave(const surface_slave& inbdry, tri_hp_ins &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			surface_slave* create(tri_hp& xin, edge_bdry &bin) const {return new surface_slave(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in);

			/* FOR COUPLED DYNAMIC BOUNDARIES */
			void tadvance();
			void rsdl(int stage);
			void update(int stage);

			void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride) {base.vloadbuff(boundary::all,pdata,0,x.NV-2,vrtstride*x.NV);}
			void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,0,x.NV-2, x.NV*vrtstride);}
			void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride); 
			void smatchsolution_rcv(FLT *sdata, int bgnmode, int endmode, int modestride);
	};


	class surface : public surface_slave {
		protected:
			Array<FLT,1> ksprg;
			Array<TinyVector<FLT,tri_mesh::ND>,1> vug_frst;
			Array<TinyVector<FLT,tri_mesh::ND>,2> vdres; //!< Driving term for multigrid (log2p, pnts)
			Array<TinyVector<FLT,tri_mesh::ND>,3> sdres; //!< Driving term for multigrid (log2p, side, order)
			const surface *fine, *coarse;

		public:
			struct global {                
				bool is_loop;
				/* FLUID PROPERTIES */
				FLT sigma,rho2,mu2;

				/* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
				Array<TinyVector<FLT,tri_mesh::ND>,1> vug0;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sug0;

				/* RESIDUALS */
				Array<TinyVector<FLT,tri_mesh::ND>,1> vres;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sres;
				Array<TinyVector<FLT,tri_mesh::ND>,1> vres0;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sres0;
#ifdef DROP
				FLT penalty,vflux;
				Array<FLT,1> vvolumeflux;
				Array<FLT,2> svolumeflux;
#endif

				/* PRECONDITIONER */
				Array<TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND>,1> vdt;
				Array<TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND>,1> sdt;
				Array<FLT,1> meshc;

#ifdef DETAILED_MINV
				Array<TinyMatrix<FLT,2*MAXP,2*MAXP>,1> ms;
				Array<FLT,5> vms;
				Array<TinyVector<int,2*MAXP>,1> ipiv;
#endif
				TinyVector<FLT,tri_mesh::ND> fadd;
				TinyMatrix<FLT,tri_mesh::ND,MAXP> cfl;
				FLT adis;
			} *gbl; 

		public:
			void* create_global_structure() {return new global;}
			surface(tri_hp_ins &xin, edge_bdry &bin) : surface_slave(xin,bin) {mytype = "surface";}
			surface(const surface& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : surface_slave(inbdry,xin,bin) {
				gbl = inbdry.gbl;
				ksprg.resize(base.maxseg);
				vug_frst.resize(base.maxseg+1);
				vdres.resize(1,base.maxseg+1);
				fine = &inbdry;
			};
			surface* create(tri_hp& xin, edge_bdry &bin) const {return new surface(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			void init(input_map& input,void* gbl_in); 

			/* FOR COUPLED DYNAMIC BOUNDARIES */
			void tadvance();
			void rsdl(int stage);
			void element_rsdl(int sind, Array<FLT,2> lf);  // FIXME Not really compatible need to make all consistent
			void maxres();
			void setup_preconditioner();
			void minvrt();
			void update(int stage);
			void mg_restrict(); 
			void element_jacobian(int indx, Array<FLT,2>& K);
#ifdef petsc
			void petsc_jacobian();
			void petsc_jacobian_dirichlet() {} // Override r_mesh constraint on mesh positions
			int petsc_rsdl(Array<FLT,1> res);
#endif
	};
	
	class actuator_disc : public surface_slave {
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {

				/* CONTINUITY */
				flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));

				// FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
				/* X&Y MOMENTUM */
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n) = norm(n);  //FIXME
					
				/* EVERYTHING ELSE */
				for (int n=tri_mesh::ND;n<x.NV-1;++n)
					flx(n) = flx(x.NV-1)*u(n);

				return;
			}
		public:
			actuator_disc(tri_hp_ins &xin, edge_bdry &bin) : surface_slave(xin,bin) {mytype = "actuator_disc";}
			actuator_disc(const actuator_disc& inbdry, tri_hp_ins &xin, edge_bdry &bin) : surface_slave(inbdry,xin,bin) {}
			actuator_disc* create(tri_hp& xin, edge_bdry &bin) const {return new actuator_disc(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			
			/* To read in data */
			void init(input_map& inmap,void* gbl_in) {
				/* LOAD K */
				
				if (!base.is_frst()) {
					// SET K TO 0.0
				}
			}

			/* Reset to normal stuff */
			void tadvance() {neumann::tadvance();}
			void rsdl(int stage) {neumann::rsdl(stage);}
			void update(int stage) {neumann::update(stage);}
	};

	/*******************************/
	/* VERTEX BOUNDARY CONDITIONS */
	/******************************/
	class surface_fixed_pt : public hp_vrtx_bdry {
		/* INTERSECTING BOUNDARY CONTAINING END POINT MUST HAVE GEOMETRY NOT BE DEFINED SOLELY BY MESH */
		protected:
			tri_hp_ins &x;
			surface *surf;
			int surfbdry;
			bool fix_norm;

		public:
			surface_fixed_pt(tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin), fix_norm(1) {mytype = "surface_fixed_pt";}
			surface_fixed_pt(const surface_fixed_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin), 
				surfbdry(inbdry.surfbdry),fix_norm(inbdry.fix_norm) {
				if (!(surf = dynamic_cast<surface *>(x.hp_ebdry(base.ebdry(surfbdry))))) {
					*x.gbl->log << "something's wrong can't find surface boundary" << std::endl;
					exit(1);
				}
			}
			surface_fixed_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_fixed_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}

			void init(input_map& inmap,void* gbl_in) {
				std::string keyword,val;
				std::istringstream data;
				std::string filename;

				hp_vrtx_bdry::init(inmap,gbl_in);

				if ((surf = dynamic_cast<surface *>(x.hp_ebdry(base.ebdry(0))))) {
					surfbdry = 0;
				}
				else if ((surf = dynamic_cast<surface *>(x.hp_ebdry(base.ebdry(1))))) {
					surfbdry = 1;
				}
				else {
					*x.gbl->log << "something's wrong neither side is a surface boundary" << std::endl;
					exit(1);
				}
				fix_norm = 1;  // THIS IS ONLY TO BE CHANGE BY INHERITED CLASSES
			}

			void rsdl(int stage) {
				if (surfbdry == 0) {
					/* SET TANGENT RESIDUAL TO ZERO */
					surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(0) = 0.0;
					if (fix_norm) {
#ifndef petsc
						/* POST-REMOVE ADDED MASS FLUX TERM FOR FIXED POINT */
						x.gbl->res.v(base.pnt,x.NV-1) += surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(1)*x.gbl->rho;
#endif
						/* AND ZERO RESIDUAL */
						surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(1) = 0.0;
					}
				}
				else {
					/* SET TANGENT RESIDUAL TO ZERO */
					surf->gbl->vres(0)(0) = 0.0;
					if (fix_norm) {
#ifndef petsc
						/* POST-REMOVE ADDED MASS FLUX TERM FOR FIXED POINT */
						x.gbl->res.v(base.pnt,x.NV-1) += surf->gbl->vres(0)(1)*x.gbl->rho;
#endif
						/* AND ZERO RESIDUAL */
						surf->gbl->vres(0)(1) = 0.0;
					}
				}
				return;
			}
			
			void vdirichlet() {
				if (surfbdry == 0) {
					for(int n=0;n<=fix_norm;++n) 
						surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(n) = 0.0;
				}
				else {
					for(int n=0;n<=fix_norm;++n) 
						surf->gbl->vres(0)(n) = 0.0;
				}
			}

			void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
				if (surfbdry == 0) {
					x.ebdry(base.ebdry(1))->mvpttobdry(0,-1.0,pt);
				}
				else {
					x.ebdry(base.ebdry(0))->mvpttobdry(x.ebdry(base.ebdry(0))->nseg-1,1.0,pt);
				}
			}
			
#ifdef petsc
			void petsc_jacobian_dirichlet() {
				/* BOTH X & Y ARE FIXED */
				if (fix_norm) {
					Array<int,1> rows(tri_mesh::ND);
					for(int n=0;n<tri_mesh::ND;++n) 
						rows(n) = (x.NV+tri_mesh::ND)*base.pnt +x.NV +n;
#ifdef MY_SPARSE
					x.my_zero_rows(tri_mesh::ND,rows);
					x.my_set_diag(tri_mesh::ND,rows,1.0);
#else				
					MatZeroRows(x.petsc_J,tri_mesh::ND,rows.data(),1.0);
#endif
				}
			}
#endif
	};
	
	class surface_outflow : public surface_fixed_pt {
		/* For a surface point sliding on a vertical or horizontal surface */
		/* For periodic wall have tri_mesh vertex type be comm */
		protected:
			FLT position;
			enum {vertical, horizontal, angled} wall_type;
			enum {prdc, free_angle, fixed_angle} contact_type;
			FLT contact_angle;  // ONLY USED FOR FIXED_ANGLE CONTACT TYPE
			TinyVector<FLT,tri_mesh::ND> wall_normal; 

		public:
			surface_outflow(tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(xin,bin) {mytype = "surface_outflow";}
			surface_outflow(const surface_outflow& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(inbdry,xin,bin), position(inbdry.position), 
				wall_type(inbdry.wall_type), contact_type(inbdry.contact_type), contact_angle(inbdry.contact_angle), wall_normal(inbdry.wall_normal) {}
			surface_outflow* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_outflow(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}

			void init(input_map& inmap,void* gbl_in) {                
				surface_fixed_pt::init(inmap,gbl_in);
				fix_norm = 0;

				std::string input;
				if (!inmap.get(base.idprefix + "_contact_type",input))
					contact_type = prdc;
				else if (input == "free_angle")
					contact_type = free_angle;
				else if (input == "fixed_angle") {
					contact_type = fixed_angle;
					inmap.getwdefault(base.idprefix +"_contact_angle",contact_angle,90.0);
					contact_angle *= M_PI/180.0;
				}
				else {
					*x.gbl->log << "Unrecognized contact type" << std::endl;
				}
				
				if (!inmap.get(base.idprefix + "_wall_type",input)) {
					wall_type = vertical;
					position = x.pnts(base.pnt)(0);
				}
				else if (input == "horizontal") {
					wall_type = horizontal;
					position = x.pnts(base.pnt)(1);
				}
				else if (input == "angled")
					wall_type = angled;
				else {
					*x.gbl->log << "Unrecognized wall type" << std::endl;
				}
				
				/* Find tangent to wall and use to constrain motion */
				int bnumwall = base.ebdry(1-surfbdry);
				TinyVector<FLT,tri_mesh::ND> rp;
				if (surfbdry == 0) {
					int sindwall = x.ebdry(bnumwall)->seg(0);
					x.crdtocht1d(sindwall);
					basis::tri(x.log2p)->ptprobe1d(2,&rp(0),&wall_normal(0),-1.0,&x.cht(0,0),MXTM);
				}
				else {
					int sindwall = x.ebdry(bnumwall)->seg(x.ebdry(bnumwall)->nseg-1);
					x.crdtocht1d(sindwall);
					basis::tri(x.log2p)->ptprobe1d(2,&rp(0),&wall_normal(0),1.0,&x.cht(0,0),MXTM);
				}
				FLT length = sqrt(wall_normal(0)*wall_normal(0) +wall_normal(1)*wall_normal(1));
				wall_normal /= length;
				FLT temp = wall_normal(0);
				wall_normal(0) = wall_normal(1);
				wall_normal(1) = -temp;
			}

			void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
				switch(wall_type) {
					case(vertical):
						pt(0) = position;
						break;
					case(horizontal):
						pt(1) = position;
						break;
					case(angled):
						surface_fixed_pt::mvpttobdry(pt);
						break;
				}
			}
			
			/* Routine to add surface tension stress */
			/* also zero's tangent residual in no petsc */
			void rsdl(int stage);
						
#ifdef petsc
			void petsc_jacobian(); 
			void petsc_jacobian_dirichlet() {}
#endif
	};
		
	class inflow_pt : public hp_vrtx_bdry {
		protected:
			tri_hp_ins &x;

		public:
			inflow_pt(tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "inflow_pt";}
			inflow_pt(const inflow_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
			inflow_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new inflow_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}

			void tadvance() { 
				for(int n=0;n<x.NV-1;++n)
					x.ug.v(base.pnt,n) = x.gbl->ibc->f(n,x.pnts(base.pnt),x.gbl->time);  
				return;
			}

			void vdirichlet2d() {
				x.gbl->res.v(base.pnt,Range(0,x.NV-2)) = 0.0;
			}
	};

	class hybrid_slave_pt : public hp_vrtx_bdry {
		protected:
			tri_hp_ins &x;

		public:
			hybrid_slave_pt(tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "hybrid_slave_pt";}
			hybrid_slave_pt(const hybrid_slave_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
			hybrid_slave_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new hybrid_slave_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}

			void vdirichlet2d() {
				x.gbl->res.v(base.pnt,Range(0,x.NV-1)) = 0.0;
			}

			void update(int stage);
	};

	class hybrid_pt : public surface_outflow {
		public:
			hybrid_pt(tri_hp_ins &xin, vrtx_bdry &bin) : surface_outflow(xin,bin) {mytype = "hybrid_pt"; contact_type = prdc;}
			hybrid_pt(const hybrid_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_outflow(inbdry,xin,bin) {}
			hybrid_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new hybrid_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			
			void init(input_map& inmap,void* gbl_in) {                
				surface_outflow::init(inmap,gbl_in);
				contact_type = prdc;
			}

			void vdirichlet2d() {
				x.gbl->res.v(base.pnt,Range(0,x.NV-1)) = 0.0;
			}

			void rsdl(int stage); 			
			void update(int stage);
	};

}
#endif
