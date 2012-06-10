/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include "../hp_boundary.h"
#include "myblas.h"
#include <symbolic_function.h>


namespace bdry_cd {
	class generic : public hp_edge_bdry {
		protected:
			tri_hp_cd &x;
			FLT diff_total,conv_total,circumference;
			
		public:
			generic(tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "generic";}
			generic(const generic& inbdry, tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin) {}
			generic* create(tri_hp& xin, edge_bdry &bin) const {return new generic(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);
	};

	class dirichlet : public generic {
		public:
			dirichlet(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin) {
				mytype = "dirichlet";
				type(0) = essential;
				essential_indices.push_back(0);
			}
			dirichlet(const dirichlet& inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			dirichlet* create(tri_hp& xin, edge_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
	};

	class characteristic : public generic {
		public:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
				FLT vel;

#ifdef CONST_A
				vel =  (x.gbl->ax-mv(0))*norm(0) +(x.gbl->ay -mv(1))*norm(1);        
#else
				vel =  (x.gbl->a->f(0,pt,x.gbl->time)-mv(0))*norm(0) +(x.gbl->a->f(1,pt,x.gbl->time) -mv(1))*norm(1);
#endif

				if (vel > 0.0)
					flx(0) = vel*u(0);
				else
					flx(0) = ibc->f(0, xpt, x.gbl->time)*vel;
			}
			characteristic(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic &inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
	};

	class mixed : public generic {
		public:
			vector_function flux_eq;

			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
				FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
				Array<FLT,1> axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND), au(1);
				axpt(0) = xpt(0); axpt(1) = xpt(1);
				amv(0) = mv(0); amv(1) = mv(1);
				anorm(0)= norm(0)/length; anorm(1) = norm(1)/length;
				flx(0) = flux_eq.Eval(u,axpt,amv,anorm,x.gbl->time)*length;
			}
		
			mixed(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "mixed";
				Array<string,1> names(4);
				Array<int,1> dims(4);
				dims = x.ND;
				names(0) = "u";
				dims(0) = x.NV;
				names(1) = "x";
				names(2) = "xt";
				names(3) = "n";
				flux_eq.set_arguments(4,dims,names);
			}
			mixed(const mixed& inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin), flux_eq(inbdry.flux_eq) {}
			mixed* create(tri_hp& xin, edge_bdry &bin) const {return new mixed(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}

			void init(input_map& inmap, void* gbl_in) {
				std::string keyword;
				std::istringstream data;
				std::string val;

				generic::init(inmap,gbl_in);
				
				if (inmap.find(base.idprefix +"_flux") != inmap.end()) {
					flux_eq.init(inmap,base.idprefix +"_flux");
				}
				else {
					*x.gbl->log << "couldn't find flux function " << base.idprefix +"_flux" << std::endl;
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
				return;
			}
	};
	
	
	class melt : public dirichlet {
		public:
			melt(tri_hp_cd &xin, edge_bdry &bin) : dirichlet(xin,bin) {mytype = "melt";}
			melt(const melt& inbdry, tri_hp_cd &xin, edge_bdry &bin) : dirichlet(inbdry,xin,bin) {}
			melt* create(tri_hp& xin, edge_bdry &bin) const {return new melt(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
						
			/* FOR COUPLED DYNAMIC BOUNDARIES */
			void init(input_map& inmap,void* gbl_in);
			void tadvance();
			void rsdl(int stage);
			void update(int stage);
#ifdef petsc
			void petsc_matchjacobian_snd();
			void petsc_matchjacobian_rcv(int phase);
			int petsc_make_1D_rsdl_vector(Array<FLT,1> res);
			void petsc_jacobian();
			void non_sparse(Array<int,1> &nnzero);
			void non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
			void non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
#endif
	};
	
	class kellerman : public melt {
		public:
			kellerman(tri_hp_cd &xin, edge_bdry &bin) : melt(xin,bin) {mytype = "kellerman";}
			kellerman(const kellerman& inbdry, tri_hp_cd &xin, edge_bdry &bin) : melt(inbdry,xin,bin) {}
			kellerman* create(tri_hp& xin, edge_bdry &bin) const {return new kellerman(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			
			/* FOR COUPLED DYNAMIC BOUNDARIES */
			void init(input_map& inmap,void* gbl_in);
#ifdef petsc
			void vdirichlet() {
				/* temperature equality constraint */
				essential_indices.push_back(0);
				dirichlet::vdirichlet();
				essential_indices.clear();
			}
			
			void sdirichlet(int m) {
				/* temperature equality constraint */
				essential_indices.push_back(0);
				dirichlet::sdirichlet(m);
				essential_indices.clear();
			}
#endif
	};
		
	
	
	class melt_end_pt : public hp_vrtx_bdry {
		/* INTERSECTING BOUNDARY CONTAINING END POINT MUST HAVE GEOMETRY NOT BE DEFINED SOLELY BY MESH */
		protected:
			tri_hp_cd& x;
			FLT res;
			FLT pnt0;
			FLT dt,dx;
			FLT diag_addition;
			
		public:
			melt_end_pt(tri_hp_cd &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "melt_end_pt";}
			melt_end_pt(const melt_end_pt& inbdry, tri_hp_cd &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
			melt_end_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new melt_end_pt(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			
			void init(input_map& inmap, void* gbl_in) {
				hp_vrtx_bdry::init(inmap,gbl_in);
				inmap.getwdefault(base.idprefix +"_diagonal_addition",diag_addition,0.0);
			}
			void rsdl(int stage);
			void element_rsdl();
			void setup_preconditioner();			
			void update(int stage);
#ifdef petsc
			void petsc_jacobian();
#endif
	};
		
		
}

