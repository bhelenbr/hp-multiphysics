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
		tri_hp_cd &x;
		bool report_flag;
		FLT diff_total,conv_total,circumference;
		
	public:
		generic(tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin), report_flag(false) {mytype = "generic";}
		generic(const generic& inbdry, tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {}
		generic* create(tri_hp& xin, edge_bdry &bin) const {return new generic(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
		void init(input_map& input,void* gbl_in) {
			hp_edge_bdry::init(input,gbl_in);
			std::string keyword = base.idprefix +"_report";
			input.getwdefault(keyword,report_flag,false);       
		}
		void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);
	};

	class dirichlet : public generic {
		tri_hp_cd &x;

		public:
			dirichlet(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin), x(xin) {mytype = "dirichlet";}
			dirichlet(const dirichlet& inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin), x(xin) {}
			dirichlet* create(tri_hp& xin, edge_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			void vdirichlet() {
				int sind=-2,v0;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					x.gbl->res.v(v0,0) = 0.0;
				}
				v0 = x.seg(sind).pnt(1);
				x.gbl->res.v(v0,0) = 0.0;
			}

			void sdirichlet(int mode) {
				int sind;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					x.gbl->res.s(sind,mode,0) = 0.0;
				}
			}
			
#ifdef petsc
			void petsc_jacobian_dirichlet() {
				hp_edge_bdry::petsc_jacobian_dirichlet();  // Apply deforming mesh stuff
				
				int sm=basis::tri(x.log2p)->sm();
				Array<int,1> indices((base.nseg+1)*(x.NV) +base.nseg*sm*(x.NV));
				
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
					for(int n=0;n<x.NV;++n) {						
						indices(counter++)=gind+n;
					}
				} while (++j < base.nseg);
				v0 = x.seg(sind).pnt(1);
				gind = v0*vdofs;
				for(int n=0;n<x.NV;++n) {
					indices(counter++)=gind+n;
				}
				
				for(int i=0;i<base.nseg;++i) {
					gind = x.npnt*vdofs+base.seg(i)*sm*x.NV;
					for(int m=0; m<sm; ++m) {
						for(int n=0;n<x.NV;++n) {
							indices(counter++)=gind+m*x.NV+n;
						}
					}
				}	
			
#ifdef MY_SPARSE
				x.J.zero_rows(counter,indices);
				x.J_mpi.zero_rows(counter,indices);
				x.J.set_diag(counter,indices,1.0);
#else
				MatZeroRows(x.petsc_J,counter,indices.data(),1.0);
#endif
			}
#endif

			void tadvance(); 
	};

	class neumann : public generic {
		protected:
			tri_hp_cd &x;
			virtual FLT flux(FLT u, TinyVector<FLT,tri_mesh::ND> x, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm) {return(0.0);}

		public:
			neumann(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin), x(xin) {mytype = "neumann";}
			neumann(const neumann& inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin), x(xin) {}
			neumann* create(tri_hp& xin, edge_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			void element_rsdl(int eind,int stage);
	};


	class characteristic : public neumann {
		public:
			FLT flux(FLT u, TinyVector<FLT,tri_mesh::ND> pt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm) {
				FLT vel;

#ifdef CONST_A
				vel =  (x.gbl->ax-mv(0))*norm(0) +(x.gbl->ay -mv(1))*norm(1);        
#else
				vel =  (x.gbl->a->f(0,pt,x.gbl->time)-mv(0))*norm(0) +(x.gbl->a->f(1,pt,x.gbl->time) -mv(1))*norm(1);
#endif

				if (vel > 0.0)
					return(vel*u);

				return(ibc->f(0, pt, x.gbl->time)*vel);
			}
			characteristic(tri_hp_cd &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic &inbdry, tri_hp_cd &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
	};

	class mixed : public neumann {
		public:
			vector_function flux_eq;

			FLT flux(FLT u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm) {
				FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
				Array<FLT,1> axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND), au(1);
				axpt(0) = xpt(0); axpt(1) = xpt(1);
				amv(0) = mv(0); amv(1) = mv(1);
				anorm(0)= norm(0)/length; anorm(1) = norm(1)/length;
				au(0) = u;
				return(flux_eq.Eval(au,axpt,amv,anorm,x.gbl->time)*length);
			}
		
			mixed(tri_hp_cd &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "mixed";
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
			mixed(const mixed& inbdry, tri_hp_cd &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), flux_eq(inbdry.flux_eq) {}
			mixed* create(tri_hp& xin, edge_bdry &bin) const {return new mixed(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}

			void init(input_map& inmap, void* gbl_in) {
				std::string keyword;
				std::istringstream data;
				std::string val;

				neumann::init(inmap,gbl_in);
				
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

