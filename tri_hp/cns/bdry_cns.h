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

#ifndef _bdry_cns_h_
#define _bdry_cns_h_


#include "tri_hp_cns.h"
#include "../hp_boundary.h"
#include <myblas.h>
#include <blitz/array.h>
#include <symbolic_function.h>


using namespace blitz;

namespace bdry_cns {

	class generic : public hp_edge_bdry {
		protected:
			tri_hp_cns &x;
			Array<FLT,1> total_flux,diff_flux,conv_flux;
			FLT circumference,moment,convect,circulation;
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx);
		
		public:
			generic(tri_hp_cns &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "generic";}
			generic(const generic& inbdry, tri_hp_cns &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin) {
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);  
			}
			generic* create(tri_hp& xin, edge_bdry &bin) const {return new generic(*this,dynamic_cast<tri_hp_cns&>(xin),bin);}
			void init(input_map& input,void* gbl_in) {
				hp_edge_bdry::init(input,gbl_in);
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);            
			}
			void output(const std::string& fname, tri_hp::filetype typ,int tlvl = 0);
	};

	class inflow : public generic {  
		public:
			inflow(tri_hp_cns &xin, edge_bdry &bin) : generic(xin,bin) {
				mytype = "inflow";
				for (int n=1;n<x.NV;++n)
					essential_indices.push_back(n);
			}
			inflow(const inflow& inbdry, tri_hp_cns &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			inflow* create(tri_hp& xin, edge_bdry &bin) const {return new inflow(*this,dynamic_cast<tri_hp_cns&>(xin),bin);}
			void modify_boundary_residual(); 		
	};
	
	
	class adiabatic : public generic {  
		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx);
		public:
			adiabatic(tri_hp_cns &xin, edge_bdry &bin) : generic(xin,bin) {
				mytype = "adiabatic";
				for (int n=1;n<x.NV-1;++n)
					essential_indices.push_back(n);
			}
			adiabatic(const adiabatic& inbdry, tri_hp_cns &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			adiabatic* create(tri_hp& xin, edge_bdry &bin) const {return new adiabatic(*this,dynamic_cast<tri_hp_cns&>(xin),bin);}
			void modify_boundary_residual(); 		
	};

	class characteristic : public generic {
		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx);
		public:
			characteristic(tri_hp_cns &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic& inbdry, tri_hp_cns &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_cns&>(xin),bin);}
	};
	
	class applied_stress : public generic {
		Array<symbolic_function<2>,1> stress;

		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx); 				

		public:
			applied_stress(tri_hp_cns &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "applied_stress";}
			applied_stress(const applied_stress& inbdry, tri_hp_cns &xin, edge_bdry &bin) : generic(inbdry,xin,bin), stress(inbdry.stress) {}
			applied_stress* create(tri_hp& xin, edge_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tri_hp_cns&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in);
	};
	
	class symmetry : public generic {
		int dir;
		
		public:
			symmetry(tri_hp_cns &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "symmetry";}
			symmetry(const symmetry& inbdry, tri_hp_cns &xin, edge_bdry &bin) : generic(inbdry,xin,bin), dir(inbdry.dir) {}
			symmetry* create(tri_hp& xin, edge_bdry &bin) const {return new symmetry(*this,dynamic_cast<tri_hp_cns&>(xin),bin);}
			void init(input_map& input,void* gbl_in) {
				generic::init(input,gbl_in);
				std::string keyword = base.idprefix +"_dir";
				input.getwdefault(keyword,dir,0);
				essential_indices.push_back(dir+1);
				type[dir] = essential;
			}
			void tadvance();
	};



	class inflow_pt : public hp_vrtx_bdry {
		protected:
			tri_hp_cns &x;
			
		public:
			inflow_pt(tri_hp_cns &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "inflow_pt";}
			inflow_pt(const inflow_pt& inbdry, tri_hp_cns &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
			inflow_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new inflow_pt(*this,dynamic_cast<tri_hp_cns&>(xin),bin);}
			
			void tadvance() { 
				for(int n=1;n<x.NV;++n)
					x.ug.v(base.pnt,n) = x.gbl->ibc->f(n,x.pnts(base.pnt),x.gbl->time);  
				return;
			}
			
			void vdirichlet() {
				x.gbl->res.v(base.pnt,Range(1,x.NV-1)) = 0.0;
			}
		};


}
#endif
