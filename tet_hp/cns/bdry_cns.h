/*
 *  cns_bdry.h
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


#include "tet_hp_cns.h"
#include "../hp_boundary.h"
#include <myblas.h>
//#include <blitz/tinyvec-et.h>
#include <blitz/array.h>
#include <symbolic_function.h>

using namespace blitz;


//#define DETAILED_DT
//#define DETAILED_MINV

//#define L2_ERROR


namespace bdry_cns {

	class generic : public hp_face_bdry {
	protected:
		tet_hp_cns &x;
		
	public:
		generic(tet_hp_cns &xin, face_bdry &bin) : hp_face_bdry(xin,bin), x(xin) {mytype = "generic";}
		generic(const generic& inbdry, tet_hp_cns &xin, face_bdry &bin) : hp_face_bdry(inbdry,xin,bin), x(xin) {}
		generic* create(tet_hp& xin, face_bdry &bin) const {return new generic(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		void init(input_map& input,void* gbl_in) {
			hp_face_bdry::init(input,gbl_in);
		}
	};
		

	class neumann : public generic {
		protected:
			virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
		
		public:
			neumann(tet_hp_cns &xin, face_bdry &bin) : generic(xin,bin) {mytype = "neumann";}
			neumann(const neumann& inbdry, tet_hp_cns &xin, face_bdry &bin) : generic(inbdry,xin,bin) {}
			neumann* create(tet_hp& xin, face_bdry &bin) const {return new neumann(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
			void rsdl(int stage);
			void element_rsdl(int find,int stage);
			FLT beta_squared;
	};



	class inflow : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx);
		
	public:
		inflow(tet_hp_cns &xin, face_bdry &bin) : neumann(xin,bin) {
			mytype = "inflow";
			ndirichlets = x.NV-1;
			dirichlets.resize(ndirichlets);
			for (int n=1;n<x.NV;++n)
				dirichlets(n-1) = n;
		}
		inflow(const inflow& inbdry, tet_hp_cns &xin, face_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		inflow* create(tet_hp& xin, face_bdry &bin) const {return new inflow(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
		void vdirichlet();
		void edirichlet();	
		void fdirichlet();	
		
		void tadvance() {
			hp_face_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		};
		//void apply_sparse_dirichlet(bool compressed_column);
		//void modify_boundary_residual();
	};
	
	
	class adiabatic : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx);
		
	public:
		adiabatic(tet_hp_cns &xin, face_bdry &bin) : neumann(xin,bin) {
			mytype = "adiabatic";
			ndirichlets = x.NV-2;
			dirichlets.resize(ndirichlets);
			for (int n=1;n<x.NV-1;++n)
				dirichlets(n-1) = n;
		}
		adiabatic(const adiabatic& inbdry, tet_hp_cns &xin, face_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		adiabatic* create(tet_hp& xin, face_bdry &bin) const {return new adiabatic(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
		void vdirichlet();
		void edirichlet();	
		void fdirichlet();
		
		void tadvance() {
			hp_face_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		};
		//void apply_sparse_dirichlet(bool compressed_column);
		//void modify_boundary_residual();
	};
	
	class characteristic : public neumann {
	protected:
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
	public:
		characteristic(tet_hp_cns &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
		characteristic(const characteristic& inbdry, tet_hp_cns &xin, face_bdry &bin) : neumann(inbdry,xin,bin) {}
		characteristic* create(tet_hp& xin, face_bdry &bin) const {return new characteristic(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
	};
	
	class applied_stress : public neumann {
		Array<symbolic_function<3>,1> stress;

		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
		public:
			applied_stress(tet_hp_cns &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "applied_stress";}
			applied_stress(const applied_stress& inbdry, tet_hp_cns &xin, face_bdry &bin) : neumann(inbdry,xin,bin), stress(inbdry.stress) {}
			applied_stress* create(tet_hp& xin, face_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in);
	};
	
	class pure_slip : public neumann {
	protected:
//		Array<int,1> dirichlets;
//		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
	public:
		pure_slip(tet_hp_cns &xin, face_bdry &bin) : neumann(xin,bin) {
			mytype = "pure_slip";
//			ndirichlets = 1; // set to 1 if you want dirichlet
//			dirichlets.resize(ndirichlets);
//			dirichlets(0) = 1; // set x-velocity to zero		
		}
		pure_slip(const pure_slip& inbdry, tet_hp_cns &xin, face_bdry &bin) : neumann(inbdry,xin,bin) {}
		pure_slip* create(tet_hp& xin, face_bdry &bin) const {return new pure_slip(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
//		void vdirichlet();
//		void edirichlet();	
//		void fdirichlet();
//		
//		void tadvance() {
//			hp_face_bdry::tadvance();
//			setvalues(ibc,dirichlets,ndirichlets);
//		};
	};
	
	class generic_edge : public hp_edge_bdry {
	protected:
		tet_hp_cns &x;
		
	public:
		generic_edge(tet_hp_cns &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "generic_edge";}
		generic_edge(const generic_edge& inbdry, tet_hp_cns &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin) {}
		generic_edge* create(tet_hp& xin, edge_bdry &bin) const {return new generic_edge(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		void init(input_map& input,void* gbl_in) {
			hp_edge_bdry::init(input,gbl_in);
		}
	};
	
	class neumann_edge : public generic_edge {
	protected:
		virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
		
	public:
		neumann_edge(tet_hp_cns &xin, edge_bdry &bin) : generic_edge(xin,bin) {mytype = "neumann_edge";}
		neumann_edge(const neumann_edge& inbdry, tet_hp_cns &xin, edge_bdry &bin) : generic_edge(inbdry,xin,bin) {}
		neumann_edge* create(tet_hp& xin, edge_bdry &bin) const {return new neumann_edge(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		void rsdl(int stage);
		void element_rsdl(int eind,int stage);
	};
	
	
	
	class inflow_edge : public neumann_edge {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx);
		
	public:
		inflow_edge(tet_hp_cns &xin, edge_bdry &bin) : neumann_edge(xin,bin) {
			mytype = "inflow_edge";
			ndirichlets = x.NV-1;
			dirichlets.resize(ndirichlets);
			for (int n=1;n<x.NV;++n)
				dirichlets(n-1) = n;
		}
		inflow_edge(const inflow_edge& inbdry, tet_hp_cns &xin, edge_bdry &bin) : neumann_edge(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		inflow_edge* create(tet_hp& xin, edge_bdry &bin) const {return new inflow_edge(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
		void vdirichlet3d();
		void edirichlet3d();	
		
		void tadvance() {
			hp_edge_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		};
	};
	
	
	class adiabatic_edge : public neumann_edge {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx);
		
	public:
		adiabatic_edge(tet_hp_cns &xin, edge_bdry &bin) : neumann_edge(xin,bin) {
			mytype = "adiabatic_edge";
			ndirichlets = x.NV-2;
			dirichlets.resize(ndirichlets);
			for (int n=1;n<x.NV-1;++n)
				dirichlets(n-1) = n;
		}
		adiabatic_edge(const adiabatic_edge& inbdry, tet_hp_cns &xin, edge_bdry &bin) : neumann_edge(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		adiabatic_edge* create(tet_hp& xin, edge_bdry &bin) const {return new adiabatic_edge(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
		void vdirichlet3d();
		void edirichlet3d();	
		
		void tadvance() {
			hp_edge_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		};
	};
	
	class inflow_pt : public hp_vrtx_bdry {
	protected:
		tet_hp_cns &x;
		Array<int,1> dirichlets;
		int ndirichlets;
		
	public:
		inflow_pt(tet_hp_cns &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {
			mytype = "inflow_pt";
			ndirichlets = x.NV-1;
			dirichlets.resize(ndirichlets);
			for (int n=1;n<x.NV;++n)
				dirichlets(n-1) = n;
		}
		inflow_pt(const inflow_pt& inbdry, tet_hp_cns &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		inflow_pt* create(tet_hp& xin, vrtx_bdry &bin) const {return new inflow_pt(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
		void tadvance() {
			hp_vrtx_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		};

		void vdirichlet3d() {
			for(int n=0; n<ndirichlets; ++n) {
				x.gbl->res.v(base.pnt,dirichlets(n)) = 0.0;
			}
		}
		
	};
	
	class adiabatic_pt : public hp_vrtx_bdry {
	protected:
		tet_hp_cns &x;
		
	public:
		adiabatic_pt(tet_hp_cns &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "adiabatic_pt";}
		adiabatic_pt(const adiabatic_pt& inbdry, tet_hp_cns &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
		adiabatic_pt* create(tet_hp& xin, vrtx_bdry &bin) const {return new adiabatic_pt(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
		void tadvance() { 
			for(int n=1;n<x.NV-1;++n)
				x.ug.v(base.pnt,n) = x.gbl->ibc->f(n,x.pnts(base.pnt),x.gbl->time);  
			return;
		}
		
		void vdirichlet3d() {
			x.gbl->res.v(base.pnt,Range(1,x.NV-2)) = 0.0;
		}
		
	};
	
	class outflow_pt : public hp_vrtx_bdry {
	protected:
		tet_hp_cns &x;
		
	public:
		outflow_pt(tet_hp_cns &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "outflow_pt";}
		outflow_pt(const outflow_pt& inbdry, tet_hp_cns &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
		outflow_pt* create(tet_hp& xin, vrtx_bdry &bin) const {return new outflow_pt(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
		void tadvance() { 
			x.ug.v(base.pnt,0) = x.gbl->ibc->f(0,x.pnts(base.pnt),x.gbl->time);  
			return;
		}
		
		void vdirichlet3d() {
			x.gbl->res.v(base.pnt,0) = 0.0;
		}
		
	};
	
}
#endif
