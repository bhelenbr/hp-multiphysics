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

#ifndef _bdry_cns_explicit_h_
#define _bdry_cns_explicit_h_


#include "tet_hp_cns_explicit.h"
#include "../hp_boundary.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array.h>
#include <symbolic_function.h>

using namespace blitz;


//#define DETAILED_DT
//#define DETAILED_MINV

//#define L2_ERROR


namespace bdry_cns_explicit {

	class generic : public hp_face_bdry {
	protected:
		tet_hp_cns_explicit &x;
		bool report_flag;
#ifdef L2_ERROR
		symbolic_function<2> l2norm;
#endif
		enum bctypes {ess, nat, mix};
		
	public:
		Array<FLT,1> total_flux,diff_flux,conv_flux;
		FLT circumference,moment,convect,circulation;
		
	public:
		generic(tet_hp_cns_explicit &xin, face_bdry &bin) : hp_face_bdry(xin,bin), x(xin) {mytype = "generic";}
		generic(const generic& inbdry, tet_hp_cns_explicit &xin, face_bdry &bin) : hp_face_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {
			if (report_flag) {
#ifdef L2_ERROR
				l2norm = inbdry.l2norm;
#endif
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);  
			}
		}
		generic* create(tet_hp& xin, face_bdry &bin) const {return new generic(*this,dynamic_cast<tet_hp_cns_explicit&>(xin),bin);}
		void init(input_map& input,void* gbl_in) {
			hp_face_bdry::init(input,gbl_in);
			std::string keyword = base.idprefix +"_report";
			input.getwdefault(keyword,report_flag,false);
			
			if (report_flag) {
#ifdef L2_ERROR
				l2norm.init(input,base.idprefix+"_norm");
#endif
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);            
			}
		}
		void output(std::ostream& fout, tet_hp::filetype typ,int tlvl = 0);
	};
		

	class neumann : public generic {
		protected:
			virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
		
		public:
			neumann(tet_hp_cns_explicit &xin, face_bdry &bin) : generic(xin,bin) {mytype = "neumann";}
			neumann(const neumann& inbdry, tet_hp_cns_explicit &xin, face_bdry &bin) : generic(inbdry,xin,bin) {}
			neumann* create(tet_hp& xin, face_bdry &bin) const {return new neumann(*this,dynamic_cast<tet_hp_cns_explicit&>(xin),bin);}
			void rsdl(int stage);
			void element_rsdl(int find,int stage);
	};



	class inflow : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx);
		
	public:
		inflow(tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(xin,bin) {
			mytype = "inflow";
			ndirichlets = x.NV;// not correct just for testing bug fix me temp
			dirichlets.resize(ndirichlets);
			for (int n=0;n<x.NV;++n)
				dirichlets(n) = n;
		}
		inflow(const inflow& inbdry, tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		inflow* create(tet_hp& xin, face_bdry &bin) const {return new inflow(*this,dynamic_cast<tet_hp_cns_explicit&>(xin),bin);}
		
		void vdirichlet();
		void edirichlet();	
		void fdirichlet();	
		
		void tadvance() {
			hp_face_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		};
		
	};
	
	
	class adiabatic : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx);
		
	public:
		adiabatic(tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(xin,bin) {
			mytype = "adiabatic";
			ndirichlets = x.NV-2;
			dirichlets.resize(ndirichlets);
			for (int n=1;n<x.NV-1;++n)
				dirichlets(n-1) = n;
		}
		adiabatic(const adiabatic& inbdry, tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		adiabatic* create(tet_hp& xin, face_bdry &bin) const {return new adiabatic(*this,dynamic_cast<tet_hp_cns_explicit&>(xin),bin);}
		
		void vdirichlet();
		void edirichlet();	
		void fdirichlet();
		
		void tadvance() {
			hp_face_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		};
		
	};
	
	class characteristic : public neumann {
	protected:
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
	public:
		characteristic(tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
		characteristic(const characteristic& inbdry, tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(inbdry,xin,bin) {}
		characteristic* create(tet_hp& xin, face_bdry &bin) const {return new characteristic(*this,dynamic_cast<tet_hp_cns_explicit&>(xin),bin);}
	};
	
	class applied_stress : public neumann {
		Array<symbolic_function<3>,1> stress;

		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
		public:
			applied_stress(tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "applied_stress";}
			applied_stress(const applied_stress& inbdry, tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(inbdry,xin,bin), stress(inbdry.stress) {}
			applied_stress* create(tet_hp& xin, face_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tet_hp_cns_explicit&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in);
	};
	
	class specified_flux : public neumann {
		Array<symbolic_function<3>,1> stress;
		
	protected:
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
	public:
		specified_flux(tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "specified_flux";}
		specified_flux(const specified_flux& inbdry, tet_hp_cns_explicit &xin, face_bdry &bin) : neumann(inbdry,xin,bin), stress(inbdry.stress) {}
		specified_flux* create(tet_hp& xin, face_bdry &bin) const {return new specified_flux(*this,dynamic_cast<tet_hp_cns_explicit&>(xin),bin);}
		void init(input_map& inmap,void* gbl_in);
	};
	
}
#endif
