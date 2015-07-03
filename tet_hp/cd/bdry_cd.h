/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cd.h"
#include "../hp_boundary.h"
#include "myblas.h"

namespace bdry_cd {
	
	class generic : public hp_face_bdry {
		tet_hp_cd &x;
		bool report_flag;
		
	public:
		generic(tet_hp_cd &xin, face_bdry &bin) : hp_face_bdry(xin,bin), x(xin), report_flag(false) {mytype = "generic";}
		generic(const generic& inbdry, tet_hp_cd &xin, face_bdry &bin) : hp_face_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {}
		generic* create(tet_hp& xin, face_bdry &bin) const {return new generic(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void init(input_map& inmap,void* gbl_in) {
			hp_face_bdry::init(inmap,gbl_in);
			std::string keyword = base.idprefix +"_report";
			inmap.getwdefault(keyword,report_flag,false);       
		}
		//void output(std::ostream& fout, tet_hp::filetype typ,int tlvl = 0);
	};
	
	class dirichlet : public generic {
		tet_hp_cd &x;
		
		public:
		dirichlet(tet_hp_cd &xin, face_bdry &bin) : generic(xin,bin), x(xin) {mytype = "dirichlet";}
		dirichlet(const dirichlet& inbdry, tet_hp_cd &xin, face_bdry &bin) : generic(inbdry,xin,bin), x(xin) {}
		dirichlet* create(tet_hp& xin, face_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void vdirichlet() { 
			int v0; 											
			for(int j=0;j<base.npnt;++j) {
				v0 = base.pnt(j).gindx;
					x.gbl->res.v(v0,0) = 0.0;
			}			
		}

		void edirichlet() {
			int sind;
			if (basis::tet(x.log2p).em > 0) {
				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j).gindx;
					x.gbl->res.e(sind,Range(0,basis::tet(x.log2p).em-1),0) = 0.0;
				}
			}
		}
			
		void fdirichlet() {
			int find;

			if (basis::tet(x.log2p).fm > 0) {
				for(int j=0;j<base.ntri;++j) {
					find = base.tri(j).gindx;
					x.gbl->res.f(find,Range(0,basis::tet(x.log2p).fm-1),0) = 0.0;
				}
			}
		}
//			void apply_sparse_dirichlet(bool compressed_column) {
//				int gind;
//				int em=basis::tet(x.log2p).em;
//				int fm=basis::tet(x.log2p).fm;
//
//				for(int i=0;i<base.npnt;++i)
//					x.sparse_dirichlet(base.pnt(i).gindx,compressed_column);
//				
//				for(int i=0;i<base.nseg;++i){
//					gind = x.npnt+base.seg(i).gindx*em;
//					for(int m=0; m<em; ++m)
//						x.sparse_dirichlet(gind+m,compressed_column);
//				}
//				
//				for(int i=0;i<base.ntri;++i){
//					gind = x.npnt+x.nseg*em+base.tri(i).gindx*fm;
//					for(int m=0; m<fm; ++m)
//						x.sparse_dirichlet(gind+m,compressed_column);
//				}				
//			}

			void tadvance(); 
	};

	class neumann : public generic {
		protected:
			tet_hp_cd &x;
			virtual FLT flux(FLT u, TinyVector<FLT,tet_mesh::ND> x, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm) {return(0.0);}
		
		public:
		neumann(tet_hp_cd &xin, face_bdry &bin) : generic(xin,bin), x(xin) {mytype = "neumann";}
		neumann(const neumann& inbdry, tet_hp_cd &xin, face_bdry &bin) : generic(inbdry,xin,bin), x(xin) {}
		neumann* create(tet_hp& xin, face_bdry &bin) const {return new neumann(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void rsdl(int stage);
		void element_rsdl(int find,int stage);
	};

	
	class generic_edge : public hp_edge_bdry {
		tet_hp_cd &x;
		bool report_flag;
		
	public:
		generic_edge(tet_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin), report_flag(false) {mytype = "generic_edge";}
		generic_edge(const generic_edge& inbdry, tet_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {}
		generic_edge* create(tet_hp& xin, edge_bdry &bin) const {return new generic_edge(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void init(input_map& inmap,void* gbl_in) {
			hp_edge_bdry::init(inmap,gbl_in);
			std::string keyword = base.idprefix +"_report";
			inmap.getwdefault(keyword,report_flag,false);       
		}
		//void output(std::ostream& fout, tet_hp::filetype typ,int tlvl = 0);
	};
	
	class dirichlet_edge : public generic_edge {
		tet_hp_cd &x;
		
	public:
		dirichlet_edge(tet_hp_cd &xin, edge_bdry &bin) : generic_edge(xin,bin), x(xin) {mytype = "dirichlet_edge";}
		dirichlet_edge(const dirichlet_edge& inbdry, tet_hp_cd &xin, edge_bdry &bin) : generic_edge(inbdry,xin,bin), x(xin) {}
		dirichlet_edge* create(tet_hp& xin, edge_bdry &bin) const {return new dirichlet_edge(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void vdirichlet3d() {
			int sind=-2,v0;
			
			for(int j=0;j<base.nseg;++j) {
				sind = base.seg(j).gindx;
				v0 = x.seg(sind).pnt(0);
				x.gbl->res.v(v0,0) = 0.0;
			}
			v0 = x.seg(sind).pnt(1);
			x.gbl->res.v(v0,0) = 0.0;
		}
		
		void edirichlet3d() {			
			int sind;
			if (basis::tet(x.log2p).em > 0) {
				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j).gindx;
					x.gbl->res.e(sind,Range(0,basis::tet(x.log2p).em-1),0) = 0.0;
				}
			}
		}
		
		void tadvance(); 
	};
	
	class neumann_edge : public generic_edge {
	protected:
		tet_hp_cd &x;
		virtual FLT flux(FLT u, TinyVector<FLT,tet_mesh::ND> x, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm) {return(0.0);}
		
	public:
		neumann_edge(tet_hp_cd &xin, edge_bdry &bin) : generic_edge(xin,bin), x(xin) {mytype = "neumann_edge";}
		neumann_edge(const neumann_edge& inbdry, tet_hp_cd &xin, edge_bdry &bin) : generic_edge(inbdry,xin,bin), x(xin) {}
		neumann_edge* create(tet_hp& xin, edge_bdry &bin) const {return new neumann_edge(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void rsdl(int stage);
		void element_rsdl(int eind,int stage);
	};
	
	
	class generic_pt : public hp_vrtx_bdry {
		tet_hp_cd &x;
		bool report_flag;
		
	public:
		generic_pt(tet_hp_cd &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin), report_flag(false) {mytype = "generic_pt";}
		generic_pt(const generic_pt& inbdry, tet_hp_cd &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {}
		generic_pt* create(tet_hp& xin, vrtx_bdry &bin) const {return new generic_pt(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void init(input_map& inmap,void* gbl_in) {
			hp_vrtx_bdry::init(inmap,gbl_in);
			std::string keyword = base.idprefix +"_report";
			inmap.getwdefault(keyword,report_flag,false);       
		}
		//void output(std::ostream& fout, tet_hp::filetype typ,int tlvl = 0);
	};
	
	class dirichlet_pt : public generic_pt {
		tet_hp_cd &x;
		
	public:
		dirichlet_pt(tet_hp_cd &xin, vrtx_bdry &bin) : generic_pt(xin,bin), x(xin) {mytype = "dirichlet_pt";}
		dirichlet_pt(const dirichlet_pt& inbdry, tet_hp_cd &xin, vrtx_bdry &bin) : generic_pt(inbdry,xin,bin), x(xin) {}
		dirichlet_pt* create(tet_hp& xin, vrtx_bdry &bin) const {return new dirichlet_pt(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void vdirichlet3d() { 
			x.gbl->res.v(base.pnt,0) = 0.0;					
		}
		
		void tadvance(); 
	};
	
	class neumann_pt : public generic_pt {
	protected:
		tet_hp_cd &x;
		virtual FLT flux(FLT u, TinyVector<FLT,tet_mesh::ND> x, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm) {return(0.0);}
		
	public:
		neumann_pt(tet_hp_cd &xin, vrtx_bdry &bin) : generic_pt(xin,bin), x(xin) {mytype = "neumann_pt";}
		neumann_pt(const neumann_pt& inbdry, tet_hp_cd &xin, vrtx_bdry &bin) : generic_pt(inbdry,xin,bin), x(xin) {}
		neumann_pt* create(tet_hp& xin, vrtx_bdry &bin) const {return new neumann_pt(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
		void rsdl(int stage);
		void element_rsdl(int stage);
	};
	
	
	
}

