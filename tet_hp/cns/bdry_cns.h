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
#include <blitz/tinyvec-et.h>
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
		bool report_flag;
#ifdef L2_ERROR
		symbolic_function<2> l2norm;
#endif
		enum bctypes {ess, nat, mix};
		
	public:
		Array<FLT,1> total_flux,diff_flux,conv_flux;
		FLT circumference,moment,convect,circulation;
		
	public:
		generic(tet_hp_cns &xin, face_bdry &bin) : hp_face_bdry(xin,bin), x(xin) {mytype = "generic";}
		generic(const generic& inbdry, tet_hp_cns &xin, face_bdry &bin) : hp_face_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {
			if (report_flag) {
#ifdef L2_ERROR
				l2norm = inbdry.l2norm;
#endif
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);  
			}
		}
		generic* create(tet_hp& xin, face_bdry &bin) const {return new generic(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
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
			virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {

				/* CONTINUITY */
				flx(0) = ibc->f(0, xpt, x.gbl->time)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1)+(u(3) -mv(2))*norm(2));
				
				/* X&Y MOMENTUM */
#ifdef INERTIALESS
				for (int n=1;n<tet_mesh::ND+1;++n)
					flx(n) = ibc->f(0, xpt, x.gbl->time)*norm(n-1);
#else
				for (int n=1;n<tet_mesh::ND+1;++n)
					flx(n) = flx(0)*u(n) +ibc->f(0, xpt, x.gbl->time)*norm(n-1);
#endif
				
				/* ENERGY EQUATION */
				double h = x.gbl->gamma/(x.gbl->gamma-1.0)*u(x.NV-1) +0.5*(u(1)*u(1)+u(2)*u(2)+u(3)*u(3));				
				flx(x.NV-1) = h*flx(0);

				return;
			}
		
		public:
			neumann(tet_hp_cns &xin, face_bdry &bin) : generic(xin,bin) {mytype = "neumann";}
			neumann(const neumann& inbdry, tet_hp_cns &xin, face_bdry &bin) : generic(inbdry,xin,bin) {}
			neumann* create(tet_hp& xin, face_bdry &bin) const {return new neumann(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
			void rsdl(int stage);
			void element_rsdl(int find,int stage);
	};



	class inflow : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
			
			/* CONTINUITY */
			flx(0) = ibc->f(0, xpt, x.gbl->time)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1)+(u(3) -mv(2))*norm(2));
			
			/* EVERYTHING ELSE DOESN'T MATTER */
			for (int n=1;n<x.NV;++n)
				flx(n) = 0.0;
			
			return;
		}
		
	public:
		inflow(tet_hp_cns &xin, face_bdry &bin) : neumann(xin,bin) {
			mytype = "inflow";
			ndirichlets = x.NV-1;
			dirichlets.resize(x.NV-1);
			for (int n=1;n<x.NV;++n)
				dirichlets(n-1) = n;
		}
		inflow(const inflow& inbdry, tet_hp_cns &xin, face_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		inflow* create(tet_hp& xin, face_bdry &bin) const {return new inflow(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
		
		void vdirichlet() {			
			for(int j=0;j<base.npnt;++j) {
				int v0 = base.pnt(j).gindx;
				x.gbl->res.v(v0,Range(1,x.NV-1)) = 0.0;
			}
		}
		
//		void edirichlet(int mode) {
//			int sind;
//			
//			for(int j=0;j<base.nseg;++j) {
//				sind = base.seg(j).gindx;
//				x.gbl->res.e(sind,mode,Range(0,x.NV-2)) = 0.0;
//			}
//		}
		void edirichlet() {
			if (basis::tet(x.log2p).em > 0) {
				for(int j=0;j<base.nseg;++j) {
					int sind = base.seg(j).gindx;
					x.gbl->res.e(sind,Range::all(),Range(1,x.NV-1)) = 0.0;
				}
			}
		}		
//		void fdirichlet(int mode) {
//			int find;
//			
//			for(int j=0;j<base.ntri;++j) {
//				find = base.tri(j).gindx;
//				x.gbl->res.f(find,mode,Range(0,x.NV-2)) = 0.0;
//			}
//		}
		void fdirichlet() {
			if (basis::tet(x.log2p).fm > 0) {
				for(int j=0;j<base.ntri;++j) {
					int find = base.tri(j).gindx;
					x.gbl->res.f(find,Range::all(),Range(1,x.NV-1)) = 0.0;
				}
			}
		}	
		
		void apply_sparse_dirichlet(bool compressed_column) {
			int gind;
			int em=basis::tet(x.log2p).em;
			int fm=basis::tet(x.log2p).fm;
			
			for(int i=0;i<base.npnt;++i){
				gind = base.pnt(i).gindx*x.NV;
				for(int n=1;n<x.NV;++n)
					x.sparse_dirichlet(gind+n,compressed_column);
			}
			
			for(int i=0;i<base.nseg;++i){
				gind = x.npnt*x.NV+base.seg(i).gindx*em*x.NV;
				for(int m=0; m<em; ++m)
					for(int n=1;n<x.NV;++n)
						x.sparse_dirichlet(gind+m*x.NV+n,compressed_column);
			}
			
			for(int i=0;i<base.ntri;++i){
				gind = x.npnt*x.NV+x.nseg*em*x.NV+base.tri(i).gindx*fm*x.NV;
				for(int m=0; m<fm; ++m)
					for(int n=1;n<x.NV;++n)
						x.sparse_dirichlet(gind+m*x.NV+n,compressed_column);
			}			
		}
		
		void tadvance() {
			hp_face_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		};
		
		void modify_boundary_residual() {
			int j,k,m,n,v0,v1,sind,info;
			TinyVector<FLT,tet_mesh::ND> pt;
			TinyVector<double,MXGP> res1d;
			TinyVector<double,MXTM> rescoef;
			char uplo[] = "U";
			
			FLT ogm1 = 1.0/(x.gbl->gamma-1.0);
			
			for(int j=0;j<base.npnt;++j) {
				int v0 = base.pnt(j).gindx;
				FLT KE = 0.5*(ibc->f(1, x.pnts(v0), x.gbl->time)*ibc->f(1, x.pnts(v0), x.gbl->time)+ibc->f(2, x.pnts(v0), x.gbl->time)*ibc->f(2, x.pnts(v0), x.gbl->time)+ibc->f(3, x.pnts(v0), x.gbl->time)*ibc->f(3, x.pnts(v0), x.gbl->time));
				x.gbl->res.v(v0,1) = x.gbl->res.v(v0,0)*ibc->f(1, x.pnts(v0), x.gbl->time);
				x.gbl->res.v(v0,2) = x.gbl->res.v(v0,0)*ibc->f(2, x.pnts(v0), x.gbl->time);
				x.gbl->res.v(v0,3) = x.gbl->res.v(v0,0)*ibc->f(3, x.pnts(v0), x.gbl->time);
				x.gbl->res.v(v0,4) = x.gbl->res.v(v0,0)*(ibc->f(4, x.pnts(v0), x.gbl->time)*ogm1+KE);
			}

			if(basis::tet(x.log2p).em) {
				for(j=0;j<base.nseg;++j) {
					sind = base.seg(j).gindx;
					v0 = x.seg(sind).pnt(0);
					v1 = x.seg(sind).pnt(1);
					
					if (is_curved()) {
						x.crdtocht1d(sind);
						for(n=0;n<tet_mesh::ND;++n)
							basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0));
					}
					else {
						for(n=0;n<tet_mesh::ND;++n) {
							basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd1d(n)(0));
							
							for(k=0;k<basis::tet(x.log2p).gpx;++k)
								x.dcrd1d(n)(k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
						}
					}
					
					/* take global coefficients and put into local vector */
					/*  only need res_rho */
					for (m=0; m<2; ++m) 
						rescoef(m) = x.gbl->res.v(x.seg(sind).pnt(m),0);					
					
					for (m=0;m<basis::tet(x.log2p).em;++m) 
						rescoef(m+2) = x.gbl->res.e(sind,m,0);					
					
					basis::tet(x.log2p).proj1d(&rescoef(0),&res1d(0));
					
					for(n=1;n<x.NV;++n)
						basis::tet(x.log2p).proj1d(x.gbl->res.v(v0,n),x.gbl->res.v(v1,n),&x.res1d(n)(0));
					
					for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
						pt(0) = x.crd1d(0)(k);
						pt(1) = x.crd1d(1)(k);
						
						FLT KE = 0.5*(ibc->f(1, pt, x.gbl->time)*ibc->f(1, pt, x.gbl->time)+ibc->f(2, pt, x.gbl->time)*ibc->f(2, pt, x.gbl->time)+ibc->f(3, pt, x.gbl->time)*ibc->f(3, pt, x.gbl->time));
						x.res1d(1)(k) -= res1d(k)*ibc->f(1, pt, x.gbl->time);
						x.res1d(2)(k) -= res1d(k)*ibc->f(2, pt, x.gbl->time);
						x.res1d(3)(k) -= res1d(k)*ibc->f(3, pt, x.gbl->time);
						x.res1d(4)(k) -= res1d(k)*(ibc->f(4, pt, x.gbl->time)*ogm1+KE);
						
					}
					
					for(n=1;n<x.NV;++n){
						basis::tet(x.log2p).intgrt1d(&x.lf(n)(0),&x.res1d(n)(0));
						
						for(m=0;m<basis::tet(x.log2p).em;++m) 
							x.gbl->res.e(sind,m,n) = -x.lf(n)(2+m)*basis::tet(x.log2p).diag1d(m);					
						
					}								
				}
			}
						
			return;
		}
	};
	
	
	
	class applied_stress : public neumann {
		Array<symbolic_function<3>,1> stress;

		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {

				
				cout << "doesn't work yet "<< endl;
//				/* CONTINUITY */
//				flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1)+(u(2) -mv(2))*norm(2));
//
//				FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1) +norm(2)*norm(2));
//				/* XYZ MOMENTUM */
//
//#ifdef INERTIALESS
//				for (int n=0;n<tet_mesh::ND;++n)
//					flx(n) = -stress(n).Eval(xpt,x.gbl->time)*length +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
//#else
//				for (int n=0;n<tet_mesh::ND;++n)
//					flx(n) = flx(x.NV-1)*u(n) -stress(n).Eval(xpt,x.gbl->time)*length +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
//#endif
//
//				/* EVERYTHING ELSE */
//				for (int n=tet_mesh::ND;n<x.NV-1;++n)
//					flx(n) = flx(x.NV-1)*u(n);

				return;
			}
		public:
			applied_stress(tet_hp_cns &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "applied_stress";}
			applied_stress(const applied_stress& inbdry, tet_hp_cns &xin, face_bdry &bin) : neumann(inbdry,xin,bin), stress(inbdry.stress) {}
			applied_stress* create(tet_hp& xin, face_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tet_hp_cns&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in);
	};
	
}
#endif
