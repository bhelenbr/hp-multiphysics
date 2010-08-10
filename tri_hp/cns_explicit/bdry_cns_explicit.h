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

#ifndef _bdry_cns_explicit_h_
#define _bdry_cns_explicit_h_


#include "tri_hp_cns_explicit.h"
#include "../hp_boundary.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array.h>
#include <symbolic_function.h>


using namespace blitz;

namespace bdry_cns_explicit {

	class generic : public hp_edge_bdry {
	protected:
		tri_hp_cns_explicit &x;
		bool report_flag;
#ifdef L2_ERROR
		symbolic_function<2> l2norm;
#endif
		enum bctypes {ess, nat, mix};
		
	public:
		Array<FLT,1> total_flux,diff_flux,conv_flux;
		FLT circumference,moment,convect,circulation;
		
	public:
		generic(tri_hp_cns_explicit &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "generic";}
		generic(const generic& inbdry, tri_hp_cns_explicit &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {
			if (report_flag) {
#ifdef L2_ERROR
				l2norm = inbdry.l2norm;
#endif
				total_flux.resize(x.NV);
				diff_flux.resize(x.NV);
				conv_flux.resize(x.NV);  
			}
		}
		generic* create(tri_hp& xin, edge_bdry &bin) const {return new generic(*this,dynamic_cast<tri_hp_cns_explicit&>(xin),bin);}
		void init(input_map& input,void* gbl_in) {
			hp_edge_bdry::init(input,gbl_in);
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
		void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);
	};
	
	
	class neumann : public generic {
	protected:
		virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
			
			// switch u from conservative to primitive 
			Array<FLT,1> cvu(x.NV);
			cvu=u;
			
			double pr = (x.gbl->gamma-1.0)*(ibc->f(3, xpt, x.gbl->time)-0.5/ibc->f(0, xpt, x.gbl->time)*(ibc->f(1, xpt, x.gbl->time)*ibc->f(1, xpt, x.gbl->time)+ibc->f(2, xpt, x.gbl->time)*ibc->f(2, xpt, x.gbl->time)));
			
			u(0) = (x.gbl->gamma-1.0)*(cvu(3)-0.5/cvu(0)*(cvu(1)*cvu(1)+cvu(2)*cvu(2)));
			u(1) = cvu(1)/cvu(0);
			u(2) = cvu(2)/cvu(0);
			u(3) = u(0)/cvu(0);
			//u(3) = pr/cvu(0);

			
			/* CONTINUITY */
			flx(0) = pr/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));
			
			/* X&Y MOMENTUM */
#ifdef INERTIALESS

			for (int n=1;n<tri_mesh::ND+1;++n)
					flx(n) = pr*norm(n-1);
#else

			for (int n=1;n<tri_mesh::ND+1;++n)
				flx(n) = flx(0)*u(n) +pr*norm(n-1);
#endif
			
			/* ENERGY EQUATION */
			double h = x.gbl->gamma/(x.gbl->gamma-1.0)*u(x.NV-1) +0.5*(u(1)*u(1)+u(2)*u(2));
			flx(x.NV-1) = h*flx(0);


			return;
		}
		
	public:
		neumann(tri_hp_cns_explicit &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "neumann";}
		neumann(const neumann& inbdry, tri_hp_cns_explicit &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
		neumann* create(tri_hp& xin, edge_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_cns_explicit&>(xin),bin);}
		void element_rsdl(int eind,int stage);
		
	};


	class inflow : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {
			
			// switch u from conservative to primitive 
			Array<FLT,1> cvu(x.NV);
			cvu=u;
			
			double pr = (x.gbl->gamma-1.0)*(ibc->f(3, xpt, x.gbl->time)-0.5/ibc->f(0, xpt, x.gbl->time)*(ibc->f(1, xpt, x.gbl->time)*ibc->f(1, xpt, x.gbl->time)+ibc->f(2, xpt, x.gbl->time)*ibc->f(2, xpt, x.gbl->time)));

			u(0) = (x.gbl->gamma-1.0)*(cvu(3)-0.5/cvu(0)*(cvu(1)*cvu(1)+cvu(2)*cvu(2)));
			u(1) = cvu(1)/cvu(0);
			u(2) = cvu(2)/cvu(0);
			u(3) = u(0)/cvu(0);
			//u(3) = pr/cvu(0);

			/* CONTINUITY */
			//flx(0) = pr/u(3)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));
			//flx(0) = cvu(0)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));
			//flx(0) = ibc->f(1, xpt, x.gbl->time)*norm(0)+ibc->f(2, xpt, x.gbl->time)*norm(1);
			//flx(0) = ibc->f(0, xpt, x.gbl->time)*(u(1)*norm(0)+u(2)*norm(1));
			//flx(0) = pr/u(3)*(u(1)*norm(0) +u(2)*norm(1));

			flx(0) = pr/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));

			/* EVERYTHING ELSE DOESN'T MATTER */
			for (int n=1;n<x.NV;++n)
				flx(n) = 0.0;

			return;
		}
		
	public:
		inflow(tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(xin,bin) {
			mytype = "inflow";
			ndirichlets = x.NV-1;
			dirichlets.resize(x.NV-1);
			for (int n=1;n<x.NV;++n)
				dirichlets(n-1) = n;
		}
		inflow(const inflow& inbdry, tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		inflow* create(tri_hp& xin, edge_bdry &bin) const {return new inflow(*this,dynamic_cast<tri_hp_cns_explicit&>(xin),bin);}
		
		void vdirichlet() {
			int sind,j,v0;
			j = 0;
			do {
				sind = base.seg(j);
				v0 = x.seg(sind).pnt(0);
				x.gbl->res.v(v0,Range(1,x.NV-1)) = 0.0;
			} while (++j < base.nseg);
			v0 = x.seg(sind).pnt(1);
			x.gbl->res.v(v0,Range(1,x.NV-1)) = 0.0;
		}
		
		void sdirichlet(int mode) {
			int sind;
			
			for(int j=0;j<base.nseg;++j) {
				sind = base.seg(j);
				x.gbl->res.s(sind,mode,Range(1,x.NV-1)) = 0.0;
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
			
			int gind,v0,sind;
			int counter = 0;
			
			int j = 0;
			do {
				sind = base.seg(j);
				v0 = x.seg(sind).pnt(0);
				gind = v0*vdofs;
				for(int n=1;n<x.NV;++n) {						
					indices(counter++)=gind+n;
				}
			} while (++j < base.nseg);
			v0 = x.seg(sind).pnt(1);
			gind = v0*vdofs;
			for(int n=1;n<x.NV;++n) {
				indices(counter++)=gind+n;
			}
			
			for(int i=0;i<base.nseg;++i) {
				gind = x.npnt*vdofs+base.seg(i)*sm*x.NV;
				for(int m=0; m<sm; ++m) {
					for(int n=1;n<x.NV;++n) {
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
		
		void tadvance() {
			hp_edge_bdry::tadvance();
			setvalues(ibc,dirichlets,ndirichlets);
		}
	};
	
	


    class characteristic : public neumann {
		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx);
		public:
			characteristic(tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic& inbdry, tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_cns_explicit&>(xin),bin);}
	};

	class applied_stress : public neumann {
		Array<symbolic_function<2>,1> stress;

		protected:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {

				// switch u from conservative to primitive 
				Array<FLT,1> cvu(x.NV);
				cvu=u;
				
				double pr = (x.gbl->gamma-1.0)*(ibc->f(3, xpt, x.gbl->time)-0.5/ibc->f(0, xpt, x.gbl->time)*(ibc->f(1, xpt, x.gbl->time)*ibc->f(1, xpt, x.gbl->time)+ibc->f(2, xpt, x.gbl->time)*ibc->f(2, xpt, x.gbl->time)));
				
				
				u(0) = (x.gbl->gamma-1.0)*(cvu(3)-0.5/cvu(0)*(cvu(1)*cvu(1)+cvu(2)*cvu(2)));
				u(1) = cvu(1)/cvu(0);
				u(2) = cvu(2)/cvu(0);
				u(3) = u(0)/cvu(0);
				
				/* CONTINUITY */
				flx(0) = pr/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));

				FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
				/* X&Y MOMENTUM */
#ifdef INERTIALESS
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n+1) = -stress(n).Eval(xpt,x.gbl->time)*length +pr*norm(n);
#else
				for (int n=0;n<tri_mesh::ND;++n)
					flx(n+1) = flx(0)*u(n+1) -stress(n).Eval(xpt,x.gbl->time)*length +pr*norm(n);
#endif

				/* ENERGY EQUATION */
				double h = x.gbl->gamma/(x.gbl->gamma-1.0)*u(x.NV-1) +0.5*(u(1)*u(1)+u(2)*u(2));
				flx(x.NV-1) = h*flx(0)-stress(2).Eval(xpt,x.gbl->time)*length;				
				
				return;
			}
		public:
			applied_stress(tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "applied_stress";}
			applied_stress(const applied_stress& inbdry, tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), stress(inbdry.stress) {}
			applied_stress* create(tri_hp& xin, edge_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tri_hp_cns_explicit&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in);
	};


	class inflow_pt : public hp_vrtx_bdry {
	protected:
		tri_hp_cns_explicit &x;
		
	public:
		inflow_pt(tri_hp_cns_explicit &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "inflow_pt";}
		inflow_pt(const inflow_pt& inbdry, tri_hp_cns_explicit &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
		inflow_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new inflow_pt(*this,dynamic_cast<tri_hp_cns_explicit&>(xin),bin);}
		
		void tadvance() { 
			for(int n=1;n<x.NV;++n)
				x.ug.v(base.pnt,n) = x.gbl->ibc->f(n,x.pnts(base.pnt),x.gbl->time);  
			return;
		}
		
		void vdirichlet2d() {
			x.gbl->res.v(base.pnt,Range(1,x.NV-1)) = 0.0;
		}
	};


}
#endif
