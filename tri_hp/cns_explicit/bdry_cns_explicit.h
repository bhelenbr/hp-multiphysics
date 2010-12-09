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

#define FIX_DENSITY

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

#ifdef FIX_DENSITY
	
	class inflow : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {
			
			FLT KE = .5*(u(1)*u(1)+u(2)*u(2))/(u(0)*u(0));
			
			flx(0) = 0.0;
			flx(1) = 0.0;
			flx(2) = 0.0; 
			// doesn't work yet
			flx(3) = u(1)/u(0)*(u(3)+(x.gbl->gamma-1.0)*(u(3)-KE))*norm(0)+u(2)/u(0)*(u(3)+(x.gbl->gamma-1.0)*(u(3)-KE))*norm(1);

			flx = 0.0;
			
			return;
		}
		
	public:
		inflow(tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(xin,bin) {
			mytype = "inflow";
			ndirichlets = x.NV-1;
			dirichlets.resize(x.NV-1);
			for (int n=0;n<x.NV-1;++n)
				dirichlets(n) = n;
		}
		inflow(const inflow& inbdry, tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		inflow* create(tri_hp& xin, edge_bdry &bin) const {return new inflow(*this,dynamic_cast<tri_hp_cns_explicit&>(xin),bin);}
		
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

		void update(int stage) {
			int j,k,m,n,v0,v1,sind,indx,info;
			double KE,rho,u,v,RT;
			TinyVector<FLT,tri_mesh::ND> pt;
			
			TinyVector<double,MXGP> u1d;
			TinyVector<double,MXTM> ucoef;
			
			char uplo[] = "U";
			
			j = 0;
			do {
				sind = base.seg(j);
				v0 = x.seg(sind).pnt(0);
				
				u = ibc->f(1, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
				v = ibc->f(2, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
				KE = 0.5*(u*u+v*v);
				RT = (ibc->f(3, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time)-KE)*(x.gbl->gamma-1.0);
							
				rho = x.ug.v(v0,3)/(RT/(x.gbl->gamma-1.0)+KE);
				
				x.ug.v(v0,0) = rho;
				x.ug.v(v0,1) = rho*u;
				x.ug.v(v0,2) = rho*v;

			} while (++j < base.nseg);
			v0 = x.seg(sind).pnt(1);
			
			u = ibc->f(1, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
			v = ibc->f(2, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
			KE = 0.5*(u*u+v*v);
			RT = (ibc->f(3, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time)-KE)*(x.gbl->gamma-1.0);
			
			rho = x.ug.v(v0,3)/(RT/(x.gbl->gamma-1.0)+KE);
			
			x.ug.v(v0,0) = rho;
			x.ug.v(v0,1) = rho*u;
			x.ug.v(v0,2) = rho*v;
			
			if (basis::tri(x.log2p)->sm()) {
				for(j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					v1 = x.seg(sind).pnt(1);
					
					if (is_curved()) {
						x.crdtocht1d(sind);
						for(n=0;n<tri_mesh::ND;++n)
							basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
					}
					else {
						for(n=0;n<tri_mesh::ND;++n) {
							basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
							
							for(k=0;k<basis::tri(x.log2p)->gpx();++k)
								x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
						}
					}
					
					/* take global coefficients and put into local vector */
					/*  only need rhoE */
					for (m=0; m<2; ++m) 
						ucoef(m) = x.ug.v(x.seg(sind).pnt(m),3);					
					
					for (m=0;m<basis::tri(x.log2p)->sm();++m) 
						ucoef(m+2) = x.ug.s(sind,m,3);					
					
					basis::tri(x.log2p)->proj1d(&ucoef(0),&u1d(0));
					
					for(n=0;n<x.NV-1;++n)
						basis::tri(x.log2p)->proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
					
					for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
						pt(0) = x.crd(0)(0,k);
						pt(1) = x.crd(1)(0,k);
						
						u = ibc->f(1, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
						v = ibc->f(2, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
						KE = 0.5*(u*u+v*v);
						RT = (ibc->f(3, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time)-KE)*(x.gbl->gamma-1.0);
						
						rho = u1d(k)/(RT/(x.gbl->gamma-1.0)+KE);
						
						x.res(0)(0,k) -= rho;
						x.res(1)(0,k) -= rho*u;
						x.res(2)(0,k) -= rho*v;
						
						
					}
					for(n=0;n<x.NV-1;++n)
						basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
					
					indx = sind*x.sm0;
					for(n=0;n<x.NV-1;++n) {
						PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(n)(2),basis::tri(x.log2p)->sm(),info);
						for(m=0;m<basis::tri(x.log2p)->sm();++m) 
							x.ug.s(sind,m,n) = -x.lf(n)(2+m);
					}
					
				}
			}
		}
	};
#else
	class inflow : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {
			
			FLT KE = .5*(u(1)*u(1)+u(2)*u(2))/(u(0)*u(0));
			
			flx(0) = 0.0;
			flx(1) = 0.0;
			flx(2) = 0.0; 
			// doesn't work yet
			flx(3) = u(1)/u(0)*(u(3)+(x.gbl->gamma-1.0)*(u(3)-KE))*norm(0)+u(2)/u(0)*(u(3)+(x.gbl->gamma-1.0)*(u(3)-KE))*norm(1);
			
			flx = 0.0;
			
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
		
		void update(int stage) {
			int j,k,m,n,v0,v1,sind,indx,info;
			double KE,rho,u,v,RT;
			TinyVector<FLT,tri_mesh::ND> pt;
			
			TinyVector<double,MXGP> u1d;
			TinyVector<double,MXTM> ucoef;
			
			char uplo[] = "U";
			
			j = 0;
			do {
				sind = base.seg(j);
				v0 = x.seg(sind).pnt(0);
				
				u = ibc->f(1, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
				v = ibc->f(2, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
				KE = 0.5*(u*u+v*v);
				RT = (ibc->f(3, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time)-KE)*(x.gbl->gamma-1.0);
				
				rho = x.ug.v(v0,0);
				
				x.ug.v(v0,1) = rho*u;
				x.ug.v(v0,2) = rho*v;
				x.ug.v(v0,3) = rho*(RT/(x.gbl->gamma-1.0)+KE);
				
			} while (++j < base.nseg);
			v0 = x.seg(sind).pnt(1);
			
			u = ibc->f(1, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
			v = ibc->f(2, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
			KE = 0.5*(u*u+v*v);
			RT = (ibc->f(3, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time)-KE)*(x.gbl->gamma-1.0);
			
			rho = x.ug.v(v0,0);
			
			x.ug.v(v0,1) = rho*u;
			x.ug.v(v0,2) = rho*v;
			x.ug.v(v0,3) = rho*(RT/(x.gbl->gamma-1.0)+KE);
			
			if(basis::tri(x.log2p)->sm()){
				for(j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					v1 = x.seg(sind).pnt(1);
					
					if (is_curved()) {
						x.crdtocht1d(sind);
						for(n=0;n<tri_mesh::ND;++n)
							basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
					}
					else {
						for(n=0;n<tri_mesh::ND;++n) {
							basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
							
							for(k=0;k<basis::tri(x.log2p)->gpx();++k)
								x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
						}
					}
					
					/* take global coefficients and put into local vector */
					/*  only need rhoE */
					for (m=0; m<2; ++m) 
						ucoef(m) = x.ug.v(x.seg(sind).pnt(m),0);					
					
					for (m=0;m<basis::tri(x.log2p)->sm();++m) 
						ucoef(m+2) = x.ug.s(sind,m,0);					
					
					basis::tri(x.log2p)->proj1d(&ucoef(0),&u1d(0));
					
					for(n=1;n<x.NV;++n)
						basis::tri(x.log2p)->proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
					
					for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
						pt(0) = x.crd(0)(0,k);
						pt(1) = x.crd(1)(0,k);
						
						u = ibc->f(1, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
						v = ibc->f(2, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
						KE = 0.5*(u*u+v*v);
						RT = (ibc->f(3, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time)-KE)*(x.gbl->gamma-1.0);
						
						rho = u1d(k);
						
						x.res(1)(0,k) -= rho*u;
						x.res(2)(0,k) -= rho*v;
						x.res(3)(0,k) -= rho*(RT/(x.gbl->gamma-1.0)+KE);
						
						
					}
					for(n=1;n<x.NV;++n)
						basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
					
					indx = sind*x.sm0;
					for(n=1;n<x.NV;++n) {
						PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(n)(2),basis::tri(x.log2p)->sm(),info);
						for(m=0;m<basis::tri(x.log2p)->sm();++m) 
							x.ug.s(sind,m,n) = -x.lf(n)(2+m);
					}
					
				}
			}					
		}
	};

#endif	
	
	class adiabatic : public neumann {  
	protected:
		Array<int,1> dirichlets;
		int ndirichlets;
		void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm,  Array<FLT,1>& flx) {
			
			FLT KE = .5*(u(1)*u(1)+u(2)*u(2))/(u(0)*u(0));
			
			flx(0) = u(1)*norm(0)+u(2)*norm(1);//doesnt work
			flx(1) = 0.0;
			flx(2) = 0.0; 
			flx(3) = u(1)/u(0)*(u(3)+(x.gbl->gamma-1.0)*(u(3)-KE))*norm(0)+u(2)/u(0)*(u(3)+(x.gbl->gamma-1.0)*(u(3)-KE))*norm(1);
			
			return;
		}
		
	public:
		adiabatic(tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(xin,bin) {
			mytype = "adiabatic";
			ndirichlets = x.NV-2;
			dirichlets.resize(ndirichlets);
			for (int n=1;n<x.NV-1;++n)
				dirichlets(n-1) = n;
		}
		adiabatic(const adiabatic& inbdry, tri_hp_cns_explicit &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), ndirichlets(inbdry.ndirichlets) {dirichlets.resize(ndirichlets), dirichlets=inbdry.dirichlets;}
		adiabatic* create(tri_hp& xin, edge_bdry &bin) const {return new adiabatic(*this,dynamic_cast<tri_hp_cns_explicit&>(xin),bin);}
		
		void vdirichlet() {
			int sind,j,v0;
			j = 0;
			do {
				sind = base.seg(j);
				v0 = x.seg(sind).pnt(0);
				x.gbl->res.v(v0,Range(1,x.NV-2)) = 0.0;
			} while (++j < base.nseg);
			v0 = x.seg(sind).pnt(1);
			x.gbl->res.v(v0,Range(1,x.NV-2)) = 0.0;
		}
		
		void sdirichlet(int mode) {
			int sind;
			
			for(int j=0;j<base.nseg;++j) {
				sind = base.seg(j);
				x.gbl->res.s(sind,mode,Range(1,x.NV-2)) = 0.0;
			}
		}
		
#ifdef petsc			
		void petsc_jacobian_dirichlet() {
			hp_edge_bdry::petsc_jacobian_dirichlet();  // Apply deforming mesh stuff
			
			int sm=basis::tri(x.log2p)->sm();
			Array<int,1> indices((base.nseg+1)*(x.NV-2) +base.nseg*sm*(x.NV-2));
			
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
				for(int n=1;n<x.NV-1;++n) {						
					indices(counter++)=gind+n;
				}
			} while (++j < base.nseg);
			v0 = x.seg(sind).pnt(1);
			gind = v0*vdofs;
			for(int n=1;n<x.NV-1;++n) {
				indices(counter++)=gind+n;
			}
			
			for(int i=0;i<base.nseg;++i) {
				gind = x.npnt*vdofs+base.seg(i)*sm*x.NV;
				for(int m=0; m<sm; ++m) {
					for(int n=1;n<x.NV-1;++n) {
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
		
		void update(int stage) {
			int j,k,m,n,v0,v1,sind,indx,info;
			TinyVector<FLT,tri_mesh::ND> pt;
			double u,v;
			TinyVector<double,MXGP> u1d;
			TinyVector<double,MXTM> ucoef;
			
			char uplo[] = "U";
			
			j = 0;
			do {
				sind = base.seg(j);
				v0 = x.seg(sind).pnt(0);
				u=ibc->f(1, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
				v=ibc->f(2, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
				x.ug.v(v0,1) = x.ug.v(v0,0)*u;
				x.ug.v(v0,2) = x.ug.v(v0,0)*v;
				
			} while (++j < base.nseg);
			v0 = x.seg(sind).pnt(1);
			u=ibc->f(1, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
			v=ibc->f(2, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
			x.ug.v(v0,1) = x.ug.v(v0,0)*u;
			x.ug.v(v0,2) = x.ug.v(v0,0)*v;
			
			if(basis::tri(x.log2p)->sm()) {
				for(j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					v1 = x.seg(sind).pnt(1);
					
					if (is_curved()) {
						x.crdtocht1d(sind);
						for(n=0;n<tri_mesh::ND;++n)
							basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
					}
					else {
						for(n=0;n<tri_mesh::ND;++n) {
							basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
							
							for(k=0;k<basis::tri(x.log2p)->gpx();++k)
								x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
						}
					}
					
					/* take global coefficients and put into local vector */
					/*  only need rho */
					for (m=0; m<2; ++m) 
						ucoef(m) = x.ug.v(x.seg(sind).pnt(m),0);					
					
					for (m=0;m<basis::tri(x.log2p)->sm();++m) 
						ucoef(m+2) = x.ug.s(sind,m,0);					
					
					basis::tri(x.log2p)->proj1d(&ucoef(0),&u1d(0));
					
					for(n=1;n<x.NV-1;++n)
						basis::tri(x.log2p)->proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
					
					for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
						pt(0) = x.crd(0)(0,k);
						pt(1) = x.crd(1)(0,k);
						u=ibc->f(1, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
						v=ibc->f(2, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
						x.res(1)(0,k) -= u1d(k)*u;
						x.res(2)(0,k) -= u1d(k)*v;
					}
					
					for(n=1;n<x.NV-1;++n)
						basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
					
					indx = sind*x.sm0;
					for(n=1;n<x.NV-1;++n) {
						PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(n)(2),basis::tri(x.log2p)->sm(),info);
						for(m=0;m<basis::tri(x.log2p)->sm();++m) 
							x.ug.s(sind,m,n) = -x.lf(n)(2+m);
					}
					
				}
			}
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
			x.gbl->res.v(base.pnt,Range(1,x.NV-2)) = 0.0;
		}
	};


}
#endif
