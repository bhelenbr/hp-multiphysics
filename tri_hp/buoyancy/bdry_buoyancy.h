/*
 *  bdry_buoyancy.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 12/13/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */
 
#ifndef _bdry_buoyancy_h_
#define _bdry_buoyancy_h_


#include "tri_hp_buoyancy.h"
#include "../ins/bdry_ins.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array.h>
#include <symbolic_function.h>


using namespace blitz;

namespace bdry_buoyancy {
	
	class surface : public bdry_ins::surface {
		Array<vector_function,1> fluxes;
		symbolic_function<1> sigma_vs_T;
		
		public:
			surface(tri_hp_ins &xin, edge_bdry &bin) : bdry_ins::surface(xin,bin) {mytype = "surface"; fluxes.resize(x.NV);
				Array<string,1> names(4);
				Array<int,1> dims(4);
				dims = x.ND;
				names(0) = "u";
				dims(0) = x.NV;
				names(1) = "x";
				names(2) = "xt";
				names(3) = "n";
				for(int n=0;n<x.NV;++n) {
					fluxes(n).set_arguments(4,dims,names);
				}
			}
			surface(const surface& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : bdry_ins::surface(inbdry,xin,bin), fluxes(inbdry.fluxes) {}
			surface* create(tri_hp& xin, edge_bdry &bin) const {return new surface(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			
			void init(input_map& input,void* gbl_in); 
			void element_rsdl(int sind, Array<FLT,2> lf);  // FIXME Not really compatible need to make all consistent
	};
	
	class melt : public bdry_ins::flexible {	
		protected:
			tri_hp_buoyancy &x;
			Array<FLT,1> ksprg;
			Array<TinyVector<FLT,tri_mesh::ND>,1> vug_frst;
			Array<TinyVector<FLT,tri_mesh::ND>,2> vdres; //!< Driving term for multigrid (log2p, pnts)
			Array<TinyVector<FLT,tri_mesh::ND>,3> sdres; //!< Driving term for multigrid (log2p, side, order)
			const melt *fine, *coarse;
			
		public:
			struct global {                
				bool is_loop;
				/* PROPERTIES */
				FLT Lf, rho_s, cp_s;
				
				/* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
				Array<TinyVector<FLT,tri_mesh::ND>,1> vug0;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sug0;
				
				/* RESIDUALS */
				Array<TinyVector<FLT,tri_mesh::ND>,1> vres;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sres;
				Array<TinyVector<FLT,tri_mesh::ND>,1> vres0;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sres0;
								
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
			melt(tri_hp_buoyancy &xin, edge_bdry &bin) : flexible(xin,bin), x(xin) {
				mytype = "melt";
			}
			melt(const melt& inbdry, tri_hp_buoyancy &xin, edge_bdry &bin)  : flexible(inbdry,xin,bin), x(xin) {
				gbl = inbdry.gbl;
				ksprg.resize(base.maxseg);
				vug_frst.resize(base.maxseg+1);
				vdres.resize(1,base.maxseg+1);
				fine = &inbdry;
			}
			melt* create(tri_hp& xin, edge_bdry &bin) const {return new melt(*this,dynamic_cast<tri_hp_buoyancy&>(xin),bin);}
			
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
			void vdirichlet();
#ifdef petsc
//			void petsc_matchjacobian_snd();
//			void petsc_matchjacobian_rcv(int phase);
			int petsc_rsdl(Array<FLT,1> res);
			void petsc_jacobian();
			void petsc_jacobian_dirichlet();
			void non_sparse(Array<int,1> &nnzero);
//			void non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
//			void non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
#endif
		
			/* For matching with solid phase not working yet */
			//void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride) {base.vloadbuff(boundary::all,pdata,0,x.NV-2,vrtstride*x.NV);}
//			void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,0,x.NV-2, x.NV*vrtstride);}
//			void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride); 
//			void smatchsolution_rcv(FLT *sdata, int bgnmode, int endmode, int modestride);

		
	};
	
	
	/*******************************/
	/* VERTEX BOUNDARY CONDITIONS */
	/******************************/
	class melt_end_pt : public hp_vrtx_bdry {
		/* INTERSECTING BOUNDARY CONTAINING END POINT MUST HAVE GEOMETRY NOT BE DEFINED SOLELY BY MESH */
	protected:
		tri_hp_buoyancy &x;
		melt *surf;
		int surfbdry;
		FLT position;
		enum {vertical, horizontal, angled, curved} wall_type;
		TinyVector<FLT,tri_mesh::ND> wall_normal; 
		
	public:
		melt_end_pt(tri_hp_buoyancy &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "melt_end_pt";}
		melt_end_pt(const melt_end_pt& inbdry, tri_hp_buoyancy &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin), surfbdry(inbdry.surfbdry),
			position(inbdry.position), wall_type(inbdry.wall_type), wall_normal(inbdry.wall_normal) {
			if (!(surf = dynamic_cast<melt *>(x.hp_ebdry(base.ebdry(surfbdry))))) {
				*x.gbl->log << "something's wrong can't find surface boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
		}
		melt_end_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new melt_end_pt(*this,dynamic_cast<tri_hp_buoyancy&>(xin),bin);}
		
		void init(input_map& inmap,void* gbl_in) {
			std::string keyword,val;
			std::istringstream data;
			std::string filename,input;
			
			hp_vrtx_bdry::init(inmap,gbl_in);
			
			if ((surf = dynamic_cast<melt *>(x.hp_ebdry(base.ebdry(0))))) {
				surfbdry = 0;
			}
			else if ((surf = dynamic_cast<melt *>(x.hp_ebdry(base.ebdry(1))))) {
				surfbdry = 1;
			}
			else {
				*x.gbl->log << "something's wrong neither side is a surface boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
					
			if (inmap.get(base.idprefix + "_wall_type",input)) {
				if (input == "vertical") {
					wall_type = vertical;
					position = x.pnts(base.pnt)(0);
				}
				else if (input == "horizontal") {
					wall_type = horizontal;
					position = x.pnts(base.pnt)(1);
				}
				else if (input == "angled")
					wall_type = angled;
				else if (input == "curved")
					wall_type = curved;
				else {
					*x.gbl->log << "Unrecognized wall type" << std::endl;
				}
			}
			else {
				wall_type = vertical;
				position = x.pnts(base.pnt)(0);
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

		void rsdl(int stage) {
			/* Fix me: Need to do more than this to get a good approximation to new point location */
			if (surfbdry == 0) {
				/* SET TANGENT RESIDUAL TO ZERO */
				surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(0) = 0.0;
			}
			else {
				/* SET TANGENT RESIDUAL TO ZERO */
				surf->gbl->vres(0)(0) = 0.0;
			}
			return;
		}
		
		void vdirichlet() {}
		
		void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
			switch(wall_type) {
				case(vertical):
					pt(0) = position;
					break;
				case(horizontal):
					pt(1) = position;
					break;
				case(angled):case(curved):
					if (surfbdry == 0) {
						x.ebdry(base.ebdry(1))->mvpttobdry(0,-1.0,pt);
					}
					else {
						x.ebdry(base.ebdry(0))->mvpttobdry(x.ebdry(base.ebdry(0))->nseg-1,1.0,pt);
					}
					break;
			}
		}

#ifdef petsc
		void petsc_jacobian() {} 
		void petsc_jacobian_dirichlet();
#endif
	};
	
	/*******************************/
	/* VERTEX BOUNDARY CONDITIONS */
	/******************************/
	class melt_inflow_pt : public hp_vrtx_bdry {
		/* INTERSECTING BOUNDARY CONTAINING END POINT MUST HAVE GEOMETRY NOT BE DEFINED SOLELY BY MESH */
		protected:
			tri_hp_buoyancy &x;
			melt *surf;
			int surfbdry;
			TinyVector<FLT,tri_mesh::ND> position; 
			
		public:
			melt_inflow_pt(tri_hp_buoyancy &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "melt_inflow_pt";}
			melt_inflow_pt(const melt_inflow_pt& inbdry, tri_hp_buoyancy &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin), surfbdry(inbdry.surfbdry),
			position(inbdry.position) {
				if (!(surf = dynamic_cast<melt *>(x.hp_ebdry(base.ebdry(surfbdry))))) {
					*x.gbl->log << "something's wrong can't find surface boundary" << std::endl;
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
			}
			melt_inflow_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new melt_inflow_pt(*this,dynamic_cast<tri_hp_buoyancy&>(xin),bin);}
			
			void init(input_map& inmap,void* gbl_in) {
				std::string keyword,val;
				std::istringstream data;
				std::string filename,input;
				
				hp_vrtx_bdry::init(inmap,gbl_in);
				
				if ((surf = dynamic_cast<melt *>(x.hp_ebdry(base.ebdry(0))))) {
					surfbdry = 0;
				}
				else if ((surf = dynamic_cast<melt *>(x.hp_ebdry(base.ebdry(1))))) {
					surfbdry = 1;
				}
				else {
					*x.gbl->log << "something's wrong neither side is a surface boundary" << std::endl;
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
				position = x.pnts(base.pnt);
			}
			
			void rsdl(int stage) {
				/* Fix me: Need to do more than this to get a good approximation to new point location */
				if (surfbdry == 0) {
					/* SET TANGENT RESIDUAL TO ZERO */
					surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(0) = 0.0;
				}
				else {
					/* SET TANGENT RESIDUAL TO ZERO */
					surf->gbl->vres(0)(0) = 0.0;
				}
				return;
			}
			
			void vdirichlet() {
				if (surfbdry == 0) {
					for(int n=0;n<tri_mesh::ND;++n) 
						surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(n) = 0.0;
				}
				else {
					for(int n=0;n<tri_mesh::ND;++n) 
						surf->gbl->vres(0)(n) = 0.0;
				}
#ifdef petsc
				r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
				r_gbl->res(base.pnt)(0) = 0.0;
				r_gbl->res(base.pnt)(1) = 0.0;
#endif	
			}
			
			void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
				pt = position;
			}
			
	#ifdef petsc
			void petsc_jacobian() {} 
			void petsc_jacobian_dirichlet() {
				Array<int,1> rows(tri_mesh::ND);
				for(int n=0;n<tri_mesh::ND;++n) 
					rows(n) = (x.NV+tri_mesh::ND)*base.pnt +x.NV +n; 
#ifdef MY_SPARSE
				x.J.zero_rows(tri_mesh::ND,rows);
				x.J_mpi.zero_rows(tri_mesh::ND,rows);
				x.J.set_diag(tri_mesh::ND,rows,1.0);
#else				
				MatZeroRows(x.petsc_J,tri_mesh::ND,rows.data(),1.0);
#endif
			}
				
	#endif
		};
	
}
#endif

