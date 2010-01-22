/*
 *  hp_boundary.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 9/3/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _hp_boundary_h_
#define _hp_boundary_h_

#include "tri_hp.h"

class hp_edge_bdry;

class hp_vrtx_bdry : public vgeometry_interface<2> {
	protected:
		std::string mytype;
		tri_hp& x;
		vrtx_bdry& base;

	public:
		hp_vrtx_bdry(tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin) {mytype = "plain";}
		hp_vrtx_bdry(const hp_vrtx_bdry &inbdry,tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin), mytype(inbdry.mytype) {}
		virtual void* create_global_structure() {return 0;}
		virtual hp_vrtx_bdry* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_vrtx_bdry(*this,xin,bin);}
		virtual void init(input_map& input,void* gbl_in) {} /**< This is to read definition data only (not solution data) */
		virtual void copy(const hp_vrtx_bdry& tgt) {}
		virtual ~hp_vrtx_bdry() {}

		/* input output functions */
		virtual void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
			switch(typ) {
				case(tri_hp::text):
					fout << base.x.gbl->idprefix << " " << mytype << std::endl;
					break;
				default:
					break;
			}
			return;
		}
		/** This is to read solution data **/
		virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0) {
			std::string idin,mytypein;

			switch(typ) {
				case(tri_hp::text):
					fin >> idin >> mytypein;
					break;
				default:
					break;
			}
			return;
		}


		/* BOUNDARY CONDITION FUNCTIONS */
		virtual void vdirichlet() {}
		virtual void vdirichlet2d() {} //!< SPECIAL CASE OF POINT BOUNDARY CONDITION FOR 2D FIELD
		virtual void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride) {base.vloadbuff(boundary::all,pdata,0,x.NV-1,x.NV*vrtstride);}
		virtual void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,0,x.NV-1,x.NV*vrtstride);}

		/* FOR COUPLED DYNAMIC BOUNDARIES */
		virtual void setup_preconditioner() {}
		virtual void tadvance() {
			int pnt = base.pnt;
			base.mvpttobdry(x.pnts(pnt));
		}
		virtual void calculate_unsteady_sources() {}
		virtual void rsdl(int stage) {}
		virtual void update(int stage) {}
		virtual void mg_restrict() {} 
		virtual void mg_prolongate() {}    
#ifdef petsc
		virtual void petsc_jacobian() {}
		virtual void petsc_jacobian_dirichlet() {}
#endif
};


class hp_edge_bdry : public egeometry_interface<2> {
	public:
		std::string mytype;
		tri_hp& x;
		edge_bdry &base;
		const hp_edge_bdry *adapt_storage;
		init_bdry_cndtn *ibc;
		bool curved, coupled;
		Array<TinyVector<FLT,tri_mesh::ND>,2> crv;
		Array<Array<TinyVector<FLT,tri_mesh::ND>,2>,1> crvbd;
		Array<TinyMatrix<FLT,tri_mesh::ND,MXGP>,2> dxdt;

	public:
		hp_edge_bdry(tri_hp& xin, edge_bdry &bin) : x(xin), base(bin), curved(false), coupled(false) {mytype = "plain"; ibc=x.gbl->ibc;}
		hp_edge_bdry(const hp_edge_bdry &inbdry, tri_hp& xin, edge_bdry &bin) : mytype(inbdry.mytype), x(xin), base(bin), adapt_storage(inbdry.adapt_storage), ibc(inbdry.ibc), curved(inbdry.curved), coupled(inbdry.coupled) {
			if (curved && !x.coarse_level) {
				crv.resize(base.maxseg,x.sm0);
				crvbd.resize(x.gbl->nhist+1);
				for(int i=1;i<x.gbl->nhist+1;++i)
					crvbd(i).resize(base.maxseg,x.sm0);
				crvbd(0).reference(crv);
			}
			dxdt.resize(x.log2pmax+1,base.maxseg);
			base.resize_buffers(base.maxseg*(x.sm0+2)*x.NV);    
		}
		virtual hp_edge_bdry* create(tri_hp& xin, edge_bdry &bin) const {return(new hp_edge_bdry(*this,xin,bin));}
		virtual void* create_global_structure() {return 0;}
		virtual void init(input_map& input,void* gbl_in); 
		virtual void copy(const hp_edge_bdry& tgt);
		virtual ~hp_edge_bdry() {}

		/* input output functions */
		virtual void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);
		/** This is to read solution data **/
		virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0); 
		void setvalues(init_bdry_cndtn *ibc, Array<int,1>& dirichlets, int ndirichlets);

		/* CURVATURE FUNCTIONS */
		bool is_curved() {return(curved);}
		FLT& crds(int ind, int mode, int dir) {return(crv(ind,mode)(dir));}
		FLT& crdsbd(int tlvl, int ind, int mode, int dir) {return(crvbd(tlvl)(ind,mode)(dir));}
		void curv_init(int tlvl = 0);

		/* BOUNDARY CONDITION FUNCTIONS */
		virtual void maxres() {}
		virtual void vdirichlet() {}
		virtual void sdirichlet(int mode) {}
		virtual void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride) {base.vloadbuff(boundary::all,pdata,0,x.NV-1,x.NV*vrtstride);}
		virtual void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,0,x.NV-1,x.NV*vrtstride);}
		virtual void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride); 
		virtual void smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride);

		/* FOR COUPLED DYNAMIC BOUNDARIES */
		virtual void setup_preconditioner() {}
		virtual void tadvance();
		virtual void calculate_unsteady_sources();
		virtual void element_rsdl(int sind,int stage) {}
		virtual void rsdl(int stage);
		virtual void element_jacobian(int sind, Array<FLT,2>& K);
		virtual void jacobian() {}
#ifdef petsc
		virtual void petsc_jacobian();
		virtual void petsc_jacobian_dirichlet(); 
		virtual int petsc_rsdl(Array<FLT,1> res) {}
#endif
		virtual void update(int stage) {}
		virtual void mg_restrict() {} 
		virtual void mg_prolongate() {} 

		/* ADAPTATION FUNCTIONS */
		virtual void updatepdata_bdry(int bel,int endpt,hp_edge_bdry *bin) {}
		virtual void movepdata_bdry(int bel,int endpt,hp_edge_bdry *bin) {}
		virtual void updatesdata_bdry(int bel,hp_edge_bdry *bin) {}
		virtual void movesdata_bdry(int bel,hp_edge_bdry *tgt) {
			int sind,tgtel,step,m,n;

			if (!curved || !x.sm0) return;

			sind = base.seg(bel);
			tgtel = tgt->x.getbdryseg(tgt->x.seg(sind).tri(1));

			for(step=0;step<x.gbl->nadapt;++step) {
				for(m=0;m<x.sm0;++m) {
					for(n=0;n<x.ND;++n) {
						crdsbd(step,bel,m,n) = tgt->crdsbd(step,tgtel,m,n);
					}
				}
			}
			return;
		}

		/* SEARCH FUNCTIONS */
		virtual void findandmovebdrypt(TinyVector<FLT,2>& xp,int &bel,FLT &psi) const;
		virtual void mvpttobdry(int nseg, FLT psi, TinyVector<FLT,tri_mesh::ND> &pt);

		/* SOME UTILITIES */
		void findmax(FLT (*fxy)(TinyVector<FLT,2> &x));
		void findintercept(FLT (*fxy)(TinyVector<FLT,2> &x));
};

#endif
