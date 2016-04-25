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
#include <symbolic_function.h>

//#define DEBUG_JAC

#ifdef DEBUG_JAC
const FLT eps_r = 0.0, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
const FLT eps_r = 1.0e-8, eps_a = 1.0e-8;  /*<< constants for accurate numerical determination of jacobians */
#endif

class hp_edge_bdry;

class hp_vrtx_bdry : public vgeometry_interface<2> {
protected:
	std::string mytype;
	tri_hp& x;
	vrtx_bdry& base;
	const hp_vrtx_bdry *adapt_storage;
	init_bdry_cndtn *ibc;
	bool coupled, frozen, report_flag;
	int jacobian_start;
	enum bctypes {essential, natural};
	std::vector<bctypes> type;
	std::vector<int> essential_indices, c0_indices, c0_indices_xy; //<! Indices of essential b.c. vars and continuous variables (for communication routines)
	virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {flx = 0.0;} //FIXME: Never been finished for symbolic default
	
public:
	hp_vrtx_bdry(tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin), ibc(x.gbl->ibc), coupled(false), frozen(false), report_flag(false), jacobian_start(0) {mytype = "plain"; type.resize(x.NV,natural);}
	hp_vrtx_bdry(const hp_vrtx_bdry &inbdry,tri_hp& xin, vrtx_bdry &bin) : mytype(inbdry.mytype), x(xin), base(bin), adapt_storage(inbdry.adapt_storage), ibc(inbdry.ibc), coupled(inbdry.coupled), frozen(inbdry.frozen),report_flag(inbdry.report_flag),
	type(inbdry.type), essential_indices(inbdry.essential_indices), c0_indices(inbdry.c0_indices), c0_indices_xy(inbdry.c0_indices_xy) {
#ifdef petsc
		base.resize_buffers((x.NV+x.ND)*60*(3 +3*x.sm0+x.im0));  // Allows for 4 elements of jacobian entries to be sent
#endif
	}
	virtual void* create_global_structure() {return 0;}
	virtual hp_vrtx_bdry* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_vrtx_bdry(*this,xin,bin);}
	virtual void init(input_map& inmap,void* gbl_in); /**< This is to read definition data only (not solution data) */
	virtual void copy(const hp_vrtx_bdry& tgt) {}
	virtual ~hp_vrtx_bdry() {}
	
	/* input output functions */
	virtual void output(const std::string& filename, tri_hp::filetype typ,int tlvl = 0) {
		switch(typ) {
			case(tri_hp::text): {
				std::string fname;
				fname = filename +"_" +x.gbl->idprefix +".txt";
				ofstream fout;
				fout.open(fname.c_str(),std::ofstream::out | std::ofstream::app);
				fout << base.idprefix << " " << mytype << std::endl;
				fout.close();
				break;
			}
			case(tri_hp::tecplot): {
				if (report_flag) {
					streamsize oldprecision = (*x.gbl->log).precision(10);
					*x.gbl->log << base.idprefix << " position: " << x.pnts(base.pnt) << std::endl;
					*x.gbl->log << base.idprefix << " value: " << x.ug.v(base.pnt,Range::all()) << std::endl;
					(*x.gbl->log).precision(oldprecision);
				}
				break;
			}
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
	virtual void vdirichlet() {
		for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n)
		x.gbl->res.v(base.pnt,*n)= 0.0;
	}
	virtual void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride) {base.vloadbuff(boundary::all,pdata,c0_indices.front(),c0_indices.back(),vrtstride*x.NV);}
	virtual void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,c0_indices.front(),c0_indices.back(), x.NV*vrtstride);}
	
	/* FOR COUPLED DYNAMIC BOUNDARIES */
	virtual void setup_preconditioner() {}
	virtual void tadvance() {
		int pnt = base.pnt;
		base.mvpttobdry(x.pnts(pnt));
		if (!frozen) {
			for(std::vector<int>::const_iterator n=essential_indices.begin();n != essential_indices.end();++n)
				x.ug.v(base.pnt,*n) = ibc->f(*n,x.pnts(base.pnt),x.gbl->time);
		}
	}
	virtual void calculate_unsteady_sources() {}
	virtual void element_rsdl(Array<FLT,1> lf) {lf = 0.0;}
	virtual void rsdl(int stage);
	virtual void update(int stage) {}
	virtual void mg_restrict() {}
	virtual void mg_prolongate() {}
	virtual void mg_source() {}
	virtual void element_jacobian(Array<FLT,2>& K);
	
#ifdef petsc
	virtual void non_sparse(Array<int,1> &nnzero) {}
	virtual void non_sparse_snd(Array<int,1> &nnzero,Array<int,1> &nnzero_mpi);
	virtual int non_sparse_rcv(Array<int,1> &nnzero,Array<int,1> &nnzero_mpi);
	virtual void petsc_jacobian();
	virtual void petsc_matchjacobian_snd();
	virtual int petsc_matchjacobian_rcv(int phase);
	virtual void petsc_jacobian_dirichlet();
	virtual void petsc_premultiply_jacobian() {}
	virtual int petsc_to_ug(PetscScalar *array) {return 0;}
	virtual void ug_to_petsc(int& ind) {}
	virtual void petsc_make_1D_rsdl_vector(Array<FLT,1> res) {}
#endif
};


class hp_edge_bdry : public egeometry_interface<2> {
public:
	
	
	std::string mytype;										/**< Class name */
	tri_hp& x;														/**< Reference to parent */
	edge_bdry &base;											/**< Reference to mesh boundary */
	const hp_edge_bdry *adapt_storage;		/**< mesh adapt storage */
	init_bdry_cndtn *ibc; /**< pointer to initial boundary condition function */
	bool curved, coupled, frozen, report_flag;  /**< Various flags */
	int jacobian_start;  /**< Index for rows of extra degrees of freedom (coupled) */
	enum bctypes {essential, natural};
	std::vector<bctypes> type;
	std::vector<int> essential_indices, c0_indices, c0_indices_xy; //<! Indices of essential b.c. vars and continuous variables (for communication routines)
	std::vector<vector_function *> fluxes, derivative_fluxes;
	virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx);
	Array<TinyVector<FLT,tri_mesh::ND>,2> crv;
	Array<Array<TinyVector<FLT,tri_mesh::ND>,2>,1> crvbd;
	Array<TinyMatrix<FLT,tri_mesh::ND,MXGP>,2> dxdt;
	symbolic_function<2> l2norm;
	
public:
	hp_edge_bdry(tri_hp& xin, edge_bdry &bin) : x(xin), base(bin), ibc(x.gbl->ibc), curved(false), coupled(false), frozen(false), report_flag(false) {mytype = "plain"; type.resize(x.NV,natural);}
	hp_edge_bdry(const hp_edge_bdry &inbdry, tri_hp& xin, edge_bdry &bin) : mytype(inbdry.mytype), x(xin), base(bin), adapt_storage(inbdry.adapt_storage), ibc(inbdry.ibc),
	curved(inbdry.curved), coupled(inbdry.coupled), frozen(inbdry.frozen), report_flag(inbdry.report_flag), type(inbdry.type), essential_indices(inbdry.essential_indices), c0_indices(inbdry.c0_indices), c0_indices_xy(inbdry.c0_indices_xy), fluxes(inbdry.fluxes), derivative_fluxes(inbdry.derivative_fluxes) {
		if (curved && !x.coarse_level) {
			crv.resize(base.maxseg,x.sm0);
			crvbd.resize(x.gbl->nhist+1);
			for(int i=1;i<x.gbl->nhist+1;++i)
				crvbd(i).resize(base.maxseg,x.sm0);
			crvbd(0).reference(crv);
		}
		if (report_flag)
			l2norm = inbdry.l2norm;
		
		dxdt.resize(x.log2pmax+1,base.maxseg);
#ifndef petsc
		base.resize_buffers(base.maxseg*(x.sm0+2)*x.NV);
#else
		base.resize_buffers(base.maxseg*(x.sm0+2)*(x.NV+x.ND)*16*(3 +3*x.sm0+x.im0));  // Allows for 4 elements of jacobian entries to be sent
#endif
	}
	virtual hp_edge_bdry* create(tri_hp& xin, edge_bdry &bin) const {return(new hp_edge_bdry(*this,xin,bin));}
	virtual void* create_global_structure() {return 0;}
	virtual void init(input_map& inmap,void* gbl_in);
	void find_matching_boundary_name(input_map& inmap, std::string& blockname, std::string& sidename);
	virtual void copy(const hp_edge_bdry& tgt);
	virtual ~hp_edge_bdry() {}
	
	/* input output functions */
	virtual void output(const std::string& fname, tri_hp::filetype typ,int tlvl = 0);
	/** This is to read solution data **/
	virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0);
	void setvalues(init_bdry_cndtn *ibc, const std::vector<int>& indices);
	
	/* CURVATURE FUNCTIONS */
	bool is_curved() {return(curved);}
	FLT& crds(int ind, int mode, int dir) {return(crv(ind,mode)(dir));}
	FLT& crdsbd(int tlvl, int ind, int mode, int dir) {return(crvbd(tlvl)(ind,mode)(dir));}
	void curv_init(int tlvl = 0);
	
	/* BOUNDARY CONDITION FUNCTIONS */
	virtual void maxres() {}
	virtual void vdirichlet();
	virtual void sdirichlet(int mode);
	virtual void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride) {base.vloadbuff(boundary::all,pdata,c0_indices.front(),c0_indices.back(),vrtstride*x.NV);}
	virtual void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,c0_indices.front(),c0_indices.back(), x.NV*vrtstride);}
	virtual void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride);
	virtual int smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride);
	
	/* FOR COUPLED DYNAMIC BOUNDARIES */
	virtual void setup_preconditioner() {}
	virtual void tadvance();
	virtual void calculate_unsteady_sources();
	virtual void element_rsdl(int eind, Array<TinyVector<FLT,MXTM>,1> lf);
	virtual void rsdl(int stage);
	virtual void element_jacobian(int sind, Array<FLT,2>& K);
	virtual int dofs(int start) {
		jacobian_start = start;
		if (curved && coupled)
			return(tri_mesh::ND*x.sm0*base.nseg);
		return 0;
	}
	virtual void jacobian() {}
#ifdef petsc
	virtual void non_sparse(Array<int,1> &nnzero) {}
	virtual void non_sparse_snd(Array<int,1> &nnzero,Array<int,1> &nnzero_mpi);
	virtual int non_sparse_rcv(Array<int,1> &nnzero,Array<int,1> &nnzero_mpi);
	virtual void petsc_jacobian();
	virtual void petsc_matchjacobian_snd();
	virtual int petsc_matchjacobian_rcv(int phase);
	virtual void petsc_jacobian_dirichlet();
	virtual void petsc_premultiply_jacobian() {}
	virtual int petsc_to_ug(PetscScalar *array);
	virtual void ug_to_petsc(int& ind);
	virtual void petsc_make_1D_rsdl_vector(Array<FLT,1> res) {}
#endif
	virtual void update(int stage) {}
	virtual void modify_boundary_residual() {}
	virtual void mg_restrict() {}
	virtual void mg_prolongate() {}
	virtual void mg_source() {}
	
	
	/* ADAPTATION FUNCTIONS */
	virtual void updatepdata_bdry(int bel,int endpt,hp_edge_bdry *bin) {}
	virtual void movepdata_bdry(int bel,int endpt,hp_edge_bdry *bin) {}
	virtual void updatesdata_bdry(int bel,hp_edge_bdry *bin) {}
	virtual void movesdata_bdry(int bel,hp_edge_bdry *tgt, int tgtel = -1) {
		int step,m,n;
		
		if (!curved || !x.sm0) return;
		
		if (tgtel < 0) {
			/* Assumes that sind's are the same */
			tgtel = tgt->x.getbdryseg(tgt->x.seg(base.seg(bel)).tri(1));
		}
		
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

class symbolic_with_integration_by_parts : public hp_edge_bdry {
	public:
		symbolic_with_integration_by_parts(tri_hp &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin) {mytype = "symbolic_with_integration_by_parts"; derivative_fluxes.resize(x.NV);
			Array<string,1> names(4);
			Array<int,1> dims(4);
			dims = x.ND;
			names(0) = "u";
			dims(0) = x.NV;
			names(1) = "x";
			names(2) = "xt";
			names(3) = "n";
			for(int n=0;n<x.NV;++n) {
				derivative_fluxes[n] = new vector_function(4,dims,names);
			}
		}
		symbolic_with_integration_by_parts(const symbolic_with_integration_by_parts& inbdry, tri_hp &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin) {}
		symbolic_with_integration_by_parts* create(tri_hp& xin, edge_bdry &bin) const {return new symbolic_with_integration_by_parts(*this,xin,bin);}
		void init(input_map& inmap,void* gbl_in);
		void element_rsdl(int eind, Array<TinyVector<FLT,MXTM>,1> lf);
	};
#endif

