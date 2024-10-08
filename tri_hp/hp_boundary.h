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
#include <tri_boundary.h>

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
	
public:
	hp_vrtx_bdry(tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin), ibc(x.hp_gbl->ibc), coupled(false), frozen(false), report_flag(false), jacobian_start(0) {mytype = "plain";}
	hp_vrtx_bdry(const hp_vrtx_bdry &inbdry,tri_hp& xin, vrtx_bdry &bin) : mytype(inbdry.mytype), x(xin), base(bin), adapt_storage(inbdry.adapt_storage), ibc(inbdry.ibc), coupled(inbdry.coupled), frozen(inbdry.frozen),report_flag(inbdry.report_flag),
	type(inbdry.type), essential_indices(inbdry.essential_indices), c0_indices(inbdry.c0_indices), c0_indices_xy(inbdry.c0_indices_xy) {
#ifdef petsc
		base.resize_buffers((x.NV+x.ND)*60*(3 +3*x.sm0+x.im0));  // Allows for 4 elements of jacobian entries to be sent
#endif
	}
	virtual hp_vrtx_bdry* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_vrtx_bdry(*this,xin,bin);}
	virtual void init(input_map& inmap); /**< This is to read definition data only (not solution data) */
	virtual void copy(const hp_vrtx_bdry& tgt) {}
	virtual ~hp_vrtx_bdry() {}
	
	/** This is to read solution data **/
	virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0);
    virtual void input(const std::string& fname,tri_hp::filetype typ,int tlvl = 0);

	/* output functions */
	virtual void output(const std::string& filename, tri_hp::filetype typ,int tlvl = 0);

	/* BOUNDARY CONDITION FUNCTIONS */
	virtual void vdirichlet() {
		for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n)
		x.hp_gbl->res.v(base.pnt,*n)= 0.0;
	}
	virtual void pmatchsolution_snd(FLT *pdata, int vrtstride) {base.vloadbuff(boundary::all,pdata,c0_indices.front(),c0_indices.back(),vrtstride*x.NV);}
	virtual void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,c0_indices.front(),c0_indices.back(), x.NV*vrtstride);}
	
	/* FOR COUPLED DYNAMIC BOUNDARIES */
    virtual int setup_preconditioner() {return(0);}
	virtual void tadvance() {
		int pnt = base.pnt;
		base.mvpttobdry(x.pnts(pnt));
		if (!frozen) {
			for(std::vector<int>::const_iterator n=essential_indices.begin();n != essential_indices.end();++n)
				x.ug.v(base.pnt,*n) = ibc->f(*n,x.pnts(base.pnt),x.gbl->time);
		}
	}
	virtual void calculate_unsteady_sources() {}
	virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {flx = 0.0;} //FIXME: Never been finished for symbolic default
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
	virtual int non_sparse_rcv(int phase,Array<int,1> &nnzero,Array<int,1> &nnzero_mpi);
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
	/* Non-shared data */
	std::string mytype;										/**< Class name */
	tri_hp& x;														/**< Reference to parent */
	edge_bdry &base;											/**< Reference to mesh boundary */
	bool shared_owner, curved, coupled, frozen, report_flag;  /**< Various flags */
	int jacobian_start;  /**< Index for rows of extra degrees of freedom (coupled) */
	Array<TinyVector<FLT,tri_mesh::ND>,2> crv;
	Array<Array<TinyVector<FLT,tri_mesh::ND>,2>,1> crvbd;
	Array<TinyMatrix<FLT,tri_mesh::ND,MXGP>,2> dxdt;
	
	/* Shared data */
	enum bctypes {essential, natural};
	std::vector<bctypes> type;
	std::vector<int> essential_indices, c0_indices, c0_indices_xy; //<! Indices of essential b.c. vars and continuous variables (for communication routines)
	std::vector<vector_function *> fluxes;
	symbolic_function<2> *l2norm;
	bool ibc_owner; /**< Fixme: this is stupid, but I need to be done */
	init_bdry_cndtn *ibc; /**< pointer to initial boundary condition function */
	const hp_edge_bdry *adapt_storage;		/**< mesh adapt storage */
	
public:
	hp_edge_bdry(tri_hp& xin, edge_bdry &bin) : x(xin), base(bin), shared_owner(false), curved(false), coupled(false), frozen(false), report_flag(false), ibc_owner(false), ibc(x.hp_gbl->ibc), adapt_storage(NULL) {mytype = "plain";}
	hp_edge_bdry(const hp_edge_bdry &inbdry, tri_hp& xin, edge_bdry &bin) : mytype(inbdry.mytype), x(xin), base(bin), shared_owner(false), curved(inbdry.curved), coupled(inbdry.coupled), frozen(inbdry.frozen), report_flag(inbdry.report_flag), type(inbdry.type), essential_indices(inbdry.essential_indices), c0_indices(inbdry.c0_indices), c0_indices_xy(inbdry.c0_indices_xy), fluxes(inbdry.fluxes), l2norm(inbdry.l2norm), ibc_owner(false), ibc(inbdry.ibc), adapt_storage(inbdry.adapt_storage) {
		
		if (curved && !x.coarse_level) {
			crv.resize(base.maxseg,x.sm0);
			crvbd.resize(x.gbl->nhist+1);
			for(int i=1;i<x.gbl->nhist+1;++i)
				crvbd(i).resize(base.maxseg,x.sm0);
			crvbd(0).reference(crv);
		}
		
		dxdt.resize(x.log2pmax+1,base.maxseg);
#ifndef petsc
		base.resize_buffers(base.maxseg*(x.sm0+2)*x.NV);
#else
		base.resize_buffers(base.maxseg*(x.sm0+2)*(x.NV+x.ND)*16*(3 +3*x.sm0+x.im0));  // Allows for 4 elements of jacobian entries to be sent
#endif
	}
	virtual hp_edge_bdry* create(tri_hp& xin, edge_bdry &bin) const {return(new hp_edge_bdry(*this,xin,bin));}
	virtual void init(input_map& inmap);
	void find_matching_boundary_name(input_map& inmap, std::string& blockname, std::string& sidename);
	virtual void copy(const hp_edge_bdry& tgt);
	virtual ~hp_edge_bdry() {
		if (shared_owner) {
			for(int n=0;n<x.NV;++n) {
				delete fluxes[n];
			}
			if (ibc_owner) {
				delete ibc;
			}
			delete l2norm;
		}
	}
	
	/** This is to read solution data **/
	virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0);
	virtual void input(const std::string& fname,tri_hp::filetype typ,int tlvl = 0);
	/* input output functions */
	virtual void output(const std::string& fname, tri_hp::filetype typ,int tlvl = 0);
    virtual void output_msh(const std::string& fname, int count_pass = 0);
	virtual void setvalues(init_bdry_cndtn *ibc, const std::vector<int>& indices);
	
	/* CURVATURE FUNCTIONS */
	bool is_curved() {return(curved);}
	FLT& crds(int ind, int mode, int dir) {return(crv(ind,mode)(dir));}
	FLT& crdsbd(int tlvl, int ind, int mode, int dir) {return(crvbd(tlvl)(ind,mode)(dir));}
	void curv_init(int tlvl = 0);
	
	/* BOUNDARY CONDITION FUNCTIONS */
	virtual void maxres() {}
	virtual void vdirichlet();
	virtual void sdirichlet(int mode);
	virtual void pmatchsolution_snd(FLT *pdata, int vrtstride) {
        if (c0_indices.size())
            base.vloadbuff(boundary::all,pdata,c0_indices.front(),c0_indices.back(),vrtstride*x.NV);
    }
	virtual void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {
        if (c0_indices.size())
            base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,c0_indices.front(),c0_indices.back(), x.NV*vrtstride);
    }
	virtual void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride);
	virtual int smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride);
	
	/* FOR COUPLED DYNAMIC BOUNDARIES */
    virtual int setup_preconditioner() {return(0);}
	virtual void tadvance();
	virtual void calculate_unsteady_sources();
    virtual void reset_timestep();
	virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx);
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
	virtual int non_sparse_rcv(int phase,Array<int,1> &nnzero,Array<int,1> &nnzero_mpi);
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
	virtual void updatepdata_bdry(int bel,int endpt,const hp_edge_bdry *bin) {}
	virtual void movepdata_bdry(int bel,int endpt,const hp_edge_bdry *bin) {}
	virtual void updatesdata_bdry(int bel,const hp_edge_bdry *bin) {}
	virtual void movesdata_bdry(int bel,const hp_edge_bdry *tgt, int tgtel = -1) {
		int step,m,n;
		
		if (!curved || !x.sm0 || tgt == NULL) return;
		
		if (tgtel < 0) {
			/* Assumes that sind's are the same */
			tgtel = tgt->x.getbdryseg(tgt->x.seg(base.seg(bel)).tri(1));
		}
		
		for(step=0;step<x.gbl->nadapt;++step) {
			for(m=0;m<x.sm0;++m) {
				for(n=0;n<x.ND;++n) {
					crdsbd(step,bel,m,n) = tgt->crvbd(step)(tgtel,m)(n); // crdsbd(step,tgtel,m,n); 
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
	std::vector<vector_function *> derivative_fluxes;
	public:
		symbolic_with_integration_by_parts(tri_hp &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin) {mytype = "symbolic_with_integration_by_parts";}
	symbolic_with_integration_by_parts(const symbolic_with_integration_by_parts& inbdry, tri_hp &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), derivative_fluxes(inbdry.derivative_fluxes) {}
		symbolic_with_integration_by_parts* create(tri_hp& xin, edge_bdry &bin) const {return new symbolic_with_integration_by_parts(*this,xin,bin);}
		void init(input_map& inmap);
		void element_rsdl(int eind, Array<TinyVector<FLT,MXTM>,1> lf);
		~symbolic_with_integration_by_parts() {
			if (shared_owner) {
				for (int n=0;n<x.NV;++n)
					delete derivative_fluxes[n];
			}
		}
	};

class hp_partition : public hp_edge_bdry {
public:
	/** Array for time history information */
	Array<tri_hp::vsi,1> ugbd;
	Array<Array<TinyVector<FLT,tri_mesh::ND>,1>,1> vrtxbd; //!< Highest level contains pre-summed unsteady mesh velocity source
	Array<tri_hp::vsi,1> dres; //!< Driving term for multigrid
	Array<FLT,2> vug_frst; //!< Solution on first entry to coarse mesh
	
	hp_partition(tri_hp &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin) {
		mytype = "hp_partition";
	}
	hp_partition(const hp_partition& inbdry, tri_hp &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin) {
		epartition& tgt = dynamic_cast<epartition&>(base);
		/** Arrays for time history information for adaptation */
		ugbd.resize(x.gbl->nadapt);;
		vrtxbd.resize(x.gbl->nadapt);; //!< Highest level contains pre-summed unsteady mesh velocity source
		for(int i=0;i<x.gbl->nadapt;++i) {
			ugbd(i).v.resize(tgt.remote_halo.maxpst,x.NV);
			ugbd(i).s.resize(tgt.remote_halo.maxpst,x.sm0,x.NV);
			ugbd(i).i.resize(tgt.remote_halo.maxpst,x.im0,x.NV);
			vrtxbd(i).resize(tgt.remote_halo.maxpst);
		}
#ifndef PETSC
		vug_frst.resize(tgt.remote_halo.maxpst,x.NV);
		dres.resize(1);
		dres(0).v.resize(tgt.remote_halo.maxpst,x.NV);
#endif
	}
	hp_partition* create(tri_hp& xin, edge_bdry &bin) const {return new hp_partition(*this,xin,bin);}
	void init(input_map& inmap);
	void copy(const hp_edge_bdry& tgt);
	void snd_solution();
	void rcv_solution();
};

/* Special communication point at the intersection of blocks with different physics */
class multi_physics_pnt : public hp_vrtx_bdry {
	/* Matching constraints for each match */
	/* Going to simplify to 1-1 matches only for now */
	std::vector<std::vector<std::pair<int,int> > > match_pairs; //<! Indices of continuous variables (for communication routines)
	Array<FLT,1> denom;
public:
	multi_physics_pnt(tri_hp &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin) {mytype = "smulti_physics_pnt";}
	multi_physics_pnt(const multi_physics_pnt& inbdry, tri_hp &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), match_pairs(inbdry.match_pairs) {}
	multi_physics_pnt* create(tri_hp& xin, vrtx_bdry &bin) const {return new multi_physics_pnt(*this,xin,bin);}
	void init(input_map& inmap);
	void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride);
	int non_sparse_rcv(int phase,Array<int,1> &nnzero,Array<int,1> &nnzero_mpi);
	int petsc_matchjacobian_rcv(int phase);
};
#endif

