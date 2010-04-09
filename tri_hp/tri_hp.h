/*
 *  tri_hp.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_h_
#define _tri_hp_h_

#include <r_tri_mesh.h>
#include <float.h>
#include <tri_basis.h>
#include <blocks.h>

#ifdef petsc
#include <petscksp.h>
#endif

#ifdef AXISYMMETRIC
#define RAD(r) (r)
#else
#define RAD(r) 1
#endif

#define MAXP 4
#define MXGP MAXP+2
#define MXTM (MAXP+1)*(MAXP+2)/2 

class hp_vrtx_bdry;
class hp_edge_bdry;

class init_bdry_cndtn {
	public:
		virtual FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) = 0;
		virtual void input(input_map &blkdata, std::string idnty) {};
		virtual ~init_bdry_cndtn() {};
};

class tri_hp_helper;

/** This class is just the data storage and nothing for multigrid */
class tri_hp : public r_tri_mesh  {
	public:
		int NV;
		int p0, sm0, im0;  /**> Initialization values */
		int log2p; /**> index of basis to use in global basis::tri array */
		int log2pmax; /**> Initialization value of log2p */
		enum movementtype {fixed,uncoupled_rigid,coupled_rigid,uncoupled_deformable,coupled_deformable} mmovement;

		/* WORK ARRAYS */
		Array<TinyMatrix<FLT,MXGP,MXGP>,1> u,res;
		Array<TinyMatrix<FLT,MXGP,MXGP>,2> du;
		TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> crd;
		TinyMatrix<FLT,MXGP,MXGP> cjcb;
		TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,ND,ND> dcrd;
		Array<TinyVector<FLT,MXTM>,1> uht,lf;
		TinyMatrix<FLT,ND,MXTM> cht, cf;
		TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> mvel; // for local mesh velocity info
		Array<TinyMatrix<FLT,MXGP,MXGP>,2> bdwk;

		/** Stores vertex, side and interior coefficients of solution */
		struct vsi {
			Array<FLT,2> v;
			Array<FLT,3> s;
			Array<FLT,3> i;
		} ug;

		/** vertex boundary information */
		Array<hp_vrtx_bdry *,1> hp_vbdry;
		virtual hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
		/** edge boundary information */
		Array<hp_edge_bdry *,1> hp_ebdry;
		virtual hp_edge_bdry* getnewsideobject(int bnum, input_map &bdrydata); 
		/** object to perform rigid mesh movement */
		tri_hp_helper *helper;

		/** Array for time history information */
		Array<vsi,1> ugbd;
		Array<Array<TinyVector<FLT,ND>,1>,1> vrtxbd; //!< Highest level contains pre-summed unsteady mesh velocity source
		Array<TinyMatrix<FLT,MXGP,MXGP>,3> dugdt; //!< Precalculated unsteady sources at Gauss points
		Array<TinyMatrix<FLT,MXGP,MXGP>,3> dxdt; //!< Precalculated mesh velocity sources at Gauss points

		/* Multigrid stuff needed on each mesh */
		bool isfrst; // FLAG TO SET ON FIRST ENTRY TO COARSE MESH
		bool coarse_flag;   // Flag to indicate coarse level
		Array<vsi,1> dres; //!< Driving term for multigrid
		Array<FLT,2> vug_frst; //!< Solution on first entry to coarse mesh
		FLT fadd; //!< Controls addition of residuals on coarse mesh

		/* THESE THINGS ARE SHARED BY MESHES OF THE SAME BLOCK */
		struct global : public r_tri_mesh::global {

			/**< Pointer to adaptation solution storage 
			* Also used for backwards difference storage in tadvance 
			* could be used for ug0 res and res_r as well? 
			*/
			tri_hp *pstr;  
			FLT curvature_sensitivity;  /**< sensitivity to boundary curvature  */

			/* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
			vsi ug0;

			/** Residual storage for equations */
			vsi res; 

			/* REAL PART FOR RESIDUAL STORAGE */
			vsi res_r;  

			/* RESIDUAL STORAGE FOR ENTRY TO MULTIGRID */
			vsi res0;

			/* PRECONDITIONER  */
			bool diagonal_preconditioner;
			Array<FLT,2> vprcn, sprcn, tprcn;  // Diagonal preconditioner
			Array<FLT,3> vprcn_ut, sprcn_ut, tprcn_ut; // Lower triangle preconditioner

			/* INITIALIZATION AND BOUNDARY CONDITION FUNCTION */
			init_bdry_cndtn *ibc;

			/* Pointers to block storage objects for edge boundary conditions */
			Array<void *,1> ebdry_gbls;

			/* Pointers to block storage objects for pnts boundary conditions */
			Array<void *,1> vbdry_gbls;

			/* Time step factor for different polynomial degree */
			TinyVector<FLT,MXGP> cfl;

		} *gbl;
		virtual init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);
		virtual tri_hp_helper* getnewhelper(input_map& inmap);

		/* FUNCTIONS FOR MOVING GLOBAL TO LOCAL */
		void ugtouht(int tind);
		void ugtouht(int tind,int nhist);
		void ugtouht_bdry(int tind);
		void ugtouht_bdry(int tind, int nhist);
		void ugtouht1d(int sind);
		void ugtouht1d(int sind, int nhist);
		void crdtocht(int tind);
		void crdtocht(int tind, int nhist);
		void crdtocht1d(int sind);
		void crdtocht1d(int sind, int nhist);
		void restouht_bdry(int tind); // USED IN MINVRT

		/* THIS FUNCTION ADDS LF TO GLOBAL VECTORS */
		void lftog(int tind, vsi gvect);

		/* SETUP V/S/T INFO */
		void setinfo();

	public:
		tri_hp() : r_tri_mesh() {}
		virtual tri_hp* create() {return new tri_hp;}
		void* create_global_structure() {return new global;}
		void init(input_map& input, void *gin);
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void tobasis(init_bdry_cndtn *ibc, int tlvl = 0);
		void curvinit();

		/* Input / Output functions */
		enum filetype {tecplot, text, binary, adapt_diagnostic, auxiliary};
		TinyVector<filetype,3> output_type;
		filetype reload_type;
		void input(const std::string &name);
		void input(const std::string &name, filetype type, int tlvl = 0);

		/** Outputs solution in various filetypes */
		void output(const std::string &name, block::output_purpose why);
		void output(const std::string &name, filetype type = tecplot, int tlvl = 0);

		/** Shift to next implicit time step */
		void tadvance();
		virtual void calculate_unsteady_sources();

		/** Makes sure vertex positions on boundaries coinside */
		void matchboundaries();

		/** Setup preconditioner */
		void setup_preconditioner();

		/** Calculate residuals */
		void rsdl() {rsdl(gbl->nstage);}
		virtual void rsdl(int stage); 
		virtual void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im) {
			*gbl->log << "I shouldn't be in generic element_rsdl" << std::endl;
		}
		void jacobian() {} // Not filled this in yet
		virtual void element_jacobian(int tind, Array<FLT,2>&);

		/** Relax solution */  
		void update();
		virtual void minvrt();
		void minvrt_test();

		/** Multigrid cycle */
		void mg_restrict();
		void mg_prolongate();

		/** Print errors */
		FLT maxres();

		/* FUNCTIONS FOR ADAPTION */ 
		void length() {*gbl->log << "using generic length\n";}
		void adapt();
		void copy(const tri_hp &tgt);
		void movepdata(int frm, int to);
		void movepdata_bdry(int bnum,int bel,int endpt);
		void updatepdata(int v);
		void updatepdata_bdry(int bnum,int bel,int endpt);
		void movesdata(int frm, int to);
		void movesdata_bdry(int bnum,int bel);
		void updatesdata(int s);
		void updatesdata_bdry(int bnum,int bel);
		void movetdata(int frm, int to);
		void updatetdata(int t);
		bool findinteriorpt(TinyVector<FLT,2> pt, int &tind, FLT &r, FLT &s);
		void findandmvptincurved(TinyVector<FLT,2>& pt,int &tind, FLT &r, FLT &s);
		bool ptprobe(TinyVector<FLT,ND> xp, Array<FLT,1> uout, int& tind, int tlvl = 0);
		void ptprobe_bdry(int bnum, TinyVector<FLT,ND> xp, Array<FLT,1> uout, int tlvl);

		/* MESSAGE PASSING ROUTINES SPECIALIZED FOR SOLUTION CONTINUITY */
		void vc0load(int phase, FLT *pdata, int vrtstride=1);
		int vc0wait_rcv(int phase,FLT *pdata, int vrtsride=1);
		int vc0rcv(int phase,FLT *pdata, int vrtstride=1);
		void sc0load(FLT *sdata, int bgnmode, int endmode, int modestride);
		int sc0wait_rcv(FLT *sdata, int bgnmode, int endmode, int modestride);
		int sc0rcv(FLT *sdata, int bgnmode, int endmode, int modestride);

		/* Some other utilities */
		void l2error(init_bdry_cndtn *toCompare);
		void findmax(int bnum, FLT (*fxy)(TinyVector<FLT,ND> &x));
		void findintercept(int bnum, FLT (*fxy)(TinyVector<FLT,ND> &x));
		void integrated_averages(Array<FLT,1> a);
		
#ifdef petsc
		/* Sparse stuff */
		void petsc_jacobian();
		void petsc_rsdl();
		void petsc_update();
		void petsc_setup_preconditioner();
		
		void petsc_initialize();
		void petsc_finalize();
		void r_jacobian_dirichlet(Array<int,1> points, int dstart, int dstop);  // Interface to Jacobian for r_tri_mesh (for now)
		void petsc_to_ug();
		void ug_to_petsc();
		void petsc_make_1D_rsdl_vector(Array<FLT,1>);
		
		int jacobian_size;
		Mat  petsc_J;           /* Jacobian matrix */
		Vec  petsc_u,petsc_f;   /* solution,residual */
		KSP  ksp;               /* linear solver context */
		PC   pc;                 /* preconditioner */
#endif	

		virtual ~tri_hp();
};

/* THIS CLASS IS TO ALLOW SPECIAL THINGS LIKE RIGIDLY MOVING MESHES OR PARAMETER CHANGING ETC.. */
class tri_hp_helper {
	protected:
		tri_hp &x;
		bool post_process;

	public:
		tri_hp_helper(tri_hp& xin) : x(xin), post_process(false) {}
		tri_hp_helper(const tri_hp_helper &in_help, tri_hp& xin) : x(xin), post_process(in_help.post_process) {}
		virtual tri_hp_helper* create(tri_hp& xin) { return new tri_hp_helper(*this,xin); }
		virtual ~tri_hp_helper() {};
		virtual void init(input_map& input, std::string idnty) {
			if (!input.get(idnty+"_post_process",post_process)) {
				input.getwdefault("post_process",post_process,false);
			}
		}
		virtual void tadvance() {
			if (post_process) {
				std::ostringstream nstr;
				std::string fname;
				nstr << x.gbl->tstep << std::flush;
				fname = "rstrt" +nstr.str() +"_" +x.gbl->idprefix;
				x.input(fname);
			} 
		}
		virtual void setup_preconditioner() {}
		virtual void rsdl(int stage) {}
		virtual void update(int stage) {}
		virtual void mg_restrict() {}
		virtual void output() {};
		
		/* Stuff to make jacobians for coupled problems */
		virtual int dofs(int dofs) { return 0;}
		virtual void non_sparse(Array<int,1> &nnzero) {}
		virtual void jacobian() {};
		virtual void jacobian_dirichlet() {};
};

namespace basis {
#ifdef AXISYMMETRIC
	extern tri_basis_array<1> tri;
#else
	extern tri_basis_array<0> tri;
#endif
}

#endif
