/*
 *  tri_hp.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#ifndef _tet_hp_h_
#define _tet_hp_h_

//#define NODAL

#ifdef NODAL
#include <tet_nodal_basis.h>
#else
#include <tet_basis.h>
#endif

#include <tet_mesh.h>
#include <float.h>
#include <blocks.h>
#include "slu_ddefs.h"

#ifdef petsc
#include <petscksp.h>
#endif

class hp_vrtx_bdry;
class hp_edge_bdry;
class hp_face_bdry;

class init_bdry_cndtn {
	public:
		virtual FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) = 0;
		virtual void input(input_map &blkdata, std::string idnty) {};
};

class tet_hp_helper;

/** This class is just the data storage and nothing for multigrid */
class tet_hp : public tet_mesh  {
	public:
		int NV;
		int p0, em0, fm0, im0;  /**> Initialization values */
		int log2p; /**> index of basis to use in global basis::tri array */
		int log2pmax; /**> Initialization value of log2p */
		enum movementtype {fixed,uncoupled_rigid,coupled_rigid,uncoupled_deformable,coupled_deformable} mmovement;
		
		/* STATIC WORK ARRAYS */
		Array<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,2> du;
		TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP> cjcb;
		TinyVector<TinyVector<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,ND>,ND> dcrd;
		TinyMatrix<FLT,ND,MXTM>  cf;
		TinyVector<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,ND> mvel; // for local mesh velocity info
		Array<TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP>,3> bdwk;	
		
		TinyVector<TinyVector<FLT,MXTM>,ND> cht; // used in crd to cht routines
		Array<TinyVector<FLT,MXTM>,1> uht; // used in ug to uht routines
		Array<TinyVector<FLT,MXTM>,1> lf;  // used in intgrt routines
		TinyVector<TinyVector<double,MXGP>,ND> crd1d,dcrd1d;
		Array<TinyVector<double,MXGP>,1> u1d,res1d;   
		TinyVector<TinyVector<TinyVector<double,MXGP>,MXGP>,ND> crd2d;
		TinyVector<TinyVector<TinyVector<TinyVector<double,MXGP>,MXGP>,2>,ND> dcrd2d;

		Array<TinyVector<TinyVector<double,MXGP>,MXGP>,1> u2d,res2d;   
		TinyVector<TinyVector<TinyVector<TinyVector<double,MXGP>,MXGP>,MXGP>,ND> crd;
		Array<TinyVector<TinyVector<TinyVector<double,MXGP>,MXGP>,MXGP>,1> u,res;  

		Array<Array<FLT,2>,1> spkmass; // mini mass matrices for every vertex
		Array<FLT,1> spkres; // spoke residual size maxnspk
		Array<Array<int,1>,1> spkpiv;
		Array<Array<TinyVector<int,2>,1>,1> spklink; // spoke list for each vertex
		Array<TinyVector<FLT,2>,3> wkseg;


		/** Stores vertex, edge, face, and interior coefficients of solution */
		struct vefi {
			Array<FLT,2> v;
			Array<FLT,3> e;
			Array<FLT,3> f;
			Array<FLT,3> i;
		} ug;
		
		/** vertex boundary information */
		Array<hp_vrtx_bdry *,1> hp_vbdry;
		virtual hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
		/** side boundary information */
		Array<hp_edge_bdry *,1> hp_ebdry;
		virtual hp_edge_bdry* getnewedgeobject(int bnum, input_map &bdrydata); 
		/** face boundary information */
		Array<hp_face_bdry *,1> hp_fbdry;
		virtual hp_face_bdry* getnewfaceobject(int bnum, input_map &bdrydata);
		/** object to perform rigid mesh movement */
		tet_hp_helper *helper;
		
		/** Array for time history information */
		Array<vefi,1> ugbd;
		Array<Array<TinyVector<FLT,ND>,1>,1> vrtxbd; //!< Highest level contains pre-summed unsteady mesh velocity source
		Array<TinyVector<TinyVector<TinyVector<double,MXGP>,MXGP>,MXGP>,3> dugdt; //!< Precalculated unsteady sources at Gauss points
		Array<TinyVector<TinyVector<TinyVector<double,MXGP>,MXGP>,MXGP>,3> dxdt; //!< Precalculated mesh velocity sources at Gauss points
		
		/* Multigrid stuff needed on each mesh */
		bool isfrst; // FLAG TO SET ON FIRST ENTRY TO COARSE MESH
		bool coarse_flag;   // Flag to indicate coarse level
		Array<vefi,1> dres; //!< Driving term for multigrid
		Array<FLT,2> vug_frst; //!< Solution on first entry to coarse mesh
		FLT fadd; //!< Controls addition of residuals on coarse mesh
		
		/* THESE THINGS ARE SHARED BY MESHES OF THE SAME BLOCK */
		struct global : public tet_mesh::global {
			
			/**< Pointer to adaptation solution storage 
			* Also used for backwards difference storage in tadvance 
			* could be used for ug0 res and res_r as well? 
			*/
			tet_hp *pstr;  
			FLT curvature_sensitivity;  /**< sensitivity to boundary curvature  

			/* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
			vefi ug0;

			/** Residual storage for equations */
			vefi res; 

			/* REAL PART FOR RESIDUAL STORAGE */
			vefi res_r;  
			
			/* RESIDUAL STORAGE FOR ENTRY TO MULTIGRID */
			vefi res0;

			/* PRECONDITIONER  */
			bool diagonal_preconditioner;
			Array<FLT,2> vprcn,eprcn,fprcn,iprcn;  // Diagonal preconditioner
			Array<FLT,3> vprcn_ut, eprcn_ut,fprcn_ut,iprcn_ut; // Lower triangle preconditioner
			
			/* INITIALIZATION AND BOUNDARY CONDITION FUNCTION */
			init_bdry_cndtn *ibc;

			/* Pointers to block storage objects for face boundary conditions */
			Array<void *,1> fbdry_gbls;
		 
			/* Pointers to block storage objects for edge boundary conditions */
			Array<void *,1> ebdry_gbls;
			
			/* Pointers to block storage objects for vrtx boundary conditions */
			Array<void *,1> vbdry_gbls;
			
			/* Time step factor for different polynomial degree */
			TinyVector<FLT,MXGP> cfl;
			
		} *gbl;
		virtual init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);
		virtual tet_hp_helper* getnewhelper(input_map& inmap);

		/* FUNCTIONS FOR MOVING GLOBAL TO LOCAL */
		void ugtouht(int tind);
		void ugtouht(int tind,int nhist);
		void ugtouht_bdry(int tind);
		void ugtouht_bdry(int tind, int nhist);
		void ugtouht2d(int find);
		void ugtouht2d(int find, int nhist);
		void ugtouht2d_bdry(int find);
		void ugtouht2d_bdry(int find, int nhist);
		void ugtouht1d(int eind);
		void ugtouht1d(int eind, int nhist);
		void crdtocht(int tind);
		void crdtocht(int tind, int nhist);
		void crdtocht2d(int find);
		void crdtocht2d(int find, int nhist);
		void crdtocht1d(int eind);
		void crdtocht1d(int eind, int nhist);
		void restouht_bdry(int tind); // USED IN MINVRT

		/* THIS FUNCTION ADDS LF TO GLOBAL VECTORS */
		void lftog(int tind, vefi gvect);

		/* SETUP V/S/T INFO */
		void setinfo();
		
	public:
		tet_hp() : tet_mesh() {}
		virtual tet_hp* create() {return new tet_hp;}
		void* create_global_structure() {return new global;}
		void init(input_map& input, void *gin);
		void init_post_findmatch();
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void copy(const tet_hp& tgt);

		void tobasis(init_bdry_cndtn *ibc, int tlvl = 0);
		void curvinit();

		/* Input / Output functions */
		enum filetype {tecplot, text, binary, adapt_diagnostic, auxiliary, trunc_error};
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
		

		/** Relax solution */  
		void update();
		void minvrt();
		void minvrt_test();
		void minvrt_iter();
		void minvrt_direct();
		void test();

		void spoke();
		void particle();
		void getudx(const TinyVector<FLT,3> &xp, int &tind,TinyVector<FLT,4> &uout,TinyVector<FLT,4> &dudx,TinyVector<FLT,4> &dudy,TinyVector<FLT,4> &dudz);
		void getu(const TinyVector<FLT,3> &xp, int &tind,TinyVector<FLT,4> &uout);

		/** Multigrid cycle */
		void mg_restrict();
		void mg_prolongate();

		/** Print errors */
		FLT maxres();
	
		/** Sparse stuff */
		void create_jacobian(bool jac_tran = true);
		void create_local_jacobian_matrix(int tind, Array<FLT,2> &K);
		void create_rsdl();
		void create_local_rsdl(int tind, Array<FLT,1> &lclres);

#ifdef petsc
		void petsc_initialize();
		void petsc_solve();
		void petsc_finalize();
		void petsc_to_ug();
		void ug_to_petsc();
		Mat  petsc_J;           /* Jacobian matrix */
		Vec  petsc_u,petsc_f;   /* solution,residual */
		KSP  ksp;               /* linear solver context */
		PC   pc;                 /* preconditioner */

#endif	

#ifndef petsc
	
		bool sparse_resized;	
		Array<FLT,1> res_vec;//residual vector
		Array<FLT,1> ug_vec;// solution vector
		//Array<int,1> ija; //sparse matrix integer storage
		//Array<FLT,1> sa; //sparse matrix element storage
	
		Array<int,1> sparse_ind; //sparse matrix index
		Array<int,1> sparse_ptr; //sparse matrix pointer
		Array<FLT,1> sparse_val; //sparse matrix element storage
	
		int number_sparse_elements; 
		int size_sparse_matrix;
	
		void insert_sparse(int row, int col, FLT value, bool compressed_column=true);//default to compressed column storage
		void sparse_dirichlet(int ind, bool compressed_col=true);//default to compressed column storage
		void apply_neumman(bool jac_tran=true);
		void find_sparse_bandwidth();
		void zero_sparse();
		void initialize_sparse();
		void vec_to_ug();
		void ug_to_vec();
		
		void superlu();
		void fgmres(int n, SuperMatrix *A, SuperMatrix *L, SuperMatrix *U, 
					double *sol, double *rhs,int *perm_c, int *perm_r,
					FLT tol, int im, int *itmax);
//		SuperMatrix A;
//		SuperMatrix L;
//		SuperMatrix U;
//		SuperMatrix B;
	
//		void lsolver();
//		struct jacobian_matrix{
//			jacobian_matrix() {}
//		};	
//		void mult(const jacobian_matrix &J, const double *vec_in, double *vec_out);

//		int dfgmr(int n, void (*dmatvec) (double, double[], double, double[]), void (*dpsolve) (int, double[], double[]),double *rhs, double *sol, double tol, int im, int *itmax, FILE * fits)
//		void ilu_gmres();/* incomplete lu with gmres solve */
//		int GLOBAL_PERM_C;
//		int GLOBAL_PERM_R;
//		SuperMatrix GLOBAL_A;
//		SuperMatrix GLOBAL_L;
//		SuperMatrix GLOBAL_U;
//		SuperLUStat_t GLOBAL_STAT;
//		void dmatvec_mult(double alpha, double x[], double beta, double y[]);/* matrix vector multiplication */
//		void dpsolve(int n, double x[], double y[]);

#endif
	



		
		/* FUNCTIONS FOR ADAPTION */ 
        void length() {*gbl->log << "using generic length\n";}
//      void adapt();
//      void copy(const tri_hp &tgt);
//      void movevdata(int frm, int to);
//      void movevdata_bdry(int bnum,int bel,int endpt);
//      void updatevdata(int v);
//      void updatevdata_bdry(int bnum,int bel,int endpt);
//      void movesdata(int frm, int to);
//      void movesdata_bdry(int bnum,int bel);
//      void updatesdata(int s);
//      void updatesdata_bdry(int bnum,int bel);
//      void movetdata(int frm, int to);
//      void updatetdata(int t);
		int findinteriorpt(TinyVector<FLT,3> pt, int &tind, FLT &r, FLT &s, FLT &t);
//      void findandmvptincurved(TinyVector<FLT,2>& pt,int &tind, FLT &r, FLT &s);
		void ptprobe(TinyVector<FLT,ND> xp, Array<FLT,1> uout, int tlvl);
//      void ptprobe_bdry(int bnum, TinyVector<FLT,ND> xp, Array<FLT,1> uout, int tlvl);

		/* MESSAGE PASSING ROUTINES SPECIALIZED FOR SOLUTION CONTINUITY */
		void pc0load(int phase, FLT *vdata, int vrtstride=1);
		int pc0wait_rcv(int phase,FLT *vdata, int vrtsride=1);
		int pc0rcv(int phase,FLT *vdata, int vrtstride=1);
		void sc0load(int phase, FLT *sdata, int bgnmode, int endmode, int modestride);
		int sc0wait_rcv(int phase, FLT *sdata, int bgnmode, int endmode, int modestride);
		int sc0rcv(int phase, FLT *sdata, int bgnmode, int endmode, int modestride);
		void tc0load(FLT *fdata, int bgnmode, int endmode, int modestride);
		int tc0wait_rcv(FLT *fdata, int bgnmode, int endmode, int modestride);
		int tc0rcv(FLT *fdata, int bgnmode, int endmode, int modestride);

		/* Some other utilities */
		void l2error(init_bdry_cndtn *toCompare);
		void findmax(int bnum, FLT (*fxy)(TinyVector<FLT,ND> &x));
		void findintercept(int bnum, FLT (*fxy)(TinyVector<FLT,ND> &x));
		//void integrated_averages(Array<FLT,1> a);
		
		virtual ~tet_hp();
};

/* THIS CLASS IS TO ALLOW SPECIAL THINGS LIKE RIGIDLY MOVING MESHES OR PARAMETER CHANGING ETC.. */
class tet_hp_helper {
	protected:
		tet_hp &x;
	public:
		tet_hp_helper(tet_hp& xin) : x(xin) {}
		tet_hp_helper(const tet_hp_helper &in_help, tet_hp& xin) : x(xin) {}
		virtual tet_hp_helper* create(tet_hp& xin) { return new tet_hp_helper(*this,xin); }
		virtual void init(input_map& input, std::string idnty) {}
		virtual void tadvance() {}
		virtual void setup_preconditioner() {}
		virtual void rsdl(int stage) {}
		virtual void update(int stage) {}
		virtual void mg_restrict() {}
		virtual void output() {};
};


#endif
