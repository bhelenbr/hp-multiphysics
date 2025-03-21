#ifndef _block_h_
#define _block_h_

#include <map>
#include <string>
#include <sstream>
#include <blitz/array.h>
#include <input_map.h>
#include <symbolic_function.h>

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

class boundary;
class multigrid_interface;
using namespace blitz;

//#define OLDRECONNECT

/** \brief Global variables for simulation
 *
 * This namespace contains global varaibles and constants accessible to all
 */
namespace sim {
	/** @name DIRK3 scheme constants
	*  Constants used to define the DIRK3 scheme
	*/
	//@{
	const FLT GRK3 = 0.43586652150845899941601945;
	const FLT C2RK3 = (2.-9*GRK3+6.*GRK3*GRK3)/(3*(1-4*GRK3+2*GRK3*GRK3));
	const FLT B2RK3 = -3*(1-4*GRK3+2*GRK3*GRK3)*(1-4*GRK3+2*GRK3*GRK3)/(4*(-1+6*GRK3-9*GRK3*GRK3+3*GRK3*GRK3*GRK3));
	//@}

	/** @name DIRK4 scheme constants
	*  Constants used to define the DIRK4 scheme
	*/
	//@{
	// const FLT GRK4 0.5  // FOR ERROR PREDICTION
	const FLT GRK4 = 0.43586652150845899941601945; // FOR L-STABILITY
	const FLT C3RK4 = 2*GRK4*(GRK4-1./4.)*(GRK4-1.)/((GRK4-0.5)*(GRK4-0.5)-1./12.);
	const FLT A32RK4 = (C3RK4*C3RK4-2*C3RK4*GRK4)/(4.*GRK4);
	const FLT B1RK4 = (1./6. +GRK4*GRK4-GRK4*GRK4*C3RK4+3./2.*GRK4*C3RK4-GRK4-1./4.*C3RK4)/(GRK4*C3RK4);
	const FLT B2RK4 = (1./3.-GRK4-1./2.*C3RK4+GRK4*C3RK4)/(2.*GRK4*(2.*GRK4-C3RK4));
	const FLT B3RK4 =( 1./3.-2*GRK4+2*GRK4*GRK4)/(C3RK4*(C3RK4-2*GRK4));
	//@}
}

/* THIS STRUCTURE STORES ALL OF THE GLOBAL INFORMATION FOR A BLOCK */
/* COMPONENTS OF BLOCK INHERIT THIS AND THEN ALLOCATE ONE BIG GLOBAL STRUCTURE */
struct block_global {
	int idnum; /**< global block number */
	std::string idprefix; /**< "b" + block # */

	FLT dti; /**< Inverse time step */
	FLT dti_prev; /**< Inverse time step for prior step (allows changes in time step) */
	FLT time; /**< Simulation time */
	int tstep; /**< Simulation time step */
	int substep; /**< For schemes requiring multiple solves per step */
	int iteration; /**< Iterative counter within substep */
    int auto_timestep_tries; /**< number of times to time step to be reduced before failure */
    FLT auto_timestep_ratio; /**< factor to increase or decrease time step */
    FLT auto_dti_min, auto_dti_max; /**< minimum & maximum inverse time step for auto_timestep */
    FLT dti_prev_store; /**< storage for dti_prev in case time step needs to be reset */
    FLT auto_timestep_maxtime; /**< maximum time step before time stepping stops */
    int recursive_timestep_levels; /**< maximum number of refinement levels for recursive timestepping */
    int recursive_level; /**< current refinement  level being used for a timestep */
    int recursive_nsuccesses; /**< number of successes required to go up a level*/
    int recursive_successes; /**< number of consecutive successes at  current_level */
    int recursive_fraction; /**< fractional progress at current level */
    
    
	FLT g;  /**< gravity */
	blitz::TinyVector<FLT,2> body; /**< General way for body forces */
	std::ostream *log; /**< log file stream */

	/** Time stepping data for simulation */
    enum time_schemes {DIRK1=1, DIRK2, DIRK3, DIRK4, AM1, BD1, BD2, BD3} time_scheme;
	int nhist; /**< number of backwards difference steps */
	int nadapt; /**< number of solutions that require adaptation */
	int stepsolves; /**< Number of implicit solutions required per timestep */
	Array<FLT,1> bd;  /**< backwards difference or diagonal DIRK constants */
	double extrapolate; /**< Extrapolate guess for next time step (0,1) */

	/** @name DIRK variables
	*  These are arrays for diagonally implicit RK Timestepping
	*/
	//@{
	Array<FLT,2> adirk; /**< ``a'' coefficient matrix */
	Array<FLT,1> cdirk; /**< ``c'' coefficient matrix */
	bool esdirk; /**< Flag to be set when using an explicit 1'st stage */
	//@}

	/** @name Multistage Explicit Scheme
	*  Constants for multistage iteration scheme
	*/
	//@{
	int nstage; /**< Number of stages */
	Array<FLT,1> alpha;  /**< Multistage time step constants (imaginary) */
	Array<FLT,1> beta; /**< Multistage time step constants (real) */
	//@}

	/** @name Adaptation parameters
	*  constants controlling adaptation
	*/
	//@{
	int adapt_interval; /**< 0 Completely turns adaptation off otherwise adapts every N timesteps)  */
    bool adaptable; /**< turns off adaptation individually by block */
	bool adapt_output; /**< Flag to tell whether to give detailed adaptation data */
	FLT tolerance; /**< Tolerance for mesh adaptation scheme */
	FLT error_target; /**< Error target for mesh adaptation scheme */
	FLT min_length; /**< minumum allowable edge length */
	FLT max_length; /**< maximum allowable edge length */
	int length_smoothing_steps; /**< number of iterations for smoothing of target length function */
	//@}
	
	/** @name Debugging flags
	 *  constants controlling debugging options
	 */
	//@{
	int rsdl_debug; /**< Output residuals in log file and exit for debugging */
	int jac_debug; /**< Output jacobian in log file and exit for debugging */
	//@}
};

/* THIS IS A BLOCK THAT HAS THE CAPABILITIES OF DRIVING A MULTIGRID CYCLE OR AN EXPLICT TIME ADVANCEMENT LOOP */
class block {
	protected:
		shared_ptr<block_global> gbl;  /**< Pointer to block globals */
		multigrid_interface* getnewlevel(input_map& blockdata);  /**< Allocates multigrid levels */
		int ntstep;  /**< Number of time steps to perform */
		int nstart; /**< Starting step (for restart from file */

		/** @name Output parameters
		*  constants determining when and what to output
		*/
		//@{
		int out_intrvl; /**< Number of time-steps between data outputs */
		int rstrt_intrvl; /**< Number of output intervals between restart files */
		int debug_output; /**< Output file every iteration */
		//@}

		/** @name Multigrid parameters
		*  constants defining multigrid iteratoin
		*/
		//@{
		Array<multigrid_interface *,1> grd;
		friend class multigrid_interface;
		int mglvls; /**< Total number of levels of multigrid */
		int ngrid; /**< Number of grids (could be more or less than mglvls) */
		int extra_finest_levels; /**< Number of extra levels to included on finest grid */
		int extra_coarsest_levels; /**< Number of extra levels to be include on coarsest grid */
		int ncycle; /**< Number of iterations per timestep */
		FLT absolute_tolerance; /**< Absolute error tolerance for iterative loop (negative means don't use) */
		FLT relative_tolerance; /**< Relative error tolerance for iterative loop (negative means don't use) */
		/** This for running a two-level iteration on performing extra iterations on coarsest mesh **/
		int error_control_level; /**< Level of multigrid to repeat iterations until convergence (negative means don't use) */
		FLT error_control_tolerance; /**< Relative error tolerance for error_control_level */
		/** Number of cycles between re-evaluation of preconditioner.
		*  negative mean reevaluate for both refinement and coarsening sweeps
		*/
		int prcndtn_intrvl;
		int itercrsn; /**< Number of iterations to performing on coarsening */
		int iterrfne; /**< Number of iteration to perform on moving to finer mesh */
		int vw; /**< V-cycle = 1, W-cycle = 2 */
		//@}

	public:
		std::string idprefix;
		int idnum;
		block(int idin) : idnum(idin) {
			std::string keyname;
            ostringstream nstr;
            nstr << "b" << idnum << std::flush;
			idprefix = nstr.str();
		}

		/** Function to start thread */
		void go(input_map input);
        void time_step();
        void recursive_time_step(int level = 0);

		/** Initialization function */
		void init(input_map& input);

		/** Outputs solution in various filetypes */
		enum output_purpose {display, restart, debug};
		void output(const std::string &filename, output_purpose why, int level = 0);

		/** Shift to next implicit time step */
		void tadvance();
        void reset_timestep();

		/** Iterate on all blocks */
		void iterate(int mglvl, int niter);

		/** Multigrid cycle */
		FLT cycle(int vw, int lvl, bool evaluate_preconditioner=true);
		FLT maxres(int lvl = 0);

		/** Mesh adaptation routines */
		void adapt();

		virtual ~block();
};

class multigrid_interface {
	public:
		/** Initialization functions */
		virtual void init(input_map& input, shared_ptr<block_global> gbl_in) {}
		enum init_purpose {duplicate, multigrid, adapt_storage, user_defined};
		virtual void init(const multigrid_interface& fine, init_purpose why=duplicate, FLT sizereduce1d=1.0) {}

		/** Outputs solution in various filetypes */
		virtual void output(const std::string &filename, block::output_purpose why) {}

		/** Shift to next implicit time step */
		virtual void tadvance() {}
        virtual void reset_timestep() {}

		void findmatch(shared_ptr<block_global> gbl, int grdlvl); /**< Sets-up parallel communications, called by init */
		class comm_info;  /**< Utility class for figuring out communication */
		virtual void calculate_halo() {} /**< calculate halo mesh structures in partition boundaries */
		virtual void matchboundaries() {} /**< Makes sure data on boundaries coinside */

		/** Setup preconditioner */
        virtual int setup_preconditioner() {return(0);}

		/** Calculate residuals */
		virtual void rsdl() {}

		/** Relax solution */
		virtual void update() {}

		/** Multigrid cycle */
		virtual void mg_restrict() {}
		virtual void mg_prolongate() {}
		virtual void connect(multigrid_interface &fine) {}

		/** Print errors */
		virtual FLT maxres() {return(0.0);}

		/** Mesh adaptation routines */
		virtual void adapt() {}
        virtual void refineby2() {}
        virtual void offset_geometry(input_map& input) {}

        /** Mesh partition routines */
		virtual void setpartition(int nparts, blitz::Array<int,1> part_list) {}
		virtual void partition(multigrid_interface& src, int part, int maxenum, int maxvnum) {}

		/* STUFF TO MATCH COMMUNICATION BOUNDARIES */
		virtual int comm_entity_size() {return 0;}
		virtual int comm_entity_list(blitz::Array<int,1>& list) {return 0;}
		virtual boundary* getvbdry(int num) {return 0;}
		virtual boundary* getebdry(int num) {return 0;}
		virtual boundary* getfbdry(int num) {return 0;}
		virtual void match_bdry_numbering() {}

		virtual ~multigrid_interface() {}
};

/** \brief Generic interface for a b.c.
 *
 *  \ingroup boundary
 *  Contains all the basic functions for parallel communications
 */
class boundary {
	public:
		int idnum;
		std::string idprefix;
		std::string mytype;

		boundary(int idin) : idnum(idin) {
			std::string keyname;
            ostringstream nstr;
            nstr << idnum << std::flush;
            idprefix = nstr.str();
			mytype = "boundary";
		}
		virtual void init(input_map& bdrydata) {}
		virtual void alloc(int n) {}
	
		/* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
		enum msg_type {flt_msg, int_msg};
		enum groups {all,all_phased,partitions,manifolds};
		enum comm_type {symmetric,master_slave,slave_master};
		enum operation {average,sum,maximum,minimum,replace};
		union  {
			bool bdum;
			int idum;
			FLT fdum;
			msg_type mdum;
		} dummy;
		virtual bool is_comm() {return(false);}
		virtual bool is_local(int matchnum) {return(true);}
		virtual bool& is_frst() {return(dummy.bdum=true);}
		virtual int& group() {return(dummy.idum=1);}
		virtual bool in_group(int group) {return(false);}
		virtual int local_cnnct(boundary *bin, int snd_tag, int rcv_tag) {return 1;}
#ifdef MPISRC
		virtual int mpi_cnnct(int proc_tgt, int snd_tag, int rcv_tag, std::string name) {return 1;}
#endif
		virtual int& matchphase(boundary::groups group, int matchnum) {return(dummy.idum=0);}
		virtual void resize_buffers(int size) {}
		virtual void *psndbuf() {return(&dummy);}
		virtual int& nmatches() {return(dummy.idum=0);}
		virtual int& msg_phase(int grp, int match) {return(dummy.idum=0);}
		virtual int& isndbuf(int indx) {return(dummy.idum);}
		virtual FLT& fsndbuf(int indx) {return(dummy.fdum);}
		virtual int& isndbuf(int m, int indx) {return(dummy.idum);}
		virtual FLT& fsndbuf(int m, int indx) {return(dummy.fdum);}
		virtual int& ircvbuf(int m,int indx) {return(dummy.idum);}
		virtual FLT& frcvbuf(int m,int indx) {return(dummy.fdum);}
		virtual int& sndsize() {return(dummy.idum=0);}
		virtual bool& one_send_buf() {return(dummy.bdum=true);}
		virtual void match_name(int m, std::string& name) {name = std::string("");}
		virtual boundary::msg_type& sndtype() {return(dummy.mdum);}
		virtual void comm_prepare(boundary::groups group, int phase, comm_type type) {}
		virtual void comm_exchange(boundary::groups group, int phase, comm_type type) {}
		virtual int comm_wait(boundary::groups group, int phase, comm_type type) {return 1;}
		virtual int comm_nowait(boundary::groups group, int phase, comm_type type) {return 1;}
		virtual bool comm_finish(boundary::groups group, int phase, comm_type type, operation op) {return(false);}
		virtual ~boundary() {}
};


/** \brief Interface & template to make a boundary type that can be used provide geometric information
 *
 */
template<int ND> class vgeometry_interface {
	public:
		virtual void mvpttobdry(TinyVector<FLT,ND> &pt) {}
		virtual ~vgeometry_interface() {}
};

template<int ND> class egeometry_interface {
	public:
		virtual void mvpttobdry(int seg, FLT psi, TinyVector<FLT,ND> &pt) {}
		virtual void bdry_normal(int seg, FLT psi, TinyVector<FLT,ND> &norm) {}
		virtual ~egeometry_interface() {}
};

template<int ND> class fgeometry_interface {
	public:
		virtual void mvpttobdry(int elem, FLT psi, FLT eta, TinyVector<FLT,ND> &pt) {}
		virtual void mvpttobdry(int seg, FLT psi, TinyVector<FLT,ND> &pt) {}
		virtual void bdry_normal(int elem, FLT psi, FLT eta, TinyVector<FLT,ND> &norm) {}
		virtual ~fgeometry_interface() {}
};
#endif
