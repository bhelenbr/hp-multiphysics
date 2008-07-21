#ifndef _block_h_
#define _block_h_

#include <map>
#include <string>
#include <sstream>
#include <blitz/array.h>
#include <input_map.h>

#define DIRK 4
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

/* THIS STRUCTURE STORES ALL OF THE GLOBAL INFORMATION FOR A BLOCK */
/* COMPONENTS OF BLOCK INHERIT THIS AND THEN ALLOCATE ONE BIG GLOBAL STRUCTURE */
struct block_global {
    std::string idprefix;
    FLT dti; /**< Inverse time step */
    FLT dti_prev; /**< Inverse time step for prior step (allows changes in time step) */
    FLT time; /**< Simulation time */
    int tstep; /**< Simulation time step */
    int substep; /**< For schemes requiring multiple solves per step */
    FLT g;  /**< gravity */
    blitz::TinyVector<FLT,2> body; /**< General way for body forces */
    std::ostream *log; /**< log file stream */
    
    /** Time stepping for simulation */
#ifdef BACKDIFF
    /** @name backdiff Backwards difference constants
    *  These are constant for backwards difference timestepping
    */
    //@{
    FLT bd[BACKDIFF+1];  /**< backwards difference constants */
    //@}
#endif
#ifdef DIRK
    /** @name DIRK variables
     *  These are arrays for diagonally implicit RK Timestepping 
     */
    //@{
    FLT bd[1]; /**< Diagonal coefficient */
    FLT adirk[DIRK][DIRK]; /**< ``a'' coefficient matrix */
    FLT cdirk[DIRK]; /**< ``c'' coefficient matrix */
    bool esdirk; /**< Flag to be set when using an explicit 1'st stage */
    //@}
#endif
    /** @name Adaptation parameters 
     *  constants controlling adaptation 
     */
    //@{
    bool adapt_flag; /**< Completely turns adaptation off (can also shut off blocks individiually)  */
    bool adapt_output; /**< Flag to tell whether to give detailed adaptation data */
    FLT tolerance; /**< Tolerance for mesh adaptation scheme */
    FLT error_target; /**< Error target for mesh adaptation scheme */
    //@}
};

/* THIS IS A BLOCK THAT HAS THE CAPABILITIES OF DRIVING A MULTIGRID CYCLE OR AN EXPLICT TIME ADVANCEMENT LOOP */
class block {
    protected:
        block_global *gbl;  /**< Pointer to block globals */
        multigrid_interface* getnewlevel(input_map& blockdata);  /**< Allocates multigrid levels */
        int ntstep;  /**< Number of time steps to perform */
        int nstart; /**< Starting step (for restart from file */
        
        /** @name Output parameters
         *  constants determining when and what to output 
         */
        //@{
        int out_intrvl; /**< Number of time-steps between data outputs */
        int rstrt_intrvl; /**< Number of output intervals between restart files */
        bool debug_output; /**< Output file every iteration */
        //@}

        /** @name Multigrid parameters 
         *  constants defining multigrid iteratoin 
        */
        //@{
        Array<multigrid_interface *,1> grd;
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
            char buffer[100];
            std::string keyname;
            sprintf(buffer,"b%d",idnum);
            idprefix = std::string(buffer);
        }
        
        /** Function to start thread */
        virtual void go(input_map input);
        
        /** Initialization function */
        virtual void init(input_map& input);
        
        /** Sets-up parallel communications, called by init */
        void findmatch(int level);
        class comm_info;  /**< Utility class for figuring out communication */

        /** Outputs solution in various filetypes */
        enum output_purpose {display, restart, debug};
        void output(const std::string &filename, output_purpose why, int level = 0);

        /** Shift to next implicit time step */
        virtual void tadvance(); 
    
        /** Iterate on all blocks */  
        void iterate(int mglvl, int niter);

        /** Multigrid cycle */
        void cycle(int vw, int lvl = 0);    
        FLT maxres(int lvl = 0);
        
        /** Mesh adaptation routines */
        void adapt();
        
        virtual ~block() {}
};

class multigrid_interface {
    public:
        /** Initialization functions */
        virtual void* create_global_structure() {return 0;}
        virtual void init(input_map& input, void *gbl_in) {}
        enum init_purpose {duplicate, multigrid, adapt_storage, user_defined};
        virtual void init(const multigrid_interface& fine, init_purpose why=duplicate, FLT sizereduce1d=1.0) {}
        
        /** Outputs solution in various filetypes */
        virtual void output(const std::string &filename, block::output_purpose why) {}

        /** Shift to next implicit time step */
        virtual void tadvance() {}
                
        /** Makes sure vertex positions on boundaries coinside */
        virtual void matchboundaries() {}
                                 
        /** Setup preconditioner */
        virtual void setup_preconditioner() {}
                
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
        
        /* STUFF TO MATCH COMMUNICATION BOUNDARIES */
        virtual int comm_entity_size() {return 0;}
        virtual int comm_entity_list(blitz::Array<int,1>& list) {return 0;}
        virtual boundary* getvbdry(int num) {return 0;}
        virtual boundary* getebdry(int num) {return 0;}
        virtual boundary* getfbdry(int num) {return 0;}
        
        virtual ~multigrid_interface() {}
};
#endif
