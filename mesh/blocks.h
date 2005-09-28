#include "block.h"
#include <fstream>
#include <blitz/array.h>
#include <utilities.h>

/* THIS IS A MULTIBLOCK MESH */
/* CAN DRIVE MULTIGRID OR STANDARD ITERATION */

#ifndef _blocks_h_
#define _blocks_h_

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

/** \brief Global variables for simulation
 *
 * This namespace contains global variables needed by all blocks of the simulation 
 */
namespace sim {
   extern FLT dti; /**< Inverse time step */
   extern FLT time; /**< Simulation time */
   extern FLT tstep; /**< Simulation time step */
   extern FLT g;  /**< gravity */
   extern std::ostream *log; /**< log file stream */
   extern sharedmem scratch; /**< Shared work memory for all blocks */
   
   /** Time stepping for simulation */
#ifdef BACKDIFF
   /** @name backdiff Backwards difference constants
   *  These are constant for backwards difference timestepping
   */
   //@{
   extern FLT bd[BACKDIFF+1];  /**< backwards difference constants */
   const int nhist = BACKDIFF; /**< number of backwards difference steps */
   //@}
#endif
#ifdef DIRK
   /** @name DIRK variables
    *  These are arrays for diagonally implicit RK Timestepping 
    */
   //@{
   extern FLT bd[1]; /**< Diagonal coefficient */
   extern FLT adirk[DIRK][DIRK]; /**< ``a'' coefficient matrix */
   extern FLT cdirk[DIRK]; /**< ``c'' coefficient matrix */
   //@}
   
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
   const int nhist = 1; 
   //@}
#endif

   /** @name Multistage Explicit Scheme
    *  Constants for multistage iteration scheme 
    */
   //@{
   const int NSTAGE = 5; /**< Number of stages */
   extern blitz::Array<FLT,1> cfl;  /**< Global array of cfl numbers (can use how you see fit) */
   const FLT alpha[NSTAGE+1] = {0.25, 1./6., .375, .5, 1.0, 1.0}; /**< Multistage time step constants (imaginary) */
   const FLT beta[NSTAGE+1] = {1.0, 0.0, 5./9., 0.0, 4./9., 1.0}; /**< Multistage time step constants (real) */
   //@}

#ifdef PV3
   /** Flag for pV3 to tell when to reload meshes */
   extern int pv3_mesh_changed;
#endif
   
   /** @name Multigrid Parameters 
    *  Constants controlling multigrid iteration 
    */
   //@{
   extern FLT fadd;   //!< Multiplicative factor for multigrid (default = 1.0)
   //@}

}

class blocks {
   protected:
      int nproc; /**< Number of processors in simulation */
      int myid;  /**< Processor number of myself for MPI */
      int nblock; /**< Number of blocks on this processor */

      /** @name Multigrid parameters 
       *  constants defining multigrid iteratoin 
      */
      //@{
      int mglvls; /**< Number of levels of multigrid */
      int ngrid; /**< Number of grids (could be more or less than mglvls) */
      int ncycle; /**< Number of iterations per timestep */
      int njacobi; /**< Number of iterations comprising a single smoothing */
      int itercrsn; /**< Number of iterations to performing on coarsening */
      int iterrfne; /**< Number of iteration to perform on moving to finer mesh */
      int vw; /**< V-cycle = 1, W-cycle = 2 */
      //@}
      int ntstep;  /**< Number of time steps to perform */
      block **blk; /**< Array containing pointers to blocks */
      
      /** This is a data structure for storing communication info
       *  only used temporarily to sort things out in findmatch
       */
      class commblocks {
         public:
            int nblock;
            struct commblock {
               int nvcomm;
               struct vid {
                  int nvbd, idnum;
               } *vcomm;
               
               int nscomm;
               struct sid {
                  int nsbd, idnum; // , vid0, vid1;
               } *scomm;
               
               int nfcomm;
               struct fid {
                  int nfbd, idnum; // , nsds;
                  // int *sids;
               } *fcomm;
            } *blkinfo;
            
            commblocks() : nblock(0) {}
            void output();             
            int unpack(blitz::Array<int,1> entitylist); 
            ~commblocks() {
               for(int i=0;i<nblock;++i) {
                  delete []blkinfo[i].vcomm;
                  delete []blkinfo[i].scomm;
                  // for(k=0;k<blkinfo[i].nfcomm;++k)
                  //    delete []blkinfo[i].fcomm[k].sids;
                  delete []blkinfo[i].fcomm;
               }
               delete []blkinfo;
            }
      };

   public:
      /** Initialize multiblock/mgrid mesh */
      blocks() : nproc(1), myid(0) {}
      
      /** Initializes blocks using data from map */
      void init(std::map<std::string,std::string> input);
      
      /** Loads map from file then initialize */
      void init(const char *infile, const char *outfile = 0);
      
      /** Makes sure vertex positions on boundaries coinside */
      void matchboundaries();
      
      /** Outputs solution in various filetypes */
      void output(char *filename, ftype::name filetype);
      
      /** Inputs solution */
      void input(char *filename) {}
      
      /** Routine that starts simulation */
      void go();
      
      /** Adapt mesh */
      void restructure(); 
      
   protected:
      /** Allocates blocks, called by init to generate blocks from initialization file */
      block* getnewblock(int idnum, std::map<std::string,std::string> *blockdata);
            
      /** Sets-up parallel communications, called by init */
      void findmatch();
 
      /** Calculate residuals */
      void rsdl(int);

      /** Iterate on all blocks */  
      void iterate(int mglvl);

      /** Multigrid cycle */
      void cycle(int vw, int lvl = 0);
      
      /** Shift to next implicit time step */
      void tadvance(int step); 
            
      /** Print errors */
      inline void maxres() {
         for(int i=0;i<nblock;++i)
            blk[i]->maxres();
      }
};
#endif




