#include"hp_mgrid.h"

/*	THIS IS A SEQUENCE OF RMESHES & STATIC STORAGE FOR A SIGNLE MULTIGRID BLOCK */
class block {
   private:
      int ngrid;
      int lg2pmax;
      class hpbasis *hpbase;
      struct hp_mgrid_glbls gbl;
      struct surface_glbls *sgbl;
      class hp_mgrid *grd;
      friend class blocks;
      
   public:
      void initialize(char *inputfile, int grds, class hpbasis *bin, int lg2p);
      void tadvance();
      void reconnect();
};

/*	THIS IS A MULTIBLOCK MESH */
class blocks {
   private:
      int lg2pmax;  // HP PARAMETER
      class hpbasis base[MXLG2P];
      
      int nblocks;  // NUMBER OF BLOCKS
      class block *blk;

      int ntstep,out_intrvl; // TIME STEP STUFF
      int mglvls,mgrids,vwcycle,ncycle;  //MULTIGRID PARAMETERS

      FLT mxr[NV];  /* STORES MAXIMUM RESIDUAL */
      
/*		5 STAGE UPDATE ON ALL BLOCKS */  
      void nstage(int grdnum, int sm, int mgrid);

/*		SET-UP COMMUNICATION BOUNDARIES */      
      void findmatch();

   public:
/*		INITIALIZE MULTIBLOCK/MGRID MESH */
      void init(int nb, int mg, int lg2p, char *filename);
      void init(char *filename);
      
/*		START SIMULATION */
      void go();
      
/*		PREPARE FOR NEW TIMESTEP */
      void tadvance();

/*		MGRID CYCLE vw = 1: v-cycle vw =2 w-cycle (CALLED RECURSIVELY) */
      void cycle(int vw, int lvl = 0);
      
/*		PRINT MAXIMUM RESIDUALS (NV) */
      inline void print_maxres() {
         for(int n=0;n<NV;++n)
            printf(" %e ",mxr[n]);
      }

/*		OUTPUT FINE MESH SOLUTION */
      void output(char *filename, FILETYPE filetype = text);
      void output(int number, FILETYPE filetype = text);

/*		MESH REFINEMENT */
};
