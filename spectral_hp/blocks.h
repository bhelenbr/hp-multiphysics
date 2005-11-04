#include"hp_mgrid.h"
#include <utilities.h>

/* THIS IS A SEQUENCE OF RMESHES & STATIC STORAGE FOR A SIGNLE MULTIGRID BLOCK */
class block {
   private:
      int ngrid;
      int lg2pmax;
      class hpbasis *hpbase;
      struct hp_mgrid_glbls gbl;
      struct surface_glbls *sgbl;
      struct r_mesh_glbls rgbl;
      class hp_mgrid *grd;
      static class hp_mgrid temp_hp; // SOLUTION STORAGE FOR ADAPTION
      friend class blocks;
      
   public:
      block() : sgbl(0) {}
      void initialize(char *inputfile, int grds, class hpbasis *bin, int lg2p);
      void tadvance(int stage = 0);
      void reconnect();
      void coarsenchk(const char *name);
};

/* THIS IS A MULTIBLOCK MESH */
class blocks {
   private:
      int lg2pmax;  // HP PARAMETER
      class hpbasis base[MXLG2P];
      
      int nblocks;  // NUMBER OF BLOCKS
      class block *blk;

      int ntstep,out_intrvl,rstrt_intrvl,readin; // TIME STEP STUFF
      int mglvls,mgrids,vwcycle,ncycle,nup,ndown;  //MULTIGRID PARAMETERS
      int adapt; //FLAG FOR ADAPTATION

      FLT mxr[NV];  /* STORES MAXIMUM RESIDUAL */
      
      /*5 STAGE UPDATE ON ALL BLOCKS */  
      void nstage(int grdnum, int sm, int mgrid);
      
      /* SETUP COMMUNICATION BOUNDARIES */
      void findmatch(int mglvl);

   public:
      /* INITIALIZE MULTIBLOCK/MGRID MESH */
      void init(char *filename, int start_sim = 1);
      void minvrt_test(int, double (*)(int, double, double));
      
      /* START SIMULATION */
      void go();
      
      /* PREPARE FOR NEW TIMESTEP */
      void tadvance(int stage = 0);

      /* MGRID CYCLE vw = 1: v-cycle vw =2 w-cycle (CALLED RECURSIVELY) */
      void cycle(int vw, int lvl = 0);
      FLT printerror(); // PRINTS DIAGNOSTIC INFO

      /* OUTPUT FINE MESH SOLUTION */
      void output(int number, FILETYPE filetype = text, int verbose = 0);
      void loadbasis() {
         for (int i=0;i<nblocks;++i)
            blk[i].grd[0].loadbasis(&base[lg2pmax]);
      }
            

      /* MESH REFINEMENT */
      void adaptation();

      /* R-MESH DEFORMATION */
      /* JACOBI ITERATION ON ALL BLOCKS */  
      void r_jacobi(int niter, int mglvl);
      /* MGRID CYCLE vw = 1: v-cycle vw =2 w-cycle */
      void r_cycle(int vw, int lvl = 0);
      /* CALCULATE DEFORMATION SOURCE TERM ON FINEST MESH */
      void r_ksrc();
      /* PRINT MAGNITUDE OF RESIDUAL */
      /* PRINT ERRORS */
      FLT r_printerror() {
         int n;
         FLT mxr[ND],biggest = 0.0;
         
         for(int i=0;i<nblocks;++i) {
            blk[i].grd[0].r_mesh::maxres(mxr);
            for(n=0;n<ND;++n) {
               printf("%.3e  ",mxr[n]);
               biggest = MAX(mxr[n],biggest);
            }
         }
         return(biggest);
      }
      
#ifdef PV3
      void viz_init(int iopt);
      void pvscal(int *key, float *v);
      void pvstruc(int& knode, int& kequiv, int& kcel1, int& kcel2, int& kcel3, int& kcel4, int& knptet, int &kptet,int& knblock,int &blocks,int &kphedra, int& ksurf,int& knsurf,int& hint);
      void pvcell(int cel1[][4], int cel2[][5], int cel3[][6], int cel4[][8], int nptet[][8], int ptet[]);
      void pvgrid(float (*xyz)[3]);
      void pvsurface(int nsurf[][3], int scon[], int scel[][4], char tsurf[][20]);
      void pvvect(int *key, float v[][3]);
#endif
};
