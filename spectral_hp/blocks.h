#include"hp_mgrid.h"

/*	THIS IS A SEQUENCE OF RMESHES & STATIC STORAGE FOR A SIGNLE MULTIGRID BLOCK */
class block {
   private:
      int ngrid;
      int lg2pmax;
      class hpbasis *hpbase;
      struct hp_mgrid_glbls gbl;
      
   public:
      class hp_mgrid *grd;
      void meshinit(int n, char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1);
      inline void setphysics(FLT rho, FLT mu, FLT sigma, FLT (*f)(int n, FLT x, FLT y)) {
         gbl.rho = rho;
         gbl.rhoi = 1.0/rho;
         gbl.mu = mu;
         gbl.nu = mu/rho;
         gbl.sigma = sigma;
         gbl.func = f;
      }
      inline void setiter(FLT fadd, FLT flowcfl[], FLT adis, int chrctr) {
         gbl.fadd = fadd;
         for (int i=0;i<MXLG2P;++i)
            gbl.flowcfl[i] = flowcfl[i];
         gbl.adis = adis;
         gbl.charyes = chrctr;
      }
      void hpinit(class hpbasis *bin, int lgpmax);
      void reconnect();
};

/*	THIS IS A MULTIBLOCK MESH */
class blocks {
   private:
      int lg2pmax;
      class hpbasis base[MXLG2P];
      int nblocks;
      int mgrids;
      class block *blk;
      
   public:
/*		INITIALIZE MULTIBLOCK/MGRID MESH */
      void init(int nb, int mg, int lg2p, char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1);
      void findmatch();

/*		5 STAGE UPDATE ON ALL BLOCKS */  
      void nstage(int mglvl);

/*		MGRID CYCLE vw = 1: v-cycle vw =2 w-cycle (CALLED RECURSIVELY) */
      void cycle(int vw, int lvl = 0);

/*		OUTPUT FINE MESH SOLUTION */
      void output(char *filename, FILETYPE filetype = easymesh);

/*		MESH REFINEMENT */      
      inline void restructure(FLT tolerance) {
         for (int i=0;i<nblocks;++i) {
            blk[i].grd[0].swap();
            blk[i].grd[0].density(0.1,0.01,10.0);
            blk[i].grd[0].yaber(1.0/tolerance);
            blk[i].grd[0].rebay(tolerance);
            blk[i].reconnect();
         }
         return;
      }
};
