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
         gbl.g = 1.0;
         gbl.func = f;
      }
      inline void setiter(FLT fadd, FLT flowcfl[], FLT adis, int chrctr) {
         gbl.fadd = fadd;
         for (int i=0;i<MXLG2P;++i)
            gbl.flowcfl[i] = flowcfl[i];
         gbl.adis = adis;
         gbl.charyes = chrctr;
      }
      inline void setsurfphysics(int snum, FLT drho, FLT rhoav, FLT muav) {
         gbl.sgbl[snum].drho = drho;
         gbl.sgbl[snum].rhoav = rhoav;
         gbl.sgbl[snum].muav = muav;
      }
      inline void setsurfiter(int snum, FLT fadd[ND], FLT cfl[MXLG2P][ND]) {
         for(int n=0;n<ND;++n)
            gbl.sgbl[snum].fadd[n] = fadd[n];
         for (int i=0;i<MXLG2P;++i)
            for (int n=0;n<ND;++n)
               gbl.sgbl[snum].cfl[i][n] = cfl[n][i];
      }
      
      
      void hpinit(class hpbasis *bin, int lgpmax);
      void reconnect();
};

/*	THIS IS A MULTIBLOCK MESH */
class blocks {
   private:
      int lg2pmax;
      int nblocks;
      int mgrids;
      int mglvls;
      class block *blk;
      class hpbasis base[MXLG2P];
      
/*		THINGS FOR ALL BLOCKS */
      FLT time, dt;
      FLT mxr[NV];  /* STORES MAXIMUM RESIDUAL */
      
/*		5 STAGE UPDATE ON ALL BLOCKS */  
      void nstage(int grdnum, int sm, int mgrid);

/*		SET-UP COMMUNICATION BOUNDARIES */      
      void findmatch();


   public:
/*		INITIALIZE MULTIBLOCK/MGRID MESH */
      void init(int nb, int mg, int lg2p, char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1);
      
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
