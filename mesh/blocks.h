#include"r_mesh.h"

/*	THIS IS A SEQUENCE OF RMESHES & STATIC STORAGE FOR A SIGNLE MULTIGRID BLOCK */
class block {
   private:
      int ngrid;
      struct r_mesh_glbls rglbl;
      
   public:
      class r_mesh *grd;
      void init(int n, char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1);
      void reconnect();
      friend class blocks;
};


/*	THIS IS A MULTIBLOCK MESH */
class blocks {
   private:
      int nblocks;
      int mglvls;
      class block *blk;
      
   public:
/*		INITIALIZE MULTIBLOCK/MGRID MESH */
      void init(int nb, int mg, char **filenames, FILETYPE filetype = easymesh, FLT grwfac = 1);
      void findmatch();

/*		JACOBI ITERATION ON ALL BLOCKS */  
      void jacobi(int niter, int mglvl);

/*		MGRID CYCLE vw = 1: v-cycle vw =2 w-cycle */
      void cycle(int vw, int lvl = 0);

/*		CALCULATE DEFORMATION SOURCE TERM ON FINEST MESH */
      void ksrc();
      
/*		OUTPUT FINE BLOCK MESH */
      void out_mesh(char **filename, FILETYPE filetype = easymesh);
      void out_mesh(char *filename, FILETYPE filetype = easymesh);
      
/*		PRINT ERRORS */
      inline void maxres() {
         for(int i=0;i<nblocks;++i)
            blk[i].grd[0].maxres();
      }
      
/*		PERTURB BLOCK MESH */
      inline void perturb() {
         for(int i=0;i<nblocks;++i)
            blk[i].grd[0].perturb();
         return;
      }
      
      void restructure(FLT tolerance) {
         int i;
         
         for (i=0;i<nblocks;++i) 
            blk[i].grd[0].length1();
            
         for(i=0;i<nblocks;++i)
            blk[i].grd[0].length_mp();
            
         for(i=0;i<nblocks;++i) {
            blk[i].grd[0].length2();
            blk[i].grd[0].swap();
            blk[i].grd[0].yaber(1.0/tolerance);
            blk[i].grd[0].treeupdate();
            blk[i].grd[0].rebay(tolerance);
            blk[i].reconnect();
         }
         return;
      }
};