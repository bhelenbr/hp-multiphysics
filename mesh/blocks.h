#include "block.h"
#include <fstream>

/* THIS IS A MULTIBLOCK MESH */
class blocks {
   private:
      int nproc, myid;  // MPI INFO
      int nblock, mglvls, ngrid;
      int ncycle, njacobi, itercrsn, iterrfne, vw;
      int ntstep;
      block **blk;
      std::ostream *log;
      std::ofstream filelog;
      
      /* private utilities */
      void localmatch(int vsf,boundary *v1, boundary *v2,int b1,int b2,int i,int j);
#ifdef MPISRC
      void mpimatch(int vsf,boundary *v1,int p,int b1,int b2,int i,int j);
#endif

   public:
      /* INITIALIZE MULTIBLOCK/MGRID MESH */
      blocks() : nproc(1), myid(0) {}
      block * getnewblock(int type);  // THIS MUST BE INSTANTIATED FOR PARTICULAR PROBLEM
      void init(std::map<std::string,std::string> input[]);
      void init(const char *infile, const char *outfile = 0);
      void findmatch();
      void matchboundaries();
      void output(char *filename, FTYPE filetype);
      void input(char *filename) {}
      void go();
      
      /* CALCULATE RESIDUALS */
      void rsdl(int);

      /* ITERATION ON ALL BLOCKS */  
      void iterate(int mglvl);

      /* MGRID CYCLE vw = 1: v-cycle vw =2 w-cycle */
      void cycle(int vw, int lvl = 0);
      
      /* MOVE BOUNDARIES */
      void tadvance(); 
      
      /* ADAPT MESH */
      void restructure(); 
      
      /* PRINT ERRORS */
      inline void maxres() {
         for(int i=0;i<nblock;++i)
            blk[i]->maxres();
      }
};




