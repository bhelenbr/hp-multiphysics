#include "block.h"
#include <fstream>

/* THIS IS A MULTIBLOCK MESH */
/* CAN DRIVE MULTIGRID OR STANDARD ITERATION */

class blocks {
   private:
      int nproc, myid;  // MPI INFO
      int nblock, mglvls, ngrid;
      int ncycle, njacobi, itercrsn, iterrfne, vw;
      int ntstep;
      block **blk;
      std::ostream *log;
      std::ofstream filelog;

   public:
      /* INITIALIZE MULTIBLOCK/MGRID MESH */
      blocks() : nproc(1), myid(0) {}
      block* getnewblock(int idnum, std::map<std::string,std::string> *blockdata);
      void init(std::map<std::string,std::string> input);
      void init(const char *infile, const char *outfile = 0);
      void findmatch();
      void matchboundaries();
      void output(char *filename, ftype::name filetype);
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




