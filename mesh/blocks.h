#include"block.h"

/* THIS IS A MULTIBLOCK MESH */
class blocks {
   private:
      int nblock, mglvls, ngrid;
      int ncycle, njacobi, itercrsn, iterrfne;
      int ntstep;
      block **blk;
      
   public:
      /* INITIALIZE MULTIBLOCK/MGRID MESH */
      virtual block * getnewblock(int type);
      void load_constants(std::map<std::string,std::string> input[]);
      void init(std::map<std::string,std::string> input[]);
      void findmatch();
      void matchboundaries();
      void output(char *filename, FILETYPE filetype);
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




