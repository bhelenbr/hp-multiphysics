#include "rblock.h"
#include "mesh.h"

class sblock : public mgrid<mesh<3> > {
      virtual FLT hgt(FLT x[3]) {return(0.0);}
      virtual FLT dhgt(int dir, FLT x[3]) {return(1.0);}
      block::ctrl tadvance(int lvl, int excpt) {
         for(int i=0;i<grd[0].nvrtx;++i)
            grd[0].vrtx[i][0] += 0.1;

         for(int i=0;i<grd[0].nsbd;++i) 
            grd[0].sbdry[i]->tadvance();
            
         for(int i=0;i<grd[0].nvbd;++i)
            grd[0].vbdry[i]->tadvance();
            
         return(stop);
      }
      block::ctrl adapt(int excpt) {return(stop);}
};




      

