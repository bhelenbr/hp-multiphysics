#include "rblock.h"
#include "mesh.h"

#ifdef CAPRI
#include <capri.h>
#endif

class capriblock : public mgrid<mesh<3> > {
   private:
      /* FOR TRANSLATION & ROTATION */
      /* EASIER / FASTER TO DO LOCALLY THAN IN CAPRI */
      FLT x,y,z,theta,phi;
      int cpri_vol, cpri_face;
      
   public:
      capriblock() : x(0.0), y(0.0), z(0.0),theta(0.0) {};
      
      void load_const(std::map <std::string,std::string>& input) {
         std::istringstream data(input["Volume"]);
         data >> cpri_vol;
         *log << "#Volume: " << cpri_vol << std::endl;
         data.clear();
         
         data.str(input["Face"]);
         data >> cpri_face;
         *log << "#Face: " << cpri_face << std::endl;
         data.clear();
         
#ifdef CAPRI
         for(int i=0;i<ngrid;++i) {
            grd[i].cpri_vol = cpri_vol;
            grd[i].cpri_face = cpri_face;
         }
#endif
      }

      block::ctrl tadvance(int lvl, int excpt) {
         int i,n;
         double pt1[3], pt2[3], uv[2]; 
         

         
         x += 0.05;
         theta += M_PI/20.0;
         z = 0.1*sin(theta);
         
         for(int i=0;i<grd[0].nvrtx;++i)
            grd[0].vrtx[i][0] += 0.05;
            
         for(i=0;i<grd[0].nvrtx;++i) {
            for(n=0;n<3;++n)
               pt1[n] = grd[0].vrtx[i][n];
            pt1[0] -= x;
            pt1[2] -= z;
#ifdef CAPRI
            uv[0] = -100.0;
            uv[1] = -100.0;
            int icode = gi_qNearestOnFace(cpri_vol, cpri_face, pt1, uv, pt2);
            
#endif
            for(n=0;n<3;++n)
               grd[0].vrtx[i][n] = pt2[n];
            grd[0].vrtx[i][0] += x;
            grd[0].vrtx[i][2] += z;
            
         }

         for(int i=0;i<grd[0].nsbd;++i) 
            grd[0].sbdry[i]->tadvance();
            
         for(int i=0;i<grd[0].nvbd;++i)
            grd[0].vbdry[i]->tadvance();
            
         return(stop);
      }
};




      

