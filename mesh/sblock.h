#include "rblock.h"
#include "mesh.h"
#include <rbd.h>

#ifdef CAPRI
#include <capri.h>
#endif

class capriblock : public mgrid<mesh<3> > {
   private:
      /* FOR TRANSLATION & ROTATION */
      /* EASIER / FASTER TO DO LOCALLY THAN IN CAPRI */
      FLT dx[3],angles[3];
      FLT deltas[4];
      int cpri_vol, cpri_face;
      int step;
      
   public:
      capriblock() : step(0) {};
      
      void init(std::map <std::string,std::string>& input, std::ostream *inlog) {
         input["ngrid"] = "1";
         mgrid<mesh<3> >::init(input,inlog);
      }

      void load_const(std::map <std::string,std::string>& input) {
         std::istringstream data(input["Volume"]);
         data >> cpri_vol;
         *log << "#Volume: " << cpri_vol << std::endl;
         data.clear();
         
         data.str(input["Face"]);
         data >> cpri_face;
         *log << "#Face: " << cpri_face << std::endl;
         data.clear();
         
         data.str(input["dx"]);
         data >> dx[0];
         data >> dx[1];
         data >> dx[2];
         *log << "#dx: " << dx[0] << ' ' << dx[1] << ' ' << dx[2] << std::endl;
         data.clear();
         
         data.str(input["angles"]);
         data >> angles[0];
         data >> angles[1];
         data >> angles[2];
         *log << "#dx: " << angles[0] << ' ' << angles[1] << ' ' << angles[2] << std::endl;
         data.clear();
         
         data.str(input["deltas"]);
         data >> deltas[0];
         data >> deltas[1];
         data >> deltas[2];
         data >> deltas[3];
         *log << "#dx: " << deltas[0] << ' ' << deltas[1] << ' ' << deltas[2] << ' ' << deltas[3] << std::endl;
         data.clear();
         
         setrbadd(dx,angles[0],angles[1],angles[2]);
         
#ifdef CAPRI
         for(int i=0;i<ngrid;++i) {
            grd[i].cpri_vol = cpri_vol;
            grd[i].cpri_face = cpri_face;
         }
#endif
      }
      
      void alloc(std::map<std::string,std::string>& input) {
         mgrid<mesh<3> >::alloc(input);
         
         for(int i=0;i<grd[0].nvrtx;++i)
            rigidbodyadd(grd[0].vrtx[i]);
      }
         


      block::ctrl tadvance(int lvl, int excpt) {
         int i,n;
         FLT xpt1[3],xpt2[3];
#ifdef CAPRI
         double uv[2] = {-100.0, -100.0};
#endif
         
         if (lvl != 0 || excpt != 0) return(stop);
         
         ++step;
         
         setrbadd(dx,angles[0],angles[1],angles[2]);
         setrbrmv();
         
         dx[0] += deltas[0]*cos(M_PI*angles[0]/180.0);
         dx[1] += deltas[0]*sin(M_PI*angles[0]/180.0);
         dx[2] = +0.5;
         
         angles[2] = deltas[2]*sin(M_PI*(step-1)/25.0);
         setrbadd(dx,angles[0],angles[1],angles[2]);

         if (cpri_face == 7 || cpri_face == 8) {
            for(i=0;i<grd[0].nvrtx;++i) {           
               for(n=0;n<3;++n)
                  xpt1[n] = grd[0].vrtx[i][n];            
               rigidbodyrmv(xpt1);
               for(n=0;n<3;++n) 
                  xpt2[n] = xpt1[n];
               rigidbodyadd(xpt2);
               xpt1[2] -= (1.0 -MIN(grd[0].vrtx[i][2]/(-2.0),1.0))*(xpt2[2]-grd[0].vrtx[i][2]);
#ifdef CAPRI
               int icode = gi_qNearestOnFace(cpri_vol, cpri_face, xpt1, uv, xpt2);
#endif	
               rigidbodyadd(xpt2);
               for(n=0;n<3;++n)
                  grd[0].vrtx[i][n] = xpt2[n];
            }
            
            setrbrmv();

            for(i=0;i<grd[0].nsbd;++i) 
               grd[0].sbdry[i]->tadvance();
            
         	for(i=0;i<grd[0].nvbd;++i)
               grd[0].vbdry[i]->tadvance();
         }
           
         else {
            for(i=0;i<grd[0].nvrtx;++i) {           
               rigidbodyrmv(grd[0].vrtx[i]);
               rigidbodyadd(grd[0].vrtx[i]);
            }
         }
            
         return(stop);
      }
      
      int comm_entity_size(int grdlvl) {
         if (grdlvl == 0) return(grd[grdlvl].comm_entity_size());
         return(0);
      }

      int comm_entity_list(int grdlvl, int *list) {
         if (grdlvl == 0)  return(grd[grdlvl].comm_entity_list(list));
         return(0);
      }

      boundary* vbdry(int grdlvl, int num) {
         if (grdlvl == 0) return(grd[grdlvl].vbdry[num]);
         return(0);
      }

      boundary* sbdry(int grdlvl, int num) {
         if (grdlvl == 0) return(grd[grdlvl].sbdry[num]);
         return(0);
      }

      block::ctrl reconnect(int grdnum, int phase) {
         return(stop);
      }

      block::ctrl matchboundaries(int lvl, int excpt) {
         if (lvl != 0) return(stop);
         switch (excpt) {
            case(0):
               mp_phase = 0;
               grd[lvl].matchboundaries1();
               return(advance);
            case(1):
               return(static_cast<ctrl>(grd[lvl].msgpass(mp_phase++)));
            case(2):
               grd[lvl].matchboundaries2();
               return(stop);
         }
         
         *log << "control flow error matchboundaries\n";
         exit(1);
         
         return(stop);
      }
};




      

