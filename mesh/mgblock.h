#include "r_mesh.h"
#include "block.h"
#include "boundary.h"

#include <map>
#include <string>
#include <sstream>
#include <iostream>
#include <ftype.h>

#ifndef _r_block_h_
#define _r_block_h_


/* GENERIC MULTIGRID BLOCK */

template<class GRD> class mgrid : public block {
   protected:
      std::ostream *log;
      int ngrid, mp_phase;
      GRD *grd;
      typedef typename GRD::transfer gtrans;
      typedef typename GRD::gbl ggbl;
      gtrans **cv_to_ft;
      gtrans **fv_to_ct;
      ggbl gstorage;
      
      bool adapt_flag;
      FLT tolerance;
   
   public:
      mgrid() : log(&std::cout), adapt_flag(0) {}
      void init(std::map <std::string,std::string>& input, std::ostream *inlog = 0);
      void output(char *filename, ftype::name filetype) {
         grd[0].output(filename, filetype);
      }
      void coarsenchk(const char *fname);
      block::ctrl reconnect(int lvl, int excpt);
      block::ctrl matchboundaries(int lvl, int excpt);
     
      /* PASS THROUGH ROUTINES FOR MATCHING BOUNDARIES*/
      int comm_entity_size(int grdlvl) {
         return(grd[grdlvl].comm_entity_size());
      }
      int comm_entity_list(int grdlvl, int *list) {
         return(grd[grdlvl].comm_entity_list(list));
      }
      boundary* vbdry(int grdlvl, int num) {
         return(grd[grdlvl].vbdry[num]);
      }
      boundary* sbdry(int grdlvl, int num) {
         return(grd[grdlvl].sbdry[num]);
      }
      boundary* fbdry(int grdlvl, int num) {
         return(0);
         // return(grd[grdlvl].fbdry[num]);
      }
      
      block::ctrl tadvance(int lvl, int excpt) {
         if (lvl == 0) {
            return(grd[lvl].tadvance(0,excpt,0,0,0));
         } else {
            return(grd[lvl].tadvance(1,excpt,fv_to_ct[lvl-1],cv_to_ft[lvl-1], &grd[lvl-1]));
         }
         return(stop);
      }
      block::ctrl mgrid<GRD>::rsdl(int lvl, int excpt) {
         return(grd[lvl].rsdl(excpt));
      }
      void maxres() {grd[0].maxres();}
      block::ctrl setup_preconditioner(int lvl, int excpt) {
         return(grd[lvl].setup_preconditioner(excpt));
      }
      block::ctrl update(int lvl, int excpt) {
         return(grd[lvl].update(excpt));
      }
      block::ctrl mg_getfres(int lvl, int excpt) {
         return(grd[lvl].mg_getfres(excpt,fv_to_ct[lvl-1],cv_to_ft[lvl-1], &grd[lvl-1]));
      }
      block::ctrl mg_getcchng(int lvl, int excpt) {
         return(grd[lvl].mg_getcchng(excpt,fv_to_ct[lvl], cv_to_ft[lvl], &grd[lvl+1]));
      }
      block::ctrl adapt(int excpt);
};

template<class GRD> void mgrid<GRD>::init(std::map <std::string,std::string>& input, std::ostream *inlog) {

   if (inlog) {
      log = inlog;
   } 

   /* LOAD NUMBER OF GRIDS */
   std::istringstream data(input["ngrid"]);
   if (!(data >> ngrid)) ngrid = 1;
   *log << "#ngrid: " << ngrid << std::endl;
   data.clear();
   
   data.str(input["adapt"]);
   if (!(data >> adapt_flag)) adapt_flag = 0;
   *log << "#adapt: " << adapt_flag << std::endl;
   data.clear();  
   
   grd = new GRD[ngrid];
   for(int i = 0; i<ngrid;++i)
      grd[i].log = log;
      
   cv_to_ft = new gtrans *[ngrid-1];
   fv_to_ct = new gtrans *[ngrid-1];
   
   FLT grwfac;
   data.str(input["growth factor"]);
   if (!(data >> grwfac)) grwfac = 2.0;
   *log << "#growth factor: " << grwfac << std::endl;
   data.clear();
	
	int filetype;
   data.str(input["filetype"]);
   if (!(data >> filetype)) filetype = ftype::grid;
   *log << "#filetype: " << filetype << std::endl;
   data.clear();
	
   std::string filename;
   data.str(input["mesh"]);
   if (!(data >> filename)) {*log << "no mesh name\n"; exit(1);}
   *log << "#mesh: " << filename << std::endl;
   data.clear();   

   std::string bdryfile;
   data.str(input["bdryfile"]);
   if (!(data >> bdryfile)) {
      bdryfile = filename +".bdry.inpt";
      input["bdryfile"] = bdryfile;
   }
   *log << "#bdryfile: " << bdryfile << std::endl;
   grd[0].mesh::input(filename.c_str(),static_cast<ftype::name>(filetype),bdryfile.c_str(),grwfac);
   data.clear();   
   
   data.str(input["tolerance"]);
   if (!(data >> tolerance)) tolerance = 0.6; 
   *log << "#tolerance: " << tolerance << std::endl;
   data.clear();
   
   grd[0].init(0,input,&gstorage);
       
#define OLDRECONNECT
   /* CREATE COARSE MESHES */
   for(int i=1;i<ngrid;++i) {
#ifdef OLDRECONNECT
      grd[i].coarsen(1.6,grd[i-1]);
#else
      FLT size_reduce = 1.0;
      if (i > 1) size_reduce = 0.25;
      grd[i].coarsen2(1.5,grd[i-1],0.25);
#endif
      grd[i].init(1,input,&gstorage);
      cv_to_ft[i-1] = new gtrans[grd[i].maxvst];
      grd[i].mgconnect(cv_to_ft[i-1],grd[i-1]);
      fv_to_ct[i-1] = new gtrans[grd[i-1].maxvst];
      grd[i-1].mgconnect(fv_to_ct[i-1],grd[i]);
   }
   
   return;
}


template<class GRD> block::ctrl mgrid<GRD>::reconnect(int grdnum, int phase) {
   
#ifdef OLDRECONNECT
   grd[grdnum].coarsen(1.6,grd[grdnum-1]);
#else
   FLT size_reduce = 1.0;
   if (i > 1) size_reduce = 0.25;
   grd[grdnum].coarsen2(1.5,grd[grdnum-1],0.25);
#endif
   grd[grdnum].mgconnect(cv_to_ft[grdnum-1],grd[grdnum-1]);
   grd[grdnum-1].mgconnect(fv_to_ct[grdnum-1],grd[grdnum]);

   return(stop);
}

template<class GRD> void mgrid<GRD>::coarsenchk(const char *fname) {
   int i;
   char name[100];

   for(i = 1; i< ngrid; ++i) {
      number_str(name,fname,i,1);
      grd[i].setbcinfo();
      grd[i].checkintegrity();
      grd[i].output(name);
      grd[i].setbcinfo();
   }
   return;
}

template<class GRD> block::ctrl mgrid<GRD>::matchboundaries(int lvl, int excpt) {
   
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

template<class GRD> block::ctrl mgrid<GRD>::adapt(int excpt) {

   if (!adapt_flag) return(stop);
   
   switch(excpt) {
      case(0):
         mp_phase = 0;
         grd[0].length1();
         return(advance);
      case(1):
         return(static_cast<ctrl>(grd[0].msgpass(mp_phase++)));
      case(2):
         grd[0].length2();
         return(advance);
      case(3):
         grd[0].adapt(tolerance);
         return(stop);
   }
   
   *log << "control flow error: adapt" << std::endl;
   exit(1);
   
   return(stop);
}

#endif
