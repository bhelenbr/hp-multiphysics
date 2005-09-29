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
      int ngrid, mp_phase;
      GRD *grd;
      typedef typename GRD::transfer gtrans;
      typedef typename GRD::gbl ggbl;
      Array<Array<gtrans,1>,1> cv_to_ft;
      Array<Array<gtrans,1>,1> fv_to_ct;
      ggbl gstorage;
      
      bool adapt_flag;
      FLT tolerance;
   
   public:
      mgrid(int idnum) : block(idnum), adapt_flag(0) {}
      void init(std::map <std::string,std::string>& input);
      void reload_scratch_pointers() {
         for(int i=0;i<ngrid;++i) {
            grd[i].reload_scratch_pointers();
         }
      }
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
      int comm_entity_list(int grdlvl, Array<int,1>& list) {
         return(grd[grdlvl].comm_entity_list(list));
      }
      boundary* vbdry(int grdlvl, int num) {
         return(grd[grdlvl].vbdry(num));
      }
      boundary* sbdry(int grdlvl, int num) {
         return(grd[grdlvl].sbdry(num));
      }
      boundary* fbdry(int grdlvl, int num) {
         return(0);
         // return(grd[grdlvl].fbdry[num]);
      }
      
      block::ctrl tadvance(int lvl, int excpt) {
         if (lvl == 0) {
            return(grd[lvl].tadvance(0,excpt,fv_to_ct(0),cv_to_ft(0),0));
         } else {
            return(grd[lvl].tadvance(1,excpt,fv_to_ct(lvl-1),cv_to_ft(lvl-1), &grd[lvl-1]));
         }
         return(stop);
      }
      block::ctrl rsdl(int lvl, int excpt) {
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
         return(grd[lvl].mg_getfres(excpt,fv_to_ct(lvl-1),cv_to_ft(lvl-1), &grd[lvl-1]));
      }
      block::ctrl mg_getcchng(int lvl, int excpt) {
         return(grd[lvl].mg_getcchng(excpt,fv_to_ct(lvl), cv_to_ft(lvl), &grd[lvl+1]));
      }
      block::ctrl adapt(int excpt);
};

template<class GRD> void mgrid<GRD>::init(std::map <std::string,std::string>& input) {
   std::string keyword;
   std::istringstream data;
   std::string filename;
   std::string bdryfile;

   /* LOAD NUMBER OF GRIDS */
   data.str(input["ngrid"]);
   if (!(data >> ngrid)) ngrid = 1;
   *sim::log << "#ngrid: " << ngrid << std::endl;
   data.clear();
   
   keyword = idprefix + ".adapt";
   data.str(input[keyword]);
   if (!(data >> adapt_flag)) {
      data.clear();
      keyword = "adapt";
      data.str(input[keyword]);
      if (!(data >> adapt_flag)) adapt_flag = 0;
   }
   *sim::log << "#adapt: " << adapt_flag << std::endl;
   data.clear();  
   
   grd = new GRD[ngrid];
      
   cv_to_ft.resize(ngrid-1);
   fv_to_ct.resize(ngrid-1);
   
   FLT grwfac;
   keyword = idprefix + ".growth factor";
   data.str(input[keyword]);
   if (!(data >> grwfac)) {
      data.clear();
      keyword = "growth factor";
      data.str(input[keyword]);
      if (!(data >> grwfac)) grwfac = 2.0;
   }
   *sim::log << "#growth factor: " << grwfac << std::endl;
   data.clear();
   
   int filetype;
   keyword = idprefix + ".filetype";
   data.str(input[keyword]);
   if (!(data >> filetype)) {
      data.clear();
      keyword = "filetype";
      data.str(input[keyword]);
      if (!(data >> filetype)) filetype = ftype::grid;
   }
   *sim::log << "#filetype: " << filetype << std::endl;
   data.clear();
   
   keyword = idprefix + ".mesh";
   data.str(input[keyword]);
   if (!(data >> filename)) {
      data.clear();
      keyword = "mesh";
      data.str(input[keyword]);
      if (data >> filename) {
         filename = filename +"_" +idprefix;
      }
      else {
         *sim::log << "no mesh name\n"; exit(1);
      }
   }
   *sim::log << "#mesh: " << filename << std::endl;
   data.clear();   

   keyword = idprefix + ".bdryfile";
   data.str(input[keyword]);
   if (!(data >> bdryfile)) {
      bdryfile = filename +"_bdry.inpt";
      input[keyword] = bdryfile;
   }
   *sim::log << "#bdryfile: " << bdryfile << std::endl;
   grd[0].mesh::input(filename.c_str(),static_cast<ftype::name>(filetype),grwfac,bdryfile.c_str());
   data.clear();   
   
   keyword = idprefix + ".tolerance";
   data.str(input[keyword]);
   if (!(data >> tolerance)) {
      data.clear();
      keyword = "tolerance";
      data.str(input[keyword]);
      if (!(data >> tolerance)) tolerance = 2.2; 
   }
   *sim::log << "#tolerance: " << tolerance << std::endl;
   data.clear();
   
   keyword = idprefix + ".coarse";
   input[keyword] = "0";
   grd[0].init(input,idprefix,&gstorage);
   
   input[keyword] = "1";
#define OLDRECONNECT
   for(int lvl=1;lvl<ngrid;++lvl) {
#ifdef OLDRECONNECT
      grd[lvl].allocate_duplicate(2.0,grd[lvl-1]);
#else
      FLT size_reduce = 1.0;
      if (lvl > 1) size_reduce = 2.0;
      grd[lvl].allocate_duplicate(size_reduce,grd[lvl-1]);
#endif
      grd[lvl].init(input,idprefix,&gstorage);
      cv_to_ft(lvl-1).resize(grd[lvl].maxvst);
      fv_to_ct(lvl-1).resize(grd[lvl-1].maxvst);
   }
   input[keyword] = "0";
   
   return;
}

template<class GRD> block::ctrl mgrid<GRD>::reconnect(int lvl, int excpt) {
   int state = block::stop;
   FLT size_reduce = 1;
   
   switch(excpt) {
      case(0):
#ifdef OLDRECONNECT
         grd[lvl].coarsen(1.6,grd[lvl-1]);
#else
         if (lvl > 1) size_reduce = 2.0;
         grd[lvl].coarsen2(1.5,grd[lvl-1],size_reduce);
#endif
      default:
         state &= grd[lvl].mgconnect(excpt,cv_to_ft(lvl-1),grd[lvl-1]);
         state &= grd[lvl-1].mgconnect(excpt,fv_to_ct(lvl-1),grd[lvl]);
   }

   return(static_cast<block::ctrl>(state));
}

template<class GRD> void mgrid<GRD>::coarsenchk(const char *fname) {
   int i;
   char name[100];

   for(i = 1; i< ngrid; ++i) {
      number_str(name,fname,i,1);
      grd[i].checkintegrity();
      grd[i].output(name,ftype::grid);
      strcat(name,"_fv_to_ct");
      grd[i-1].testconnect(name,fv_to_ct(i-1),&grd[i]);
      number_str(name,fname,i,1);
      strcat(name,"_cv_to_ft");
      grd[i].testconnect(name,cv_to_ft(i-1),&grd[i-1]);         
   }
   
   return;
}

template<class GRD> block::ctrl mgrid<GRD>::matchboundaries(int lvl, int excpt) {
   
   switch (excpt) {
      case(0):
         mp_phase = -1;
         return(block::advance);
      case(1):
         ++mp_phase;
         /* MESSAGE PASSING SEQUENCE */
         switch(mp_phase%3) {
            case(0):
               grd[lvl].matchboundaries1(mp_phase/3);
               return(stay);
            case(1):
               grd[lvl].vmsgpass(mp_phase/3);
               return(stay);
            case(2):
               return(static_cast<ctrl>(grd[lvl].matchboundaries2(mp_phase/3)));
         }
      case(2):
         return(stop);
   }
   
   *sim::log << "control flow error matchboundaries\n";
   exit(1);
   
   return(stop);
}

template<class GRD> block::ctrl mgrid<GRD>::adapt(int excpt) {
   
   if (!adapt_flag) return(stop);
   
   switch(excpt) {
      case(0):
         mp_phase =-1;
         // grd[0].length();
         return(advance);
      case(1):
         ++mp_phase;
         /* MESSAGE PASSING SEQUENCE */
         switch(mp_phase%3) {
            case(0):
               grd[0].vmsgload(mp_phase/3,grd[0].vlngth.data(),0,0,1);
               return(stay);
            case(1):
               grd[0].vmsgpass(mp_phase/3);
               return(stay);
            case(2):
               return(static_cast<ctrl>(grd[0].vmsgwait_rcv(mp_phase/3,grd[0].vlngth.data(),0,0,1)));
         }
      default:
         return(grd[0].adapt(excpt-2,tolerance));
   }
   
   *sim::log << "control flow error: adapt" << std::endl;
   exit(1);
   
   return(stop);
}

#endif
