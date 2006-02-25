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
      typedef typename GRD::filetype outtype;
      typedef typename GRD::transfer gtrans;
      typedef typename GRD::gbl ggbl;
      Array<Array<gtrans,1>,1> cv_to_ft;
      Array<Array<gtrans,1>,1> fv_to_ct;
   public:
      ggbl gstorage;
      bool adapt_flag;
      FLT tolerance;
   private:
      int excpt, excpt1;
   
   public:
      mgrid(int idnum) : block(idnum), adapt_flag(0) {}
      void init(input_map& input);
      void reload_scratch_pointers() {
         for(int i=0;i<ngrid;++i) {
            grd[i].reload_scratch_pointers();
         }
      }
      void input(const std::string &filename) {
         std::string fapp;
         fapp = idprefix +"_" +filename;
         grd[0].input(fapp);
      }
      void output(const std::string &filename, block::output_purpose why, int level = 0) {
         std::string fapp;
         fapp = idprefix +"_" +filename;
         grd[level].output(fapp,why);
      }
      block::ctrl reconnect(int lvl, block::ctrl ctrl_message);
      block::ctrl matchboundaries(int lvl, block::ctrl ctrl_message);
     
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
      
      block::ctrl tadvance(int lvl, block::ctrl ctrl_message) {
         if (lvl == 0) {
            return(grd[lvl].tadvance(0,ctrl_message,fv_to_ct(0),cv_to_ft(0),0));
         } else {
            return(grd[lvl].tadvance(1,ctrl_message,fv_to_ct(lvl-1),cv_to_ft(lvl-1), &grd[lvl-1]));
         }
         return(stop);
      }
      block::ctrl rsdl(int lvl, block::ctrl ctrl_message) {
         return(grd[lvl].rsdl(ctrl_message));
      }
      void maxres() {grd[0].maxres();}
      block::ctrl setup_preconditioner(int lvl, block::ctrl ctrl_message) {
         return(grd[lvl].setup_preconditioner(ctrl_message));
      }
      block::ctrl update(int lvl, block::ctrl ctrl_message) {
         return(grd[lvl].update(ctrl_message));
      }
      block::ctrl mg_getfres(int lvl, int lvlm, block::ctrl ctrl_message) {
         return(grd[lvl].mg_getfres(ctrl_message,fv_to_ct(lvlm),cv_to_ft(lvlm),&grd[lvlm]));
      }
      block::ctrl mg_getcchng(int lvl, int lvlp, block::ctrl ctrl_message) {
         return(grd[lvl].mg_getcchng(ctrl_message,fv_to_ct(lvl), cv_to_ft(lvl), &grd[lvlp]));
      }
      block::ctrl adapt(block::ctrl ctrl_message);
};

template<class GRD> void mgrid<GRD>::init(input_map& input) {
   std::string keyword;
   std::istringstream data;
   std::string filename;
   std::string bdryfile;

   /* LOAD NUMBER OF GRIDS */
   input.getwdefault("ngrid",ngrid,1);
   *sim::log << "#ngrid: " << ngrid << std::endl;
   
   keyword = idprefix + ".adapt";
   if (!input.get(keyword,adapt_flag)) {
      input.getwdefault("adapt",adapt_flag,false);
   }
   *sim::log << "#adapt: " << adapt_flag << std::endl;
   
   grd = new GRD[ngrid];
      
   cv_to_ft.resize(ngrid);
   fv_to_ct.resize(ngrid);
   
   FLT grwfac;
   keyword = idprefix + ".growth factor";
   if (!input.get(keyword,grwfac)) {
      input.getwdefault("growth factor",grwfac,2.0);
   }
   *sim::log << "#growth factor: " << grwfac << std::endl;
   
   int filetype;
   keyword = idprefix + ".filetype";
   if (!input.get(keyword,filetype)) {
      input.getwdefault("filetype",filetype,static_cast<int>(mesh::grid));
   }
   *sim::log << "#filetype: " << filetype << std::endl;
   
   keyword = idprefix + ".mesh";
   if (!input.get(keyword,filename)) {
      if (input.get("mesh",filename)) {
         filename = filename +"_" +idprefix;
      }
      else {
         *sim::log << "no mesh name\n";
         exit(1);
      }
   }
   *sim::log << "#mesh: " << filename << std::endl;

   keyword = idprefix + ".bdryfile";
   if (!input.get(keyword,bdryfile)) {
      bdryfile = filename +"_bdry.inpt";
      input[keyword] = bdryfile;
   }
   *sim::log << "#bdryfile: " << bdryfile << std::endl;
   grd[0].mesh::input(filename.c_str(),static_cast<mesh::filetype>(filetype),grwfac,bdryfile.c_str());
   
   keyword = idprefix + ".tolerance";
   if (!input.get(keyword,tolerance)) {
      input.getwdefault("tolerance",tolerance,2.2);
   }
   *sim::log << "#tolerance: " << tolerance << std::endl;
   
   keyword = idprefix + ".coarse";
   input[keyword] = "0";
   grd[0].idprefix = idprefix;
   grd[0].init(input,&gstorage);
   
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
      grd[lvl].idprefix = idprefix;
      grd[lvl].init(input,&gstorage);
      cv_to_ft(lvl-1).resize(grd[lvl].maxvst);
      fv_to_ct(lvl-1).resize(grd[lvl-1].maxvst);
   }
   input[keyword] = "0";
   
   return;
}


template<class GRD> block::ctrl mgrid<GRD>::reconnect(int lvl, block::ctrl ctrl_message) {
   int state;
   std::string name,fname;
   std::ostringstream nstr;
   
   if (ctrl_message == block::begin) excpt1 = 0;
   
   name = idprefix + "_coarsen";

   switch(excpt1) {
      case(0):
#ifdef OLDRECONNECT
         grd[lvl].coarsen(1.6,grd[lvl-1]);
#else
         FLT size_reduce = 1.0;
         if (lvl > 1) size_reduce = 2.0;
         grd[lvl].coarsen2(1.5,grd[lvl-1],size_reduce);
#endif
         ++excpt1;
         return(block::advance);
      case(1):
         state = grd[lvl].mgconnect(block::begin,cv_to_ft(lvl-1),grd[lvl-1]);
         ++excpt1;
         if (state != block::stop) return(state);
         else return(block::advance2);
      case(2):
         if (ctrl_message != block::advance2) {
            state = grd[lvl].mgconnect(ctrl_message,cv_to_ft(lvl-1),grd[lvl-1]);
            if (state != block::stop) return(state);
            else return(block::advance2);
         }
         else ++excpt1;
      case(3):
         state = grd[lvl-1].mgconnect(block::begin,fv_to_ct(lvl-1),grd[lvl]);
         ++excpt1;
         if (state != block::stop) return(state);
         else return(block::advance2);
      case(4):
         if (ctrl_message != block::advance2) {
            state = grd[lvl-1].mgconnect(ctrl_message,fv_to_ct(lvl-1),grd[lvl]);
            if (state != block::stop) return(state);
            else return(block::advance2);
         }
         else {
            /* THIS IS FOR DIAGNOSIS OF MULTI-GRID FAILURES */
            grd[lvl].checkintegrity();
            nstr << lvl << flush;
            fname = name +nstr.str();
            grd[lvl].mesh::output(fname,mesh::grid);
            fname = name +nstr.str() + "_ft_to_cv";
            grd[lvl-1].testconnect(fname,fv_to_ct(lvl-1),&grd[lvl]);
            fname = name +nstr.str() + "_cv_to_ft";
            grd[lvl].testconnect(fname,cv_to_ft(lvl-1),&grd[lvl-1]);      
            nstr.str("");
            grd[lvl].setinfo();
            ++excpt1;
         }
   }
   return(block::stop);
}

template<class GRD> block::ctrl mgrid<GRD>::matchboundaries(int lvl, block::ctrl ctrl_message) {

   if (ctrl_message == block::begin) excpt = 0;
   else excpt += ctrl_message;
   
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

template<class GRD> block::ctrl mgrid<GRD>::adapt(block::ctrl ctrl_message) {
   block::ctrl state;
   
   if (!adapt_flag) return(stop);
   
   if (ctrl_message == block::begin) excpt = 0;
   
   if (excpt == 0) {
      /* CALCULATE TARGET LENGTH FUNCTION */
      state = grd[0].length(ctrl_message);
      if (state != block::stop) return(state);
      ++excpt;
      ctrl_message = block::stay;
      mp_phase = -1;
   }
   
   excpt += ctrl_message;
   
   if (excpt == 1) {
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
   }
   else if(excpt == 2) {
      return(grd[0].adapt(block::begin,tolerance));
   }
   else {
      return(grd[0].adapt(ctrl_message,tolerance));
   }
   return(stop);
}

#endif
