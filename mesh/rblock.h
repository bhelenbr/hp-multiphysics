#include"r_mesh.h"
#include"block.h"

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#ifndef _r_block_h_
#define _r_block_h_


/* GENERIC MULTIGRID BLOCK */

template<class GRD> class mgrid : public block {
   protected:
      std::ostream *log;
      int ngrid, mp_phase;
      typename GRD::gbl gbl_store;
      GRD *grd;
      FLT tolerance;
   
   public:
      void init(std::map <std::string,std::string>& input, std::ostream *inlog = 0);
      void load_const(std::map <std::string,std::string>& input) {}
      void alloc(std::map<std::string,std::string>& input);
      void input(char *filename) {
         grd[0].in_mesh(filename,text);
      }
      void output(char *filename, FTYPE filetype = easymesh) {
         grd[0].setbcinfo();
         grd[0].out_mesh(filename,filetype);
      }
      void coarsenchk(const char *fname);
      block::ctrl reconnect(int lvl, int excpt);
      block::ctrl matchboundaries(int lvl, int excpt);
      int comm_entity_size(int grdlvl);
      int comm_entity_list(int grdlvl, int *list);
      boundary* vbdry(int grdlvl, int num);
      boundary* sbdry(int grdlvl, int num);
      boundary* fbdry(int grdlvl, int num);
      
      block::ctrl tadvance(int lvl, int excpt) {return(stop);}
      block::ctrl rsdl(int lvl, int excpt) {return(stop);}
      void maxres() {}
      block::ctrl vddt(int lvl, int excpt) {return(stop);}
      block::ctrl update(int lvl, int excpt) {return(stop);}
      block::ctrl mg_getfres(int lvl, int excpt) {return(stop);}
      block::ctrl mg_getcchng(int lvl, int excpt) {return(stop);}
      block::ctrl adapt(int excpt);
};

template<class GRD> void mgrid<GRD>::init(std::map <std::string,std::string>& input, std::ostream *inlog) {

   if (inlog) {
      log = inlog;
   } 
   else {
      log = &std::cout;
   }
   
   /* LOAD NUMBER OF GRIDS */
   std::istringstream data(input["ngrid"]);
   data >> ngrid;
   *log << "#ngrid: " << ngrid << std::endl;
   data.clear();

   grd = new GRD[ngrid];
   for(int i = 0; i<ngrid;++i)
      grd[i].log = log;
      
}

template<class GRD> void mgrid<GRD>::alloc(std::map<std::string,std::string>& input) {
   int i;
   
	FLT grwfac;
   std::istringstream data(input["growth factor"]);
   data >> grwfac;
   *log << "#growth factor: " << grwfac << std::endl;
   data.clear();
	
	int filetype;
   data.str(input["filetype"]);
   data >> filetype;
   *log << "#filetype: " << filetype << std::endl;
   data.clear();
	
   std::string filename;
   data.str(input["mesh"]);
   data >> filename;
   *log << "#mesh: " << filename << std::endl;
   data.clear();
   
   grd[0].in_mesh(filename.c_str(),static_cast<FTYPE>(filetype),grwfac);
   
   /* CREATE COARSE MESHES */
   for(i=1;i<ngrid;++i)
      reconnect(i,0);

   grd[0].allocate(0,&gbl_store);
   for(i=1;i<ngrid;++i)
      grd[i].allocate(1,&gbl_store);

   return;
}

template<class GRD> int mgrid<GRD>::comm_entity_size(int grdlvl) {
   return(grd[grdlvl].comm_entity_size());
}

template<class GRD> int mgrid<GRD>::comm_entity_list(int grdlvl, int *list) {
   return(grd[grdlvl].comm_entity_list(list));
}

template<class GRD> boundary* mgrid<GRD>::vbdry(int grdlvl, int num) {
   return(grd[grdlvl].vbdry[num]);
}

template<class GRD> boundary* mgrid<GRD>::sbdry(int grdlvl, int num) {
   return(grd[grdlvl].sbdry[num]);
}

template<class GRD> boundary* mgrid<GRD>::fbdry(int grdlvl, int num) {
   return(0);
}

template<class GRD> block::ctrl mgrid<GRD>::reconnect(int grdnum, int phase) {
#define OLDRECONNECT
#ifdef OLDRECONNECT
   grd[grdnum].coarsen(1.6,grd[grdnum-1]);
#else
   grd[grdnum].coarsen2(1.5,grd[grdnum-1],temp_hp);
#endif
   grd[grdnum].setbcinfo();
   grd[grdnum].setfine(grd[grdnum-1]);
   grd[grdnum-1].setcoarse(grd[grdnum]);

   return(stop);
}

template<class GRD> void mgrid<GRD>::coarsenchk(const char *fname) {
   int i;
   char name[100];

   for(i = 1; i< ngrid; ++i) {
      number_str(name,fname,i,1);
      grd[i].setbcinfo();
      grd[i].checkintegrity();
      grd[i].out_mesh(name);
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




/* THIS IS A SEQUENCE OF MESHES & STORAGE FOR A DEFORMABLE BLOCK */
template<class GRD> class rblock : public mgrid<GRD> {

   public:
      void load_const(std::map <std::string,std::string>& input);
      void maxres() { grd[0].maxres();}
      block::ctrl rsdl(int lvl, int excpt);
      block::ctrl vddt(int lvl, int excpt);
      block::ctrl update(int lvl, int excpt);
      block::ctrl mg_getfres(int lvl, int excpt);
      block::ctrl mg_getcchng(int lvl, int excpt);
      block::ctrl tadvance(int lvl, int excpt);
};



template<class GRD> void rblock<GRD>::load_const(std::map <std::string,std::string>& input) {

   std::istringstream data(input["tolerance"]);
   data >> tolerance; 
   *log << "#tolerance: " << tolerance << std::endl;
   data.clear();
   
   data.str(input["fadd"]);
   data >> GRD::fadd;
   *log << "#fadd: " << GRD::fadd << std::endl;
   data.clear();

   data.str(input["vnn"]);
   data >> GRD::vnn; 
   *log << "#vnn: " << GRD::vnn << std::endl;
   data.clear();
      
   return;
}

template<class GRD> block::ctrl rblock<GRD>::tadvance(int lvl, int execpoint) {
   
   if (lvl == 0) {
      switch (execpoint) {
         case (0):
            /* SETUP SPRING CONSTANTS  */
            mp_phase = 0;
            grd[0].rklaplace();
            return(advance);
            
         case (1):
            return(static_cast<block::ctrl>(grd[0].msgpass(mp_phase++)));
            
         case (2):
            grd[0].kvoli();
            return(advance);

#ifdef FOURTH
#define P2 2
         case (3):
            grd[0].rsdl1();
            mp_phase = 0;
            return(advance);
            
         case (4):
            return(static_cast<ctrl>(grd[0].msgpass(mp_phase++)));
#else
#define P2 0
#endif            
         case (3+P2):
            mp_phase = 0;
            grd[0].zero_source();
            grd[0].rsdl();
            return(advance);

         case (4+P2):
            return(static_cast<block::ctrl>(grd[0].msgpass(mp_phase++)));
         
         case (5+P2):
            grd[0].rsdl_finalrcv();
            grd[0].sumsrc();
            grd[0].tadvance();
            return(stop);

         default:
            *log << "error in control flow tadvance 1" << std::endl;
            exit(1);
      }
   }
   else {
#ifdef GEOMETRIC
      switch (execpoint) {
         case (0):
            /* SETUP SPRING CONSTANTS  */
            mp_phase = 0;
            grd[lvl].rklaplace();
            return(advance);
            
         case (1):
            return(static_cast<block::ctrl>(grd[lvl].msgpass(mp_phase++)));
            
         case (2):
            grd[lvl].kvoli();
            return(stop);
         
         default:
            *log << "error in control flow tadvance 2" << std::endl;
            exit(1);
      }
   }
#else
   /* USE MULTIGRID INTERPOLATION (ALGEBRAIC) */
   /* MUST BE DONE THIS WAY FOR SPRING METHOD */
   /* SETUP FIRST MESH */
      switch (phase) {
         case(0):
            mp_phase = 0;
            grd[lvl].rkmgrid();
            return(advance);
         case(1):
            return(static_cast<ctrl>(grd[lvl].msgpass(mp_phase++)));
         case(2):
            grd[lvl].kvoli();
            return(stop);
      }
   }
#endif      

   return(stop);
}

template<class GRD> block::ctrl rblock<GRD>::rsdl(int lvl, int excpt) {

   switch (excpt) {
#ifdef FOURTH
      case(0):
         mp_phase = 0;
         grd[lvl].rsdl1();
         return(advance);
      case(1):
         return(static_cast<ctrl>(grd[lvl].msgpass(mp_phase++)));
#endif
      case(0+P2):
         mp_phase = 0;
         grd[lvl].rsdl();
         return(advance);
         
      case(1+P2):
         return(static_cast<block::ctrl>(grd[lvl].msgpass(mp_phase++)));
         
      case(2+P2):
         grd[lvl].rsdl_finalrcv();
         return(stop);
         
      default:
         *log << "flow control error, rsdl" << std::endl;
         exit(1);
   }
   
   return(stop);
}

template<class GRD> block::ctrl rblock<GRD>::vddt(int lvl, int excpt) {

   switch (excpt) {
#ifdef FOURTH
      case(0):
         mp_phase = 0;
         grd[lvl].vddt1();
         return(advance);
      case(1):
         return(static_cast<ctrl>(grd[lvl].msgpass(mp_phase++)));
#endif
      case(0+P2):
         mp_phase = 0;
         grd[lvl].vddt();
         return(advance);
         
      case(1+P2):
         return(static_cast<block::ctrl>(grd[lvl].msgpass(mp_phase++)));
         
      case(2+P2):
         grd[lvl].vddti();
         return(stop);
   }
   
   *log << "flow control error: vddt" << std::endl;
   exit(1);
   
   return(stop);
}

template<class GRD> block::ctrl rblock<GRD>::update(int lvl, int excpt) {
   grd[lvl].update();
   return(stop);
}

template<class GRD> block::ctrl rblock<GRD>::mg_getfres(int lvl, int excpt) {
   grd[lvl].mg_getfres();
   return(stop);
}

template<class GRD> block::ctrl rblock<GRD>::mg_getcchng(int lvl, int excpt) {
   grd[lvl].mg_getcchng();
   return(stop);
}
#endif
