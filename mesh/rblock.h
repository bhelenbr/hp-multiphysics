#include"r_mesh.h"
#include"block.h"

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#ifndef _r_block_h_
#define _r_block_h_

/* THIS IS A SEQUENCE OF MESHES & STORAGE FOR A DEFORMABLE BLOCK */
template<class GRD> class rblock : public block {
   protected:
      int ngrid, mp_phase;
      typename GRD::gbl gbl_store;
      GRD *grd;
      FLT tolerance;
      
   public:
      void init(std::map <std::string,std::string>& input);
      void load_const(std::map <std::string,std::string>& input);
      void alloc(std::map <std::string,std::string>& input);
      void input(char *filename) {
         grd[0].in_mesh(filename,text);
      }
      void output(char *filename, FTYPE filetype = easymesh) {
         grd[0].mesh::setbcinfo();
         grd[0].out_mesh(filename,filetype);
      }
      void coarsenchk(const char *fname);
      void maxres() { grd[0].maxres();}
      control_state rsdl(int lvl, int excpt);
      control_state vddt(int lvl, int excpt);
      control_state update(int lvl, int excpt);
      control_state mg_getfres(int lvl, int excpt);
      control_state mg_getcchng(int lvl, int excpt);
      control_state reconnect(int lvl, int excpt);
      control_state tadvance(int lvl, int excpt);
      control_state adapt(int excpt);
      control_state matchboundaries(int lvl, int excpt);
      int comm_entity_size(int grdlvl);
      int comm_entity_list(int grdlvl, int *list);
      boundary* vbdry(int grdlvl, int num);
      boundary* sbdry(int grdlvl, int num);
      boundary* fbdry(int grdlvl, int num);
};

template<class GRD> void rblock<GRD>::init(std::map <std::string,std::string>& input) {

   /* LOAD NUMBER OF GRIDS */
   std::istringstream data(input["ngrid"]);
   data >> ngrid;
   std::cout << "#ngrid: " << ngrid << std::endl;
   data.clear();

   grd = new GRD[ngrid];
}

template<class GRD> void rblock<GRD>::load_const(std::map <std::string,std::string>& input) {

   std::istringstream data(input["tolerance"]);
   data >> tolerance; 
   std::cout << "#tolerance: " << tolerance << std::endl;
   data.clear();
   
   data.str(input["fadd"]);
   data >> GRD::fadd;
   std::cout << "#fadd: " << GRD::fadd << std::endl;
   data.clear();

   data.str(input["vnn"]);
   data >> GRD::vnn; 
   std::cout << "#vnn: " << GRD::vnn << std::endl;
   data.clear();
      
   return;
}

template<class GRD> void rblock<GRD>::alloc(std::map<std::string,std::string>& input) {
   int i;
   
	FLT grwfac;
   std::istringstream data(input["growth factor"]);
   data >> grwfac;
   std::cout << "#growth factor: " << grwfac << std::endl;
   data.clear();
	
	int filetype;
   data.str(input["filetype"]);
   data >> filetype;
   std::cout << "#filetype: " << filetype << std::endl;
   data.clear();
	
   std::string filename;
   data.str(input["mesh"]);
   data >> filename;
   std::cout << "#mesh: " << filename << std::endl;
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

template<class GRD> int rblock<GRD>::comm_entity_size(int grdlvl) {
   return(grd[grdlvl].comm_entity_size());
}

template<class GRD> int rblock<GRD>::comm_entity_list(int grdlvl, int *list) {
   return(grd[grdlvl].comm_entity_list(list));
}

template<class GRD> boundary* rblock<GRD>::vbdry(int grdlvl, int num) {
   return(grd[grdlvl].vbdry[num]);
}

template<class GRD> boundary* rblock<GRD>::sbdry(int grdlvl, int num) {
   return(grd[grdlvl].sbdry[num]);
}

template<class GRD> boundary* rblock<GRD>::fbdry(int grdlvl, int num) {
   return(0);
}

template<class GRD> block::control_state rblock<GRD>::reconnect(int grdnum, int phase) {
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

template<class GRD> void rblock<GRD>::coarsenchk(const char *fname) {
   int i;
   char name[100];

   for(i = 1; i< ngrid; ++i) {
      number_str(name,fname,i,1);
      grd[i].mesh::setbcinfo();
      grd[i].checkintegrity();
      grd[i].out_mesh(name);
      grd[i].setbcinfo();
   }
   return;
}

template<class GRD> block::control_state rblock<GRD>::tadvance(int lvl, int execpoint) {
   
   if (lvl == 0) {
      switch (execpoint) {
         case (0):
            /* SETUP SPRING CONSTANTS  */
            mp_phase = 0;
            grd[0].rklaplace();
            return(advance);
            
         case (1):
            return(static_cast<control_state>(grd[0].msgpass(mp_phase++)));
            
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
            return(static_cast<control_state>(grd[0].msgpass(mp_phase++)));
#else
#define P2 0
#endif            
         case (3+P2):
            mp_phase = 0;
            grd[0].zero_source();
            grd[0].rsdl();
            return(advance);

         case (4+P2):
            return(static_cast<control_state>(grd[0].msgpass(mp_phase++)));
         
         case (5+P2):
            grd[0].rsdl_finalrcv();
            grd[0].sumsrc();
            grd[0].tadvance();
            return(stop);

         default:
            printf("error in control flow tadvance 1\n");
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
            return(static_cast<control_state>(grd[lvl].msgpass(mp_phase++)));
            
         case (2):
            grd[lvl].kvoli();
            return(stop);
         
         default:
            printf("error in control flow tadvance 2\n");
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
            return(static_cast<control_state>(grd[lvl].msgpass(mp_phase++)));
         case(2):
            grd[lvl].kvoli();
            return(stop);
      }
   }
#endif      

   return(stop);
}

template<class GRD> block::control_state rblock<GRD>::adapt(int excpt) {
   
   switch(excpt) {
      case(0):
         mp_phase = 0;
         grd[0].length1();
         return(advance);
      case(1):
         return(static_cast<control_state>(grd[0].msgpass(mp_phase++)));
      case(2):
         grd[0].length2();
         return(advance);
      case(3):
         grd[0].adapt(tolerance);
         return(stop);
   }
   
   printf("control flow error: adapt\n");
   exit(1);
   
   return(stop);
}

template<class GRD> block::control_state rblock<GRD>::rsdl(int lvl, int excpt) {

   switch (excpt) {
#ifdef FOURTH
      case(0):
         mp_phase = 0;
         grd[lvl].rsdl1();
         return(advance);
      case(1):
         return(static_cast<control_state>(grd[lvl].msgpass(mp_phase++)));
#endif
      case(0+P2):
         mp_phase = 0;
         grd[lvl].rsdl();
         return(advance);
         
      case(1+P2):
         return(static_cast<control_state>(grd[lvl].msgpass(mp_phase++)));
         
      case(2+P2):
         grd[lvl].rsdl_finalrcv();
         return(stop);
         
      default:
         printf("flow control error, rsdl\n");
         exit(1);
   }
   
   return(stop);
}

template<class GRD> block::control_state rblock<GRD>::vddt(int lvl, int excpt) {

   switch (excpt) {
#ifdef FOURTH
      case(0):
         mp_phase = 0;
         grd[lvl].vddt1();
         return(advance);
      case(1):
         return(static_cast<control_state>(grd[lvl].msgpass(mp_phase++)));
#endif
      case(0+P2):
         mp_phase = 0;
         grd[lvl].vddt();
         return(advance);
         
      case(1+P2):
         return(static_cast<control_state>(grd[lvl].msgpass(mp_phase++)));
         
      case(2+P2):
         grd[lvl].vddti();
         return(stop);
   }
   
   printf("flow control error: vddt\n");
   exit(1);
   
   return(stop);
}

template<class GRD> block::control_state rblock<GRD>::update(int lvl, int excpt) {
   grd[lvl].update();
   return(stop);
}

template<class GRD> block::control_state rblock<GRD>::mg_getfres(int lvl, int excpt) {
   grd[lvl].mg_getfres();
   return(stop);
}

template<class GRD> block::control_state rblock<GRD>::mg_getcchng(int lvl, int excpt) {
   grd[lvl].mg_getcchng();
   return(stop);
}

template<class GRD> block::control_state rblock<GRD>::matchboundaries(int lvl, int excpt) {
   
   switch (excpt) {
      case(0):
         mp_phase = 0;
         grd[lvl].matchboundaries1();
         return(advance);
      case(1):
         return(static_cast<control_state>(grd[lvl].msgpass(mp_phase++)));
      case(2):
         grd[lvl].matchboundaries2();
         return(stop);
   }
   
	printf("control flow error matchboundaries\n");
   exit(1);
   
   return(stop);
}
#endif
