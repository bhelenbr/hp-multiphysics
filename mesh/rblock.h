#include"r_mesh.h"
#include"block.h"

#include<map>
#include<string>
#include<sstream>

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
      void findmatch(int grdlvl, block *match);
      void vmatch(int lvl, class vrtx_boundary *vin);
      void smatch(int lvl, class side_boundary *sin);
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
            return(static_cast<control_state>(grd[0].kvol_mp(mp_phase++)));
            
         case (2):
            grd[0].kvoli(mp_phase);
            return(advance);
            
         case (3):
            mp_phase = 0;
            grd[0].zero_source();
            grd[0].rsdl();
            grd[0].rsdl_snd(mp_phase);
            return(advance);
         
#ifdef FOURTH
#define P2 2
         case (4):
            grd[0].rsdl_rcv(mp_phase++);
            return(grd[0].rsdl1_snd(mp_phase));
            
         case (5):
            grd[0].rsdl1(mp_phase);
            mp_phase = 0;
            return(advance);
#else
#define P2 0
#endif

         case (4+P2):
            grd[0].rsdl_rcv(mp_phase++);
            return(static_cast<control_state>(grd[0].rsdl_snd(mp_phase)));
         
         case (5+P2):
            grd[0].rsdl_rcv(mp_phase);
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
            return(static_cast<control_state>(grd[lvl].kvol_mp(mp_phase++)));
            
         case (2):
            grd[lvl].kvoli(mp_phase);
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
            return(grd[lvl].rkmgrid_mp(mp_phase++));
         case(2):
            grd[lvl].rkmgridi(lastphase);
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
         return(static_cast<control_state>(grd[0].length_mp(mp_phase++)));
      case(2):
         grd[0].length2(mp_phase);
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
         return(grd[lvl].rsdl1_mp(mp_phase++));
#endif
      case(0+P2):
         mp_phase = 0;
         grd[lvl].rsdl();
         grd[lvl].rsdl_snd(mp_phase);
         return(advance);
         
      case(1+P2):
         grd[lvl].rsdl_rcv(mp_phase++);
         return(static_cast<control_state>(grd[lvl].rsdl_snd(mp_phase)));
         
      case(2+P2):
         grd[lvl].rsdl_rcv(mp_phase);
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
         return(grd[lvl].vddt1_mp(mp_phase++));
#endif
      case(0+P2):
         mp_phase = 0;
         grd[lvl].vddt();
         return(advance);
         
      case(1+P2):
         return(static_cast<control_state>(grd[lvl].vddt_mp(mp_phase++)));
         
      case(2+P2):
         grd[lvl].vddti(mp_phase);
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


template<class GRD> void rblock<GRD>::findmatch(int lvl, block *match) {
   int i,j;
   
   if (match != this) {
      for(i=0;i<grd[lvl].nvbd;++i)
         match->vmatch(lvl,grd[lvl].vbdry[i]);
         
      for(i=0;i<grd[lvl].nsbd;++i)
         match->smatch(lvl,grd[lvl].sbdry[i]);
   }
   else {
      for(i=0;i<grd[lvl].nvbd;++i) {
         for(j=0;j<grd[lvl].nvbd;++j) {
            if (i == j) continue;  // CAN'T MATCH VERTEX TO ITSELF */
            grd[lvl].vbdry[i]->match(grd[lvl].vbdry[j]);
         }
      }
      
      for(i=0;i<grd[lvl].nsbd;++i) {
         for(j=0;j<grd[lvl].nsbd;++j) {
            if (i == j) continue;  // CAN'T MATCH SIDE TO ITSELF */
            if (grd[lvl].sbdry[i]->match(grd[lvl].sbdry[j])) break;
         }
      }
   }
   
   return;
}

template<class GRD> void rblock<GRD>::vmatch(int lvl, class vrtx_boundary *vin) {
   int i;
   
   for(i=0;i<grd[lvl].nvbd;++i) 
      vin->match(grd[lvl].vbdry[i]);
   
   return;
}

template<class GRD> void rblock<GRD>::smatch(int lvl, class side_boundary *sin) {
   int i;
   
   for(i=0;i<grd[lvl].nsbd;++i) 
      sin->match(grd[lvl].sbdry[i]);
   
   return;
}

template<class GRD> block::control_state rblock<GRD>::matchboundaries(int lvl, int excpt) {
   
   switch (excpt) {
      case(0):
         mp_phase = 0;
         grd[lvl].matchboundaries1();
         return(advance);
      case(1):
         return(static_cast<control_state>(grd[lvl].matchboundaries_mp(mp_phase++)));
      case(2):
         grd[lvl].matchboundaries2(mp_phase);
         return(stop);
   }
   
	printf("control flow error matchboundaries\n");
   exit(1);
   
   return(stop);
}
#endif
