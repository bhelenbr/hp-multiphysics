#include<map>
#include<string>
#include"mesh.h"

/* THIS IS A SEQUENCE OF MESHES & STORAGE FOR A SINGLE MULTIGRID BLOCK */
template<class GRD> class block {
   protected:
      int ngrid;
      typename GRD::gbl gbl_store;
      GRD *grd;
      FLT tolerance;
      
   public:
      void init(std::map <std::string,std::string>& input);
      void load_const(std::map <std::string,std::string>& input);
      void alloc(std::map <std::string,std::string>& input);
      void output(char *filename, FILETYPE filetype = easymesh) {
         grd[0].mesh::setbcinfo();
         grd[0].out_mesh(filename,filetype);
      }
      void input(char *filename) {
         grd[0].in_mesh(filename,text);
      }
      inline GRD& lvl(int i) {return(grd[i]);}
      void reconnect();
      void coarsenchk(const char *fname);
      void tadvance();
      void adapt();
};

template<class GRD> void block<GRD>::init(std::map <std::string,std::string>& input) {

   /* LOAD NUMBER OF GRIDS */
   std::istringstream data(input["ngrid"]);
   data >> ngrid;
   std::cout << "#ngrid: " << ngrid << endl;
   data.clear();

   grd = new GRD[ngrid];
}

template<class GRD> void block<GRD>::load_const(std::map <std::string,std::string>& input) {

   std::istringstream data(input["tolerance"]);
   data >> tolerance; 
   std::cout << "#tolerance: " << tolerance << endl;
   data.clear();
   
   data.str(input["fadd"]);
   data >> GRD::fadd;
   std::cout << "#fadd: " << GRD::fadd << endl;
   data.clear();

   data.str(input["vnn"]);
   data >> GRD::vnn; 
   std::cout << "#vnn: " << GRD::vnn << endl;
   data.clear();
      
   return;
}

template<class GRD> void block<GRD>::alloc(std::map<std::string,std::string>& input) {
   int i;
   
	FLT grwfac;
   std::istringstream data(input["growth factor"]);
   data >> grwfac;
   std::cout << "#growth factor: " << grwfac << endl;
   data.clear();
	
	int filetype;
   data.str(input["filetype"]);
   data >> filetype;
   std::cout << "#filetype: " << filetype << endl;
   data.clear();
	
   std::string filename;
   data.str(input["mesh"]);
   data >> filename;
   std::cout << "#mesh: " << filename << endl;
   data.clear();
   
   grd[0].in_mesh(filename.c_str(),static_cast<FILETYPE>(filetype),grwfac);
   
   /* CREATE COARSE MESHES */
   reconnect();

   grd[0].allocate(0,&gbl_store);
   for(i=1;i<ngrid;++i)
      grd[i].allocate(1,&gbl_store);

   return;
}

#define OLDRECONNECT

template<class GRD> void block<GRD>::reconnect() {
   int i;

   for(i = 1; i< ngrid; ++i) {
#ifdef OLDRECONNECT
      grd[i].coarsen(1.6,grd[i-1]);
#else
      grd[i].coarsen2(1.5,grd[i-1],temp_hp);
#endif
      grd[i].setbcinfo();
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }

   return;
}

template<class GRD> void block<GRD>::coarsenchk(const char *fname) {
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

template<class GRD> void block<GRD>::tadvance() {
   grd[0].tadvance();
}

template<class GRD> void block<GRD>::adapt() {
   grd[0].yaber(1.0/tolerance,1,0.0);
   grd[0].setbcinfo();
   grd[0].out_mesh("coarse",grid);
   grd[0].treeupdate();
   grd[0].rebay(tolerance);
   grd[0].setbcinfo();
   grd[0].out_mesh("refine",grid);
}
   
