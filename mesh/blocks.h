#include<map>
#include<string>
#include<sstream>
#include"mesh.h"

/* THIS IS A MULTIBLOCK MESH */
template<class BLK> class blocks {
   private:
      int nblock, ngrid;
      int niter, itercrsn, iterrfne;
      int ntstep;
      static const int lastphase = 1;
      BLK *blk;
      
   public:
      /* INITIALIZE MULTIBLOCK/MGRID MESH */
      void load_constants(std::map<std::string,std::string>& input);
      void init(std::map<std::string,std::string>& input);
      void findmatch();
      void matchboundaries();
      void output(char *filename, FILETYPE filetype);
      void input(char *filename) {}
      void go();

      /* ITERATION ON ALL BLOCKS */  
      void iterate(int mglvl);

      /* MGRID CYCLE vw = 1: v-cycle vw =2 w-cycle */
      void cycle(int vw, int lvl = 0);

      /* CALCULATE DEFORMATION SOURCE TERM ON FINEST MESH */
      void ksrc();
      
      /* PRINT ERRORS */
      inline void maxres() {
         for(int i=0;i<nblock;++i)
            blk[i].lvl(0).maxres();
      }
      
      /* MOVE BOUNDARIES */
      inline void tadvance() {
         ksrc();
         
         for(int i=0;i<nblock;++i)
            blk[i].tadvance();
         return;
      }

      inline void restructure() {
         int i,phase;
         
         matchboundaries();
         
         for (i=0;i<nblock;++i) 
            blk[i].lvl(0).length1();
            
         for(phase=0;phase<lastphase;++phase)
            for(i=0;i<nblock;++i)
               blk[i].lvl(0).length_mp(phase);
            
         for(i=0;i<nblock;++i) {
            blk[i].lvl(0).length2(lastphase);
            blk[i].adapt();
            blk[i].reconnect();
         }
         
         return;
      }
};

template<class BLK> void blocks<BLK>::load_constants(std::map<std::string,std::string>& input) {

   std::istringstream data(input["itercrsn"]);   
   data >> itercrsn;
   std::cout << "#itercrsn:" << itercrsn << endl;
   data.clear();
   
   data.str(input["niter"]);   
   data >> niter;
   std::cout << "#niter:" << niter << endl;
   data.clear();
   
   data.str(input["ntstep"]);   
   data >> ntstep;
   std::cout << "#ntstep:" << ntstep << endl;
   data.clear(); 
   
   return;
}

template<class BLK> void blocks<BLK>::init(std::map<std::string,std::string>& input) {
   int i;
   
   /* LOAD NUMBER OF GRIDS */
   std::istringstream data(input["nblock"]);
   data >> nblock;
   std::cout << "#nblock:" << nblock << endl;
   data.clear();
   
   data.str(input["ngrid"]);   
   data >> ngrid;
   std::cout << "#ngrid:" << ngrid << endl;
   data.clear();

   blk = new BLK[nblock];

   for (i=0;i<nblock;++i) {
      blk[i].init(input);
      blk[i].load_const(input);
      blk[i].alloc(input);
   }
   
   findmatch();
   matchboundaries();
   
   return;
}

template<class BLK> void blocks<BLK>::findmatch() {
   int j,k,grdlvl;
   
   for(grdlvl=0;grdlvl<ngrid;++grdlvl) {
      for(j=0;j<nblock;++j) 
         for(k=0;k<nblock;++k) 
            blk[j].lvl(grdlvl).findmatch(blk[k].lvl(grdlvl));
   }
   
   return;
}

   /* MATCH BOUNDARIES */
template<class BLK> void blocks<BLK>::matchboundaries() {
   int i,j,phase;
   
   for(i=0;i<ngrid;++i) {
      for(j=0;j<nblock;++j)
         blk[j].lvl(i).matchboundaries1();
      
      for(j=0;j<nblock;++j)
         for(phase=0;phase<lastphase;++phase)
            blk[j].lvl(i).matchboundaries_mp(phase);
      
      for(j=0;j<nblock;++j)
         blk[j].lvl(i).matchboundaries2(lastphase);  
   }
}


template<class BLK> void blocks<BLK>::output(char *filename, FILETYPE filetype) {
   int i;   
   char fnmcat[80],outname[80];

   /* ASSUME FOR NOW MESHES ARE LABELED a,b,c... */
   /* I HAVEN'T FIGURED OUT HOW THIS IS GOING TO WORK IN THE TOTALLY GENERAL CASE */
   if (nblock > 1) {
      for (i=0;i<nblock;++i) {
         strcpy(fnmcat,filename);
         strcat(fnmcat,".");
         number_str(outname, fnmcat, i, 1);
         blk[i].output(outname,filetype);
      }
   }
   else {
      blk[0].output(filename,filetype);
   }
   
   return;
}

template<class BLK> void blocks<BLK>::iterate(int lvl) {
   int i,iter,phase;
      
/*****************************************/
/* JACOBI-ITERATION FOR MESH POSITION ****/
/*****************************************/

   for(i=0;i<nblock;++i)
      blk[i].lvl(lvl).vddt();

#ifdef FOURTH
   for (phase=0;phase<lastphase;++phase)
      for(i=0;i<nblock;++i)
         blk[i].lvl(lvl).vddt1_mp(phase);   
   
   for(i=0;i<nblock;++i)
      blk[i].lvl(lvl).vddt1(lastphase);
#endif
   for (phase=0;phase<lastphase;++phase)
      for(i=0;i<nblock;++i)
         blk[i].lvl(lvl).vddt_mp(phase);     
   
   for(i=0;i<nblock;++i)
      blk[i].lvl(lvl).vddti(lastphase);

   for(iter=0;iter<itercrsn;++iter) {
   
      for(i=0;i<nblock;++i)
         blk[i].lvl(lvl).rsdl();

#ifdef FOURTH
      for (phase=0;phase<lastphase;++phase)
         for(i=0;i<nblock;++i)
            blk[i].lvl(lvl).rsdl1_mp(phase);     
         
      for(i=0;i<nblock;++i)
         blk[i].lvl(lvl).rsdl1(lastphase);
#endif
   for (phase=0;phase<lastphase;++phase)
      for(i=0;i<nblock;++i)
         blk[i].lvl(lvl).rsdl_mp(phase);

      for(i=0;i<nblock;++i) 
         blk[i].lvl(lvl).update(lastphase);
   }

   return;
}

template<class BLK> void blocks<BLK>::cycle(int vw, int lvl) {
   int i,j,phase;  // DON'T MAKE THESE STATIC SCREWS UP RECURSION
   
   for (i=0;i<vw;++i) {
      iterate(lvl);
      
      if (lvl == ngrid-1) return;
      
      for(j=0;j<nblock;++j)
         blk[j].lvl(lvl).rsdl();
         
#ifdef FOURTH
      for(j=0;j<nblock;++j)
         blk[j].lvl(lvl).rsdl1_mp();      

      for(j=0;j<nblock;++j)
         blk[j].lvl(lvl).rsdl1();
#endif
      for (phase=0;phase<lastphase;++phase)
         for(j=0;j<nblock;++j)
            blk[j].lvl(lvl).rsdl_mp(phase);      
      
      for(j=0;j<nblock;++j)
         blk[j].lvl(lvl+1).mg_getfres(lastphase);
      
      cycle(vw, lvl+1);

      for(j=0;j<nblock;++j)
         blk[j].lvl(lvl).mg_getcchng();
   }

   return;
}

#define GEOMETRIC

template<class BLK> void blocks<BLK>::ksrc() {
   int i,j,phase;

#ifdef GEOMETRIC   
   /* SETUP SPRING CONSTANTS  */
   for(i=0;i<ngrid;++i) {
      for(j=0;j<nblock;++j)
         blk[j].lvl(i).rklaplace();
      
      for (phase=0;phase<lastphase;++phase)
         for(j=0;j<nblock;++j)
            blk[j].lvl(i).kvol_mp(phase);
               
      for(j=0;j<nblock;++j)
         blk[j].lvl(i).kvoli(lastphase);
   }
#else
   /* USE MULTIGRID INTERPOLATION (ALGEBRAIC) */
   /* MUST BE DONE THIS WAY FOR SPRING METHOD */
   /* SETUP FIRST MESH */
   for(j=0;j<nblock;++j) 
      blk[j].lvl(0).rklaplace();
   
   for (phase=0;phase<lastphase;++phase)
      for(j=0;j<nblock;++j)
         blk[j].lvl(0).kvol_mp(phase);
               
   for(j=0;j<nblock;++j)
      blk[j].lvl(0).kvoli(lastphase);
   
   /* SETUP COARSE GRIDS */
   for(i=1;i<ngrid;++i) {
      for(j=0;j<nblock;++j)
         blk[j].lvl(i).rkmgrid();
      
      for (phase=0;phase<lastphase;++phase)
         for(j=0;j<nblock;++j)
            blk[j].lvl(i).rkmgrid_mp(phase);
               
      for(j=0;j<nblock;++j)
         blk[j].lvl(i).rkmgridi(lastphase);
   }
#endif      

/* CALCULATE SOURCE TERM ON FINEST MESH */
   for(i=0;i<nblock;++i)
      blk[i].lvl(0).source();

#ifdef FOURTH
   for (phase=0;phase<lastphase;++phase)
      for(i=0;i<nblock;++i)
         blk[i].lvl(0).rsdl1_mp(phase);  

   for(i=0;i<nblock;++i)
      blk[i].lvl(0).rsdl1(lastphase);
#endif
   for (phase=0;phase<lastphase;++phase)
      for(i=0;i<nblock;++i)
         blk[i].lvl(0).rsdl_mp(phase);
   
   for(i=0;i<nblock;++i)
      blk[i].lvl(0).sumsrc(lastphase);

   return;
}

template<class BLK> void blocks<BLK>::go() {
   int i,step;
   char outname[100];
   
   for(step = 1;step<=ntstep;++step) {
      tadvance();
      for(i=0;i<niter;++i) {
         cycle(2);
         printf("%d ",i);
         maxres();
         printf("\n");
      }
      output("deformed",grid);
      restructure();
      number_str(outname, "end", step, 2);
      output(outname,grid);
   }
   return;
}


