/*
 *  blocks.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Thu Sep 26 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "blocks.h"
#include <utilities.h>
#include<iostream>
#include<map>
#include<string>
#include<sstream>


void blocks::load_constants(std::map<std::string,std::string> input[0]) {

   std::istringstream data(input[0]["itercrsn"]);   
   data >> itercrsn;
   std::cout << "#itercrsn:" << itercrsn << std::endl;
   data.clear();
   
   data.str(input[0]["iterrfne"]);   
   data >> iterrfne;
   std::cout << "#iterrfne:" << iterrfne << std::endl;
   data.clear();
   
   data.str(input[0]["njacobi"]);   
   data >> njacobi;
   std::cout << "#njacobi:" << njacobi << std::endl;
   data.clear();
   
   data.str(input[0]["ncycle"]);   
   data >> ncycle;
   std::cout << "#ncycle:" << ncycle << std::endl;
   data.clear();
   
   data.str(input[0]["ntstep"]);   
   data >> ntstep;
   std::cout << "#ntstep:" << ntstep << std::endl;
   data.clear(); 
   
   return;
}

void blocks::init(std::map<std::string,std::string> input[]) {
   int i,type;
   std::map<std::string,std::string>::const_iterator mi;
   std::map<std::string,std::string> merge;
   
   /* LOAD NUMBER OF GRIDS */
   std::istringstream data(input[0]["nblock"]);
   data >> nblock;
   std::cout << "#nblock:" << nblock << std::endl;
   data.clear();
   
   data.str(input[0]["ngrid"]);   
   data >> ngrid;
   std::cout << "#ngrid:" << ngrid << std::endl;
   data.clear();
   
   data.str(input[0]["mglvls"]);   
   data >> mglvls;
   std::cout << "#mglvls:" << mglvls << std::endl;
   data.clear();

   blk = new block *[nblock];

   for (i=0;i<nblock;++i) {
      merge = input[0];
      for (mi=input[i+1].begin(); mi != input[i+1].end(); ++mi)
         merge[mi->first] = mi->second;
      
      data.str(input[i+1]["blktype"]);   
      data >> type;
      std::cout << "#blktype:" << type << std::endl;
      data.clear();
      blk[i] = getnewblock(type);
      blk[i]->init(merge);
      blk[i]->load_const(merge);
      blk[i]->alloc(merge);
   }
   
   findmatch();
   matchboundaries();
   
   return;
}

void blocks::findmatch() {
   int j,k,grdlvl;
   
   for(grdlvl=0;grdlvl<ngrid;++grdlvl) {
      for(j=0;j<nblock;++j) 
         for(k=0;k<nblock;++k) 
            blk[j]->findmatch(grdlvl,blk[k]);
   }
   
   return;
}

   /* MATCH BOUNDARIES */
void blocks::matchboundaries() {
   int i,lvl,excpt;
   int state;
   
   for (lvl=0;lvl<ngrid;++lvl) {
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->matchboundaries(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
   }
}


void blocks::output(char *filename, FTYPE filetype) {
   int i;   
   char fnmcat[80],outname[80];

   /* ASSUME FOR NOW MESHES ARE LABELED a,b,c... */
   /* I HAVEN'T FIGURED OUT HOW THIS IS GOING TO WORK IN THE TOTALLY GENERAL CASE */
   if (nblock > 1) {
      for (i=0;i<nblock;++i) {
         strcpy(fnmcat,filename);
         strcat(fnmcat,".");
         number_str(outname, fnmcat, i, 1);
         blk[i]->output(outname,filetype);
      }
   }
   else {
      blk[0]->output(filename,filetype);
   }
   
   return;
}

void blocks::rsdl(int lvl) {
   int i,excpt;
   int state;
   
   excpt = 0;
   do {
      state = block::stop;
      for(i=0;i<nblock;++i)
         state &= blk[i]->rsdl(lvl,excpt);
      excpt += state;
   } while (state != block::stop);

   return;
}
      
      


void blocks::iterate(int lvl) {
/*****************************************/
/* JACOBI-ITERATION FOR MESH POSITION ****/
/*****************************************/
   int i,iter,excpt;
   int state;
   
   excpt = 0;
   do {
      state = block::stop;
      for(i=0;i<nblock;++i)
         state &= blk[i]->vddt(lvl,excpt);
      excpt += state;
   } while (state != block::stop);

   for(iter=0;iter<njacobi;++iter) {
      rsdl(lvl);
   
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->update(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
   }

   return;
}

void blocks::cycle(int vw, int lvl) {
   int i,vcount,iter,excpt;  // DON'T MAKE THESE STATIC SCREWS UP RECURSION
   int state;
   
   for (vcount=0;vcount<vw;++vcount) {
   
      for(iter=0;iter<itercrsn;++iter)
         iterate(lvl);
      
      if (lvl == mglvls-1) return;
      
      rsdl(lvl);
      
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->mg_getfres(lvl+1,excpt);
         excpt += state;
      } while (state != block::stop);
          
      cycle(vw, lvl+1);

      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->mg_getcchng(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
         
      for(iter=0;iter<iterrfne;++iter)
         iterate(lvl);
   }

   return;
}

void blocks::go() {
   int i,step;
   char outname[100];
   
   for(step = 1;step<=ntstep;++step) {
      tadvance();
      for(i=0;i<ncycle;++i) {
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

void blocks::tadvance() {
   int i,lvl,excpt;
   int state;
   
   for (lvl=0;lvl<ngrid;++lvl) {
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->tadvance(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
   }

   return;
}

void blocks::restructure() {
   int i,lvl,excpt;
   int state;
   
   matchboundaries();
   
   excpt = 0;
   do {
      state = block::stop;
      for(i=0;i<nblock;++i)
         state &= blk[i]->adapt(excpt);
      excpt += state;
   } while (state != block::stop);
   
   for(lvl=1;lvl<ngrid;++lvl) {
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->reconnect(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
   }
            
   return;
}