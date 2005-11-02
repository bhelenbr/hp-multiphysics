/*
 *  adapt.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on 12/2/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "boundary.h"

int nlst; 

block::ctrl mesh::adapt(int excpt, FLT tolsize) {
   int i;
#ifdef MPISRC
   MPI_Status status;
#endif

   switch (excpt) {
      case(0):         
         setup_for_adapt();
         return(block::advance);
         
      case(1):
         /* SWAP EDGES */
         swap(1.0e-10); 
         return(block::advance);
         
      case(2):         
         bdry_yaber(tolsize);
         return(block::advance);

      case(3):
         for(i=0;i<nsbd;++i) 
            sbdry(i)->master_slave_transmit();
         return(block::advance);
         
      case(4):
         /* COARSEN MATCHING BOUNDARIES */
         bdry_yaber1();
         return(block::advance);
         
      case(5):
         /* INTERIOR SIDE COARSENING */
         yaber(tolsize);
         return(block::advance);
         
      case(6):
         bdry_rebay(tolsize);
         return(block::advance);
         
      case(7):
         for(i=0;i<nsbd;++i) 
            sbdry(i)->master_slave_transmit();
         return(block::advance);
         
      case(8):
         bdry_rebay1();
         return(block::advance);
         
      case(9):
         rebay(tolsize);
         return(block::advance);
         
      case(10):
         cleanup_after_adapt();
         return(block::stop);
   }
 
   return(block::stop);
}


void mesh::setup_for_adapt() {
   int i, v0, v1;
   
#ifdef PV3
   sim::pv3_mesh_changed = true;
#endif

   /* SET-UP ADAPTION TRACKING STUFF */
   /* For vertices 0 = untouched, 1 = touched, 2 = deleted, 3 = special */
   /* For sides 0 = untouced, 1 = touched, 2 = deleted */
   /* For triangles 0 = untouched, 1 = touched, 2 = deleted, 3 = searched */
   for(i=0;i<maxvst;++i)
      td(i).info = 0;
      
   /* Back reference into list of entities needing refinement */
   for(i=0;i<maxvst;++i)
      vd(i).info = -1;
      
   /* MARK BEGINNING/END OF SIDE GROUPS & SPECIAL VERTEX POINTS */
   /* THESE SHOULD NOT BE DELETED */
   for(i=0;i<nsbd;++i) {
      v0 = sd(sbdry(i)->el(0)).vrtx(0);
      v1 = sd(sbdry(i)->el(sbdry(i)->nel-1)).vrtx(1);
      td(v0).info |=  VSPEC;
      td(v1).info |=  VSPEC;
   }
   
   for(i=0;i<nvbd;++i)
      td(vbdry(i)->v0).info |= VSPEC;
      
   return;
}
      
      
void mesh::cleanup_after_adapt() {
   int i,j,sind,v0;
   
   /* DELETE SIDES FROM BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i) {
      for(j=sbdry(i)->nel-1;j>=0;--j) {
         if (td(sbdry(i)->el(j)).info&SDLTE) 
            sbdry(i)->el(j) = sbdry(i)->el(--sbdry(i)->nel);
      }
      sbdry(i)->reorder();
   }
   
   /* UPDATE BOUNDARY CONDITION DATA */
   for (i=0;i<nsbd;++i) {
      /* THIS IS PROBABLY UNNECESSARY SINCE FIRST */
      /* & LAST VERTEX SHOULD NEVER BE CHANGED */
      sind = sbdry(i)->el(0);
      v0 = sd(sind).vrtx(0);
      if (td(v0).info&VTOUC) {
         updatevdata_bdry(i,0,0);
         td(v0).info &= ~VTOUC;
      }
      else movevdata_bdry(i,0,0);
      
      for(j=0;j<sbdry(i)->nel;++j) {
         sind = sbdry(i)->el(j);
         v0 = sd(sind).vrtx(1);
         if (td(v0).info&VTOUC) {
            updatevdata_bdry(i,j,1);
            td(v0).info &= ~VTOUC;
         }
         else movevdata_bdry(i,j,1);
         if (td(sind).info&STOUC) {
            updatesdata_bdry(i,j);
            td(sind).info &= ~STOUC;
         }
         else movesdata_bdry(i,j);
      }      
   }
   
   /* DELETE LEFTOVER VERTICES */
   /* VINFO > NVRTX STORES VRTX MOVEMENT HISTORY */         
   for(i=0;i<nvrtx;++i) {
      if (td(i).info&VTOUC) {
         vd(i).info = -2;
         updatevdata(i);
      }
      if (td(i).info&VDLTE) dltvrtx(i);
   }
   
   /* FIX BOUNDARY CONDITION POINTERS */
   for(i=0;i<nvbd;++i)
      if (vd(vbdry(i)->v0).info > -1) 
         vbdry(i)->v0 = vd(vbdry(i)->v0).info;  
                        
   /* CLEAN UP SIDES */
   /* SINFO WILL END UP STORING -1 UNTOUCHED, -2 TOUCHED, or INITIAL INDEX OF UNTOUCHED SIDE */
   /* SINFO > NSIDE WILL STORE MOVEMENT LOCATION */  /* TEMPORARY HAVEN"T TESTED THIS */
   for(i=0;i<nside;++i) {
      sd(i).info = -1;
      if (td(i).info&STOUC) {
         updatesdata(i);
         sd(i).info = -2;
      }
      if (td(i).info&SDLTE) dltsd(i);
   }
         
   /* FIX BOUNDARY CONDITION POINTERS */
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry(i)->nel;++j) 
         if (sbdry(i)->el(j) >= nside) 
            sbdry(i)->el(j) = sd(sbdry(i)->el(j)).info; 
      
   /* CLEAN UP DELETED TRIS */
   /* TINFO < NTRI STORES INDEX OF ORIGINAL TRI ( > 0), TINFO = 0 -> UNMOVED */
   /* TINFO > NTRI STORES TRI MOVEMENT HISTORY */
   for(i=0;i<ntri;++i) {
      if (td(i).info&TDLTE) dlttri(i);
      else if (td(i).info&TTOUC) {
         td(i).info = -2;
         updatetdata(i);
      }
   }
      
   bdrylabel();
   
   cnt_nbor();
      
   return;
}


void mesh::putinlst(int sind) {
   int i, temp, top, bot, mid;
      
   /* CREATE ORDERED LIST OF SIDES SMALLEST FSCR1 TO LARGEST */
   bot = 0;
   if (nlst > 0) {
      top = 0;
      bot = nlst-1;
      if (fscr1(sind) < fscr1(sd(top).info)) {
         bot = 0;
      }
      else if (fscr1(sind) > fscr1(sd(bot).info)) {
         bot = nlst;
      }
      else {
         while(top < bot-1) {
            mid = top + (bot -top)/2;
            if (fscr1(sind) > fscr1(sd(mid).info))
               top = mid;
            else
               bot = mid;
         }
      }
      for(i=nlst-1;i>=bot;--i) {
         temp = sd(i).info;
         sd(i+1).info = temp;
         vd(temp).info = i+1;
      }
   }
   sd(bot).info= sind;
   vd(sind).info = bot;
   ++nlst;

   assert(nlst < maxvst -1);
   
   return;
}

void mesh::tkoutlst(int sind) {
   int bgn,temp,i;
   
   bgn = vd(sind).info;
   for(i=bgn+1;i<nlst;++i) {
      temp = sd(i).info;
      sd(i-1).info = temp;
      vd(temp).info = i-1;
   }
   vd(sind).info = -1;
   --nlst;
   sd(nlst).info = -1;
   
   return;
}


      

   
   


      

