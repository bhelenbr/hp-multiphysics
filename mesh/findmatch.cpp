#include"blocks.h"
#include<cstdio>

int mesh::findmatch(class mesh& tgt) {
   int i,j,matches,found;

   for(i=0;i<nvbd;++i) {
      if (vbdry[i].type & ALLD_MP) {
         found = 0;
         for(j=0;j<tgt.nvbd;++j) {
            if (vbdry[i].type == tgt.vbdry[j].type) {
               if (&tgt == this && i == j) continue;  // CAN'T MATCH VERTEX TO ITSELF */
               vbdry[i].adjmesh = &tgt;
               vbdry[i].adjbnum = j;
               found = 1;
               break;
            }
         }
         if (found) {
            printf("#matched vertex %d to %d\n",vbdry[i].el[0],vbdry[j].el[0]);
         }
      }
   }
   
   matches = 0;
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & ALLD_MP) {
         for(j=0;j<tgt.nsbd;++j) {
            if (sbdry[i].type == tgt.sbdry[j].type) {
               if (&tgt == this && i == j) continue;  // CAN'T MATCH SIDE TO ITSELF */
               sbdry[i].adjmesh = &tgt;
               sbdry[i].adjbnum = j;
               ++matches;
               break;
            }
         }
      }
   }
   
   return(matches);
}

int mesh::alld_mp() {
   int i,matches;
   
   matches = 0;
   for(i=0;i<nsbd;++i) 
      if (sbdry[i].type & ALLD_MP) ++matches;
      
   return(matches);
}