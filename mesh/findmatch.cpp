#include "mesh.h"
#include <cstdio>
#include <stdlib.h>

void mesh::findmatch(class mesh& tgt) {
   int i,j,found;

   for(i=0;i<nvbd;++i) {
      if (vbdry[i].type & ALLD_MP) {
         found = 0;
         for(j=0;j<tgt.nvbd;++j) {
            if (vbdry[i].type == tgt.vbdry[j].type) {
               if (&tgt == this && i == j) continue;  // CAN'T MATCH VERTEX TO ITSELF */
               vbdry[i].adjmesh = &tgt;
               vbdry[i].adjbnum = j;
               vbdry[i].isfrst = 1 -tgt.vbdry[j].isfrst;
               found = 1;
               break;
            }
         }
         if (found) {
            printf("#matched vertex %d to %d\n",vbdry[i].el[0],vbdry[j].el[0]);
         }
      }
   }
   
   for(i=0;i<nsbd;++i) {
      for(j=0;j<tgt.nsbd;++j) {
         if (&tgt == this && i == j) continue;  // CAN'T MATCH SIDE TO ITSELF */
         if (sbdry[i]->match(tgt.sbdry[j])) goto NEXT;
      }
      printf("#Error: couldn't match boundary %d\n",sbdry[i]->idnty());
      exit(1);
NEXT:;
   }
   
   return;
}

