#include "mesh.h"
#include <cstdio>
#include <stdlib.h>

void mesh::findmatch(class mesh& tgt) {
   int i,j;

   for(i=0;i<nvbd;++i) {
      for(j=0;j<tgt.nvbd;++j) {
         if (&tgt == this && i == j) continue;  // CAN'T MATCH VERTEX TO ITSELF */
         vbdry[i]->match(tgt.vbdry[j]);
      }
   }
   
   for(i=0;i<nsbd;++i) {
      for(j=0;j<tgt.nsbd;++j) {
         if (&tgt == this && i == j) continue;  // CAN'T MATCH SIDE TO ITSELF */
         if (sbdry[i]->match(tgt.sbdry[j])) break;
      }
   }
   
   return;
}

