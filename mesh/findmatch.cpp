#include"blocks.h"
#include<cstdio>
#include<stdlib.h>

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
   
   matches = 0;
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & ALLD_MP) {
         for(j=0;j<tgt.nsbd;++j) {
            if (sbdry[i].type == tgt.sbdry[j].type) {
               if (&tgt == this && i == j) continue;  // CAN'T MATCH SIDE TO ITSELF */
               sbdry[i].adjmesh = &tgt;
               sbdry[i].adjbnum = j;
               sbdry[i].isfrst = 1 -tgt.sbdry[j].isfrst;
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

void mesh::zerobdryisfrst() {
   int i;
   
   for(i=0;i<nvbd;++i)
      vbdry[i].isfrst = 0;
      
   for(i=0;i<nsbd;++i)
      sbdry[i].isfrst = 0;
      
   return;
}

void mesh::matchboundaries1() {
   int i,j,v0,sind,bnum,count;
   class mesh *tgt;
   
   /* SEND POSITIONS & NUMBER TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & ALLD_MP) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         
         /* SEND NUMBER AS FLOAT? */
         tgt->sbuff[bnum][count++] = sbdry[i].num;
         
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = vrtx[v0][0];
            tgt->sbuff[bnum][count++] = vrtx[v0][1];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = vrtx[v0][0];
         tgt->sbuff[bnum][count++] = vrtx[v0][1];
      }
   }

   return;
   
}

void mesh::matchboundaries2() {
   int i,j,v0,sind,count,num;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & PRDX_MASK) {
         count = 0;
         
         /*	RECV NUMBER */
         num = static_cast<int>(sbuff[i][count++]);
         if (num != sbdry[i].num) {
            printf("non matching number of boundaries %d\n",sbdry[i].type);
            exit(1);
         }
         
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            ++count; // SKIP X DIRECTION INFO
            vrtx[v0][1] = 0.5*(vrtx[v0][1] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         ++count; // SKIP X DIRECTION INFO
         vrtx[v0][1] = 0.5*(vrtx[v0][1] +sbuff[i][count++]);
         continue;
      }

      if (sbdry[i].type & PRDY_MASK) {
         count = 0;
         
         /*	RECV NUMBER */
         num = static_cast<int>(sbuff[i][count++]);
         if (num != sbdry[i].num) {
            printf("non matching number of boundaries %d\n",sbdry[i].type);
            exit(1);
         }
         
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            vrtx[v0][0] = 0.5*(vrtx[v0][0] +sbuff[i][count++]);
            ++count; // SKIP X DIRECTION INFO

         }
         v0 = svrtx[sind][0];
         vrtx[v0][0] = 0.5*(vrtx[v0][0] +sbuff[i][count++]);
         continue;
      }

      if (sbdry[i].type & ALLD_MP) {
         count = 0;
         
         /*	RECV NUMBER */
         num = static_cast<int>(sbuff[i][count++]);
         if (num != sbdry[i].num) {
            printf("non matching number of boundaries %d\n",sbdry[i].type);
            exit(1);
         }
         
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            vrtx[v0][0] = 0.5*(vrtx[v0][0] +sbuff[i][count++]);
            vrtx[v0][1] = 0.5*(vrtx[v0][1] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         vrtx[v0][0] = 0.5*(vrtx[v0][0] +sbuff[i][count++]);
         vrtx[v0][1] = 0.5*(vrtx[v0][1] +sbuff[i][count++]);
      }      
   }
   
   return;
}
