#include "mesh.h"
#include <cstdio>
#include <stdlib.h>

int mesh::comm_entity_size() {
   int i,tsize,nvcomm,nscomm;
   
   /*	VERTEX INFO */   
   nvcomm = 0;
   for(i=0;i<nvbd;++i) 
      if (vbdry[i]->is_comm()) ++nvcomm;
      
   tsize = 1 + 2*nvcomm; // bdry number & id
   
   /* SIDE INFO */
   nscomm = 0;
   for(i=0;i<nsbd;++i)
      if (sbdry[i]->is_comm()) ++nscomm; 
      
   tsize += 1 +4*nscomm; // bdry number, id, v0id, v1id
   
   tsize += 1;  // nfcomm = 0
   
   return(tsize);
}
   


int mesh::comm_entity_list(int *list) {
   int i,j,nvcomm,nscomm,tsize,v0,v0id;
   
   /* MAKE 1D PACKED LIST OF ALL INFORMATION ON COMMUNICATION BOUNDARIES */
   tsize = 0;

   /*	VERTEX INFO */   
   nvcomm = 0;
   for(i=0;i<nvbd;++i) 
      if (vbdry[i]->is_comm()) ++nvcomm;
      
   list[tsize++] = nvcomm;
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i]->is_comm()) {
         list[tsize++] = i;
         list[tsize++] = vbdry[i]->idnty();
      }
   }
   
   /* SIDE INFO */
   nscomm = 0;
   for(i=0;i<nsbd;++i)
      if (sbdry[i]->is_comm()) ++nscomm;
      
   list[tsize++] = nscomm;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i]->is_comm()) {
         list[tsize++] = i;
         list[tsize++] = sbdry[i]->idnty();
         v0 = svrtx[sbdry[i]->sd(0)][0];
         v0id = -1;
         for(j=0;j<nvbd;++j) {
            if (vbdry[j]->v() == v0) {
               v0id = vbdry[j]->idnty();
               break;
            }
         }
         list[tsize++] = v0id;
         v0 = svrtx[sbdry[i]->sd(sbdry[i]->nsd()-1)][1];
         v0id = -1;
         for(j=0;j<nvbd;++j) {
            if (vbdry[j]->v() == v0) {
               v0id = vbdry[j]->idnty();
               break;
            }
         }
         list[tsize++] = v0id;
      }
   }
   
   /* FACE BOUNDARIES */
   list[tsize++] = 0;
   
   return(tsize);
}

int mesh::msgpass(int phase) {
   int stop=1;

   for(int i=0;i<nsbd;++i) 
      stop &= sbdry[i]->rcv(phase);
   for(int i=0;i<nvbd;++i) 
      stop &= vbdry[i]->rcv(phase);
   
   
   if (!stop) {
      for(int i=0;i<nsbd;++i)
         sbdry[i]->snd(phase+1);
      for(int i=0;i<nvbd;++i)
         vbdry[i]->snd(phase+1);
   }
   
   return(stop);
}

/*	MAKE SURE MATCHING BOUNDARIES ARE AT EXACTLY THE SAME POSITIONS */
void mesh::matchboundaries1() {
   
   for(int i=0;i<nvbd;++i)
      vbdry[i]->sndpositions();
   for(int i=0;i<nsbd;++i) 
      sbdry[i]->sndpositions();
   
   /* FIRST PHASE OF SENDING */
   for(int i=0;i<nsbd;++i)
      sbdry[i]->snd(0);
   for(int i=0;i<nvbd;++i)
      vbdry[i]->snd(0);
}

void mesh::matchboundaries2() {
   for(int i=0;i<nsbd;++i)
      sbdry[i]->rcvpositions();
   for(int i=0;i<nvbd;++i)
      vbdry[i]->rcvpositions();
}

void mesh::length1() {
   int i;
   
   setlength();
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      sbdry[i]->loadbuff(vlngth,0,0,1);
   for(i=0;i<nvbd;++i)
      vbdry[i]->loadbuff(vlngth,0,0,1);
   

   for(i=0;i<nsbd;++i)
      sbdry[i]->snd(0);
   for(i=0;i<nvbd;++i)
      vbdry[i]->snd(0);

   return;
   
}

void mesh::length2() {
   int i;
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) 
      sbdry[i]->finalrcv(vlngth,0,0,1);
   for(i=0;i<nvbd;++i)
      vbdry[i]->finalrcv(vlngth,0,0,1);

   return;
}