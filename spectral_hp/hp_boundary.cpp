/*
 *  hp_boundary.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/2/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"


void tri_hp::vmsgload(int phase, FLT *vdata) {
   int i;
      
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_snd(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_snd(phase,vdata);
   
   for(i=0;i<nsbd;++i)
      sbdry(i)->comm_prepare(phase);
   for(i=0;i<nvbd;++i)
      vbdry(i)->comm_prepare(phase);
   
   return;
}
int tri_hp::vmsgwait_rcv(int phase, FLT *vdata) {
   int stop = 1;
   int i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_wait(phase);
   for(i=0;i<nvbd;++i)
      stop &= vbdry(i)->comm_wait(phase);
      

   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_rcv(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_rcv(phase,vdata);
      
   return(stop);
}

int tri_hp::vmsgrcv(int phase, FLT *vdata) {
   int stop = 1,i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_nowait(phase);
   for(i=0;i<nvbd;++i)
      stop &= vbdry(i)->comm_nowait(phase);
      
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_rcv(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_rcv(phase,vdata);
      
   return(stop);
}

void tri_hp::smsgload(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) {
   int i;
      
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_snd(phase,sdata,bgnmode,endmode,modestride);
   
   for(i=0;i<nsbd;++i)
      sbdry(i)->comm_prepare(phase);
   
   return;
}

int tri_hp::smsgwait_rcv(int phase,FLT *sdata, int bgnmode, int endmode, int modestride) {
   int stop = 1;
   int i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_wait(phase);

   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_rcv(phase,sdata,bgnmode,endmode,modestride);
      
   return(stop);
}

int tri_hp::smsgrcv(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) {
   int stop = 1,i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_nowait(phase);
      
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_rcv(phase,sdata,bgnmode,endmode,modestride);
      
   return(stop);
}
