/*
 *  r_length.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Wed Jan 09 2002.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"r_mesh.h"

void r_mesh::length1() {
   int i;
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) 
      sbdry[i]->sendy(vlngth,0,0,1);

   return;
   
}

void r_mesh::length_mp() {
   int i;
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) 
      sbdry[i]->rcvy(vlngth,0,0,1);
  
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) 
      sbdry[i]->sendx(vlngth,0,0,1);
}


void r_mesh::length2() {
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(int i=0;i<nsbd;++i) 
      sbdry[i]->rcvx(vlngth,0,0,1);
   
   return;
}
