/*
 *  mg_alloc.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 30 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include <Carbon/Carbon.h>

FLT **hp_mgrid::cv00,**hp_mgrid::cv01,**hp_mgrid::cv10,**hp_mgrid::cv11;
FLT **hp_mgrid::e00,**hp_mgrid::e01,**hp_mgrid::e10,**hp_mgrid::e11;	
int hp_mgrid::size;

void hp_mgrid::allocate(struct hp_mgrid_glbls ginit) {
   
   if (spectral_hp::size == 0 or mesh::initialized == 0) {
      printf("must initialize mesh/spectral_hp first\n");
      exit(1);
   }
   
   if (size == 0) {
      mat_alloc(cv00,b.gpx,b.gpn,FLT);
      mat_alloc(cv01,b.gpx,b.gpn,FLT);
      mat_alloc(cv10,b.gpx,b.gpn,FLT);
      mat_alloc(cv11,b.gpx,b.gpn,FLT);
      	
      mat_alloc(e00,b.gpx,b.gpn,FLT);
      mat_alloc(e01,b.gpx,b.gpn,FLT);
      mat_alloc(e10,b.gpx,b.gpn,FLT);
      mat_alloc(e11,b.gpx,b.gpn,FLT);
      size = b.p;
   }
   else {
      if (size < b.p) {
         printf("allocate from largest too smallest\n");
         exit(1);
      }
   }

/*	INITIALIZE GLOBAL STRUCTURE */   
   gbl = ginit;
      
/*		THINGS NEEDED ON EACH HP_MGRID MESH FOR MGRID */
      FLT (*vug_frst)[NV];
      FLT (*sug_frst)[NV];
      FLT (*iug_frst)[NV];
      double (*dg)[NV],
