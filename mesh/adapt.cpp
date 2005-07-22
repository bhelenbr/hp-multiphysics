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

extern int nslst;  // (THIS NEEDS TO BE FIXED)

#define VDLTE 0x10
#define VSPEC 0x01
#define SDLTE 0x1000
#define STOUC 0x0100
#define TDLTE 0x100000
#define TTOUC 0x010000

void mesh::setup_for_adapt() {
   int i,v0,v1;
   
   /* For vertices 0 = untouched, 1 = special, 2 = deleted */
   /* For sides 0 = untouced, 1 = touched, 2 = deleted */
   /* For triangles 0 = untouched, 1 = touched, 2 = deleted */
   for(i=0;i<maxvst;++i)
      td(i).info = 0;

   /* List of entities needing refinement  */
   for(i=0;i<maxvst;++i)
      sd(i).info = -1;
      
   /* Back reference into list of entities needing refinement */
   for(i=0;i<maxvst;++i)
      vd(i).info = -1;
      
   /* MARK BEGINNING/END OF SIDE GROUPS & SPECIAL VERTEX POINTS */
   /* THESE SHOULD NOT BE DELETED */
   for(i=0;i<nsbd;++i) {
      v0 = sd(sbdry(i)->el(0)).vrtx(0);
      v1 = sd(sbdry(i)->el(sbdry(i)->nel-1).vrtx(1);
      td(v0).info |=  VSPEC;
      td(v1).info |=  VSPEC;
   }
   
   for(i=0;i<nvbd;++i)
      td(vbdry(i)->v0).info |= VSPEC;
            
   for(i=0;i<nside;++i)
      
         
   return;
}

/* NEXT STEP IS TO SWAP EDGES */

/* THEN COARSEN BOUNDARIES */
block::ctrl mesh::coarsen_bdry(int excpt,FLT tolsize) {
   int i,tind,sind,endpt;
   int ntdel,tdel[maxlst];
   int nsdel,sdel[maxlst];
   
   
   switch (excpt) {
      case(0):
         mp_phase = 0;
         
         for(int bnum=0;bnum<nsbd;++bnum) {
            if (!sbdry(bnum)->is_frst()) {
               sbdry(bnum)->master_slave_prepare();
               continue;
            }
            
            nslst = 0;
            for(int indx=0;indx<sbdry(bnum)->nel;++indx) {
               sind = sbdry(bnum)->el(indx);
               fscr1(sind) = side_lngth_ratio(sind);
               if (fscr1(sind) < tolsize) {
                  putinlst(sind);
               }
            }
            
            sbdry(bnum)->sndsize() = 0;
            while (nslst > 0) {
               // START WITH LARGEST SIDE LENGTH RATIO
               sind = i2wk(nslst-1);
               el = getbdryel(sd(sind).tri(1));
               
               /* COLLAPSE EDGE */
               /* ASSUME COLLAPSE WILL MARK SIDES AS DELETED OR AFFECTED */
               endpt = collapse(sind,ntdel,tdel,nsdel,sdel);
               
               tkoutlst(sind);
               if (endpt < 0) {
                  *log << "#Warning: boundary side collapse failed" << sind << std::endl;
                  continue;
               }
               /* STORE ADAPTATION INFO FOR COMMUNICATION BOUNDARIES */
               sbdry(bnum)->isndbuf(sbdry(bnum)->sndsize()++) = el;
               sbdry(bnum)->isndbuf(sbdry(bnum)->sndsize()++) = endpt;
               
               /* RECACULCULATE SIDE RATIO FOR AFFECTED SIDES */
               /* FIND ADJACENT NON-DELETED SIDE */
               if (endpt > 0) {
                  for(next = el+1;next<sbdry(bnum)->nel;++next) {
                     if(!(td(sbdry(bnum)->el(next)).info&SDLTE)) goto found;
                  }
                  continue;
               }
               else {
                  for(next = el-1;next>=0;--next) {
                     if(!(td(sbdry(bnum)->el(next)).info&SDLTE)) goto found;
                  }
                  continue;
               }
               
               found: 
               sind = sbdry(bnum)->el(next);
               if (i3wk(sind) > -1) tkoutlst(sind);
               fscr1(sind) = side_lngth_ratio(sind);
               if (fscr1(sind) > tolsize) putinlst(sind);
            }
         }
         return(block::advance);

      case(1):
         for(int i=0;i<nsbd;++i) 
            sbdry(i)->master_slave_transmit();
         
         return(block::advance);
         
      case(2):
         /* RECEIVE ADAPTATION INFO FOR COMM BOUNDARIES AND ADAPT SECOND SIDES */
         for(int i=0;i<nsbd;++i) 
            sbdry(i)->master_slave_wait();
            
         int j,k,m,count,offset,sind;
#ifdef MPISRC
         MPI_Status status;
#endif

         for(int bnum=0;bnum<nsbd;++bnum) {
            
            if (sbdry(bnum)->is_frst()) continue;
            
            sbdry(bnum)->master_slave_wait();
            
            
            
            /* ASSUMES REVERSE ORDERING OF SIDES */
            /* WON'T WORK IN 3D */
            
            
            
            /* RELOAD FROM BUFFER */
            /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
            /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */
            count = 0;
            for(j=0;j<nel;++j) {
               sind = el(j);
               offset = x.sd(sind).vrtx(0)*stride;
               for (k=bgn;k<=end;++k)
                  base[offset+k] = fsndbuf(count++);
            }
            offset = x.sd(sind).vrtx(1)*stride;
            for (k=bgn;k<=end;++k)
               base[offset+k] = fsndbuf(count++);
            
            for(m=0;m<nlocal_match;++m) {            
               count = 0;
               for(j=nel-1;j>=0;--j) {
                  sind = el(j);
                  offset = x.sd(sind).vrtx(1)*stride;
                  for (k=bgn;k<=end;++k) 
                     base[offset+k] += static_cast<FLT *>(local_rcv_buf[m])[count++];
               }
               offset = x.sd(sind).vrtx(0)*stride;
               for (k=bgn;k<=end;++k) 
                     base[offset+k] += static_cast<FLT *>(local_rcv_buf[m])[count++];            
            }
            

         #ifdef MPISRC
            /* MPI PASSES */
            for(m=0;m<nmpi_match;++m) {
               MPI_Wait(&mpi_rcvrqst[m], &status);
               count = 0;
               for(j=nel-1;j>=0;--j) {
                  sind = el(j);
                  offset = x.sd(sind).vrtx(1)*stride;
                  for (k=bgn;k<=end;++k) 
                     base[offset+k] += static_cast<FLT *>(mpi_rcv_buf[m])[count++];
               }
               offset = x.sd(sind).vrtx(0)*stride;
               for (k=bgn;k<=end;++k) 
                     base[offset+k] += static_cast<FLT *>(mpi_rcv_buf[m])[count++];
            }
         #endif

            count = 0;
            for(j=0;j<nel;++j) {
               sind = el(j);
               offset = x.sd(sind).vrtx(0)*stride;
               for (k=bgn;k<=end;++k)
                  base[offset+k] /= (1 +nlocal_match +nmpi_match);
            }
            offset = x.sd(sind).vrtx(1)*stride;
            for (k=bgn;k<=end;++k)
               base[offset+k] /= (1 +nlocal_match +nmpi_match);
         }

   return;
}
      
      

   
   


      

