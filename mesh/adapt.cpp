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

void mesh::setup_for_adapt() {
   int i,v0,v1;
   
   /* NEED TO INITIALIZE TO ZERO TO KEEP TRACK OF DELETED TRIS (-1) */
   /* ALSO TO DETERMINE TRI'S ON BOUNDARY OF COARSENING REGION */
   for(i=0;i<ntri;++i)
      td[i].info = 0;

   /* KEEPS TRACK OF DELETED SIDES = -3, TOUCHED SIDES =-2, UNTOUCHED SIDES =-1 */
   for(i=0;i<nside;++i)
      sd[i].info = -1;
      
   /* VINFO TO KEEP TRACK OF SPECIAL VERTICES (1) DELETED VERTICES (-1) */
   for(i=0;i<nvrtx;++i)
      vd[i].info = 0;
      
   /* MARK BEGINNING/END OF SIDE GROUPS & SPECIAL VERTEX POINTS */
   /* THESE SHOULD NOT BE DELETED */
   for(i=0;i<nsbd;++i) {
      v0 = sd[sbdry[i]->el[0]].vrtx[0];
      v1 = sd[sbdry[i]->el[sbdry[i]->nel-1]].vrtx[1];
      vd[v0].info = 1;
      vd[v1].info = 1;
   }
   
   for(i=0;i<nvbd;++i)
      vd[vbdry[i]->v0].info = 1;
      
   for(i=0;i<nside;++i)
      fwk[i] = side_lngth_ratio(i);
         
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
            
            if (!sbdry[bnum]->is_frst()) {
               sbdry[bnum]->master_slave_prepare();
               continue;
            }
            
            nslst = 0;
            for(int indx=0;indx<sbdry[bnum]->nel;++indx) {
               sind = sbdry[bnum]->el[indx];
               if (fwk[sind] < tolsize) {
                  putinlst(sind);
               }
            }
            
            sbdry[bnum]->sndsize() = 0;
            while (nslst > 0) {
               // START WITH LARGEST SIDE TO DENSITY RATIO
               sind = i2wk[nslst-1];
               
               /* COLLAPSE EDGE */
               endpt = collapse(sind,ntdel,tdel,nsdel,sdel);
               
               tkoutlst(sind);
               if (endpt < 0) {
                  *log << "#Warning: boundary side collapse failed" << sind << std::endl;
                  continue;
               }
               /* STORE ADAPTATION INFO FOR COMMUNICATION BOUNDARIES */
               sbdry[bnum]->isndbuf(sbdry[bnum]->sndsize()++) = sind;
               sbdry[bnum]->isndbuf(sbdry[bnum]->sndsize()++) = endpt;
               
               /* LOOP THROUGH AFFECTED TRIANGLES */
               for(int i=0;i<ntdel;++i) {
                  tind = tdel[i];
                  for(int j=0;j<3;++j) {
                     sind = td[tind].side[j];
                     if (i3wk[sind] > -1) tkoutlst(sind);
                     fwk[sind] = side_lngth_ratio(sind);
                     if (fwk[sind] > tolsize && (-sd[sind].tri[1])>>16 == bnum+1) {
                        putinlst(sind);
                     }
                  }
               }
            }
         }
         return(block::advance);

      case(1):
         for(int i=0;i<nsbd;++i) 
            sbdry[i]->master_slave_transmit();
         
         return(block::advance);
         
      case(2):
         /* RECEIVE ADAPTATION INFO FOR COMM BOUNDARIES AND ADAPT SECOND SIDES */
         for(int i=0;i<nsbd;++i)
            sbdry[i]->master_slave_wait();
            
         int j,k,m,count,offset,sind;
#ifdef MPISRC
         MPI_Status status;
#endif

         for(int bnum=0;bnum<nsbd;++bnum) {
            
            if (sbdry[bnum]->is_frst()) continue;
            /* ASSUMES REVERSE ORDERING OF SIDES */
            /* WON'T WORK IN 3D */
            
            /* RELOAD FROM BUFFER */
            /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
            /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */
            count = 0;
            for(j=0;j<nel;++j) {
               sind = el[j];
               offset = x.sd[sind].vrtx[0]*stride;
               for (k=bgn;k<=end;++k)
                  base[offset+k] = fsndbuf(count++);
            }
            offset = x.sd[sind].vrtx[1]*stride;
            for (k=bgn;k<=end;++k)
               base[offset+k] = fsndbuf(count++);
            
            for(m=0;m<nlocal_match;++m) {            
               count = 0;
               for(j=nel-1;j>=0;--j) {
                  sind = el[j];
                  offset = x.sd[sind].vrtx[1]*stride;
                  for (k=bgn;k<=end;++k) 
                     base[offset+k] += static_cast<FLT *>(local_rcv_buf[m])[count++];
               }
               offset = x.sd[sind].vrtx[0]*stride;
               for (k=bgn;k<=end;++k) 
                     base[offset+k] += static_cast<FLT *>(local_rcv_buf[m])[count++];            
            }
            

         #ifdef MPISRC
            /* MPI PASSES */
            for(m=0;m<nmpi_match;++m) {
               MPI_Wait(&mpi_rcvrqst[m], &status);
               count = 0;
               for(j=nel-1;j>=0;--j) {
                  sind = el[j];
                  offset = x.sd[sind].vrtx[1]*stride;
                  for (k=bgn;k<=end;++k) 
                     base[offset+k] += static_cast<FLT *>(mpi_rcv_buf[m])[count++];
               }
               offset = x.sd[sind].vrtx[0]*stride;
               for (k=bgn;k<=end;++k) 
                     base[offset+k] += static_cast<FLT *>(mpi_rcv_buf[m])[count++];
            }
         #endif

            count = 0;
            for(j=0;j<nel;++j) {
               sind = el[j];
               offset = x.sd[sind].vrtx[0]*stride;
               for (k=bgn;k<=end;++k)
                  base[offset+k] /= (1 +nlocal_match +nmpi_match);
            }
            offset = x.sd[sind].vrtx[1]*stride;
            for (k=bgn;k<=end;++k)
               base[offset+k] /= (1 +nlocal_match +nmpi_match);
         }

      
   return;
}
      
      

   
   


      

