#include "mesh.h"
#include "boundary.h"
#include <cfloat>
#include <iostream>
#include <stdlib.h>
#include <utilities.h>
#include <string.h>

 void mesh::mgconnect(mesh::transfer *cnnct, const class mesh& tgt) {
   int i,j,k,bnum,tind,sind,v0;
   double x,y,wgt[3],ainv;
   double dx,dy,a,b,c,minneg;
   int neg_count, triloc, sidloc;
      
   /* LOOP THROUGH VERTICES AND FIND SURROUNDING TRIANGLE */
   for(i=0;i<nvrtx;++i) {
      tgt.qtree.nearpt(vrtx[i],v0);
      cnnct[i].tri = tgt.findtri(vrtx[i],v0);
      tgt.getwgts(cnnct[i].wt);
   }

   /* REDO BOUNDARY SIDES TO DEAL WITH CURVATURE */
   for(bnum=0;bnum<nsbd;++bnum) {
   
      /* CHECK TO MAKE SURE THESE ARE THE SAME SIDES */
      if(sbdry[bnum]->idnum != tgt.sbdry[bnum]->idnum) {
         *log << "error: sides are not numbered the same" << std::endl;
         exit(1);
      }
         
      for(k=0;k<sbdry[bnum]->nel;++k) {
         v0 = sd[sbdry[bnum]->el[k]].vrtx[0];
         x = vrtx[v0][0];
         y = vrtx[v0][1];
         minneg = -1.0E32;
         
         /* LOOP THROUGH TARGET SIDES TO FIND TRIANGLE */
         for(i=0;i<tgt.sbdry[bnum]->nel;++i) {
            sind = tgt.sbdry[bnum]->el[i];
            tind = tgt.sd[sind].tri[0];
            if (tind < 0) {
               *log << "boundary side in wrong direction" << sind << tind << std::endl;
               exit(1);
            }
            neg_count = 0;
            for(j=0;j<3;++j) {
               wgt[j] = 
               ((tgt.vrtx[tgt.td[tind].vrtx[(j+1)%3]][0]
                  -tgt.vrtx[tgt.td[tind].vrtx[j]][0])*
               (y-tgt.vrtx[tgt.td[tind].vrtx[j]][1])-
               (tgt.vrtx[tgt.td[tind].vrtx[(j+1)%3]][1]
               -tgt.vrtx[tgt.td[tind].vrtx[j]][1])*
               (x-tgt.vrtx[tgt.td[tind].vrtx[j]][0]));
               if (wgt[j] < 0.0) ++neg_count;
            }
            if (neg_count==0) {
               sidloc = sind;
               break;
            }
            else {
               if(neg_count == 1) {
                  for(j=0;j<3;++j) {
                     if (wgt[j] < 0 && wgt[j] > minneg) {
                        minneg = wgt[j];
                        sidloc = sind;
                        triloc = tind;
                     }
                  }
               }
            }
         }
         sind = sidloc;   

         /* PROJECT LOCATION NORMAL TO CURVED FACE */
         dx = x-tgt.vrtx[tgt.sd[sind].vrtx[0]][0];
         dy = y-tgt.vrtx[tgt.sd[sind].vrtx[0]][1];
         a = sqrt(dx*dx+dy*dy);
         dx = x-tgt.vrtx[tgt.sd[sind].vrtx[1]][0];
         dy = y-tgt.vrtx[tgt.sd[sind].vrtx[1]][1];
         b = sqrt(dx*dx+dy*dy);
         dx = tgt.vrtx[tgt.sd[sind].vrtx[0]][0]
            -tgt.vrtx[tgt.sd[sind].vrtx[1]][0];
         dy = tgt.vrtx[tgt.sd[sind].vrtx[0]][1]
            -tgt.vrtx[tgt.sd[sind].vrtx[1]][1];
         c = sqrt(dx*dx+dy*dy);
         a = (b*b+c*c-a*a)/(2.*c*c);
         x = tgt.vrtx[tgt.sd[sind].vrtx[1]][0] +dx*a;
         y = tgt.vrtx[tgt.sd[sind].vrtx[1]][1] +dy*a;
         tind = tgt.sd[sind].tri[0];                     
         for(j=0;j<3;++j) {
            wgt[j] = 
            ((tgt.vrtx[tgt.td[tind].vrtx[(j+1)%3]][0]
            -tgt.vrtx[tgt.td[tind].vrtx[j]][0])*
            (y-tgt.vrtx[tgt.td[tind].vrtx[j]][1])-
            (tgt.vrtx[tgt.td[tind].vrtx[(j+1)%3]][1]
            -tgt.vrtx[tgt.td[tind].vrtx[j]][1])*
            (x-tgt.vrtx[tgt.td[tind].vrtx[j]][0]));
         }
         cnnct[v0].tri = tind;
         ainv = 1.0/(tgt.area(tind));
         for (j=0;j<3;++j) {
            cnnct[v0].wt[(j+2)%3] = wgt[j]*ainv;
            if (wgt[j]*ainv > 1.0) 
               cnnct[v0].wt[(j+2)%3] = 1.0; 
               
            if (wgt[j]*ainv < 0.0)
               cnnct[v0].wt[(j+2)%3] = 0.0;
         }
      }
   }

   return;
}


/*	 THIS ROUTINE DETERMINES THE POSITION OF COARSE VERTICES  */
/*  TO TEST USING MULTI-GRID CONNECTION */
/* THE MULTIGRID CONNECTIONS */
 void mesh::testconnect(char *fname,transfer *cnnct, mesh *cmesh) {
   int i,j,n,tind;
   FLT (*work)[ND];
   
   work = (FLT (*)[ND]) xmalloc(ND*maxvst*sizeof(FLT));
   
   if (cmesh != NULL) {
            
      /* LOOP THROUGH VERTICES TO TO CALCULATE POSITION OF COARSE VERTICES  */
      for(i=0;i<nvrtx;++i) {
         tind = cnnct[i].tri;
         
         for(n=0;n<ND;++n)
            work[i][n] = 0.0;
         
         for(j=0;j<3;++j) {
            for(n=0;n<ND;++n)
               work[i][n] += cnnct[i].wt[j]*cmesh->vrtx[cmesh->td[tind].vrtx[j]][n];
         }
      }
      FLT (*storevrtx)[ND] = vrtx;
      vrtx = work;
      output(fname, ftype::grid);
      vrtx = storevrtx;
   }

   free(work);
   
   return;
}
   
