#include"mesh.h"
#include<cstdio>
#include<cfloat>
#include<iostream>
#include<stdlib.h>
#include<utilities.h>
#include<string.h>

void mesh::setfine(class mesh& tgt) {
   fmpt = &tgt;
   if (fine == NULL) 
      fine = new struct mg_trans[maxvst];
   mgconnect(fine,tgt);
   return;
}

void mesh::setcoarse(class mesh& tgt) {
   cmpt = &tgt;
   if (coarse == NULL)
      coarse = new struct mg_trans[maxvst];
   mgconnect(coarse,tgt);
   return;
}

/* THIS IS THE NEW WAY USING THE QUADTREE DATA STRUCTURE */
void mesh::mgconnect(struct mg_trans *cnnct, const class mesh& tgt) {
   int i,j,k,bnum,tind,sind,v0;
   double x,y,wgt[3],ainv;
   double dx,dy,a,b,c,minneg;
   int neg_count, triloc, sidloc;
   
   
   /* LOOP THROUGH VERTICES AND FIND SURROUNDING TRIANGLE */
   for(i=0;i<nvrtx;++i) {
      tgt.qtree.nearpt(vrtx[i],v0);
      cnnct[i].tri = tgt.findtri(vrtx[i][0],vrtx[i][1],v0);
      tgt.getwgts(cnnct[i].wt);
   }

   /* REDO BOUNDARY SIDES TO DEAL WITH CURVATURE */
   for(bnum=0;bnum<nsbd;++bnum) {
   
      /* CHECK TO MAKE SURE THESE ARE THE SAME SIDES */
      if(sbdry[bnum]->idnty() != tgt.sbdry[bnum]->idnty()) {
         printf("error: sides are not numbered the same\n");
         exit(1);
      }
         
      for(k=0;k<sbdry[bnum]->nsd();++k) {
         v0 = svrtx[sbdry[bnum]->sd(k)][0];
         x = vrtx[v0][0];
         y = vrtx[v0][1];
         minneg = -1.0E32;
         
         /* LOOP THROUGH TARGET SIDES TO FIND TRIANGLE */
         for(i=0;i<tgt.sbdry[bnum]->nsd();++i) {
            sind = tgt.sbdry[bnum]->sd(i);
            tind = tgt.stri[sind][0];
            if (tind < 0) {
               printf("boundary side in wrong direction %d %d\n",sind,tind);
               exit(1);
            }
            neg_count = 0;
            for(j=0;j<3;++j) {
               wgt[j] = 
               ((tgt.vrtx[tgt.tvrtx[tind][(j+1)%3]][0]
                  -tgt.vrtx[tgt.tvrtx[tind][j]][0])*
               (y-tgt.vrtx[tgt.tvrtx[tind][j]][1])-
               (tgt.vrtx[tgt.tvrtx[tind][(j+1)%3]][1]
               -tgt.vrtx[tgt.tvrtx[tind][j]][1])*
               (x-tgt.vrtx[tgt.tvrtx[tind][j]][0]));
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
         dx = x-tgt.vrtx[tgt.svrtx[sind][0]][0];
         dy = y-tgt.vrtx[tgt.svrtx[sind][0]][1];
         a = sqrt(dx*dx+dy*dy);
         dx = x-tgt.vrtx[tgt.svrtx[sind][1]][0];
         dy = y-tgt.vrtx[tgt.svrtx[sind][1]][1];
         b = sqrt(dx*dx+dy*dy);
         dx = tgt.vrtx[tgt.svrtx[sind][0]][0]
            -tgt.vrtx[tgt.svrtx[sind][1]][0];
         dy = tgt.vrtx[tgt.svrtx[sind][0]][1]
            -tgt.vrtx[tgt.svrtx[sind][1]][1];
         c = sqrt(dx*dx+dy*dy);
         a = (b*b+c*c-a*a)/(2.*c*c);
         x = tgt.vrtx[tgt.svrtx[sind][1]][0] +dx*a;
         y = tgt.vrtx[tgt.svrtx[sind][1]][1] +dy*a;
         tind = tgt.stri[sind][0];                     
         for(j=0;j<3;++j) {
            wgt[j] = 
            ((tgt.vrtx[tgt.tvrtx[tind][(j+1)%3]][0]
            -tgt.vrtx[tgt.tvrtx[tind][j]][0])*
            (y-tgt.vrtx[tgt.tvrtx[tind][j]][1])-
            (tgt.vrtx[tgt.tvrtx[tind][(j+1)%3]][1]
            -tgt.vrtx[tgt.tvrtx[tind][j]][1])*
            (x-tgt.vrtx[tgt.tvrtx[tind][j]][0]));
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
void mesh::testconnect(char *fname) {
   int i,j,n,tind;
   class mesh *fmesh, *cmesh;
   char fappnd[100];
   FLT (*work)[ND];
   
   work = (FLT (*)[ND]) xmalloc(ND*maxvst*sizeof(FLT));
   
   cmesh = static_cast<class mesh *>(cmpt);

   if (cmesh != NULL) {
            
      /* LOOP THROUGH VERTICES TO TO CALCULATE POSITION OF COARSE VERTICES  */
      for(i=0;i<nvrtx;++i) {
         tind = coarse[i].tri;
         
         for(n=0;n<ND;++n)
            work[i][n] = 0.0;
         
         for(j=0;j<3;++j) {
            for(n=0;n<ND;++n)
               work[i][n] += coarse[i].wt[j]*cmesh->vrtx[cmesh->tvrtx[tind][j]][n];
         }
      }
      strcpy(fappnd,fname);
      strcat(fappnd,"_tocrse");
      out_mesh(work, fappnd, grid);
   }
   
   
   fmesh = static_cast<class mesh *>(fmpt);

   if (fmesh != NULL) {
      /* LOOP THROUGH VERTICES   */
      /* TO CALCULATE VRTX POSITION FROM FINE MESH */
      for(i=0;i<nvrtx;++i) {
         tind = fine[i].tri;
   
         for(n=0;n<ND;++n)
            work[i][n] = 0.0;
            
         for(j=0;j<3;++j) {
            for(n=0;n<ND;++n)
               work[i][n] += fine[i].wt[j]*fmesh->vrtx[fmesh->tvrtx[tind][j]][n];
         }
      }
   
      strcpy(fappnd,fname);
      strcat(fappnd,"_tofine");
      out_mesh(work, fappnd, grid);
   }

   free(work);
   
   return;
}

   
