#include"mesh.h"
#include"utilities.h"
#include<cstdio>
#include<float.h>

/* CREATE SIDELIST FROM TRIANGLE VERTEX LIST */
/* USES VINFO TO STORE FIRST SIND FROM VERTEX */
/* USES SINFO TO STORE NEXT SIND FROM SIND */
/* TVRTX MUST BE COUNTERCLOCKWISE ORDERED */
void mesh::createsideinfo(void) {
   int i,j,tind,v1,v2,vout,temp,minv,maxv,order,sind,sindprev;
   
   for(i=0;i<nvrtx;++i)
      vinfo[i] = -1;
      
   nside = 0;
   for(tind=0;tind<ntri;++tind) {
      vout = tvrtx[tind][0];
      v1 = tvrtx[tind][1];
      v2 = tvrtx[tind][2];
      for(j=0;j<3;++j) {
         /* CHECK IF SIDE HAS BEEN CREATED ALREADY */
         if (v2 > v1) {
            minv = v1;
            maxv = v2;
            order = 0;
         }
         else {
            minv = v2;
            maxv = v1;
            order = 1;
         }
         sind = vinfo[minv];
         while (sind >= 0) {
            if (maxv == svrtx[sind][order]) {
               if (stri[sind][1] >= 0) {
                  printf("Error: side %d has been matched with Triangle %d 3 times\n",sind,tind);
                  exit(1);
               }
               else {
                  stri[sind][1] = tind;
                  tside[tind].side[j] = sind;
                  tside[tind].sign[j] = -1;
                  goto NEXTTRISIDE;
               }
            }
            sindprev = sind;
            sind = sinfo[sind];
         }
         /* NEW SIDE */
         svrtx[nside][0] = v1;
         svrtx[nside][1] = v2;
         stri[nside][0] = tind;
         stri[nside][1] = -1;
         tside[tind].side[j] = nside;
         tside[tind].sign[j] = 1;
         sinfo[nside] = -1;
         if (vinfo[minv] < 0)
            vinfo[minv] = nside;
         else 
            sinfo[sindprev] = nside;
         ++nside;
NEXTTRISIDE:
         temp = vout;
         vout = v1;
         v1 = v2;
         v2 = temp;
      }
   }

   return;
}

void mesh::createtsidestri(void) {
   int i,j,tind,v1,v2,vout,temp,minv,maxv,order,sind,sindprev;
   
   for(i=0;i<nvrtx;++i)
      vinfo[i] = -1;
      
   for(i=0;i<nside;++i) {
      v1 = svrtx[i][0];
      v2 = svrtx[i][1];
      minv = (v1 < v2 ? v1 : v2);
      sind = vinfo[minv];
      while (sind >= 0) {
         sindprev = sind;
         sind = sinfo[sind];
      }
      if (vinfo[minv] < 0)
         vinfo[minv] = i;
      else 
         sinfo[sindprev] = i;
      sinfo[i] = -1;
   }

   for(i=0;i<nside;++i)
      stri[i][1] = -1;

   for(tind=0;tind<ntri;++tind) {
      vout = tvrtx[tind][0];
      v1 = tvrtx[tind][1];
      v2 = tvrtx[tind][2];
      for(j=0;j<3;++j) {
         /* CHECK IF SIDE HAS BEEN CREATED ALREADY */
         if (v2 > v1) {
            minv = v1;
            maxv = v2;
            order = 0;
         }
         else {
            minv = v2;
            maxv = v1;
            order = 1;
         }
         sind = vinfo[minv];
         while (sind >= 0) {
            if (maxv == svrtx[sind][1]) {
               stri[sind][order] = tind;
               tside[tind].side[j] = sind;
               tside[tind].sign[j] = 1 -2*order;
               goto NEXTTRISIDE;
            }
            if (maxv == svrtx[sind][0]) {
               stri[sind][1-order] = tind;
               tside[tind].side[j] = sind;
               tside[tind].sign[j] = 2*order -1;
               goto NEXTTRISIDE;
            }
            sind = sinfo[sind];
         }
         printf("didn't match side %d %d\n",v1,v2);
         exit(1);
         
NEXTTRISIDE:
         temp = vout;
         vout = v1;
         v1 = v2;
         v2 = temp;
      }
   }

   return;
}



void mesh::createvtri(void) {
   int i,tind;
   
   /* THIS ALLOWS US TO GET TO LOCAL HIGHER ENTITIES FROM VERTEX NUMBER */
   for (tind=0;tind<ntri;++tind)
      for(i=0;i<3;++i)
         vtri[tvrtx[tind][i]] = tind;
   
   return;
}

/* CALCULATE NUMBER OF NEIGHBORS TO EACH CELL */
void mesh::cnt_nbor(void) {
   int i;
   
   for (i=0;i<nvrtx;++i)
      nnbor[i] = 0;
   
   for(i=0;i<nside;++i) {
      ++nnbor[svrtx[i][0]];
      ++nnbor[svrtx[i][1]];
   }

   return;
}


/* CREATES TRIANGLE TO TRIANGLE POINTER */
void mesh::createttri(void) {
   int tind,sind,j,flip;
   
   for(tind=0;tind<ntri;++tind) {
      for(j=0;j<3;++j) {
         sind = tside[tind].side[j];
         flip = (1 +tside[tind].sign[j])/2;
         ttri[tind][j] = stri[sind][flip];
      }
   }
   
   return;
}

void mesh::treeinit() {
   int i,j,sind,v0;
   FLT x1,y1,x2,y2;
   
   x1 = vrtx[0][0];
   y1 = vrtx[0][1];
   x2 = x1;
   y2 = y1;   
   
   for (i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nsd();++j) {
         sind = sbdry[i]->sd(j);
         v0 = svrtx[sind][0];
         x1 = MIN(x1,vrtx[v0][0]);
         y1 = MIN(y1,vrtx[v0][1]);
         x2 = MAX(x2,vrtx[v0][0]);
         y2 = MAX(y2,vrtx[v0][1]);
      }
   }

   qtree.init(x1-0.5*(x2-x1),y1-0.5*(y2-y1),x2+0.5*(x2-x1),y2+0.5*(y2-y1));
      
   for(i=0;i<nvrtx;++i) 
      qtree.addpt(i);

   return;
}

/* FIX STRI TTRI TO POINT TO GROUP/SIDE ON BOUNDARY */
void mesh::bdrylabel() {
   int i,j,k,sind,tind;
   
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nsd();++j) {
         sind = sbdry[i]->sd(j);
         stri[sind][1] = -(((i+1)<<16) +j);
         tind = stri[sind][0];
         for(k=0;k<3;++k)
            if (tside[tind].side[k] == sind) break;
            
         ttri[tind][k] = stri[sind][1];
      }
   }
   
   return;
}

void mesh::initvlngth() {
   int i,v0,v1;
   FLT l;
   
   for(i=0;i<nvrtx;++i)
      vlngth[i] = 0.0;
      
   for(i=0;i<nside;++i) {
      v0 = svrtx[i][0];
      v1 = svrtx[i][1];
      l = distance(svrtx[i][0],svrtx[i][1]);
      vlngth[v0] += l;
      vlngth[v1] += l;
   }
   
   for(i=0;i<nvrtx;++i)
      vlngth[i] /= nnbor[i];

   return;
}
      
      
