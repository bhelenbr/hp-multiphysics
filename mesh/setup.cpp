#include "mesh.h"
#include "boundary.h"
#include "utilities.h"
#include <float.h>

/* CREATE SIDELIST FROM TRIANGLE VERTEX LIST */
/* USES VINFO TO STORE FIRST SIND FROM VERTEX */
/* USES SINFO TO STORE NEXT SIND FROM SIND */
/* TVRTX MUST BE COUNTERCLOCKWISE ORDERED */
void mesh::createsideinfo(void) {
   int i,j,tind,v1,v2,vout,temp,minv,maxv,order,sind,sindprev;
   
   for(i=0;i<nvrtx;++i)
      vd[i].info = -1;
      
   nside = 0;
   for(tind=0;tind<ntri;++tind) {
      vout = td[tind].vrtx[0];
      v1 = td[tind].vrtx[1];
      v2 = td[tind].vrtx[2];
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
         sind = vd[minv].info;
         while (sind >= 0) {
            if (maxv == sd[sind].vrtx[order]) {
               if (sd[sind].tri[1] >= 0) {
                  *log << "Error: side " << sind << "has been matched with Triangle" << tind << "3 times" << std::endl;                  exit(1);
               }
               else {
                  sd[sind].tri[1] = tind;
                  td[tind].side[j] = sind;
                  td[tind].sign[j] = -1;
                  goto NEXTTRISIDE;
               }
            }
            sindprev = sind;
            sind = sd[sind].info;
         }
         /* NEW SIDE */
         sd[nside].vrtx[0] = v1;
         sd[nside].vrtx[1] = v2;
         sd[nside].tri[0] = tind;
         sd[nside].tri[1] = -1;
         td[tind].side[j] = nside;
         td[tind].sign[j] = 1;
         sd[nside].info = -1;
         if (vd[minv].info < 0)
            vd[minv].info = nside;
         else 
            sd[sindprev].info = nside;
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

void mesh::createtdstri(void) {
   int i,j,tind,v1,v2,vout,temp,minv,maxv,order,sind,sindprev;
   
   for(i=0;i<nvrtx;++i)
      vd[i].info = -1;
      
   for(i=0;i<nside;++i) {
      v1 = sd[i].vrtx[0];
      v2 = sd[i].vrtx[1];
      minv = (v1 < v2 ? v1 : v2);
      sind = vd[minv].info;
      while (sind >= 0) {
         sindprev = sind;
         sind = sd[sind].info;
      }
      if (vd[minv].info < 0)
         vd[minv].info = i;
      else 
         sd[sindprev].info = i;
      sd[i].info = -1;
   }

   for(i=0;i<nside;++i)
      sd[i].tri[1] = -1;

   for(tind=0;tind<ntri;++tind) {
      vout = td[tind].vrtx[0];
      v1 = td[tind].vrtx[1];
      v2 = td[tind].vrtx[2];
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
         sind = vd[minv].info;
         while (sind >= 0) {
            if (maxv == sd[sind].vrtx[1]) {
               sd[sind].tri[order] = tind;
               td[tind].side[j] = sind;
               td[tind].sign[j] = 1 -2*order;
               goto NEXTTRISIDE;
            }
            if (maxv == sd[sind].vrtx[0]) {
               sd[sind].tri[1-order] = tind;
               td[tind].side[j] = sind;
               td[tind].sign[j] = 2*order -1;
               goto NEXTTRISIDE;
            }
            sind = sd[sind].info;
         }
         *log << "didn't match side: " << v1 << v2 << std::endl;
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
         vd[td[tind].vrtx[i]].tri = tind;
   
   return;
}

/* CALCULATE NUMBER OF NEIGHBORS TO EACH CELL */
void mesh::cnt_nbor(void) {
   int i;
   
   for (i=0;i<nvrtx;++i)
      vd[i].nnbor = 0;
   
   for(i=0;i<nside;++i) {
      ++vd[sd[i].vrtx[0]].nnbor;
      ++vd[sd[i].vrtx[1]].nnbor;
   }

   return;
}

/* CREATES TRIANGLE TO TRIANGLE POINTER */
void mesh::createttri(void) {
   int tind,sind,j,flip;
   
   for(tind=0;tind<ntri;++tind) {
      for(j=0;j<3;++j) {
         sind = td[tind].side[j];
         flip = (1 +td[tind].sign[j])/2;
         td[tind].tri[j] = sd[sind].tri[flip];
      }
   }
   
   return;
}

void mesh::treeinit() {
   int i,j,n,sind,v0;
   FLT x1[ND], x2[ND], dx;
   
   for(n=0;n<ND;++n)	{
      x1[n] = vrtx[0][n];
      x2[n] = vrtx[0][n];
   }

   
   for (i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nel;++j) {
         sind = sbdry[i]->el[j];
         v0 = sd[sind].vrtx[0];
         for(n=0;n<ND;++n) {
            x1[n] = MIN(x1[n],vrtx[v0][n]);
            x2[n] = MAX(x2[n],vrtx[v0][n]);
         }
      }
   }
   
   for(n=0;n<ND;++n) {
      dx = MAX(x2[n]-x1[n],100.0*EPSILON);
      x1[n] -= 0.25*dx;
      x2[n] += 0.25*dx;
   }
      
   // *log << nsbd << "max:" << maxvst << std::endl;
   // *log << x1[0] << "," << x1[1] << std::endl;
   // *log << x2[0] << "," << x2[1] << std::endl;
   qtree.init(x1,x2);
      
   for(i=0;i<nvrtx;++i) 
      qtree.addpt(i);

   return;
}

/* FIX STRI TTRI TO POINT TO GROUP/SIDE ON BOUNDARY */
void mesh::bdrylabel() {
   int i,j,k,sind,tind;
   
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nel;++j) {
         sind = sbdry[i]->el[j];
         sd[sind].tri[1] = -(((i+1)<<16) +j);
         tind = sd[sind].tri[0];
         for(k=0;k<3;++k)
            if (td[tind].side[k] == sind) break;
            
         td[tind].tri[k] = sd[sind].tri[1];
      }
   }
   
   return;
}

void mesh::initvlngth() {
   int i,j,v0,v1;
   FLT l;
   
   for(i=0;i<nvrtx;++i)
      vlngth[i] = 0.0;
      
   for(i=0;i<nside;++i) {
      v0 = sd[i].vrtx[0];
      v1 = sd[i].vrtx[1];
      l = distance(sd[i].vrtx[0],sd[i].vrtx[1]);
      vlngth[v0] += l;
      vlngth[v1] += l;
   }
   
   for(i=0;i<nvrtx;++i)
      vlngth[i] /= vd[i].nnbor;
      
   
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nel;++j) {
         v0 = sd[sbdry[i]->el[j]].vrtx[0];
         vlngth[v0] = 1.0e32;
      }
   }
           
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nel;++j) {
         v0 = sd[sbdry[i]->el[j]].vrtx[0];
         v1 = sd[sbdry[i]->el[j]].vrtx[1];
         l = distance(v0,v1);
         vlngth[v0] = MIN(l,vlngth[v0]);
         vlngth[v1] = MIN(l,vlngth[v1]);
      }
   }

   return;
}

void mesh::settrim() {
   int i,j,n,bsd,tin,tind,nsrch,ntdel;
   
   for(i=0;i<ntri;++i)
      td[i].info = 0;

   ntdel = 0;

   for (bsd=0;bsd<sbdry[0]->nel;++bsd) {
      tind = sd[sbdry[0]->el[bsd]].tri[0];
      if (td[tind].info > 0) continue;
      
      i1wk[0] = tind;
      td[tind].info = 1;
      nsrch = ntdel+1;
      
      /* NEED TO SEARCH SURROUNDING TRIANGLES */
      for(i=ntdel;i<nsrch;++i) {
         tin = i1wk[i];
         for (n=0;n<3;++n)
            if (fwk[td[tin].vrtx[n]] < 0.0) goto NEXT;
            
         i1wk[ntdel++] = tin;

         for(j=0;j<3;++j) {
            tind = td[tin].tri[j];
            if (tind < 0) continue;
            if (td[tind].info > 0) continue; 
            td[tind].info = 1;        
            i1wk[nsrch++] = tind;
         }
         NEXT: continue;
      }
   }
   
   for(i=0;i<ntri;++i)
      td[i].info = 0;
      
   for(i=0;i<ntdel;++i)
      td[i1wk[i]].info = 1;
      
   for(i=0;i<maxvst;++i)
      i1wk[i] = -1;

   return;
}
