#include"mesh.h"
#include"utilities.h"
#include<float.h>

/* CREATE SIDELIST FROM TRIANGLE VERTEX LIST */
/* USES VINFO TO STORE FIRST SIND FROM VERTEX */
/* USES SINFO TO STORE NEXT SIND FROM SIND */
/* TVRTX MUST BE COUNTERCLOCKWISE ORDERED */
template void mesh<2>::createsideinfo(void);
template void mesh<3>::createsideinfo(void);

template<int ND> void mesh<ND>::createsideinfo(void) {
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
                  *log << "Error: side " << sind << "has been matched with Triangle" << tind << "3 times" << std::endl;                  exit(1);
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

template void mesh<2>::createtsidestri(void);
template void mesh<3>::createtsidestri(void);

template<int ND> void mesh<ND>::createtsidestri(void) {
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


template void mesh<2>::createvtri(void);
template void mesh<3>::createvtri(void);

template<int ND> void mesh<ND>::createvtri(void) {
   int i,tind;
   
   /* THIS ALLOWS US TO GET TO LOCAL HIGHER ENTITIES FROM VERTEX NUMBER */
   for (tind=0;tind<ntri;++tind)
      for(i=0;i<3;++i)
         vtri[tvrtx[tind][i]] = tind;
   
   return;
}

/* CALCULATE NUMBER OF NEIGHBORS TO EACH CELL */
template void mesh<2>::cnt_nbor(void);
template void mesh<3>::cnt_nbor(void);

template<int ND> void mesh<ND>::cnt_nbor(void) {
   int i;
   
   for (i=0;i<nvrtx;++i)
      nnbor[i] = 0;
   
   for(i=0;i<nside;++i) {
      ++nnbor[svrtx[i][0]];
      ++nnbor[svrtx[i][1]];
   }

   return;
}

template void mesh<2>::createttri(void);
template void mesh<3>::createttri(void);

/* CREATES TRIANGLE TO TRIANGLE POINTER */
template<int ND> void mesh<ND>::createttri(void) {
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

template void mesh<2>::treeinit();
template void mesh<3>::treeinit();

template<int ND> void mesh<ND>::treeinit() {
   int i,j,n,sind,v0;
   FLT x1[ND], x2[ND], dx;
   
   for(n=0;n<ND;++n)	{
      x1[n] = vrtx[0][n];
      x2[n] = vrtx[0][n];
   }

   
   for (i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nsd();++j) {
         sind = sbdry[i]->sd(j);
         v0 = svrtx[sind][0];
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
      

   qtree.init(x1,x2);
      
   for(i=0;i<nvrtx;++i) 
      qtree.addpt(i);

   return;
}

template void mesh<2>::bdrylabel();
template void mesh<3>::bdrylabel();

/* FIX STRI TTRI TO POINT TO GROUP/SIDE ON BOUNDARY */
template<int ND> void mesh<ND>::bdrylabel() {
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

template void mesh<2>::initvlngth();
template void mesh<3>::initvlngth();

template<int ND> void mesh<ND>::initvlngth() {
   int i,j,v0,v1;
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
      
   
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nsd();++j) {
         v0 = svrtx[sbdry[i]->sd(j)][0];
         vlngth[v0] = 1.0e32;
      }
   }
           
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nsd();++j) {
         v0 = svrtx[sbdry[i]->sd(j)][0];
         v1 = svrtx[sbdry[i]->sd(j)][1];
         l = distance(v0,v1);
         vlngth[v0] = MIN(l,vlngth[v0]);
         vlngth[v1] = MIN(l,vlngth[v1]);
      }
   }

   return;
}

template void mesh<2>::settrim();
template void mesh<3>::settrim();

template<int ND> void mesh<ND>::settrim() {
   int i,j,n,bsd,tin,tind,nsrch;
   
   for(i=0;i<ntri;++i)
      tinfo[i] = 0;

   ntdel = 0;

   for (bsd=0;bsd<sbdry[0]->nsd();++bsd) {
      tind = stri[sbdry[0]->sd(bsd)][0];
      if (tinfo[tind] > 0) continue;
      
      intwk1[0] = tind;
      tinfo[tind] = 1;
      nsrch = ntdel+1;
      
      /* NEED TO SEARCH SURROUNDING TRIANGLES */
      for(i=ntdel;i<nsrch;++i) {
         tin = intwk1[i];
         for (n=0;n<3;++n)
            if (fltwk[tvrtx[tin][n]] < 0.0) goto NEXT;
            
         intwk1[ntdel++] = tin;

         for(j=0;j<3;++j) {
            tind = ttri[tin][j];
            if (tind < 0) continue;
            if (tinfo[tind] > 0) continue; 
            tinfo[tind] = 1;        
            intwk1[nsrch++] = tind;
         }
         NEXT: continue;
      }
   }
   
   for(i=0;i<ntri;++i)
      tinfo[i] = 0;
      
   for(i=0;i<ntdel;++i)
      tinfo[intwk1[i]] = 1;
      
   for(i=0;i<maxvst;++i)
      intwk1[i] = -1;

   return;
}
