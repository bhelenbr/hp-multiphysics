#include"spectral_hp.h"
#include<assert.h>

extern FLT hgt(int type, FLT x, FLT y);
extern FLT dhgtdx(int type, FLT x, FLT y);
extern FLT dhgtdy(int type, FLT x, FLT y);

void spectral_hp::ptprobe(FLT xp, FLT yp, FLT uout[NV]) {
   FLT r,s;
   int tind;
   
   tind = findinteriorpt(xp,yp,r,s);
   ugtouht(tind);  
   b.ptprobe(NV,uht,uout,r,s);
}

void spectral_hp::ptprobe1d(int typ, FLT xp, FLT yp, FLT uout[NV]) {
   FLT psi;
   int sind;
   
   sind = findbdrypt(typ,xp,yp,psi);
   ugtouht1d(sind);  
   b.ptprobe1d(NV,uht,uout,psi);

}

int spectral_hp::findandmvptincurved(FLT &xp, FLT &yp, FLT &r, FLT &s) {
   FLT x[ND],wgt[3];
   int tind,v0;

   x[0] = xp;
   x[1] = yp;
   qtree.nearpt(x,v0);
   tind = findtri(xp,yp,v0);
   getwgts(wgt);

   assert(tind > -1);

   /* TRIANGLE COORDINATES */   
   s = wgt[0]*2 -1.0;
   r = wgt[2]*2 -1.0;
   
   if (tinfo[tind] < 0) return(tind);
   
   printf("#In find and move %f %f\n",xp,yp);

   /* MOVE POINT WITH SIDE CURVATURE */
   crdtocht(tind);
   b.ptprobe_bdry(ND,cht,x,r,s);
   xp = x[0];
   yp = x[1];


   return(tind);
}


int spectral_hp::findinteriorpt(FLT xp, FLT yp, FLT &r, FLT &s) {
   FLT dr,ds,dx,dy,x[ND],ddr[3],dds[3],wgt[3],det,roundoff,dxmax[ND];
   int n,iter,sind,tind,v0;

   x[0] = xp;
   x[1] = yp;
   qtree.nearpt(x,v0);
   tind = findtri(xp,yp,v0);
   getwgts(wgt);

   if (tind == -1) {
      /* POINT IS IN CONVEX CURVED TRIANGLE NEAR BOUNDARY */
      /* THIS ASSUMES SIDE RADIUS OF CURVATURE IS LESS THAN 1/2 LENGTH OF SIDE ON BDRY */
      x[0] = xp;
      x[1] = yp;
      sind = findbdryside(x,v0,CURV_MASK,CURV_MASK);
      if (sind < 0) {
         printf("Warning: couldn't find boundary tri (%f,%f) nearpt %d neartri %d sind %d\n",xp,yp,v0,vtri[v0],sind);
      }
      tind = stri[sind][0];
      wgt[2] = 0.5;
      wgt[0] = 0.5;
   }


   /* TRIANGLE COORDINATES */   
   s = wgt[0]*2 -1.0;
   r = wgt[2]*2 -1.0;
   
   if (tinfo[tind] < 0) return(tind);

   /* DEAL WITH CURVED SIDES */
   crdtocht(tind);
   
   for(n=0;n<ND;++n)
      dxmax[n] = fabs(cht[n][0]-cht[n][1]) +fabs(cht[n][1]-cht[n][2]);
   roundoff = 10.0*EPSILON*(1.0 +(fabs(xp)*dxmax[1] +fabs(yp)*dxmax[0])/area(tind));
   
   iter = 0;
   do {
      b.ptprobe_bdry(ND,cht,x,ddr,dds,r,s);
      det = 1.0/(fabs(ddr[0]*dds[1] - ddr[1]*dds[0]) +10.0*EPSILON);
      dx = xp-x[0];
      dy = yp-x[1];
      dr =  (dds[1]*dx -dds[0]*dy)*det;
      ds = -(ddr[1]*dx -ddr[0]*dy)*det;

      r += dr;
      s += ds;
      if (iter++ > 100) {
         printf("#Warning: max iterations for curved triangle %d loc: %f,%f (r,s) %f,%f error: %e\n",tind,xp,yp,r,s,fabs(dr) +fabs(ds));
         break;
      }
   } while (fabs(dr) +fabs(ds) > roundoff);

   return(tind);
}
      
int spectral_hp::findbdrypt(int typ, FLT &x, FLT &y, FLT &psi) {
   int vnear,sind,tind,snum,snumnew,v0,v1,iter,bnum;
   FLT dpsi,xp[ND],dx,dy,ol,roundoff;
   
   /* SEARCH FOR TRI ADJACENT TO BOUNDARY NEAR POINT */
   xp[0] = x; xp[1] = y;
   qtree.nearpt(xp,vnear);
   sind = findbdryside(xp,vnear,typ);
   
   if (sind > -1) {
      bnum = (-stri[sind][1]>>16) -1;
      snum = -stri[sind][1]&0xFFFF;
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];
      dx = vrtx[v1][0] - vrtx[v0][0];
      dy = vrtx[v1][1] - vrtx[v0][1];
      ol = 2./(dx*dx +dy*dy);
      psi = ol*((x -vrtx[v0][0])*dx +(y -vrtx[v0][1])*dy) -1.;
   }
   else {
      printf("#Warning: brute force boundary locate (%f,%f) nearpt %d neartri %d type %d\n",x,y,vnear,vtri[vnear],typ);
      for(bnum=0;bnum<nsbd;++bnum)
         if (sbdry[bnum].type == typ) break;
      
      /* SEARCH FROM MIDDLE */
      sind = sbdry[bnum].el[sbdry[bnum].num/2];
      snumnew = -stri[sind][1]&0xFFFF;
      snum = snumnew;
      for(;;) {
         sind = sbdry[bnum].el[snumnew];
         v0 = svrtx[sind][0];
         v1 = svrtx[sind][1];
         dx = vrtx[v1][0] - vrtx[v0][0];
         dy = vrtx[v1][1] - vrtx[v0][1];
         ol = 2./(dx*dx +dy*dy);
         psi = ol*((x -vrtx[v0][0])*dx +(y -vrtx[v0][1])*dy) -1.;
   
         if (psi < -1.) {
            if (snumnew <= snum) {
               snum = snumnew;
               snumnew -= 1;
               if (snumnew == -1) snumnew = sbdry[bnum].num-1;
               continue;
            }
            else {
               printf("#Warning: trouble finding side point %f %f, best guess side %d psi -1.0\n",x,y,sbdry[bnum].el[snumnew]);
               psi = -1.;
            }
         }
         else if (psi >  1.) {
            if (snumnew >= snum) {
               snum = snumnew;
               snumnew += 1;
               if (snumnew == sbdry[bnum].num) snumnew = 0;
               continue;
            }
            else {
               printf("#Warning: trouble finding side point %f %f, best guess side %d psi 1.0\n",x,y,sbdry[bnum].el[snumnew]);
               psi = 1.;
            }
         }
         break;
      }
   }

   if (!(typ&CURV_MASK)) {
      x = vrtx[v0][0] +dx*(psi +1.)*.5;
      y = vrtx[v0][1] +dy*(psi +1.)*.5;
   }
   else {
      /* FIND PSI SUCH THAT TANGENTIAL POSITION ALONG LINEAR SIDE STAYS THE SAME */
      /* THIS WAY, MULTIPLE CALLS WILL NOT GIVE DIFFERENT RESULTS */ 
      iter = 0;
      crdtocht1d(sind);
      dx *= ol;
      dy *= ol;
      roundoff = 10.0*EPSILON*(1.0 +(fabs(x*dx) +fabs(y*dy)));
      do {
         b.ptprobe1d(ND,cht,xp,psi);
         dpsi = (x -xp[0])*dx +(y -xp[1])*dy;
         psi += dpsi;
         if (iter++ > 100) {
            printf("#Warning: max iterations for curved triangle in bdry_locate tri: %d type: %d loc: %f,%f error: %e\n",tind,typ,x,y,dpsi);
            break;
         }  
      } while (fabs(dpsi) > roundoff);
      x = xp[0];
      y = xp[1]; 
   }
   
   return(sind);
}

