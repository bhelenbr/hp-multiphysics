#include"spectral_hp.h"
#include<assert.h>

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

int spectral_hp::findinteriorpt(FLT xp, FLT yp, FLT &r, FLT &s) {
   FLT dr,ds,dx,dy,x[ND],ddr[3],dds[3],wgt[3],det;
   int iter,sind,tind,v0;

   qtree.nearpt(xp,yp,v0);
   tind = findtri(xp,yp,v0);
   getwgts(wgt);

   if (tind == -1) {
      /* POINT IS IN CONVEX CURVED TRIANGLE NEAR BOUNDARY */
      /* THIS ASSUMES SIDE RADIUS OF CURVATURE IS LESS THAN 1/2 LENGTH OF SIDE ON BDRY */
      x[0] = xp;
      x[1] = yp;
      sind = findbdryside(x,v0);
      if (sind < 0 || sinfo[sind] < 0) {
         printf("Warning: error finding boundary tri (%f,%f) nearpt %d sind %d\n",xp,yp,v0,sind);
      }
      tind = stri[sind][0];
      wgt[1] = 0.5;
      wgt[2] = 0.5;
   }


   /* TRIANGLE COORDINATES */   
   s = wgt[2]*2 -1.0;
   r = wgt[1]*2 -1.0;
   if (tinfo[tind] < 0) return(tind);
   

   /* DEAL WITH CURVED SIDES */
   crdtouht(tind);

   iter = 0;
   do {
      b.ptprobe_bdry(ND,uht,x,ddr,dds,r,s);
      det = 1.0/(fabs(ddr[0]*dds[1] - ddr[1]*dds[0]) +10.0*EPSILON);
      dx = xp-x[0];
      dy = yp-x[1];
      dr =  (dds[1]*dx -dds[0]*dy)*det;
      ds = -(ddr[1]*dx -ddr[0]*dy)*det;

      r += dr;
      s += ds;
      if (iter++ > 50) {
         printf("#Warning: max iterations for curved triangle %d loc: %f,%f (r,s) %f,%f error: %f\n",tind,xp,yp,r,s,fabs(dr) +fabs(ds));
         break;
      }
   } while (fabs(dr) +fabs(ds) > 100.*EPSILON);

   return(tind);
}
      
int spectral_hp::findbdrypt(int typ, FLT &x, FLT &y, FLT &psi) {
   int j,vnear,sind,tind,stoptri,vn,told,snum,snumnew,v0,v1,iter,bnum,dir;
   FLT dpsi,xp[ND],dx,dy,ol;
   
   /* SEARCH FOR TRI ADJACENT TO BOUNDARY NEAR POINT */
   qtree.nearpt(x,y,vnear);
   xp[0] = x; xp[1] = y;
   sind = findbdryside(xp,vnear);
   if (sind < 0) {
      printf("#Warning: brute force boundary locate (%f,%f) %d\n",x,y,typ);
      for(bnum=0;bnum<nsbd;++bnum)
         if (sbdry[bnum].type == typ) break;
      sind = sbdry[bnum].el[sbdry[bnum].num/2];
   }
   else {
      bnum = (-stri[sind][1]>>16) -1;
      if (sbdry[bnum].type != typ) {
         printf("#Warning: brute force boundary locate (%f,%f) sind %d, type: %d %d\n",x,y,sind,sbdry[bnum].type,typ);
         for(bnum=0;bnum<nsbd;++bnum)
            if (sbdry[bnum].type == typ) break;
         sind = sbdry[bnum].el[sbdry[bnum].num/2];
      }
   }
   
#ifdef OLDWAY
   tind = vtri[vnear];
   stoptri = tind;
   dir = 1;
   for(;;) {
      told = tind;
      for(vn=0;vn<3;++vn) 
         if (tvrtx[tind][vn] == vnear) break;
      
      assert(vn != 3);
      
      if (ttri[tind][vn] < 0) {
         bnum = (-ttri[tind][vn]>>16) -1;
         assert(bnum > -1 && bnum < nsbd);
         if (sbdry[bnum].type == typ) {
            sind = tside[tind].side[vn];
            break;
         }
      }
      
      tind = ttri[tind][(vn +dir)%3];
      if (tind < 0) {
         bnum = (-tind>>16) -1;
         assert(bnum > -1 && bnum < nsbd);
         if (sbdry[bnum].type == typ) {
            sind = tside[told].side[(vn+dir)%3];
            break; 
         }
         if (dir > 1) {
            /* DIDN'T FIND SIDE SO DO BRUTE FORCE METHOD */
            printf("#Warning: brute force boundary locate\n");
            for(bnum=0;bnum<nsbd;++bnum)
               if (sbdry[bnum].type == typ) break;
            sind = sbdry[bnum].el[0];
            break;
         }
         /* REVERSE DIRECTION AND GO BACK TO START */
         ++dir;
         tind = vtri[vnear];
         stoptri = -1;
      }
      
      if (tind == stoptri) {
         /* COULDN'T FIND SIDE DO BRUTE FORCE */
         printf("#Warning: brute force bdry locate\n");
         for(bnum=0;bnum<nsbd;++bnum)
            if (sbdry[bnum].type == typ) break;
         sind = sbdry[bnum].el[0];
         break;
      }
   }
#endif
   

   /* SEARCH AROUND THIS SIDE */
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
   };

   /* FIND PSI SUCH THAT TANGENTIAL POSITION ALONG LINEAR SIDE STAYS THE SAME */
   if (typ&CURV_MASK) {
      crdtouht1d(sind);
      dx *= ol;
      dy *= ol;
      iter = 0;
      do {
         b.ptprobe1d(ND,uht,xp,psi);
         
         dpsi = (x -xp[0])*dx +(y -xp[1])*dy;
         psi += dpsi;
         if (iter++ > 100) {
            printf("#Warning: max iterations for curved triangle in bdry_locate tri: %d type: %d loc: %f,%f\n",tind,typ,x,y);
            break;
         }
         /*
         if (psi > 1.0) {
            psi = 1.0;
            break;
         }
         if (psi < -1.0) {
            psi = 1.0;
            break;
         }   */     
      } while (fabs(dpsi) > 100.*EPSILON);
      x = xp[0];
      y = xp[1]; 
   }
   else {
      x = vrtx[v0][0] +dx*(psi +1.)*.5;
      y = vrtx[v0][1] +dy*(psi +1.)*.5;
   }
   
   return(sind);
}