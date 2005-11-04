/*
 *  length.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#include"hp_mgrid.h"
#include<math.h>
#include<utilities.h>

/* THIS FUNCTION WILL SET THE VLNGTH VALUES BASED ON THE TRUNCATION ERROR */
#ifdef DROP
extern FLT dydt;
#endif

void hp_mgrid::energy(FLT& esum, FLT& asum) {
   int tind,j,v0;
   FLT q,p,duv,um,vm,u,v;
   
   for(tind=0;tind<ntri;++tind) {
      q = 0.0;
      p = 0.0;
      duv = 0.0;
      um = ug.v[tvrtx[tind][2]][0];
      vm = ug.v[tvrtx[tind][2]][1];
#ifdef DROP
      vm -= dydt;
#endif
      for(j=0;j<3;++j) {
         v0 = tvrtx[tind][j];
         u = ug.v[v0][0];
         v = ug.v[v0][1];
#ifdef DROP
         v -= dydt;
#endif
         q += pow(u,2) +pow(v,2);
         p += fabs(ug.v[v0][2]);
         duv += fabs(u-um)+fabs(v-vm);
         um = u;
         vm = v;
      }
      esum += 1./3.*( (0.5*gbl->rho*q +p)*area(tind) +duv*gbl->mu*sqrt(area(tind)) );
      asum += area(tind);
   }
   return;
}

void hp_mgrid::length1(FLT norm) {
   int i,j,k,v0,v1,v2,indx,sind,tind,bnum,count;
   FLT sum,u,v,ruv,lgtol,lgf,ratio;
   class mesh *tgt;
   FLT dx0[2],dx1[2],dx2[2],dedpsi[2],ep[2];
   FLT length0,length1,length2,lengthept;
   FLT ang1,curved1,ang2,curved2;

   lgtol = -log(vlngth_tol);
   
   for(i=0;i<nvrtx;++i)
      fltwk[i] = 0.0;

   switch(b->p) {
      case(1):
         for(i=0;i<nside;++i) {
            v0 = svrtx[i][0];
            v1 = svrtx[i][1];
            u = fabs(ug.v[v0][0] +ug.v[v1][0]);
            v = fabs(ug.v[v0][1] +ug.v[v1][1]);
#ifdef DROP
            v -= dydt;
#endif
            ruv = gbl->rho*0.5*(u + v) +gbl->mu/distance(v0,v1);
            sum = distance2(v0,v1)*(ruv*(fabs(ug.v[v0][0] -ug.v[v1][0]) +fabs(ug.v[v0][1] -ug.v[v1][1])) +fabs(ug.v[v0][2] -ug.v[v1][2]));
            fltwk[v0] += sum;
            fltwk[v1] += sum;
         }
         break;
         
      default:
         indx = b->sm-1;
         for(i=0;i<nside;++i) {
            v0 = svrtx[i][0];
            v1 = svrtx[i][1];
            u = fabs(ug.v[v0][0] +ug.v[v1][0]);
            v = fabs(ug.v[v0][1] +ug.v[v1][1]);
#ifdef DROP
            v -= dydt;
#endif
            ruv = gbl->rho*0.5*(u + v) +gbl->mu/distance(v0,v1);
            sum = distance2(v0,v1)*(ruv*(fabs(ug.s[indx][0]) +fabs(ug.s[indx][1])) +fabs(ug.s[indx][2]));
            fltwk[v0] += sum;
            fltwk[v1] += sum;
            indx += sm0;
            
         }

        /* BOUNDARY CURVATURE */
         for(i=0;i<nsbd;++i) {
            if (!(sbdry[i].type&CURV_MASK)) continue;
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               v1 = svrtx[sind][0];
               v2 = svrtx[sind][1];
               
               crdtocht1d(sind);
                                    
               /* FIND ANGLE BETWEEN LINEAR SIDES */
               tind = stri[sind][0];
               for(k=0;k<3;++k)
                  if (tside[tind].side[k] == sind) break;
               
               v0 = tvrtx[tind][k];
               
               dx0[0] = vrtx[v2][0]-vrtx[v1][0];
               dx0[1] = vrtx[v2][1]-vrtx[v1][1];
               length0 = dx0[0]*dx0[0] +dx0[1]*dx0[1];
               
               dx1[0] = vrtx[v0][0]-vrtx[v2][0];
               dx1[1] = vrtx[v0][1]-vrtx[v2][1];
               length1 = dx1[0]*dx1[0] +dx1[1]*dx1[1];
               
               dx2[0] = vrtx[v1][0]-vrtx[v0][0];
               dx2[1] = vrtx[v1][1]-vrtx[v0][1];
               length2 = dx2[0]*dx2[0] +dx2[1]*dx2[1];
               
               b->ptprobe1d(2,ep,dedpsi,-1.0,cht[0],MXTM);
               lengthept = dedpsi[0]*dedpsi[0] +dedpsi[1]*dedpsi[1];
               
               ang1 = acos(-(dx0[0]*dx2[0] +dx0[1]*dx2[1])/sqrt(length0*length2));
               curved1 = acos((dx0[0]*dedpsi[0] +dx0[1]*dedpsi[1])/sqrt(length0*lengthept));
               
               b->ptprobe1d(2,ep,dedpsi,1.0,cht[0],MXTM);
               lengthept = dedpsi[0]*dedpsi[0] +dedpsi[1]*dedpsi[1];
               
               ang2 = acos(-(dx0[0]*dx1[0] +dx0[1]*dx1[1])/sqrt(length0*length1));
               curved2 = acos((dx0[0]*dedpsi[0] +dx0[1]*dedpsi[1])/sqrt(length0*lengthept));                     

               sum = bdrysensitivity*(curved1/ang1 +curved2/ang2);
               fltwk[v0] += sum*trncerr*norm*nnbor[v0];
               fltwk[v1] += sum*trncerr*norm*nnbor[v1];
            }
         }
         break;
   }

   for(i=0;i<nvrtx;++i) {
      fltwk[i] = pow(fltwk[i]/(norm*nnbor[i]*trncerr),1./(b->p+1+ND));
#ifndef DROP
      lgf = log(fltwk[i]);
      fltwk[i] = exp(lgtol*lgf/(lgtol +fabs(lgf)));
#endif
      vlngth[i] /= fltwk[i];
#ifdef ONELAYER
      vlngth[i] = MIN(vlngth[i],0.5);
#endif
#ifdef THREELAYER
#define TRES 0.025/THREELAYER
      if (vrtx[i][1] > 0.525) {
         vlngth[i] = MIN(vlngth[i],TRES +(vrtx[i][1]-0.525)*(9*TRES)/0.475);
      }
      else if (vrtx[i][1] < 0.475) {
         vlngth[i] = MIN(vlngth[i],TRES +(0.475 -vrtx[i][1])*(9*TRES)/0.475);
      }
      else {
         vlngth[i] = MIN(vlngth[i],TRES);
      }
#endif
#ifdef TWOLAYER
      vlngth[i] = MIN(vlngth[i],0.3333); 
#endif
   }
   
   /* AVOID HIGH ASPECT RATIOS */
   int nsweep = 0;
   do {
      count = 0;
      for(i=0;i<nside;++i) {
         v0 = svrtx[i][0];
         v1 = svrtx[i][1];
         ratio = vlngth[v1]/vlngth[v0];
         
         if (ratio > 3.0) {
            vlngth[v1] = 2.5*vlngth[v0];
            ++count;
         }
         else if (ratio < 0.333) {
            vlngth[v0] = 2.5*vlngth[v1];
            ++count;
         }
      }
      ++nsweep;
      printf("#aspect ratio fixes %d: %d\n",nsweep,count);
   } while(count > 0 && nsweep < 5);

   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMY_MASK +IFCE_MASK)) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = vlngth[v0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = vlngth[v0];
      }
   }

   return;
   
}

void hp_mgrid::length_mp() {
   int i,j,v0,sind,bnum,count;
   class mesh *tgt;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMY_MASK +IFCE_MASK)) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            vlngth[v0] = 0.5*(vlngth[v0] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         vlngth[v0] = 0.5*(vlngth[v0] +sbuff[i][count++]);
      }
   }
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = vlngth[v0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = vlngth[v0];
      }
   }
   return;
}

void hp_mgrid::length2() {
   int i,j,v0,sind,count;

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            vlngth[v0] = 0.5*(vlngth[v0] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         vlngth[v0] = 0.5*(vlngth[v0] +sbuff[i][count++]);
      }
   }
   
   return;
}

#include<string.h>

void hp_mgrid::outlength(char *name, FILETYPE type) {
   char fnmapp[100];
   int v0,v1,indx;
   FLT sum,u,v,ruv;
   FILE *out;
   int i,n,tind;

   strcpy(fnmapp,name);

   switch(type) {
      case(text):
         strcat(fnmapp,".txt");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open tecplot output file %s\n",fnmapp);
            exit(1);
         }
         for(i=0;i<nvrtx;++i)
            fprintf(out,"%e\n",vlngth[i]);
            
         break;

      case(tecplot):               
         strcat(fnmapp,".lgth.dat");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open tecplot output file %s\n",fnmapp);
            exit(1);
         }
         
         for(i=0;i<nvrtx;++i)
            fltwk[i] = 0.0;

         switch(b->p) {
            case(1):
               for(i=0;i<nside;++i) {
                  v0 = svrtx[i][0];
                  v1 = svrtx[i][1];
                  u = fabs(ug.v[v0][0] +ug.v[v1][0]);
                  v = fabs(ug.v[v0][1] +ug.v[v1][1]);
                  ruv = 0.5*gbl->rho*(u + v);
                  sum = distance2(v0,v1)*(ruv*(fabs(ug.v[v0][0] -ug.v[v1][0]) +fabs(ug.v[v0][1] -ug.v[v1][1])) +fabs(ug.v[v0][2] -ug.v[v1][2]));
                  fltwk[v0] += sum;
                  fltwk[v1] += sum;
               }
               break;
               
            default:
               indx = 0;
               for(i=0;i<nside;++i) {
                  v0 = svrtx[i][0];
                  v1 = svrtx[i][1];
                  u = fabs(ug.v[v0][0] +ug.v[v1][0]);
                  v = fabs(ug.v[v0][1] +ug.v[v1][1]);
                  ruv = 0.5*gbl->rho*(u + v);
                  sum = distance2(v0,v1)*(ruv*(fabs(ug.s[indx+b->sm -1][0]) +fabs(ug.s[indx+b->sm -1][1])) +fabs(ug.s[indx+b->sm -1][2]));
                  fltwk[v0] += sum;
                  fltwk[v1] += sum;
                  indx += sm0;
               }
               
               break;
         }
            
         
         for(i=0;i<nvrtx;++i)
            fltwk[i] = log10(fltwk[i]/(nnbor[i]*trncerr));
      
         fprintf(out,"ZONE F=FEPOINT, ET=TRIANGLE, N=%d, E=%d\n",nvrtx,ntri);
      
         /* VERTEX MODES */
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               fprintf(out,"%e ",vrtx[i][n]);
            fprintf(out,"%.6e %.6e\n",vlngth[i],fltwk[i]);               
         }
         
         /* OUTPUT CONNECTIVY INFO */
         fprintf(out,"\n#CONNECTION DATA#\n");
         
         for(tind=0;tind<ntri;++tind)
            fprintf(out,"%d %d %d\n"
               ,tvrtx[tind][0]+1,tvrtx[tind][1]+1,tvrtx[tind][2]+1);

         break;
         
      default:
         printf("Output of length function is not supported in that filetype\n");
         break;
   }
         
   fclose(out);
   
   return;
}

void hp_mgrid::inlength(char *name) {
   char fnmapp[100];
   FILE *out;
   int i;

   strcpy(fnmapp,name);
   strcat(fnmapp,".txt");
   out = fopen(fnmapp,"r");
   if (out == NULL ) {
      printf("#couldn't open vlngth input file %s. Reinitializing vlngth\n",fnmapp);
      initvlngth();
      return;
   }
   for(i=0;i<nvrtx;++i) 
      fscanf(out,"%le\n",&vlngth[i]);
               
   return;
}

