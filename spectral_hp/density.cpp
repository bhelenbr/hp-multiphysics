/*
 *  density.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"hp_mgrid.h"
#include<math.h>
#include<utilities.h>

FLT rhotemporary;

/* THIS FUNCTION WILL SET THE vlngth VALUES BASED ON THE TRUNCATION ERROR */

void spectral_hp::density1(FLT trncerr, FLT min, FLT max) {
   int i,j,v0,v1,indx,sind,bnum,count;
   FLT sum,u,v,ruv;
   class mesh *tgt;

//   for(i=0;i<nvrtx;++i)
//      vlngth[i] = 1.05*3.1415/20.0*(1.0  - 0.875*exp(-((vrtx[i][0] - center)*(vrtx[i][0] -center) + vrtx[i][1]*vrtx[i][1]) +0.5*0.5));
//
//   for(i=0;i<nvrtx;++i)
//      vlngth[i] = 0.5;

   
   for(i=0;i<nvrtx;++i)
      fltwk[i] = 0.0;

   switch(b.p) {
      case(1):
         for(i=0;i<nside;++i) {
            v0 = svrtx[i][0];
            v1 = svrtx[i][1];
            u = fabs(vug[v0][0] +vug[v1][0]);
            v = fabs(vug[v0][1] +vug[v1][1]);
            ruv = rhotemporary*0.5*(u + v);
            sum = ruv*(fabs(vug[v0][0] -vug[v1][0]) +fabs(vug[v0][1] -vug[v1][1])) +fabs(vug[v0][2] -vug[v1][2]);
            fltwk[v0] += sum;
            fltwk[v1] += sum;
         }
         break;
         
      default:
         indx = 0;
         for(i=0;i<nside;++i) {
            v0 = svrtx[i][0];
            v1 = svrtx[i][1];
            u = fabs(vug[v0][0] +vug[v1][0]);
            v = fabs(vug[v0][1] +vug[v1][1]);
            ruv = rhotemporary*0.5*(u + v);
            sum = ruv*(fabs(sug[indx+b.sm -1][0]) +fabs(sug[indx+b.sm -1][1])) +fabs(sug[indx+b.sm -1][2]);
            fltwk[v0] += sum;
            fltwk[v1] += sum;
            indx += sm0;
         }
         break;
   }
   
   for(i=0;i<nvrtx;++i) {
      fltwk[i] = pow(fltwk[i]/(nnbor[i]*trncerr),1./(b.p+1));
//      fltwk[i] = MAX(0.5,fltwk[i]);
//      fltwk[i] = MIN(2.0,fltwk[i]);
      vlngth[i] /= fltwk[i];
      vlngth[i] = MIN(vlngth[i],max);
      vlngth[i] = MAX(vlngth[i],min);
   }

/*	SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & ALLD_MP) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
/*			SEND VERTEX INFO */
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

void spectral_hp::density2() {
   int i,j,v0,sind,count;

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & ALLD_MP) {
         count = 0;
/*			RECV VERTEX INFO */
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

void spectral_hp::outdensity(char *name) {
   char fnmapp[100];
	FILE *out;
	int i,n,tind;

   strcpy(fnmapp,name);
   strcat(fnmapp,".dat");
   out = fopen(fnmapp,"w");
   if (out == NULL ) {
      printf("couldn't open tecplot output file %s\n",fnmapp);
      exit(1);
   }

   fprintf(out,"ZONE F=FEPOINT, ET=TRIANGLE, N=%d, E=%d\n",nvrtx,ntri);

/*	VERTEX MODES */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n)
         fprintf(out,"%e ",vrtx[i][n]);
      fprintf(out,"%.6e %.6e",vlngth[i],fltwk[i]);					
      fprintf(out,"\n");
   }
   
/*	OUTPUT CONNECTIVY INFO */
   fprintf(out,"\n#CONNECTION DATA#\n");
   
   for(tind=0;tind<ntri;++tind)
      fprintf(out,"%d %d %d\n"
         ,tvrtx[tind][0]+1,tvrtx[tind][1]+1,tvrtx[tind][2]+1);
         
   fclose(out);
   
   return;
}