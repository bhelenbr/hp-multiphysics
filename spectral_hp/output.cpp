/*
 *  output.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 15 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"
#include<cstring>

void spectral_hp::output(char *name, FILETYPE typ = tecplot) {
   char fnmapp[100];
	FILE *out;
	int i,j,k,n,v0,v1,sind,tind,indx,sgn;
   static int ijind[MXTM][MXTM];

   switch (typ) {
      case (text):
         strcpy(fnmapp,name);
         strcat(fnmapp,".txt");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file\n");
            exit(1);
         }
         
/*			HEADER INFORMATION */
         fprintf(out,"sm = %d\n",sm0);
         fprintf(out,"nvrtx = %d, nside = %d, ntri = %d\n",nvrtx,nside,ntri);
         fprintf(out,"nsbd = %d\n",nsbd);
         for(i=0;i<nsbd;++i)
            fprintf(out," type = %d, number = %d\n",sbdry[i].type,sbdry[i].num);
/*			END HEADER */
         fprintf(out,"END OF HEADER\n");

         for(i=0;i<nvrtx;++i) {
            for(n=0;n<NV;++n)
               fprintf(out,"%15.8e ",vug[i][n]);
            fprintf(out,"\n");
         }
         
         for(i=0;i<nside*sm0;++i) {
            for(n=0;n<NV;++n)
               fprintf(out,"%15.8e ",sug[i][n]);
            fprintf(out,"\n");
         }
         
         for(i=0;i<ntri*im0;++i) {
            for(n=0;n<NV;++n)
               fprintf(out,"%15.8e ",iug[i][n]);
            fprintf(out,"\n");
         }

/*			BOUNDARY INFO */
         for(i=0;i<nsbd;++i) {
            fprintf(out,"Boundary %d, type %d, num %d\n",i,sbdry[i].type,sbdry[i].num);
            for(j=0;j<sbdry[i].num*sm0;++j)
               fprintf(out,"%15.8e %15.8e %15.8e\n",binfo[i][j].curv[0],binfo[i][j].curv[1],binfo[i][j].sfct);
         }
         fclose(out);
         break;
      
      case(tecplot):
         strcpy(fnmapp,name);
         strcat(fnmapp,".dat");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file\n");
            exit(1);
         }
      
         fprintf(out,"ZONE F=FEPOINT, ET=TRIANGLE, N=%d, E=%d\n",nvrtx+b.sm*nside+b.im*ntri,ntri*(b.sm+1)*(b.sm+1));

/*			VERTEX MODES */
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               fprintf(out,"%e ",vrtx[i][n]);
            for(n=0;n<NV;++n)
               fprintf(out,"%.6e ",vug[i][n]);					
            fprintf(out,"\n");
         }
         
         if (b.p > 1) {
/*				SIDE MODES */
            for(sind=0;sind<nside;++sind) {
               if (sinfo[sind] < 0) {
                  v0 = svrtx[sind][0];
                  v1 = svrtx[sind][1];
                  for(n=0;n<ND;++n)
                     b.proj1d_leg(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
               }
               else {
                  crdtouht1d(sind);
                  for(n=0;n<ND;++n)
                     b.proj1d_leg(uht[n],crd[n][0]);
               }
               
               ugtouht1d(sind);
               for(n=0;n<NV;++n)
                  b.proj1d_leg(uht[n],u[n][0]);

               for(i=1;i<b.sm+1;++i) {
                  for(n=0;n<ND;++n)
                     fprintf(out,"%e ",crd[n][0][i]);
                  for(n=0;n<NV;++n)
                     fprintf(out,"%.6e ",u[n][0][i]);					
                  fprintf(out,"\n");
               }
            }
   
/*				INTERIOR MODES */
            if (b.p > 2) {
               for(tind = 0; tind < ntri; ++tind) {
                  ugtouht(tind);
                  for(n=0;n<NV;++n)
                     b.proj_leg(uht[n],u[n]);
                     
                  if (tinfo[tind] < 0) {
                     for(n=0;n<ND;++n)
                        b.proj_leg(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);
                  }
                  else {
                     crdtouht(tind);
                     for(n=0;n<ND;++n)
                        b.proj_bdry_leg(uht[n],crd[n]);
                  }
                  
                  for(i=1;i<b.sm;++i) {
                     for(j=1;j<b.sm-(i-1);++j) {
                        for(n=0;n<ND;++n)
                           fprintf(out,"%e ",crd[n][i][j]);
                        for(n=0;n<NV;++n)
                           fprintf(out,"%.6e ",u[n][i][j]);					
                        fprintf(out,"\n");
                     }
                  }
               }
            }
         }
      
/*			OUTPUT CONNECTIVY INFO */
         fprintf(out,"\n#CONNECTION DATA#\n");
         
         for(tind=0;tind<ntri;++tind) {

/*				VERTICES */
            ijind[0][0] = tvrtx[tind][0];
            ijind[b.sm+1][0] = tvrtx[tind][1];
            ijind[0][b.sm+1] = tvrtx[tind][2];
   
/*				SIDES */		   
            indx = tside[tind].side[0];
            sgn = tside[tind].sign[0];
            if (sgn > 0) {
               for(i=0;i<b.sm;++i)
                  ijind[b.sm-i][i+1] = nvrtx +indx*b.sm +i;
            }
            else {
               for(i=0;i<b.sm;++i)
                  ijind[b.sm-i][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
            }

            indx = tside[tind].side[1];
            sgn = tside[tind].sign[1];
            if (sgn > 0) {
               for(i=0;i<b.sm;++i)
                  ijind[0][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
            }
            else {
               for(i=0;i<b.sm;++i)
                  ijind[0][i+1] = nvrtx +indx*b.sm +i;
            }
            
            indx = tside[tind].side[2];
            sgn = tside[tind].sign[2];
            if (sgn < 0) {
               for(i=0;i<b.sm;++i)
                  ijind[i+1][0] = nvrtx +(indx+1)*b.sm -(i+1);
            }
            else {
               for(i=0;i<b.sm;++i)
                  ijind[i+1][0] = nvrtx +indx*b.sm +i;
            }
   
/*				INTERIOR VERTICES */
            k = 0;
            for(i=1;i<b.sm;++i) {
               for(j=1;j<b.sm-(i-1);++j) {
                  ijind[i][j] = nvrtx +nside*b.sm +tind*b.im +k;
                  ++k;
               }
            }
   
/*				OUTPUT CONNECTION LIST */		
            for(i=0;i<b.sm+1;++i) {
               for(j=0;j<b.sm-i;++j) {
                  fprintf(out,"%d %d %d\n"
                  ,ijind[i][j]+1,ijind[i+1][j]+1,ijind[i][j+1]+1);
                  fprintf(out,"%d %d %d\n"
                  ,ijind[i+1][j]+1,ijind[i+1][j+1]+1,ijind[i][j+1]+1);
               }
               fprintf(out,"%d %d %d\n"
               ,ijind[i][b.sm-i]+1,ijind[i+1][b.sm-i]+1,ijind[i][b.sm+1-i]+1);
            }
         }
         break;
         
      case(easymesh):
         out_mesh(name,easymesh);
         break;
      
      case(gambit):
         out_mesh(name,gambit);
         break;
  	}
   
 	return;
}


