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
#include<assert.h>
#include<stdlib.h>

void spectral_hp::output(struct vsi g, FLT (*vin)[ND], struct bistruct **bin, char *name, FILETYPE typ = tecplot) {
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
            printf("couldn't open text output file %s\n",fnmapp);
            exit(1);
         }
         
/*			HEADER INFORMATION */
         fprintf(out,"p0 = %d\n",p0);
         fprintf(out,"nvrtx = %d, nside = %d, ntri = %d\n",nvrtx,nside,ntri);
         fprintf(out,"nsbd = %d\n",nsbd);
         for(i=0;i<nsbd;++i)
            fprintf(out," type = %d, number = %d\n",sbdry[i].type,sbdry[i].num);
/*			END HEADER */
         fprintf(out,"END OF HEADER\n");

         for(i=0;i<nvrtx;++i) {
            for(n=0;n<NV;++n)
               fprintf(out,"%15.8e ",g.v[i][n]);
            fprintf(out,"\n");
         }
         
         for(i=0;i<nside*sm0;++i) {
            for(n=0;n<NV;++n)
               fprintf(out,"%15.8e ",g.s[i][n]);
            fprintf(out,"\n");
         }
         
         for(i=0;i<ntri*im0;++i) {
            for(n=0;n<NV;++n)
               fprintf(out,"%15.8e ",g.i[i][n]);
            fprintf(out,"\n");
         }

/*			BOUNDARY INFO */
         for(i=0;i<nsbd;++i) {
            fprintf(out,"Boundary %d, type %d, num %d, frstvrtx %d\n",i,sbdry[i].type,sbdry[i].num,svrtx[sbdry[i].el[0]][0]);
            for(j=0;j<sbdry[i].num*sm0;++j)
               fprintf(out,"%15.8e %15.8e\n",bin[i][j].curv[0],bin[i][j].curv[1]);
         }
         fclose(out);
         break;
      
      case(tecplot):
         strcpy(fnmapp,name);
         strcat(fnmapp,".dat");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open tecplot output file %s\n",fnmapp);
            exit(1);
         }
      
         fprintf(out,"ZONE F=FEPOINT, ET=TRIANGLE, N=%d, E=%d\n",nvrtx+b.sm*nside+b.im*ntri,ntri*(b.sm+1)*(b.sm+1));

/*			VERTEX MODES */
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               fprintf(out,"%e ",vin[i][n]);
            for(n=0;n<NV;++n)
               fprintf(out,"%.6e ",g.v[i][n]);					
            fprintf(out,"\n");
         }
         
         if (b.p > 1) {
/*				SIDE MODES */
            for(sind=0;sind<nside;++sind) {
               if (sinfo[sind] < 0) {
                  v0 = svrtx[sind][0];
                  v1 = svrtx[sind][1];
                  for(n=0;n<ND;++n)
                     b.proj1d_leg(vin[v0][n],vin[v1][n],crd[n][0]);
               }
               else {
                  crdtouht1d(sind,vin,bin);
                  for(n=0;n<ND;++n)
                     b.proj1d_leg(uht[n],crd[n][0]);
               }
               
               ugtouht1d(sind,g);
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
                  ugtouht(tind,g);
                  for(n=0;n<NV;++n)
                     b.proj_leg(uht[n],u[n]);
                     
                  if (tinfo[tind] < 0) {
                     for(n=0;n<ND;++n)
                        b.proj_leg(vin[tvrtx[tind][0]][n],vin[tvrtx[tind][1]][n],vin[tvrtx[tind][2]][n],crd[n]);
                  }
                  else {
                     crdtouht(tind,vin,bin);
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
         fclose(out);
         break;
         
      default:
         printf("That filetype is not supported for spectral_hp output\n");
         break;
  	}
   
 	return;
}

#define NEW

void spectral_hp::input(struct vsi g, FLT (*vin)[ND], struct bistruct **bin, char *name, FILETYPE typ = text) {
   int i,j,k,m,n,pin,indx,bind;
   int bnum,innum,intyp,ierr,v0;
   FILE *in;
   char buffer[200];
   
   switch(typ) {
   
#ifdef NEW
      case (text):
         strcpy(buffer,name);
         strcat(buffer,".txt");
         in = fopen(buffer,"r");
         if (in == NULL ) {
            printf("couldn't open text input file %s\n",name);
            exit(1);
         }
         
/*			HEADER INFORMATION */
/*			INPUT # OF SIDE MODES (ONLY THING THAT CAN BE DIFFERENT) THEN SKIP THE REST */
         fscanf(in,"p0 = %d\n",&pin);
         
         do {
            fscanf(in,"%[^\n]",buffer);
            fscanf(in,"\n");
         } while (buffer[0] != 'E');  // END OF HEADER

         for(i=0;i<nvrtx;++i) {
            for(n=0;n<NV;++n)
               fscanf(in,"%le ",&g.v[i][n]);
            fscanf(in,"\n");
         }
         
         indx = 0;
         for(i=0;i<nside;++i) {
            for(m=0;m<(pin-1);++m) {
               for(n=0;n<NV;++n)
                  fscanf(in,"%le ",&g.s[indx][n]);
               fscanf(in,"\n");
               ++indx;
            }
            indx += p0 -pin;
         }
         
         indx = 0;
         for(i=0;i<ntri;++i) {
            for(m=1;m<pin-1;++m) {
               for(k=0;k<pin-1-m;++k) {
                  for(n=0;n<NV;++n) 
                     fscanf(in,"%le ",&g.i[indx][n]);
                  fscanf(in,"\n");
                  ++indx;
               }
               indx += p0 -pin;
            }
         }
         

/*			BOUNDARY INFO */
         for(i=0;i<nsbd;++i) {
            fscanf(in,"Boundary %*d, type %d, num %d, frstvrtx %d\n",&intyp,&innum,&v0);
/*				FIND MATCHING BOUNDARY (CAN CHANGE NUMBER/ORDER) */
            for(bnum=0;bnum<nsbd;++bnum)
               if (svrtx[sbdry[bnum].el[0]][0] == v0) 
                  break;
            
            if (bnum == nsbd || innum != sbdry[bnum].num) {
               printf("Trouble reading boundary info %d %d %d %d\n",innum,sbdry[bnum].num,sbdry[bnum].type,intyp);
               exit(1);
            }

            indx = 0;
            for(j=0;j<sbdry[bnum].num;++j) {
               for(m=0;m<pin -1;++m) {
                  fscanf(in,"%le %le\n",&bin[bnum][indx].curv[0],&bin[bnum][indx].curv[1]);
                  ++indx;
               }
               indx += p0 -pin;
            }
         }
         fclose(in);
         break;
#else
      case (text):
         in = fopen(name,"r");
         if (in == NULL ) {
            printf("couldn't open text input file %s\n",name);
            exit(1);
         }
         
/*			SKIP FIRST LINE */
         fscanf(in,"%*[^\n]");
         
         for(i=0;i<nvrtx;++i) {
            ierr = fscanf(in,"%le %le %le %le %le %*e %*e\n",
            &vin[i][0],&vin[i][1],&g.v[i][0],&g.v[i][1],&g.v[i][2]);
            if(ierr != 5) {
               printf("error in read file %d\n",i);
               exit(1);
            }
         }
      
         if (b.sm > 0) {
            for(i=0;i<nside;++i) {
               indx = b.sm*i;
               if (sinfo[i] > -1) {
                  bnum = (-stri[i][1]>>16) -1;
                  bind = (-stri[i][1]&0xFFFF)*sm0;
                  assert(bnum > -1 && bnum < nsbd);
                  assert(indx > -1 && indx <= sbdry[bnum].num*sm0);
                  
                  for(m=0;m<b.sm;++m) {
                     ierr = fscanf(in,"%le %le %le %le %le %*e %*e\n",
                     &bin[bnum][bind+m].curv[0],&bin[bnum][bind+m].curv[1],
                     &g.s[indx+m][0],&g.s[indx+m][1],&g.s[indx+m][2]);
                     if(ierr != 5) {
                        printf("error in reading curved side %d\n",i);
                        exit(1);
                     }
                  }
               }
               else {
                  for(m=0;m<b.sm;++m) {
                     ierr = fscanf(in,"%le %le %le\n",
                     &g.s[indx+m][0],&g.s[indx+m][1],&g.s[indx+m][2]);
                     if(ierr != 3) {
                        printf("error in reading straight side %d\n",i);
                        exit(1);
                     }
                  }
               }
            }
         
            if (b.im > 0) {
               for(i=0;i<b.im*ntri;++i) {
                  ierr = fscanf(in,"%le %le %le\n",
                  &g.i[i][0],&g.i[i][1],&g.i[i][2]);
                  if(ierr != 3) {
                     printf("error in read file6\n");
                     exit(1);
                  }
               }
            }
         }
         
         break;
#endif
      default:
         printf("Spectral_hp input of that filetype is not supported\n");
         exit(1);
         break;
   }
   
   return;
}
