/*
 *  lengthcd.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on Thu May 29 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#include"hp_mgrid.h"
#include<math.h>
#include<utilities.h>

/* THIS FUNCTION WILL SET THE VLNGTH VALUES BASED ON THE TRUNCATION ERROR */

void hp_mgrid::energy(FLT& esum, FLT& asum) {
   int tind,j,v0;
   FLT q;
   
   for(tind=0;tind<ntri;++tind) {
      q = 0.0;
      for(j=0;j<3;++j) {
         v0 = td(tind).vrtx(j);
         q += pow(ug.v(v0,0),2);
      }
      esum += hp_gbl->rho*q/3.*area(tind);
      asum += area(tind);
   }
   return;
}

void hp_mgrid::length1(FLT norm) {
   int i,j,v0,v1,indx,sind,bnum,count;
   FLT sum,lgtol,lgf,ratio;
   class mesh *tgt;
   
   lgtol = -log(vlngth_tol);
   
   for(i=0;i<nvrtx;++i)
      fscr1(i) = 0.0;

   switch(basis::tri(log2p).p) {
      case(1):
         for(i=0;i<nside;++i) {
            v0 = sd(i).vrtx(0);
            v1 = sd(i).vrtx(1);
            sum = fabs(ug.v(v0,0) -ug.v(v1,0));
            fscr1(v0) += sum;
            fscr1(v1) += sum;
         }
         
         /* BOUNDARY CURVATURE (FOR LINEAR THIS IS DIFFICULT TO TEST)
         for(i=0;i<nsbd;++i) {
            if (!(sbdry[i].type&CURV_MASK)) continue;
<<<<<<< lengthcd.cpp
            for(j=0;j<sbdry(i)->nel;++j) {
               sind = sbdry(i)->el(j);
               v0 = sd(sind).vrtx(0);
               v1 = sd(sind).vrtx(1);
               sum = trncerr*bdrysensitivity*(fabs(vrtx(v0)(0) -vrtx(v1)(0)) +fabs(vrtx(v0)(1) -vrtx(v1)(1)));
               fscr1(v0) += sum;
               fscr1(v1) += sum;
=======
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               v0 = svrtx[sind][0];
               v1 = svrtx[sind][1];
               sum = trncerr*bdrysensitivity*(fabs(vrtx[v0][0] -vrtx[v1][0]) +fabs(vrtx[v0][1] -vrtx[v1][1]));
               fltwk[v0] += sum;
               fltwk[v1] += sum;
>>>>>>> 1.4
            }
         }
         */
                     
         break;
         
      default:
         indx = basis::tri(log2p).sm-1;
         for(i=0;i<nside;++i) {
            v0 = sd(i).vrtx(0);
            v1 = sd(i).vrtx(1);
            sum = fabs(ug.s(indx)(0));
            fscr1(v0) += sum;
            fscr1(v1) += sum;
            indx += sm0;
            
         }

         /* BOUNDARY CURVATURE? */
         for(i=0;i<nsbd;++i) {
            if (!(sbdry[i].type&CURV_MASK)) continue;
<<<<<<< lengthcd.cpp
            indx = basis::tri(log2p).sm-1;
            for(j=0;j<sbdry(i)->nel;++j) {
               sind = sbdry(i)->el(j);
               v0 = sd(sind).vrtx(0);
               v1 = sd(sind).vrtx(1);
               /* THIS LIMITS BOUNDARY CURVATURE TO 1/bdrysensitivity VARIATION */
               sum = pow(bdrysensitivity*0.5*(fabs(binfo[i][indx].curv[0]) +fabs(binfo[i][indx].curv[1]))/distance(v0,v1),(basis::tri(log2p).p+1)/2.0);
               fscr1(v0) += sum*trncerr*nnbor[v0];
               fscr1(v1) += sum*trncerr*nnbor[v1];
=======
            indx = b->sm-1;
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               v0 = svrtx[sind][0];
               v1 = svrtx[sind][1];
               /* THIS LIMITS BOUNDARY CURVATURE TO 1/bdrysensitivity VARIATION */
               sum = pow(bdrysensitivity*0.5*(fabs(binfo[i][indx].curv[0]) +fabs(binfo[i][indx].curv[1]))/distance(v0,v1),(b->p+1)/2.0);
               fltwk[v0] += sum*trncerr*nnbor[v0];
               fltwk[v1] += sum*trncerr*nnbor[v1];
>>>>>>> 1.4
               indx += sm0;
            }
         }
         break;
   }
         
   
   for(i=0;i<nvrtx;++i) {
      fscr1(i) = pow(fscr1(i)/(norm*nnbor[i]*trncerr),1./(basis::tri(log2p).p+1));
      lgf = log(fscr1(i)+EPSILON);
      fscr1(i) = exp(lgtol*lgf/(lgtol +fabs(lgf)));
      vlngth(i) /= fscr1(i);
   }
   
   /* AVOID HIGH ASPECT RATIOS */
   do {
      count = 0;
      for(i=0;i<nside;++i) {
         v0 = sd(i).vrtx(0);
         v1 = sd(i).vrtx(1);
         ratio = vlngth(v1)/vlngth(v0);
         
         if (ratio > 3.0) {
            vlngth(v1) = 2.5*vlngth(v0);
            ++count;
         }
         else if (ratio < 0.333) {
            vlngth(v0) = 2.5*vlngth(v1);
            ++count;
         }
      }
      printf("#aspect ratio fixes %d\n",count);
   } while(count > 0);

   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMY_MASK +IFCE_MASK)) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            tgt->sbuff[bnum][count++] = vlngth(v0);
         }
         v0 = sd(sind).vrtx(1);
         tgt->sbuff[bnum][count++] = vlngth(v0);
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
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            vlngth(v0) = 0.5*(vlngth(v0) +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         vlngth(v0) = 0.5*(vlngth(v0) +sbuff[i][count++]);
      }
   }
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            tgt->sbuff[bnum][count++] = vlngth(v0);
         }
         v0 = sd(sind).vrtx(1);
         tgt->sbuff[bnum][count++] = vlngth(v0);
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
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            vlngth(v0) = 0.5*(vlngth(v0) +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         vlngth(v0) = 0.5*(vlngth(v0) +sbuff[i][count++]);
      }
   }
   
   return;
}

#include<string.h>

void hp_mgrid::outlength(char *name, ftype::name type) {
   char fnmapp[100];
   int j,v0,v1,indx,sind;
   FLT sum;
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
            fprintf(out,"%e\n",vlngth(i));
            
         break;

      case(tecplot):               
         strcat(fnmapp,".dat");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open tecplot output file %s\n",fnmapp);
            exit(1);
         }
         
         for(i=0;i<nvrtx;++i)
            fscr1(i) = 0.0;

         switch(basis::tri(log2p).p) {
            case(1):
               for(i=0;i<nside;++i) {
                  v0 = sd(i).vrtx(0);
                  v1 = sd(i).vrtx(1);
                  sum = fabs(ug.v(v0,0) -ug.v(v1,0));
                  fscr1(v0) += sum;
                  fscr1(v1) += sum;
               }
               
               /* BOUNDARY CURVATURE? */
               for(i=0;i<nsbd;++i) {
                  if (!(sbdry[i].type&CURV_MASK)) continue;
<<<<<<< lengthcd.cpp
                  for(j=0;j<sbdry(i)->nel;++j) {
                     sind = sbdry(i)->el(j);
                     v0 = sd(sind).vrtx(0);
                     v1 = sd(sind).vrtx(1);
                     sum = trncerr*bdrysensitivity*(fabs(vrtx(v0)(0) -vrtx(v1)(0)) +fabs(vrtx(v0)(1) -vrtx(v1)(1)));
                     fscr1(v0) += sum;
                     fscr1(v1) += sum;
=======
                  for(j=0;j<sbdry[i].num;++j) {
                     sind = sbdry[i].el[j];
                     v0 = svrtx[sind][0];
                     v1 = svrtx[sind][1];
                     sum = trncerr*bdrysensitivity*(fabs(vrtx[v0][0] -vrtx[v1][0]) +fabs(vrtx[v0][1] -vrtx[v1][1]));
                     fltwk[v0] += sum;
                     fltwk[v1] += sum;
>>>>>>> 1.4
                  }
               }
               break;
               
            default:
               indx = 0;
               for(i=0;i<nside;++i) {
                  v0 = sd(i).vrtx(0);
                  v1 = sd(i).vrtx(1);
                  sum = fabs(ug.s(indx+basis::tri(log2p).sm -1)(0));
                  fscr1(v0) += sum;
                  fscr1(v1) += sum;
                  indx += sm0;
               }
               
               /* BOUNDARY CURVATURE? */
               for(i=0;i<nsbd;++i) {
                  if (!(sbdry[i].type&CURV_MASK)) continue;
                  indx = 0;
<<<<<<< lengthcd.cpp
                  for(j=0;j<sbdry(i)->nel;++j) {
                     sind = sbdry(i)->el(j);
                     v0 = sd(sind).vrtx(0);
                     v1 = sd(sind).vrtx(1);
                     sum = trncerr*bdrysensitivity*(fabs(binfo[i][indx+basis::tri(log2p).sm-1].curv[0]) +fabs(binfo[i][indx+basis::tri(log2p).sm-1].curv[1]));
                     fscr1(v0) += sum;
                     fscr1(v1) += sum;
=======
                  for(j=0;j<sbdry[i].num;++j) {
                     sind = sbdry[i].el[j];
                     v0 = svrtx[sind][0];
                     v1 = svrtx[sind][1];
                     sum = trncerr*bdrysensitivity*(fabs(binfo[i][indx+b->sm-1].curv[0]) +fabs(binfo[i][indx+b->sm-1].curv[1]));
                     fltwk[v0] += sum;
                     fltwk[v1] += sum;
>>>>>>> 1.4
                     indx += sm0;
                  }
               }
               break;
         }
            
         
         for(i=0;i<nvrtx;++i)
            fscr1(i) = log10(fscr1(i)/nnbor[i]);
      
         fprintf(out,"ZONE F=FEPOINT, ET=TRIANGLE, N=%d, E=%d\n",nvrtx,ntri);
      
         /* VERTEX MODES */
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               fprintf(out,"%e ",vrtx(i)(n));
            fprintf(out,"%.6e %.6e\n",vlngth(i),fscr1(i));               
         }
         
         /* OUTPUT CONNECTIVY INFO */
         fprintf(out,"\n#CONNECTION DATA#\n");
         
         for(tind=0;tind<ntri;++tind)
            fprintf(out,"%d %d %d\n"
               ,td(tind).vrtx(0)+1,td(tind).vrtx(1)+1,td(tind).vrtx(2)+1);

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
      printf("couldn't open vlgth input file %s\n",fnmapp);
      return;
   }
   for(i=0;i<nvrtx;++i)
      fscanf(out,"%lf\n",&vlngth(i));
               
   return;
}

