#include"hpbasis.h"
#include<stdio.h>
#include<myblas.h>
#include<utilities.h>
#include <time.h>

#define NOROTATE
#define DEBUG

#define P 4

int main (int argc, const char * argv[]) {
   clock_t cpu_time;
   hpbasis b;
   
   b.initialize(P);
   return 0;

   
#ifdef ROTATE
   FLT *uht[2],**rot;
   FLT **u,**crd[2];
   FLT f;
   int i,j,m,indx,indx1,info;
   char uplo[] = "U";

   vect_alloc(uht[0],b.tm,FLT);
   vect_alloc(uht[1],b.tm,FLT);
   mat_alloc(rot,b.tm,b.tm,FLT);
   mat_alloc(u,b.gpx,b.gpn,FLT);
   mat_alloc(crd[0],b.gpx,b.gpn,FLT);
   mat_alloc(crd[1],b.gpx,b.gpn,FLT);
   
   /* GENERATE ROTATION MATRIX */
   for(i=0;i<b.tm;++i)
      for(j=0;j<b.tm;++j)
         rot[i][j] = 0.0;
         
   /* VERTEX AND SIDE IS EASY */
   for(i=0;i<3;++i)
      rot[i][(i+2)%3] = 1.0;
   
   for(i=0;i<3;++i) {
      indx = 3 +i*b.sm;
      indx1 = 3 +((i+2)%3)*b.sm;
      for(m=0;m<b.sm;++m)
         rot[indx+m][indx1+m] = 1.0;
   }
   
   /* INTERIOR MODES ARE HARD */
   for(indx=0;indx<b.tm;++indx) {
      for(i=0;i<b.tm;++i) {
         uht[0][i] = 0.0;
         uht[1][i] = 0.0;
      }
      uht[0][indx] = 1.0;
      
      /* ROTATE SIDE COEFFICIENT */
      for(i=0;i<b.bm;++i)
         for(j=0;j<b.bm;++j)
            uht[1][i] += rot[i][j]*uht[0][j];
            
      b.proj_bdry(uht[1],u);
      
      /* ROTATE COORDINATES */
      b.proj(1.0,-1.0,-1.0,crd[0]);
      b.proj(-1.0,1.0,-1.0,crd[1]);

      for (i=0;i<b.gpx;++i) {
         for(j=0;j<b.gpn;++j) {
            b.ptprobe(1, uht, &f, crd[0][i][j], crd[1][i][j]);
            u[i][j] -= f;
         }
      }
      
      b.intgrt(u,uht[0]);
      DPBTRS(uplo,b.im,b.ibwth,1,b.idiag[0],b.ibwth+1,&uht[0][b.bm],b.im,info);
      for(i=0;i<b.im;++i)
         rot[i+b.bm][indx] = -uht[0][b.bm+i];
   }
   
   FILE *fp;
   fp = fopen("rotate","w");
   if (fp == NULL)
      printf("trouble\n");
   
   for(i=0;i<b.tm;++i) {
      for(j=0;j<b.tm;++j) {
         fprintf(fp,"%+18.15e ",rot[i][j]);
      }
      fprintf(fp,"\n");
   }
   fclose(fp);
#endif

      

   return 0;
}
