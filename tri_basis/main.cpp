#include"hpbasis.h"
#include<stdio.h>
#include<myblas.h>
#include<utilities.h>

/* LAGRANGIAN1D / MASS */
#define MASS

#define P 12
int main (int argc, const char * argv[]) {
   int i,j,info;
   hpbasis b;
      
   b.initialize(P);      
     
#ifdef LAGRANGIAN1D
   double l1d[P+1][P+1], m[P+1][P+1], k[P+1][P+1],work[P+1][P+1];
   double c[P+1],d[P+1],e[P+1],f[P+1];
   int ipiv[2*(P+1)];
   char trans[] = "T";
   
   for(i=0;i<P+1;++i) {
      for(j=0;j<P+1;++j) {
         l1d[i][j] = 0.0;
         m[i][j] = 0.0;
         k[i][j] = 0.0;
      }
   }
   
   /* VERTEX COMPONENTS */
   l1d[0][0] = 1;   
   l1d[1][P] = 1;
   
   /* SIDE COMPONENTS */
   for(i=0;i<P+1;++i) {
      for(j=0;j<P+1;++j)
         c[j] = 0.0;
      c[i] = 1.0;
      b.proj1d_leg(c,d);
      for(j=1;j<P;++j)
         l1d[i][j] = d[j];
   }

/*   
   for(i=0;i<P+1;++i) {
      for(j=0;j<P+1;++j)
         printf("%3.2e ",l1d[j][i]);
      printf("\n");
   }
*/
   
   for(i=0;i<P+1;++i) {
      for(j=0;j<P+1;++j) {
         c[j] = 0.0;
         f[j] = 0.0;
      }
         
      c[i] = 1.0;
      
      b.proj1d(c,d,e);
      b.intgrt1d(d,c);
      b.intgrtx1d(e,f);
      for(j=0;j<P+1;++j) {
         m[i][j] = 0.5*c[j];
         k[i][j] = 2.0*f[j];
      }
   }

   
   printf("\nMASS\n");
   for(i=0;i<P+1;++i) {
      for(j=0;j<P+1;++j)
         printf("%+3.2f ",m[j][i]);
      printf("\n");
   }
   return 0;
   
/*
   
   printf("\nSTIFFNESS\n");
   for(i=0;i<P+1;++i) {
      for(j=0;j<P+1;++j)
         printf("%3.2e ",k[j][i]);
      printf("\n");
   }
*/
   
   
   /* FORM INVERSE OF L1D */
   GETRF(P+1,P+1,&l1d[0][0],P+1,ipiv,info);
   if (info != 0) {
      printf("DGETRF FAILED\n");
      exit(1);
   }
   GETRS(trans,P+1,P+1,&l1d[0][0],P+1,ipiv,&m[0][0],P+1,info);
   
   for(i=0;i<P+1;++i)
      for(j=0;j<P+1;++j)
         work[i][j] = m[j][i];
         
   GETRS(trans,P+1,P+1,&l1d[0][0],P+1,ipiv,&work[0][0],P+1,info);
   
   for(i=0;i<P+1;++i)
      for(j=0;j<P+1;++j)
         m[i][j] = work[j][i];
         
   printf("\nMASS\n");
   for(i=0;i<P+1;++i) {
      for(j=0;j<P+1;++j)
         printf("%7.6e ",m[j][i]);
      printf("\n");
   }   
   
   GETRS(trans,P+1,P+1,&l1d[0][0],P+1,ipiv,&k[0][0],P+1,info);
   
   for(i=0;i<P+1;++i)
      for(j=0;j<P+1;++j)
         work[i][j] = k[j][i];
         
   GETRS(trans,P+1,P+1,&l1d[0][0],P+1,ipiv,&work[0][0],P+1,info);
   
   for(i=0;i<P+1;++i)
      for(j=0;j<P+1;++j)
         k[i][j] = work[j][i];
         
   printf("\nSTIFFNESS\n");
   for(i=0;i<P+1;++i) {
      for(j=0;j<P+1;++j)
         printf("%7.6e ",k[j][i]);
      printf("\n");
   }   
#endif

#ifdef MASS
   int one = 1;
   int k,cnt,ipiv[2*P*P];
   double m[P*P][P*P],mlmp[P*P][P*P],c[P*P],d[P*P];
   double **gp;
   char trans[] = "T";

   mat_alloc(gp,P*P,P*P,double);
   
   for(i=0;i<b.tm;++i) {
      for(j=0;j<b.tm;++j)
         c[j] = 0.0;
         
      c[i] = 1.0;
      
      b.proj(c,gp);
      b.intgrt(gp,c);
      for(j=0;j<b.tm;++j)
         m[i][j] = c[j];
   }
   

   

   for(i=0;i<b.tm;++i) {
      for(j=0;j<b.tm;++j)
         c[j] = 0.0;
         
      c[i] = 1.0;
      
      if (i < b.bm) {
         for(j=b.bm;j<b.tm;++j)
            c[j] = -1.0/m[j][j]*m[i][j];
            
         cnt = b.bm;
         for(j=2;j<P;++j) {
            for(k=1;k<P-j;++k)
               ++cnt;
            c[cnt++] = 0.0;
         }
      }
      
      b.proj(c,gp);
      b.intgrt(gp,c);
      for(j=0;j<b.tm;++j)
         mlmp[i][j] = c[j];
   }
   
   printf("\nMLMP\n");
   for(i=0;i<b.tm;++i) {
      for(j=0;j<b.tm;++j) {
         if (fabs(mlmp[i][j]) > 1.0e-10) 
            printf("* ");
         else 
            printf(". "); 
      }
      printf("\n");
   }
   
   FILE *fp;
   fp = fopen("mlmp","w");
   if (fp == NULL)
      printf("trouble\n");
   
   for(i=0;i<b.bm;++i) {
      for(j=0;j<b.bm;++j) {
        // if (fabs(m[i][j]) > 1.0e-10) 
        //    printf("* ");
        // else 
        //    printf(". ");
         fprintf(fp,"%+14.10e ",mlmp[i][j]);
      }
      fprintf(fp,"\n");
   }
   
   int ind = b.bm;
   for(i=2;i<P;++i) {
      for(int n=1;n<P-i;++n)
         ++ind;
      for(j=0;j<b.bm;++j) {
        // if (fabs(m[i][j]) > 1.0e-10) 
        //    printf("* ");
        // else 
        //    printf(". ");
         fprintf(fp,"%+14.10e ",mlmp[ind][j]);
      }
      fprintf(fp,"\n");
      ++ind;
   }
   
   
   return(0);
   
   printf("\nSIDE 0 LUMPED\n");
   
   /* SIDE 0 LUMPED INVERSION */
   for(i=0;i<b.sm-1;++i) {
      cnt = b.bm;
      for(j=2;j<i+2;++j) 
         for(k=1;k<=P-j;++k)
            ++cnt;
            
      for(j=i;j<b.sm-1;++j) {
         for(k=i;k<b.sm-1;++k)
            mlmp[j-i][k-i] = m[cnt+k-i][3+b.sm+j];
         c[j-i] = m[3+i][3+b.sm+j];
      }
      /* FORM INVERSE */
      GETRF(b.sm-1-i,b.sm-1-i,&mlmp[0][0],P*P,ipiv,info);
      if (info != 0) {
         printf("DGETRF FAILED\n");
         exit(1);
      }
      GETRS(trans,b.sm-1-i,1,&mlmp[0][0],P*P,ipiv,&c[0],P*P,info);
      
      /* CALCULATE NEW LUMPED MATRIX */
      for(j=0;j<b.tm;++j)
         d[j] = 0.0;
         
      d[3+i] = 1.0;
      for(k=i;k<b.sm-1;++k)
         d[cnt+k-i] = -c[k-i];
      
      b.proj(d,gp);
      b.intgrt(gp,c);
      for(j=0;j<b.tm;++j) {
         if (fabs(c[j]) > 1.0e-10) 
            printf("* ");
         else 
            printf(". ");
      }
      printf("\n");
   }
   
#endif
   
   return 0;
}
