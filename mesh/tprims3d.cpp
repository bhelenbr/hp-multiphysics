#include"mesh.h"
#include<cmath>

template FLT mesh<3>::area(int v0, int v1, int v2) const;

template<int ND> FLT mesh<ND>::area(int v0, int v1, int v2) const {
   int n;
   FLT dx0[ND],dx1[ND];
   FLT at2;
   
   for(n=0;n<ND;++n) {
      dx0[n] =  (vrtx[v1][n]-vrtx[v0][n]);
      dx1[n] =  (vrtx[v2][n]-vrtx[v0][n]);
   }

   at2  = pow(dx0[1]*dx1[2] -dx0[2]*dx1[1],2);
   at2 += pow(dx0[2]*dx1[0] -dx0[0]*dx1[2],2);
   at2 += pow(dx0[0]*dx1[1] -dx0[1]*dx1[0],2);

   return(sqrt(at2));
}

static FLT a[3];

template FLT mesh<3>::intri(int tind, FLT x[3]) const;

template<int ND> FLT mesh<ND>::intri(int tind, FLT x[ND]) const {
   int n,v0,v1,v2;
   FLT dx0[ND],dx1[ND],dx2[ND],x1[2];
   FLT dist1,tdist,dist2;

   v0 = tvrtx[tind][0];
   v1 = tvrtx[tind][1];
   v2 = tvrtx[tind][2];
   
   dist1 = 0.0;
   for(n=0;n<ND;++n) {
      dx0[n] =  (vrtx[v1][n]-vrtx[v0][n]);
      dx1[n] =  (vrtx[v2][n]-vrtx[v0][n]);
      dx2[n] =  (x[n] -vrtx[v0][n]);
      dist1 += pow(dx0[n],2);
   }
   
   /* FIND ORTHOGONAL VECTOR TO DX0 IN PLANE */
   tdist = 0.0;
   for(n=0;n<ND;++n) 
      tdist += dx1[n]*dx0[n];
   tdist /= dist1;
   
   dist2 = 0.0;
   for(n=0;n<ND;++n) {
      dx1[n] -= tdist*dx0[n];
      dist2 += pow(dx1[n],2);
   }
      
   x1[0] = 0.0;
   x1[1] = 0.0;
   for(n=0;n<ND;++n) {
      x1[0] += dx2[n]*dx0[n];
      x1[1] += dx2[n]*dx1[n];
   }
   x1[0] = x1[0]/dist1 -tdist;
   x1[1] /= dist2;
   
   /* VERTICES ARE (0,0) (1,0) (0,1) (NON-ORTHOGONAL VECTORS) */
   a[0] = -x1[0] -x1[1] +1.0;
   a[1] = x1[0];
   a[2] = x1[1];
   
   return(fabs(a[0]) +fabs(a[1]) +fabs(a[2]) - (a[0] +a[1] +a[2]));
}

template void mesh<3>::getwgts(FLT *wt) const;

template<int ND> void mesh<ND>::getwgts(FLT *wt) const {
   int i;
   FLT sum;
   
   sum = a[0] +a[1] +a[2];
   for(i=0;i<3;++i) 
      wt[i] = a[i]/sum;
   
   return;
}
