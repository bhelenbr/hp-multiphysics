#include"mesh.h"
#include<cmath>

template FLT mesh<2>::incircle(int tind, FLT *a) const;
template FLT mesh<3>::incircle(int tind, FLT *a) const;

template<int ND> FLT mesh<ND>::incircle(int tind, FLT *a) const {
   int i;
   FLT pt[3][2];
   FLT l2[3];
   FLT determ;
   
   for(i=0;i<3;++i) {
      pt[i][0] = vrtx[tvrtx[tind][i]][0]-a[0];
      pt[i][1] = vrtx[tvrtx[tind][i]][1]-a[1];
   }

   for(i=0;i<3;++i)
      l2[i] = pt[i][0]*pt[i][0] +pt[i][1]*pt[i][1];

   determ = 0.0;
   determ +=  pt[0][0]*(pt[1][1]*l2[2] -pt[2][1]*l2[1]);
   determ += -pt[0][1]*(pt[1][0]*l2[2] -pt[2][0]*l2[1]);
   determ +=  l2[0]*(pt[1][0]*pt[2][1] -pt[1][1]*pt[2][0]);
   
   return(determ);
}

template FLT mesh<2>::insidecircle(int sind, FLT *a) const;
template FLT mesh<3>::insidecircle(int sind, FLT *a) const;

template<int ND> FLT mesh<ND>::insidecircle(int sind, FLT *a) const {
   int v0,v1;
   FLT ctr[2],dist2;
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   ctr[0] = 0.5*(vrtx[v0][0] +vrtx[v1][0]);
   ctr[1] = 0.5*(vrtx[v0][1] +vrtx[v1][1]);
   dist2 = (a[0]-ctr[0])*(a[0]-ctr[0]) +(a[1]-ctr[1])*(a[1]-ctr[1]);
   return(0.25*distance2(v0,v1) -dist2);
}


template FLT mesh<2>::area(int v0, int v1, int v2) const;

template<int ND> FLT mesh<ND>::area(int v0, int v1, int v2) const {
   FLT dx1,dy1,dx2,dy2;
   
   dx1 =  (vrtx[v0][0]-vrtx[v2][0]);
   dy1 =  (vrtx[v0][1]-vrtx[v2][1]);
   dx2 =  (vrtx[v1][0]-vrtx[v0][0]);
   dy2 =  (vrtx[v1][1]-vrtx[v0][1]);

   return(dx1*dy2 -dy1*dx2);
}

template FLT mesh<2>::area(int snum, int v2) const;
template FLT mesh<3>::area(int snum, int v2) const;

template<int ND> FLT mesh<ND>::area(int snum, int v2) const {
   FLT dx1,dy1,dx2,dy2;
   int v0, v1;
   
   v0 = svrtx[snum][0];
   v1 = svrtx[snum][1];
   
   dx1 =  (vrtx[v0][0]-vrtx[v2][0]);
   dy1 =  (vrtx[v0][1]-vrtx[v2][1]);
   dx2 =  (vrtx[v1][0]-vrtx[v0][0]);
   dy2 =  (vrtx[v1][1]-vrtx[v0][1]);

   return(dx1*dy2 -dy1*dx2);
}

template FLT mesh<2>::area(int tind) const;
template FLT mesh<3>::area(int tind) const;

template<int ND> FLT mesh<ND>::area(int tind) const {
   FLT dx1,dy1,dx2,dy2;
   int v0, v1, v2;
   
   v0 = tvrtx[tind][0];
   v1 = tvrtx[tind][1];
   v2 = tvrtx[tind][2];
   
   dx1 =  (vrtx[v0][0]-vrtx[v2][0]);
   dy1 =  (vrtx[v0][1]-vrtx[v2][1]);
   dx2 =  (vrtx[v1][0]-vrtx[v0][0]);
   dy2 =  (vrtx[v1][1]-vrtx[v0][1]);

   return(dx1*dy2 -dy1*dx2);
}

static FLT a[3];

template FLT mesh<2>::intri(int tind, FLT x[2]) const;

template<int ND> FLT mesh<ND>::intri(int tind, FLT x[ND]) const {
   int v0,v1,v2;
   FLT dx0,dy0,dx1,dy1,dx2,dy2;

   v0 = tvrtx[tind][0];
   v1 = tvrtx[tind][1];
   v2 = tvrtx[tind][2];
   
   dx0 =  (x[0] -vrtx[v0][0]);
   dy0 =  (x[1] -vrtx[v0][1]); 
   dx1 =  (x[0] -vrtx[v1][0]);
   dy1 =  (x[1] -vrtx[v1][1]);
   dx2 =  (x[0] -vrtx[v2][0]);
   dy2 =  (x[1] -vrtx[v2][1]);
   
   a[0] = (dy2*dx1 -dx2*dy1);
   a[1] = (dy0*dx2 -dx0*dy2);
   a[2] = (dy1*dx0 -dx1*dy0);
   
   return(fabs(a[0]) +fabs(a[1]) +fabs(a[2]) - (a[0] +a[1] +a[2]));
}

/* RETURNS WEIGHTS FROM INTRI FUNCTION */
template void mesh<2>::getwgts(FLT *wt) const;

template<int ND> void mesh<ND>::getwgts(FLT *wt) const {
   int i;
   FLT sum;
   
   sum = a[0] +a[1] +a[2];
   for(i=0;i<3;++i) 
      wt[i] = a[i]/sum;
   
   return;
}
      
template FLT mesh<2>::minangle(int v0, int v1, int v2) const;
template FLT mesh<3>::minangle(int v0, int v1, int v2) const;

template<int ND> FLT mesh<ND>::minangle(int v0, int v1, int v2) const {
   int i, i1, i2;
   FLT l[3];
   FLT dx[3], dy[3];
   FLT crossprod;
   
   dx[0] = vrtx[v2][0] -vrtx[v1][0];
   dy[0] = vrtx[v2][1] -vrtx[v1][1];
   l[0] = dx[0]*dx[0] +dy[0]*dy[0];   

   dx[1] = vrtx[v0][0] -vrtx[v2][0];
   dy[1] = vrtx[v0][1] -vrtx[v2][1];
   l[1] = dx[1]*dx[1] +dy[1]*dy[1];
   
   dx[2] = vrtx[v1][0] -vrtx[v0][0];
   dy[2] = vrtx[v1][1] -vrtx[v0][1];
   l[2] = dx[2]*dx[2] +dy[2]*dy[2];
      
   i = (l[0] < l[1] ? 0 : 1);
   i = (l[i] < l[2] ? i : 2);
   i1 = (i+1)%3;
   i2 = (i+2)%3;
   
   crossprod = -dx[i2]*dy[i1] +dy[i2]*dx[i1];

   return(crossprod/sqrt(l[i1]*l[i2]));
   
}
   
template FLT mesh<2>::angle(int v0, int v1, int v2) const;
template FLT mesh<3>::angle(int v0, int v1, int v2) const;

template<int ND> FLT mesh<ND>::angle(int v0, int v1, int v2) const {
   FLT l[3];
   FLT dx, dy;
   
   dx = vrtx[v1][0] -vrtx[v0][0];
   dy = vrtx[v1][1] -vrtx[v0][1];
   l[0] = dx*dx +dy*dy;

   dx = vrtx[v2][0] -vrtx[v1][0];
   dy = vrtx[v2][1] -vrtx[v1][1];
   l[1] = dx*dx +dy*dy;   

   dx = vrtx[v0][0] -vrtx[v2][0];
   dy = vrtx[v0][1] -vrtx[v2][1];
   l[2] = dx*dx +dy*dy;  
         
   return((l[0] +l[1] -l[2])/(2.*sqrt(l[0]*l[1])));
   
}

template FLT mesh<2>::tradius(int tind) const;
template FLT mesh<3>::tradius(int tind) const;

template<int ND> FLT mesh<ND>::tradius(int tind) const {
   FLT alpha,beta;
   FLT xmid1,ymid1,xmid2,ymid2,xcen,ycen;
   FLT dx1,dy1,dx2,dy2,area;
   int v0, v1, v2;
   
   v0 = tvrtx[tind][0];
   v1 = tvrtx[tind][1];
   v2 = tvrtx[tind][2];
   
   dx1 =  (vrtx[v0][0]-vrtx[v2][0]);
   dy1 =  (vrtx[v0][1]-vrtx[v2][1]);
   dx2 =  (vrtx[v1][0]-vrtx[v0][0]);
   dy2 =  (vrtx[v1][1]-vrtx[v0][1]);

   /* RELATIVE TO POINT 0 TO AVOID ROUNDOFF */
   xmid1 = 0.5*(vrtx[v2][0] -vrtx[v0][0]);
   ymid1 = 0.5*(vrtx[v2][1] -vrtx[v0][1]);   
   xmid2 = 0.5*(vrtx[v1][0] -vrtx[v0][0]);
   ymid2 = 0.5*(vrtx[v1][1] -vrtx[v0][1]);
      
   area        = 1.0/(dx1*dy2 -dy1*dx2);
   alpha       = dx2*xmid2 +dy2*ymid2;
   beta        = dx1*xmid1 +dy1*ymid1;
   xcen = area*(beta*dy2 -alpha*dy1);
   ycen = area*(alpha*dx1 -beta*dx2);

   return(sqrt(xcen*xcen +ycen*ycen));
}

template void mesh<2>::tcenter(int tind, FLT x[2]) const;
template void mesh<3>::tcenter(int tind, FLT x[3]) const;

template<int ND> void mesh<ND>::tcenter(int tind, FLT x[ND]) const {
   FLT alpha,beta;
   FLT xmid1,ymid1,xmid2,ymid2;
   FLT dx1,dy1,dx2,dy2,area;
   int v0, v1, v2;
   
   v0 = tvrtx[tind][0];
   v1 = tvrtx[tind][1];
   v2 = tvrtx[tind][2];
   
   dx1 =  (vrtx[v0][0]-vrtx[v2][0]);
   dy1 =  (vrtx[v0][1]-vrtx[v2][1]);
   dx2 =  (vrtx[v1][0]-vrtx[v0][0]);
   dy2 =  (vrtx[v1][1]-vrtx[v0][1]);

   /* RELATIVE TO POINT V0 */
   xmid1 = 0.5*(vrtx[v2][0] -vrtx[v0][0]);
   ymid1 = 0.5*(vrtx[v2][1] -vrtx[v0][1]);   
   xmid2 = 0.5*(vrtx[v1][0] -vrtx[v0][0]);
   ymid2 = 0.5*(vrtx[v1][1] -vrtx[v0][1]);
      
   area        = 1.0/(dx1*dy2 -dy1*dx2);
   alpha       = dx2*xmid2 +dy2*ymid2;
   beta        = dx1*xmid1 +dy1*ymid1;
   x[0] = area*(beta*dy2 -alpha*dy1) +vrtx[v0][0];
   x[1] = area*(alpha*dx1 -beta*dx2) +vrtx[v0][1];
   
   return;
}

template FLT mesh<2>::aspect(int tind) const;
template FLT mesh<3>::aspect(int tind) const;

template<int ND> FLT mesh<ND>::aspect(int tind) const {
   int v0,v1,v2;
   FLT dx1,dy1,dx2,dy2,area,perim;
   
   v0 = tvrtx[tind][0];
   v1 = tvrtx[tind][1];
   v2 = tvrtx[tind][2];
   
   dx1 =  (vrtx[v0][0]-vrtx[v2][0]);
   dy1 =  (vrtx[v0][1]-vrtx[v2][1]);
   dx2 =  (vrtx[v1][0]-vrtx[v0][0]);
   dy2 =  (vrtx[v1][1]-vrtx[v0][1]);
   
   area = (dx1*dy2 -dy1*dx2);
   perim = sqrt(dx1*dx1 +dy1*dy1) +sqrt(dx2*dx2 +dy2*dy2)
      +sqrt((dx1+dx2)*(dx1+dx2) +(dy1+dy2)*(dy1+dy2));
      
   return(area/(perim*perim));
}

