#include"tri_mesh.h"
#include<cmath>

FLT tri_mesh::incircle(int tind, const TinyVector<FLT,ND> &a) const {
    int i;
    TinyMatrix<FLT,3,2> pt;
    TinyVector<FLT,3> l2;
    FLT determ;
    
    for(i=0;i<3;++i) {
        pt(i,0) = vrtx(td(tind).vrtx(i))(0)-a(0);
        pt(i,1) = vrtx(td(tind).vrtx(i))(1)-a(1);
    }

    for(i=0;i<3;++i)
        l2(i) = pt(i,0)*pt(i,0) +pt(i,1)*pt(i,1);

    determ = 0.0;
    determ +=  pt(0,0)*(pt(1,1)*l2(2) -pt(2,1)*l2(1));
    determ += -pt(0,1)*(pt(1,0)*l2(2) -pt(2,0)*l2(1));
    determ +=  l2(0)*(pt(1,0)*pt(2,1) -pt(1,1)*pt(2,0));
    
    return(determ);
}

FLT tri_mesh::insidecircle(int sind, const TinyVector<FLT,ND> &a) const {
    int v0,v1;
    TinyVector<FLT,2> ctr;
    FLT dist2;
    
    v0 = sd(sind).vrtx(0);
    v1 = sd(sind).vrtx(1);
    ctr(0) = 0.5*(vrtx(v0)(0) +vrtx(v1)(0));
    ctr(1) = 0.5*(vrtx(v0)(1) +vrtx(v1)(1));
    dist2 = (a(0)-ctr(0))*(a(0)-ctr(0)) +(a(1)-ctr(1))*(a(1)-ctr(1));
    return(0.25*distance2(v0,v1) -dist2);
}


FLT tri_mesh::area(int v0, int v1, int v2) const {
    FLT dx1,dy1,dx2,dy2;
    
    dx1 =  (vrtx(v0)(0)-vrtx(v2)(0));
    dy1 =  (vrtx(v0)(1)-vrtx(v2)(1));
    dx2 =  (vrtx(v1)(0)-vrtx(v0)(0));
    dy2 =  (vrtx(v1)(1)-vrtx(v0)(1));

    return(dx1*dy2 -dy1*dx2);
}

FLT tri_mesh::area(int snum, int v2) const {
    FLT dx1,dy1,dx2,dy2;
    int v0, v1;
    
    v0 = sd(snum).vrtx(0);
    v1 = sd(snum).vrtx(1);
    
    dx1 =  (vrtx(v0)(0)-vrtx(v2)(0));
    dy1 =  (vrtx(v0)(1)-vrtx(v2)(1));
    dx2 =  (vrtx(v1)(0)-vrtx(v0)(0));
    dy2 =  (vrtx(v1)(1)-vrtx(v0)(1));

    return(dx1*dy2 -dy1*dx2);
}

FLT tri_mesh::area(int tind) const {
    FLT dx1,dy1,dx2,dy2;
    int v0, v1, v2;
    
    v0 = td(tind).vrtx(0);
    v1 = td(tind).vrtx(1);
    v2 = td(tind).vrtx(2);
    
    dx1 =  (vrtx(v0)(0)-vrtx(v2)(0));
    dy1 =  (vrtx(v0)(1)-vrtx(v2)(1));
    dx2 =  (vrtx(v1)(0)-vrtx(v0)(0));
    dy2 =  (vrtx(v1)(1)-vrtx(v0)(1));

    return(dx1*dy2 -dy1*dx2);
}

FLT tri_mesh::intri(int tind, const TinyVector<FLT,2> &x) {
    int v0,v1,v2;
    FLT dx0,dy0,dx1,dy1,dx2,dy2;

    v0 = td(tind).vrtx(0);
    v1 = td(tind).vrtx(1);
    v2 = td(tind).vrtx(2);
    
    dx0 =  (x(0) -vrtx(v0)(0));
    dy0 =  (x(1) -vrtx(v0)(1)); 
    dx1 =  (x(0) -vrtx(v1)(0));
    dy1 =  (x(1) -vrtx(v1)(1));
    dx2 =  (x(0) -vrtx(v2)(0));
    dy2 =  (x(1) -vrtx(v2)(1));
    
    tri_wgt(0) = (dy2*dx1 -dx2*dy1);
    tri_wgt(1) = (dy0*dx2 -dx0*dy2);
    tri_wgt(2) = (dy1*dx0 -dx1*dy0);
    
    return(fabs(tri_wgt(0)) +fabs(tri_wgt(1)) +fabs(tri_wgt(2)) - (tri_wgt(0) +tri_wgt(1) +tri_wgt(2)));
}

/* RETURNS WEIGHTS FROM INTRI FUNCTION */
void tri_mesh::getwgts(TinyVector<FLT,3> &wt) const {
    int i;
    FLT sum;
    
    sum = tri_wgt(0) +tri_wgt(1) +tri_wgt(2);
    for(i=0;i<3;++i) 
        wt(i) = tri_wgt(i)/sum;
    
    return;
}
        
FLT tri_mesh::minangle(int v0, int v1, int v2) const {
    int i, i1, i2;
    TinyVector<FLT,3> l,dx,dy;
    FLT crossprod;
    
    dx(0) = vrtx(v2)(0) -vrtx(v1)(0);
    dy(0) = vrtx(v2)(1) -vrtx(v1)(1);
    l(0) = dx(0)*dx(0) +dy(0)*dy(0);    

    dx(1) = vrtx(v0)(0) -vrtx(v2)(0);
    dy(1) = vrtx(v0)(1) -vrtx(v2)(1);
    l(1) = dx(1)*dx(1) +dy(1)*dy(1);
    
    dx(2) = vrtx(v1)(0) -vrtx(v0)(0);
    dy(2) = vrtx(v1)(1) -vrtx(v0)(1);
    l(2) = dx(2)*dx(2) +dy(2)*dy(2);
        
    i = (l(0) < l(1) ? 0 : 1);
    i = (l(i) < l(2) ? i : 2);
    i1 = (i+1)%3;
    i2 = (i+2)%3;
    
    crossprod = -dx(i2)*dy(i1) +dy(i2)*dx(i1);

    return(crossprod/sqrt(l(i1)*l(i2)));
    
}
    
FLT tri_mesh::angle(int v0, int v1, int v2) const {
    TinyVector<FLT,3> l;
    FLT dx, dy;
    
    dx = vrtx(v1)(0) -vrtx(v0)(0);
    dy = vrtx(v1)(1) -vrtx(v0)(1);
    l(0) = dx*dx +dy*dy;

    dx = vrtx(v2)(0) -vrtx(v1)(0);
    dy = vrtx(v2)(1) -vrtx(v1)(1);
    l(1) = dx*dx +dy*dy;    

    dx = vrtx(v0)(0) -vrtx(v2)(0);
    dy = vrtx(v0)(1) -vrtx(v2)(1);
    l(2) = dx*dx +dy*dy;  
            
    return((l(0) +l(1) -l(2))/(2.*sqrt(l(0)*l(1))));
    
}

FLT tri_mesh::circumradius(int tind) const {
    FLT alpha,beta;
    FLT xmid1,ymid1,xmid2,ymid2,xcen,ycen;
    FLT dx1,dy1,dx2,dy2,area;
    int v0, v1, v2;
    
    v0 = td(tind).vrtx(0);
    v1 = td(tind).vrtx(1);
    v2 = td(tind).vrtx(2);
    
    dx1 =  (vrtx(v0)(0)-vrtx(v2)(0));
    dy1 =  (vrtx(v0)(1)-vrtx(v2)(1));
    dx2 =  (vrtx(v1)(0)-vrtx(v0)(0));
    dy2 =  (vrtx(v1)(1)-vrtx(v0)(1));

    /* RELATIVE TO POINT 0 TO AVOID ROUNDOFF */
    xmid1 = 0.5*(vrtx(v2)(0) -vrtx(v0)(0));
    ymid1 = 0.5*(vrtx(v2)(1) -vrtx(v0)(1));    
    xmid2 = 0.5*(vrtx(v1)(0) -vrtx(v0)(0));
    ymid2 = 0.5*(vrtx(v1)(1) -vrtx(v0)(1));
        
    area          = 1.0/(dx1*dy2 -dy1*dx2);
    alpha         = dx2*xmid2 +dy2*ymid2;
    beta          = dx1*xmid1 +dy1*ymid1;
    xcen = area*(beta*dy2 -alpha*dy1);
    ycen = area*(alpha*dx1 -beta*dx2);

    return(sqrt(xcen*xcen +ycen*ycen));
}

void tri_mesh::circumcenter(int tind, TinyVector<FLT,2> &x) const {
    FLT alpha,beta;
    FLT xmid1,ymid1,xmid2,ymid2;
    FLT dx1,dy1,dx2,dy2,area;
    int v0, v1, v2;
    
    v0 = td(tind).vrtx(0);
    v1 = td(tind).vrtx(1);
    v2 = td(tind).vrtx(2);
    
    dx1 =  (vrtx(v0)(0)-vrtx(v2)(0));
    dy1 =  (vrtx(v0)(1)-vrtx(v2)(1));
    dx2 =  (vrtx(v1)(0)-vrtx(v0)(0));
    dy2 =  (vrtx(v1)(1)-vrtx(v0)(1));

    /* RELATIVE TO POINT V0 */
    xmid1 = 0.5*(vrtx(v2)(0) -vrtx(v0)(0));
    ymid1 = 0.5*(vrtx(v2)(1) -vrtx(v0)(1));    
    xmid2 = 0.5*(vrtx(v1)(0) -vrtx(v0)(0));
    ymid2 = 0.5*(vrtx(v1)(1) -vrtx(v0)(1));
        
    area          = 1.0/(dx1*dy2 -dy1*dx2);
    alpha         = dx2*xmid2 +dy2*ymid2;
    beta          = dx1*xmid1 +dy1*ymid1;
    x(0) = area*(beta*dy2 -alpha*dy1) +vrtx(v0)(0);
    x(1) = area*(alpha*dx1 -beta*dx2) +vrtx(v0)(1);
    
    return;
}

FLT tri_mesh::inscribedradius(int tind) const {
    int v0,v1,v2;
    FLT dx1,dy1,dx2,dy2,area,perim;
    
    v0 = td(tind).vrtx(0);
    v1 = td(tind).vrtx(1);
    v2 = td(tind).vrtx(2);
    
    dx1 =  (vrtx(v0)(0)-vrtx(v2)(0));
    dy1 =  (vrtx(v0)(1)-vrtx(v2)(1));
    dx2 =  (vrtx(v1)(0)-vrtx(v0)(0));
    dy2 =  (vrtx(v1)(1)-vrtx(v0)(1));
    
    area = (dx1*dy2 -dy1*dx2);
    perim = sqrt(dx1*dx1 +dy1*dy1) +sqrt(dx2*dx2 +dy2*dy2)
        +sqrt((dx1+dx2)*(dx1+dx2) +(dy1+dy2)*(dy1+dy2));
        
    return(area/perim);
}


FLT tri_mesh::aspect(int tind) const {
    int v0,v1,v2;
    FLT dx1,dy1,dx2,dy2,area,perim;
    
    v0 = td(tind).vrtx(0);
    v1 = td(tind).vrtx(1);
    v2 = td(tind).vrtx(2);
    
    dx1 =  (vrtx(v0)(0)-vrtx(v2)(0));
    dy1 =  (vrtx(v0)(1)-vrtx(v2)(1));
    dx2 =  (vrtx(v1)(0)-vrtx(v0)(0));
    dy2 =  (vrtx(v1)(1)-vrtx(v0)(1));
    
    area = (dx1*dy2 -dy1*dx2);
    perim = sqrt(dx1*dx1 +dy1*dy1) +sqrt(dx2*dx2 +dy2*dy2)
        +sqrt((dx1+dx2)*(dx1+dx2) +(dy1+dy2)*(dy1+dy2));
        
    return(area/(perim*perim));
}

