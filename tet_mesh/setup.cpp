/*
 *  setup.cpp
 *  mesh3d
 *
 *  Created by Michael Brazell on 7/30/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_mesh.h"
#include <utilities.h>
#include <float.h>
#include <math.h>


void tet_mesh::fixvertexinfo(void){
    /* orient elements so that v0 has the lowest global index
    and v1 has the second lowest global index.*/
    int tind,j,ind,temp;
    double rnx,rny,rnz,volume;
    TinyVector<int,4> v;
    
    for(tind = 0; tind < ntet; ++tind){
        
        v=tet(tind).pnt;
        // move v0 to first index
        ind=0;
        for(j = 1; j < 4; ++j){
            if (v(j) < v(ind)){
                ind = j;                
            }
        }
        
        temp = v(ind);
        v(ind) = v(0);
        v(0) = temp;
        
        // move v1 to second index
        ind = 1;
        for(j = 2; j < 4; ++j){
            if (v(j) < v(ind)){
                ind = j;                
            }
        }
        temp = v(ind);
        v(ind) = v(1);
        v(1) = temp;

        
        /* check orientation of v2 & v3    by calculating volume and checking sign */
        rnx = (pnts(v(1))(1)-pnts(v(0))(1))*(pnts(v(2))(2)-pnts(v(0))(2))-(pnts(v(2))(1)-pnts(v(0))(1))*(pnts(v(1))(2)-pnts(v(0))(2));
        rny = (pnts(v(1))(2)-pnts(v(0))(2))*(pnts(v(2))(0)-pnts(v(0))(0))-(pnts(v(2))(2)-pnts(v(0))(2))*(pnts(v(1))(0)-pnts(v(0))(0));
        rnz = (pnts(v(1))(0)-pnts(v(0))(0))*(pnts(v(2))(1)-pnts(v(0))(1))-(pnts(v(2))(0)-pnts(v(0))(0))*(pnts(v(1))(1)-pnts(v(0))(1));
        volume=-rnx*(pnts(v(3))(0)-pnts(v(0))(0))-rny*(pnts(v(3))(1)-pnts(v(0))(1))-rnz*(pnts(v(3))(2)-pnts(v(0))(2));
        
        //cout << "volume = " << volume << endl;
        // temporary may have to flip sign
        if(volume < 0.0){
            temp = v(2);
            v(2) = v(3);
            v(3) = temp;
        }
        tet(tind).pnt = v;
        
        /* calculates volume */
        tet(tind).vol=fabs(volume);            
            
    }

    return;
}

void tet_mesh::feedinvertexinfo(void){
    /* orient elements so that v0 has the lowest pnt.info index
    and v1 has the second lowest pnt.info index */
    int tind,j,ind,temp;
    double rnx,rny,rnz,volume;
    TinyVector<int,4> v,vinfo;
    
    for(tind = 0; tind < ntet; ++tind){
        
        v=tet(tind).pnt;
        
        for(j = 0; j < 4; ++j)
            vinfo(j)=pnt(v(j)).info;
            
        // move v0 to first index
        ind=0;
        for(j = 1; j < 4; ++j){
            if (vinfo(j) < vinfo(ind)){
                ind = j;                
            }
        }
        
        temp = v(ind);
        v(ind) = v(0);
        v(0) = temp;
        
        temp = vinfo(ind);
        vinfo(ind) = vinfo(0);
        vinfo(0) = temp;
        
        // move v1 to second index
        ind = 1;
        for(j = 2; j < 4; ++j){
            if (vinfo(j) < vinfo(ind)){
                ind = j;                
            }
        }
        temp = v(ind);
        v(ind) = v(1);
        v(1) = temp;

        
        /* check orientation of v2 & v3    by calculating volume and checking sign */
        rnx = (pnts(v(1))(1)-pnts(v(0))(1))*(pnts(v(2))(2)-pnts(v(0))(2))-(pnts(v(2))(1)-pnts(v(0))(1))*(pnts(v(1))(2)-pnts(v(0))(2));
        rny = (pnts(v(1))(2)-pnts(v(0))(2))*(pnts(v(2))(0)-pnts(v(0))(0))-(pnts(v(2))(2)-pnts(v(0))(2))*(pnts(v(1))(0)-pnts(v(0))(0));
        rnz = (pnts(v(1))(0)-pnts(v(0))(0))*(pnts(v(2))(1)-pnts(v(0))(1))-(pnts(v(2))(0)-pnts(v(0))(0))*(pnts(v(1))(1)-pnts(v(0))(1));
        volume=-rnx*(pnts(v(3))(0)-pnts(v(0))(0))-rny*(pnts(v(3))(1)-pnts(v(0))(1))-rnz*(pnts(v(3))(2)-pnts(v(0))(2));
        
        /* if negative flip vertex */
        if(volume < 0.0){
            temp = v(2);
            v(2) = v(3);
            v(3) = temp;
        }
        tet(tind).pnt = v;
        
        /* calculates volume */
        tet(tind).vol=fabs(volume);            
            
    }

    return;
}



void tet_mesh::vertexnnbor(void){
    /* counts the number of neigboring tet's connected to a pnt*/

    TinyVector<int,4> v;
    
    for(int i=0;i<npnt;++i){
        pnt(i).nnbor = 0;
        pnt(i).nspk = 0;
    }
    
    for(int tind = 0; tind < ntet; ++tind){    
        for(int j = 0; j < 4; ++j)
            v(j)=tet(tind).pnt(j);        
        for(int j = 0; j < 4; ++j){        
            ++pnt(v(j)).nnbor;    
            pnt(v(j)).tet=tind;
        }                
    }
    
    for(int eind = 0; eind < nseg; ++eind){    
        for(int j = 0; j < 2; ++j)
            v(j)=seg(eind).pnt(j);        
        for(int j = 0; j < 2; ++j){    
            ++pnt(v(j)).nspk;
            pnt(v(j)).seg=eind;
        }                    
    }

    return;
}


void tet_mesh::createedgeinfo(void) {
    int i,tind,minv,eind,eindprev,ne,lcl0;
    long lcl1,lcl2;
    TinyMatrix<int,6,2> vs;
    TinyVector<int,2> v,a;

    for(i=0;i<npnt;++i){
        pnt(i).info = -1;
    }
    for(i=0;i<nseg;++i){
        seg(i).nnbor = 0;
    }

    vs(0,0)=2, vs(0,1)=3; 
    vs(1,0)=3, vs(1,1)=1; 
    vs(2,0)=2, vs(2,1)=1;
    vs(3,0)=1, vs(3,1)=0;
    vs(4,0)=2, vs(4,1)=0; 
    vs(5,0)=3, vs(5,1)=0;


    ne=0;
    for(tind = 0; tind < ntet; ++tind){
        for(i = 0; i < 6; ++i){            
            v(0)=tet(tind).pnt(vs(i,0));
            v(1)=tet(tind).pnt(vs(i,1));            
            minv=min(v);
            eind=pnt(minv).info;
            /* search for edges already created */
            lcl0=v(0)+v(1);
            lcl1=(v(0)+1)*(v(1)+1);
            while(eind >= 0){                
                a(0)=seg(eind).pnt(0);    
                a(1)=seg(eind).pnt(1);            
                lcl2=fabs(lcl0-a(0)-a(1));
                lcl2+=fabs(lcl1-(a(0)+1)*(a(1)+1));
                if(lcl2 == 0) {
                    tet(tind).seg(i)=eind;   // 6 edges attached to tet
                    ++seg(eind).nnbor;        // number of neighbors surrounding an seg
                    seg(eind).tet = tind;       // one tet connected to an seg
                    if(v(0) == a(1))
                        tet(tind).sgn(i) = -1;
                    else
                        tet(tind).sgn(i) = 1;
                    goto NEXTEDGE; 
                }            
                eindprev=eind;
                eind=seg(eind).info;            
            }        

            /* create new seg    */        
            seg(ne).tet = tind;       // one tet connected to an seg
            tet(tind).seg(i)=ne;     // 6 edges connected to a tet
            ++seg(ne).nnbor;          // number of neighbors surrounding an seg
            seg(ne).pnt=v;           // two vertex attached to seg
            seg(ne).info=-1;
            tet(tind).sgn(i) = 1;
            if (pnt(minv).info < 0){
                pnt(minv).info=ne;
            }
            else {
                seg(eindprev).info=ne;
            }
            ne=ne+1;
NEXTEDGE:;

        }    
    }
    nseg = ne;

    return;
}


void tet_mesh::edgeinfo(void) {
    int i,tind,minv,eind,eindprev,lcl0;
    long lcl1,lcl2;
    TinyMatrix<int,6,2> vs;
    TinyVector<int,2> v,a;

    for(i=0;i<npnt;++i){
        pnt(i).info = -1;
    }
    for(i=0;i<nseg;++i){
        seg(i).nnbor = 0;
    }
    
    for(i=0;i<ntet;++i){        
        tet(i).sgn = 1;
    }
    
    vs(0,0)=2, vs(0,1)=3; 
    vs(1,0)=3, vs(1,1)=1; 
    vs(2,0)=2, vs(2,1)=1;
    vs(3,0)=1, vs(3,1)=0;
    vs(4,0)=2, vs(4,1)=0; 
    vs(5,0)=3, vs(5,1)=0;


    for(int j = 0; j < nseg; ++j){
        v(0) = seg(j).pnt(0);
        v(1) = seg(j).pnt(1);
        minv = min(v);
        eind = pnt(minv).info;
        while(eind >= 0){
            eindprev = eind;
            eind = seg(eind).info;        
        }
        if(pnt(minv).info < 0)
            pnt(minv).info = j;
        else
            seg(eindprev).info = j;
        seg(j).info = -1;    
    }
    
    for(tind = 0; tind < ntet; ++tind){
        for(i = 0; i < 6; ++i){            
            v(0)=tet(tind).pnt(vs(i,0));
            v(1)=tet(tind).pnt(vs(i,1));            
            minv=min(v);
            eind=pnt(minv).info;
            /* search for edges already created */
            lcl0=v(0)+v(1);
            lcl1=(v(0)+1)*(v(1)+1);
            while(eind >= 0){                
                a(0)=seg(eind).pnt(0);    
                a(1)=seg(eind).pnt(1);            
                lcl2=fabs(lcl0-a(0)-a(1));
                lcl2+=fabs(lcl1-(a(0)+1)*(a(1)+1));
                if(lcl2 == 0) {
                    tet(tind).seg(i)=eind;   // 6 edges attached to tet
                    ++seg(eind).nnbor;        // number of neighbors surrounding an seg
                    seg(eind).tet = tind;     // one tet connected to an seg
                    if(v(0) == a(1))
                        tet(tind).sgn(i) = -1;            
                    goto NEXTEDGE; 
                }            
                eind=seg(eind).info;            
            }        

NEXTEDGE:;

        }    
    }

    return;
}



void tet_mesh::createfaceinfo(void) {
    int i,j,tind,minv,find,findprev,nf,lcl0;
    long lcl1,lcl2;
    TinyMatrix<int,4,3> vf;
    TinyVector<int,3> v,a;

    vf(0,0)=1, vf(0,1)=2, vf(0,2)=3;
    vf(1,0)=0, vf(1,1)=3, vf(1,2)=2;
    vf(2,0)=0, vf(2,1)=1, vf(2,2)=3;
    vf(3,0)=0, vf(3,1)=2, vf(3,2)=1; 
    
    for(i=0;i<npnt;++i)
        pnt(i).info = -1;
    
    nf=0;
    for(tind = 0; tind < ntet; ++tind){
        for(i = 0; i < 4; ++i){
            for(j = 0; j < 3; ++j)
                v(j)=tet(tind).pnt(vf(i,j));
            minv=min(v);
            find=pnt(minv).info;
            //search for face already created
            lcl0=v(0)+v(1)+v(2);
            lcl1=(v(0)+1)*(v(1)+1)*(v(2)+1);
            while(find >= 0){        
                
                /* Proof by Michael Felland Clarkson University 8/1/07
                Face matching algorithm: The only way to match is if both sets are identical. 
                Let a be the common value--call the sets {a,b,c} and {a,d,e}.  You want a*b*c=a*d*e 
                and a+b+c=a+d+e; a is not 0.  This means b*c=d*e and b+c=d+e.  Since the sums are the 
                same, one of the two pairs of values must lie between the other--say b<d<e<c. If m
                is the common mean value, m=(b+c)/2=(d+e)/2, then there are x and y so that
                d=m-x,e=m+x,b=m-y,c=m+y with x<y.  d*e = m^2-x^2 > m^2-y^2 = b*c (this is if
                all are positive or all negative--cases where some are negative can be also handled). */
                
                a(0)=tri(find).pnt(0);    
                a(1)=tri(find).pnt(1);
                a(2)=tri(find).pnt(2);
                lcl2=fabs(lcl0-a(0)-a(1)-a(2));
                lcl2+=fabs(lcl1-(a(0)+1)*(a(1)+1)*(a(2)+1));
                if(lcl2 == 0 ){        
                    tri(find).tet(1)=tind;    // 2nd tet connected to a face
                    tet(tind).tri(i)=find;      // ith face on to a tet
                    tet(tind).rot(i) = -1;        // cw orientation 
                    goto NEXTFACE; 
                }            
                findprev=find;
                find=tri(find).info;            
            }        

            //create new face
            tri(nf).tet(0) = tind;       // 1st tet connected to a face
            tri(nf).tet(1) = -1;
            tet(tind).tri(i)=nf;        // 4 faces on a tet
            tri(nf).pnt=v;              // 3 vertices on a face
            tri(nf).info=-1; 
            tet(tind).rot(i) = 1;        // ccw standard orientation
            if (pnt(minv).info < 0) {
                pnt(minv).info=nf;
            }
            else {
                tri(findprev).info=nf;
            }
    
            nf=nf+1;
NEXTFACE:;
            
    
        }
            
    }
    ntri = nf;
    return;
}

void tet_mesh::faceinfo(void) {
    int i,j,tind,minv,find,findprev,lcl0;
    long lcl1,lcl2;
    TinyMatrix<int,4,3> vf;
    TinyVector<int,3> v,a;

    // may need rearrangment to correspond to my standard element
    vf(0,0)=1, vf(0,1)=2, vf(0,2)=3;
    vf(1,0)=0, vf(1,1)=3, vf(1,2)=2;
    vf(2,0)=0, vf(2,1)=1, vf(2,2)=3;
    vf(3,0)=0, vf(3,1)=2, vf(3,2)=1;  
    
    for(i=0;i<npnt;++i)
        pnt(i).info = -1;
    for(i=0;i<ntri;++i){
        tri(i).tet(0) = -1;
        tri(i).tet(1) = -1;
    }
        
    for(int j = 0; j < ntri; ++j){
        v(0) = tri(j).pnt(0);
        v(1) = tri(j).pnt(1);
        v(2) = tri(j).pnt(2);
        minv = min(v);
        find = pnt(minv).info;
        while(find >= 0){
            findprev = find;
            find = tri(find).info;        
        }
        if(pnt(minv).info < 0)
            pnt(minv).info = j;
        else
            tri(findprev).info = j;
        tri(j).info = -1;    
    }
    
    for(tind = 0; tind < ntet; ++tind){
        for(i = 0; i < 4; ++i){
            for(j = 0; j < 3; ++j)
                v(j)=tet(tind).pnt(vf(i,j));            
            minv=min(v);
            find=pnt(minv).info;
            //search for face already created
            lcl0=v(0)+v(1)+v(2);
            lcl1=(v(0)+1)*(v(1)+1)*(v(2)+1);
            while(find >= 0){        
                            
                /* Proof by Michael Felland Clarkson University 8/1/07
                Face matching algorithm: The only way to match is if both sets are identical. 
                Let a be the common value--call the sets {a,b,c} and {a,d,e}.  You want a*b*c=a*d*e 
                and a+b+c=a+d+e; a is not 0.  This means b*c=d*e and b+c=d+e.  Since the sums are the 
                same, one of the two pairs of values must lie between the other--say b<d<e<c. If m
                is the common mean value, m=(b+c)/2=(d+e)/2, then there are x and y so that
                d=m-x,e=m+x,b=m-y,c=m+y with x<y.  d*e = m^2-x^2 > m^2-y^2 = b*c (this is if
                all are positive or all negative--cases where some are negative can be also handled). */
                
                a(0)=tri(find).pnt(0);    
                a(1)=tri(find).pnt(1);
                a(2)=tri(find).pnt(2);
                lcl2=fabs(lcl0-a(0)-a(1)-a(2));
                lcl2+=fabs(lcl1-(a(0)+1)*(a(1)+1)*(a(2)+1));
                if(lcl2 == 0 ){        
                    if(a(1) == v(2))
                        tet(tind).rot(i) = -1;        // cw orientation            
                    else 
                        tet(tind).rot(i) = 1;        // ccw orientation
                                            
                    if (tri(find).tet(0) >= 0)
                        tri(find).tet(1)=tind;    // 2nd tet connected to a face
                    else 
                        tri(find).tet(0)=tind;     // 1st tet connected to a face
                    
                    tet(tind).tri(i)=find;          // 4 face on to a tet
                    goto NEXTFACE;
                }            
                find=tri(find).info;            
            }            

NEXTFACE:;
    
        }    
    }    
    
    return;
}

void tet_mesh::morefaceinfo(void){
    int i,tind,find;
    TinyMatrix<int,4,3> sf;
    TinyMatrix<int,3,2> vs;
    TinyVector<int,3> v;

    sf(0,0)=0, sf(0,1)=1, sf(0,2)=2;
    sf(1,0)=0, sf(1,1)=4, sf(1,2)=5;
    sf(2,0)=1, sf(2,1)=5, sf(2,2)=3;
    sf(3,0)=2, sf(3,1)=3, sf(3,2)=4;

    vs(0,0)=1, vs(0,1)=2;
    vs(1,0)=2, vs(1,1)=0;
    vs(2,0)=0, vs(2,1)=1;
    
    int findtri=-1;
    int ind;
    
    for(find = 0; find < ntri; ++find){
        tind = tri(find).tet(0);
        ind = 0;
        for(i = 0; i < 4; ++i){
            if(tet(tind).tri(i) == find){
                findtri = ind;
                //break;
            }
            ++ind;
        }
        for(i = 0; i < 3; ++i){
            tri(find).seg(i)=tet(tind).seg(sf(findtri,i));    // 3 edges on a face
        }
        v=tri(find).pnt;
        for(int i = 0; i < 3; ++i){
            if(v(vs(i,0)) == seg(tri(find).seg(i)).pnt(0) && v(vs(i,1)) == seg(tri(find).seg(i)).pnt(1)){
                tri(find).sgn(i) = 1;
            }
            else if (v(vs(i,1)) == seg(tri(find).seg(i)).pnt(0) && v(vs(i,0)) == seg(tri(find).seg(i)).pnt(1)){
                tri(find).sgn(i) = -1;
            }
            else{
                cout << "bad sign on tri: " << find << endl;
            }
        }
    }
    return;
}


void tet_mesh::createfaceorientation(void) {
    int i,j,tind,find;
    TinyMatrix<int,4,3> vf;
    TinyVector<int,3> v,a;

    // may need rearrangment to correspond to my standard element
    vf(0,0)=1, vf(0,1)=2, vf(0,2)=3;
    vf(1,0)=0, vf(1,1)=3, vf(1,2)=2;
    vf(2,0)=0, vf(2,1)=1, vf(2,2)=3;
    vf(3,0)=0, vf(3,1)=2, vf(3,2)=1;  
    
    for(tind = 0; tind < ntet; ++tind){
        for(i = 0; i < 4; ++i){
            tet(tind).rot(i) = 1;
            find = tet(tind).tri(i);
            for(j = 0; j < 3; ++j){
                v(j)=tet(tind).pnt(vf(i,j));
                a(j)=tri(find).pnt(j);
            }
            
            if(a(0) == v(0)){
                if(a(1) == v(2)){
                    tet(tind).rot(i) = -1;
                }
            }
            
            if(a(0) == v(1)){
                if(a(1) == v(2)){
                    tet(tind).rot(i) = 2;
                }
                else {
                    tet(tind).rot(i) = -2;
                }
            }
            
            if(a(0) == v(2)){
                if(a(1) == v(0)){
                    tet(tind).rot(i) = 3;
                }
                else {
                    tet(tind).rot(i) = -3;
                }
            }

        }
        cout << tet(tind).rot << endl;

    }
    
    return;
}


void tet_mesh::createtetinfo(void) {
    int tind,find,j;
//    int nbf = 0;
    for(tind = 0; tind < ntet; ++tind) {
        for(j = 0; j < 4; ++j){
            tet(tind).tet(j)=-1;
            find=tet(tind).tri(j);
            if(tri(find).tet(0) == tind) {
                tet(tind).tet(j)=tri(find).tet(1);
            }
            else {
                tet(tind).tet(j)=tri(find).tet(0);
            }
        }        
    }
    return;
}

void tet_mesh::treeinit() {
    int i,j,n,v0;
    FLT x1[ND], x2[ND], dx;
    
    for(n=0;n<ND;++n)    {
        x1[n] = pnts(0)(n);
        x2[n] = pnts(0)(n);
    }

    
    for (i=0;i<nfbd;++i) {
        for(j=0;j<fbdry(i)->npnt;++j) {
            v0 = fbdry(i)->pnt(j).gindx;
            for(n=0;n<ND;++n) {
                x1[n] = MIN(x1[n],pnts(v0)(n));
                x2[n] = MAX(x2[n],pnts(v0)(n));
            }
        }
    }
    
    for(n=0;n<ND;++n) {
        dx = MAX(x2[n]-x1[n],100.0*EPSILON);
        x1[n] -= 0.25*dx;
        x2[n] += 0.25*dx;
    }
    
    treeinit(x1,x2);

    return;
}

void tet_mesh::treeinit(FLT x1[ND], FLT x2[ND]) {
    
    otree.init(x1,x2);
        
    for(int i=0;i<npnt;++i) 
        otree.addpt(i);
    
    return;
}

void tet_mesh::initlngth() {
    int i,p0,p1;
    FLT l;
    
    for(i=0;i<npnt;++i)
        lngth(i) = 0.0;
        
    for(i=0;i<nseg;++i) {
        p0 = seg(i).pnt(0);
        p1 = seg(i).pnt(1);
        l = distance(seg(i).pnt(0),seg(i).pnt(1));
        lngth(p0) += l;
        lngth(p1) += l;
    }
    
    for(i=0;i<npnt;++i)
        lngth(i) /= pnt(i).nspk;

    return;
}

/* FIX tri.tet and tet.tet TO POINT TO GROUP/SIDE ON BOUNDARY */
void tet_mesh::bdrylabel() {
    int i,j,k,tind,ttind;
    
    for(i=0;i<nfbd;++i) {
        for(j=0;j<fbdry(i)->ntri;++j) {
            tind = fbdry(i)->tri(j).gindx;
            tri(tind).tet(1) = tetnumatbdry(i,j);
            ttind = tri(tind).tet(0);
            for(k=0;k<4;++k)
                if (tet(ttind).tri(k) == tind) break;
                
            tet(ttind).tet(k) = tri(tind).tet(1);
        }
    }
    
    return;
}

//
//void tet_mesh::settrim() {
//    int i,j,n,bsd,tin,tind,nsrch,ntdel;
//    
//    /* ASSUMES fscr1 HAS BEEN SET WITH VALUES TO DETERMINE HOW MUCH TO TRIM OFF OF BOUNDARIES */
//    
//    for(i=0;i<ntri;++i)
//        tet(i).info = 0;
//
//    ntdel = 0;
//
//    for (bsd=0;bsd<sbdry(0)->nseg;++bsd) {
//        tind = seg(sbdry(0)->el(bsd)).tri(0);
//        if (tet(tind).info > 0) continue;
//        
//        i1wk(0) = tind;
//        tet(tind).info = 1;
//        nsrch = ntdel+1;
//        
//        /* NEED TO SEARCH SURROUNDING TRIANGLES */
//        for(i=ntdel;i<nsrch;++i) {
//            tin = i1wk(i);
//            for (n=0;n<3;++n)
//                if (fscr1(tet(tin).pnt(n)) < 0.0) goto NEXT;
//                
//            i1wk(ntdel++) = tin;
//
//            for(j=0;j<3;++j) {
//                tind = tet(tin).tri(j);
//                if (tind < 0) continue;
//                if (tet(tind).info > 0) continue; 
//                tet(tind).info = 1;          
//                i1wk(nsrch++) = tind;
//            }
//            NEXT: continue;
//        }
//    }
//    
//    for(i=0;i<ntri;++i)
//        tet(i).info = 0;
//        
//    for(i=0;i<ntdel;++i)
//        tet(i1wk(i)).info = 1;
//        
//    for(i=0;i<maxvst;++i)
//        i1wk(i) = -1;
//        
//    return;
//}


