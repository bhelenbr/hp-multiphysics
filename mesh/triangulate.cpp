/*
 *  triangulate.cpp
 *  mblock
 *
 *  Created by helenbrk on Mon Aug 13 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"
#include <utilities.h>
#include<float.h>
#include<assert.h>

void tri_mesh::triangulate(int nsd) {
    int i,j,n,vcnt,sck,dirck,stest,nv,vtry;
    int sind2;
    int bgn,end;
    int sind,dir;
    int sind1,sindprev;
    int nsidebefore,ntest;
    int ngood;
    int minv,maxv,itemp;
    TinyVector<int,2> v;
    int v2,v3;
    TinyVector<FLT,2> xmid,xcen;
    FLT xm1dx1,xm2dx2;
    FLT hmin,height;
    FLT temp;
    FLT ds1,ds2;
    TinyVector<FLT,ND> xmax,xmin,xmax1,xmin1;
    TinyVector<FLT,ND> dx1,dx2,dx3,dx4;
    FLT det,det13,det14,det34,h1,h3;
    
    /* CREATE VERTEX LIST */
    nv = 0;
    for(i=0;i<nsd;++i) {
        sind = abs(gbl->i2wk_lst1(i)) -1;
        dir = (1 -SIGN(gbl->i2wk_lst1(i)))/2;
        sd(sind).tri(dir) = -1;
        gbl->i2wk_lst2(nv++) = sd(sind).vrtx(dir);
    }
    if (nv > maxvst -2) {
            *gbl->log << gbl->idprefix << " coarse mesh is not big enough " << nv << ' ' << maxvst << std::endl;
            exit(1);
    }
    
    /* SETUP SIDE POINTER INFO */
    for(i=0;i<nv;++i)
        vd(gbl->i2wk_lst2(i)).info = -1;
        
    for(i=0;i<nsd;++i) {
        sind = abs(gbl->i2wk_lst1(i)) -1;
        v(0) = sd(sind).vrtx(0);
        v(1) = sd(sind).vrtx(1);
        if (v(1) > v(0)) {
            minv = v(0);
            maxv = v(1);
        }
        else {
            minv = v(1);
            maxv = v(0);
        }
        sind1 = vd(minv).info;
        while (sind1 >= 0) {
            sindprev = sind1;
            sind1 = sd(sind1).info;
        }
        sd(sind).info = -1;
        if (vd(minv).info < 0)
            vd(minv).info = sind;
        else 
            sd(sindprev).info = sind;
    }
    
    bgn = 0;
    end = nsd;
    ntest = end;
    while(bgn < end) {
        nsidebefore = nside;
        for(sind2=bgn;sind2<end;++sind2) {
            sind = abs(gbl->i2wk_lst1(sind2)) -1;
            dir =  (1 -SIGN(gbl->i2wk_lst1(sind2)))/2;
            
            if (sd(sind).tri(dir) > -1) continue; // SIDE HAS ALREADY BEEN MATCHED 
            
            v(0) = sd(sind).vrtx(dir);
            v(1) = sd(sind).vrtx(1 -dir);

            /* SEARCH FOR GOOD POINTS */
            for(n=0;n<ND;++n) {
                dx2(n) = vrtx(v(1))(n) -vrtx(v(0))(n);
                xmid(n) = 0.5*(vrtx(v(1))(n) -vrtx(v(0))(n));
            }
            hmin = 1.0e99;
            
            /* FIND NODES WHICH MAKE POSITIVE TRIANGLE WITH SIDE */
            for(i=0;i<nv;++i) {
                vtry = gbl->i2wk_lst2(i);
                if (vtry == v(0) || vtry == v(1)) continue;
                
        
                for(n=0;n<ND;++n)
                    dx1(n) = vrtx(v(0))(n) -vrtx(vtry)(n);
                det          = dx1(0)*dx2(1) -dx1(1)*dx2(0);
                if (det <= 0.0) continue;
                
                /* CIRCUMCENTER IS AT INTERSECTION OF NORMAL TO SIDES THROUGH MIDPOINT */
                det = 1./det;
                xm1dx1         = -0.5*(dx1(0)*dx1(0) +dx1(1)*dx1(1));
                xm2dx2         =  0.5*(dx2(0)*dx2(0) +dx2(1)*dx2(1));
                xcen(0) = det*(xm1dx1*dx2(1) -xm2dx2*dx1(1));
                xcen(1) = det*(xm2dx2*dx1(0) -xm1dx1*dx2(0));
                        
                /* FIND TRIANGLE FOR WHICH THE HEIGHT OF THE CIRCUMCENTER */
                /* ABOVE THE EDGE MID-POINT IS MINIMIZED (MINIMIZES RADIUS) */
                height = dx2(0)*(xcen(1) -xmid(1)) -dx2(1)*(xcen(0) -xmid(0));
                
                if (height > hmin) continue;
    
                /* CHECK FOR INTERSECTION OF TWO CREATED SIDES */
                /* WITH ALL OTHER BOUNDARY SIDES */
                for(vcnt=0;vcnt<2;++vcnt) {
                    minv = MIN(vtry,v(vcnt));
                    maxv = MAX(vtry,v(vcnt));
                    
                    /* LOOK THROUGH ALL SIDES CONNECTED TO MINV FOR DUPLICATE */
                    /* IF DUPLICATE THEN SIDE IS OK - NO NEED TO CHECK */
                    sind1 = vd(minv).info;
                    while (sind1 >= 0) {
                        if (maxv == sd(sind1).vrtx(0) || maxv == sd(sind1).vrtx(1)) {
                            goto next_vrt;
                        }
                        sind1 = sd(sind1).info;
                    }
                    
                    /* FIND BOUNDING BOX OF POTENTIAL SIDE */
                    for(n=0;n<ND;++n) {
                        xmin(n)        = MIN(vrtx(vtry)(n),vrtx(v(vcnt))(n));
                        xmax(n)        = MAX(vrtx(vtry)(n),vrtx(v(vcnt))(n));
                    }
                    
                    for(sck=0;sck<ntest;++sck) {
                        stest = abs(gbl->i2wk_lst1(sck))-1;
                        if (stest == sind) continue;
                        dirck =  (1 -SIGN(gbl->i2wk_lst1(sck)))/2;
                        v2 = sd(stest).vrtx(dirck);
                        v3 = sd(stest).vrtx(1-dirck);
                        
                        /* NO NEED TO CHECK FOR INTERSECTIONS IF CONNECTED TO ENDPOINT */
                        if (v2 == minv || v3 == minv) continue;
                        if (v2 == maxv || v3 == maxv) continue;
                        
                        /* FIND BOUNDING BOX OF SIDE TO CHECK AGAINST */
                        for(n=0;n<ND;++n) {
                            xmin1(n)        = MIN(vrtx(v2)(n),vrtx(v3)(n));
                            xmax1(n)        = MAX(vrtx(v2)(n),vrtx(v3)(n));
                        }
                        /* IF BOUNDING BOXES DON'T OVERLAP THEN NO INTERSECTION */
                        for(n=0;n<ND;++n)
                            if (xmax(n) < xmin1(n) || xmin(n) > xmax1(n)) goto next_bdry_side;
                        
                        /* CHECK FOR INTERSECTION OF SIDES */
                        for(n=0;n<ND;++n) {
                            dx1(n) = vrtx(maxv)(n) -vrtx(minv)(n);
                            dx3(n) = vrtx(v3)(n)-vrtx(v2)(n);
                            dx4(n) = vrtx(v2)(n)-vrtx(minv)(n);
                        }
                        /* DETERMINANT IS POSITIVE IF CCW FROM VECTOR DX1 TO DX3 */
                        /* AREA OF TRIANGLE FORMED BY JOINING BASE OF VECTORS */
                        det13 = dx1(0)*dx3(1) -dx1(1)*dx3(0);
                        /* CHECK FOR PARALLELISM */
                        if (fabs(det13) < EPSILON*100.0*(fabs(xmax(0))+fabs(xmax(1)))) continue;
                        
                        det14 = dx1(0)*dx4(1) -dx4(0)*dx1(1);
                        det34 = dx3(0)*dx4(1) -dx4(0)*dx3(1);
                      
                        det13 = 1./det13;
                        /* Height ratio relative to side 1 frame of reference*/
                        h1 = det14*det13;
                        /* Height ratio relative to side 3 frame of reference */
                        h3 = det34*det13;
             
                        /* FIRST & SECOND PART CHECKS WHETHER HEIGHTS ARE IN SAME DIRECTION */
                        /* SECOND & THIRD CHECKS WHETHER HEIGHT MOVING AWAY IS GREATER THAN RETURNING */     
                        if (h1 > 0.0 || h3 > 0.0 || h1 < -1.0 || h3 < -1.0) continue;
                        
                        goto vtry_failed;
                        
next_bdry_side:    continue;
                    }
next_vrt:        continue;
                }
        
                /* CHECK IF DEGENERATE */
                if (height > hmin-200.0*EPSILON*sqrt(xcen(0)*xcen(0)+xcen(1)*xcen(1))) {
                    gbl->i2wk_lst3(ngood++) = vtry;
                    continue;
                }
                
                ngood = 0;
                gbl->i2wk_lst3(ngood++) = vtry;
                hmin = height+100.0*EPSILON*sqrt(xcen(0)*xcen(0)+xcen(1)*xcen(1));
vtry_failed:continue;
            }
            
            if (ngood > 1) {
                /* ORDER COCIRCULAR POINTS */    
                /* CALCULATE SIDE ANGLE */
                ds2 = 1./sqrt(dx2(0)*dx2(0) +dx2(1)*dx2(1));
                for(i=0;i<ngood;++i) {
                    vtry = gbl->i2wk_lst3(i);
                    dx1(0) = vrtx(v(0))(0) -vrtx(vtry)(0);
                    dx1(1) = vrtx(v(0))(1) -vrtx(vtry)(1);
                    ds1 = 1./sqrt(dx1(0)*dx1(0) +dx1(1)*dx1(1));
                    gbl->fltwk(i) = -(dx2(0)*dx1(0)  +dx2(1)*dx1(1))*ds2*ds1;
                }
            
                /* ORDER POINTS BY ANGLE */        
                for(i=0;i<ngood-1;++i) {
                    for(j=i+1;j<ngood;++j) {
        
                        /* TO ELIMINATE POSSIBILITY OF REPEATED VERTICES IN gbl->i2wk_lst2 */
                        if (gbl->i2wk_lst3(i) == gbl->i2wk_lst3(j)) {
                            gbl->i2wk_lst3(j) = gbl->i2wk_lst3(ngood-1);
                            --ngood;
                        }
        
                        /* ORDER BY ANGLE */
                        if(gbl->fltwk(i) > gbl->fltwk(j)) {
                            temp = gbl->fltwk(i);
                            gbl->fltwk(i) = gbl->fltwk(j);
                            gbl->fltwk(j) = temp;
                            itemp = gbl->i2wk_lst3(i);
                            gbl->i2wk_lst3(i) = gbl->i2wk_lst3(j);
                            gbl->i2wk_lst3(j) = itemp;
                        }
                    }
                    // *gbl->log << "degenerate case" << v(0) << ' ' << v(1) << std::endl;
                }
            }
            addtri(v(0),v(1),gbl->i2wk_lst3(0),sind,dir);
            /* ADD ANY DEGENERATE TRIANGLES */              
            for(i=1;i<ngood;++i)
                addtri(gbl->i2wk_lst3(i-1),v(1),gbl->i2wk_lst3(i),-1,-1);
        }
      
        bgn = end;
        end += nside -nsidebefore;
        for(i=nsidebefore;i<nside;++i)
            gbl->i2wk_lst1(bgn+i-nsidebefore) = -(i + 1);
    }
            
    return;

}

void tri_mesh::addtri(int v0,int v1, int v2, int sind, int dir) {
    int i,j,k,end,sind1,tind;
    int minv,maxv,order,sindprev,temp;
        
    /* ADD NEW TRIANGLE */
    td(ntri).vrtx(0) = v0;
    td(ntri).vrtx(1) = v1;
    td(ntri).vrtx(2) = v2;
                    
    vd(v0).tri = ntri;
    vd(v1).tri = ntri;
    vd(v2).tri = ntri;
    
    end = 3;
    if (sind > -1) {
        /* SIDE 2 INFO IS KNOWN ALREADY */
        td(ntri).side(2) = sind;
        td(ntri).sign(2) = 1 -2*dir;
        sd(sind).tri(dir) = ntri;
        tind = sd(sind).tri(1-dir);
        td(ntri).tri(2) = tind;
        if (tind > -1) {
            for(i=0;i<3;++i) {
                if (td(tind).side(i) == sind) {
                    td(tind).tri(i) = ntri;
                    break;
                }
            }
        }
        end = 2;
    }

    /* LOOP THROUGH SIDES */
    for(k=0;k<end;++k) {
        if (v2 > v1) {
            minv = v1;
            maxv = v2;
            order = 1;
        }
        else {
            minv = v2;
            maxv = v1;
            order = 0;
        }
        
        sind1 = vd(minv).info;
        while (sind1 >= 0) {
            if (maxv == sd(sind1).vrtx(order)) {
                /* SIDE IN SAME DIRECTION */
                if (sd(sind1).tri(0) >= 0) {
                    *gbl->log << "1:side already matched?" << sind1 << ' ' << v1 << ' ' << v2 << std::endl;
                    output("error",tecplot);
                    output("error",grid);
                    exit(1);
                }
                sd(sind1).tri(0) = ntri;
                td(ntri).side(k) = sind1;
                td(ntri).sign(k) = 1;
                tind = sd(sind1).tri(1);
                td(ntri).tri(k) = tind;
                if (tind > -1) {
                    for(j=0;j<3;++j) {
                        if (td(tind).side(j) == sind1) {
                            td(tind).tri(j) = ntri;
                            break;
                        }
                    }
                }
                goto NEXTTRISIDE;
            }
            else if(maxv == sd(sind1).vrtx(1-order)) {
                /* SIDE IN OPPOSITE DIRECTION */
                if (sd(sind1).tri(1) >= 0) {
                    *gbl->log << "2:side already matched?" << sind1 << ' ' << v1 << ' ' << v2 << std::endl;
                    output("error",tecplot);
                    output("error",grid);
                    exit(1);
                }
                sd(sind1).tri(1) = ntri;
                td(ntri).side(k) = sind1;
                td(ntri).sign(k) = -1;
                tind = sd(sind1).tri(0);
                td(ntri).tri(k) = tind;
                if (tind > -1) {
                    for(j=0;j<3;++j) {
                        if (td(tind).side(j) == sind1) {
                            td(tind).tri(j) = ntri;
                            break;
                        }
                    }
                }
                goto NEXTTRISIDE;      
            }
            sindprev = sind1;
            sind1 = sd(sind1).info;
        }
        /* NEW SIDE */
        sd(nside).vrtx(0) = v1;
        sd(nside).vrtx(1) = v2;
        sd(nside).tri(0) = ntri;
        sd(nside).tri(1) = -1;
        td(ntri).side(k) = nside;
        td(ntri).sign(k) = 1;
        sd(nside).info = -1;
        if (vd(minv).info < 0)
            vd(minv).info = nside;
        else 
            sd(sindprev).info = nside;
        ++nside;
        assert(nside < maxvst -1);
            
NEXTTRISIDE:
        temp = v0;
        v0 = v1;
        v1 = v2;
        v2 = temp;
    }
    ++ntri;
    
    assert(ntri < maxvst -1);
        
    return;
}

