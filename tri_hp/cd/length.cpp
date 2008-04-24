/*
 *  lengthcd.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on Thu May 29 2003.
 *  Copyright (c) 2003 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include "../hp_boundary.h"
#include <math.h>
#include <utilities.h>

void tri_hp_cd::length() {
    int i,j,k,v0,v1,v2,sind,tind,count;
    TinyVector<FLT,2> dx0,dx1,dx2,ep,dedpsi;
    FLT sum,ratio;
    FLT length0,length1,length2,lengthept;
    FLT ang1,curved1,ang2,curved2;
    
    gbl->fltwk(Range(0,npnt-1)) = 0.0;

    switch(basis::tri(log2p).p) {
        case(1): {
            for(i=0;i<nseg;++i) {
                v0 = seg(i).pnt(0);
                v1 = seg(i).pnt(1);
                sum = distance2(v0,v1)*fabs(ug.v(v0,0) -ug.v(v1,0));
                gbl->fltwk(v0) += sum;
                gbl->fltwk(v1) += sum;
            }
            break;
        }
                  
        default: {
            for(i=0;i<nseg;++i) {
                v0 = seg(i).pnt(0);
                v1 = seg(i).pnt(1);
                sum = distance2(v0,v1)*fabs(ug.s(i,sm0-1,0));
                gbl->fltwk(v0) += sum;
                gbl->fltwk(v1) += sum;
                /* UNCOMMENT FOR SHOCK DETECTION 
                sum = abs(ug.s(i,sm0-1,0)/(abs(ug.s(i,0,0))+1.0e-5));
                gbl->fltwk(v0) = MAX(sum,gbl->fltwk(v0));
                gbl->fltwk(v1) = MAX(sum,gbl->fltwk(v1)); */
            }

          /* BOUNDARY CURVATURE */
            for(i=0;i<nebd;++i) {
                if (!(hp_ebdry(i)->is_curved())) continue;
                
                for(j=0;j<ebdry(i)->nseg;++j) {
                    sind = ebdry(i)->seg(j);
                    v1 = seg(sind).pnt(0);
                    v2 = seg(sind).pnt(1);
                    
                    crdtocht1d(sind);
                                                
                    /* FIND ANGLE BETWEEN LINEAR SIDES */
                    tind = seg(sind).tri(0);
                    for(k=0;k<3;++k)
                        if (tri(tind).seg(k) == sind) break;
                    
                    v0 = tri(tind).pnt(k);
                    
                    dx0(0) = pnts(v2)(0)-pnts(v1)(0);
                    dx0(1) = pnts(v2)(1)-pnts(v1)(1);
                    length0 = dx0(0)*dx0(0) +dx0(1)*dx0(1);
                    
                    dx1(0) = pnts(v0)(0)-pnts(v2)(0);
                    dx1(1) = pnts(v0)(1)-pnts(v2)(1);
                    length1 = dx1(0)*dx1(0) +dx1(1)*dx1(1);
                    
                    dx2(0) = pnts(v1)(0)-pnts(v0)(0);
                    dx2(1) = pnts(v1)(1)-pnts(v0)(1);
                    length2 = dx2(0)*dx2(0) +dx2(1)*dx2(1);
                    
                    basis::tri(log2p).ptprobe1d(2,&ep(0),&dedpsi(0),-1.0,&cht(0,0),MXTM);
                    lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);
                    
                    ang1 = acos(-(dx0(0)*dx2(0) +dx0(1)*dx2(1))/sqrt(length0*length2));
                    curved1 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));
                    
                    basis::tri(log2p).ptprobe1d(2,&ep(0),&dedpsi(0),1.0,&cht(0,0),MXTM);
                    lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);
                    
                    ang2 = acos(-(dx0(0)*dx1(0) +dx0(1)*dx1(1))/sqrt(length0*length1));
                    curved2 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));                            

                    sum = gbl->curvature_sensitivity*(curved1/ang1 +curved2/ang2);
                    gbl->fltwk(v0) += sum*gbl->error_target*pnt(v0).nnbor;
                    gbl->fltwk(v1) += sum*gbl->error_target*pnt(v1).nnbor;
                }
            }
            break;
        }
    }
    
    // output_error(); FOR SHOCK DETECTION

    for(i=0;i<npnt;++i) {
        gbl->fltwk(i) = pow(gbl->fltwk(i)/(pnt(i).nnbor*gbl->error_target),1./(basis::tri(log2p).p+1+ND));
        lngth(i) /= gbl->fltwk(i);  
        lngth(i) = MAX(lngth(i),gbl->minlngth);
        lngth(i) = MIN(lngth(i),gbl->maxlngth);
    }

    /* AVOID HIGH ASPECT RATIOS */
    int nsweep = 0;
    do {
        count = 0;
        for(i=0;i<nseg;++i) {
            v0 = seg(i).pnt(0);
            v1 = seg(i).pnt(1);
            ratio = lngth(v1)/lngth(v0);
            
            if (ratio > 3.0) {
                lngth(v1) = 2.5*lngth(v0);
                ++count;
            }
            else if (ratio < 0.333) {
                lngth(v0) = 2.5*lngth(v1);
                ++count;
            }
        }
        ++nsweep;
        *gbl->log << "#aspect ratio fixes " << nsweep << ' ' << count << std::endl;
    } while(count > 0 && nsweep < 5);

    return;
}
