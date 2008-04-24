#include <math.h>

#include "tri_hp_swe.h"
#include "../hp_boundary.h"

void tri_hp_swe::setup_preconditioner() {
    int tind,i,j,side,v0;
    FLT jcb,hmax,q,qmax,umax,vmax,c,c2,pre,rtpre,fmax,cflnow,alpha,alpha2,sigma;
    FLT dx, dxmax, lambdamax;
    TinyVector<int,3> v;
    
    cflnow = 0.0;

    /***************************************/
    /** DETERMINE FLOW PSEUDO-TIME STEP ****/
    /***************************************/
    gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
    if (basis::tri(log2p).sm > 0) {
        gbl->sprcn(Range(0,nseg-1),Range::all()) = 0.0;
    }

    for(tind = 0; tind < ntri; ++tind) {
        jcb = 0.25*area(tind);  // area is 2 x triangle area
        v = tri(tind).pnt;
        dxmax = 0.0;
        for(j=0;j<3;++j) {
            dx = pow(pnts(v(j))(0) -pnts(v((j+1)%3))(0),2.0) + 
            pow(pnts(v(j))(1) -pnts(v((j+1)%3))(1),2.0);
            dxmax = (dx > dxmax ? dx : dxmax);
        }
        dxmax = sqrt(dxmax);
        dx = 4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*dxmax);
        dxmax /= 0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1);
                
        qmax = 0.0;
        hmax = 0.0;
        umax = 0.0;
        vmax = 0.0;
        fmax = 0.0;
        for(j=0;j<3;++j) {
            v0 = v(j);
            q = pow(ug.v(v0,0)/ug.v(v0,NV-1) -gbl->bd[0]*(pnts(v0)(0) -vrtxbd(1)(v0)(0)),2.0) 
                +pow(ug.v(v0,1)/ug.v(v0,NV-1) -gbl->bd[0]*(pnts(v0)(1) -vrtxbd(1)(v0)(1)),2.0);  
            qmax = MAX(qmax,q);
            hmax = MAX(hmax,ug.v(v0,NV-1));
            umax = MAX(umax,ug.v(v0,0)/ug.v(v0,NV-1));
            vmax = MAX(vmax,ug.v(v0,1)/ug.v(v0,NV-1));
            fmax = MAX(fmax,fabs(gbl->f0 +gbl->beta*pnts(v0)(1)));
        }
        if (!(jcb > 0.0) || !(hmax > 0.0)) {  // THIS CATCHES NAN'S TOO
            *gbl->log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
            *gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1) << std::endl;
            tri_mesh::output("negative",grid);
            exit(1);
        }
        
        q = sqrt(qmax);
        c2 = gbl->g*hmax;
        c = sqrt(c2);
        sigma = MAX((qmax -c2)/qmax,0);
        alpha = gbl->cd*dxmax/(2*hmax);
        alpha2 = alpha*alpha;
//            pre = gbl->ptest*(pow(dxmax*(gbl->bd[0]+fmax),2.0) +(3.+alpha2)*qmax)/(dxmax*dxmax*(gbl->bd[0]*gbl->bd[0] +sigma*fmax*fmax) +c2 +(3.+sigma*alpha2)*qmax);
        pre = gbl->ptest*(pow(dxmax*(gbl->bd[0]*.5+fmax),2.0) +(1.+alpha2)*qmax)/(dxmax*dxmax*(gbl->bd[0]*gbl->bd[0]*.25 +sigma*fmax*fmax) +c2 +(1.+sigma*alpha2)*qmax);

        rtpre = sqrt(pre);
        lambdamax = ((1. +(2.+sqrt(7)-sqrt(3))/2.*MAX(1 -rtpre,0))*q +rtpre*c +2*alpha*q)/dx +gbl->bd[0] +fmax;
        
        /* SET UP DISSIPATIVE COEFFICIENTS */
        gbl->tau(tind,0) = adis/(lambdamax*jcb);
        gbl->tau(tind,1) = gbl->tau(tind,0);
        gbl->tau(tind,NV-1) = pre*gbl->tau(tind,0);
             
        /* SET UP DIAGONAL PRECONDITIONER */
        jcb *=  lambdamax;
        
        cflnow = MAX((q/dx +fmax +c/dx)/gbl->bd[0],cflnow);
        
        gbl->tprcn(tind,0) = jcb;    
        gbl->tprcn(tind,1) = jcb;      
        gbl->tprcn(tind,2) =  jcb/pre;
        for(i=0;i<3;++i) {
            gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
            if (basis::tri(log2p).sm > 0) {
                side = tri(tind).seg(i);
                gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
            }
        }
    }
    // if (!coarse && log2p == log2pmax) *gbl->log << "#cfl is " << cflnow << std::endl;

    tri_hp::setup_preconditioner();  

    return;
}
