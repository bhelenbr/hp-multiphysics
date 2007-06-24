#include <math.h>

#include "tri_hp_swe.h"
#include "../hp_boundary.h"

block::ctrl tri_hp_swe::setup_preconditioner(block::ctrl ctrl_message) {
   int tind,i,j,side,v0;
   FLT jcb,hmax,q,qmax,umax,vmax,c,c2,pre,rtpre,fmax,cflnow,alpha,alpha2,sigma;
   FLT dx, dxmax, lambdamax;
   TinyVector<int,3> v;
   
   if (ctrl_message == block::begin) excpt = 0;
   if (ctrl_message == block::advance1) ++excpt;   
   
   cflnow = 0.0;
   
   if (excpt == 3) {
      ++excpt;
   	        
      /***************************************/
      /** DETERMINE FLOW PSEUDO-TIME STEP ****/
      /***************************************/
      gbl_ptr->vprcn(Range(0,nvrtx-1),Range::all()) = 0.0;
      if (basis::tri(log2p).sm > 0) {
         gbl_ptr->sprcn(Range(0,nside-1),Range::all()) = 0.0;
      }

      for(tind = 0; tind < ntri; ++tind) {
         jcb = 0.25*area(tind);  // area is 2 x triangle area
         v = td(tind).vrtx;
         dxmax = 0.0;
         for(j=0;j<3;++j) {
            dx = pow(vrtx(v(j))(0) -vrtx(v((j+1)%3))(0),2.0) + 
            pow(vrtx(v(j))(1) -vrtx(v((j+1)%3))(1),2.0);
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
            q = pow(ug.v(v0,0)/ug.v(v0,NV-1) -sim::bd[0]*(vrtx(v0)(0) -vrtxbd(1)(v0)(0)),2.0) 
               +pow(ug.v(v0,1)/ug.v(v0,NV-1) -sim::bd[0]*(vrtx(v0)(1) -vrtxbd(1)(v0)(1)),2.0);  
            qmax = MAX(qmax,q);
            hmax = MAX(hmax,ug.v(v0,NV-1));
            umax = MAX(umax,ug.v(v0,0)/ug.v(v0,NV-1));
            vmax = MAX(vmax,ug.v(v0,1)/ug.v(v0,NV-1));
            fmax = MAX(fmax,fabs(gbl_ptr->f0 +gbl_ptr->beta*vrtx(v0)(1)));
         }
         if (!(jcb > 0.0) || !(hmax > 0.0)) {  // THIS CATCHES NAN'S TOO
            *sim::log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
            *sim::log << "approximate location: " << vrtx(v(0))(0) << ' ' << vrtx(v(0))(1) << std::endl;
            mesh::output("negative",grid);
            exit(1);
         }
         
         q = sqrt(qmax);
         c2 = sim::g*hmax;
         c = sqrt(c2);
         sigma = MAX((qmax -c2)/qmax,0);
         alpha = gbl_ptr->cd*dxmax/(2*hmax);
         alpha2 = alpha*alpha;
//         pre = gbl_ptr->ptest*(pow(dxmax*(sim::bd[0]+fmax),2.0) +(3.+alpha2)*qmax)/(dxmax*dxmax*(sim::bd[0]*sim::bd[0] +sigma*fmax*fmax) +c2 +(3.+sigma*alpha2)*qmax);
         pre = gbl_ptr->ptest*(pow(dxmax*(sim::bd[0]*.5+fmax),2.0) +(1.+alpha2)*qmax)/(dxmax*dxmax*(sim::bd[0]*sim::bd[0]*.25 +sigma*fmax*fmax) +c2 +(1.+sigma*alpha2)*qmax);

         rtpre = sqrt(pre);
         lambdamax = ((1. +(2.+sqrt(7)-sqrt(3))/2.*MAX(1 -rtpre,0))*q +rtpre*c +2*alpha*q)/dx +sim::bd[0] +fmax;
         
         /* SET UP DISSIPATIVE COEFFICIENTS */
         gbl_ptr->tau(tind,0) = adis/(lambdamax*jcb);
         gbl_ptr->tau(tind,1) = gbl_ptr->tau(tind,0);
         gbl_ptr->tau(tind,NV-1) = pre*gbl_ptr->tau(tind,0);
             
         /* SET UP DIAGONAL PRECONDITIONER */
         jcb *=  lambdamax;
         
         cflnow = MAX((q/dx +fmax +c/dx)/sim::bd[0],cflnow);
         
         gbl_ptr->tprcn(tind,0) = jcb;   
         gbl_ptr->tprcn(tind,1) = jcb;     
         gbl_ptr->tprcn(tind,2) =  jcb/pre;
         for(i=0;i<3;++i) {
            gbl_ptr->vprcn(v(i),Range::all())  += gbl_ptr->tprcn(tind,Range::all());
            if (basis::tri(log2p).sm > 0) {
               side = td(tind).side(i);
               gbl_ptr->sprcn(side,Range::all()) += gbl_ptr->tprcn(tind,Range::all());
            }
         }
      }
      // if (!coarse && log2p == log2pmax) *sim::log << "#cfl is " << cflnow << std::endl;
   }
   else {
      ctrl_message = tri_hp::setup_preconditioner(ctrl_message);  
   }

   return(ctrl_message);
}
