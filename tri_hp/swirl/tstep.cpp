#include "tri_hp_swirl.h"
#include "hp_boundary.h"
#include <math.h>

block::ctrl tri_hp_swirl::setup_preconditioner(block::ctrl ctrl_message) {
   int tind,i,j,side,v0;
   FLT jcb,h,hmax,q,qmax,lam1,gam;
   TinyVector<int,3> v;

   if (ctrl_message == block::begin) excpt = 0;
   ctrl_message = tri_hp::setup_preconditioner(ctrl_message);
   
   if (ctrl_message == block::advance1) ++excpt;
   
   if (excpt == 3) {
      ++excpt;

      /***************************************/
      /** DETERMINE FLOW PSEUDO-TIME STEP ****/
      /***************************************/
      swirl_gbl->vprcn(Range(0,nvrtx-1),Range::all()) = 0.0;
      if (basis::tri(log2p).sm > 0) {
         hp_gbl->sprcn(Range(0,nside-1),Range::all()) = 0.0;
      }

      for(tind = 0; tind < ntri; ++tind) {
         jcb = 0.25*area(tind);  // area is 2 x triangle area
         v = td(tind).vrtx;
         hmax = 0.0;
         for(j=0;j<3;++j) {
            h = pow(vrtx(v(j))(0) -vrtx(v((j+1)%3))(0),2.0) + 
            pow(vrtx(v(j))(1) -vrtx(v((j+1)%3))(1),2.0);
            hmax = (h > hmax ? h : hmax);
         }
         hmax = sqrt(hmax);
         
         if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
            *sim::log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
            *sim::log << "approximate location: " << vrtx(v(0))(0) << ' ' << vrtx(v(0))(1) << std::endl;
            mesh::output("negative",grid);
            exit(1);
         }
         h = 4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
         hmax = hmax/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1));
      
         qmax = 0.0;
         for(j=0;j<3;++j) {
            v0 = v(j);
            q = pow(ug.v(v0,0)-0.5*(sim::bd[0]*(vrtx(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
               +pow(ug.v(v0,1)-0.5*(sim::bd[0]*(vrtx(v0)(1) -vrtxbd(1)(v0)(1))),2.0);
				q += pow(ug.v(v0,2),2.0);
            qmax = MAX(qmax,q);
         }
         gam = 3.0*qmax +(0.5*hmax*sim::bd[0] +2.*swirl_gbl->nu/hmax)*(0.5*hmax*sim::bd[0] +2.*swirl_gbl->nu/hmax);
         if (swirl_gbl->mu + sim::bd[0] == 0.0) gam = MAX(gam,0.01);
         q = sqrt(qmax);
         lam1 = q + sqrt(qmax +gam);
         
         /* SET UP DISSIPATIVE COEFFICIENTS */
         swirl_gbl->tau(tind) = adis*h/(jcb*sqrt(gam));
         swirl_gbl->delt(tind) = qmax*swirl_gbl->tau(tind);
         
         /* SET UP DIAGONAL PRECONDITIONER */
         // jcb *= 8.*hp_gbl->nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +sim::bd[0];
         jcb *= 2.*swirl_gbl->nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
         jcb *= (vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.;

         swirl_gbl->tprcn(tind,0) = swirl_gbl->rho*jcb;   
         swirl_gbl->tprcn(tind,1) = swirl_gbl->rho*jcb;  
         swirl_gbl->tprcn(tind,2) = swirl_gbl->rho*jcb;     
         swirl_gbl->tprcn(tind,3) =  jcb/gam;
         for(i=0;i<3;++i) {
            hp_gbl->vprcn(v(i),Range::all())  += hp_gbl->tprcn(tind,Range::all());
            if (basis::tri(log2p).sm > 0) {
               side = td(tind).side(i);
               hp_gbl->sprcn(side,Range::all()) += hp_gbl->tprcn(tind,Range::all());
            }
         }
      }
   }

   return(ctrl_message);
}