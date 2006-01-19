#include "tri_hp_cd.h"
#include <math.h>
#include <utilities.h>
#include "hp_boundary.h"


block::ctrl tri_hp_cd::setup_preconditioner(int excpt) {
   int tind,i,j,side,v0;
   FLT jcb,h,hmax,q,qmax,lam1;
   TinyVector<int,3> v;
   
   switch(excpt) {
      case(0): {
         /***************************************/
         /** DETERMINE FLOW PSEUDO-TIME STEP ****/
         /***************************************/
         hp_gbl->vprcn(Range(0,nvrtx-1),Range::all()) = 0.0;
         if (basis::tri(log2p).sm > 0) {
            hp_gbl->sprcn(Range(0,nside-1),Range::all()) = 0.0;
         }
         
#ifdef TIMEACCURATE
         FLT dtstari = 0.0;
#endif

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
               q = pow(cd_gbl->ax -(sim::bd[0]*(vrtx(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
                  +pow(cd_gbl->ay -(sim::bd[0]*(vrtx(v0)(1) -vrtxbd(1)(v0)(1))),2.0);
               qmax = MAX(qmax,q);
            }
            q = sqrt(qmax);
            
            lam1  = (q +1.5*cd_gbl->nu/h +h*sim::bd[0]);
                        
            /* SET UP DISSIPATIVE COEFFICIENTS */
            cd_gbl->tau(tind)  = adis*h/(jcb*lam1);
            
            jcb *= lam1/h;

      
            /* SET UP DIAGONAL PRECONDITIONER */
#ifdef TIMEACCURATE
            dtstari = MAX(lam1/h,dtstari);
         }
         printf("#iterative to physical time step ratio: %f\n",sim::bd[0]/dtstari);
            
         for(tind=0;tind<ntri;++tind) {
            v = td(tind).vrtx;
            jcb = 0.25*area(tind)*dtstari;
#endif
#ifdef AXISYMMETRIC
            jcb *= (vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.;
#endif
            hp_gbl->tprcn(tind,0) = jcb;   
            for(i=0;i<3;++i) {
               hp_gbl->vprcn(v(i),Range::all())  += hp_gbl->tprcn(tind,Range::all());
               if (basis::tri(log2p).sm > 0) {
                  side = td(tind).side(i);
                  hp_gbl->sprcn(side,Range::all()) += hp_gbl->tprcn(tind,Range::all());
               }
            }
         }
            
         /* SET UP TSTEP FOR ACTIVE BOUNDARIES */   
         for(i=0;i<nsbd;++i)
            hp_sbdry(i)->setup_preconditioner(0);
   
         mp_phase = -1;
         return(block::advance);
      }
   }
   
   return(tri_hp::setup_preconditioner(excpt));
}