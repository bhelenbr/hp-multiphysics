#include "tri_hp_ps.h"
#include "hp_boundary.h"
#include <math.h>

block::ctrl tri_hp_ps::setup_preconditioner(block::ctrl ctrl_message) {
   int tind,i,j,side;
   FLT jcb,h,hmax,lam1,gam,gami;
   TinyVector<int,3> v;

   if (ctrl_message == block::begin) excpt = 0;
   if (ctrl_message == block::advance1) ++excpt;

   if (excpt == 3) {
      ++excpt;

      /***************************************/
      /** DETERMINE PSEUDO-TIME STEP ****/
      /***************************************/
      gbl_ptr->vprcn(Range(0,nvrtx-1),Range::all()) = 0.0;
      if (basis::tri(log2p).sm > 0) {
         gbl_ptr->sprcn(Range(0,nside-1),Range::all()) = 0.0;
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

         gami = pow(hmax/(2.*gbl_ptr->mu),2.0) +gbl_ptr->lami*hmax*hmax/gbl_ptr->mu;
         gam = 1./gami; 
         lam1 = sqrt(gam);
         
         /* SET UP DISSIPATIVE COEFFICIENTS */
         gbl_ptr->tau(tind)  = adis*h/(jcb*lam1);
         jcb *= (8.*gbl_ptr->mu*(1./(hmax*hmax) +1./(h*h))) ;

         gbl_ptr->tprcn(tind,0) = jcb;   
         gbl_ptr->tprcn(tind,1) = jcb;     
         gbl_ptr->tprcn(tind,2) =  jcb/gam;
         for(i=0;i<3;++i) {
            gbl_ptr->vprcn(v(i),Range::all())  += gbl_ptr->tprcn(tind,Range::all());
            if (basis::tri(log2p).sm > 0) {
               side = td(tind).side(i);
               gbl_ptr->sprcn(side,Range::all()) += gbl_ptr->tprcn(tind,Range::all());
            }
         }
      }
   }
   else {
      ctrl_message = tri_hp::setup_preconditioner(ctrl_message);  
   }

   return(ctrl_message);
}
