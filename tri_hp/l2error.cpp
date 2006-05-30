/*
 *  l2error.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include <utilities.h>
#include <boundary.h>
#include "hp_boundary.h"

 void tri_hp::l2error(FLT (*func)(int, TinyVector<FLT,2> &x)) {
	int i,j,n,tind,loc[NV];
	FLT err,mxr[NV],l2r[NV];
   TinyVector<FLT,2> pt;
   
#ifdef TAYLOR
   extern FLT ppipi;
//   
//   ptprobe(0.5,0.5,l2r);
//   ppipi = l2r[2];
   
/* MATCH PRESSURE AT ONE POINT */
   ppipi = 0.0;
   ppipi = -(*func)(2,vrtx(0)(0),vrtx(0)(1))+ug.v(0,2);
#endif
   
	for(n=0;n<NV;++n) {
		mxr[n] = 0.0;
		l2r[n] = 0.0;
	}
	
	for(tind=0;tind<ntri;++tind) {
      
      if (td(tind).info > -1) {
         crdtocht(tind);
         for(n=0;n<ND;++n)
            basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
      }
      else {
         for(n=0;n<ND;++n)
            basis::tri(log2p).proj(vrtx(td(tind).vrtx(0))(n),vrtx(td(tind).vrtx(1))(n),vrtx(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);

         for(i=0;i<basis::tri(log2p).gpx;++i) {
            for(j=0;j<basis::tri(log2p).gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd(n,0)(i,j) = 0.5*(vrtx(td(tind).vrtx(1))(n) -vrtx(td(tind).vrtx(0))(n));
                  dcrd(n,1)(i,j) = 0.5*(vrtx(td(tind).vrtx(2))(n) -vrtx(td(tind).vrtx(0))(n));
               }
            }
         }
      }

      ugtouht(tind);
		for(n=0;n<NV;++n)
			basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
		
 		for (i=0;i<basis::tri(log2p).gpx;++i) {	
			for (j=0;j<basis::tri(log2p).gpn;++j) {
            cjcb(i,j) = (dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
            pt(0) = crd(0)(i,j);
            pt(1) = crd(1)(i,j);
            for(n=0;n<NV;++n) {
               err =  fabs(u(n)(i,j)-func(n,pt));
               if (err > mxr[n]) {
                  mxr[n] = err;
                  loc[n] = tind;
               }
               l2r[n] += err*err*basis::tri(log2p).wtx(i)*basis::tri(log2p).wtn(j)*cjcb(i,j);
            }
         }
      }	
	}
   
	for(n=0;n<NV;++n) {
		l2r[n] = sqrt(l2r[n]); 
		printf("#L_2: %.3e L_inf %.3e %4d ",l2r[n],mxr[n],loc[n]);
	}
	printf("\n");
		
	return;
}

 block::ctrl tri_hp::findmax(block::ctrl ctrl_message, int bnum, FLT (*fxy)(TinyVector<FLT,2> &x)) {
   FLT ddpsi1, ddpsi2, psil, psir;
   TinyVector<FLT,2> xp, dx, maxloc, minloc;
   FLT max,min;
   int v0, sind;
   
   if (ctrl_message == block::begin) excpt = 0;
   else ++excpt;
   
   switch (excpt) {
      case(0):
         /* CALCULATE SLOPE AT ENDPOINT & TRANSMIT TO NEXT SURFACE */
         sind = sbdry(bnum)->el(sbdry(bnum)->nel-1);
         crdtocht1d(sind);
         basis::tri(log2p).ptprobe1d(ND,&xp(0),&dx(0),1.0,&cht(0,0),MXTM);
         ddpsi2 = (*fxy)(dx);
         vbdry(sbdry(bnum)->vbdry(1))->vloadbuff(boundary::manifolds,&ddpsi2,0,1,1);
         vbdry(sbdry(bnum)->vbdry(0))->comm_prepare(boundary::manifolds,0,boundary::master_slave);
         return(block::advance);
      
      case(1): 
         vbdry(sbdry(bnum)->vbdry(1))->comm_transmit(boundary::manifolds,0,boundary::master_slave);
         return(block::advance);
      
      case(2):
         vbdry(sbdry(bnum)->vbdry(0))->comm_wait(boundary::manifolds,0,boundary::master_slave);
         if (vbdry(sbdry(bnum)->vbdry(0))->is_comm()) 
            ddpsi2 = vbdry(sbdry(bnum)->vbdry(0))->frcvbuf(0,0);
         else
            ddpsi2 = 0.0;
            
         max = -1.0e99;
         min = 1.0e99;
         for(int indx=0;indx<sbdry(bnum)->nel;++indx) {
            sind = sbdry(bnum)->el(indx);
            crdtocht1d(sind);
            basis::tri(log2p).ptprobe1d(ND, &xp(0), &dx(0), -1.0, &cht(0,0), MXTM);
            ddpsi1 = (*fxy)(dx);
            if (ddpsi1 * ddpsi2 <= 0.0) {
               v0 = sd(sbdry(bnum)->el(indx)).vrtx(0);
               if ((*fxy)(vrtx(v0)) > max) {
                  maxloc[0] = vrtx(v0)(0);
                  maxloc[1] = vrtx(v0)(1);
                  max = (*fxy)(vrtx(v0));
               }
               if ((*fxy)(vrtx(v0)) < min) {
                  minloc[0] = vrtx(v0)(0);
                  minloc[1] = vrtx(v0)(1);
                  min = (*fxy)(vrtx(v0));
               }
               *sim::log << "#LOCAL EXTREMA: " << vrtx(v0)(0) << ' ' << vrtx(v0)(1) << ' ' <<(*fxy)(vrtx(v0)) << std::endl;
            }
            basis::tri(log2p).ptprobe1d(ND, &xp(0), &dx(0), 1.0, &cht(0,0), MXTM);
            ddpsi2 = (*fxy)(dx);
            if (ddpsi1 *ddpsi2 <= 0.0) {
               /* INTERIOR MAXIMUM */
               psil = -1.0;
               psir = 1.0;
               while (psir-psil > 1.0e-10) {
                  basis::tri(log2p).ptprobe1d(ND, &xp(0), &dx(0), 0.5*(psil +psir), &cht(0,0), MXTM);
                  if ((*fxy)(dx)*ddpsi1 < 0.0) 
                     psir = 0.5*(psil+psir);
                  else
                     psil = 0.5*(psil+psir);
               }
               if ((*fxy)(xp) > max) {
                  maxloc[0] = xp[0];
                  maxloc[1] = xp[1];
                  max = (*fxy)(xp);
               }
               if ((*fxy)(xp) < min) {
                  minloc[0] = xp[0];
                  minloc[1] = xp[1];
                  min = (*fxy)(xp);
               }
               *sim::log << "#LOCAL EXTREMA: " << xp[0] << ' ' << xp[1] << ' ' << (*fxy)(xp) << std::endl;
            }  
      }
      *sim::log << "#MAX EXTREMA: " << maxloc[0] << ' ' << maxloc[1] << ' ' << max << std::endl;
      *sim::log << "#MIN EXTREMA: " << minloc[0] << ' ' << minloc[1] << ' ' << min << std::endl;
      return(block::stop);
   }
   
   *sim::log << "flow control error" << std::endl;
   return(block::stop);
}

 void tri_hp::findintercept(int bnum, FLT (*fxy)(TinyVector<FLT,2> &x)) {
   FLT psil, psir;
   TinyVector<FLT,2> xp, dx;
   int v0, sind;
   FLT vl, vr;
   
   sind = sbdry(bnum)->el(0);
   crdtocht1d(sind);
   v0 = sd(sind).vrtx(0);
   vl = (*fxy)(vrtx(v0));

   for(int indx=0;indx<sbdry(bnum)->nel;++indx) {
      sind = sbdry(bnum)->el(indx);
      crdtocht1d(sind);
      v0 = sd(sind).vrtx(1);
      vr = (*fxy)(vrtx(v0));

      if (vl*vr <= 0.0) {
         /* INTERIOR INTERCEPT */
         psil = -1.0;
         psir = 1.0;
         while (psir-psil > 1.0e-10) {
            basis::tri(log2p).ptprobe1d(ND,&xp(0),&dx(0),0.5*(psil+psir),&cht(0,0),MXTM);
            if ((*fxy)(xp)*vl < 0.0) 
               psir = 0.5*(psil+psir);
            else
               psil = 0.5*(psil+psir);
         }
         *sim::log << "#INTERSECTION: " << xp[0] << ' ' << xp[1] << ' ' << (*fxy)(xp) << std::endl;
      }
      vl = vr; 
   }
   
   return;
}
