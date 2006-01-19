/*
 *  length.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_ins.h"
#include "hp_boundary.h"
#include <math.h>
#include<utilities.h>

/* THIS FUNCTION WILL SET THE VLNGTH VALUES BASED ON THE TRUNCATION ERROR */

#ifdef DROP
extern FLT dydt;
#endif

block::ctrl tri_hp_ins::length(int excpt) {
   int i,j,k,v0,v1,v2,indx,sind,tind,count;
   TinyVector<FLT,2> dx0,dx1,dx2,ep,dedpsi;
   FLT q,p,duv,um,vm,u,v;
   FLT sum,ruv,ratio;
   FLT length0,length1,length2,lengthept;
   FLT ang1,curved1,ang2,curved2;
   FLT norm;
   
   switch(excpt) {
      case(0): {
         for(tind=0;tind<ntri;++tind) {
            q = 0.0;
            p = 0.0;
            duv = 0.0;
            um = ug.v(td(tind).vrtx(2),0);
            vm = ug.v(td(tind).vrtx(2),1);
#ifdef DROP
            vm -= dydt;
#endif
            for(j=0;j<3;++j) {
               v0 = td(tind).vrtx(j);
               u = ug.v(v0,0);
               v = ug.v(v0,1);
#ifdef DROP
               v -= dydt;
#endif
               q += pow(u,2) +pow(v,2);
               p += fabs(ug.v(v0,2));
               duv += fabs(u-um)+fabs(v-vm);
               um = u;
               vm = v;
            }
            ins_gbl->eanda(0) += 1./3.*( (0.5*ins_gbl->rho*q +p)*area(tind) +duv*ins_gbl->mu*sqrt(area(tind)) );
            ins_gbl->eanda(1) += area(tind);
         }
         sim::blks.allreduce1(ins_gbl->eanda.data(),ins_gbl->eanda_recv.data());
         return(block::advance);
      }
      
      case(1): {
         sim::blks.allreduce2(2,blocks::flt_msg,blocks::sum);
         return(block::advance);
      }
      
      case(2): {
         norm = ins_gbl->eanda_recv(0)/ins_gbl->eanda_recv(1);
         fscr1(Range(0,nvrtx-1)) = 0.0;

         switch(basis::tri(log2p).p) {
            case(1): {
               for(i=0;i<nside;++i) {
                  v0 = sd(i).vrtx(0);
                  v1 = sd(i).vrtx(1);
                  u = fabs(ug.v(v0,0) +ug.v(v1,0));
                  v = fabs(ug.v(v0,1) +ug.v(v1,1));
#ifdef DROP
                  v -= dydt;
#endif
                  ruv = ins_gbl->rho*0.5*(u + v) +ins_gbl->mu/distance(v0,v1);
                  sum = distance2(v0,v1)*(ruv*(fabs(ug.v(v0,0) -ug.v(v1,0)) +fabs(ug.v(v0,1) -ug.v(v1,1))) +fabs(ug.v(v0,2) -ug.v(v1,2)));
                  fscr1(v0) += sum;
                  fscr1(v1) += sum;
               }                     
               break;
            }
         
            default: {
               indx = basis::tri(log2p).sm-1;
               for(i=0;i<nside;++i) {
                  v0 = sd(i).vrtx(0);
                  v1 = sd(i).vrtx(1);
                  u = fabs(ug.v(v0,0) +ug.v(v1,0));
                  v = fabs(ug.v(v0,1) +ug.v(v1,1));
#ifdef DROP
                  v -= dydt;
#endif
                  ruv = ins_gbl->rho*0.5*(u + v) +ins_gbl->mu/distance(v0,v1);
                  sum = distance2(v0,v1)*(ruv*(fabs(ug.s(i,indx,0)) +fabs(ug.s(i,indx,1))) +fabs(ug.s(i,indx,2)));
                  fscr1(v0) += sum;
                  fscr1(v1) += sum;
               }

              /* BOUNDARY CURVATURE */
               for(i=0;i<nsbd;++i) {
                  if (!(hp_sbdry(i)->is_curved())) continue;
                  
                  for(j=0;j<sbdry(i)->nel;++j) {
                     sind = sbdry(i)->el(j);
                     v1 = sd(sind).vrtx(0);
                     v2 = sd(sind).vrtx(1);
                     
                     crdtocht1d(sind);
                                          
                     /* FIND ANGLE BETWEEN LINEAR SIDES */
                     tind = sd(sind).tri(0);
                     for(k=0;k<3;++k)
                        if (td(tind).side(k) == sind) break;
                     
                     v0 = td(tind).vrtx(k);
                     
                     dx0(0) = vrtx(v2)(0)-vrtx(v1)(0);
                     dx0(1) = vrtx(v2)(1)-vrtx(v1)(1);
                     length0 = dx0(0)*dx0(0) +dx0(1)*dx0(1);
                     
                     dx1(0) = vrtx(v0)(0)-vrtx(v2)(0);
                     dx1(1) = vrtx(v0)(1)-vrtx(v2)(1);
                     length1 = dx1(0)*dx1(0) +dx1(1)*dx1(1);
                     
                     dx2(0) = vrtx(v1)(0)-vrtx(v0)(0);
                     dx2(1) = vrtx(v1)(1)-vrtx(v0)(1);
                     length2 = dx2(0)*dx2(0) +dx2(1)*dx2(1);
                     
                     basis::tri(log2p).ptprobe1d(2,&ep(0),&dedpsi(0),-1.0,&cht(0,0),MXTM);
                     lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);
                     
                     ang1 = acos(-(dx0(0)*dx2(0) +dx0(1)*dx2(1))/sqrt(length0*length2));
                     curved1 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));
                     
                     basis::tri(log2p).ptprobe1d(2,&ep(0),&dedpsi(0),1.0,&cht(0,0),MXTM);
                     lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);
                     
                     ang2 = acos(-(dx0(0)*dx1(0) +dx0(1)*dx1(1))/sqrt(length0*length1));
                     curved2 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));                     

                     sum = bdrysensitivity*(curved1/ang1 +curved2/ang2);
                     fscr1(v0) += sum*trncerr*norm*vd(v0).nnbor;
                     fscr1(v1) += sum*trncerr*norm*vd(v1).nnbor;
                  }
               }
               break;
            }
         }

         for(i=0;i<nvrtx;++i) {
            fscr1(i) = pow(fscr1(i)/(norm*vd(i).nnbor*trncerr),1./(basis::tri(log2p).p+1+ND));
            vlngth(i) /= fscr1(i);      

#ifdef THREELAYER
#define TRES 0.025/THREELAYER
            if (vrtx(i)(1) > 0.525) {
               vlngth(i) = MIN(vlngth(i),TRES +(vrtx(i)(1)-0.525)*(9*TRES)/0.475);
            }
            else if (vrtx(i)(1) < 0.475) {
               vlngth(i) = MIN(vlngth(i),TRES +(0.475 -vrtx(i)(1))*(9*TRES)/0.475);
            }
            else {
               vlngth(i) = MIN(vlngth(i),TRES);
            }
#endif
#ifdef TWOLAYER
            vlngth(i) = MIN(vlngth(i),0.3333); 
#endif
         }
   
         /* AVOID HIGH ASPECT RATIOS */
         int nsweep = 0;
         do {
            count = 0;
            for(i=0;i<nside;++i) {
               v0 = sd(i).vrtx(0);
               v1 = sd(i).vrtx(1);
               ratio = vlngth(v1)/vlngth(v0);
               
               if (ratio > 3.0) {
                  vlngth(v1) = 2.5*vlngth(v0);
                  ++count;
               }
               else if (ratio < 0.333) {
                  vlngth(v0) = 2.5*vlngth(v1);
                  ++count;
               }
            }
            ++nsweep;
            printf("#aspect ratio fixes %d: %d\n",nsweep,count);
         } while(count > 0 && nsweep < 5);
      }
   }
   return(block::stop);
}
