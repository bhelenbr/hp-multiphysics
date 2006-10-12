/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_lvlset.h"
#include "hp_boundary.h"

void tri_hp_lvlset::init(input_map& input, gbl *gin) {
   bool coarse, adapt_storage;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   keyword = idprefix + ".nvariable";
   input[keyword] = "4";
   
   tri_hp_ins::init(input,gin);
   
   /* Load pointer to block stuff */
   gbl_ptr = gin;
   
   keyword = idprefix + ".adapt_storage";
   input.getwdefault(keyword,adapt_storage,false);
   if (adapt_storage) return;
   
   keyword = idprefix + ".coarse";
   input.getwdefault(keyword,coarse,false);
   if (coarse) return;

   gbl_ptr->tau_lvlset.resize(maxvst);
  
   keyword = idprefix + ".rho2";
   input.getwdefault(keyword,gbl_ptr->rho,1.0);

   keyword = idprefix + ".mu2";
   input.getwdefault(keyword,gbl_ptr->mu,0.0);
   
   keyword = idprefix + ".sigma";
   input.getwdefault(keyword,gbl_ptr->sigma,0.0);

   keyword = idprefix + ".width";
   input.getwdefault(keyword,gbl_ptr->width,0.02);
   
   return;
}

void tri_hp_lvlset::calculate_unsteady_sources(bool coarse) {
   int i,j,n,tind;
   FLT rho;
   
   for (log2p=0;log2p<=log2pmax;++log2p) {
      for(tind=0;tind<ntri;++tind) {
         if (td(tind).info > -1) {
            crdtocht(tind,1);
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
         }
         else {
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj(vrtxbd(1)(td(tind).vrtx(0))(n),vrtxbd(1)(td(tind).vrtx(1))(n),vrtxbd(1)(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);

            for(i=0;i<basis::tri(log2p).gpx;++i) {
               for(j=0;j<basis::tri(log2p).gpn;++j) {
                  for(n=0;n<ND;++n) {
                     dcrd(n,0)(i,j) = 0.5*(vrtxbd(1)(td(tind).vrtx(1))(n) -vrtxbd(1)(td(tind).vrtx(0))(n));
                     dcrd(n,1)(i,j) = 0.5*(vrtxbd(1)(td(tind).vrtx(2))(n) -vrtxbd(1)(td(tind).vrtx(0))(n));
                  }
               }
            }
         }
         
         ugtouht(tind,1);
         for(n=0;n<NV;++n)
            basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);

         for(i=0;i<basis::tri(log2p).gpx;++i) {
            for(j=0;j<basis::tri(log2p).gpn;++j) {   
               if (u(2)(i,j) < -gbl_ptr->width) rho = gbl_ptr->rho;
               else if (u(2)(i,j) > gbl_ptr->width) rho = gbl_ptr->rho2;
               else rho = 0.5*(gbl_ptr->rho +gbl_ptr->rho2 +(gbl_ptr->rho2 -gbl_ptr->rho)*sin(M_PI*u(2)(i,j)/gbl_ptr->width));
               cjcb(i,j) = -sim::bd[0]*rho*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
               for(n=0;n<NV-1;++n)
                  dugdt(log2p,tind,n)(i,j) = u(n)(i,j)*cjcb(i,j);
               dugdt(log2p,tind,NV-1)(i,j) = cjcb(i,j);

               for(n=0;n<ND;++n)
                  dxdt(log2p,tind,n)(i,j) = crd(n)(i,j);
            }            
         }
      }
   }
   log2p=log2pmax;
   
   return;
}