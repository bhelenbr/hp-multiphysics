/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ins.h"
#include "hp_boundary.h"

#ifdef DROP
TinyVector<FLT,mesh::ND> tri_hp_ins::mesh_ref_vel = 0.0;
#endif

 void tri_hp_ins::init(input_map& input, gbl *gin) {
   bool coarse, adapt_storage;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   keyword = idprefix + ".nvariable";
   input[keyword] = "3";
   
   tri_hp::init(input,gin);
   
   /* Load pointer to block stuff */
   ins_gbl = gin;
  
   keyword = idprefix + ".adapt_storage";
   input.getwdefault(keyword,adapt_storage,false);
   if (adapt_storage) return;

   keyword = idprefix + ".coarse";
   input.getwdefault(keyword,coarse,false);
   
   keyword = idprefix + ".dissipation";
   input.getwdefault(keyword,adis,1.0);
   
   if (coarse) return;
   
   ins_gbl->tau.resize(maxvst);
   ins_gbl->delt.resize(maxvst);
  
   keyword = idprefix + ".rho";
   input.getwdefault(keyword,ins_gbl->rho,1.0);

   keyword = idprefix + ".mu";
   input.getwdefault(keyword,ins_gbl->mu,0.0);

   ins_gbl->nu = ins_gbl->mu/ins_gbl->rho;
      
   return;
}

/* OVERRIDE VIRTUAL FUNCTION FOR INCOMPRESSIBLE FLOW */
void tri_hp_ins::calculate_unsteady_sources(bool coarse) {
   int i,j,n,tind;
      
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
               cjcb(i,j) = -sim::bd[0]*ins_gbl->rho*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
               for(n=0;n<ND;++n)
                  dugdt(log2p,tind,n)(i,j) = u(n)(i,j)*cjcb(i,j);
               dugdt(log2p,tind,ND)(i,j) = cjcb(i,j);

               for(n=0;n<ND;++n)
                  dxdt(log2p,tind,n)(i,j) = crd(n)(i,j);
            }            
         }
      }
   }
   log2p=log2pmax;
   
   return;
}