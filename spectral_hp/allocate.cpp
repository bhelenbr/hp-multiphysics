/*
 *  allocate.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 11 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include <utilities.h>
#include <input_map.h>
#include <boundary.h>
#include "hp_boundary.h"

/* STATIC WORK ARRAYS */
Array<TinyMatrix<FLT,MXGP,MXGP>,1> tri_hp::u,tri_hp::res;
Array<TinyMatrix<FLT,MXGP,MXGP>,2> tri_hp::du;
TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_hp::ND> tri_hp::crd;
TinyMatrix<FLT,MXGP,MXGP> tri_hp::cjcb;
TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,tri_hp::ND,tri_hp::ND> tri_hp::dcrd;
Array<TinyVector<FLT,MXTM>,1> tri_hp::uht,tri_hp::lf;
TinyMatrix<FLT,tri_hp::ND,MXTM> tri_hp::cht, tri_hp::cf;
TinyVector<TinyMatrix<FLT,MXGP,MXGP>,tri_hp::ND> tri_hp::mvel;
Array<TinyMatrix<FLT,MXGP,MXGP>,2> tri_hp::bdwk;

 void tri_hp::init(std::map <std::string,std::string>& input, std::string prefix, tri_hp::gbl *hp_in) {
   int i,p,coarse,adapt_storage;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   keyword = prefix + ".coarse";
   data.str(input[keyword]);
   if (!(data >> coarse)) {
      coarse = 0;
   }
   data.clear();
   
   keyword = prefix + ".adapt_storage";
   data.str(input[keyword]);
   if (!(data >> adapt_storage)) {
      adapt_storage = 0;
   }
   data.clear();

   /* Initialize stuff for r_mesh */
   if (!adapt_storage && ((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable))) r_mesh::init(input,prefix,hp_in);
   
   keyword = prefix + ".nvariable";
   data.str(input[keyword]);
   if (!(data >> NV)) {
      NV = 1;
   }
   data.clear();
   
   keyword = prefix + ".log2p";
   data.str(input[keyword]);
   if (!(data >> log2p)) {
      data.clear();
      keyword = "log2p";
      data.str(input[keyword]);
      if (!(data >> log2p)) log2p = 0;
   }
   *sim::log << "#" << prefix << ".log2p: " << log2p << std::endl;
   log2pmax = log2p;
   data.clear();
   
   int npts;
   keyword = prefix + ".ngaussx";
   data.str(input[keyword]);
   if (!(data >> npts)) {
      data.clear();
      keyword = "ngaussx";
      data.str(input[keyword]);
      if (!(data >> npts)) npts = 0;
   }
   *sim::log << "#" << prefix << ".ngaussx: " << npts << std::endl;
   data.clear();

   /* Check that global basis has been allocated & allocate if necessary */
   if (basis::tri.extent(firstDim) < log2p +1) {
      basis::tri.resize(log2p+1);
      
      p = 1;
      for(i=0;i<log2p;++i) {
         basis::tri(i).initialize(p,p+1+npts);
         p = p<<1;
      }
   }
   
   if (coarse) {
      log2p = 0;
      p0 = 1;
      sm0 = 0;
      im0 = 0;
   }
   else {
      p0 = basis::tri(log2p).p;
      sm0 = basis::tri(log2p).sm;
      im0 = basis::tri(log2p).im;
   }
      
   /* Check that static work arrays are big enough */
   if (u.extent(firstDim) < NV) {
      u.resize(NV);
      res.resize(NV);
      du.resize(NV,ND);
      uht.resize(NV);
      lf.resize(NV);
      bdwk.resize(sim::nhist+1,NV);
   }
   
    
   /* Allocate solution vector storage */
   ug.v.resize(maxvst,NV);
   ug.s.resize(maxvst,sm0,NV);
   ug.i.resize(maxvst,im0,NV);
   
   /* Find bdryfile name in input map */
   std::string bdryfile;
   std::map<std::string,std::string> bdrymap;
   keyword = prefix + ".bdryfile";
   data.str(input[keyword]);
   if (data >> bdryfile) {
      if (!(strncmp("${HOME}",bdryfile.c_str(), 7))) {
         filename = getenv("HOME");
         filename = filename + (bdryfile.c_str()+7);
      }
      else 
         filename = bdryfile;
   }
   else {
      keyword = prefix + ".mesh";
      data.str(input[keyword]);
      filename = filename +"_bdry.inpt";
   }
   input_map(bdrymap,filename.c_str());
   data.clear();
      
   hp_sbdry.resize(nsbd);
   for(i=0;i<nsbd;++i) {
      hp_sbdry(i) = getnewsideobject(i,&bdrymap);
      hp_sbdry(i)->init(input,prefix);
   }
         
   setbcinfo();
      
   /* Load pointer to block stuff */
   hp_gbl = hp_in;
   
   if (!coarse) {
      /* Allocate time history stuff */
      /* For ease of access have level 0 reference ug */
      ugbd(0).v.reference(ug.v);
      ugbd(0).s.reference(ug.s);
      ugbd(0).i.reference(ug.i);
      vrtxbd(0).reference(vrtx);
      
      if (!adapt_storage) {
         for(i=1;i<sim::nhist+1;++i) {
            ugbd(i).v.resize(maxvst,NV);
            ugbd(i).s.resize(maxvst,sm0,NV);
            ugbd(i).i.resize(maxvst,im0,NV);
            vrtxbd(i).resize(maxvst);
         }
      }
      else {
         for(i=1;i<sim::nadapt+1;++i) {
            ugbd(i).v.resize(maxvst,NV);
            ugbd(i).s.resize(maxvst,sm0,NV);
            ugbd(i).i.resize(maxvst,im0,NV);
            vrtxbd(i).resize(maxvst);
         }
         return;
      }

#ifdef PV3
      /** Variables to understand iterative convergence using pV3 */
      ugpv3.v.resize(maxvst,NV);
      ugpv3.s.resize(maxvst,sm0,NV);
      ugpv3.i.resize(maxvst,im0,NV);
      vrtxpv3.resize(maxvst);
#endif

      /* Multigrid Storage all except highest order (log2p+1)*/
      dres.resize(log2p);
      for(i=0;i<log2p;++i) {
         dres(i).v.resize(maxvst,NV);
         dres(i).s.resize(maxvst,basis::tri(i).sm,NV);
         dres(i).i.resize(maxvst,basis::tri(i).im,NV);
      }
            
      /* Allocate block stuff */
      hp_gbl->ug0.v.resize(maxvst,NV);
      hp_gbl->ug0.s.resize(maxvst,sm0,NV);
      hp_gbl->ug0.i.resize(maxvst,im0,NV);
      
      hp_gbl->res.v.resize(maxvst,NV);
      hp_gbl->res.s.resize(maxvst,sm0,NV);
      hp_gbl->res.i.resize(maxvst,im0,NV);
      
      hp_gbl->res_r.v.resize(maxvst,NV);
      hp_gbl->res_r.s.resize(maxvst,sm0,NV);
      hp_gbl->res_r.i.resize(maxvst,im0,NV);

      hp_gbl->res0.v.resize(maxvst,NV);
      hp_gbl->res0.s.resize(maxvst,basis::tri(log2p).sm,NV);
      hp_gbl->res0.i.resize(maxvst,basis::tri(log2p).im,NV); 
      
      
#ifndef MATRIX_PRECONDITIONER
      hp_gbl->vprcn.resize(maxvst,NV);
      hp_gbl->sprcn.resize(maxvst,NV);
      hp_gbl->tprcn.resize(maxvst,NV);
#else
      hp_gbl->vprcn.resize(maxvst,NV,NV);
      hp_gbl->sprcn.resize(maxvst,NV,NV);
      hp_gbl->tprcn.resize(maxvst,NV,NV);
#endif
      
      /* Use global scratch for adaptation work */
      /* resize if necessary */
      /* Total scratch size needed for adaptation */
#ifdef USE_SIMSCRATCH
      if (sim::scratch.size() < needed_scratch_size()) {
         sim::scratch.resize(size);
      }
#endif
      reload_scratch_pointers();
   }
   else {
      ugbd(0).v.reference(ug.v);
      ugbd(0).s.reference(ug.s);
      ugbd(0).i.reference(ug.i);
      vrtxbd(0).reference(vrtx);
      
      /* HERE I MAKE 1 REFERENCE 0 BECAUSE I DON'T NEED TO STORE SOURCE TERMS ON COARSE MESHES */
      /* I JUST USE BD(1) TO CALCULATE UNSTEADY SOURCES (DUGDT,DVRTDT) IN TADVANCE */
      ugbd(1).v.reference(ug.v);
      ugbd(1).s.reference(ug.s);
      ugbd(1).i.reference(ug.i);
      vrtxbd(1).reference(vrtx);
      
      dres.resize(1);
      dres(0).v.resize(maxvst,NV);
      vrtx_frst.resize(maxvst);
   }
   
   /* UNSTEADY SOURCE TERMS */
   dugdt.resize(log2p+1,maxvst,NV);
   dxdt.resize(log2p+1,maxvst,ND);
   
      
   return;
}

 tri_hp::~tri_hp() {
   for(int i=0;i<nsbd;++i)
      delete hp_sbdry(i);
}

 void tri_hp::setbcinfo() {
   int i,j,sind;
   
   /* SET UP VRTX BC INFORMATION FOR OUTPUT */
   for(i=0;i<nvrtx;++i)
      vd(i).info = -1;

   for(i=0;i<nvbd;++i)
      vd(vbdry(i)->v0).info = 0;

   /* SET UP SIDE BC INFORMATION FOR CURVED SIDES OUTPUT */
   for(i=0;i<nside;++i)
      sd(i).info = -1;
   
   for(i=0;i<ntri;++i)
      td(i).info = -1;
   
   for(i=0;i<nsbd;++i) {
      if (hp_sbdry(i)->is_curved()) {
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            sd(sind).info = 0;
            td(sd(sind).tri(0)).info = 0;
         }
      } 
   }
   
   return;
}

 FLT tri_hp::maxres(FLT *mxr) {
   int i,n;
   FLT emax = 0.0;
   
   for(n=0;n<NV;++n)
      mxr[n] = 0.0;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         mxr[n] = MAX(mxr[n],fabs(hp_gbl->res.v(i,n)));

   for(n=0;n<NV;++n)
      emax = MAX(emax,mxr[n]);
      
   return(emax);
}   


