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

const int nmovetypes = 5;
const char movetypes[nmovetypes][80] = {"fixed","uncoupled_rigid","coupled_rigid","uncoupled_deformable","coupled_deformable"};


 void tri_hp::init(input_map& inmap, tri_hp::gbl *hp_in) {
   int i,ival,p,adapt_storage;
   std::string keyword, line;
   std::istringstream data;
   std::string filename;
   
   keyword = idprefix + ".coarse";
   inmap.getwdefault(keyword,coarse,false);
   *sim::log << "#" << keyword << ": " << coarse << std::endl;
   
   keyword = idprefix + ".adapt_storage";
   inmap.getwdefault(keyword,adapt_storage,0);
      
   keyword = idprefix + ".mesh_movement";
   if (!inmap.get(keyword,line)) {
      keyword = "mesh_movement";
      inmap.getwdefault(keyword,line,std::string("fixed"));
   }
   
   for (i=0;i<nmovetypes;++i)
      if (line == movetypes[i]) break;
      
   if (i == nmovetypes) 
      *sim::log << "unrecognized mesh movement type" << std::endl;
   
   mmovement = static_cast<movementtype>(i);
   *sim::log << "#" << keyword << ": " << line << std::endl;

   /* Initialize stuff for r_mesh */
   if (!adapt_storage && ((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable))) r_mesh::init(inmap,hp_in);
   
   keyword = idprefix + ".nvariable";
   inmap.getwdefault(keyword,NV,1);
   *sim::log << "#" << keyword << ": " << NV << std::endl;

   keyword = idprefix + ".log2p";
   if (!inmap.get(keyword,log2p)) {
      inmap.getwdefault("log2p",log2p,0);
   }
   *sim::log << "#" << keyword << ": " << log2p << std::endl;
   
#ifdef AXISYMMETRIC
   int npts = 1;
#else
   int npts = 0;
#endif

   /* Check that global basis has been allocated & allocate if necessary */
   if (basis::tri.extent(firstDim) < log2p +1) {
      basis::tri.resize(log2p+1);
      
      p = 1;
      for(i=0;i<log2p+1;++i) {
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
   log2pmax = log2p;
   
   TinyVector<std::string,3> output_purposes;
   TinyVector<filetype,3> defaults;
   output_purposes(0) = "display_type";
   defaults(0) = tri_hp::tecplot;
   output_purposes(1) = "restart_type";
   defaults(1) = tri_hp::text;
   output_purposes(2) = "debug_type";
   defaults(2) = tri_hp::tecplot;

   for(int i=0 ;i<3;++i) {
      keyword = idprefix + "." + output_purposes(i);
      if (inmap.get(keyword,ival)) {
         output_type(i) = static_cast<filetype>(ival);
      }
      else {
         if (inmap.get(output_purposes(i),ival)) {
            output_type(i) = static_cast<filetype>(ival);
         }
         else {
            output_type(i) = defaults(i);
         }
      }
   }

         
   /* Check that static work arrays are big enough */
   if (u.extent(firstDim) < NV) {
      u.resize(NV);
      res.resize(NV);
      du.resize(NV,ND);
      uht.resize(NV);
      lf.resize(MAX(NV,ND));
      bdwk.resize(sim::nhist+1,MAX(NV,ND));
   }
   
    
   /* Allocate solution vector storage */
   ug.v.resize(maxvst,NV);
   ug.s.resize(maxvst,sm0,NV);
   ug.i.resize(maxvst,im0,NV);
   
   /* Find bdryfile name in inmap map */
   std::string bdryfile;
   input_map bdrymap;
   keyword = idprefix + ".bdryfile";
   if (!inmap.get(keyword,bdryfile)) {
      keyword = idprefix + ".mesh";
      inmap.get(keyword,bdryfile);
      bdryfile += "_bdry.inpt";
   }
   if (bdryfile.substr(0,7) == "${HOME}") {
      filename  = getenv("HOME") +bdryfile.substr(7,bdryfile.length());
   }
   else {
      filename = bdryfile;
   }
   bdrymap.input(filename);
   
   /* Load pointer to block stuff */
   hp_gbl = hp_in;
   
   /* ALLOCATE BOUNDARY CONDITION STUFF */
   if (!coarse) {
      hp_gbl->sbdry_gbls.resize(nsbd);
      hp_gbl->vbdry_gbls.resize(nvbd);
   }
   hp_sbdry.resize(nsbd);
   hp_vbdry.resize(nvbd);
   for(i=0;i<nsbd;++i) hp_sbdry(i) = getnewsideobject(i,bdrymap);
   for(i=0;i<nvbd;++i) hp_vbdry(i) = getnewvrtxobject(i,bdrymap);
   for(i=0;i<nsbd;++i) hp_sbdry(i)->init(bdrymap,hp_gbl->sbdry_gbls(i));
   for(i=0;i<nvbd;++i) hp_vbdry(i)->init(bdrymap,hp_gbl->vbdry_gbls(i));
         
   if (!coarse) {
      setinfo();
      
      /* Allocate time history stuff */
      /* For ease of access have level 0 reference ug */
      ugbd(0).v.reference(ug.v);
      ugbd(0).s.reference(ug.s);
      ugbd(0).i.reference(ug.i);
      vrtxbd(0).reference(vrtx);
      
      /* GET INITIAL CONDITION FUNCTION */
      hp_gbl->ibc = getnewibc(inmap);
      hp_gbl->ibc->input(inmap,idprefix);
      
      /* GET MESH MOVEMENT FUNCTION */
      hp_gbl->mover = getnewmesh_mover(inmap);
      hp_gbl->mover->init(inmap,idprefix);
      
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

      keyword = idprefix +"cfl";
      if (!inmap.getline(keyword,line)) {
         if (!inmap.getline("cfl",line)) {
            line = "2.5 1.5 1.0";
         }
      }
      data.str(line);
      for(i=0;i<log2pmax+1;++i) {
         data >> hp_gbl->cfl(i);
      }
      data.clear();
 

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
      
      vug_frst.resize(maxvst,NV);     
      dres.resize(1);
      dres(0).v.resize(maxvst,NV);
   }
   
   /* UNSTEADY SOURCE TERMS */
   dugdt.resize(log2p+1,maxvst,NV);
   dxdt.resize(log2p+1,maxvst,ND);
   
   if (!coarse && !adapt_storage) {
      keyword = idprefix + ".adapt";
      if (!inmap.get(keyword,adapt_flag)) {
         inmap.getwdefault("adapt",adapt_flag,false);
      }
      *sim::log << "#adapt: " << adapt_flag << std::endl;
   
      if (adapt_flag) {
         inmap.getwdefault("error",trncerr,1.0e-2);
         inmap.getwdefault("bdryangle",bdrysensitivity,5.0);
         bdrysensitivity *= M_PI/180.0;
         inmap.getwdefault("length_tol", vlngth_tol,0.25);

         /* NOW ALLOCATE A COPY SO CAN PERFORM ADAPTATION */
         hp_gbl->pstr = create();
         hp_gbl->pstr->idprefix = idprefix;
         hp_gbl->pstr->mesh::copy(*this);
         keyword = idprefix + ".adapt_storage";
         inmap[keyword] = "1";
         (*hp_gbl->pstr).init(inmap, hp_in);
         inmap[keyword] = "0";
      }
   }
   
   /* RESTART SEQUENCE OR INITIAL CONDITION SEQUENCE */
   if (!coarse) {
      int restartfile;
      if (inmap.get("restart",restartfile)) {
         std::ostringstream nstr;
         std::string fname;
         nstr << restartfile << std::flush;
         fname = idprefix +"rstrt" +nstr.str();
         input(fname);
      } 
      else {
         for(i=0;i<nsbd;++i)
            hp_sbdry(i)->curv_init();
            
         /* USE TOBASIS TO INITALIZE SOLUTION */
         tobasis(hp_gbl->ibc);
      }
   }
   
   return;
}

 tri_hp::~tri_hp() {
   for(int i=0;i<nsbd;++i)
      delete hp_sbdry(i);
}

 void tri_hp::setinfo() {
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

 void tri_hp::maxres() {
   int i,n;
   Array<FLT,1> mxr(NV);
   
   if (mmovement == coupled_deformable) r_mesh::maxres();
   
   for(n=0;n<NV;++n)
      mxr(n) = 0.0;

   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         mxr(n) = MAX(mxr(n),fabs(hp_gbl->res.v(i,n)));
         
   for(n=0;n<NV;++n)
      *sim::log << ' ' << mxr(n) << ' ';
}


