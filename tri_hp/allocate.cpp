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
   int i,ival,p;
   bool adapt_storage, fine_mesh;
   std::string keyword, line;
   std::istringstream data;
   std::string filename;
   
   keyword = idprefix + ".coarse";
   inmap.getwdefault(keyword,coarse,false);
   
   keyword = idprefix + ".adapt_storage";
   inmap.getwdefault(keyword,adapt_storage,false);
   
   fine_mesh = false;
   if (!adapt_storage && !coarse) fine_mesh = true;
      
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
   
   keyword = idprefix + ".nvariable";
   inmap.getwdefault(keyword,NV,1);
   
   keyword = idprefix + ".log2p";
   if (!inmap.get(keyword,log2p)) {
      inmap.getwdefault("log2p",log2p,0);
   }
   if (coarse) log2p = 0;
   
#ifdef AXISYMMETRIC
   *sim::log << "#AXISYMMETRIC is defined" << std::endl;
   int npts = 1;
#else
   *sim::log << "#AXISYMMETRIC is NOT defined" << std::endl;
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
   
   p0 = basis::tri(log2p).p;
   sm0 = basis::tri(log2p).sm;
   im0 = basis::tri(log2p).im;
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
      
   /* Load pointer to block stuff */
   gbl_ptr = hp_in;

   /* ALLOCATE BOUNDARY CONDITION STUFF */
   if (fine_mesh) {
      gbl_ptr->sbdry_gbls.resize(nsbd);
      gbl_ptr->vbdry_gbls.resize(nvbd);
   }
   

   hp_sbdry.resize(nsbd);
   hp_vbdry.resize(nvbd);
   for(i=0;i<nsbd;++i) hp_sbdry(i) = getnewsideobject(i,inmap);
   for(i=0;i<nvbd;++i) hp_vbdry(i) = getnewvrtxobject(i,inmap);
   for(i=0;i<nsbd;++i) hp_sbdry(i)->init(inmap,gbl_ptr->sbdry_gbls(i));
   for(i=0;i<nvbd;++i) hp_vbdry(i)->init(inmap,gbl_ptr->vbdry_gbls(i));

   /* For ease of access have level 0 in time history reference ug */
   ugbd(0).v.reference(ug.v);
   ugbd(0).s.reference(ug.s);
   ugbd(0).i.reference(ug.i);
   vrtxbd(0).reference(vrtx); 
   
   /***************************************/
   /* ALLOCATE JUST ENOUGH FOR ADAPTATION */
   /* NO MORE ADAPT STORAGE AFTER HERE */
   /************************************/
   if (adapt_storage) {
      for(i=1;i<sim::nadapt+1;++i) {
         ugbd(i).v.resize(maxvst,NV);
         ugbd(i).s.resize(maxvst,sm0,NV);
         ugbd(i).i.resize(maxvst,im0,NV);
         vrtxbd(i).resize(maxvst);
      }
      return;
   }

   inmap.getwdefault("hp_fadd",fadd,1.0);
      
   /* Initialize stuff for r_mesh */
   if ((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable)) r_mesh::init(inmap,hp_in);
             
   /* GET MESH MOVEMENT FUNCTION */
   mover = getnewmesh_mover(inmap);
   mover->init(inmap,idprefix);
   
   /* UNSTEADY SOURCE TERMS */
   dugdt.resize(log2p+1,maxvst,NV);
   dxdt.resize(log2p+1,maxvst,ND);
   
   /***************************************/
   /* ALLOCATE FOR COARSE MESH */
   /* NO MORE COARSE MESH AFTER HERE */
   /************************************/
   if (coarse) {
      ugbd(0).v.reference(ug.v);
      ugbd(0).s.reference(ug.s);
      ugbd(0).i.reference(ug.i);
      vrtxbd(0).reference(vrtx);
      
      /* THESE STORE BACKWARDS DIFFERENCE SOURCE TERM AT VERTICES */
      ugbd(1).v.resize(maxvst,NV);
      vrtxbd(1).resize(maxvst);
      
      vug_frst.resize(maxvst,NV);     
      dres.resize(1);
      dres(0).v.resize(maxvst,NV);
      return;
   }
      
   /****************************/
   /* STUFF FOR FINE MESH ONLY */
   /***************************/
   /* GET INITIAL CONDITION FUNCTION */
   gbl_ptr->ibc = getnewibc(inmap);
   gbl_ptr->ibc->input(inmap,idprefix);

   for(i=1;i<sim::nhist+1;++i) {
      ugbd(i).v.resize(maxvst,NV);
      ugbd(i).s.resize(maxvst,sm0,NV);
      ugbd(i).i.resize(maxvst,im0,NV);
      vrtxbd(i).resize(maxvst);
   }

   /* Multigrid Storage all except highest order (log2p+1)*/
   dres.resize(log2p);
   for(i=0;i<log2p;++i) {
      dres(i).v.resize(maxvst,NV);
      dres(i).s.resize(maxvst,basis::tri(i).sm,NV);
      dres(i).i.resize(maxvst,basis::tri(i).im,NV);
   }
         
   /* Allocate block stuff */
   gbl_ptr->ug0.v.resize(maxvst,NV);
   gbl_ptr->ug0.s.resize(maxvst,sm0,NV);
   gbl_ptr->ug0.i.resize(maxvst,im0,NV);
          
   gbl_ptr->res.v.resize(maxvst,NV);
   gbl_ptr->res.s.resize(maxvst,sm0,NV);
   gbl_ptr->res.i.resize(maxvst,im0,NV);
   
   gbl_ptr->res_r.v.resize(maxvst,NV);
   gbl_ptr->res_r.s.resize(maxvst,sm0,NV);
   gbl_ptr->res_r.i.resize(maxvst,im0,NV);

   gbl_ptr->res0.v.resize(maxvst,NV);
   gbl_ptr->res0.s.resize(maxvst,basis::tri(log2p).sm,NV);
   gbl_ptr->res0.i.resize(maxvst,basis::tri(log2p).im,NV); 
   
#ifndef MATRIX_PRECONDITIONER
   gbl_ptr->vprcn.resize(maxvst,NV);
   gbl_ptr->sprcn.resize(maxvst,NV);
   gbl_ptr->tprcn.resize(maxvst,NV);
#else
   gbl_ptr->vprcn.resize(maxvst,NV,NV);
   gbl_ptr->sprcn.resize(maxvst,NV,NV);
   gbl_ptr->tprcn.resize(maxvst,NV,NV);
#endif

   if (!inmap.getline(idprefix +".cfl",line)) inmap.getlinewdefault("cfl",line,"2.5 1.5 1.0"); 
   data.str(line);
   for(i=0;i<log2pmax+1;++i) {
      data >> gbl_ptr->cfl(i);
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
   
   /*********************************/
   /* ALLOCATE ADAPTATION STORAGE   */
   /* BY CALLING THIS ROUTINE AGAIN */
   /*********************************/
   keyword = idprefix + ".adapt";
   if (!inmap.get(keyword,adapt_flag)) {
      inmap.getwdefault("adapt",adapt_flag,false);
   }
   if (adapt_flag) {
      inmap.getwdefault("error",trncerr,1.0e-2);
      inmap.getwdefault("bdryangle",bdrysensitivity,0.1);
      bdrysensitivity *= 180.0/M_PI;
      inmap.getwdefault("length_tol", vlngth_tol,0.25);

      /* NOW ALLOCATE A COPY SO CAN PERFORM ADAPTATION */
      gbl_ptr->pstr = create();
      gbl_ptr->pstr->idprefix = idprefix;
      gbl_ptr->pstr->mesh::copy(*this);
      keyword = idprefix + ".adapt_storage";
      inmap[keyword] = "1";
      gbl_ptr->pstr->tri_hp::init(inmap, hp_in);
      inmap[keyword] = "0";
      
      /* LET EACH BOUNDARY CONDITION DIRECTLY FIND ITS ADAPTATION STORAGE */
      for(i=0;i<nsbd;++i)
         hp_sbdry(i)->adapt_storage = gbl_ptr->pstr->hp_sbdry(i);
         
   }

   /***************************************************/
   /* RESTART SEQUENCE OR INITIAL CONDITION SEQUENCE */
   /**************************************************/
   int restartfile;
   if (inmap.get("restart",restartfile)) {
      std::ostringstream nstr;
      std::string fname;
      nstr << restartfile << std::flush;
      fname = idprefix +"_rstrt" +nstr.str();
      input(fname);
   } 
   else {
      for(i=0;i<nsbd;++i)
         hp_sbdry(i)->curv_init();  /* TEMPO WILL NEED TO CHANGE THIS TO "tobasis" */
         
      /* USE TOBASIS TO INITALIZE SOLUTION */
      tobasis(gbl_ptr->ibc);
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
   
   if (log2p > 0) {
      for(i=0;i<nsbd;++i) {
         if (hp_sbdry(i)->is_curved()) {
            for(j=0;j<sbdry(i)->nel;++j) {
               sind = sbdry(i)->el(j);
               sd(sind).info = 0;
               td(sd(sind).tri(0)).info = 0;
            }
         } 
      }
   }
   
   return;
}

FLT tri_hp::maxres() {
   int i,n;
   Array<FLT,1> mxr(NV);
   FLT mesherror, flowerror;
   
   /* THIS ROUTINE WILL HAVE TO BE OVERWRITTEN TO GET CORRECT DIMENSIONAL NORM FOR EACH SYSTEM */
   if (mmovement == coupled_deformable) mesherror = r_mesh::maxres();
   
   for(n=0;n<NV;++n)
      mxr(n) = 0.0;

   for(i=0;i<nvrtx;++i) {
      for(n=0;n<NV;++n) {
         mxr(n) = MAX(mxr(n),fabs(gbl_ptr->res.v(i,n)));
      }
   }
         
   for(n=0;n<NV;++n)
      *sim::log << ' ' << mxr(n) << ' ';
      
   for(i=0;i<nsbd;++i)
      hp_sbdry(i)->maxres();
      
   flowerror = 0.0;
   for(n=0;n<NV;++n)
      flowerror = MAX(flowerror,mxr(n));
      
   return(flowerror);
   
}


