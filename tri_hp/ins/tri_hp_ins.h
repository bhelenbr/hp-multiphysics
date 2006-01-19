/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_ins_h_
#define _tri_hp_ins_h_

#include "tri_hp.h"
#include <blocks.h>

#ifdef AXISYMMETRIC
#define RAD(I,J) crd(0)(I,J)
#define RAD1D(I) crd(0)(0,I)
#else
#define RAD(I,J) 1
#define RAD1D(I) 1
#endif

class tri_hp_ins : public tri_hp {
   public:
      /* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
      struct gbl : public tri_hp::gbl {
         /* STABILIZATION */
         Array<FLT,1> tau,delt;
            
         /* PHYSICAL CONSTANTS */
         FLT rho, mu, nu;
         
         /* STORAGE FOR CALCULATION OF ENERGY AND AREA */
         TinyVector<FLT,2> eanda, eanda_recv;

      } *ins_gbl;
      
      FLT adis; // DISSIPATION CONSTANT
      
      hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map *bdrydata);
      hp_side_bdry* getnewsideobject(int bnum, input_map *bdrydata);
      
   public:
      void init(input_map& input, gbl *gin); 
      tri_hp_ins* create() { return new tri_hp_ins(); }
      block::ctrl length(int excpt);
      block::ctrl setup_preconditioner(int excpt);
      block::ctrl rsdl(int excpt, int stage=sim::NSTAGE);
      void calculate_unsteady_sources(bool coarse);
      
      /* BOUNDARY CONDITION ROUTINES */
      // void drag(int bdry_id);
      // void integrated_averages(FLT a[]);
      // void addbflux();
      // block::ctrl tadvance(bool coarse,int execpoint,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp_ins *fmesh);

};
#endif