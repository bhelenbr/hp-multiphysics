/*
 *  getnewins_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_ins.h"
#include "bdry_ins.h"
#include <input_map.h>

using namespace bdry_ins;

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_ins_vtype {
   public:
      static const int ntypes = 5;
      enum ids {unknown=-1,surface_inflow,surface_periodic,surface_outflow,surface_outflow_planar,inflow};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i);
         return(unknown);
      }
};

const char tri_hp_ins_vtype::names[ntypes][40] = {"surface_inflow","surface_periodic","surface_outflow","surface_outflow_planar","inflow"};

hp_vrtx_bdry* tri_hp_ins::getnewvrtxobject(int bnum, input_map &bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   int type;        
   hp_vrtx_bdry *temp;  
   
   keyword = vbdry(bnum)->idprefix + "_ins_type";
   if (bdrydata.get(keyword,val)) {
      type = tri_hp_ins_vtype::getid(val.c_str());
      if (type == tri_hp_ins_vtype::unknown)  {
         *sim::log << "unknown vertex type:" << val << std::endl;
         exit(1);
      }
   }
   else {
      type = tri_hp_ins_vtype::unknown;
   }
   
   
   switch(type) {
      case tri_hp_ins_vtype::surface_inflow: {
         temp = new surface_fixed_pt(*this,*vbdry(bnum));
         break;
      }
      case tri_hp_ins_vtype::surface_periodic: {
         temp = new surface_periodic_pt(*this,*vbdry(bnum));
         break;
      }
      case tri_hp_ins_vtype::surface_outflow: {
         temp = new surface_outflow_endpt(*this,*vbdry(bnum));
         break;
      }
      case tri_hp_ins_vtype::surface_outflow_planar: {
         temp = new surface_outflow_planar(*this,*vbdry(bnum));
         break;
      }
      case tri_hp_ins_vtype::inflow: {
         temp = new inflow_pt(*this,*vbdry(bnum));
         break;
      }
      default: {
         temp = tri_hp::getnewvrtxobject(bnum,bdrydata);
         break;
      }
   } 
   return(temp);
}


/** \brief Helper object for side_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_ins_stype {
   public:
      static const int ntypes = 10;
      enum ids {unknown=-1,plain,inflow,outflow,characteristic,euler,
         symmetry,applied_stress,surface,surface_slave,hybrid_surface_levelset};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i);
         return(-1);
      }
};

const char tri_hp_ins_stype::names[ntypes][40] = {"plain","inflow","outflow","characteristic","euler",
   "symmetry","applied_stress","surface","surface_slave","hybrid_surface_levelset"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_side_bdry* tri_hp_ins::getnewsideobject(int bnum, input_map &bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   int type;        
   hp_side_bdry *temp;  
   

   keyword =  sbdry(bnum)->idprefix + "_ins_type";
   if (bdrydata.get(keyword,val)) {
      type = tri_hp_ins_stype::getid(val.c_str());
      if (type == tri_hp_ins_stype::unknown)  {
         *sim::log << "unknown side type:" << val << std::endl;
         exit(1);
      }
   }
   else {
      type = tri_hp_ins_stype::unknown;
   }

   switch(type) {
      case tri_hp_ins_stype::plain: {
         temp = new generic(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::inflow: {
         temp = new inflow(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::outflow: {
         temp = new neumann(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::characteristic: {
         temp = new characteristic(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::euler: {
         temp = new euler(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::symmetry: {
         temp = new symmetry(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::applied_stress: {
         temp = new applied_stress(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_ins_stype::surface: {
         temp = new surface(*this,*sbdry(bnum));
         dynamic_cast<sgeometry_pointer *>(sbdry(bnum))->solution_data = temp;
         break;
      }
      case tri_hp_ins_stype::surface_slave: {
         temp = new surface_slave(*this,*sbdry(bnum));
         dynamic_cast<sgeometry_pointer *>(sbdry(bnum))->solution_data = temp;
         break;
      }
      case tri_hp_ins_stype::hybrid_surface_levelset: {
         temp = new hybrid_surface_levelset(*this,*sbdry(bnum));
         break;
      }
      default: {
         temp = tri_hp::getnewsideobject(bnum,bdrydata);
         break;
      }
   }   
   return(temp);
}


