#include "tri_hp_cd.h"
#include "cd_bdry.h"
#include "hp_boundary.h"
#include <myblas.h>
#include <input_map.h>

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_cd_vtype {
   public:
      static const int ntypes = 1;
      enum ids {plain=1};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char tri_hp_cd_vtype::names[ntypes][40] = {"plain"};

hp_vrtx_bdry* tri_hp_cd::getnewvrtxobject(int bnum, input_map *bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   char idntystring[10];
   int type;        
   hp_vrtx_bdry *temp;  
   int idnum = vbdry(bnum)->idnum;
   
   type = idnum&0xffff;

   if (bdrydata) {
      sprintf(idntystring,"v%d",idnum);
      keyword = std::string(idntystring) + ".cd_type";
      
      if ((*bdrydata).get(keyword,val)) {
         type = tri_hp_cd_vtype::getid(val.c_str());
         if (type < 0)  {
            *sim::log << "unknown vertex type:" << val << std::endl;
            exit(1);
         }
      }
   }
   
   switch(type) {
      case tri_hp_cd_vtype::plain: {
         temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
         break;
      }
      default: {
         std::cout << "unrecognizable vrtx type: " <<  type << " idnum: " << idnum << std::endl;
         temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
         break;
      }
   } 
   
   if (bdrydata) temp->input(*bdrydata);
   
   return(temp);
}


/** \brief Helper object for side_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_cd_stype {
   public:
      static const int ntypes = 6;
      enum ids {dirichlet=1,curved_dirichlet,adiabatic,curved_adiabatic,characteristic,curved_characteristic};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char tri_hp_cd_stype::names[ntypes][40] = {"dirichlet","curved_dirichlet","adiabatic","curved_adiabatic","characteristic","curved_characteristic"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_side_bdry* tri_hp_cd::getnewsideobject(int bnum, input_map *bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   char idntystring[10];
   int type;        
   hp_side_bdry *temp;  
   int idnum = sbdry(bnum)->idnum;

   type = idnum&0xff;
   
   if (bdrydata) {
      sprintf(idntystring,"s%d",idnum);
      keyword = std::string(idntystring) + ".cd_type";
      
      if ((*bdrydata).get(keyword,val)) {
         type = tri_hp_cd_stype::getid(val.c_str());
         if (type < 0)  {
            *sim::log << "unknown side type:" << val << std::endl;
            exit(1);
         }
      }
      else {
         *sim::log << "#couldn't find " << keyword << std::endl;
         type = tri_hp_cd_stype::dirichlet;
      }
   }

   switch(type) {
      case tri_hp_cd_stype::dirichlet: {
         temp = new dirichlet_cd<hp_side_bdry>(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_cd_stype::curved_dirichlet: {
         temp = new dirichlet_cd<hp_curved>(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_cd_stype::adiabatic: {
         temp = new neumann_cd<hp_side_bdry>(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_cd_stype::curved_adiabatic: {
         temp = new neumann_cd<hp_curved>(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_cd_stype::characteristic: {
         temp = new char_cd<hp_side_bdry>(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_cd_stype::curved_characteristic: {
         temp = new char_cd<hp_curved>(*this,*sbdry(bnum));
         break;
      }
   }
   
   if (bdrydata) temp->input(*bdrydata);
   
   return(temp);
}
