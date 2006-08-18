/*
 *  mvpttobdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "boundaries.h"

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class vtype {
   public:
      static const int ntypes = 3;
      enum ids {plain=1,comm,prdc};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char vtype::names[ntypes][40] = {"plain","comm","prdc"};

vrtx_bdry* mesh::getnewvrtxobject(int idnum, input_map& bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   ostringstream nstr;
   int type;        
   vrtx_bdry *temp;  

   nstr.str("");
   nstr << idnum << std::flush;      
   keyword = idprefix +".v" +nstr.str() + ".type";
   if (bdrydata.get(keyword,val)) {
      type = vtype::getid(val.c_str());
      if (type < 0)  {
         *sim::log << "unknown vertex type:" << val << std::endl;
         exit(1);
      }
   }
   else {
      type = vtype::plain;
   }
      
   switch(type) {
      case vtype::plain: {
         temp = new vrtx_bdry(idnum,*this);
         break;
      }
      case vtype::comm: {
         temp = new vcomm(idnum,*this);
         break;
      }
      case vtype::prdc: {
         temp = new vprdc(idnum,*this);
         break;
      }
      default: {
         std::cout << "unrecognizable vrtx type: " <<  type << " idnum: " << idnum << std::endl;
         temp = new vrtx_bdry(idnum,*this);
         break;
      }
   } 
   
   temp->input(bdrydata);
   
   return(temp);
}


/** \brief Helper object for side_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class stype {
   public:
      static const int ntypes = 18;
      enum ids {plain=1, comm, prdc, sinewave, circle, spline, partition, naca, gaussian,parabola,hyperbola,coupled_sinewave,coupled_circle,coupled_sinewave_comm,coupled_circle_comm,coupled_parabola,coupled_hyperbola,coupled_hyperbola_comm};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char stype::names[ntypes][40] = {"plain", "comm", "prdc", "sinewave", "circle", "spline",
   "partition","naca","gaussian","parabola","hyperbola","coupled_sinewave","coupled_circle",
	"coupled_sinewave_comm","coupled_circle_comm","coupled_parabola","coupled_hyperbola",
	"coupled_hyperbola_comm"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
side_bdry* mesh::getnewsideobject(int idnum, input_map& bdrydata) {
   std::string keyword,val;
   std::istringstream data;
   ostringstream nstr;
   int type;        
   side_bdry *temp;  
   
   type = stype::plain;

   nstr.str("");
   nstr << idnum << std::flush;      
   keyword = idprefix +".s" +nstr.str() + ".type";
   if (bdrydata.get(keyword,val)) {
      type = stype::getid(val.c_str());
      if (type < 0)  {
         *sim::log << "unknown side type:" << val << std::endl;
         exit(1);
      }
   }
   else {
      type = stype::plain;
   }

   switch(type) {
      case stype::plain: {
         temp = new side_bdry(idnum,*this);
         break;
      }
      case stype::comm: {
         temp = new scomm(idnum,*this); 
         break;
      }
      case stype::prdc: {
         temp = new sprdc(idnum,*this);
         break;
      }
      case stype::sinewave: {
         temp = new analytic_geometry<side_bdry,sinewave>(idnum,*this);
         break;
      }
      case stype::circle: {
         temp = new analytic_geometry<side_bdry,circle>(idnum,*this);
         break;
      }
      case stype::spline: {
    //     temp = new spline(idnum,*this);
         break;
      }
      case stype::partition: {
         temp = new spartition(idnum,*this);
         break;
      }
      case stype::naca: {
         temp = new analytic_geometry<side_bdry,naca>(idnum,*this);
         break;
      }
      case stype::gaussian: {
         temp = new analytic_geometry<side_bdry,gaussian>(idnum,*this);
         break;
      }
		case stype::parabola: {
         temp = new analytic_geometry<side_bdry,parabola>(idnum,*this);
         break;
      }
		case stype::hyperbola: {
         temp = new analytic_geometry<side_bdry,hyperbola>(idnum,*this);
         break;
      }
      case stype::coupled_sinewave: {
         temp = new ssolution_geometry<analytic_geometry<side_bdry,sinewave> >(idnum,*this);
         break;
      }
      case stype::coupled_circle: {
         temp = new ssolution_geometry<analytic_geometry<side_bdry,circle> >(idnum,*this);
         break;
      }
      case stype::coupled_sinewave_comm: {
         temp = new ssolution_geometry<analytic_geometry<scomm,sinewave> >(idnum,*this);
         break;
      }
      case stype::coupled_circle_comm: {
         temp = new ssolution_geometry<analytic_geometry<scomm,circle> >(idnum,*this);
         break;
      }
		case stype::coupled_parabola: {
         temp = new ssolution_geometry<analytic_geometry<side_bdry,parabola> >(idnum,*this);
         break;
      }
		case stype::coupled_hyperbola: {
         temp = new ssolution_geometry<analytic_geometry<side_bdry,hyperbola> >(idnum,*this);
         break;
      }
		case stype::coupled_hyperbola_comm: {
         temp = new ssolution_geometry<analytic_geometry<scomm,hyperbola> >(idnum,*this);
         break;
      }		
      
      default: {
         temp = new side_bdry(idnum,*this);
         std::cout << "unrecognizable side type: " << idnum << "type " << type << std::endl;
         break;
      }
   }
   
   temp->input(bdrydata);
   
   return(temp);
}

   
   
