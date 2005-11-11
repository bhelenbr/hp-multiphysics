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

vrtx_bdry* mesh::getnewvrtxobject(int idnum, std::map<std::string,std::string> *bdrydata) {
   std::string keyword;
   std::istringstream data;
   std::map<std::string,std::string>::const_iterator mi;
   char idntystring[10];
   int type;        
   vrtx_bdry *temp;  
   
   type = idnum&0xffff;

   if (bdrydata) {
      sprintf(idntystring,"v%d",idnum);
      keyword = std::string(idntystring) + ".type";
      mi = (*bdrydata).find(keyword);
      if (mi != (*bdrydata).end()) {
         type = vtype::getid((*mi).second.c_str());
         if (type < 0)  {
            *sim::log << "unknown vertex type:" << (*mi).second << std::endl;
            exit(1);
         }
      }
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
   
   if (bdrydata) temp->input(*bdrydata);
   
   temp->output(*sim::log);

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
      static const int ntypes = 9;
      enum ids {plain=1, comm, prdc, sinewave, circle, spline, partition, naca, gaussian};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char stype::names[ntypes][40] = {"plain", "comm", "prdc", "sinewave", "circle", "spline", "partition","naca","gaussian"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
side_bdry* mesh::getnewsideobject(int idnum, std::map<std::string,std::string> *bdrydata) {
   std::string keyword;
   std::istringstream data;
   std::map<std::string,std::string>::const_iterator mi;
   char idntystring[10];
   int type;        
   side_bdry *temp;  
   
   type = idnum&0xff;

   if (bdrydata) {
      sprintf(idntystring,"s%d",idnum);
      keyword = std::string(idntystring) + ".type";
      mi = (*bdrydata).find(keyword);
      if (mi != (*bdrydata).end()) {
         type = stype::getid((*mi).second.c_str());
         if (type < 0)  {
            *sim::log << "unknown side type:" << (*mi).second << std::endl;
            exit(1);
         }
      }
      else {
         *sim::log << "#couldn't find type for side: " << idnum << std::endl;
         type = stype::plain;
      }
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
         temp = new sinewave(idnum,*this);
         break;
      }
      case stype::circle: {
         temp = new circle(idnum,*this);
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
         temp = new naca(idnum,*this);
         break;
      }
      case stype::gaussian: {
         temp = new gaussian(idnum,*this);
         break;
      }
      default: {
         temp = new side_bdry(idnum,*this);
         std::cout << "unrecognizable side type: " << idnum << "type " << type << std::endl;
         break;
      }
   }
   
   if (bdrydata) temp->input(*bdrydata);
   
   temp->output(*sim::log);

   return(temp);
}

   
   
