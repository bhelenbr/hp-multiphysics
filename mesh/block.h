#include <map>
#include <string>
#include <sstream>
#include <ftype.h>
#include <blitz/array.h>
#ifndef _block_h_
#define _block_h_

class boundary;

class block {
   public:
      std::string idprefix; 
      int idnum;
      enum ctrl {stay = 0, advance = 1, stop = 3};
      block(int idin) : idnum(idin) {
         char buffer[100];
         std::string keyname;
         sprintf(buffer,"b%d",idnum);
         idprefix = std::string(buffer);
      }
      virtual void init(std::map <std::string,std::string>& input) = 0;
      virtual void reload_scratch_pointers() = 0;
      virtual void output(char *filename, ftype::name filetype = ftype::easymesh) = 0;
      virtual ctrl matchboundaries(int lvl, int excpt) = 0;
      virtual ctrl rsdl(int lvl, int excpt) = 0;
      virtual ctrl setup_preconditioner(int lvl, int excpt) = 0;
      virtual ctrl update(int lvl, int excpt) = 0;
      virtual void maxres() = 0;
      virtual ctrl mg_getfres(int lvl, int excpt) = 0;
      virtual ctrl mg_getcchng(int lvl, int excpt) = 0;
      virtual ctrl reconnect(int lvl, int excpt) = 0;
      virtual void coarsenchk(const char *fname) = 0;
      virtual ctrl tadvance(int lvl, int excpt) = 0;
      virtual ctrl adapt(int excpt) = 0;
      
      /* STUFF TO MATCH COMMUNICATION BOUNDARIES */
      virtual int comm_entity_size(int grdlvl) = 0;
      virtual int comm_entity_list(int grdlvl, blitz::Array<int,1>& list) = 0;
      virtual boundary* vbdry(int grdlvl,int num) = 0;
      virtual boundary* sbdry(int grdlvl,int num) = 0;
      virtual boundary* fbdry(int grdlvl,int num) = 0;
      virtual ~block() {}
};

class btype {
   public:
      enum ids {plain=1};
};
#endif
