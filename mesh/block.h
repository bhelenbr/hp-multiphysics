#ifndef _block_h_
#define _block_h_

#include <map>
#include <string>
#include <sstream>
#include <blitz/array.h>
#include <input_map.h>

class boundary;

class block {
   public:
      std::string idprefix; 
      int idnum;
      typedef int ctrl;
      static const ctrl begin = -1, stay = 0, advance = 1, advance1 = 3, advance2 = 7, stop = 15;
      block(int idin) : idnum(idin) {
         char buffer[100];
         std::string keyname;
         sprintf(buffer,"b%d",idnum);
         idprefix = std::string(buffer);
      }
      virtual void init(input_map& input) = 0;
      virtual void reload_scratch_pointers() = 0;
      enum output_purpose {display, restart, debug};
      virtual void output(const std::string &filename, output_purpose why, int lvl = 0) = 0;
      virtual ctrl matchboundaries(int lvl, block::ctrl ctrl_message) = 0;
      virtual ctrl rsdl(int lvl, block::ctrl ctrl_message) = 0;
      virtual ctrl setup_preconditioner(int lvl, block::ctrl ctrl_message) = 0;
      virtual ctrl update(int lvl, block::ctrl ctrl_message) = 0;
      virtual double maxres(int lvl = 0) = 0;
      virtual ctrl mg_getfres(int lvl, int lvlm, block::ctrl ctrl_message) = 0;
      virtual ctrl mg_getcchng(int lvl, int lvlp, block::ctrl ctrl_message) = 0;
      virtual ctrl reconnect(int lvl, block::ctrl ctrl_message) = 0;
      virtual ctrl tadvance(int lvl, block::ctrl ctrl_message) = 0;
      virtual ctrl adapt(block::ctrl ctrl_message) = 0;
      
      /* STUFF TO MATCH COMMUNICATION BOUNDARIES */
      virtual int comm_entity_size(int grdlvl) = 0;
      virtual int comm_entity_list(int grdlvl, blitz::Array<int,1>& list) = 0;
      virtual boundary* vbdry(int grdlvl,int num) = 0;
      virtual boundary* sbdry(int grdlvl,int num) = 0;
      virtual boundary* fbdry(int grdlvl,int num) = 0;
      virtual ~block() {}
};
#endif
