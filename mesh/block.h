#include <map>
#include <string>
#include <sstream>
#include <ftype.h>
#ifndef _block_h_
#define _block_h_

class boundary;

class block {
   public:
      enum ctrl {stay = 0, advance = 1, stop = 3};
      virtual void init(std::map <std::string,std::string>& input, std::ostream *inlog = 0) = 0;
      virtual void load_const(std::map <std::string,std::string>& input) = 0;
      virtual void alloc(std::map <std::string,std::string>& input) = 0;
      virtual void input(char *filename) = 0;
      virtual void output(char *filename, ftype::name filetype = ftype::easymesh) = 0;
      virtual int comm_entity_size(int grdlvl) = 0;
      virtual int comm_entity_list(int grdlvl, int *list) = 0;
      virtual boundary* vbdry(int grdlvl,int num) = 0;
      virtual boundary* sbdry(int grdlvl,int num) = 0;
      virtual boundary* fbdry(int grdlvl,int num) = 0;
      virtual void maxres() = 0;
      virtual ctrl matchboundaries(int lvl, int excpt) = 0;
      virtual ctrl rsdl(int lvl, int excpt) = 0;
      virtual ctrl vddt(int lvl, int excpt) = 0;
      virtual ctrl update(int lvl, int excpt) = 0;
      virtual ctrl mg_getfres(int lvl, int excpt) = 0;
      virtual ctrl mg_getcchng(int lvl, int excpt) = 0;
      virtual ctrl reconnect(int lvl, int excpt) = 0;
      virtual ctrl tadvance(int lvl, int excpt) = 0;
      virtual ctrl adapt(int excpt) = 0;
};
#endif
