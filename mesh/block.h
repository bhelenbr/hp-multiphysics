#include <map>
#include <string>
#include <sstream>
#include <ftype.h>
#ifndef _block_h_
#define _block_h_

class boundary;

class block {
   public:
      enum control_state {stay = 0, advance = 1, stop = 3};
      virtual void init(std::map <std::string,std::string>& input) = 0;
      virtual void load_const(std::map <std::string,std::string>& input) = 0;
      virtual void alloc(std::map <std::string,std::string>& input) = 0;
      virtual void input(char *filename) = 0;
      virtual void output(char *filename, FTYPE filetype = easymesh) = 0;
      virtual int comm_entity_size(int grdlvl) = 0;
      virtual int comm_entity_list(int grdlvl, int *list) = 0;
      virtual boundary* vbdry(int grdlvl,int num) = 0;
      virtual boundary* sbdry(int grdlvl,int num) = 0;
      virtual boundary* fbdry(int grdlvl,int num) = 0;
      virtual void maxres() = 0;
      virtual control_state matchboundaries(int lvl, int excpt) = 0;
      virtual control_state rsdl(int lvl, int excpt) = 0;
      virtual control_state vddt(int lvl, int excpt) = 0;
      virtual control_state update(int lvl, int excpt) = 0;
      virtual control_state mg_getfres(int lvl, int excpt) = 0;
      virtual control_state mg_getcchng(int lvl, int excpt) = 0;
      virtual control_state reconnect(int lvl, int excpt) = 0;
      virtual control_state tadvance(int lvl, int excpt) = 0;
      virtual control_state adapt(int excpt) = 0;
};
#endif
