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
        block(int idin) : idnum(idin) {
            char buffer[100];
            std::string keyname;
            sprintf(buffer,"b%d",idnum);
            idprefix = std::string(buffer);
        }
        virtual void init(input_map& input) = 0;
        enum output_purpose {display, restart, debug};
        virtual void output(const std::string &filename, output_purpose why, int lvl = 0) = 0;
        virtual void matchboundaries(int lvl) = 0;
        virtual void rsdl(int lvl) = 0;
        virtual void setup_preconditioner(int lvl) = 0;
        virtual void update(int lvl) = 0;
        virtual double maxres(int lvl = 0) = 0;
        virtual void mg_getfres(int lvl, int lvlm) = 0;
        virtual void mg_getcchng(int lvl, int lvlp) = 0;
        virtual void reconnect(int lvl) = 0;
        virtual void tadvance(int lvl) = 0;
        virtual void adapt() = 0;
        
        /* STUFF TO MATCH COMMUNICATION BOUNDARIES */
        virtual int comm_entity_size(int grdlvl) = 0;
        virtual int comm_entity_list(int grdlvl, blitz::Array<int,1>& list) = 0;
        virtual boundary* vbdry(int grdlvl,int num) = 0;
        virtual boundary* sbdry(int grdlvl,int num) = 0;
        virtual boundary* fbdry(int grdlvl,int num) = 0;
        virtual ~block() {}
};
#endif
