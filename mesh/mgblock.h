#include "r_mesh.h"
#include "block.h"

#include <map>
#include <string>
#include <sstream>
#include <iostream>

#ifndef _r_block_h_
#define _r_block_h_


/* GENERIC MULTIGRID BLOCK */

template<class GRD> class mgrid : public block {
    protected:
        int ngrid;
        GRD *grd;
        typedef typename GRD::filetype outtype;
        typedef typename GRD::transfer gtrans;
        typedef typename GRD::gbl ggbl;
        Array<Array<gtrans,1>,1> cv_to_ft;
        Array<Array<gtrans,1>,1> fv_to_ct;
    public:
        ggbl gstorage;
        bool adapt_flag;
        FLT tolerance;
    
    public:
        mgrid(int idnum) : block(idnum), adapt_flag(0) {}
        void init(input_map& input);
        void output(const std::string &filename, block::output_purpose why, int level = 0) {
            std::string fapp;
            fapp = idprefix +"_" +filename;
            grd[level].output(fapp,why);
        }
        void reconnect(int lvl);
        void matchboundaries(int lvl) {return(grd[lvl].matchboundaries());}
      
        /* PASS THROUGH ROUTINES FOR MATCHING BOUNDARIES*/
        int comm_entity_size(int grdlvl) {
            return(grd[grdlvl].comm_entity_size());
        }
        int comm_entity_list(int grdlvl, Array<int,1>& list) {
            return(grd[grdlvl].comm_entity_list(list));
        }
        boundary* vbdry(int grdlvl, int num) {
            return(grd[grdlvl].vbdry(num));
        }
        boundary* sbdry(int grdlvl, int num) {
            return(grd[grdlvl].sbdry(num));
        }
        boundary* fbdry(int grdlvl, int num) {
            return(0);
            // return(grd[grdlvl].fbdry[num]);
        }
        
        void tadvance(int lvl) {
            if (lvl == 0) {
                grd[lvl].tadvance(0,fv_to_ct(0),cv_to_ft(0),0);
            } else {
                grd[lvl].tadvance(1,fv_to_ct(lvl-1),cv_to_ft(lvl-1), &grd[lvl-1]);
            }
        }
        void rsdl(int lvl) {
            grd[lvl].rsdl();
        }
        FLT maxres(int lvl = 0) {return(grd[lvl].maxres());}
        void setup_preconditioner(int lvl) {
            grd[lvl].setup_preconditioner();
        }
        void update(int lvl) {
            grd[lvl].update();
        }
        void mg_getfres(int lvl, int lvlm) {
            grd[lvl].mg_getfres(fv_to_ct(lvlm),cv_to_ft(lvlm),&grd[lvlm]);
        }
        void mg_getcchng(int lvl, int lvlp) {
            grd[lvl].mg_getcchng(fv_to_ct(lvl), cv_to_ft(lvl), &grd[lvlp]);
        }
        void adapt();
};

template<class GRD> void mgrid<GRD>::init(input_map& input) {
    std::string keyword;
    std::istringstream data;
    std::string filename;
    std::string bdryfile;

    /* LOAD NUMBER OF GRIDS */
    input.getwdefault("ngrid",ngrid,1);
    cv_to_ft.resize(ngrid);
    fv_to_ct.resize(ngrid);
    grd = new GRD[ngrid];

    keyword = idprefix + "_adapt";
    if (!input.get(keyword,adapt_flag)) {
        input.getwdefault("adapt",adapt_flag,false);
    }
    
    keyword = idprefix + "_tolerance";
    if (!input.get(keyword,tolerance)) {
        input.getwdefault("tolerance",tolerance,1.9);
    }
    
    grd[0].idprefix = idprefix;
    keyword = idprefix + "_coarse";
    input[keyword] = "0";
    grd[0].init(input,&gstorage);
    grd[0].setinfo();
    
    input[keyword] = "1";
#define OLDRECONNECT
    for(int lvl=1;lvl<ngrid;++lvl) {
        grd[lvl].idprefix = idprefix;
#ifdef OLDRECONNECT
        grd[lvl].mesh::init(grd[lvl-1],2.0);
#else
        FLT size_reduce = 1.0;
        if (lvl > 1) size_reduce = 2.0;
        grd[lvl].mesh::init(grd[lvl-1],size_reduce);
#endif
        grd[lvl].init(input,&gstorage);
        cv_to_ft(lvl-1).resize(grd[lvl].maxvst);
        fv_to_ct(lvl-1).resize(grd[lvl-1].maxvst);
    }
    input[keyword] = "0";
    
    return;
}


template<class GRD> void mgrid<GRD>::reconnect(int lvl) {
    std::string name,fname;
    std::ostringstream nstr;
        
    name = idprefix + "_coarsen";
    
#ifdef OLDRECONNECT
    grd[lvl].coarsen(1.6,grd[lvl-1]);
#else
    FLT size_reduce = 1.0;
    if (lvl > 1) size_reduce = 2.0;
    grd[lvl].coarsen2(1.5,grd[lvl-1],size_reduce);
#endif

    grd[lvl].mgconnect(cv_to_ft(lvl-1),grd[lvl-1]);
    grd[lvl-1].mgconnect(fv_to_ct(lvl-1),grd[lvl]);
    
    /* THIS IS FOR DIAGNOSIS OF MULTI-GRID FAILURES */
    grd[lvl].checkintegrity();
    if (sim::adapt_output) {
        std::string adapt_file;
        std::ostringstream nstr;
        nstr.str("");
        nstr << sim::tstep << std::flush;
        name = idprefix +"_coarsen" +nstr.str() +".";
        nstr.str("");
        nstr << lvl << flush;
        fname = name +nstr.str();
        grd[lvl].mesh::output(fname,mesh::grid);
        fname = name +nstr.str() + "_ft_to_cv";
        grd[lvl-1].testconnect(fname,fv_to_ct(lvl-1),&grd[lvl]);
        fname = name +nstr.str() + "_cv_to_ft";
        grd[lvl].testconnect(fname,cv_to_ft(lvl-1),&grd[lvl-1]);
    }

    return;
}

template<class GRD> void mgrid<GRD>::adapt() {
    int last_phase, mp_phase;
    if (!adapt_flag) return;
 
    grd[0].length();
    
    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        grd[0].vmsgload(boundary::all_phased,mp_phase, boundary::symmetric,grd[0].vlngth.data(),0,0,1);
        grd[0].vmsgpass(boundary::all_phased,mp_phase, boundary::symmetric);
        last_phase = true;
        last_phase &= grd[0].vmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average, grd[0].vlngth.data(),0,0,1);
    }
    
    grd[0].adapt(tolerance);
    
    return;
}

#endif
