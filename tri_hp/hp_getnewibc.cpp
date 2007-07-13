/*
 *  hp_initial_conditions.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include <mathclass.h>

class symbolic_ibc : public init_bdry_cndtn {
    private:
        Array<symbolic_function<2>,1> fcn;
    public:
        FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            return(fcn(n).Eval(x,sim::time));
        }
        void input(input_map &inmap,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
            std::ostringstream nstr;
            int nvar;
            
            keyword = idnty +"_nvariable";

            if (!inmap.get(keyword,nvar))
                inmap.getwdefault("nvariable",nvar,1);
            
            fcn.resize(nvar);

            for(int n=0;n<nvar;++n) {
                nstr.str("");
                nstr << idnty << "_ibc" << n << std::flush;
                if (inmap.find(nstr.str() +"_expression") != inmap.end()) {
                    fcn(n).init(inmap,nstr.str());
                }
                else {
                    nstr.str("");
                    nstr << "ibc" << n << std::flush;
                    if (inmap.find(nstr.str() +"_expression") != inmap.end()) {
                        fcn(n).init(inmap,nstr.str());
                    }
                    else {
                        *sim::log << "couldn't find initial condition function\n";
                        exit(1);
                    }
                }
            }
        }
};

class ibc_type {
    public:
        const static int ntypes = 1;
        enum ids {symbolic};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i);
            return(-1);
        }
};
const char ibc_type::names[ntypes][40] = {"symbolic"};



init_bdry_cndtn *tri_hp::getnewibc(input_map& inmap) {
    std::string keyword,ibcname;
    int type;

    /* FIND INITIAL CONDITION TYPE */
    keyword = std::string(idprefix) + "_ibc";
    if (!inmap.get(keyword,ibcname)) {
        if (!inmap.get("ibc",ibcname)) {
            *sim::log << "couldn't find initial condition type" << std::endl;
        }
    }
    
    type = ibc_type::getid(ibcname.c_str());
        
    switch(type) {
        case ibc_type::symbolic: {
            init_bdry_cndtn *temp = new symbolic_ibc;
            return(temp);
        }
        default: {
            *sim::log << "couldn't find initial condition function " << ibcname << std::endl;
            exit(1);
        }
    }
}
        
        
class translating : public mesh_mover {
    public: 
        tri_hp &x;
        TinyVector<FLT,2> velocity;
        translating(tri_hp &xin) :mesh_mover(xin), x(xin) {}

        virtual void init(input_map& input, std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
            std::ostringstream nstr;
            
            FLT vdflt[2] = {0.0, 0.0};
            if (!input.get(idnty +"_velocity",velocity.data(),2)) input.getwdefault("velocity",velocity.data(),2,vdflt);
        }
        
        block::ctrl tadvance(block::ctrl ctrl_message) {
            
            if (x.coarse) return(block::stop);
            
            if (ctrl_message == block::begin) {
                if (sim::substep == 0) x.l2error(x.gbl_ptr->ibc);
                
                TinyVector<FLT,2> dx;
                for (int n=0;n<mesh::ND;++n)
                    dx(n) = sim::cdirk[sim::substep]/sim::dti*velocity(n);
                
                for(int i=0;i<x.nvrtx;++i)
                    x.vrtx(i) += dx;
            }
            
            return(block::stop);
        }
};

class gcl_test : public mesh_mover {
    public: 
        tri_hp &x;
        Array<symbolic_function<2>,1> vel;
        gcl_test(tri_hp &xin) :mesh_mover(xin), x(xin) {}

        virtual void init(input_map& input, std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
            std::ostringstream nstr;
            
            vel.resize(mesh::ND);
            
            for(int n=0;n<mesh::ND;++n) {
                nstr.str("");
                nstr << idnty << "_gcl_velocity" << n << std::flush;
                if (input.find(nstr.str() +"_expression") != input.end()) {
                    vel(n).init(input,nstr.str());
                }
                else {
                    nstr.str("");
                    nstr << "gcl_velocity" << n << std::flush;
                    if (input.find(nstr.str() +"_expression") != input.end()) {
                        vel(n).init(input,nstr.str());
                    }
                    else {
                        *sim::log << "couldn't find mesh movement function\n";
                        exit(1);
                    }
                }
            }
        }
        
        block::ctrl tadvance(block::ctrl ctrl_message) {
            
            if (x.coarse) return(block::stop);
            
            if (ctrl_message == block::begin) {
                if (sim::substep == 0) x.l2error(x.gbl_ptr->ibc);
                
                TinyVector<FLT,2> dx;
                
                for(int i=0;i<x.nvrtx;++i) {
                    for (int n=0;n<mesh::ND;++n)
                        dx(n) = sim::cdirk[sim::substep]/sim::dti*vel(n).Eval(x.vrtx(i),sim::time);
                    x.vrtx(i) += dx;
                }
            }
            
            return(block::stop);
        }
};

class l2_error : public mesh_mover {
    public: 
        tri_hp &x;
        l2_error(tri_hp &xin) :mesh_mover(xin), x(xin) {}
        block::ctrl tadvance(block::ctrl ctrl_message) {
            
            if (x.coarse) return(block::stop);
            
            if (ctrl_message == block::begin && sim::substep == 0) {
                x.l2error(x.gbl_ptr->ibc);
            }
            return(block::stop);
        }
};

class print_averages : public mesh_mover {
    public: 
        tri_hp &x;
        print_averages(tri_hp &xin) :mesh_mover(xin), x(xin) {}
        block::ctrl tadvance(block::ctrl ctrl_message) {
            
            if (x.coarse) return(block::stop);
            
            if (ctrl_message == block::begin && sim::substep == 0) {
                Array<FLT,1> av(x.NV+3);
                x.integrated_averages(av);
                *sim::log << "# area: " << av(0) << " xbar: " << av(1) << " ybar: " << av(2);
                for(int n=0;n<x.NV;++n) 
                    *sim::log << " V" << n << ": " << av(n+3);
                *sim::log << std::endl;
            }
            
            return(block::stop);
        }
};

class mesh_mover_type {
    public:
        const static int ntypes = 4;
        enum ids {translating,print_averages,l2error,gcl_test};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i);
            return(-1);
        }
};
const char mesh_mover_type::names[ntypes][40] = {"translating","print_averages","l2error","gcl_test"};


mesh_mover *tri_hp::getnewmesh_mover(input_map& inmap) {
    std::string keyword,movername;
    int type;
    
    /* FIND INITIAL CONDITION TYPE */
    keyword = std::string(idprefix) + "_mesh_mover";
    if (!inmap.get(keyword,movername)) {
        if (!inmap.get("mesh_mover",movername)) {
            type = -1;
        }
    }
    
    type = mesh_mover_type::getid(movername.c_str());
        
    switch(type) {
        case mesh_mover_type::translating: {
            mesh_mover *temp = new translating(*this);
            return(temp);
        }
        case mesh_mover_type::print_averages: {
            mesh_mover *temp = new print_averages(*this);
            return(temp);
        }
        case mesh_mover_type::l2error: {
            mesh_mover *temp = new l2_error(*this);
            return(temp);
        }
        case mesh_mover_type::gcl_test: {
            mesh_mover *temp = new gcl_test(*this);
            return(temp);
        }
        default: {
            mesh_mover *temp = new mesh_mover(*this);
            return(temp);
        }
    }
}

