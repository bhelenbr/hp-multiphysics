/*
 *  hp_initial_conditions.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include <symbolic_function.h>

class symbolic_ibc : public init_bdry_cndtn {
    private:
        Array<symbolic_function<2>,1> fcn;
    public:
        FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
            return(fcn(n).Eval(x,time));
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
                nstr << idnty << n << std::flush;
                if (inmap.find(nstr.str()) != inmap.end()) {
                    fcn(n).init(inmap,nstr.str());
                }
                else {
					std:cerr << "couldn't find initial condition function\n";
					exit(1);
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



init_bdry_cndtn *tri_hp::getnewibc(std::string suffix, input_map& inmap) {
    std::string keyword,ibcname;
	init_bdry_cndtn *temp;
    int type;

    /* FIND INITIAL CONDITION TYPE */
	keyword = gbl->idprefix + "_" +suffix;
	if (!inmap.get(keyword,ibcname)) {
		keyword = suffix;
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find initial condition type" << std::endl;
		}
	}
    
    type = ibc_type::getid(ibcname.c_str());
        
    switch(type) {
        case ibc_type::symbolic: {
            temp = new symbolic_ibc;
			break;
        }
        default: {
            *gbl->log << "couldn't find initial condition function " << ibcname << std::endl;
            exit(1);
        }
    }
	temp->input(inmap,keyword);
	return(temp);

}
        
        
class translating : public tri_hp_helper {
    public: 
        tri_hp &x;
        TinyVector<FLT,2> velocity;
        translating(tri_hp &xin) :tri_hp_helper(xin), x(xin) {}

        virtual void init(input_map& input, std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
            std::ostringstream nstr;
            
            FLT vdflt[2] = {0.0, 0.0};
            if (!input.get(idnty +"_velocity",velocity.data(),2)) input.getwdefault("velocity",velocity.data(),2,vdflt);
        }
        
        void tadvance() {
            
            if (x.coarse_level) return;
            
            if (x.gbl->substep == 0) x.l2error(x.gbl->ibc);
            
            TinyVector<FLT,2> dx;
            for (int n=0;n<tri_mesh::ND;++n)
                dx(n) = x.gbl->cdirk[x.gbl->substep]/x.gbl->dti*velocity(n);
            
            for(int i=0;i<x.npnt;++i)
                x.pnts(i) += dx;
            
            return;
        }
};

class gcl_test : public tri_hp_helper {
    public: 
        tri_hp &x;
        Array<symbolic_function<2>,1> vel;
        gcl_test(tri_hp &xin) : tri_hp_helper(xin), x(xin) {}

        virtual void init(input_map& input, std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
            std::ostringstream nstr;
            
            vel.resize(tri_mesh::ND);
            
            for(int n=0;n<tri_mesh::ND;++n) {
                nstr.str("");
                nstr << idnty << "_gcl_velocity" << n << std::flush;
                if (input.find(nstr.str()) != input.end()) {
                    vel(n).init(input,nstr.str());
                }
                else {
                    nstr.str("");
                    nstr << "gcl_velocity" << n << std::flush;
                    if (input.find(nstr.str()) != input.end()) {
                        vel(n).init(input,nstr.str());
                    }
                    else {
                        *x.gbl->log << "couldn't find mesh movement function\n";
                        exit(1);
                    }
                }
            }
        }
        
        void tadvance() {
            if (x.coarse_level) return;
            
            if (x.gbl->substep == 0) x.l2error(x.gbl->ibc);
                
            TinyVector<FLT,2> dx;
            
            for(int i=0;i<x.npnt;++i) {
                for (int n=0;n<tri_mesh::ND;++n)
                    dx(n) = x.gbl->cdirk[x.gbl->substep]/x.gbl->dti*vel(n).Eval(x.pnts(i),x.gbl->time);
                x.pnts(i) += dx;
            }            
            return;
        }
};

class l2_error : public tri_hp_helper {
    public: 
        tri_hp &x;
        l2_error(tri_hp &xin) :tri_hp_helper(xin), x(xin) {}
        void tadvance() {
            if (x.coarse_level) return;
            
            if (x.gbl->substep == 0) {
                x.l2error(x.gbl->ibc);
            }
            return;
        }
};

class print_averages : public tri_hp_helper {
    public: 
        tri_hp &x;
        print_averages(tri_hp &xin) :tri_hp_helper(xin), x(xin) {}
        void tadvance() {
            
            if (x.coarse_level) return;
            
            if (x.gbl->substep == 0) {
                Array<FLT,1> av(x.NV+3);
                x.integrated_averages(av);
                *x.gbl->log << "# area: " << av(0) << " xbar: " << av(1) << " ybar: " << av(2);
                for(int n=0;n<x.NV;++n) 
                    *x.gbl->log << " V" << n << ": " << av(n+3);
                *x.gbl->log << std::endl;
            }
            
            return;
        }
};

class helper_type {
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
const char helper_type::names[ntypes][40] = {"translating","print_averages","l2error","gcl_test"};


tri_hp_helper *tri_hp::getnewhelper(input_map& inmap) {
    std::string movername;
    int type;
    
    /* FIND INITIAL CONDITION TYPE */
    if (!inmap.get(gbl->idprefix + "_helper",movername))
		inmap.getwdefault("tri_hp_helper",movername,std::string("default"));

	type = helper_type::getid(movername.c_str());
        
    switch(type) {
        case helper_type::translating: {
            tri_hp_helper *temp = new translating(*this);
            return(temp);
        }
        case helper_type::print_averages: {
            tri_hp_helper *temp = new print_averages(*this);
            return(temp);
        }
        case helper_type::l2error: {
            tri_hp_helper *temp = new l2_error(*this);
            return(temp);
        }
        case helper_type::gcl_test: {
            tri_hp_helper *temp = new gcl_test(*this);
            return(temp);
        }
        default: {
            tri_hp_helper *temp = new tri_hp_helper(*this);
            return(temp);
        }
    }
}

