/*
 *  hp_initial_conditions.cpp
 *  tet_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"
#include <symbolic_function.h>

class symbolic_ibc : public init_bdry_cndtn {
	private:
		Array<symbolic_function<tet_mesh::ND>,1> fcn;
	public:
		FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
			return(fcn(n).Eval(x,time));
		}
		void init(input_map &inmap, std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
			std::ostringstream nstr;
			int nvar;

			/* Figure out how many variables */
			for(nvar=0, nstr.str(""), nstr << idnty << nvar << std::flush; inmap.find(nstr.str()) != inmap.end(); ++nvar, nstr.str(""), nstr << idnty << nvar << std::flush);
		
			fcn.resize(nvar);
			
			for(int n=0;n<nvar;++n) {
				nstr.str("");
				nstr << idnty << n << std::flush;
				if (inmap.find(nstr.str()) != inmap.end()) {
						fcn(n).init(inmap,nstr.str());
				}
				else {
					std::cerr << "couldn't find initial condition function\n";
					exit(1);
				}
			}
		}
};

class ibc_type {
	public:
		const static int ntypes = 1;
		enum ids {unknown=-1,symbolic};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};
const char ibc_type::names[ntypes][40] = {"symbolic"};



init_bdry_cndtn *tet_hp::getnewibc(std::string ibcname) {
	init_bdry_cndtn *temp = NULL;
	int type;

	type = ibc_type::getid(ibcname.c_str());
	if (type == ibc_type::unknown) {
		*gbl->log << "unknown ibc type:" << ibcname << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	switch(type) {
		case ibc_type::symbolic: {
			temp = new symbolic_ibc;
			break;
		}
		default: {
			*gbl->log << "couldn't find initial condition function " << ibcname << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	return(temp);

}
		
		
class translating : public tet_hp_helper {
	public: 
		tet_hp &x;
		TinyVector<FLT,tet_mesh::ND> velocity;
		translating(tet_hp &xin) :tet_hp_helper(xin), x(xin) {}

		virtual void init(input_map& inmap, std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
			std::ostringstream nstr;

			FLT vdflt[tet_mesh::ND] = {0.0, 0.0, 0.0};
			if (!inmap.get(idnty +"_velocity",velocity.data(),3)) inmap.getwdefault("velocity",velocity.data(),2,vdflt);
		}
		
		void tadvance() {

			if (x.coarse_level) return;

			if (x.gbl->substep == 0) x.l2error(x.hp_gbl->ibc);

			TinyVector<FLT,tet_mesh::ND> dx;
#ifdef DIRK
			for (int n=0;n<tet_mesh::ND;++n)
				dx(n) = x.gbl->cdirk(x.gbl->substep)/x.gbl->dti*velocity(n);
#else
			for (int n=0;n<tet_mesh::ND;++n)
				dx(n) = velocity(n)/x.gbl->dti;
#endif

			for(int i=0;i<x.npnt;++i)
				x.pnts(i) += dx;

			return;
		}
};

class gcl_test : public tet_hp_helper {
	public: 
		tet_hp &x;
		Array<symbolic_function<tet_mesh::ND>,1> vel;
		gcl_test(tet_hp &xin) : tet_hp_helper(xin), x(xin) {}

		virtual void init(input_map& inmap, std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
			std::ostringstream nstr;

			vel.resize(tet_mesh::ND);

			for(int n=0;n<tet_mesh::ND;++n) {
				nstr.str("");
				nstr << idnty << "_gcl_velocity" << n << std::flush;
				if (inmap.find(nstr.str()) != inmap.end()) {
						vel(n).init(inmap,nstr.str());
				}
				else {
						nstr.str("");
						nstr << "gcl_velocity" << n << std::flush;
						if (inmap.find(nstr.str()) != inmap.end()) {
							vel(n).init(inmap,nstr.str());
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

			if (x.gbl->substep == 0) x.l2error(x.hp_gbl->ibc);

			FLT dt;
			dt = x.gbl->cdirk(x.gbl->substep)/x.gbl->dti;          
			TinyVector<FLT,tet_mesh::ND> dx;
			for(int i=0;i<x.npnt;++i) {
				for (int n=0;n<tet_mesh::ND;++n)
						dx(n) = dt*vel(n).Eval(x.pnts(i),x.gbl->time);
				x.pnts(i) += dx;
			}			
			return;
		}
};

class l2_error : public tet_hp_helper {
	public: 
		tet_hp &x;
		l2_error(tet_hp &xin) :tet_hp_helper(xin), x(xin) {}
		void output() {
			x.l2error(x.hp_gbl->ibc);
		}
};

class helper_type {
	public:
		const static int ntypes = 4;
		enum ids {unknown=-1,plain,translating,l2error,gcl_test};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};
const char helper_type::names[ntypes][40] = {"plain","translating","l2error","gcl_test"};


tet_hp_helper *tet_hp::getnewhelper(std::string helpername) {
	std::string movername;
	int type;
	
	type = helper_type::getid(helpername.c_str());
	if (type == helper_type::unknown) {
		*gbl->log << "unknown helper type:" << helpername << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
		
	switch(type) {
		case helper_type::translating: {
			tet_hp_helper *temp = new translating(*this);
			return(temp);
		}
		case helper_type::l2error: {
			tet_hp_helper *temp = new l2_error(*this);
			return(temp);
		}
		case helper_type::gcl_test: {
			tet_hp_helper *temp = new gcl_test(*this);
			return(temp);
		}
		default: {
			tet_hp_helper *temp = new tet_hp_helper(*this);
			return(temp);
		}
	}
}

