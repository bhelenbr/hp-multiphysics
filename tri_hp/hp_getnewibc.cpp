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
			*gbl->log << "couldn't find initial condition type " << keyword << std::endl;
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
#ifdef DIRK
			for (int n=0;n<tri_mesh::ND;++n)
				dx(n) = x.gbl->cdirk(x.gbl->substep)/x.gbl->dti*velocity(n);
#else
			for (int n=0;n<tri_mesh::ND;++n)
				dx(n) = velocity(n)/x.gbl->dti;
#endif

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

			FLT dt;
#ifdef DIRK
			dt = x.gbl->cdirk(x.gbl->substep)/x.gbl->dti;
#else
			dt = 1./x.gbl->dti;
#endif                
			TinyVector<FLT,tri_mesh::ND> dx;
			for(int i=0;i<x.npnt;++i) {
				for (int n=0;n<tri_mesh::ND;++n)
					dx(n) = dt*vel(n).Eval(x.pnts(i),x.gbl->time);
				x.pnts(i) += dx;
			}			
			return;
		}
};

class l2_error : public tri_hp_helper {
	public: 
		tri_hp &x;
		l2_error(tri_hp &xin) :tri_hp_helper(xin), x(xin) {}
		void output() {
			x.l2error(x.gbl->ibc);
		}
};

class print_averages : public tri_hp_helper {
	public: 
		tri_hp &x;
		print_averages(tri_hp &xin) :tri_hp_helper(xin), x(xin) {}
		void output() {            
			Array<FLT,1> av(x.NV+3);
			x.integrated_averages(av);
			*x.gbl->log << "# area: " << av(0) << " xbar: " << av(1) << " ybar: " << av(2);
			for(int n=0;n<x.NV;++n) 
				*x.gbl->log << " V" << n << ": " << av(n+3);
			*x.gbl->log << std::endl;

			return;
		}
};

class output_contour : public tri_hp_helper {
	protected:
		int var;
		FLT c;
		symbolic_function<2> norm;


	public: 
		tri_hp &x;
		output_contour(tri_hp &xin) :tri_hp_helper(xin), x(xin) {}

		virtual void init(input_map& input, std::string idnty) {
			input.getwdefault(idnty+"_var",var,0);
			input.getwdefault(idnty+"_level",c,0.0);
			norm.init(input,idnty+"_norm");
		}

		void output() {


			Array<FLT,1> u(x.NV), du(x.NV);
			TinyVector<FLT,2> xpt;
			FLT psi,dpsi;
			int iter,sind,v0,v1;
			int nintsct = 0;
			FLT norm_sum = 0.0;

			std::ofstream out;
			std:ostringstream filename;

			filename << "contour" << x.gbl->tstep << ".dat";
			out.open(filename.str().c_str());


			/* FIND CONTOUR EDGE INTERSECTION POINTS */
			for (sind=0;sind < x.nseg;++sind) {
				v0 = x.seg(sind).pnt(0);
				v1 = x.seg(sind).pnt(1);
				if ( (x.ug.v(v0,var)-c)*(x.ug.v(v1,var)-c) < 0.0) {
					/* Have intersection */
					psi = -1. +2*(c -x.ug.v(v0,var))/(x.ug.v(v1,var)-x.ug.v(v0,var));
					x.ugtouht1d(sind);
					for (iter = 0, dpsi = 1.0; iter < 100 && fabs(dpsi) > 10*FLT_EPSILON; ++iter) {
						basis::tri(x.log2p).ptprobe1d(x.NV,u.data(),du.data(),psi,&x.uht(0)(0),MXTM);
						dpsi = -(u(var)-c)/du(var);
						psi += dpsi;
					}
					if (fabs(psi) > 1.0 || iter > 99) continue;

					x.crdtocht1d(sind);
					for(int n=0;n<x.ND;++n) {
						basis::tri(x.log2p).ptprobe1d(x.ND,xpt.data(),&x.cht(0,0),MXTM);
						out << xpt(n) << ' ';
					}
					out << '\n';
					norm_sum += norm.Eval(xpt,x.gbl->time);
					++nintsct;
				}

			}
			out.close();
			*x.gbl->log << "Contour norm: " << sqrt(norm_sum)/nintsct << '\n';            
			return;
		}
};



class helper_type {
	public:
		const static int ntypes = 5;
		enum ids {translating,print_averages,l2error,output_contour,gcl_test};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};
const char helper_type::names[ntypes][40] = {"translating","print_averages","l2error","output_contour","gcl_test"};


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
		case helper_type::output_contour: {
			tri_hp_helper *temp = new output_contour(*this);
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

