#ifndef _mathclass_h_
#define _mathclass_h_

#include <muParser.h>
#include <blitz/array.h>
#include <input_map.h>
#include <math.h>

template<int N> class symbolic_function {
	protected:
		mu::Parser p;
		blitz::TinyVector<double,N+1> x; 
		int nchildren;
		blitz::Array<symbolic_function<N> *, 1> children;
		blitz::Array<double,1> child_values;
      
	public:		
	symbolic_function() {
		std::ostringstream varname;
		
		p.DefineFun("erf", erf, false);
		p.DefineFun("erfc", erfc, false);
		
		if (N > 1) {
			for (int n=0;n<N;++n) {
				varname.str("");
				varname << "x" << n << std::flush;
				p.DefineVar(varname.str(), &x(n));
				varname.clear();
			}
		}
		else if (N == 1) {
			p.DefineVar("x", &x(0));
		}
      p.DefineVar("t",&x(N));
	}
	
	void init(input_map& input, std::string idprefix);
	
	 void copy(const symbolic_function<N>& in) {
	     p = in.p;
	 }
     
    void copy_consts(symbolic_function<N>& in) {
     mu::valmap_type cmap = in.p.GetConst();
     if (cmap.size()) {
         mu::valmap_type::const_iterator item = cmap.begin();
         for (; item!=cmap.end(); ++item) {
             p.DefineConst(item->first ,item->second);
         }
     }
    }
     
     
		
	
	double Eval(blitz::TinyVector<double,N> xin, double t = 0.0) {
 		double rslt;
		for (int n = 0; n < N; ++n)
			x(n) = xin(n);
		
		x(N) = t;

		try {
			for (int n=0;n<nchildren;++n) {
				child_values(n) = children(n)->Eval(xin,t);
			}
			rslt = p.Eval();
		}
		catch (mu::Parser::exception_type &e) {
			std::cout << "Message:  " << e.GetMsg() << std::endl;
			std::cout << "Formula:  " << e.GetExpr() << std::endl;
			std::cout << "Token:    " << e.GetToken() << std::endl;
			std::cout << "Position: " << e.GetPos() << std::endl;
			std::cout << "Errc:     " << e.GetCode() << std::endl;
		}
		return(rslt);
	}
	
	/* Special case for one-spatial variable function */
	double Eval(double xin, double t) {
		double rslt;

		x(0) = xin;
		x(N) = t;
		try {
			for (int n=0;n<nchildren;++n) {
				child_values(n) = children(n)->Eval(xin,t);
			}
			rslt = p.Eval();
		}
		catch (mu::Parser::exception_type &e) {
			std::cout << "Message:  " << e.GetMsg() << std::endl;
			std::cout << "Formula:  " << e.GetExpr() << std::endl;
			std::cout << "Token:    " << e.GetToken() << std::endl;
			std::cout << "Position: " << e.GetPos() << std::endl;
			std::cout << "Errc:     " << e.GetCode() << std::endl;
		}
		return(rslt);
	}

	/* Special case for zero-variable functions */
	double Eval(double t) {
		double rslt;

		x(0) = t;
		try {
			for (int n=0;n<nchildren;++n) {
				child_values(n) = children(n)->Eval(t);
			}
			rslt = p.Eval();
		}
		catch (mu::Parser::exception_type &e) {
			std::cout << "Message:  " << e.GetMsg() << std::endl;
			std::cout << "Formula:  " << e.GetExpr() << std::endl;
			std::cout << "Token:    " << e.GetToken() << std::endl;
			std::cout << "Position: " << e.GetPos() << std::endl;
			std::cout << "Errc:     " << e.GetCode() << std::endl;
		}
		return(rslt);
	}
		
};	

template<int N> void symbolic_function<N>::init(input_map& input, std::string idprefix) {
	std::ostringstream conststring, varname;
	std::istringstream constname;
	std::string buffer,name,keyword;
	double value;
	
	/* LOAD CONSTANTS IN FORMULA */
	mu::Parser ptemp;
	ptemp.DefineFun("erf", erf, false);
	ptemp.DefineFun("erfc", erfc, false);
  
	if (!input.getline(idprefix,buffer)) {
	   std::cout << "couldn't find expression" << idprefix << '\n';
	}
	ptemp.SetExpr(buffer);
	mu::varmap_type variables = ptemp.GetUsedVar();
	children.resize(variables.size());  // This is too big
	child_values.resize(variables.size());
	nchildren = 0;
	
	for (mu::varmap_type::const_iterator item = variables.begin(); item!=variables.end(); ++item) {
		if (N > 1) {
			for (int n=0;n<N;++n) {
				varname.str("");
				varname << "x" << n << std::flush;
				if (item->first == varname.str()) goto NEXT;
				varname.clear();
			}
		}
		else {
			if (item->first == "x") goto NEXT;
		}
		if (item->first == "t") goto NEXT;
		
		if (!input.get(item->first,value)) {
			/* Check if it is not there or if it is also defined as a formula */
			std::map<std::string,std::string>::const_iterator mi;
    		mi = input.find(item->first);
			if (mi != input.end()) {
				/* It is there */
				children(nchildren) = new symbolic_function<N>;
				children(nchildren)->init(input,item->first);
				p.DefineVar(item->first, &child_values(nchildren));
				++nchildren;
				goto NEXT;
			}
			else {
				std::cout << "couldn't find expression " << item->first << '\n';
				exit(1);
			}
		}
		p.DefineConst(item->first,value);
		
		NEXT: continue;
	 }
	 p.SetExpr(buffer);
	
	 return;
 }
	 
	 

	
#endif
