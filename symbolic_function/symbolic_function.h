#ifndef _symbolic_function_h_
#define _symbolic_function_h_

#include <muParser/muParser.h>
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
	
	symbolic_function(const symbolic_function& tgt) : p(tgt.p), nchildren(tgt.nchildren) {
		std::ostringstream varname;

		/* Reassociate Variables */
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
      	
      	children.resize(nchildren);
      	child_values.resize(nchildren);
      	

		mu::varmap_type variables = p.GetUsedVar();

		/* Reassociate Children */
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
			
			/* This is a child variable */
			children(nchildren) = new symbolic_function<N>(*tgt.children(nchildren));
			p.DefineVar(item->first, &child_values(nchildren));
			++nchildren;
			NEXT: continue;
		 }
	}
	
	void init(input_map& input, std::string idprefix);
		
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
			rslt = 0.0;
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
			rslt = 0.0;
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
			rslt = 0.0;
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
	
	try {
		ptemp.SetExpr(buffer);
	}
	catch (mu::Parser::exception_type &e) {
		std::cout << "Message:  " << e.GetMsg() << std::endl;
		std::cout << "Formula:  " << e.GetExpr() << std::endl;
		std::cout << "Token:     " << e.GetToken() << std::endl;
		std::cout << "Position: " << e.GetPos() << std::endl;
		std::cout << "Errc:      " << e.GetCode() << std::endl;
	}
	mu::varmap_type variables = ptemp.GetUsedVar();



	/* Find out how many children are needed */
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
				++nchildren;
				goto NEXT;
			}
			else {
				std::cout << "couldn't find expression " << item->first << '\n';
				exit(1);
			}
		}
		
		NEXT: continue;
	 }

	children.resize(nchildren);
	child_values.resize(nchildren);
	/* REPEAT EXCEPT THIS TIME ALLOCATE CHILDREN */
	nchildren = 0;
	for (mu::varmap_type::const_iterator item = variables.begin(); item!=variables.end(); ++item) {
		if (N > 1) {
			for (int n=0;n<N;++n) {
				varname.str("");
				varname << "x" << n << std::flush;
				if (item->first == varname.str()) goto NEXT1;
				varname.clear();
			}
		}
		else {
			if (item->first == "x") goto NEXT1;
		}
		if (item->first == "t") goto NEXT1;
		
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
				goto NEXT1;
			}
			else {
				std::cout << "couldn't find expression " << item->first << '\n';
				exit(1);
			}
		}
		p.DefineConst(item->first,value);
		
		NEXT1: continue;
	 }
	 p.SetExpr(buffer);
	
	 return;
 }
	 
	 

	
#endif
