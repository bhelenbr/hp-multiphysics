#ifndef _symbolic_function_h_
#define _symbolic_function_h_

#include <muParser/muParser.h>
#include <blitz/array.h>
#include <input_map.h>
#include <math.h>

template<int N> class symbolic_function {
	protected:
		mu::Parser p;
		mutable blitz::TinyVector<double,N+1> x; 
		int nchildren;
		blitz::Array<symbolic_function<N> *, 1> children;
		mutable blitz::Array<double,1> child_values;
      
	public:		
	symbolic_function() : nchildren(0) {
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
		mu::varmap_type variables = p.GetVar();

		/* Reassociate Children */
		int nchild_temp = 0;
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
			
			/* This is a child variable? */
			children(nchild_temp) = new symbolic_function<N>(*tgt.children(nchild_temp));
			p.DefineVar(item->first, &child_values(nchild_temp));
			++nchild_temp;
			NEXT: continue;
		}
		assert(nchild_temp == nchildren);
	}
	
	void init(input_map& input, std::string idprefix);
		
	double Eval(blitz::TinyVector<double,N> xin, double t = 0.0) const {
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
	double Eval(double xin, double t) const {
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
	double Eval(double t) const {
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
	}
	catch (mu::Parser::exception_type &e) {
		std::cout << "Message:  " << e.GetMsg() << std::endl;
		std::cout << "Formula:  " << e.GetExpr() << std::endl;
		std::cout << "Token:     " << e.GetToken() << std::endl;
		std::cout << "Position: " << e.GetPos() << std::endl;
		std::cout << "Errc:      " << e.GetCode() << std::endl;
	}
	
	return;
}
	 
class vector_function {
protected:
	mu::Parser p;
	int nargs;
	blitz::Array<int, 1> dims;
	blitz::Array<std::string, 1> names;
	blitz::Array<double,2> xargs;
	double time;
	int nchildren;
	blitz::Array<vector_function *, 1> children;
	blitz::Array<double,1> child_values;
	
public:		
	vector_function() : p(), nargs(0), nchildren(0) {
		p.DefineFun("erf", erf, false);
		p.DefineFun("erfc", erfc, false);
		p.DefineVar("t", &time);
	}
	vector_function(const vector_function& tgt) : p(tgt.p), nargs(tgt.nargs), dims(tgt.dims), names(tgt.names), xargs(tgt.xargs), nchildren(tgt.nchildren) {

		/* Reasociate Variables */
		std::ostringstream varname;
		for (int n=0;n<nargs;++n) {
			if (dims(n) > 1) {
				for(int m=0;m<dims(n);++m) {
					varname.str("");
					varname << names(n) << m << std::flush;
					p.DefineVar(varname.str(), &xargs(n,m));
					varname.clear();
				}
			}
			else {
				p.DefineVar(names(n), &xargs(n,0));
			}
		}
		p.DefineVar("t", &time);

		
		/* Reassociate Children */
		children.resize(nchildren);
		child_values.resize(nchildren);
		mu::varmap_type variables = p.GetUsedVar();
		nchildren = 0;
		for (mu::varmap_type::const_iterator item = variables.begin(); item!=variables.end(); ++item) {
			
			for (int n=0;n<nargs;++n) {
				if (dims(n) > 1) {
					for(int m=0;m<dims(n);++m) {
						varname.str("");
						varname << names(n) << m << std::flush;
						if (item->first == varname.str()) goto NEXT;
						varname.clear();
					}
				}
				else {
					if (item->first == names(n)) goto NEXT;
				}
			}
			if (item->first == "t") goto NEXT;
			
			/* This is a child variable */
			children(nchildren) = new vector_function(*tgt.children(nchildren));
			p.DefineVar(item->first, &child_values(nchildren));
			++nchildren;
			NEXT: continue;
		}
	}

	void set_arguments(int nargs_,blitz::Array<int,1> dims_, blitz::Array<std::string,1> names_) {
		nargs = nargs_;
		dims.resize(nargs);
		dims = dims_;
		names.resize(nargs);
		names = names_;
		int maxdim = 1;
		for (int n=0;n<nargs;++n) {
			maxdim = (dims(n)>maxdim ? dims(n) : maxdim);
		}
		xargs.resize(nargs,maxdim);
		std::ostringstream varname;
		for (int n=0;n<nargs;++n) {
			if (dims(n) > 1) {
				for(int m=0;m<dims(n);++m) {
					varname.str("");
					varname << names(n) << m << std::flush;
					p.DefineVar(varname.str(), &xargs(n,m));
					varname.clear();
				}
			}
			else {
				p.DefineVar(names(n), &xargs(n,0));
			}
		}
	}
	vector_function(int nargs_,blitz::Array<int,1> dims_, blitz::Array<std::string,1> names_) {
		p.DefineFun("erf", erf, false);
		p.DefineFun("erfc", erfc, false);
		p.DefineVar("t", &time);
		// vector_function(); // Why doesn't this work????
		set_arguments(nargs_,dims_,names_);
	}
	
	double Eval(blitz::Array<double,1> a0, double t = 0.0) {
		assert(nargs == 1);
		
 		double rslt;
		time = t;
		for (int m = 0; m < dims(0); ++m)
			xargs(0,m) = a0(m);
				
		try {
			for (int n=0;n<nchildren;++n) {
				child_values(n) = children(n)->Eval(a0,t);
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
	
	double Eval(blitz::Array<double,1> a0, blitz::Array<double,1> a1, double t = 0.0) {
		assert(nargs == 2);
		
 		double rslt;
		time = t;
		for (int m = 0; m < dims(0); ++m)
			xargs(0,m) = a0(m);
			
		for (int m = 0; m < dims(1); ++m)
			xargs(1,m) = a1(m);
		
		try {
			for (int n=0;n<nchildren;++n) {
				child_values(n) = children(n)->Eval(a0,a1,t);
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
	
	double Eval(blitz::Array<double,1> a0, blitz::Array<double,1> a1, blitz::Array<double,1> a2, double t = 0.0) {
		assert(nargs == 3);
		
 		double rslt;
		time = t;
		for (int m = 0; m < dims(0); ++m)
			xargs(0,m) = a0(m);
		
		for (int m = 0; m < dims(1); ++m)
			xargs(1,m) = a1(m);
		
		for (int m = 0; m < dims(2); ++m)
			xargs(2,m) = a2(m);
		
		try {
			for (int n=0;n<nchildren;++n) {
				child_values(n) = children(n)->Eval(a0,a1,a2,t);
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
	
	double Eval(blitz::Array<double,1> a0, blitz::Array<double,1> a1, blitz::Array<double,1> a2, blitz::Array<double,1> a3, double t = 0.0) {
		assert(nargs == 4);
		
 		double rslt;
		time = t;
		for (int m = 0; m < dims(0); ++m)
			xargs(0,m) = a0(m);
		
		for (int m = 0; m < dims(1); ++m)
			xargs(1,m) = a1(m);
		
		for (int m = 0; m < dims(2); ++m)
			xargs(2,m) = a2(m);
		
		for (int m = 0; m < dims(3); ++m)
			xargs(3,m) = a3(m);
					
		try {
			for (int n=0;n<nchildren;++n) {
				child_values(n) = children(n)->Eval(a0,a1,a2,a3,t);
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

	double Eval(blitz::Array<double,1> a0, blitz::Array<double,1> a1, blitz::Array<double,1> a2, blitz::Array<double,1> a3, blitz::Array<double,1> a4, double t = 0.0) {
		assert(nargs == 5);
		
 		double rslt;
		time = t;
		for (int m = 0; m < dims(0); ++m)
			xargs(0,m) = a0(m);
		
		for (int m = 0; m < dims(1); ++m)
			xargs(1,m) = a1(m);
		
		for (int m = 0; m < dims(2); ++m)
			xargs(2,m) = a2(m);
		
		for (int m = 0; m < dims(3); ++m)
			xargs(3,m) = a3(m);
			
		for (int m = 0; m < dims(4); ++m)
			xargs(4,m) = a4(m);
		
		try {
			for (int n=0;n<nchildren;++n) {
				child_values(n) = children(n)->Eval(a0,a1,a2,a3,a4,t);
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
	
	void init(input_map& input, std::string idprefix) {
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
			mu::varmap_type variables = ptemp.GetUsedVar();
			
			
			/* Reassociate Children */
			nchildren = 0;
			for (mu::varmap_type::const_iterator item = variables.begin(); item!=variables.end(); ++item) {
				for (int n=0;n<nargs;++n) {
					if (dims(n) > 1) {
						for(int m=0;m<dims(n);++m) {
							varname.str("");
							varname << names(n) << m << std::flush;
							if (item->first == varname.str()) goto NEXT;
							varname.clear();
						}
					}
					else {
						if (item->first == names(n)) goto NEXT;
					}
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
				for (int n=0;n<nargs;++n) {
					if (dims(n) > 1) {
						for(int m=0;m<dims(n);++m) {
							varname.str("");
							varname << names(n) << m << std::flush;
							if (item->first == varname.str()) goto NEXT1;
							varname.clear();
						}
					}
					else {
						if (item->first == names(n)) goto NEXT1;
					}
				}
				if (item->first == "t") goto NEXT1;

				if (!input.get(item->first,value)) {
					/* Check if it is not there or if it is also defined as a formula */
					std::map<std::string,std::string>::const_iterator mi;
					mi = input.find(item->first);
					if (mi != input.end()) {
						/* It is there */
						children(nchildren) = new vector_function(nargs,dims,names);
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
			
			/* TESTING */
		}
		catch (mu::Parser::exception_type &e) {
			std::cout << "Message:  " << e.GetMsg() << std::endl;
			std::cout << "Formula:  " << e.GetExpr() << std::endl;
			std::cout << "Token:     " << e.GetToken() << std::endl;
			std::cout << "Position: " << e.GetPos() << std::endl;
			std::cout << "Errc:      " << e.GetCode() << std::endl;
		}
		
		return;
	}
};

	
#endif
