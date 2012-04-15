#include "symbolic_function.h"
#include <input_map.h>

#define TESTING

int main (int argc, char *argv[]) {
#ifdef TESTING
	input_map inmap;
	std::string prefix;
	symbolic_function<1> f; 
	symbolic_function<2> e1 ,e2, e4, dhdx0;
	blitz::TinyVector<double,2>  x2d(1.0, 0.75), x2da(0.5, 0.25);
		
	inmap["Re"] = "100.0";
	inmap["l1"] = "37.0";
	inmap["f"] = "exp(-x/l1*Re)";
	inmap["e1"] = "sin(-x0/l1*Re+x1)";
	inmap["e2"] = "x0*e1";
	inmap["dhdx0"] = "2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k-2*s*(x0-ex0)/denom^2*exp(-(x0-ex0)^2/denom^2)*(1+(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)^2)^(1/2)-s*(1-exp(-(x0-ex0)^2/denom^2))/(1+(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)^2)^(1/2)*(2*h*ke^2*x0*exp(-ke^2*x0^2)*cos(k*x0-w*t)-h*(1-exp(-ke^2*x0^2))*sin(k*x0-w*t)*k)*(2*h*ke^2*exp(-ke^2*x0^2)*cos(k*x0-w*t)-4*h*ke^4*x0^2*exp(-ke^2*x0^2)*cos(k*x0-w*t)-4*h*ke^2*x0*exp(-ke^2*x0^2)*sin(k*x0-w*t)*k-h*(1-exp(-ke^2*x0^2))*cos(k*x0-w*t)*k^2)";
	inmap["e4"] = "cos(2*_pi*x1)";
	std::cout << inmap << std::endl;
	inmap.echo = true;

	f.init(inmap,"f");
	std::cout << "f " << f.Eval(1.0) << ' ' << exp(-1.0/37.0*100.0) << std::endl;
 
	e1.init(inmap,"e1");
	std::cout << " e1 eval " << e1.Eval(x2d) << ' ' << sin(-x2d(0)/37.0*100.0 +x2d(1)) << std::endl;
	
	e2.init(inmap,"e2");
	std::cout << " e2 eval " << e2.Eval(x2d) << ' ' << x2d(0)*sin(-x2d(0)/37.0*100.0 +x2d(1)) << std::endl;
	
	symbolic_function<2> e3(e2);
	
	std::cout << " e3 eval " << e3.Eval(x2da) << ' ' << x2da(0)*sin(-x2da(0)/37.0*100.0 +x2da(1)) << std::endl;
	
	e1 = e2;
	std::cout << " = eval " << e1.Eval(x2d) << ' ' << x2d(0)*sin(-x2d(0)/37.0*100.0 +x2d(1)) << std::endl;
	
	
//	dhdx0.init(inmap,"dhdx0");
//	std::cout << "dhdx0 " << dhdx0.Eval(x2d,1.0) << std::endl;
	
	
	/* Vector Function Testing */
	blitz::Array<std::string,1> names;
	blitz::Array<int,1> dims;
	names.resize(3);
	dims.resize(3);
	names(0) = "x";
	names(1) = "u";
	names(2) = "n";
	dims(0) = 2;
	dims(1) = 4;
	dims(2) = 2;
	vector_function vf(3,dims,names);
	
	inmap["Vfunc"] = "V3*u2";
	inmap["V2"] = "u3*x0";
	inmap["V3"] = "n1*t*V2";
	// t*u2*u3*x0
	
	vf.init(inmap,"Vfunc");
	blitz::Array<double,1> x(2);
	x(0) = 0.1;
	x(1) = 0.3;
	double t = 0.7;
	blitz::Array<double,1> u(4);
	u(0) = 1.9;
	u(1) = 2.9;
	u(2) = 3.34;
	u(3) = 5;
	blitz::Array<double,1> n(2);
	n(0) = -100;
	n(1) = -200;
	std::cout << vf.Eval(x,u,n,t) << ' ' << t*u(2)*u(3)*x(0)*n(1) << std::endl;
	

#else
    input_map mymap;
	  
    mymap.input(argv[1]);
	
	symbolic_function<3> f;
	f.init(mymap,"formula");
	
	blitz::TinyVector<double,3> x;
	mymap.get("position",x.data(),3);
	
	double t;
	mymap.getwdefault("time",t,0.0);
	
	std::cout << f.Eval(x) << std::endl;
	
#endif

	return(0);
}
