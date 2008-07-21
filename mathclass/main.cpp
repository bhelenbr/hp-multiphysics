#include "mathclass.h"
#include <input_map.h>

int main() {
	input_map inmap;
	std::string prefix;
	symbolic_function<1> f;
	symbolic_function<2> g;
	symbolic_function<2> g1;
    double x = 1.0;
	blitz::TinyVector<double,2>  x2d(1.0, 0.0);
		
	inmap["Re"] = "100.0";
	inmap["l1"] = "37.0";
	inmap["E1"] = "exp(-x/l1*Re)";
	inmap["E2"] = "sin(-x0/l1*Re+x1)";
	inmap["E3"] = "x0*E2";

	std::cout << inmap << std::endl;

	inmap.echo = true;

	f.init(inmap,"E1");
	std::cout << f.Eval(1.0) << std::endl;

	g.init(inmap,"E2");
	std::cout << g.Eval(x2d) << std::endl;
	
	g1.init(inmap,"E3");
	std::cout << g1.Eval(x2d) << std::endl;
	

//	g1.copy_consts(g);
	

	return(0);
}
