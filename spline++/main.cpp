/*
 *  main.cpp
 *  spline++
 *
 *  Created by Brian Helenbrook on 4/17/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "spline.h"

int main() {
	spline<2> myspline;
	std::string file("naca.spl");
	
	myspline.read(file);
	
	TinyVector<double,2> x;
	
	for (double s=0.0;s<2.0;s+=0.001) {
		myspline.interpolate(s,x);
		std::cout << s << ' ' << x(0) << ' ' << x(1) << std::endl;
	}
	
	myspline.interpolate(1.2,x);
	
	double s;
	myspline.find(s,x);
	std::cout << s << ' ' << x << std::endl;	
		
	return(0);
}