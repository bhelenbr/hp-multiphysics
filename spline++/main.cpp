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
	spline3<2> myspline3;
	
	std::string file("naca.spl");
	
	myspline.read(file);
	myspline3.read(file);
	
	TinyVector<double,2> x;
	
	for (double s=0.0;s<5.0;s+=0.1) {
		myspline.interpolate(s,x);
		std::cout << s << ' ' << x(0) << ' ' << x(1) << std::endl;
	}
	
	for (double s=0.0;s<4.32;s+=0.1) {
		myspline3.interpolate(s,x);
		std::cout << s << ' ' << x(0) << ' ' << x(1) << std::endl;
	}
	
	
	double s;
	myspline.find(s,x);
	std::cout << s << ' ' << x << std::endl;	
		
	return(0);
}

