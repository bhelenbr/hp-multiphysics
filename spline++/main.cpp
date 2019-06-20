/*
 *  main.cpp
 *  spline++
 *
 *  Created by Brian Helenbrook on 4/17/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "spline.h"

#define TESTING

int main(int argc, char** argv) {

	
	if (argc != 4) {
		std::cerr << "usage: spline <filename> <ND> <npts>, outputs uniform interpolation of spline" << std::endl;
	}
	std::string file(argv[1]);
	
	int nd;
	sscanf(argv[2],"%d",&nd);
	
	int npts;
	sscanf(argv[3],"%d",&npts);
	
	switch(nd) {
		case 1: {
			spline<1> myspline;
			spline3<1> myspline3;
			
			myspline.read(file);
			myspline3.read(file);
						
			TinyVector<double,1> x;
			
			for (double s=myspline.start();s<myspline.stop();s+=(myspline.stop()-myspline.start())/npts) {
				myspline.interpolate(s,x);
				std::cout << s << ' ' << x(0);
				myspline3.interpolate(s,x);
				std::cout << ' ' << x(0) << std::endl;
			}
			break;
		}
		case 2: {
			spline<2> myspline;
			spline3<2> myspline3;
			
			myspline.read(file);
			myspline3.read(file);
			
			TinyVector<double,2> x,tan,curv;
			
			for (double s=myspline.start();s<myspline.stop();s+=(myspline.stop()-myspline.start())/npts) {
				myspline.interpolate(s,x);
                myspline.tangent(s, tan);
                myspline.curvature(s, curv);
                std::cout << s << ' ' << x(0) << ' ' << x(1) << ' ' << tan(0) << ' ' << tan(1) << ' ' << curv(0) << ' ' << curv(1) << std::endl;
				myspline3.interpolate(s,x);
                myspline3.tangent(s, tan);
                myspline3.curvature(s, curv);
                std::cout << s << ' ' << x(0) << ' ' << x(1) << ' ' << tan(0) << ' ' << tan(1) << ' ' << curv(0) << ' ' << curv(1) << std::endl;
			}

#ifdef TESTING
			double s;
            x(0) = 0.3;
            x(1) = 0.2;
			myspline.find(s,x);
			std::cout << s << ' ' << x << std::endl;
            myspline3.find(s,x);
            std::cout << s << ' ' << x << std::endl;
#endif
            
			break;
		}
		case 3: {
			spline<3> myspline;
			spline3<3> myspline3;
			
			myspline.read(file);
			myspline3.read(file);
			
			TinyVector<double,3> x;
			
			for (double s=myspline.start();s<myspline.stop();s+=(myspline.stop()-myspline.start())/npts) {
				myspline.interpolate(s,x);
				std::cout << s << ' ' << x(0) << ' ' << x(1) << ' ' << x(2);
				myspline3.interpolate(s,x);
				std::cout << ' ' << x(0) << ' ' << x(1) << ' ' << x(2) << std::endl;
			}
			break;
		}
	}
		
	return(0);
}

