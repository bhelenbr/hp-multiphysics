/*
 *  spline.h
 *  spline++
 *
 *  Created by Brian Helenbrook on 4/16/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _spline_h_
#define _spline_h_

#include <blitz/array.h>
using namespace blitz;

template<int ND> class spline {
	private:
		int npts;
		Array<double,1> x;
		Array<TinyVector<double,ND>,1> y;
		Array<TinyVector<TinyVector<double,ND>,6>,1> c;
		
	public:
		spline() : npts(0) {}
		spline(const spline<ND>& inpt) : npts(inpt.npts) {
			x.resize(npts);
			y.resize(npts);
			c.resize(npts);
			x = inpt.x;
			y = inpt.y;
			c = inpt.c;
		}
		int read(std::string filename);
		int interpolate(const double spt, TinyVector<double,ND>& loc);
		int find(double &spt, TinyVector<double,ND>& loc);
};

template<int ND> class spline3 {
	private:
		int npts;
		Array<double,1> x;
		Array<TinyVector<double,ND>,1> y;
		Array<TinyVector<TinyVector<double,ND>,4>,1> c;
		
	public:
		spline3() : npts(0) {}
		spline3(const spline3<ND>& inpt) : npts(inpt.npts) {
			x.resize(npts);
			y.resize(npts);
			c.resize(npts);
			x = inpt.x;
			y = inpt.y;
			c = inpt.c;
		}
		int read(std::string filename);
		int interpolate(const double spt, TinyVector<double,ND>& loc);
		int find(double &spt, TinyVector<double,ND>& loc);
};

#endif
