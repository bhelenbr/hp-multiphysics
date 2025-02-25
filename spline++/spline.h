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


template<int ND> class spline {
private:
    int npts;
    blitz::Array<double,1> x;
    blitz::Array<blitz::TinyVector<double,ND>,1> y;
    blitz::Array<blitz::TinyVector<blitz::TinyVector<double,ND>,6>,1> c;
    
public:
    spline() : npts(0) {
#ifdef BZ_DEBUG
    std::cerr << "#spline: BZ_DEBUG is set\n";
#endif
#ifdef DEBUG
    std::cerr << "#spline: Running in Xcode's DEBUG Mode\n";
#endif
    }
    spline(const spline<ND>& inpt) : npts(inpt.npts) {
        x.resize(npts);
        y.resize(npts);
        c.resize(npts);
        x = inpt.x;
        y = inpt.y;
        c = inpt.c;
    }
    int read(std::string filename);
    int interpolate(const double spt, blitz::TinyVector<double,ND>& loc) const;
    int offset(const double spt, const double distance, blitz::TinyVector<double,ND>& loc) const;
    int tangent(const double spt, blitz::TinyVector<double,ND>& tan) const;
    int curvature(const double spt, blitz::TinyVector<double,ND>& curv) const;
    int find(double &spt, blitz::TinyVector<double,ND>& loc) const;
    int find_with_guess(double &spt, blitz::TinyVector<double,ND>& loc) const;
    double start() const {return(x(0));}
    double stop() const {return(x(npts-1));}
    int size() const {return(npts);}
};

template<int ND> class spline3 {
private:
    int npts;
    blitz::Array<double,1> x;
    blitz::Array<blitz::TinyVector<double,ND>,1> y;
    blitz::Array<blitz::TinyVector<blitz::TinyVector<double,ND>,4>,1> c;
    
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
    int interpolate(const double spt, blitz::TinyVector<double,ND>& loc) const;
    int offset(const double spt, const double distance, blitz::TinyVector<double,ND>& loc) const;
    int tangent(const double spt, blitz::TinyVector<double,ND>& tan) const;
    int curvature(const double spt, blitz::TinyVector<double,ND>& curv) const;
    int find(double &spt, blitz::TinyVector<double,ND>& loc) const;
    int find_with_guess(double &spt, blitz::TinyVector<double,ND>& loc) const;
    double start() const {return(x(0));}
    double stop() const {return(x(npts-1));}
    int size() const {return(npts);}
};

namespace spline_functions2D {
void transform2D(blitz::TinyVector<double,2>& xpt, const double size, const double angle, const blitz::TinyVector<double,2> offset);
void transform2Di(blitz::TinyVector<double,2>& xpt, const double size, const double angle, const blitz::TinyVector<double,2> offset);
void interpolate(blitz::TinyVector<double,2>& loc, blitz::TinyVector<double,2>& tan, blitz::TinyVector<double,2>& curv, const spline<2>& myspline, double s, const double size, const double angle, const blitz::TinyVector<double,2> offset, double norm_dist);
int find(const blitz::TinyVector<double,2>& loc, const spline<2>& myspline, double& s, const double size, const double angle,  const blitz::TinyVector<double,2> offset, double &norm_dist);
int find_with_guess(const blitz::TinyVector<double,2>& loc, const spline<2>& myspline, double& s, const double size, const double angle,  const blitz::TinyVector<double,2> offset, double &norm_dist);
}



#endif
