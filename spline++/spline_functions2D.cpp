//
//  spline_functions2D.cpp
//  spline++
//
//  Created by Brian Helenbrook on 7/14/19.
//

#include <stdio.h>
#include "spline.h"

using namespace spline_functions2D;
using namespace blitz;

void spline_functions2D::transform2D(TinyVector<double,2>& xpt, const double size, const double angle, const TinyVector<double,2> offset) {
    // xpt is assumed in physical space and being moved to spline space;
    const int ND = 2;
    xpt -= offset;
    TinyVector<double,ND> temp;
    temp(0) = cos(-angle)*xpt(0) -sin(-angle)*xpt(1);
    temp(1) = sin(-angle)*xpt(0) +cos(-angle)*xpt(1);
    xpt = temp/size;
}

void spline_functions2D::transform2Di(TinyVector<double,2>& xpt, const double size, const double angle, const TinyVector<double,2> offset) {
    // xpt is assumed in spline space and being moved to physical space;
    const int ND = 2;
    xpt = xpt*size;
    TinyVector<double,ND> temp;
    temp(0) = cos(angle)*xpt(0) -sin(angle)*xpt(1);
    temp(1) = sin(angle)*xpt(0) +cos(angle)*xpt(1);
    xpt = temp +offset;
}

void spline_functions2D::interpolate(const spline<2>& myspline, double s, const double size, const double angle, const TinyVector<double,2> dx, double norm_dist) {
    const int ND = 2;
    TinyVector<double,ND> x, tan, curv, zero = 0.0;
    myspline.interpolate(s,x);
    transform2Di(x,size,angle,dx);
    myspline.tangent(s, tan);
    transform2Di(tan,size,angle,zero);
    TinyVector<double,ND> norm;
    norm(0) = tan(1);
    norm(1) = -tan(0);
    norm /= sqrt(dot(norm,norm));
    myspline.curvature(s, curv);
    transform2Di(curv,size,angle,zero);
    std::cout << std::setprecision(10) << s << ' ' << x(0) +norm(0)*norm_dist << ' ' << x(1) +norm(1)*norm_dist << ' ' << tan(0) << ' ' << tan(1) << ' ' << curv(0) << ' ' << curv(1) << std::endl;
}
