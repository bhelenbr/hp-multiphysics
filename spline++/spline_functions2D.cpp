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

const int ND = 2;

void spline_functions2D::transform2D(TinyVector<double,ND>& loc, const double size, const double angle, const TinyVector<double,ND> offset) {
    // xpt is assumed in physical space and being moved to spline space;
    loc -= offset;
    TinyVector<double,ND> temp;
    temp(0) = cos(-angle)*loc(0) -sin(-angle)*loc(1);
    temp(1) = sin(-angle)*loc(0) +cos(-angle)*loc(1);
    loc = temp/size;
}

void spline_functions2D::transform2Di(TinyVector<double,ND>& loc, const double size, const double angle, const TinyVector<double,ND> offset) {
    // xpt is assumed in spline space and being moved to physical space;
    loc = loc*size;
    TinyVector<double,ND> temp;
    temp(0) = cos(angle)*loc(0) -sin(angle)*loc(1);
    temp(1) = sin(angle)*loc(0) +cos(angle)*loc(1);
    loc = temp +offset;
}

/* p = x(s) +n x k *norm_dist */
/* dp/ds = dx/ds +curv tan * norm_dist */
/* dp/dn = norm */
void spline_functions2D::interpolate(TinyVector<double,ND>& loc, TinyVector<double,ND>& tan, TinyVector<double,ND>& curv, const spline<ND>& myspline, double s, const double size, const double angle, const TinyVector<double,ND> offset, double norm_dist) {
    TinyVector<double,ND> zero = 0.0;
    myspline.offset(s,norm_dist/size,loc);
    transform2Di(loc,size,angle,offset);
    myspline.tangent(s, tan);
    transform2Di(tan,size,angle,zero);
    myspline.curvature(s, curv);
    transform2Di(curv,size,angle,zero);
}

void spline_functions2D::find(TinyVector<double,ND>& loc, const spline<ND>& myspline, double& s, const double size, const double angle,  const TinyVector<double,ND> offset, const double norm_dist) {
    
    TinyVector<double,ND> tan, curv;
    transform2D(loc, size, angle, offset);
    myspline.find(s, loc);
    interpolate(loc, tan, curv, myspline, s, size, angle, offset, norm_dist);
}
