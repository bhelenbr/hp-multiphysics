//
//  mapped.cpp
//  tri_mesh
//
//  Created by Brian Helenbrook on 6/6/24.
//

#include "mapped_mesh.h"

void spline_mapping::init(input_map& input, std::string idprefix, std::ostream *log) {
    trsfm.init(input,idprefix);
    std::string line;
    if (!input.get(idprefix+"_spline",line)) {
        *log << "Couldn't fine spline file name in input file\n";
        sim::abort(__LINE__,__FILE__,log);
    }
    my_spline.read(line);
    input.getwdefault(idprefix+"_scale",scale,1.0);
}

void spline_mapping::to_physical_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    TinyVector<FLT,2> tan, curv;
    spline_functions2D::interpolate(to, tan, curv, my_spline, from(0), scale, trsfm.theta,trsfm.pos, -from(1));
}

void spline_mapping::calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {
    TinyVector<FLT,2> pnt, tan, curv;
    spline_functions2D::interpolate(pnt, tan, curv, my_spline, loc(0), scale, trsfm.theta,trsfm.pos, -loc(1));
    /* p = x(s) +n*norm_dist */
    /* dp/ds = dx/ds +curv * norm_dist */
    /* dp/dn = norm */
    
    /* Derivatives with respect to s*/
    jacobian(0,0) = tan(0) -curv(0)*loc(1);
    jacobian(1,0) = tan(1) -curv(1)*loc(1);
    /* Derivaties with respect to norm_dist */
    jacobian(0,1) = -tan(1);
    jacobian(1,1) = +tan(0);
}


void polar_mapping::init(input_map& input, std::string idprefix, std::ostream *log) {
    if (!input.get(idprefix+"_pnt",pnt.data(),2)) {
        *log << "Couldn't read location of polar pnt " << idprefix+"_pnt" << std::endl;;
        sim::abort(__LINE__,__FILE__,log);
    }
    if (!input.get(idprefix+"_theta_length",theta_length)) {
        *log << "Couldn't read length to scale theta " << idprefix+"_theta_length" << std::endl;;
        sim::abort(__LINE__,__FILE__,log);
    }
}

void polar_mapping::to_physical_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    const FLT r = from(1);
    const FLT theta = -from(0)/theta_length;
    to(0) = pnt(0) +r*cos(theta);
    to(1) = pnt(1) +r*sin(theta);
}

void polar_mapping::calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {
    const FLT r = loc(1);
    const FLT theta = -loc(0)/theta_length;
    
    /* Derivatives with respect to theta*length */
    jacobian(0,0) = +r*sin(theta)/theta_length;  // dx/dt
    jacobian(1,0) = -r*cos(theta)/theta_length; // dy/dt
    /* Derivaties with respect to r */
    jacobian(0,1) = cos(theta); // dx/dr
    jacobian(1,1) = sin(theta); // dy/dr
}

void polar_log_mapping::init(input_map& input, std::string idprefix, std::ostream *log) {
    polar_mapping::init(input,idprefix,log);
    
    if (!input.get(idprefix+"_r0",r0)) {
        *log << "Couldn't read r0 " << idprefix+"_r0" << std::endl;;
        sim::abort(__LINE__,__FILE__,log);
    }
    input.getwdefault(idprefix+"_r_eps",r_eps,DBL_EPSILON);
    r_eps = r_eps*r0;
}

void polar_log_mapping::to_physical_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    const FLT r = exp(from(1))*(r0+r_eps) -r_eps;
    const FLT theta = -from(0)/theta_length;
    to(0) = pnt(0) +r*cos(theta);
    to(1) = pnt(1) +r*sin(theta);
}

void polar_log_mapping::calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {
    const FLT r = exp(loc(1))*(r0+r_eps) -r_eps;
    const FLT theta = -loc(0)/theta_length;
    
    const FLT drdlogr = exp(loc(1))*(r0+r_eps);
    
    /* Derivatives with respect to theta*length */
    jacobian(0,0) = +r*sin(theta)/theta_length;  // dx/dt
    jacobian(1,0) = -r*cos(theta)/theta_length; // dy/dt
    /* Derivaties with respect to logr */
    jacobian(0,1) = cos(theta)*drdlogr; // dx/dlogr
    jacobian(1,1) = sin(theta)*drdlogr; // dy/dlogr
}



void mapped_mesh::init(input_map& input, void *gbl_in) {
    r_tri_mesh::init(input,gbl_in);
    map->init(input,gbl->idprefix,gbl->log);
    mapped_pnts.resize(maxpst);
    map_pnts();
}

void mapped_mesh::init(const multigrid_interface& mgin, init_purpose why, FLT sizereduce1d) {
    r_tri_mesh::init(mgin,why,sizereduce1d);
    const mapped_mesh& smm = dynamic_cast<const mapped_mesh&>(mgin);
    map = smm.map;
    mapped_pnts.resize(maxpst);
    map_pnts();
}

void mapped_mesh::copy(const mapped_mesh& tgt) {
    r_tri_mesh::copy(tgt);
    map = tgt.map;
}
/** Outputs solution in various filetypes */
void mapped_mesh::output(const std::string &outname,block::output_purpose why) {
    if (why == block::display) {
        map_pnts();
        Array<TinyVector<FLT,ND>,1> temp; /**< Physical location of the points in the mesh */
        temp.reference(pnts);
        pnts.reference(mapped_pnts);
        tri_mesh::output(outname,output_type);
        pnts.reference(temp);
    }
    else {
        tri_mesh::output(outname,output_type);
    }
}

void mapped_mesh::map_pnts() {
    for (int i = 0; i < npnt; ++i) {
        map->to_physical_frame(pnts(i), mapped_pnts(i));
    }
}

