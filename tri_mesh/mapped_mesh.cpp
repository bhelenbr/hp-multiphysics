//
//  mapped.cpp
//  tri_mesh
//
//  Created by Brian Helenbrook on 6/6/24.
//

#include "mapped_mesh.h"

shared_ptr<mapping> getnewmapping(input_map& inmap, std::string mapname) {
    
    if (mapname == "none") {
        return(make_shared<no_mapping>());
    }
    
    std::string maptype;
    if (inmap.get(mapname+"_type",maptype)) {
        if (maptype == "polar") {
            return(make_shared<polar_mapping>());
        }
        else if (maptype == "polar_log") {
            return(make_shared<polar_log_mapping>());
        }
        else if (maptype == "polar_chi") {
            return(make_shared<polar_chi_mapping>());
        }
        else if (maptype == "spline") {
            return(make_shared<spline_mapping>());
        }
        else if (maptype == "spline_log") {
            return(make_shared<spline_log_mapping>());
        }
    }
    std::cerr << "Unrecognized mapping " << mapname << ' ' << maptype << std::endl;
    exit(1);
}

void spline_mapping::init(input_map& input, std::string idprefix, std::ostream *log) {
    trsfm.init(input,idprefix);
    std::string line;
    if (!input.get(idprefix+"_spline",line)) {
        *log << "Couldn't fine spline file name in input file " << idprefix +"_spline" <<std::endl;
        sim::abort(__LINE__,__FILE__,log);
    }
    my_spline.read(line);
    input.getwdefault(idprefix+"_scale",scale,1.0);
}

int spline_mapping::to_physical_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    TinyVector<FLT,2> tan, curv;
    spline_functions2D::interpolate(to, tan, curv, my_spline, from(0), scale, trsfm.theta,trsfm.pos, -from(1));
    return(0);
}

int spline_mapping::to_parametric_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    TinyVector<FLT,2> tan, curv;
    to(1) *= -1;
    int err = spline_functions2D::find_with_guess(from, my_spline, to(0), scale, trsfm.theta,trsfm.pos, to(1));
    to(1) *= -1;
    if (to(1) < -FLT_EPSILON || err) {
        err = spline_functions2D::find(from, my_spline, to(0), scale, trsfm.theta,trsfm.pos, to(1));
        to(1) *= -1;
    }
    return(err);
}

int spline_mapping::calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {
    TinyVector<FLT,2> pnt, tan, curv;
    spline_functions2D::interpolate(pnt, tan, curv, my_spline, loc(0), scale, trsfm.theta,trsfm.pos, -loc(1));
    /* p = x(s) +n*norm_dist */
    /* p = x(loc(0)) +loc(1)*(-tan(1),tan(0))*/
    /* dp/ds = dx/ds +curv * norm_dist */
    /* dp/dn = norm */
    
    /* Derivatives with respect to s*/
    jacobian(0,0) = tan(0) -curv(1)*loc(1);
    jacobian(1,0) = tan(1) +curv(0)*loc(1);
    /* Derivaties with respect to norm_dist */
    jacobian(0,1) = -tan(1);
    jacobian(1,1) = +tan(0);
    return(0);
}


void spline_log_mapping::init(input_map& input, std::string idprefix, std::ostream *log) {
    spline_mapping::init(input,idprefix,log);
    
    if (!input.get(idprefix+"_r0",r0)) {
        *log << "Couldn't read r0 " << idprefix+"_r0" << std::endl;;
        sim::abort(__LINE__,__FILE__,log);
    }
    input.getwdefault(idprefix+"_r_eps",r_eps,DBL_EPSILON);
    r_eps = r_eps*r0;
}

int spline_log_mapping::to_physical_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    TinyVector<FLT,2> from2;
    from2(0) = from(0);
    from2(1) = exp(from(1))*(r0+r_eps) -r_eps;
    int err = spline_mapping::to_physical_frame(from2, to);
    return(err);
}

int spline_log_mapping::to_parametric_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    int err = spline_mapping::to_parametric_frame(from, to);
    to(1) = log((to(1)+r_eps)/(r0+r_eps));
    return(err);
}

int spline_log_mapping::calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {
    const FLT r = exp(loc(1))*(r0+r_eps) -r_eps;
    const FLT drdlogr = exp(loc(1))*(r0+r_eps);

    TinyVector<FLT,2> loc2;
    loc2(0) = loc(0);
    loc2(1) = r;
    int err = spline_mapping::calc_metrics(loc2, jacobian);
    jacobian(0,1) *= drdlogr;
    jacobian(1,1) *= drdlogr;
    return(err);
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
    
    input.getwdefault(idprefix+"_theta0", theta0, 0.0);
   
}

int polar_mapping::to_physical_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    const FLT r = from(1);
    const FLT theta = -from(0)/theta_length;
    to(0) = pnt(0) +r*cos(theta);
    to(1) = pnt(1) +r*sin(theta);
    return(0);
}

int polar_mapping::to_parametric_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    to = from-pnt;
    const FLT r = sqrt(to(0)*to(0) +to(1)*to(1));
    FLT alpha = atan2(to(1),to(0)) -theta0;
    alpha += (alpha < -M_PI ? 2.*M_PI : 0.0) +theta0;
    to(0) = -alpha*theta_length;
    to(1) = r;
    return(0);
}


int polar_mapping::calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {
    const FLT r = loc(1);
    const FLT theta = -loc(0)/theta_length;
    
    /* Derivatives with respect to theta*length */
    jacobian(0,0) = +r*sin(theta)/theta_length;  // dx/dt
    jacobian(1,0) = -r*cos(theta)/theta_length; // dy/dt
    /* Derivaties with respect to r */
    jacobian(0,1) = cos(theta); // dx/dr
    jacobian(1,1) = sin(theta); // dy/dr
    return(0);
}

void polar_log_mapping::init(input_map& input, std::string idprefix, std::ostream *log) {
    polar_mapping::init(input,idprefix,log);
    
    if (!input.get(idprefix+"_r0",r0)) {
        *log << "Couldn't read r0 " << idprefix+"_r0" << std::endl;;
        sim::abort(__LINE__,__FILE__,log);
    }
    input.getwdefault(idprefix+"_eps",r_eps,1.0e-8);
    r_eps *= r0;
}

int polar_log_mapping::to_physical_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    const FLT r = exp(from(1)/r0)*(r0+r_eps) -r_eps;
    const FLT theta = -from(0)/theta_length;
    to(0) = pnt(0) +r*cos(theta);
    to(1) = pnt(1) +r*sin(theta);
    return(0);
}

int polar_log_mapping::to_parametric_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    to = from-pnt;
    const FLT r = sqrt(to(0)*to(0) +to(1)*to(1));
    FLT alpha = atan2(to(1),to(0)) -theta0;
    alpha += (alpha < -M_PI ? 2.*M_PI : 0.0) +theta0;
    to(0) = -alpha*theta_length;
    to(1) = r0*log((r+r_eps)/(r0+r_eps));
    return(0);
}

int polar_log_mapping::calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {
    const FLT r = exp(loc(1)/r0)*(r0+r_eps) -r_eps;
    const FLT theta = -loc(0)/theta_length;
    
    const FLT drdlogr = exp(loc(1)/r0)*(r0+r_eps)/r0;
    
    /* Derivatives with respect to theta*length */
    jacobian(0,0) = +r*sin(theta)/theta_length;  // dx/dt
    jacobian(1,0) = -r*cos(theta)/theta_length; // dy/dt
    /* Derivaties with respect to logr */
    jacobian(0,1) = cos(theta)*drdlogr; // dx/dlogr
    jacobian(1,1) = sin(theta)*drdlogr; // dy/dlogr
    return(0);
}

void polar_chi_mapping::init(input_map& input, std::string idprefix, std::ostream *log) {
    polar_mapping::init(input,idprefix,log);
    if (!input.get(idprefix+"_r0",r0)) {
        *log << "Couldn't read r0 " << idprefix+"_r0" << std::endl;;
        sim::abort(__LINE__,__FILE__,log);
    }
    input.getwdefault(idprefix+"_m0",m0,5.0); // May be better to set it as p+1
    input.getwdefault(idprefix+"_eta_s",eta_s,0.95);

}

int polar_chi_mapping::to_physical_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {

    const FLT denom = eta_s*(1.0 -exp(1.0 +1.0/(pow(eta_s,2.0*m0-1.0) -1.0)));
    const FLT r = r0*from(1)*(1.0 -exp(1.0 +1.0/(pow(from(1),2.0*m0-1.0) -1.0)))/denom;
    const FLT theta = -from(0)/theta_length;
    to(0) = pnt(0) +r*cos(theta);
    to(1) = pnt(1) +r*sin(theta);
    return(0);
}

int polar_chi_mapping::to_parametric_frame(const TinyVector<double, 2> &from, TinyVector<double, 2> &to) {
    to = from-pnt;
    const FLT r = sqrt(to(0)*to(0) +to(1)*to(1));
    FLT alpha = atan2(to(1),to(0)) -theta0;
    alpha += (alpha < -M_PI ? 2.*M_PI : 0.0) +theta0;
    to(0) = -alpha*theta_length;
    
    // Newton-Raphson solver for eta
    const int max_iter = 50;
    const FLT tol = 1e-12;
    
    const FLT k = 2.0*m0 - 1.0;
    const FLT eta_ks = pow(eta_s,k) -1.0;
    const FLT Es = exp(1.0 +1.0/eta_ks);
    const FLT r_s = eta_s*(1.0 -Es);
    
    FLT eta = r / r0; // initial guess
    
    for(int i=0; i<max_iter; ++i)  {
        FLT denom = pow(eta, k) - 1.0;
        
        // Avoid division blow-up near eta=1
        if(abs(denom) < 1e-14)
            denom = (denom < 0 ? -1e-14 : 1e-14);
        
        FLT E = exp(1.0 + 1.0/denom);
        
        // F(eta)
        FLT F = r0*eta*(1.0 -E)/r_s -r;
        
        // derivative
        FLT dF = r0*((1.0 -E) +k*pow(eta,k)*E/pow(denom,2.0))/r_s;
        
        FLT delta = F/dF;
        eta -= F/dF;
        
        // Keep eta in valid range
        eta = min(max(eta, FLT(1e-8)), eta_s);
        
        if(abs(delta) < tol)
            break;}
    cout << "to_parametric_called!" << endl;
    to(1) = eta;
    return(0);
}

int polar_chi_mapping::calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {
    
    const FLT theta = -loc(0)/theta_length;
    const FLT k =2.0*m0 -1.0;
    const FLT denom = pow(loc(1),k) -1.0;
    const FLT E = exp(1.0 +1.0/denom);
    const FLT eta_ks = pow(eta_s,k) -1.0;
    const FLT Es = exp(1.0 +1.0/eta_ks);
    const FLT r_s = eta_s*(1.0 -Es);
    const FLT r = r0*loc(1)*(1.0 -E)/r_s;
    const FLT drdeta = r0*((1.0 -E) +k*pow(loc(1),k)*E/pow(denom,2.0))/r_s;
    
   // std::cout << drdeta << std::endl;
    
    /* Derivatives with respect to theta*length */
    jacobian(0,0) = +r*sin(theta)/theta_length;  // dx/dt
    jacobian(1,0) = -r*cos(theta)/theta_length; // dy/dt
    /* Derivaties with respect to eta */
    jacobian(0,1) = cos(theta)*drdeta; // dx/deta
    jacobian(1,1) = sin(theta)*drdeta; // dy/deta
    return(0);
}



void mapped_mesh::init(input_map& input, shared_ptr<block_global> gbl_in) {
    r_tri_mesh::init(input,gbl_in);

    std::string mapname;
    if (input.get(gbl->idprefix+"_mapping",mapname)) {
        map = getnewmapping(input,mapname);
    }
    else {
        /* No mapping */
        map = getnewmapping(input,"none");
    }
    map->init(input,mapname,gbl->log);
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

