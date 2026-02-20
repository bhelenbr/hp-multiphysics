//
//  mappings.h
//  tri_mesh
//
//  Created by Brian Helenbrook on 2/20/25.
//
#include "block.h"
#include <spline.h>
#include "boundary.h"

#ifndef _mappings_h
#define _mappings_h

/* Generic interface to allow mapping between coordinate systems */
class mapping {
public:
    virtual void init(input_map& inmap,std::string idprefix,std::ostream *log) {}
    virtual int to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) {return(0);}
    virtual int to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) {return(0);}
    virtual int calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {return(0);}
    virtual ~mapping() {}
};

class no_mapping : public mapping {
public:
    virtual int to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override {to = from; return(0);}
    virtual int to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override {to = from; return(0);}
    virtual int calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override {
        jacobian(0,0) = 1.0;
        jacobian(1,0) = 0.0;
        jacobian(0,1) = 0.0;
        jacobian(1,1) = 1.0;
        return(0);
    }
    virtual ~no_mapping() {}
};

class spline_mapping : public mapping {
protected:
    SPLINE<2> my_spline;
    FLT scale;
public:
    rigid_movement_interface2D trsfm;
    void init(input_map& inmap,std::string idprefix,std::ostream *log) override;
    int to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override;
};

class spline_log_mapping : public spline_mapping {
protected:
    FLT r0, r_eps;
public:
    void init(input_map& inmap,std::string idprefix,std::ostream *log) override;
    int to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override;
};

class polar_mapping : public mapping {
protected:
    TinyVector<FLT,2> pnt;
    FLT theta_length, theta0;
    void init(input_map& inmap,std::string idprefix,std::ostream *log) override;
    int to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override;
};

class polar_log_mapping : public polar_mapping {
    FLT r0, r_eps;
    void init(input_map& inmap,std::string idprefix,std::ostream *log) override;
    int to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override;
};

class polar_chi_mapping : public polar_mapping {
    FLT r0, m0, eta_s;
    void init(input_map& inmap,std::string idprefix,std::ostream *log) override;
    int to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    int calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override;
};

shared_ptr<mapping> getnewmapping(input_map& inmap, std::string mapname);
#endif
