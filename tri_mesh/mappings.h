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
    virtual void to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) {}
    virtual void to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) {}
    virtual void calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) {}
    virtual ~mapping() {}
};

class spline_mapping : public mapping {
protected:
    spline<2> my_spline;
    FLT scale;
public:
    rigid_movement_interface2D trsfm;
    void init(input_map& inmap,std::string idprefix,std::ostream *log) override;
    void to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    void to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    void calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override;
};

class polar_mapping : public mapping {
protected:
    TinyVector<FLT,2> pnt;
    FLT theta_length;
    void init(input_map& inmap,std::string idprefix,std::ostream *log) override;
    void to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    void to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    void calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override;
};

class polar_log_mapping : public polar_mapping {
    FLT r0, r_eps;
    void init(input_map& inmap,std::string idprefix,std::ostream *log) override;
    void to_parametric_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    void to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) override;
    void calc_metrics(const TinyVector<FLT,2> loc, TinyMatrix<FLT,2,2>& jacobian) override;
};
#endif
