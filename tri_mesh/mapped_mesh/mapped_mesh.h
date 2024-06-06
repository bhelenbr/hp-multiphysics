//
//  mapped.hpp
//  tri_mesh
//
//  Created by Brian Helenbrook on 6/6/24.
//

#ifndef _mapped_mesh_h
#define _mapped_mesh_h

#include <stdio.h>
#include "r_tri_mesh.h"
#include "block.h"
#include <spline.h>
#include "boundary.h"

/* AN R-DEFORMABLE MAPPED MULTI-GRID MESH OBJECT */

/* Generic interface to allow mapping between coordinate systems */
class mapping {
public:
    virtual void init(input_map& inmap,std::string idprefix,std::ostream *log) {}
    virtual void to_geometry_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) {}
    virtual void to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to) {}
    virtual ~mapping() {}
};

class spline_mapping : public mapping {
    spline<tri_mesh::ND> my_spline;
    FLT scale;
    rigid_movement_interface2D trsfm;
    void init(input_map& inmap,std::string idprefix,std::ostream *log);
    //void to_geometry_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to);
    void to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to);
};

class polar_mapping : public mapping {
    TinyVector<FLT,2> pnt;
    FLT theta_length;
    void init(input_map& inmap,std::string idprefix,std::ostream *log);
    //void to_geometry_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to);
    void to_physical_frame(const TinyVector<FLT,2>& from, TinyVector<FLT,2>& to);
};

class mapped_mesh : public r_tri_mesh {
public:
    mapping *map;
    Array<TinyVector<FLT,ND>,1> mapped_pnts; /**< Physical location of the points in the mesh */

    void init(input_map& input, void *gbl_in);
    /** Routine to initialze from another mesh with option of increasing or decreasing storage (compatible with block.h) */
    void init(const multigrid_interface& mgin, init_purpose why=duplicate, FLT sizereduce1d=1.0);
    /** Routine to copy */
    void copy(const mapped_mesh& tgt);
    /** Outputs solution in various filetypes */
    void output(const std::string &outname,block::output_purpose why);
    /** update physical location of points */
    void map_pnts();
    ~mapped_mesh() {delete map;}
};
#endif /* mapped_hpp */
