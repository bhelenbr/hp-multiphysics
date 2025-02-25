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
#include "mappings.h"

/* AN R-DEFORMABLE MAPPED MULTI-GRID MESH OBJECT */
class mapped_mesh : public r_tri_mesh {
public:
    shared_ptr<mapping> map;
    Array<TinyVector<FLT,ND>,1> mapped_pnts; /**< Physical location of the points in the mesh */

    void init(input_map& input, shared_ptr<block_global>) override;
    /** Routine to initialze from another mesh with option of increasing or decreasing storage (compatible with block.h) */
    void init(const multigrid_interface& mgin, init_purpose why=duplicate, FLT sizereduce1d=1.0) override;
    /** Routine to copy */
    void copy(const mapped_mesh& tgt);
    /** Outputs solution in various filetypes */
    void output(const std::string &outname,block::output_purpose why) override;
    /** Update mapped locations after mesh adaptation */
    void adapt() override {r_tri_mesh::adapt(); map_pnts();}
    /** update physical location of points */
    void map_pnts();
};
#endif /* mapped_hpp */
