//
//  cd_mapped.hpp
//  tri_hp
//
//  Created by Brian Helenbrook on 6/8/24.
//

#ifndef _cd_mapped_h_
#define _cd_mapped_h_

#include <stdio.h>
#include "tri_hp_cd.h"
#include <mapped_mesh.h>

class cd_mapped : public tri_hp_cd {
public:
    shared_ptr<mapping> map;
    cd_mapped* create() override { return new cd_mapped(); }
    void init(input_map& inmap, void *gin) override;
    void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0) override;
    void calc_metrics(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND>& crd, TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,ND,ND>& dcrd, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND>& mvel) const;
    void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im) override;
    void l2error(init_bdry_cndtn *comparison) override;
    void output(const std::string& fname, block::output_purpose why) override;
};
#endif /* _cd_mapped_h_ */
