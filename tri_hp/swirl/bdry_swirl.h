/*
 *  swirl_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swirl.h"
#include "../ins/bdry_ins.h"
#include "../hp_boundary.h"

namespace bdry_swirl {

class symmetry : public bdry_ins::generic {
    void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
        /* THESE DON'T GET USED */
        flx(Range(0,x.NV-1)) = 0.0;
        return;
    }
    
public:
    /* Fixme: Replace init with a constructor that accepts an input_map? */
    symmetry(tri_hp_swirl &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "symmetry";}
    symmetry(const symmetry& inbdry, tri_hp_swirl &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
    symmetry* create(tri_hp& xin, edge_bdry &bin) const {return new symmetry(*this,dynamic_cast<tri_hp_swirl&>(xin),bin);}
    void init(input_map& inmap) {
        generic::init(inmap);
        type[0] = essential;
        essential_indices.push_back(0);
        type[2] = essential;
        essential_indices.push_back(2);
    }
    void tadvance();
};

class actuator_disc: public bdry_ins::actuator_disc {
protected:
    symbolic_function<3> dq; // differential torque
    
    void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
        /* CONTINUITY */
        flx(x.NV-1) = x.hp_ins_gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
        
        if (base.is_frst()) {
            FLT norm_vel = flx(x.NV-1)/x.hp_ins_gbl->rho;
            TinyVector<FLT,3> inpt(xpt(0),xpt(1),norm_vel);
            FLT delta_p = dp.Eval(inpt,x.gbl->time);
            FLT delta_q = dq.Eval(inpt,x.gbl->time);
            
            /* X&Y MOMENTUM */
            for (int n=0;n<tri_mesh::ND;++n)
                flx(n) = norm(n)*delta_p;
            
            flx(tri_mesh::ND) = norm(1)*delta_q;
        }
        else {
            for (int n=0;n<x.NV-1;++n)
                flx(n) = 0.0;
        }
        
        /* EVERYTHING ELSE */
        for (int n=tri_mesh::ND+1;n<x.NV-1;++n)
            flx(n) = flx(x.NV-1)*u(n);
        
        return;
    }
    
public:
    actuator_disc(tri_hp_ins &xin, edge_bdry &bin) : bdry_ins::actuator_disc(xin,bin) {}
    actuator_disc(const actuator_disc& inbdry, tri_hp_ins &xin, edge_bdry &bin) : bdry_ins::actuator_disc(inbdry,xin,bin), dq(inbdry.dq) {}
    actuator_disc* create(tri_hp& xin, edge_bdry &bin) const {return new actuator_disc(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
    
    /* To read in data */
    void init(input_map& inmap) {
        
        bdry_ins::actuator_disc::init(inmap);
        
        /* LOAD PRESSURE JUMP FUNCTION */
        dq.init(inmap,base.idprefix+"_shear");
    }
    // void output(const std::string& fname, tri_hp::filetype typ,int tlvl = 0);
};
}
