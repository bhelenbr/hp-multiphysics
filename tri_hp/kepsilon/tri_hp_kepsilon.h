/*
 *  tri_hp_kepsilon.h
 *  planar++
 *
 *  Standard (Fluent) k-epsilon turbulence model.
 *  Derived by converting the komega module (see ../komega).
 *
 */

#ifndef _tri_hp_kepsilon_h_
#define _tri_hp_kepsilon_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>
#include <symbolic_function.h>

class tri_hp_kepsilon : public tri_hp_ins {
public:
    struct hp_kepsilon_global {
        /* Physical Constants */
        FLT linf,uinf;

        /* Free-stream / reference turbulence values (used for ibc, MMS) */
        FLT kinf, epsinf;

        /* Floors used to protect divisions by k and epsilon */
        FLT kmin, epsmin;

        /* Standard k-epsilon model constants (Fluent defaults) */
        FLT Cmu, C1eps, C2eps, sigmak, sigmaeps;

        /* Whether the 2/3 rho k term is included in the momentum equation */
        int kmom_on;

        /* SOURCE FUNCTION FOR MMS */
        init_bdry_cndtn *src;
    };
    shared_ptr<hp_kepsilon_global> hp_kepsilon_gbl;


public:
    tri_hp_kepsilon* create() { return new tri_hp_kepsilon(); }

    void init(input_map& inmap,shared_ptr<block_global> gin);
    void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

    void error_estimator();
    void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
};
#endif
