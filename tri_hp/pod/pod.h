/*
 *  pod.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/26/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

template<class BASE> class pod_generate : public BASE {
    protected:
        int nsnapshots;
        int nmodes;
        Array<FLT,1> scaling;
        Array<FLT,1> coeffs;
        typedef typename BASE::vsi vsi;
        Array<vsi,1> modes;
        Array<FLT,1> psimatrix,psimatrix_recv;
    public:
        void init(input_map& input, void *gin); 
        pod_generate<BASE>* create() { return new pod_generate<BASE>();}
        void tadvance();
};

template<class BASE> class pod_simulate : public BASE {
    protected:
        int nmodes;
        Array<FLT,1> coeffs, rsdls, rsdls_recv;
        typedef typename BASE::vsi vsi;
        Array<vsi,1> modes;
        Array<FLT,2> jacobian;
        Array<int,1> ipiv;
        
    public:
        void init(input_map& input, void *gin); 
        pod_simulate<BASE>* create() { return new pod_simulate<BASE>();}
        void rsdl(int stage);
        void setup_preconditioner();
        void update();
        FLT maxres();        
};




#include "pod.cpp"
