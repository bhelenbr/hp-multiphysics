/*
 *  pod.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/26/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#define LOWNOISE

template<class BASE> class pod_gen_edge_bdry;
template<class BASE> class pod_gen_vrtx_bdry;

#ifdef LOWNOISE
template<class BASE> class pod_generate : public BASE {
    protected:
        int nsnapshots;
        int nmodes;
        Array<FLT,1> scaling;
        Array<FLT,1> coeffs;
		Array<pod_gen_edge_bdry<BASE> *, 1> pod_ebdry;
		Array<pod_gen_vrtx_bdry<BASE> *, 1> pod_vbdry;
	
    public:
        void init(input_map& input, void *gin); 
        pod_generate<BASE>* create() { return new pod_generate<BASE>();}
        void tadvance();
};
#else
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
#endif

template<class BASE> class pod_gen_edge_bdry {
	protected:
		pod_generate<BASE> &x;
		edge_bdry &base;
		int nmodes;
		bool active;
	
	public:
		pod_gen_edge_bdry(pod_generate<BASE>& xin, edge_bdry &bin) : x(xin), base(bin) {}
		void init(input_map& input);
		void zero_bdry();
		void calculate_modes();
		void output();
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
