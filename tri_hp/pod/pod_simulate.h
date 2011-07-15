/*
 *  pod.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/26/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _POD_SIMULATE_H_
#define _POD_SIMULATE_H_


#ifdef POD_BDRY
template<class BASE> class pod_sim_edge_bdry;
template<class BASE> class pod_sim_vrtx_bdry;
#endif

template<class BASE> class pod_simulate : public BASE {
	protected:
		int pod_id;
		int nmodes;
		int tmodes;
		Array<FLT,1> coeffs, rsdls, rsdls_recv, rsdls0;
		typedef typename BASE::vsi vsi;
		Array<vsi,1> modes;
#ifdef POD_BDRY
		Array<FLT,1> multiplicity;
		Array<pod_sim_edge_bdry<BASE> *, 1> pod_ebdry;
		Array<pod_sim_vrtx_bdry<BASE> *, 1> pod_vbdry;
#endif
		Array<FLT,2> jacobian;
		Array<int,1> ipiv;
#ifdef POD_BDRY
		friend class pod_sim_edge_bdry<BASE>;
#endif

	public:
		void init(input_map& input, void *gin); 
		pod_simulate<BASE>* create() { return new pod_simulate<BASE>();}
		void output(const std::string& fname, block::output_purpose why);
		void tadvance();
		void rsdl(int stage);
		void setup_preconditioner();
		void update();
		FLT maxres();

#ifdef POD_BDRY
		/* communication for boundary modes */
		void sc0load();
		int sc0wait_rcv();
#endif
};

#ifdef POD_BRY
template<class BASE> class pod_sim_edge_bdry {
	protected:
		pod_simulate<BASE> &x;
		edge_bdry &base;
		int nmodes;
		int pod_id;
		int bindex; 
		bool active;

		struct vs {
			Array<FLT,2> v;
			Array<FLT,3> s;
		} ug;
		Array<vs,1> modes;
		friend class pod_simulate<BASE>;

	public:
		pod_sim_edge_bdry(pod_simulate<BASE>& xin, edge_bdry &bin) : x(xin), base(bin) {}
		void init(input_map& input);
		void rsdl();
		void addto2Dsolution(struct tri_hp::vsi ug);
		void addto2Dsolution(struct tri_hp::vsi ug, int mode, FLT coeff);
		void update();
		void loadbuff(Array<FLT,1>& sdata);
		void finalrcv(Array<FLT,1>& sdata);
};

//template<class BASE> class pod_coefficients : public pod_simulate<BASE> {
//	public:
//		pod_coefficents<BASE>* create() { return new pod_coefficients<BASE>();}
//		void tadvance();
//}

#endif

#include "pod_simulate.cpp"

#endif
