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
		Array<FLT,1> scaling;
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
		tri_hp_helper* getnewhelper(input_map& inmap);
		void output(const std::string& fname, block::output_purpose why);
		void calc_coeffs();
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

#ifdef POD_BDRY
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
#endif

class svv_ins : public pod_simulate<tri_hp_ins> {
	public:
		int cutoff; // Mode number to turn on additional viscosity
		FLT alpha; // Additional viscosity beyond cut-off mode
		bool remove_pressure; // Flag to remove pressure from equations based on incompressibility
	
		svv_ins* create() { return new svv_ins();}
		void init(input_map& input, void *gin) { 
			pod_simulate<tri_hp_ins>::init(input,gin);
			if (!input.get("svv_cutoff", cutoff)) {
				*gbl->log << "Failed to find svv cutoff number" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			input.getwdefault("svv_alpha", alpha,0.0);	
			input.getwdefault("svv_remove_pressure",remove_pressure,false);
			
			if (remove_pressure) {
				for(int m=0;m<tmodes;++m) {
					modes(m).v(Range::all(),NV-1) = 0.0;
					modes(m).s(Range::all(),Range::all(),NV-1) = 0.0;
					modes(m).i(Range::all(),Range::all(),NV-1) = 0.0;
				}
			}
		}
	
	
		void rsdl(int stage){
			
			// calculate residaul with normal viscosity
			pod_simulate<tri_hp_ins>::rsdl(stage);
			// Store this POD resdiual in rsdls0
			rsdls0 = rsdls_recv;

			// Add spetral vanishing viscosity
			gbl->mu += alpha;
			// Recalculate residuals 
			pod_simulate<tri_hp_ins>::rsdl(stage);
			
			// combine previously stored rsdls_recv with this new one.
			rsdls_recv(Range(0,cutoff-1)) = rsdls0(Range(0,cutoff-1));
			
			// Remove spectral vansishing viscosity
			gbl->mu -= alpha;
		}
};

#include "pod_simulate.cpp"

#endif
