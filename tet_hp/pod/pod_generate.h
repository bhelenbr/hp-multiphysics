/*
 *  pod_generate.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 1/18/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#ifndef _POD_GENERATE_H_
#define _POD_GENERATE_H_

//#define LOWNOISE

template<class BASE> class pod_gen_face_bdry;
template<class BASE> class pod_gen_edge_bdry;
template<class BASE> class pod_gen_vrtx_bdry;

#ifdef LOWNOISE
template<class BASE> class pod_generate : public BASE {
	public:
		int nsnapshots;
		int nmodes;
		int pod_id;
		Array<FLT,1> scaling;
		Array<FLT,1> coeffs;
		Array<pod_gen_face_bdry<BASE> *, 1> pod_fbdry;
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
		typedef typename BASE::vefi vefi;
		Array<vsi,1> modes;
		Array<FLT,1> psimatrix,psimatrix_recv;

	public:
		void init(input_map& input, void *gin); 
		pod_generate<BASE>* create() { return new pod_generate<BASE>();}
		void tadvance();
};
#endif

template<class BASE> class pod_gen_face_bdry {
	protected:
		pod_generate<BASE> &x;
		face_bdry &base;
		int nmodes;
		int pod_id;
		bool active;

	public:
		pod_gen_face_bdry(pod_generate<BASE>& xin, face_bdry &bin) : x(xin), base(bin) {}
		void init(input_map& input);
		void zero_bdry(tet_hp::vefi ug);
		void calculate_modes();
		void output();
};

template<class BASE> class pod_gen_edge_bdry {
protected:
	pod_generate<BASE> &x;
	edge_bdry &base;
	int nmodes;
	int pod_id;
	bool active;
	
public:
	pod_gen_edge_bdry(pod_generate<BASE>& xin, edge_bdry &bin) : x(xin), base(bin) {}
	void init(input_map& input);
	void zero_bdry(tet_hp::vefi ug);
	void calculate_modes();
	void output();
};

#include "pod_generate.cpp"

#endif