/*
 *  pod_generate.h
 *  tet_hp
 *
 *  Created by Brian Helenbrook on 1/18/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */

#ifndef _POD_GENERATE_H_
#define _POD_GENERATE_H_

//#define LOWNOISE

#ifdef POD_BDRY
template<class BASE> class pod_gen_face_bdry;
template<class BASE> class pod_gen_edge_bdry;
template<class BASE> class pod_gen_vrtx_bdry;
#endif

template<class BASE> class pod_generate : public BASE {
	public:
		int nsnapshots;
		int restart_interval;
		int nmodes;
		int restartfile;
		Array<FLT,1> scaling;
#ifdef LOWNOISE
		int pod_id;
#else
		typedef typename BASE::vefi vefi;
		Array<vefi,1> modes;
		Array<FLT,1> psimatrix,psimatrix_recv;
#endif
#ifdef POD_BDRY
		Array<pod_gen_face_bdry<BASE> *, 1> pod_fbdry;
		Array<pod_gen_edge_bdry<BASE> *, 1> pod_ebdry;
		Array<pod_gen_vrtx_bdry<BASE> *, 1> pod_vbdry;
#endif

	public:
		void init(input_map& inmap, void *gin); 
		pod_generate<BASE>* create() { return new pod_generate<BASE>();}
		void tadvance();
};

#ifdef POD_BDRY
template<class BASE> class pod_gen_face_bdry {
	protected:
		pod_generate<BASE> &x;
		face_bdry &base;
		int nmodes;
		int pod_id;
		bool active;

	public:
		pod_gen_face_bdry(pod_generate<BASE>& xin, face_bdry &bin) : x(xin), base(bin) {}
		void init(input_map& inmap);
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
	void init(input_map& inmap);
	void zero_bdry(tet_hp::vefi ug);
	void calculate_modes();
	void output();
};
#endif

#include "pod_generate.cpp"

#endif
