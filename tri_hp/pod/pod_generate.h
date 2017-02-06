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

//#define DIRECT_METHOD // This creates modes using the two-point correlation tensor not the snapshot method
#define USING_MASS_MATRIX // This uses the mass matrix to form inner products rather than numerically evaluate them

#define POD_BDRY

#if (defined(DIRECT_METHOD)  && !defined(USING_MASS_MATRIX))
error: DIRECT_METHOD only works with USING_MASS_MATRIX defined.
#endif

#ifdef POD_BDRY
template<class BASE> class pod_gen_edge_bdry;
template<class BASE> class pod_gen_vrtx_bdry;
#endif

template<class BASE> class pod_generate : public BASE {
	public:
		int nsnapshots; // number of snapshots used to make modes
		int nmodes;  // nmodes to output
		int restartfile; // number of first snapshot
		int restart_interval; // interval between numbering of snapshots
		int coefficient_start,coefficient_interval,coefficient_end; // restart numbering to output coefficient files
		int pod_id;  // For problems with separate PODs on groups of blocks
		int ndeflation, nmodes_per_deflation;  // To use the deflation method
		int nsets, nsnapshots_per_set; // For combining modes each with their own vector of eigenvalue weights
		
		Array<double,1> set_eigenvalues;
		Array<FLT,1> scaling;
    
#ifdef USING_MASS_MATRIX
		sparse_row_major mass;
		Array<FLT,1> mass_times_snapshot;
#endif
		typedef typename BASE::vsi vsi;
		Array<vsi,1> modes;
		
#ifdef POD_BDRY
		Array<pod_gen_edge_bdry<BASE> *, 1> pod_ebdry;
		Array<pod_gen_vrtx_bdry<BASE> *, 1> pod_vbdry;
#endif

	public:
		void init(input_map& inmap, void *gin); 
		pod_generate<BASE>* create() { return new pod_generate<BASE>();}
		void tadvance();
		void create_mass_matrix(sparse_row_major& mass);
		void project_to_gauss(vsi& target);
		void time_average(vsi& target, Array<FLT,2> correlation_matrix);
		FLT inner_product_with_projection(vsi& target);
		FLT norm2(vsi& target);
		void test_orthogonality();
		void output(vsi& target,std::string filename);
};
/*template<class BASE> class wdpod_generate : public pod_generate<BASE>{
	public:
		int nsets;
		int ntruncation;
		Array<double,1> set_eigenvalues ;
}*/

template<class BASE> class gram_schmidt : public pod_generate<BASE> {
	public:
		gram_schmidt<BASE>* create() { return new gram_schmidt<BASE>();}
		void tadvance();
};


template<class BASE> class pod_generate_with_r_mesh : public pod_generate<BASE> {
	public:
		void init(input_map& inmap, void *gin); 
		pod_generate<BASE>* create() { return new pod_generate<BASE>();}
		void tadvance();
};


#ifdef POD_BDRY
template<class BASE> class pod_gen_edge_bdry {
	protected:
		pod_generate<BASE> &x;
		edge_bdry &base;
		int nmodes;
		int pod_id;
		bool active;
		typedef typename BASE::vsi vsi;

	public:
		pod_gen_edge_bdry(pod_generate<BASE>& xin, edge_bdry &bin) : x(xin), base(bin) {}
		void init(input_map& inmap);
		void zero_bdry(tri_hp::vsi& ug);
		void calculate_modes();
		void output(vsi& target,std::string filename,typename BASE::filetype typ);
};
#endif

#include "pod_generate.cpp"

#endif
