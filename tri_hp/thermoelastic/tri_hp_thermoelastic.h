/*
 *  tri_hp_thermoelastic.h
 *  tri_hp
 *
 *  Plane strain elasticity driven by thermal strain, coupled to a
 *  convection-diffusion equation for the temperature.
 *
 *  Variables are (u, v, T, p): x/y displacement, temperature, and
 *  p = -lambda*div(u) +(3*lambda+2*mu)*alpha*(T-T_ref), so that the
 *  in-plane stress is sigma = -p I +mu*(grad(u) +grad(u)^T).
 *
 */

#ifndef _tri_hp_thermoelastic_h_
#define _tri_hp_thermoelastic_h_

#include "../planestrain/tri_hp_ps.h"
#include <blocks.h>

class tri_hp_thermoelastic : public tri_hp_ps {
	public:
		/* THINGS SHARED BY ALL tri_hp_thermoelastic in same multigrid block */
		struct hp_thermoelastic_global {
			/* STABILIZATION FOR TEMPERATURE EQUATION */
			Array<FLT,1> tau;

			/* PHYSICAL CONSTANTS */
			FLT rhocv,kcond,ax,ay;
			FLT alpha;  // thermal expansion coefficient
			FLT T_ref;  // temperature at which there is no thermal strain
			FLT E;      // Young's modulus (for scaling the adaptation measure)

			/* SOURCE FUNCTION */
			init_bdry_cndtn *src;
		};
		shared_ptr<hp_thermoelastic_global> hp_te_gbl;

		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);
		init_bdry_cndtn* getnewibc(std::string name);

	public:
		tri_hp_thermoelastic* create() { return new tri_hp_thermoelastic(); }

		void init(input_map& inmap,shared_ptr<block_global> gin);
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void length();
		void error_estimator();
		int setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		void calculate_unsteady_sources();
};
#endif
