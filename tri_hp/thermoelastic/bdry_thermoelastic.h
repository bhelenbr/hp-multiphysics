/*
 *  bdry_thermoelastic.h
 *  tri_hp
 *
 *  Boundary conditions for thermoelastic plane strain.  The plane
 *  strain conditions (bdry_ps) are reused: their continuity flux goes
 *  to the last variable (p), and traction (_flux0, _flux1) and heat
 *  flux (_flux2) are natural conditions (zero unless set).
 *
 */

#ifndef _bdry_thermoelastic_h_
#define _bdry_thermoelastic_h_

#include "tri_hp_thermoelastic.h"
#include "../planestrain/bdry_ps.h"

namespace bdry_thermoelastic {

	/* DISPLACEMENT AND TEMPERATURE SET */
	class dirichlet : public bdry_ps::neumann {
		public:
			dirichlet(tri_hp_thermoelastic &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "dirichlet";}
			dirichlet(const dirichlet& inbdry, tri_hp_thermoelastic &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			dirichlet* create(tri_hp& xin, edge_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_thermoelastic&>(xin),bin);}
			void init(input_map& inmap) {
				neumann::init(inmap);
				for(int n=0;n<x.NV-1;++n) {
					if (type[n] != essential) essential_indices.push_back(n);
					type[n] = essential;
				}
			}
	};

	/* TEMPERATURE SET, TRACTION SPECIFIED */
	class temperature : public bdry_ps::neumann {
		public:
			temperature(tri_hp_thermoelastic &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "temperature";}
			temperature(const temperature& inbdry, tri_hp_thermoelastic &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			temperature* create(tri_hp& xin, edge_bdry &bin) const {return new temperature(*this,dynamic_cast<tri_hp_thermoelastic&>(xin),bin);}
			void init(input_map& inmap) {
				neumann::init(inmap);
				if (type[2] != essential) essential_indices.push_back(2);
				type[2] = essential;
			}
	};

	/* UPWINDED CONVECTIVE HEAT FLUX (SAME AS cd characteristic), TRACTION SPECIFIED */
	class characteristic : public bdry_ps::neumann {
		protected:
			tri_hp_thermoelastic &x;
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) override {
				neumann::flux(u,xpt,mv,norm,side_length,flx);

				FLT vel = (x.hp_te_gbl->ax -mv(0))*norm(0) +(x.hp_te_gbl->ay -mv(1))*norm(1);
				if (vel > 0.0)
					flx(2) += x.hp_te_gbl->rhocv*vel*u(2);
				else
					flx(2) += x.hp_te_gbl->rhocv*vel*ibc->f(2, xpt, x.gbl->time);
			}
		public:
			characteristic(tri_hp_thermoelastic &xin, edge_bdry &bin) : neumann(xin,bin), x(xin) {mytype = "characteristic";}
			characteristic(const characteristic& inbdry, tri_hp_thermoelastic &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), x(xin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_thermoelastic&>(xin),bin);}
	};
}

#endif
