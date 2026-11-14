/*
 *  melt_thermoelastic.h
 *  tri_hp
 *
 *  Solid side of a melt (solidification) interface for the thermoelastic
 *  module.  The temperature coupling, interface motion and kinetics are
 *  handled exactly as for a cd solid (bdry_cd::melt_cd, with the liquid
 *  as master).  The displacement equations get the same natural conditions
 *  as a plane strain neumann boundary: traction from _flux0/_flux1 (zero
 *  unless set) and the u.n continuity flux for p.
 *
 */

#ifndef _melt_thermoelastic_h_
#define _melt_thermoelastic_h_

#include "tri_hp_thermoelastic.h"
#include "../cd/melt_cd.h"

namespace bdry_thermoelastic {

	class melt : public bdry_cd::melt_cd {
		protected:
			tri_hp_thermoelastic &x;
			bool elastic_only; // When true element_rsdl gives only the displacement & continuity fluxes

			void solid_velocity(TinyVector<FLT,tri_mesh::ND>& us) override {
				us(0) = x.hp_te_gbl->ax;
				us(1) = x.hp_te_gbl->ay;
			}

		public:
			melt(tri_hp_thermoelastic &xin, edge_bdry &bin) : melt_cd(xin,bin), x(xin), elastic_only(false) {mytype = "melt"; solid_side = true;}
			melt(const melt& inbdry, tri_hp_thermoelastic &xin, edge_bdry &bin) : melt_cd(inbdry,xin,bin), x(xin), elastic_only(false) {}
			melt* create(tri_hp& xin, edge_bdry &bin) const override {return new melt(*this,dynamic_cast<tri_hp_thermoelastic&>(xin),bin);}

			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) override {
				/* TRACTIONS & HEAT FLUX (ZERO UNLESS SET WITH _flux0/_flux1/_flux2) */
				hp_edge_bdry::flux(u,xpt,mv,norm,side_length,flx);
				/* CONTINUITY (PRESSURE IS THE LAST VARIABLE) */
				flx(x.NV-1) = u(0)*norm(0) +u(1)*norm(1);
			}

			void element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) override {
				if (elastic_only) {
					hp_edge_bdry::element_rsdl(indx,lf);
					lf(2) = 0.0;  // Temperature flux is handled by the coupled (liquid) side
				}
				else {
					melt_cd::element_rsdl(indx,lf);
				}
			}

			void rsdl(int stage) override {
				melt_cd::rsdl(stage);
				elastic_only = true;
				hp_edge_bdry::rsdl(stage);
				elastic_only = false;
			}

#ifdef petsc
			void petsc_jacobian() override;
#endif
	};
}
#endif
