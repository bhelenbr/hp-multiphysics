/*
 *  bdry_ps.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ps.h"
#include "../hp_boundary.h"
#include <myblas.h>

namespace bdry_ps {

	class neumann : public hp_edge_bdry {
		protected:
			tri_hp_ps &x;
			virtual void flux(TinyVector<FLT,3> u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> norm, TinyVector<FLT,3>& flx) {
				flx(2) = u(0)*norm(0) +u(1)*norm(1);
				flx(0) = 0.0;
				flx(1) = 0.0;

				return;
			}

		public:
			neumann(tri_hp_ps &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "neumann";}
			neumann(const neumann& inbdry, tri_hp_ps &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin) {}
			neumann* create(tri_hp& xin, edge_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_ps&>(xin),bin);}
	};



	class dirichlet : public neumann {        
		public:
			dirichlet(tri_hp_ps &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "dirichlet";}
			dirichlet(const dirichlet& inbdry, tri_hp_ps &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			dirichlet* create(tri_hp& xin, edge_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_ps&>(xin),bin);}
			void init(input_map& inmap) {
				neumann::init(inmap);
				for(int n=0;n<x.ND;++n) {
					essential_indices.push_back(n);
					type[n] = essential;
				}
			}
	};

	class friction_wall : public hp_edge_bdry {
		protected:
			tri_hp_ps &x;
			int dir;
			FLT muwall;

			virtual void flux(TinyVector<FLT,3> u, TinyVector<FLT,2> stress, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> norm, TinyVector<FLT,3>& flx) {
				flx(2) = u(0)*norm(0) +u(1)*norm(1);
				flx(dir) = 0.0;  // not used;
				if (u(1-dir) > 1.0e-4) 
					flx(1-dir) = -muwall*fabs(stress(dir));
				else if (u(1-dir) < -1.0e-4)
					flx(1-dir) =  muwall*fabs(stress(dir));
				else
					flx(1-dir) = 0.0;

				flx(1-dir) = -muwall*fabs(stress(dir));

				return;
			}

		public:
			friction_wall(tri_hp_ps &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "friction_wall";}
			friction_wall(const friction_wall& inbdry, tri_hp_ps &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin), dir(inbdry.dir), muwall(inbdry.muwall) {}
			friction_wall* create(tri_hp& xin, edge_bdry &bin) const {return new friction_wall(*this,dynamic_cast<tri_hp_ps&>(xin),bin);}
			void init(input_map& inmap, void* gbl_in) {
				std::string keyword;

				hp_edge_bdry::init(inmap);

				keyword = base.idprefix + "_friction";
				inmap.getwdefault(keyword,muwall,0.2);

				keyword = base.idprefix + "_dir";
				inmap.getwdefault(keyword,dir,0);

				return;
			}

			void rsdl(int stage);

			void vdirichlet() {
				int v0,sind=-1;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					x.hp_gbl->res.v(v0,dir) = 0.0;
				}
				v0 = x.seg(sind).pnt(1);
				x.hp_gbl->res.v(v0,dir) = 0.0;
			}

			void sdirichlet(int mode) {
				int sind;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					x.hp_gbl->res.s(sind,mode,dir) = 0.0;
				}
			}

	};

    class curve_edges : public dirichlet {
        public:
            curve_edges(tri_hp_ps &xin, edge_bdry &bin) : dirichlet(xin,bin) {mytype = "curve_edges";}
            curve_edges(const curve_edges& inbdry, tri_hp_ps &xin, edge_bdry &bin) : dirichlet(inbdry,xin,bin) {}
            curve_edges* create(tri_hp& xin, edge_bdry &bin) const {return new curve_edges(*this,dynamic_cast<tri_hp_ps&>(xin),bin);}
            void init(input_map& inmap, void* gbl_in) {
                inmap[base.idprefix+"_curved"] = "0";
                dirichlet::init(inmap);
                return;
            }
            void setvalues(init_bdry_cndtn *ibc, const std::vector<int>& indices);
    };
}
