/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include "../hp_boundary.h"
#include "myblas.h"

namespace bdry_cd {
    class dirichlet : public hp_edge_bdry {
		tri_hp_cd &x;

		public:
			dirichlet(tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "dirichlet";}
			dirichlet(const dirichlet& inbdry, tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin) {}
			dirichlet* create(tri_hp& xin, edge_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			void vdirichlet() {
				int sind,v0;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					v0 = x.seg(sind).pnt(0);
					x.gbl->res.v(v0,0) = 0.0;
				}
				v0 = x.seg(sind).pnt(1);
				x.gbl->res.v(v0,0) = 0.0;
			}

			void sdirichlet(int mode) {
				int sind;

				for(int j=0;j<base.nseg;++j) {
					sind = base.seg(j);
					x.gbl->res.s(sind,mode,0) = 0.0;
				}
			}

			void tadvance(); 
    };

    class neumann : public hp_edge_bdry {
		protected:
			tri_hp_cd &x;
			virtual FLT flux(FLT u, TinyVector<FLT,tri_mesh::ND> x, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm) {return(0.0);}

		public:
			neumann(tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "neumann";}
			neumann(const neumann& inbdry, tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin) {}
			neumann* create(tri_hp& xin, edge_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			void rsdl(int stage);
    };


    class characteristic : public neumann {
		public:
			FLT flux(FLT u, TinyVector<FLT,tri_mesh::ND> pt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm) {
				FLT vel;

				vel =  (x.gbl->ax-mv(0))*norm(0) +(x.gbl->ay -mv(1))*norm(1);        


				if (vel > 0.0)
					return(vel*u);

				return(ibc->f(0, pt, x.gbl->time)*vel);
			}
			characteristic(tri_hp_cd &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic &inbdry, tri_hp_cd &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
    };

    class mixed : public neumann {
		public:
			TinyVector<FLT,5> c;

			FLT flux(FLT u, TinyVector<FLT,tri_mesh::ND> pt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm) {

				FLT fout = 0.0;
				for(int i=0;i<5;++i)
					fout += c(i)*pow(u,i);

				return(fout);
			}
			mixed(tri_hp_cd &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "mixed";}
			mixed(const mixed& inbdry, tri_hp_cd &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), c(inbdry.c) {}
			mixed* create(tri_hp& xin, edge_bdry &bin) const {return new mixed(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}

			void init(input_map& inmap, void* gbl_in) {
				std::string keyword;
				std::istringstream data;
				std::string val;

				neumann::init(inmap,gbl_in);

				keyword = base.idprefix + "_cd_mixed_coefficients";

				if (inmap.getline(keyword,val)) {
						data.str(val);
						data >> c(0) >> c(1) >> c(2) >> c(3) >> c(4);  
						data.clear(); 
				}
				else {
					*x.gbl->log << "couldn't find coefficients" << std::endl;
				}
				return;
			}
    };
}

