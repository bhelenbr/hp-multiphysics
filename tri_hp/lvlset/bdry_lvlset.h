/*
 *  lvlset_bdry.h
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

#ifndef _bdry_lvlset_h_
#define _bdry_lvlset_h_


#include "tri_hp_lvlset.h"
#include "../hp_boundary.h"
#include "../ins/bdry_ins.h"
#include <myblas.h>

namespace bdry_lvlset {

	template <class BASE> class characteristic : public BASE {
		protected:
			tri_hp_lvlset &x;
		public:
			characteristic(tri_hp_lvlset &xin, edge_bdry &bin) : BASE(xin,bin), x(xin) {BASE::mytype="characteristic";}
			characteristic(const characteristic<BASE>& inbdry, tri_hp_lvlset &xin, edge_bdry &bin) : BASE(inbdry,xin,bin), x(xin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic<BASE>(*this,dynamic_cast<tri_hp_lvlset&>(xin),bin);}

			void vdirichlet() {
				int sind,v0,v1;
				TinyVector<FLT,2> nrm, vel;
				FLT normvel;

				BASE::vdirichlet();

				for(int j=0;j<BASE::base.nseg;++j) {
					sind = BASE::base.seg(j);
					v0 = x.seg(sind).pnt(0);
					v1 = x.seg(sind).pnt(1);

					nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
					nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));

					vel(0) = 0.5*(x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0))) +
							  x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0))));
					vel(1) = 0.5*(x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1))) +
							  x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1))));
#ifdef MESH_REF_VEL
					vel -= x.gbl->mesh_ref_vel;
#endif
					normvel = vel(0)*nrm(0)+vel(1)*nrm(1);

					/* normvel is defined positive outward */
					if (normvel < 0.0) {
						x.gbl->res.v(v0,2) = 0.0;
						x.gbl->res.v(v1,2) = 0.0;
					}
				}
			}

			void sdirichlet(int mode) {
				int sind,v0,v1;
				TinyVector<FLT,2> nrm, vel;
				FLT normvel;

				BASE::sdirichlet(mode);

				for(int j=0;j<BASE::base.nseg;++j) {
					sind = BASE::base.seg(j);
					v0 = x.seg(sind).pnt(0);
					v1 = x.seg(sind).pnt(1);

					nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
					nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));

					vel(0) = 0.5*(x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0))) +
							  x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0))));
					vel(1) = 0.5*(x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1))) +
							  x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1))));
#ifdef MESH_REF_VEL
					vel -= x.gbl->mesh_ref_vel;
#endif
					normvel = vel(0)*nrm(0)+vel(1)*nrm(1);

					if (normvel < 0.0) {
						x.gbl->res.s(sind,mode,2) = 0.0;
					}
				}
			}
	};

	class hybrid : public characteristic<bdry_ins::generic> {                             
		public:
			hybrid(tri_hp_lvlset &xin, edge_bdry &bin) : characteristic<bdry_ins::generic>(xin,bin) {mytype = "hybrid";}
			hybrid(const hybrid& inbdry, tri_hp_lvlset &xin, edge_bdry &bin) : characteristic<bdry_ins::generic>(inbdry,xin,bin) {}
			hybrid* create(tri_hp& xin, edge_bdry &bin) const {return new hybrid(*this,dynamic_cast<tri_hp_lvlset&>(xin),bin);}

			void rsdl(int stage) {};
			void update(int stage);

			void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride);			
			void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride); 
			void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride); 
			void smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride);

	};

	class hybrid_pt : public hp_vrtx_bdry {
		protected:
			tri_hp_lvlset &x;
			bool tangset;
			TinyVector<FLT,tri_mesh::ND> tang;
			
		public:
			hybrid_pt(tri_hp_lvlset &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin), tangset(false) {mytype = "hybrid_pt";}
			hybrid_pt(const hybrid_pt& inbdry, tri_hp_lvlset &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin), tangset(false) {}			
			hybrid_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new hybrid_pt(*this,dynamic_cast<tri_hp_lvlset&>(xin),bin);}

			void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride); 
			void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride);
			
			void vdirichlet() {
				x.gbl->res.v(base.pnt,Range(0,x.ND-1)) = 0.0;
				x.gbl->res.v(base.pnt,x.NV-1) = 0.0;  // Needs to be fixed
			}
			void update(int stage);
			friend class hybrid;
	};


	}
#endif
