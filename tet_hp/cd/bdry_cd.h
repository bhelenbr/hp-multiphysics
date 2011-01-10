/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cd.h"
#include "../hp_boundary.h"
#include "myblas.h"

namespace bdry_cd {
	class dirichlet : public hp_face_bdry {
		tet_hp_cd &x;
		
		public:
			dirichlet(tet_hp_cd &xin, face_bdry &bin) : hp_face_bdry(xin,bin), x(xin) {mytype = "dirichlet";}
			dirichlet(const dirichlet& inbdry, tet_hp_cd &xin, face_bdry &bin) : hp_face_bdry(inbdry,xin,bin), x(xin) {}
			dirichlet* create(tet_hp& xin, face_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
			void vdirichlet() { 
				int v0; 											
				for(int j=0;j<base.npnt;++j) {
					v0 = base.pnt(j).gindx;
						x.gbl->res.v(v0,0) = 0.0;
				}			
			}

			void edirichlet() {
				int sind;
				if (basis::tet(x.log2p).em > 0) {
					for(int j=0;j<base.nseg;++j) {
						sind = base.seg(j).gindx;
						x.gbl->res.e(sind,Range(0,basis::tet(x.log2p).em-1),0) = 0.0;
					}
				}
			}
			
			void fdirichlet() {
				int find;

				if (basis::tet(x.log2p).fm > 0) {
					for(int j=0;j<base.ntri;++j) {
						find = base.tri(j).gindx;
						x.gbl->res.f(find,Range(0,basis::tet(x.log2p).fm-1),0) = 0.0;
					}
				}
			}
//			void apply_sparse_dirichlet(bool compressed_column) {
//				int gind;
//				int em=basis::tet(x.log2p).em;
//				int fm=basis::tet(x.log2p).fm;
//
//				for(int i=0;i<base.npnt;++i)
//					x.sparse_dirichlet(base.pnt(i).gindx,compressed_column);
//				
//				for(int i=0;i<base.nseg;++i){
//					gind = x.npnt+base.seg(i).gindx*em;
//					for(int m=0; m<em; ++m)
//						x.sparse_dirichlet(gind+m,compressed_column);
//				}
//				
//				for(int i=0;i<base.ntri;++i){
//					gind = x.npnt+x.nseg*em+base.tri(i).gindx*fm;
//					for(int m=0; m<fm; ++m)
//						x.sparse_dirichlet(gind+m,compressed_column);
//				}				
//			}

			void tadvance(); 
	};

	class neumann : public hp_face_bdry {
		protected:
			tet_hp_cd &x;
			virtual FLT flux(FLT u, TinyVector<FLT,tet_mesh::ND> x, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm) {return(0.0);}
		
		public:
			neumann(tet_hp_cd &xin, face_bdry &bin) : hp_face_bdry(xin,bin), x(xin) {mytype = "neumann";}
			neumann(const neumann& inbdry, tet_hp_cd &xin, face_bdry &bin) : hp_face_bdry(inbdry,xin,bin), x(xin) {}
			neumann* create(tet_hp& xin, face_bdry &bin) const {return new neumann(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
//            void rsdl(int stage);
	};


//   class characteristic : public neumann {
//        public:
//            FLT flux(FLT u, TinyVector<FLT,tet_mesh::ND> pt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm) {
//                FLT vel;
//
//                vel =  (x.gbl->ax-mv(0))*norm(0) +(x.gbl->ay -mv(1))*norm(1);        
//
//
//                if (vel > 0.0)
//                    return(vel*u);
//
//                return(x.gbl->ibc->f(0, pt, x.gbl->time)*vel);
//            }
//            characteristic(tet_hp_cd &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
//            characteristic(const characteristic &inbdry, tet_hp_cd &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
//            characteristic* create(tet_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
//    };
//
//    class mixed : public neumann {
//        public:
//            TinyVector<FLT,5> c;
//            
//            FLT flux(FLT u, TinyVector<FLT,tet_mesh::ND> pt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm) {
//                
//                FLT fout = 0.0;
//                for(int i=0;i<5;++i)
//                    fout += c(i)*pow(u,i);
//                
//                return(fout);
//            }
//            mixed(tet_hp_cd &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "mixed";}
//            mixed(const mixed& inbdry, tet_hp_cd &xin, edge_bdry &bin) : neumann(inbdry,xin,bin), c(inbdry.c) {}
//            mixed* create(tet_hp& xin, edge_bdry &bin) const {return new mixed(*this,dynamic_cast<tet_hp_cd&>(xin),bin);}
//            
//            void init(input_map& inmap, void* gbl_in) {
//                std::string keyword;
//                std::istringstream data;
//                std::string val;
//                
//                neumann::init(inmap,gbl_in);
//
//                keyword = base.idprefix + "_cd_mixed_coefficients";
//                
//                if (inmap.getline(keyword,val)) {
//                        data.str(val);
//                        data >> c(0) >> c(1) >> c(2) >> c(3) >> c(4);  
//                        data.clear(); 
//                }
//                else {
//                    *x.gbl->log << "couldn't find coefficients" << std::endl;
//                }
//                return;
//            }
//    };
}

