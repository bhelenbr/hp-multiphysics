/*
 *  ins_bdry.h
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
 
#ifndef _bdry_ins_h_
#define _bdry_ins_h_


#include "tet_hp_ins.h"
#include "../hp_boundary.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>
#include <symbolic_function.h>

//#define DETAILED_DT
//#define DETAILED_MINV


namespace bdry_ins {

    class generic : public hp_face_bdry {
        protected:
            tet_hp_ins &x;
            bool report_flag;
        public:
            Array<FLT,1> total_flux,diff_flux,conv_flux;
            FLT circumference,moment,convect,circulation;

        public:
            generic(tet_hp_ins &xin, edge_bdry &bin) : hp_face_bdry(xin,bin), x(xin) {mytype = "generic";}
            generic(const generic& inbdry, tet_hp_ins &xin, edge_bdry &bin) : hp_face_bdry(inbdry,xin,bin), x(xin), report_flag(inbdry.report_flag) {
                if (report_flag) {
                    total_flux.resize(x.NV);
                    diff_flux.resize(x.NV);
                    conv_flux.resize(x.NV);  
                }
            }
            generic* create(tet_hp& xin, edge_bdry &bin) const {return new generic(*this,dynamic_cast<tet_hp_ins&>(xin),bin);}
            void init(input_map& input,void* gbl_in) {
                hp_face_bdry::init(input,gbl_in);
                std::string keyword = base.idprefix +"_report";
                input.getwdefault(keyword,report_flag,false);
                
                if (report_flag) {
                    total_flux.resize(x.NV);
                    diff_flux.resize(x.NV);
                    conv_flux.resize(x.NV);            
                }
            }
            void output(std::ostream& fout, tet_hp::filetype typ,int tlvl = 0);
};
        

    class neumann : public generic {
        protected:
            virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
                
                /* CONTINUITY */
                flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1))+(u(2) -mv(2))*norm(2));
              
                /* X&Y MOMENTUM */
#ifdef INERTIALESS
                for (int n=0;n<tet_mesh::ND;++n)
                    flx(n) = ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#else
                for (int n=0;n<tet_mesh::ND;++n)
                    flx(n) = flx(x.NV-1)*u(n) +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#endif

              
                /* EVERYTHING ELSE */
                 for (int n=tet_mesh::ND;n<x.NV-1;++n)
                    flx(n) = flx(x.NV-1)*u(n);
                    
                return;
            }
        
        public:
            neumann(tet_hp_ins &xin, face_bdry &bin) : generic(xin,bin) {mytype = "neumann";}
            neumann(const neumann& inbdry, tet_hp_ins &xin, face_bdry &bin) : generic(inbdry,xin,bin) {}
            neumann* create(tet_hp& xin, face_bdry &bin) const {return new neumann(*this,dynamic_cast<tet_hp_ins&>(xin),bin);}
            void rsdl(int stage);
    };



    class inflow : public neumann {        
        void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {

            /* CONTINUITY */
            flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1))+(u(2) -mv(2))*norm(2));

            /* EVERYTHING ELSE DOESN'T MATTER */
             for (int n=0;n<x.NV-1;++n)
                flx(n) = 0.0;
                
            return;
        }
        
        public:
            inflow(tet_hp_ins &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "inflow";}
            inflow(const inflow& inbdry, tet_hp_ins &xin, face_bdry &bin) : neumann(inbdry,xin,bin) {}
            inflow* create(tet_hp& xin, face_bdry &bin) const {return new inflow(*this,dynamic_cast<tet_hp_ins&>(xin),bin);}
			
            void vdirichlet() {
                int sind,v0;
                            
                for(int j=0;j<base.nseg;++j) {
                    sind = base.seg(j);
                    v0 = x.seg(sind).pnt(0);
                    x.gbl->res.v(v0,Range(0,x.NV-2)) = 0.0;
                }
                v0 = x.seg(sind).pnt(1);
                x.gbl->res.v(v0,Range(0,x.NV-2)) = 0.0;
            }
            
            void edirichlet() {
                int sind;

                for(int j=0;j<base.nseg;++j) {
                    sind = base.seg(j);
                    x.gbl->res.s(sind,Range(0,x.tet_basis(log2p).em-1),Range(0,x.NV-2)) = 0.0;
                }
            }
			
			void fdirichlet(int mode) {
                int find;

                for(int j=0;j<base.ntri;++j) {
                    find = base.tri(j);
                    x.gbl->res.f(find,mode,Range(0,x.NV-2)) = 0.0;
                }
            }
			
			void setvalues(init_bdry_cndtn *ibc);
			void tadvance() {
				hp_face_bdry::tadvance();
				setvalues(ibc);
			};
    };
	
	
	class rigid : public init_bdry_cndtn {
		public:
//			FLT omega;
//			TinyVector<FLT,2> vel, ctr;
//			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
//				FLT r,cost,sint,v; 
//				TinyVector<FLT,2> dx;
//				
//				dx = x-ctr;
//				r = sqrt(dx(0)*dx(0) +dx(1)*dx(1));
//				cost = dx(0)/r;
//				sint = dx(1)/r;
//				if(n == 0)
//					v = -r*omega*sint +vel(0);
//				else
//					v = r*omega*cost +vel(1);
//					
//				return(v);
//			}
//			rigid() : omega(0.0), vel(0.0), ctr(0.0) {}
	};		
    
	class force_coupling : public inflow {
		rigid my_ibc;
		
		void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
             for (int n=0;n<x.NV;++n)
                flx(n) = 0.0;
            return;
        }
		
		public:
			force_coupling(tet_hp_ins &xin, face_bdry &bin) : inflow(xin,bin) {mytype = "force_coupling";}
			force_coupling(const force_coupling& inbdry, tet_hp_ins &xin, face_bdry &bin) : inflow(inbdry,xin,bin) {}
			force_coupling* create(tet_hp& xin, face_bdry &bin) const {return new force_coupling(*this,dynamic_cast<tet_hp_ins&>(xin),bin);}
			void tadvance() {hp_face_bdry::tadvance();}
			void update(int stage) {
				if (!x.coarse_flag) setvalues(&my_ibc);
			}
			void set_omega(FLT dwdt) {my_ibc.omega = dwdt;}
			void set_ctr_rot(TinyVector<FLT,3> c) {my_ibc.ctr = c;}
			void set_vel(TinyVector<FLT,3> v) {my_ibc.vel = v;}
			void input(ifstream& fin,tet_hp::filetype typ,int tlvl) {
				inflow::input(fin,typ,tlvl);
				if (typ == tet_hp::text) {
					fin >> my_ibc.omega >> my_ibc.vel >> my_ibc.ctr;
				}
			}
			void output(std::ostream& fout, tet_hp::filetype typ,int tlvl = 0) {
				inflow::output(fout,typ,tlvl);
				if (typ == tet_hp::text) {
					fout << my_ibc.omega << ' ' << my_ibc.vel << ' ' << my_ibc.ctr << '\n';
				}
			}

	};
			
	
		
		
    class euler : public neumann {
        protected:
            void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
                Array<FLT,1> ub(x.NV);
                
                for(int n=0;n<x.NV;++n)
                    ub(n) = ibc->f(n,xpt,x.gbl->time);
                
                flx(x.NV-1) = x.gbl->rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1)+(ub(2) -mv(2))*norm(2));
                
                /* X&Y MOMENTUM */
                for (int n=0;n<tet_mesh::ND;++n)
                    flx(n) = flx(x.NV-1)*ub(n) +u(x.NV-1)*norm(n);
              
                /* EVERYTHING ELSE */
                 for (int n=tet_mesh::ND;n<x.NV-1;++n)
                    flx(n) = flx(x.NV-1)*ub(n);
                
                return;
            }
        public:
            euler(tet_hp_ins &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "euler";}
            euler(const euler& inbdry, tet_hp_ins &xin, face_bdry &bin) : neumann(inbdry,xin,bin) {}
            euler* create(tet_hp& xin, face_bdry &bin) const {return new euler(*this,dynamic_cast<tet_hp_ins&>(xin),bin);}
    };
    
    class characteristic : public neumann {
        protected:
            void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx);
        public:
            characteristic(tet_hp_ins &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
            characteristic(const characteristic& inbdry, tet_hp_ins &xin, face_bdry &bin) : neumann(inbdry,xin,bin) {}
            characteristic* create(tet_hp& xin, face_bdry &bin) const {return new characteristic(*this,dynamic_cast<tet_hp_ins&>(xin),bin);}
    };
    
    class symmetry : public generic {
        int dir;
        
        public:
            symmetry(tet_hp_ins &xin, face_bdry &bin) : generic(xin,bin) {mytype = "symmetry";}
            symmetry(const symmetry& inbdry, tet_hp_ins &xin, face_bdry &bin) : generic(inbdry,xin,bin), dir(inbdry.dir) {}
            symmetry* create(tet_hp& xin, face_bdry &bin) const {return new symmetry(*this,dynamic_cast<tet_hp_ins&>(xin),bin);}
            void init(input_map& input,void* gbl_in) {
                generic::init(input,gbl_in);
                std::string keyword = base.idprefix +"_dir";
                input.getwdefault(keyword,dir,0);
            }

            void vdirichlet() {
                int sind,v0;
                            
                for(int j=0;j<base.nseg;++j) {
                    sind = base.seg(j);
                    v0 = x.seg(sind).pnt(0);
                    x.gbl->res.v(v0,dir) = 0.0;
                }
                v0 = x.seg(sind).pnt(1);
                x.gbl->res.v(v0,dir) = 0.0;
            }
            
            void edirichlet() {
                int sind;

                for(int j=0;j<base.nseg;++j) {
                    sind = base.seg(j);
                    x.gbl->res.s(sind,Range(0,x.tet_basis(log2p).em-1),dir) = 0.0;
                }
            }
			
			void fdirichlet() {
                int find;

                for(int j=0;j<base.ntri;++j) {
                    find = base.tri(j);
                    x.gbl->res.f(find,Range(0,x.tet_basis(log2p).fm-1),dir) = 0.0;
                }
            }
                
            void tadvance();
    };
    
    class applied_stress : public neumann {
        Array<symbolic_function<3>,1> stress;

        protected:
            void flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
                
                /* CONTINUITY */
                flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1)+(u(2) -mv(2))*norm(2));
              
                FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1)+norm(2)*norm(2));
                /* X&Y MOMENTUM */
#ifdef INERTIALESS
                for (int n=0;n<tet_mesh::ND;++n)
                    flx(n) = -stress(n).Eval(xpt,x.gbl->time)*length +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#else
                for (int n=0;n<tet_mesh::ND;++n)
                    flx(n) = flx(x.NV-1)*u(n) -stress(n).Eval(xpt,x.gbl->time)*length +ibc->f(x.NV-1, xpt, x.gbl->time)*norm(n);
#endif
              
                /* EVERYTHING ELSE */
                 for (int n=tet_mesh::ND;n<x.NV-1;++n)
                    flx(n) = flx(x.NV-1)*u(n);
                    
                return;
            }
        public:
            applied_stress(tet_hp_ins &xin, face_bdry &bin) : neumann(xin,bin) {mytype = "applied_stress";}
            applied_stress(const applied_stress& inbdry, tet_hp_ins &xin, face_bdry &bin) : neumann(inbdry,xin,bin), stress(inbdry.stress) {}
            applied_stress* create(tet_hp& xin, face_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tet_hp_ins&>(xin),bin);}
            void init(input_map& inmap,void* gbl_in);
    };
    
    class surface_slave : public neumann {
        public:
//            surface_slave(tri_hp_ins &xin, edge_bdry &bin) : neumann(xin,bin) {mytype = "surface_slave";}
//            surface_slave(const surface_slave& inbdry, tri_hp_ins &xin, edge_bdry &bin) : neumann(inbdry,xin,bin) {}
//            surface_slave* create(tri_hp& xin, edge_bdry &bin) const {return new surface_slave(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//            void init(input_map& inmap,void* gbl_in);
//
//            /* FOR COUPLED DYNAMIC BOUNDARIES */
//            void tadvance();
//            void rsdl(int stage);
//            void update(int stage);
//                      
//            void pmatchsolution_snd(int phase, FLT *pdata, int vrtstride=1) {base.ploadbuff(boundary::all,pdata,0,x.NV-2,vrtstride*x.NV);}
//            void pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride=1) {base.pfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,0,x.NV-2, x.NV*vrtstride);}
//            void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride); 
//            void smatchsolution_rcv(FLT *sdata, int bgnmode, int endmode, int modestride);
    };
    

    class surface : public surface_slave {
//        protected:
//            Array<FLT,1> ksprg;
//            Array<TinyVector<FLT,tri_mesh::ND>,1> vug_frst;
//            Array<TinyVector<FLT,tri_mesh::ND>,2> vdres; //!< Driving term for multigrid (log2p, pnts)
//            Array<TinyVector<FLT,tri_mesh::ND>,3> sdres; //!< Driving term for multigrid (log2p, side, order)
//            const surface *fine, *coarse;
//            
//        public:
//            struct global {                
//                bool is_loop;
//                /* FLUID PROPERTIES */
//                FLT sigma,rho2,mu2;
//
//                /* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
//                Array<TinyVector<FLT,tri_mesh::ND>,1> vug0;
//                Array<TinyVector<FLT,tri_mesh::ND>,2> sug0;
//                
//                /* RESIDUALS */
//                Array<TinyVector<FLT,tri_mesh::ND>,1> vres;
//                Array<TinyVector<FLT,tri_mesh::ND>,2> sres;
//                Array<TinyVector<FLT,tri_mesh::ND>,1> vres0;
//                Array<TinyVector<FLT,tri_mesh::ND>,2> sres0;
//#ifdef DROP
//                FLT penalty,vflux;
//                Array<FLT,1> vvolumeflux;
//                Array<FLT,2> svolumeflux;
//#endif
//                
//                /* PRECONDITIONER */
//                Array<TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND>,1> vdt;
//                Array<TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND>,1> sdt;
//                Array<FLT,1> meshc;
//
//#ifdef DETAILED_MINV
//                Array<TinyMatrix<FLT,2*MAXP,2*MAXP>,1> ms;
//                Array<FLT,5> vms;
//                Array<TinyVector<int,2*MAXP>,1> ipiv;
//#endif
//                TinyVector<FLT,tri_mesh::ND> fadd;
//                TinyMatrix<FLT,tri_mesh::ND,MAXP> cfl;
//                FLT adis;
//            } *gbl; 
//
//        public:
//            void* create_global_structure() {return new global;}
//            surface(tri_hp_ins &xin, edge_bdry &bin) : surface_slave(xin,bin) {mytype = "surface";}
//            surface(const surface& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : surface_slave(inbdry,xin,bin) {
//                gbl = inbdry.gbl;
//                ksprg.resize(base.maxseg);
//                vug_frst.resize(base.maxseg+1);
//                vdres.resize(1,base.maxseg+1);
//                fine = &inbdry;
//            };
//            surface* create(tri_hp& xin, edge_bdry &bin) const {return new surface(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//            void init(input_map& input,void* gbl_in); 
//
//            /* FOR COUPLED DYNAMIC BOUNDARIES */
//            void tadvance();
//            void rsdl(int stage);
//            void maxres();
//            void setup_preconditioner();
//            void minvrt();
//            void update(int stage);
//            void mg_restrict(); 
    };
    
    class hybrid_surface_levelset : public surface_slave { 
//        protected:
//            virtual void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
//                
//                /* CONTINUITY */
//                flx(x.NV-1) = x.gbl->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
//              
//                /* EVERYTHING ELSE */
//                 for (int n=0;n<x.NV-1;++n)
//                    flx(n) = 0.0;
//                    
//                return;
//            }      
//   
//        public:
//            hybrid_surface_levelset(tri_hp_ins &xin, edge_bdry &bin) : surface_slave(xin,bin) {mytype = "hybrid_surface_levelset";}
//            hybrid_surface_levelset(const hybrid_surface_levelset& inbdry, tri_hp_ins &xin, edge_bdry &bin) : surface_slave(inbdry,xin,bin) {}
//            hybrid_surface_levelset* create(tri_hp& xin, edge_bdry &bin) const {return new hybrid_surface_levelset(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//            void init(input_map& inmap,void* gbl_in) {
//                std::string keyword;
//                
//                keyword = base.idprefix + "_curved";
//                inmap[keyword] = "0";
//
//                neumann::init(inmap,gbl_in);
//                
//                return;
//            }
//            
//            void update(int stage) {return;}
//            void rsdl(int stage) {neumann::rsdl(stage);}

    };



    /*******************************/
    /* VERTEX BOUNDARY CONDITIONS */
    /******************************/
    class surface_fixed_pt : public hp_vrtx_bdry {
//        protected:
//            tri_hp_ins &x;
//            surface *surf;
//            int surfbdry;
//            int dirstart,dirstop;
//            
//         public:
//            surface_fixed_pt(tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "surface_fixed_pt";}
//            surface_fixed_pt(const surface_fixed_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin), 
//                surfbdry(inbdry.surfbdry), dirstart(inbdry.dirstart), dirstop(inbdry.dirstop) {
//                if (!(surf = dynamic_cast<surface *>(x.hp_ebdry(base.ebdry(surfbdry))))) {
//                    *x.gbl->log << "something's wrong can't find surface boundary" << std::endl;
//                    exit(1);
//                }
//            }
//            surface_fixed_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_fixed_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//            
//            void init(input_map& inmap,void* gbl_in) {
//                std::string keyword,val;
//                std::istringstream data;
//                std::string filename;
//
//                hp_vrtx_bdry::init(inmap,gbl_in);
//                
//                if (surf = dynamic_cast<surface *>(x.hp_ebdry(base.ebdry(0)))) {
//                    surfbdry = 0;
//                }
//                else if (surf = dynamic_cast<surface *>(x.hp_ebdry(base.ebdry(1)))) {
//                    surfbdry = 1;
//                }
//                else {
//                    *x.gbl->log << "something's wrong neither side is a surface boundary" << std::endl;
//                    exit(1);
//                }
//                
//                keyword = base.idprefix + "_dirstart";
//                inmap.getwdefault(keyword,dirstart,0);
//                
//                keyword = base.idprefix + "_dirstop";
//                inmap.getwdefault(keyword,dirstop,1);
//            }
//            
//            void rsdl(int stage) {
//                if (surfbdry == 0) {
//                    /* SET TANGENT RESIDUAL TO ZERO */
//                    surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(0) = 0.0;
//                    if (dirstop > 0) {
//                        /* POST-REMOVE ADDED MASS FLUX TERM FOR FIXED POINT */
//                        x.gbl->res.v(base.pnt,x.NV-1) += surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(1)*x.gbl->rho;
//                        /* AND ZERO RESIDUAL */
//                        surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(1) = 0.0;
//                    }
//                }
//                else {
//                    /* SET TANGENT RESIDUAL TO ZERO */
//                    surf->gbl->vres(0)(0) = 0.0;
//                    if (dirstop > 0) {
//                        /* POST-REMOVE ADDED MASS FLUX TERM FOR FIXED POINT */
//                        x.gbl->res.v(base.pnt,x.NV-1) += surf->gbl->vres(0)(1)*x.gbl->rho;
//                        /* AND ZERO RESIDUAL */
//                        surf->gbl->vres(0)(1) = 0.0;
//                    }
//                }
//                return;
//            }
//            
//            void vdirichlet() {
//                if (surfbdry == 0) {
//                    for(int n=dirstart;n<=dirstop;++n) 
//                        surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg)(n) = 0.0;
//                }
//                else {
//                    for(int n=dirstart;n<=dirstop;++n) 
//                        surf->gbl->vres(0)(n) = 0.0;
//                }
//            }
//            
//            void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
//                if (surfbdry == 0) {
//                    x.ebdry(base.ebdry(1))->mvpttobdry(0,-1.0,pt);
//                }
//                else {
//                    x.ebdry(base.ebdry(0))->mvpttobdry(x.ebdry(base.ebdry(0))->nseg-1,1.0,pt);
//                }
//            }
    };
    
    class surface_periodic_pt : public surface_fixed_pt {
//        protected:
//            FLT position;
//            bool vertical;
//            
//        public:
//            surface_periodic_pt(tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(xin,bin) {mytype = "surface_periodic_pt";}
//            surface_periodic_pt(const surface_periodic_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(inbdry,xin,bin), position(inbdry.position), vertical(inbdry.vertical) {}
//            surface_periodic_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_periodic_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//            
//            void init(input_map& inmap,void* gbl_in) {                
//                surface_fixed_pt::init(inmap,gbl_in);
//                dirstart = 0;
//                dirstop = 0;
//                inmap.getwdefault(base.idprefix + "_vertical",vertical,true); 
//                position = x.pnts(base.pnt)(1-vertical);
//            }
//            
//            void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
//                pt(1-vertical) = position;
//            }
    };
    

    class surface_outflow_endpt : public surface_fixed_pt {
//         public:
//            surface_outflow_endpt(tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(xin,bin) {mytype = "surface_outflow_endpt";}
//            surface_outflow_endpt(const surface_outflow_endpt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(inbdry,xin,bin) {}
//            surface_outflow_endpt* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_outflow_endpt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//            void init(input_map& input,void* gbl_in) {
//                surface_fixed_pt::init(input,gbl_in);
//                dirstart = 0;
//                dirstop = 0;
//            }
//                
//            /* FOR COUPLED DYNAMIC BOUNDARIES */
//            void rsdl(int stage);
    };
    
    class surface_outflow_planar : public surface_outflow_endpt {
//        protected:
//            bool vertical;
//            FLT position;
//            
//        public:
//            surface_outflow_planar(tri_hp_ins &xin, vrtx_bdry &bin) : surface_outflow_endpt(xin,bin) {mytype = "surface_outflow_planar";}
//            surface_outflow_planar(const surface_outflow_planar& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_outflow_endpt(inbdry,xin,bin), vertical(inbdry.vertical), position(inbdry.position) {}
//            surface_outflow_planar* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_outflow_planar(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//            void init(input_map& inmap,void* gbl_in) {
//                std::string keyword;
//                std::ostringstream nstr;
//                
//                surface_outflow_endpt::init(inmap,gbl_in);
//
//                keyword = base.idprefix + "_vertical";
//                inmap.getwdefault(keyword,vertical,true);
//                     
//                keyword = base.idprefix +"_position";
//                if (!inmap.get(keyword,position)) {
//                    position = x.pnts(base.pnt)(1-vertical);
//                    nstr.str("");
//                    nstr << position;
//                    inmap[keyword] = nstr.str();
//                }
//            }
//                
//            void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
//                pt(1-vertical) = position;
//            }
    };

    class inflow_pt : public hp_vrtx_bdry {
//        protected:
//            tri_hp_ins &x;
//            
//         public:
//            inflow_pt(tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "inflow_pt";}
//            inflow_pt(const inflow_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
//            inflow_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new inflow_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
//
//            void tadvance() { 
//                for(int n=0;n<x.NV-1;++n)
//                    x.ug.v(base.pnt,n) = x.gbl->ibc->f(n,x.pnts(base.pnt),x.gbl->time);  
//                return;
//            }
//                                    
//            void vdirichlet2d() {
//                x.gbl->res.v(base.pnt,Range(0,x.NV-2)) = 0.0;
//            }
    };

    
    
    
    
/*    clase surface_vertex_comm : public hp_vrtx_bdry { */
                          /* THIS IS VERY COMPLICATED TO DO THIS WAY */
            /* FOR NOW USE PHASED SWEEPING OF BOUNDARY SIDES (NO COMMUNICATION THROUGH VERTICES) */
            /*
            void pmatchsolution_snd(int phase, FLT *pdata) {base.ploadbuff(boundary::all,pdata,0,x.NV-2,x.NV);}
            void pmatchsolution_rcv(int phase, FLT *pdata) {base.pfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,pdata,0,x.NV-2,x.NV);}
            */
/*    };*/

}
#endif
