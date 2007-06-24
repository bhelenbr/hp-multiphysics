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


#include "tri_hp_ins.h"
#include "../hp_boundary.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>
#include <mathclass.h>


namespace bdry_ins {

    class generic : public hp_side_bdry {
        protected:
            tri_hp_ins &x;
            Array<FLT,1> fluxstorage;
            bool report_flag;
        
        public:
            generic(tri_hp_ins &xin, side_bdry &bin) : hp_side_bdry(xin,bin), x(xin) {mytype = "generic";}
            generic(const generic& inbdry, tri_hp_ins &xin, side_bdry &bin) : hp_side_bdry(inbdry,xin,bin), x(xin) {}
            generic* create(tri_hp& xin, side_bdry &bin) const {return new generic(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            void init(input_map& input,void* &gbl_in) {
                hp_side_bdry::init(input,gbl_in);
                std::string keyword = base.idprefix +"_report";
                input.getwdefault(keyword,report_flag,false);
                
                if (report_flag) {
                    fluxstorage.resize(x.NV);
                }
            }
            void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);
};
        

    class neumann : public generic {
        protected:
            virtual void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx) {
                
                /* CONTINUITY */
                flx(x.NV-1) = x.gbl_ptr->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
              
                /* X&Y MOMENTUM */
#ifdef INERTIALESS
                for (int n=0;n<mesh::ND;++n)
                    flx(n) = x.gbl_ptr->ibc->f(x.NV-1, xpt)*norm(n);
#else
                for (int n=0;n<mesh::ND;++n)
                    flx(n) = flx(x.NV-1)*u(n) +x.gbl_ptr->ibc->f(x.NV-1, xpt)*norm(n);
#endif

              
                /* EVERYTHING ELSE */
                 for (int n=mesh::ND;n<x.NV-1;++n)
                    flx(n) = flx(x.NV-1)*u(n);
                    
                return;
            }
        
        public:
            neumann(tri_hp_ins &xin, side_bdry &bin) : generic(xin,bin) {mytype = "neumann";}
            neumann(const neumann& inbdry, tri_hp_ins &xin, side_bdry &bin) : generic(inbdry,xin,bin) {}
            neumann* create(tri_hp& xin, side_bdry &bin) const {return new neumann(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            block::ctrl rsdl(block::ctrl ctrl_message);
    };



    class inflow : public neumann {        
        void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm,  Array<FLT,1>& flx) {
            
            /* CONTINUITY */
            flx(x.NV-1) = x.gbl_ptr->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));

            /* EVERYTHING ELSE */
             for (int n=0;n<x.NV-1;++n)
                flx(n) = 0.0;
                
            return;
        }
        
        public:
            inflow(tri_hp_ins &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "inflow";}
            inflow(const inflow& inbdry, tri_hp_ins &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
            inflow* create(tri_hp& xin, side_bdry &bin) const {return new inflow(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            void vdirichlet() {
                int sind,v0;
                            
                for(int j=0;j<base.nel;++j) {
                    sind = base.el(j);
                    v0 = x.sd(sind).vrtx(0);
                    x.gbl_ptr->res.v(v0,Range(0,x.NV-2)) = 0.0;
                }
                v0 = x.sd(sind).vrtx(1);
                x.gbl_ptr->res.v(v0,Range(0,x.NV-2)) = 0.0;
            }
            
            void sdirichlet(int mode) {
                int sind;

                for(int j=0;j<base.nel;++j) {
                    sind = base.el(j);
                    x.gbl_ptr->res.s(sind,mode,Range(0,x.NV-2)) = 0.0;
                }
            }
                
            block::ctrl tadvance(bool coarse, block::ctrl ctrl_message);
    };
    
    class euler : public neumann {
        protected:
            void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm,  Array<FLT,1>& flx) {
                Array<FLT,1> ub(x.NV);
                
                for(int n=0;n<x.NV;++n)
                    ub(n) = x.gbl_ptr->ibc->f(n,xpt);
                
                flx(x.NV-1) = x.gbl_ptr->rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1));
                
                /* X&Y MOMENTUM */
                for (int n=0;n<mesh::ND;++n)
                    flx(n) = flx(x.NV-1)*ub(n) +u(x.NV-1)*norm(n);
              
                /* EVERYTHING ELSE */
                 for (int n=mesh::ND;n<x.NV-1;++n)
                    flx(n) = flx(x.NV-1)*ub(n);
                
                return;
            }
        public:
            euler(tri_hp_ins &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "euler";}
            euler(const euler& inbdry, tri_hp_ins &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
            euler* create(tri_hp& xin, side_bdry &bin) const {return new euler(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
    };
    
    class characteristic : public neumann {
        protected:
            void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx);
        public:
            characteristic(tri_hp_ins &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "characteristic";}
            characteristic(const characteristic& inbdry, tri_hp_ins &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
            characteristic* create(tri_hp& xin, side_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
    };
    
    class symmetry : public generic {
        public:
            symmetry(tri_hp_ins &xin, side_bdry &bin) : generic(xin,bin) {mytype = "symmetry";}
            symmetry(const symmetry& inbdry, tri_hp_ins &xin, side_bdry &bin) : generic(inbdry,xin,bin) {}
            symmetry* create(tri_hp& xin, side_bdry &bin) const {return new symmetry(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}

            void vdirichlet() {
                int sind,v0;
                            
                for(int j=0;j<base.nel;++j) {
                    sind = base.el(j);
                    v0 = x.sd(sind).vrtx(0);
                    x.gbl_ptr->res.v(v0,0) = 0.0;
                }
                v0 = x.sd(sind).vrtx(1);
                x.gbl_ptr->res.v(v0,0) = 0.0;
            }
            
            void sdirichlet(int mode) {
                int sind;

                for(int j=0;j<base.nel;++j) {
                    sind = base.el(j);
                    x.gbl_ptr->res.s(sind,mode,0) = 0.0;
                }
            }
                
            block::ctrl tadvance(bool coarse, block::ctrl ctrl_message);
    };
    
    class applied_stress : public neumann {
        Array<symbolic_function<2>,1> stress;

        protected:
            void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx) {
                
                /* CONTINUITY */
                flx(x.NV-1) = x.gbl_ptr->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
              
                FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
                /* X&Y MOMENTUM */
#ifdef INERTIALESS
                for (int n=0;n<mesh::ND;++n)
                    flx(n) = -stress(n).Eval(xpt,sim::time)*length +x.gbl_ptr->ibc->f(x.NV-1, xpt)*norm(n);
#else
                for (int n=0;n<mesh::ND;++n)
                    flx(n) = flx(x.NV-1)*u(n) -stress(n).Eval(xpt,sim::time)*length +x.gbl_ptr->ibc->f(x.NV-1, xpt)*norm(n);
#endif
              
                /* EVERYTHING ELSE */
                 for (int n=mesh::ND;n<x.NV-1;++n)
                    flx(n) = flx(x.NV-1)*u(n);
                    
                return;
            }
        public:
            applied_stress(tri_hp_ins &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "applied_stress";}
            applied_stress(const applied_stress& inbdry, tri_hp_ins &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
            applied_stress* create(tri_hp& xin, side_bdry &bin) const {return new applied_stress(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            void init(input_map& inmap,void* &gbl_in);
    };
    
    class surface_slave : public neumann {
        private:
            int excpt, excpt1, stage;
            
        public:
            surface_slave(tri_hp_ins &xin, side_bdry &bin) : neumann(xin,bin) {mytype = "surface_slave";}
            surface_slave(const surface_slave& inbdry, tri_hp_ins &xin, side_bdry &bin) : neumann(inbdry,xin,bin) {}
            surface_slave* create(tri_hp& xin, side_bdry &bin) const {return new surface_slave(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            void init(input_map& inmap,void* &gbl_in);

            /* FOR COUPLED DYNAMIC BOUNDARIES */
            block::ctrl rsdl(block::ctrl ctrl_message);
            block::ctrl update(block::ctrl ctrl_message);
                      
            void vmatchsolution_snd(int phase, FLT *vdata, int vrtstride=1) {base.vloadbuff(boundary::all,vdata,0,x.NV-2,vrtstride*x.NV);}
            void vmatchsolution_rcv(int phase, FLT *vdata, int vrtstride=1) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,vdata,0,x.NV-2, x.NV*vrtstride);}
            void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride); 
            void smatchsolution_rcv(FLT *sdata, int bgnmode, int endmode, int modestride);
    };
    

    class surface : public surface_slave {
        protected:
            Array<FLT,1> ksprg;
            Array<TinyVector<FLT,mesh::ND>,1> vug_frst;
            Array<TinyVector<FLT,mesh::ND>,2> vdres; //!< Driving term for multigrid (log2p, vrtx)
            Array<TinyVector<FLT,mesh::ND>,3> sdres; //!< Driving term for multigrid (log2p, side, order)
        
        private:
            int mp_phase,excpt,excpt1,stage;
                     
        public:
            struct gbl {                
                bool is_loop;
                /* FLUID PROPERTIES */
                FLT sigma,rho2,mu2;

                /* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
                Array<TinyVector<FLT,mesh::ND>,1> vug0;
                Array<TinyVector<FLT,mesh::ND>,2> sug0;
                
                /* RESIDUALS */
                Array<TinyVector<FLT,mesh::ND>,1> vres;
                Array<TinyVector<FLT,mesh::ND>,2> sres;
                Array<TinyVector<FLT,mesh::ND>,1> vres0;
                Array<TinyVector<FLT,mesh::ND>,2> sres0;
#ifdef DROP
                FLT penalty,vflux;
                Array<FLT,1> vvolumeflux;
                Array<FLT,2> svolumeflux;
#endif
                
                /* PRECONDITIONER */
                Array<TinyMatrix<FLT,mesh::ND,mesh::ND>,1> vdt;
                Array<TinyMatrix<FLT,mesh::ND,mesh::ND>,1> sdt;
                Array<FLT,1> meshc;
                
                TinyVector<FLT,mesh::ND> fadd;
                TinyMatrix<FLT,mesh::ND,MAXP> cfl;
                FLT adis;
            } *surf_gbl; 

        public:
            surface(tri_hp_ins &xin, side_bdry &bin) : surface_slave(xin,bin) {mytype = "surface";}
            surface(const surface& inbdry, tri_hp_ins &xin, side_bdry &bin) : surface_slave(inbdry,xin,bin) {}
            surface* create(tri_hp& xin, side_bdry &bin) const {return new surface(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}

            void init(input_map& input,void* &gbl_in); 

            /* FOR COUPLED DYNAMIC BOUNDARIES */
            block::ctrl tadvance(bool coarse,block::ctrl ctrl_message);
            block::ctrl rsdl(block::ctrl ctrl_message);
            void maxres();
            block::ctrl setup_preconditioner(block::ctrl ctrl_message);
            block::ctrl minvrt(block::ctrl ctrl_message);
            block::ctrl update(block::ctrl ctrl_message);
            block::ctrl mg_getfres(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh, int bnum); 
    };
    
    class hybrid_surface_levelset : public surface_slave { 
        protected:
            virtual void flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx) {
                
                /* CONTINUITY */
                flx(x.NV-1) = x.gbl_ptr->rho*((u(0) -mv(0))*norm(0) +(u(1) -mv(1))*norm(1));
              
                /* EVERYTHING ELSE */
                 for (int n=0;n<x.NV-1;++n)
                    flx(n) = 0.0;
                    
                return;
            }      
        private:
            int mp_phase,excpt,excpt1,stage;
                     
        public:
            hybrid_surface_levelset(tri_hp_ins &xin, side_bdry &bin) : surface_slave(xin,bin) {mytype = "hybrid_surface_levelset";}
            hybrid_surface_levelset(const hybrid_surface_levelset& inbdry, tri_hp_ins &xin, side_bdry &bin) : surface_slave(inbdry,xin,bin) {}
            hybrid_surface_levelset* create(tri_hp& xin, side_bdry &bin) const {return new hybrid_surface_levelset(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            void init(input_map& inmap,void* &gbl_in) {
                std::string keyword;
                
                keyword = base.idprefix + "_curved";
                inmap[keyword] = "0";

                neumann::init(inmap,gbl_in);
                
                return;
            }
            
            block::ctrl update(block::ctrl ctrl_message) {return(block::stop);}
            block::ctrl rsdl(block::ctrl ctrl_message) {return(neumann::rsdl(ctrl_message));}

    };



    /*******************************/
    /* VERTEX BOUNDARY CONDITIONS */
    /******************************/
    class surface_fixed_pt : public hp_vrtx_bdry {
        protected:
            tri_hp_ins &x;
            surface *surf;
            int surfbdry;
            int dirstart,dirstop;
            
         public:
            surface_fixed_pt(tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "surface_fixed_pt";}
            surface_fixed_pt(const surface_fixed_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin), 
                surfbdry(inbdry.surfbdry), dirstart(inbdry.dirstart), dirstop(inbdry.dirstop) {}
            surface_fixed_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_fixed_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            void init(input_map& inmap,void* &gbl_in) {
                std::string keyword,val;
                std::istringstream data;
                std::string filename;

                hp_vrtx_bdry::init(inmap,gbl_in);
                
                if (surf = dynamic_cast<surface *>(x.hp_sbdry(base.sbdry(0)))) {
                    surfbdry = 0;
                }
                else if (surf = dynamic_cast<surface *>(x.hp_sbdry(base.sbdry(1)))) {
                    surfbdry = 1;
                }
                else {
                    *sim::log << "something's wrong neither side is a surface boundary" << std::endl;
                    exit(1);
                }
                
                keyword = base.idprefix + "_dirstart";
                inmap.getwdefault(keyword,dirstart,0);
                
                keyword = base.idprefix + "_dirstop";
                inmap.getwdefault(keyword,dirstop,1);
            }
            
            block::ctrl rsdl(block::ctrl ctrl_message) {
                if (ctrl_message == block::begin) {
                    if (surfbdry == 0) {
                        /* SET TANGENT RESIDUAL TO ZERO */
                        surf->surf_gbl->vres(x.sbdry(base.sbdry(0))->nel)(0) = 0.0;
                        if (dirstop > 0) {
                            /* POST-REMOVE ADDED MASS FLUX TERM FOR FIXED POINT */
                            x.gbl_ptr->res.v(base.v0,x.NV-1) += surf->surf_gbl->vres(x.sbdry(base.sbdry(0))->nel)(1)*x.gbl_ptr->rho;
                            /* AND ZERO RESIDUAL */
                            surf->surf_gbl->vres(x.sbdry(base.sbdry(0))->nel)(1) = 0.0;
                        }
                    }
                    else {
                        /* SET TANGENT RESIDUAL TO ZERO */
                        surf->surf_gbl->vres(0)(0) = 0.0;
                        if (dirstop > 0) {
                            /* POST-REMOVE ADDED MASS FLUX TERM FOR FIXED POINT */
                            x.gbl_ptr->res.v(base.v0,x.NV-1) += surf->surf_gbl->vres(0)(1)*x.gbl_ptr->rho;
                            /* AND ZERO RESIDUAL */
                            surf->surf_gbl->vres(0)(1) = 0.0;
                        }
                    }
                }
                return(block::stop);
            }
            
            void vdirichlet() {
                if (surfbdry == 0) {
                    for(int n=dirstart;n<=dirstop;++n) 
                        surf->surf_gbl->vres(x.sbdry(base.sbdry(0))->nel)(n) = 0.0;
                }
                else {
                    for(int n=dirstart;n<=dirstop;++n) 
                        surf->surf_gbl->vres(0)(n) = 0.0;
                }
            }
            
            void mvpttobdry(TinyVector<FLT,mesh::ND> &pt) {
                if (surfbdry == 0) {
                    x.sbdry(base.sbdry(1))->mvpttobdry(0,-1.0,pt);
                }
                else {
                    x.sbdry(base.sbdry(0))->mvpttobdry(x.sbdry(base.sbdry(0))->nel-1,1.0,pt);
                }
            }
    };
    
    class surface_periodic_pt : public surface_fixed_pt {
        protected:
            FLT position;
            bool vertical;
            
        public:
            surface_periodic_pt(tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(xin,bin) {mytype = "surface_periodic_pt";}
            surface_periodic_pt(const surface_periodic_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(xin,bin), position(inbdry.position), vertical(inbdry.vertical) {*sim::log << position << std::endl;}
            surface_periodic_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_periodic_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            
            void init(input_map& inmap,void* &gbl_in) {
                bool coarse,first;
                ostringstream nstr;
                
                surface_fixed_pt::init(inmap,gbl_in);
                dirstart = 0;
                dirstop = 0;
                inmap.getwdefault(base.idprefix + "_vertical",vertical,true); 
                inmap.getwdefault(x.idprefix + "_coarse",coarse,false);
                if (!coarse) {
                    position = x.vrtx(base.v0)(1-vertical);
                    nstr.str("");
                    nstr << position << std::flush;
                    inmap.getwdefault(base.idprefix+"_first",first,true);
                    if (first) {
                        inmap[base.idprefix+"_first"] = "0";
                        inmap[base.idprefix+"_prdc_pos_first"] = nstr.str();
                    }
                    else {
                        inmap[base.idprefix+"_first"] = "1";
                        inmap[base.idprefix+"_prdc_pos_second"] = nstr.str();
                    }
                }
                else {
                    if (!inmap.get(base.idprefix+"_first",first)) {
                        *sim::log << "Something went wrong finding periodic point position" << std::endl;
                        exit(1);
                    }
                    if (first) {
                        inmap[base.idprefix+"_first"] = "0";
                        inmap.get(base.idprefix+"_prdc_pos_first",position);
                    }
                    else {
                        inmap[base.idprefix+"_first"] = "1";
                        inmap.get(base.idprefix+"_prdc_pos_second",position);
                    }
                }
            }
            
            void mvpttobdry(TinyVector<FLT,mesh::ND> &pt) {
                pt(1-vertical) = position;
            }
    };
    

    class surface_outflow_endpt : public surface_fixed_pt {
         public:
            surface_outflow_endpt(tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(xin,bin) {mytype = "surface_outflow_endpt";}
            surface_outflow_endpt(const surface_outflow_endpt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt(xin,bin) {}
            surface_outflow_endpt* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_outflow_endpt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            void init(input_map& input,void* &gbl_in) {
                surface_fixed_pt::init(input,gbl_in);
                dirstart = 0;
                dirstop = 0;
            }
                
            /* FOR COUPLED DYNAMIC BOUNDARIES */
            block::ctrl rsdl(block::ctrl ctrl_message);
    };
    
    class surface_outflow_planar : public surface_outflow_endpt {
        protected:
            bool vertical;
            FLT position;
            
        public:
            surface_outflow_planar(tri_hp_ins &xin, vrtx_bdry &bin) : surface_outflow_endpt(xin,bin) {mytype = "surface_outflow_planar";}
            surface_outflow_planar(const surface_outflow_planar& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_outflow_endpt(xin,bin), vertical(inbdry.vertical), position(inbdry.position) {}
            surface_outflow_planar* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_outflow_planar(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
            void init(input_map& inmap,void* &gbl_in) {
                std::string keyword;
                std::ostringstream nstr;
                
                surface_outflow_endpt::init(inmap,gbl_in);

                keyword = base.idprefix + "_vertical";
                inmap.getwdefault(keyword,vertical,true);
                     
                keyword = base.idprefix +"_position";
                if (!inmap.get(keyword,position)) {
                    position = x.vrtx(base.v0)(1-vertical);
                    nstr.str("");
                    nstr << position;
                    inmap[keyword] = nstr.str();
                }
            }
                
            void mvpttobdry(TinyVector<FLT,mesh::ND> &pt) {
                pt(1-vertical) = position;
            }
    };

    class inflow_pt : public hp_vrtx_bdry {
        protected:
            tri_hp_ins &x;
            
         public:
            inflow_pt(tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin) {mytype = "inflow_pt";}
            inflow_pt(const inflow_pt& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin) {}
            inflow_pt* create(tri_hp& xin, vrtx_bdry &bin) const {return new inflow_pt(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}

            block::ctrl tadvance(block::ctrl ctrl_message) { 
                if (ctrl_message == block::begin) {
                    for(int n=0;n<x.NV-1;++n)
                        x.ug.v(base.v0,n) = x.gbl_ptr->ibc->f(n,x.vrtx(base.v0));  
                }
                return(block::stop);
            }
                                    
            void vdirichlet2d() {
                x.gbl_ptr->res.v(base.v0,Range(0,x.NV-2)) = 0.0;
            }
    };

    
    
    
    
/*    clase surface_vertex_comm : public hp_vrtx_bdry { */
                          /* THIS IS VERY COMPLICATED TO DO THIS WAY */
            /* FOR NOW USE PHASED SWEEPING OF BOUNDARY SIDES (NO COMMUNICATION THROUGH VERTICES) */
            /*
            void vmatchsolution_snd(int phase, FLT *vdata) {base.vloadbuff(boundary::all,vdata,0,x.NV-2,x.NV);}
            void vmatchsolution_rcv(int phase, FLT *vdata) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,vdata,0,x.NV-2,x.NV);}
            */
/*    };*/

}
#endif
