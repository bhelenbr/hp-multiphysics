/*
 *  hp_boundary.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 9/3/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _hp_boundary_h_
#define _hp_boundary_h_
 
#include <boundary.h>
#include <hpbasis.h>
#include <myblas.h>

class hp_side_bdry;

class hp_vrtx_bdry : public vgeometry_interface {
    protected:
        std::string mytype;
        tri_hp& x;
        vrtx_bdry& base;
        static FLT dummy;
        int excpt,excpt1;
        
    public:
        hp_vrtx_bdry(tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin) {mytype = "plain";}
        hp_vrtx_bdry(const hp_vrtx_bdry &inbdry,tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin), mytype(inbdry.mytype) {}
        virtual hp_vrtx_bdry* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_vrtx_bdry(*this,xin,bin);}
        virtual void init(input_map& input,void* &gbl_in) {} /**< This is to read definition data only (not solution data) */
        virtual void copy(const hp_vrtx_bdry& tgt) {}
        virtual ~hp_vrtx_bdry() {}
        
        /* input output functions */
        virtual void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
            switch(typ) {
                case(tri_hp::text):
                    fout << base.idprefix << " " << mytype << std::endl;
                    break;
                default:
                    break;
            }
            return;
        }
        /** This is to read solution data **/
        virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0) {
            std::string idin,mytypein;

            switch(typ) {
                case(tri_hp::text):
                    fin >> idin >> mytypein;
                    break;
                default:
                    break;
            }
            return;
        }
        
                
        /* BOUNDARY CONDITION FUNCTIONS */
        virtual void vdirichlet() {}
        virtual void vdirichlet2d() {} //!< SPECIAL CASE OF POINT BOUNDARY CONDITION FOR 2D FIELD
        virtual void vmatchsolution_snd(int phase, FLT *vdata, int vrtstride=1) {base.vloadbuff(boundary::all,vdata,0,x.NV-1,x.NV*vrtstride);}
        virtual void vmatchsolution_rcv(int phase, FLT *vdata, int vrtstride=1) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,vdata,0,x.NV-1,x.NV*vrtstride);}
                
        /* FOR COUPLED DYNAMIC BOUNDARIES */
        virtual block::ctrl setup_preconditioner(block::ctrl ctrl_message) {return(block::stop);}
        virtual block::ctrl tadvance(block::ctrl ctrl_message) {return(block::stop);}
        virtual void calculate_unsteady_sources(bool coarse) {}
        virtual block::ctrl rsdl(block::ctrl ctrl_message) {return(block::stop);}
        virtual block::ctrl update(block::ctrl ctrl_message) {return(block::stop);}
        virtual block::ctrl mg_getfres(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh, int bnum) {return(block::stop);} 
        virtual block::ctrl mg_getcchng(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh, int bnum) {return(block::stop);}    
};


class hp_side_bdry : public sgeometry_interface {
    public:
        std::string mytype;
        tri_hp& x;
        side_bdry &base;
        const hp_side_bdry *adapt_storage;
        bool curved, coupled;
        Array<TinyVector<FLT,mesh::ND>,2> crv;
        TinyVector<Array<TinyVector<FLT,mesh::ND>,2>,sim::nhist+1> crvbd;
        Array<TinyMatrix<FLT,mesh::ND,MXGP>,2> dxdt;
        int excpt,excpt1;

    public:
        hp_side_bdry(tri_hp& xin, side_bdry &bin) : x(xin), base(bin), curved(false), coupled(false) {mytype = "plain";}
        hp_side_bdry(const hp_side_bdry &inbdry, tri_hp& xin, side_bdry &bin) : mytype(inbdry.mytype), x(xin), base(bin), curved(inbdry.curved), coupled(inbdry.coupled) {}
        virtual hp_side_bdry* create(tri_hp& xin, side_bdry &bin) const {return(new hp_side_bdry(*this,xin,bin));}
        virtual void init(input_map& input,void* &gbl_in); 
        virtual void copy(const hp_side_bdry& tgt);
        virtual ~hp_side_bdry() {}
                
        /* input output functions */
        virtual void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);
        /** This is to read solution data **/
        virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0); 
        
        /* CURVATURE FUNCTIONS */
        bool is_curved() {return(curved);}
        FLT& crds(int ind, int mode, int dir) {return(crv(ind,mode)(dir));}
        FLT& crdsbd(int tlvl, int ind, int mode, int dir) {return(crvbd(tlvl)(ind,mode)(dir));}
        void curv_init(int tlvl = 0);
                        
        /* BOUNDARY CONDITION FUNCTIONS */
        virtual block::ctrl rsdl(block::ctrl ctrl_message) {return(block::stop);}
        virtual void maxres() {}
        virtual void vdirichlet() {}
        virtual void sdirichlet(int mode) {}
        virtual void vmatchsolution_snd(int phase, FLT *vdata, int vrtstride=1) {base.vloadbuff(boundary::all,vdata,0,x.NV-1,x.NV*vrtstride);}
        virtual void vmatchsolution_rcv(int phase, FLT *vdata, int vrtstride=1) {base.vfinalrcv(boundary::all_phased,phase,boundary::symmetric,boundary::average,vdata,0,x.NV-1,x.NV*vrtstride);}
        virtual void smatchsolution_snd(FLT *sdata, int bgnmode, int endmode, int modestride) {
            base.sloadbuff(boundary::all,sdata,bgnmode*x.NV,(endmode+1)*x.NV-1,x.NV*modestride);
            return;
        }
        virtual void smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride);
    
                
        /* FOR COUPLED DYNAMIC BOUNDARIES */
        virtual block::ctrl setup_preconditioner(block::ctrl ctrl_message) {return(block::stop);}
        virtual block::ctrl tadvance(bool coarse, block::ctrl ctrl_message);
        virtual void calculate_unsteady_sources(bool coarse);
        virtual block::ctrl update(block::ctrl ctrl_message) {return(block::stop);}
        virtual block::ctrl mg_getfres(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh, int bnum) {return(block::stop);} 
        virtual block::ctrl mg_getcchng(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh, int bnum) {return(block::stop);}    
        
        /* ADAPTATION FUNCTIONS */
        virtual void updatevdata_bdry(int bel,int endpt,hp_side_bdry *bin) {}
        virtual void movevdata_bdry(int bel,int endpt,hp_side_bdry *bin) {}
        virtual void updatesdata_bdry(int bel,hp_side_bdry *bin) {}
        virtual void movesdata_bdry(int bel,hp_side_bdry *tgt) {
            int sind,tgtel,step,m,n;
            
            if (!curved || !x.sm0) return;
    
            sind = base.el(bel);
            tgtel = tgt->x.getbdryel(tgt->x.sd(sind).tri(1));
                        
            for(step=0;step<sim::nadapt;++step) {
                for(m=0;m<x.sm0;++m) {
                    for(n=0;n<x.ND;++n) {
                        crdsbd(step,bel,m,n) = tgt->crdsbd(step,tgtel,m,n);
                    }
                }
            }
            return;
        }
        
        /* SEARCH FUNCTIONS */
        virtual void findandmovebdrypt(TinyVector<FLT,2>& xp,int &bel,FLT &psi) const;
        virtual void mvpttobdry(int nel, FLT psi, TinyVector<FLT,mesh::ND> &pt);
        
        /* SOME UTILITIES */
        block::ctrl findmax(block::ctrl ctrl_message, FLT (*fxy)(TinyVector<FLT,2> &x));
        void findintercept(FLT (*fxy)(TinyVector<FLT,2> &x));
};

#endif
