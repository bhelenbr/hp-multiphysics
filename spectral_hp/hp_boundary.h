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

class hp_vrtx_bdry {
   protected:
      std::string mytype;
      tri_hp& x;
      vrtx_bdry& base;
      static FLT dummy;

   public:
      hp_vrtx_bdry(tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin) {mytype = "plain";}
      hp_vrtx_bdry(const hp_vrtx_bdry &inbdry,tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin), mytype(inbdry.mytype) {}
      virtual hp_vrtx_bdry* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_vrtx_bdry(*this,xin,bin);}
      virtual void init(input_map& input) {} /**< This is to read definition data only (not solution data) */
      virtual void copy_data(const hp_vrtx_bdry& tgt) {}
      virtual ~hp_vrtx_bdry() {}
      
      /* input output functions */
      virtual void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
         switch(typ) {
            case(text):
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
            case(text):
               fin >> idin >> mytypein;
               break;
            default:
               break;
         }
         return;
      }
      
            
      /* BOUNDARY CONDITION FUNCTIONS */
      virtual void vdirichlet() {}
      virtual void vmatchsolution_snd(int phase, FLT *vdata) {base.vloadbuff(vdata,0,x.NV-1,x.NV);}
      virtual void vmatchsolution_rcv(int phase, FLT *vdata) {base.vfinalrcv(phase,vdata,0,x.NV-1,x.NV);}
            
      /* FOR COUPLED DYNAMIC BOUNDARIES */
      virtual block::ctrl setup_preconditioner(int excpt) {return(block::stop);}
      virtual block::ctrl tadvance(int excpt) {return(block::stop);}
      virtual block::ctrl rsdl(int excpt) {return(block::stop);}
      virtual block::ctrl update(int excpt) {return(block::stop);}
      virtual block::ctrl minvrt(int excpt) {return(block::stop);}
      virtual void mg_getfres(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {} 
      virtual void mg_getcchng(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {}   
};

class hp_side_bdry {
   protected:
      std::string mytype;
      tri_hp& x;
      side_bdry &base;
      bool curved, coupled;
      Array<TinyVector<FLT,mesh::ND>,2> crv;
      TinyVector<Array<TinyVector<FLT,mesh::ND>,2>,sim::nhist+1> crvbd;
      Array<TinyMatrix<FLT,mesh::ND,MXGP>,2> dxdt;

   public:
      hp_side_bdry(tri_hp& xin, side_bdry &bin) : x(xin), base(bin), curved(false), coupled(false) {mytype = "plain";}
      hp_side_bdry(const hp_side_bdry &inbdry, tri_hp& xin, side_bdry &bin) : mytype(inbdry.mytype), x(xin), base(bin), curved(inbdry.curved), coupled(inbdry.coupled) {}
      virtual hp_side_bdry* create(tri_hp& xin, side_bdry &bin) const {return(new hp_side_bdry(*this,xin,bin));}
      virtual void init(input_map& input); 
      virtual void copy_data(const hp_side_bdry& tgt);
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
      virtual void addbflux() {}
      virtual void vdirichlet() {}
      virtual void sdirichlet(int mode) {}
      void vmatchsolution_snd(int phase, FLT *vdata) {base.vloadbuff(vdata,0,x.NV-1,x.NV);}
      void vmatchsolution_rcv(int phase, FLT *vdata) {base.vfinalrcv(phase,vdata,0,x.NV-1,x.NV);}
      void smatchsolution_snd(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) {
         base.sloadbuff(sdata,bgnmode*x.NV,(endmode+1)*x.NV-1,x.NV*modestride);
         return;
      }
      void smatchsolution_rcv(int phi, FLT *sdata, int bgn, int end, int stride);
   
            
      /* FOR COUPLED DYNAMIC BOUNDARIES */
      virtual block::ctrl setup_preconditioner(int excpt) {return(block::stop);}
      virtual block::ctrl tadvance(int excpt);
      virtual block::ctrl rsdl(int excpt) {return(block::stop);}
      virtual block::ctrl update(int excpt) {return(block::stop);}
      virtual block::ctrl minvrt(int excpt) {return(block::stop);}
      virtual void mg_getfres(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {} 
      virtual void mg_getcchng(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {}   
      
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
      virtual void findbdrypt(const TinyVector<FLT,2> xp,int &bel,FLT &psi);            

        
      /* ENPOINT COMMUNICATION ROUTINES */
      virtual void findendpts(class mesh& x) {}
      virtual void msgload_vrtx(int endpt,FLT *base,int bgn,int end,int stride) {}
      virtual void msgpass_vrtx(int endpt) {}
      virtual void msgwait_rcv_vrtx(int endpt,FLT *base,int bgn,int end,int stride) {}
};
#endif