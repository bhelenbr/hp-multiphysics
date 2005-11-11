/*
 *  hp_boundary.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 9/3/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
 
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
      virtual void init(std::map <std::string,std::string>& input, std::string prefix) {}
      virtual hp_vrtx_bdry* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_vrtx_bdry(xin,bin);}
      virtual void copy_data(const hp_vrtx_bdry& tgt) {}
      virtual ~hp_vrtx_bdry() {}
            
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

      virtual void input(std::map<std::string,std::string>& inmap) {         
         int tlvl;
         std::istringstream data(inmap[base.idprefix+".hp_tlvl"]);
         if (!(data >> tlvl)) tlvl = 0;
         data.clear();
         
         std::string filename;
         bool readfile = true;
         data.str(inmap[base.idprefix+".hp_datafile"]);
         if (!(data >> filename)) readfile = false;
         data.clear();
         
         int typ;
         data.str(inmap["hp_filetype"]);
         if (!(data >> typ)) typ = tri_hp::text;
         data.clear();
         
         if (readfile) {
            ifstream fin;
            fin.open(filename.c_str());
            input(fin,static_cast<tri_hp::filetype>(typ),tlvl);
            fin.close();
         }
         
         return;
      }
      virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0) {} 
};

class hp_side_bdry {
   protected:
      std::string mytype;
      tri_hp& x;
      side_bdry &base;
      Array<TinyMatrix<FLT,mesh::ND,MXGP>,2> dxdt;
      static FLT dummy;

   public:
      hp_side_bdry(tri_hp& xin, side_bdry &bin) : x(xin), base(bin) {mytype = "plain";}
      virtual void init(std::map <std::string,std::string>& input, std::string prefix) {
         dxdt.resize(x.log2pmax,base.maxel);
      }
      virtual hp_side_bdry* create(tri_hp& xin, side_bdry &bin) const {return(new hp_side_bdry(xin,bin));}
      virtual void copy_data(const hp_side_bdry& tgt) {}
      virtual ~hp_side_bdry() {}
      
      /* CURVATURE FUNCTIONS */
      virtual bool is_curved() {return(false);}
      virtual FLT& crds(int ind, int mode, int dir) {return(dummy=0);}
      virtual FLT& crdsbd(int ind, int mode, int dir, int tlvl) {return(dummy=0);}
      virtual void curvinit(int tlvl) {}
      
      /* BOUNDARY CONDITION FUNCTIONS */
      virtual void addbflux(bool mgrid) {}
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
      virtual block::ctrl tadvance(int excpt) {return(block::stop);}
      virtual block::ctrl rsdl(int excpt) {return(block::stop);}
      virtual block::ctrl update(int excpt) {return(block::stop);}
      virtual block::ctrl minvrt(int excpt) {return(block::stop);}
      virtual void mg_getfres(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {} 
      virtual void mg_getcchng(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {}   
      
      /* ADAPTATION FUNCTIONS */
      virtual void updatevdata_bdry(int bel,int endpt,hp_side_bdry *bin) {}
      virtual void movevdata_bdry(int bel,int endpt,hp_side_bdry *bin) {}
      virtual void updatesdata_bdry(int bel,hp_side_bdry *bin) {}
      virtual void movesdata_bdry(int bel,hp_side_bdry *bin) {}
      
      /* SEARCH FUNCTIONS */
      virtual void findbdrypt(TinyVector<FLT,2> xp,int &bel,FLT &psi) {
         base.findbdrypt(xp,bel,psi);
         basis::tri(x.log2p).ptvalues1d(psi);
      }
      
      /* input output functions */
      virtual void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);
      virtual void input(std::map<std::string,std::string>& inmap);
      virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0) {}

        
      /* ENPOINT COMMUNICATION ROUTINES */
      virtual void findendpts(class mesh& x) {}
      virtual void msgload_vrtx(int endpt,FLT *base,int bgn,int end,int stride) {}
      virtual void msgpass_vrtx(int endpt) {}
      virtual void msgwait_rcv_vrtx(int endpt,FLT *base,int bgn,int end,int stride) {}
};


class hp_curved : public hp_side_bdry {
   protected:
      bool coupled;
      Array<TinyVector<FLT,mesh::ND>,2> crv;
      TinyVector<Array<TinyVector<FLT,mesh::ND>,2>,sim::nhist+1> crvbd;
   public:

      hp_curved(tri_hp& hpin, side_bdry &bin) : hp_side_bdry(hpin,bin) { 
         mytype = "curved";
      }
      hp_curved* create(tri_hp& xin, side_bdry &bin) const {return(new hp_curved(xin,bin));}
      void copy_data(const hp_side_bdry &tgt);

      /* CURVATURE FUNCTIONS */
      bool is_curved() {return(true);}
      FLT& crds(int ind, int mode, int dir) {return(crv(ind,mode)(dir));}
      FLT& crdsbd(int tlvl, int ind, int mode, int dir) {return(crvbd(tlvl)(ind,mode)(dir));}
      
      void init(std::map <std::string,std::string>& input, std::string prefix);      
      void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0);         
      virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0);      
      void curv_init(int tlvl = 0);
            
      void findbdrypt(TinyVector<FLT,2> xp,int &bel,FLT &psi);      
      block::ctrl tadvance(int excpt);
};

