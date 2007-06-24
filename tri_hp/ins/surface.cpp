#include "bdry_ins.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>

//#define CTRL_DEBUG
//#define MPDEBUG
//#define DEBUG

#define BODYFORCE

using namespace bdry_ins;

// extern FLT body[ND];

void surface_slave::init(input_map& inmap,void* &gbl_in) {
   std::string keyword;
   
   keyword = base.idprefix + "_curved";
   inmap[keyword] = "1";
   
   keyword = base.idprefix + "_coupled";
   inmap[keyword] = "1";

   neumann::init(inmap,gbl_in);
   
   return;
}
         
void surface::init(input_map& inmap,void* &gbl_in) {
   std::string keyword,val;
   std::istringstream data;
   std::string filename;
   bool adapt_storage;

   surface_slave::init(inmap,gbl_in);
   
   keyword = x.idprefix + "_adapt_storage";
   inmap.getwdefault(keyword,adapt_storage,false);
   if (adapt_storage) return;
   
   ksprg.resize(base.maxel);
   if (!x.coarse) {
      vdres.resize(x.log2pmax,base.maxel);
      sdres.resize(x.log2pmax,base.maxel,x.sm0);
   
      surf_gbl = new gbl;      
      gbl_in = surf_gbl;
      
      keyword = base.idprefix + "_sigma";
      inmap.getwdefault(keyword,surf_gbl->sigma,0.0);
      
      keyword = base.idprefix + "_matching_block";
      if (!inmap.get(keyword,val)) {
         surf_gbl->mu2 = 0.0;
         surf_gbl->rho2 = 0.0;
      }
      else {
         keyword = val +"_mu";
         if (!inmap.get(keyword,surf_gbl->mu2)) {
            *sim::log << "couldn't find matching blocks viscosity" << std::endl;
            exit(1);
         }
         
         keyword = val +"_rho";
         if (!inmap.get(keyword,surf_gbl->rho2)) {
            *sim::log << "couldn't find matching blocks density" << std::endl;
            exit(1);
         }
      }
      
        
      if (x.sd(base.el(0)).vrtx(0) == x.sd(base.el(base.nel-1)).vrtx(1)) surf_gbl->is_loop = true;
      else surf_gbl->is_loop = false;
   
      surf_gbl->vug0.resize(base.maxel+1);
      surf_gbl->sug0.resize(base.maxel,x.sm0);
      
      surf_gbl->vres.resize(base.maxel+1);
      surf_gbl->sres.resize(base.maxel,x.sm0); 
      surf_gbl->vres0.resize(base.maxel+1);
      surf_gbl->sres0.resize(base.maxel,x.sm0); 
      
#ifdef DROP      
      surf_gbl->vvolumeflux.resize(base.maxel+1);
      surf_gbl->svolumeflux.resize(base.maxel,x.sm0);
#endif

           
      surf_gbl->vdt.resize(base.maxel+1);
      surf_gbl->sdt.resize(base.maxel);   
      
      surf_gbl->meshc.resize(base.maxel);
      
      /* Multigrid Storage all except highest order (log2p+1)*/
      vdres.resize(x.log2p+1,base.maxel+1);
      sdres.resize(x.log2p+1,base.maxel,x.sm0);
      
      keyword = base.idprefix + "_fadd";
      inmap.getlinewdefault(keyword,val,"1.0 1.0");
      data.str(val);
      data >> surf_gbl->fadd(0) >> surf_gbl->fadd(1);  
      data.clear(); 
      
      double CFLtdflt[3] = {2.5, 1.5, 1.0};
      inmap.getwdefault(base.idprefix + "_cfltangent",&surf_gbl->cfl(0,0),3,CFLtdflt); 

      double CFLndflt[3] = {2.0, 1.25, 0.75};
      inmap.getwdefault(base.idprefix + "_cflnormal",&surf_gbl->cfl(1,0),3,CFLndflt); 
       
      surf_gbl->adis = 1.0;  /* TEMPORARY */
   }
   else {
      surf_gbl = reinterpret_cast<gbl *>(gbl_in);
      vug_frst.resize(base.maxel+1);
      vdres.resize(1,base.maxel+1);
   }
   return;
}

block::ctrl surface::tadvance(bool coarse, block::ctrl ctrl_message) {
   int i,j,m,sind;
   block::ctrl state;
   
   
   if (ctrl_message == block::begin) {
      excpt = 0;
   }
   
   switch(excpt) {
      case 0: {
         if (ctrl_message != block::advance1) {
            state = hp_side_bdry::tadvance(coarse,ctrl_message);
   
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         else 
            ++excpt;
      }
         
      case 1: {
         if (sim::substep == 0) {
            /* SET SPRING CONSTANTS */
            for(j=0;j<base.nel;++j) {
               sind = base.el(j);
               ksprg(j) = 1.0/x.distance(x.sd(sind).vrtx(0),x.sd(sind).vrtx(1));
            }
            
            /* CALCULATE TANGENT SOURCE TERM FOR FINE MESH */
            /* ZERO TANGENTIAL MESH MOVEMENT SOURCE */   
            if (!coarse) {
               for(i=0;i<base.nel+1;++i)
                  vdres(x.log2p,i)(0) = 0.0;

               for(i=0;i<base.nel;++i) 
                  for(m=0;m<basis::tri(x.log2p).sm;++m)
                     sdres(x.log2p,i,m)(0) = 0.0;
                     
               rsdl(block::begin);

               for(i=0;i<base.nel+1;++i)
                  vdres(x.log2p,i)(0) = -surf_gbl->vres(i)(0);

               for(i=0;i<base.nel;++i) 
                  for(m=0;m<basis::tri(x.log2p).sm;++m)
                     sdres(x.log2p,i,m)(0) = -0.0*surf_gbl->sres(i,m)(0); // TEMPO FOR UNIFORM SPACING OF HIGHER MODES ALONG SIDE
            }
         }
         ++excpt;
      }
   }
   
   return(block::stop);
}

block::ctrl surface::rsdl(block::ctrl ctrl_message) {
   int i,j,m,n,sind,indx,count,v0,v1;
   TinyVector<FLT,mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV);
   FLT jcb;
   Array<TinyVector<FLT,MXGP>,1> u(x.NV);
   TinyMatrix<FLT,mesh::ND,MXGP> crd, dcrd, mvel;
   TinyMatrix<FLT,8,MXGP> res;
   TinyMatrix<FLT,4,MXGP> lf,lf1;
   block::ctrl state;   
   
   if (ctrl_message == block::begin) {
      excpt = 0;
   }
   
   switch(excpt) {
      case(0): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In surface::rsdl step 0 with ctrl_message " << ctrl_message << std::endl;
#endif
         /**************************************************/
         /* DETERMINE MESH RESIDUALS & SURFACE TENSION     */
         /**************************************************/
         for(n=0;n<mesh::ND;++n)
            surf_gbl->vres(0)(n) = 0.0;
            
#ifdef DROP
         surf_gbl->vvolumeflux(0) = 0.0;
#endif

          for(indx=0;indx<base.nel;++indx) {
            sind = base.el(indx);
            v0 = x.sd(sind).vrtx(0);
            v1 = x.sd(sind).vrtx(1);
            
            x.crdtocht1d(sind);
            for(n=0;n<mesh::ND;++n)
               basis::tri(x.log2p).proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
               
            x.ugtouht1d(sind);
            for(n=0;n<mesh::ND;++n)
               basis::tri(x.log2p).proj1d(&x.uht(n)(0),&u(n)(0));   
            
            for(i=0;i<basis::tri(x.log2p).gpx;++i) {
               norm(0) =  dcrd(1,i);
               norm(1) = -dcrd(0,i);
               jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
               
               /* RELATIVE VELOCITY STORED IN MVEL(N)*/
               for(n=0;n<mesh::ND;++n) {
                  mvel(n,i) = u(n)(i) -(sim::bd[0]*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef DROP
                  mvel(n,i) -= tri_hp_ins::mesh_ref_vel(n);
#endif   
               }
               /* TANGENTIAL SPACING */            
               res(0,i) = -ksprg(indx)*jcb;
               /* NORMAL FLUX */
               res(1,i) = -RAD(crd(0,i))*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1));    
               /* UPWINDING BASED ON TANGENTIAL VELOCITY */
               res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*surf_gbl->meshc(indx);
               
#ifdef DROP
               res(3,i) = +RAD(crd(0,i))*surf_gbl->vflux*jcb;
#endif 
          
               /* SURFACE TENSION SOURCE TERM X-DIRECTION */ 
               res(4,i) = +RAD(crd(0,i))*(x.gbl_ptr->rho -surf_gbl->rho2)*sim::g*crd(1,i)*norm(0);
#ifdef AXISYMMETRIC
               res(4,i) += surf_gbl->sigma*jcb;
#endif
               /* AND INTEGRATION BY PARTS TERM */
               res(5,i) = +RAD(crd(0,i))*surf_gbl->sigma*norm(1)/jcb;


               /* SURFACE TENSION SOURCE TERM Y-DIRECTION */
               res(6,i) = +RAD(crd(0,i))*(x.gbl_ptr->rho -surf_gbl->rho2)*sim::g*crd(1,i)*norm(1);            
               /* AND INTEGRATION BY PARTS TERM */
               res(7,i) = -RAD(crd(0,i))*surf_gbl->sigma*norm(0)/jcb;
            }
				
            for(m=0;m<basis::tri(x.log2p).sm+2;++m)
               lf(0,m) = 0.0;

            /* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */               
            basis::tri(x.log2p).intgrtx1d(&lf(0,0),&res(0,0));
            basis::tri(x.log2p).intgrt1d(&lf(1,0),&res(1,0));
            basis::tri(x.log2p).intgrtx1d(&lf(1,0),&res(2,0));
#ifdef DROP
            basis::tri(x.log2p).intgrt1d(&lf(2,0),&res(3,0));
#endif
            
            /* TO LEAVE TANGENTIAL POSITION TOTALLY FREE */
            /* for(m=0;m<basis::tri(x.log2p).sm+2;++m)
               cf(0,m) = 0.0; */
            
            /* STORE IN RES */
            for(n=0;n<mesh::ND;++n) {
               surf_gbl->vres(indx)(n) += lf(n,0);
               surf_gbl->vres(indx+1)(n) = lf(n,1);
               for(m=0;m<basis::tri(x.log2p).sm;++m)
                  surf_gbl->sres(indx,m)(n) = lf(n,m+2);
            }
                        
#ifdef DROP
            surf_gbl->vres(indx)(1) += lf(2,0);
            surf_gbl->vres(indx+1)(1) += lf(2,1);
            for(m=0;m<basis::tri(x.log2p).sm;++m)
               surf_gbl->sres(indx,m)(1) += lf(2,m+2);      
               
            surf_gbl->vvolumeflux(indx) += lf(2,0);
            surf_gbl->vvolumeflux(indx+1) = lf(2,1);
            for(m=0;m<basis::tri(x.log2p).sm;++m)
               surf_gbl->svolumeflux(indx,m) = lf(2,m+2);               
#endif
            
            /* INTEGRATE & STORE SURFACE TENSION SOURCE TERM */
            basis::tri(x.log2p).intgrt1d(&lf1(0,0),&res(4,0));
            basis::tri(x.log2p).intgrtx1d(&lf1(0,0),&res(5,0));
            basis::tri(x.log2p).intgrt1d(&lf1(1,0),&res(6,0));
            basis::tri(x.log2p).intgrtx1d(&lf1(1,0),&res(7,0));
				
				/* MASS FLUX PRECONDITIONER */				
				for(n=2;n<x.NV-1;++n) /* For swirling case */
					for(m=0;m<basis::tri(x.log2p).sm+2;++m)
						lf1(n,m) = 0.0;

            for(m=0;m<basis::tri(x.log2p).sm+2;++m)
               lf1(x.NV-1,m) = -x.gbl_ptr->rho*lf(1,m); 

#ifndef INERTIALESS
            for (n=0;n<x.NV-1;++n) 
               ubar(n) = 0.5*(x.uht(n)(0) +x.uht(n)(1));
					               
            for (n=0;n<x.NV-1;++n) {
               lf1(n,0) -= x.uht(n)(0)*(x.gbl_ptr->rho -surf_gbl->rho2)*lf(1,0);
               lf1(n,1) -= x.uht(n)(1)*(x.gbl_ptr->rho -surf_gbl->rho2)*lf(1,1);
               for(m=0;m<basis::tri(x.log2p).sm;++m)
                  lf1(n,m+2) -= ubar(n)*(x.gbl_ptr->rho -surf_gbl->rho2)*lf(1,m+2);
            }
#endif
            /* ADD TO RESIDUAL */
            for(n=0;n<x.NV;++n)
               x.gbl_ptr->res.v(v0,n) += lf1(n,0);

            for(n=0;n<x.NV;++n)
               x.gbl_ptr->res.v(v1,n) += lf1(n,1);
            
            for(m=0;m<basis::tri(x.log2p).sm;++m) {
               for(n=0;n<x.NV;++n)
                  x.gbl_ptr->res.s(sind,m,n) += lf1(n,m+2);
            }  
         }

          /* CALL VERTEX RESIDUAL HERE */
         for(i=0;i<2;++i)
            state &= x.hp_vbdry(base.vbdry(i))->rsdl(block::begin);
         
         /************************************************/
         /* MODIFY SURFACE RESIDUALS ON COARSER MESHES   */
         /************************************************/   
         if(x.coarse) {
            if (x.isfrst) {
               for(i=0;i<base.nel+1;++i) 
                  for(n=0;n<mesh::ND;++n)
                     vdres(x.log2p,i)(n) = surf_gbl->fadd(n)*surf_gbl->vres0(i)(n) -surf_gbl->vres(i)(n);
               
               for(i=0;i<base.nel;++i) 
                  for(m=0;m<basis::tri(x.log2p).sm;++m) 
                     for(n=0;n<mesh::ND;++n)
                        sdres(x.log2p,i,m)(n) = surf_gbl->fadd(n)*surf_gbl->sres0(i,m)(n) -surf_gbl->sres(i,m)(n);
               
            }
            for(i=0;i<base.nel+1;++i) 
               for(n=0;n<mesh::ND;++n)
                  surf_gbl->vres(i)(n) += vdres(x.log2p,i)(n);
            
            for(i=0;i<base.nel;++i) 
               for(m=0;m<basis::tri(x.log2p).sm;++m) 
                  for(n=0;n<mesh::ND;++n)
                     surf_gbl->sres(i,m)(n) += sdres(x.log2p,i,m)(n);
         }
         else {
            /* ADD TANGENTIAL MESH MOVEMENT SOURCE */   
            for(i=0;i<base.nel+1;++i)
               surf_gbl->vres(i)(0) += vdres(x.log2p,i)(0);

            for(i=0;i<base.nel;++i)
               for(m=0;m<basis::tri(x.log2p).sm;++m) 
                  surf_gbl->sres(i,m)(0) += sdres(x.log2p,i,m)(0);
         }

         if (base.is_comm()) { 
            count = 0;
            for(j=0;j<base.nel+1;++j) {
               base.fsndbuf(count++) = surf_gbl->vres(j)(1)*surf_gbl->rho2;
#ifdef MPDEBUG 
               *sim::log << surf_gbl->vres(j)(1)*surf_gbl->rho2 << '\n';
#endif
#ifdef DROP
               base.fsndbuf(count-1) -= surf_gbl->vvolumeflux(j)*surf_gbl->rho2;
#endif
            }
            for(j=0;j<base.nel;++j) {
               for(m=0;m<basis::tri(x.log2p).sm;++m) {
                  base.fsndbuf(count++) = surf_gbl->sres(j,m)(1)*surf_gbl->rho2;
#ifdef DROP
                  base.fsndbuf(count-1) -= surf_gbl->svolumeflux(j,m)*surf_gbl->rho2;
#endif
               }
            }
            base.sndsize() = count;
            base.sndtype() = boundary::flt_msg;
            base.comm_prepare(boundary::all,0,boundary::master_slave);
         }         
         ++excpt;
         return(block::advance);
      }
      case(1): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In surface::rsdl step 1 with ctrl_message " << ctrl_message << std::endl;
#endif
         base.comm_exchange(boundary::all,0,boundary::master_slave);
         ++excpt;
         return(block::advance);
      }
      case(2): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In surface::rsdl step 2 with ctrl_message " << ctrl_message << std::endl;
#endif
         base.comm_wait(boundary::all,0,boundary::master_slave);
         ++excpt;
         return(block::advance);
      }

   }
   return(block::stop);
}

block::ctrl surface_slave::rsdl(block::ctrl ctrl_message) {
   int i,m,msgn,count,sind,v0;
   
   if (ctrl_message == block::begin) {
      excpt = 0;
   }
   
   switch(excpt) {
      case(0): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 0 of surface_slave::rsdl with ctrl_message: " << ctrl_message << std::endl;
#endif
         base.comm_prepare(boundary::all,0,boundary::master_slave);
         ++excpt;
         return(block::advance);
      }
      case(1): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 1 of surface_slave::rsdl with ctrl_message: " << ctrl_message << std::endl;
#endif
         base.comm_exchange(boundary::all,0,boundary::master_slave);
         ++excpt;
         return(block::advance);
      }
      case(2): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 2 of surface_slave::rsdl with ctrl_message: " << ctrl_message << std::endl;
#endif
         base.comm_wait(boundary::all,0,boundary::master_slave);        
         count = 0;
         for(i=base.nel-1;i>=0;--i) {
            sind = base.el(i);
            v0 = x.sd(sind).vrtx(1);
#ifdef MPDEBUG
            *sim::log << x.gbl_ptr->res.v(v0,x.NV-1) << ' ' << base.frcvbuf(0,count) << '\n';
#endif
            x.gbl_ptr->res.v(v0,x.NV-1) += base.frcvbuf(0,count++);
         }
         v0 = x.sd(sind).vrtx(0);
#ifdef MPDEBUG
         *sim::log << x.gbl_ptr->res.v(v0,x.NV-1) << ' ' << base.frcvbuf(0,count) << '\n';
#endif    
         x.gbl_ptr->res.v(v0,x.NV-1) += base.frcvbuf(0,count++);
    
         for(i=base.nel-1;i>=0;--i) {
            sind = base.el(i);
            msgn = 1;
            for(m=0;m<basis::tri(x.log2p).sm;++m) {
               x.gbl_ptr->res.s(sind,m,x.NV-1) += msgn*base.frcvbuf(0,count++);
               msgn *= -1;
            }
         }

         ++excpt;
      }
   }
   return(block::stop);
}

block::ctrl surface_outflow_endpt::rsdl(block::ctrl ctrl_message) {
   int bnum,sind;
   TinyVector<FLT,mesh::ND> ubar, tangent, rp;
   FLT jcb;
   
   if (ctrl_message == block::begin) {
      /* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
      /* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
      bnum = base.sbdry(surfbdry);
      if (surfbdry == 0) {
         sind = x.sbdry(bnum)->el(x.sbdry(bnum)->nel-1);
         x.crdtocht1d(sind);
         basis::tri(x.log2p).ptprobe1d(2,&rp(0),&tangent(0),1.0,&x.cht(0,0),MXTM);
         jcb = sqrt(tangent(0)*tangent(0) +tangent(1)*tangent(1));
         x.gbl_ptr->res.v(base.v0,0) += -RAD(rp(0))*surf->surf_gbl->sigma*tangent(0)/jcb;
         x.gbl_ptr->res.v(base.v0,1) += -RAD(rp(0))*surf->surf_gbl->sigma*tangent(1)/jcb;
      }
      else {
         sind = x.sbdry(bnum)->el(0);
         x.crdtocht1d(sind);
         basis::tri(x.log2p).ptprobe1d(2,&rp(0),&tangent(0),-1.0,&x.cht(0,0),MXTM);
         jcb = sqrt(tangent(0)*tangent(0) +tangent(1)*tangent(1));
         x.gbl_ptr->res.v(base.v0,0) -= -RAD(rp(0))*surf->surf_gbl->sigma*tangent(0)/jcb;
         x.gbl_ptr->res.v(base.v0,1) -= -RAD(rp(0))*surf->surf_gbl->sigma*tangent(1)/jcb;
      }
   }
   return(block::stop);
}

block::ctrl surface::minvrt(block::ctrl ctrl_message) {
   int i,m,n,indx;
   FLT temp;
   
   if (ctrl_message == block::begin) excpt = 0;
   else excpt += ctrl_message;
   
   switch(excpt) {
      case(0): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 0 of surface::minvrt with ctrl_message: " << ctrl_message << std::endl;
#endif
         /* INVERT MASS MATRIX */
         /* LOOP THROUGH SIDES */
         if (basis::tri(x.log2p).sm > 0) {
            for(indx = 0; indx<base.nel; ++indx) {
               /* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */         
               for (m=0; m <basis::tri(x.log2p).sm; ++m) {
                  for(n=0;n<mesh::ND;++n)
                     surf_gbl->vres(indx)(n) -= basis::tri(x.log2p).sfmv1d(0,m)*surf_gbl->sres(indx,m)(n);
                  for(n=0;n<mesh::ND;++n)
                     surf_gbl->vres(indx+1)(n) -= basis::tri(x.log2p).sfmv1d(1,m)*surf_gbl->sres(indx,m)(n);
               }
            }
         }
         mp_phase = -1;
         ++excpt;
         ctrl_message = block::stay;
      }
      case(1): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 1 of surface::minvrt with ctrl_message: " << ctrl_message << std::endl;
#endif
         ++mp_phase;
         switch(mp_phase%3) {
            case(0):
               x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&surf_gbl->vres(0)(0),0,1,0);
               x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&surf_gbl->vres(base.nel)(0),0,1,0);
               x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase/3,boundary::symmetric);
               x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase/3,boundary::symmetric);
               return(block::stay);
            case(1):
               x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase/3,boundary::symmetric);
               x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase/3,boundary::symmetric);
               return(block::stay);
            case(2):
               i = x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase/3,boundary::symmetric);
               i &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase/3,boundary::symmetric);
               x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase/3,boundary::symmetric,boundary::average,&surf_gbl->vres(0)(0),0,1,0);
               x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase/3,boundary::symmetric,boundary::average,&surf_gbl->vres(base.nel)(0),0,1,0);
               return(static_cast<block::ctrl>(i));
         }
      }
      
      case(2): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 2 of surface::minvrt with ctrl_message: " << ctrl_message << std::endl;
#endif
         if (surf_gbl->is_loop) {
            surf_gbl->vres(0)(1) = 0.5*(surf_gbl->vres(0)(1) +surf_gbl->vres(base.nel+1)(1));
            surf_gbl->vres(base.nel+1)(1) = surf_gbl->vres(0)(1);
            surf_gbl->vres(0)(0) = 0.0;
            surf_gbl->vres(base.nel+1)(0) = 0.0;
         }
         x.hp_vbdry(base.vbdry(0))->vdirichlet();
         x.hp_vbdry(base.vbdry(1))->vdirichlet();
         

         /* SOLVE FOR VERTEX MODES */
         for(i=0;i<base.nel+1;++i) {
            temp                = surf_gbl->vres(i)(0)*surf_gbl->vdt(i)(0,0) +surf_gbl->vres(i)(1)*surf_gbl->vdt(i)(0,1);
            surf_gbl->vres(i)(1) = surf_gbl->vres(i)(0)*surf_gbl->vdt(i)(1,0) +surf_gbl->vres(i)(1)*surf_gbl->vdt(i)(1,1);
            surf_gbl->vres(i)(0) = temp;
         }
         
         /* SOLVE FOR SIDE MODES */
         if (basis::tri(x.log2p).sm > 0) {
            for(indx = 0; indx<base.nel; ++indx) {
               
               /* INVERT SIDE MODES */
               DPBTRSNU2(&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,&(surf_gbl->sres(indx,0)(0)),mesh::ND);
               for(m=0;m<basis::tri(x.log2p).sm;++m) {
                  temp                      = surf_gbl->sres(indx,m)(0)*surf_gbl->sdt(indx)(0,0) +surf_gbl->sres(indx,m)(1)*surf_gbl->sdt(indx)(0,1);
                  surf_gbl->sres(indx,m)(1) = surf_gbl->sres(indx,m)(0)*surf_gbl->sdt(indx)(1,0) +surf_gbl->sres(indx,m)(1)*surf_gbl->sdt(indx)(1,1);       
                  surf_gbl->sres(indx,m)(0) = temp;
               }

               for(m=0;m<basis::tri(x.log2p).sm;++m) {
                  for(n=0;n<mesh::ND;++n) {
                     surf_gbl->sres(indx,m)(n) -= basis::tri(x.log2p).vfms1d(0,m)*surf_gbl->vres(indx)(n);
                     surf_gbl->sres(indx,m)(n) -= basis::tri(x.log2p).vfms1d(1,m)*surf_gbl->vres(indx+1)(n);
                  }
               }
            }
         }
      }
   }
 
   return(block::stop);
}

block::ctrl surface::setup_preconditioner(block::ctrl ctrl_message) {
   int i,indx,m,n,sind,v0,v1;
   TinyVector<FLT,mesh::ND> nrm;
   FLT h, hsm;
   FLT dttang, dtnorm;
   FLT uvel, vvel, vslp, strss;
   FLT drho, srho, smu;
   FLT nu1, nu2;
   FLT qmax, gam1, gam2;
   
   if (ctrl_message == block::begin) excpt = 0;
   else excpt += ctrl_message;

   drho = x.gbl_ptr->rho -surf_gbl->rho2;
   srho = x.gbl_ptr->rho +surf_gbl->rho2;
   smu = x.gbl_ptr->mu +surf_gbl->mu2;
   nu1 = x.gbl_ptr->mu/x.gbl_ptr->rho;
   if (surf_gbl->rho2 > 0.0) 
      nu2 = surf_gbl->mu2/surf_gbl->rho2;
   else
      nu2 = 0.0;

   switch(excpt) {
      case(0): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In surface::setup_preconditioner step 0 with ctrl_message " << ctrl_message << std::endl;
#endif
      
         /**************************************************/
         /* DETERMINE SURFACE MOVEMENT TIME STEP           */
         /**************************************************/
         surf_gbl->vdt(0) = 0.0;
                  
         for(indx=0; indx < base.nel; ++indx) {
            sind = base.el(indx);
            v0 = x.sd(sind).vrtx(0);
            v1 = x.sd(sind).vrtx(1);

            nrm(0) =  (x.vrtx(v1)(1) -x.vrtx(v0)(1));
            nrm(1) = -(x.vrtx(v1)(0) -x.vrtx(v0)(0));
            h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
            
            uvel = x.ug.v(v0,0)-(sim::bd[0]*(x.vrtx(v0)(0) -x.vrtxbd(1)(v0)(0)));
            vvel = x.ug.v(v0,1)-(sim::bd[0]*(x.vrtx(v0)(1) -x.vrtxbd(1)(v0)(1)));
#ifdef DROP
            uvel -= tri_hp_ins::mesh_ref_vel(0);
            vvel -= tri_hp_ins::mesh_ref_vel(1);
#endif

            qmax = uvel*uvel+vvel*vvel;
            vslp = fabs(-uvel*nrm(1)/h +vvel*nrm(0)/h);

            uvel = x.ug.v(v1,0)-(sim::bd[0]*(x.vrtx(v1)(0) -x.vrtxbd(1)(v1)(0)));
            vvel = x.ug.v(v1,1)-(sim::bd[0]*(x.vrtx(v1)(1) -x.vrtxbd(1)(v1)(1)));
#ifdef DROP
            uvel -= tri_hp_ins::mesh_ref_vel(0);
            vvel -= tri_hp_ins::mesh_ref_vel(1);
#endif
            qmax = MAX(qmax,uvel*uvel+vvel*vvel);
            vslp = MAX(vslp,fabs(-uvel*nrm(1)/h +vvel*nrm(0)/h));
            
            hsm = h/(.25*(basis::tri(x.log2p).p+1)*(basis::tri(x.log2p).p+1));
                  
            dttang = 2.*ksprg(indx)*(.25*(basis::tri(x.log2p).p+1)*(basis::tri(x.log2p).p+1))/hsm;
#ifndef BODYFORCE
            strss =  4.*surf_gbl->sigma/(hsm*hsm) +fabs(drho*sim::g*nrm(1)/h);
#else
            strss =  4.*surf_gbl->sigma/(hsm*hsm) +fabs(drho*(-sim::body(0)*nrm(0) +(sim::g-sim::body(1))*nrm(1))/h);
#endif

            gam1 = 3.0*qmax +(0.5*hsm*sim::bd[0] + 2.*nu1/hsm)*(0.5*hsm*sim::bd[0] + 2.*nu1/hsm);
            gam2 = 3.0*qmax +(0.5*hsm*sim::bd[0] + 2.*nu2/hsm)*(0.5*hsm*sim::bd[0] + 2.*nu2/hsm);

            if (sim::bd[0] + x.gbl_ptr->mu == 0.0) gam1 = MAX(gam1,0.1);

#ifdef INERTIALESS
            gam1 = (2.*nu1/hsm)*(2.*nu1/hsm);
            gam2 = (2.*nu2/hsm)*(2.*nu2/hsm);
#endif
            dtnorm = 2.*vslp/hsm +sim::bd[0] +1.*strss/(x.gbl_ptr->rho*sqrt(qmax +gam1) +surf_gbl->rho2*sqrt(qmax +gam2));            
            
            /* SET UP DISSIPATIVE COEFFICIENT */
            /* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
            /* RESIDUAL HAS DX/2 WEIGHTING */
            /* |a| dx/2 dv/dx  dx/2 dpsi */
            /* |a| dx/2 2/dx dv/dpsi  dpsi */
            /* |a| dv/dpsi  dpsi */
            // surf_gbl->meshc(indx) = surf_gbl->adis/(h*dtnorm*0.5);
            surf_gbl->meshc(indx) = surf_gbl->adis/(h*(vslp/hsm +sim::bd[0]));
            
            dtnorm *= RAD(0.5*(x.vrtx(v0)(0) +x.vrtx(v1)(0)));
            
            nrm *= 0.5;
                  
            surf_gbl->vdt(indx)(0,0) += -dttang*nrm(1)*basis::tri(x.log2p).vdiag1d;
            surf_gbl->vdt(indx)(0,1) +=  dttang*nrm(0)*basis::tri(x.log2p).vdiag1d;
            surf_gbl->vdt(indx)(1,0) +=  dtnorm*nrm(0)*basis::tri(x.log2p).vdiag1d;
            surf_gbl->vdt(indx)(1,1) +=  dtnorm*nrm(1)*basis::tri(x.log2p).vdiag1d;
            surf_gbl->vdt(indx+1)(0,0) = -dttang*nrm(1)*basis::tri(x.log2p).vdiag1d;
            surf_gbl->vdt(indx+1)(0,1) =  dttang*nrm(0)*basis::tri(x.log2p).vdiag1d;
            surf_gbl->vdt(indx+1)(1,0) =  dtnorm*nrm(0)*basis::tri(x.log2p).vdiag1d;
            surf_gbl->vdt(indx+1)(1,1) =  dtnorm*nrm(1)*basis::tri(x.log2p).vdiag1d;
                     
            if (basis::tri(x.log2p).sm) {
               surf_gbl->sdt(indx)(0,0) = -dttang*nrm(1);
               surf_gbl->sdt(indx)(0,1) =  dttang*nrm(0);
               surf_gbl->sdt(indx)(1,0) =  dtnorm*nrm(0);
               surf_gbl->sdt(indx)(1,1) =  dtnorm*nrm(1);
            }
         }
         mp_phase = -1;
         ++excpt;
         ctrl_message = block::stay;
      }
      case(1): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In surface::setup_preconditioner step 1 with ctrl_message " << ctrl_message << std::endl;
#endif
         ++mp_phase;
         switch(mp_phase%3) {
            case(0):
               x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&surf_gbl->vdt(0)(0,0),0,3,0);
               x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&surf_gbl->vdt(base.nel)(0,0),0,3,0);
               x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase/3,boundary::symmetric);
               x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase/3,boundary::symmetric);
               return(block::stay);
            case(1):
               x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase/3,boundary::symmetric);
               x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase/3,boundary::symmetric);
               return(block::stay);
            case(2):
               i = x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase/3,boundary::symmetric);
               i &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase/3,boundary::symmetric);
               x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase/3,boundary::symmetric,boundary::average,&surf_gbl->vdt(0)(0,0),0,3,0);
               x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase/3,boundary::symmetric,boundary::average,&surf_gbl->vdt(base.nel)(0,0),0,3,0);
               return(static_cast<block::ctrl>(i));
         }
      }
      
      case(2): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In surface::setup_preconditioner step 2 with ctrl_message " << ctrl_message << std::endl;
#endif
         if (surf_gbl->is_loop) {
            for(m=0;m<mesh::ND;++m)
               for(n=0;n<mesh::ND;++n)
                  surf_gbl->vdt(0)(m,n) = 0.5*(surf_gbl->vdt(0)(m,n) +surf_gbl->vdt(base.nel+1)(m,n));
            surf_gbl->vdt(base.nel+1) = surf_gbl->vdt(0);
         }

         FLT jcbi,temp;
         for(indx=0;indx<base.nel+1;++indx) {   
            /* INVERT VERTEX MATRIX */
            jcbi = 1.0/(surf_gbl->vdt(indx)(0,0)*surf_gbl->vdt(indx)(1,1) 
                       -surf_gbl->vdt(indx)(0,1)*surf_gbl->vdt(indx)(1,0));
                       
            temp = surf_gbl->vdt(indx)(0,0)*jcbi*surf_gbl->cfl(1,x.log2p);
            surf_gbl->vdt(indx)(0,0) = surf_gbl->vdt(indx)(1,1)*jcbi*surf_gbl->cfl(0,x.log2p);
            surf_gbl->vdt(indx)(1,1) = temp;
            surf_gbl->vdt(indx)(0,1) *= -jcbi*surf_gbl->cfl(1,x.log2p);
            surf_gbl->vdt(indx)(1,0) *= -jcbi*surf_gbl->cfl(0,x.log2p);
         }

         /* INVERT SIDE MATRIX */   
         if (basis::tri(x.log2p).sm > 0) {
            for(indx=0;indx<base.nel;++indx) {
               /* INVERT SIDE MVDT MATRIX */
               jcbi = 1.0/(surf_gbl->sdt(indx)(0,0)*surf_gbl->sdt(indx)(1,1)
                          -surf_gbl->sdt(indx)(0,1)*surf_gbl->sdt(indx)(1,0));

               temp = surf_gbl->sdt(indx)(0,0)*jcbi*surf_gbl->cfl(1,x.log2p);
               surf_gbl->sdt(indx)(0,0) = surf_gbl->sdt(indx)(1,1)*jcbi*surf_gbl->cfl(0,x.log2p);
               surf_gbl->sdt(indx)(1,1) = temp;
               surf_gbl->sdt(indx)(0,1) *= -jcbi*surf_gbl->cfl(1,x.log2p);
               surf_gbl->sdt(indx)(1,0) *= -jcbi*surf_gbl->cfl(0,x.log2p);
            }
         }
      }
   }
   
   return(block::stop);
}

block::ctrl surface::update(block::ctrl ctrl_message) {
   int i,m,n,count,sind,indx,v0;
   block::ctrl state;

   if (ctrl_message == block::begin) excpt1 = 0;

   switch(excpt1) {
      case(0): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 0 of surface::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         indx = 0;
         for(i=0;i<base.nel;++i) {
            sind = base.el(i);
            v0 = x.sd(sind).vrtx(0);
            surf_gbl->vug0(i) = x.vrtx(v0);
         }
         v0 = x.sd(sind).vrtx(1);
         surf_gbl->vug0(base.nel) = x.vrtx(v0);
            
         if (basis::tri(x.log2p).sm > 0) surf_gbl->sug0(Range(0,base.nel-1),Range(0,basis::tri(x.log2p).sm-1)) = crv(Range(0,base.nel-1),Range(0,basis::tri(x.log2p).sm-1));

         ++excpt1;
         stage = 0;
         return(block::stop);
      }
      
      case(1): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 1 of surface::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         if (ctrl_message == block::advance) {
            ++excpt1;
            ctrl_message = block::begin;
         }
         else return(block::stop);
      }
      
      case(2): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 2 of surface::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         if (ctrl_message != block::advance1) {
            state = minvrt(ctrl_message);
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         ++excpt1;
         ctrl_message = block::begin;
      }
      
      case(3): {       
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 3 of surface::update with ctrl_message: " << ctrl_message << std::endl;
#endif

#ifdef DEBUG
      for(i=0;i<base.nel+1;++i)
         printf("vdt: %d %8.4e %8.4e %8.4e %8.4e\n",i,surf_gbl->vdt(i)(0,0),surf_gbl->vdt(i)(0,1),surf_gbl->vdt(i)(1,0),surf_gbl->vdt(i)(1,1));
            
      for(i=0;i<base.nel;++i)
         printf("sdt: %d %8.4e %8.4e %8.4e %8.4e\n",i,surf_gbl->sdt(i)(0,0),surf_gbl->sdt(i)(0,1),surf_gbl->sdt(i)(1,0),surf_gbl->sdt(i)(1,1));

      for(i=0;i<base.nel+1;++i) {
         printf("vres: %d ",i);
         for(n=0;n<mesh::ND;++n) {
            if (fabs(surf_gbl->vres(i)(n)) > 1.0e-9) printf("%8.4e ",surf_gbl->vres(i)(n));
            else printf("%8.4e ",0.0);
         }
         printf("\n");
      }
         
      for(i=0;i<base.nel;++i) {
         for(m=0;m<basis::tri(x.log2p).sm;++m) {
            printf("sres: %d ",i);
            for(n=0;n<mesh::ND;++n) {
               if (fabs(surf_gbl->sres(i,m)(n)) > 1.0e-9) printf("%8.4e ",surf_gbl->sres(i,m)(n));
               else printf("%8.4e ",0.0);
            }
            printf("\n");
         }
      }
                  
      for(i=0;i<base.nel;++i) {
         sind = base.el(i);
         v0 = x.sd(sind).vrtx(0);
         printf("vertex positions %d %8.4e %8.4e\n",v0,x.vrtx(v0)(0),x.vrtx(v0)(1));
      }
      v0 = x.sd(sind).vrtx(1);
      printf("vertex positions %d %8.4e %8.4e\n",v0,x.vrtx(v0)(0),x.vrtx(v0)(1));
      
      for(i=0;i<base.nel;++i)
         for(m=0;m<basis::tri(x.log2p).sm;++m)
            printf("spos: %d %d %8.4e %8.4e\n",i,m,crv(i,m)(0),crv(i,m)(1));

#endif
         
         for(i=0;i<base.nel;++i) {
            sind = base.el(i);
            v0 = x.sd(sind).vrtx(0);
            x.vrtx(v0) = surf_gbl->vug0(i) -sim::alpha[stage]*surf_gbl->vres(i);
         }
         v0 = x.sd(sind).vrtx(1);
         x.vrtx(v0) = surf_gbl->vug0(base.nel) -sim::alpha[stage]*surf_gbl->vres(base.nel);
         
         if (basis::tri(x.log2p).sm > 0) crv(Range(0,base.nel-1),Range(0,basis::tri(x.log2p).sm-1)) = surf_gbl->sug0(Range(0,base.nel-1),Range(0,basis::tri(x.log2p).sm-1)) -sim::alpha[stage]*surf_gbl->sres(Range(0,base.nel-1),Range(0,basis::tri(x.log2p).sm-1));
      
         /* FIX POINTS THAT SLIDE ON CURVE */
         x.hp_vbdry(base.vbdry(1))->mvpttobdry(x.vrtx(v0));
         sind = base.el(0);
         v0 = x.sd(sind).vrtx(0);
         x.hp_vbdry(base.vbdry(0))->mvpttobdry(x.vrtx(v0));
         
#ifdef DEBUG
         for(i=0;i<base.nel;++i) {
            sind = base.el(i);
            v0 = x.sd(sind).vrtx(0);
            printf("vertex positions %d %e %e\n",v0,x.vrtx(v0)(0),x.vrtx(v0)(1));
         }
         v0 = x.sd(sind).vrtx(1);
         printf("vertex positions %d %e %e\n",v0,x.vrtx(v0)(0),x.vrtx(v0)(1));
         
         for(i=0;i<base.nel;++i)
            for(m=0;m<basis::tri(x.log2p).sm;++m)
               printf("spos: %d %d %e %e\n",i,m,crv(i,m)(0),crv(i,m)(1));
#endif
         
         if (base.is_comm()) {            
            count = 0;
            for(i=0;i<base.nel;++i) {
               sind = base.el(i);
               v0 = x.sd(sind).vrtx(0);
               for(n=0;n<mesh::ND;++n)
                  base.fsndbuf(count++) = x.vrtx(v0)(n);
            }
            v0 = x.sd(sind).vrtx(1);
            for(n=0;n<mesh::ND;++n)
               base.fsndbuf(count++) = x.vrtx(v0)(n);            
            
            for(i=0;i<base.nel;++i) {
               for(m=0;m<basis::tri(x.log2p).sm;++m) {
                  for(n=0;n<mesh::ND;++n)
                     base.fsndbuf(count++) = crv(i,m)(n);
               }
            }
            base.sndsize() = count;
            base.sndtype() = boundary::flt_msg;
            base.comm_prepare(boundary::all,0,boundary::master_slave);
         }
         ++excpt1;
         return(block::advance);
      }
      
      case(4): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 4 of surface::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         base.comm_exchange(boundary::all,0,boundary::master_slave);
         ++excpt1;
         return(block::advance);
      }
      
      case(5): {
         base.comm_wait(boundary::all,0,boundary::master_slave);
         ++stage;
         ++excpt1;
         if (stage < sim::NSTAGE) excpt1 = 1;
         return(block::stop);
      }
   }
   
   return(block::stop);
}

block::ctrl surface_slave::update(block::ctrl ctrl_message) {
   int i,m,n,msgn,count,sind,v0;

   if (ctrl_message == block::begin) excpt1 = 0;

   switch(excpt1) {
      case(0): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 0 of surface_slave::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         ++excpt1;
         stage = 0;
         return(block::stop);
      }
      
      case(1): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 1 of surface_slave::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         if (ctrl_message == block::advance) {
            ++excpt1;
            ctrl_message = block::begin;
         }
         else return(block::stop);
      }
      
      case(2): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 2 of surface_slave::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         if (ctrl_message != block::advance1) {
            return(block::advance1);
         }
         ++excpt1;
         ctrl_message = block::begin;
      }
      
      case(3): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 3 of surface_slave::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         base.comm_prepare(boundary::all,0,boundary::master_slave);
         ++excpt1;
         return(block::advance);
      }
      case(4): {
#ifdef CTRL_DEBUG
         *sim::log << base.idprefix << " In step 4 of surface_slave::update with ctrl_message: " << ctrl_message << std::endl;
#endif
         base.comm_exchange(boundary::all,0,boundary::master_slave);
         ++excpt1;
         return(block::advance);
      }
      case(5): {

         base.comm_wait(boundary::all,0,boundary::master_slave);
         
         count = 0;
         for(i=base.nel-1;i>=0;--i) {
            sind = base.el(i);
            v0 = x.sd(sind).vrtx(1);
            for(n=0;n<mesh::ND;++n) 
               x.vrtx(v0)(n) = base.frcvbuf(0,count++);
         }
         v0 = x.sd(sind).vrtx(0);
         for(n=0;n<mesh::ND;++n)
            x.vrtx(v0)(n) = base.frcvbuf(0,count++);
         
         if (basis::tri(x.log2p).sm > 0) {
            for(i=base.nel-1;i>=0;--i) {
               sind = base.el(i);
               msgn = 1;
               for(m=0;m<basis::tri(x.log2p).sm;++m) {
                  for(n=0;n<mesh::ND;++n)
                     crds(i,m,n) = msgn*base.frcvbuf(0,count++);
                  msgn *= -1;
               }
            }
         }
         ++stage;
         ++excpt1;
         if (stage < sim::NSTAGE) excpt1 = 1;
         return(block::stop);
      }
   }
   return(block::stop);
}

block::ctrl surface::mg_getfres(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh, int bnum) {
   int i,indx,tind,v0,snum,sind;
   surface* fbdry = dynamic_cast<surface *>(fmesh->hp_sbdry(bnum));

   if (ctrl_message == block::begin) excpt = 0;
         
   switch(excpt) {
      case(0): {
         if(x.p0 > 1) {
            /* TRANSFER IS ON FINEST MESH */
            surf_gbl->vres0(Range(0,base.nel)) = surf_gbl->vres(Range(0,base.nel));
            if (basis::tri(x.log2p).sm > 0) surf_gbl->sres0(Range(0,base.nel-1),Range(0,basis::tri(x.log2p).sm-1)) = surf_gbl->sres(Range(0,base.nel-1),Range(0,basis::tri(x.log2p).sm-1));
            return(block::stop);
         }
         else {
            /* TRANSFER IS BETWEEN DIFFERENT MESHES */
            surf_gbl->vres0(Range(0,base.nel)) = 0.0;
               
            /* CALCULATE COARSE RESIDUALS */
            /* DO ENDPOINTS FIRST */
            surf_gbl->vres0(0) = surf_gbl->vres(0);
            surf_gbl->vres0(base.nel) = surf_gbl->vres(fbdry->base.nel);
               
            for(i=1;i<fbdry->base.nel;++i) {
               sind = fbdry->base.el(i);
               v0 = fmesh->sd(sind).vrtx(0);
               tind = fv_to_ct(v0).tri;
               for(snum=0;snum<3;++snum) 
                  if (x.getbdrynum(x.td(tind).tri(snum))  == bnum) break;
               assert(snum != 3);
               indx = x.getbdryel(x.td(tind).tri(snum));
               surf_gbl->vres0(indx) += fv_to_ct(v0).wt((snum+1)%3)*surf_gbl->vres(i);
               surf_gbl->vres0(indx+1) += fv_to_ct(v0).wt((snum+2)%3)*surf_gbl->vres(i);
            }
         }
         /* MESH POSITIONS AND VRTX_FIRST STORAGE WILL BE HANDLED BY R_MESH */
         ++excpt;
      }
   }
   
   return(block::stop);
}

void surface_slave::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
   int j,k,n,countup,offset;
   
   if (!base.is_comm()) return;
   
#ifdef MPDEBUG
      *sim::log << base.idprefix << " In surface_snd"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
   
   countup = 0;
   for(j=0;j<base.nel;++j) {
      offset = base.el(j)*stride*x.NV;
      for(k=bgn;k<=end;++k) {
         for(n=0;n<x.NV-1;++n) {
            base.fsndbuf(countup++) = sdata[offset +k*x.NV +n];
#ifdef MPDEBUG
               *sim::log << "\t" << sdata[offset +k*x.NV +n] << std::endl;
#endif
         }
      }
   }
   base.sndsize() = countup;
   base.sndtype() = boundary::flt_msg;
   return;
}


void surface_slave::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
   
   if (!base.is_comm()) return;
      
   /* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
   int j,k,m,n,count,countdn,countup,offset,sind,sign;
   FLT mtchinv;
   
   /* ASSUMES REVERSE ORDERING OF SIDES */
   /* WON'T WORK IN 3D */
   
   int matches = 1;
   
   int bgnsign = (bgn % 2 ? -1 : 1);
   
   /* RELOAD FROM BUFFER */
   /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
   /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */   
   for(m=0;m<base.matches();++m) {         
      ++matches;
      
      int ebp1 = end-bgn+1;

      countdn = (base.nel-1)*ebp1*(x.NV-1);
      countup = 0;
      for(j=0;j<base.nel;++j) {
         sign = bgnsign;
         for(k=0;k<ebp1;++k) {
            for(n=0;n<x.NV-1;++n) {
               base.fsndbuf(countup++) += sign*base.frcvbuf(m,countdn++);
            }
            sign *= -1;
         }
         countdn -= 2*ebp1*(x.NV-1);
      }
   }
   
   if (matches > 1) {
      mtchinv = 1./matches;

#ifdef MPDEBUG
      *sim::log << base.idprefix << " In surface_rcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
      count = 0;
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         offset = sind*stride*x.NV;
         for (k=bgn;k<=end;++k) {
            for(n=0;n<x.NV-1;++n) {
               sdata[offset +k*x.NV +n] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
               *sim::log << "\t" << sdata[offset +k*x.NV +n] << std::endl;
#endif
            }
         }

      }
   }
   return;
}
   

void surface::maxres() {
   int i,n;
   TinyVector<FLT,mesh::ND> mxr;

   mxr = 0.0;
   
   for(i=0;i<base.nel+1;++i)
      for(n=0;n<mesh::ND;++n)
         mxr(n) = MAX(fabs(surf_gbl->vres(i)(n)),mxr(n));

   for(n=0;n<mesh::ND;++n)
      *sim::log << ' ' << mxr(n) << ' ';
   
   return;
}
