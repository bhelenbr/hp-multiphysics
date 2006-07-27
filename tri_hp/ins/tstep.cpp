#include "tri_hp_ins.h"
#include "hp_boundary.h"
#include <math.h>

#ifndef MATRIX_PRECONDITIONER
block::ctrl tri_hp_ins::setup_preconditioner(block::ctrl ctrl_message) {
   int tind,i,j,side,v0;
   FLT jcb,h,hmax,q,qmax,lam1,gam;
   TinyVector<int,3> v;
   
   if (ctrl_message == block::begin) excpt = 0;
   ctrl_message = tri_hp::setup_preconditioner(ctrl_message);
   
   if (ctrl_message == block::advance1) ++excpt;
   
   if (excpt == 3) {
      ++excpt;
   
      FLT nu = gbl_ptr->mu/gbl_ptr->rho;
      
      /***************************************/
      /** DETERMINE FLOW PSEUDO-TIME STEP ****/
      /***************************************/
      gbl_ptr->vprcn(Range(0,nvrtx-1),Range::all()) = 0.0;
      if (basis::tri(log2p).sm > 0) {
         gbl_ptr->sprcn(Range(0,nside-1),Range::all()) = 0.0;
      }
      
#ifdef TIMEACCURATE
      gam = 10.0;
      FLT dtstari = 0.0;
#endif

      for(tind = 0; tind < ntri; ++tind) {
         jcb = 0.25*area(tind);  // area is 2 x triangle area
         v = td(tind).vrtx;
         hmax = 0.0;
         for(j=0;j<3;++j) {
            h = pow(vrtx(v(j))(0) -vrtx(v((j+1)%3))(0),2.0) + 
            pow(vrtx(v(j))(1) -vrtx(v((j+1)%3))(1),2.0);
            hmax = (h > hmax ? h : hmax);
         }
         hmax = sqrt(hmax);
         
         if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
            *sim::log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
            *sim::log << "approximate location: " << vrtx(v(0))(0) << ' ' << vrtx(v(0))(1) << std::endl;
            mesh::output("negative",grid);
            exit(1);
         }
         h = 4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
         hmax = hmax/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1));
      
         qmax = 0.0;
         for(j=0;j<3;++j) {
            v0 = v(j);
#ifdef DROP
            q = pow(ug.v(v0,0)-0.5*(sim::bd[0]*(vrtx(v0)(0) -vrtxbd(1)(v0)(0)) +mesh_ref_vel(0)),2.0) 
               +pow(ug.v(v0,1)-0.5*(sim::bd[0]*(vrtx(v0)(1) -vrtxbd(1)(v0)(1)) +mesh_ref_vel(1)),2.0);
#else
            q = pow(ug.v(v0,0)-0.5*(sim::bd[0]*(vrtx(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
               +pow(ug.v(v0,1)-0.5*(sim::bd[0]*(vrtx(v0)(1) -vrtxbd(1)(v0)(1))),2.0);  
#endif
            qmax = MAX(qmax,q);
         }

#ifndef TIMEACCURATE
         gam = 3.0*qmax +(0.5*hmax*sim::bd[0] +2.*nu/hmax)*(0.5*hmax*sim::bd[0] +2.*nu/hmax);
         if (gbl_ptr->mu + sim::bd[0] == 0.0) gam = MAX(gam,0.1);
#endif
         q = sqrt(qmax);
         lam1 = q + sqrt(qmax +gam);
         
         /* SET UP DISSIPATIVE COEFFICIENTS */
         gbl_ptr->tau(tind) = adis*h/(jcb*sqrt(gam));
         gbl_ptr->delt(tind) = qmax*gbl_ptr->tau(tind);
         
         /* SET UP DIAGONAL PRECONDITIONER */
         // jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +sim::bd[0];
         jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned


#ifdef INERTIALESS
         gam = pow(2.*nu/hmax,2); 
         lam1 = sqrt(gam);
         
         /* SET UP DISSIPATIVE COEFFICIENTS */
         gbl_ptr->tau(tind)  = adis*h/(jcb*sqrt(gam));
         gbl_ptr->delt(tind) = 0.0;

         jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax;
#endif

#ifdef TIMEACCURATE
         dtstari = MAX((nu/(h*h) +lam1/h +sim::bd[0]),dtstari);

      }
      printf("#iterative to physical time step ratio: %f\n",sim::bd[0]/dtstari);
         
      for(tind=0;tind<ntri;++tind) {
         v = td(tind).vrtx;
         jcb = 0.25*area(tind)*dtstari;
#endif

         jcb *= RAD((vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.);

         gbl_ptr->tprcn(tind,0) = gbl_ptr->rho*jcb;   
         gbl_ptr->tprcn(tind,1) = gbl_ptr->rho*jcb;     
         gbl_ptr->tprcn(tind,2) =  jcb/gam;
         for(i=0;i<3;++i) {
            gbl_ptr->vprcn(v(i),Range::all())  += gbl_ptr->tprcn(tind,Range::all());
            if (basis::tri(log2p).sm > 0) {
               side = td(tind).side(i);
               gbl_ptr->sprcn(side,Range::all()) += gbl_ptr->tprcn(tind,Range::all());
            }
         }
      }
   }
   
   return(ctrl_message);
}

#else
void hp_mgrid::tstep1(void) {
    int tind,i,j,n,m,sind,count,bnum,side,v0,v1,*v;
    FLT jcb,h,hmax,q,qmax,lam1,c,u[ND],dd,gam,dtstari;
   class mesh *tgt;

   /***************************************/
   /** DETERMINE FLOW PSEUDO-TIME STEP ****/
   /***************************************/
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         for(m=0;m<NV;++m)
            gbl_ptr->vprcn[i][n][m] = 0.0;

   if (basis::tri(log2p).sm > 0) {
      for(i=0;i<nside;++i)
         for(n=0;n<NV;++n)
            for(m=0;m<NV;++m)
               gbl_ptr->sprcn[i][n][m] = 0.0;
   }

   for(tind = 0; tind < ntri; ++tind) {
      jcb = 0.25*area(tind);

      /* CALCULATE SOME MEAN / MAX QUANTITIES ON TRIANGLE */      
      v = tvrtx[tind];
      hmax = 0.0;
      qmax = 0.0;
      u[0] = 0.0;
      u[1] = 0.0;

      v1 = v[2];
      for(j=0;j<3;++j) {
         v0 = v1;
         v1 = v[j];
         dd = vrtx(v1)(0) -vrtx(v0)(0);
         h = dd*dd;
         dd = vrtx(v1)(1) -vrtx(v0)(1);
         h += dd*dd;
         hmax = (h > hmax ? h : hmax);

         dd = ug.v(v1,0)-0.5*(sim::bd[0]*vrtx(v1)(0) +dvrtdt[v1][0]);
         q = dd*dd;
         dd = ug.v(v1,1)-0.5*(sim::bd[0]*vrtx(v1)(1) +dvrtdt[v1][1]);
         q += dd*dd;
         u[0] += ug.v(v1,0);
         u[1] += ug.v(v1,1);
         dd = ug.v(v1,0)-(sim::bd[0]*vrtx(v1)(0) +dvrtdt[v1][0]);
         q = dd*dd;
         dd = ug.v(v1,1)-(sim::bd[0]*vrtx(v1)(1) +dvrtdt[v1][1]);
         q += dd*dd;
         qmax = MAX(qmax,q);
      }
      hmax = sqrt(hmax);
      h = 2.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
      u[0] *= 1./3.;
      u[1] *= 1./3.;

      /* THIS IS GAMMA (DIAGONAL PRECONDITIONER FOR CONTINUITY) */
      gam = qmax +(0.25*h*sim::bd[0] + nu/h)*(0.25*h*sim::bd[0] +nu/h); 
      q = sqrt(qmax);
      c = sqrt(qmax+gam);
      lam1  = (q+c);

      /* SET UP DISSIPATIVE COEFFICIENTS */
      gbl_ptr->tau(tind)  = adis*h/(jcb*sqrt(gam));
      gbl_ptr->delt(tind) = qmax*gbl_ptr->tau(tind);

      /* STORE PRECONDITIONER (THIS IS TO DRIVE ITERATION USING NONCONSERVATIVE SYSTEM) */
      dtstari = jcb*(nu/(h*h) +lam1/h +sim::bd[0])*RAD((vrtx[v[0]][0] +vrtx[v[1]][0] +vrtx[v[2]][0])/3.);

      gbl_ptr->tprcn[tind][0][0]  = dtstari*gbl_ptr->rho;
   	/* gbl_ptr->tprcn[tind][1][1] = gbl_ptr->tprcn[tind][0][0]; */
      gbl_ptr->tprcn[tind][NV-1][NV-1] =  dtstari/gam;
      gbl_ptr->tprcn[tind][0][NV-1] = u[0]/gam*dtstari;
      gbl_ptr->tprcn[tind][1][NV-1] = u[1]/gam*dtstari;
      
      for(i=0;i<3;++i) {
         /* ASSEMBLE VERTEX PRECONDITIONER */
         gbl_ptr->vprcn[v[i]][0][0]  += gbl_ptr->tprcn[tind][0][0];
         /* gbl_ptr->vprcn[v[i]][1][1]  += gbl_ptr->tprcn[tind][0][0];  */
         gbl_ptr->vprcn[v[i]][NV-1][NV-1]  += gbl_ptr->tprcn[tind][NV-1][NV-1];
         gbl_ptr->vprcn[v[i]][0][NV-1]  += gbl_ptr->tprcn[tind][0][NV-1];
         gbl_ptr->vprcn[v[i]][1][NV-1]  += gbl_ptr->tprcn[tind][1][NV-1];

         if (basis::tri(log2p).sm > 0) {
            /* ASSEMBLE SIDE PRECONDITIONER */
            side = td(tind).side(i);
            gbl_ptr->sprcn[side][0][0]  += gbl_ptr->tprcn[tind][0][0];
            /* gbl_ptr->sprcn[side][1][1]  += gbl_ptr->tprcn[tind][0][0];  */
            gbl_ptr->sprcn[side][NV-1][NV-1]  += gbl_ptr->tprcn[tind][NV-1][NV-1];
            gbl_ptr->sprcn[side][0][NV-1]  += gbl_ptr->tprcn[tind][0][NV-1];
            gbl_ptr->sprcn[side][1][NV-1]  += gbl_ptr->tprcn[tind][1][NV-1];
         }
      }
   }
   
   /* SEND Y-DIRECTION BOUNDARY INFORMATION */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][0];
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][NV-1];
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][1][NV-1];
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][NV-1][NV-1];
         }
         v0 = sd(sind).vrtx(1);
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][0];
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][NV-1];
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][1][NV-1];
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][NV-1][NV-1];

         /* SEND SIDE INFO */
         if (basis::tri(log2p).sm) {
            for(j=0;j<sbdry(i)->nel;++j) {
               sind = sbdry(i)->el(j);
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][0][0];
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][0][NV-1];
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][1][NV-1];
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][NV-1][NV-1];
            }
         }
      }
      
      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][0];
         }
         v0 = sd(sind).vrtx(1);
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][0];

         /* SEND SIDE INFO */
         if (basis::tri(log2p).sm) {
            for(j=0;j<sbdry(i)->nel;++j) {
               sind = sbdry(i)->el(j);
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][0][0];
            }
         }
      }
   }
   
   /* SET UP TSTEP FOR ACTIVE BOUNDARIES */   
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfdt1(i);
   
   return;
}

void hp_mgrid::tstep_mp() {
   int i,j,sind,v0,count,bnum;
   class mesh *tgt;
   
   /* RECEIVE COMMUNCATION PACKETS FROM OTHER MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            gbl_ptr->vprcn[v0][0][0] = 0.5*(gbl_ptr->vprcn[v0][0][0] +sbuff[i][count++]);
            gbl_ptr->vprcn[v0][0][NV-1] = 0.5*(gbl_ptr->vprcn[v0][0][NV-1] +sbuff[i][count++]);
            gbl_ptr->vprcn[v0][1][NV-1] = 0.5*(gbl_ptr->vprcn[v0][1][NV-1] +sbuff[i][count++]);
            gbl_ptr->vprcn[v0][NV-1][NV-1] = 0.5*(gbl_ptr->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         gbl_ptr->vprcn[v0][0][0] = 0.5*(gbl_ptr->vprcn[v0][0][0] +sbuff[i][count++]);
         gbl_ptr->vprcn[v0][0][NV-1] = 0.5*(gbl_ptr->vprcn[v0][0][NV-1] +sbuff[i][count++]);
         gbl_ptr->vprcn[v0][1][NV-1] = 0.5*(gbl_ptr->vprcn[v0][1][NV-1] +sbuff[i][count++]);
         gbl_ptr->vprcn[v0][NV-1][NV-1] = 0.5*(gbl_ptr->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (basis::tri(log2p).sm > 0) {
            for(j=sbdry(i)->nel-1;j>=0;--j) {
               sind = sbdry(i)->el(j);
               gbl_ptr->sprcn[sind][0][0] = 0.5*(gbl_ptr->sprcn[sind][0][0] +sbuff[i][count++]);
               gbl_ptr->sprcn[sind][0][NV-1] = 0.5*(gbl_ptr->sprcn[sind][0][NV-1] +sbuff[i][count++]);
               gbl_ptr->sprcn[sind][1][NV-1] = 0.5*(gbl_ptr->sprcn[sind][1][NV-1] +sbuff[i][count++]);
               gbl_ptr->sprcn[sind][NV-1][NV-1] = 0.5*(gbl_ptr->sprcn[sind][NV-1][NV-1] +sbuff[i][count++]);
            }
         }
      }

      /* ONLY SEND & RECEIVE DIAGV'S FOR INTERFACE */
      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            gbl_ptr->vprcn[v0][0][0] = 0.5*(gbl_ptr->vprcn[v0][0][0] +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         gbl_ptr->vprcn[v0][0][0] = 0.5*(gbl_ptr->vprcn[v0][0][0] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (basis::tri(log2p).sm > 0) {
            for(j=sbdry(i)->nel-1;j>=0;--j) {
               sind = sbdry(i)->el(j);
               gbl_ptr->sprcn[sind][0][0] = 0.5*(gbl_ptr->sprcn[sind][0][0] +sbuff[i][count++]);
            }
         }
      }
   }
   
   /* SEND X-DIRECTION BOUNDARY INFORMATION */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][0];
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][NV-1];
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][1][NV-1];
            tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][NV-1][NV-1];
         }
         v0 = sd(sind).vrtx(1);
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][0];
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][0][NV-1];
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][1][NV-1];
         tgt->sbuff[bnum][count++] = gbl_ptr->vprcn[v0][NV-1][NV-1];

         /* SEND SIDE INFO */
         if (basis::tri(log2p).sm) {
            for(j=0;j<sbdry(i)->nel;++j) {
               sind = sbdry(i)->el(j);
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][0][0];
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][0][NV-1];
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][1][NV-1];
               tgt->sbuff[bnum][count++] = gbl_ptr->sprcn[sind][NV-1][NV-1];
            }
         }
      }
   }
   
   return;
}

void hp_mgrid::tstep2(void) {
   int i,j,sind,v0,count;
   
   /* RECEIVE COMMUNCATION PACKETS FROM OTHER MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            gbl_ptr->vprcn[v0][0][0] = 0.5*(gbl_ptr->vprcn[v0][0][0] +sbuff[i][count++]);
            gbl_ptr->vprcn[v0][0][NV-1] = 0.5*(gbl_ptr->vprcn[v0][0][NV-1] +sbuff[i][count++]);
            gbl_ptr->vprcn[v0][1][NV-1] = 0.5*(gbl_ptr->vprcn[v0][1][NV-1] +sbuff[i][count++]);
            gbl_ptr->vprcn[v0][NV-1][NV-1] = 0.5*(gbl_ptr->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         gbl_ptr->vprcn[v0][0][0] = 0.5*(gbl_ptr->vprcn[v0][0][0] +sbuff[i][count++]);
         gbl_ptr->vprcn[v0][0][NV-1] = 0.5*(gbl_ptr->vprcn[v0][0][NV-1] +sbuff[i][count++]);
         gbl_ptr->vprcn[v0][1][NV-1] = 0.5*(gbl_ptr->vprcn[v0][1][NV-1] +sbuff[i][count++]);
         gbl_ptr->vprcn[v0][NV-1][NV-1] = 0.5*(gbl_ptr->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (basis::tri(log2p).sm > 0) {
            for(j=sbdry(i)->nel-1;j>=0;--j) {
               sind = sbdry(i)->el(j);
               gbl_ptr->sprcn[sind][0][0] = 0.5*(gbl_ptr->sprcn[sind][0][0] +sbuff[i][count++]);
               gbl_ptr->sprcn[sind][0][NV-1] = 0.5*(gbl_ptr->sprcn[sind][0][NV-1] +sbuff[i][count++]);
               gbl_ptr->sprcn[sind][1][NV-1] = 0.5*(gbl_ptr->sprcn[sind][1][NV-1] +sbuff[i][count++]);
               gbl_ptr->sprcn[sind][NV-1][NV-1] = 0.5*(gbl_ptr->sprcn[sind][NV-1][NV-1] +sbuff[i][count++]);
            }
         }
      }
   }
   
   /* INVERT PRECONDITIONER FOR VERTICES */
   for(i=0;i<nvrtx;++i) {
      gbl_ptr->vprcn[i][0][0]  = 1.0/(basis::tri(log2p).vdiag*gbl_ptr->vprcn[i][0][0]);
      gbl_ptr->vprcn[i][NV-1][NV-1]  = 1.0/(basis::tri(log2p).vdiag*gbl_ptr->vprcn[i][NV-1][NV-1]);
      gbl_ptr->vprcn[i][0][NV-1] = -basis::tri(log2p).vdiag*gbl_ptr->vprcn[i][0][NV-1]*gbl_ptr->vprcn[i][0][0]*gbl_ptr->vprcn[i][NV-1][NV-1];
      gbl_ptr->vprcn[i][1][NV-1] = -basis::tri(log2p).vdiag*gbl_ptr->vprcn[i][1][NV-1]*gbl_ptr->vprcn[i][0][0]*gbl_ptr->vprcn[i][NV-1][NV-1];      
   }
   
   if (basis::tri(log2p).sm > 0) {
      /* INVERT PRECONDITIONER FOR SIDES */            
      for(i=0;i<nside;++i) {
         gbl_ptr->sprcn[i][0][0] = 1.0/gbl_ptr->sprcn[i][0][0];
         gbl_ptr->sprcn[i][NV-1][NV-1] = 1.0/gbl_ptr->sprcn[i][NV-1][NV-1];
         gbl_ptr->sprcn[i][0][NV-1] = -gbl_ptr->sprcn[i][0][NV-1]*gbl_ptr->sprcn[i][0][0]*gbl_ptr->sprcn[i][NV-1][NV-1];
         gbl_ptr->sprcn[i][1][NV-1] = -gbl_ptr->sprcn[i][1][NV-1]*gbl_ptr->sprcn[i][0][0]*gbl_ptr->sprcn[i][NV-1][NV-1];
      }
   }
   
   /* SET UP TSTEP FOR ACTIVE BOUNDARIES */   
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfdt2(i);
   
   return;
}
#endif
