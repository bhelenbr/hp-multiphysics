#include"hp_mgrid.h"
#include<math.h>
#include<utilities.h>

#ifdef CONSERV
void hp_mgrid::tstep1(void) {
   int tind,i,j,n,sind,count,bnum,side,v0,*v;
   FLT jcb,h,hmax,q,qmax,lam1,gam;
   class mesh *tgt;
   
   /***************************************/
   /** DETERMINE FLOW PSEUDO-TIME STEP ****/
   /***************************************/
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->vprcn[i][n][n] = 0.0;

   if (b.sm > 0) {
      for(i=0;i<nside;++i) 
         for(n=0;n<NV;++n)
            gbl->sprcn[i][n][n] = 0.0;
   }
   
#ifdef TIMEACCURATE
   gam = 10.0;
   FLT dtstari = 0.0;
#endif

   for(tind = 0; tind < ntri; ++tind) {
      jcb = 0.25*area(tind);
      v = tvrtx[tind];
      hmax = 0.0;
      for(j=0;j<3;++j) {
         h = pow(vrtx[v[j]][1] -vrtx[v[(j+1)%3]][1],2.0) + 
         pow(vrtx[v[j]][0] -vrtx[v[(j+1)%3]][0],2.0);
         hmax = (h > hmax ? h : hmax);
      }
      hmax = sqrt(hmax);
      
      if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
         printf("negative triangle area caught in tstep %d %d\n",nvrtx,ntri);
         output("negative",tecplot);
         out_mesh("negative",grid);
         exit(1);
      }
      // h = 2.*jcb/(0.25*(b.p +1)*(b.p+1)*hmax); THIS IS WRONG BY FACTOR OF 2
      h = 4.*jcb/(0.25*(b.p +1)*(b.p+1)*hmax);
      hmax = hmax/(0.25*(b.p +1)*(b.p+1));
      
      qmax = 0.0;
      for(j=0;j<3;++j) {
         v0 = v[j];
         q = pow(ug.v[v0][0]-0.5*(bd[0]*vrtx[v0][0] +dvrtdt[v0][0]),2.0) 
            +pow(ug.v[v0][1]-0.5*(bd[0]*vrtx[v0][1] +dvrtdt[v0][1]),2.0);
         qmax = MAX(qmax,q);
      }
#ifndef TIMEACCURATE
      gam = 3.0*qmax +(0.5*hmax*bd[0] +2.*gbl->nu/hmax)*(0.5*hmax*bd[0] +2.*gbl->nu/hmax);
      // gam = MAX(gam,0.01); // TEMPORARY FOR INVISCID ONLY
#endif

      q = sqrt(qmax);
      lam1 = q + sqrt(qmax +gam);
      
      /* SET UP DISSIPATIVE COEFFICIENTS */
      gbl->tau[tind]  = adis*h/(jcb*sqrt(gam));
      gbl->delt[tind] = qmax*gbl->tau[tind];
      
#ifdef INERTIALESS
      gam = pow(2.*gbl->nu/hmax,2); 
      lam1 = sqrt(gam);
      
      /* SET UP DISSIPATIVE COEFFICIENTS */
      gbl->tau[tind]  = adis*h/(jcb*sqrt(gam));
      gbl->delt[tind] = 0.0;
#endif
      
      /* SET UP DIAGONAL PRECONDITIONER */
#ifdef TIMEACCURATE
      dtstari = MAX((gbl->nu/(h*h) +lam1/h +bd[0]),dtstari);
#else

#ifndef INERTIALESS
      jcb *= 8.*gbl->nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +bd[0];
#else
      jcb *= 8.*gbl->nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax;
#endif
      
#ifdef AXISYMMETRIC
      jcb *= (vrtx[v[0]][0] +vrtx[v[1]][0] +vrtx[v[2]][0])/3.;
#endif
#endif

#ifdef TIMEACCURATE
   }
   printf("#iterative to physical time step ratio: %f\n",bd[0]/dtstari);
      
   for(tind=0;tind<ntri;++tind) {
      v = tvrtx[tind];
      jcb = 0.25*area(tind)*dtstari;
#ifdef AXISYMMETRIC
      jcb *= (vrtx[v[0]][0] +vrtx[v[1]][0] +vrtx[v[2]][0])/3.;
#endif
#endif
      gbl->tprcn[tind][0][0] = gbl->rho*jcb;   
      gbl->tprcn[tind][1][1] = gbl->rho*jcb;     
      gbl->tprcn[tind][NV-1][NV-1] =  jcb/gam;
      for(i=0;i<3;++i) {
         gbl->vprcn[v[i]][0][0]  += gbl->tprcn[tind][0][0];
         gbl->vprcn[v[i]][NV-1][NV-1]  += gbl->tprcn[tind][NV-1][NV-1];
         if (b.sm > 0) {
            side = tside[tind].side[i];
            gbl->sprcn[side][0][0] += gbl->tprcn[tind][0][0];
            gbl->sprcn[side][NV-1][NV-1] += gbl->tprcn[tind][NV-1][NV-1];
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
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];

         /* SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][0];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][NV-1][NV-1];
            }
         }
      }
      
      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];

         /* SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][0];
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
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
            gbl->vprcn[v0][NV-1][NV-1] = 0.5*(gbl->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
         gbl->vprcn[v0][NV-1][NV-1] = 0.5*(gbl->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (b.sm > 0) {
            for(j=sbdry[i].num-1;j>=0;--j) {
               sind = sbdry[i].el[j];
               gbl->sprcn[sind][0][0] = 0.5*(gbl->sprcn[sind][0][0] +sbuff[i][count++]);
               gbl->sprcn[sind][NV-1][NV-1] = 0.5*(gbl->sprcn[sind][NV-1][NV-1] +sbuff[i][count++]);
            }
         }
      }

      /* ONLY SEND & RECEIVE DIAGV'S FOR INTERFACE */
      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (b.sm > 0) {
            for(j=sbdry[i].num-1;j>=0;--j) {
               sind = sbdry[i].el[j];
               gbl->sprcn[sind][0][0] = 0.5*(gbl->sprcn[sind][0][0] +sbuff[i][count++]);
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
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];

         /* SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][0];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][NV-1][NV-1];
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
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
            gbl->vprcn[v0][NV-1][NV-1] = 0.5*(gbl->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
         gbl->vprcn[v0][NV-1][NV-1] = 0.5*(gbl->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (b.sm > 0) {
            for(j=sbdry[i].num-1;j>=0;--j) {
               sind = sbdry[i].el[j];
               gbl->sprcn[sind][0][0] = 0.5*(gbl->sprcn[sind][0][0] +sbuff[i][count++]);
               gbl->sprcn[sind][NV-1][NV-1] = 0.5*(gbl->sprcn[sind][NV-1][NV-1] +sbuff[i][count++]);
            }
         }
      }
   }
      
   /* FORM DIAGANOL PRECONDITIONER FOR VERTICES */
   for(i=0;i<nvrtx;++i) {
      gbl->vprcn[i][0][0] = 1.0/(b.vdiag*gbl->vprcn[i][0][0]);
      gbl->vprcn[i][1][1] = gbl->vprcn[i][0][0];
      gbl->vprcn[i][NV-1][NV-1] = 1.0/(b.vdiag*gbl->vprcn[i][NV-1][NV-1]);
   }
      
   if (b.sm > 0) {
      /* FORM DIAGANOL PRECONDITIONER FOR SIDES */            
      for(i=0;i<nside;++i) {
         gbl->sprcn[i][0][0] = 1.0/gbl->sprcn[i][0][0];
         gbl->sprcn[i][1][1] = gbl->sprcn[i][0][0];
         gbl->sprcn[i][NV-1][NV-1] = 1.0/gbl->sprcn[i][NV-1][NV-1];
      }
   }
   
   /* SET UP TSTEP FOR ACTIVE BOUNDARIES */   
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfdt2(i);
   
   return;
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
            gbl->vprcn[i][n][m] = 0.0;

   if (b.sm > 0) {
      for(i=0;i<nside;++i)
         for(n=0;n<NV;++n)
            for(m=0;m<NV;++m)
               gbl->sprcn[i][n][m] = 0.0;
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
         dd = vrtx[v1][0] -vrtx[v0][0];
         h = dd*dd;
         dd = vrtx[v1][1] -vrtx[v0][1];
         h += dd*dd;
         hmax = (h > hmax ? h : hmax);

         dd = ug.v[v1][0]-0.5*(bd[0]*vrtx[v1][0] +dvrtdt[v1][0]);
         q = dd*dd;
         dd = ug.v[v1][1]-0.5*(bd[0]*vrtx[v1][1] +dvrtdt[v1][1]);
         q += dd*dd;
         u[0] += ug.v[v1][0];
         u[1] += ug.v[v1][1];
         dd = ug.v[v1][0]-(bd[0]*vrtx[v1][0] +dvrtdt[v1][0]);
         q = dd*dd;
         dd = ug.v[v1][1]-(bd[0]*vrtx[v1][1] +dvrtdt[v1][1]);
         q += dd*dd;
         qmax = MAX(qmax,q);
      }
      hmax = sqrt(hmax);
      h = 2.*jcb/(0.25*(b.p +1)*(b.p+1)*hmax);
      u[0] *= 1./3.;
      u[1] *= 1./3.;

      /* THIS IS GAMMA (DIAGONAL PRECONDITIONER FOR CONTINUITY) */
      gam = qmax +(0.25*h*bd[0] + gbl->nu/h)*(0.25*h*bd[0] + gbl->nu/h); 
      q = sqrt(qmax);
      c = sqrt(qmax+gam);
      lam1  = (q+c);

      /* SET UP DISSIPATIVE COEFFICIENTS */
      gbl->tau[tind]  = adis*h/(jcb*sqrt(gam));
      gbl->delt[tind] = qmax*gbl->tau[tind];

      /* STORE PRECONDITIONER (THIS IS TO DRIVE ITERATION USING NONCONSERVATIVE SYSTEM) */
      dtstari = jcb*(gbl->nu/(h*h) +lam1/h +bd[0]);
#ifdef AXISYMMETRIC
      dtstari *= (vrtx[v[0]][0] +vrtx[v[1]][0] +vrtx[v[2]][0])/3.;
#endif
      gbl->tprcn[tind][0][0]  = dtstari*gbl->rho;
   	/* gbl->tprcn[tind][1][1] = gbl->tprcn[tind][0][0]; */
      gbl->tprcn[tind][NV-1][NV-1] =  dtstari/gam;
      gbl->tprcn[tind][0][NV-1] = u[0]/gam*dtstari;
      gbl->tprcn[tind][1][NV-1] = u[1]/gam*dtstari;
      
      for(i=0;i<3;++i) {
         /* ASSEMBLE VERTEX PRECONDITIONER */
         gbl->vprcn[v[i]][0][0]  += gbl->tprcn[tind][0][0];
         /* gbl->vprcn[v[i]][1][1]  += gbl->tprcn[tind][0][0];  */
         gbl->vprcn[v[i]][NV-1][NV-1]  += gbl->tprcn[tind][NV-1][NV-1];
         gbl->vprcn[v[i]][0][NV-1]  += gbl->tprcn[tind][0][NV-1];
         gbl->vprcn[v[i]][1][NV-1]  += gbl->tprcn[tind][1][NV-1];

         if (b.sm > 0) {
            /* ASSEMBLE SIDE PRECONDITIONER */
            side = tside[tind].side[i];
            gbl->sprcn[side][0][0]  += gbl->tprcn[tind][0][0];
            /* gbl->sprcn[side][1][1]  += gbl->tprcn[tind][0][0];  */
            gbl->sprcn[side][NV-1][NV-1]  += gbl->tprcn[tind][NV-1][NV-1];
            gbl->sprcn[side][0][NV-1]  += gbl->tprcn[tind][0][NV-1];
            gbl->sprcn[side][1][NV-1]  += gbl->tprcn[tind][1][NV-1];
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
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][NV-1];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][1][NV-1];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][NV-1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][1][NV-1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];

         /* SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][0];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][NV-1];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][1][NV-1];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][NV-1][NV-1];
            }
         }
      }
      
      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];

         /* SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][0];
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
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
            gbl->vprcn[v0][0][NV-1] = 0.5*(gbl->vprcn[v0][0][NV-1] +sbuff[i][count++]);
            gbl->vprcn[v0][1][NV-1] = 0.5*(gbl->vprcn[v0][1][NV-1] +sbuff[i][count++]);
            gbl->vprcn[v0][NV-1][NV-1] = 0.5*(gbl->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
         gbl->vprcn[v0][0][NV-1] = 0.5*(gbl->vprcn[v0][0][NV-1] +sbuff[i][count++]);
         gbl->vprcn[v0][1][NV-1] = 0.5*(gbl->vprcn[v0][1][NV-1] +sbuff[i][count++]);
         gbl->vprcn[v0][NV-1][NV-1] = 0.5*(gbl->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (b.sm > 0) {
            for(j=sbdry[i].num-1;j>=0;--j) {
               sind = sbdry[i].el[j];
               gbl->sprcn[sind][0][0] = 0.5*(gbl->sprcn[sind][0][0] +sbuff[i][count++]);
               gbl->sprcn[sind][0][NV-1] = 0.5*(gbl->sprcn[sind][0][NV-1] +sbuff[i][count++]);
               gbl->sprcn[sind][1][NV-1] = 0.5*(gbl->sprcn[sind][1][NV-1] +sbuff[i][count++]);
               gbl->sprcn[sind][NV-1][NV-1] = 0.5*(gbl->sprcn[sind][NV-1][NV-1] +sbuff[i][count++]);
            }
         }
      }

      /* ONLY SEND & RECEIVE DIAGV'S FOR INTERFACE */
      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (b.sm > 0) {
            for(j=sbdry[i].num-1;j>=0;--j) {
               sind = sbdry[i].el[j];
               gbl->sprcn[sind][0][0] = 0.5*(gbl->sprcn[sind][0][0] +sbuff[i][count++]);
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
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][NV-1];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][1][NV-1];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][NV-1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][1][NV-1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];

         /* SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][0];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][NV-1];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][1][NV-1];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][NV-1][NV-1];
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
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
            gbl->vprcn[v0][0][NV-1] = 0.5*(gbl->vprcn[v0][0][NV-1] +sbuff[i][count++]);
            gbl->vprcn[v0][1][NV-1] = 0.5*(gbl->vprcn[v0][1][NV-1] +sbuff[i][count++]);
            gbl->vprcn[v0][NV-1][NV-1] = 0.5*(gbl->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         gbl->vprcn[v0][0][0] = 0.5*(gbl->vprcn[v0][0][0] +sbuff[i][count++]);
         gbl->vprcn[v0][0][NV-1] = 0.5*(gbl->vprcn[v0][0][NV-1] +sbuff[i][count++]);
         gbl->vprcn[v0][1][NV-1] = 0.5*(gbl->vprcn[v0][1][NV-1] +sbuff[i][count++]);
         gbl->vprcn[v0][NV-1][NV-1] = 0.5*(gbl->vprcn[v0][NV-1][NV-1] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (b.sm > 0) {
            for(j=sbdry[i].num-1;j>=0;--j) {
               sind = sbdry[i].el[j];
               gbl->sprcn[sind][0][0] = 0.5*(gbl->sprcn[sind][0][0] +sbuff[i][count++]);
               gbl->sprcn[sind][0][NV-1] = 0.5*(gbl->sprcn[sind][0][NV-1] +sbuff[i][count++]);
               gbl->sprcn[sind][1][NV-1] = 0.5*(gbl->sprcn[sind][1][NV-1] +sbuff[i][count++]);
               gbl->sprcn[sind][NV-1][NV-1] = 0.5*(gbl->sprcn[sind][NV-1][NV-1] +sbuff[i][count++]);
            }
         }
      }
   }
   
   /* INVERT PRECONDITIONER FOR VERTICES */
   for(i=0;i<nvrtx;++i) {
      gbl->vprcn[i][0][0]  = 1.0/(b.vdiag*gbl->vprcn[i][0][0]);
      gbl->vprcn[i][NV-1][NV-1]  = 1.0/(b.vdiag*gbl->vprcn[i][NV-1][NV-1]);
      gbl->vprcn[i][0][NV-1] = -b.vdiag*gbl->vprcn[i][0][NV-1]*gbl->vprcn[i][0][0]*gbl->vprcn[i][NV-1][NV-1];
      gbl->vprcn[i][1][NV-1] = -b.vdiag*gbl->vprcn[i][1][NV-1]*gbl->vprcn[i][0][0]*gbl->vprcn[i][NV-1][NV-1];      
   }
   
   if (b.sm > 0) {
      /* INVERT PRECONDITIONER FOR SIDES */            
      for(i=0;i<nside;++i) {
         gbl->sprcn[i][0][0] = 1.0/gbl->sprcn[i][0][0];
         gbl->sprcn[i][NV-1][NV-1] = 1.0/gbl->sprcn[i][NV-1][NV-1];
         gbl->sprcn[i][0][NV-1] = -gbl->sprcn[i][0][NV-1]*gbl->sprcn[i][0][0]*gbl->sprcn[i][NV-1][NV-1];
         gbl->sprcn[i][1][NV-1] = -gbl->sprcn[i][1][NV-1]*gbl->sprcn[i][0][0]*gbl->sprcn[i][NV-1][NV-1];
      }
   }
   
   /* SET UP TSTEP FOR ACTIVE BOUNDARIES */   
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfdt2(i);
   
   return;
}
#endif
