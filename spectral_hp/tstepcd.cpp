#include"hp_mgrid.h"
#include<math.h>
#include<utilities.h>

extern FLT axext, ayext, nuext;

void hp_mgrid::tstep1(void) {
   int tind,i,j,sind,count,bnum,side,v0,*v;
   FLT jcb,h,hmax,lam1,dtstari;
   class mesh *tgt;

   /***************************************/
   /** DETERMINE FLOW PSEUDO-TIME STEP ****/
   /***************************************/
   for(i=0;i<nvrtx;++i) {
      hp_gbl->vprcn[i][0][0] = 0.0;
   }

   if (basis::tri(log2p).sm > 0) {
      for(i=0;i<nside;++i) {
         hp_gbl->sprcn[i][0][0] = 0.0;
      }
   }
   
#ifdef TIME_ACCURATE
   dtstari = 0.0;
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
         printf("negative triangle area\n");
         output("negative",tecplot);
         exit(1);
      }
      
      h   =  4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
      lam1  = (sqrt(axext*axext +ayext*ayext) +1.5*nuext/h +h*sim::bd[0]);

      /* SET UP DISSIPATIVE COEFFICIENTS */
      hp_gbl->tau(tind)  = adis*h/(jcb*lam1);
      
      /* SET UP DIAGONAL PRECONDITIONER */
      dtstari = MAX(lam1/h,dtstari);
   }
   
   printf("#iterative to physical time step ratio: %f %e %e\n",sim::bd[0]/dtstari,sim::bd[0],dtstari);
      
   for(tind=0;tind<ntri;++tind) {
      v = tvrtx[tind];
      jcb = 0.25*area(tind)*dtstari;
#ifdef AXISYMMETRIC
      jcb *= (vrtx[v[0]][0] +vrtx[v[1]][0] +vrtx[v[2]][0])/3.;
#endif
      hp_gbl->tprcn[tind][0][0] = jcb;      
      for(i=0;i<3;++i) {
         hp_gbl->vprcn[v[i]][0][0]  += hp_gbl->tprcn[tind][0][0];
         if (basis::tri(log2p).sm > 0) {
            side = td(tind).side(i);
            hp_gbl->sprcn[side][0][0] += hp_gbl->tprcn[tind][0][0];
         }
      }
   } 
#else
   for(tind = 0; tind < ntri; ++tind) {
   	/* THIS IS THE JACOBIAN WEIGHTING FACTOR FOR */
      /* THE STANDARD TRIANGLE */
      /* FACTOR OF 4 BECAUSE OF LENGTH 2 SIDES ON ST */
      jcb = 0.25*area(tind);
      if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
         printf("negative triangle area\n");
         output("negative",tecplot);
         exit(1);
      }
          
      v = tvrtx[tind];
      hmax = 0.0;
      for(j=0;j<3;++j) {
         h = pow(vrtx[v[j]][1] -vrtx[v[(j+1)%3]][1],2.0) + 
         pow(vrtx[v[j]][0] -vrtx[v[(j+1)%3]][0],2.0);
         hmax = (h > hmax ? h : hmax);
      }
      hmax = sqrt(hmax);

      /* SMALLEST LENGTH FOR CALCULATING EIGENVALUES */
      /* MOST RESTRICTIVE */
      h   =  4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
      
      /* 1.5 IS BASED ON COMPARISON TO R_MESH FOR LAPLACE */
      lam1  = (sqrt(axext*axext +ayext*ayext) +1.5*nuext/h +h*sim::bd[0]);

      /* SET UP DISSIPATIVE COEFFICIENTS */
      hp_gbl->tau(tind)  = adis*h/(jcb*lam1);
      
      /* SET UP DIAGONAL PRECONDITIONER */
      dtstari = jcb*lam1/h;
         
#ifdef AXISYMMETRIC
      dtstari *= (vrtx[v[0]][0] +vrtx[v[1]][0] +vrtx[v[2]][0])/3.;
#endif
      hp_gbl->tprcn[tind][0][0] = dtstari;      
      for(i=0;i<3;++i) {
         hp_gbl->vprcn[v[i]][0][0]  += hp_gbl->tprcn[tind][0][0];

         if (basis::tri(log2p).sm > 0) {
            side = td(tind).side(i);
            hp_gbl->sprcn[side][0][0] += hp_gbl->tprcn[tind][0][0];
         }
      }
   }
#endif
   
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
            tgt->sbuff[bnum][count++] = hp_gbl->vprcn[v0][0][0];
         }
         v0 = sd(sind).vrtx(1);
         tgt->sbuff[bnum][count++] = hp_gbl->vprcn[v0][0][0];

         /* SEND SIDE INFO */
         if (basis::tri(log2p).sm) {
            for(j=0;j<sbdry(i)->nel;++j) {
               sind = sbdry(i)->el(j);
               tgt->sbuff[bnum][count++] = hp_gbl->sprcn[sind][0][0];
            }
         }
      }
   }
   
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
            hp_gbl->vprcn[v0][0][0] = 0.5*(hp_gbl->vprcn[v0][0][0] +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         hp_gbl->vprcn[v0][0][0] = 0.5*(hp_gbl->vprcn[v0][0][0] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (basis::tri(log2p).sm > 0) {
            for(j=sbdry(i)->nel-1;j>=0;--j) {
               sind = sbdry(i)->el(j);
               hp_gbl->sprcn[sind][0][0] = 0.5*(hp_gbl->sprcn[sind][0][0] +sbuff[i][count++]);
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
            tgt->sbuff[bnum][count++] = hp_gbl->vprcn[v0][0][0];
         }
         v0 = sd(sind).vrtx(1);
         tgt->sbuff[bnum][count++] = hp_gbl->vprcn[v0][0][0];

         /* SEND SIDE INFO */
         if (basis::tri(log2p).sm) {
            for(j=0;j<sbdry(i)->nel;++j) {
               sind = sbdry(i)->el(j);
               tgt->sbuff[bnum][count++] = hp_gbl->sprcn[sind][0][0];
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
            hp_gbl->vprcn[v0][0][0] = 0.5*(hp_gbl->vprcn[v0][0][0] +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         hp_gbl->vprcn[v0][0][0] = 0.5*(hp_gbl->vprcn[v0][0][0] +sbuff[i][count++]);

         /* RECV SIDE INFO */
         if (basis::tri(log2p).sm > 0) {
            for(j=sbdry(i)->nel-1;j>=0;--j) {
               sind = sbdry(i)->el(j);
               hp_gbl->sprcn[sind][0][0] = 0.5*(hp_gbl->sprcn[sind][0][0] +sbuff[i][count++]);
            }
         }
      }
   }
      
   /* FORM DIAGANOL PRECONDITIONER FOR VERTICES */
   for(i=0;i<nvrtx;++i) {
      hp_gbl->vprcn[i][0][0]  = 1.0/(basis::tri(log2p).vdiag*hp_gbl->vprcn[i][0][0]);
   }

   if (basis::tri(log2p).sm > 0) {
      /* FORM DIAGANOL PRECONDITIONER FOR SIDES */            
      for(i=0;i<nside;++i) {
         hp_gbl->sprcn[i][0][0] = 1.0/hp_gbl->sprcn[i][0][0];
      }
   }
   
   return;
}