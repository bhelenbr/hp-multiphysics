#include"hp_mgrid.h"
#include<math.h>
#include<utilities.h>

void hp_mgrid::tstep1(void) {
	static int tind,i,j,sind,count,bnum,side,v0,*v;
	static FLT jcb,h,hmax,q,qmax,lam1,c;
   class mesh *tgt;

/***************************************/
/** DETERMINE FLOW PSEUDO-TIME STEP ****/
/***************************************/
	for(i=0;i<nvrtx;++i) {
		gbl->vdiagv[i] = 0.0;
		gbl->vdiagp[i] = 0.0;
	}

	if (b.sm > 0) {
		for(i=0;i<nside;++i) {
			gbl->sdiagv[i] = 0.0;
			gbl->sdiagp[i] = 0.0;
		}
	}

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
		h = 2.*jcb/(0.25*(b.p +1)*(b.p+1)*hmax);
		
		qmax = 0.0;
		for(j=0;j<3;++j) {
			v0 = v[j];
			q = pow(ug.v[v0][0]-0.5*(bd[0]*vrtx[v0][0] -dvrtdt[v0][0]),2.0) 
				+pow(ug.v[v0][1]-0.5*(bd[0]*vrtx[v0][1] -dvrtdt[v0][1]),2.0);
			qmax = MAX(qmax,q);
		}
		gbl->gam[tind] = qmax +(0.25*h*bd[0] + gbl->nu/h)*(0.25*h*bd[0] + gbl->nu/h);  
		
		q = sqrt(qmax);
		c = sqrt(qmax+gbl->gam[tind]);
		lam1  = (q+c);
		
		gbl->dtstar[tind]  = gbl->nu/(h*h) +lam1/h +bd[0];
		gbl->dtstar[tind] *= jcb*gbl->rho;

/*		SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau[tind]  = adis*h/(jcb*sqrt(gbl->gam[tind]));
		gbl->delt[tind] = qmax*gbl->tau[tind];
		gbl->gam[tind] =  1.0/gbl->gam[tind];
		for(i=0;i<3;++i) {
			gbl->vdiagv[v[i]]  += gbl->dtstar[tind]*b.vdiag;
			gbl->vdiagp[v[i]]  += gbl->dtstar[tind]*gbl->gam[tind]*b.vdiag*gbl->rhoi;
			if (b.sm > 0) {
				side = tside[tind].side[i];
				gbl->sdiagv[side] += gbl->dtstar[tind];
				gbl->sdiagp[side] += gbl->dtstar[tind]*gbl->gam[tind]*gbl->rhoi;
			}
		}
	}
   
/*	SEND BOUNDARY INFORMATION */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & HP_MGRID_MP) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
/*			SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vdiagv[v0];
            tgt->sbuff[bnum][count++] = gbl->vdiagp[v0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vdiagv[v0];
         tgt->sbuff[bnum][count++] = gbl->vdiagp[v0];

/*			SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sdiagv[sind];
               tgt->sbuff[bnum][count++] = gbl->sdiagp[sind];
            }
         }
      }
      
      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
/*			SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vdiagv[v0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vdiagv[v0];

/*			SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sdiagv[sind];
            }
         }
      }
   }
   
/*	SET UP TSTEP FOR ACTIVE BOUNDARIES */   
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfdt1(i);
   
   return;
}

void hp_mgrid::tstep2(void) {
   static int i,j,sind,v0,count;
   
/*	RECEIVE COMMUNCATION PACKETS FROM OTHER MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & HP_MGRID_MP) {
         count = 0;
/*			RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            gbl->vdiagv[v0] = 0.5*(gbl->vdiagv[v0] +sbuff[i][count++]);
            gbl->vdiagp[v0] = 0.5*(gbl->vdiagp[v0] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         gbl->vdiagv[v0] = 0.5*(gbl->vdiagv[v0] +sbuff[i][count++]);
         gbl->vdiagp[v0] = 0.5*(gbl->vdiagp[v0] +sbuff[i][count++]);

/*			RECV SIDE INFO */
         if (b.sm > 0) {
            for(j=sbdry[i].num-1;j>=0;--j) {
               sind = sbdry[i].el[j];
               gbl->sdiagv[sind] = 0.5*(gbl->sdiagv[sind] +sbuff[i][count++]);
               gbl->sdiagp[sind] = 0.5*(gbl->sdiagp[sind] +sbuff[i][count++]);
            }
         }
      }

/*		ONLY SEND & RECEIVE DIAGV'S FOR INTERFACE */
      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
/*			RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            gbl->vdiagv[v0] = 0.5*(gbl->vdiagv[v0] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         gbl->vdiagv[v0] = 0.5*(gbl->vdiagv[v0] +sbuff[i][count++]);

/*			RECV SIDE INFO */
         if (b.sm > 0) {
            for(j=sbdry[i].num-1;j>=0;--j) {
               sind = sbdry[i].el[j];
               gbl->sdiagv[sind] = 0.5*(gbl->sdiagv[sind] +sbuff[i][count++]);
            }
         }
      }
   }
   
/*	FORM DIAGANOL PRECONDITIONER FOR VERTICES */
   for(i=0;i<nvrtx;++i) {
      gbl->vdiagv[i]  = 1.0/gbl->vdiagv[i];
      gbl->vdiagp[i]  = 1.0/gbl->vdiagp[i];
   }
      
   if (b.sm > 0) {
/*		FORM DIAGANOL PRECONDITIONER FOR SIDES */				
		for(i=0;i<nside;++i) {
			gbl->sdiagv[i] = 1.0/gbl->sdiagv[i];
			gbl->sdiagp[i] = 1.0/gbl->sdiagp[i];
		}
	}
   
/*	SET UP TSTEP FOR ACTIVE BOUNDARIES */   
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfdt2(i);
   
	return;
}