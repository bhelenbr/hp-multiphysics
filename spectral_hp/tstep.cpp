#include"hp_mgrid.h"
#include<math.h>
#include<utilities.h>

void hp_mgrid::tstep1(void) {
	static int tind,i,j,side,v0,*v;
	static FLT jcb,h,hmax,q,qmax,lam1,c;

/***************************************/
/** DETERMINE FLOW PSEUDO-TIME STEP ****/
/***************************************/
	for(i=0;i<nvrtx;++i) {
		gbl.vdiagv[i] = 0.0;
		gbl.vdiagp[i] = 0.0;
	}

	if (b.sm > 0) {
		for(i=0;i<nside;++i) {
			gbl.sdiagv[i] = 0.0;
			gbl.sdiagp[i] = 0.0;
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
#ifdef MOVING_MESH
		for(j=0;j<3;++j) {
			v0 = v[j];
			q = pow(vug[v0][0]-0.5*mv[v0][0],2.0) 
				+pow(vug[v0][1]-0.5*mv[v0][1],2.0);
			qmax = MAX(qmax,q);
		}
#else
		for(j=0;j<3;++j) {
			v0 = v[j];
			q = vug[v0][0]*vug[v0][0] +vug[v0][1]*vug[v0][1];
			qmax = MAX(qmax,q);
		}
#endif
		gbl.gam[tind] = qmax +(0.25*h*dt0 + gbl.nu/h)*(0.25*h*dt0 + gbl.nu/h);  
		
		q = sqrt(qmax);
		c = sqrt(qmax+gbl.gam[tind]);
		lam1  = (q+c);
		
		gbl.dtstar[tind]  = gbl.nu/(h*h) +lam1/h +dt0;
		gbl.dtstar[tind] *= jcb*gbl.rho;

/*		SET UP DISSIPATIVE COEFFICIENTS */
		gbl.tau[tind]  = gbl.adis*h/(jcb*sqrt(gbl.gam[tind]));
		gbl.delt[tind] = qmax*gbl.tau[tind];
		gbl.gam[tind] =  1.0/gbl.gam[tind];
		for(i=0;i<3;++i) {
			gbl.vdiagv[v[i]]  += gbl.dtstar[tind]*b.vdiag;
			gbl.vdiagp[v[i]]  += gbl.dtstar[tind]*gbl.gam[tind]*b.vdiag*gbl.rhoi;
			if (b.sm > 0) {
				side = tside[tind].side[i];
				gbl.sdiagv[side] += gbl.dtstar[tind];
				gbl.sdiagp[side] += gbl.dtstar[tind]*gbl.gam[tind]*gbl.rhoi;
			}
		}
	}

/*	SEND COMMUNICATION PACKETS TO OTHER MESHES HERE */
   
   return;
}

void hp_mgrid::tstep2(void) {
   static int i;
   
/*	RECEIVE COMMUNCATION PACKETS FROM OTHER MESHES */

/*	FORM DIAGANOL PRECONDITIONER FOR VERTICES */
   for(i=0;i<nvrtx;++i) {
      gbl.vdiagv[i]  = 1.0/gbl.vdiagv[i];
      gbl.vdiagp[i]  = 1.0/gbl.vdiagp[i];
   }
      
   if (b.sm > 0) {
/*		FORM DIAGANOL PRECONDITIONER FOR SIDES */				
		for(i=0;i<nside;++i) {
			gbl.sdiagv[i] = 1.0/gbl.sdiagv[i];
			gbl.sdiagp[i] = 1.0/gbl.sdiagp[i];
		}
	}
   
	return;
}
