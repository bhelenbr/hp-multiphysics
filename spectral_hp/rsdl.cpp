/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "hp_mgrid.h"
   
void hp_mgrid::rsdl(int stage, int mgrid) {
   static int i,j,n,tind;
   static FLT fluxx,fluxy;
   static FLT visc0[ND][ND][ND], tres[NV];
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl.vres[i][n] = 0.0;

    for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         gbl.sres[i][n] = 0.0;
                 
    for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         gbl.ires[i][n] = 0.0;  

	for(tind = 0; tind<ntri;++tind) {
   
      if (tinfo[tind] > -1) {
         crdtouht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(uht[n], crd[n], dcrd[n][0], dcrd[n][1]);

         ugtouht(tind);
         if (beta[stage] > 0.0) {
            b.proj(uht[0],u[0],du[0][0],du[0][1]);
            b.proj(uht[1],u[1],du[1][0],du[1][1]);
            b.proj(uht[2],u[2]);
         }
         else {
            b.proj(uht[0],u[0]);
            b.proj(uht[1],u[1]);
            b.proj(uht[2],u[2]);
         }

#ifdef MOVING_MESH
         for(i=0;i<b.gpx;++i)
            for(j=0;j<b.gpn;++i)
               for(n=0;n<ND;++n)
                  crd[n][i][j] = dt0*crd[n][i][j] +xdt[tind][n][i][j];
#endif			

         for(n=0;n<NV;++n)
            for(i=0;i<b.tm;++i)
               lf[n][i] = 0.0;

/*			CONVECTIVE TERMS (IMAGINARY FIRST)*/
         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {
#ifdef MOVING_MESH
               fluxx = gbl.rho*(u[0][i][j] -crd[0][i][j]);
               fluxy = gbl.rho*(u[1][i][j] -crd[1][i][j]);
#else
               fluxx = gbl.rho*u[0][i][j];
               fluxy = gbl.rho*u[1][i][j];
#endif
               du[2][0][i][j] = -dcrd[1][1][i][j]*fluxx +dcrd[0][1][i][j]*fluxy;
               du[2][1][i][j] = +dcrd[1][0][i][j]*fluxx -dcrd[0][0][i][j]*fluxy;
   
               cv00[i][j] =  u[0][i][j]*du[2][0][i][j] -dcrd[1][1][i][j]*u[2][i][j];
               cv01[i][j] =  u[0][i][j]*du[2][1][i][j] +dcrd[1][0][i][j]*u[2][i][j];
               cv10[i][j] =  u[1][i][j]*du[2][0][i][j] +dcrd[0][1][i][j]*u[2][i][j];
               cv11[i][j] =  u[1][i][j]*du[2][1][i][j] -dcrd[0][0][i][j]*u[2][i][j];
            }
         }
         b.intgrtrs(cv00,cv01,lf[0]);
         b.intgrtrs(cv10,cv11,lf[1]);
         b.intgrtrs(du[2][0],du[2][1],lf[2]);

/*			ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
         lftog(tind,gbl.vres,gbl.sres,gbl.ires);

/*			NEGATIVE REAL TERMS */
         if (beta[stage] > 0.0) {
         		
/*				TIME DERIVATIVE TERMS */
            if (!mgrid) {			
               for(i=0;i<b.gpx;++i) {
                  for(j=0;j<b.gpn;++j) {
                     cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
                     res[0][i][j] = gbl.rho*dt0*u[0][i][j]*cjcb[i][j]; // +udt[tind][0][i][j];
                     res[1][i][j] = gbl.rho*dt0*u[1][i][j]*cjcb[i][j];// +udt[tind][1][i][j];
                     res[2][i][j] = gbl.rho*dt0*cjcb[i][j];// +udt[tind][2][i][j];
                  }
               }
            }
            else {
               for(i=0;i<b.gpx;++i) {
                  for(j=0;j<b.gpn;++j) {
                     cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
                     res[0][i][j] = gbl.rho*dt0*u[0][i][j]*cjcb[i][j];
                     res[1][i][j] = gbl.rho*dt0*u[1][i][j]*cjcb[i][j];
                     res[2][i][j] = gbl.rho*dt0*cjcb[i][j];
                  }
               }
            }
            b.intgrt(res[0],lf[0]);
            b.intgrt(res[1],lf[1]);
            b.intgrt(res[2],lf[2]);

/*				VISCOUS TERMS  */
            for(i=0;i<b.gpx;++i) {
               for(j=0;j<b.gpn;++j) {

                  cjcb[i][j] = gbl.mu/cjcb[i][j];
/*						BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
/*						INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: DERIVATIVE (R OR S) 4: DERIVATIVE (R OR S)*/
                  visc0[0][0][0] =  cjcb[i][j]*(2.*dcrd[1][1][i][j]*dcrd[1][1][i][j] +dcrd[0][1][i][j]*dcrd[0][1][i][j]);
                  visc0[0][1][1] =  cjcb[i][j]*(2.*dcrd[1][0][i][j]*dcrd[1][0][i][j] +dcrd[0][0][i][j]*dcrd[0][0][i][j]);
                  visc0[0][0][1] = -cjcb[i][j]*(2.*dcrd[1][1][i][j]*dcrd[1][0][i][j] +dcrd[0][1][i][j]*dcrd[0][0][i][j]);
/*    	         visc0[0][1][0] = visc0[0][0][1];  */
                  
                  visc0[1][0][0] = -cjcb[i][j]*dcrd[0][1][i][j]*dcrd[1][1][i][j];
                  visc0[1][1][1] = -cjcb[i][j]*dcrd[0][0][i][j]*dcrd[1][0][i][j];
                  visc0[1][0][1] =  cjcb[i][j]*dcrd[0][1][i][j]*dcrd[1][0][i][j];
                  visc0[1][1][0] =  cjcb[i][j]*dcrd[0][0][i][j]*dcrd[0][1][i][j];
   
/*						OTHER SYMMETRIES                
                  visc1[0][0][0] = visc0[1][0][0];
                  visc1[0][1][1] = visc0[1][1][1];
                  visc1[0][0][1] = visc0[1][1][0];
                  visc1[0][1][0] = visc0[1][0][1];
                  
                  visc1[1][0][0] = visc0[0][0][0];
                  visc1[1][1][1] = visc0[0][1][1];  
                  visc1[1][0][1] = visc0[0][0][1];
                  visc1[1][1][0] = visc0[0][1][0] = visc00[0][1];
*/
   
                  e00[i][j] = +visc0[0][0][0]*du[0][0][i][j] +visc0[1][0][0]*du[1][0][i][j]
                              +visc0[0][0][1]*du[0][1][i][j] +visc0[1][0][1]*du[1][1][i][j];
   
                  e01[i][j] = +visc0[0][0][1]*du[0][0][i][j] +visc0[1][1][0]*du[1][0][i][j]
                              +visc0[0][1][1]*du[0][1][i][j] +visc0[1][1][1]*du[1][1][i][j];
   
                  e10[i][j] = +visc0[1][0][0]*du[0][0][i][j] +visc0[0][0][0]*du[1][0][i][j]
                              +visc0[1][1][0]*du[0][1][i][j] +visc0[0][0][1]*du[1][1][i][j];
   
                  e11[i][j] = +visc0[1][0][1]*du[0][0][i][j] +visc0[0][0][1]*du[1][0][i][j]
                              +visc0[1][1][1]*du[0][1][i][j] +visc0[0][1][1]*du[1][1][i][j];
   
   
                  cv00[i][j] = -cv00[i][j] -e00[i][j];
                  cv01[i][j] = -cv01[i][j] -e01[i][j];
                  cv10[i][j] = -cv10[i][j] -e10[i][j];
                  cv11[i][j] = -cv11[i][j] -e11[i][j];
   
#ifdef MOVING_MESH
                  res[2][i][j] += rho*
                                 (dcrd[1][1][i][j]*(du[0][0][i][j] -dlmvx[0][i][j])
                                 -dcrd[0][1][i][j]*(du[1][0][i][j] -dlmvx[1][i][j])
                                 -dcrd[1][0][i][j]*(du[0][1][i][j] -dlmvy[0][i][j])
                                 +dcrd[0][0][i][j]*(du[1][1][i][j] -dlmvy[1][i][j]));
#else					
                  res[2][i][j] += gbl.rho*(dcrd[1][1][i][j]*du[0][0][i][j] -dcrd[0][1][i][j]*du[1][0][i][j]
                                    -dcrd[1][0][i][j]*du[0][1][i][j] +dcrd[0][0][i][j]*du[1][1][i][j]);
#endif
   
               }
            }
            b.derivr(cv00,res[0]);
            b.derivs(cv01,res[0]);
            b.derivr(cv10,res[1]);
            b.derivs(cv11,res[1]);

/*				THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
            for(i=0;i<b.gpx;++i) {
               for(j=0;j<b.gpn;++j) {
                  tres[0] = gbl.tau[tind]*res[0][i][j];
                  tres[1] = gbl.tau[tind]*res[1][i][j];
                  tres[2] = gbl.delt[tind]*res[2][i][j];

#ifdef MOVING_MESH							
                  e00[i][j] += (dcrd[1][1][i][j]*(2*u[0][i][j]-crd[0][i][j])
                              -dcrd[0][1][i][j]*(u[1][i][j]-crd[1][i][j]))*tres[0]
                              -dcrd[0][1][i][j]*u[0][i][j]*tres[1]
                              +dcrd[1][1][i][j]*tres[2];
                  e01[i][j] += (-dcrd[1][0][i][j]*(2*u[0][i][j]-crd[0][i][j])
                              +dcrd[0][0][i][j]*(u[1][i][j]-crd[1][i][j]))*tres[0]
                              +dcrd[0][0][i][j]*u[0][i][j]*tres[1]
                              -dcrd[1][0][i][j]*tres[2];
                  e10[i][j] += +dcrd[1][1][i][j]*u[1][i][j]*tres[0]
                              +(dcrd[1][1][i][j]*(u[0][i][j]-crd[0][i][j])
                              -dcrd[0][1][i][j]*(2.*u[1][i][j]-crd[1][i][j]))*tres[1]
                              -dcrd[0][1][i][j]*tres[2];
                  e11[i][j] += -dcrd[1][0][i][j]*u[1][i][j]*tres[0]
                              +(-dcrd[1][0][i][j]*(u[0][i][j]-crd[0][i][j])
                              +dcrd[0][0][i][j]*(2.*u[1][i][j]-crd[1][i][j]))*tres[1]
                              +dcrd[0][0][i][j]*tres[2];
#else
                  e00[i][j] += (dcrd[1][1][i][j]*(2*u[0][i][j])
                              -dcrd[0][1][i][j]*(u[1][i][j]))*tres[0]
                              -dcrd[0][1][i][j]*u[0][i][j]*tres[1]
                              +dcrd[1][1][i][j]*tres[2];
                  e01[i][j] += (-dcrd[1][0][i][j]*(2*u[0][i][j])
                              +dcrd[0][0][i][j]*(u[1][i][j]))*tres[0]
                              +dcrd[0][0][i][j]*u[0][i][j]*tres[1]
                              -dcrd[1][0][i][j]*tres[2];
                  e10[i][j] += +dcrd[1][1][i][j]*u[1][i][j]*tres[0]
                              +(dcrd[1][1][i][j]*(u[0][i][j])
                              -dcrd[0][1][i][j]*(2.*u[1][i][j]))*tres[1]
                              -dcrd[0][1][i][j]*tres[2];
                  e11[i][j] += -dcrd[1][0][i][j]*u[1][i][j]*tres[0]
                              +(-dcrd[1][0][i][j]*(u[0][i][j])
                              +dcrd[0][0][i][j]*(2.*u[1][i][j]))*tres[1]
                              +dcrd[0][0][i][j]*tres[2];
#endif
                              
                  du[2][0][i][j] = (dcrd[1][1][i][j]*tres[0] -dcrd[0][1][i][j]*tres[1]);
                  du[2][1][i][j] = (-dcrd[1][0][i][j]*tres[0] +dcrd[0][0][i][j]*tres[1]);
               }
            }
            b.intgrtrs(e00,e01,lf[0]);
            b.intgrtrs(e10,e11,lf[1]);
            b.intgrtrs(du[2][0],du[2][1],lf[2]);
   
            for(n=0;n<NV;++n)
               for(i=0;i<b.tm;++i)
                  lf[n][i] *= beta[stage];
                  
            lftog(tind,gbl.vvf,gbl.svf,gbl.ivf);
         }
      }
	}
   
/*	ADD IN VISCOUS/DISSIPATIVE FLUX */
	for(i=0;i<nvrtx;++i)
		for(n=0;n<NV;++n)
			gbl.vres[i][n] += gbl.vvf[i][n];

 	for(i=0;i<nside*b.sm;++i)
		for(n=0;n<NV;++n)
			gbl.sres[i][n] += gbl.svf[i][n];        

	for(i=0;i<ntri*b.im;++i)
		for(n=0;n<NV;++n)
			gbl.ires[i][n] += gbl.ivf[i][n];         

/*	ADD IN BOUNDARY SOURCES */
//	flowsrc();
	
/*********************************************/
/*	MODIFY RESIDUALS ON COARSER MESHES			*/
/*********************************************/	
	if(mgrid) {
/* 	CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
		if(isfrst) {
         for(i=0;i<nvrtx;++i)
            for(n=0;n<NV;++n)
               vdres[log2p][i][n] = gbl.fadd*gbl.vres0[i][n] - gbl.vres[i][n];
      
         for(i=0;i<nside*b.sm;++i)
            for(n=0;n<NV;++n)
               sdres[log2p][i][n] = gbl.fadd*gbl.sres0[i][n] - gbl.sres[i][n];        
      
         for(i=0;i<ntri*b.im;++i)
            for(n=0;n<NV;++n)
               idres[log2p][i][n] = gbl.fadd*gbl.sres0[i][n] - gbl.sres[i][n];
		}

		for(i=0;i<nvrtx;++i) 
			for(n=0;n<NV;++n)		
				gbl.vres[i][n] += vdres[log2p][i][n];	
            
		for(i=0;i<nside*b.sm;++i) 
			for(n=0;n<NV;++n)		
				gbl.sres[i][n] += sdres[log2p][i][n];
            
		for(i=0;i<ntri*b.im;++i) 
			for(n=0;n<NV;++n)		
				gbl.sres[i][n] += sdres[log2p][i][n];      
	}

/*	APPLY BOUNDARY CONDITIONS */
//	flowbnds(); 

/* SEND COMMUNICATION PACKETS */   
//   send(XDIR_MP, (FLT *) res, 0, 1, 2);
   
	return;
}
