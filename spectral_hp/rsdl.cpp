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
   static int i,j,n,sind,v0,v1;
   static FLT dx,dy;
   static FLT visc0[ND][ND][ND];
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         vres[i][n] = 0.0;

    for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         sres[i][n] = 0.0;
                 
    for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         ires[i][n] = 0.0;  

	for(tind = 0; tind<ntri;++tind) {
   
      if (tinfo[tind] > -1) {
         crdtouht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(uht[n], crd[n], dcrd[n][0], dcrd[n][1]);

         ugtouht(tind);
         if (beta[stage] > 0.0) {
            proj(uht[0],u[0],du[0][0],du[0][1]);
            proj(uht[1],u[1],du[1][0],du[1][1]);
            proj2(uht[2],u[2]);
         }
         else {
            proj2(uht[0],u[0]);
            proj2(uht[1],u[1]);
            proj2(uht[2],u[2]);
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
               fluxx = rho*(u[0][i][j] -crd[0][i][j]);
               fluxy = rho*(u[1][i][j] -crd[1][i][j]);
#else
               fluxx = rho*u[0][i][j];
               fluxy = rho*u[1][i][j];
#endif
               dur[2][i][j] = -dcrd[1][1][i][j]*fluxx +dcrd[0][1][i][j]*fluxy;
               dus[2][i][j] = +dcrd[1][0][i][j]*fluxx -dcrd[0][0][i][j]*fluxy;
   
               cv00[i][j] =  u[0][i][j]*dur[2][i][j] -dcrd[1][1][i][j]*u[2][i][j];
               cv01[i][j] =  u[0][i][j]*dus[2][i][j] +dcrd[1][0][i][j]*u[2][i][j];
               cv10[i][j] =  u[1][i][j]*dur[2][i][j] +dcrd[0][1][i][j]*u[2][i][j];
               cv11[i][j] =  u[1][i][j]*dus[2][i][j] -dcrd[0][0][i][j]*u[2][i][j];
            }
         }
         b.intgrtrs(cv00,cv01,lf[0]);
         b.intgrtrs(cv10,cv11,lf[1]);
         b.intgrtrs(dur[2],dus[2],lf[2]);

/*			ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
         lftog(tind,vgf,sgf,igf);

/*			NEGATIVE REAL TERMS */
         if (beta[stage] > 0.0) {
         		
/*				TIME DERIVATIVE TERMS */
            if (!mgrid) {			
               for(i=0;i<b.gpx;++i) {
                  for(j=0;j<b.gpn;++j) {
                     cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
                     res[0][i][j] = rho*dt0*u[0][i][j]*cjcb[i][j] +udt[tind][0][i][j];
                     res[1][i][j] = rho*dt0*u[1][i][j]*cjcb[i][j] +udt[tind][1][i][j];
                     res[2][i][j] = rho*dt0*cjcb[i][j] +udt[tind][2][i][j];
                  }
               }
            }
            else {
               for(i=0;i<b.gpx;++i) {
                  for(j=0;j<b.gpn;++j) {
                     cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
                     res[0][i][j] = rho*dt0*u[0][i][j]*cjcb[i][j];
                     res[1][i][j] = rho*dt0*u[1][i][j]*cjcb[i][j];
                     res[2][i][j] = rho*dt0*cjcb[i][j];
                  }
               }
            }
            b.intgrt(res[0],lf[0]);
            b.intgrt(res[1],lf[1]);
            b.intgrt(res[2],lf[2]);

/*				VISCOUS TERMS  */
            for(i=0;i<b.gpx;++i) {
               for(j=0;j<b.gpy;++j) {

               cjcb[i][j] = mu/cjcb[i][j];
/*					BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
/*					INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: DERIVATIVE (R OR S) 4: DERIVATIVE (R OR S)*/
               visc0[0][0][0] =  cjcb[i][j]*(2.*dcrd[1][1][i][j]*dcrd[1][1][i][j] +dcrd[0][1][i][j]*dcrd[0][1][i][j]);
               visc0[0][1][1] =  cjcb[i][j]*(2.*dcrd[1][0][i][j]*dcrd[1][0][i][j] +dcrd[0][0][i][j]*dcrd[0][0][i][j]);
               visc0[0][0][1] = -cjcb[i][j]*(2.*dcrd[1][1][i][j]*dcrd[1][0][i][j] +dcrd[0][1][i][j]*dcrd[0][0][i][j]);
/*             visc0[0][1][0] = visc0[0][0][1];  */
               
               visc0[1][0][0] = -cjcb[i][j]*dcrd[0][1][i][j]*dcrd[1][1][i][j];
               visc0[1][1][1] = -cjcb[i][j]*dcrd[0][0][i][j]*dcrd[1][0][i][j];
               visc0[1][0][1] =  cjcb[i][j]*dcrd[0][1][i][j]*dcrd[1][0][i][j];
               visc0[1][1][0] =  cjcb[i][j]*dcrd[0][0][i][j]*dcrd[0][1][i][j];

/*					OTHER SYMMETRIES                
               visc1[0][0][0] = visc0[1][0][0];
               visc1[0][1][1] = visc0[1][1][1];
               visc1[0][0][1] = visc0[1][1][0];
               visc1[0][1][0] = visc0[1][0][1];
               
               visc1[1][0][0] = visc0[0][0][0];
               visc1[1][1][1] = visc0[0][1][1];  
               visc1[1][0][1] = visc0[0][0][1];
               visc1[1][1][0] = visc0[0][1][0] = visc00[0][1];
*/

					e00[i][j] = +visc0[0][0][0]*dur[0][i][j] +visc0[1][0][0]*dur[1][i][j]
									+visc0[0][0][1]*dus[0][i][j] +visc0[1][0][1]*dus[1][i][j];

					e01[i][j] = +visc0[0][0][1]*dur[0][i][j] +visc0[1][1][0]*dur[1][i][j]
									+visc0[0][1][1]*dus[0][i][j] +visc0[1][1][1]*dus[1][i][j];

					e10[i][j] = +visc0[1][0][0]*dur[0][i][j] +visc0[0][0][0]*dur[1][i][j]
									+visc0[1][1][0]*dus[0][i][j] +visc0[0][0][1]*dus[1][i][j];

					e11[i][j] = +visc0[1][0][1]*dur[0][i][j] +visc0[0][0][1]*dur[1][i][j]
									+visc0[1][1][1]*dus[0][i][j] +visc0[0][1][1]*dus[1][i][j];


					cv00[i][j] = -cv00[i][j] -e00[i][j];
					cv01[i][j] = -cv01[i][j] -e01[i][j];
					cv10[i][j] = -cv10[i][j] -e10[i][j];
					cv11[i][j] = -cv11[i][j] -e11[i][j];

#ifdef MOVING_MESH
					res[2][i][j] += rho*
						 				(dcrd[1][1][i][j]*(dur[0][i][j] -dlmvx[0][i][j])
									   -dcrd[0][1][i][j]*(dur[1][i][j] -dlmvx[1][i][j])
									   -dcrd[1][0][i][j]*(dus[0][i][j] -dlmvy[0][i][j])
									   +dcrd[0][0][i][j]*(dus[1][i][j] -dlmvy[1][i][j]));
#else					
               res[2][i][j] += rho*(dcrd[1][1][i][j]*dur[0][i][j] -dcrd[0][1][i][j]*dur[1][i][j]
                                   -dcrd[1][0][i][j]*dus[0][i][j] +dcrd[0][0][i][j]*dus[1][i][j]);
#endif
   
				}
			}
			b.derivr(cv00,res[0]);
			b.derivs(cv01,res[0]);
			b.derivr(cv10,res[1]);
			b.derivs(cv11,res[1]);

/*			THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(i=0;i<b.gpx;++i) {
				for(j=0;j<b.gpn;++j) {
					tres[0] = tau[tind]*res[0][i][j];
					tres[1] = tau[tind]*res[1][i][j];
					tres[2] = delt[tind]*res[2][i][j];

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
								  	 
					dur[2][i][j] = (dcrd[1][1][i][j]*tres[0] -dcrd[0][1][i][j]*tres[1]);
					dus[2][i][j] = (-dcrd[1][0][i][j]*tres[0] +dcrd[0][0][i][j]*tres[1]);
				}
			}

			b.intgrtrs(e00,e01,lf[0]);
			b.intgrtrs(e10,e11,lf[1]);
			b.intgrtrs(dur[2],dus[2],lf[2]);

			for(n=0;n<NV;++n)
				for(i=0;i<tm;++i)
					lf[n][i] *= beta[stage];
               
         lftog(tind,vvfwk,svfwk,ivfwk);
		}				
	}
   
/*	ADD IN VISCOUS/DISSIPATIVE FLUX */
	for(i=0;i<nvrtx;++i)
		for(n=0;n<NV;++n)
			vgf[i][n] += vvfwk[i][n];

 	for(i=0;i<nside*b.sm;++i)
		for(n=0;n<NV;++n)
			sgf[i][n] += svfwk[i][n];        

	for(i=0;i<ntri*b.im;++i)
		for(n=0;n<NV;++n)
			igf[i][n] += ivfwk[i][n];         

/*	ADD IN BOUNDARY SOURCES */
	flowsrc();
	
/*********************************************/
/*	MODIFY RESIDUALS ON COARSER MESHES			*/
/*********************************************/	
	if(mgrid) {
/* 	CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
		if(isfrst) {
         for(i=0;i<nvrtx;++i)
            for(n=0;n<NV;++n)
               vdg[i][n] = ffadd*vgf0[i][n] - vgf[i][n];
      
         for(i=0;i<nside*b.sm;++i)
            for(n=0;n<NV;++n)
               sdg[i][n] = ffadd*sgf0[i][n] - sgf[i][n];        
      
         for(i=0;i<ntri*b.im;++i)
            for(n=0;n<NV;++n)
               idg[i][n] = ffadd*sgf0[i][n] - sgf[i][n];
		}

		for(i=0;i<nvrtx;++i) 
			for(n=0;n<NV;++n)		
				vgf[i][n] += vdg[i][n];	
            
		for(i=0;i<nside*b.sm;++i) 
			for(n=0;n<NV;++n)		
				sgf[i][n] += sdg[i][n];
            
		for(i=0;i<ntri*b.im;++i) 
			for(n=0;n<NV;++n)		
				sgf[i][n] += sdg[i][n];      
	}

/*	APPLY BOUNDARY CONDITIONS */
	flowbnds(); 
		
	return;
}
   
/* APPLY DIRICHLET BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & FIXX_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            res[svrtx[sind][0]][0] = 0.0;
            res[svrtx[sind][1]][0] = 0.0;
         }
      }
      if (sbdry[i].type & FIXY_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            res[svrtx[sind][0]][1] = 0.0;
            res[svrtx[sind][1]][1] = 0.0;
         }
      }
   }

/* SEND COMMUNICATION PACKETS */   
   send(XDIR_MP, (FLT *) res, 0, 1, 2);
   
   return;
}


