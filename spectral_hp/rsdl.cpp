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
   int i,j,n,tind;
   FLT rho;
   FLT fluxx,fluxy;
   FLT visc[ND][ND][ND][ND], tres[NV];
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->res.v[i][n] = 0.0;

    for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         gbl->res.s[i][n] = 0.0;
                 
    for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         gbl->res.i[i][n] = 0.0;  
         
    for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->vf.v[i][n] *= (1. -beta[stage]);

    for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         gbl->vf.s[i][n]  *= (1. -beta[stage]);
                 
    for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         gbl->vf.i[i][n]  *= (1. -beta[stage]);          

   for(tind = 0; tind<ntri;++tind) {
   
      if (tinfo[tind] > -5) {
         crdtouht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(uht[n], crd[n], dcrd[n][0], dcrd[n][1]);
            
         crdtouht(tind,dvrtdt,gbl->dbinfodt);
         for(n=0;n<ND;++n)
            b.proj_bdry(uht[n],u[n]);
      }
      else {
         for(n=0;n<ND;++n)
            b.proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);
            
         for(n=0;n<ND;++n)
            b.proj(dvrtdt[tvrtx[tind][0]][n],dvrtdt[tvrtx[tind][1]][n],dvrtdt[tvrtx[tind][2]][n],u[n]);
            
         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd[n][0][i][j] = 0.5*(vrtx[tvrtx[tind][1]][n] -vrtx[tvrtx[tind][0]][n]);
                  dcrd[n][1][i][j] = 0.5*(vrtx[tvrtx[tind][2]][n] -vrtx[tvrtx[tind][0]][n]);
               }
            }
         }
      }

      /* CALCULATE MESH VELOCITY */
      for(i=0;i<b.gpx;++i)
         for(j=0;j<b.gpn;++j)
            for(n=0;n<ND;++n)
               crd[n][i][j] = bd[0]*crd[n][i][j] +u[n][i][j];

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

      for(n=0;n<NV;++n)
         for(i=0;i<b.tm;++i)
            lf[n][i] = 0.0;

      /* CONVECTIVE TERMS (IMAGINARY FIRST)*/
      for(i=0;i<b.gpx;++i) {
         for(j=0;j<b.gpn;++j) {

            fluxx = gbl->rho*(u[0][i][j] -crd[0][i][j]);
            fluxy = gbl->rho*(u[1][i][j] -crd[1][i][j]);

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

      /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
      lftog(tind,gbl->res);

      /* NEGATIVE REAL TERMS */
      if (beta[stage] > 0.0) {
            
         /* TIME DERIVATIVE TERMS */
         if (!mgrid) {         
            for(i=0;i<b.gpx;++i) {
               for(j=0;j<b.gpn;++j) {
                  cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
                  rho = gbl->rho*(1 +0.1*u[2][i][j]);
                  res[0][i][j] = gbl->rho*bd[0]*u[0][i][j]*cjcb[i][j] +gbl->dugdt[0][tind][i][j];
                  res[1][i][j] = gbl->rho*bd[0]*u[1][i][j]*cjcb[i][j] +gbl->dugdt[1][tind][i][j];
                  res[2][i][j] = rho*bd[0]*cjcb[i][j] +gbl->dugdt[2][tind][i][j];
               }
            }
         }
         else {
            for(i=0;i<b.gpx;++i) {
               for(j=0;j<b.gpn;++j) {
                  cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
                  res[0][i][j] = gbl->rho*bd[0]*u[0][i][j]*cjcb[i][j];
                  res[1][i][j] = gbl->rho*bd[0]*u[1][i][j]*cjcb[i][j];
                  res[2][i][j] = gbl->rho*bd[0]*cjcb[i][j];
               }
            }
         }
         b.intgrt(res[0],lf[0]);
         b.intgrt(res[1],lf[1]);
         b.intgrt(res[2],lf[2]);

         /* VISCOUS TERMS  */
         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {

               cjcb[i][j] = gbl->mu/cjcb[i][j];
               
               /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
               /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
               visc[0][0][0][0] =  cjcb[i][j]*(2.*dcrd[1][1][i][j]*dcrd[1][1][i][j] +dcrd[0][1][i][j]*dcrd[0][1][i][j]);
               visc[0][0][1][1] =  cjcb[i][j]*(2.*dcrd[1][0][i][j]*dcrd[1][0][i][j] +dcrd[0][0][i][j]*dcrd[0][0][i][j]);
               visc[0][0][0][1] = -cjcb[i][j]*(2.*dcrd[1][1][i][j]*dcrd[1][0][i][j] +dcrd[0][1][i][j]*dcrd[0][0][i][j]);
#define           viscI0II0II1II0I visc[0][0][0][1]

               visc[1][1][0][0] =  cjcb[i][j]*(dcrd[1][1][i][j]*dcrd[1][1][i][j] +2.*dcrd[0][1][i][j]*dcrd[0][1][i][j]);
               visc[1][1][1][1] =  cjcb[i][j]*(dcrd[1][0][i][j]*dcrd[1][0][i][j] +2.*dcrd[0][0][i][j]*dcrd[0][0][i][j]);
               visc[1][1][0][1] = -cjcb[i][j]*(dcrd[1][1][i][j]*dcrd[1][0][i][j] +2.*dcrd[0][1][i][j]*dcrd[0][0][i][j]);
#define           viscI1II1II1II0I visc[1][1][0][1]
               
               visc[0][1][0][0] = -cjcb[i][j]*dcrd[0][1][i][j]*dcrd[1][1][i][j];
               visc[0][1][1][1] = -cjcb[i][j]*dcrd[0][0][i][j]*dcrd[1][0][i][j];
               visc[0][1][0][1] =  cjcb[i][j]*dcrd[0][1][i][j]*dcrd[1][0][i][j];
               visc[0][1][1][0] =  cjcb[i][j]*dcrd[0][0][i][j]*dcrd[1][1][i][j];

                  /* OTHER SYMMETRIES    */            
#define           viscI1II0II0II0I visc[0][1][0][0]
#define           viscI1II0II1II1I visc[0][1][1][1]
#define           viscI1II0II0II1I visc[0][1][1][0]
#define           viscI1II0II1II0I visc[0][1][0][1]
               


               e00[i][j] = +visc[0][0][0][0]*du[0][0][i][j] +visc[0][1][0][0]*du[1][0][i][j]
                           +visc[0][0][0][1]*du[0][1][i][j] +visc[0][1][0][1]*du[1][1][i][j];

               e01[i][j] = +viscI0II0II1II0I*du[0][0][i][j] +visc[0][1][1][0]*du[1][0][i][j]
                           +visc[0][0][1][1]*du[0][1][i][j] +visc[0][1][1][1]*du[1][1][i][j];

               e10[i][j] = +viscI1II0II0II0I*du[0][0][i][j] +visc[1][1][0][0]*du[1][0][i][j]
                           +viscI1II0II0II1I*du[0][1][i][j] +visc[1][1][0][1]*du[1][1][i][j];

               e11[i][j] = +viscI1II0II1II0I*du[0][0][i][j] +viscI1II1II1II0I*du[1][0][i][j]
                           +viscI1II0II1II1I*du[0][1][i][j] +visc[1][1][1][1]*du[1][1][i][j];
                           
               cv00[i][j] = -cv00[i][j] -e00[i][j];
               cv01[i][j] = -cv01[i][j] -e01[i][j];
               cv10[i][j] = -cv10[i][j] -e10[i][j];
               cv11[i][j] = -cv11[i][j] -e11[i][j];
               
               du[2][0][i][j] *= -1;
               du[2][1][i][j] *= -1;
               
#ifdef OLDWAYS
   /* THESE DON'T WORK FOR MOVING MESH RIGHT NOW */
                  res[2][i][j] += gbl->rho*
                                 (dcrd[1][1][i][j]*(du[0][0][i][j] -dlmvx[0][i][j])
                                 -dcrd[0][1][i][j]*(du[1][0][i][j] -dlmvx[1][i][j])
                                 -dcrd[1][0][i][j]*(du[0][1][i][j] -dlmvy[0][i][j])
                                 +dcrd[0][0][i][j]*(du[1][1][i][j] -dlmvy[1][i][j]));

                 res[2][i][j] += gbl->rho*(dcrd[1][1][i][j]*du[0][0][i][j] -dcrd[0][1][i][j]*du[1][0][i][j]
                                          -dcrd[1][0][i][j]*du[0][1][i][j] +dcrd[0][0][i][j]*du[1][1][i][j]);
#endif
            
            }
         }
         b.derivr(cv00,res[0]);
         b.derivs(cv01,res[0]);
         b.derivr(cv10,res[1]);
         b.derivs(cv11,res[1]);
         b.derivr(du[2][0],res[2]);
         b.derivs(du[2][1],res[2]);
         

         /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {

               tres[0] = gbl->tau[tind]*res[0][i][j];
               tres[1] = gbl->tau[tind]*res[1][i][j];
               tres[2] = gbl->delt[tind]*res[2][i][j];

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
               
         lftog(tind,gbl->vf);
      }
   }

   /*
   for(i=0;i<nvrtx;++i)
      printf("%d %f %f %f\n",i,gbl->res.v[i][0],gbl->res.v[i][1],gbl->res.v[i][2]);
      
   for(i=0;i<nside*b.sm;++i)
      printf("%d %f %f %f\n",i,gbl->res.s[i][0],gbl->res.s[i][1],gbl->res.s[i][2]);

   for(i=0;i<ntri*b.im;++i)
      printf("%d %f %f %f\n",i,gbl->res.i[i][0],gbl->res.i[i][1],gbl->res.i[i][2]);      
   */      

   /* ADD IN VISCOUS/DISSIPATIVE FLUX */
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->res.v[i][n] += gbl->vf.v[i][n];

    for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         gbl->res.s[i][n] += gbl->vf.s[i][n];        

   for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         gbl->res.i[i][n] += gbl->vf.i[i][n];         
      
   /* ADD IN BOUNDARY FLUXES */
   addbflux(mgrid);

/*********************************************/
   /* MODIFY RESIDUALS ON COARSER MESHES         */
/*********************************************/   
   if(mgrid) {
   /* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
      if(isfrst) {
         for(i=0;i<nvrtx;++i)
            for(n=0;n<NV;++n)
               dres[log2p].v[i][n] = fadd*gbl->res0.v[i][n] -gbl->res.v[i][n];
      
         for(i=0;i<nside*b.sm;++i)
            for(n=0;n<NV;++n)
               dres[log2p].s[i][n] = fadd*gbl->res0.s[i][n] -gbl->res.s[i][n];        
      
         for(i=0;i<ntri*b.im;++i)
            for(n=0;n<NV;++n)
               dres[log2p].i[i][n] = fadd*gbl->res0.i[i][n] -gbl->res.s[i][n];
               
         isfrst = false;
      }

      for(i=0;i<nvrtx;++i) 
         for(n=0;n<NV;++n)      
            gbl->res.v[i][n] += dres[log2p].v[i][n];   
            
      for(i=0;i<nside*b.sm;++i) 
         for(n=0;n<NV;++n)      
            gbl->res.s[i][n] += dres[log2p].s[i][n];
            
      for(i=0;i<ntri*b.im;++i) 
         for(n=0;n<NV;++n)      
            gbl->res.s[i][n] += dres[log2p].i[i][n];  
   }

   /*
   for(i=0;i<nvrtx;++i)
      printf("v: %d %f %f %f\n",i,gbl->res.v[i][0],gbl->res.v[i][1],gbl->res.v[i][2]);
      
   for(i=0;i<nside*b.sm;++i)
      printf("s: %d %f %f %f\n",i,gbl->res.s[i][0],gbl->res.s[i][1],gbl->res.s[i][2]);

   for(i=0;i<ntri*b.im;++i)
      printf("i: %d %f %f %f\n",i,gbl->res.i[i][0],gbl->res.i[i][1],gbl->res.i[i][2]); 
        
   exit(1);
   */

   return;
}
