/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "hp_mgrid.h"
   
extern FLT forcing(FLT x,FLT y);
extern FLT axext, ayext, nuext;

void hp_mgrid::rsdl(int stage, int mgrid) {
   int i,j,n,tind;
   FLT fluxx,fluxy;
   FLT visc[ND][ND][ND][ND], tres[NV];
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->res.v[i][n] = 0.0;

    for(i=0;i<nside*b->sm;++i)
      for(n=0;n<NV;++n)
         gbl->res.s[i][n] = 0.0;
                 
    for(i=0;i<ntri*b->im;++i)
      for(n=0;n<NV;++n)
         gbl->res.i[i][n] = 0.0;
         
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->vf.v[i][n] *= (1. -beta[stage]);

    for(i=0;i<nside*b->sm;++i)
      for(n=0;n<NV;++n)
         gbl->vf.s[i][n]  *= (1. -beta[stage]);
                 
    for(i=0;i<ntri*b->im;++i)
      for(n=0;n<NV;++n)
         gbl->vf.i[i][n]  *= (1. -beta[stage]);          

   for(tind = 0; tind<ntri;++tind) {
   
      if (tinfo[tind] > -1) {
         crdtocht(tind);
         for(n=0;n<ND;++n)
            b->proj_bdry(cht[n], crd[n], dcrd[n][0], dcrd[n][1],MXGP);
            
         crdtocht(tind,dvrtdt,gbl->dbinfodt);
         for(n=0;n<ND;++n)
            b->proj_bdry(cht[n],u[n],MXGP);
      }
      else {
         for(n=0;n<ND;++n)
            b->proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n],MXGP);
            
         for(n=0;n<ND;++n)
            b->proj(dvrtdt[tvrtx[tind][0]][n],dvrtdt[tvrtx[tind][1]][n],dvrtdt[tvrtx[tind][2]][n],u[n],MXGP);
            
         for(i=0;i<b->gpx;++i) {
            for(j=0;j<b->gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd[n][0][i][j] = 0.5*(vrtx[tvrtx[tind][2]][n] -vrtx[tvrtx[tind][1]][n]);
                  dcrd[n][1][i][j] = 0.5*(vrtx[tvrtx[tind][0]][n] -vrtx[tvrtx[tind][1]][n]);
               }
            }
         }
      }

      /* CALCULATE MESH VELOCITY */
      for(i=0;i<b->gpx;++i)
         for(j=0;j<b->gpn;++j)
            for(n=0;n<ND;++n)
               mvel[n][i][j] = bd[0]*crd[n][i][j] +u[n][i][j];

      ugtouht(tind);
      b->proj(uht[0],u[0],du[0][0],du[0][1],MXGP);


      for(n=0;n<NV;++n)
         for(i=0;i<b->tm;++i)
            lf[n][i] = 0.0;

      /* CONVECTION */
      for(i=0;i<b->gpx;++i) {
         for(j=0;j<b->gpn;++j) {

            fluxx = RAD(i,j)*(axext -mvel[0][i][j])*u[0][i][j];
            fluxy = RAD(i,j)*(ayext -mvel[1][i][j])*u[0][i][j];

            cv00[i][j] = +dcrd[1][1][i][j]*fluxx -dcrd[0][1][i][j]*fluxy;
            cv01[i][j] = -dcrd[1][0][i][j]*fluxx +dcrd[0][0][i][j]*fluxy;
         }
      }
      b->intgrtrs(lf[0],cv00,cv01,MXGP);

      /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
      lftog(tind,gbl->res);

      /* NEGATIVE REAL TERMS */
      if (beta[stage] > 0.0) {
            
         /* TIME DERIVATIVE TERMS */
         if (!mgrid) {         
            for(i=0;i<b->gpx;++i) {
               for(j=0;j<b->gpn;++j) {
                  cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
                  res[0][i][j] = RAD(i,j)*bd[0]*u[0][i][j]*cjcb[i][j] +gbl->dugdt[0][tind][i][j];
                  res[0][i][j] += RAD(i,j)*cjcb[i][j]*forcing(crd[0][i][j],crd[1][i][j]);
               }
            }
         }
         else {
            for(i=0;i<b->gpx;++i) {
               for(j=0;j<b->gpn;++j) {
                  cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
                  res[0][i][j] = RAD(i,j)*bd[0]*u[0][i][j]*cjcb[i][j];
               }
            }
         }
         b->intgrt(lf[0],res[0],MXGP);

         /* DIFFUSIVE TERMS  */
         for(i=0;i<b->gpx;++i) {
            for(j=0;j<b->gpn;++j) {

               cjcb[i][j] = nuext*RAD(i,j)/cjcb[i][j];
               
               /* DIFFUSION TENSOR (LOTS OF SYMMETRY THOUGH)*/
               /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
               visc[0][0][0][0] = -cjcb[i][j]*(dcrd[1][1][i][j]*dcrd[1][1][i][j] +dcrd[0][1][i][j]*dcrd[0][1][i][j]);
               visc[0][0][1][1] = -cjcb[i][j]*(dcrd[1][0][i][j]*dcrd[1][0][i][j] +dcrd[0][0][i][j]*dcrd[0][0][i][j]);
               visc[0][0][0][1] =  cjcb[i][j]*(dcrd[1][1][i][j]*dcrd[1][0][i][j] +dcrd[0][1][i][j]*dcrd[0][0][i][j]);
#define        viscI0II0II1II0I visc[0][0][0][1]

               e00[i][j] = +visc[0][0][0][0]*du[0][0][i][j] +visc[0][0][0][1]*du[0][1][i][j];
               e01[i][j] = +viscI0II0II1II0I*du[0][0][i][j] +visc[0][0][1][1]*du[0][1][i][j];
                           
               cv00[i][j] += e00[i][j];
               cv01[i][j] += e01[i][j];
             }
         }
         b->derivr(cv00,res[0],MXGP);
         b->derivs(cv01,res[0],MXGP);

         /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
         for(i=0;i<b->gpx;++i) {
            for(j=0;j<b->gpn;++j) {

               tres[0] = gbl->tau[tind]*res[0][i][j];

               e00[i][j] -= (dcrd[1][1][i][j]*(axext-mvel[0][i][j])
                            -dcrd[0][1][i][j]*(ayext-mvel[1][i][j]))*tres[0];
               e01[i][j] -= (-dcrd[1][0][i][j]*(axext-mvel[0][i][j])
                             +dcrd[0][0][i][j]*(ayext-mvel[1][i][j]))*tres[0];
           }
         }
         b->intgrtrs(lf[0],e00,e01,MXGP); 

         for(n=0;n<NV;++n)
            for(i=0;i<b->tm;++i)
               lf[n][i] *= beta[stage];
               
         lftog(tind,gbl->vf);
      }
   }

   /* ADD IN VISCOUS/DISSIPATIVE FLUX */
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->res.v[i][n] += gbl->vf.v[i][n];

    for(i=0;i<nside*b->sm;++i)
      for(n=0;n<NV;++n)
         gbl->res.s[i][n] += gbl->vf.s[i][n];        

   for(i=0;i<ntri*b->im;++i)
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
      
         for(i=0;i<nside*b->sm;++i)
            for(n=0;n<NV;++n)
               dres[log2p].s[i][n] = fadd*gbl->res0.s[i][n] -gbl->res.s[i][n];        
      
         for(i=0;i<ntri*b->im;++i)
            for(n=0;n<NV;++n)
               dres[log2p].i[i][n] = fadd*gbl->res0.i[i][n] -gbl->res.s[i][n];
               
         isfrst = false;
      }

      for(i=0;i<nvrtx;++i) 
         for(n=0;n<NV;++n)      
            gbl->res.v[i][n] += dres[log2p].v[i][n];   
            
      for(i=0;i<nside*b->sm;++i) 
         for(n=0;n<NV;++n)      
            gbl->res.s[i][n] += dres[log2p].s[i][n];
            
      for(i=0;i<ntri*b->im;++i) 
         for(n=0;n<NV;++n)      
            gbl->res.s[i][n] += dres[log2p].i[i][n];  
   }

   return;
}