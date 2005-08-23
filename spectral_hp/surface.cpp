#include"hp_mgrid.h"
#include<assert.h>
#include<utilities.h>
#include<myblas.h>
#include<stdio.h>

FLT surface::fadd[ND];
FLT surface::cfl[MXLG2P][ND];

extern FLT body[ND];

void surface::alloc(int maxside, int log2p, int mgrid, int fmesh, struct surface_glbls *store) {
   int i,p;
   
   p = 1;
   for (i=0;i<log2p;++i)
      p = p<<1;

   if (!mgrid) {
      gbl_alloc(maxside, p, store);
      vect_alloc(ksprg,maxside,FLT);
      vug = (FLT (*)[ND]) xmalloc((maxside+1)*ND*sizeof(FLT));
      sug = (FLT (*)[ND]) xmalloc(maxside*(p-1)*ND*sizeof(FLT));
   }
   vdres[log2p] = (FLT (*)[ND]) xmalloc(ND*(maxside+1)*sizeof(FLT));
   sdres[log2p] = (FLT (*)[ND]) xmalloc(ND*maxside*(p-1)*sizeof(FLT)); 
   gbl = store;         
   
   if (!fmesh) {
      vect_alloc(ksprg,maxside,FLT);
      vug = (FLT (*)[ND]) xmalloc((maxside+1)*ND*sizeof(FLT));
      vug_frst = (FLT (*)[ND]) xmalloc(ND*(maxside+1)*sizeof(FLT));
   }

   return;
}

void surface::gbl_alloc(int maxside, int p, struct surface_glbls *store) {
   
   /* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   store->vug0 = (FLT (*)[ND]) xmalloc((maxside+1)*ND*sizeof(FLT));
   store->sug0 = (FLT (*)[ND]) xmalloc(maxside*(p-1)*ND*sizeof(FLT));

   /* RESIDUALS */
   store->vres = (FLT (*)[ND]) xmalloc((maxside+1)*ND*sizeof(FLT));
   store->sres = (FLT (*)[ND]) xmalloc(maxside*(p-1)*ND*sizeof(FLT));
   store->vres0 = (FLT (*)[ND]) xmalloc((maxside+1)*ND*sizeof(FLT));
   store->sres0 = (FLT (*)[ND]) xmalloc(maxside*(p-1)*ND*sizeof(FLT));

   /* PSEUDO TIME ITERATION */
   store->vdt = (FLT (*)[ND][ND]) xmalloc((maxside+1)*ND*ND*sizeof(FLT));
   store->sdt = (FLT (*)[ND][ND]) xmalloc(maxside*ND*ND*sizeof(FLT));
   vect_alloc(store->dtfnrm,maxside,FLT);
   vect_alloc(store->normc,maxside,FLT);
   vect_alloc(store->meshc,maxside,FLT);
   
   return;
}

void hp_mgrid::setksprg1d() {
   int i,j,sind;
   class surface *srf;

   for(i=0;i<nsbd;++i) {
      if(sbdry[i].type&(FSRF_MASK+IFCE_MASK)) {
         srf = static_cast<class surface *>(sbdry[i].misc);
         if (srf != NULL) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               srf->ksprg[j] = 1.0/distance(svrtx[sind][0],svrtx[sind][1]);
            }
         }
      }
   }
   
   
   return;
}
void hp_mgrid::surfksrc1d() {
   int i,bnum;
   class surface *srf;   

   for(bnum=0;bnum<nsbd;++bnum) {
      if(sbdry[bnum].type&(FSRF_MASK+IFCE_MASK)) {
         srf = static_cast<class surface *>(sbdry[bnum].misc);
         if (srf != NULL) {
            /* ZERO TANGENTIAL MESH MOVEMENT SOURCE */   
            for(i=0;i<sbdry[bnum].num+1;++i)
               srf->vdres[log2p][i][0] = 0.0;

            for(i=0;i<sbdry[bnum].num*b->sm;++i) 
               srf->sdres[log2p][i][0] = 0.0;
         }
      }
   }
   surfrsdl(0);
   
   for(bnum=0;bnum<nsbd;++bnum) {
      if(sbdry[bnum].type&(FSRF_MASK+IFCE_MASK)) {
         srf = static_cast<class surface *>(sbdry[bnum].misc);
         if (srf != NULL) {         
            /* TANGENTIAL MESH MOVEMENT SOURCE */   
            for(i=0;i<sbdry[bnum].num+1;++i) {
               srf->vdres[log2p][i][0] = -srf->gbl->vres[i][0];
            }

            for(i=0;i<sbdry[bnum].num*b->sm;++i) {
               srf->sdres[log2p][i][0] = -0.0*srf->gbl->sres[i][0];  // TEMPO FOR UNIFORM SPACING OF HIGHER MODES ALONG SIDE
            }
         }
      }
   }
   
   
   return;
}

#ifdef DROP
extern FLT amp;
FLT dydt = 0.0;
#endif

void hp_mgrid::surfrsdl(int mgrid) {
   int i,m,n,sind,indx,indx1,count,v0,bnum;
   FLT ubar[ND], norm[ND], rp[ND], jcb;
   FLT sigor, drhor;
   class surface *srf;


   for(bnum=0;bnum<nsbd;++bnum) {
      if (!(sbdry[bnum].type&(FSRF_MASK +IFCE_MASK))) continue;

      srf = static_cast<class surface *>(sbdry[bnum].misc);
      if (srf == NULL) continue;
      
      count = 0;
      sigor = srf->gbl->sigma/(gbl->rho +srf->gbl->rho2);
      drhor = (gbl->rho -srf->gbl->rho2)/(gbl->rho +srf->gbl->rho2);
      

#ifdef DROP
      /* DETRMINE DX CORRECTION TO CONSERVE AREA */
      /* IMPORTANT FOR STEADY SOLUTIONS */
      /* SINCE THERE ARE MULTIPLE STEADY-STATES */
      /* TO ENSURE GET CORRECT VOLUME */
      FLT avg[5],rbar, vflux, kc; 
      kc = srf->gbl->sigma/(gbl->mu +srf->gbl->mu2);
      integrated_averages(avg);
      rbar  = pow(3.*0.5*avg[0],1.0/3.0);
      vflux =  amp*kc*(rbar -0.5);
      dydt = amp*kc*avg[2] +avg[4];
      /* C_D TO G CONVERSION REMINDER 
      re = 1.0/srf->gbl->mu2;
      cd = 24./re*(1 +0.1935*pow(re,0.6305));
      cd /= 16.0; // (1/2 rho u^2 * Pi r^2 / 2 pi);
      g = amp*(avg[2] +avg[4]) +12.*cd/(gbl->rho -srf->gbl->rho2);
      */
#endif
            
      /**************************************************/
      /* DETERMINE MESH RESIDUALS & SURFACE TENSION     */
      /**************************************************/
      for(n=0;n<ND;++n)
         srf->gbl->vres[0][n] = 0.0;
    
      for(n=0;n<NV;++n) 
         binfo[bnum][0].flx[n] = 0.0;

      for(indx=0;indx<sbdry[bnum].num;++indx) {
         sind = sbdry[bnum].el[indx];

         /* UG STORES THE CURRENT BOUNDARY POSITION */
         for(n=0;n<ND;++n) {
            uht[n][0] = srf->vug[indx][n];
            uht[n][1] = srf->vug[indx+1][n];
         }
         if (b->sm > 0) {
            indx1 = indx*sm0;
            for(m=0;m<b->sm;++m)
               for(n=0;n<ND;++n)
                  uht[n][m+2] = srf->sug[indx1 +m][n];
         }
         for(n=0;n<ND;++n)
            b->proj1d(uht[n],crd[n][0],dcrd[n][0][0]);
            
         crdtocht1d(sind,dvrtdt,gbl->dbinfodt);
         for(n=0;n<ND;++n)
            b->proj1d(cht[n],mvel[n][0]);

         ugtouht1d(sind);
         for(n=0;n<ND;++n)
            b->proj1d(uht[n],u[n][0]);   
         
         for(i=0;i<b->gpx;++i) {
            norm[0] =  dcrd[1][0][0][i];
            norm[1] = -dcrd[0][0][0][i];
            jcb = sqrt(norm[0]*norm[0] +norm[1]*norm[1]);
            
            /* RELATIVE VELOCITY STORED IN MVEL[N][0]*/
            for(n=0;n<ND;++n)
               mvel[n][0][i] = u[n][0][i] -(bd[0]*crd[n][0][i] +mvel[n][0][i]);
               
#ifdef DROP
            mvel[1][0][i] -= dydt;
#endif
            /* TANGENTIAL SPACING & NORMAL FLUX */            
            res[0][0][i] = -srf->ksprg[indx]*jcb;
            res[1][0][i] = -RAD1D(i)*(mvel[0][0][i]*norm[0] +mvel[1][0][i]*norm[1]);
            res[1][1][i] = -res[1][0][i]*(-norm[1]*mvel[0][0][i] +norm[0]*mvel[1][0][i])/jcb*srf->gbl->meshc[indx];
#ifdef DROP
            res[2][0][i] = +RAD1D(i)*vflux*jcb;
#endif
            
            /* SURFACE TENSION SOURCE TERM */ 
            u[0][0][i] = +RAD1D(i)*srf->gbl->sigma*norm[1]/jcb;
            u[0][1][i] = +RAD1D(i)*(gbl->rho -srf->gbl->rho2)*g*crd[1][0][i]*norm[0];
#ifdef AXISYMMETRIC
            u[0][1][i] += srf->gbl->sigma*jcb;
#endif
            u[1][0][i] = -RAD1D(i)*srf->gbl->sigma*norm[0]/jcb;
            u[1][1][i] = +RAD1D(i)*(gbl->rho -srf->gbl->rho2)*g*crd[1][0][i]*norm[1];            
         }
         
         for(m=0;m<b->sm+2;++m)
            cf[0][m] = 0.0;

         /* INTEGRATE & STORE SURFACE RESIDUALS */               
         b->intgrtx1d(cf[0],res[0][0]);
         b->intgrt1d(cf[1],res[1][0]);
         b->intgrtx1d(cf[1],res[1][1]);
         
#ifdef DROP
         b->intgrt1d(lf[2],res[2][0]);
#endif
         
         /* TO LEAVE TANGENTIAL POSITION TOTALLY FREE */
         /* for(m=0;m<b->sm+2;++m)
            cf[0][m] = 0.0; */
         
         /* STORE IN RES */
         for(n=0;n<ND;++n) {
            srf->gbl->vres[indx][n] += cf[n][0];
            srf->gbl->vres[indx+1][n] = cf[n][1];
            for(m=0;m<b->sm;++m)
               srf->gbl->sres[b->sm*indx +m][n] = cf[n][m+2];
         }
         
#ifdef DROP
         srf->gbl->vres[indx][1] += lf[2][0];
         srf->gbl->vres[indx+1][1] += lf[2][1];
         for(m=0;m<b->sm;++m)
            srf->gbl->sres[b->sm*indx +m][1] += lf[2][m+2];
#endif
         
         /* INTEGRATE & STORE SURFACE TENSION SOURCE TERM */
         b->intgrt1d(lf[0],u[0][1]);
         b->intgrtx1d(lf[0],u[0][0]);
         b->intgrt1d(lf[1],u[1][1]);
         b->intgrtx1d(lf[1],u[1][0]);

         /* MASS FLUX PRECONDITIONER */
         for(m=0;m<b->sm+2;++m)
            lf[2][m] = -gbl->rho*cf[1][m]; 

#ifndef INERTIALESS
         for (n=0;n<ND;++n) 
            ubar[n] = 0.5*(uht[n][0] +uht[n][1]);
            
         for (n=0;n<ND;++n) {
            lf[n][0] -= uht[n][0]*(gbl->rho -srf->gbl->rho2)*cf[1][0];
            lf[n][1] -= uht[n][1]*(gbl->rho -srf->gbl->rho2)*cf[1][1];
            for(m=0;m<b->sm;++m)
               lf[n][m+2] -= ubar[n]*(gbl->rho -srf->gbl->rho2)*cf[1][m+2];
         }
#endif

         /* STORE IN BINFO.FLUX */
         for(n=0;n<NV;++n) 
            binfo[bnum][count].flx[n] += lf[n][0];
         ++count;
         
         for(m=0;m<b->sm;++m) {
            for(n=0;n<NV;++n)
               binfo[bnum][count].flx[n] = lf[n][m+2];
            ++count;
         }
         
         for(n=0;n<NV;++n)
            binfo[bnum][count].flx[n] = lf[n][1];
      }
      
      /* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
      v0 = svrtx[sbdry[bnum].el[0]][0];
      for(i=0;i<nvbd;++i) {
         if (vbdry[i].el[0] != v0) continue;
         
          if  (vbdry[i].type&OUTF_MASK) {

            /* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
            crdtocht1d(sbdry[bnum].el[0]);
            b->ptprobe1d(2,rp,norm,-1.0,cht[0],MXTM);
            jcb = sqrt(norm[0]*norm[0] +norm[1]*norm[1]);
#ifdef AXISYMMETRIC
            binfo[bnum][0].flx[0] -= -rp[0]*srf->gbl->sigma*norm[0]/jcb;
            binfo[bnum][0].flx[1] -= -rp[0]*srf->gbl->sigma*norm[1]/jcb;
#else
            binfo[bnum][0].flx[0] -= -srf->gbl->sigma*norm[0]/jcb;
            binfo[bnum][0].flx[1] -= -srf->gbl->sigma*norm[1]/jcb;
#endif
         }
         if (vbdry[i].type&PRDX_MASK && vbdry[i].type&PRDY_MASK) {
            binfo[bnum][0].flx[2] = 0.0;
         }
      }
      
      /* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
      v0 = svrtx[sbdry[bnum].el[sbdry[bnum].num -1]][1];
      for(i=0;i<nvbd;++i) {
         if (vbdry[i].el[0] != v0) continue;
         
         if  (vbdry[i].type&OUTF_MASK) {

            /* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
            crdtocht1d(sbdry[bnum].el[sbdry[bnum].num -1]);
            b->ptprobe1d(2,rp,norm,1.0,cht[0],MXTM);
            jcb = sqrt(norm[0]*norm[0] +norm[1]*norm[1]);
#ifdef AXISYMMETRIC
            binfo[bnum][count].flx[0] += -rp[0]*srf->gbl->sigma*norm[0]/jcb;
            binfo[bnum][count].flx[1] += -rp[0]*srf->gbl->sigma*norm[1]/jcb;
#else
            binfo[bnum][count].flx[0] += -srf->gbl->sigma*norm[0]/jcb;
            binfo[bnum][count].flx[1] += -srf->gbl->sigma*norm[1]/jcb;
#endif
         }
         if (vbdry[i].type&PRDX_MASK && vbdry[i].type&PRDY_MASK) {
            binfo[bnum][count].flx[2] = 0.0;
         }
      }   
         
      /************************************************/
      /* MODIFY SURFACE RESIDUALS ON COARSER MESHES   */
      /************************************************/   
      if(mgrid) {
         if (isfrst) {
            for(i=0;i<sbdry[bnum].num+1;++i) {
               for(n=0;n<ND;++n)
                  srf->vdres[log2p][i][n] = surface::fadd[n]*srf->gbl->vres0[i][n] -srf->gbl->vres[i][n];
            }
            for(i=0;i<sbdry[bnum].num*b->sm;++i) {
               for(n=0;n<ND;++n)
                  srf->sdres[log2p][i][n] = surface::fadd[n]*srf->gbl->sres0[i][n] -srf->gbl->sres[i][n];
            }
         }
         for(i=0;i<sbdry[bnum].num+1;++i) {
            for(n=0;n<ND;++n)
               srf->gbl->vres[i][n] += srf->vdres[log2p][i][n];
         }
         for(i=0;i<sbdry[bnum].num*b->sm;++i) {
            for(n=0;n<ND;++n)      
               srf->gbl->sres[i][n] += srf->sdres[log2p][i][n];
         }
      }
      else {
         /* ADD TANGENTIAL MESH MOVEMENT SOURCE */   
         for(i=0;i<sbdry[bnum].num+1;++i)
            srf->gbl->vres[i][0] += srf->vdres[log2p][i][0];

         for(i=0;i<sbdry[bnum].num*b->sm;++i) 
            srf->gbl->sres[i][0] += srf->sdres[log2p][i][0];
      }
      
      
      if (sbdry[bnum].type&IFCE_MASK) {
      /* THIS IS TO SEND CONTINUITY FLUX FOR INTERFACE */
         class spectral_hp *tgt = static_cast<spectral_hp *>(sbdry[bnum].adjmesh);
         int vbnum = sbdry[bnum].adjbnum;
         int msgn;
                  
         int count1 = 0;
         for(i=0;i<sbdry[bnum].num;++i) {
            tgt->binfo[vbnum][count1++].flx[2] = -srf->gbl->rho2*binfo[bnum][count].flx[2]/gbl->rho;
                     
            count -= b->sm;
            msgn = 1;
            for(m=0;m<b->sm;++m) {
               tgt->binfo[vbnum][count1++].flx[2] = -srf->gbl->rho2*msgn*binfo[bnum][count++].flx[2]/gbl->rho;
               msgn *= -1;
            }
            count -= (b->sm+1);
            
         }
         tgt->binfo[vbnum][count1++].flx[2] = -srf->gbl->rho2*binfo[bnum][0].flx[2]/gbl->rho;
         
         for(i=0;i<count1;++i)
            for (n=0;n<ND;++n)
               tgt->binfo[vbnum][i].flx[n] = 0.0;
      }

//      for(i=0;i<sbdry[bnum].num+1;++i)
//         printf("vres: %d %e %e\n",i,srf->gbl->vres[i][0],srf->gbl->vres[i][1]);
//            
//      for(i=0;i<sbdry[bnum].num*b->sm;++i)
//         printf("sres: %f %f\n",srf->gbl->sres[i][0],srf->gbl->sres[i][1]);
         
   }

   return;
}

void hp_mgrid::surfinvrt1(int bnum) {
   int i,m,n,v0,v1,indx,indx1,end,vbnum;
   class surface *srf;
   class mesh *tgt;

/* POINTER TO STUFF NEEDED FOR SURFACES IS STORED IN MISCELLANEOUS */      
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   
   if (srf == NULL) return;
   
   /* INVERT MASS MATRIX */
   /* LOOP THROUGH SIDES */
   if (b->sm > 0) {
      indx1 = 0;
      for(indx = 0; indx<sbdry[bnum].num; ++indx) {
         /* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */         
         for (m=0; m <b->sm; ++m) {
            for(n=0;n<ND;++n)
               srf->gbl->vres[indx][n] -= b->sfmv1d(0,m)*srf->gbl->sres[indx1][n];
            for(n=0;n<ND;++n)
               srf->gbl->vres[indx+1][n] -= b->sfmv1d(1,m)*srf->gbl->sres[indx1][n];
            ++indx1;
         }
      }
   }

   /* SEND COMMUNICATION INFO FOR ENDPOINTS */
   end = sbdry[bnum].num;
   v0 = svrtx[sbdry[bnum].el[0]][0];
   v1 = svrtx[sbdry[bnum].el[end-1]][1];
   if (v0 == v1) {
      /* SURFACE IS A LOOP ON SINGLE BLOCK */
      srf->gbl->vres[0][1] = 0.5*(srf->gbl->vres[0][1] +srf->gbl->vres[end][1]);
      srf->gbl->vres[end][1] = srf->gbl->vres[0][1];
      srf->gbl->vres[0][0] = 0.0;
      srf->gbl->vres[end][0] = 0.0;
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      if (vbdry[i].type&(COMX_MASK+COMY_MASK)) {
         vbnum = vbdry[i].adjbnum;
         tgt = vbdry[i].adjmesh;
         tgt->vbuff[vbnum][0] = srf->gbl->vres[0][0];
         tgt->vbuff[vbnum][1] = srf->gbl->vres[0][1];
      }
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      
      if (vbdry[i].type&(COMX_MASK+COMY_MASK)) {
         vbnum = vbdry[i].adjbnum;
         tgt = vbdry[i].adjmesh;
         tgt->vbuff[vbnum][0] = srf->gbl->vres[end][0];
         tgt->vbuff[vbnum][1] = srf->gbl->vres[end][1];
      }
   }

   return;
}

void hp_mgrid::surfinvrt2(int bnum) {
   int i,m,n,v0,v1,indx,indx1,end;
   FLT temp;
   class surface *srf;

/* POINTER TO STUFF NEEDED FOR SURFACES IS STORED IN MISCELLANEOUS */      
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;
   
/* RECEIVE MESSAGES AND ZERO VERTEX MODES */
   end = sbdry[bnum].num;
   v0 = svrtx[sbdry[bnum].el[0]][0];
   v1 = svrtx[sbdry[bnum].el[end-1]][1];
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      
      if (vbdry[i].type&(COMX_MASK +COMY_MASK)) {
         srf->gbl->vres[0][0] = 0.5*(srf->gbl->vres[0][0] +vbuff[i][0]);
         srf->gbl->vres[0][1] = 0.5*(srf->gbl->vres[0][1] +vbuff[i][1]);
      }
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
            
      if (vbdry[i].type&(PRDX_MASK +SYMM_MASK +PRDY_MASK +CURV_MASK)) {
         srf->gbl->vres[0][0] = 0.0;
      }
      if (vbdry[i].type&PRDX_MASK && vbdry[i].type&PRDY_MASK) {  // TOTALLY FIXED POINT
         srf->gbl->vres[0][1] = 0.0;
      }
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      
      if (vbdry[i].type&(COMX_MASK +COMY_MASK)) {
         srf->gbl->vres[end][0] = 0.5*(srf->gbl->vres[end][0] +vbuff[i][0]);
         srf->gbl->vres[end][1] = 0.5*(srf->gbl->vres[end][1] +vbuff[i][1]);
      }
   }
   
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
            
      if (vbdry[i].type&(PRDX_MASK +PRDY_MASK +SYMM_MASK +CURV_MASK)) {
         srf->gbl->vres[end][0] = 0.0;
      }
      if (vbdry[i].type&PRDX_MASK && vbdry[i].type&PRDY_MASK) {  // TOTALLY FIXED POINT
         srf->gbl->vres[end][1] = 0.0;
      }
   }

   /* SOLVE FOR VERTEX MODES */
   for(i=0;i<=end;++i) {
      temp                = srf->gbl->vres[i][0]*srf->gbl->vdt[i][0][0] +srf->gbl->vres[i][1]*srf->gbl->vdt[i][0][1];
      srf->gbl->vres[i][1] = srf->gbl->vres[i][0]*srf->gbl->vdt[i][1][0] +srf->gbl->vres[i][1]*srf->gbl->vdt[i][1][1];
      srf->gbl->vres[i][0] = temp;
   }
   
/* HAVE TO CORRECT AT PRDC BOUNDARIES SO POINT DOESN'T MOVE OFF LINE */
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      if (vbdry[i].type&(PRDX_MASK +SYMM_MASK)) {
         srf->gbl->vres[0][0] = 0.0;
      }
      if (vbdry[i].type&PRDY_MASK)
         srf->gbl->vres[0][1] = 0.0;
   }       
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      if (vbdry[i].type&(PRDX_MASK +SYMM_MASK)) {
         srf->gbl->vres[end][0] = 0.0;
      }
      if (vbdry[i].type&PRDY_MASK)
         srf->gbl->vres[end][1] = 0.0;
   }   
   
   /* SOLVE FOR SIDE MODES */
   if (b->sm > 0) {
      indx1 = 0;
      for(indx = 0; indx<sbdry[bnum].num; ++indx) {
         
         /* INVERT SIDE MODES */
         DPBTRSNU2(&b->sdiag1d(0,0),b->sbwth+1,b->sm,b->sbwth,&(srf->gbl->sres[indx1][0]),ND);
         for(m=0;m<b->sm;++m) {
            temp                      = srf->gbl->sres[indx1][0]*srf->gbl->sdt[indx][0][0] +srf->gbl->sres[indx1][1]*srf->gbl->sdt[indx][0][1];
            srf->gbl->sres[indx1][1] = srf->gbl->sres[indx1][0]*srf->gbl->sdt[indx][1][0] +srf->gbl->sres[indx1][1]*srf->gbl->sdt[indx][1][1];       
            srf->gbl->sres[indx1][0] = temp;
            ++indx1;
         }
         
         indx1 -= b->sm;
         for(m=0;m<b->sm;++m) {
            for(n=0;n<ND;++n) {
               srf->gbl->sres[indx1][n] -= b->vfms1d(0,m)*srf->gbl->vres[indx][n];
               srf->gbl->sres[indx1][n] -= b->vfms1d(1,m)*srf->gbl->vres[indx+1][n];
            }
            ++indx1;
         }
      }
   }
   
//   for(i=0;i<sbdry[bnum].num*b->sm;++i) 
//      printf("s: %f %f\n",srf->gbl->sres[i][0],srf->gbl->sres[i][1]);
   
         
   return;
}

void hp_mgrid::surfdt1(int bnum) {
   int i,m,n,indx,end,vbnum,count,sind,v0,v1;
   FLT nrm[ND], h, hsm;
   FLT dttang, dtnorm;
   FLT uvel, vvel, vslp, strss;
   class surface *srf;
   class mesh *tgt;
   FLT drho, srho, smu;
   FLT nu1, nu2;
   FLT qmax, gam1, gam2;
   
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;
      
   drho = gbl->rho -srf->gbl->rho2;
   srho = gbl->rho +srf->gbl->rho2;
   smu = gbl->mu +srf->gbl->mu2;
   nu1 = gbl->mu/gbl->rho;
   if (srf->gbl->rho2 > 0.0) 
      nu2 = srf->gbl->mu2/srf->gbl->rho2;
   else
      nu2 = 0.0;

   /**************************************************/
   /* DETERMINE SURFACE MOVEMENT TIME STEP           */
   /**************************************************/
   for(n=0;n<ND;++n)
      for(m=0;m<ND;++m)
         srf->gbl->vdt[0][m][n] = 0.0;
            
   for(indx=0; indx < sbdry[bnum].num; ++indx) {
      sind = sbdry[bnum].el[indx];
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];

      nrm[0] =  (vrtx[v1][1] -vrtx[v0][1]);
      nrm[1] = -(vrtx[v1][0] -vrtx[v0][0]);
      h = sqrt(nrm[0]*nrm[0] +nrm[1]*nrm[1]);
      
      
      uvel = ug.v[v0][0]-(bd[0]*vrtx[v0][0] +dvrtdt[v0][0]);
      vvel = ug.v[v0][1]-(bd[0]*vrtx[v0][1] +dvrtdt[v0][1]);  
#ifdef DROP
      vvel -= dydt;
#endif
      qmax = uvel*uvel+vvel*vvel;
      vslp = fabs(-uvel*nrm[1]/h +vvel*nrm[0]/h);

      uvel = ug.v[v1][0]-(bd[0]*vrtx[v1][0] +dvrtdt[v1][0]);
      vvel = ug.v[v1][1]-(bd[0]*vrtx[v1][1] +dvrtdt[v1][1]);
#ifdef DROP
      vvel  -= dydt;
#endif
      qmax = MAX(qmax,uvel*uvel+vvel*vvel);
      vslp = MAX(vslp,fabs(-uvel*nrm[1]/h +vvel*nrm[0]/h));
      
      hsm = h/(.25*(b->p+1)*(b->p+1));
            
      dttang = 2.*srf->ksprg[indx]*(.25*(b->p+1)*(b->p+1))/hsm;
#ifndef BODY
      strss =  4.*srf->gbl->sigma/(hsm*hsm) +fabs(drho*g*nrm[1]/h);
#else
      strss =  4.*srf->gbl->sigma/(hsm*hsm) +fabs(drho*(-body[0]*nrm[0] +(g-body[1])*nrm[1])/h);
#endif

      gam1 = 3.0*qmax +(0.5*hsm*bd[0] + 2.*nu1/hsm)*(0.5*hsm*bd[0] + 2.*nu1/hsm);
      gam2 = 3.0*qmax +(0.5*hsm*bd[0] + 2.*nu2/hsm)*(0.5*hsm*bd[0] + 2.*nu2/hsm);
#ifdef INERTIALESS
      gam1 = (2.*nu1/hsm)*(2.*nu1/hsm);
      gam2 = (2.*nu2/hsm)*(2.*nu2/hsm);
#endif
      dtnorm = 2.*vslp/hsm +bd[0] +1.*strss/(gbl->rho*sqrt(qmax +gam1) +srf->gbl->rho2*sqrt(qmax +gam2));
      /* SET UP DISSIPATIVE COEFFICIENT */
      /* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
      /* RESIDUAL HAS DX/2 WEIGHTING */
      /* |a| dx/2 dv/dx [a du/dx] dx/2 dpsi */
      /* |a| dx/2 2/dx dv/dpsi [a du/dx dx/2] dpsi */
      /* |a| dv/dpsi [a du/dx dx/2] dpsi */
      srf->gbl->meshc[indx] = adis/(h*dtnorm*0.5);
      /* srf->gbl->meshc[indx] = adis/(h*vslp/hsm); */
         
#ifdef AXISYMMETRIC
      dtnorm *= 0.5*(vrtx[v0][0] +vrtx[v1][0]);
#endif

      nrm[0] *= 0.5;
      nrm[1] *= 0.5;
            
      srf->gbl->vdt[indx][0][0] += -dttang*nrm[1]*b->vdiag1d;
      srf->gbl->vdt[indx][0][1] +=  dttang*nrm[0]*b->vdiag1d;
      srf->gbl->vdt[indx][1][0] +=  dtnorm*nrm[0]*b->vdiag1d;
      srf->gbl->vdt[indx][1][1] +=  dtnorm*nrm[1]*b->vdiag1d;
      srf->gbl->vdt[indx+1][0][0] = -dttang*nrm[1]*b->vdiag1d;
      srf->gbl->vdt[indx+1][0][1] =  dttang*nrm[0]*b->vdiag1d;
      srf->gbl->vdt[indx+1][1][0] =  dtnorm*nrm[0]*b->vdiag1d;
      srf->gbl->vdt[indx+1][1][1] =  dtnorm*nrm[1]*b->vdiag1d;
               
      if (b->sm) {
         srf->gbl->sdt[indx][0][0] = -dttang*nrm[1];
         srf->gbl->sdt[indx][0][1] =  dttang*nrm[0];
         srf->gbl->sdt[indx][1][0] =  dtnorm*nrm[0];
         srf->gbl->sdt[indx][1][1] =  dtnorm*nrm[1];
      }
   }
   
   // printf("# Surface Max Time Step: %e\n",dtmax);
   
   /* SEND COMMUNICATION INFO FOR ENDPOINTS */
   end = sbdry[bnum].num;
   v0 = svrtx[sbdry[bnum].el[0]][0];
   v1 = svrtx[sbdry[bnum].el[end-1]][1];
   if (v0 == v1) {
      /* SURFACE IS A LOOP ON SINGLE BLOCK */
      for(m=0;m<ND;++m) {
         for(n=0;n<ND;++n) {
            srf->gbl->vdt[0][m][n] = 0.5*(srf->gbl->vdt[0][m][n] +srf->gbl->vdt[end][m][n]);
            srf->gbl->vdt[end][m][n] = srf->gbl->vdt[0][m][n];
         }
      }
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      
      if (vbdry[i].type&(COMX_MASK+COMY_MASK)) {
         vbnum = vbdry[i].adjbnum;
         tgt = vbdry[i].adjmesh;
         count = 0;
         for(m=0;m<ND;++m)
            for(n=0;n<ND;++n) 
               tgt->vbuff[vbnum][count++] = srf->gbl->vdt[0][m][n];
      }
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      
      if (vbdry[i].type&(COMX_MASK+COMY_MASK)) {
         vbnum = vbdry[i].adjbnum;
         tgt = vbdry[i].adjmesh; 
         count = 0;
         for(m=0;m<ND;++m)
            for(n=0;n<ND;++n) 
               tgt->vbuff[vbnum][count++] = srf->gbl->vdt[end][m][n];
      }
   }
   
   return;
}

void hp_mgrid::surfdt2(int bnum) {
   int i,m,n,indx,v0,v1,end,count;
   FLT jcbi,temp;
   class surface *srf;
   
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;

   /* RECEIVE COMMUNICATION PACKETS */
   /* RECEIVE MESSAGES AND ZERO VERTEX MODES */
   end = sbdry[bnum].num;
   v0 = svrtx[sbdry[bnum].el[0]][0];
   v1 = svrtx[sbdry[bnum].el[end-1]][1];
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      
      if (vbdry[i].type&(COMX_MASK+COMY_MASK)) {
         count = 0;
         for(m=0;m<ND;++m)
            for(n=0;n<ND;++n) 
               srf->gbl->vdt[0][m][n] = 0.5*(srf->gbl->vdt[0][m][n] +vbuff[i][count++]);
      }
   }

   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      
      if (vbdry[i].type&(COMX_MASK+COMY_MASK)) {
         count = 0;
         for(m=0;m<ND;++m)
            for(n=0;n<ND;++n) 
               srf->gbl->vdt[end][m][n] = 0.5*(srf->gbl->vdt[end][m][n] +vbuff[i][count++]);
      }
   } 
   
   for(indx=0;indx<sbdry[bnum].num+1;++indx) {   
      /* INVERT VERTEX MATRIX */
      jcbi = 1.0/(srf->gbl->vdt[indx][0][0]*srf->gbl->vdt[indx][1][1] 
                 -srf->gbl->vdt[indx][0][1]*srf->gbl->vdt[indx][1][0]);
                 
      temp = srf->gbl->vdt[indx][0][0]*jcbi*surface::cfl[log2p][1];
      srf->gbl->vdt[indx][0][0] = srf->gbl->vdt[indx][1][1]*jcbi*surface::cfl[log2p][0];
      srf->gbl->vdt[indx][1][1] = temp;
      srf->gbl->vdt[indx][0][1] *= -jcbi*surface::cfl[log2p][1];
      srf->gbl->vdt[indx][1][0] *= -jcbi*surface::cfl[log2p][0];
   }

   /* INVERT SIDE MATRIX */   
   if (b->sm > 0) {
      for(indx=0;indx<sbdry[bnum].num;++indx) {
         /* INVERT SIDE MVDT MATRIX */
         jcbi = 1.0/(srf->gbl->sdt[indx][0][0]*srf->gbl->sdt[indx][1][1] 
                    -srf->gbl->sdt[indx][0][1]*srf->gbl->sdt[indx][1][0]);

         temp = srf->gbl->sdt[indx][0][0]*jcbi*surface::cfl[log2p][1];
         srf->gbl->sdt[indx][0][0] = srf->gbl->sdt[indx][1][1]*jcbi*surface::cfl[log2p][0];
         srf->gbl->sdt[indx][1][1] = temp;
         srf->gbl->sdt[indx][0][1] *= -jcbi*surface::cfl[log2p][1];
         srf->gbl->sdt[indx][1][0] *= -jcbi*surface::cfl[log2p][0];
      }
   }
   
   return;
}

void hp_mgrid::surfnstage1(int bnum) {
   int i,n;
   class surface *srf;
   
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;

   for(i=0;i<sbdry[bnum].num+1;++i)
      for(n=0;n<ND;++n)
         srf->gbl->vug0[i][n] = srf->vug[i][n];
         
   if (b->sm > 0) {
      for(i=0;i<sbdry[bnum].num*sm0;++i)
         for(n=0;n<ND;++n)
            srf->gbl->sug0[i][n] = srf->sug[i][n];
   }
   
   return;
}

void hp_mgrid::surfnstage2(int bnum, int stage) {
   int i,m,n,indx,indx1,end,v0,v1;
   class surface *srf;
   
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;
   
   for(i=0;i<sbdry[bnum].num+1;++i) {
      for(n=0;n<ND;++n) {
         srf->vug[i][n] = srf->gbl->vug0[i][n] -alpha[stage]*srf->gbl->vres[i][n];
      }
   }
   
   /* FIX POINTS THAT SLIDE ON CURVE HERE? */
   end = sbdry[bnum].num;
   v0 = svrtx[sbdry[bnum].el[0]][0];
   v1 = svrtx[sbdry[bnum].el[end-1]][1];
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      if (vbdry[i].type&(CURV_MASK)) {
         mvpttobdry(vbdry[i].type, srf->vug[0][0], srf->vug[0][1]);
      }
   }  
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      if (vbdry[i].type&(CURV_MASK)) {
         mvpttobdry(vbdry[i].type, srf->vug[end][0], srf->vug[end][1]);
      }
   }       
   

   if (b->sm > 0) {
      indx = 0;
      indx1 = 0;
      for(i=0;i<sbdry[bnum].num;++i) {
         for(m=0;m<b->sm;++m) {
            for(n=0;n<ND;++n)
               srf->sug[indx1][n] = srf->gbl->sug0[indx1][n] -alpha[stage]*srf->gbl->sres[indx][n];
            ++indx;
            ++indx1;
         }
         indx1 += sm0 -b->sm;
      }
   }
   
   return;
}

void hp_mgrid::surfugtovrt1() {
   int i,m,n,sind,indx,v0,bnum,vbnum,count;
   class surface *srf;
   class mesh *tgt;
   
   for (bnum=0;bnum<nsbd;++bnum) {
      if (!(sbdry[bnum].type&(FSRF_MASK +IFCE_MASK))) continue;
   
      srf = static_cast<class surface *>(sbdry[bnum].misc);
      if (srf == NULL) continue;
      
      indx = 0;
      for(i=0;i<sbdry[bnum].num;++i) {
         sind = sbdry[bnum].el[i];
         v0 = svrtx[sind][0];
         for(n=0;n<ND;++n)
            vrtx[v0][n] = srf->vug[i][n];
         
         for(m=0;m<b->sm;++m) {
            for(n=0;n<ND;++n)
               binfo[bnum][indx].curv[n] = srf->sug[indx][n];
            ++indx;
         }
         indx += sm0 -b->sm;
      }
      v0 = svrtx[sind][1];
      for(n=0;n<ND;++n)
         vrtx[v0][n] = srf->vug[sbdry[bnum].num][n];
            
      /* SEND MESSAGES TO MATCHING INTERFACE */      
      if (sbdry[bnum].type&IFCE_MASK) {
         tgt = sbdry[bnum].adjmesh;
         vbnum = sbdry[bnum].adjbnum;
         count = 0;
         for(i=0;i<sbdry[bnum].num;++i) {
            for(n=0;n<ND;++n)
               tgt->sbuff[vbnum][count++] = srf->vug[i][n];
            for(m=0;m<b->sm;++m)
               for(n=0;n<ND;++n)
                  tgt->sbuff[vbnum][count++] = srf->sug[i*sm0 +m][n];
         }
         for(n=0;n<ND;++n)
            tgt->sbuff[vbnum][count++] = srf->vug[sbdry[bnum].num][n];
      }
   }
   
   return;
}

void hp_mgrid::surfugtovrt2() {
   int j,m,n,sind,indx,msgn,v0,count,bnum;
   class surface *srf;
   
   for (bnum=0;bnum<nsbd;++bnum) {
      if (!(sbdry[bnum].type&IFCE_MASK)) continue;
      
      srf = static_cast<class surface *>(sbdry[bnum].misc);
      if (srf != NULL) return;
   
      /* RECEIVE MESSAGES FROM MATCHING INTERFACE */      
      count = 0;
      for(j=sbdry[bnum].num-1;j>=0;--j) {
         sind = sbdry[bnum].el[j];
         v0 = svrtx[sind][1];
         for(n=0;n<ND;++n)
            vrtx[v0][n] = sbuff[bnum][count++];
         
         indx = j*sm0;
         msgn = 1;
         for(m=0;m<b->sm;++m) {
            for(n=0;n<ND;++n)
               binfo[bnum][indx+m].curv[n] = msgn*sbuff[bnum][count++];
            msgn *= -1;
         }      
      }
      v0 = svrtx[sind][0];
      for(n=0;n<ND;++n)
         vrtx[v0][n] = sbuff[bnum][count++];
   }
   
   return;
}

void hp_mgrid::surfvrttoug() {
   int i,m,n,indx,sind,v0,bnum;
   class surface *srf;
   
   for (bnum=0;bnum<nsbd;++bnum) {
      if (!(sbdry[bnum].type&(FSRF_MASK +IFCE_MASK))) continue;
      
      srf = static_cast<class surface *>(sbdry[bnum].misc);
      if (srf == NULL) continue;
            
      indx = 0;
      for(i=0;i<sbdry[bnum].num;++i) {
         sind = sbdry[bnum].el[i];
         v0 = svrtx[sind][0];
         for(n=0;n<ND;++n)
            srf->vug[i][n] = vrtx[v0][n];
         
         for(m=0;m<b->sm;++m) {
            for(n=0;n<ND;++n)
               srf->sug[indx][n] = binfo[bnum][indx].curv[n];
            ++indx;
         }
      }
      v0 = svrtx[sind][1];
      for(n=0;n<ND;++n)
         srf->vug[sbdry[bnum].num][n] = vrtx[v0][n];
   }
         
   return;
}

void hp_mgrid::surfgetfres(int bnum) {
   int i,j,n,indx,indx1,tind,v0,snum,sind;
   class surface *srf, *finesrf;
   class hp_mgrid *fmesh;
   
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;
   
   /* TRANSFER INTERFACE RESIDUALS */

   if(p0 > 1) {
      /* TRANSFER IS ON FINEST MESH */
      for(i=0;i<sbdry[bnum].num+1;++i)
         for(n=0;n<ND;++n)
            srf->gbl->vres0[i][n] = srf->gbl->vres[i][n];

      if (b->sm > 0) {
         indx = 0;
         indx1 = 0;
         for(i=0;i<sbdry[bnum].num;++i) {
            for (j=0;j<b->sm;++j) {
               for(n=0;n<ND;++n)
                     srf->gbl->sres0[indx][n] = srf->gbl->sres[indx1][n];
               ++indx;
               ++indx1;
            }
            indx1 += b->p;
         }
      }
      return;
   }
   
   fmesh = static_cast<class hp_mgrid *>(fmpt);
   
   /* TRANSFER IS BETWEEN DIFFERENT MESHES */
   for(i=0;i<sbdry[bnum].num +1;++i)
      for(n=0;n<ND;++n)
         srf->gbl->vres0[i][n] = 0.0;
         
   /* CALCULATE COARSE RESIDUALS */
   /* DO ENDPOINTS FIRST */
   for(n=0;n<ND;++n)
      srf->gbl->vres0[0][n] = srf->gbl->vres[0][n];

   for(n=0;n<ND;++n)
      srf->gbl->vres0[sbdry[bnum].num][n] = srf->gbl->vres[fmesh->sbdry[bnum].num][n];
      
   for(i=1;i<fmesh->sbdry[bnum].num;++i) {
      sind = fmesh->sbdry[bnum].el[i];
      v0 = fmesh->svrtx[sind][0];
      tind = fmesh->coarse[v0].tri;
      for(snum=0;snum<3;++snum) 
         if ((-ttri[tind][snum]>>16) -1  == bnum) break;
      assert(snum != 3);
      indx = -ttri[tind][snum]&0xFFFF;
      for(n=0;n<ND;++n) {
         srf->gbl->vres0[indx][n] += fmesh->coarse[v0].wt[(snum+1)%3]*srf->gbl->vres[i][n];
         srf->gbl->vres0[indx+1][n] += fmesh->coarse[v0].wt[(snum+2)%3]*srf->gbl->vres[i][n];
      }
   }

   /* CALCULATE VALUES OF SOLUTION ON COARSE MESH */
   finesrf = static_cast<class surface *>(fmesh->sbdry[bnum].misc);
   
   for(i=0;i<sbdry[bnum].num +1;++i)
      for(n=0;n<ND;++n)
         srf->vug[i][n] = 0.0;
         
   /* DO ENDPOINTS FIRST */
   for(n=0;n<ND;++n)
      srf->vug[0][n] = finesrf->vug[0][n];

   for(n=0;n<ND;++n)
      srf->vug[sbdry[bnum].num][n] = finesrf->vug[fmesh->sbdry[bnum].num][n];

   /* NOW CALCULATE INTERIOR POINTS */
   for(i=1;i<sbdry[bnum].num;++i) {
      sind = sbdry[bnum].el[i];
      v0 = svrtx[sind][0];
      tind = fine[v0].tri;
      for(snum=0;snum<3;++snum) 
         if ((-fmesh->ttri[tind][snum]>>16) -1  == bnum) break;
      assert(snum != 3);
      indx = -fmesh->ttri[tind][snum]&0xFFFF;
          
      for(n=0;n<ND;++n) {
         srf->vug[i][n] += fine[v0].wt[(snum+1)%3]*finesrf->vug[indx][n];
         srf->vug[i][n] += fine[v0].wt[(snum+2)%3]*finesrf->vug[indx+1][n];
      }
   }
    
   /* STORE UG_FRST */
   for(i=0;i<sbdry[bnum].num+1;++i)
      for(n=0;n<ND;++n)
         srf->vug_frst[i][n] = srf->vug[i][n];
   
   return;
}

void hp_mgrid::surfgetcchng(int bnum) {
   int i,n,tind,sind,v0,snum,indx;
   class surface *srf, *coarsesrf;
   class hp_mgrid *cmesh;

    if(b->p > 1) {
      return;
   } 
    
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;
   
   cmesh = static_cast<class hp_mgrid *>(cmpt);
   coarsesrf = static_cast<class surface *>(cmesh->sbdry[bnum].misc);
   
   /* DETERMINE CORRECTIONS ON COARSE MESH   */   
   for(i=0;i<cmesh->sbdry[bnum].num+1;++i)
      for(n=0;n<ND;++n)
         coarsesrf->vug_frst[i][n] -= coarsesrf->vug[i][n];
   
   /* LOOP THROUGH FINE VERTICES   */
   /* TO DETERMINE CHANGE IN SOLUTION */   
   /* DO ENDPOINTS FIRST */
   for(n=0;n<ND;++n)
      srf->gbl->vres[0][n] = -coarsesrf->vug_frst[0][n];

   for(n=0;n<ND;++n)
      srf->gbl->vres[sbdry[bnum].num][n] = -coarsesrf->vug_frst[cmesh->sbdry[bnum].num][n];
      
   for(i=1;i<sbdry[bnum].num;++i) {
      
      for(n=0;n<ND;++n)
         srf->gbl->vres[i][n] = 0.0;
         
      sind = sbdry[bnum].el[i];
      v0 = svrtx[sind][0];
      tind = coarse[v0].tri;
      for(snum=0;snum<3;++snum) 
         if ((-cmesh->ttri[tind][snum]>>16) -1  == bnum) break;
      assert(snum != 3);
      indx = -cmesh->ttri[tind][snum]&0xFFFF;
         
      for(n=0;n<ND;++n) {
         srf->gbl->vres[i][n] -= coarse[v0].wt[(snum+1)%3]*coarsesrf->vug_frst[indx][n];
         srf->gbl->vres[i][n] -= coarse[v0].wt[(snum+2)%3]*coarsesrf->vug_frst[indx+1][n];
      }
   }
   
   for(i=0;i<sbdry[bnum].num+1;++i)
      for(n=0;n<ND;++n) 
         srf->vug[i][n] += srf->gbl->vres[i][n];

   return;
}
    
void hp_mgrid::surfmaxres() {
   int i,n,bnum;
   class surface *srf;
   FLT mxr[ND];
   
   for (bnum=0;bnum<nsbd;++bnum) {
      if (!(sbdry[bnum].type&(FSRF_MASK +IFCE_MASK))) continue;
   
      srf = static_cast<class surface *>(sbdry[bnum].misc);
      if (srf == NULL) continue;
      

      for(n=0;n<ND;++n)
         mxr[n] = 0.0;
      
      for(i=0;i<sbdry[bnum].num+1;++i)
         for(n=0;n<ND;++n)
            mxr[n] = MAX(fabs(srf->gbl->vres[i][n]),mxr[n]);

      for(n=0;n<ND;++n)
         printf("%.3e  ",mxr[n]);
   }
   
   return;
}


/* CALCULATE AREA/CIRCUMFERENCE/YBAR */
void hp_mgrid::integrated_averages(FLT a[]) {
   int i,j,n,tind;

   /* a[0] = area */
   /* a[1] = xbar */
   /* a[2] = ybar */
   /* a[3] = ubar */
   /* a[4] = vbar */
   for (i = 0; i <5; ++i)
      a[i] = 0.0;
   
   for(tind=0;tind<ntri;++tind) {
      if (tinfo[tind] > -1) {
         crdtocht(tind);
         for(n=0;n<ND;++n)
            b->proj_bdry(cht[n], crd[n][0], dcrd[n][0][0], dcrd[n][1][0],MXGP);
      }
      else {
         for(n=0;n<ND;++n)
            b->proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n][0],MXGP);

         for(i=0;i<b->gpx;++i) {
            for(j=0;j<b->gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd[n][0][i][j] = 0.5*(vrtx[tvrtx[tind][2]][n] -vrtx[tvrtx[tind][1]][n]);
                  dcrd[n][1][i][j] = 0.5*(vrtx[tvrtx[tind][0]][n] -vrtx[tvrtx[tind][1]][n]);
               }
            }
         }
      }
      
      ugtouht(tind);
      b->proj(uht[0],u[0][0],MXGP);
      b->proj(uht[1],u[1][0],MXGP);
      for(i=0;i<b->gpx;++i) {
         for(j=0;j<b->gpn;++j) {
            cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
            a[0] += RAD(i,j)*b->wtx(i)*b->wtn(j)*cjcb[i][j];
            a[1] += crd[0][i][j]*RAD(i,j)*b->wtx(i)*b->wtn(j)*cjcb[i][j];
            a[2] += crd[1][i][j]*RAD(i,j)*b->wtx(i)*b->wtn(j)*cjcb[i][j];
            a[3] += u[0][i][j]*RAD(i,j)*b->wtx(i)*b->wtn(j)*cjcb[i][j];
            a[4] += u[1][i][j]*RAD(i,j)*b->wtx(i)*b->wtn(j)*cjcb[i][j];
         }
      }
   }
   for (i=1;i<5;++i)
      a[i] /= a[0];
      
//   /* CALCULATE CIRCUMFERENCE */
//   FLT circ = 0.0;
//   for(indx=0;indx<sbdry[bnum].num;++indx) {
//      sind = sbdry[bnum].el[indx];
//      
//      /* UG STORES THE CURRENT BOUNDARY POSITION */
//      for(n=0;n<ND;++n) {
//         uht[n][0] = srf->vug[indx][n];
//         uht[n][1] = srf->vug[indx+1][n];
//      }
//      if (b->sm > 0) {
//         indx1 = indx*sm0;
//         for(m=0;m<b->sm;++m)
//            for(n=0;n<ND;++n)
//               uht[n][m+2] = srf->sug[indx1 +m][n];
//      }
//      
//      for(n=0;n<ND;++n)
//         b->proj1d(uht[n],crd[n][0],dcrd[n][0][0]);
//      
//      for(i=0;i<b->gpx;++i) {
//         norm[0] =  dcrd[1][0][0][i];
//         norm[1] = -dcrd[0][0][0][i];
//         jcb = RAD1D(i)*sqrt(norm[0]*norm[0] +norm[1]*norm[1]);
//         circ += jcb*b->wtx[i];
//      }
//   }
   
   return;
}

