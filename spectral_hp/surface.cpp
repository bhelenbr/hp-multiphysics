#include"hp_mgrid.h"
#include<assert.h>
#include<utilities.h>
#include<myblas.h>

FLT surface::fadd[ND];
FLT surface::cfl[MXLG2P][ND];

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
   vect_alloc(store->normc,maxside,FLT);
   vect_alloc(store->meshc,maxside,FLT);
   
   return;
}

void hp_mgrid::setksprg1d() {
   int i,j,sind;
   class surface *srf;

   for(i=0;i<nsbd;++i) {
      if(sbdry[i].type&(CURV_MASK+IFCE_MASK)) {
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
      if(sbdry[bnum].type&(CURV_MASK+IFCE_MASK)) {
         srf = static_cast<class surface *>(sbdry[bnum].misc);
         if (srf != NULL) {
            /* ZERO TANGENTIAL MESH MOVEMENT SOURCE */   
            for(i=0;i<sbdry[bnum].num+1;++i)
               srf->vdres[log2p][i][0] = 0.0;

            for(i=0;i<sbdry[bnum].num*b.sm;++i) 
               srf->sdres[log2p][i][0] = 0.0;

            surfrsdl(bnum, 0);
            
            /* TANGENTIAL MESH MOVEMENT SOURCE */   
            for(i=0;i<sbdry[bnum].num+1;++i) {
               srf->vdres[log2p][i][0] = -srf->gbl->vres[i][0];
            }

            for(i=0;i<sbdry[bnum].num*b.sm;++i) {
               srf->sdres[log2p][i][0] = -0.0*srf->gbl->sres[i][0];  // TEMPORARY FOR UNIFORM SPACING OF HIGHER MODES ALONG SIDE
            }
         }
      }
   }
   
   
   return;
}



void hp_mgrid::surfrsdl(int bnum, int mgrid) {
   int i,m,n,sind,indx,indx1,count,v0;
   FLT norm[ND], rp[ND], jcb, tau, tabs;
   FLT dnormdt = 0.0, hsm;
   FLT sigor, drhor;
   class surface *srf;

   /* DETRMINE DX CORRECTION TO CONSERVE AREA */
   /* IMPORTANT FOR STEADY SOLUTIONS */
   /* SINCE THERE ARE MULTIPLE STEADY-STATES */
   if (!(sbdry[bnum].type&(IFCE_MASK +FSRF_MASK))) {
      printf("error shouldn't be in surfrsdl\n");
      exit(1);
   }

   srf = static_cast<class surface *>(sbdry[bnum].misc);
   
/* POINTER TO STUFF NEEDED FOR SURFACES IS STORED IN MISCELLANEOUS */      
   if (srf == NULL) return;
   
   count = 0;
   sigor = srf->gbl->sigma/(gbl->rho +srf->gbl->rho2);
   drhor = (gbl->rho -srf->gbl->rho2)/(gbl->rho +srf->gbl->rho2);
      
   /* CONSERVE AREA FOR STEADY CLOSED BDRY PROBLEMS */   
//   if (bd[0] == 0.0) dnormdt = srf->gbl->lamv*cnsrvarea(bnum);
//   else dnormdt = 0.0;

   /**************************************************/
   /* DETERMINE MESH RESIDUALS & SURFACE TENSION     */
   /**************************************************/
   for(n=0;n<ND;++n)
      srf->gbl->vres[0][n] = 0.0;
 
   for(n=0;n<ND;++n) 
      binfo[bnum][0].flx[n] = 0.0;

   for(indx=0;indx<sbdry[bnum].num;++indx) {
      sind = sbdry[bnum].el[indx];

      /* THIS IS FOR LINEARIZED SURFACES */      
//      crdtocht1d(sind);
//      for(n=0;n<ND;++n)
//         b.proj1d(cht[n],crd[n][0],dcrd[n][0][0]);

      /* UG STORES THE CURRENT BOUNDARY POSITION */
      for(n=0;n<ND;++n) {
         uht[n][0] = srf->vug[indx][n];
         uht[n][1] = srf->vug[indx+1][n];
      }
      if (b.sm > 0) {
         indx1 = indx*sm0;
         for(m=0;m<b.sm;++m)
            for(n=0;n<ND;++n)
               uht[n][m+2] = srf->sug[indx1 +m][n];
      }

      for(n=0;n<ND;++n)
         b.proj1d(uht[n],crd[n][0],dcrd[n][0][0]);
         
      crdtocht1d(sind,dvrtdt,gbl->dbinfodt);
      for(n=0;n<ND;++n)
         b.proj1d(cht[n],mvel[n][0]);

      ugtouht1d(sind);
      for(n=0;n<ND;++n)
         b.proj1d(uht[n],u[n][0]);   
         
      for(i=0;i<b.gpx;++i) {
         norm[0] =  dcrd[1][0][0][i];
         norm[1] = -dcrd[0][0][0][i];
         jcb = sqrt(norm[0]*norm[0] +norm[1]*norm[1]);
         
         /* RELATIVE VELOCITY STORED IN MVEL[N][0]*/
         for(n=0;n<ND;++n)
            mvel[n][0][i] = u[n][0][i] -(bd[0]*crd[n][0][i] +mvel[n][0][i] +dnormdt*norm[n]/jcb); 

         hsm = jcb/(.25*(b.p+1)*(b.p+1));
         tau = (mvel[0][0][i]*dcrd[0][0][0][i] +mvel[1][0][i]*dcrd[1][0][0][i])/jcb;
         tabs = fabs(tau) + 10.*EPSILON;
         tau = tau/(jcb*(tabs/hsm +bd[0] +(sigor/(hsm*hsm) +fabs(drhor*g*norm[1]/jcb))/tabs));

         /* TANGENTIAL SPACING & NORMAL FLUX */            
         res[0][0][i] = -srf->ksprg[indx]*jcb;
         res[1][0][i] = -RAD1D(i)*(mvel[0][0][i]*norm[0] +mvel[1][0][i]*norm[1]);
         res[1][1][i] = -res[1][0][i]*tau;
         
         /* SURFACE TENSION SOURCE TERM */ 
         u[0][0][i] = +RAD1D(i)*srf->gbl->sigma*norm[1]/jcb;
         u[0][1][i] = +RAD1D(i)*(gbl->rho -srf->gbl->rho2)*g*crd[1][0][i]*norm[0];
#ifdef AXISYMMETRIC
         u[0][1][i] += srf->gbl->sigma*jcb;
#endif
         u[1][0][i] = -RAD1D(i)*srf->gbl->sigma*norm[0]/jcb;
         u[1][1][i] = +RAD1D(i)*(gbl->rho -srf->gbl->rho2)*g*crd[1][0][i]*norm[1];            
      }
      
      for(m=0;m<b.sm+2;++m)
         lf[0][m] = 0.0;

      /* INTEGRATE & STORE SURFACE RESIDUALS */               
      b.intgrtx1d(res[0][0],lf[0]);
      b.intgrt1d(res[1][0],lf[1]);
      b.intgrtx1d(res[1][1],lf[1]);

      /* STORE IN RES */
      for(n=0;n<ND;++n) {
         srf->gbl->vres[indx][n] += lf[n][0];
         srf->gbl->vres[indx+1][n] = lf[n][1];
         for(m=0;m<b.sm;++m)
            srf->gbl->sres[b.sm*indx +m][n] = lf[n][m+2];
      }
      
      /* INTEGRATE & STORE SURFACE TENSION SOURCE TERM */
      b.intgrt1d(u[0][1],lf[0]);
      b.intgrtx1d(u[0][0],lf[0]);
      b.intgrt1d(u[1][1],lf[1]);
      b.intgrtx1d(u[1][0],lf[1]);
      
      /* STORE IN BINFO.FLUX */
      for(n=0;n<ND;++n) 
         binfo[bnum][count].flx[n] += lf[n][0];
      ++count;
      
      for(m=0;m<b.sm;++m) {
         for(n=0;n<ND;++n)
            binfo[bnum][count].flx[n] = lf[n][m+2];
         ++count;
      }
      
      for(n=0;n<ND;++n)
         binfo[bnum][count].flx[n] = lf[n][1];
   }
   
   /* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
   v0 = svrtx[sbdry[bnum].el[0]][0];
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0 || !(vbdry[i].type&OUTF_MASK)) continue;

      /* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
      crdtocht1d(sbdry[bnum].el[0]);
      b.ptprobe1d(2,cht,rp,norm,-1.0);
      jcb = sqrt(norm[0]*norm[0] +norm[1]*norm[1]);
#ifdef AXISYMMETRIC
      binfo[bnum][0].flx[0] -= -rp[0]*srf->gbl->sigma*norm[0]/jcb;
      binfo[bnum][0].flx[1] -= -rp[0]*srf->gbl->sigma*norm[1]/jcb;
#else
      binfo[bnum][0].flx[0] -= -srf->gbl->sigma*norm[0]/jcb;
      binfo[bnum][0].flx[1] -= -srf->gbl->sigma*norm[1]/jcb;
#endif
   }
   
   /* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
   v0 = svrtx[sbdry[bnum].el[sbdry[bnum].num -1]][1];
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0 || !(vbdry[i].type&OUTF_MASK)) continue;

      /* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
      crdtocht1d(sbdry[bnum].el[sbdry[bnum].num -1]);
      b.ptprobe1d(2,cht,rp,norm,1.0);
      jcb = sqrt(norm[0]*norm[0] +norm[1]*norm[1]);
#ifdef AXISYMMETRIC
      binfo[bnum][count].flx[0] += -rp[0]*srf->gbl->sigma*norm[0]/jcb;
      binfo[bnum][count].flx[1] += -rp[0]*srf->gbl->sigma*norm[1]/jcb;
#else
      binfo[bnum][count].flx[0] += -srf->gbl->sigma*norm[0]/jcb;
      binfo[bnum][count].flx[1] += -srf->gbl->sigma*norm[1]/jcb;
#endif
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
         for(i=0;i<sbdry[bnum].num*b.sm;++i) {
            for(n=0;n<ND;++n)
               srf->sdres[log2p][i][n] = surface::fadd[n]*srf->gbl->sres0[i][n] -srf->gbl->sres[i][n];
         }
      }
      for(i=0;i<sbdry[bnum].num+1;++i) {
         for(n=0;n<ND;++n)
            srf->gbl->vres[i][n] += srf->vdres[log2p][i][n];
      }
      for(i=0;i<sbdry[bnum].num*b.sm;++i) {
         for(n=0;n<ND;++n)      
            srf->gbl->sres[i][n] += srf->sdres[log2p][i][n];
      }
   }
   else {
      /* ADD TANGENTIAL MESH MOVEMENT SOURCE */   
      for(i=0;i<sbdry[bnum].num+1;++i)
         srf->gbl->vres[i][0] += srf->vdres[log2p][i][0];

      for(i=0;i<sbdry[bnum].num*b.sm;++i) 
         srf->gbl->sres[i][0] += srf->sdres[log2p][i][0];
   }

   /*
      for(i=0;i<sbdry[bnum].num+1;++i)
            printf("vres: %f %f\n",srf->gbl->vres[i][0],srf->gbl->vres[i][1]);
            
      for(i=0;i<sbdry[bnum].num*b.sm;++i)
         printf("sres: %f %f\n",srf->gbl->sres[i][0],srf->gbl->sres[i][1]);
   */
   
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
   if (b.sm > 0) {
      indx1 = 0;
      for(indx = 0; indx<sbdry[bnum].num; ++indx) {
         /* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */         
         for (m=0; m <b.sm; ++m) {
            for(n=0;n<ND;++n)
               srf->gbl->vres[indx][n] -= b.sfmv1d[0][m]*srf->gbl->sres[indx1][n];
            for(n=0;n<ND;++n)
               srf->gbl->vres[indx+1][n] -= b.sfmv1d[1][m]*srf->gbl->sres[indx1][n];
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
            
      if (vbdry[i].type&(PRDX_MASK +SYMM_MASK +PRDY_MASK)) {
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
            
      if (vbdry[i].type&(PRDX_MASK +PRDY_MASK +SYMM_MASK)) {
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
      if (vbdry[i].type&(PRDX_MASK +SYMM_MASK))
         srf->gbl->vres[0][0] = 0.0;
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
   if (b.sm > 0) {
      indx1 = 0;
      for(indx = 0; indx<sbdry[bnum].num; ++indx) {
         
         /* INVERT SIDE MODES */
         DPBTRSNU(b.sdiag1d,b.sm,b.sbwth,&(srf->gbl->sres[indx1][0]),ND);
         for(m=0;m<b.sm;++m) {
            temp                      = srf->gbl->sres[indx1][0]*srf->gbl->sdt[indx][0][0] +srf->gbl->sres[indx1][1]*srf->gbl->sdt[indx][0][1];
            srf->gbl->sres[indx1][1] = srf->gbl->sres[indx1][0]*srf->gbl->sdt[indx][1][0] +srf->gbl->sres[indx1][1]*srf->gbl->sdt[indx][1][1];       
            srf->gbl->sres[indx1][0] = temp;
            ++indx1;
         }
         
         indx1 -= b.sm;
         for(m=0;m<b.sm;++m) {
            for(n=0;n<ND;++n) {
               srf->gbl->sres[indx1][n] -= b.vfms1d[0][m]*srf->gbl->vres[indx][n];
               srf->gbl->sres[indx1][n] -= b.vfms1d[1][m]*srf->gbl->vres[indx+1][n];
            }
            ++indx1;
         }
      }
   }
   
//   for(i=0;i<sbdry[bnum].num*b.sm;++i) 
//      printf("s: %f %f\n",srf->gbl->sres[i][0],srf->gbl->sres[i][1]);
   
         
   return;
}

void hp_mgrid::surfdt1(int bnum) {
   int i,m,n,indx,end,vbnum,count,sind,v0,v1;
   FLT nrm[ND], h, hsm;
   FLT dttang, dtnorm;
   FLT vslp, strss, cnvct, dtfli;
   class surface *srf;
   class mesh *tgt;
   FLT drho, srho, smu;
   
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;
      
   drho = gbl->rho -srf->gbl->rho2;
   srho = gbl->rho +srf->gbl->rho2;
   smu = gbl->mu +srf->gbl->mu2;

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

      nrm[0] =  0.5*(vrtx[v1][1] -vrtx[v0][1]);
      nrm[1] = -0.5*(vrtx[v1][0] -vrtx[v0][0]);
      h = 2.0*sqrt(nrm[0]*nrm[0] +nrm[1]*nrm[1]);
 
      vslp = fabs(-((ug.v[v0][0]-(bd[0]*vrtx[v0][0] +dvrtdt[v0][0]))+(ug.v[v1][0]-(bd[0]*vrtx[v1][0] +dvrtdt[v1][0])))*nrm[1]/h
               +((ug.v[v0][1]-(bd[0]*vrtx[v0][1] +dvrtdt[v0][1]))+(ug.v[v1][1]-(bd[0]*vrtx[v1][1] +dvrtdt[v1][1])))*nrm[0]/h);
      
      hsm = h/(.25*(b.p+1)*(b.p+1));

      strss = srf->gbl->sigma/(hsm*hsm) +2.*fabs(drho*g*nrm[1]/h);
      cnvct = bd[0] + vslp/hsm;
      dtfli = smu/(srho*hsm*hsm) +vslp/hsm +bd[0];
      dttang  = 2.*srf->ksprg[indx]*(.25*(b.p+1)*(b.p+1))/hsm;
      dtnorm  = srho*hsm*cnvct*dtfli +strss;
      
      srf->gbl->normc[indx] = srho*hsm*cnvct*dtfli/dtnorm;
      srf->gbl->meshc[indx] = srho*hsm*dtfli;
      dtnorm = dtnorm/srf->gbl->meshc[indx];
#ifdef AXISYMMETRIC
      dtnorm *= 0.5*(vrtx[v0][0] +vrtx[v1][0]);
#endif

      srf->gbl->vdt[indx][0][0] += -dttang*nrm[1]*b.vdiag1d;
      srf->gbl->vdt[indx][0][1] +=  dttang*nrm[0]*b.vdiag1d;
      srf->gbl->vdt[indx][1][0] +=  dtnorm*nrm[0]*b.vdiag1d;
      srf->gbl->vdt[indx][1][1] +=  dtnorm*nrm[1]*b.vdiag1d;
      srf->gbl->vdt[indx+1][0][0] = -dttang*nrm[1]*b.vdiag1d;
      srf->gbl->vdt[indx+1][0][1] =  dttang*nrm[0]*b.vdiag1d;
      srf->gbl->vdt[indx+1][1][0] =  dtnorm*nrm[0]*b.vdiag1d;
      srf->gbl->vdt[indx+1][1][1] =  dtnorm*nrm[1]*b.vdiag1d;
               
      if (b.sm) {
         srf->gbl->sdt[indx][0][0] = -dttang*nrm[1];
         srf->gbl->sdt[indx][0][1] =  dttang*nrm[0];
         srf->gbl->sdt[indx][1][0] =  dtnorm*nrm[0];
         srf->gbl->sdt[indx][1][1] =  dtnorm*nrm[1];
      }
   }
   
   

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
   if (b.sm > 0) {
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

#ifdef SKIP
void hp_mgrid::surf_preconditioner(int bnum) {
   int tind,sind,i,j,k,m,n,indx,indx1,indx2,*v,v0,v1,ind;
   int sign[3],msgn,sgn,side[3];
   double nrm[ND],nrmo[ND],olength;
   double normco, meshco, nrmres, tanres;
   double w[3];

   /*********************************************************************/   
   /* LINEAR COMBINATIONS OF NORMAL FLOW CORRECTION & MESH CORRECTION   */
   /* THIS IMPROVES CONDITIONING OF THE SYSTEM                            */
   /* SOMEDAY I'M GOING TO FIX THIS ONCE & FOR ALL                        */
   /*********************************************************************/

   sind = sbdry[bnum].el[0];
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   nrmo[0]  =   (vrtx[v1][1] -vrtx[v0][1]);
   nrmo[1]  =  -(vrtx[v1][0] -vrtx[v0][0]);
   olength = 1/sqrt(nrmo[0]*nrmo[0] +nrmo[1]*nrmo[1]);
   nrmo[0]  *=  olength;
   nrmo[1]  *=  olength;
   meshco = meshc[0];
   normco = normc[0];
   for(indx=0;indx<sbdry[bnum].num;++indx) {
      sind = sbdry[bnum].el[indx];
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];
      nrm[0]  =   (vrtx[v1][1] -vrtx[v0][1]);
      nrm[1]  =  -(vrtx[v1][0] -vrtx[v0][0]);
      olength = 1/sqrt(nrm[0]*nrm[0] +nrm[1]*nrm[1]);
      nrm[0]  *=  olength;
      nrm[1]  *=  olength;

      /* LET ERROR FLUX OUT OF DOMAIN */
      gf[v0][2] += srf->gbl->fo*srf->gbl->rhoav*srf->gbl->res[indx][2] +surfgf[sind][1]*lamv2;
         
      /* FORM LINEAR COMBINATIONS AT VERTEX POINT */
      nrmres = 0.5*( (nrm[0]+nrmo[0])*gf[v0][0] +(nrm[1]+nrmo[1])*gf[v0][1]);
      tanres = 0.5*(-(nrm[1]+nrmo[1])*gf[v0][0] +(nrm[0]+nrmo[0])*gf[v0][1]);
      surfgf[sind][1] += -0.5*(meshc[sind]+meshco)*surfgf[sind][2]+nr*nrmres;

      /* REDUCE FLOW NORMAL RESIDUAL */
      nrmres *= 0.5*(normc[sind] +normco);
      gf[v0][0] = -0.5*(nrm[1]+nrmo[1])*tanres +0.5*(nrm[0]+nrmo[0])*nrmres;
      gf[v0][1] =  0.5*(nrm[0]+nrmo[0])*tanres +0.5*(nrm[1]+nrmo[1])*nrmres;
      

      /* FORM LINEAR COMBINATIONS FOR SIDE */
      for(m=0;m<sm-2;++m) {
         ind = eifce+1+sind*(sm-2)+m;
                     
         /* LET ERROR FLUX OUT OF DOMAIN */            
         gf[indx+m][2] += fo*rho[0]*surfgf[ind][2] +surfgf[ind][1]*lamv2;
         
         nrmres =  nrm[0]*gf[indx+m][0] +nrm[1]*gf[indx+m][1];
         tanres = -nrm[1]*gf[indx+m][0] +nrm[0]*gf[indx+m][1];
         surfgf[ind][1] += -meshc[sind]*surfgf[ind][2] +nr*nrmres;

         /* REDUCE FLOW NORMAL RESIDUAL */            
         nrmres *= normc[sind];
         gf[indx+m][0]   = -nrm[1]*tanres +nrm[0]*nrmres;
         gf[indx+m][1]   =  nrm[0]*tanres +nrm[1]*nrmres;
      }
      
      nrmo[0] = nrm[0];
      nrmo[1] = nrm[1];
      meshco = meshc[sind];
      normco = normc[sind];
   }
   
   /* LET ERROR FLUX OUT OF DOMAIN */         
   gf[v0][2] += fo*rho[0]*surfgf[eifce][2] +surfgf[eifce][1]*lamv2;
   nrmres =  nrm[0]*gf[v1][0] +nrm[1]*gf[v1][1];
   tanres = -nrm[1]*gf[v1][0] +nrm[0]*gf[v1][1];
   surfgf[eifce][1] += -meshc[eifce-1]*surfgf[eifce][2] +nr*nrmres;
   
   nrmres *= normc[eifce -1];
   gf[v1][0]   = -nrm[1]*tanres +nrm[0]*nrmres;
   gf[v1][1]   =  nrm[0]*tanres +nrm[1]*nrmres;   
}
#endif

void hp_mgrid::surfnstage1(int bnum) {
   int i,n;
   class surface *srf;
   
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;
   
   for(i=0;i<sbdry[bnum].num+1;++i)
      for(n=0;n<ND;++n)
         srf->gbl->vug0[i][n] = srf->vug[i][n];
         
   if (b.sm > 0) {
      for(i=0;i<sbdry[bnum].num*sm0;++i)
         for(n=0;n<ND;++n)
            srf->gbl->sug0[i][n] = srf->sug[i][n];
   }
   
   return;
}

void hp_mgrid::surfnstage2(int bnum, int stage) {
   int i,m,n,indx,indx1;
   class surface *srf;
   
   srf = static_cast<class surface *>(sbdry[bnum].misc);
   if (srf == NULL) return;
      
   for(i=0;i<sbdry[bnum].num+1;++i) {
      for(n=0;n<ND;++n)
         srf->vug[i][n] = srf->gbl->vug0[i][n] -alpha[stage]*srf->gbl->vres[i][n];
   }

   if (b.sm > 0) {
      indx = 0;
      indx1 = 0;
      for(i=0;i<sbdry[bnum].num;++i) {
         for(m=0;m<b.sm;++m) {
            for(n=0;n<ND;++n)
               srf->sug[indx1][n] = srf->gbl->sug0[indx1][n] -alpha[stage]*srf->gbl->sres[indx][n];
            ++indx;
            ++indx1;
         }
         indx1 += sm0 -b.sm;
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
      if (srf == NULL) return;
      
      indx = 0;
      for(i=0;i<sbdry[bnum].num;++i) {
         sind = sbdry[bnum].el[i];
         v0 = svrtx[sind][0];
         for(n=0;n<ND;++n)
            vrtx[v0][n] = srf->vug[i][n];
         
         for(m=0;m<b.sm;++m) {
            for(n=0;n<ND;++n)
               binfo[bnum][indx].curv[n] = srf->sug[indx][n];
            ++indx;
         }
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
            for(m=0;m<b.sm;++m)
               for(n=0;n<ND;++n)
                  tgt->sbuff[vbnum][count++] = srf->sug[i*b.sm +m][n];
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
         for(m=0;m<b.sm;++m) {
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
      if (srf == NULL) return;
      
      indx = 0;
      for(i=0;i<sbdry[bnum].num;++i) {
         sind = sbdry[bnum].el[i];
         v0 = svrtx[sind][0];
         for(n=0;n<ND;++n)
            srf->vug[i][n] = vrtx[v0][n];
         
         for(m=0;m<b.sm;++m) {
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

      if (b.sm > 0) {
         indx = 0;
         indx1 = 0;
         for(i=0;i<sbdry[bnum].num;++i) {
            for (j=0;j<b.sm;++j) {
               for(n=0;n<ND;++n)
                     srf->gbl->sres0[indx][n] = srf->gbl->sres[indx1][n];
               ++indx;
               ++indx1;
            }
            indx1 += b.p;
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

    if(b.p > 1) {
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
      if (srf == NULL) return;
      
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

   
