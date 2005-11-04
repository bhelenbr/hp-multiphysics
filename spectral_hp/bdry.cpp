#include"hp_mgrid.h"
#include<myblas.h>

      for(i=0;i<nsbd;++i)
         binfopv3[i] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
        
        /* UNSTEADY SOURCE TERMS (NEEDED ON FINE MESH ONLY) */
   for(i=0;i<TMSTORE;++i) {
      for(j=0;j<nsbd;++j)
         store->binfobd[i][j] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
   }
   
   for(j=0;j<nsbd;++j)
      store->dbinfodt[j] = new struct bistruct[maxsbel+1 +maxsbel*sm0];   
      
   /* ALLOCATE STORAGE FOR BOUNDARIES */
   for(i=0;i<nsbd;++i)
      binfo[i] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
      
struct bistruct *tri_hp::binfowk[TMADAPT][MAXSB]; // STORAGE FOR UNSTEADY ADAPTATION BOUNDARY BD INFO
 FLT tri_hp::bdwk(TMADAPT,NV)(MXGP,MXGP); // WORK FOR ADAPTATION
     

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION    */
/*************************************************/
void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]);

#ifdef OUTF_STRESS
double df1d(int n, double x, double y, int dir);
#endif

#ifdef DROP
extern FLT dydt;
#endif

void hp_mgrid::setinflow() {
    int i,j,k,m,n,indx,v0,v1,info;
    FLT x,y,mvel[ND];
    int sind;
   char uplo[] = "U";

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         /* INFLOW BOUNDARIES */
         /* SET VERTEX VALUES OF U,V */   
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
   
            x = vrtx(v0)(0);      
            y = vrtx(v0)(1);
            ug.v(v0,0) = (*(hp_gbl->func))(0,x,y);
            ug.v(v0,1) = (*(hp_gbl->func))(1,x,y);
         }
         v0 = sd(sind).vrtx(1);
         x = vrtx(v0)(0);      
         y = vrtx(v0)(1);
         ug.v(v0,0) = (*(hp_gbl->func))(0,x,y);
         ug.v(v0,1) = (*(hp_gbl->func))(1,x,y);
         
         /**********************************/   
         /* SET SIDE VALUES & FLUXES */
         /**********************************/
         /* ZERO FLUX FOR FIRST VERTEX */
         binfo[i][0].flx[2] = 0.0;
            
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            v1 = sd(sind).vrtx(1);
            
            if (sbdry[i].type&CURV_MASK) {
               crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(0,0),&dcrd(n,0)(0,0));
               
               crdtocht1d(sind,dvrtdt,hp_gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(1,0));
            }
            else {
               for(n=0;n<ND;++n) {
                  basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));
                  
                  for(k=0;k<basis::tri(log2p).gpx;++k)
                     dcrd(n,0)(0,k) = 0.5*(vrtx(v1)(n)-vrtx(v0)(n));
               
                  basis::tri(log2p).proj1d(dvrtdt[v0][n],dvrtdt[v1][n],&crd(n)(1,0));
               }
            }

            if (basis::tri(log2p).sm) {
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(ug.v(v0,n),ug.v(v1,n),&res(n)(0,0));
         
               for(k=0;k<basis::tri(log2p).gpx; ++k)
                  for(n=0;n<ND;++n)
                     res(n)(0,k) -= (*(hp_gbl->func))(n,crd(0)(0,k),crd(1)(0,k));
                     
               for(n=0;n<ND;++n)
                  basis::tri(log2p).intgrt1d(&lf(n)(0),&res(n)(0,0));
         
               indx = sind*sm0;
               for(n=0;n<ND;++n) {
                  PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
                  for(m=0;m<basis::tri(log2p).sm;++m) 
                     ug.s(indx+m)(n) = -lf(n)(2+m);
               }
            }
            
            /* NOW SET FLUXES */
            ugtouht1d(sind);
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj1d(&uht(n)(0),&u(n)(0,0));

            for(k=0;k<basis::tri(log2p).gpx;++k) {
               for(n=0;n<ND;++n)
                  mvel[n] = sim::bd[0]*crd(n)(0,k) +crd(n)(1,k);
#ifdef DROP
               mvel[1] += dydt;
#endif
               
               res(2)(0,k) = hp_gbl->rho*RAD1D(k)*((u(0)(0,k) -mvel[0])*dcrd(1,0)(0,k) -(u(1)(0,k) -mvel[1])*dcrd(0,0)(0,k));
            }
            
            basis::tri(log2p).intgrt1d(&lf(0)(0),&res(2)(0,0));
            
            indx = j*(basis::tri(log2p).sm +1);
            binfo[i][indx++].flx[2] += lf(0)(0);
            for(m=0;m<basis::tri(log2p).sm;++m)
               binfo[i][indx++].flx[2] = lf(0)(m+2);
            binfo[i][indx].flx[2] = lf(0)(1);
         }
      }
      
#ifdef OUTF_STRESS
      double nrm[ND];
      
      if (sbdry[i].type&OUTF_MASK) {
      
         for(n=0;n<ND;++n)
            binfo[i][0].flx[n] = 0.0;
            
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            
            if (sbdry[i].type&CURV_MASK) {
               crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(0,0),&dcrd(n,0)(0,0));
            }
            else {
               v0 = sd(sind).vrtx(0);
               v1 = sd(sind).vrtx(1); 
               for(n=0;n<ND;++n) {
                  basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));
                  
                  for(k=0;k<basis::tri(log2p).gpx;++k)
                     dcrd(n,0)(0,k) = 0.5*(vrtx(v1)(n)-vrtx(v0)(n));
               }
            }
            
            for(k=0;k<basis::tri(log2p).gpx;++k) {
               nrm[0] = dcrd(1,0)(0,k);
               nrm[1] = -dcrd(0,0)(0,k);
               
               res(0)(0,k) = -hp_gbl->mu*RAD1D(k)*(nrm[0]*(2.*df1d(0,crd(0)(0,k),crd(1)(0,k),0))
                             +nrm[1]*(df1d(0,crd(0)(0,k),crd(1)(0,k),1) +df1d(1,crd(0)(0,k),crd(1)(0,k),0)));
               res(1)(0,k) = -hp_gbl->mu*RAD1D(k)*(nrm[1]*(2.*df1d(1,crd(0)(0,k),crd(1)(0,k),1))
                             +nrm[0]*(df1d(0,crd(0)(0,k),crd(1)(0,k),1) +df1d(1,crd(0)(0,k),crd(1)(0,k),0)));
            }
            
            for(n=0;n<ND;++n)
               basis::tri(log2p).intgrt1d(&lf(n)(0),&res(n)(0,0));
               
            indx = j*(basis::tri(log2p).sm +1);
            for(n=0;n<ND;++n)
               binfo[i][indx].flx[n] += lf(n)(0);
            ++indx;
            
            for(m=0;m<basis::tri(log2p).sm;++m) {
               for(n=0;n<ND;++n) {
                  binfo[i][indx].flx[n] = lf(n)(m+2);
               }
               ++indx;
            }
            for(n=0;n<ND;++n) 
               binfo[i][indx].flx[n] = lf(n)(1);
         }
      }
#endif
   }
   
   return;
}

void hp_mgrid::addbflux(int mgrid) {
   int i,j,k,n,indx,indx1;
   int sind,v0,v1;
   FLT gam, nrm[ND], wl[NV], wr[NV];
   FLT mvel[ND] = {0.0, 0.0};
   

   setinflow();
   
   /* THESE ARE SOURCE TERMS WHICH CHANGE WITH THE SOLUTION */
   /* MUST BE UPDATED DURING MGRID FOR GOOD CONVERGENCE */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         indx = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            sind=sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            hp_gbl->res.v(v0,2) += binfo[i][indx++].flx[2];
            indx1 = sind*basis::tri(log2p).sm;
            for(k=0;k<basis::tri(log2p).sm;++k)
               hp_gbl->res.s(indx1++)(2) += binfo[i][indx++].flx[2];
         }
         v0 = sd(sind).vrtx(1);
         hp_gbl->res.v(v0,2) += binfo[i][indx].flx[2];
      }
      
#ifdef OUTF_STRESS
      if (sbdry[i].type&OUTF_MASK) {
         /* ADD SURFACE TENSION SOURCE TERM */
         indx = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            sind=sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for(n=0;n<ND;++n)
               hp_gbl->res.v(v0,n) += binfo[i][indx].flx[n];
            ++indx;
            indx1 = sind*basis::tri(log2p).sm;
            for(k=0;k<basis::tri(log2p).sm;++k) {
               for(n=0;n<ND;++n)
                  hp_gbl->res.s(indx1)(n) += binfo[i][indx].flx[n];
               ++indx;
               ++indx1;
            }
         }
         v0 = sd(sind).vrtx(1);
         for(n=0;n<ND;++n)
            hp_gbl->res.v(v0,n) += binfo[i][indx].flx[n];
      }
#endif

      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK)) {
         /* ADD SURFACE TENSION SOURCE TERM */
         indx = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            sind=sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) += binfo[i][indx].flx[n];
            ++indx;
            indx1 = sind*basis::tri(log2p).sm;
            for(k=0;k<basis::tri(log2p).sm;++k) {
               for(n=0;n<NV;++n)
                  hp_gbl->res.s(indx1)(n) += binfo[i][indx].flx[n];
               ++indx;
               ++indx1;
            }
         }
         v0 = sd(sind).vrtx(1);
         for(n=0;n<NV;++n)
            hp_gbl->res.v(v0,n) += binfo[i][indx].flx[n];
      }


      /* OUTFLOW BOUNDARY CONDITION    */
      if (sbdry[i].type&OUTF_MASK) {
         indx = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            v1 = sd(sind).vrtx(1);
            
            if (sbdry[i].type&CURV_MASK) {
               crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(0,0),&dcrd(n,0)(0,0));
               
               crdtocht1d(sind,dvrtdt,hp_gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(1,0));
            }
            else {
               for(n=0;n<ND;++n) {
                  basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));
                  
                  for(k=0;k<basis::tri(log2p).gpx;++k)
                     dcrd(n,0)(0,k) = 0.5*(vrtx(v1)(n)-vrtx(v0)(n));
               
                  basis::tri(log2p).proj1d(dvrtdt[v0][n],dvrtdt[v1][n],&crd(n)(1,0));
               }
            }
            
            ugtouht1d(sind);
            for(n=0;n<NV;++n)
               basis::tri(log2p).proj1d(&uht(n)(0),&u(n)(0,0));
            
            gam = hp_gbl->rhoi*hp_gbl->tprcn[sd(sind).tri(0)][0][0]/hp_gbl->tprcn[sd(sind).tri(0)][NV-1][NV-1];
            for(k=0;k<basis::tri(log2p).gpx;++k) {
               for(n=0;n<NV;++n) {
                  wl[n] = u(n)(0,k);
                  wr[n] = (hp_gbl->func)(n,crd(0)(0,k),crd(1)(0,k));
               }
               nrm[0] = dcrd(1,0)(0,k);
               nrm[1] = -dcrd(0,0)(0,k);

               for(n=0;n<ND;++n)
                  mvel[n] = sim::bd[0]*crd(n)(0,k) +crd(n)(1,k);
#ifdef DROP
               mvel[1] += dydt;
#endif
                  
               if (!charyes)
                  wl[2] = wr[2];
               else 
                  chrctr(hp_gbl->rho,gam,wl,wr,nrm,mvel);
               
               res(2)(0,k) = hp_gbl->rho*RAD1D(k)*((wl[0] -mvel[0])*nrm[0] +(wl[1] -mvel[1])*nrm[1]);
#ifndef INERTIALESS
               res(0)(0,k) = res(2)(0,k)*wl[0] +wl[2]*RAD1D(k)*nrm[0];
               res(1)(0,k) = res(2)(0,k)*wl[1] +wl[2]*RAD1D(k)*nrm[1];
#else
               res(0)(0,k) = wl[2]*RAD1D(k)*nrm[0];
               res(1)(0,k) = wl[2]*RAD1D(k)*nrm[1]; 
#endif
            }
            
            for(n=0;n<NV;++n)
               basis::tri(log2p).intgrt1d(&lf(n)(0),&res(n)(0,0));
            
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) += lf(n)(0);

            for(n=0;n<NV;++n)
               hp_gbl->res.v(v1,n) += lf(n)(1);
            
            indx1 = sind*basis::tri(log2p).sm;
            indx = 2;
            for(k=0;k<basis::tri(log2p).sm;++k) {
               for(n=0;n<NV;++n)
                  hp_gbl->res.s(indx1)(n) += lf(n)(indx);
               ++indx1;
               ++indx;
            }
         }
      }

      /* INVISCID WALL BOUNDARY CONDITION */
      if (sbdry[i].type&EULR_MASK) {
         indx = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            v1 = sd(sind).vrtx(1);
            
            if (sbdry[i].type&CURV_MASK) {
               crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(0,0),&dcrd(n,0)(0,0));
               
               crdtocht1d(sind,dvrtdt,hp_gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(1,0));
            }
            else {
               for(n=0;n<ND;++n) {
                  basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));
                  
                  for(k=0;k<basis::tri(log2p).gpx;++k)
                     dcrd(n,0)(0,k) = 0.5*(vrtx(v1)(n)-vrtx(v0)(n));
               
                  basis::tri(log2p).proj1d(dvrtdt[v0][n],dvrtdt[v1][n],&crd(n)(1,0));
               }
            }
            
            ugtouht1d(sind);
            basis::tri(log2p).proj1d(&uht(2)(0)&u(2)(0,0));
               
            for(k=0;k<basis::tri(log2p).gpx;++k) {
               wl[0] = (hp_gbl->func)(0,crd(0)(0,k),crd(1)(0,k));
               wl[1] = (hp_gbl->func)(1,crd(0)(0,k),crd(1)(0,k));
               for(n=0;n<ND;++n)
                  mvel[n] = sim::bd[0]*crd(n)(0,k) +crd(n)(1,k);
               res(2)(0,k) = hp_gbl->rho*RAD1D(k)*((wl[0] -mvel[0])*dcrd(1,0)(0,k) -(wl[1] -mvel[1])*dcrd(0,0)(0,k));
               res(0)(0,k) = res(2)(0,k)*wl[0] +u(2)(0,k)*RAD1D(k)*dcrd(1,0)(0,k);
               res(1)(0,k) = res(2)(0,k)*wl[1] -u(2)(0,k)*RAD1D(k)*dcrd(0,0)(0,k);
            }
            
            for(n=0;n<NV;++n)
               basis::tri(log2p).intgrt1d(&lf(n)(0),&res(n)(0,0));
            
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) += lf(n)(0);

            for(n=0;n<NV;++n)
               hp_gbl->res.v(v1,n) += lf(n)(1);
            
            indx1 = sind*basis::tri(log2p).sm;
            indx = 2;
            for(k=0;k<basis::tri(log2p).sm;++k) {
               for(n=0;n<NV;++n)
                  hp_gbl->res.s(indx1)(n) += lf(n)(indx);
               ++indx1;
               ++indx;
            }
         }
      }
   }
   
   return;
}

void hp_mgrid::bdry_vsnd() {
   int i,j,n,sind,count,v0,bnum;
   class mesh *tgt;
   
   /* SEND VERTEX INFO FOR Y_DIR*/
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
         }
         v0 = sd(sind).vrtx(1);
         for (n=0;n<NV;++n) 
            tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
      }
      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for (n=0;n<ND;++n) 
               tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
         }
         v0 = sd(sind).vrtx(1);
         for (n=0;n<ND;++n) 
            tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
      }
   }
   
   return;
}

void hp_mgrid::bdry_mp() {
   int i,j,n,sind,count,v0,bnum;
   class mesh *tgt;
   
   /* THIS PART IS TO RECEIVE AND ZERO FOR VERTICES */
   /* RECEIVE VRTX MESSAGES */
   /* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         for(n=0;n<NV;++n)
            hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
      }

      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<ND;++n)
               hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         for(n=0;n<ND;++n)
            hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
      }         
   }
   
   /* SEND VERTEX INFO FOR X_DIR*/
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
         }
         v0 = sd(sind).vrtx(1);
         for (n=0;n<NV;++n) 
            tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
      }
   }
   
   return;
}


void hp_mgrid::bdry_vrcvandzero() {
    int i,j,n;
    int sind,v0,count;
   
   /* THIS PART IS TO RECEIVE AND ZERO FOR VERTICES */
   /* RECEIVE VRTX MESSAGES */
   /* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         for(n=0;n<NV;++n)
            hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
      }         
   }

   /* APPLY VRTX DIRICHLET CONDITIONS TO RES */
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].type&INFL_MASK) {
         for(j=0;j<vbdry[i].num;++j) {
            v0 = vbdry[i].el[j];
            for(n=0;n<ND;++n)
               hp_gbl->res.v(v0,n) = 0.0;
         }
      }
      
      if (vbdry[i].type&SYMM_MASK) {
         for(j=0;j<vbdry[i].num;++j) {
            v0 = vbdry[i].el[j];
            hp_gbl->res.v(v0,0) = 0.0;
         }
      }
   }

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for(n=0;n<ND;++n)
               hp_gbl->res.v(v0,n) = 0.0;
         }
         v0 = sd(sind).vrtx(1);
         for(n=0;n<ND;++n)
            hp_gbl->res.v(v0,n) = 0.0;
      }
      
      if (sbdry[i].type&SYMM_MASK) {
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            hp_gbl->res.v(v0,0) = 0.0;
         }
         v0 = sd(sind).vrtx(1);
         hp_gbl->res.v(v0,0) = 0.0;
      }
   }
   
   return;
}

void hp_mgrid::bdry_ssnd(int mode) {
   int i,j,n,count,indx,bnum;
   class mesh *tgt;
   
   /* SEND SIDE INFO */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMX_MASK +COMY_MASK)) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            indx = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = hp_gbl->res.s(indx)(n);
         }
      } 

      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            indx = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            for (n=0;n<ND;++n) 
               tgt->sbuff[bnum][count++] = hp_gbl->res.s(indx)(n);
//            printf("Sent %d %f %f\n",count,tgt->sbuff[bnum][count-2],tgt->sbuff[bnum][count-1]);
         }
      }         
   }
   
   return;
}

   
void hp_mgrid::bdry_srcvandzero(int mode) {
    int i,j,n;
    int sind,count,indx,sign;
   
   sign = (mode % 2 ? -1 : 1);
   
   /* THIS PART TO RECIEVE AND ZERO FOR SIDES */
   /* RECEIVE P'TH SIDE MODE MESSAGES */
   /* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMX_MASK +COMY_MASK)) {
         count = 0;
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            indx = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            for(n=0;n<NV;++n)
               hp_gbl->res.s(indx)(n) = 0.5*(hp_gbl->res.s(indx)(n) +sign*sbuff[i][count++]);
         }
      }
      
      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            indx = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            for(n=0;n<ND;++n)
               hp_gbl->res.s(indx)(n) = 0.5*(hp_gbl->res.s(indx)(n) +sign*sbuff[i][count++]);
//            printf("Recieved %d %f %f\n",count,sbuff[i][count-2],sbuff[i][count-1]);

         }
      }         
   }

   /* APPLY SIDE DIRICHLET CONDITIONS TO MODE */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            for(n=0;n<ND;++n)
               hp_gbl->res.s(sind)(n) = 0.0;
         }
      }
      
      if (sbdry[i].type&SYMM_MASK) {
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            hp_gbl->res.s(sind)(0) = 0.0;
         }
      }
   }
   
   return;
}


void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]) {
   FLT ul,vl,ur,vr,pl,pr,cl,cr,rhoi;
   FLT u,um,v,c,den,lam0,lam1,lam2,uvp[3],mag;
   
   rhoi = 1./rho;

   /* CHARACTERISTIC FAR-FIELD B.C. */   
   mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
   
   norm[0] /= mag;
   norm[1] /= mag;
   
   ul =  wl[0]*norm[0] +wl[1]*norm[1];
   vl = -wl[0]*norm[1] +wl[1]*norm[0];
   pl =  wl[2];
      
   /* DEPENDENT ON FREESTREAM CONDITIONS */
   ur =  wr[0]*norm[0] +wr[1]*norm[1];
   vr = -wr[0]*norm[1] +wr[1]*norm[0];
   pr =  wr[2];
      
   um = mv[0]*norm[0] +mv[1]*norm[1];
   
   cl = sqrt((ul-.5*um)*(ul-.5*um) +gam);
   cr = sqrt((ur-.5*um)*(ur-.5*um) +gam);
   c = 0.5*(cl+cr);
   u = 0.5*(ul+ur);
   v = 0.5*(vl+vr);
   
   den = 1./(2*c);
   lam0 = u -um;
   lam1 = u-.5*um +c; /* always positive */
   lam2 = u-.5*um -c; /* always negative */
      
   /* PERFORM CHARACTERISTIC SWAP */
   /* BASED ON LINEARIZATION AROUND UL,VL,PL */
   uvp[0] = ((pl-pr)*rhoi +(ul*lam1 -ur*lam2))*den;
   if (lam0 > 0.0)
      uvp[1] = v*((pr-pl)*rhoi +lam2*(ur-ul))*den/(lam0-lam2) +vl;
   else
      uvp[1] = v*((pr-pl)*rhoi +lam1*(ur-ul))*den/(lam0-lam1) +vr;
   uvp[2] = (rho*(ul -ur)*gam - lam2*pl +lam1*pr)*den;
   
   /* CHANGE BACK TO X,Y COORDINATES */
   wl[0] =  uvp[0]*norm[0] -uvp[1]*norm[1];
   wl[1] =  uvp[0]*norm[1] +uvp[1]*norm[0];
   wl[2] =  uvp[2];

   /* SHOULDN'T CHANGE NORM */   
   norm[0] *= mag;
   norm[1] *= mag;
      
   return;
 
}


void hp_mgrid::matchboundaries1() {
   int i,j,m,n,v0,sind,bnum,count,indx,indx1;
   class mesh *tgt;
   
   /* SEND POSITIONS, NUMBER, AND SOLUTION TO ADJACENT MESHES */
   /* SENDS EVERYTHING FOR ALL BOUNDARIES */
   /* SORT OUT WHAT TO DO AFTER RECEIVE */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & ALLD_MP) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         
         /* SEND NUMBER AS FLOAT? */
         tgt->sbuff[bnum][count++] = sbdry(i)->nel;
         
         for(j=0;j<sbdry(i)->nel;++j) {
            /* SEND VERTEX INFO */
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for(n=0;n<ND;++n)
               tgt->sbuff[bnum][count++] = vrtx(v0)(n);
            for(n=0;n<NV;++n)
               tgt->sbuff[bnum][count++] = ug.v(v0,n);
               
            /* SEND SIDE INFO */
            indx = sind*basis::tri(log2p).sm;
            indx1 = j*basis::tri(log2p).sm;
            for(m=0;m<basis::tri(log2p).sm;++m) {
               for (n=0;n<ND;++n) 
                  tgt->sbuff[bnum][count++] = binfo[i][indx1+m].curv[n];
               for (n=0;n<NV;++n) 
                  tgt->sbuff[bnum][count++] = ug.s(indx+m)(n);
            }
         }
         v0 = sd(sind).vrtx(1);
         for(n=0;n<ND;++n)
            tgt->sbuff[bnum][count++] = vrtx(v0)(n);
         for(n=0;n<NV;++n)
            tgt->sbuff[bnum][count++] = ug.v(v0,n);
      } 
   }

   return;
   
}

void hp_mgrid::matchboundaries2() {
   int i,j,m,n,v0,num,sind,bnum,count,indx,indx1,msgn;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & PRDX_MASK) {
         count = 0;
         
         /*	RECV NUMBER */
         num = static_cast<int>(sbuff[i][count++]);
         if (num != sbdry(i)->nel) {
            printf("non matching number of boundaries %d\n",sbdry[i].type);
            exit(1);
         }
         
         /* RECV INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            ++count; // SKIP X
            for(n=1;n<ND;++n)
               vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]);
            for(n=0;n<NV;++n)
               ug.v(v0,n) = 0.5*(ug.v(v0,n) +sbuff[i][count++]);
               
            
            indx = sind*basis::tri(log2p).sm;
            indx1 = j*basis::tri(log2p).sm;
            msgn = 1;
            for(m=0;m<basis::tri(log2p).sm;++m) {
               ++count; // SKIP X
               for(n=1;n<ND;++n)
                  binfo[i][indx1+m].curv[n] = 0.5*(binfo[i][indx1+m].curv[n] +msgn*sbuff[i][count++]);
               for(n=0;n<NV;++n)
                  ug.s(indx+m)(n) = 0.5*(ug.s(indx+m)(n) +msgn*sbuff[i][count++]);
               msgn *= -1;
            } 
         }
         v0 = sd(sind).vrtx(0);
         ++count; // SKIP X
         for(n=1;n<ND;++n)
            vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]);
         for(n=0;n<NV;++n)
            ug.v(v0,n) = 0.5*(ug.v(v0,n) +sbuff[i][count++]);
            
         continue;
      }
      if (sbdry[i].type & PRDY_MASK) {
         count = 0;
         
         /*	RECV NUMBER */
         num = static_cast<int>(sbuff[i][count++]);
         if (num != sbdry(i)->nel) {
            printf("non matching number of boundaries %d\n",sbdry[i].type);
            exit(1);
         }
         
         /* RECV INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<1;++n)
               vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]);
            ++count; // SKIP Y
            for(n=0;n<NV;++n)
               ug.v(v0,n) = 0.5*(ug.v(v0,n) +sbuff[i][count++]);
               
            
            indx = sind*basis::tri(log2p).sm;
            indx1 = j*basis::tri(log2p).sm;
            msgn = 1;
            for(m=0;m<basis::tri(log2p).sm;++m) {
               for(n=0;n<1;++n)
                  binfo[i][indx1+m].curv[n] = 0.5*(binfo[i][indx1+m].curv[n] +msgn*sbuff[i][count++]);
               ++count; // SKIP Y
               for(n=0;n<NV;++n)
                  ug.s(indx+m)(n) = 0.5*(ug.s(indx+m)(n) +msgn*sbuff[i][count++]);
               msgn *= -1;
            } 
         }
         v0 = sd(sind).vrtx(0);
         for(n=0;n<1;++n)
            vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]);
         ++count; // SKIP Y
         for(n=0;n<NV;++n)
            ug.v(v0,n) = 0.5*(ug.v(v0,n) +sbuff[i][count++]);
            
         continue;
      }
      
      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
         
         /*	RECV NUMBER */
         num = static_cast<int>(sbuff[i][count++]);
         if (num != sbdry(i)->nel) {
            printf("non matching number of boundaries %d\n",sbdry[i].type);
            exit(1);
         }
         
         /* RECV INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<ND;++n)
               vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]);
            for(n=0;n<ND;++n)
               ug.v(v0,n) = 0.5*(ug.v(v0,n) +sbuff[i][count++]);
            ++count; // SKIP P
            
            indx = sind*basis::tri(log2p).sm;
            indx1 = j*basis::tri(log2p).sm;
            msgn = 1;
            for(m=0;m<basis::tri(log2p).sm;++m) {
               for(n=0;n<ND;++n)
                  binfo[i][indx1+m].curv[n] = 0.5*(binfo[i][indx1+m].curv[n] +msgn*sbuff[i][count++]);
               for(n=0;n<ND;++n)
                  ug.s(indx+m)(n) = 0.5*(ug.s(indx+m)(n) +msgn*sbuff[i][count++]);
               ++count; // SKIP P
               msgn *= -1;
            } 
         }
         v0 = sd(sind).vrtx(0);
         for(n=0;n<ND;++n)
            vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]);
         for(n=0;n<ND;++n)
            ug.v(v0,n) = 0.5*(ug.v(v0,n) +sbuff[i][count++]);
         ++count; // SKIP P
            
         continue;
      }
      
      if (sbdry[i].type & ALLD_MP) {
         count = 0;
         
         /*	RECV NUMBER */
         num = static_cast<int>(sbuff[i][count++]);
         if (num != sbdry(i)->nel) {
            printf("non matching number of boundaries %d\n",sbdry[i].type);
            exit(1);
         }
         
         /* RECV INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<ND;++n)
               vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]);
            for(n=0;n<NV;++n)
               ug.v(v0,n) = 0.5*(ug.v(v0,n) +sbuff[i][count++]);
               
            
            indx = sind*basis::tri(log2p).sm;
            indx1 = j*basis::tri(log2p).sm;
            msgn = 1;
            for(m=0;m<basis::tri(log2p).sm;++m) {
               for(n=0;n<ND;++n)
                  binfo[bnum][indx1+m].curv[n] = 0.5*(binfo[i][indx1+m].curv[n] +msgn*sbuff[i][count++]);
               for(n=0;n<NV;++n)
                  ug.s(indx+m)(n) = 0.5*(ug.s(indx+m)(n) +msgn*sbuff[i][count++]);
               msgn *= -1;
            } 
         }
         v0 = sd(sind).vrtx(0);
         for(n=0;n<ND;++n)
            vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]);
         for(n=0;n<NV;++n)
            ug.v(v0,n) = 0.5*(ug.v(v0,n) +sbuff[i][count++]);
            
         continue;
      }
   }
   return;
}
