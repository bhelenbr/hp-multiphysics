#include"hp_mgrid.h"
#include<myblas.h>

extern FLT axext, ayext;
/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION    */
/*************************************************/
extern FLT df1d(int,FLT,FLT);
void chrctr(FLT ax, FLT ay, double wl[NV], double wr[NV], double norm[ND], double mv[ND]);

void hp_mgrid::setinflow() {
    int i,j,k,m,n,indx,v0,v1,info;
    FLT x,y;
    int sind;
   char uplo[] = "U";

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         /* INFLOW BOUNDARIES */
         /* SET VERTEX VALUES OF U,V */   
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
   
            x = vrtx[v0][0];      
            y = vrtx[v0][1];
            ug.v[v0][0] = (*(gbl->func))(0,x,y);
         }
         v0 = svrtx[sind][1];
         x = vrtx[v0][0];      
         y = vrtx[v0][1];
         ug.v[v0][0] = (*(gbl->func))(0,x,y);
         
         /**********************************/   
         /* SET SIDE VALUES & FLUXES */
         /**********************************/
         /* ZERO FLUX FOR FIRST VERTEX */
         for(n=0;n<NV;++n)
            binfo[i][0].flx[n] = 0.0;
            
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            
            if (sbdry[i].type&CURV_MASK) {
               crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  b->proj1d(cht[n],crd[n][0],dcrd[n][0][0]);
               
               crdtocht1d(sind,dvrtdt,gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  b->proj1d(cht[n],crd[n][1]);
            }
            else {
               for(n=0;n<ND;++n) {
                  b->proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
                  
                  for(k=0;k<b->gpx;++k)
                     dcrd[n][0][0][k] = 0.5*(vrtx[v1][n]-vrtx[v0][n]);
               
                  b->proj1d(dvrtdt[v0][n],dvrtdt[v1][n],crd[n][1]);
               }
            }

            if (b->sm) {
               for(n=0;n<NV;++n)
                  b->proj1d(ug.v[v0][n],ug.v[v1][n],res[n][0]);
         
               for(k=0;k<b->gpx; ++k)
                  for(n=0;n<NV;++n)
                     res[n][0][k] -= (*(gbl->func))(n,crd[0][0][k],crd[1][0][k]);
                     
               for(n=0;n<NV;++n)
                  b->intgrt1d(lf[n],res[n][0]);
         
               indx = sind*sm0;
               for(n=0;n<NV;++n) {
                  PBTRS(uplo,b->sm,b->sbwth,1,b->sdiag1d[0],b->sbwth+1,&lf[n][2],b->sm,info);
                  for(m=0;m<b->sm;++m) 
                     ug.s[indx+m][n] = -lf[n][2+m];
               }
            }
         }
      }
      
      if (sbdry[i].type&OUTF_MASK) {
         for(n=0;n<NV;++n)
            binfo[i][0].flx[n] = 0.0;
            
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            
           
            for(n=0;n<ND;++n) {
               b->proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
               
               for(k=0;k<b->gpx;++k)
                  dcrd[n][0][0][k] = 0.5*(vrtx[v1][n]-vrtx[v0][n]);
            }
            
            /* NOW SET FLUXES */
            for(k=0;k<b->gpx;++k) 
               res[0][0][k] = -gbl->mu*RAD1D(k)*(df1d(0,crd[0][0][k],crd[1][0][k])*dcrd[1][0][0][k]);
            
            b->intgrt1d(lf[0],res[0][0]);
            
            indx = j*(b->sm +1);
            binfo[i][indx++].flx[0] += lf[0][0];
            for(m=0;m<b->sm;++m)
               binfo[i][indx++].flx[0] = lf[0][m+2];
            binfo[i][indx].flx[0] = lf[0][1];
         }
      }
   }
   
   return;
}

void hp_mgrid::addbflux(int mgrid) {
    int i,j,k,n,indx,indx1;
    int sind,v0,v1;
    FLT nrm[ND], wl[NV], wr[NV];
   FLT mvel[ND] = {0.0, 0.0};
   
   /***********************************/
   /* ADD SOURCE TERMS ON FINEST MESH */
   /***********************************/
   if(!mgrid) {
//      setinflow();  //TEMPORARY

      for(i=0;i<nsbd;++i) {         
         if (sbdry[i].type&OUTF_MASK) {
            /* ALLOWS FOR APPLIED STRESS ON BOUNDARY */
            indx = 0;
            for(j=0;j<sbdry[i].num;++j) {
               sind=sbdry[i].el[j];
               v0 = svrtx[sind][0];
               indx1 = sind*b->sm;
               for(n=0;n<NV;++n)
                  gbl->res.v[v0][n] += binfo[i][indx].flx[n];
               ++indx;
               for(k=0;k<b->sm;++k) {
                  for(n=0;n<NV;++n)
                     gbl->res.s[indx1][n] += binfo[i][indx].flx[n];
                  ++indx;
                  ++indx1;
               }
            }
            v0 = svrtx[sind][1];
            for(n=0;n<NV;++n)
               gbl->res.v[v0][n] += binfo[i][indx].flx[n];
         }
      }
   }
   
   
   /* THESE ARE SOURCE TERMS WHICH CHANGE WITH THE SOLUTION */
   /* MUST BE UPDATED DURING MGRID FOR GOOD CONVERGENCE */
   for(i=0;i<nsbd;++i) {

      /* OUTFLOW BOUNDARY CONDITION    */
      if (sbdry[i].type&OUTF_MASK) {
         indx = 0;
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            
            if (sbdry[i].type&CURV_MASK) {
               crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  b->proj1d(cht[n],crd[n][0],dcrd[n][0][0]);
               
               crdtocht1d(sind,dvrtdt,gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  b->proj1d(cht[n],crd[n][1]);
            }
            else {
               for(n=0;n<ND;++n) {
                  b->proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
                  
                  for(k=0;k<b->gpx;++k)
                     dcrd[n][0][0][k] = 0.5*(vrtx[v1][n]-vrtx[v0][n]);
               
                  b->proj1d(dvrtdt[v0][n],dvrtdt[v1][n],crd[n][1]);
               }
            }
            
            ugtouht1d(sind);
            for(n=0;n<NV;++n)
               b->proj1d(uht[n],u[n][0]);
            
            for(k=0;k<b->gpx;++k) {
               for(n=0;n<NV;++n) {
                  wl[n] = u[n][0][k];
                  wr[n] = (gbl->func)(n,crd[0][0][k],crd[1][0][k]);
               }
               nrm[0] = dcrd[1][0][0][k];
               nrm[1] = -dcrd[0][0][0][k];

               for(n=0;n<ND;++n)
                  mvel[n] = bd[0]*crd[n][0][k] +crd[n][1][k];
                  
               chrctr(axext,ayext,wl,wr,nrm,mvel);
                                 
               res[0][0][k] = wl[0]*RAD1D(k)*((axext -mvel[0])*nrm[0] +(ayext -mvel[1])*nrm[1]);
            }
            
            for(n=0;n<NV;++n)
               b->intgrt1d(lf[n],res[n][0]);
            
            for(n=0;n<NV;++n)
               gbl->res.v[v0][n] += lf[n][0];

            for(n=0;n<NV;++n)
               gbl->res.v[v1][n] += lf[n][1];
            
            indx1 = sind*b->sm;
            indx = 2;
            for(k=0;k<b->sm;++k) {
               for(n=0;n<NV;++n)
                  gbl->res.s[indx1][n] += lf[n][indx];
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
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = gbl->res.v[v0][n];
         }
         v0 = svrtx[sind][1];
         for (n=0;n<NV;++n) 
            tgt->sbuff[bnum][count++] = gbl->res.v[v0][n];
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
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            for(n=0;n<NV;++n)
               gbl->res.v[v0][n] = 0.5*(gbl->res.v[v0][n] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         for(n=0;n<NV;++n)
            gbl->res.v[v0][n] = 0.5*(gbl->res.v[v0][n] +sbuff[i][count++]);
      } 
   }
   
   /* SEND VERTEX INFO FOR X_DIR*/
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = gbl->res.v[v0][n];
         }
         v0 = svrtx[sind][1];
         for (n=0;n<NV;++n) 
            tgt->sbuff[bnum][count++] = gbl->res.v[v0][n];
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
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            for(n=0;n<NV;++n)
               gbl->res.v[v0][n] = 0.5*(gbl->res.v[v0][n] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         for(n=0;n<NV;++n)
            gbl->res.v[v0][n] = 0.5*(gbl->res.v[v0][n] +sbuff[i][count++]);
      }         
   }

   /* APPLY VRTX DIRICHLET CONDITIONS TO RES */
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].type&INFL_MASK) {
         for(j=0;j<vbdry[i].num;++j) {
            v0 = vbdry[i].el[j];
            for(n=0;n<NV;++n)
               gbl->res.v[v0][n] = 0.0;
         }
      }
   }

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            for(n=0;n<NV;++n)
               gbl->res.v[v0][n] = 0.0;
         }
         v0 = svrtx[sind][1];
         for(n=0;n<NV;++n)
            gbl->res.v[v0][n] = 0.0;
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
         for(j=0;j<sbdry[i].num;++j) {
            indx = sbdry[i].el[j]*b->sm +mode;
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = gbl->res.s[indx][n];
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
         for(j=sbdry[i].num-1;j>=0;--j) {
            indx = sbdry[i].el[j]*b->sm +mode;
            for(n=0;n<NV;++n)
               gbl->res.s[indx][n] = 0.5*(gbl->res.s[indx][n] +sign*sbuff[i][count++]);
         }
      }
   }

   /* APPLY SIDE DIRICHLET CONDITIONS TO MODE */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j]*b->sm +mode;
            for(n=0;n<NV;++n)
               gbl->res.s[sind][n] = 0.0;
         }
      }
   }
   
   return;
}


void chrctr(FLT ax, FLT ay,double wl[NV], double wr[NV], double norm[ND], double mv[ND]) {
   FLT ul;
   FLT um,lam0,mag;
   
   /* CHARACTERISTIC FAR-FIELD B.C. */   
   mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
   
   norm[0] /= mag;
   norm[1] /= mag;
   
   ul =  ax*norm[0] +ay*norm[1];      
   um = mv[0]*norm[0] +mv[1]*norm[1];

   lam0 = ul-um;

   if (lam0 > 0.0)
      wl[0] = wl[0];
   else
      wl[0] = wr[0];

   /* SHOULDN'T CHANGE NORM */   
   norm[0] *= mag;
   norm[1] *= mag;
   
   return;
 
}

void hp_mgrid::bdrycheck1() {
   int i,j,n,v0,count,bnum,sind;
   class mesh *tgt;
   
   /* SEND SIDE INFO */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMX_MASK +COMY_MASK +IFCE_MASK)) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            for (n=0;n<ND;++n) 
               tgt->sbuff[bnum][count++] = vrtx[v0][n];
         }
         v0 = svrtx[sind][1];
         for (n=0;n<ND;++n) 
            tgt->sbuff[bnum][count++] = vrtx[v0][n]; 
      }    
   }
   
   return;
}

void hp_mgrid::bdrycheck2() {
    int i,j,n,v0;
    int sind,count;
      
   /* THIS PART TO RECIEVE AND ZERO FOR SIDES */
   /* RECEIVE P'TH SIDE MODE MESSAGES */
   /* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & PRDX_MASK) {
         count = 1;
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            vrtx[v0][1] = 0.5*(vrtx[v0][1] +sbuff[i][count++]);
            ++count;
         }
         v0 = svrtx[sind][0];
         vrtx[v0][1] = 0.5*(vrtx[v0][1] +sbuff[i][count++]); 
         continue;   
      }
      
      if (sbdry[i].type & PRDY_MASK) {
         count = 0;
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            vrtx[v0][0] = 0.5*(vrtx[v0][0] +sbuff[i][count++]);
            ++count;
         }
         v0 = svrtx[sind][0];
         vrtx[v0][0] = 0.5*(vrtx[v0][0] +sbuff[i][count++]); 
         continue;   
      }
      
      if (sbdry[i].type & (COMX_MASK +COMY_MASK +IFCE_MASK)) {
         count = 0;
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            for(n=0;n<ND;++n)
               vrtx[v0][n] = 0.5*(vrtx[v0][n] +sbuff[i][count++]); 
         }
      }
   }
   
   return;
}




