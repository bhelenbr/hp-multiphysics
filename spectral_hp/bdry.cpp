#include"hp_mgrid.h"
#include<myblas.h>

/*************************************************/
/*	SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/*	(THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION	 */
/*************************************************/
void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]);

void hp_mgrid::setinflow() {
	static int i,j,k,m,n,indx,v0,v1,info;
	static FLT x,y,mvel[ND];
	static int sind;
   char uplo[] = "U";

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
/*			INFLOW BOUNDARIES */
/*			SET VERTEX VALUES OF U,V */	
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
   
            x = vrtx[v0][0];		
            y = vrtx[v0][1];
            ug.v[v0][0] = (*(gbl->func))(0,x,y);
            ug.v[v0][1] = (*(gbl->func))(1,x,y);
         }
         v0 = svrtx[sind][1];
         x = vrtx[v0][0];		
         y = vrtx[v0][1];
         ug.v[v0][0] = (*(gbl->func))(0,x,y);
         ug.v[v0][1] = (*(gbl->func))(1,x,y);
         
/**********************************/   
/*			SET SIDE VALUES & FLUXES */
/**********************************/
/*			ZERO FLUX FOR FIRST VERTEX */
         for(n=0;n<NV;++n)
            binfo[i][0].flx[n] = 0.0;
            
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            
            if (sbdry[i].type&CURV_MASK) {
               crdtouht1d(sind);
               for(n=0;n<ND;++n)
                  b.proj1d(uht[n],crd[n][0],dcrd[n][0][0]);
               
               crdtouht1d(sind,dvrtdt,gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  b.proj1d(uht[n],crd[n][1]);
            }
            else {
               for(n=0;n<ND;++n) {
                  b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
                  
                  for(k=0;k<b.gpx;++k)
                     dcrd[n][0][0][k] = 0.5*(vrtx[v1][n]-vrtx[v0][n]);
               
                  b.proj1d(dvrtdt[v0][n],dvrtdt[v1][n],crd[n][1]);
               }
            }

            if (b.sm) {
               for(n=0;n<ND;++n)
                  b.proj1d(ug.v[v0][n],ug.v[v1][n],res[n][0]);
         
               for(k=0;k<b.gpx; ++k)
                  for(n=0;n<ND;++n)
                     res[n][0][k] -= (*(gbl->func))(n,crd[0][0][k],crd[1][0][k]);
                     
               for(n=0;n<ND;++n)
                  b.intgrt1d(res[n][0],lf[n]);
         
               indx = sind*sm0;
               for(n=0;n<ND;++n) {
                  PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
                  for(m=0;m<b.sm;++m) 
                     ug.s[indx+m][n] = -lf[n][2+m];
               }
            }
            
/*				NOW SET FLUXES */
            ugtouht1d(sind);
            for(n=0;n<ND;++n)
               b.proj1d(uht[n],u[n][0]);

            for(k=0;k<b.gpx;++k) {
               for(n=0;n<ND;++n)
                  mvel[n] = bd[0]*crd[n][0][k] +crd[n][1][k];
               
               res[2][0][k] = gbl->rho*((u[0][0][k] -mvel[0])*dcrd[1][0][0][k] -(u[1][0][k] -mvel[1])*dcrd[0][0][0][k]);
            }
            
            b.intgrt1d(res[2][0],lf[0]);
            
            indx = j*(b.sm +1);
            binfo[i][indx++].flx[2] += lf[0][0];
            for(m=0;m<b.sm;++m)
               binfo[i][indx++].flx[2] = lf[0][m+2];
            binfo[i][indx].flx[2] = lf[0][1];
         }
      }
	}
   
	return;
}

void hp_mgrid::addbflux(int mgrid) {
	static int i,j,k,n,indx,indx1;
	static int sind,v0,v1;
	static FLT gam, nrm[ND], wl[NV], wr[NV];
   static FLT mvel[ND] = {0.0, 0.0};
	
/***********************************/
/*	ADD SOURCE TERMS ON FINEST MESH */
/***********************************/
	if(!mgrid) {
      for(i=0;i<nsbd;++i) {
         if (sbdry[i].type&INFL_MASK) {
            indx = 0;
            for(j=0;j<sbdry[i].num;++j) {
               sind=sbdry[i].el[j];
               v0 = svrtx[sind][0];
               gbl->res.v[v0][2] += binfo[i][indx++].flx[2];
               indx1 = sind*b.sm;
               for(k=0;k<b.sm;++k)
                  gbl->res.s[indx1++][2] += binfo[i][indx++].flx[2];
            }
            v0 = svrtx[sind][1];
            gbl->res.v[v0][2] += binfo[i][indx].flx[2];
         }
         
         if (sbdry[i].type&OUTF_MASK) {
/*				ALLOWS FOR APPLIED STRESS ON BOUNDARY */
            indx = 0;
            for(j=0;j<sbdry[i].num;++j) {
               sind=sbdry[i].el[j];
               v0 = svrtx[sind][0];
               indx1 = sind*b.sm;
               for(n=0;n<ND;++n)
                  gbl->res.v[v0][n] += binfo[i][indx].flx[n];
               ++indx;
               for(k=0;k<b.sm;++k) {
                  for(n=0;n<ND;++n)
                     gbl->res.s[indx1][n] += binfo[i][indx].flx[n];
                  ++indx;
                  ++indx1;
               }
            }
            v0 = svrtx[sind][1];
            for(n=0;n<ND;++n)
               gbl->res.v[v0][n] += binfo[i][indx].flx[n];
         }
      }
   }
   
/*	THESE ARE SOURCE TERMS WHICH CHANGE WITH THE SOLUTION */
/* MUST BE UPDATED DURING MGRID FOR GOOD CONVERGENCE */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK)) {

/*			CALCULATE RESIDUAL / SURFACE TENSION TERMS */
         surfrsdl(i,mgrid);

/*			ADD SURFACE TENSION SOURCE TERM */
         indx = 0;
         for(j=0;j<sbdry[i].num;++j) {
            sind=sbdry[i].el[j];
            v0 = svrtx[sind][0];
            for(n=0;n<ND;++n)
               gbl->res.v[v0][n] += binfo[i][indx].flx[n];
            ++indx;
            indx1 = sind*b.sm;
            for(k=0;k<b.sm;++k) {
               for(n=0;n<ND;++n)
                  gbl->res.s[indx1][n] += binfo[i][indx].flx[n];
               ++indx;
               ++indx1;
            }
         }
         v0 = svrtx[sind][1];
         for(n=0;n<ND;++n)
            gbl->res.v[v0][n] += binfo[i][indx].flx[n];
      }


/* 	OUTFLOW BOUNDARY CONDITION    */
      if (sbdry[i].type&OUTF_MASK) {
         indx = 0;
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            
            if (sbdry[i].type&CURV_MASK) {
               crdtouht1d(sind);
               for(n=0;n<ND;++n)
                  b.proj1d(uht[n],crd[n][0],dcrd[n][0][0]);
               
               crdtouht1d(sind,dvrtdt,gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  b.proj1d(uht[n],crd[n][1]);
            }
            else {
               for(n=0;n<ND;++n) {
                  b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
                  
                  for(k=0;k<b.gpx;++k)
                     dcrd[n][0][0][k] = 0.5*(vrtx[v1][n]-vrtx[v0][n]);
               
                  b.proj1d(dvrtdt[v0][n],dvrtdt[v1][n],crd[n][1]);
               }
            }
            
            ugtouht1d(sind);
            for(n=0;n<NV;++n)
               b.proj1d(uht[n],u[n][0]);
            
            gam = 1.0/gbl->gam[stri[sind][0]];
            for(k=0;k<b.gpx;++k) {
               for(n=0;n<NV;++n) {
                  wl[n] = u[n][0][k];
                  wr[n] = (gbl->func)(n,crd[0][0][k],crd[1][0][k]);
               }
               nrm[0] = dcrd[1][0][0][k];
               nrm[1] = -dcrd[0][0][0][k];

               for(n=0;n<ND;++n)
                  mvel[n] = bd[0]*crd[n][0][k] +crd[n][1][k];
                  
               if (!charyes)
                  wl[2] = wr[2];
               else 
                  chrctr(gbl->rho,gam,wl,wr,nrm,mvel);
                  
               res[2][0][k] = gbl->rho*((wl[0] -mvel[0])*nrm[0] +(wl[1] -mvel[1])*nrm[1]);
               res[0][0][k] = res[2][0][k]*wl[0] +wl[2]*nrm[0];
               res[1][0][k] = res[2][0][k]*wl[1] +wl[2]*nrm[1];
            }
            
            for(n=0;n<NV;++n)
               b.intgrt1d(res[n][0],lf[n]);
            
            for(n=0;n<NV;++n)
               gbl->res.v[v0][n] += lf[n][0];

            for(n=0;n<NV;++n)
               gbl->res.v[v1][n] += lf[n][1];
            
            indx1 = sind*b.sm;
            indx = 2;
            for(k=0;k<b.sm;++k) {
               for(n=0;n<NV;++n)
                  gbl->res.s[indx1][n] += lf[n][indx];
               ++indx1;
               ++indx;
            }
         }
      }

/* 	INVISCID WALL BOUNDARY CONDITION */
      if (sbdry[i].type&EULR_MASK) {
         indx = 0;
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            
            if (sbdry[i].type&CURV_MASK) {
               crdtouht1d(sind);
               for(n=0;n<ND;++n)
                  b.proj1d(uht[n],crd[n][0],dcrd[n][0][0]);
               
               crdtouht1d(sind,dvrtdt,gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  b.proj1d(uht[n],crd[n][1]);
            }
            else {
               for(n=0;n<ND;++n) {
                  b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
                  
                  for(k=0;k<b.gpx;++k)
                     dcrd[n][0][0][k] = 0.5*(vrtx[v1][n]-vrtx[v0][n]);
               
                  b.proj1d(dvrtdt[v0][n],dvrtdt[v1][n],crd[n][1]);
               }
            }
            
            ugtouht1d(sind);
            b.proj1d(uht[2],u[2][0]);
               
            for(k=0;k<b.gpx;++k) {
               wl[0] = (gbl->func)(0,crd[0][0][k],crd[1][0][k]);
               wl[1] = (gbl->func)(1,crd[0][0][k],crd[1][0][k]);
               for(n=0;n<ND;++n)
                  mvel[n] = bd[0]*crd[n][0][k] +crd[n][1][k];
               res[2][0][k] = gbl->rho*((wl[0] -mvel[0])*dcrd[1][0][0][k] -(wl[1] -mvel[1])*dcrd[0][0][0][k]);
               res[0][0][k] = res[2][0][k]*wl[0] +u[2][0][k]*dcrd[1][0][0][k];
               res[1][0][k] = res[2][0][k]*wl[1] -u[2][0][k]*dcrd[0][0][0][k];
            }
            
            for(n=0;n<NV;++n)
               b.intgrt1d(res[n][0],lf[n]);
            
            for(n=0;n<NV;++n)
               gbl->res.v[v0][n] += lf[n][0];

            for(n=0;n<NV;++n)
               gbl->res.v[v1][n] += lf[n][1];
            
            indx1 = sind*b.sm;
            indx = 2;
            for(k=0;k<b.sm;++k) {
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
   
/*	SEND VERTEX INFO FOR Y_DIR*/
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
/*			SEND VERTEX INFO */
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
      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
/*			SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            for (n=0;n<ND;++n) 
               tgt->sbuff[bnum][count++] = gbl->res.v[v0][n];
         }
         v0 = svrtx[sind][1];
         for (n=0;n<ND;++n) 
            tgt->sbuff[bnum][count++] = gbl->res.v[v0][n];
      }
   }
   
   return;
}

void hp_mgrid::bdry_mp() {
   int i,j,n,sind,count,v0,bnum;
   class mesh *tgt;
   
/*	THIS PART IS TO RECEIVE AND ZERO FOR VERTICES */
/*	RECEIVE VRTX MESSAGES */
/* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         count = 0;
/*   		RECV VERTEX INFO */
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

      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
/*   		RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            for(n=0;n<ND;++n)
               gbl->res.v[v0][n] = 0.5*(gbl->res.v[v0][n] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         for(n=0;n<ND;++n)
            gbl->res.v[v0][n] = 0.5*(gbl->res.v[v0][n] +sbuff[i][count++]);
      }         
   }
   
/*	SEND VERTEX INFO FOR X_DIR*/
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
/*			SEND VERTEX INFO */
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
	static int i,j,n;
	static int sind,v0,count;
   
/*	THIS PART IS TO RECEIVE AND ZERO FOR VERTICES */
/*	RECEIVE VRTX MESSAGES */
/* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         count = 0;
/*   		RECV VERTEX INFO */
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

/*	APPLY VRTX DIRICHLET CONDITIONS TO RES */
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].type&INFL_MASK) {
         for(j=0;j<vbdry[i].num;++j) {
            v0 = vbdry[i].el[j];
            for(n=0;n<ND;++n)
               gbl->res.v[v0][n] = 0.0;
         }
      }
      
      if (vbdry[i].type&SYMM_MASK) {
         for(j=0;j<vbdry[i].num;++j) {
            v0 = vbdry[i].el[j];
            gbl->res.v[v0][0] = 0.0;
         }
      }
   }

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            for(n=0;n<ND;++n)
               gbl->res.v[v0][n] = 0.0;
         }
         v0 = svrtx[sind][1];
         for(n=0;n<ND;++n)
            gbl->res.v[v0][n] = 0.0;
      }
      
      if (sbdry[i].type&SYMM_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            gbl->res.v[v0][0] = 0.0;
         }
         v0 = svrtx[sind][1];
         gbl->res.v[v0][0] = 0.0;
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
            indx = sbdry[i].el[j]*b.sm +mode;
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = gbl->res.s[indx][n];
         }
      } 

      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         for(j=0;j<sbdry[i].num;++j) {
            indx = sbdry[i].el[j]*b.sm +mode;
            for (n=0;n<ND;++n) 
               tgt->sbuff[bnum][count++] = gbl->res.s[indx][n];
         }
      }         
   }
   
   return;
}

	
void hp_mgrid::bdry_srcvandzero(int mode) {
	static int i,j,n;
	static int sind,count,indx;
   
/*	THIS PART TO RECIEVE AND ZERO FOR SIDES */
/*	RECEIVE P'TH SIDE MODE MESSAGES */
/* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMX_MASK +COMY_MASK)) {
         count = 0;
         for(j=sbdry[i].num-1;j>=0;--j) {
            indx = sbdry[i].el[j]*b.sm +mode;
            for(n=0;n<NV;++n)
               gbl->res.s[indx][n] = 0.5*(gbl->res.s[indx][n] +sbuff[i][count++]);
         }
      }
      if (sbdry[i].type & IFCE_MASK) {
         count = 0;
         for(j=sbdry[i].num-1;j>=0;--j) {
            indx = sbdry[i].el[j]*b.sm +mode;
            for(n=0;n<ND;++n)
               gbl->res.s[indx][n] = 0.5*(gbl->res.s[indx][n] +sbuff[i][count++]);
         }
      }         
   }

/*	APPLY SIDE DIRICHLET CONDITIONS TO MODE */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j]*b.sm +mode;
            for(n=0;n<ND;++n)
               gbl->res.s[sind][n] = 0.0;
         }
      }
      
      if (sbdry[i].type&SYMM_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j]*b.sm +mode;
            gbl->res.s[sind][0] = 0.0;
         }
      }
   }
   
	return;
}


void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]) {
	static FLT ul,vl,ur,vr,pl,pr,rhoi;
	static FLT um,c,den,lam0,lam1,lam2,v[3],mag;
	
   rhoi = 1./rho;

/*	CHARACTERISTIC FAR-FIELD B.C. */	
	mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
	
	norm[0] /= mag;
	norm[1] /= mag;
	
	ul =  wl[0]*norm[0] +wl[1]*norm[1];
	vl = -wl[0]*norm[1] +wl[1]*norm[0];
	pl =  wl[2];
		
/*	DEPENDENT ON FREESTREAM CONDITIONS */
	ur =  wr[0]*norm[0] +wr[1]*norm[1];
	vr = -wr[0]*norm[1] +wr[1]*norm[0];
	pr =  wr[3];
		
	um = mv[0]*norm[0] +mv[1]*norm[1];
	
	c = sqrt((ul-.5*um)*(ul-.5*um) +gam);
	
	den = 1./(2*c);
	
	lam0 = ul-um;
	lam1 = ul-.5*um +c; /* always positive */
	lam2 = ul-.5*um -c; /* always negative */
	
/*	PERFORM CHARACTERISTIC SWAP */
/* BASED ON LINEARIZATION AROUND UL,VL,PL */
	v[0] = ((pl-pr)*rhoi +(ul*lam1 -ur*lam2))*den;

	if (lam0 > 0.0)
		v[1] = vl*((pl-pr)*rhoi +lam2*(ul-ur))*den/(ul-lam1) +vl;
	else
		v[1] = vl*((pl-pr)*rhoi +lam1*(ul-ur))*den/(ul-lam2) +vr;

	v[2] = (rho*(ul -ur)*gam - lam2*pl +lam1*pr)*den;

/*	CHANGE BACK TO X,Y COORDINATES */
	wl[0] =  v[0]*norm[0] -v[1]*norm[1];
	wl[1] =  v[0]*norm[1] +v[1]*norm[0];
	wl[2] =  v[2];

/*	SHOULDN'T CHANGE NORM */   
   norm[0] *= mag;
   norm[1] *= mag;
   
	return;
 
}

