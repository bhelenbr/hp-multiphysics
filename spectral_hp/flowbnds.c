#include"hp_mgrid.h"
#include<myblas.h>

/*************************************************/
/*	SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/*	(THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION	 */
/*************************************************/

void hp_mgrid::setinflow(FLT (*func)(int, FLT, FLT)) {
	static int i,k,m,indx,indx1,indm,v0,v1,info;
	static double x,y;
	static int sind,sgn,ind;
	static double nrm[ND], flx[ND];

   for(i=0;i<nsbd;++i) {

      if (sbdry[i].type&INFL_MASK) {
/*			INFLOW BOUNDARIES */
/*			SET VERTEX VALUES OF U,V */	
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
   
            x = vrtx[v0][0];		
            y = vrtx[v0][1];
            ug[v0][0] = (*func)(0,x,y);
            ug[v0][1] = (*func)(1,x,y);
         }
         v0 = svrtx[sind][1];
         x = vrtx[v0][0];		
         y = vrtx[v0][1];
         ug[v0][0] = (*func)(0,x,y);
         ug[v0][1] = (*func)(1,x,y);
         
         if (!b.sm) continue;

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
            
            if (!sbdry[i].type&CURV_MASK) {
               for(n=0;n<ND;++n)
                  b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
                  for(k=0;k<b.gpx;++k)
                     dcrd[n][0][0][k] = 0.5*(vrtx[v1][n]-vrtx[v0][n]);
            }
            else {
               crdtouht1d(sind);
               for(n=0;n<ND;++n)
                  b.proj1d(uht[n],crd[n][0],dcrd[n][0][0]);
            }
            
            for(n=0;n<ND;++n)
               b.proj1d(vug[v0][n],vug[v1][n],res[n][0]);
      
            for(k=0;k<b.gpx; ++k)
               for(n=0;n<ND;++n)
                  res[n][0][k] -= (*func)(n,crd[0][0][k],crd[1][0][k]);
                  
            for(n=0;n<ND;++n)
               b.intgrt1d(res[n][0],lf[n]);
      
            indx = sind*sm0;
            for(n=0;n<ND;++n) {
               PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
               for(m=0;m<b.sm;++m) 
                  sug[indx+m][n] = -lf[n][2+m];
            }
            
/*				NOW SET FLUXES */
            ugtouht1d(sind);
            for(n=0;n<ND;++n)
               b.proj1d(uht[n],u[n][0]);

/*				NEED TO FIGURE OUT WHAT I'M DOING HERE FOR MOVING MESH */
               
            for(k=0;k<b.gpx;++k)
               res[2][0][k] = gbl.rho*(u[0][0][k]*dcrd[1][0][0][k] -u[1][0][k]*dcrd[0][0][0][k]);
            
            b.intgrt1d(res[2][0],lf[0]);
            
            indx = j*(b.sm +1);
            binfo[i][indx].flx[2] += lf[0][0];
            for(m=0;m<b.sm;++m)
               binfo[i][indx+m].flx[2] = lf[0][m+2];
            binfo[i][indx +b.sm+1].flx[2] = lf[0][1];
         }
      }
	}
   
	return;
}

void hp_mgrid::addbflux(int mgrid) {
	static int i,m,n,indx,indx1,indm, sgn;
	static int sind,indx1,v0,v1;
	static double fx[MXSM][NV], x[ND], nrm[ND], w[NV+ND];
	
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
               indx1 = sind*b.sm;
               vres[v0][2] += binfo[i][indx++].flx[2];
               for(k=0;k<b.sm;++k)
                  sres[indx1++][2] += binfo[i][indx++].flx[2];
            }
            v0 = svrtx[sind][1];
            vres[v0][2] += binfo[i].indx.flx[2];
         }
         
         if (sbdry[i].type&OUTF_MASK) {
/*				ALLOWS FOR APPLIED STRESS ON BOUNDARY */
            indx = 0;
            for(j=0;j<sbdry[i].num;++j) {
               sind=sbdry[i].el[j];
               v0 = svrtx[sind][0];
               indx1 = sind*b.sm;
               for(n=0;n<ND;++n)
                  vres[v0][n] += binfo[i][indx].flx[n];
               ++indx;
               for(k=0;k<b.sm;++k) {
                  for(n=0;n<ND;++n)
                     sres[indx1][n] += binfo[i][indx].flx[n];
                  ++indx;
                  ++indx1;
               }
            }
            v0 = svrtx[sind][1];
            for(n=0;n<ND;++n)
               vres[v0][n] += binfo[i].indx.flx[n];
         }
      }
   }
   
/*	THESE ARE SOURCE TERMS WHICH CHANGE WITH THE SOLUTION */
/* MUST BE UPDATED DURING MGRID FOR GOOD CONVERGENCE */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK)) {
/*			SURFACE TENSION SOURCE TERM */
         indx = 0;
         for(j=0;j<sbdry[i].num;++j) {
            sind=sbdry[i].el[j];
            v0 = svrtx[sind][0];
            indx1 = sind*b.sm;
            for(n=0;n<ND;++n)
               vres[v0][n] += binfo[i][indx].flx[n];
            ++indx;
            for(k=0;k<b.sm;++k) {
               for(n=0;n<ND;++n)
                  sres[indx1][n] += binfo[i][indx].flx[n];
               ++indx;
               ++indx1;
            }
         }
         v0 = svrtx[sind][1];
         for(n=0;n<ND;++n)
            vres[v0][n] += binfo[i].indx.flx[n];
      }


/* 	OUTFLOW BOUNDARY CONDITION    */
      if (sbdry[i].type&OUTF_MASK) {
         indx = 0;
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            
            if (!sbdry[i].type&CURV_MASK) {
               for(n=0;n<ND;++n)
                  b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
                  for(k=0;k<b.gpx;++k)
                     dcrd[n][0][0][k] = 0.5*(vrtx[v1][n]-vrtx[v0][n]);
            }
            else {
               crdtouht1d(sind);
               for(n=0;n<ND;++n)
                  b.proj1d(uht[n],crd[n][0],dcrd[n][0][0]);
            }
            
            ugtouht1d(sind);
            for(n=0;n<ND;++n)
               b.proj1d(uht[n],u[n][0]);

/*				NEED TO FIGURE OUT WHAT I'M DOING HERE FOR MOVING MESH */
               
            for(k=0;k<b.gpx;++k) {
               res[2][0][k] = gbl.rho*(u[0][0][k]*dcrd[1][0][0][k] -u[1][0][k]*dcrd[0][0][0][k]);
               res[0][0][k] = res[2][0][k]*u[0][0][k];
               res[1][0][k] = res[2][0][k]*u[1][0][k];
            }
            
            b.intgrt1d(res[2][0],lf[0]);
            
            indx = j*(b.sm +1);
            binfo[i][indx].flx[2] += lf[0][0];
            for(m=0;m<b.sm;++m)
               binfo[i][indx+m].flx[2] = lf[0][m+2];
            binfo[i][indx +b.sm+1].flx[2] = lf[0][1];
         }
         
         
		indx = sind*(mg_sm-2) + nvrtx;
		v0 = svrtx[sind][0];
		v1 = svrtx[sind][1];
		ind = tfluid[stri[sind][0]];
		nrm[0]  =  .5*(vrtx[v1][1] -vrtx[v0][1]);
		nrm[1]  = -.5*(vrtx[v1][0] -vrtx[v0][0]);

/*		CALCULATE NORMAL AND FLUX ON SIDE */
		for(i=0;i<sm;++i) {
			x[0] = vrtx[v0][0]*gx[0][i] +vrtx[v1][0]*gx[1][i];
			x[1] = vrtx[v0][1]*gx[0][i] +vrtx[v1][1]*gx[1][i];
			w[0] = ug[v0][0]*gx[0][i] +ug[v1][0]*gx[1][i];
			w[1] = ug[v0][1]*gx[0][i] +ug[v1][1]*gx[1][i];
			w[2] = ug[v0][2]*gx[0][i] +ug[v1][2]*gx[1][i];
			w[3] = mv[v0][0]*gx[0][i] +mv[v1][0]*gx[1][i];
			w[4] = mv[v0][1]*gx[0][i] +mv[v1][1]*gx[1][i];
			for(m=0;m<sm-2;++m) {
				w[0] += ug[indx+m][0]*gx[m+3][i];
				w[1] += ug[indx+m][1]*gx[m+3][i];
				w[2] += ug[indx+m][2]*gx[m+3][i];
			}
			chrctr(stri[sind][0],w,nrm,x);
			fx[i][2] = rho[ind]*((w[0]-w[3])*nrm[0] +(w[1]-w[4])*nrm[1]);
			fx[i][0] = fx[i][2]*w[0] +w[2]*nrm[0];
			fx[i][1] = fx[i][2]*w[1] +w[2]*nrm[1];
		}

/*		INTEGRATE FLUX */
		indx1 = sind*(sm-2) + nvrtx;
		for(n=0;n<NV;++n) {			
			for(m=0;m<sm;++m) {
				indm = m + (m > 1 ? 1 : 0);
				lf[n][m] = 0.0;
				for(i=0;i<sm;++i)
					lf[n][m] += wtx[i]*gx[indm][i]*fx[i][n];
			}

			gf[v0][n] += lf[n][0];
			gf[v1][n] += lf[n][1];
			for(m=0;m<sm-2;++m)
				gf[indx1+m][n] += lf[n][m+2];
		}
	}	
	
/* INVISCID WALL BOUNDARY CONDITION CURVED */
	for(sind=binvc; sind < einvc; ++sind) {
		indx = sind*(mg_sm-2) + nvrtx;
		v0 = svrtx[sind][0];
		v1 = svrtx[sind][1];
		ind = tfluid[stri[sind][0]];

/*		CALCULATE NORMAL AND FLUX ON SIDE */
		for(i=0;i<sm;++i) {
			nrm[0]  =  .5*(vrtx[v1][1] -vrtx[v0][1]);
			nrm[1]  = -.5*(vrtx[v1][0] -vrtx[v0][0]);
			x[0] = vrtx[v0][0]*gx[0][i] +vrtx[v1][0]*gx[1][i];
			x[1] = vrtx[v0][1]*gx[0][i] +vrtx[v1][1]*gx[1][i];
			w[2] = ug[v0][2]*gx[0][i] +ug[v1][2]*gx[1][i];
			w[3] = mv[v0][0]*gx[0][i] +mv[v1][0]*gx[1][i];
			w[4] = mv[v0][1]*gx[0][i] +mv[v1][1]*gx[1][i];
			for(m=0;m<sm-2;++m) {
				x[0] = vrtx[indx+m][0]*gx[m+3][i];
				x[1] = vrtx[indx+m][1]*gx[m+3][i];
				w[3] = mv[indx+m][0]*gx[m+3][i];
				w[4] = mv[indx+m][1]*gx[m+3][i];
				nrm[0] += vrtx[indx+m][1]*dgx[m+3][i];
				nrm[1] -= vrtx[indx+m][0]*dgx[m+3][i];
				w[2]   += ug[indx+m][2]*gx[m+3][i];
			}
			fx[i][2] = rho[ind]*((f1(x[0],x[1],ind)-w[3])*nrm[0]
									  +(f2(x[0],x[1],ind)-w[4])*nrm[1]);
			fx[i][0] = fx[i][2]*f1(x[0],x[1],ind) +w[2]*nrm[0];
			fx[i][1] = fx[i][2]*f2(x[0],x[1],ind) +w[2]*nrm[1];
		}

/*		INTEGRATE FLUX */
		indx1 = sind*(sm-2) + nvrtx;
		for(n=0;n<NV;++n) {			
			for(m=0;m<sm;++m) {
				indm = m + (m > 1 ? 1 : 0);
				lf[n][m] = 0.0;
				for(i=0;i<sm;++i)
					lf[n][m] += wtx[i]*gx[indm][i]*fx[i][n];
			}

			gf[v0][n] += lf[n][0];
			gf[v1][n] += lf[n][1];
			for(m=0;m<sm-2;++m)
				gf[indx1+m][n] += lf[n][m+2];
		}
	}
	
/* INVISCID WALL BOUNDARY CONDITION STRAIGHT */
	for(sind=binvs; sind < einvs; ++sind) {
		indx = sind*(mg_sm-2) + nvrtx;
		v0 = svrtx[sind][0];
		v1 = svrtx[sind][1];
		ind = tfluid[stri[sind][0]];
		nrm[0]  =  .5*(vrtx[v1][1] -vrtx[v0][1]);
		nrm[1]  = -.5*(vrtx[v1][0] -vrtx[v0][0]);

/*		CALCULATE NORMAL AND FLUX ON SIDE */
		for(i=0;i<sm;++i) {
			x[0] = vrtx[v0][0]*gx[0][i] +vrtx[v1][0]*gx[1][i];
			x[1] = vrtx[v0][1]*gx[0][i] +vrtx[v1][1]*gx[1][i];
			w[2] = ug[v0][2]*gx[0][i] +ug[v1][2]*gx[1][i];
			w[3] = mv[v0][0]*gx[0][i] +mv[v1][0]*gx[1][i];
			w[4] = mv[v0][1]*gx[0][i] +mv[v1][1]*gx[1][i];
			for(m=0;m<sm-2;++m) {
				w[2] += ug[indx+m][2]*gx[m+3][i];
			}
			fx[i][2] = rho[ind]*((f1(x[0],x[1],ind)-w[3])*nrm[0]
									  +(f2(x[0],x[1],ind)-w[4])*nrm[1]);
			fx[i][0] = fx[i][2]*f1(x[0],x[1],ind) +w[2]*nrm[0];
			fx[i][1] = fx[i][2]*f2(x[0],x[1],ind) +w[2]*nrm[1];
		}

/*		INTEGRATE FLUX */
		indx1 = sind*(sm-2) + nvrtx;
		for(n=0;n<NV;++n) {			
			for(m=0;m<sm;++m) {
				indm = m + (m > 1 ? 1 : 0);
				lf[n][m] = 0.0;
				for(i=0;i<sm;++i)
					lf[n][m] += wtx[i]*gx[indm][i]*fx[i][n];
			}

			gf[v0][n] += lf[n][0];
			gf[v1][n] += lf[n][1];
			for(m=0;m<sm-2;++m)
				gf[indx1+m][n] += lf[n][m+2];
		}
	}
	
	return;
}
	
void flowbnds(void) {
	static int i,m,n,ind,indx,indm, sgn;
	static int sind,indx1,v0,v1;

/***********************************/
/*	APPLY BOUNDARY CONDITIONS TO GF */
/***********************************/
/*	INFC BOUNDARY */
	for(i=0; i< vinfc; ++i) {
		for(n=0;n<2;++n)
			gf[vinfc_lst[i]][n] = 0.0;
	}

/*	INFS BOUNDARY */
	for(i=0; i< vinfs; ++i) {
		for(n=0;n<2;++n)
			gf[vinfs_lst[i]][n] = 0.0;
	}

/*	SYMMETRY BOUNDARY */
	for(i=0; i< vsymc; ++i)
		gf[vsymc_lst[i]][0] = 0.0;
	
/*	SUM INTERFACE RESIDUALS */
	for(i=0; i < vifce; ++i) {
		indx = vifce_lst[i];
		indx1 = vifc2_lst[vifce-1-i];
		gf[indx][0] += gf[indx1][0];
		gf[indx][1] += gf[indx1][1];
		gf[indx1][0] = 0;
		gf[indx1][1] = 0;
	}

/*	PERIODIC BOUNDARY  */
	for(i=0; i < vpdx1; ++i) {
		indx = vpdx1_lst[i];
		indx1 = vpdx2_lst[vpdx1-1-i];
		gf[indx][0] += gf[indx1][0];
		gf[indx][1] += gf[indx1][1];
		gf[indx][2] += gf[indx1][2];
		gf[indx1][0] = 0.0;
		gf[indx1][1] = 0.0;
		gf[indx1][2] = 0.0;		
	}

/*	PERIODIC BOUNDARY  */
	for(i=0; i < vpdy1; ++i) {
		indx = vpdy1_lst[i];
		indx1 = vpdy2_lst[vpdy1-1-i];
		gf[indx][0] += gf[indx1][0];
		gf[indx][1] += gf[indx1][1];
		gf[indx][2] += gf[indx1][2];
		gf[indx1][0] = 0.0;
		gf[indx1][1] = 0.0;
		gf[indx1][2] = 0.0;
	}
	
			
/*	SIDES */
	if (sm > 2) {
/*		INFC & INFS BOUNDARY */
		for(i = binfc; i < einfs; ++i) {
			indx = nvrtx+(sm-2)*i;
			for(m=0;m<sm-2;++m)
				for(n=0;n<2;++n) 
					gf[indx+m][n] = 0.0;
		}

/*		SYMMETRY BOUNDARY */		
		for(i = bsymc; i < esymc; ++i) {
			indx = nvrtx+(sm-2)*i;
			for(m=0;m<sm-2;++m)
				gf[indx+m][0] = 0.0;
		}

/*		INTERFACE BOUNDARY */
		for(i=0; i< sifce; ++i) {
			indx  = nvrtx + (i+bifce)*(mg_sm-2);
			indx1 = nvrtx + (eifc2-1-i)*(mg_sm-2);	
			sgn = 1;	
			for (m=0;m<sm-2;++m) {
				gf[indx +m][0] += sgn*gf[indx1 +m][0];		
				gf[indx +m][1] += sgn*gf[indx1 +m][1];
				gf[indx1 +m][0] = 0.0;
				gf[indx1 +m][1] = 0.0;
				sgn *= -1;
			}
		}	

/*		PERIODIC BOUNDARY X */
		for(i=0; i< spdx1; ++i) {
			indx  = nvrtx + (i+bpdx1)*(mg_sm-2);
			indx1 = nvrtx + (epdx2-1-i)*(mg_sm-2);
			sgn = 1;
			for (m=0;m<sm-2;++m) {
				gf[indx +m][0] += sgn*gf[indx1 +m][0];		
				gf[indx +m][1] += sgn*gf[indx1 +m][1];
				gf[indx +m][2] += sgn*gf[indx1 +m][2];
				gf[indx1 +m][0] = 0.0;
				gf[indx1 +m][1] = 0.0;
				gf[indx1 +m][2] = 0.0;
				sgn *= -1;
			}
		}
		
/*		PERIODIC BOUNDARY Y */
		for(i=0; i< spdy1; ++i) {
			indx  = nvrtx + (i+bpdy1)*(mg_sm-2);
			indx1 = nvrtx + (epdy2-1-i)*(mg_sm-2);
			sgn = 1;
			for (m=0;m<sm-2;++m) {
				gf[indx +m][0] += sgn*gf[indx1 +m][0];		
				gf[indx +m][1] += sgn*gf[indx1 +m][1];
				gf[indx +m][2] += sgn*gf[indx1 +m][2];
				gf[indx1 +m][0] = 0.0;
				gf[indx1 +m][1] = 0.0;
				gf[indx1 +m][2] = 0.0;
				sgn *= -1;
			}
		}
		
/*		SYMMETRIC BOUNDARY V = 0 */
		for(i=nvrtx+(mg_sm-2)*bsymc; i< nvrtx+(mg_sm-2)*esymc; ++i) 
			gf[i][0] = 0.0;
	}
	
	return;
}

void chrctr(FLT igam, double wl[3], double wr[3], double norm[2], double vnrm) {
	static int ind;
	static double ul,vl,ur,vr,pl,pr;
	static double um,c,den,lam0,lam1,lam2,v[3],mag;
	
	ind = tfluid[tind];

/*	PRESSURE FIXED B.C. */
	if (!charyes) {
		wl[2] = wr[2];
		return;
	}

/*	CHARACTERISTIC FAR-FIELD B.C. */	
	mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
	
	norm[0] /= mag;
	norm[1] /= mag;
	
	ul =  wl[0]*norm[0] +wl[1]*norm[1];
	vl = -wl[0]*norm[1] +wl[1]*norm[0];
	pl =  wl[2];
		
/*	DEPENDENT ON FREESTREAM CONDITIONS */
	ur =  f1(x[0],x[1],ind)*norm[0] +f2(x[0],x[1],ind)*norm[1];
	vr = -f1(x[0],x[1],ind)*norm[1] +f2(x[0],x[1],ind)*norm[0];
	pr =  f3(x[0],x[1],ind);
		
	um = wl[3]*norm[0] +wl[4]*norm[1];
	
	c = sqrt((ul-.5*um)*(ul-.5*um) + 1.0/gam[tind]);
	
	den = 1./(2*c);
	
	lam0 = ul-um;
	lam1 = ul-.5*um +c; /* always positive */
	lam2 = ul-.5*um -c; /* always negative */
	
/*	PERFORM CHARACTERISTIC SWAP */
/* BASED ON LINEARIZATION AROUND UL,VL,PL */
	v[0] = ((pl-pr)*rhoi[ind] +(ul*lam1 -ur*lam2))*den;

	if (lam0 > 0.0)
		v[1] = vl*((pl-pr)*rhoi[ind] +lam2*(ul-ur))*den/(ul-lam1) +vl;
	else
		v[1] = vl*((pl-pr)*rhoi[ind] +lam1*(ul-ur))*den/(ul-lam2) +vr;

	v[2] = (rho[ind]*(ul -ur)/gam[tind] - lam2*pl +lam1*pr)*den;

/*	CHANGE BACK TO X,Y COORDINATES */
	wl[0] =  v[0]*norm[0] -v[1]*norm[1];
	wl[1] =  v[0]*norm[1] +v[1]*norm[0];
	wl[2] =  v[2];
		
	norm[0] *= mag;
	norm[1] *= mag;
		
	return;
 
}

