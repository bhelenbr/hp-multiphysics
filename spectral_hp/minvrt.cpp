#include"hp_mgrid.h"


/* THIS REQUIRES A FAIR AMOUNT OF MESSAGE PASSING (NOT NEARLY AS MUCH AS AN ACTUAL INVERSION)*/
/* ANYWAY TO DO IT QUICKER? */

void hp_mgrid::minvrt1(void) {
	static int tind,sind,i,j,k,m,n,indx,indx1,indx2,*v,v0,v1,ind;
	static int sign[3],msgn,sgn,side[3];
	static double nrm[ND],nrmo[ND],olength;
	static double normco, meshco, nrmres, tanres;

	
/************************************************/
/**********		INVERT MASS MATRIX		**********/
/************************************************/
/*	LOOP THROUGH SIDES */
	if (b.sm > 0) {
      indx = 0;
		for(sind = 0; sind<nside;++sind) {
/*			SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */			
			for (i=0; i<2; ++i) {
				v0 = svrtx[sind][i];
				for (k=0; k <b.sm; ++k) {
					for(n=0;n<NV;++n)
						gbl.vres[v0][n] -= sfmv[i][k]*gbl.sres[indx][n];
               ++indx;
            }
			}
		}
			
		if (b.im > 0) {
/*			SUBTRACT INTERIORS */
         indx = 0;
			for(tind = 0; tind<ntri;++tind) {
            indx2 = 3;
				for (i=0; i<3; ++i) {
					v0 = tvrtx[tind][i];
					for (k=0;k<b.im;++k)
						for(n=0;n<NV;++n)
							gbl.vres[v0][n] -= ifmb[i][k]*gbl.ires[indx +k][n];

					indx1 = tside[tind].side[i]*sm;
					sgn = tside[tind].sign[i];
					msgn = 1;
					for (j=0;j<b.sm;++j) {
						for (k=0;k<b.im;++k) 
							for(n=0;n<NV;++n)
								gbl.sres[indx1][n] -= msgn*ifmb[indx2][k]*gbl.res[indx+k][n];
						msgn *= sgn;
                  ++indx1;
                  ++indx2;
					}
				}

/*				MULTIPLY INTERIOR MATRIX * LOCAL FORCING FOR INTERIOR MODES */		
				PBSLN(b.idiag,b.im,b.ibwth,&gbl.ires[indx][0],NV);
				for(i=0;i<im;++i) {
					gbl.res[indx][0] /= gbl.dtstar[tind];
					gbl.res[indx][1] /= gbl.dtstar[tind];
					gbl.res[indx][2] /= gbl.gam[tind]*gbl.dtstar[tind]*gbl.rhoi;
               ++indx;
				}
			}
		}
	}
   
   return;

/********************************/
/*	MESSAGE PASSING MUST GO HERE */	
/********************************/

void hp_mgrid::minvrt2(void)
/* SOLVE FOR VERTEX MODES */
	for(i=0;i<nvrtx;++i) {
		gbl.vres[i][0] *= gbl.vdiagv[i];
		gbl.vres[i][1] *= gbl.vdiagv[i];
		gbl.vres[i][2] *= gbl.vdiagp[i];
	}

/*	REMOVE VERTEX CONTRIBUTION FROM INTERIOR & SIDE MODES */
	if (b.sm > 0) {
/*		SOLVE FOR SIDE MODES */
/*		PART 1 REMOVE VERTEX CONTRIBUTIONS */
		for(tind=0;tind<ntri;++tind) {
			for(i=0;i<3;++i) {
            v0 = tvrtx[tind][i];
				uht[0][i] = gbl.vres[v0][0]*gbl.dtstar[tind];
				uht[1][i] = gbl.vres[v0][1]*gbl.dtstar[tind];
				uht[2][i] = gbl.vres[v0][2]*gbl.dtstar[tind]*gbl.gam[tind]*gbl.rhoi;
			}

			for(i=0;i<3;++i) {
				indx = tside[tind].side[i]*b.sm;
				sgn  = tside[tind].sign[i];
				for(j=0;j<3;++j) {
					indx1 = (i+j)%3;
					msgn = 1;
					for(k=0;k<b.sm;++k) {
						gbl.sres[indx][0] -= msgn*vfms[j][k]*uht[0][indx1];
						gbl.sres[indx][1] -= msgn*vfms[j][k]*uht[1][indx1];
						gbl.sres[indx][2] -= msgn*vfms[j][k]*uht[2][indx1];
						msgn *= sgn;
                  ++indx;
					}
				}
			}
		}
   }
   return;
}

/********************************/
/*	MESSAGE PASSING MUST GO HERE */	
/********************************/

void hp_mgrid::minvrt3(int mode) {  

/*	SOLVE FOR MODE */
   indx = mode;
   for(sind = 0; sind < nside; ++sind) {
      gbl.res[indx][0] *= sdiagv[sind]*sdiag[mode];
      gbl.res[indx][1] *= sdiagv[sind]*sdiag[mode];
      gbl.res[indx][2] *= sdiagp[sind]*sdiag[mode];
      indx += b.sm;
   }

/*	REMOVE MODE FROM HIGHER MODES */
   for(tind=0;tind<ntri;++tind) {
      for(i=0;i<3;++i) {
         side[i] = tside[tind].side[i]*b.sm;
         sign[i] = tside[tind].sign[i];
         sgn     = (mode % 2 ? sign[i] : 1);
         uht[0][i] = sgn*gbl.res[side[i]+mode][0]*gbl.dtstar[tind];
         uht[1][i] = sgn*gbl.res[side[i]+mode][1]*gbl.dtstar[tind];
         uht[2][i] = sgn*gbl.res[side[i]+mode][2]*gbl.dtstar[tind]*gbl.gam[tind]*gbl.rhoi;
      }
      
/*		REMOVE MODES J,K FROM MODE I,M */
      for(i=0;i<3;++i) {
         msgn = (mode +1 % 2 ? sign[i] : 1);
         for(m=mode+1;m<b.sm;++m) {
            for(j=0;j<3;++j) {
               indx = (i+j)%3;
               for(n=0;n<NV;++n) 
                  gbl.res[side[i]+m][n] -= msgn*sfms[mode][m][j]*uht[n][indx];
            }
            msgn *= sign[i];
         }
      }
   }
      
   return;
}

/********************************/
/*	MESSAGE PASSING MUST GO HERE */	
/********************************/

void hp_mgrid::minvrt4() {  
   int sind;
			
   indx = b.sm-1;
   for(sind = 0; sind < nside; ++sind) {
      gbl.res[indx][0] *= gbl.sdiagv[sind]*b.sdiag[sm-1];
      gbl.res[indx][1] *= gbl.sdiagv[sind]*b.sdiag[sm-1];
      gbl.res[indx][2] *= gbl.sdiagp[sind]*b.sdiag[sm-1];
      indx += b.sm;
   }					

/*		SOLVE FOR INTERIOR MODES */			
   if (im > 0) {
      for(tind = 0; tind < ntri; ++tind) {
         indx = gbl_e + tind*im;
         gtol_b(tind);
         for(k=0; k <im; ++k)
            for (i = 0; i<bm; ++i)
               for(n=0;n<NV;++n) 
                  gbl.res[indx +k][n] -= bfmi[i][k]*uht[n][i];
      }
   }
}

/**************************************/
/*	RE-APPLY BOUNDARY CONDITIONS TO GF */
/**************************************/	
/*	INFC BOUNDARY */
	for(i=0; i< vinfc; ++i) {
		for(n=0;n<2;++n)
			gbl.res[vinfc_lst[i]][n] = 0.0;
	}

/*	INFS BOUNDARY */
	for(i=0; i< vinfs; ++i) {
		for(n=0;n<2;++n)
			gbl.res[vinfs_lst[i]][n] = 0.0;
	}

/*	SYMMETRY BOUNDARY */
	for(i=0; i< vsymc; ++i)
		gbl.res[vsymc_lst[i]][0] = 0.0;
	
	if (sm > 2) {
/*		INFC & INFS BOUNDARY */
		for(i = binfc; i < einfs; ++i) {
			indx = nvrtx+(sm-2)*i;
			for(m=0;m<sm-2;++m)
				for(n=0;n<2;++n) 
					gbl.res[indx+m][n] = 0.0;
		}

/*		SYMMETRY BOUNDARY */		
		for(i = bsymc; i < esymc; ++i) {
			indx = nvrtx+(sm-2)*i;
			for(m=0;m<sm-2;++m)
				gbl.res[indx+m][0] = 0.0;
		}
	}

	return;
}
