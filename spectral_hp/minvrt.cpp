#include"hp_mgrid.h"
#include<myblas.h>

/* THIS REQUIRES MESSAGE PASSING FOR EACH P SIDE MODE */
void hp_mgrid::minvrt1(void) {
	static int tind,sind,i,j,k,n,indx,indx1,indx2,v0,sgn,msgn;
	
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
						gbl.vres[v0][n] -= b.sfmv[i][k]*gbl.sres[indx][n];
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
							gbl.vres[v0][n] -= b.ifmb[i][k]*gbl.ires[indx +k][n];

					indx1 = tside[tind].side[i]*b.sm;
					sgn = tside[tind].sign[i];
					msgn = 1;
					for (j=0;j<b.sm;++j) {
						for (k=0;k<b.im;++k) 
							for(n=0;n<NV;++n)
								gbl.sres[indx1][n] -= msgn*b.ifmb[indx2][k]*gbl.ires[indx+k][n];
						msgn *= sgn;
                  ++indx1;
                  ++indx2;
					}
				}

/*				MULTIPLY INTERIOR MATRIX * LOCAL FORCING FOR INTERIOR MODES */		
				DPBSLN(b.idiag,b.im,b.ibwth,&gbl.ires[indx][0],NV);
				for(i=0;i<b.im;++i) {
					gbl.ires[indx][0] /= gbl.dtstar[tind];
					gbl.ires[indx][1] /= gbl.dtstar[tind];
					gbl.ires[indx][2] /= gbl.gam[tind]*gbl.dtstar[tind]*gbl.rhoi;
               ++indx;
				}
			}
		}
	}
   
/*********************************************/
/*	SEND MESSAGES FOR VERTICES 					*/	
/*********************************************/

   return;
}

void hp_mgrid::minvrt2(void) {
   static int i,k,tind,v0,indx,j,indx1,sgn,msgn;
   
/**********************************/
/*  RECEIVE MESSAGES FOR VERTICES */
/* APPLY DIRCHLET B.C.S TO VERTICES */
/**********************************/
   bdry_rcvandzero(-1);

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
						gbl.sres[indx][0] -= msgn*b.vfms[j][k]*uht[0][indx1];
						gbl.sres[indx][1] -= msgn*b.vfms[j][k]*uht[1][indx1];
						gbl.sres[indx][2] -= msgn*b.vfms[j][k]*uht[2][indx1];
						msgn *= sgn;
                  ++indx;
					}
				}
			}
		}
   }
   
/* SEND MESSAGE FOR LOWEST ORDER MODE */
   
   return;
}

/********************************/
/*	MESSAGE PASSING MUST GO HERE */	
/********************************/

void hp_mgrid::minvrt3(int mode) {  
   static int i,j,m,n,indx,sind,tind;
   static int sign[3],msgn,sgn,side[3];
   
/* RECEIVE MESSAGE FOR MODE */
/* APPLY DIRCHLET B.C.S TO MODE */
   bdry_rcvandzero(mode);

/*	SOLVE FOR MODE */
   indx = mode;
   for(sind = 0; sind < nside; ++sind) {
      gbl.sres[indx][0] *= gbl.sdiagv[sind]*b.sdiag[mode];
      gbl.sres[indx][1] *= gbl.sdiagv[sind]*b.sdiag[mode];
      gbl.sres[indx][2] *= gbl.sdiagp[sind]*b.sdiag[mode];
      indx += b.sm;
   }

/*	REMOVE MODE FROM HIGHER MODES */
   for(tind=0;tind<ntri;++tind) {
      for(i=0;i<3;++i) {
         side[i] = tside[tind].side[i]*b.sm;
         sign[i] = tside[tind].sign[i];
         sgn     = (mode % 2 ? sign[i] : 1);
         uht[0][i] = sgn*gbl.sres[side[i]+mode][0]*gbl.dtstar[tind];
         uht[1][i] = sgn*gbl.sres[side[i]+mode][1]*gbl.dtstar[tind];
         uht[2][i] = sgn*gbl.sres[side[i]+mode][2]*gbl.dtstar[tind]*gbl.gam[tind]*gbl.rhoi;
      }
      
/*		REMOVE MODES J,K FROM MODE I,M */
      for(i=0;i<3;++i) {
         msgn = (mode +1 % 2 ? sign[i] : 1);
         for(m=mode+1;m<b.sm;++m) {
            for(j=0;j<3;++j) {
               indx = (i+j)%3;
               for(n=0;n<NV;++n) 
                  gbl.sres[side[i]+m][n] -= msgn*b.sfms[mode][m][j]*uht[n][indx];
            }
            msgn *= sign[i];
         }
      }
   }
   
/* SEND MESSAGES FOR NEXT MODE */
      
   return;
}

void hp_mgrid::minvrt4() {  
   int i,k,n,sind,indx,tind;
   
/* RECEIVE MESSAGE FOR LAST MODE */
/* APPLY DIRICHLET B.C.'S */
   indx = b.sm-1;
   bdry_rcvandzero(indx);
			
   for(sind = 0; sind < nside; ++sind) {
      gbl.sres[indx][0] *= gbl.sdiagv[sind]*b.sdiag[b.sm-1];
      gbl.sres[indx][1] *= gbl.sdiagv[sind]*b.sdiag[b.sm-1];
      gbl.sres[indx][2] *= gbl.sdiagp[sind]*b.sdiag[b.sm-1];
      indx += b.sm;
   }					

/*	SOLVE FOR INTERIOR MODES */			
   if (b.im > 0) {
      indx = 0;
      for(tind = 0; tind < ntri; ++tind) {
         ugtouht_bdry(tind);
         for(k=0;k<b.im;++k) {
            for (i=0;i<b.bm;++i)
               for(n=0;n<NV;++n) 
                  gbl.ires[indx][n] -= b.bfmi[i][k]*uht[n][i];
            ++indx;
         }
      }
   }
   
   return;
}
