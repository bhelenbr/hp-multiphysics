#include"surface.h"

 
void hp_mgrid::surfrsdl(int bnum, int mgrid) {
	static int i,m,n,sind,indx,count;
	static FLT norm[ND], jcb, tau, tabs;
	static FLT dnormdt, hsm;
   FLT sigor, drhor;
   class surface *srf;

/*	DETRMINE DX CORRECTION TO CONSERVE AREA */
/*	IMPORTANT FOR STEADY SOLUTIONS */
/*	SINCE THERE ARE MULTIPLE STEADY-STATES */
   if (!sbdry[bnum].type&(IFCE_MASK +FSRF_MASK)) {
      printf("error shouldn't be in surfrsdl\n");
      exit(1);
   }

   srf = static_cast<class surface *>(sbdry[i].misc);
   
/* POINTER TO STUFF NEEDED FOR SURFACES IS STORED IN MISCELLANEOUS */      
   if (srf->gbl.first == 0) return;
   
   count = 0;
   sigor = gbl.sigma/srf->gbl.rhoav;
   drhor = srf->gbl.drho/srf->gbl.rhoav;
      
/*	CONSERVE AREA FOR CLOSED BDRY PROBLEMS */	
//   if (dt0 == 0.0 && lamv > 0.0) dnormdt = lamv*cnsrvarea(bnum);
//   else dnormdt = 0.0;

/**************************************************/
/*	DETERMINE MESH RESIDUALS & SURFACE TENSION	  */
/**************************************************/
   for(n=0;n<ND;++n)
      srf->gbl.res[0][n] = 0.0;

   for(indx=0;indx<sbdry[bnum].num;++indx) {
      sind = sbdry[bnum].el[indx];
      
      crdtouht1d(sind);
      for(n=0;n<ND;++n)
         b.proj1d(uht[n],crd[n][0],dcrd[n][0][0]);

      ugtouht1d(sind);
      for(n=0;n<ND;++n)
         b.proj1d(uht[n],u[n][0]);

      for(i=0;i<b.gpx;++i) {
         norm[0] =  dcrd[1][0][0][i];
         norm[1] = -dcrd[0][0][0][i];
         jcb = sqrt(norm[0]*norm[0] +norm[1]*norm[1]);
/*			FIGURE OUT WHAT TO DO HERE FOR RELATIVE VELOCITY STORED IN CRD[N][1]*/
         for(n=0;n<ND;++n)
            crd[n][1][i] = u[n][0][i] -(dt0*crd[n][0][i] +dnormdt*norm[n]/jcb); // PLUS DX/DT TERMS!!?? 
            
         hsm = jcb/(.25*(b.p+1)*(b.p+1));
         tau = (crd[0][1][i]*dcrd[0][0][0][i] +crd[1][1][i]*dcrd[1][0][0][i])/jcb;
         tabs = fabs(tau) + EPSILON;
         tau = tau/(jcb*(tabs*tabs/hsm +dt0 +(sigor/(hsm*hsm) +drhor*gbl.g*fabs(norm[1]/jcb))/tabs));

/*			TANGENTIAL SPACING & NORMAL FLUX */            
         res[0][0][i] = srf->ksprg[indx]*jcb;
         res[1][0][i] = crd[0][1][i]*norm[0] +crd[1][1][i]*norm[1];
         res[1][1][i] = res[1][0][i]*tau;
         
/*			SURFACE TENSION SOURCE TERM */
         u[0][0][i] = -gbl.sigma*norm[1]/jcb;
         u[0][1][i] = +srf->gbl.drho*gbl.g*crd[1][0][i]*norm[0];
         u[1][0][i] = +gbl.sigma*norm[0]/jcb;
         u[1][1][i] = +srf->gbl.drho*gbl.g*crd[1][0][i]*norm[0];            
      }
      
      for(m=0;m<b.sm+2;++m)
         lf[n][0] = 0.0;

/*		INTEGRATE & STORE SURFACE RESIDUALS */               
      b.intgrtx1d(res[0][0],lf[0]);
      b.intgrt1d(res[1][0],lf[1]);
      b.intgrtx1d(res[1][1],lf[1]);

/*		STORE IN RES */
      for(n=0;n<ND;++n) {
         srf->gbl.res[count][n] += lf[n][0];
         for(m=0;m<b.sm;++m)
            srf->gbl.res[count+m+1][n] = lf[n][m+2];
         srf->gbl.res[count+b.sm+1][n] = lf[n][1];
      }
      
/*		INTEGRATE & STORE SURFACE TENSION SOURCE TERM */
      b.intgrt1d(u[0][1],lf[0]);
      b.intgrtx1d(u[0][0],lf[0]);
      b.intgrt1d(u[1][1],lf[1]);
      b.intgrt1d(u[1][0],lf[1]);

/*		STORE IN BINFO.FLUX */
      for(n=0;n<ND;++n) {
         binfo[bnum][count].flx[n] += lf[n][0];
         for(m=0;m<b.sm;++m)
            binfo[bnum][count+m+1].flx[n] = lf[n][m+2];
         binfo[bnum][count+b.sm+1].flx[n] = lf[n][1];
      }
      count += b.sm +1;
   }

/************************************************/
/*	MODIFY SURFACE RESIDUALS ON COARSER MESHES	*/
/************************************************/	
   if(mgrid) {
      if (isfrst) {
         for(i=0;i<count;++i) {
            srf->dres[log2p][i][0] = srf->gbl.fadd0*srf->gbl.res0[i][0] -srf->gbl.res[i][0];
            srf->dres[log2p][i][1] = srf->gbl.fadd1*srf->gbl.res0[i][1] -srf->gbl.res[i][1];
         }
      }
      for(i=0;i<count;++i)
         for(n=0;n<ND;++n)		
            srf->gbl.res[i][n] += srf->dres[log2p][i][n];
   }
   else {
/*		ADD TANGENTIAL MESH MOVEMENT SOURCE */
      for(i=0;i<count;++i)
         srf->gbl.res[i][0] += srf->dres[log2p][i][0];
   }

   return;
}

#ifdef SKIP
void surfinvrt(void) {
	static int i,m,n,sind,v0,v1,indx,indx1,sgn;
	static double temp;
	
/*	INVERT MASS MATRIX */
/*	LOOP THROUGH SIDES */
	if (sm > 2) {
		for(sind = 0; sind<eifce; ++sind) {
			indx = eifce+1 +sind*(sm-2);
/*			SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */			
			for (i=0; i<2; ++i) {
				for (m=0; m <sm-2; ++m)
					for(n=0;n<ND+1;++n)
						surfgf[sind+i][n] -= sfmv1d[i][m]*surfgf[indx+m][n];
			}
		}
	}
	
	if (spdx1 || surfclsd) {
		surfgf[0][0] += surfgf[eifce][0];
		surfgf[0][1] += surfgf[eifce][1];
		surfgf[0][2] += surfgf[eifce][2];
		surfgf[eifce][0] = surfgf[0][0];
		surfgf[eifce][1] = surfgf[0][1];
		surfgf[eifce][2] = surfgf[0][2];
	}

/*	SOLVE FOR VERTEX MODES */
	for(i=0;i<eifce;++i) {
		v0 = svrtx[i][0];
		temp 			= surfgf[i][0]*vsurfdt[i][0][0]
						 +surfgf[i][1]*vsurfdt[i][0][1];
		mvgf[v0][1] = surfgf[i][0]*vsurfdt[i][1][0]
						 +surfgf[i][1]*vsurfdt[i][1][1];
		mvgf[v0][0] = temp;
		
		surfgf[i][0]   = surfgf[i][2]*vsurfdt[i][0][2];
		surfgf[i][1]   = surfgf[i][2]*vsurfdt[i][1][2];
		gf[v0][0] += surfgf[i][0];
		gf[v0][1] += surfgf[i][1];
	}
	v0 = svrtx[eifce-1][1];
	temp 			= surfgf[eifce][0]*vsurfdt[eifce][0][0]
					 +surfgf[eifce][1]*vsurfdt[eifce][0][1];
	mvgf[v0][1] = surfgf[eifce][0]*vsurfdt[eifce][1][0]
					 +surfgf[eifce][1]*vsurfdt[eifce][1][1];
	mvgf[v0][0] = temp;
	
	surfgf[eifce][0]   = surfgf[eifce][2]*vsurfdt[eifce][0][2];
	surfgf[eifce][1]   = surfgf[eifce][2]*vsurfdt[eifce][1][2];	
	gf[v0][0] += surfgf[eifce][0];
	gf[v0][1] += surfgf[eifce][1];
	
/*	SOLVE FOR SIDE MODES */
	if (sm > 2) {
		for(sind = 0; sind<eifce; ++sind) {
			indx = eifce+1 +sind*(sm-2);
			indx1 = nvrtx +sind*(sm-2);

/*			INVERT SIDE MODES */
			DPBSLN(sdiag1d,sm-2,SBWTH,&surfgf[indx][0],ND+1);		
			for(m=0;m<sm-2;++m) {
				temp = 				  surfgf[indx+m][0]*ssurfdt[sind][0][0]
				 						 +surfgf[indx+m][1]*ssurfdt[sind][0][1];
				mvgf[indx1+m][1] = surfgf[indx+m][0]*ssurfdt[sind][1][0]
				 						 +surfgf[indx+m][1]*ssurfdt[sind][1][1]; 		
				mvgf[indx1+m][0] = temp;
				
				gf[indx1+m][0]   += surfgf[indx+m][2]*ssurfdt[sind][0][2];
				gf[indx1+m][1]   += surfgf[indx+m][2]*ssurfdt[sind][1][2];							}
			
			for(i=0;i<2;++i) {
				v0 = svrtx[0+sind][i];
				for(m=0;m<sm-2;++m) {
					for(n=0;n<ND;++n) {
						mvgf[indx1+m][n] -= vfms1d[i][m]*mvgf[v0][n];
						gf[indx1+m][n] -= vfms1d[i][m]*surfgf[sind+i][n];
					}
				}
			}
		}
	}

	v0 = svrtx[0][0];
	v1 = svrtx[eifce-1][1];
			
	if (spdx1) {
		mvgf[v0][0] = 0.0;
		mvgf[v1][0] = 0.0;
	}
		
	if (spdy1)	{
		mvgf[v0][1] = 0.0;
		mvgf[v1][1] = 0.0;
	}
	
	if (ssymc) {
		mvgf[v0][0] = 0.0;
		mvgf[v1][0] = 0.0;
		gf[v0][0] = 0.0;
		gf[v1][0] = 0.0;
	}

#ifdef FIXEND
	mvgf[v0][0] = 0.0;
	mvgf[v1][0] = 0.0;
	mvgf[v1][1] = 0.0;
	gf[v1][0] = 0.0;
	gf[v1][1] = 0.0;
#endif


/*	EQUATE INTERFACE RESIDUALS */
	if (sifce) {
		for(i=0; i < vifce; ++i) {
			indx = vifce_lst[i];
			indx1 = vifc2_lst[vifce-1-i];
			gf[indx1][0] = gf[indx][0];
			gf[indx1][1] = gf[indx][1];
		}

		for(sind=0; sind < sifce; ++sind) {
			indx  = nvrtx + (sind+bifce)*(sm-2);
			indx1 = nvrtx + (eifc2-1-sind)*(sm-2);
			for(m=0;m<sm-2;++m) {
				sgn     = (m % 2 ? -1 : 1);		
				gf[indx1+m][0] = sgn*gf[indx+m][0];
				gf[indx1+m][1] = sgn*gf[indx+m][1];
			}		
		}
	}
	
	return;
}
#endif

