#include"surface.h"
#include<myblas.h>

 
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
      srf->gbl.vres[0][n] = 0.0;

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
         srf->gbl.vres[indx][n] += lf[n][0];
         srf->gbl.vres[indx+1][n] = lf[n][1];
         for(m=0;m<b.sm;++m)
            srf->gbl.sres[b.sm*indx +m][n] = lf[n][m+2];
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
         for(i=0;i<sbdry[bnum].num+1;++i) {
            srf->vdres[log2p][i][0] = srf->gbl.fadd0*srf->gbl.vres0[i][0] -srf->gbl.vres[i][0];
            srf->vdres[log2p][i][1] = srf->gbl.fadd1*srf->gbl.vres0[i][1] -srf->gbl.vres[i][1];
         }
         for(i=0;i<sbdry[bnum].num*b.sm;++i) {
            srf->sdres[log2p][i][0] = srf->gbl.fadd0*srf->gbl.sres0[i][0] -srf->gbl.sres[i][0];
            srf->sdres[log2p][i][1] = srf->gbl.fadd1*srf->gbl.sres0[i][1] -srf->gbl.sres[i][1];
         }
      }
      for(i=0;i<sbdry[bnum].num+1;++i) {
         for(n=0;n<ND;++n)		
            srf->gbl.vres[i][n] += srf->vdres[log2p][i][n];
      }
      for(i=0;i<sbdry[bnum].num*b.sm;++i) {
         for(n=0;n<ND;++n)		
            srf->gbl.sres[i][n] += srf->sdres[log2p][i][n];
      }
   }
   else {
/*		ADD TANGENTIAL MESH MOVEMENT SOURCE */
      for(i=0;i<sbdry[bnum].num+1;++i) 
         srf->gbl.vres[i][0] += srf->vdres[log2p][i][0];

      for(i=0;i<sbdry[bnum].num*b.sm;++i) 
         srf->gbl.sres[i][0] += srf->sdres[log2p][i][0];
   }

   return;
}

void hp_mgrid::surfinvrt1(int bnum) {
	static int i,m,n,v0,v1,indx,indx1,end,vbnum;
   class surface *srf;
   class mesh *tgt;

/* POINTER TO STUFF NEEDED FOR SURFACES IS STORED IN MISCELLANEOUS */      
   srf = static_cast<class surface *>(sbdry[i].misc);
	
   if (srf->gbl.first == 0) return;
   
/*	INVERT MASS MATRIX */
/*	LOOP THROUGH SIDES */
	if (b.sm > 0) {
      indx1 = 0;
		for(indx = 0; indx<sbdry[bnum].num; ++indx) {
/*			SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */			
         for (m=0; m <b.sm; ++m) {
            for(n=0;n<ND;++n)
               srf->gbl.vres[indx][n] -= b.sfmv1d[0][m]*srf->gbl.sres[indx1][n];
            for(n=0;n<ND;++n)
               srf->gbl.vres[indx+1][n] -= b.sfmv1d[1][m]*srf->gbl.sres[indx1][n];
            ++indx1;
			}
      }
	}

/*	SEND COMMUNICATION INFO FOR ENDPOINTS */
   end = sbdry[bnum].num;
   v0 = svrtx[sbdry[bnum].el[0]][0];
   v1 = svrtx[sbdry[bnum].el[end-1]][1];
   if (v0 == v1) {
/*		SURFACE IS A LOOP ON SINGLE BLOCK */
      srf->gbl.vres[0][1] = 0.5*(srf->gbl.vres[0][1] +srf->gbl.vres[end][1]);
      srf->gbl.vres[end][1] = srf->gbl.vres[0][1];
      srf->gbl.vres[0][0] = 0.0;
      srf->gbl.vres[end][0] = 0.0;
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      
      if (vbdry[i].type&(HP_MGRID_MP)) {
         vbnum = vbdry[i].adjbnum;
         tgt = vbdry[i].adjmesh;
         tgt->vbuff[vbnum][0] = srf->gbl.vres[0][0];
         tgt->vbuff[vbnum][1] = srf->gbl.vres[0][1];
      }
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      
      if (vbdry[i].type&(HP_MGRID_MP)) {
         vbnum = vbdry[i].adjbnum;
         tgt = vbdry[i].adjmesh;
         tgt->vbuff[vbnum][0] = srf->gbl.vres[end][0];
         tgt->vbuff[vbnum][1] = srf->gbl.vres[end][1];
      }
   }

   return;
}

void hp_mgrid::surfinvrt2(int bnum) {
	static int i,m,n,v0,v1,indx,indx1,end,vbnum;
   FLT temp;
   class surface *srf;

/* POINTER TO STUFF NEEDED FOR SURFACES IS STORED IN MISCELLANEOUS */      
   srf = static_cast<class surface *>(sbdry[i].misc);
   if (srf->gbl.first == 0) return;
   
/* RECEIVE MESSAGES AND ZERO VERTEX MODES */
   end = sbdry[bnum].num;
   v0 = svrtx[sbdry[bnum].el[0]][0];
   v1 = svrtx[sbdry[bnum].el[end-1]][1];
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      
      if (vbdry[i].type&(COMX_MASK +COMY_MASK)) {
         srf->gbl.vres[0][0] = 0.5*(srf->gbl.vres[0][0] +vbuff[vbnum][0]);
         srf->gbl.vres[0][1] = 0.5*(srf->gbl.vres[0][1] +vbuff[vbnum][1]);
      }
      
      if (vbdry[i].type&(PRDX_MASK +PRDY_MASK +SYMM_MASK)) {
         srf->gbl.vres[0][0] = 0.0;
         srf->gbl.vres[0][1] = 0.5*(srf->gbl.vres[0][1] +vbuff[vbnum][1]);
      }
   }
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      
      if (vbdry[i].type&(COMX_MASK +COMY_MASK)) {
         srf->gbl.vres[end][0] = 0.5*(srf->gbl.vres[end][0] +vbuff[vbnum][0]);
         srf->gbl.vres[end][1] = 0.5*(srf->gbl.vres[end][1] +vbuff[vbnum][1]);
      }
      
      if (vbdry[i].type&(PRDX_MASK +PRDY_MASK +SYMM_MASK)) {
         srf->gbl.vres[end][0] = 0.0;
         srf->gbl.vres[end][1] = 0.5*(srf->gbl.vres[end][1] +vbuff[vbnum][1]);
      }
   }
   
/*	SOLVE FOR VERTEX MODES */
	for(i=0;i<=end;++i) {
		temp                = srf->gbl.vres[i][0]*srf->gbl.vdt[i][0][0] +srf->gbl.vres[i][1]*srf->gbl.vdt[i][0][1];
		srf->gbl.vres[i][1] = srf->gbl.vres[i][0]*srf->gbl.vdt[i][1][0] +srf->gbl.vres[i][1]*srf->gbl.vdt[i][1][1];
		srf->gbl.vres[i][2] = temp;
	}
   
/* HAVE TO CORRECT AT PRDC BOUNDARIES SO POINT DOESN'T MOVE OFF LINE */
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v0) continue;
      if (vbdry[i].type&(PRDX_MASK +SYMM_MASK))
         srf->gbl.vres[0][0] = 0.0;
      if (vbdry[i].type&PRDY_MASK)
         srf->gbl.vres[0][1] = 0.0;
   }       
 
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].el[0] != v1) continue;
      if (vbdry[i].type&(PRDX_MASK +SYMM_MASK))
         srf->gbl.vres[end][0] = 0.0;
      if (vbdry[i].type&PRDY_MASK)
         srf->gbl.vres[end][1] = 0.0;
   }   
	
/*	SOLVE FOR SIDE MODES */
	if (b.sm > 2) {
      indx1 = 0;
		for(indx = 0; indx<sbdry[bnum].num; ++indx) {
         
/*			INVERT SIDE MODES */
			DPBSLN(b.sdiag1d,b.sm,b.sbwth,&srf->gbl.sres[indx1][0],ND);		
			for(m=0;m<b.sm;++m) {
				temp 							= srf->gbl.sres[indx1][0]*srf->gbl.sdt[indx][0][0] +srf->gbl.sres[indx1][1]*srf->gbl.sdt[indx][0][1];
				srf->gbl.sres[indx1][1] = srf->gbl.sres[indx1][0]*srf->gbl.sdt[indx][1][0] +srf->gbl.sres[indx1][1]*srf->gbl.sdt[indx][1][1]; 		
				srf->gbl.sres[indx1][0] = temp;
            ++indx1;
         }
			
         indx1 -= b.sm;
         for(m=0;m<b.sm;++m) {
            for(n=0;n<ND;++n) {
               gbl.sres[indx1][n] -= b.vfms1d[0][m]*srf->gbl.vres[indx][n];
               gbl.sres[indx1][n] -= b.vfms1d[1][m]*srf->gbl.vres[indx+1][n];
            }
            ++indx1;
         }
		}
	}
   
   return;
}

