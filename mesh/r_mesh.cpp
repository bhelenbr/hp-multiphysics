#include"r_mesh.h"
#include"utilities.h"
#include<iostream>
#include<cstdio>
#include<cmath>

FLT r_mesh::fadd, r_mesh::vnn;
int r_mesh::fixx_mask = (0xFF -PRDY_MASK);
int r_mesh::fixy_mask = (0xFF -PRDX_MASK -SYMM_MASK);
int r_mesh::fix2x_mask = (0xFF -FSRF_MASK -IFCE_MASK -PRDX_MASK -PRDY_MASK -EULR_MASK);
int r_mesh::fix2y_mask = (0xFF -FSRF_MASK -IFCE_MASK -PRDX_MASK -PRDY_MASK -EULR_MASK -SYMM_MASK);

void r_mesh::allocate(bool coarse, struct r_mesh_glbls *rginit) {

/*	global mgrid arrays */
	rg = rginit;
   if (!coarse) gbl_alloc(rginit);
   
/*	local storage */   
   ksprg = new FLT[maxvst];
   kvol = new FLT[maxvst];
   src = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   if (coarse) vrtx_frst = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   isfrst = false;
}

void r_mesh::gbl_alloc(struct r_mesh_glbls *store) {
   store->work = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   store->res = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   store->diag = new FLT[maxvst];
}

int r_mesh::setfine(class r_mesh& tgt) {
   fmpt = &tgt;
   if (fine == NULL)
      fine = new struct mg_trans[maxvst];
   return(mgconnect(fine,tgt));
}

int r_mesh::setcoarse(class r_mesh& tgt) {
   cmpt = &tgt;
   if (coarse == NULL)
      coarse = new struct mg_trans[maxvst];
   return(mgconnect(coarse,tgt));
}  
   

void r_mesh::rklaplace(void) {
   /* static */int sind,tind,v0,v1,k;
   /* static */FLT dx,dy,l;
   
   for(sind=0;sind<nside;++sind)
      ksprg[sind] = 0.0;

/* COEFFICIENTS FOR LAPLACE EQUATION */
/* THIS REQUIRES 2 EVALUATIONS OF SIDE LENGTH FOR EACH SIDE */
/* BUT IS LOGISTICALLY SIMPLE        */         
   for(tind=0;tind<ntri;++tind) {
      for(k=0;k<3;++k) {
         sind = tside[tind].side[k];
         v0 = svrtx[sind][0];
         v1 = svrtx[sind][1];
         dx = vrtx[v1][0] -vrtx[v0][0];
         dy = vrtx[v1][1] -vrtx[v0][1];      
         l  = (dx*dx +dy*dy)/area(tind);

         ksprg[sind] -= l;
         sind = tside[tind].side[(k+1)%3];
         ksprg[sind] += l;
         sind = tside[tind].side[(k+2)%3];
         ksprg[sind] += l;
      }
   }      

	calc_kvol();
   
   return;
}


void r_mesh::rksprg(void) {
   /* static */int sind,v0,v1;
   /* static */double dx,dy;

/* 2D SPRING CONSTANTS FINE MESH*/
   for(sind=0;sind<nside;++sind) {
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];
      dx = vrtx[v1][0] -vrtx[v0][0];
      dy = vrtx[v1][1] -vrtx[v0][1];
      ksprg[sind] = 1.0/sqrt(dx*dx +dy*dy);
   }

   calc_kvol();

   return;
}

void r_mesh::calc_kvol() {
   /* static */int i,tind;
   
   for(i=0;i<nvrtx;++i) 
      kvol[i] = 0.0;

   for(tind=0;tind<ntri;++tind) 
      for(i=0;i<3;++i) 
         kvol[tvrtx[tind][i]] += area(tind);
                  
/*	SEND COMMUNICATION PACKETS IN XFIRST */   
   send(XDIR_MP,(FLT *) kvol,0,0,1);
   
   return;
}

void r_mesh::kvol_mp() {
   rcv(XDIR_MP,(FLT *) kvol,0,0,1);
   send(YDIR_MP,(FLT *) kvol,0,0,1);
   return;
}

void r_mesh::kvoli() {
   /* static */int i;
   
   rcv(YDIR_MP,(FLT *) kvol,0,0,1);
   
   for(i=0;i<nvrtx;++i)
      kvol[i] = 1./kvol[i];
}

void r_mesh::rkmgrid(void) {
   /* static */int i,j,sind,tind,tind0,tind1,v0,v1;
   /* static */class r_mesh *fmesh;
   
   fmesh = static_cast<class r_mesh *>(fmpt);

/* TEMPORARILY USE DIAG TO STORE DIAGONAL SUM */
   for(i=0;i<fmesh->nvrtx;++i)
      rg->diag[i] = 0.0;

/* FORM KIJ SUM AT VERTICES */
   for(sind=0;sind<fmesh->nside;++sind) {
      v0 = fmesh->svrtx[sind][0];
      v1 = fmesh->svrtx[sind][1];
      rg->diag[v0] += fmesh->ksprg[sind];
      rg->diag[v1] += fmesh->ksprg[sind];
   }

   for(i=0;i<nside;++i)
      ksprg[i] = 0.0;

/* LOOP THROUGH FINE VERTICES   */
/* TO CALCULATE KSPRG ON COARSE MESH */   
   for(i=0;i<fmesh->nvrtx;++i) {
      tind = fmesh->coarse[i].tri;
      for(j=0;j<3;++j) {
         sind = tside[tind].side[j];
         ksprg[sind] -= fmesh->coarse[i].wt[j]*fmesh->coarse[i].wt[(j+1)%3]*rg->diag[i];
      }
   }

/* LOOP THROUGH FINE SIDES */
   for(i=0;i<fmesh->nside;++i) {
      v0 = fmesh->svrtx[i][0];
      v1 = fmesh->svrtx[i][1];
      tind0 = fmesh->coarse[v0].tri;
      tind1 = fmesh->coarse[v1].tri;
               
/*    TEMPORARILY STORE WEIGHTS FOR FINE POINTS (0,1) */
/*    FOR EACH COARSE VERTEX */
      for(j=0;j<3;++j)  {
         rg->work[tvrtx[tind1][j]][0] = 0.0;
         rg->work[tvrtx[tind0][j]][1] = 0.0;
      }

      for(j=0;j<3;++j)  {
         rg->work[tvrtx[tind0][j]][0] = fmesh->coarse[v0].wt[j];
         rg->work[tvrtx[tind1][j]][1] = fmesh->coarse[v1].wt[j];
      }
      
/*    LOOP THROUGH COARSE TRIANGLE 0 SIDES */
      for(j=0;j<3;++j) {
         sind = tside[tind0].side[j];
         ksprg[sind] += fmesh->ksprg[i]*
            (rg->work[svrtx[sind][0]][0]*rg->work[svrtx[sind][1]][1]
            +rg->work[svrtx[sind][1]][0]*rg->work[svrtx[sind][0]][1]);
      }

      if (tind0 != tind1) {
         for(j=0;j<3;++j) {
            sind = tside[tind1].side[j];
            if (stri[sind][0] +stri[sind][1] != tind0 +tind1) {
               ksprg[sind] += fmesh->ksprg[i]*
                  (rg->work[svrtx[sind][0]][0]*rg->work[svrtx[sind][1]][1]
                  +rg->work[svrtx[sind][1]][0]*rg->work[svrtx[sind][0]][1]);

            }
         }
      }
   }
   
/*	CALCULATE KVOL USING MULTIGRID I VOL I^T*/
/* LOOP THROUGH FINE VERTICES TO CALCULATE COARSE VOLUME  */
   for(i=0;i<nvrtx;++i)
      kvol[i] = 0.0;

   for(i=0;i<fmesh->nvrtx;++i) {
      tind = fmesh->coarse[i].tri;
      for(j=0;j<3;++j) {
         v0 = tvrtx[tind][j];
         kvol[v0] += fmesh->coarse[i].wt[j]/fmesh->kvol[i];
      }
   }

/*	SEND COMMUNICATION PACKETS IN XFIRST */   
   send(XDIR_MP,(FLT *) kvol,0,0,1);
   
   return;
}

void r_mesh::rkmgrid_mp() {
   rcv(XDIR_MP,(FLT *) kvol,0,0,1);
   send(YDIR_MP,(FLT *) kvol,0,0,1);
   return;
}

void r_mesh::rkmgridi() {
   /* static */int i;
   
   rcv(YDIR_MP,(FLT *) kvol,0,0,1);
   
   for(i=0;i<nvrtx;++i)
      kvol[i] = 1./kvol[i];
   
   return;
}   

/*************************************/ 
/* SECOND ORDER MESH MOVEMENT SCHEME */
/*************************************/
#ifndef FOURTH
void r_mesh::rsdl() {
   /* static */int i,j,n,sind,v0,v1;
   /* static */FLT dx,dy;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         rg->res[i][n] = 0.0;

   for(sind=0;sind<nside;++sind) {
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];

      dx = ksprg[sind]*(vrtx[v1][0]-vrtx[v0][0]);
      dy = ksprg[sind]*(vrtx[v1][1]-vrtx[v0][1]);

      rg->res[v0][0] -= dx;
      rg->res[v0][1] -= dy;

      rg->res[v1][0] += dx;
      rg->res[v1][1] += dy;
   }
   
/* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
   if (isfrst) {
      for(i=0;i<nvrtx;++i) 
         for(n=0;n<ND;++n)
            src[i][n] -= rg->res[i][n];

      isfrst = false;
   }

/* ADD IN MULTIGRID SOURCE OR FMESH SOURCE */
   for(i=0;i<nvrtx;++i) 
      for(n=0;n<ND;++n)
         rg->res[i][n] += src[i][n];
   
/* APPLY DIRICHLET BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & fixx_mask) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            rg->res[svrtx[sind][0]][0] = 0.0;
            rg->res[svrtx[sind][1]][0] = 0.0;
         }
      }
      if (sbdry[i].type & fixy_mask) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            rg->res[svrtx[sind][0]][1] = 0.0;
            rg->res[svrtx[sind][1]][1] = 0.0;
         }
      }
   }

/* SEND COMMUNICATION PACKETS */   
   send(XDIR_MP, (FLT *) rg->res, 0, 1, 2);
   
   return;
}

void r_mesh::vddt(void)
{
   /* static */int i,v0,v1,sind;

/**************************************************/
/*   DETERMINE MESH MOVEMENT TIME STEP           */
/**************************************************/
   for(i=0;i<nvrtx;++i)
      rg->diag[i] = 0.0;

/* FORM TIME STEP FOR MV_UPDATE */
   for(sind=0;sind<nside;++sind) {
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];
      rg->diag[v0] += fabs(ksprg[sind]);
      rg->diag[v1] += fabs(ksprg[sind]);
   }

/*	SEND COMMUNICATION PACKETS */   
   send(XDIR_MP, (FLT *) rg->diag,0,0,1);
   
   return;
}

#else
/*************************************/ 
/* FOURTH ORDER MESH MOVEMENT SCHEME */
/*	THIS HAS TO BE DONE IN TWO PARTS  */
/* TO CACLUATE LAPLACIAN 				 */
/*************************************/
void r_mesh::rsdl() {
   /* static */int i,j,n,v0,v1,sind;
   /* static */FLT dx,dy;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         rg->work[i][n] = 0.0;

   for (sind = 0; sind < nside; ++sind) {
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];

      dx = ksprg[sind]*(vrtx[v1][0]-vrtx[v0][0]);
      dy = ksprg[sind]*(vrtx[v1][1]-vrtx[v0][1]);

      rg->work[v0][0] -= dx;
      rg->work[v0][1] -= dy;

      rg->work[v1][0] += dx;
      rg->work[v1][1] += dy;      
   }

/* APPLY ZERO SECOND DERIVATIVE BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & fix2x_mask) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            rg->work[svrtx[sind][0]][0] = 0.0;
            rg->work[svrtx[sind][1]][0] = 0.0;
         }
      }
      if (sbdry[i].type & fix2y_mask) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            rg->work[svrtx[sind][0]][1] = 0.0;
            rg->work[svrtx[sind][1]][1] = 0.0;
         }
      }
   }
   
   send(XDIR_MP,(FLT *) rg->work,0,1,2);
   
   return;
}

void r_mesh::rsdl1_mp() {
   rcv(XDIR_MP,(FLT *) rg->work,0,1,2);
   send(YDIR_MP,(FLT *) rg->work,0,1,2);
   return;
}

void r_mesh::rsdl1() {
   /* static */int i,j,n,v0,v1,sind;
   /* static */FLT dx,dy;
   
   rcv(YDIR_MP,(FLT *) rg->work, 0,1,2);
   
/* DIVIDE BY VOLUME FOR AN APPROXIMATION TO D^2/DX^2 */
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n) 
         rg->work[i][n] *= kvol[i];
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         rg->res[i][n] = 0.0;

   for (sind = 0; sind < nside; ++sind) {
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];

      dx = ksprg[sind]*(rg->work[v1][0]-rg->work[v0][0]);
      dy = ksprg[sind]*(rg->work[v1][1]-rg->work[v0][1]);

      rg->res[v0][0] -= dx;
      rg->res[v0][1] -= dy;

      rg->res[v1][0] += dx;
      rg->res[v1][1] += dy;      
   }

/* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
   if (isfrst) {
      for(i=0;i<nvrtx;++i) 
         for(n=0;n<ND;++n)
            src[i][n] -= rg->res[i][n];

      isfrst = false;
   }
   
/* ADD IN MULTIGRID SOURCE OR FMESH SOURCE */
   for(i=0;i<nvrtx;++i) 
      for(n=0;n<ND;++n)
         rg->res[i][n] += src[i][n];  

/* APPLY BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & fixx_mask) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            rg->res[svrtx[sind][0]][0] = 0.0;
            rg->res[svrtx[sind][1]][0] = 0.0;
         }
      }
      if (sbdry[i].type & fixy_mask) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            rg->res[svrtx[sind][0]][1] = 0.0;
            rg->res[svrtx[sind][1]][1] = 0.0;
         }
      }
   }
   
   send(XDIR_MP, (FLT *) rg->res,0,1,2);

   return;
}

void r_mesh::vddt(void)
{
   /* static */int i,v0,v1,sind;

/**************************************************/
/*   DETERMINE MESH MOVEMENT TIME STEP           */
/**************************************************/
   for(i=0;i<nvrtx;++i)
      rg->diag[i] = 0.0;
      
   for(sind=0;sind<nside;++sind) {
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];
      rg->diag[v0] += fabs(ksprg[sind]);
      rg->diag[v1] += fabs(ksprg[sind]);
   }
   
   send(XDIR_MP,(FLT *) rg->diag,0,0,1);
   
   return;
}

void r_mesh::vddt1_mp() {
   rcv(XDIR_MP,(FLT *) rg->diag,0,0,1);
   send(YDIR_MP,(FLT *) rg->diag,0,0,1);
   return;
}

void r_mesh::vddt1(void) {
   /* static */int i,v0,v1,sind;

   rcv(YDIR_MP,(FLT *) rg->diag, 0,0,1);
   
   for(i=0;i<nvrtx;++i)
      rg->work[i][0] = rg->diag[i]*kvol[i];
      
   for(i=0;i<nvrtx;++i)
      rg->diag[i] = 0.0;
   
   for(sind=0;sind<nside;++sind) {
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];
      rg->diag[v0] += fabs(ksprg[sind])*(rg->work[v0][0] +fabs(ksprg[sind])*kvol[v1]);
      rg->diag[v1] += fabs(ksprg[sind])*(rg->work[v1][0] +fabs(ksprg[sind])*kvol[v0]);
   }
   
   send(XDIR_MP, (FLT *) rg->diag,0,0,1);

   return;
}
#endif

void r_mesh::rsdl_mp() {
   rcv(XDIR_MP, (FLT *) rg->res, 0, 1, 2);
   send(YDIR_MP, (FLT *) rg->res, 0, 1, 2);
   return;
}

void r_mesh::vddt_mp() {
   rcv(XDIR_MP, (FLT *) rg->diag,0,0,1);
   send(YDIR_MP, (FLT *) rg->diag,0,0,1);
   return;
}

void r_mesh::vddti(void) {
   /* static */int i;
   
   rcv(YDIR_MP,(FLT *) rg->diag, 0,0,1);
   
   for(i=0;i<nvrtx;++i)
      rg->diag[i] = vnn/rg->diag[i];
}

void r_mesh::update() {
   /* static */int i,n;

   rcv(YDIR_MP,(FLT *) rg->res, 0, 1, 2);
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         vrtx[i][n] -= rg->diag[i]*rg->res[i][n];
   
   return;
}

void r_mesh::source() {
   /* static */int i,n;
   
   for(i=0;i<nvrtx;++i) 
      for(n=0;n<ND;++n) 
         src[i][n] = 0.0;
         
   rsdl();
   
   return;
}

void r_mesh::sumsrc() {
   /* static */int i,n;
   
   rcv(YDIR_MP,(FLT *) rg->res, 0, 1, 2);

   for(i=0;i<nvrtx;++i) 
      for(n=0;n<ND;++n)
         src[i][n] = -1.0*rg->res[i][n];   
   return;
}

void r_mesh::mg_getfres() {
   /* static */int i,j,n,tind,v0;
   /* static */class r_mesh *fmesh;
   
   fmesh = static_cast<class r_mesh *>(fmpt);
      
   fmesh->rcv(YDIR_MP,(FLT *) rg->res,0,1,2); 
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         src[i][n] = 0.0;
         
/* LOOP THROUGH FINE VERTICES TO CALCULATE RESIDUAL  */
   for(i=0;i<fmesh->nvrtx;++i) {
      tind = fmesh->coarse[i].tri;
      for(j=0;j<3;++j) {
         v0 = tvrtx[tind][j];
         for(n=0;n<ND;++n)
            src[v0][n] += fadd*fmesh->coarse[i].wt[j]*rg->res[i][n];
      }
   }
   
/* LOOP THROUGH COARSE VERTICES   */
/* TO CALCULATE VRTX ON COARSE MESH */
   for(i=0;i<nvrtx;++i) {
      tind = fine[i].tri;

      for(n=0;n<ND;++n)
         vrtx[i][n] = 0.0;
         
      for(j=0;j<3;++j) {
         for(n=0;n<ND;++n)
            vrtx[i][n] += fine[i].wt[j]*fmesh->vrtx[fmesh->tvrtx[tind][j]][n];
      }
      
      for(n=0;n<ND;++n)
         vrtx_frst[i][n] = vrtx[i][n];
   }

   isfrst = true;

   return;
}

void r_mesh::mg_getcchng() {
   /* static */int i,j,n,ind,tind;
   /* static */class r_mesh *cmesh;
   
   cmesh = static_cast<class r_mesh *>(cmpt);

/* DETERMINE CORRECTIONS ON COARSE MESH   */   
   for(i=0;i<cmesh->nvrtx;++i)
      for(n=0;n<ND;++n) 
         cmesh->vrtx_frst[i][n] -= cmesh->vrtx[i][n];

/* LOOP THROUGH FINE VERTICES   */
/* TO DETERMINE CHANGE IN SOLUTION */   
   for(i=0;i<nvrtx;++i) {
      
      for(n=0;n<ND;++n)
         rg->res[i][n] = 0.0;
      
      tind = coarse[i].tri;
      
      for(j=0;j<3;++j) {
         ind = cmesh->tvrtx[tind][j];
         for(n=0;n<ND;++n) 
            rg->res[i][n] -= coarse[i].wt[j]*cmesh->vrtx_frst[ind][n];
      }
   }
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n) 
         vrtx[i][n] += rg->res[i][n];

   return;
}

void r_mesh::send(int mask, FLT *base,int bgn,int end, int stride) {
   /* static */int i,j,k,sind,count,offset,bnum;
   /* static */class mesh *tgt;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & mask) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            offset = svrtx[sind][0]*stride;
            for (k=bgn;k<=end;++k) 
               tgt->sbuff[bnum][count++] = base[offset+k];
         }
         offset = svrtx[sind][1]*stride;
         for (k=bgn;k<=end;++k) 
            tgt->sbuff[bnum][count++] = base[offset+k];
      }
   }
   
   return;
}

void r_mesh::rcv(int mask, FLT *base,int bgn,int end, int stride) {
   /* static */int i,j,k,sind,count,offset;
      
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & mask) {
         count = 0;
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            offset = svrtx[sind][1]*stride;
            for (k=bgn;k<=end;++k) 
               base[offset+k] = 0.5*(base[offset+k] +sbuff[i][count++]);
         }
         offset = svrtx[sind][0]*stride;
         for (k=bgn;k<=end;++k) 
               base[offset+k] = 0.5*(base[offset+k] +sbuff[i][count++]);
      }
   }
   
   return;
}

void r_mesh::maxres() {
   /* static */int i,n;
   /* static */FLT mxr[ND];

   for(n=0;n<ND;++n)
      mxr[n] = 0.0;

   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         mxr[n] = MAX(mxr[n],fabs(rg->res[i][n]));
         
   for(n=0;n<ND;++n)
      printf("%.3e  ",mxr[n]);
         
   return;
}

