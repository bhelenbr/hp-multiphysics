/*
 *  adapt.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 23 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"hp_mgrid.h"
#include<myblas.h>
#include"utilities.h"
#include<assert.h>

/* THIS IS USED IN THE MVPTTOBDRY FUNCTION */
extern class spectral_hp *tgt;

void hp_mgrid::adapt(class hp_mgrid& str, FLT tolerance) {
   int i,j,m,n,v0,v1,sind,stgt,info,indx,indx1,indx2,nvrt0,tind,touchd,snum,step;
   FLT r,s,x,y,psi,upt[NV];
   char uplo[] = "U";
   
/*	COPY SOLUTION & MESH TO BEGIN */
   str.spectral_hp::copy(*this);
   
/*	COPY UNSTEADY SOURCE TERMS TOO */
   for(step=0;step<nstep-1;++step) {
      for(i=0;i<nvrtx;++i)
         for(n=0;n<NV;++n)
            ugstr[step].v[i][n] = gbl->ugbd[step].v[i][n];
      
      for(i=0;i<nside*b.sm;++i)
         for(n=0;n<NV;++n)
            ugstr[step].s[i][n] = gbl->ugbd[step].s[i][n];

      for(i=0;i<ntri*b.im;++i)
         for(n=0;n<NV;++n)
            ugstr[step].i[i][n] = gbl->ugbd[step].i[i][n];            
               
      for(i=0;i<nvrtx;++i)
         for(n=0;n<ND;++n)
            vrtxstr[step][i][n] = gbl->vrtxbd[step][i][n];
            
      for(i=0;i<nsbd;++i)
         if (sbdry[i].type&CURV_MASK) 
            for (j=0;j<sbdry[i].num;++j)
               binfostr[step][i][j] = gbl->binfobd[step][i][j];
   }

/*	SET TARGET POINTER: USED IN EXTERNAL FUNCTION MVPTTOBDRY PROVIDED TO MESH */
   tgt = &str;

/*	REDUCE maxsrch (KEEP FAILED TRIANGLE SEARCHES LOCAL) */
   mesh::maxsrch = 15;
   
/* BEGIN ADAPTION PROCEDURE */
   swap();

/*	CALCULATE SIDE LENGTH RATIO */
   for(i=0;i<nside;++i)
      fltwk[i] = (vlngth[svrtx[i][0]] +vlngth[svrtx[i][1]])/
         (2.*distance(svrtx[i][0],svrtx[i][1]));

/* COARSEN */ 
   nvrt0 = nvrtx;
   yaber(1.0/tolerance);
      
/*	MOVE KEPT VERTEX VALUES TO NEW POSITIONS */
   for(i=nvrt0-1;i>=nvrtx;--i) {
      if (vinfo[i] < 0) continue;
      v0 = vinfo[i];
      for(n=0;n<NV;++n)
         ug.v[v0][n] = ug.v[i][n];
         
      for(step=0;step<nstep-1;++step)
         for(n=0;n<NV;++n)
            gbl->ugbd[step].v[v0][n] = gbl->ugbd[step].v[i][n];
   }
 
/*	REFINE */
   nvrt0 = nvrtx;
   rebay(tolerance);
   
/*	TO SEE MESH MANIPULATION */
   for(i=0;i<nside;++i)
      sinfo[i] += 2;
   out_mesh("adapt");
   for(i=0;i<nside;++i)
      sinfo[i] -= 2;

/* MARK BOUNDARY VERTICES */
   for(i=nvrt0;i<nvrtx;++i)
      vinfo[i] = -1;

   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i].num;++j)
         vinfo[svrtx[sbdry[i].el[j]][0]] = 0;
      
/*	ASSIGN NEW VALUES */
   for(i=nvrt0;i<nvrtx;++i) {
      if (vinfo[i] < 0) {
         tind = str.findinteriorpt(vrtx[i][0],vrtx[i][1],r,s);
         str.ugtouht(tind);  
         str.b.ptprobe(NV,uht,ug.v[i],r,s);
         
         for(step=0;step<nstep-1;++step) {
            str.ugtouht(tind,ugstr[step]);
            str.b.ptprobe(NV,uht,gbl->ugbd[step].v[i]);
         }
      }
   }
         

/* ASSIGN NEW BOUNDARY VERTEX VALUES */
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i].num;++j) {
         v0 = svrtx[sbdry[i].el[j]][0];
         if (v0 >= nvrt0) {
            sind = findbdrypt(sbdry[i].type,vrtx[v0][0],vrtx[v0][1],psi);
            str.ugtouht1d(sind);  
            str.b.ptprobe1d(NV,uht,ug.v[v0],psi);
            
            for(step=0;step<nstep-1;++step) {
               str.ugtouht1d(tind,ugstr[step]);
               str.b.ptprobe1d(NV,uht,gbl->ugbd[step].v[i]);
            }
         }
      }
   }
   
   if (b.sm <= 0) {
      setbcinfo();
      return;
   }
         
      
/* ASSIGN NEW SIDE VALUES */
   indx = 0;
   indx1 = 0;
   for(sind=0;sind<nside;++sind) {
      switch (sinfo[sind]) {
         case(-1): 
/*				UNTOUCHED/UNMOVED */            
            break;

         case(-2):
/*				TOUCHED */
            if (stri[sind][1] < 0) break; // DO BOUNDARY SIDES SEPARATELY

            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];

            for(n=0;n<ND;++n)
               b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
 
            for(n=0;n<NV;++n)
               b.proj1d(ug.v[v0][n],ug.v[v1][n],res[n][0]);
               
            for(step=0;step<nstep-1;++step)
               for(n=0;n<NV;++n)
                  b.proj1d(ugstr[step].v[v0][n],ugstr[step].v[v1][n],bdwk[step][n][0]);
      
            for(i=0;i<b.gpx;++i) {
               tind = str.findinteriorpt(crd[0][0][i],crd[1][0][i],r,s);
               str.ugtouht(tind);  
               str.b.ptprobe(NV,uht,upt,r,s);
               for(n=0;n<NV;++n)
                  res[n][0][i] -= upt[n];
               
               for(step=0;step<nstep-1;++step) {
                  str.ugtouht(tind,ugstr[step]);
                  str.b.ptprobe(NV,uht,upt);
                  for(n=0;n<NV;++n)   
                     bdwk[step][n][0][i] -= upt[n];
               }
            }            
                  
            for(n=0;n<NV;++n)
               b.intgrt1d(res[n][0],lf[n]);
         
            for(n=0;n<NV;++n) {
               PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
               for(m=0;m<b.sm;++m) 
                  ug.s[indx+m][n] = -lf[n][2+m];
            }
            
            for(step=0;step<nstep-1;++step) {
               for(n=0;n<NV;++n)
                  b.intgrt1d(bdwk[step][n][0],lf[n]);
         
               for(n=0;n<NV;++n) {
                  PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
                  for(m=0;m<b.sm;++m) 
                     gbl->ugbd[step].s[indx+m][n] = -lf[n][2+m];
               }
            }
            
            break;

         default:
/*				SIDE INTACT BUT MOVED FROM IT'S ORIGINAL LOCATION */
            assert(sinfo[sind] > -1);
            
            indx2 = sinfo[sind]*str.sm0;
            for(m=0;m<b.sm;++m)
               for(n=0;n<NV;++n)
                  ug.s[indx+m][n] = str.ug.s[indx2+m][n];
                  
            for(step=0;step<nstep-1;++step)
               for(m=0;m<b.sm;++m)
                  for(n=0;n<NV;++n)
                     gbl->ugbd[step].s[indx+m][n] = gbl->ugbd[step].s[indx+m][n];
                  
            break;
      }     
      indx1 += str.sm0;
      indx += sm0;
   }
   
/*	UPDATE BOUNDARY SIDES */
   for(i=0;i<nsbd;++i) {
      indx = 0;
      for(j=0;j<sbdry[i].num;++j) {
         sind = sbdry[i].el[j];
         
         switch(sinfo[sind]) {
            case(-1): // UNTOUCHED
               indx1 = ((-str.stri[sind][1])%maxsbel)*str.sm0;
               for(m=0;m<b.sm;++m)
                  binfo[i][indx+m] = str.binfo[i][indx1+m];
               
               for(step=0;step<nstep-1;++step)
                  for(m=0;m<b.sm;++m)
                     gbl->binfobd[step][i][indx+m] = binfostr[step][i][indx1+m];   
               
               break;
               
            case(-2): // TOUCHED
               v0 = svrtx[sind][0];
               v1 = svrtx[sind][1];
   
               for(n=0;n<ND;++n)
                  b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
   
               for(n=0;n<NV;++n)
                  b.proj1d(ug.v[v0][n],ug.v[v1][n],res[n][0]);
                  
               for(step=0;step<nstep-1;++step)
                  for(n=0;n<NV;++n)
                     b.proj1d(ugstr[step].v[v0][n],ugstr[step].v[v1][n],bdwk[step][n][0]);
                     
               for(m=0;m<b.gpx;++m) {
                  x = crd[0][0][m];
                  y = crd[1][0][m];

/*						THIS MOVES X,Y PT TO BOUNDRY */            
                  stgt = str.findbdrypt(sbdry[i].type,x,y,psi);
                  crd[0][0][m] -= x;
                  crd[1][0][m] -= y;

/*						CALCULATE VALUE OF SOLUTION AT POINT */
                  str.ugtouht1d(stgt);
                  str.b.ptprobe1d(NV,uht,upt,psi);
                  for(n=0;n<NV;++n)
                  	res[n][0][m] -= upt[n]; 
                  
                  for(step=0;step<nstep-1;++step) {
                     str.ugtouht1d(stgt,ugstr[step]);
                     str.b.ptprobe(NV,uht,upt);
                     for(n=0;n<NV;++n)   
                        bdwk[step][n][0][i] -= upt[n];
                  }   
               }     
               
               if (sbdry[i].type&CURV_MASK) {
                  for(n=0;n<ND;++n) {
                     b.intgrt1d(crd[n][0],lf[n]);
                     PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
                  
                     for(m=0;m<b.sm;++m)
                        binfo[i][indx+m].curv[n] = -lf[n][m+2];
                  }
               }
               else {
                  for(n=0;n<ND;++n)
                     for(m=0;m<b.sm;++m)
                        binfo[i][indx+m].curv[n] = 0.0;
               }
               
               indx1 = sind*sm0;      
               for(n=0;n<NV;++n)
                  b.intgrt1d(res[n][0],lf[n]);
            
               for(n=0;n<NV;++n) {
                  PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
                  for(m=0;m<b.sm;++m) 
                     ug.s[indx1+m][n] = -lf[n][2+m];
               }
               
               for(step=0;step<nstep-1;++step) {
                  for(n=0;n<NV;++n)
                     b.intgrt1d(bdwk[step][n][0],lf[n]);
            
                  for(n=0;n<NV;++n) {
                     PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
                     for(m=0;m<b.sm;++m) 
                        gbl->ugbd[step].s[indx1+m][n] = -lf[n][2+m];
                  }
               }
               break;
               
            default:
/*					SIDE INTACT BUT MOVED FROM IT'S ORIGINAL LOCATION */
/*					ALREADY MOVED UG INFO JUST NEED TO MOVE BINFO */
               assert(sinfo[sind] > -1);
               indx1 = ((-str.stri[sinfo[sind]][1])%maxsbel)*str.sm0;
               for(m=0;m<b.sm;++m)
                  binfo[i][indx+m] = str.binfo[i][indx1+m];
                  
               for(step=0;step<nstep-1;++step) 
                  for(m=0;m<b.sm;++m)
                     gbl->binfobd[step][i][indx+m] = binfostr[step][i][indx1+m];
                     
               break;
         }
         indx += sm0;
      }
   }
   
   if (b.im <= 0) {
      setbcinfo();
      return;
   }

/*	FIGURE OUT WHICH TRIANGLES HAVE BEEN TOUCHED */
   for(tind=0;tind<ntri;++tind) {
      touchd = 0;
      for(snum=0;snum<3;++snum)
         touchd = MIN(touchd,sinfo[tside[tind].side[snum]]);
         
      if (touchd == -2) intwk2[tind] = -2;
      else intwk2[i] = tinfo[tind];
   }
   
/* SET UP BDRY CONDITION INFO */
   setbcinfo();
      
/*	RESET INTERIOR VALUES */
   indx = 0;
   for(tind=0;tind<ntri;++tind) {
      switch (intwk2[tind]) {
         case(0): // UNTOUCHED & UNMOVED
            break;
            
         case(-2): // TOUCHED
            ugtouht_bdry(tind);
            for(n=0;n<NV;++n)
               b.proj_bdry(uht[n],u[n]);
               
            for(step=0;step<nstep-1;++step) {
               ugtouht_bdry(tind,ugstr[step]);
               for(n=0;n<NV;++n)
                  b.proj_bdry(uht[n],bdwk[step][n]);
            }
               
            if (tinfo[tind] < 0) {
               for(n=0;n<ND;++n)
                  b.proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);
            }
            else {
               crdtouht(tind);
               for(n=0;n<ND;++n)
                  b.proj_bdry(uht[n],crd[n]);
            }
               
            for (i=0; i < b.gpx; ++i ) {
               for (j=0; j < b.gpn; ++j ) {
                  tind = str.findinteriorpt(crd[0][i][j],crd[1][i][j],r,s);
                  str.ugtouht(tind);
                  str.b.ptprobe(NV,uht,upt,r,s);
                  for(n=0;n<NV;++n)
                     u[n][i][j] -= upt[n];
                  
                  for(step=0;step<nstep-1;++step) {
                     str.ugtouht(tind,ugstr[step]);
                     str.b.ptprobe(NV,uht,upt);
                     for(n=0;n<NV;++n)
                        bdwk[step][n][i][j] -= upt[n];
                  }
               }
            }
                           
            for(n=0;n<NV;++n) {
               b.intgrt(u[n],lf[n]);
               PBTRS(uplo,b.im,b.ibwth,1,b.idiag[0],b.ibwth+1,&lf[n][b.bm],b.im,info);
               for(i=0;i<b.im;++i)
                  ug.i[indx+i][n] = -lf[n][b.bm+i];
            }
            
            for(step=0;step<nstep-1;++step) {
               for(n=0;n<NV;++n) {
                  b.intgrt(bdwk[step][n],lf[n]);
                  PBTRS(uplo,b.im,b.ibwth,1,b.idiag[0],b.ibwth+1,&lf[n][b.bm],b.im,info);
                  for(i=0;i<b.im;++i)
                     gbl->ugbd[step].i[indx+i][n] = -lf[n][b.bm+i];
               }
            }
            
            break;

         default:  // MOVED BUT NOT TOUCHED
            
            indx1 = intwk2[tind]*str.im0;
            
            for(i=0;i<b.im;++i)
               for(n=0;n<NV;++n)
                  ug.i[indx+i][n] = str.ug.i[indx1+m][n];
            
            for(step=0;step<nstep-1;++step) {
               for(n=0;n<NV;++n)
                  for(i=0;i<b.im;++i)
                     gbl->ugbd[step].i[indx+i][n] = ugstr[step].i[indx+i][n];
            }
            
            break;
      }
      indx += im0;
   }
   
/*	RESET intwk2 */
   for(i=0;i<ntri;++i)
      intwk2[i] = -1;
      
/*	RESTORE maxsrch in findtri */
   mesh::maxsrch = 3*MAXLST/4;
            
   return;
}
   
   
   
   
   

   

