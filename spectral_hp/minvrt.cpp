#include"hp_mgrid.h"
#include<myblas.h>

/* THIS REQUIRES MESSAGE PASSING FOR EACH P SIDE MODE */
void hp_mgrid::minvrt1(void) {
    int tind,sind,i,j,k,n,indx,indx1,indx2,v0,sgn,msgn;
         
   /************************************************/
   /**********      INVERT MASS MATRIX      **********/
   /************************************************/
   /* LOOP THROUGH SIDES */
   if (b.sm > 0) {
      indx = 0;
      for(sind = 0; sind<nside;++sind) {
         /* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */         
         for (k=0; k <b.sm; ++k) {
            for (i=0; i<2; ++i) {
               v0 = svrtx[sind][i];
               for(n=0;n<NV;++n)
                  gbl->res.v[v0][n] -= b.sfmv[i][k]*gbl->res.s[indx][n];
            }
            ++indx;
         }
      }
         
      if (b.im > 0) {
         /* SUBTRACT INTERIORS */
         indx = 0;
         for(tind = 0; tind<ntri;++tind) {
            indx2 = 3;
            for (i=0; i<3; ++i) {
               v0 = tvrtx[tind][i];
               for (k=0;k<b.im;++k)
                  for(n=0;n<NV;++n)
                     gbl->res.v[v0][n] -= b.ifmb[i][k]*gbl->res.i[indx +k][n];

               indx1 = tside[tind].side[i]*b.sm;
               sgn = tside[tind].sign[i];
               msgn = 1;
               for (j=0;j<b.sm;++j) {
                  for (k=0;k<b.im;++k)
                     for(n=0;n<NV;++n)
                        gbl->res.s[indx1][n] -= msgn*b.ifmb[indx2][k]*gbl->res.i[indx+k][n];
                  msgn *= sgn;
                  ++indx1;
                  ++indx2;
               }
            }
           indx += b.im;
         }
      }
   }
   
   /* SOLVE FOR VERTEX MODES */
#ifdef CONSERV
   for(i=0;i<nvrtx;++i) {
      gbl->res.v[i][0] *= gbl->vprcn[i][0][0];
      gbl->res.v[i][1] *= gbl->vprcn[i][0][0];
      gbl->res.v[i][2] *= gbl->vprcn[i][NV-1][NV-1];
   }
#else
   for(i=0;i<nvrtx;++i) {
      gbl->res.v[i][0] = gbl->res.v[i][0]*gbl->vprcn[i][0][0] +gbl->res.v[i][2]*gbl->vprcn[i][0][NV-1];
      gbl->res.v[i][1] = gbl->res.v[i][1]*gbl->vprcn[i][0][0] +gbl->res.v[i][2]*gbl->vprcn[i][1][NV-1];
      gbl->res.v[i][0] *= gbl->vprcn[i][NV-1][NV-1];
   }
#endif
      
   /* INVERT MATRICES FOR COUPLED BOUNDARY EQUATIONS */
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK+IFCE_MASK))
         surfinvrt1(i);
   
   /*********************************************/
   /* SEND MESSAGES FOR VERTICES                */   
   /*********************************************/
   bdry_vsnd();
   
   return;
}

   /* CALL bdry_mp() here */

void hp_mgrid::minvrt2(void) {
   int i,k,sind,tind,v0,indx,j,indx1,sgn,msgn;
   
   /**********************************/
   /*  RECEIVE MESSAGES FOR VERTICES */
   /* APPLY DIRCHLET B.C.S TO VERTICES */
   /**********************************/
   bdry_vrcvandzero();

   /* FINISH INVERSION FOR COUPLED BOUNDARY EQUATIONS */
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK+IFCE_MASK))
         surfinvrt2(i);
         
   if(b.sm == 0) return;
   
   /* REMOVE VERTEX CONTRIBUTION FROM SIDE MODES */
   /* SOLVE FOR SIDE MODES */
   /* PART 1 REMOVE VERTEX CONTRIBUTIONS */
   for(tind=0;tind<ntri;++tind) {
   
#ifdef CONSERV
      for(i=0;i<3;++i) {
         v0 = tvrtx[tind][i];
         uht[0][i] = gbl->res.v[v0][0]*gbl->tprcn[tind][0][0];
         uht[1][i] = gbl->res.v[v0][1]*gbl->tprcn[tind][0][0];
         uht[2][i] = gbl->res.v[v0][2]*gbl->tprcn[tind][NV-1][NV-1];
      }

      for(i=0;i<3;++i) {
         indx = tside[tind].side[i]*b.sm;
         sgn  = tside[tind].sign[i];
         for(j=0;j<3;++j) {
            indx1 = (i+j)%3;
            msgn = 1;
            for(k=0;k<b.sm;++k) {
               gbl->res.s[indx +k][0] -= msgn*b.vfms[j][k]*uht[0][indx1];
               gbl->res.s[indx +k][1] -= msgn*b.vfms[j][k]*uht[1][indx1];
               gbl->res.s[indx +k][2] -= msgn*b.vfms[j][k]*uht[2][indx1];
               msgn *= sgn;
            }
         }
      }
#else
      /* THIS IS TO USE THE NONCONSERVATIVE FORM */
      for(i=0;i<3;++i) {
         v0 = tvrtx[tind][i];
         uht[0][i] = gbl->res.v[v0][0]*gbl->tprcn[tind][0][0];
         uht[1][i] = gbl->res.v[v0][1]*gbl->tprcn[tind][0][0];
         uht[2][i] = gbl->res.v[v0][2];
      }

      for(i=0;i<3;++i) {
         indx = tside[tind].side[i]*b.sm;
         sgn  = tside[tind].sign[i];
         for(j=0;j<3;++j) {
            indx1 = (i+j)%3;
            msgn = 1;
            for(k=0;k<b.sm;++k) {
               gbl->res.s[indx +k][0] -= msgn*b.vfms[j][k]*(uht[0][indx1] +gbl->tprcn[tind][0][NV-1]*uht[2][indx1]);
               gbl->res.s[indx +k][1] -= msgn*b.vfms[j][k]*(uht[1][indx1] +gbl->tprcn[tind][1][NV-1]*uht[2][indx1]);
               gbl->res.s[indx +k][2] -= msgn*b.vfms[j][k]*(uht[2][indx1]*gbl->tprcn[tind][NV-1][NV-1]);
               msgn *= sgn;
            }
         }
      } 
#endif
   }
   
   /* SOLVE FOR LOWEST ORDER MODE */
#ifdef CONSERV
   indx = 0;
   for(sind = 0; sind < nside; ++sind) {
      gbl->res.s[indx][0] *= gbl->sprcn[sind][0][0]*b.sdiag[0];
      gbl->res.s[indx][1] *= gbl->sprcn[sind][0][0]*b.sdiag[0];
      gbl->res.s[indx][2] *= gbl->sprcn[sind][NV-1][NV-1]*b.sdiag[0];
      indx += b.sm;
   }
#else
   indx = 0;
   for(sind = 0; sind < nside; ++sind) {
      gbl->res.s[indx][0] = b.sdiag[0]*(gbl->res.s[indx][0]*gbl->sprcn[sind][0][0] +gbl->res.s[indx][2]*gbl->sprcn[sind][0][NV-1]);
      gbl->res.s[indx][1] = b.sdiag[0]*(gbl->res.s[indx][1]*gbl->sprcn[sind][0][0] +gbl->res.s[indx][2]*gbl->sprcn[sind][1][NV-1]);
      gbl->res.s[indx][2] *= b.sdiag[0]*gbl->sprcn[sind][NV-1][NV-1];
      indx += b.sm;
   }
#endif
   
   /* SEND MESSAGE FOR LOWEST ORDER MODE */
   bdry_ssnd(0);

   return;
}

/* call inline void hp_mgrid::minvrt3_mp() here */

void hp_mgrid::minvrt3(int mode) {  
   int i,j,m,n,indx,sind,tind;
   int sign[3],msgn,sgn,side[3];
   
   /* RECEIVE MESSAGE FOR MODE */
   /* APPLY DIRCHLET B.C.S TO MODE */
   bdry_srcvandzero(mode);

   /* REMOVE MODE FROM HIGHER MODES */
   for(tind=0;tind<ntri;++tind) {

#ifdef CONSERV
      for(i=0;i<3;++i) {
         side[i] = tside[tind].side[i]*b.sm;
         sign[i] = tside[tind].sign[i];
         sgn     = (mode % 2 ? sign[i] : 1);
         uht[0][i] = sgn*gbl->res.s[side[i]+mode][0]*gbl->tprcn[tind][0][0];
         uht[1][i] = sgn*gbl->res.s[side[i]+mode][1]*gbl->tprcn[tind][0][0];
         uht[2][i] = sgn*gbl->res.s[side[i]+mode][2]*gbl->tprcn[tind][NV-1][NV-1];
      }
      
      /* REMOVE MODES J,K FROM MODE I,M */
      for(i=0;i<3;++i) {
         msgn = ( (mode +1) % 2 ? sign[i] : 1);
         for(m=mode+1;m<b.sm;++m) {
            for(j=0;j<3;++j) {
               indx = (i+j)%3;
               for(n=0;n<NV;++n) {
                  gbl->res.s[side[i]+m][n] -= msgn*b.sfms[mode][m][j]*uht[n][indx];
               }
            }
            msgn *= sign[i];
         }
      }
#else
      for(i=0;i<3;++i) {
         side[i] = tside[tind].side[i]*b.sm;
         sign[i] = tside[tind].sign[i];
         sgn     = (mode % 2 ? sign[i] : 1);
         uht[0][i] = sgn*gbl->res.s[side[i]+mode][0]*gbl->tprcn[tind][0][0];
         uht[1][i] = sgn*gbl->res.s[side[i]+mode][1]*gbl->tprcn[tind][0][0];
         uht[2][i] = sgn*gbl->res.s[side[i]+mode][2];
      }
      
      /* REMOVE MODES J,K FROM MODE I,M */
      for(i=0;i<3;++i) {
         msgn = ( (mode +1) % 2 ? sign[i] : 1);
         for(m=mode+1;m<b.sm;++m) {
            for(j=0;j<3;++j) {
               indx = (i+j)%3;
               gbl->res.s[side[i]+m][0] -= msgn*b.sfms[mode][m][j]*(uht[0][indx] +gbl->tprcn[tind][0][NV-1]*uht[2][indx]);
               gbl->res.s[side[i]+m][1] -= msgn*b.sfms[mode][m][j]*(uht[1][indx] +gbl->tprcn[tind][1][NV-1]*uht[2][indx]);
               gbl->res.s[side[i]+m][1] -= msgn*b.sfms[mode][m][j]*(uht[2][indx]*gbl->tprcn[tind][NV-1][NV-1]);            
            }
            msgn *= sign[i];
         }
      }
#endif
   }
   
      /* SOLVE FOR NEXT MODE */
#ifdef CONSERV
      indx = mode +1;
      for(sind = 0; sind < nside; ++sind) {
         gbl->res.s[indx][0] *= gbl->sprcn[sind][0][0]*b.sdiag[mode+1];
         gbl->res.s[indx][1] *= gbl->sprcn[sind][0][0]*b.sdiag[mode+1];
         gbl->res.s[indx][2] *= gbl->sprcn[sind][NV-1][NV-1]*b.sdiag[mode+1];
         indx += b.sm;
      }
#else
      indx = mode +1;
      for(sind = 0; sind < nside; ++sind) {
         gbl->res.s[indx][0] = b.sdiag[mode+1]*(gbl->res.s[indx][0]*gbl->sprcn[sind][0][0] +gbl->res.s[indx][2]*gbl->sprcn[sind][0][NV-1]);
         gbl->res.s[indx][1] = b.sdiag[mode+1]*(gbl->res.s[indx][1]*gbl->sprcn[sind][0][0] +gbl->res.s[indx][2]*gbl->sprcn[sind][1][NV-1]);
         gbl->res.s[indx][2] *= b.sdiag[mode+1]*gbl->sprcn[sind][NV-1][NV-1];
         indx += b.sm;
      }
#endif

   return;
}

/* call inline void hp_mgrid::minvrt3_mp() here */


void hp_mgrid::minvrt4() {  
   int i,k,n,indx,tind;
   
   /* RECEIVE MESSAGE FOR LAST MODE */
   /* APPLY DIRICHLET B.C.'S */
   bdry_srcvandzero(b.sm-1);

   /* SOLVE FOR INTERIOR MODES */
   if (b.im > 0) {
      indx = 0;
      for(tind = 0; tind < ntri; ++tind) {
         DPBTRSNU(b.idiag,b.im,b.ibwth,&(gbl->res.i[indx][0]),NV);
         restouht_bdry(tind);
#ifdef CONSERV
         for(k=0;k<b.im;++k) {
            gbl->res.i[indx][0] /= gbl->tprcn[tind][0][0];
            gbl->res.i[indx][1] /= gbl->tprcn[tind][0][0];
            gbl->res.i[indx][2] /= gbl->tprcn[tind][NV-1][NV-1];
            
            for (i=0;i<b.bm;++i)
               for(n=0;n<NV;++n) 
                  gbl->res.i[indx][n] -= b.bfmi[i][k]*uht[n][i];
           
             ++indx;            
         }
#else      
         for(k=0;k<b.im;++k) {
            /* SUBTRACT BOUNDARY MODES (bfmi is multipled by interior inverse matrix so do this after DPBSLN) */
            for (i=0;i<b.bm;++i) {
                  gbl->res.i[indx][0] -= b.bfmi[i][k]*(uht[0][i]*gbl->tprcn[tind][0][0] +gbl->tprcn[tind][0][NV-1]*uht[2][i]);
                  gbl->res.i[indx][1] -= b.bfmi[i][k]*(uht[1][i]*gbl->tprcn[tind][0][0] +gbl->tprcn[tind][1][NV-1]*uht[2][i]);
                  gbl->res.i[indx][2] -= b.bfmi[i][k]*(uht[2][i]*gbl->tprcn[tind][NV-1][NV-1]);
            }
            
            /* INVERT PRECONDITIONER (Warning: tprcn is not preinverted like sprcn and vprcn) */
            gbl->res.i[indx][0] = (gbl->res.i[indx][0] -gbl->res.i[indx][2]*gbl->tprcn[tind][0][NV-1]/gbl->tprcn[tind][NV-1][NV-1])/gbl->tprcn[tind][0][0];
            gbl->res.i[indx][1] = (gbl->res.i[indx][1] -gbl->res.i[indx][2]*gbl->tprcn[tind][1][NV-1]/gbl->tprcn[tind][NV-1][NV-1])/gbl->tprcn[tind][0][0];
            gbl->res.i[indx][2] /= gbl->tprcn[tind][NV-1][NV-1];
            ++indx;            
         }
#endif
      }
   }

   return;
}

void hp_mgrid::restouht_bdry(int tind) {
    int i,m,n,indx,cnt;
    int sign, msgn;
   
   for (i=0; i<3; ++i) {
      indx = tvrtx[tind][i];
      for(n=0; n<NV; ++n)
         uht[n][i] = gbl->res.v[indx][n];
   }

   cnt = 3;
   for(i=0;i<3;++i) {
      indx = tside[tind].side[i]*b.sm;
      sign = tside[tind].sign[i];
      msgn = 1;
      for (m = 0; m < b.sm; ++m) {
         for(n=0; n<NV; ++n)
            uht[n][cnt] = msgn*gbl->res.s[indx +m][n];
         msgn *= sign;
         ++cnt;
      }
   }
     
   return;
}

void hp_mgrid::minvrt_test_bgn(FLT (*func)(int, FLT, FLT)) {
   int i,j,n,tind;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->res.v[i][n] = 0.0;

    for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         gbl->res.s[i][n] = 0.0;
                 
    for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         gbl->res.i[i][n] = 0.0;  
     

   for(tind = 0; tind<ntri;++tind) {
   
      if (tinfo[tind] > -1) {
         crdtocht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(cht[n], crd[n], dcrd[n][0], dcrd[n][1]);

      }
      else {
         for(n=0;n<ND;++n)
            b.proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);

         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd[n][0][i][j] = 0.5*(vrtx[tvrtx[tind][1]][n] -vrtx[tvrtx[tind][0]][n]);
                  dcrd[n][1][i][j] = 0.5*(vrtx[tvrtx[tind][2]][n] -vrtx[tvrtx[tind][0]][n]);
               }
            }
         }
      }
      
      ugtouht(tind);
      for(n=0;n<NV;++n)
         b.proj(uht[n],u[n]);
       
      for(n=0;n<NV;++n)
         for(i=0;i<b.tm;++i)
            lf[n][i] = 0.0;

      for(i=0;i<b.gpx;++i) {
         for(j=0;j<b.gpn;++j) {
            cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
            for(n=0;n<NV;++n)
               res[n][i][j] = RAD(i,j)*(u[n][i][j] -(*func)(n,crd[0][i][j],crd[1][i][j]))*cjcb[i][j];
         }
      }
      for(n=0;n<NV;++n)
         b.intgrt(res[n],lf[n]);
                    
      lftog(tind,gbl->res);
   }
   
   return;
}

void hp_mgrid::minvrt_test_tstep() {
   int tind,i,j,n,sind,count,bnum,side,v0,*v;
   FLT jcb,dtstari;
   class mesh *tgt;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         ug.v[i][n] = 0.0;

    for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         ug.s[i][n] = 0.0;
                 
    for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         ug.i[i][n] = 0.0;
         
   setinflow();
         
   /*	SET TIME STEP TO BE 1 */
   for(tind = 0; tind < ntri; ++tind) {
      jcb = 0.25*area(tind);
      v = tvrtx[tind];

      /* SET UP DIAGONAL PRECONDITIONER */
      dtstari = jcb*1;
#ifdef AXISYMMETRIC
      dtstari *= (vrtx[v[0]][0] +vrtx[v[1]][0] +vrtx[v[2]][0])/3.;
#endif
      gbl->tprcn[tind][0][0] = dtstari;      
      gbl->tprcn[tind][NV-1][NV-1] =  dtstari;
      for(i=0;i<3;++i) {
         gbl->vprcn[v[i]][0][0]  += gbl->tprcn[tind][0][0];
         gbl->vprcn[v[i]][NV-1][NV-1]  += gbl->tprcn[tind][NV-1][NV-1];
         if (b.sm > 0) {
            side = tside[tind].side[i];
            gbl->sprcn[side][0][0] += gbl->tprcn[tind][0][0];
            gbl->sprcn[side][NV-1][NV-1] += gbl->tprcn[tind][NV-1][NV-1];
         }
      }
   }
   
   /* SEND Y-DIRECTION BOUNDARY INFORMATION */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][NV-1][NV-1];

         /* SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][0];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][NV-1][NV-1];
            }
         }
      }
      
      if (sbdry[i].type & IFCE_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = gbl->vprcn[v0][0][0];

         /* SEND SIDE INFO */
         if (b.sm) {
            for(j=0;j<sbdry[i].num;++j) {
               sind = sbdry[i].el[j];
               tgt->sbuff[bnum][count++] = gbl->sprcn[sind][0][0];
            }
         }
      }
   }

   return;
}

void hp_mgrid::minvrt_test_end()  {
   int i,m,k,n,indx,indx1;
   FLT cflalpha;

   cflalpha = 1.0;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         ug.v[i][n] -= cflalpha*gbl->res.v[i][n];

   if (b.sm > 0) {
      indx = 0;
      indx1 = 0;
      for(i=0;i<nside;++i) {
         for (m=0;m<b.sm;++m) {
            for(n=0;n<NV;++n)
               ug.s[indx1][n] -= cflalpha*gbl->res.s[indx][n];
            ++indx;
            ++indx1;
         }
         indx1 += sm0 -b.sm;
      }         

      if (b.im > 0) {
         indx = 0;
         indx1 = 0;
         for(i=0;i<ntri;++i) {
            for(m=1;m<b.sm;++m) {
               for(k=0;k<b.sm-m;++k) {
                  for(n=0;n<NV;++n) {
                     ug.i[indx1][n] -= cflalpha*gbl->res.i[indx][n];
                  }
                  ++indx; ++indx1;
               }
               indx1 += sm0 -b.sm;
            }
         }
      }
   }
   
   return;
}


