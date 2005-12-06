#include "tri_hp.h"
#include "hp_boundary.h"
#include <myblas.h>

/************************************************/
/**********      INVERT MASS MATRIX    **********/
/************************************************/
block::ctrl tri_hp::minvrt(int excpt) {
   int i,j,k,m,n,tind,sind,v0,indx,indx1,indx2,sgn,msgn;
   TinyVector<int,3> sign,side;
   static int mode;
   Array<FLT,2> tinv(NV,NV);
   Array<FLT,1> temp(NV);
   
   if(basis::tri(log2p).sm == 0 && excpt > 2) return(block::stop);

   switch(excpt) {
      case(0): {
         /* LOOP THROUGH SIDES */
         if (basis::tri(log2p).sm > 0) {
            indx = 0;
            for(sind = 0; sind<nside;++sind) {
               /* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */         
               for (k=0; k <basis::tri(log2p).sm; ++k) {
                  for (i=0; i<2; ++i) {
                     v0 = sd(sind).vrtx(i);
                     for(n=0;n<NV;++n)
                        hp_gbl->res.v(v0,n) -= basis::tri(log2p).sfmv(i,k)*hp_gbl->res.s(sind,k,n);
                  }
                  ++indx;
               }
            }
               
            if (basis::tri(log2p).im > 0) {
               /* SUBTRACT INTERIORS */
               indx = 0;
               for(tind = 0; tind<ntri;++tind) {
                  indx2 = 3;
                  for (i=0; i<3; ++i) {
                     v0 = td(tind).vrtx(i);
                     for (k=0;k<basis::tri(log2p).im;++k)
                        for(n=0;n<NV;++n)
                           hp_gbl->res.v(v0,n) -= basis::tri(log2p).ifmb(i,k)*hp_gbl->res.i(tind,k,n);

                     sind = td(tind).side(i);
                     sgn = td(tind).sign(i);
                     msgn = 1;
                     for (j=0;j<basis::tri(log2p).sm;++j) {
                        for (k=0;k<basis::tri(log2p).im;++k)
                           for(n=0;n<NV;++n)
                              hp_gbl->res.s(sind,j,n) -= msgn*basis::tri(log2p).ifmb(indx2,k)*hp_gbl->res.i(tind,k,n);
                        msgn *= sgn;
                        ++indx2;
                     }
                  }
                 indx += basis::tri(log2p).im;
               }
            }
         }
         
         /* SOLVE FOR VERTEX MODES */
#ifndef MATRIX_PRECONDITIONER
         hp_gbl->res.v(Range(0,nvrtx-1),Range::all()) *= hp_gbl->vprcn(Range(0,nvrtx-1),Range::all());
#else
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<NV;++n) {
               temp(n) = hp_gbl->vprcn(i,n,0)*hp_gbl->res.v(i,0);
               for(m=1;m<NV;++m) {
                  temp(n) += hp_gbl->vprcn(i,n,m)*hp_gbl->res.v(i,m);
               }
            }
            hp_gbl->res.v(i,Range::all()) = temp(Range::all());
         }
#endif
            
         
         /* START MASS MATRIX INVERSION FOR COUPLED BOUNDARY EQUATIONS */
         /* TEMPORARY NEED TO UNDERSTAND HOW ENDPOINT MESSAGE PASSING WILL WORK */
         for(i=0;i<nsbd;++i)
            hp_sbdry(i)->minvrt(0);

         /* PREPARE MESSAGE PASSING */
         mp_phase = -1;
         return(block::advance);
      }
      
      case(1): {
         ++mp_phase;
         switch(mp_phase%3) {
            case(0):
               vc0load(mp_phase/3,hp_gbl->res.v.data());
               return(block::stay);
            case(1):
               vmsgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(vc0wait_rcv(mp_phase/3,hp_gbl->res.v.data())));
         }
      }
      
      case(2): {
         /* APPLY VERTEX DIRICHLET B.C.'S */
         for(i=0;i<nsbd;++i)
            hp_sbdry(i)->vdirichlet();
            
         /* FINISH INVERSION FOR COUPLED BOUNDARY EQUATIONS */
         /* TEMPORARY NEED TO UNDERSTAND HOW ENDPOINT MESSAGE PASSING WILL WORK */
         for(i=0;i<nsbd;++i)
            hp_sbdry(i)->minvrt(1);
               
         if(basis::tri(log2p).sm == 0) return(block::stop);
         
         /* REMOVE VERTEX CONTRIBUTION FROM SIDE MODES */
         /* SOLVE FOR SIDE MODES */
         /* PART 1 REMOVE VERTEX CONTRIBUTIONS */
         for(tind=0;tind<ntri;++tind) {
         
#ifndef MATRIX_PRECONDITIONER
            for(i=0;i<3;++i) {
               v0 = td(tind).vrtx(i);
               for(n=0;n<NV;++n)
                  uht(n)(i) = hp_gbl->res.v(v0,n)*hp_gbl->tprcn(tind,n);
            }

            for(i=0;i<3;++i) {
               sind = td(tind).side(i);
               sgn  = td(tind).sign(i);
               for(j=0;j<3;++j) {
                  indx1 = (i+j)%3;
                  msgn = 1;
                  for(k=0;k<basis::tri(log2p).sm;++k) {
                     for(n=0;n<NV;++n)
                        hp_gbl->res.s(sind,k,n) -= msgn*basis::tri(log2p).vfms(j,k)*uht(n)(indx1);
                     msgn *= sgn;
                  }
               }
            }
#else
            /* THIS IS TO USE A MATRIX PRECONDITIONER */
            for(i=0;i<3;++i) {
               v0 = td(tind).vrtx(i);
               for(n=0;n<NV;++n)
                  uht(n)(i) = hp_gbl->res.v(v0,n);
            }

            for(i=0;i<3;++i) {
               indx = td(tind).side(i);
               sgn  = td(tind).sign(i);
               for(j=0;j<3;++j) {
                  indx1 = (i+j)%3;
                  msgn = 1;
                  for(k=0;k<basis::tri(log2p).sm;++k) {
                     for(n=0;n<NV;++n) {
                        for(m=0;m<NV;++m) {
                           hp_gbl->res.s(indx,k,n) -= msgn*basis::tri(log2p).vfms(j,k)*hp_gbl->tprcn(tind,n,m)*uht(m)(indx1);
                        }
                     }
                     msgn *= sgn;
                  }
               }
            } 
#endif
         }
         
         /* SOLVE FOR LOWEST ORDER MODE */
#ifndef MATRIXPRECONDITIONER
         hp_gbl->res.s(Range(0,nside-1),0,Range::all()) *= hp_gbl->sprcn(Range(0,nside-1),Range::all())*basis::tri(log2p).sdiag(0);
#else
         for(sind = 0; sind < nside; ++sind) {
            for(n=0;n<NV;++n) {
               temp(n) = hp_gbl->sprcn(sind,n,0)*hp_gbl->res.s(sind,0,0);
               for(m=1;m<NV;++m) {
                  temp(n) += hp_gbl->sprcn(sind,n,m)*hp_gbl->res.s(sind,0,m);
               }
            }
            hp_gbl->res.s(sind,0,Range::all()) = temp(Range::all());
         }
#endif
         mode = 0;
         return(block::advance);
      }
   }
   
   /* INTERNALLY RESET TO ZERO FOR EASIER THINKING */
   excpt -= 2;
   
   /* THIS PART MUST BE REPEATED FOR EACH MODE EXCEPT FOR PHASE 2 OF LAST */
   if (excpt < 2*(basis::tri(log2p).sm-1)+1) {
      excpt = excpt%2;
      switch(excpt) {  
         case(0): {
            ++mp_phase;
            switch(mp_phase%3) {
               case(0):
                  sc0load(mp_phase/3,hp_gbl->res.s.data(),mode,mode,hp_gbl->res.s.extent(secondDim));
                  return(block::stay);
               case(1):
                  smsgpass(mp_phase/3);
                  return(block::stay);
               case(2):
                  return(static_cast<block::ctrl>(sc0wait_rcv(mp_phase/3,hp_gbl->res.s.data(),mode,mode,hp_gbl->res.s.extent(secondDim))));
            }
         }
         
         case(1): {
         
            /* APPLY DIRCHLET B.C.S TO MODE */
            for(i=0;i<nsbd;++i)
               hp_sbdry(i)->sdirichlet(mode);

            /* REMOVE MODE FROM HIGHER MODES */
            for(tind=0;tind<ntri;++tind) {

#ifndef MATRIX_PRECONDITIONER
               for(i=0;i<3;++i) {
                  side(i) = td(tind).side(i);
                  sign(i) = td(tind).sign(i);
                  sgn     = (mode % 2 ? sign(i) : 1);
                  for(n=0;n<NV;++n)
                     uht(n)(i) = sgn*hp_gbl->res.s(side(i),mode,n)*hp_gbl->tprcn(tind,n);
               }
               
               /* REMOVE MODES J,K FROM MODE I,M */
               for(i=0;i<3;++i) {
                  msgn = ( (mode +1) % 2 ? sign(i) : 1);
                  for(m=mode+1;m<basis::tri(log2p).sm;++m) {
                     for(j=0;j<3;++j) {
                        indx = (i+j)%3;
                        for(n=0;n<NV;++n) {
                           hp_gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p).sfms(mode,m,j)*uht(n)(indx);
                        }
                     }
                     msgn *= sign(i);
                  }
               }
#else
               for(i=0;i<3;++i) {
                  side(i) = td(tind).side(i);
                  sign(i) = td(tind).sign(i);
                  sgn     = (mode % 2 ? sign[i] : 1);
                  for(n=0;n<NV;++n)
                     uht(n)(i) = sgn*hp_gbl->res.s(side(i),mode,n)
               }
               
               /* REMOVE MODES J,K FROM MODE I,M */
               for(i=0;i<3;++i) {
                  msgn = ( (mode +1) % 2 ? sign(i) : 1);
                  for(m=mode+1;m<basis::tri(log2p).sm;++m) {
                     for(j=0;j<3;++j) {
                        indx = (i+j)%3;
                        for(n=0;n<NV;++n) {
                           hp_gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p).sfms(mode,m,j)*hp_gbl->tprcn(tind,n,0)*uht(0)(indx);
                           for(k=1;k<NV;++k) {
                              hp_gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p).sfms(mode,m,j)*hp_gbl->tprcn(tind,n,k)*uht(k)(indx);
                           }
                        }
                     }
                     msgn *= sign(i);
                  }
               }
#endif
            }
            
            /* SOLVE FOR NEXT MODE */
#ifndef MATRIX_PRECONDITIONER
            for(sind = 0;sind < nside;++sind)
               for (n=0;n<NV;++n)
                  hp_gbl->res.s(sind,mode+1,n) *= hp_gbl->sprcn(sind,n)*basis::tri(log2p).sdiag(mode+1);
#else
            for(sind = 0; sind < nside; ++sind) {
               for(n=0;n<NV;++n) {
                  temp(n) = hp_gbl->sprcn(sind,n,0)*hp_gbl->res.s(sind,mode+1,0);
                  for(m=1;m<NV;++m) {
                     temp(n) += hp_gbl->sprcn(sind,n,m)*hp_gbl->res.s(sind,mode+1,m);
                  }
               }
               hp_gbl->res.s(sind,mode+1,Range::all()) = temp(Range::all());
            }
#endif
            ++mode;
            return(block::advance);
         }
      }
   }
   
   /* IF FALL THROUGH TO HERE THEN MUST BE TIME TO SOLVE FOR INTERIOR MODES */
   if (excpt == 2*(basis::tri(log2p).sm-1)+1) {
      /* APPLY DIRCHLET B.C.S TO MODE */
      for(i=0;i<nsbd;++i)
         hp_sbdry(i)->sdirichlet(mode);
         
         /* SOLVE FOR INTERIOR MODES */
         if (basis::tri(log2p).im > 0) {
            for(tind = 0; tind < ntri; ++tind) {
               DPBTRSNU2(&basis::tri(log2p).idiag(0,0),basis::tri(log2p).ibwth+1,basis::tri(log2p).im,basis::tri(log2p).ibwth,&(hp_gbl->res.i(tind,0,0)),NV);
               restouht_bdry(tind);
#ifndef MATRIX_PRECONDITIONER
               for(k=0;k<basis::tri(log2p).im;++k) {
                  hp_gbl->res.i(tind,k,Range::all()) /= hp_gbl->tprcn(tind,Range::all());
                  
                  for (i=0;i<basis::tri(log2p).bm;++i)
                     for(n=0;n<NV;++n) 
                        hp_gbl->res.i(tind,k,n) -= basis::tri(log2p).bfmi(i,k)*uht(n)(i);
               }
#else      
               /* INVERT PRECONDITIONER (tprcn is not preinverted like sprcn and vprcn) */
               tinv = tprcn(tind,Range::all(),Range::all())
               GETRF(NV,NV,tinv.data(),NV,ipiv,info);
                  
               for(k=0;k<basis::tri(log2p).im;++k) {
                  /* SUBTRACT BOUNDARY MODES (bfmi is multipled by interior inverse matrix so do this after DPBSLN) */
                  for (i=0;i<basis::tri(log2p).bm;++i) {
                     for(n=0;n<NV;++n) {
                        for(m=0;m<NV;++m) {
                           hp_gbl->res.i(tind,k,n) -= basis::tri(log2p).bfmi(i,k)*uht(m)(i)*hp_gbl->tprcn(tind,n,m);
                        }
                     }
                  }
                  GETRS(trans,NV,1,tinv.data(),NV,ipiv,&hp_gbl->res.i(tind,k,0),NV,info);
               }
#endif
            }
         }
      }
      
      return(block::stop);
}

void tri_hp::restouht_bdry(int tind) {
    int i,m,n,indx,cnt;
    int sign, msgn;
   
   for (i=0; i<3; ++i) {
      indx = td(tind).vrtx(i);
      for(n=0; n<NV; ++n)
         uht(n)(i) = hp_gbl->res.v(indx,n);
   }

   cnt = 3;
   for(i=0;i<3;++i) {
      indx = td(tind).side(i);
      sign = td(tind).sign(i);
      msgn = 1;
      for (m = 0; m < basis::tri(log2p).sm; ++m) {
         for(n=0; n<NV; ++n)
            uht(n)(cnt) = msgn*hp_gbl->res.s(indx,m,n);
         msgn *= sign;
         ++cnt;
      }
   }
     
   return;
}

block::ctrl tri_hp::setup_preconditioner(int excpt) {
   int i,tind,side;
   TinyVector<int,3> v;
   FLT jcb,dtstari;

   switch (excpt) {
      case(0): {
         /*	SET TIME STEP TO BE 1 */
         
         for(tind = 0; tind < ntri; ++tind) {
            jcb = 0.25*area(tind);
            v = td(tind).vrtx;

            /* SET UP DIAGONAL PRECONDITIONER */
            dtstari = jcb*1;
#ifdef AXISYMMETRIC
            dtstari *= (vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.;
#endif
            hp_gbl->tprcn(tind,Range::all()) = dtstari;      
            
            for(i=0;i<3;++i) {
               hp_gbl->vprcn(v(i),Range::all())  += hp_gbl->tprcn(tind,Range::all());
               if (basis::tri(log2p).sm > 0) {
                  side = td(tind).side(i);
                  hp_gbl->sprcn(side,Range::all()) += hp_gbl->tprcn(tind,Range::all());
               }
            }
         }
         mp_phase = -1;
         return(block::advance);
      }
      
      case(1): {
         ++mp_phase;
         switch(mp_phase%3) {
            case(0):
               vc0load(mp_phase/3,hp_gbl->vprcn.data());
               return(block::stay);
            case(1):
               vmsgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(vc0wait_rcv(mp_phase/3,hp_gbl->vprcn.data())));
         }
      }
      
      case(2): {
         ++mp_phase;
         switch(mp_phase%3) {
            case(0):
               sc0load(mp_phase/3,hp_gbl->sprcn.data(),0,0,1);
               return(block::stay);
            case(1):
               smsgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(sc0wait_rcv(mp_phase/3,hp_gbl->sprcn.data(),0,0,1)));
         }
      }
      
      case(3): 
#ifndef MATRIX_PRECONDITIONER
         /* PREINVERT PRECONDITIONER FOR VERTICES */
         hp_gbl->vprcn(Range(0,nvrtx-1),Range::all()) = 1.0/(basis::tri(log2p).vdiag*hp_gbl->vprcn(Range(0,nvrtx-1),Range::all()));
        
         if (basis::tri(log2p).sm > 0) {
            /* INVERT DIAGANOL PRECONDITIONER FOR SIDES */            
            hp_gbl->sprcn(Range(0,nside-1),Range::all()) = 1.0/hp_gbl->sprcn(Range(0,nside-1),Range::all());
         }
#else
         /* NEED MATRIX INVERSION HERE */
#endif
   }
   return(block::stop);
}

   

block::ctrl tri_hp::minvrt_test(int excpt) {
   int i,j,k,m,n,tind,indx,indx1;
   TinyVector<int,3> v;
   block::ctrl step;
   TinyVector<FLT,mesh::ND> pt;
   

   if (excpt < 3) {
      return(setup_preconditioner(excpt));
   }
   
   switch(excpt) {
      case(3): {
         hp_gbl->res.v(Range(0,nvrtx-1),Range::all()) = 0.0;
         hp_gbl->res.s(Range(0,nside-1),Range::all(),Range::all()) = 0.0;
         hp_gbl->res.i(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
         
         ug.v(Range(0,nvrtx-1),Range::all()) = 0.0;
         ug.s(Range(0,nside-1),Range::all(),Range::all()) = 0.0;
         ug.i(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;

         for(tind = 0; tind<ntri;++tind) {
         
            if (td(tind).info > -1) {
               crdtocht(tind);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);

            }
            else {
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj(vrtx(td(tind).vrtx(0))(n),vrtx(td(tind).vrtx(1))(n),vrtx(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);

               for(i=0;i<basis::tri(log2p).gpx;++i) {
                  for(j=0;j<basis::tri(log2p).gpn;++j) {
                     for(n=0;n<ND;++n) {
                        dcrd(n,0)(i,j) = 0.5*(vrtx(td(tind).vrtx(1))(n) -vrtx(td(tind).vrtx(0))(n));
                        dcrd(n,1)(i,j) = 0.5*(vrtx(td(tind).vrtx(2))(n) -vrtx(td(tind).vrtx(0))(n));
                     }
                  }
               }
            }
             
            for(n=0;n<NV;++n)
               for(i=0;i<basis::tri(log2p).tm;++i)
                  lf(n)(i) = 0.0;

            for(i=0;i<basis::tri(log2p).gpx;++i) {
               for(j=0;j<basis::tri(log2p).gpn;++j) {
                  pt(0) = crd(0)(i,j);
                  pt(1) = crd(1)(i,j);
                  cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                  for(n=0;n<NV;++n)
                     res(n)(i,j) = RAD(i,j)*(*hp_gbl->initfunc)(n,pt)*cjcb(i,j);
               }
            }
            for(n=0;n<NV;++n)
               basis::tri(log2p).intgrt(&lf(n)(0),&res(n)(0,0),MXGP);
                          
            lftog(tind,hp_gbl->res);
         }
      }
   }
   
   
   step = minvrt(excpt-4);
   if(step != block::stop) return(step);
   
   /* Inversion finished */
   ug.v(Range(0,nvrtx-1),Range::all()) = hp_gbl->res.v(Range(0,nvrtx-1),Range::all());

   if (basis::tri(log2p).sm > 0) {
      ug.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = hp_gbl->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());
 
      if (basis::tri(log2p).im > 0) {

         for(i=0;i<ntri;++i) {
            indx = 0;
            indx1 = 0;
            for(m=1;m<basis::tri(log2p).sm;++m) {
               for(k=0;k<basis::tri(log2p).sm-m;++k) {
                  for(n=0;n<NV;++n) {
                     ug.i(i,indx1,n) = hp_gbl->res.i(i,indx,n);
                  }
                  ++indx; ++indx1;
               }
               indx1 += sm0 -basis::tri(log2p).sm;
            }
         }
      }
   }
   
   return(block::stop);
}


