#include"hp_mgrid.h"
#include<assert.h>
 
void hp_mgrid::drag(int type) {
	int i,m,n,bnum,ind,sind,tind,sd;
	FLT drg[2];
	FLT visc[ND][ND][ND][ND];
   FLT circ;
   
   for(bnum=0;bnum<nsbd;++bnum) {
   
      if (sbdry[bnum].type != type) continue;

      drg[0] = 0.0;
      drg[1] = 0.0;
      circ = 0.0;
      
      for(ind=0; ind < sbdry[bnum].num; ++ind) {
         sind = sbdry[bnum].el[ind];
         tind = stri[sind][0];
         
         for(sd=0;sd<3;++sd)
            if (tside[tind].side[sd] == sind) break;
         assert(sd != 3);
         sd = (sd+1)%3; // BASIS SIDE ORDERING & MESH SIDE ORDERING ARE DIFFERENT
         
         crdtouht(tind);
         for(m=b.bm;m<b.tm;++m)
            for(n=0;n<ND;++n)
               uht[n][m] = 0.0;
               
         for(n=0;n<ND;++n)
            b.proj_side(sd,uht[n], crd[n][0], dcrd[n][0][0], dcrd[n][1][0]);

         ugtouht(tind);
         for(n=0;n<NV;++n)
            b.proj_side(sd,uht[n],u[n][0],du[n][0][0],du[n][1][0]);
         
         for (i=0;i<b.gpx;++i) {
            circ += b.wtx[i]*sqrt(dcrd[0][0][0][i]*dcrd[0][0][0][i] +dcrd[1][0][0][i]*dcrd[1][0][0][i]);
            cjcb[0][i] = gbl->mu/(dcrd[0][0][0][i]*dcrd[1][1][0][i] -dcrd[1][0][0][i]*dcrd[0][1][0][i]);
            
            printf("%d %f\n",sd,cjcb[0][i]/gbl->mu*(du[0][0][0][i]*dcrd[1][1][0][i] -du[0][1][0][i]*dcrd[1][0][0][i]));
//            printf("%d %f %f\n",sd,gbl->mu/cjcb[0][i],area(tind)/4.0);
//            printf("%d %f %f\n",sd,dcrd[0][1][0][i],0.5*(vrtx[tvrtx[tind][(sd+1)%3]][0]-vrtx[tvrtx[tind][(sd+2)%3]][0]));

            /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
            /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
            visc[0][0][0][0] =  cjcb[0][i]*(2.*dcrd[1][1][0][i]*dcrd[1][1][0][i] +dcrd[0][1][0][i]*dcrd[0][1][0][i]);
            visc[0][0][1][1] =  cjcb[0][i]*(2.*dcrd[1][0][0][i]*dcrd[1][0][0][i] +dcrd[0][0][0][i]*dcrd[0][0][0][i]);
            visc[0][0][0][1] = -cjcb[0][i]*(2.*dcrd[1][1][0][i]*dcrd[1][0][0][i] +dcrd[0][1][0][i]*dcrd[0][0][0][i]);
#define     viscI0II0II1II0I visc[0][0][0][1]

            visc[1][1][0][0] =  cjcb[0][i]*(dcrd[1][1][0][i]*dcrd[1][1][0][i] +2.*dcrd[0][1][0][i]*dcrd[0][1][0][i]);
            visc[1][1][1][1] =  cjcb[0][i]*(dcrd[1][0][0][i]*dcrd[1][0][0][i] +2.*dcrd[0][0][0][i]*dcrd[0][0][0][i]);
            visc[1][1][0][1] = -cjcb[0][i]*(dcrd[1][1][0][i]*dcrd[1][0][0][i] +2.*dcrd[0][1][0][i]*dcrd[0][0][0][i]);
#define     viscI1II1II1II0I visc[1][1][0][1]
            
            visc[0][1][0][0] = -cjcb[0][i]*dcrd[0][1][0][i]*dcrd[1][1][0][i];
            visc[0][1][1][1] = -cjcb[0][i]*dcrd[0][0][0][i]*dcrd[1][0][0][i];
            visc[0][1][0][1] =  cjcb[0][i]*dcrd[0][1][0][i]*dcrd[1][0][0][i];
            visc[0][1][1][0] =  cjcb[0][i]*dcrd[0][0][0][i]*dcrd[1][1][0][i];

            /* OTHER SYMMETRIES    */            
#define     viscI1II0II0II0I visc[0][1][0][0]
#define     viscI1II0II1II1I visc[0][1][1][1]
#define     viscI1II0II0II1I visc[0][1][1][0]
#define     viscI1II0II1II0I visc[0][1][0][1]


 //           printf("%d %d %f %f\n",i,sd,du[0][1][0][i],du[1][1][0][i]);

#ifdef SKIP
            printf("%d %d %f %f\n",i,sd,-viscI0II0II1II0I*du[0][0][0][i] -visc[0][1][1][0]*du[1][0][0][i]
                        -visc[0][0][1][1]*du[0][1][0][i] -visc[0][1][1][1]*du[1][1][0][i],-viscI1II0II1II0I*du[0][0][0][i] -viscI1II1II1II0I*du[1][0][0][i]
                        -viscI1II0II1II1I*du[0][1][0][i] -visc[1][1][1][1]*du[1][1][0][i]);

            printf("%d %d %f %f\n",i,sd,(8.-64.*pow(crd[0][0][i],2.0))*dcrd[1][0][0][i] -64.*crd[0][0][i]*crd[1][0][i]*dcrd[0][0][0][i],
            64.*crd[0][0][i]*crd[1][0][i]*dcrd[1][0][0][i] -(8.-64.*pow(crd[1][0][i],2.0))*dcrd[0][0][0][i]);
#endif

           
            drg[0] -=   b.wtx[i]*(-u[2][0][i]*dcrd[1][0][0][i] 
                        -viscI0II0II1II0I*du[0][0][0][i] -visc[0][1][1][0]*du[1][0][0][i]
                        -visc[0][0][1][1]*du[0][1][0][i] -visc[0][1][1][1]*du[1][1][0][i]);															
            drg[1] -=   b.wtx[i]*(-u[2][0][i]*dcrd[0][0][0][i] 
                        -viscI1II0II1II0I*du[0][0][0][i] -viscI1II1II1II0I*du[1][0][0][i]
                        -viscI1II0II1II1I*du[0][1][0][i] -visc[1][1][1][1]*du[1][1][0][i]);
         }				
      }
      
      printf("#drag: %e %e %e\n",drg[0],drg[1],circ);
   }
   
	return;
}
