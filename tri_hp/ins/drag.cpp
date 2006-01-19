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
      
      for(ind=0; ind < sbdry(bnum)->nel; ++ind) {
         sind = sbdry(bnum)->el(ind);
         tind = sd(sind).tri(0);      
         
         for(sd=0;sd<3;++sd)
            if (td(tind).side(sd) == sind) break;
         assert(sd != 3);
         
         crdtocht(tind);
         for(m=basis::tri(log2p).bm;m<basis::tri(log2p).tm;++m)
            for(n=0;n<ND;++n)
               cht(n,m) = 0.0;
               
         for(n=0;n<ND;++n)
            basis::tri(log2p).proj_side(sd,&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0));

         ugtouht(tind);
 
         for(n=0;n<NV;++n)
            basis::tri(log2p).proj_side(sd,&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0));

         for (i=0;i<basis::tri(log2p).gpx;++i) {
            circ += basis::tri(log2p).wtx(i)*sqrt(dcrd(0,0)(0,i)*dcrd(0,0)(0,i) +dcrd(1,0)(0,i)*dcrd(1,0)(0,i));
            cjcb(0,i) = hp_gbl->mu*RAD1D(i)/(dcrd(0,0)(0,i)*dcrd(1,1)(0,i) -dcrd(1,0)(0,i)*dcrd(0,1)(0,i));
            
            /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
            /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
            visc[0][0][0][0] =  cjcb(0,i)*(2.*dcrd(1,1)(0,i)*dcrd(1,1)(0,i) +dcrd(0,1)(0,i)*dcrd(0,1)(0,i));
            visc[0][0][1][1] =  cjcb(0,i)*(2.*dcrd(1,0)(0,i)*dcrd(1,0)(0,i) +dcrd(0,0)(0,i)*dcrd(0,0)(0,i));
            visc[0][0][0][1] = -cjcb(0,i)*(2.*dcrd(1,1)(0,i)*dcrd(1,0)(0,i) +dcrd(0,1)(0,i)*dcrd(0,0)(0,i));
#define     viscI0II0II1II0I visc[0][0][0][1]

            visc[1][1][0][0] =  cjcb(0,i)*(dcrd(1,1)(0,i)*dcrd(1,1)(0,i) +2.*dcrd(0,1)(0,i)*dcrd(0,1)(0,i));
            visc[1][1][1][1] =  cjcb(0,i)*(dcrd(1,0)(0,i)*dcrd(1,0)(0,i) +2.*dcrd(0,0)(0,i)*dcrd(0,0)(0,i));
            visc[1][1][0][1] = -cjcb(0,i)*(dcrd(1,1)(0,i)*dcrd(1,0)(0,i) +2.*dcrd(0,1)(0,i)*dcrd(0,0)(0,i));
#define     viscI1II1II1II0I visc[1][1][0][1]
            
            visc[0][1][0][0] = -cjcb(0,i)*dcrd(0,1)(0,i)*dcrd(1,1)(0,i);
            visc[0][1][1][1] = -cjcb(0,i)*dcrd(0,0)(0,i)*dcrd(1,0)(0,i);
            visc[0][1][0][1] =  cjcb(0,i)*dcrd(0,1)(0,i)*dcrd(1,0)(0,i);
            visc[0][1][1][0] =  cjcb(0,i)*dcrd(0,0)(0,i)*dcrd(1,1)(0,i);

            /* OTHER SYMMETRIES    */            
#define     viscI1II0II0II0I visc[0][1][0][0]
#define     viscI1II0II1II1I visc[0][1][1][1]
#define     viscI1II0II0II1I visc[0][1][1][0]
#define     viscI1II0II1II0I visc[0][1][0][1]

            drg[0] -=   basis::tri(log2p).wtx(i)*(-u(2)(0,i)*RAD1D(i)*dcrd(1,0)(0,i) 
                        -viscI0II0II1II0I*du(0,0)(0,i) -visc[0][1][1][0]*du(1,0)(0,i)
                        -visc[0][0][1][1]*du(0,1)(0,i) -visc[0][1][1][1]*du(1,1)(0,i));															
            drg[1] -=   basis::tri(log2p).wtx(i)*( u(2)(0,i)*RAD1D(i)*dcrd(0,0)(0,i)
                        -viscI1II0II1II0I*du(0,0)(0,i) -viscI1II1II1II0I*du(1,0)(0,i)
                        -viscI1II0II1II1I*du(0,1)(0,i) -visc[1][1][1][1]*du(1,1)(0,i));
         }				
      }
      
      printf("#circumference, drag_x, drag_y: %e %e %e\n",circ,drg[0],drg[1]);
   }
   
	return;
}
