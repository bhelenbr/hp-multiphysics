#include"mesh.h"
#include<cstdio>
#include<cfloat>
#include<iostream>

int mesh::setfine(const class mesh& tgt) {
   if (fine == NULL) 
      fine = new struct mg_trans[maxvst];
   return(mgconnect(fine,tgt));
}

int mesh::setcoarse(const class mesh& tgt) {
   if (coarse == NULL)
      coarse = new struct mg_trans[maxvst];
   return(mgconnect(coarse,tgt));
}

/*	THIS IS THE NEW WAY USING THE QUADTREE DATA STRUCTURE */
int mesh::mgconnect(struct mg_trans *cnnct, const class mesh& tgt) {
	int i,j,k,bnum,tind,sind,v0,error=0;
	double x,y,wgt[3],ainv;
	double dx,dy,a,b,c,minneg;
	int neg_count, triloc, sidloc;
   
/*	LOOP THROUGH VERTICES AND FIND SURROUNDING TRIANGLE */
   for(i=0;i<nvrtx;++i) {
      x = vrtx[i][0];
      y = vrtx[i][1];
      
      tgt.qtree.nearpt(x,y,v0);
      cnnct[i].tri = tgt.findtri(x,y,v0);
      tgt.getwgts(cnnct[i].wt);
   }

/*	REDO BOUNDARY SIDES TO DEAL WITH CURVATURE */
   for(bnum=0;bnum<nsbd;++bnum) {
   
/*		CHECK TO MAKE SURE THESE ARE THE SAME SIDES */
      if(sbdry[bnum].type != tgt.sbdry[bnum].type) {
         printf("error: sides are not numbered the same\n");
         exit(1);
      }
         
      for(k=0;k<sbdry[bnum].num;++k) {
         v0 = svrtx[sbdry[bnum].el[k]][0];
         x = vrtx[v0][0];
         y = vrtx[v0][1];
         minneg = -1.0E32;
         
   /* 	LOOP THROUGH TARGET SIDES TO FIND TRIANGLE */
         for(i=0;i<tgt.sbdry[bnum].num;++i) {
            sind = tgt.sbdry[bnum].el[i];
            tind = tgt.stri[sind][0];
            if (tind < 0) {
               printf("boundary side in wrong direction %d %d\n",sind,tind);
               exit(1);
            }
            neg_count = 0;
            for(j=0;j<3;++j) {
               wgt[j] = 
               ((tgt.vrtx[tgt.tvrtx[tind][(j+1)%3]][0]
                  -tgt.vrtx[tgt.tvrtx[tind][j]][0])*
               (y-tgt.vrtx[tgt.tvrtx[tind][j]][1])-
               (tgt.vrtx[tgt.tvrtx[tind][(j+1)%3]][1]
               -tgt.vrtx[tgt.tvrtx[tind][j]][1])*
               (x-tgt.vrtx[tgt.tvrtx[tind][j]][0]));
               if (wgt[j] < 0.0) ++neg_count;
            }
            if (neg_count==0) {
               sidloc = sind;
               break;
            }
            else {
               if(neg_count == 1) {
                  for(j=0;j<3;++j) {
                     if (wgt[j] < 0 && wgt[j] > minneg) {
                        minneg = wgt[j];
                        sidloc = sind;
                        triloc = tind;
                     }
                  }
               }
            }
         }
         sind = sidloc;	

   /*		PROJECT LOCATION NORMAL TO CURVED FACE */
         dx = x-tgt.vrtx[tgt.svrtx[sind][0]][0];
         dy = y-tgt.vrtx[tgt.svrtx[sind][0]][1];
         a = sqrt(dx*dx+dy*dy);
         dx = x-tgt.vrtx[tgt.svrtx[sind][1]][0];
         dy = y-tgt.vrtx[tgt.svrtx[sind][1]][1];
         b = sqrt(dx*dx+dy*dy);
         dx = tgt.vrtx[tgt.svrtx[sind][0]][0]
            -tgt.vrtx[tgt.svrtx[sind][1]][0];
         dy = tgt.vrtx[tgt.svrtx[sind][0]][1]
            -tgt.vrtx[tgt.svrtx[sind][1]][1];
         c = sqrt(dx*dx+dy*dy);
         a = (b*b+c*c-a*a)/(2.*c*c);
         x = tgt.vrtx[tgt.svrtx[sind][1]][0] +dx*a;
         y = tgt.vrtx[tgt.svrtx[sind][1]][1] +dy*a;
         tind = tgt.stri[sind][0];							
         for(j=0;j<3;++j) {
            wgt[j] = 
            ((tgt.vrtx[tgt.tvrtx[tind][(j+1)%3]][0]
            -tgt.vrtx[tgt.tvrtx[tind][j]][0])*
            (y-tgt.vrtx[tgt.tvrtx[tind][j]][1])-
            (tgt.vrtx[tgt.tvrtx[tind][(j+1)%3]][1]
            -tgt.vrtx[tgt.tvrtx[tind][j]][1])*
            (x-tgt.vrtx[tgt.tvrtx[tind][j]][0]));
         }
         cnnct[v0].tri = tind;
         ainv = 1.0/(tgt.area(tind));
         for (j=0;j<3;++j) {
            cnnct[v0].wt[(j+2)%3] = wgt[j]*ainv;
            if (wgt[j]*ainv > 1.0) 
               cnnct[v0].wt[(j+2)%3] = 1.0; 
               
            if (wgt[j]*ainv < 0.0)
               cnnct[v0].wt[(j+2)%3] = 0.0;
         }
      }
   }

	return(error);
}