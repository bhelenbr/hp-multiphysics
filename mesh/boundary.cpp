/*
 *  boundary.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Fri Jun 07 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include <float.h>
#include <assert.h>
#include <utilities.h>
#include "mesh.h"

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
side_boundary* mesh::getnewsideobject(int type) {
   if (type & PRDX_MASK)
      return(new prdx_boundary(*this,type));
   if (type & PRDY_MASK)
      return(new prdy_boundary(*this,type));
   if (type & (COMX_MASK+COMY_MASK))
      return(new comm_boundary(*this,type));    
   
   return(new side_boundary(*this,type));
}


void side_boundary::mvpttobdry(int indx, FLT psi, FLT pt[2]) {
   /* FOR A LINEAR SIDE */
   pt[0] = (1. -psi)*x.vrtx[x.svrtx[sd(indx)][0]][0] +psi*x.vrtx[x.svrtx[sd(indx)][1]][0];
   pt[1] = (1. -psi)*x.vrtx[x.svrtx[sd(indx)][0]][1] +psi*x.vrtx[x.svrtx[sd(indx)][1]][1];
   
   return;
}

void side_boundary::getgeometryfrommesh() {
   FLT length,xn,yn,x0,y0;
   int i;
   
   x0 = x.vrtx[x.svrtx[sd(0)][0]][0];
   y0 = x.vrtx[x.svrtx[sd(0)][0]][1];
   length = 0;
   
   for(i=0;i<nsd();++i) {
      s[i][0] = length;
      xn = b().vrtx[b().svrtx[sd(i)][1]][0];
      yn = b().vrtx[b().svrtx[sd(i)][1]][1];
      length += sqrt(pow(xn-x0,2) +pow(yn-y0,2));
      s[i][1] = length;
      x0 = xn;
      y0 = yn;
   }
   
   return;
}

void side_boundary::alloc(int n) {
   maxel = n;
   vect_alloc(el,n,int);
   s = (FLT (*)[ND]) xmalloc(2*n*sizeof(FLT));
}

     
void side_boundary::copy(const side_boundary &bin) {
   int i,n;
   
   if (!maxel) alloc(bin.maxel);
	else assert(bin.nel < maxel);
   assert(idnum == bin.idnum);
   
   nel = bin.nel;
   
   for(i=0;i<nel;++i)
      el[i] = bin.el[i];
      
   for(i=0;i<nel;++i)
      for(n=0;n<2;++n)
         s[i][n] = bin.s[i][n];
      
   return;
}

void side_boundary::output(FILE *out) {
	int i;
	
	fprintf(out,"idnty: %d\nnumber: %d\n",idnum,nel);
	
	for(i=0;i<nel;++i)
		fprintf(out,"%d %17.10e %17.10e\n",el[i],s[i][0],s[i][1]);
		
	return;
}

void side_boundary::input(FILE *in, FLT grwfac) {
	int i;
	
	fscanf(in,"number: %d\n",&nel);
	
	if (!maxel) alloc(static_cast<int>(grwfac*nel));
	else assert(nel < maxel);
	
	for(i=0;i<nel;++i)
		fscanf(in,"%d %lf %lf\n",&el[i],&s[i][0],&s[i][1]);
		
	return;
}

void side_boundary::findbdrypt(const class side_boundary *tgt,int ntgt,FLT psitgt,int *nout, FLT *psiout) {
   int top,bot;
   FLT s0;

   s0 = psitgt*tgt->s[ntgt][1] +(1.-psitgt)*tgt->s[ntgt][0];
   
   /* SEARCH BOUNDARY */
   bot = 0;
   top = nel-1;
   do {
      *nout = (top +bot)/2;
      *psiout = (s0-s[*nout][0])/(s[*nout][1] -s[*nout][0]);
      if (*psiout > 1.0 +10.*EPSILON)
         bot = *nout+1;
      else if (*psiout < 0.0 -10.*EPSILON)
         top = *nout-1;
      else
         break;
   } while (top >= bot);
   
   return;
}

/* SWAP ELEMENTS IN LIST */
void side_boundary::swap(int s1, int s2) {
   FLT temp[2];
   int ind;
   
   temp[0] = s[s1][0];
   temp[1] = s[s1][1];
   s[s1][0] = s[s2][0];
   s[s1][1] = s[s2][1];
   s[s2][0] = temp[0];
   s[s2][1] = temp[1];
   
   ind = el[s1];
   el[s1] = el[s2];
   el[s2] = ind;

   return;
}

  

/* REORDERS BOUNDARIES TO BE SEQUENTIAL */
/* USES INTWK1 & INTWK2 AS WORK ARRAYS */
side_boundary* side_boundary::reorder() {
   int i,count,total,sind,minv,first;

   total = nsd();
   
   /* STORE SIDE INDICES BY VERTEX NUMBER */
   for(i=0; i < nsd(); ++i) {
      sind = sd(i);
      x.intwk1[x.svrtx[sind][1]] = i;
      x.intwk2[x.svrtx[sind][0]] = i;
   }

   /* FIND FIRST SIDE */   
   first = -1;
   for(i=0;i<nsd();++i) {
      sind = sd(i);
      if (x.intwk1[x.svrtx[sind][0]] == -1) {
         first = i;
         break;
      }
   }
   
   /* SPECIAL CONSTRAINT IF LOOP */
   /* THIS IS TO ELIMINATE ANY INDEFINITENESS ABOUT SIDE ORDERING FOR LOOP */
   if (first < 0) {
      minv = x.nvrtx;
      for(i=0;i<nsd();++i) {
         sind = sd(i);
         if (x.svrtx[sind][1] < minv) {
            first = i;
            minv = x.svrtx[sind][1];
         }
      }
   }
   
   /* SWAP FIRST SIDE */
   count = 0;
   swap(count,first);
   x.intwk1[x.svrtx[sd(first)][1]] = first;
   x.intwk2[x.svrtx[sd(first)][0]] = first;
   x.intwk1[x.svrtx[sd(count)][1]] = count;
   x.intwk2[x.svrtx[sd(count)][0]] = -1;  // TO MAKE SURE LOOP STOPS

   /* REORDER LIST */
   while ((first = x.intwk2[x.svrtx[sd(count++)][1]]) >= 0) {
      swap(count,first);
      x.intwk1[x.svrtx[sd(first)][1]] = first;
      x.intwk2[x.svrtx[sd(first)][0]] = first;
      x.intwk1[x.svrtx[sd(count)][1]] = count;
      x.intwk2[x.svrtx[sd(count)][0]] = count;
   }
   
   /* RESET INTWK TO -1 */
   for(i=0; i <total; ++i) {
      sind = sd(i);
      x.intwk1[x.svrtx[sind][1]] = -1;
      x.intwk2[x.svrtx[sind][0]] = -1;
   }
   
   if (count < total) {
      side_boundary *temp = x.getnewsideobject(idnty());
      temp->copy(*this);
      nsd() = count;

      for(i=0;i<total-nsd();++i)
         temp->swap(i,i+nsd());
      temp->nsd() = total-nsd();
      return temp;
   }
   
   return 0;
}


/* GENERIC VERTEX COMMUNICATIONS */
void comm_boundary::send (FLT *base,int bgn,int end, int stride) {
   int j,k,sind,count,offset;
   
   count = 0;
   for(j=0;j<nsd();++j) {
      sind = sd(j);
      offset = b().svrtx[sind][0]*stride;
      for (k=bgn;k<=end;++k) 
         bdrymatch->sbuff[count++] = base[offset+k];
   }
   offset = b().svrtx[sind][1]*stride;
   for (k=bgn;k<=end;++k) 
      bdrymatch->sbuff[count++] = base[offset+k]; 
      
   bdrymatch->msgsize = nsd();
      
   return;
}

void comm_boundary::rcv (FLT *base,int bgn,int end, int stride) {
   int j,k,sind,count,offset;
   
   /*	RECV NUMBER */
   if (nsd() != msgsize) {
      printf("non matching number of comm boundaries %d: %d %d\n",idnty(),nsd(),msgsize);
      exit(1);
   }
   
   count = 0;
   for(j=nsd()-1;j>=0;--j) {
      sind = sd(j);
      offset = b().svrtx[sind][1]*stride;
      for (k=bgn;k<=end;++k) 
         base[offset+k] = 0.5*(base[offset+k] +sbuff[count++]);
   }
   offset = b().svrtx[sind][0]*stride;
   for (k=bgn;k<=end;++k) 
         base[offset+k] = 0.5*(base[offset+k] +sbuff[count++]);
}

void curv_boundary::mvpttobdry(int indx, FLT psi, FLT pt[2]) {
   int iter,n;
   FLT mag, delt_dist;
      
   /* FOR AN ANALYTIC DEFINED SIDE */
   for (n=0;n<ND;++n)
      pt[n] = (1. -psi)*b().vrtx[b().svrtx[sd(indx)][0]][n] +psi*b().vrtx[b().svrtx[sd(indx)][1]][n];

   iter = 0;
   do {
      mag = 0.0;
      for(n=0;n<ND;++n)
         mag += pow(dhgt(n,pt),2);
      mag = sqrt(mag);
      delt_dist = -hgt(pt)/mag;
      for(n=0;n<ND;++n)
         pt[n] += delt_dist*dhgt(n,pt)/mag;
      if (++iter > 100) {
         printf("iterations exceeded curved boundary %d %f %f\n",idnty(),pt[0],pt[1]);
         exit(1);
      }
   } while (fabs(delt_dist) > 10.*EPSILON);
   
   return;
}



