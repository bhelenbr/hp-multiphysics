#include "boundaries.h"
#include <assert.h>
#include <float.h>

#ifdef MPISRC
#include <mpi.h>
#endif

/********************/
/* VERTEX FUNCTIONS */
/********************/
   
/* GENERIC VERTEX COMMUNICATIONS */
void vcomm::loadbuff(FLT *base,int bgn,int end, int stride) {
   int i,offset;
   
   sndsize()=end-bgn+1;
   sndtype()=flt_msg;
   
   /* LOAD SEND BUFFER */   
   offset = v0*stride +bgn;
   for (i=0;i<end-bgn+1;++i) 
      fsndalias[i] = base[offset+i];
}

void vcomm::finalrcv(FLT *base,int bgn,int end, int stride) {
   int i,m,offset;

   /* REMINDER DON'T OVERWRITE SEND BUFFER */
#ifdef MPISRC
   MPI_Status status;
   /* MPI PASSES */
   for(m=0;m<nmpi_match;++m)
      MPI_Wait(&mpi_rcvrqst[m], &status);         
#endif
   
   offset = v0*stride +bgn;
   for(i=0;i<end-bgn+1;++i) {
      base[offset] = fsndalias[i];  // ELIMINATES UNDESIRABLE SIDE/VRTX INTERACTIONS
      for(m=0;m<nlocal_match;++m)
         base[offset] += static_cast<FLT *>(local_rcv_buf[m])[i];

#ifdef MPISRC
      for(m=0;m<nmpi_match;++m)
         base[offset] += static_cast<FLT *>(mpi_rcv_buf[m])[i];
#endif
      base[offset++] /= (1 +nlocal_match +nmpi_match);
   }
}


/**************************************/
/* GENERIC FUNCTIONS FOR SIDES        */
/**************************************/
void side_bdry::alloc(int n) {
   maxel = n;
   el = new int[n];
}
     
void side_bdry::copy(const boundary& b) {
   int i;
   
   const side_bdry& bin = dynamic_cast<const side_bdry&>(b);
      
   if (!maxel) alloc(bin.maxel);
	else assert(bin.nel < maxel);
   
   nel = bin.nel;
   
   for(i=0;i<nel;++i)
      el[i] = bin.el[i];
      
   return;
}

void side_bdry::mvpttobdry(int indx, FLT psi, FLT pt[mesh::DIM]) {
   /* FOR A LINEAR SIDE */
   int n;
   
   for (n=0;n<mesh::DIM;++n)
      pt[n] = (1. -psi)*x.vrtx[x.sd[el[indx]].vrtx[0]][n] +psi*x.vrtx[x.sd[el[indx]].vrtx[1]][n];
   
   return;
}

#ifdef SKIP
void side_bdry::setupcoordinates() {
   int i,n;
   FLT length,x1[mesh::DIM],x0[mesh::DIM];
   
   for(n=0;n<mesh::DIM;++n)
      x0[n] = x.vrtx[x.sd[el[0]].vrtx[0]][n];
   
   s[0] = 0.0;
   for(i=0;i<nel;++i) {
      length = 0.0;
      for(n=0;n<mesh::DIM;++n) {
         x1[n] = x.vrtx[x.sd[el[i]].vrtx[1]][n];
         length += pow(x1[n]-x0[n],2);
         x0[n] = x1[n];
      }
      s[i+1] = s[0] +sqrt(length);
   }
   
   return;
}
#endif


#ifdef SKIP
void side_bdry::findbdrypt(const boundary *bin,int ntgt,FLT psitgt,int *nout, FLT *psiout) {
   int top,bot;
   FLT s0;

   const side_bdry *tgt = dynamic_cast<const side_bdry*>(bin);
   s0 = psitgt*tgt->s[ntgt+1] +(1.-psitgt)*tgt->s[ntgt];
   
   /* SEARCH BOUNDARY */
   bot = 0;
   top = nel-1;
   do {
      *nout = (top +bot)/2;
      *psiout = (s0-s[*nout])/(s[*nout +1] -s[*nout]);
      if (*psiout > 1.0 +10.*EPSILON)
         bot = *nout+1;
      else if (*psiout < 0.0 -10.*EPSILON)
         top = *nout-1;
      else
         break;
   } while (top >= bot);
   
   return;
}
#endif

/* SWAP ELEMENTS IN LIST */
void side_bdry::swap(int s1, int s2) {
   int ind;
   
   /* TEMPORARY NOT SURE HOW TO SWAP S VALUES */   
   ind = el[s1];
   el[s1] = el[s2];
   el[s2] = ind;

   return;
}


/* REORDERS BOUNDARIES TO BE SEQUENTIAL */
/* USES i1wk & i2wk AS WORK ARRAYS */
void side_bdry::reorder() {
   int i,count,total,sind,minv,first;

   total = nel;
   
   /* STORE SIDE INDICES BY VERTEX NUMBER */
   for(i=0; i < nel; ++i) {
      sind = el[i];
      x.i1wk[x.sd[sind].vrtx[1]] = i;
      x.i2wk[x.sd[sind].vrtx[0]] = i;
   }

   /* FIND FIRST SIDE */   
   first = -1;
   for(i=0;i<nel;++i) {
      sind = el[i];
      if (x.i1wk[x.sd[sind].vrtx[0]] == -1) {
         first = i;
         break;
      }
   }
   
   /* SPECIAL CONSTRAINT IF LOOP */
   /* THIS IS TO ELIMINATE ANY INDEFINITENESS ABOUT SIDE ORDERING FOR LOOP */
   if (first < 0) {
      minv = x.nvrtx;
      for(i=0;i<nel;++i) {
         sind = el[i];
         if (x.sd[sind].vrtx[1] < minv) {
            first = i;
            minv = x.sd[sind].vrtx[1];
         }
      }
   }
   
   /* SWAP FIRST SIDE */
   count = 0;
   swap(count,first);
   x.i1wk[x.sd[el[first]].vrtx[1]] = first;
   x.i2wk[x.sd[el[first]].vrtx[0]] = first;
   x.i1wk[x.sd[el[count]].vrtx[1]] = count;
   x.i2wk[x.sd[el[count]].vrtx[0]] = -1;  // TO MAKE SURE LOOP STOPS

   /* REORDER LIST */
   while ((first = x.i2wk[x.sd[el[count++]].vrtx[1]]) >= 0) {
      swap(count,first);
      x.i1wk[x.sd[el[first]].vrtx[1]] = first;
      x.i2wk[x.sd[el[first]].vrtx[0]] = first;
      x.i1wk[x.sd[el[count]].vrtx[1]] = count;
      x.i2wk[x.sd[el[count]].vrtx[0]] = count;
   }
   
   /* RESET INTWK TO -1 */
   for(i=0; i <total; ++i) {
      sind = el[i];
      x.i1wk[x.sd[sind].vrtx[1]] = -1;
      x.i2wk[x.sd[sind].vrtx[0]] = -1;
   }
   
   if (count < total) {

      ++x.nsbd;
      x.sbdry[x.nsbd-1] = create(x);
      x.sbdry[x.nsbd-1]->copy(*this);
      nel = count;

      for(i=0;i<total-nel;++i)
         x.sbdry[x.nsbd-1]->swap(i,i+nel);
      x.sbdry[x.nsbd-1]->nel = total-nel;
      *x.log << "#creating new boundary: " << idnum << " num: " << x.sbdry[x.nsbd-1]->nel << std::endl;
      return;
   }
   
   return;
}

void scomm::loadbuff(FLT *base,int bgn,int end, int stride) {
   int j,k,count,sind,offset;

   count = 0;
   for(j=0;j<nel;++j) {
      sind = el[j];
      offset = x.sd[sind].vrtx[0]*stride;
      for (k=bgn;k<=end;++k) 
         fsndalias[count++] = base[offset+k];
   }
   offset = x.sd[sind].vrtx[1]*stride;
   for (k=bgn;k<=end;++k) 
      fsndalias[count++] = base[offset+k]; 
      
   sndsize() = count;
   sndtype() = boundary::flt_msg;
}

void scomm::finalrcv(FLT *base,int bgn,int end, int stride) {
   int j,k,m,count,offset,sind;
#ifdef MPISRC
   MPI_Status status;
#endif
   /* ASSUMES REVERSE ORDERING OF SIDES */
   /* WON'T WORK IN 3D */
   
   /* RELOAD FROM BUFFER */
   /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
   /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */
   count = 0;
   for(j=0;j<nel;++j) {
      sind = el[j];
      offset = x.sd[sind].vrtx[0]*stride;
      for (k=bgn;k<=end;++k)
         base[offset+k] = fsndalias[count++];
   }
   offset = x.sd[sind].vrtx[1]*stride;
   for (k=bgn;k<=end;++k)
      base[offset+k] = fsndalias[count++];
   
   for(m=0;m<nlocal_match;++m) {            
      count = 0;
      for(j=nel-1;j>=0;--j) {
         sind = el[j];
         offset = x.sd[sind].vrtx[1]*stride;
         for (k=bgn;k<=end;++k) 
            base[offset+k] += static_cast<FLT *>(local_rcv_buf[m])[count++];
      }
      offset = x.sd[sind].vrtx[0]*stride;
      for (k=bgn;k<=end;++k) 
            base[offset+k] += static_cast<FLT *>(local_rcv_buf[m])[count++];            
   }
   

#ifdef MPISRC
   /* MPI PASSES */
   for(m=0;m<nmpi_match;++m) {
      MPI_Wait(&mpi_rcvrqst[m], &status);
      count = 0;
      for(j=nel-1;j>=0;--j) {
         sind = el[j];
         offset = x.sd[sind].vrtx[1]*stride;
         for (k=bgn;k<=end;++k) 
            base[offset+k] += static_cast<FLT *>(mpi_rcv_buf[m])[count++];
      }
      offset = x.sd[sind].vrtx[0]*stride;
      for (k=bgn;k<=end;++k) 
            base[offset+k] += static_cast<FLT *>(mpi_rcv_buf[m])[count++];
   }
#endif

   count = 0;
   for(j=0;j<nel;++j) {
      sind = el[j];
      offset = x.sd[sind].vrtx[0]*stride;
      for (k=bgn;k<=end;++k)
         base[offset+k] /= (1 +nlocal_match +nmpi_match);
   }
   offset = x.sd[sind].vrtx[1]*stride;
   for (k=bgn;k<=end;++k)
      base[offset+k] /= (1 +nlocal_match +nmpi_match);
}

void curved_analytic::mvpttobdry(int indx, FLT psi, FLT pt[mesh::DIM]) {
   
   /* GET LINEAR APPROXIMATION */
   side_bdry::mvpttobdry(indx,psi,pt);
   
   int iter,n;
   FLT mag, delt_dist;
      
   /* FOR AN ANALYTIC SURFACE */
   iter = 0;
   do {
      mag = 0.0;
      for(n=0;n<mesh::DIM;++n)
         mag += pow(dhgt(n,pt),2);
      mag = sqrt(mag);
      delt_dist = -hgt(pt)/mag;
      for(n=0;n<mesh::DIM;++n)
         pt[n] += delt_dist*dhgt(n,pt)/mag;
      if (++iter > 100) {
         *x.log << "iterations exceeded curved boundary " << idnum << ' ' << pt[0] << ' ' << pt[1] << '\n';
         exit(1);
      }
   } while (fabs(delt_dist) > 10.*EPSILON);
   
   return;
}