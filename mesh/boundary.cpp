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
      fsndbuf(i) = base[offset+i];
}

void vcomm::finalrcv(FLT *base,int bgn,int end, int stride) {
   int i,m,offset;

   /* REMINDER DON'T OVERWRITE SEND BUFFER */
   offset = v0*stride +bgn;
   for(i=0;i<end-bgn+1;++i) {
      base[offset] = fsndbuf(i);  // ELIMINATES UNDESIRABLE SIDE/VRTX INTERACTIONS
      for(m=0;m<nmatch;++m)
         base[offset] += frcvbuf(m,i);

      base[offset++] /= (1 +nmatch);
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

block::ctrl side_bdry::mgconnect(int excpt, mesh::transfer *cnnct, const class mesh& tgt, int bnum) {
   int i,k,sind,tind,v0,neg_count,j,triloc,sidloc;
   FLT xp[2],dx,dy,wgt[3],ainv,a,b,c,minneg;
   
   if (excpt != 0) return(block::stop);
   
   /* BOUNDARY IS A PHYSICAL BOUNDARY (ALIGNED CURVES) */
   for(k=0;k<nel;++k) {
      v0 = x.sd[el[k]].vrtx[0];
      xp[0] = x.vrtx[v0][0];
      xp[1] = x.vrtx[v0][1];
      minneg = -1.0E32;
      
      /* LOOP THROUGH TARGET SIDES TO FIND TRIANGLE */
      for(i=0;i<tgt.sbdry[bnum]->nel;++i) {
         sind = tgt.sbdry[bnum]->el[i];
         tind = tgt.sd[sind].tri[0];
         if (tind < 0) {
            *x.log << "boundary side in wrong direction" << sind << tind << std::endl;
            exit(1);
         }
         neg_count = 0;
         for(j=0;j<3;++j) {
            wgt[j] = 
            ((tgt.vrtx[tgt.td[tind].vrtx[(j+1)%3]][0]
               -tgt.vrtx[tgt.td[tind].vrtx[j]][0])*
            (xp[1]-tgt.vrtx[tgt.td[tind].vrtx[j]][1])-
            (tgt.vrtx[tgt.td[tind].vrtx[(j+1)%3]][1]
            -tgt.vrtx[tgt.td[tind].vrtx[j]][1])*
            (xp[0]-tgt.vrtx[tgt.td[tind].vrtx[j]][0]));
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

      /* PROJECT LOCATION NORMAL TO CURVED FACE */
      dx = xp[0]-tgt.vrtx[tgt.sd[sind].vrtx[0]][0];
      dy = xp[1]-tgt.vrtx[tgt.sd[sind].vrtx[0]][1];
      a = sqrt(dx*dx+dy*dy);
      dx = xp[0]-tgt.vrtx[tgt.sd[sind].vrtx[1]][0];
      dy = xp[1]-tgt.vrtx[tgt.sd[sind].vrtx[1]][1];
      b = sqrt(dx*dx+dy*dy);
      dx = tgt.vrtx[tgt.sd[sind].vrtx[0]][0]
         -tgt.vrtx[tgt.sd[sind].vrtx[1]][0];
      dy = tgt.vrtx[tgt.sd[sind].vrtx[0]][1]
         -tgt.vrtx[tgt.sd[sind].vrtx[1]][1];
      c = sqrt(dx*dx+dy*dy);
      a = (b*b+c*c-a*a)/(2.*c*c);
      xp[0] = tgt.vrtx[tgt.sd[sind].vrtx[1]][0] +dx*a;
      xp[1] = tgt.vrtx[tgt.sd[sind].vrtx[1]][1] +dy*a;
      tind = tgt.sd[sind].tri[0];                     
      for(j=0;j<3;++j) {
         wgt[j] = 
         ((tgt.vrtx[tgt.td[tind].vrtx[(j+1)%3]][0]
         -tgt.vrtx[tgt.td[tind].vrtx[j]][0])*
         (xp[1]-tgt.vrtx[tgt.td[tind].vrtx[j]][1])-
         (tgt.vrtx[tgt.td[tind].vrtx[(j+1)%3]][1]
         -tgt.vrtx[tgt.td[tind].vrtx[j]][1])*
         (xp[0]-tgt.vrtx[tgt.td[tind].vrtx[j]][0]));
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
 
   return(block::stop);
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
      for (k=bgn;k<=end;++k) {
         fsndbuf(count++) = base[offset+k];
      }
   }
   offset = x.sd[sind].vrtx[1]*stride;
   for (k=bgn;k<=end;++k) 
      fsndbuf(count++) = base[offset+k]; 
      
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
         base[offset+k] = fsndbuf(count++);
   }
   offset = x.sd[sind].vrtx[1]*stride;
   for (k=bgn;k<=end;++k)
      base[offset+k] = fsndbuf(count++);
   
   for(m=0;m<nmatch;++m) {            
      count = 0;
      for(j=nel-1;j>=0;--j) {
         sind = el[j];
         offset = x.sd[sind].vrtx[1]*stride;
         for (k=bgn;k<=end;++k) {
            base[offset+k] += frcvbuf(m,count++);
         }
      }
      offset = x.sd[sind].vrtx[0]*stride;
      for (k=bgn;k<=end;++k) 
            base[offset+k] += frcvbuf(m,count++);            
   }
   
   count = 0;
   for(j=0;j<nel;++j) {
      sind = el[j];
      offset = x.sd[sind].vrtx[0]*stride;
      for (k=bgn;k<=end;++k)
         base[offset+k] /= (1 +nmatch);
   }
   offset = x.sd[sind].vrtx[1]*stride;
   for (k=bgn;k<=end;++k)
      base[offset+k] /= (1 +nmatch);
}

block::ctrl spartition::mgconnect(int excpt, mesh::transfer *cnnct, const class mesh& tgt, int bnum) {
   int i,j,k,v0;
   
   switch(excpt) {
      case(0):
         /* BOUNDARY IS AN INTERNAL PARTITION BOUNDARY */
         /* MAKE SURE ENDPOINTS ARE OK */
         i = x.sd[el[0]].vrtx[0];
         if (cnnct[i].tri < 0) {
            tgt.qtree.nearpt(x.vrtx[i],v0);
            cnnct[i].tri=tgt.vd[v0].tri;
            for(j=0;j<3;++j) {
               cnnct[i].wt[j] = 0.0;
               if (tgt.td[cnnct[i].tri].vrtx[j] == v0) cnnct[i].wt[j] = 1.0;
            }
         }
         i = x.sd[el[nel-1]].vrtx[1];
         if (cnnct[i].tri < 0) {
            tgt.qtree.nearpt(x.vrtx[i],v0);
            cnnct[i].tri=tgt.vd[v0].tri;
            for(j=0;j<3;++j) {
               cnnct[i].wt[j] = 0.0;
               if (tgt.td[cnnct[i].tri].vrtx[j] == v0) cnnct[i].wt[j] = 1.0;
            }
         }
         
         if (first) {
            sndsize() = 0;
            sndtype() = int_msg;
            for(k=1;k<nel;++k) {
               v0 = x.sd[el[k]].vrtx[0];
               if (cnnct[v0].tri > 0) {
                  isndbuf(sndsize()++) = -1;
               }
               else {
                  isndbuf(sndsize()++) = +1;
                  cnnct[v0].tri = 0;
                  for(j=0;j<3;++j)
                     cnnct[v0].wt[j] = 0.0;
               }
            }
            slave_master_prepare(); 
         }
         return(block::advance);                       
      case(1):
         slave_master_transmit();
         return(block::advance);
      case(2):
         slave_master_wait();
         if (!first) {
            i = 0;
            for(k=nel-1;k>0;--k) {
               v0 = x.sd[el[k]].vrtx[1];
               if (ircvbuf(0,i) < 0) {
                  cnnct[v0].tri = 0;
                  for(j=0;j<3;++j)
                     cnnct[v0].wt[j] = 0.0;
               }
            }
         }               
   }
   return(block::stop);
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