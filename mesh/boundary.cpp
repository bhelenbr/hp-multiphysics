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

void vcomm::finalrcv(int phi, FLT *base,int bgn,int end, int stride) {
   int i,m,offset;
   int matches = 1;
   
   for(m=0;m<nmatch;++m) {
      if (phase[m] != phi) continue;
      ++matches;
      
      for(i=0;i<end-bgn+1;++i) 
         fsndbuf(i) += frcvbuf(m,i);
   }
   
   if (matches > 1) {
      offset = v0*stride +bgn;
      for(i=0;i<end-bgn+1;++i) 
         base[offset++] = fsndbuf(i)/matches;
   }

   return;
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

void side_bdry::findbdryside(FLT *xpt, int &sidloc, FLT &psiloc) const {
   int k,sind,v0,v1,sidlocprev;
   FLT dx,dy,ol,psi,normdist;
   FLT psiprev,normdistprev;
   FLT mindist = 1.0e32;
      
   if (x.sd[el[0]].vrtx[0] == x.sd[el[nel-1]].vrtx[1]) {
      /* BOUNDARY IS A LOOP */
      sind = el[nel-1];
      v0 = x.sd[sind].vrtx[0];
      v1 = x.sd[sind].vrtx[1];
      dx = x.vrtx[v1][0] - x.vrtx[v0][0];
      dy = x.vrtx[v1][1] - x.vrtx[v0][1];
      ol = 2./(dx*dx +dy*dy);
      psi = ol*((xpt[0] -x.vrtx[v0][0])*dx +(xpt[1] -x.vrtx[v0][1])*dy) -1.;
      normdist = dx*(xpt[1]-x.vrtx[v0][1])-dy*(xpt[0]-x.vrtx[v1][0]);
      normdist *= sqrt(ol/2.);
      psiprev = psi;
      normdistprev = normdist;
      sidlocprev = nel-1;
   } 
   else {
      psiprev = -1.0;
   }
   
   for(k=0;k<nel;++k) {
      sind = el[k];
      v0 = x.sd[sind].vrtx[0];
      v1 = x.sd[sind].vrtx[1];
      dx = x.vrtx[v1][0] - x.vrtx[v0][0];
      dy = x.vrtx[v1][1] - x.vrtx[v0][1];
      ol = 2./(dx*dx +dy*dy);
      psi = ol*((xpt[0] -x.vrtx[v0][0])*dx +(xpt[1] -x.vrtx[v0][1])*dy) -1.;
      normdist = dx*(xpt[1]-x.vrtx[v0][1])-dy*(xpt[0]-x.vrtx[v1][0]);
      normdist *= sqrt(ol/2.);
      
      if (psi <= -1.0 && psiprev >= 1.0) {
         /* PREVIOUS & THIS SIDE ARE POTENTIAL MATCHES */
         if (fabs(normdist) < mindist) {
            mindist = fabs(normdist);
            sidloc = k;
            psiloc = -1.0;
         }
         if (fabs(normdistprev) < mindist) {
            mindist = fabs(normdistprev);
            sidloc = sidlocprev;
            psiloc = 1.0;
         }
      }
      else if (psi >= -1.0 && psi <= 1.0) {
         /* POTENTIAL SIDE */
         if (fabs(normdist) < mindist) {
            mindist = fabs(normdist);
            sidloc = k;
            psiloc = psi;
         }
      }
      psiprev = psi;
      normdistprev = normdist;
      sidlocprev = k;
   }
   
   return;
}



block::ctrl side_bdry::mgconnect(int excpt, mesh::transfer *cnnct, const class mesh& tgt, int bnum) {
   int j,k,sind,tind,v0,sidloc;
   FLT psiloc;
   
   if (excpt != 0) return(block::stop);
   
   for(k=1;k<nel;++k) {
      v0 = x.sd[el[k]].vrtx[0];
      tgt.sbdry[bnum]->findbdryside(x.vrtx[v0], sidloc, psiloc);
      sind = tgt.sbdry[bnum]->el[sidloc];
      tind = tgt.sd[sind].tri[0];                     
      cnnct[v0].tri = tind;
      for (j=0;j<3;++j) 
         if (tgt.td[tind].side[j] == sind) break;
      assert(j < 3);
      cnnct[v0].wt[j] = 0.0;
      cnnct[v0].wt[(j+1)%3] = 0.5*(1.-psiloc);
      cnnct[v0].wt[(j+2)%3] = 0.5*(1.+psiloc);
   } 
   
   return(block::stop);
}

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

void scomm::finalrcv(int phi, FLT *base,int bgn,int end, int stride) {
   int j,k,m,count,countdn,countup,offset,sind;
   FLT mtchinv;
   /* ASSUMES REVERSE ORDERING OF SIDES */
   /* WON'T WORK IN 3D */
   
   int matches = 1;
   
   /* RELOAD FROM BUFFER */
   /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
   /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */   
   for(m=0;m<nmatch;++m) {   
      if (phase[m] != phi) continue;
      
      ++matches;
      
      int ebp1 = end-bgn+1;
      countdn = nel*ebp1;
      countup = 0;
      for(j=0;j<nel+1;++j) {
         for(k=0;k<ebp1;++k)
            fsndbuf(countup +k) += frcvbuf(m,countdn +k);
         countup += ebp1;
         countdn -= ebp1;
      }
   }
   
   if (matches > 1) {
      mtchinv = 1./matches;

#ifdef MPDEBUG
      std::cout << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
      count = 0;
      for(j=0;j<nel;++j) {
         sind = el[j];
         offset = x.sd[sind].vrtx[0]*stride;
         for (k=bgn;k<=end;++k) {
            base[offset+k] = fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
            std::cout << "\t" << base[offset+k] << std::endl;
#endif
         }

      }
      offset = x.sd[sind].vrtx[1]*stride;
      for (k=bgn;k<=end;++k) {
         base[offset+k] = fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
         std::cout << "\t" << base[offset+k] << std::endl;
#endif
      }
   }
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

