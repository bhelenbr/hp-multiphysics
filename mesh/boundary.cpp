#include "boundary.h"
#include <assert.h>
#include <float.h>

#ifdef MPISRC
#include <mpi.h>
#endif

/**************************************/
/* GENERIC FUNCTION FOR COMMUNICATION */
/**************************************/

template<class BASE, class MESH> void comm_boundary<BASE,MESH>::alloc(int size) {
   buffsize = size; 
   sndbuff = xmalloc(buffsize*sizeof(FLT)); 
   isndalias = static_cast<int *>(sndbuff);
   fsndalias = static_cast<FLT *>(sndbuff);        
}

template<class BASE, class MESH> void comm_boundary<BASE,MESH>::setphase(int phase, int msg_tag) {
   int i;
   for(i=0;i<nlocal_match;++i) {
      if (local_tags[i] == msg_tag) {
         local_phase[i] = phase;
         maxphase = MAX(maxphase,phase);
      }
   }

#ifdef MPISRC
   for(i=0;i<nmpi_match;++i) {
      if (mpi_tags[i] == msg_tag) {
         mpi_phase[i] = phase;
         maxphase = MAX(maxphase,phase);
      }
   }
#endif

   return;
}

/* MATCH BOUNDARIES */
template<class BASE, class MESH> int comm_boundary<BASE,MESH>::local_cnnct(boundary *bin, int msg_tag) {
   if (bin->idnty() == idnty()) {
      local_match[nlocal_match] = dynamic_cast<comm_boundary<BASE,MESH> *>(bin);
      local_tags[nlocal_match] = msg_tag;
      local_phase[nlocal_match] = 0; // DEFAULT ALL MESSAGE PASSING IN ONE PHASE
      local_rcv_buf[nlocal_match] = xmalloc(buffsize*sizeof(FLT));
      ++nlocal_match;
      return(0);
   }
   *b().log << "error: not local match" << idnty() << bin->idnty() << std::endl;
   return(1);
}

#ifdef MPISRC
template<class BASE, class MESH> int comm_boundary<BASE,MESH>::mpi_cnnct(int nproc, int msg_tag) {
   mpi_match[nmpi_match] = nproc;
   mpi_tags[nmpi_match] = msg_tag;
   mpi_phase[nmpi_match] = 0; // DEFAULT ALL MESSAGE PASSING IN ONE PHASE
   mpi_rcv_buf[nmpi_match] = xmalloc(buffsize*sizeof(FLT));
   ++nmpi_match;
   return(0);
}
#endif

/* GENERIC MECHANISM FOR SENDING */
template<class BASE, class MESH> int comm_boundary<BASE,MESH>::rcv(int phase) {
   int i,m;
   
   switch(msgtype) {
      case(boundary::flt_msg):
         /* LOCAL PASSES */
         for(m=0;m<nlocal_match;++m) {
            if (phase != local_phase[m]) continue;
            
            for(i=0;i<msgsize;++i)
               static_cast<FLT *>(local_rcv_buf[m])[i] = local_match[m]->fsndbuf(i);
         }

#ifdef MPISRC
         /* MPI PASSES */
         for(m=0;m<nmpi_match;++m) {
            if (phase != mpi_phase[m]) continue;
#ifdef SINGLE
            MPI_Isend(fsndalias, msgsize, MPI_FLOAT, 
               mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
#else
            MPI_Isend(fsndalias, msgsize, MPI_DOUBLE, 
               mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);             
#endif
         }
#endif
         break;
         
      case(boundary::int_msg):
         /* LOCAL PASSES */
         for(m=0;m<nlocal_match;++m) {
            if (phase != local_phase[m]) continue;
            
            for(i=0;i<msgsize;++i)
               static_cast<int *>(local_rcv_buf[m])[i] = local_match[m]->isndbuf(i);
         }

#ifdef MPISRC
         /* MPI PASSES */
         for(m=0;m<nmpi_match;++m) {
            if (phase != mpi_phase[m]) continue;
            
            MPI_Isend(isndalias, msgsize, MPI_INT, 
               mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_sndrqst[m]);
         }
#endif

         break;              
   }
   
   /* ONE MEANS FINISHED 0 MEANS MORE TO DO */
   return((phase-maxphase >= 0 ? 1 : 0));
}

template<class BASE, class MESH> void comm_boundary<BASE,MESH>::snd(int phase) {
#ifdef MPISRC
   int m;
#endif
   /* NOTHING TO DO FOR LOCAL PASSES (BUFFER ALREADY LOADED) */

#ifdef MPISRC
   switch(msgtype) {
      case(flt_msg):
         /* MPI POST RECEIVES FIRST */
         for(m=0;m<nmpi_match;++m) {
            if (phase != mpi_phase[m]) continue;
#ifdef SINGLE
            MPI_Irecv(static_cast<FLT *>(mpi_rcv_buf[m]), msgsize, MPI_FLOAT, 
               mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);
#else
            MPI_Irecv(static_cast<FLT *>(mpi_rcv_buf[m]), msgsize, MPI_DOUBLE, 
               mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);             
#endif
         }
         break;
         
      case(int_msg):
         /* MPI POST RECEIVES FIRST */
         for(m=0;m<nmpi_match;++m) {
            if (phase != mpi_phase[m]) continue;
            
            MPI_Irecv(static_cast<int *>(mpi_rcv_buf[m]), msgsize, MPI_INT, 
               mpi_match[m], mpi_tags[m],MPI_COMM_WORLD, &mpi_rcvrqst[m]);
         }
         break;              
   }
#endif

}


/**************************************/
/* GENERIC FUNCTIONS FOR SIDES        */
/**************************************/

template<class MESH> void side_template<MESH>::mvpttobdry(int indx, FLT psi, FLT pt[MESH::ND]) {
   /* FOR A LINEAR SIDE */
   int n;
   
   for (n=0;n<MESH::ND;++n)
      pt[n] = (1. -psi)*b().vrtx[b().svrtx[sd(indx)][0]][n] +psi*b().vrtx[b().svrtx[sd(indx)][1]][n];
   
   return;
}

template<class MESH> void side_template<MESH>::getgeometryfrommesh() {
   int i,n;
   FLT length,x1[MESH::ND],x0[MESH::ND];
   
   for(n=0;n<MESH::ND;++n)
      x0[n] = b().vrtx[b().svrtx[sd(0)][0]][n];
   
   s[0] = 0.0;
   for(i=0;i<nsd();++i) {
      length = 0.0;
      for(n=0;n<MESH::ND;++n) {
         x1[n] = b().vrtx[b().svrtx[sd(i)][1]][n];
         length += pow(x1[n]-x0[n],2);
         x0[n] = x1[n];
      }
      s[i+1] = s[0] +sqrt(length);
   }
   
   return;
}

template<class MESH> void side_template<MESH>::alloc(int n) {
   maxel = n;
   el = new int[n];
   s = new FLT[maxel+1];
}

     
template<class MESH> void side_template<MESH>::copy(const boundary& bin) {
   int i;
   
   const side_template<MESH>& temp = dynamic_cast<const side_template<MESH>&>(bin);
   
   if (!maxel) alloc(temp.maxel);
	else assert(temp.nel < maxel);
   
   nel = temp.nel;
   
   for(i=0;i<nel;++i)
      el[i] = temp.el[i];
      
   for(i=0;i<nel+1;++i)
      s[i] = temp.s[i];
      
   return;
}

template<class MESH> void side_template<MESH>::output(std::ostream& out) {
	int i;
	
   boundary::output(out);
   
	out << "number: " << nel << '\n';
	for(i=0;i<nel;++i)
		out << el[i] << ' ' << s[i] << ' ' << s[i+1] << '\n';
		
	return;
}

template<class MESH> void side_template<MESH>::input(FILE *in, FLT grwfac) {
	int i;
	
	fscanf(in,"number: %d\n",&nel);
	
	if (!maxel) alloc(static_cast<int>(grwfac*nel));
	else assert(nel < maxel);
	
	for(i=0;i<nel;++i)
		fscanf(in,"%d %lf %lf\n",&el[i],&s[i],&s[i+1]);
		
	return;
}

template<class MESH> void side_template<MESH>::findbdrypt(const boundary *bin,int ntgt,FLT psitgt,int *nout, FLT *psiout) {
   int top,bot;
   FLT s0;

   const side_template<MESH> *tgt = dynamic_cast<const side_template<MESH> *>(bin);
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

/* SWAP ELEMENTS IN LIST */
template<class MESH> void side_template<MESH>::swap(int s1, int s2) {
   int ind;
   
   /* TEMPORARY NOT SURE HOW TO SWAP S VALUES */   
   ind = el[s1];
   el[s1] = el[s2];
   el[s2] = ind;

   return;
}

  

/* REORDERS BOUNDARIES TO BE SEQUENTIAL */
/* USES INTWK1 & INTWK2 AS WORK ARRAYS */
template<class MESH> void side_template<MESH>::reorder() {
   int i,count,total,sind,minv,first;

   total = nsd();
   
   /* STORE SIDE INDICES BY VERTEX NUMBER */
   for(i=0; i < nsd(); ++i) {
      sind = sd(i);
      b().intwk1[b().svrtx[sind][1]] = i;
      b().intwk2[b().svrtx[sind][0]] = i;
   }

   /* FIND FIRST SIDE */   
   first = -1;
   for(i=0;i<nsd();++i) {
      sind = sd(i);
      if (b().intwk1[b().svrtx[sind][0]] == -1) {
         first = i;
         break;
      }
   }
   
   /* SPECIAL CONSTRAINT IF LOOP */
   /* THIS IS TO ELIMINATE ANY INDEFINITENESS ABOUT SIDE ORDERING FOR LOOP */
   if (first < 0) {
      minv = b().nvrtx;
      for(i=0;i<nsd();++i) {
         sind = sd(i);
         if (b().svrtx[sind][1] < minv) {
            first = i;
            minv = b().svrtx[sind][1];
         }
      }
   }
   
   /* SWAP FIRST SIDE */
   count = 0;
   swap(count,first);
   b().intwk1[b().svrtx[sd(first)][1]] = first;
   b().intwk2[b().svrtx[sd(first)][0]] = first;
   b().intwk1[b().svrtx[sd(count)][1]] = count;
   b().intwk2[b().svrtx[sd(count)][0]] = -1;  // TO MAKE SURE LOOP STOPS

   /* REORDER LIST */
   while ((first = b().intwk2[b().svrtx[sd(count++)][1]]) >= 0) {
      swap(count,first);
      b().intwk1[b().svrtx[sd(first)][1]] = first;
      b().intwk2[b().svrtx[sd(first)][0]] = first;
      b().intwk1[b().svrtx[sd(count)][1]] = count;
      b().intwk2[b().svrtx[sd(count)][0]] = count;
   }
   
   /* RESET INTWK TO -1 */
   for(i=0; i <total; ++i) {
      sind = sd(i);
      b().intwk1[b().svrtx[sind][1]] = -1;
      b().intwk2[b().svrtx[sind][0]] = -1;
   }
   
   if (count < total) {

      b().getnewsideobject(b().nsbd,idnty());
      ++b().nsbd;
      b().sbdry[b().nsbd-1]->copy(*this);
      nsd() = count;

      for(i=0;i<total-nsd();++i)
         b().sbdry[b().nsbd-1]->swap(i,i+nsd());
      b().sbdry[b().nsbd-1]->nsd() = total-nsd();
      *b().log << "#creating new boundary: " << idnty() << " num: " << b().sbdry[b().nsbd-1]->nsd() << std::endl;
      return;
   }
   
   return;
}



