#include "mesh.h"
#include "boundary.h"
#include <stdlib.h>

int mesh::comm_entity_size() {
   int i,tsize,nvcomm,nscomm;
   
   /*	VERTEX INFO */   
   nvcomm = 0;
   for(i=0;i<nvbd;++i) 
      if (vbdry[i]->is_comm()) ++nvcomm;
      
   tsize = 1 + 2*nvcomm; // bdry number & id
   
   /* SIDE INFO */
   nscomm = 0;
   for(i=0;i<nsbd;++i)
      if (sbdry[i]->is_comm()) ++nscomm; 
   
   // tsize += 1 +4*nscomm; // bdry number, id, v0id, v1id
   tsize += 1 +2*nscomm; // bdry number, id

   tsize += 1;  // nfcomm = 0
   
   return(tsize);
}

int mesh::comm_entity_list(int *list) {
   int i,nvcomm,nscomm,tsize;
   
   /* MAKE 1D PACKED LIST OF ALL INFORMATION ON COMMUNICATION BOUNDARIES */
   tsize = 0;

   /*	VERTEX INFO */   
   nvcomm = 0;
   for(i=0;i<nvbd;++i) 
      if (vbdry[i]->is_comm()) ++nvcomm;
      
   list[tsize++] = nvcomm;
   
   for(i=0;i<nvbd;++i) {
      if (vbdry[i]->is_comm()) {
         list[tsize++] = i;
         list[tsize++] = vbdry[i]->idnum;
      }
   }
   
   /* SIDE INFO */
   nscomm = 0;
   for(i=0;i<nsbd;++i)
      if (sbdry[i]->is_comm()) ++nscomm;
      
   list[tsize++] = nscomm;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i]->is_comm()) {
         list[tsize++] = i;
         list[tsize++] = sbdry[i]->idnum;
#ifdef SKIP
         v0 = sd[sbdry[i]->el[0]].vrtx[0];
         v0id = -1;
         for(j=0;j<nvbd;++j) {
            if (vbdry[j]->v0 == v0) {
               v0id = vbdry[j]->idnum;
               break;
            }
         }
         list[tsize++] = v0id;
         v0 = sd[sbdry[i]->el[sbdry[i]->nel-1]].vrtx[1];
         v0id = -1;
         for(j=0;j<nvbd;++j) {
            if (vbdry[j]->v0 == v0) {
               v0id = vbdry[j]->idnum;
               break;
            }
         }
         list[tsize++] = v0id;
#endif
      }
   }
   
   /* FACE BOUNDARIES */
   list[tsize++] = 0;
   
   return(tsize);
}

void mesh::msgload(int phase,FLT *base,int bgn, int end, int stride) {
   int i;
      
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      sbdry[i]->loadbuff(base,bgn,end,stride);
   for(i=0;i<nvbd;++i)
      vbdry[i]->loadbuff(base,bgn,end,stride);
   
   for(i=0;i<nsbd;++i)
      sbdry[i]->comm_prepare(phase);
   for(i=0;i<nvbd;++i)
      vbdry[i]->comm_prepare(phase);
   
   return;
}

void mesh::msgpass(int phase) {

   for(int i=0;i<nsbd;++i) 
      sbdry[i]->comm_transmit(phase);
   for(int i=0;i<nvbd;++i) 
      vbdry[i]->comm_transmit(phase);

   return;
}

int mesh::msgwait_rcv(int phase,FLT *base,int bgn, int end, int stride) {
   int stop = 1;
   int i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry[i]->comm_wait(phase);
   for(i=0;i<nvbd;++i)
      stop &= vbdry[i]->comm_wait(phase);
      

   for(i=0;i<nsbd;++i) 
      sbdry[i]->finalrcv(phase,base,bgn,end,stride);
   for(i=0;i<nvbd;++i)
      vbdry[i]->finalrcv(phase,base,bgn,end,stride);
      
   return(stop);
}

int mesh::msgrcv(int phase,FLT *base,int bgn, int end, int stride) {
   int stop = 1,i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry[i]->comm_nowait(phase);
   for(i=0;i<nvbd;++i)
      stop &= vbdry[i]->comm_nowait(phase);
      
   for(i=0;i<nsbd;++i) 
      sbdry[i]->finalrcv(phase,base,bgn,end,stride);
   for(i=0;i<nvbd;++i)
      vbdry[i]->finalrcv(phase,base,bgn,end,stride);
      
   return(stop);
}
   


/*	MAKE SURE MATCHING BOUNDARIES ARE AT EXACTLY THE SAME POSITIONS */
void mesh::matchboundaries1(int phase) {
   
   /* LOAD POSITIONS INTO BUFFERS */
   for(int i=0;i<nvbd;++i)
      vbdry[i]->loadpositions();
   for(int i=0;i<nsbd;++i) 
      sbdry[i]->loadpositions();
   
   /* FIRST PHASE OF SENDING, POST ALL RECEIVES */
   for(int i=0;i<nsbd;++i)
      sbdry[i]->comm_prepare(phase);
   for(int i=0;i<nvbd;++i)
      vbdry[i]->comm_prepare(phase);
}

int mesh::matchboundaries2(int phase) {
   int stop=1;

   /* FINAL PHASE OF SENDING */
   for(int i=0;i<nsbd;++i)
      stop &= sbdry[i]->comm_wait(phase);
   for(int i=0;i<nvbd;++i)
      stop &= vbdry[i]->comm_wait(phase);
      
   for(int i=0;i<nsbd;++i)
      sbdry[i]->rcvpositions(phase);
   for(int i=0;i<nvbd;++i)
      vbdry[i]->rcvpositions(phase);
            
   return(stop);
}

#ifdef METIS

extern "C" void METIS_PartMeshNodal(int *ne, int *nn, int *elmnts, int *etype, int *numflag, int *nparts, int *edgecut,
int *epart, int *npart);

void mesh::setpartition(int nparts) {
   int i,n;
   int etype = 1;
   int numflag = 0;
   int edgecut;
   
   /* CREATE ISOLATED TVRTX ARRAY */
   int (*tvrtx)[3] = static_cast<int (*)[3]>(scratch->p);
   for(i=0;i<ntri;++i)
      for(n=0;n<3;++n)
         tvrtx[i][n] = td[i].vrtx[n];
         
   METIS_PartMeshNodal(&ntri, &nvrtx, static_cast<int *>(&tvrtx[0][0]), &etype, &numflag, &nparts, &edgecut,i1wk,i2wk);
   delete []tvrtx;
   
   for(i=0;i<ntri;++i)
      td[i].info = i1wk[i];
            
   return;
}
#endif

void mesh::partition(const class mesh& xin, int npart) {
   int i,j,n,tind,v0,indx;
   int bcntr[MAXSB];
   int bnum,bel,match;
   
   for(i=0;i<xin.nvrtx;++i)
      xin.vd[i].info = -1;
      
   ntri = 0;
   for(i=0;i<xin.ntri;++i) {
      if (xin.td[i].info == npart) {
         ++ntri;
         for(n=0;n<3;++n)
            xin.vd[xin.td[i].vrtx[n]].info = npart;
      }
   }
   
   printf("New mesh with %d of %d tris\n",ntri,xin.ntri);
   
   if (!initialized) {
      maxvst = 3*ntri;
      allocate(maxvst,xin.scratch);
   }
   else if (3*ntri > maxvst) {
      *log << "mesh is too large" << std::endl;
      exit(1);
   }

   nvrtx = 0;
   for(i=0;i<xin.nvrtx;++i) {
      if (xin.vd[i].info == npart) {
         for(n=0;n<ND;++n)
            vrtx[nvrtx][n] = xin.vrtx[i][n];
         i1wk[i] = nvrtx;
         ++nvrtx;
      }
   }

   ntri = 0;
   for(i=0;i<xin.ntri;++i) {
      if (xin.td[i].info == npart) {
         for(n=0;n<3;++n)
            td[ntri].vrtx[n] = i1wk[xin.td[i].vrtx[n]];
         i2wk[ntri] = i;
         ++ntri;
      }
   }

   createsideinfo();
   
   /* MOVE BOUNDARY INFO */
   nvbd = 0;
   for(i=0;i<xin.nvbd;++i) {
      if (xin.vd[vbdry[i]->v0].info == npart) {
         vbdry[nvbd] = xin.vbdry[i]->create(*this);
         vbdry[nvbd]->alloc(1);
         vbdry[nvbd]->v0 = xin.vd[vbdry[i]->v0].info;
         ++nvbd;
         if (nvbd >= MAXVB) {
            *log << "Too many vertex boundary conditions: increase MAXSB: " << nvbd << std::endl;
            exit(1);
         }
      }
   }
   
   
   nsbd = 0;
   for(i=0;i<nside;++i) {
      sd[i].info = -1;
      if (sd[i].tri[1] < 0) {
         tind = i2wk[sd[i].tri[0]];

         v0 = sd[i].vrtx[0];
         for(n=0;n<3;++n)
            if (i1wk[xin.td[tind].vrtx[n]] == v0) break;
         if (n==3) printf("error in partitioning\n");
         n = (n+2)%3;

         indx = xin.td[tind].tri[n];
         if (indx < 0) {
            /* BOUNDARY SIDE */
            bnum = (-indx>>16) -1;
            bel = -indx&0xFFFF;

            for (j = 0; j <nsbd;++j) {
               if (xin.sbdry[bnum]->idnum == sbdry[j]->idnum) {
                  ++bcntr[j];
                  sd[i].info = j;
                  goto next1;
               }
            }
            /* NEW SIDE */
            sbdry[nsbd] = xin.sbdry[bnum]->create(*this);
            sd[i].info = nsbd;
            bcntr[nsbd++] = 1;

            if (nsbd > MAXSB) {
               *log << "error: too many different side boundaries: increase MAXSB" << std::endl;
               exit(1);
            }
         }
         else {
            /* PARTITION SIDE */
            match = xin.td[indx].info;
            if (match < npart) bnum = (match<<16) + (npart << 24) + stype::partition;
            else bnum = (npart<<16) + (match << 24) + stype::partition;
            for (j = 0; j <nsbd;++j) {
               if (sbdry[j]->idnum == bnum) {
                  ++bcntr[j];
                  sd[i].info = j;
                  goto next1;
               }
            }
            /* NEW SIDE */
            sbdry[nsbd] = getnewsideobject(bnum,0);
            sd[i].info = nsbd;
            bcntr[nsbd++] = 1;
            if (nsbd > MAXSB) {
               *log << "error: too many different side boundaries: increase MAXSB" << std::endl;
               exit(1);
            }
         }
      }
      next1: continue;
   }

   for(i=0;i<nsbd;++i) {
      sbdry[i]->alloc(static_cast<int>(bcntr[i]*2));
      sbdry[i]->nel = 0;
   }       
   
   
   for(i=0;i<nside;++i) {
      if (sd[i].info > -1) 
         sbdry[sd[i].info]->el[sbdry[sd[i].info]->nel++] = i;
   }
   
   /* i1wk,2 SHOULD ALWAYS BE RESET TO NEGATIVE 1 AFTER USE */
   for(i=0;i<xin.maxvst;++i)
      i1wk[i] = -1;
   for(i=0;i<xin.maxvst;++i)
      i2wk[i] = -1;
    
   for(i=0;i<nsbd;++i) {
      /* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
      sbdry[i]->reorder();
      sbdry[i]->setupcoordinates();
   }
   
   bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT

   createttri();
   createvtri();
   cnt_nbor();
   FLT xmin[ND], xmax[ND];
   for(n=0;n<ND;++n) {
      xmin[n] = xin.qtree.xmin(n);
      xmax[n] = xin.qtree.xmax(n);
   }
   treeinit(xmin,xmax);

   initialized = 1;

   return;
}
            
   
   
   
         

