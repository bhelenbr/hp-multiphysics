#include "mesh.h"
#include <stdlib.h>

template int mesh<2>::comm_entity_size();
template int mesh<3>::comm_entity_size();

template<int ND> int mesh<ND>::comm_entity_size() {
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
   
   tsize += 1 +4*nscomm; // bdry number, id, v0id, v1id
   
   tsize += 1;  // nfcomm = 0
   
   return(tsize);
}
   
template int mesh<2>::comm_entity_list(int *list);
template int mesh<3>::comm_entity_list(int *list);

template<int ND> int mesh<ND>::comm_entity_list(int *list) {
   int i,j,nvcomm,nscomm,tsize,v0,v0id;
   
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
         list[tsize++] = vbdry[i]->idnty();
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
         list[tsize++] = sbdry[i]->idnty();
         v0 = svrtx[sbdry[i]->sd(0)][0];
         v0id = -1;
         for(j=0;j<nvbd;++j) {
            if (vbdry[j]->v() == v0) {
               v0id = vbdry[j]->idnty();
               break;
            }
         }
         list[tsize++] = v0id;
         v0 = svrtx[sbdry[i]->sd(sbdry[i]->nsd()-1)][1];
         v0id = -1;
         for(j=0;j<nvbd;++j) {
            if (vbdry[j]->v() == v0) {
               v0id = vbdry[j]->idnty();
               break;
            }
         }
         list[tsize++] = v0id;
      }
   }
   
   /* FACE BOUNDARIES */
   list[tsize++] = 0;
   
   return(tsize);
}

template int mesh<2>::msgpass(int phase);
template int mesh<3>::msgpass(int phase);

template<int ND> int mesh<ND>::msgpass(int phase) {
   int stop=1;

   for(int i=0;i<nsbd;++i) 
      stop &= sbdry[i]->rcv(phase);
   for(int i=0;i<nvbd;++i) 
      stop &= vbdry[i]->rcv(phase);
   
   
   if (!stop) {
      for(int i=0;i<nsbd;++i)
         sbdry[i]->snd(phase+1);
      for(int i=0;i<nvbd;++i)
         vbdry[i]->snd(phase+1);
   }
   
   return(stop);
}

template void mesh<2>::matchboundaries1();
template void mesh<3>::matchboundaries1();

/*	MAKE SURE MATCHING BOUNDARIES ARE AT EXACTLY THE SAME POSITIONS */
template<int ND> void mesh<ND>::matchboundaries1() {
   
   for(int i=0;i<nvbd;++i)
      vbdry[i]->sndpositions();
   for(int i=0;i<nsbd;++i) 
      sbdry[i]->sndpositions();
   
   /* FIRST PHASE OF SENDING */
   for(int i=0;i<nsbd;++i)
      sbdry[i]->snd(0);
   for(int i=0;i<nvbd;++i)
      vbdry[i]->snd(0);
}

template void mesh<2>::matchboundaries2();
template void mesh<3>::matchboundaries2();

template<int ND> void mesh<ND>::matchboundaries2() {
   for(int i=0;i<nsbd;++i)
      sbdry[i]->rcvpositions();
   for(int i=0;i<nvbd;++i)
      vbdry[i]->rcvpositions();
}

template void mesh<2>::length1();
template void mesh<3>::length1();

template<int ND> void mesh<ND>::length1() {
   int i;
   
   setlength();
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      sbdry[i]->loadbuff(vlngth,0,0,1);
   for(i=0;i<nvbd;++i)
      vbdry[i]->loadbuff(vlngth,0,0,1);
   

   for(i=0;i<nsbd;++i)
      sbdry[i]->snd(0);
   for(i=0;i<nvbd;++i)
      vbdry[i]->snd(0);

   return;
   
}

template void mesh<2>::length2();
template void mesh<3>::length2();

template<int ND> void mesh<ND>::length2() {
   int i;
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) 
      sbdry[i]->finalrcv(vlngth,0,0,1);
   for(i=0;i<nvbd;++i)
      vbdry[i]->finalrcv(vlngth,0,0,1);

   return;
}

#ifdef METIS
//#include <metis.h>

extern "C" void METIS_PartMeshNodal(int *ne, int *nn, int *elmnts, int *etype, int *numflag, int *nparts, int *edgecut,
int *epart, int *npart);

template void mesh<2>::setpartition(int nparts);
template void mesh<3>::setpartition(int nparts);

template<int ND> void mesh<ND>::setpartition(int nparts) {
   int etype = 1;
   int numflag = 0;
   int edgecut;
   METIS_PartMeshNodal(&ntri, &nvrtx, static_cast<int *>(&tvrtx[0][0]), &etype, &numflag, &nparts, &edgecut,tinfo,vinfo);
   return;
}
#endif

template void mesh<2>::partition(const class mesh& xmesh, int npart);
template void mesh<3>::partition(const class mesh& xmesh, int npart);

template<int ND> void mesh<ND>::partition(const class mesh& xin, int npart) {
   int i,j,n,tind,v0,indx;
   int bcntr[MAXSB];
   int bnum,bel,match;
   
   ntri = 0;
   for(i=0;i<xin.ntri;++i) {
      if (xin.tinfo[i] == npart) {
         ++ntri;
         for(n=0;n<3;++n)
            xin.vinfo[xin.tvrtx[i][n]] = npart;
      }
   }
   
   printf("New mesh with %d of %d tris\n",ntri,xin.ntri);
   
   if (!initialized) {
      maxvst = 3*ntri;
      allocate(maxvst);
   }
   else if (3*ntri > maxvst) {
      *log << "mesh is too large" << std::endl;
      exit(1);
   }

   nvrtx = 0;
   for(i=0;i<xin.nvrtx;++i) {
      if (xin.vinfo[i] == npart) {
         for(n=0;n<ND;++n)
            vrtx[nvrtx][n] = xin.vrtx[i][n];
         intwk1[i] = nvrtx;
         ++nvrtx;
      }
   }

   ntri = 0;
   for(i=0;i<xin.ntri;++i) {
      if (xin.tinfo[i] == npart) {
         for(n=0;n<3;++n)
            tvrtx[ntri][n] = intwk1[xin.tvrtx[i][n]];
         intwk2[ntri] = i;
         ++ntri;
      }
   }

   createsideinfo();
   
   /* MOVE BOUNDARY INFO */
   nvbd = 0;
   for(i=0;i<xin.nvbd;++i) {
      if (xin.vinfo[vbdry[i]->v()] == npart) {
         getnewvrtxobject(nvbd,xin.vbdry[i]->idnty());
         vbdry[nvbd]->alloc(1);
         vbdry[nvbd]->v() = xin.vinfo[vbdry[i]->v()];
         ++nvbd;
         if (nvbd >= MAXVB) {
            *log << "Too many vertex boundary conditions: increase MAXSB: " << nvbd << std::endl;
            exit(1);
         }
      }
   }
   
   
   nsbd = 0;
   for(i=0;i<nside;++i) {
      sinfo[i] = -1;
      if (stri[i][1] < 0) {
         tind = intwk2[stri[i][0]];

         v0 = svrtx[i][0];
         for(n=0;n<3;++n)
            if (intwk1[xin.tvrtx[tind][n]] == v0) break;
         if (n==3) printf("error in partitioning\n");
         n = (n+2)%3;

         indx = xin.ttri[tind][n];
         if (indx < 0) {
            /* BOUNDARY SIDE */
            bnum = (-indx>>16) -1;
            bel = -indx&0xFFFF;
            for (j = 0; j <nsbd;++j) {
               if (xin.sbdry[bnum]->idnty() == sbdry[j]->idnty()) {
                  ++bcntr[j];
                  sinfo[i] = j;
                  goto next1;
               }
            }
            /* NEW SIDE */
            getnewsideobject(nsbd,xin.sbdry[bnum]->idnty());
            sinfo[i] = nsbd;
            bcntr[nsbd++] = 1;

            if (nsbd > MAXSB) {
               *log << "error: too many different side boundaries: increase MAXSB" << std::endl;
               exit(1);
            }
         }
         else {
            /* PARTITION SIDE */
            match = tinfo[indx];
            if (match < npart) bnum = (match<<16) + (npart << 24) + (1<<8);
            else bnum = (npart<<16) + (match << 24) + (1<<8);
            for (j = 0; j <nsbd;++j) {
               if (sbdry[j]->idnty() == bnum) {
                  ++bcntr[j];
                  sinfo[i] = j;
                  goto next1;
               }
            }
            /* NEW SIDE */
            getnewsideobject(nsbd,bnum);
            sinfo[i] = nsbd;
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
      sbdry[i]->nsd() = 0;
   }       
   
   
   for(i=0;i<nside;++i) {
      if (sinfo[i] > -1) 
         sbdry[sinfo[i]]->sd(sbdry[sinfo[i]]->nsd()++) = i;
   }
   
   /* INTWK1,2 SHOULD ALWAYS BE RESET TO NEGATIVE 1 AFTER USE */
   for(i=0;i<xin.maxvst;++i)
      intwk1[i] = -1;
   for(i=0;i<xin.maxvst;++i)
      intwk2[i] = -1;
    
   for(i=0;i<nsbd;++i) {
      /* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
      sbdry[i]->reorder();
      sbdry[i]->getgeometryfrommesh();
   }
   
   bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT

   *log << "#Boundaries" << std::endl;
   for(i=0;i<nsbd;++i)
      sbdry[i]->summarize(*log);

   createttri();
   createvtri();
   cnt_nbor();
   if (!initialized) qtree.allocate(vrtx,maxvst);
   treeinit();
   initvlngth();

   initialized = 1;

   return;
}
            
   
   
   
         

