/*
 *  blocks.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Thu Sep 26 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "blocks.h"
#include "boundary.h"
#include <time.h>
#include <input_map.h>
#include <utilities.h>
#include <iostream>
#include <map>
#include <string>
#include <sstream>

#ifdef MPISRC
#include <mpi.h>
#endif

#ifdef CAPRI
#include <capri.h>
#endif

void blocks::init(const char *infile, const char *outfile) {
   int i,nb,total;
   char ctemp[100];
   std::map<std::string,std::string> maptemp;
   
#ifdef MPISRC
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif   

   input_map(maptemp,infile);
   
   /* EXTRACT NBLOCKS FOR MYID */
   /* SPACE DELIMITED ARRAY OF NBLOCKS FOR EACH PROCESSOR */
   /* TO AVOID BUGS IN GCC 2.97 */
   total = 0;
   std::istringstream data(maptemp["nblock"]);
   for (i=0;i<myid;++i) {
      std::cout << i << std::endl;
      if (!(data >> nb)) {
         std::cout << "error reading blocks\n"; 
         exit(1);
      }
      total += nb;
   }
   if (!(data >> ctemp)) {
      std::cout << "error reading blocks\n"; 
      exit(1);
   }
   maptemp["nblock"] = ctemp;  // NUMBER OF BLOCKS FOR THIS PROCESSOR */
   data.clear();
   maptemp["bstart"] = total; // BLOCK NUMBER TO START FROM FOR THIS PROCESSOR

      
   if (outfile) {
      maptemp["logfile"] = outfile;
   }

   init(maptemp);
   
   return;
}
   


void blocks::init(std::map<std::string,std::string> input) {
   int i;
   char outfile[100];
   std::map<std::string,std::string>::const_iterator mi;
   std::map<std::string,std::string> merge;
   
#ifdef MPISRC
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
   
   
   std::istringstream data;

   /* OPEN LOGFILES FOR EACH PROCESSOR */
   mi = input.find("logfile");
   if (mi != input.end()) {
      data.str(mi->second);
      data >> outfile;
      data.clear();
      strcat(outfile,".");
      number_str(outfile,outfile,myid,1);
      strcat(outfile,".log");
      filelog.setf(std::ios::scientific, std::ios::floatfield);
      filelog.precision(3);
      filelog.open(outfile);
      log = &filelog;
   }
   else {
      std::cout.setf(std::ios::scientific, std::ios::floatfield);
      std::cout.precision(3);
      log = &std::cout;
   }
#ifdef CAPRI
   int status = gi_uStart();
   *log << "gi_uStart status = ", status, "\n";
   if (status != CAPRI_SUCCESS) exit(1);
   
   mi = input.find("BRep");
   if (mi != input.end()) {
      char temp[100];
      strcpy(temp,mi->second.c_str());
      status = gi_uLoadPart(temp);
      *log << mi->second << ": gi_uLoadPart status = " << status << std::endl;
      if (status != CAPRI_SUCCESS) exit(1);
   }
#endif   
   
   /* LOAD BASIC CONSTANTS FOR MULTIGRID */
   data.str(input["itercrsn"]);   
   if (!(data >> itercrsn)) itercrsn = 1;
   *log << "#itercrsn: " << itercrsn << std::endl;
   data.clear();
   
   data.str(input["iterrfne"]);   
   if (!(data >> iterrfne)) iterrfne = 0;
   *log << "#iterrfne: " << iterrfne << std::endl;
   data.clear();
      
   data.str(input["njacobi"]);   
   if (!(data >> njacobi)) njacobi = 1;
   *log << "#njacobi: " << njacobi << std::endl;
   data.clear();
   
   data.str(input["ncycle"]);   
   if (!(data >> ncycle)) ncycle = 20;
   *log << "#ncycle: " << ncycle << std::endl;
   data.clear();
   
   data.str(input["vwcycle"]);   
   if (!(data >> vw)) vw = 2;
   *log << "#vwcycle: " << vw << std::endl;
   data.clear();
   
   data.str(input["ntstep"]);   
   if (!(data >> ntstep)) ntstep = 1;
   *log << "#ntstep: " << ntstep << std::endl;
   data.clear();
   
   /* LOAD NUMBER OF GRIDS */
   data.str(input["nblock"]);
   if (!(data >> nblock)) nblock = 1;
   *log << "#nblock: " << nblock << std::endl;
   data.clear();
   
   data.str(input["ngrid"]);   
   if (!(data >> ngrid)) ngrid = 1;
   *log << "#ngrid: " << ngrid << std::endl;
   data.clear();
   
   data.str(input["mglvls"]);   
   if (!(data >> mglvls)) mglvls = 1;
   *log << "#mglvls: " << mglvls << std::endl;
   data.clear();
   
   int bstart;
   data.str(input["bstart"]);   
   if (!(data >> bstart)) bstart = 0;
   data.clear();

   blk = new block *[nblock];

   int lvl = 0;
   int excpt = 0;
   int state = block::stop;
   for (i=0;i<nblock;++i) {
      blk[i] = getnewblock(bstart+i,&input);
      blk[i]->init(input,log);
   }
   
   findmatch();
   matchboundaries();
   for(lvl=1;lvl<ngrid;++lvl) {
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i) {
            state &= blk[i]->reconnect(lvl,excpt);
         }
         excpt += state;
      } while (state != block::stop);
   }
   output("matched",ftype::grid);
   
   for(i=0;i<nblock;++i) {
      number_str(outfile,"coarsenchk",i,1);
      strcat(outfile,".");
      blk[i]->coarsenchk(outfile);
   }
   
   return;
}

static int maxblock,maxbdry,pdig,bdig,ndig;

int tagid(int vsf,int p1,int p2,int b1,int b2,int n1,int n2) {
   int tag = vsf;
   int shift =2;
   tag += (p1+p2)<<shift;
   shift += pdig;
   tag += (b1+b2)<<shift;
   shift += bdig;
   tag += (n1+n2)<<shift;
   shift += ndig;
   tag += abs(p1-p2)<<shift;
   shift += pdig;
   tag += abs(b1-b2)<<shift;
   shift += bdig;
   tag += abs(n1-n2)<<shift;
   return(tag);
}

void blocks::findmatch() {
   int b,b1,b2,i,j,k,grdlvl,p,count,tsize;
   int *entitylist, *sndentitylist, *size, *sndsize;
      
   struct commblocks {
      int nblock;
      struct commblock {
         int nvcomm;
         struct vid {
            int nvbd, idnum;
         } *vcomm;
         
         int nscomm;
         struct sid {
            int nsbd, idnum, vid0, vid1;
         } *scomm;
         
         int nfcomm;
         struct fid {
            int nfbd, idnum, nsds;
            int *sids;
         } *fcomm;
      } *blkinfo;
   } *blksinfo;
   commblocks::commblock *bp, *bp1, *bp2;
   blksinfo = new commblocks[nproc];
   sndsize = new int[nproc];
   for(i=0;i<nproc;++i)
      sndsize[i] = 0;

   /* FIRST DETERMINE TOTAL SIZE OF LIST */
   for(grdlvl=0;grdlvl<ngrid;++grdlvl) {
      
      sndsize[myid] = 1;  // 1 INT FOR # OF BLOCKS 
      for(b=0;b<nblock;++b) {
         sndsize[myid] += blk[b]->comm_entity_size(grdlvl);
      }
#ifdef MPISRC
      size = new int[nproc];
      for(i=0;i<nproc;++i)
         size[i] = 0;
      MPI_Allreduce(sndsize,size,nproc,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      delete []sndsize;
#else
      size = sndsize;
#endif

      tsize = 0;
      for(p=0;p<nproc;++p)
         tsize += size[p];
         
      sndentitylist = new int[tsize];
      for(i=0;i<tsize;++i)
         sndentitylist[i] = 0;

      /* BEGIN ASSEMBLING COMMUNICATION GRAPH LIST*/
      count = 0;
      for(p=0;p<myid;++p)
         count += size[p];
         
      sndentitylist[count++] = nblock;
      for(b=0;b<nblock;++b)
         count += blk[b]->comm_entity_list(grdlvl,&sndentitylist[count]);

#ifdef MPISRC   
      entitylist = new int[tsize];
      for(i=0;i<tsize;++i)
         entitylist[i] = 0;
        
      MPI_Allreduce(sndentitylist,entitylist,tsize,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
      delete []sndentitylist;
#else
      entitylist = sndentitylist;
#endif

      /* UNPACK ENTITY LIST INTO MORE USABLE FORM */
      maxblock = 0; // SO I CAN GENERATE UNIQUE TAGS
      maxbdry = 0; // SO I CAN GENERATE UNIQUE TAGS
      count = 0;
      for(p=0;p<nproc;++p) {
         blksinfo[p].nblock = entitylist[count++];
         maxblock = MAX(maxblock,blksinfo[p].nblock);
         blksinfo[p].blkinfo = new commblocks::commblock[blksinfo[p].nblock];
         for(b=0;b<blksinfo[p].nblock;++b) {
            bp = &blksinfo[p].blkinfo[b];
            bp->nvcomm = entitylist[count++];
            maxbdry = MAX(maxbdry,bp->nvcomm);
            bp->vcomm = new commblocks::commblock::vid[bp->nvcomm];
            for(i=0;i<bp->nvcomm;++i) {
               bp->vcomm[i].nvbd = entitylist[count++];
               bp->vcomm[i].idnum = entitylist[count++];
            }
            
            bp->nscomm = entitylist[count++];
            maxbdry = MAX(maxbdry,bp->nscomm);
            bp->scomm = new commblocks::commblock::sid[bp->nscomm];
            for(i=0;i<bp->nscomm;++i) {
               bp->scomm[i].nsbd = entitylist[count++];
               bp->scomm[i].idnum = entitylist[count++];
               bp->scomm[i].vid0 = entitylist[count++];
               bp->scomm[i].vid1 = entitylist[count++];
            }
            
            bp->nfcomm = entitylist[count++];
            maxbdry = MAX(maxbdry,bp->nfcomm);
            bp->fcomm = new commblocks::commblock::fid[bp->nfcomm];
            for(i=0;i<bp->nfcomm;++i) {
               bp->fcomm[i].nfbd = entitylist[count++];
               bp->fcomm[i].idnum = entitylist[count++];
               bp->fcomm[i].nsds = entitylist[count++];
               bp->fcomm[i].sids = new int[bp->fcomm[i].nsds];
               for (j=0;j<bp->fcomm[i].nsds;++j)
                  bp->fcomm[i].sids[j] = entitylist[count++];
            }
         }
      }
      
      /* CALCULATE NUMBER OF BINARY DIGITS NEEDED TO MAKE TAGS */
      pdig = 1;
      while ((nproc>>pdig) > 0) ++pdig;
      bdig = 1;
      while ((maxblock>>bdig) > 0) ++bdig;
      ndig = 1;
      while ((maxbdry>>ndig) > 0) ++ndig;
      if (2*(pdig+bdig+ndig)+2 > 32) *log << "can't guarantee unique tags\n";
      
      /* CAN NOW START TO LOOK FOR MATCHES */
      /* FOR NOW ALL COMMUNICATION IN 1 PHASE */
      /* LATER DO SOMETHING FANCY WITH ALL THIS INFO ? */
      
      for(b1=0;b1<blksinfo[myid].nblock;++b1) {
         bp1 = &blksinfo[myid].blkinfo[b1];

         for (p=0;p<nproc;++p) {
            for(b2=0;b2<blksinfo[p].nblock;++b2) {
               
               bp2 = &blksinfo[p].blkinfo[b2];
               /* LOOK FOR VERTEX MATCHES */
               for(i=0;i<bp1->nvcomm;++i) {
                  bool first_found = false;
                  for(j=0;j<bp2->nvcomm;++j) {
                     if (bp1->vcomm[i].idnum == bp2->vcomm[j].idnum) {
                        if (p == myid) {
                           if (bp1 == bp2 && i == j) {
                              if (!first_found) first_found = true;  // BOUNDARY'S FIRST FLAG SHOULD ALREADY BE SET
                              continue;  // CAN"T MATCH TO MYSELF
                           }
                           boundary *v1 = blk[b1]->vbdry(grdlvl,bp1->vcomm[i].nvbd);
                           boundary *v2 = blk[b2]->vbdry(grdlvl,bp2->vcomm[j].nvbd);
                           v1->local_cnnct(v2,tagid(1,myid,myid,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false; 
                              first_found = true;
                           }
                           *log << "#mgrid " << grdlvl << "vrtx local match id: " << v1->idnum << " b1: " << b1 << " b2: " << b2 << " n1: " << bp1->vcomm[i].nvbd << " n2: " << bp2->vcomm[j].nvbd << "\n";
                        }
#ifdef MPISRC
                        else {
                           boundary *v1 = blk[b1]->vbdry(grdlvl,bp1->vcomm[i].nvbd);
                           v1->mpi_cnnct(p,tagid(1,myid,p,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false;
                              first_found = true;
                           }
                           *log << "#mgrid " << grdlvl << "vrtx mpi match " << "p: " << p << " id: " << v1->idnum << " b1: " << b1 << " b2: " << b2 << " n1: " << bp1->vcomm[i].nvbd << " n2: " << bp2->vcomm[j].nvbd << "\n";
                        }
#endif
                        
                     }
                  }
               }
               
               /* LOOK FOR SIDE MATCHES */
               for(i=0;i<bp1->nscomm;++i) {
                  bool first_found = false;
                  for(j=0;j<bp2->nscomm;++j) {
                     if (bp1->scomm[i].idnum == bp2->scomm[j].idnum) {
                        if (p == myid) {
                           if (bp1 == bp2 && i == j) {
                              if (!first_found) first_found = true;  // BOUNDARY'S FIRST FLAG SHOULD ALREADY BE SET
                              continue;  // CAN"T MATCH TO MYSELF
                           }
                           boundary *v1 = blk[b1]->sbdry(grdlvl,bp1->scomm[i].nsbd);
                           boundary *v2 = blk[b2]->sbdry(grdlvl,bp2->scomm[j].nsbd);
                           v1->local_cnnct(v2,tagid(2,myid,myid,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false;
                              first_found = true;
                           }
                           *log << "#mgrid " << grdlvl << "side local match id: " << v1->idnum << " b1: " << b1 << " b2: " << b2 << " n1: " << bp1->vcomm[i].nvbd << " n2: " << bp2->vcomm[j].nvbd << "\n";

                        }
#ifdef MPISRC
                        else {
                           boundary *v1 = blk[b1]->sbdry(grdlvl,bp1->scomm[i].nsbd);
                           v1->mpi_cnnct(p,tagid(2,myid,p,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false;
                              first_found = true;
                           }
                           *log << "#mgrid " << grdlvl << "side mpi match " << "p: " << p << " id: " << v1->idnum << " b1: " << b1 << " b2: " << b2 << " n1: " << bp1->vcomm[i].nvbd << " n2: " << bp2->vcomm[j].nvbd << "\n";

                        }
#endif
                     }
                  }
               }
               
               /* LOOK FOR FACE MATCHES */
               for(i=0;i<bp1->nfcomm;++i) {
                  bool first_found = false;
                  for(j=0;j<bp2->nfcomm;++j) {
                     if (bp1->fcomm[i].idnum == bp2->fcomm[j].idnum) {
                        if (p == myid) {
                           if (bp1 == bp2 && i == j) {
                              if (!first_found) first_found = true;  // BOUNDARY'S FIRST FLAG SHOULD ALREADY BE SET
                              continue;  // CAN"T MATCH TO MYSELF
                           }
                           boundary *v1 = blk[b1]->fbdry(grdlvl,bp1->fcomm[i].nfbd);
                           boundary *v2 = blk[b2]->fbdry(grdlvl,bp2->fcomm[j].nfbd);
                           v1->local_cnnct(v2,tagid(3,myid,myid,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false;
                              first_found = true;
                           }
                           *log << "#mgrid " << grdlvl << "face local match id: " << v1->idnum << " b1: " << b1 << " b2: " << b2 << " n1: " << bp1->vcomm[i].nvbd << " n2: " << bp2->vcomm[j].nvbd << "\n";

                        }
#ifdef MPISRC
                        else {
                           boundary *v1 = blk[b1]->fbdry(grdlvl,bp1->fcomm[i].nfbd);
                           v1->mpi_cnnct(p,tagid(1,myid,p,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false;
                              first_found = true;
                           }
                           *log << "#mgrid " << grdlvl << "face mpi match " << "p: " << p << " id: " << v1->idnum << " b1: " << b1 << " b2: " << b2 << " n1: " << bp1->vcomm[i].nvbd << " n2: " << bp2->vcomm[j].nvbd << "\n";

                        }
#endif
                     }
                  }
               }
            }
         }
      }
               
      /* DELETE DATA STRUCTURE */
      for(i=0;i<nproc;++i) {
         for(j=0;j<blksinfo[i].nblock;++j) {
            bp = &blksinfo[i].blkinfo[j];
            delete []bp->vcomm;
            delete []bp->scomm;
            for(k=0;k<bp->nfcomm;++k)
               delete []bp->fcomm[k].sids;
            delete []bp->fcomm;
         }
         delete []blksinfo[i].blkinfo;
      }
      delete []entitylist;
   }
   
   delete []blksinfo;
   delete []size;

   return;
}

            
            
   /* MATCH BOUNDARIES */
void blocks::matchboundaries() {
   int i,lvl,excpt;
   int state;
   
   for (lvl=0;lvl<1;++lvl) {
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->matchboundaries(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
   }
}


void blocks::output(char *filename, ftype::name filetype) {
   int i;   
   char fnmcat[80], fnmcat1[80];
   
   strcpy(fnmcat,filename);
#ifdef MPISRC
   strcat(fnmcat,".");
   number_str(fnmcat, fnmcat, myid, 1);
#endif

   /* ASSUME FOR NOW MESHES ARE LABELED a,b,c... */
   /* I HAVEN'T FIGURED OUT HOW THIS IS GOING TO WORK IN THE TOTALLY GENERAL CASE */
   if (nblock > 1) {
      strcat(fnmcat,".");
      for (i=0;i<nblock;++i) {
         number_str(fnmcat1, fnmcat, i, 1);
         blk[i]->output(fnmcat1,filetype);
      }
   }
   else {
      blk[0]->output(fnmcat,filetype);
   }
   
   return;
}

void blocks::rsdl(int lvl) {
   int i,excpt;
   int state;
   
   excpt = 0;
   do {
      state = block::stop;
      for(i=0;i<nblock;++i)
         state &= blk[i]->rsdl(lvl,excpt);
      excpt += state;
   } while (state != block::stop);

   return;
}
      
      


void blocks::iterate(int lvl) {
/*****************************************/
/* JACOBI-ITERATION FOR MESH POSITION ****/
/*****************************************/
   int i,iter,excpt;
   int state;
   
   excpt = 0;
   do {
      state = block::stop;
      for(i=0;i<nblock;++i)
         state &= blk[i]->setup_preconditioner(lvl,excpt);
      excpt += state;
   } while (state != block::stop);

   for(iter=0;iter<njacobi;++iter) {
      rsdl(lvl);
   
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->update(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
   }

   return;
}

void blocks::cycle(int vw, int lvl) {
   int i,vcount,iter,excpt;  // DON'T MAKE THESE STATIC SCREWS UP RECURSION
   int state;
   
   for (vcount=0;vcount<vw;++vcount) {
   
      for(iter=0;iter<itercrsn;++iter)
         iterate(lvl);
      
      if (lvl == mglvls-1) return;
      
      rsdl(lvl);
      
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->mg_getfres(lvl+1,excpt);
         excpt += state;
      } while (state != block::stop);
          
      cycle(vw, lvl+1);

      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->mg_getcchng(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
         
      for(iter=0;iter<iterrfne;++iter)
         iterate(lvl);
   }

   return;
}

void blocks::go() {
   int i,step;
   char outname[100];
   clock_t cpu_time;

   clock();
   for(step = 1;step<=ntstep;++step) {
      tadvance();
      matchboundaries();

      for(i=0;i<ncycle;++i) {
         cycle(vw);
         *log << i << ' ';
         maxres();
         *log << '\n';
      }
      number_str(outname, "end", step, 3);
      output(outname,ftype::tecplot);
      restructure();
   }
   cpu_time = clock();
   *log << "that took " << cpu_time << " cpu time" << std::endl;
   
   return;
}

void blocks::tadvance() {
   int i,lvl,excpt;
   int state;
   
   for (lvl=0;lvl<ngrid;++lvl) {
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i)
            state &= blk[i]->tadvance(lvl,excpt);
         excpt += state;
      } while (state != block::stop);
   }
   
   matchboundaries();

   return;
}

void blocks::restructure() {
   int i,lvl,excpt;
   int state;
   
   matchboundaries();
   
   excpt = 0;
   do {
      state = block::stop;
      for(i=0;i<nblock;++i)
         state &= blk[i]->adapt(excpt);
      excpt += state;
   } while (state != block::stop);
   
   for(lvl=1;lvl<ngrid;++lvl) {
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i) {
            state &= blk[i]->reconnect(lvl,excpt);
         }
         excpt += state;
      } while (state != block::stop);
   }
            
   return;
}

