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
#include <utilities.h>
#ifdef MPISRC
#include <mpi.h>
#endif

#ifdef CAPRI
#include <capri.h>
#endif

std::ostream *sim::log = &std::cout;
double sim::time, sim::dti, sim::g;
sharedmem sim::scratch;

#ifdef BACKDIFF
FLT sim::bd[BACKDIFF+1];
#endif
#ifdef DIRK
FLT sim::bd[1];
#if (DIRK == 3)
/* THIS IS THE STANDARD FORM */
// FLT sim::adirk[DIRK][DIRK] = {{GRK3,0.0,0.0},{C2RK3-GRK3,GRK3,0.0},{1-B2RK3-GRK3,B2RK3,GRK3}} 
// FLT sim::bdirk[DIRK] = {1-B2RK3-GRK3,B2RK3,GRK3};
// FLT sim::cdirk[DIRK] = {GRK3,C2RK3,1.0};
/* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
FLT sim::adirk[DIRK][DIRK] = { {1./GRK3,0.0,0.0}, {C2RK3-GRK3,1./GRK3,0.0}, {1.-B2RK3-C2RK3,B2RK3,1./GRK3} };
FLT sim::cdirk[DIRK] = {GRK3,C2RK3-GRK3,1.0-C2RK3};
#else
/* THIS IS THE STANDARD FORM */
// FLT sim::adirk[DIRK][DIRK] = {{0.0,0.0,0.0,0.0},{GRK4,GRK4,0.0,0.0},{C3RK4-A32RK4-GRK4,A32RK4,GRK4,0.0},{B1RK4,B2RK4,B3RK4,GRK4}} 
// FLT sim::bdirk[DIRK] = {B1RK4,B2RK4,B3RK4,GRK4};
// FLT sim::cdirk[DIRK] = {0.0,2.*GRK4,C3RK4,1.0};
/* THIS IS THE INCREMENTAL FORM WITH DIAGONAL TERM IS INVERTED */
FLT sim::adirk[DIRK][DIRK] = {{1./GRK4,0.0,0.0,0.0},{GRK4,1./GRK4,0.0,0.0},{C3RK4-A32RK4-2.*GRK4,A32RK4,1./GRK4,0.0},{B1RK4-(C3RK4-A32RK4-GRK4),B2RK4-A32RK4,B3RK4,1./GRK4}}; 
FLT sim::cdirk[DIRK] = {2.*GRK4,C3RK4-2.*GRK4,1.0-C3RK4,0.0};
#endif
#endif

blitz::Array<FLT,1> sim::cfl;  

#ifdef PV3
int sim::pv3_mesh_changed;
#endif
   
FLT sim::fadd;  

void blocks::init(const char *infile, const char *outfile) {
   std::map<std::string,std::string> maptemp;
   char name[100];
   
   strcpy(name,infile);
   strcat(name,".inpt");
   input_map(maptemp,name);
   
   if (outfile) {
      maptemp["logfile"] = outfile;
   }
   
   init(maptemp);
   
   return;
}
   


void blocks::init(std::map<std::string,std::string> input) {
   int i,nb;
   char outfile[100];

   std::map<std::string,std::string>::const_iterator mi;
   std::map<std::string,std::string> merge;
   std::istringstream data;
   
#ifdef MPISRC
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif

   /* OPEN LOGFILES FOR EACH PROCESSOR */
   mi = input.find("logfile");
   if (mi != input.end()) {
      data.str(mi->second);
      data >> outfile;
      data.clear();
      strcat(outfile,".");
      number_str(outfile,outfile,myid,1);
      strcat(outfile,".log");
      std::ofstream *filelog = new std::ofstream;
      filelog->setf(std::ios::scientific, std::ios::floatfield);
      filelog->precision(3);
      filelog->open(outfile);
      sim::log = filelog;
   }
   else {
      std::cout.setf(std::ios::scientific, std::ios::floatfield);
      std::cout.precision(3);
      sim::log = &std::cout;
   }
#ifdef CAPRI
   int status = gi_uStart();
   *sim::log << "gi_uStart status = ", status, "\n";
   if (status != CAPRI_SUCCESS) exit(1);
   
   mi = input.find("BRep");
   if (mi != input.end()) {
      char temp[100];
      strcpy(temp,mi->second.c_str());
      status = gi_uLoadPart(temp);
      *sim::log << mi->second << ": gi_uLoadPart status = " << status << std::endl;
      if (status != CAPRI_SUCCESS) exit(1);
   }
#endif 

   /* EXTRACT NBLOCKS FOR MYID */
   /* SPACE DELIMITED ARRAY OF NBLOCKS FOR EACH PROCESSOR */
   int bstart = 0;
   data.str(input["nblock"]);
   for (i=0;i<myid;++i) {
      if (!(data >> nb)) {
         std::cerr << "error reading blocks\n"; 
         exit(1);
      }
      bstart += nb;
   }
   if (!(data >> nblock)) {
      std::cerr << "error reading blocks\n"; 
      exit(1);
   }
   *sim::log << "#starting block index: " << bstart << std::endl;
   *sim::log << "#nblocks for this processor: " << nblock << std::endl;
   
   /* LOAD TIME STEPPING INFO */
   sim::time = 0.0;  // Simulation starts at t = 0
   data.str(input["dti"]);   
   if (!(data >> sim::dti)) sim::dti=1.0;
   *sim::log << "#dti: " << sim::dti << std::endl;
   data.clear();
   
   /* OTHER UNIVERSAL CONSTANTS */
   data.str(input["gravity"]);   
   if (!(data >> sim::g)) sim::g = 0.0;
   *sim::log << "#gravity: " << sim::g << std::endl;
   data.clear();
   
   data.str(input["iterrfne"]);   
   if (!(data >> iterrfne)) iterrfne = 0;
   *sim::log << "#iterrfne: " << iterrfne << std::endl;
   data.clear();
   
   /* LOAD BASIC CONSTANTS FOR MULTIGRID */
   data.str(input["itercrsn"]);   
   if (!(data >> itercrsn)) itercrsn = 1;
   *sim::log << "#itercrsn: " << itercrsn << std::endl;
   data.clear();
   
   data.str(input["iterrfne"]);   
   if (!(data >> iterrfne)) iterrfne = 0;
   *sim::log << "#iterrfne: " << iterrfne << std::endl;
   data.clear();
      
   data.str(input["njacobi"]);   
   if (!(data >> njacobi)) njacobi = 1;
   *sim::log << "#njacobi: " << njacobi << std::endl;
   data.clear();
   
   data.str(input["ncycle"]);   
   if (!(data >> ncycle)) ncycle = 20;
   *sim::log << "#ncycle: " << ncycle << std::endl;
   data.clear();
   
   data.str(input["vwcycle"]);   
   if (!(data >> vw)) vw = 2;
   *sim::log << "#vwcycle: " << vw << std::endl;
   data.clear();
   
   data.str(input["ntstep"]);   
   if (!(data >> ntstep)) ntstep = 1;
   *sim::log << "#ntstep: " << ntstep << std::endl;
   data.clear();
   
   /* LOAD NUMBER OF GRIDS */
   data.str(input["nblock"]);
   if (!(data >> nblock)) nblock = 1;
   *sim::log << "#nblock: " << nblock << std::endl;
   data.clear();
   
   data.str(input["ngrid"]);   
   if (!(data >> ngrid)) ngrid = 1;
   *sim::log << "#ngrid: " << ngrid << std::endl;
   data.clear();
   
   data.str(input["mglvls"]);   
   if (!(data >> mglvls)) mglvls = 1;
   *sim::log << "#mglvls: " << mglvls << std::endl;
   data.clear();

   blk = new block *[nblock];

   int lvl = 0;
   int excpt = 0;
   int state = block::stop;
   for (i=0;i<nblock;++i) {
      blk[i] = getnewblock(bstart+i,&input);
      blk[i]->init(input);
   }
   
   for (i=0;i<nblock;++i)
      blk[i]->reload_scratch_pointers();
   
   findmatch();
   matchboundaries();  // ONLY DOES FIRST LEVEL
   output("matched",ftype::grid);

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
   
   for(i=0;i<nblock;++i) {
      number_str(outfile,"coarsenchk",i,1);
      strcat(outfile,".");
      blk[i]->coarsenchk(outfile);
   }
      
   return;
}

static int maxblock,maxbdry,pdig,bdig,ndig;

int tagid(int vsf,int p1,int p2,int b1,int b2,int n1,int n2) {
   int tag,lefttag,righttag;
   
   lefttag = (p1<<(bdig+ndig)) +(b1<<(ndig)) +n1;
   righttag = (p2<<(bdig+ndig)) +(b2<<(ndig)) +n2;

   if (lefttag > righttag) {
      tag = (vsf<<(2*(bdig+ndig+pdig))) +(lefttag<<(bdig+ndig+pdig)) +righttag;
   }
   else {
      tag = (vsf<<(2*(bdig+ndig+pdig)))+(righttag<<(bdig+ndig+pdig)) +lefttag;
   }
   
   return(tag);
}

void blocks::commblocks::output() {
   *sim::log << "# number of blocks: " << nblock << std::endl;
   for (int i=0;i<nblock;++i) {
      *sim::log << "#\tblock: " << i << std::endl;
      
      *sim::log << "#\t\tnvcomm: " << blkinfo[i].nvcomm << std::endl;
      for (int v=0;v<blkinfo[i].nvcomm;++v)
         *sim::log << "#\t\t\tnvbd: " << blkinfo[i].vcomm[v].nvbd << " idnum: " << blkinfo[i].vcomm[v].idnum << std::endl;
      
      *sim::log << "#\t\tnscomm: " << blkinfo[i].nscomm << std::endl;
      for (int v=0;v<blkinfo[i].nscomm;++v)
         *sim::log << "#\t\t\tnsbd: " << blkinfo[i].scomm[v].nsbd << " idnum: " << blkinfo[i].scomm[v].idnum << std::endl;
      
      *sim::log << "#\t\tnfcomm: " << blkinfo[i].nfcomm << std::endl;
      for (int v=0;v<blkinfo[i].nfcomm;++v)
         *sim::log << "#\t\t\tnfbd: " << blkinfo[i].fcomm[v].nfbd << " idnum: " << blkinfo[i].fcomm[v].idnum << std::endl;
   }
}

int blocks::commblocks::unpack(blitz::Array<int,1> entitylist) {
   int count = 0;
   nblock = entitylist(count++);
   blkinfo = new commblock[nblock];
   for(int b=0;b<nblock;++b) {
      blkinfo[b].nvcomm = entitylist(count++);
      blkinfo[b].vcomm = new commblock::vid[blkinfo[b].nvcomm];
      for(int i=0;i<blkinfo[b].nvcomm;++i) {
         blkinfo[b].vcomm[i].nvbd = entitylist(count++);
         blkinfo[b].vcomm[i].idnum = entitylist(count++);
      }
      
      blkinfo[b].nscomm = entitylist(count++);
      blkinfo[b].scomm = new commblock::sid[blkinfo[b].nscomm];
      for(int i=0;i<blkinfo[b].nscomm;++i) {
         blkinfo[b].scomm[i].nsbd = entitylist(count++);
         blkinfo[b].scomm[i].idnum = entitylist(count++);
         // blkinfo[b].scomm[i].vid0 = entitylist(count++);
         // blkinfo[b].scomm[i].vid1 = entitylist(count++);
      }
      
      blkinfo[b].nfcomm = entitylist(count++);
      blkinfo[b].fcomm = new commblock::fid[blkinfo[b].nfcomm];
      for(int i=0;i<blkinfo[b].nfcomm;++i) {
         blkinfo[b].fcomm[i].nfbd = entitylist(count++);
         blkinfo[b].fcomm[i].idnum = entitylist(count++);
         // bp->fcomm[i].nsds = entitylist(count++);
         // bp->fcomm[i].sids = new int[bp->fcomm[i].nsds];
         // for (j=0;j<bp->fcomm[i].nsds;++j)
            // bp->fcomm[i].sids[j] = entitylist(count++);
      }
   }
   return(count);
}


void blocks::findmatch() {
   int b,b1,b2,i,j,grdlvl,p,count,tsize;
   Array<int,1> entitylist, sndentitylist, sublist;
   int *size, *sndsize;
   commblocks *blksinfo;
   commblocks::commblock *bp1, *bp2;
   
   /* FIRST DETERMINE TOTAL SIZE OF LIST */
   for(grdlvl=0;grdlvl<ngrid;++grdlvl) {
      blksinfo = new commblocks[nproc];

      *sim::log << "# finding matches at multigrid level " << grdlvl << " for processor " << myid << std::endl;
      
      sndsize = new int[nproc];
      for(i=0;i<nproc;++i)
         sndsize[i] = 0;
         
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
         
      sndentitylist.resize(tsize);
      sndentitylist = 0;

      /* BEGIN ASSEMBLING COMMUNICATION GRAPH LIST*/
      count = 0;
      for(p=0;p<myid;++p)
         count += size[p];
         
      sndentitylist(count++) = nblock;
      for(b=0;b<nblock;++b) {
         sublist.reference(sndentitylist(Range(count,toEnd)));
         count += blk[b]->comm_entity_list(grdlvl,sublist);
      }
      ~sublist;

#ifdef MPISRC   
      entitylist.resize(tsize);
      entitylist = 0;
        
      MPI_Allreduce(sndentitylist.data(),entitylist.data(),tsize,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#else
      entitylist.reference(sndentitylist);
#endif
      ~sndentitylist;

      /* UNPACK ENTITY LIST INTO MORE USABLE FORM */
      maxblock = 0; // SO I CAN GENERATE UNIQUE TAGS
      maxbdry = 0; // SO I CAN GENERATE UNIQUE TAGS
      count = 0;
      for(p=0;p<nproc;++p) {
         count += blksinfo[p].unpack(entitylist(Range(count,toEnd)));
         maxblock = MAX(maxblock,blksinfo[p].nblock);
         for(int b=0;b<blksinfo[p].nblock;++b) {
            maxbdry = MAX(maxbdry,blksinfo[p].blkinfo[b].nvcomm);
            maxbdry = MAX(maxbdry,blksinfo[p].blkinfo[b].nscomm);
            maxbdry = MAX(maxbdry,blksinfo[p].blkinfo[b].nfcomm);
         }
         /* OUTPUT LIST FOR DEBUGGING */
         *sim::log << "# processor " << p << " block data" << std::endl;
         blksinfo[p].output();
      }
      ~entitylist;
      
      
      /* CALCULATE NUMBER OF BINARY DIGITS NEEDED TO MAKE TAGS */
      pdig = 1;
      while ((nproc>>pdig) > 0) ++pdig;
      bdig = 1;
      while ((maxblock>>bdig) > 0) ++bdig;
      ndig = 1;
      while ((maxbdry>>ndig) > 0) ++ndig;
      if (2*(pdig+bdig+ndig)+2 > 32) *sim::log << "can't guarantee unique tags\n"; 

      /* CAN NOW START TO LOOK FOR MATCHES */
      /* FOR NOW ALL COMMUNICATION IN 1 PHASE */
      /* LATER DO SOMETHING FANCY WITH ALL THIS INFO ? */
      
      for(b1=0;b1<blksinfo[myid].nblock;++b1) {
         bp1 = &blksinfo[myid].blkinfo[b1];
         *sim::log << "# finding matches for block " << b1 << std::endl;

         /* LOOK FOR VERTEX MATCHES */
         for(i=0;i<bp1->nvcomm;++i) {
            bool first_found = false;
            *sim::log << "#\tvertex " << bp1->vcomm[i].nvbd << " idnum: " << bp1->vcomm[i].idnum << std::endl;
            for (p=0;p<nproc;++p) {
               for(b2=0;b2<blksinfo[p].nblock;++b2) {
                  bp2 = &blksinfo[p].blkinfo[b2];
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
                           *sim::log <<  "#\t\tlocal match to processor " << p << " block: " << b2 << " vrtx: " << bp2->vcomm[j].nvbd << " tag: " << tagid(1,myid,myid,b1,b2,i,j) << std::endl;
                        }
#ifdef MPISRC
                        else {
                           boundary *v1 = blk[b1]->vbdry(grdlvl,bp1->vcomm[i].nvbd);
                           v1->mpi_cnnct(p,tagid(1,myid,p,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false;
                              first_found = true;
                           }
                           *sim::log <<  "#\t\tmpi match to processor " << p << " block: " << b2 << " vrtx: " << bp2->vcomm[j].nvbd << " tag: " << tagid(1,myid,p,b1,b2,i,j) << std::endl;
                        }
#endif
                        
                     }
                  }
               }
            }
         }
               
         /* LOOK FOR SIDE MATCHES */
         for(i=0;i<bp1->nscomm;++i) {
            bool first_found = false;
            *sim::log << "#\tside " << bp1->scomm[i].nsbd << " idnum: " << bp1->scomm[i].idnum << std::endl;
            for (p=0;p<nproc;++p) {
               for(b2=0;b2<blksinfo[p].nblock;++b2) {
                  bp2 = &blksinfo[p].blkinfo[b2];
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
                           *sim::log << "#\t\tlocal match to processor " << p << " block: " << b2 << " side: " << bp2->scomm[j].nsbd << " tag: " << tagid(2,myid,myid,b1,b2,i,j) << std::endl;
                        }
#ifdef MPISRC
                        else {
                           boundary *v1 = blk[b1]->sbdry(grdlvl,bp1->scomm[i].nsbd);
                           v1->mpi_cnnct(p,tagid(2,myid,p,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false;
                              first_found = true;
                           }
                           *sim::log << "#\t\tmpi match to processor " << p << " block: " << b2 << " side: " << bp2->scomm[j].nsbd << " tag: " << tagid(2,myid,p,b1,b2,i,j) << std::endl;
                        }
#endif
                     }
                  }
               }
            }
         }
         
         /* LOOK FOR FACE MATCHES */
         for(i=0;i<bp1->nfcomm;++i) {
            bool first_found = false;
            *sim::log << "#\tface " << bp1->fcomm[i].nfbd << " idnum: " << bp1->fcomm[i].idnum << std::endl;
            for (p=0;p<nproc;++p) {
               for(b2=0;b2<blksinfo[p].nblock;++b2) {
                  bp2 = &blksinfo[p].blkinfo[b2];
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
                           *sim::log <<  "#\t\tlocal match to processor " << p << " block: " << b2 << " face: " << bp2->fcomm[j].nfbd << " tag: " << tagid(3,myid,myid,b1,b2,i,j) << std::endl;

                        }
#ifdef MPISRC
                        else {
                           boundary *v1 = blk[b1]->fbdry(grdlvl,bp1->fcomm[i].nfbd);
                           v1->mpi_cnnct(p,tagid(1,myid,p,b1,b2,i,j));
                           if (!first_found) {
                              v1->is_frst() = false;
                              first_found = true;
                           }
                           *sim::log << "#\t\t mpi match to processor " << p << " block: " << b2 << " face: " << bp2->fcomm[j].nfbd << " tag: " << tagid(3,myid,p,b1,b2,i,j) << std::endl;

                        }
#endif
                     }
                  }
               }
            }
         }
      }
               
      /* DELETE DATA STRUCTURE */
      delete []blksinfo;
   }
   delete []size;

   return;
}

            
            
   /* MATCH BOUNDARIES */
void blocks::matchboundaries() {
   int i,lvl,excpt;
   int state;
   
   for (lvl=0;lvl<1;++lvl) {  // TEMPORARY
      excpt = 0;
      do {
         state = block::stop;
         for(i=0;i<nblock;++i) {
            state &= blk[i]->matchboundaries(lvl,excpt);
         }
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
      tadvance(step);
      matchboundaries();

      for(i=0;i<ncycle;++i) {
         cycle(vw);
         *sim::log << i << ' ';
         maxres();
         *sim::log << '\n';
      }
      number_str(outname, "end", step, 3);
      output(outname,ftype::grid);
      restructure();
   }
   cpu_time = clock();
   *sim::log << "that took " << cpu_time << " cpu time" << std::endl;
   
   return;
}

void blocks::tadvance(int step) {
   int i,lvl,excpt;
   int state;
   
   if (sim::dti > 0.0) sim::time += 1./sim::dti;

#ifdef BACKDIFF
   for(i=0;i<BACKDIFF+1;++i)
      sim::bd[i] = 0.0;
   
   switch(step) {
      case(1):
         sim::bd[0] =  sim::dti;
         sim::bd[1] = -sim::dti;
         break;
#if (BACKDIFF > 1)
      case(2):
         sim::bd[0] =  1.5*sim::dti;
         sim::bd[1] = -2.0*sim::dti;
         sim::bd[2] =  0.5*sim::dti;
         break;
#endif
#if (BACKDIFF > 2)
      case(3):
         sim::bd[0] = 11./6*sim::dti;
         sim::bd[1] = -3.*sim::dti;
         sim::bd[2] = 1.5*sim::dti;
         sim::bd[3] = -1./3.*sim::dti;
         break;
#endif
   }
#endif

#if (DIRK == 4)
   /* STARTUP SEQUENCE */
   if (step == 1) {
      sim::adirk[0][0] = 0.0; sim::adirk[0][1] = 0.0;            sim::adirk[0][2] = 0.0;     sim::adirk[0][3] = 0.0;
      sim::adirk[1][0] = 0.0; sim::adirk[1][1] = 1./sim::GRK3;   sim::adirk[1][2] = 0.0;     sim::adirk[1][3] = 0.0;
      sim::adirk[2][0] = 0.0; sim::adirk[2][1] = sim::C2RK3-sim::GRK3;     sim::adirk[2][2] = 1./sim::GRK3; sim::adirk[2][3] = 0.0;
      sim::adirk[3][0] = 0.0; sim::adirk[3][1] = 1.-sim::B2RK3-sim::C2RK3; sim::adirk[3][2] = sim::B2RK3;   sim::adirk[3][3] = 1./sim::GRK3;
      sim::cdirk[0] = sim::GRK3; sim::cdirk[1] = sim::C2RK3-sim::GRK3; sim::cdirk[2] = 1.0-sim::C2RK3;
   } 
   else if (step == 2) {
      sim::adirk[0][0] = 1./sim::GRK3;                   sim::adirk[0][1] = 0.0;          sim::adirk[0][2] = 0.0;     sim::adirk[0][3] = 0.0;
      sim::adirk[1][0] = sim::GRK4;                      sim::adirk[1][1] = 1./sim::GRK4;      sim::adirk[1][2] = 0.0;     sim::adirk[1][3] = 0.0;
      sim::adirk[2][0] = sim::C3RK4-sim::A32RK4-2.*sim::GRK4;      sim::adirk[2][1] = sim::A32RK4;       sim::adirk[2][2] = 1./sim::GRK4; sim::adirk[2][3] = 0.0;
      sim::adirk[3][0] = sim::B1RK4-(sim::C3RK4-sim::A32RK4-sim::GRK4); sim::adirk[3][1] = sim::B2RK4-sim::A32RK4; sim::adirk[3][2] = sim::B3RK4;   sim::adirk[3][3] = 1./sim::GRK4; 
      sim::cdirk[0] = 2.*sim::GRK4; sim::cdirk[1] = sim::C3RK4-2.*sim::GRK4; sim::cdirk[2] = 1.0-sim::C3RK4;
   }
   else {
      sim::adirk[0][0] = 1./sim::GRK4;                   sim::adirk[0][1] = 0.0;          sim::adirk[0][2] = 0.0;     sim::adirk[0][3] = 0.0;
      sim::adirk[1][0] = sim::GRK4;                      sim::adirk[1][1] = 1./sim::GRK4;      sim::adirk[1][2] = 0.0;     sim::adirk[1][3] = 0.0;
      sim::adirk[2][0] = sim::C3RK4-sim::A32RK4-2.*sim::GRK4;      sim::adirk[2][1] = sim::A32RK4;       sim::adirk[2][2] = 1./sim::GRK4; sim::adirk[2][3] = 0.0;
      sim::adirk[3][0] = sim::B1RK4-(sim::C3RK4-sim::A32RK4-sim::GRK4); sim::adirk[3][1] = sim::B2RK4-sim::A32RK4; sim::adirk[3][2] = sim::B3RK4;   sim::adirk[3][3] = 1./sim::GRK4; 
      sim::cdirk[0] = 2.*sim::GRK4; sim::cdirk[1] = sim::C3RK4-2.*sim::GRK4; sim::cdirk[2] = 1.0-sim::C3RK4;
   }
#endif

   
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
   char outfile[120];
   
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
   
   for(i=0;i<nblock;++i) {
      number_str(outfile,"coarsenchk",i,1);
      strcat(outfile,".");
      blk[i]->coarsenchk(outfile);
   }
            
   return;
}

