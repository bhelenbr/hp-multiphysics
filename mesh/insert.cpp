/*
 *  insert.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Aug 31 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "boundary.h"
#include <math.h>
#include <float.h>   
#include <assert.h>
#include <stdlib.h>

int mesh::insert(FLT x[ND]) {
   int n,tind,vnear,err;
   int ntdel, tdel[maxlst];
   int nsdel, sdel[maxlst];
   
   for(n=0;n<ND;++n)
      vrtx[nvrtx][n] = x[n];
   
   qtree.addpt(nvrtx);
   qtree.nearpt(nvrtx,vnear);

   /* FIND TRIANGLE CONTAINING POINT */      
   tind = findtri(x,vnear);
   if (tind < 0) {
      std::cerr << "couldn't find triangle for point: " << x[0] << ' ' << x[1] << " vnear: " << vnear << std::endl;
      std::cerr << "maxsrch: " << maxsrch << "vtri: " << vd[vnear].tri << std::endl;
      output("error",ftype::grid);
      exit(1);
   }    
   if (nvrtx >= maxvst) {
      *log << "need to use larger growth factor: too many vertices" << std::endl;
      exit(1);
   }
   err = insert(tind,nvrtx,0.0,ntdel,tdel,nsdel,sdel);
   nvrtx += 1 -err;
   
   return(err);
}

int mesh::insert(int tind, int vnum, FLT theta, int &ntdel, int *tdel, int &nsdel, int *sdel) {
   int i,j,tin,v0,dir;
   int sct;
   int nskeep, skeep[maxlst];
   int sind, sind1;
   
   if (tind < 0) {
      *log << "Warning: trying to insert point outd domain" << std::endl;
      return(1);
   }
   
   /* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
   ntdel = 0;
   tdel[ntdel++] = tind;
   i1wk[tind] = 0;
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = td[tin].tri[j];
         if (tind < 0) continue;
         if (i1wk[tind] == 0) continue;

         if (incircle(tind,vrtx[vnum]) > 0.0) {
            tdel[ntdel++] = tind;
            i1wk[tind] = 0;
            assert(ntdel < maxlst -1);
         }
      }
   }

   /* CREATE LIST OF SIDES TO BE DELETED AND BOUNDARY SIDES */
   nsdel = 0;
   nskeep = 0;
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = td[tin].tri[j];
         if (tind < 0 || i1wk[tind] != 0) {
            /* KEEP SIDE */
            skeep[nskeep++] = td[tin].side[j];
            assert(nskeep < maxlst);
            continue;
         }
         else {
            /* DELETE SIDE */
            if (td[tin].sign[j] > 0) {
               sind = td[tin].side[j];
               sdel[nsdel++] = sind;
               assert(nsdel < maxlst);
            }
         }
      }
   }
   
   assert(nskeep == nsdel+3);
   assert(nskeep == ntdel+2);
   
   /* RESET i1wk */
   for(i=0;i<ntdel;++i) {
      i1wk[tdel[i]] = -1;
   }
   
   /*	CHECK THAT WE AREN'T INSERTING POINT VERY CLOSE TO BOUNDARY */
   for(i=0;i<nskeep;++i) {
      sind = skeep[i];
      if(fabs(minangle(vnum, sd[sind].vrtx[0] , sd[sind].vrtx[1])) < theta*M_PI/180.0) {
         *log << "#Warning: inserting too close to boundary" << std::endl;
         return(1);
      }
   }
   
   /* ADD NEW INDICES TO DEL LIST */
   for(i=nsdel;i<nskeep;++i)
      sdel[i] = nside +(i-nsdel);
   
   for(i=ntdel;i<nskeep;++i)
      tdel[i] = ntri +(i-ntdel);
      
   ntri += 2;
   nside += 3;
   
   if (ntri > maxvst || nside > maxvst) {
      *log << "need to use bigger growth factor: too many sides/tris" << nside << ntri << std::endl;
      exit(1);
   }
   
   sct = 0;
   for(i=0;i<nskeep;++i) {
      sind = skeep[i];
      tind = tdel[i];
      /* CHECK SIDE DIRECTION  */
      if (area(sind,vnum) > 0.0) 
         dir = 0;
      else 
         dir = 1;
         
      /* CREATE NEW INFO */
      v0 = sd[sind].vrtx[dir];
      td[tind].vrtx[0] = v0;
      vd[v0].tri = tind;
      v0 = sd[sind].vrtx[1-dir];
      td[tind].vrtx[1] = v0;
      vd[v0].tri = tind;
      td[tind].vrtx[2] = vnum;
      vd[vnum].tri = tind;

      /* SIDE 2 INFO */      
      td[tind].side[2] = sind;
      td[tind].sign[2] = 1 -2*dir;
      
      sd[sind].tri[dir] = tind;
      tin = sd[sind].tri[1-dir];
      td[tind].tri[2] = tin;
      if (tin > -1) {
         for(j=0;j<3;++j) {
            if (td[tin].side[j] == sind) {
               td[tin].tri[j] = tind;
               break;
            }
         }
      }

      /* CREATE SIDE 1 INFO */
      v0 = sd[sind].vrtx[dir];
      sind1 = i1wk[v0];
      if (sind1 > -1) {
         td[tind].side[1] = sind1;
         td[tind].sign[1] = -1;
         
         sd[sind1].tri[1] = tind;
         tin = sd[sind1].tri[0];
         td[tind].tri[1] = tin;
         if (tin > -1) {
            for(j=0;j<3;++j) {
               if (td[tin].side[j] == sind1) {
                  td[tin].tri[j] = tind;
                  break;
               }
            }
         }
      }
      else {
         sind1 = sdel[sct++];
         td[tind].side[1] = sind1;
         td[tind].sign[1] = -1;
         
         sd[sind1].tri[1] = tind;
         sd[sind1].vrtx[0] = v0;
         sd[sind1].vrtx[1] = vnum;
         i1wk[v0] = sind1;
      }
      

      /* CREATE SIDE 0 INFO */
      v0 = sd[sind].vrtx[1-dir];
      sind1 = i1wk[v0];
      if (sind1 > -1) {
         td[tind].side[0] = sind1;
         td[tind].sign[0] = 1;
         
         sd[sind1].tri[0] = tind;
         tin = sd[sind1].tri[1];
         td[tind].tri[0] = tin;
         if (tin > -1) {
            for(j=0;j<3;++j) {
               if (td[tin].side[j] == sind1) {
                  td[tin].tri[j] = tind;
                  break;
               }
            }
         }
      }
      else {
         sind1 = sdel[sct++];
         td[tind].side[0] = sind1;
         td[tind].sign[0] = 1;
         
         sd[sind1].tri[0] = tind;
         
         sd[sind1].vrtx[0] = v0;
         sd[sind1].vrtx[1] = vnum;
         i1wk[v0] = sind1;
      }
   }
   
   /* RESET i1wk */
   for(i=0;i<nskeep;++i) {
      sind = skeep[i];
      i1wk[sd[sind].vrtx[0]] = -1;
      i1wk[sd[sind].vrtx[1]] = -1;
   }
      
   return(0);
}

void mesh::bdry_insert(int tind, int snum, int vnum, int &ntdel, int *tdel, int &nsdel, int *sdel) {
   int i,j,tin,v0,dir,tbdry;
   int nskeep, skeep[maxlst];
   int sct;
   int sind, sind1;
   
   tbdry = tind;
      
   /* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
   ntdel = 0;
   tdel[ntdel++] = tind;
     i1wk[tind] = 0;
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = td[tin].tri[j];
         if (tind < 0) continue;
         if (i1wk[tind] == 0) continue;

         if (incircle(tind,vrtx[vnum]) > 0.0) {
            tdel[ntdel++] = tind;
            i1wk[tind] = 0;
            assert(ntdel < maxlst);
         }
      }
   }

   /* CREATE LIST OF SIDES TO BE DELETED AND BOUNDARY SIDES */
   nsdel = 0;
   sdel[nsdel++] = td[tbdry].side[snum];
   nskeep = 0;
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = td[tin].tri[j];
         if (tind < 0) {
            if (tin == tbdry && j == snum) continue;               
            /* KEEP SIDE */
            skeep[nskeep++] = td[tin].side[j];
            assert(nskeep < maxlst);
            continue;
         }
         else if (i1wk[tind] != 0) {
            /* KEEP SIDE */
            skeep[nskeep++] = td[tin].side[j];
            assert(nskeep < maxlst);
            continue; 
         }
         else {
            /* DELETE SIDE */
            if (td[tin].sign[j] > 0) {
               sind = td[tin].side[j];
               sdel[nsdel++] = sind;
               assert(nsdel < maxlst);
            }
         }
      }
   }
   
   /* ALTER OLD BOUNDARY SIDE & CREATE NEW SIDE */
   sind = td[tbdry].side[snum];
   sd[nside].vrtx[0] = vnum;
   sd[nside].vrtx[1] = sd[sind].vrtx[1];
   sd[sind].vrtx[1] = vnum;
   /* ADD NEW SIDE TO BOUNDARY GROUP */
   /* NEED TO REORDER WHEN FINISHED */
   i = (-sd[sind].tri[1]>>16) -1;
   sd[nside].tri[1] = -(((i+1)<<16) +sbdry[i]->nel);
   sbdry[i]->el[sbdry[i]->nel++] = nside;
   if (sbdry[i]->nel > sbdry[i]->maxel) {
      *log << "too many boundary elements" <<  sbdry[i]->idnum  << ' ' << sbdry[i]->maxel << std::endl;
      exit(1);
   }
   ++nside;

   assert(nskeep == nsdel+1);
   assert(nskeep == ntdel+1);
   
   /* RESET i1wk */
   for(i=0;i<nsdel;++i)
      i1wk[tdel[i]] = -1;
   
   /* ADD NEW INDICES TO DEL LIST */
   for(i=nsdel;i<nskeep;++i)
      sdel[i] = nside +(i-nsdel);
      
   /* ADD THIS ONE TOO (NEVER GETS ACCESSED HERE)*/
   /* BUT FOR ACCOUNTING PURPOSES IN REBAY */
   sdel[nskeep] = nside -1;
   
   for(i=ntdel;i<nskeep;++i)
      tdel[i] = ntri +(i-ntdel);

   /* USE i1wk TO FIND SIDES */
   sind = td[tbdry].side[snum];
   i1wk[sd[sind].vrtx[0]] = sind;
   i1wk[sd[nside-1].vrtx[1]] = nside-1;
      
   ntri += 1;
   nside += 1;  
   
   if (ntri > maxvst || nside > maxvst) {
      *log << "need to use bigger growth factor: too many sides/tris:" << nside << ntri << std::endl;
      exit(1);
   }

   sct = 1;   // SKIP FIRST SIDE 
   for(i=0;i<nskeep;++i) {
      sind = skeep[i];
      tind = tdel[i];
      /* CHECK SIDE DIRECTION  */
      if (area(sind,vnum) > 0.0) 
         dir = 0;
      else 
         dir = 1;

      /* CREATE NEW INFO */
      v0 = sd[sind].vrtx[dir];
      td[tind].vrtx[0] = v0;
      vd[v0].tri = tind;
      v0 = sd[sind].vrtx[1-dir];
      td[tind].vrtx[1] = v0;
      vd[v0].tri = tind;
      td[tind].vrtx[2] = vnum;
      vd[vnum].tri = tind;

      /* SIDE 2 INFO */      
      td[tind].side[2] = sind;
      td[tind].sign[2] = 1 -2*dir;
      
      sd[sind].tri[dir] = tind;
      tin = sd[sind].tri[1-dir];
      td[tind].tri[2] = tin;
      if (tin > -1) {
         for(j=0;j<3;++j) {
            if (td[tin].side[j] == sind) {
               td[tin].tri[j] = tind;
               break;
            }
         }
      }

      /* CREATE SIDE 1 INFO */
      /* CREATE IN OPPOSITE DIRECTION */
      /* IF CREATED IT IS IN SAME DIRECTION */
      v0 = sd[sind].vrtx[dir];
      sind1 = i1wk[v0];
      if (sind1 > -1) {
         td[tind].side[1] = sind1;
         td[tind].sign[1] = 1;
         
         sd[sind1].tri[0] = tind;
         tin = sd[sind1].tri[1];
         td[tind].tri[1] = tin;
         if (tin > -1) {
            for(j=0;j<3;++j) {
               if (td[tin].side[j] == sind1) {
                  td[tin].tri[j] = tind;
                  break;
               }
            }
         }
      }
      else {
         sind1 = sdel[sct++];
         td[tind].side[1] = sind1;
         td[tind].sign[1] = -1;
         
         sd[sind1].tri[1] = tind;
         sd[sind1].vrtx[0] = v0;
         sd[sind1].vrtx[1] = vnum;
         i1wk[v0] = sind1;
      }
      

      /* CREATE SIDE 0 INFO */
      /* CREATE IN OPPOSITE DIRECTION */
      /* IF CREATED IT IS IN SAME DIRECTION */
      v0 = sd[sind].vrtx[1-dir];
      sind1 = i1wk[v0];
      if (sind1 > -1) {
         td[tind].side[0] = sind1;
         td[tind].sign[0] = 1;
         
         sd[sind1].tri[0] = tind;
         tin = sd[sind1].tri[1];
         td[tind].tri[0] = tin;
         if (tin > -1) {
            for(j=0;j<3;++j) {
               if (td[tin].side[j] == sind1) {
                  td[tin].tri[j] = tind;
                  break;
               }
            }
         }
      }
      else {
         sind1 = sdel[sct++];
         td[tind].side[0] = sind1;
         td[tind].sign[0] = -1;
         
         sd[sind1].tri[1] = tind;
         
         sd[sind1].vrtx[0] = vnum;
         sd[sind1].vrtx[1] = v0;
         i1wk[v0] = sind1;
      }
   }

   /* RESET i1wk */   
   for(i=0;i<nskeep;++i) {
      sind = skeep[i];
      i1wk[sd[sind].vrtx[0]] = -1;
      i1wk[sd[sind].vrtx[1]] = -1;
   } 
   
   sind = td[tbdry].side[snum];
   i1wk[sd[sind].vrtx[0]] = -1;
   i1wk[sd[nside-1].vrtx[1]] = -1;

   return;
}

int mesh::findtri(FLT x[ND], int vnear) const {
   int i,j,vn,dir,stoptri,tin,tind;
   int ntdel, tdel[maxlst];
   int tclose,nsurround;
   FLT minclosest,closest;
   
   /* HERE WE USE i1wk THIS MUST BE -1 BEFORE USING */
   tind = vd[vnear].tri;
   stoptri = tind;
   dir = 1;
   ntdel = 0;
   do {
      if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND;
      i1wk[tind] = 0;
      tdel[ntdel++] = tind;
      assert(ntdel < maxlst -1);
   
      for(vn=0;vn<3;++vn) 
         if (td[tind].vrtx[vn] == vnear) break;
      assert(vn != 3);
      
      tind = td[tind].tri[(vn +dir)%3];
      if (tind < 0) {
         if (dir > 1) break;
         /* REVERSE DIRECTION AND GO BACK TO START */
         ++dir;
         tind = vd[vnear].tri;
         for(vn=0;vn<3;++vn) 
            if (td[tind].vrtx[vn] == vnear) break;
         assert(vn != 3);
         tind = td[tind].tri[(vn +dir)%3];
         if (tind < 0) break;
         stoptri = -1;
      }
   } while(tind != stoptri); 
   
   nsurround = ntdel;
      
   /* DIDN'T FIND TRIANGLE */
   /* NEED TO SEARCH SURROUNDING TRIANGLES */
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = td[tin].tri[j];
         if (tind < 0) continue;
         if (i1wk[tind] == 0) continue;
         i1wk[tind] = 0;
         tdel[ntdel++] = tind;         
         if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND;
      }
      if (ntdel >= maxsrch-4) break;
   }
   std::cerr << "couldn't find tri for point " << x[0] << ' ' << x[1] << ' ' << vnear << std::endl;
   tind = tdel[0];
   minclosest = intri(tind,x)/area(tind);
   tclose = tind;
   for (i=1;i<nsurround;++i) {
      tind = tdel[i];
      if ((closest = intri(tind,x)/area(tind)) < minclosest) {
         minclosest = closest;
         tclose = tind;
      }
   }
   intri(tclose,x);
   tind = -tclose;
      
FOUND:
   /* RESET INTWKW1 TO -1 */
   for(i=0;i<ntdel;++i)
      i1wk[tdel[i]] = -1;

   return(tind);
}

