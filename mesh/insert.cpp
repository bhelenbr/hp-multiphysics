/*
 *  insert.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Aug 31 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include"mesh.h"
#include<math.h>
#include<float.h>   
#include<assert.h>
#include<stdlib.h>

int ntdel, tdel[MAXLST+1];
int nsdel, sdel[MAXLST+1];

int mesh::insert(FLT x, FLT y) {
   int tind,vnear,err;
   
   vrtx[nvrtx][0] = x;
   vrtx[nvrtx][1] = y;
   
   qtree.addpt(nvrtx);
   qtree.nearpt(nvrtx,vnear);

   /* FIND TRIANGLE CONTAINING POINT */      
   tind = findtri(x,y,vnear);
   assert(tind > -1);  
   if (nvrtx >= maxvst) {
      printf("need to use larger growth factor: too many vertices\n");
      exit(1);
   }
   err = insert(tind,nvrtx);
   nvrtx += 1 -err;
   
   return(err);
}

int mesh::insert(int tind, int vnum, FLT theta = 0.0) {
   int i,j,tin,v0,dir;
   int sct, nskeep, skeep[MAXLST+1];
   int sind, sind1;
   
   if (tind < 0) {
      printf("Warning: trying to insert point outside domain\n");
      return(1);
   }
   
   /* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
   ntdel = 0;
   tdel[ntdel++] = tind;
   intwk1[tind] = 0;
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = ttri[tin][j];
         if (tind < 0) continue;
         if (intwk1[tind] == 0) continue;

         if (incircle(tind,vrtx[vnum]) > 0.0) {
            tdel[ntdel++] = tind;
            intwk1[tind] = 0;
            assert(ntdel < MAXLST);
         }
      }
   }

   /* CREATE LIST OF SIDES TO BE DELETED AND BOUNDARY SIDES */
   nsdel = 0;
   nskeep = 0;
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = ttri[tin][j];
         if (tind < 0 || intwk1[tind] != 0) {
            /* KEEP SIDE */
            skeep[nskeep++] = tside[tin].side[j];
            assert(nskeep < MAXLST);
            continue;
         }
         else {
            /* DELETE SIDE */
            if (tside[tin].sign[j] > 0) {
               sind = tside[tin].side[j];
               sdel[nsdel++] = sind;
               assert(nsdel < MAXLST);
            }
         }
      }
   }
   
   assert(nskeep == nsdel+3);
   assert(nskeep == ntdel+2);
   
   /* RESET INTWK1 */
   for(i=0;i<nskeep;++i)
      intwk1[tdel[i]] = -1;
      
   /*	CHECK THAT WE AREN'T INSERTING POINT VERY CLOSE TO BOUNDARY */
   for(i=0;i<nskeep;++i) {
      sind = skeep[i];
      if(fabs(minangle(vnum, svrtx[sind][0] , svrtx[sind][1])) < theta*M_PI/180.0) {
         printf("#Warning: inserting too close to boundary\n");
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
      printf("need to use bigger growth factor: too many sides/tris %d/%d\n",nside,ntri);
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
      v0 = svrtx[sind][dir];
      tvrtx[tind][0] = v0;
      vtri[v0] = tind;
      v0 = svrtx[sind][1-dir];
      tvrtx[tind][1] = v0;
      vtri[v0] = tind;
      tvrtx[tind][2] = vnum;
      vtri[vnum] = tind;

      /* SIDE 2 INFO */      
      tside[tind].side[2] = sind;
      tside[tind].sign[2] = 1 -2*dir;
      
      stri[sind][dir] = tind;
      tin = stri[sind][1-dir];
      ttri[tind][2] = tin;
      if (tin > -1) {
         for(j=0;j<3;++j) {
            if (tside[tin].side[j] == sind) {
               ttri[tin][j] = tind;
               break;
            }
         }
      }

      /* CREATE SIDE 1 INFO */
      v0 = svrtx[sind][dir];
      sind1 = intwk1[v0];
      if (sind1 > -1) {
         tside[tind].side[1] = sind1;
         tside[tind].sign[1] = -1;
         
         stri[sind1][1] = tind;
         tin = stri[sind1][0];
         ttri[tind][1] = tin;
         if (tin > -1) {
            for(j=0;j<3;++j) {
               if (tside[tin].side[j] == sind1) {
                  ttri[tin][j] = tind;
                  break;
               }
            }
         }
      }
      else {
         sind1 = sdel[sct++];
         tside[tind].side[1] = sind1;
         tside[tind].sign[1] = -1;
         
         stri[sind1][1] = tind;
         svrtx[sind1][0] = v0;
         svrtx[sind1][1] = vnum;
         intwk1[v0] = sind1;
      }
      

      /* CREATE SIDE 0 INFO */
      v0 = svrtx[sind][1-dir];
      sind1 = intwk1[v0];
      if (sind1 > -1) {
         tside[tind].side[0] = sind1;
         tside[tind].sign[0] = 1;
         
         stri[sind1][0] = tind;
         tin = stri[sind1][1];
         ttri[tind][0] = tin;
         if (tin > -1) {
            for(j=0;j<3;++j) {
               if (tside[tin].side[j] == sind1) {
                  ttri[tin][j] = tind;
                  break;
               }
            }
         }
      }
      else {
         sind1 = sdel[sct++];
         tside[tind].side[0] = sind1;
         tside[tind].sign[0] = 1;
         
         stri[sind1][0] = tind;
         
         svrtx[sind1][0] = v0;
         svrtx[sind1][1] = vnum;
         intwk1[v0] = sind1;
      }
   }
   
   /* RESET INTWK1 */
   for(i=0;i<nskeep;++i) {
      sind = skeep[i];
      intwk1[svrtx[sind][0]] = -1;
      intwk1[svrtx[sind][1]] = -1;
   }
      
   return(0);
}

void mesh::bdry_insert(int tind, int snum, int vnum) {
   int i,j,tin,v0,dir,tbdry;
   int sct, nskeep, skeep[MAXLST+1];
   int sind, sind1;
   
   tbdry = tind;
      
   /* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
   ntdel = 0;
   tdel[ntdel++] = tind;
     intwk1[tind] = 0;
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = ttri[tin][j];
         if (tind < 0) continue;
         if (intwk1[tind] == 0) continue;

         if (incircle(tind,vrtx[vnum]) > 0.0) {
            tdel[ntdel++] = tind;
            intwk1[tind] = 0;
            assert(ntdel < MAXLST);
         }
      }
   }

   /* CREATE LIST OF SIDES TO BE DELETED AND BOUNDARY SIDES */
   nsdel = 0;
   sdel[nsdel++] = tside[tbdry].side[snum];
   nskeep = 0;
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = ttri[tin][j];
         if (tind < 0) {
            if (tin == tbdry && j == snum) continue;               
            /* KEEP SIDE */
            skeep[nskeep++] = tside[tin].side[j];
            assert(nskeep < MAXLST);
            continue;
         }
         else if (intwk1[tind] != 0) {
            /* KEEP SIDE */
            skeep[nskeep++] = tside[tin].side[j];
            assert(nskeep < MAXLST);
            continue; 
         }
         else {
            /* DELETE SIDE */
            if (tside[tin].sign[j] > 0) {
               sind = tside[tin].side[j];
               sdel[nsdel++] = sind;
               assert(nsdel < MAXLST);
            }
         }
      }
   }
   
   /* ALTER OLD BOUNDARY SIDE & CREATE NEW SIDE */
   sind = tside[tbdry].side[snum];
   svrtx[nside][0] = vnum;
   svrtx[nside][1] = svrtx[sind][1];
   svrtx[sind][1] = vnum;
   /* ADD NEW SIDE TO BOUNDARY GROUP */
   /* NEED TO REORDER WHEN FINISHED */
   i = (-stri[sind][1]>>16) -1;
   stri[nside][1] = -(((i+1)<<16) +sbdry[i].num);
   sbdry[i].el[sbdry[i].num++] = nside;
   if (sbdry[i].num > maxsbel) {
      printf("too many boundary elements\n");
      exit(1);
   }
   ++nside;

   assert(nskeep == nsdel+1);
   assert(nskeep == ntdel+1);
   
   /* RESET INTWK1 */
   for(i=0;i<nskeep;++i)
      intwk1[tdel[i]] = -1;
   
   /* ADD NEW INDICES TO DEL LIST */
   for(i=nsdel;i<nskeep;++i)
      sdel[i] = nside +(i-nsdel);
      
   /* ADD THIS ONE TOO (NEVER GETS ACCESSED HERE)*/
   /* BUT FOR ACCOUNTING PURPOSES IN REBAY */
   sdel[nskeep] = nside -1;
   
   for(i=ntdel;i<nskeep;++i)
      tdel[i] = ntri +(i-ntdel);

   /* USE INTWK1 TO FIND SIDES */
   sind = tside[tbdry].side[snum];
   intwk1[svrtx[sind][0]] = sind;
   intwk1[svrtx[nside-1][1]] = nside-1;
      
   ntri += 1;
   nside += 1;  
   
   if (ntri > maxvst || nside > maxvst) {
      printf("need to use bigger growth factor: too many sides/tris %d/%d\n",nside,ntri);
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
      v0 = svrtx[sind][dir];
      tvrtx[tind][0] = v0;
      vtri[v0] = tind;
      v0 = svrtx[sind][1-dir];
      tvrtx[tind][1] = v0;
      vtri[v0] = tind;
      tvrtx[tind][2] = vnum;
      vtri[vnum] = tind;

      /* SIDE 2 INFO */      
      tside[tind].side[2] = sind;
      tside[tind].sign[2] = 1 -2*dir;
      
      stri[sind][dir] = tind;
      tin = stri[sind][1-dir];
      ttri[tind][2] = tin;
      if (tin > -1) {
         for(j=0;j<3;++j) {
            if (tside[tin].side[j] == sind) {
               ttri[tin][j] = tind;
               break;
            }
         }
      }

      /* CREATE SIDE 1 INFO */
      /* CREATE IN OPPOSITE DIRECTION */
      /* IF CREATED IT IS IN SAME DIRECTION */
      v0 = svrtx[sind][dir];
      sind1 = intwk1[v0];
      if (sind1 > -1) {
         tside[tind].side[1] = sind1;
         tside[tind].sign[1] = 1;
         
         stri[sind1][0] = tind;
         tin = stri[sind1][1];
         ttri[tind][1] = tin;
         if (tin > -1) {
            for(j=0;j<3;++j) {
               if (tside[tin].side[j] == sind1) {
                  ttri[tin][j] = tind;
                  break;
               }
            }
         }
      }
      else {
         sind1 = sdel[sct++];
         tside[tind].side[1] = sind1;
         tside[tind].sign[1] = -1;
         
         stri[sind1][1] = tind;
         svrtx[sind1][0] = v0;
         svrtx[sind1][1] = vnum;
         intwk1[v0] = sind1;
      }
      

      /* CREATE SIDE 0 INFO */
      /* CREATE IN OPPOSITE DIRECTION */
      /* IF CREATED IT IS IN SAME DIRECTION */
      v0 = svrtx[sind][1-dir];
      sind1 = intwk1[v0];
      if (sind1 > -1) {
         tside[tind].side[0] = sind1;
         tside[tind].sign[0] = 1;
         
         stri[sind1][0] = tind;
         tin = stri[sind1][1];
         ttri[tind][0] = tin;
         if (tin > -1) {
            for(j=0;j<3;++j) {
               if (tside[tin].side[j] == sind1) {
                  ttri[tin][j] = tind;
                  break;
               }
            }
         }
      }
      else {
         sind1 = sdel[sct++];
         tside[tind].side[0] = sind1;
         tside[tind].sign[0] = -1;
         
         stri[sind1][1] = tind;
         
         svrtx[sind1][0] = vnum;
         svrtx[sind1][1] = v0;
         intwk1[v0] = sind1;
      }
   }

   /* RESET INTWK1 */   
   for(i=0;i<nskeep;++i) {
      sind = skeep[i];
      intwk1[svrtx[sind][0]] = -1;
      intwk1[svrtx[sind][1]] = -1;
   } 
   
   sind = tside[tbdry].side[snum];
   intwk1[svrtx[sind][0]] = -1;
   intwk1[svrtx[nside-1][1]] = -1;

   return;
}

int mesh::maxsrch = MAXLST*3/4;

int mesh::findtri(FLT x, FLT y, int vnear) const {
   int i,j,vn,dir,stoptri,tin,tind;

   /* HERE WE USE INTWK1 THIS MUST BE -1 BEFORE USING */
   tind = vtri[vnear];
   stoptri = tind;
   dir = 1;
   ntdel = 0;
   do {
      for(vn=0;vn<3;++vn) 
         if (tvrtx[tind][vn] == vnear) break;
      
      assert(vn != 3);
      
      tind = ttri[tind][(vn +dir)%3];
      if (tind < 0) {
         if (dir > 1) break;
         /* REVERSE DIRECTION AND GO BACK TO START */
         ++dir;
         tind = vtri[vnear];
         stoptri = -1;
      }
      
      intwk1[tind] = 0;
      tdel[ntdel++] = tind;
      assert(ntdel < MAXLST -1);

      if (intri(tind,x,y) < EPSILON) goto FOUND;
         
   } while(tind != stoptri); 
   
   /* DIDN'T FIND TRIANGLE */
   /* NEED TO SEARCH SURROUNDING TRIANGLES */
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = ttri[tin][j];
         if (tind < 0) continue;
         if (intwk1[tind] == 0) continue;
         intwk1[tind] = 0;
         tdel[ntdel++] = tind;         
         if (intri(tind,x,y) < EPSILON) goto FOUND;
      }
      if (i >= maxsrch) break;
   }
   tind = -1;
   
FOUND:

   /* RESET INTWKW1 TO -1 */
   for(i=0;i<ntdel;++i)
      intwk1[tdel[i]] = -1;

   return(tind);
}


int mesh::findbdryside(FLT *x, int vnear) const {
   int i,j,tin,tind,sind;
   
   /* HERE WE USE INTWK1 THIS MUST BE -1 BEFORE USING */
   tind = vtri[vnear];
   tdel[0] = tind;
   ntdel = 1;
   /* SEARCH TRIANGLES FOR BOUNDARY SIDE CONTAINING POINT */
   for(i=0;i<ntdel;++i) {
      tin = tdel[i];
      for(j=0;j<3;++j) {
         tind = ttri[tin][j];
         if (tind < 0) {
            sind = tside[tin].side[j];
            if (insidecircle(sind,x) > -EPSILON) goto FOUND;
            continue;
         }
         if (intwk1[tind] == 0) continue;
         intwk1[tind] = 0;
         tdel[ntdel++] = tind; 
      }
      if (i >= maxsrch) break;
   }
   sind = -1;
   
FOUND:

   /* RESET INTWKW1 TO -1 */
   for(i=0;i<ntdel;++i)
      intwk1[tdel[i]] = -1;

   return(sind);
}