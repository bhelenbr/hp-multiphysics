#include<stdio.h>
#include<math.h>
#include"blocks.h"
#include"pV3.h"
#include<utilities.h>

/** WARNING THIS ONLY WORKS FOR ONE BLOCK RIGHT NOW */

#define NKEYS 6
#define MNOD 20000
#define LEN_TKEYS 40 /* MUST BE LONGER THAN 32 */
#define KNBLOCK 1

class blocks *thisblocks;

void flotov(struct vsi,int nvar, float *v);
void logflotov(struct vsi,int nvar, float *v);

int strcpn(char *str1,char *str2);
int strcpb(char *str1,char *str2, int len);

void blocks::viz_init(int iopt) {
   /**********************************************/
   /*  The important declarations needed for pV3 */
   /**********************************************/
   char data[5][80], titl[80];
   static char blank[] = " ";
   static int len_blank = 1;
   static char ttl[] = "# pV3 Initialization";
   static int len_titl;
   static INT num_keys = NKEYS;
   static INT len_tkeys = LEN_TKEYS;
   static char tkeys[NKEYS][LEN_TKEYS] =
   {{"u-Velocity                             "},
    {"v-Velocity                             "},
    {"Pressure                               "},
    {"Residual-u                             "},
    {"Residual-v                             "},
    {"Residual-p                             "}};
#ifdef SKIP
    {"x-Velocity                             "},
    {"y-Velocity                             "},     
    {"Residual-x                             "},
    {"Residual-y                             "},
    {"Flow vectors                           "}};
#endif
    
   static INT zero = 0, cid = 1;
   static INT fkeys[NKEYS] = { 1,1,1,1,1,1}; // ,1,1,1,1,2};
   static INT ikeys[NKEYS] = { 'u','v','p','U','V','P'}; // ,'x','y','X','Y','a'};
   static FLOAT flims[NKEYS][2] = {{0.0, 0.0},  {0.0, 0.0}, {0.0, 0.0}, {-16.0, 1.0},
                            {-16.0, 1.0},  {-16.0, 1.0}}; //, {0.0, 0.0},  {0.0, 0.0}, 
                        //    {0.0, 0.0},  {0.0, 0.0}, {1.0, 0.0}};

   INT istat;
   INT mirr;

/*   
   TITL:i: C Title (up to 80 characters used)
   CID:i: I This informs pV3 of the unique client-id for this client within
   a multi-client application. The index must be in the range 1 to
   the total number of clients in this discipline.
   CNAME:i: C*(*) This string identies the partition (client) with a name (up to 20 characters are used).
   *DNAME:i: C*(*) The name for the discipline. Required to be non-blank for multidiscipline cases (up to 20 characters used).
   IOPT:i: I Unsteady control parameter
   IOPT=-3 structure unsteady with connectivity supplied
   IOPT=-2 unsteady grid/data with connectivity supplied
   IOPT=-1 steady grid and unsteady data with connectivity supplied
   IOPT=0 steady grid and data
   IOPT=1 steady grid and unsteady data
   IOPT=2 unsteady grid and data
   *NPGCUT:i: I Number of programmer-defined cuts
   *TPGCUT:i: C(NPGCUT) Title for each cut (up to 32 characters used)
   *NKEYS:i: I Number of active keyboard keys
   *IKEYS:i: I(NKEYS) X-keypress return code for each key
   *TKEYS:i: C(NKEYS) Title for each key (up to 32 characters used)
   *FKEYS:i: I(NKEYS) Type of function controlled by each key:
   FKEYS()=1 Scalar
   FKEYS()=2 Vector
   FKEYS()=3 Surface scalar
   FKEYS()=4 Surface vector
   FKEYS()=5 Threshold
   FLIMS:i: R(2,NKEYS) Function limits/scales
   FKEYS()=1,3,5 Min and max values of function
   FKEYS()=2,4 Arrow/tuft scaling (only the first element is used)
   MIRROR:i: I Mirror/Replication flag:
   MIRROR=-nrep Replicate the data nrep times using REPMAT
   MIRROR=0 No mirroring or replication
   MIRROR=1 Mirror about the plane X=0.0
   MIRROR=2 Mirror about the plane Y=0.0
   MIRROR=3 Mirror about the plane Z=0.0
   REPMAT:i: R(16) The Replication matrix (required for negative MATRIX values only).
   MAXBLK:i: I The maximum number of structured blocks used during the session.
   ISTAT:i/o: I On input this sets the startup/terminate state (the following are additive):
   0 do not wait, do not terminate and Time Accurate
   1 wait for the server to startup (the first time)
   2 terminate with the server
   *4 Non-Time Accurate mode
   ISTAT = 3 is the only valid condition for steady-state cases
   (IOPT = 0).
   On output any non-zero value is the indication of a startup error
   and the task is not included in the pV3 client pool. See the
   Appendix for a list of the error codes.
   
*/
   thisblocks = this;   
   printf("%s",ttl);

   mirr   = 0;
   istat  = 0;
   len_titl = strcpn(titl,ttl);
   cid = 1;
      
   /* NOTE: string lengths are appended onto the stack */
   pV_INIT(titl, &cid, blank, blank, &iopt, &zero, &data[0][0], &num_keys,
           ikeys, &tkeys[0][0], fkeys, &flims[0][0], &mirr, NULL, &zero,
           &istat, len_titl, len_blank, len_blank, LEN_TKEYS, LEN_TKEYS);

   printf(" = %d with iopt: %d\n",istat,iopt);
   
}

/****************************************************************************
  pVStruc subroutine which supplies basic mesh information
*****************************************************************************/

void pVSTRUC (INT *knode,  INT *kequiv,  INT *kcel1,
				 INT *kcel2,  INT *kcel3,   INT *kcel4,
				 INT *knptet, INT *kptet,   INT *knblock,
				 INT *blocks, INT *kphedra, INT *ksurf, 
				 INT *knsurf, INT *hint ) {
                
   thisblocks->pvstruc(*knode, *kequiv, *kcel1, *kcel2, *kcel3, *kcel4,*knptet,*kptet,*knblock,*blocks,*kphedra,*ksurf,*knsurf,*hint);
   
   return;
}

void blocks::pvstruc(int& knode, int& kequiv, int& kcel1, int& kcel2, int& kcel3, int& kcel4, int& knptet, int
&kptet,int& knblock,int &blocks,int &kphedra, int& ksurf,int& knsurf,int& hint) {
   int i;
   static int first = 1;
   
   kequiv = 0;
   knode = 0;
   kcel1  = 0;
   kcel2  = 0;
   kcel3 = 0;
   kcel4  = 0;
   knptet = 0;
   kptet  = 0;
   knblock = 0;
   kphedra = 0; 
   knsurf = nblocks;
   ksurf = 0;
   hint = 0;

/*   if (first) {
      for(i=0;i<nblocks;++i) {
         knode += 2*blk[i].grd[0].maxvst*(1 +blk[i].grd[0].b.sm +blk[i].grd[0].b.im);
         kcel3 += blk[i].grd[0].maxvst*(blk[i].grd[0].b.sm+1)*(blk[i].grd[0].b.sm+1);
         ksurf += blk[i].grd[0].maxvst*(blk[i].grd[0].b.sm+1)*(blk[i].grd[0].b.sm+1);;
      }
      first = 0;
   } 
   else { */
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].pvstruc(knode, kequiv, kcel1, kcel2, kcel3, kcel4, knptet,kptet,knblock,blocks,kphedra,ksurf,knsurf,hint);
   
   return;
}


void hp_mgrid::pvstruc(int& knode, int& kequiv, int& kcel1, int& kcel2, int& kcel3, int& kcel4, int& knptet, int
&kptet,int& knblock,int &blocks,int &kphedra, int& ksurf,int& knsurf,int& hint) {
   int i;
   
   if (changed) {
      knode  += 2*(nvrtx +b.sm*nside+b.im*ntri);
   }
   else {
      knode = -1;
      return;
   }

   kcel3 += ntri*(b.sm+1)*(b.sm+1);
   ksurf += 2*ntri*(b.sm+1)*(b.sm+1);
   for(i=0;i<nsbd;++i)
      ksurf += sbdry[i].num*(b.sm+1);
   
   return;
}



/****************************************************************************
  pVCell subroutine which supplies connection information
*****************************************************************************/

void pVCELL(INT *cel1, INT *cel2,  INT *cel3, INT *cel4, INT *nptet, INT *ptet) {
   thisblocks->pvcell((int (*)[4]) cel1,(int (*)[5]) cel2, (int (*)[6]) cel3,(int (*)[8]) cel4,(int (*)[8]) nptet,ptet);
   return;
}

void blocks::pvcell(int cel1[][4], int cel2[][5], int cel3[][6], int cel4[][8], int nptet[][8], int ptet[]) {
   int i,kn,kpoffset;

   kn = 0;
   kpoffset = 0;
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].pvcell(kn,kpoffset,cel1,cel2,cel3,cel4,nptet,ptet);
      
   return;
}
      
   
void hp_mgrid::pvcell(int &kn,int &kpoffset, int cel1[][4], int cel2[][5], int cel3[][6], int cel4[][8], int nptet[][8], int ptet[]) {
   static int i, j, k, ijind[30][30], dkp, tind, sgn, indx;

   dkp = (nvrtx +b.sm*nside +b.im*ntri);
   
   /* FRONT FACE OF 2D MESH */
   for(tind=0;tind<ntri;++tind) {



      /* VERTICES */
      ijind[0][b.sm+1] = tvrtx[tind][0];
      ijind[0][0] = tvrtx[tind][1];
      ijind[b.sm+1][0] = tvrtx[tind][2];

      /* SIDES */
      indx = tside[tind].side[0];
      sgn = tside[tind].sign[0];
      if (sgn < 0) {
         for(i=0;i<b.sm;++i)
            ijind[i+1][0] = nvrtx +(indx+1)*b.sm -(i+1);
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[i+1][0] = nvrtx +indx*b.sm +i;
      }

      indx = tside[tind].side[1];
      sgn = tside[tind].sign[1];
      if (sgn > 0) {
         for(i=0;i<b.sm;++i)
            ijind[b.sm-i][i+1] = nvrtx +indx*b.sm +i;
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[b.sm-i][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
      }

      indx = tside[tind].side[2];
      sgn = tside[tind].sign[2];
      if (sgn > 0) {
         for(i=0;i<b.sm;++i)
            ijind[0][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[0][i+1] = nvrtx +indx*b.sm +i;
      }

      /* INTERIOR VERTICES */
      k = 0;
      for(i=1;i<b.sm;++i) {
         for(j=1;j<b.sm-(i-1);++j) {
            ijind[i][j] = nvrtx +nside*b.sm +tind*b.im +k;
            ++k;
         }
      }
      /* OUTPUT CONNECTION LIST */      
      for(i=0;i<b.sm+1;++i) {
         for(j=0;j<b.sm-i;++j) {
            cel3[kn][0] = ijind[i][j]+1 +kpoffset;
            cel3[kn][1] = ijind[i][j]+1 +dkp +kpoffset;
            cel3[kn][2] = ijind[i+1][j]+1 +dkp +kpoffset;
            cel3[kn][3] = ijind[i+1][j]+1 +kpoffset;
            cel3[kn][4] = ijind[i][j+1]+1 +dkp +kpoffset;
            cel3[kn++][5] = ijind[i][j+1]+1 +kpoffset;

            cel3[kn][0] = ijind[i+1][j]+1 +kpoffset;
            cel3[kn][1] = ijind[i+1][j]+1 +dkp +kpoffset;
            cel3[kn][2] = ijind[i+1][j+1]+1 +dkp +kpoffset;
            cel3[kn][3] = ijind[i+1][j+1]+1 +kpoffset;
            cel3[kn][4] = ijind[i][j+1]+1 +dkp +kpoffset;
            cel3[kn++][5] = ijind[i][j+1]+1 +kpoffset;   
         }
         cel3[kn][0] = ijind[i][b.sm-i]+1 +kpoffset;
         cel3[kn][1] = ijind[i][b.sm-i]+1 +dkp +kpoffset;
         cel3[kn][2] = ijind[i+1][b.sm-i]+1 +dkp +kpoffset;
         cel3[kn][3] = ijind[i+1][b.sm-i]+1 +kpoffset;
         cel3[kn][4] = ijind[i][b.sm+1-i]+1 +dkp +kpoffset;
         cel3[kn++][5] = ijind[i][b.sm+1-i]+1 +kpoffset;

      }      
   }
   
   kpoffset += 2*dkp;

   return;   
}


/****************************************************************************
  pVSURFACE subroutine which supplies surface group data
*****************************************************************************/
void pVSURFACE(INT  *nsurf, INT *scon,   INT *scel,
				   char *tsurf, int tsurfLEN) {
   
   thisblocks->pvsurface((int (*)[3]) nsurf,scon,(int (*)[4]) scel,(char (*)[20]) tsurf);
}

void blocks::pvsurface(int nsurf[][3], int scon[], int scel[][4], char tsurf[][20]) {
   int i,vrtxoffset;
   
   nsurf[0][0] = 0;
   vrtxoffset = 0;
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].pvsurface(i,vrtxoffset,nsurf,scon,scel,tsurf);

   return;
}

void hp_mgrid::pvsurface(int snum, int &offset, int nsurf[][3], int scon[], int scel[][4], char tsurf[][20]) {
   static int i, j, k, ijind[30][30], kn, tind, sgn, indx;
   char sname[30];
   
   kn = nsurf[MAX(snum-1,0)][0];  // END OF PREVIOUS SURFACE OR ZERO
      
   /* FRONT FACE OF 2D MESH */
   for(tind=0;tind<ntri;++tind) {

      /* VERTICES */
      ijind[0][b.sm+1] = tvrtx[tind][0];
      ijind[0][0] = tvrtx[tind][1];
      ijind[b.sm+1][0] = tvrtx[tind][2];

      /* SIDES */
      indx = tside[tind].side[0];
      sgn = tside[tind].sign[0];
      if (sgn < 0) {
         for(i=0;i<b.sm;++i)
            ijind[i+1][0] = nvrtx +(indx+1)*b.sm -(i+1);
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[i+1][0] = nvrtx +indx*b.sm +i;
      }

      indx = tside[tind].side[1];
      sgn = tside[tind].sign[1];
      if (sgn > 0) {
         for(i=0;i<b.sm;++i)
            ijind[b.sm-i][i+1] = nvrtx +indx*b.sm +i;
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[b.sm-i][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
      }

      indx = tside[tind].side[2];
      sgn = tside[tind].sign[2];
      if (sgn > 0) {
         for(i=0;i<b.sm;++i)
            ijind[0][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[0][i+1] = nvrtx +indx*b.sm +i;
      }

      /* INTERIOR VERTICES */
      k = 0;
      for(i=1;i<b.sm;++i) {
         for(j=1;j<b.sm-(i-1);++j) {
            ijind[i][j] = nvrtx +nside*b.sm +tind*b.im +k;
            ++k;
         }
      }
      /* OUTPUT CONNECTION LIST */      
      for(i=0;i<b.sm+1;++i) {
         for(j=0;j<b.sm-i;++j) {
            scel[kn][0] = ijind[i][j]+1 +offset;
            scel[kn][1] = ijind[i+1][j]+1 +offset;
            scel[kn][2] = ijind[i][j+1]+1 +offset;
            scel[kn][3] = 0;
            scon[kn++] = 0;

            scel[kn][0] = ijind[i+1][j]+1 +offset;
            scel[kn][1] = ijind[i+1][j+1]+1 +offset;
            scel[kn][2] = ijind[i][j+1]+1 +offset;
            scel[kn][3] = 0;
            scon[kn++] = 0;   
         }
         scel[kn][0] = ijind[i][b.sm-i]+1 +offset;
         scel[kn][1] = ijind[i+1][b.sm-i]+1 +offset;
         scel[kn][2] = ijind[i][b.sm+1-i]+1 +offset;
         scel[kn][3] = 0;
         scon[kn++] = 0;   
      }
   }

   nsurf[snum][0] = kn;
   nsurf[snum][1] = 16;
   nsurf[snum][2] = 1000  +snum;
   number_str(sname,"block",snum,1);
   strcpb(tsurf[snum],sname,20);
   
   offset += 2*(nvrtx +b.sm*nside +b.im*ntri);
   
   return;
}


/****************************************************************************
  pVGRID subroutine which passes the grid coordinates to pV3
*****************************************************************************/
void pVGRID(float *xyz) {
   thisblocks->pvgrid((float (*)[3]) xyz);
}

void blocks::pvgrid(float (*xyz)[3]) {
   int i,kn;
   
   kn = 0;
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].pvgrid(kn,xyz);
      
   return;
}

void hp_mgrid::pvgrid(int &kn, float (*xyz)[3]) {
   static int i,j,n,tind,sind,knstart;
   static double zplane;
   static int v0, v1;
   
   knstart = kn;
   zplane = -0.5;
   
   /* VERTEX MODES */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n)
         xyz[kn][n] = vrtx[i][n];
      xyz[kn++][2] = zplane;
   }
         
   if (b.p > 1) {
      /* SIDE MODES */
      for(sind=0;sind<nside;++sind) {
         if (sinfo[sind] < 0) {
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            for(n=0;n<ND;++n)
               b.proj1d_leg(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
         }
         else {
            crdtocht1d(sind);
            for(n=0;n<ND;++n)
               b.proj1d_leg(cht[n],crd[n][0]);
         }
         for(i=1;i<b.sm+1;++i) {
            for(n=0;n<ND;++n)
               xyz[kn][n] = crd[n][0][i];
            xyz[kn++][2] = zplane;
         }
      }

      /* INTERIOR MODES */
      if (b.p > 2) {
         for(tind = 0; tind < ntri; ++tind) {
         
            if (tinfo[tind] < 0) {
               for(n=0;n<ND;++n)
                  b.proj_leg(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);
            }
            else {
               crdtocht(tind);
               for(n=0;n<ND;++n)
                  b.proj_bdry_leg(cht[n],crd[n]);
            }
            
            for(i=1;i<b.sm;++i) {
               for(j=1;j<b.sm-(i-1);++j) {
                  for(n=0;n<ND;++n)
                     xyz[kn][n] = crd[n][i][j];
                  xyz[kn++][2] = zplane;
               }
            }
         }
      }
   }
   
   for(i=knstart;i<kn;++i) {
      for(n=0;n<ND;++n)
         xyz[i+kn-knstart][n] = xyz[i][n];
      xyz[i+kn-knstart][2] = 0.5;
   }
   
   kn += kn-knstart;

   return;
}

/****************************************************************************
  pV3SCAL subroutine which provides all the scalar data to pV3
*****************************************************************************/
void pVSCAL(int *key, float *v) {
   thisblocks->pvscal(key,v);
}

void blocks::pvscal(int *key, float *v) {
   int i,offset;
   
   offset = 0;
   for(i=0;i<nblocks;++i) {
      switch (*key) {
         case 1: /* u-velocity */
            blk[i].grd[0].flotov(offset,blk[i].grd[0].ug,0,v);
            break;
         case 2: /* v-velocity */
            blk[i].grd[0].flotov(offset,blk[i].grd[0].ug,1,v);
            break;
         case 3: /* pressure */
            blk[i].grd[0].flotov(offset,blk[i].grd[0].ug,2,v);
            break;
   
         /* RESIDUALS */
         case 4: /* u-velocity */
            blk[i].grd[0].logflotov(offset,blk[i].gbl.res,0,v);
            break;
         case 5: /* v-velocity */
            blk[i].grd[0].logflotov(offset,blk[i].gbl.res,1,v);
            break;
         case 6: /* pressure */
            blk[i].grd[0].logflotov(offset,blk[i].gbl.res,2,v);
            break;         
            
#ifdef SKIP         
         case 7: /* x-velocity */
   //         mvtov(mv,0,v);
            break;
         case 8: /* y-velocity */
   //         mvtov(mv,1,v);
            break;         
         case 9: /* x-velocity */
   //         mvgftov(mvgf,0,v);
            break;
         case 10: /* y-velocity */
   //         mvgftov(mvgf,1,v);
            break;
#endif
      }  
   }

   return;
}

void hp_mgrid::flotov(int &kn, struct vsi flo,int nvar, float *v) {
   static int i,j,tind,sind,knstart;
   
   knstart = kn;
   /* VERTEX MODES */
   for(i=0;i<nvrtx;++i)
      v[kn++] = flo.v[i][nvar];
   
   if (b.p > 1) {
      /* SIDE MODES */
      for(sind=0;sind<nside;++sind) {
         ugtouht1d(sind,flo);
         b.proj1d_leg(uht[nvar],u[nvar][0]);

         for(i=1;i<b.sm+1;++i)
            v[kn++] = u[nvar][0][i];   
      }

      /* INTERIOR MODES */
      if (b.p > 2) {
         for(tind = 0; tind < ntri; ++tind) {
            ugtouht(tind,flo);
            b.proj_leg(uht[nvar],u[nvar]);
               
            for(i=1;i<b.sm;++i)
               for(j=1;j<b.sm-(i-1);++j)
                  v[kn++] = u[nvar][i][j];               
         }
      }
   }
   
   for(i=knstart;i<kn;++i)
      v[i+kn-knstart] = v[i];

   kn += kn-knstart;
   
   return;
}

void hp_mgrid::logflotov(int &kn, struct vsi flo,int nvar, float *v) {
   static int i,j,tind,sind,knstart;
   
    knstart = kn;
   /* VERTEX MODES */
   for(i=0;i<nvrtx;++i)
      v[kn++] = log10(fabs(flo.v[i][nvar]) +EPSILON);
   
   if (b.p > 1) {
      /* SIDE MODES */
      for(sind=0;sind<nside;++sind) {
         ugtouht1d(sind,flo);
         b.proj1d_leg(uht[nvar],u[nvar][0]);

         for(i=1;i<b.sm+1;++i)
            v[kn++] = log10(fabs(u[nvar][0][i]) +EPSILON);   
      }

      /* INTERIOR MODES */
      if (b.p > 2) {
         for(tind = 0; tind < ntri; ++tind) {
            ugtouht(tind,flo);
            b.proj_leg(uht[nvar],u[nvar]);
               
            for(i=1;i<b.sm;++i)
               for(j=1;j<b.sm-(i-1);++j)
                  v[kn++] = log10(fabs(u[nvar][i][j])+EPSILON);               
         }
      }
   }
   
   for(i=knstart;i<kn;++i)
      v[i+kn-knstart] = v[i];

   kn += kn-knstart;

   return;
}


#ifdef SKIP
void hp_mgrid::pvvect(int *key,FLOAT v[][3]) {
   static int i,j,n,tind,sind,kn;

    kn = 0;
   /* VERTEX MODES */
   for(i=0;i<nvrtx;++i) {
      v[kn][0] = ug.v[i][0];
        v[kn][1] = ug.v[i][1];
        v[kn++][2] = 0.0;
   }
   
   if (b.p > 1) {
      /* SIDE MODES */
      for(sind=0;sind<nside;++sind) {
         ugtouht1d(sind);
         for(n=0;n<ND;++n)
            b.proj1d_leg(uht[n],u[n][0]);

         for(i=1;i<b.sm+1;++i) {
            for(n=0;n<ND;++n)
               v[kn][n] = u[n][0][i];   
            v[kn++][2] = 0.0;
         }
      }

      /* INTERIOR MODES */
      if (b.p > 2) {
         for(tind = 0; tind < ntri; ++tind) {
            ugtouht(tind);
            for(n=0;n<ND;++n)
               b.proj_leg(uht[n],u[n]);
               
            for(i=1;i<b.sm;++i) {
               for(j=1;j<b.sm-(i-1);++j) {
                  for(n=0;n<ND;++n)
                     v[kn][n] = u[n][i][j];               
                  v[kn++][2] = 0.0;
               }
            }
         }
      }
   }

    for(i=0;i<kn;++i)
      for(j=0;j<3;++j)
         v[i+kn][j] = v[i][j];  
   
   return;
}
#endif

/* strcpb copies str2 into str1 up to, but not including the first null and
   then pads with blanks to the len character */
   
int strcpb(char *str1,char *str2, int len) {
    register int i,j;

    i = 0;
    while (str2[i] != '\0') {
      str1[i] = str2[i];
      i++;
    }
    for (j=i; j<len; j++) str1[j] = ' ';

    return i;
}

/* strcpn copies str2 into str1 up to, but not including the first null */

int strcpn(char *str1,char *str2) {
    register int i;

    i = 0;
    while (str2[i] != '\0') {
      str1[i] = str2[i];
      i++;
    }

    return i;
}
