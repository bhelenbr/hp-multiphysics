#include<stdio.h>
#include<math.h>
#include"blocks.h"
#include"osdepend.h"

/** WARNING THIS ONLY WORKS FOR ONE BLOCK RIGHT NOW */

#define NKEYS 6
#define MNOD 20000
#define LEN_TKEYS 40 /* MUST BE LONGER THAN 32 */
#define KNBLOCK 1

extern "C" void pVCELL (int cel1[][4], int cel2[][5], int cel3[][6], int cel4[][8], int nptet[][8], int ptet[]);
extern "C" void pVSURFACE(int nsurf[][3], int scon[], int scel[][4], char tsurf[][20]);
extern "C" void pVGRID(float (*xyz)[3]);
extern "C" void pVSCAL(int *key, float *v);
extern "C" void pVSCAL(int *key, float *v);
extern "C" void pVVECT(int *key,FLOAT v[][3]);
extern "C" void pVSTRUC(int& knode, int& kequiv, int& kcel1, int& kcel2, int& kcel3, int& kcel4, int& knptet, int
&kptet,int& knblock,int blocks[][3],int& ksurf,int& knsurf,int& hint);

int npvgrds;
class hp_mgrid **pvgrds;

void flotov(struct vsi,int nvar, float *v);

int strcpn(char *str1,char *str2);
int strcpb(char *str1,char *str2, int len);

void blocks::viz_init(void) {
/**********************************************/
/*  The important declarations needed for pV3 */
/**********************************************/
   int i;
	static char ttl[] = "# pV3 Initialization";
	static char titl[80];
	static int len_titl;
	int num_keys = NKEYS;
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
    
	static INT zero = 0;
	static INT fkeys[NKEYS] = { 1,1,1,1,1,1}; // ,1,1,1,1,2};
	static INT ikeys[NKEYS] = { 'u','v','p','U','V','P'}; // ,'x','y','X','Y','a'};
	static FLOAT flims[NKEYS][2] = {{0.0, 0.0},  {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
									 {0.0, 0.0},  {0.0, 0.0}}; //, {0.0, 0.0},  {0.0, 0.0}, 
								//	 {0.0, 0.0},  {0.0, 0.0}, {1.0, 0.0}};

	INT iopt, mirr, knode, kequiv, kcel1, kcel2, kcel3, kcel4, knptet, kptet;
	INT knblock, iblock[KNBLOCK][3], ksurf, knsurf, istat;
   
   pvgrds = new (hp_mgrid *)[nblocks];
   for(i=0;i<nblocks;++i)
      pvgrds[i] = &blk[i].grd[0];
   npvgrds = nblocks;
   
	printf("%s",ttl);

/*	SET UP VISUAL3 CONSTANTS 
	IOPT   = 0   STEADY GRID AND DATA
	IOPT   = 1   STEADY GRID AND UNSTEADY DATA
	IOPT   = 2   UNSTEADY GRID AND DATA

	MIRR   = 0   NO MIRRORING
	MIRR   = 1   MIRRORING ABOUT X=0.0 PLANE
	MIRR   = 1   MIRRORING ABOUT Y=0.0 PLANE
	MIRR   = 1   MIRRORING ABOUT K=0.0 PLANE

	KNODE        NUMBER OF NON-STRUCTURED BLOCK NODES
	KEQUIV       NUMBER OF EQUIVALENCE PAIRS
	KCEL1        NUMBER OF TETRAHEDRA
	KCEL2        NUMBER OF PYRAMIDS
	KCEL3        NUMBER OF PRISMS
	KCEL4        NUMBER OF HEXAHEDRA
	KNPTET       NUMBER OF POLY-TETRAHEDRAL STRIPS
	KPTET        NUMBER OF TETRAHEDRAL CELLS IN ALL STRIPS
	
	KNBLOCK      NUMBER OF STRUCTURED BLOCKS
	IBLOCKS(3,*) BLOCK DIMENSIONS

	KSURF        NUMBER OF SURFACE FACES
	KNSURF       NUMBER OF SURFACE GROUPS
*/

	iopt   = -3;
	mirr   = 0;
   knode  = 2*blk[0].grd[0].maxvst*(1 +blk[0].grd[0].b.sm +blk[0].grd[0].b.im);
	kequiv = 0;
	kcel1  = 0;
	kcel2  = 0;
   kcel3  = blk[0].grd[0].maxvst*(blk[0].grd[0].b.sm+1)*(blk[0].grd[0].b.sm+1);
	kcel4  = 0;
	knptet = 0;
	kptet  = 0;
	knblock= 0;
	 
	knsurf = 1;
	ksurf  = kcel3;
		
	istat  = 0;
   len_titl = strcpn(titl,ttl);
	
	pV_INIT (titl, &iopt, &zero, &tkeys[0][0], &num_keys, ikeys, &tkeys[0][0],
			   fkeys, &flims[0][0], &mirr, &knode, &kequiv, &kcel1, &kcel2, &kcel3,
				&kcel4, &knptet, &kptet, &knblock, &iblock[0][0], &ksurf, &knsurf,
				&istat, len_titl, LEN_TKEYS, LEN_TKEYS);

	printf(" = %d\n",istat);
	
}


/****************************************************************************
  pVCell subroutine which supplies connection information
*****************************************************************************/

void pVCELL (int cel1[][4], int cel2[][5], int cel3[][6], int cel4[][8], int nptet[][8], int ptet[]) {
   pvgrds[0]->pvcell(cel1,cel2,cel3,cel4,nptet,ptet);
   return;
}
   
void hp_mgrid::pvcell(int cel1[][4], int cel2[][5], int cel3[][6], int cel4[][8], int nptet[][8], int ptet[]) {
	static int i, j, k, ijind[30][30], kp, kn, tind, sgn, indx;

   kp = (nvrtx +b.sm*nside +b.im*ntri);
	kn = 0;
	
/*	FRONT FACE OF 2D MESH */
   for(tind=0;tind<ntri;++tind) {

/*		VERTICES */
      ijind[0][0] = tvrtx[tind][0];
      ijind[b.sm+1][0] = tvrtx[tind][1];
      ijind[0][b.sm+1] = tvrtx[tind][2];

/*		SIDES */		   
      indx = tside[tind].side[0];
      sgn = tside[tind].sign[0];
      if (sgn > 0) {
         for(i=0;i<b.sm;++i)
            ijind[b.sm-i][i+1] = nvrtx +indx*b.sm +i;
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[b.sm-i][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
      }
 
      indx = tside[tind].side[1];
      sgn = tside[tind].sign[1];
      if (sgn > 0) {
         for(i=0;i<b.sm;++i)
            ijind[0][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[0][i+1] = nvrtx +indx*b.sm +i;
      }
      
      indx = tside[tind].side[2];
      sgn = tside[tind].sign[2];
      if (sgn < 0) {
         for(i=0;i<b.sm;++i)
            ijind[i+1][0] = nvrtx +(indx+1)*b.sm -(i+1);
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[i+1][0] = nvrtx +indx*b.sm +i;
      }

/*		INTERIOR VERTICES */
      k = 0;
      for(i=1;i<b.sm;++i) {
         for(j=1;j<b.sm-(i-1);++j) {
            ijind[i][j] = nvrtx +nside*b.sm +tind*b.im +k;
            ++k;
         }
      }
		
/*		OUTPUT CONNECTION LIST */		
   	for(i=0;i<b.sm+1;++i) {
      	for(j=0;j<b.sm-i;++j) {
				cel3[kn][0] = ijind[i][j]+1;
				cel3[kn][1] = ijind[i][j]+1 +kp;
				cel3[kn][2] = ijind[i+1][j]+1 +kp;
				cel3[kn][3] = ijind[i+1][j]+1;
				cel3[kn][4] = ijind[i][j+1]+1 +kp;
				cel3[kn++][5] = ijind[i][j+1]+1;

				cel3[kn][0] = ijind[i+1][j]+1;
				cel3[kn][1] = ijind[i+1][j]+1 +kp;
				cel3[kn][2] = ijind[i+1][j+1]+1 +kp;
				cel3[kn][3] = ijind[i+1][j+1]+1;
				cel3[kn][4] = ijind[i][j+1]+1 +kp;
				cel3[kn++][5] = ijind[i][j+1]+1;	
      	}
			cel3[kn][0] = ijind[i][b.sm-i]+1;
			cel3[kn][1] = ijind[i][b.sm-i]+1 +kp;
			cel3[kn][2] = ijind[i+1][b.sm-i]+1 +kp;
			cel3[kn][3] = ijind[i+1][b.sm-i]+1;
			cel3[kn][4] = ijind[i][b.sm+1-i]+1 +kp;
			cel3[kn++][5] = ijind[i][b.sm+1-i]+1;

   	}		
	}

	return;	
}

void pVSTRUC(int& knode, int& kequiv, int& kcel1, int& kcel2, int& kcel3, int& kcel4, int& knptet, int
&kptet,int& knblock,int blocks[][3],int& ksurf,int& knsurf,int& hint) {

	pvgrds[0]->pvstruc(knode, kequiv, kcel1, kcel2, kcel3, kcel4,knptet,kptet,knblock,blocks,ksurf,knsurf,hint);
	
	return;
}

void hp_mgrid::pvstruc(int& knode, int& kequiv, int& kcel1, int& kcel2, int& kcel3, int& kcel4, int& knptet, int
&kptet,int& knblock,int blocks[][3],int& ksurf,int& knsurf,int& hint) {

	if (changed) {
   	knode  = 2*(nvrtx +b.sm*nside+b.im*ntri);
	}
	else {
		knode = -1;
		return;
	}
	kequiv = 0;
	kcel1  = 0;
	kcel2  = 0;
   kcel3 = 0;
   kcel3  = ntri*(b.sm+1)*(b.sm+1);
	kcel4  = 0;
	knptet = 0;
	kptet  = 0;
	knblock= 0;
	 
	knsurf = 1;
	ksurf  = kcel3;
	
	hint = 0;
	return;
}



/****************************************************************************
  pVSURFACE subroutine which supplies surface group data
*****************************************************************************/
void pVSURFACE(int nsurf[][3], int scon[], int scel[][4], char tsurf[][20]) {
   pvgrds[0]->pvsurface(nsurf,scon,scel,tsurf);
}

void hp_mgrid::pvsurface(int nsurf[][3], int scon[], int scel[][4], char tsurf[][20]) {
	static int i, j, k, ijind[30][30], kp, kn, tind, sgn, indx;
   
   kp = (nvrtx +b.sm*nside +b.im*ntri);
	kn = 0;
	
/*	FRONT FACE OF 2D MESH */
   for(tind=0;tind<ntri;++tind) {

/*		VERTICES */
      ijind[0][0] = tvrtx[tind][0];
      ijind[b.sm+1][0] = tvrtx[tind][1];
      ijind[0][b.sm+1] = tvrtx[tind][2];

/*		SIDES */		   
      indx = tside[tind].side[0];
      sgn = tside[tind].sign[0];
      if (sgn > 0) {
         for(i=0;i<b.sm;++i)
            ijind[b.sm-i][i+1] = nvrtx +indx*b.sm +i;
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[b.sm-i][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
      }

      indx = tside[tind].side[1];
      sgn = tside[tind].sign[1];
      if (sgn > 0) {
         for(i=0;i<b.sm;++i)
            ijind[0][i+1] = nvrtx +(indx+1)*b.sm -(i+1);
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[0][i+1] = nvrtx +indx*b.sm +i;
      }
      
      indx = tside[tind].side[2];
      sgn = tside[tind].sign[2];
      if (sgn < 0) {
         for(i=0;i<b.sm;++i)
            ijind[i+1][0] = nvrtx +(indx+1)*b.sm -(i+1);
      }
      else {
         for(i=0;i<b.sm;++i)
            ijind[i+1][0] = nvrtx +indx*b.sm +i;
      }

/*		INTERIOR VERTICES */
      k = 0;
      for(i=1;i<b.sm;++i) {
         for(j=1;j<b.sm-(i-1);++j) {
            ijind[i][j] = nvrtx +nside*b.sm +tind*b.im +k;
            ++k;
         }
      }

/*		OUTPUT CONNECTION LIST */		
      for(i=0;i<b.sm+1;++i) {
         for(j=0;j<b.sm-i;++j) {
				scel[kn][0] = ijind[i][j]+1;
				scel[kn][1] = ijind[i+1][j]+1;
				scel[kn][2] = ijind[i][j+1]+1;
				scel[kn][3] = 0;
				scon[kn++] = 0;

				scel[kn][0] = ijind[i+1][j]+1;
				scel[kn][1] = ijind[i+1][j+1]+1;
				scel[kn][2] = ijind[i][j+1]+1;
				scel[kn][3] = 0;
				scon[kn++] = 0;	
         }
         scel[kn][0] = ijind[i][b.sm-i]+1;
         scel[kn][1] = ijind[i+1][b.sm-i]+1;
         scel[kn][2] = ijind[i][b.sm+1-i]+1;
         scel[kn][3] = 0;
         scon[kn++] = 0;	
      }
   }
   nsurf[0][0] = kn;
   nsurf[0][1] = 16;
   nsurf[0][2] = 1000;
   strcpb(tsurf[0],"2D MESH",20);
	
   
   return;
}


/****************************************************************************
  pVGRID subroutine which passes the grid coordinates to pV3
*****************************************************************************/
void pVGRID(float (*xyz)[3]) {
   pvgrds[0]->pvgrid(xyz);
}

void hp_mgrid::pvgrid(float (*xyz)[3]) {
	static int i,j,n,tind,sind,kn;
	static double zplane;
	static int v0, v1;
   
	kn = 0;
	
	zplane = -0.5;
	
/*	VERTEX MODES */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n)
         xyz[kn][n] = vrtx[i][n];
      xyz[kn++][2] = zplane;
   }
         
   if (b.p > 1) {
/*		SIDE MODES */
      for(sind=0;sind<nside;++sind) {
         if (sinfo[sind] < 0) {
            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];
            for(n=0;n<ND;++n)
               b.proj1d_leg(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
         }
         else {
            crdtouht1d(sind);
            for(n=0;n<ND;++n)
               b.proj1d_leg(uht[n],crd[n][0]);
         }
         for(i=1;i<b.sm+1;++i) {
            for(n=0;n<ND;++n)
               xyz[kn][n] = crd[n][0][i];
            xyz[kn++][2] = zplane;
         }
      }

/*		INTERIOR MODES */
      if (b.p > 2) {
         for(tind = 0; tind < ntri; ++tind) {
         
            if (tinfo[tind] < 0) {
               for(n=0;n<ND;++n)
                  b.proj_leg(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);
            }
            else {
               crdtouht(tind);
               for(n=0;n<ND;++n)
                  b.proj_bdry_leg(uht[n],crd[n]);
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
	
	for(i=0;i<kn;++i) {
		for(n=0;n<ND;++n)
			xyz[i+kn][n] = xyz[i][n];
		xyz[i+kn][2] = 0.5;
	}


	return;
}

/****************************************************************************
  pV3SCAL subroutine which provides all the scalar data to pV3
*****************************************************************************/
void pVSCAL(int *key, float *v) {
   pvgrds[0]->pvscal(key,v);
}

void hp_mgrid::pvscal(int *key, float *v) {

	switch (*key) {
		case 1: /* u-velocity */
			flotov(ug,0,v);
			break;
		case 2: /* v-velocity */
			flotov(ug,1,v);
			break;
		case 3: /* pressure */
			flotov(ug,2,v);
			break;

/*		RESIDUALS */
		case 4: /* u-velocity */
			flotov(gbl->res,0,v);
			break;
		case 5: /* v-velocity */
			flotov(gbl->res,1,v);
			break;
		case 6: /* pressure */
			flotov(gbl->res,2,v);
			break;			
			
#ifdef SKIP			
		case 7: /* x-velocity */
//			mvtov(mv,0,v);
			break;
		case 8: /* y-velocity */
//			mvtov(mv,1,v);
			break;			
		case 9: /* x-velocity */
//			mvgftov(mvgf,0,v);
			break;
		case 10: /* y-velocity */
//			mvgftov(mvgf,1,v);
			break;
#endif
	}				

	return;
}

void pVVECT(int *key,FLOAT v[][3]) {
   pvgrds[0]->pvvect(key,v);
}

void hp_mgrid::pvvect(int *key,FLOAT v[][3]) {
	static int i,j,n,tind,sind,kn;

 	kn = 0;
/*	VERTEX MODES */
   for(i=0;i<nvrtx;++i) {
      v[kn][0] = ug.v[i][0];
  		v[kn][1] = ug.v[i][1];
  		v[kn++][2] = 0.0;
   }
   
   if (b.p > 1) {
/*		SIDE MODES */
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

/*		INTERIOR MODES */
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


void hp_mgrid::flotov(struct vsi flo,int nvar, float *v) {
	static int i,j,tind,sind,kn;
   
 	kn = 0;
/*	VERTEX MODES */
   for(i=0;i<nvrtx;++i)
      v[kn++] = flo.v[i][nvar];
   
   if (b.p > 1) {
/*		SIDE MODES */
      for(sind=0;sind<nside;++sind) {
         ugtouht1d(sind,flo);
         b.proj1d_leg(uht[nvar],u[nvar][0]);

         for(i=1;i<b.sm+1;++i)
            v[kn++] = u[nvar][0][i];	
      }

/*		INTERIOR MODES */
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
	
	for(i=0;i<kn;++i)
		v[i+kn] = v[i];

	return;
}

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
