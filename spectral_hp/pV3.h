/*
 *	pV3 Rev 2 Client Include
 *
 *	Allows OS/machine dependent stuff to work with FORTRAN
 *
 *	Copyright 1994 - 2001, Massachusetts Institute of Technology.
 */

/*
 *      pV3.h
 *
 *      Note: if using MPI, have this include specified after mpi.h
 *
 *
$Log$
Revision 1.1  2002/05/22 13:54:08  helenbrk
Upgrading pV3.  Added ptprobe routine for insert into curved triangles.  Increased number of max iterations for findcurvedpt.

 */

#ifndef _PV3_H_

#define _PV3_H_

#ifdef ABSOFT
#undef WIN32
#endif


/* main program entry point definition */

#undef MAINPROG

#ifdef DOUBUNS
#define MAINPROG   MAIN__
#endif

#ifdef NOUNS
#define MAINPROG   main
#endif

#ifdef LCUNS
#define MAINPROG   main_
#endif

#ifdef CRAY
#define MAINPROG   main
#endif

#ifndef MAINPROG
#define MAINPROG   MAIN_
#endif


/* routine definitions */

#ifdef CRAY
#define pV_INIT     PV_INIT
#define pV_UPDATE   PV_UPDATE
#define pV_STAT     PV_STAT
#define pV_CONSOLE  PV_CONSOLE
#define pV_TERMIN   PV_TERMIN
#define pV_GETSTRUC PV_GETSTRUC
#define pV_FIELD    PV_FIELD
#define pV_SENDXI   PV_SENDXI
#define pV_SENDXR   PV_SENDXR
#define pV_GSERDATA PV_GETSERVERDATA
#define pV_MPISTOP  PV_MPISTOP
#define pVCELL      PVCELL
#define pVSURFACE   PVSURFACE
#define pVEQUIV     PVEQUIV
#define pVBLANK     PVBLANK
#define pVGRID      PVGRID
#define pVSCAL      PVSCAL
#define pVTHRES     PVTHRES
#define pVVECT      PVVECT
#define pVSTRUC     PVSTRUC
#define pVLOCATE    PVLOCATE
#define pVCONNECT   PVCONNECT
#define pVPOLYHEDRA PVPOLYHEDRA
#define pVPGRID     PVPGRID
#define pVPSCAL     PVPSCAL
#define pVPTHRES    PVPTHRES
#define pVPVECT     PVPVECT
#define pVZPRIME    PVZPRIME
#define pVXYPRIME   PVXYPRIME
#define pVPZPRIME   PVPZPRIME
#define pVPXYPRIME  PVPXYPRIME
#define pVSURF      PVSURF
#define pVXYSURF    PVXYSURF
#define pVPXYSURF   PVPXYSURF
#define pVSSURF     PVSSURF
#define pVVSURF     PVVSURF
#define pVPSSURF    PVPSSURF
#define pVPVSURF    PVPVSURF
#define pVSTRING    PVSTRING
#define pVCATCH     PVCATCH
#define pVEXTRACT   PVEXTRACT
#define pVSTATE     PVSTATE
#define pVPSTATE    PVPSTATE
#define pVSSINIT    PVSSINIT
#define pVUPDATE    PVUPDATE
#define pVNOTIFY    PVNOTIFY
#define pVLOCAL     PVLOCAL
#define pVBOUNCE    PVBOUNCE
#define pVSETSTRUC  PVSETSTRUC
#define pREFIX      void
#define pREFI       int
#endif

#ifdef WIN32
#define pV_INIT     PV_INIT
#define pV_UPDATE   PV_UPDATE
#define pV_STAT     PV_STAT
#define pV_CONSOLE  PV_CONSOLE
#define pV_TERMIN   PV_TERMIN
#define pV_GETSTRUC PV_GETSTRUC
#define pV_FIELD    PV_FIELD
#define pV_SENDXI   PV_SENDXI
#define pV_SENDXR   PV_SENDXR
#define pV_GSERDATA PV_GETSERVERDATA
#define pV_MPISTOP  PV_MPISTOP
#define pVCELL      PVCELL
#define pVSURFACE   PVSURFACE
#define pVEQUIV     PVEQUIV
#define pVBLANK     PVBLANK
#define pVGRID      PVGRID
#define pVSCAL      PVSCAL
#define pVTHRES     PVTHRES
#define pVVECT      PVVECT
#define pVSTRUC     PVSTRUC
#define pVLOCATE    PVLOCATE
#define pVCONNECT   PVCONNECT
#define pVPOLYHEDRA PVPOLYHEDRA
#define pVPGRID     PVPGRID
#define pVPSCAL     PVPSCAL
#define pVPTHRES    PVPTHRES
#define pVPVECT     PVPVECT
#define pVZPRIME    PVZPRIME
#define pVXYPRIME   PVXYPRIME
#define pVPZPRIME   PVPZPRIME
#define pVPXYPRIME  PVPXYPRIME
#define pVSURF      PVSURF
#define pVXYSURF    PVXYSURF
#define pVPXYSURF   PVPXYSURF
#define pVSSURF     PVSSURF
#define pVVSURF     PVVSURF
#define pVPSSURF    PVPSSURF
#define pVPVSURF    PVPVSURF
#define pVSTRING    PVSTRING
#define pVCATCH     PVCATCH
#define pVEXTRACT   PVEXTRACT
#define pVSTATE     PVSTATE
#define pVPSTATE    PVPSTATE
#define pVSSINIT    PVSSINIT
#define pVUPDATE    PVUPDATE
#define pVNOTIFY    PVNOTIFY
#define pVLOCAL     PVLOCAL
#define pVBOUNCE    PVBOUNCE
#define pVSETSTRUC  PVSETSTRUC
#define pREFIX      void __stdcall
#define pREFI       int  __stdcall
#endif

#ifdef ABSOFT
#define pV_INIT     pv_init
#define pV_UPDATE   pv_update
#define pV_STAT     pv_stat
#define pV_CONSOLE  pv_console
#define pV_TERMIN   pv_termin
#define pV_GETSTRUC pv_getstruc
#define pV_FIELD    pv_field
#define pV_SENDXI   pv_sendxi
#define pV_SENDXR   pv_sendxr
#define pV_GSERDATA pv_getserverdata
#define pV_MPISTOP  pv_mpistop
#define pVCELL      pvcell
#define pVSURFACE   pvsurface
#define pVEQUIV     pvequiv
#define pVBLANK     pvblank
#define pVGRID      pvgrid
#define pVSCAL      pvscal
#define pVTHRES     pvthres
#define pVVECT      pvvect
#define pVSTRUC     pvstruc
#define pVLOCATE    pvlocate
#define pVPOLYHEDRA pvpolyhedra
#define pVPGRID     pvpgrid
#define pVPSCAL     pvpscal
#define pVPTHRES    pvpthres
#define pVPVECT     pvpvect
#define pVCONNECT   pvconnect
#define pVZPRIME    pvzprime
#define pVXYPRIME   pvxyprime
#define pVPZPRIME   pvpzprime
#define pVPXYPRIME  pvpxyprime
#define pVSURF      pvsurf
#define pVXYSURF    pvxysurf
#define pVPXYSURF   pvpxysurf
#define pVSSURF     pvssurf
#define pVVSURF     pvvsurf
#define pVPSSURF    pvpssurf
#define pVPVSURF    pvpvsurf
#define pVSTRING    pvstring
#define pVCATCH     pvcatch
#define pVEXTRACT   pvextract
#define pVSTATE     pvstate
#define pVPSTATE    pvpstate
#define pVSSINIT    pvssinit
#define pVUPDATE    pvupdate
#define pVNOTIFY    pvnotify
#define pVLOCAL     pvlocal
#define pVBOUNCE    pvbounce
#define pVSETSTRUC  pvsetstruc
#define pREFIX      void __cdecl
#define pREFI       int  __cdecl
#endif

#ifdef UNDERSCORE
#define pV_INIT     pv_init_
#define pV_UPDATE   pv_update_
#define pV_STAT     pv_stat_
#define pV_CONSOLE  pv_console_
#define pV_TERMIN   pv_termin_
#define pV_GETSTRUC pv_getstruc_
#define pV_FIELD    pv_field_
#define pV_SENDXI   pv_sendxi_
#define pV_SENDXR   pv_sendxr_
#define pV_GSERDATA pv_getserverdata_
#define pV_MPISTOP  pv_mpistop_
#define pVCELL      pvcell_
#define pVSURFACE   pvsurface_
#define pVEQUIV     pvequiv_
#define pVBLANK     pvblank_
#define pVGRID      pvgrid_
#define pVSCAL      pvscal_
#define pVTHRES     pvthres_
#define pVVECT      pvvect_
#define pVSTRUC     pvstruc_
#define pVLOCATE    pvlocate_
#define pVCONNECT   pvconnect_
#define pVPOLYHEDRA pvpolyhedra_
#define pVPGRID     pvpgrid_
#define pVPSCAL     pvpscal_
#define pVPTHRES    pvpthres_
#define pVPVECT     pvpvect_
#define pVZPRIME    pvzprime_
#define pVXYPRIME   pvxyprime_
#define pVPZPRIME   pvpzprime_
#define pVPXYPRIME  pvpxyprime_
#define pVSURF      pvsurf_
#define pVXYSURF    pvxysurf_
#define pVPXYSURF   pvpxysurf_
#define pVSSURF     pvssurf_
#define pVVSURF     pvvsurf_
#define pVPSSURF    pvpssurf_
#define pVPVSURF    pvpvsurf_
#define pVSTRING    pvstring_
#define pVCATCH     pvcatch_
#define pVEXTRACT   pvextract_
#define pVSTATE     pvstate_
#define pVPSTATE    pvpstate_
#define pVSSINIT    pvssinit_
#define pVUPDATE    pvupdate_
#define pVNOTIFY    pvnotify_
#define pVLOCAL     pvlocal_
#define pVBOUNCE    pvbounce_
#define pVSETSTRUC  pvsetstruc_
#define pREFIX      void
#define pREFI       int
#endif

#ifdef NOUNDERSCORE
#define pV_INIT     pv_init
#define pV_UPDATE   pv_update
#define pV_STAT     pv_stat
#define pV_CONSOLE  pv_console
#define pV_TERMIN   pv_termin
#define pV_GETSTRUC pv_getstruc
#define pV_FIELD    pv_field
#define pV_SENDXI   pv_sendxi
#define pV_SENDXR   pv_sendxr
#define pV_GSERDATA pv_getserverdata
#define pV_MPISTOP  pv_mpistop
#define pVCELL      pvcell
#define pVSURFACE   pvsurface
#define pVEQUIV     pvequiv
#define pVBLANK     pvblank
#define pVGRID      pvgrid
#define pVSCAL      pvscal
#define pVTHRES     pvthres
#define pVVECT      pvvect
#define pVSTRUC     pvstruc
#define pVLOCATE    pvlocate
#define pVPOLYHEDRA pvpolyhedra
#define pVPGRID     pvpgrid
#define pVPSCAL     pvpscal
#define pVPTHRES    pvpthres
#define pVPVECT     pvpvect
#define pVCONNECT   pvconnect
#define pVZPRIME    pvzprime
#define pVXYPRIME   pvxyprime
#define pVPZPRIME   pvpzprime
#define pVPXYPRIME  pvpxyprime
#define pVSURF      pvsurf
#define pVXYSURF    pvxysurf
#define pVPXYSURF   pvpxysurf
#define pVSSURF     pvssurf
#define pVVSURF     pvvsurf
#define pVPSSURF    pvpssurf
#define pVPVSURF    pvpvsurf
#define pVSTRING    pvstring
#define pVCATCH     pvcatch
#define pVEXTRACT   pvextract
#define pVSTATE     pvstate
#define pVPSTATE    pvpstate
#define pVSSINIT    pvssinit
#define pVUPDATE    pvupdate
#define pVNOTIFY    pvnotify
#define pVLOCAL     pvlocal
#define pVBOUNCE    pvbounce
#define pVSETSTRUC  pvsetstruc
#define pREFIX      void
#define pREFI       int
#endif


/* fortran integer definitions */

#ifdef INTEGERLONG
#define INT long
#else
#define INT int
#endif


/* fortran real definitions */

#ifdef FLOATDOUBLE
#define FLOAT double
#else
#define FLOAT float
#endif


/* routine prototypes */

#ifdef __ProtoGlarp__
#undef __ProtoGlarp__
#endif
#if defined(__STDC__) || defined(__cplusplus) || defined(WIN32)
#define __ProtoGlarp__(x) x
#else
#define __ProtoGlarp__(x) ()
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef CRAY

#ifdef T3E
pREFIX	pV_INIT __ProtoGlarp__(( char *title,  int titleLEN,  INT *cid, 
				 char *cname,  int cnameLEN,  char *dname, 
				 int dnameLEN, INT *iopt,     INT *npgcut,
				 char *tpgcut, int tpgcutLEN, INT *nkeys, 
				 INT *ikeys,   char *tkeys,   int tkeysLEN,
				 INT *fkeys,   FLOAT *flims,  INT *mirror,
				 FLOAT *rmat,  INT *maxblk,   INT *istat ));
pREFIX	pV_CONSOLE __ProtoGlarp__(( char *string, int stringLEN ));
#else
/* title, tpgcut and tkeys are pointers packed with the string len */
pREFIX	pV_INIT __ProtoGlarp__(( int title,   int *cid,     int cname,
				 int dname,   int *iopt,    int *npgcut,
				 int tpgcut,  int *nkeys,   int *ikeys,
				 int tkeys,   int *fkeys,   float *flims,
				 int *mirror, float *rmat,  int *maxblk,
				 int *istat ));
/* string is a char pointer packed with the string length */
pREFIX	pV_CONSOLE __ProtoGlarp__(( int string ));
#endif

#else

#ifdef WIN32
pREFIX  pV_INIT __ProtoGlarp__(( char *title,  int titleLEN,  INT *cid,
				 char *cname,  int cnameLEN,  char *dname, 
				 int dnameLEN, INT *iopt,     INT *npgcut,
				 char *tpgcut, int tpgcutLEN, INT *nkeys, 
				 INT *ikeys,   char *tkeys,   int tkeysLEN,
				 INT *fkeys,   FLOAT *flims,  INT *mirror,
				 FLOAT *rmat,  INT *maxblk,   INT *istat ));
#else
pREFIX  pV_INIT __ProtoGlarp__(( char *title,   INT *cid,      char *cname,
				 char *dname,   INT *iopt,     INT *npgcut,
				 char *tpgcut,  INT *nkeys,    INT *ikeys,
				 char *tkeys,   INT *fkeys,    FLOAT *flims,
				 INT *mirror,   FLOAT *rmat,   INT *maxblk,
				 INT *istat,    int titleLEN,  int cnameLEN,
				 int dnameLEN,  int tpgcutLEN, int tkeysLEN  ));
#endif
pREFIX	pV_CONSOLE __ProtoGlarp__(( char *string, int stringLEN ));
#endif

pREFIX	pV_UPDATE __ProtoGlarp__(( FLOAT *time ));
pREFIX	pV_STAT  __ProtoGlarp__(( INT *istate ));
pREFIX	pV_TERMIN __ProtoGlarp__(( void ));
pREFIX	pV_GETSTRUC __ProtoGlarp__(( INT *opt, void **ptr, INT *len ));
pREFIX	pV_FIELD __ProtoGlarp__(( INT *ls, INT *lv, INT *lt, INT *lz ));
pREFIX	pV_SENDXI __ProtoGlarp__(( INT *index,   INT *exnum, INT *subindex,
                                   INT *subsize, INT *len,   INT *ivec ));
pREFIX	pV_SENDXR __ProtoGlarp__(( INT *index,   INT *exnum, INT *subindex,
                                   INT *subsize, INT *len,   FLOAT *rvec ));
pREFIX	pV_GSERDATA __ProtoGlarp__(( INT *handle, INT *mid,    INT *nint,
                                     INT **ints,  INT *nfloat, FLOAT **floats,
                                     INT *nchar,  char **chars ));

#ifdef _MPI_INCLUDE
int	pV_MPIStart __ProtoGlarp__(( MPI_Comm cin, int rank, int ncl,
	                             int bsize, MPI_Comm *cout ));
void	pV_MPISTOP __ProtoGlarp__(( void ));
#endif

pREFIX	pVCELL __ProtoGlarp__(( INT *cel1, INT *cel2,  INT *cel3,
				INT *cel4, INT *nptet, INT *ptet ));
pREFIX	pVSURFACE __ProtoGlarp__(( INT  *nsurf, INT *scon,   INT *scel,
				   char *tsurf, int tsurfLEN ));
pREFIX	pVEQUIV __ProtoGlarp__(( INT *listeq ));
pREFIX	pVBLANK __ProtoGlarp__(( INT *iblank, INT *tidcon ));
pREFIX	pVGRID __ProtoGlarp__(( FLOAT *xyz ));
pREFIX	pVSCAL __ProtoGlarp__(( INT *jkey, FLOAT *s ));
pREFIX	pVTHRES __ProtoGlarp__(( INT *jkey, FLOAT *xyz, FLOAT *t ));
pREFIX	pVVECT __ProtoGlarp__(( INT *jkey, FLOAT *v ));
pREFIX	pVSTRUC __ProtoGlarp__(( INT *knode,  INT *kequiv,  INT *kcel1,
				 INT *kcel2,  INT *kcel3,   INT *kcel4,
				 INT *knptet, INT *kptet,   INT *knblock,
				 INT *blocks, INT *kphedra, INT *ksurf, 
				 INT *knsurf, INT *hint ));
pREFIX	pVLOCATE __ProtoGlarp__(( FLOAT *pxyz, INT *kcold, INT *kcnew ));
pREFIX	pVCONNECT __ProtoGlarp__(( INT *kcout, INT *kfout, INT *kcin,
				   INT *idtin ));
pREFIX	pVPOLYHEDRA __ProtoGlarp__(( INT *index, INT *nlnode, INT *ntets,
                                     INT *nfacet ));
pREFIX	pVPGRID __ProtoGlarp__(( INT *index, FLOAT *xyz, INT *tets, INT *dsg,
				 INT *scn ));
pREFIX	pVPSCAL __ProtoGlarp__(( INT *jkey, INT *index, FLOAT *s ));
pREFIX	pVPTHRES __ProtoGlarp__(( INT *jkey, INT *index, FLOAT *xyz,
                                  FLOAT *t ));
pREFIX	pVPVECT __ProtoGlarp__(( INT *jkey, INT *index, FLOAT *v ));
pREFIX	pVZPRIME __ProtoGlarp__(( INT *idcut, FLOAT *xyz,    INT *nnode,
				  FLOAT *zp,  FLOAT *zprime, FLOAT *xpc,
                                  FLOAT *ypc, FLOAT *halfw ));
pREFIX	pVXYPRIME __ProtoGlarp__(( FLOAT *zprime, INT *kn, FLOAT *xyz,
                                   INT *n, FLOAT *xyp ));
pREFIX	pVPZPRIME __ProtoGlarp__(( FLOAT *zprime, INT *index, FLOAT *xyz,
                                   FLOAT *zp ));
pREFIX	pVPXYPRIME __ProtoGlarp__(( FLOAT *zprime, INT *index, FLOAT *xyz,
                                    FLOAT *xyp ));
pREFIX	pVSURF __ProtoGlarp__(( INT *isurf, FLOAT *xpc, FLOAT *ypc,
				FLOAT *halfw ));
pREFIX	pVXYSURF __ProtoGlarp__(( INT *kn, FLOAT *xyz, INT *n, FLOAT *xyp ));
pREFIX	pVPXYSURF __ProtoGlarp__(( INT *index, INT *kn, FLOAT *xyz, INT *n, 
                                   FLOAT *xyp ));
pREFIX	pVSSURF __ProtoGlarp__(( INT *jkey, INT *kn, FLOAT *xyz,
				 INT *n, FLOAT *s ));
pREFIX	pVVSURF __ProtoGlarp__(( INT *jkey, INT *kn, FLOAT *xyz,
				 INT *n, FLOAT *v ));
pREFIX	pVPSSURF __ProtoGlarp__(( INT *jkey, INT *index, INT *kn, FLOAT *xyz,
				  INT *n, FLOAT *s ));
pREFIX	pVPVSURF __ProtoGlarp__(( INT *jkey, INT *index, INT *kn, FLOAT *xyz,
				  INT *n, FLOAT *v ));
pREFIX	pVSTRING __ProtoGlarp__(( char *string ));
pREFIX	pVCATCH __ProtoGlarp__(( char *string ));
pREFIX	pVEXTRACT __ProtoGlarp__(( INT *index, INT *exnum, INT *reqmask,
                                   INT *ival, FLOAT *rvec ));
pREFIX	pVSTATE __ProtoGlarp__(( INT *rank, INT *kn, FLOAT *sv ));
pREFIX	pVPSTATE __ProtoGlarp__(( INT *rank, INT *index, INT *kn, FLOAT *sv ));
pREFIX	pVSSINIT __ProtoGlarp__(( void ));
pREFIX	pVUPDATE __ProtoGlarp__(( INT *stat ));
pREFIX	pVNOTIFY __ProtoGlarp__(( INT *handle, INT *messid ));
pREFIX	pVLOCAL __ProtoGlarp__(( INT *kb, INT *ir, INT *jr, INT *kr ));
pREFIX	pVBOUNCE __ProtoGlarp__(( INT *kc, FLOAT *x0, FLOAT *u0, FLOAT *dt, FLOAT *x ));
pREFI 	pVSETSTRUC __ProtoGlarp__(( INT *type, INT *len, void **ptr ));

#ifdef __cplusplus
}
#endif

#endif  /*_PV3_H_*/
