/*
 *	pV3 Client Include
 *
 *	Allows OS/machine dependent stuff to work with FORTRAN
 *
 *	Copyright 1994 - 1997, Massachusetts Institute of Technology.
 */

/*
 *      osdepend.h
 *
 *      Note: if using MPI have this include specified after mpi.h
 *
 *
$Log$
Revision 1.1  2002/02/21 23:13:38  helenbrk
Added pV3 and many changes for compilation in LINUX

 */

#ifndef _PV3_H_

#define _PV3_H_


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

#ifndef CRAY
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
#define pV_MPISTOP  pv_mpistop_
#define pVCLIENT    pvclient_
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
#define pVZPRIME    pvzprime_
#define pVXYPRIME   pvxyprime_
#define pVSURF      pvsurf_
#define pVXYSURF    pvxysurf_
#define pVSSURF     pvssurf_
#define pVVSURF     pvvsurf_
#define pVSTRING    pvstring_
#define pVCATCH     pvcatch_
#define pVEXTRACT   pvextract_
#define pVSTATE     pvstate_

#else

#define pV_INIT     pv_init
#define pV_UPDATE   pv_update
#define pV_STAT     pv_stat
#define pV_CONSOLE  pv_console
#define pV_TERMIN   pv_termin
#define pV_GETSTRUC pv_getstruc
#define pV_FIELD    pv_field
#define pV_SENDXI   pv_sendxi
#define pV_SENDXR   pv_sendxr
#define pV_MPISTOP  pv_mpistop
#define pVCLIENT    pvclient
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
#define pVCONNECT   pvconnect
#define pVZPRIME    pvzprime
#define pVXYPRIME   pvxyprime
#define pVSURF      pvsurf
#define pVXYSURF    pvxysurf
#define pVSSURF     pvssurf
#define pVVSURF     pvvsurf
#define pVSTRING    pvstring
#define pVCATCH     pvcatch
#define pVEXTRACT   pvextract
#define pVSTATE     pvstate

#endif

#else

#define pV_INIT     PV_INIT
#define pV_UPDATE   PV_UPDATE
#define pV_STAT     PV_STAT
#define pV_CONSOLE  PV_CONSOLE
#define pV_TERMIN   PV_TERMIN
#define pV_GETSTRUC PV_GETSTRUC
#define pV_FIELD    PV_FIELD
#define pV_SENDXI   PV_SENDXI
#define pV_SENDXR   PV_SENDXR
#define pV_MPISTOP  PV_MPISTOP
#define pVCLIENT    PVCLIENT
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
#define pVZPRIME    PVZPRIME
#define pVXYPRIME   PVXYPRIME
#define pVSURF      PVSURF
#define pVXYSURF    PVXYSURF
#define pVSSURF     PVSSURF
#define pVVSURF     PVVSURF
#define pVSTRING    PVSTRING
#define pVCATCH     PVCATCH
#define pVEXTRACT   PVEXTRACT
#define pVSTATE     PVSTATE

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
#if defined(__STDC__) || defined(__cplusplus)
#define __ProtoGlarp__(x) x
#else
#define __ProtoGlarp__(x) ()
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef CRAY

#ifdef T3E
void	pV_INIT __ProtoGlarp__(( char *title, int title_len, INT *iopt,
				 INT *npgcut, char *tpgcut, int tpgcut_len,
				 INT *nkeys, INT *ikeys, char *tkeys, 
				 int tkeys_len, INT *fkeys, FLOAT *flims,
				 INT *mirror, INT *knode, INT *kequiv,
				 INT *kcel1, INT *kcel2, INT *kcel3,
				 INT *kcel4, INT *knptet, INT *kptet,
				 INT *knblock, INT *blocks, INT *ksurf,
				 INT *knsurf, INT *istat ));
void	pV_CONSOLE __ProtoGlarp__(( char *string, int string_len ));
#else
/* title, tpgcut and tkeys are pointers packed with the string len */
void	pV_INIT __ProtoGlarp__(( int title, int *iopt, int *npgcut,
				 int tpgcut, int *nkeys, int *ikeys,
				 int tkeys, int *fkeys, float *flims,
				 int *mirror, int *knode, int *kequiv,
				 int *kcel1, int *kcel2, int *kcel3,
				 int *kcel4, int *knptet, int *kptet,
				 int *knblock, int *blocks, int *ksurf,
				 int *knsurf, int *istat));
/* string is a char pointer packed with the string length */
void	pV_CONSOLE __ProtoGlarp__(( int string ));
#endif

#else
void	pV_INIT __ProtoGlarp__(( char *title, INT *iopt, INT *npgcut,
				 char *tpgcut, INT *nkeys, INT *ikeys,
				 char *tkeys, INT *fkeys, FLOAT *flims,
				 INT *mirror, INT *knode, INT *kequiv,
				 INT *kcel1, INT *kcel2, INT *kcel3,
				 INT *kcel4, INT *knptet, INT *kptet,
				 INT *knblock, INT *blocks, INT *ksurf,
				 INT *knsurf, INT *istat, int title_len,
				 int tpgcut_len, int tkeys_len ));
void	pV_CONSOLE __ProtoGlarp__(( char *string, int string_len ));
#endif
void	pV_UPDATE __ProtoGlarp__(( FLOAT *time ));
void	pV_STAT  __ProtoGlarp__(( INT *istate ));
void	pV_TERMIN __ProtoGlarp__(( void ));
void	pV_GETSTRUC __ProtoGlarp__(( INT *opt, void **ptr, INT *len ));
void	pV_FIELD __ProtoGlarp__(( INT *ls, INT *lv, INT *lt, INT *lz ));
void	pV_SENDXI __ProtoGlarp__(( INT *index, INT *exnum, INT *subindex,
                                   INT *subsize, INT *len, INT *ivec ));
void	pV_SENDXR __ProtoGlarp__(( INT *index, INT *exnum, INT *subindex,
                                   INT *subsize, INT *len, FLOAT *rvec ));

#ifdef _MPI_INCLUDE
int	pV_MPIStart __ProtoGlarp__(( MPI_Comm cin, int rank, int ncl,
	                             int bsize, MPI_Comm *cout ));
void	pV_MPISTOP __ProtoGlarp__(( void ));
#endif
/*
void	pVCLIENT __ProtoGlarp__(( INT *cid, char *cname, char *dname ));
void	pVCELL __ProtoGlarp__(( INT *cel1, INT *cel2, INT *cel3,
				INT *cel4, INT *nptet, INT *ptet ));
void	pVSURFACE __ProtoGlarp__(( INT *nsurf, INT *scon, INT *scel,
				   char *tsurf ));
void	pVEQUIV __ProtoGlarp__(( INT *listeq ));
void	pVBLANK __ProtoGlarp__(( INT *iblank, INT *tidcon ));
void	pVGRID __ProtoGlarp__(( FLOAT *xyz ));
void	pVSCAL __ProtoGlarp__(( INT *jkey, FLOAT *s ));
void	pVTHRES __ProtoGlarp__(( INT *jkey, FLOAT *xyz, FLOAT *t ));
void	pVVECT __ProtoGlarp__(( INT *jkey, FLOAT *v ));
void	pVSTRUC __ProtoGlarp__(( INT *knode, INT *kequiv,
				 INT *kcel1, INT *kcel2, INT *kcel3,
				 INT *kcel4, INT *knptet, INT *kptet,
				 INT *knblock, INT *blocks, INT *ksurf,
				 INT *knsurf, INT *hint ));
void	pVLOCATE __ProtoGlarp__(( FLOAT *pxyz, INT *kcold, INT *kcnew ));
void	pVCONNECT __ProtoGlarp__(( INT *kcout, INT *kfout, INT *kcin,
				   INT *idtin ));
void	pVZPRIME __ProtoGlarp__(( INT *idcut, FLOAT *xyz, INT *nnode,
				  FLOAT *zprime, FLOAT *xpc, FLOAT *ypc,
				  FLOAT *halfw ));
void	pVXYPRIME __ProtoGlarp__(( FLOAT *zprime, INT *kn, FLOAT *xyz,
                                   INT *n, FLOAT *xyp ));
void	pVSURF __ProtoGlarp__(( INT *isurf, FLOAT *xpc, FLOAT *ypc,
				FLOAT *halfw ));
void	pVXYSURF __ProtoGlarp__(( INT *kn, FLOAT *xyz, INT *n, FLOAT *xyp ));
void	pVSSURF __ProtoGlarp__(( INT *jkey, INT *kn, FLOAT *xyz,
				 INT *n, FLOAT *s ));
void	pVVSURF __ProtoGlarp__(( INT *jkey, INT *kn, FLOAT *xyz,
				 INT *n, FLOAT *v ));
void	pVSTRING __ProtoGlarp__(( char *string ));
void	pVCATCH __ProtoGlarp__(( char *string ));
void	pVEXTRACT __ProtoGlarp__(( INT *index, INT *exnum, INT *reqmask,
                                   INT *ival, FLOAT *rvec ));
void	pVSTATE __ProtoGlarp__(( INT *rank, INT *kn, INT *n, FLOAT *sv ));
*/

#ifdef __cplusplus
}
#endif

#endif  /*_PV3_H_*/
