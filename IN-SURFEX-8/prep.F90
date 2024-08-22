!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!------------------------
PROGRAM PREP
!------------------------
!!
!!    PURPOSE
!!    -------
!!   This program prepares the initial file for offline run
!!
!!    METHOD
!!    ------
!!   
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    P. LeMoigne                  Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     22/04/04
!!
!----------------------------------------------------------------------------
!
USE MODE_POS_SURF
!
USE MODD_SURFEX_OMP, ONLY : NWORK, NWORK2, XWORK, XWORK2, XWORK3, NBLOCKTOT, &
                             NWORK_FULL, NWORK2_FULL, XWORK_FULL, XWORK2_FULL
!
USE MODD_IO_SURF_ASC
USE MODD_IO_SURF_FA
USE MODD_IO_SURF_LFI
USE MODD_IO_SURF_NC
USE MODD_SURF_PAR
USE MODD_SURF_CONF, ONLY : CSOFTWARE
!
USE MODD_SFX_OASIS, ONLY : LOASIS
!
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_READ_ALL_NAMELISTS
USE MODI_GET_LUOUT
!
USE MODI_INIT_PGD_SURF_ATM
USE MODI_IO_BUFF_CLEAN
USE MODI_PREP_SURF_ATM
USE MODI_WRITE_DIAG_SURF_ATM_n
USE MODI_WRITE_HEADER_MNH
USE MODI_WRITE_SURF_ATM_n
!
USE MODI_GET_LONLAT_n
USE MODI_FLAG_UPDATE
USE MODI_ABOR1_SFX
!
USE MODI_SFX_OASIS_INIT
USE MODI_SFX_OASIS_READ_NAM
USE MODI_SFX_OASIS_PREP
USE MODI_SFX_OASIS_END
!
USE MODI_INIT_OUTPUT_NC_n
!
USE MODN_IO_OFFLINE
!------------------------------------------------------------------------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODD_OFF_SURFEX_n
USE MODE_MODELN_SURFEX_HANDLER
!
IMPLICIT NONE
!
#ifndef AIX64
INCLUDE 'omp_lib.h'
#endif
!
#ifdef SFX_MPI
INCLUDE 'mpif.h'
#endif
!
!*    0.     Declaration of local variables
!            ------------------------------
!
INTEGER            :: ILUOUT
INTEGER            :: ILUNAM
INTEGER            :: IYEAR, IMONTH, IDAY
REAL               :: ZTIME
LOGICAL            :: GFOUND

REAL, DIMENSION(0) :: ZZS
 CHARACTER(LEN=28)  :: YATMFILE  ='                            '  ! name of the Atmospheric file
 CHARACTER(LEN=6)   :: YATMFILETYPE ='      '                     ! type of the Atmospheric file
 CHARACTER(LEN=28)  :: YPGDFILE  ='                            '  ! name of the pgd file
 CHARACTER(LEN=6)   :: YPGDFILETYPE ='      '                     ! type of the pgd file
 CHARACTER(LEN=28)  :: YLUOUT    ='LISTING_PREP                '  ! name of listing
!
INTEGER, DIMENSION(11)  :: IDATEF
!
INTEGER :: JNW, INW
INTEGER :: IRET, INB
#ifdef CPLOASIS
INTEGER :: ILOCAL_COMM, INFOMPI, INPROC
#endif
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!
!*    1.      Set default names and parallelized I/O
!             --------------------------------------
!
#ifdef CPLOASIS
!Must be call before DRHOOK !
 CALL SFX_OASIS_INIT(CNAMELIST,ILOCAL_COMM,'PRE')
#else
LOASIS = .FALSE.
#endif
!
IF (LHOOK) CALL DR_HOOK('PREP',0,ZHOOK_HANDLE)
!
#ifdef CPLOASIS
IF(LOASIS)THEN
  CALL MPI_COMM_SIZE(ILOCAL_COMM,INPROC,INFOMPI)
  IF(INPROC>1)THEN
    CALL ABOR1_SFX('PREP: FOR PREP"WITH OASIS ONLY 1 PROC MUST BE USED')
  ENDIF
ENDIF
#endif
!
!    Allocations of Surfex Types
 CALL SURFEX_ALLOC_LIST(1)
!
 CSOFTWARE='PREP'
!
!     1.1     initializations
!             ---------------
!
IYEAR    = NUNDEF
IMONTH   = NUNDEF
IDAY     = NUNDEF
ZTIME    = XUNDEF
!
LPREP    = .TRUE.
!
!     1.2     output listing
!             --------------
 CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
 CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN (UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
!     1.3     output file name read in namelist
!             ---------------------------------
!
 CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
!
 CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
!
 CFILEPGD     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
 CFILEIN      = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')      ! output of PGD program
 CFILEIN_SAVE = CFILEIN
 CFILEOUT     = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
!
 CFILEPGD_FA  = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
 CFILEIN_FA   = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
 CFILEIN_FA_SAVE = CFILEIN_FA
 CFILEOUT_FA  = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
!
 CFILEPGD_LFI = CPGDFILE
 CFILEIN_LFI  = CPGDFILE
 CFILEIN_LFI_SAVE = CFILEIN_LFI
 CFILEOUT_LFI = CPREPFILE
!
 CFILEPGD_NC = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
 CFILEIN_NC  = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
 CFILEIN_NC_SAVE  = CFILEIN_NC
 CFILEOUT_NC = ADJUSTL(ADJUSTR(CPREPFILE)//'.nc')
!
 CALL CLOSE_NAMELIST('ASCII ',ILUNAM)
!
! Reading all namelist (also assimilation)
 YSURF_CUR => YSURF_LIST(1)
 CALL READ_ALL_NAMELISTS(YSURF_CUR, &
                        CSURF_FILETYPE,'PRE',.FALSE.)
!
!*      1.4.   Reads SFX - OASIS coupling namelists
!              ------------------------------------
!
 CALL SFX_OASIS_READ_NAM(CSURF_FILETYPE,XTSTEP_SURF,'PRE')
!
!*      1.5.   Goto model of Surfex Types
!              ---------------------------
!
  ICURRENT_MODEL = 1
!
!*    2.      Preparation of surface physiographic fields
!             -------------------------------------------
!
!$OMP PARALLEL
!$ NBLOCKTOT = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
!
 CALL IO_BUFF_CLEAN
 CALL INIT_PGD_SURF_ATM(YSURF_CUR, &
                              CSURF_FILETYPE,'PRE',YATMFILE,YATMFILETYPE, &
                         IYEAR, IMONTH, IDAY, ZTIME            ) 
!
 CALL IO_BUFF_CLEAN
 CALL PREP_SURF_ATM(YSURF_CUR, &
                    CSURF_FILETYPE,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
!*    3.      Preparation of SFX-OASIS grid, mask, areas files
!             ------------------------------------------------
!
IF(LOASIS)THEN
  CALL SFX_OASIS_PREP(YSURF_CUR%IM%I, YSURF_CUR%UG, YSURF_CUR%U, &
                      CSURF_FILETYPE)
ENDIF
!
!*    4.      Store of surface physiographic fields
!             -------------------------------------
!
 CALL FLAG_UPDATE(YSURF_CUR%IM%DGI, YSURF_CUR%DGU, &
                 .FALSE.,.TRUE.,.FALSE.,.FALSE.)
!
!* opens the file
IF (CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
  LFANOCOMPACT = .TRUE.
  IDATEF(1)=YSURF_CUR%U%TTIME%TDATE%YEAR
  IDATEF(2)=YSURF_CUR%U%TTIME%TDATE%MONTH
  IDATEF(3)=YSURF_CUR%U%TTIME%TDATE%DAY
  IDATEF(4)=NINT(YSURF_CUR%U%TTIME%TIME/3600.) 
  IDATEF(5)=NINT(YSURF_CUR%U%TTIME%TIME/60.) - IDATEF(4) * 60 
  IDATEF(6)=1 
  IDATEF(7:11)=0  
  CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'NEW',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
  CALL FANDAR(IRET,NUNIT_FA,IDATEF)
#endif
END IF
!
LDEF = .TRUE.
!
IF (CSURF_FILETYPE=="NC    ") THEN
  CALL INIT_OUTPUT_NC_n (YSURF_CUR%TM%BDD, YSURF_CUR%CHE, YSURF_CUR%CHN, YSURF_CUR%CHU, &
                         YSURF_CUR%SM%DTS, YSURF_CUR%TM%DTT, YSURF_CUR%DTZ, YSURF_CUR%IM%I, &
                         YSURF_CUR%UG, YSURF_CUR%U, YSURF_CUR%DGU)
ENDIF
!
INW = 1
IF (CSURF_FILETYPE=="NC    ") INW = 2
!
DO JNW = 1,INW
  !
  IF (LWRITE_COORD) CALL GET_LONLAT_n(YSURF_CUR, &
                                      CSURF_FILETYPE)
  !
  !* writes into the file
  CALL IO_BUFF_CLEAN
  !
  ! FLAG_UPDATE now in WRITE_PGD_SURF_ATM_n
  CALL WRITE_SURF_ATM_n(YSURF_CUR, &
                        CSURF_FILETYPE,'PRE',LLAND_USE) !no pgd field
  CALL WRITE_DIAG_SURF_ATM_n(YSURF_CUR, &
                             CSURF_FILETYPE,'ALL')
  !
  LDEF = .FALSE.
  CALL IO_BUFF_CLEAN  
  !
ENDDO
!
!* closes the file
IF (CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
  CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
#endif
END IF
!
!* add informations in the file
IF (CSURF_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH
!
!
!*    4.     Close parallelized I/O
!            ----------------------
!
WRITE(ILUOUT,*) ' '
WRITE(ILUOUT,*) '    -----------------------'
WRITE(ILUOUT,*) '    | PREP ENDS CORRECTLY |'
WRITE(ILUOUT,*) '    -----------------------'
!
WRITE(*,*) ' '
WRITE(*,*) '    -----------------------'
WRITE(*,*) '    | PREP ENDS CORRECTLY |'
WRITE(*,*) '    -----------------------'
!
 CLOSE(ILUOUT)
!
 CALL SURFEX_DEALLO_LIST
!
IF (ASSOCIATED(NWORK)) DEALLOCATE(NWORK)
IF (ASSOCIATED(XWORK)) DEALLOCATE(XWORK)
IF (ASSOCIATED(NWORK2)) DEALLOCATE(NWORK2)
IF (ASSOCIATED(XWORK2)) DEALLOCATE(XWORK2)
IF (ASSOCIATED(XWORK3)) DEALLOCATE(XWORK3)
IF (ASSOCIATED(NWORK_FULL)) DEALLOCATE(NWORK_FULL)
IF (ASSOCIATED(XWORK_FULL)) DEALLOCATE(XWORK_FULL)
IF (ASSOCIATED(NWORK2_FULL)) DEALLOCATE(NWORK2_FULL)
IF (ASSOCIATED(XWORK2_FULL)) DEALLOCATE(XWORK2_FULL)
!
IF (LHOOK) CALL DR_HOOK('PREP',1,ZHOOK_HANDLE)
!
! * OASIS must be finalized after the last DR_HOOK call
!
IF(LOASIS)THEN
  CALL SFX_OASIS_END
ENDIF
!
!-------------------------------------------------------------------------------
!
END PROGRAM PREP
