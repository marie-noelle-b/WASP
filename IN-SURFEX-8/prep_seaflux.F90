!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE PREP_SEAFLUX (DTCO, UG, U, SM, &
                         HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_SEAFLUX* - prepares variables for SEAFLUX scheme
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     S. Malardel 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      S. Riette   06/2009 PREP_SEAFLUX_SBL has no more argument
!!      Modified    07/2012, P. Le Moigne : CMO1D phasing
!!      Modified    01/2014, S. Senesi : introduce sea-ice model 
!!      Modified    01/2015, R. Séférian : introduce ocean surface albedo 
!!------------------------------------------------------------------
!
!
USE MODD_SURFEX_n, ONLY : SEAFLUX_MODEL_t
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODI_PREP_HOR_SEAFLUX_FIELD
USE MODI_PREP_VER_SEAFLUX
USE MODI_PREP_OUTPUT_GRID
USE MODI_PREP_SEAFLUX_SBL
USE MODI_PREP_SEAICE
USE MODI_GET_LUOUT
!
USE MODN_PREP_SEAFLUX
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
USE MODD_PREP,           ONLY : XZS_LS
USE MODD_SURF_ATM,       ONLY : LVERTSHIFT
USE MODD_SFX_OASIS,      ONLY : LCPL_WAVE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_CLEAN_PREP_OUTPUT_GRID
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SEAFLUX_MODEL_t), INTENT(INOUT) :: SM
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
!
INTEGER :: JMTH,INMTH
INTEGER :: ILUOUT
LOGICAL :: GFOUND         ! Return code when searching namelist
INTEGER :: ILUNAM         ! logical unit of namelist file
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------------
!
!*      0.     Default of configuration
!
!
IF (LHOOK) CALL DR_HOOK('PREP_SEAFLUX',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
 CALL PREP_OUTPUT_GRID(UG, U, &
                       ILUOUT,SM%SG%CGRID,SM%SG%XGRID_PAR,SM%SG%XLAT,SM%SG%XLON)
!
!-------------------------------------------------------------------------------------
!
!*      1.     Read namelist
!
SM%S%LSBL = LSEA_SBL
SM%O%LMERCATOR = LOCEAN_MERCATOR
SM%O%LCURRENT  = LOCEAN_CURRENT
! Relaxation-forcing parameters
SM%OR%XTAU_REL   = XTIME_REL
SM%OR%XQCORR     = XCORFLX
!
SM%OR%LREL_CUR   = LCUR_REL
SM%OR%LREL_TS    = LTS_REL
SM%OR%LFLUX_NULL = LZERO_FLUX
SM%OR%LFLX_CORR  = LCORR_FLUX
SM%OR%LDIAPYCNAL = LDIAPYC
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading and horizontal interpolations
!
!
!*      2.0    Large scale orography
!
 CALL PREP_HOR_SEAFLUX_FIELD(DTCO, UG, U, &
                            SM%DTS, SM%O, SM%OR, SM%SG, SM%S, &
                            HPROGRAM,'ZS     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.1.1    Temperature
!
 CALL PREP_HOR_SEAFLUX_FIELD(DTCO, UG, U, &
                            SM%DTS, SM%O, SM%OR, SM%SG, SM%S, &
                            HPROGRAM,'SST    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.1.2    Salinity
!

 CALL PREP_HOR_SEAFLUX_FIELD(DTCO, UG, U, &
                            SM%DTS, SM%O, SM%OR, SM%SG, SM%S, &
                            HPROGRAM,'SSS    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.1.3   Sea-ice
!
IF (CSEAICE_SCHEME /= 'NONE  ') THEN 
   CALL PREP_SEAICE(UG, &
                    DTCO, SM%DTS, SM%O, SM%OR, SM%SG, SM%S, U, &
                    HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
ENDIF
!
!
!*      2.2    Roughness
!
ALLOCATE(SM%S%XZ0(SIZE(SM%S%XSST)))
SM%S%XZ0 = 0.001
!
ALLOCATE(SM%S%XZ0H(SIZE(SM%S%XSST)))
SM%S%XZ0H = SM%S%XZ0
!
!*      2.3   Ocean Surface Albedo
!
IF(SM%S%CSEA_ALB=='RS14')THEN
  ALLOCATE(SM%S%XDIR_ALB(SIZE(SM%S%XSST)))
  ALLOCATE(SM%S%XSCA_ALB(SIZE(SM%S%XSST)))
  SM%S%XDIR_ALB = 0.065
  SM%S%XSCA_ALB = 0.065
ENDIF
!
!* Peak frequency and significant wave height
!
ALLOCATE(SM%S%XHS(SIZE(SM%S%XSST)))
ALLOCATE(SM%S%XTP(SIZE(SM%S%XSST)))
!
SM%S%XHS=0.0
SM%S%XTP=0.0
!
IF (LCPL_WAVE) THEN
  ALLOCATE(SM%S%XCHARN  (SIZE(SM%S%XSST)))
  SM%S%XCHARN=0.011
ENDIF
!
!-------------------------------------------------------------------------------------
 CALL CLEAN_PREP_OUTPUT_GRID
!-------------------------------------------------------------------------------------
!
!*      3.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_SEAFLUX(SM%S)
ENDIF
!
DEALLOCATE(XZS_LS)
!
!-------------------------------------------------------------------------------------
!
!*      4.     Preparation of optional interpolation of monthly sst
!
SM%S%LINTERPOL_SST=.FALSE.
IF(TRIM(SM%S%CINTERPOL_SST)/='NONE')THEN
!
  SM%S%LINTERPOL_SST=.TRUE.
!
! Precedent, Current, Next, and Second-next Monthly SST
  INMTH=4
!
  ALLOCATE(SM%S%XSST_MTH(SIZE(SM%S%XSST),INMTH))
  DO JMTH=1,INMTH
     SM%S%XSST_MTH(:,JMTH)=SM%S%XSST(:)
  ENDDO
!
ENDIF
!
!-------------------------------------------------------------------------------------
!
!
!*      5.     Optional preparation of interpolation of monthly Sea Surface salinity
!
SM%S%LINTERPOL_SSS=.FALSE.
IF(TRIM(SM%S%CINTERPOL_SSS)/='NONE')THEN
!
   SM%S%LINTERPOL_SSS=.TRUE.
   !
   ! Precedent, Current, Next, and Second-next Monthly SSS
   INMTH=4
   !
   ALLOCATE(SM%S%XSSS_MTH(SIZE(SM%S%XSSS),INMTH))
   DO JMTH=1,INMTH
      SM%S%XSSS_MTH(:,JMTH)=SM%S%XSSS(:)
   ENDDO
   !
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      6.     Preparation of SBL air variables
!
!
IF (SM%S%LSBL) CALL PREP_SEAFLUX_SBL(SM%SG, SM%SSB)
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_SEAFLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_SEAFLUX
