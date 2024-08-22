!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #################################################################################
SUBROUTINE DEALLOC_SEAFLUX_n (CHS, SG, S)
!     #################################################################################
!
!!****  *DEALLOC_SEAFLUX_n * - Deallocate all arrays
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
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      S. Belamari 03/2014   other _SEA_ variables
!!      S. Senesi   09/2013   introduce sea-ice-cover ans sea-surface salinity
!!      Modified    11/2014 : J. Pianezze : deallocation of wave parameters
!!-------------------------------------------------------------------------
!
!
!
!
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_SEAFLUX_GRID_n, ONLY : SEAFLUX_GRID_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
!
USE MODI_GLTOOLS_DEALLOC
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!*      0.2    declarations of local variables
!
!
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: CHS
TYPE(SEAFLUX_GRID_t), INTENT(INOUT) :: SG
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEALLOC_SEAFLUX_N',0,ZHOOK_HANDLE)
!
IF (ASSOCIATED(S%LCOVER ))   DEALLOCATE(S%LCOVER )
IF (ASSOCIATED(S%XCOVER ))   DEALLOCATE(S%XCOVER )
IF (ASSOCIATED(S%XZS    ))   DEALLOCATE(S%XZS    )
IF (ASSOCIATED(S%XSST   ))   DEALLOCATE(S%XSST   )
IF (ASSOCIATED(S%XSSS   ))   DEALLOCATE(S%XSSS   )
IF (ASSOCIATED(S%XSIC   ))   DEALLOCATE(S%XSIC   )
IF (ASSOCIATED(S%XFSIC  ))   DEALLOCATE(S%XFSIC  )
IF (ASSOCIATED(S%XFSIT  ))   DEALLOCATE(S%XFSIT  )
IF (ASSOCIATED(S%XTP    ))   DEALLOCATE(S%XTP    )
IF (ASSOCIATED(S%XHS    ))   DEALLOCATE(S%XHS    )
IF (ASSOCIATED(S%XZ0    ))   DEALLOCATE(S%XZ0    )
IF (ASSOCIATED(S%XZ0H   ))   DEALLOCATE(S%XZ0H   )
IF (ASSOCIATED(S%XSEABATHY)) DEALLOCATE(S%XSEABATHY)
IF (ASSOCIATED(S%XEMIS  ))   DEALLOCATE(S%XEMIS  )
IF (ASSOCIATED(S%XDIR_ALB))  DEALLOCATE(S%XDIR_ALB)
IF (ASSOCIATED(S%XSCA_ALB))  DEALLOCATE(S%XSCA_ALB)
!
!-------------------------------------------------------------------------------------
!
IF (ASSOCIATED(SG%XGRID_PAR )) DEALLOCATE(SG%XGRID_PAR )
IF (ASSOCIATED(SG%XLAT      )) DEALLOCATE(SG%XLAT      )
IF (ASSOCIATED(SG%XLON      )) DEALLOCATE(SG%XLON      )
IF (ASSOCIATED(SG%XMESH_SIZE)) DEALLOCATE(SG%XMESH_SIZE)
!
!-------------------------------------------------------------------------------------
!
IF(ASSOCIATED(CHS%XDEP))      DEALLOCATE(CHS%XDEP)
IF(ASSOCIATED(CHS%CCH_NAMES)) DEALLOCATE(CHS%CCH_NAMES)
IF(ASSOCIATED(CHS%SVS%CSV))       DEALLOCATE(CHS%SVS%CSV)
!
!-------------------------------------------------------------------------------------
!
IF(ASSOCIATED(S%XCPL_SEA_WIND))      DEALLOCATE(S%XCPL_SEA_WIND)
IF(ASSOCIATED(S%XCPL_SEA_FWSU))      DEALLOCATE(S%XCPL_SEA_FWSU)
IF(ASSOCIATED(S%XCPL_SEA_FWSV))      DEALLOCATE(S%XCPL_SEA_FWSV)
IF(ASSOCIATED(S%XCPL_SEA_SNET))      DEALLOCATE(S%XCPL_SEA_SNET)
IF(ASSOCIATED(S%XCPL_SEA_HEAT))      DEALLOCATE(S%XCPL_SEA_HEAT)
IF(ASSOCIATED(S%XCPL_SEA_EVAP))      DEALLOCATE(S%XCPL_SEA_EVAP)
IF(ASSOCIATED(S%XCPL_SEA_RAIN))      DEALLOCATE(S%XCPL_SEA_RAIN)
IF(ASSOCIATED(S%XCPL_SEA_SNOW))      DEALLOCATE(S%XCPL_SEA_SNOW)
!
!-------------------------------------------------------------------------------------
!
IF (ASSOCIATED(S%TGLT%bat) .AND. S%CSEAICE_SCHEME=='GELATO' ) CALL GLTOOLS_DEALLOC(S%TGLT)
!
IF (LHOOK) CALL DR_HOOK('DEALLOC_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DEALLOC_SEAFLUX_n


