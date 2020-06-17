!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
MODULE MODE_CHARNOCK_WA
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
INTERFACE CHARNOCK_WA
 MODULE PROCEDURE CHARNOCK_WAGE
END INTERFACE
!
CONTAINS
!
!---------------------------------------------------------------------------------------
!
!#######################################################################################
FUNCTION CHARNOCK_WAGE(PWIND,PWAGE) RESULT(PCHARNWA)
!#######################################################################################
!
!****  *CHARNOCK_WAGE*
!
!       PURPOSE
!       -------
! Compute the Charnock parameter using the wave age and the surface wind
! The formulation used for that depends on the surface wind range
!
!! AUTHOR
!! ------
!!
!! MODIFICATIONS
!-------------------------------------------------------------------------------
IMPLICIT NONE
!
!        0.  Declaration
!
!        0.1 declaration of arguments
!
REAL, INTENT(IN) :: PWIND ! wind
REAL, INTENT(IN) :: PWAGE ! wave age
!
REAL :: PCHARNWA ! Charnock parameter
!
!* 0.2 declarations of local variables
!-------------------------------------------------------------------------------
!
REAL                         :: ZLIMCHAR,ZLIMCHAR2,ZAA,ZLIMCHAR1,ZBB
REAL, DIMENSION(1:2)         :: ZCOEFU,ZCOEFA2,ZCOEFB2
REAL, DIMENSION(1:4)         :: ZPOLYU,ZCOEFA,ZCOEFB
REAL(KIND=JPRB)              :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODE_CHARNOCK_WA:CHARNOCK_WAGE',0,ZHOOK_HANDLE)
!
ZCOEFU = (/ 0.70, -2.52 /)
ZCOEFA2 = (/ 2.27, -6.67E-02 /)
ZCOEFB2 = (/ -2.41, 4.30E-02 /)
ZCOEFA = (/ -9.202, 2.265, -0.134, 2.35E-03 /)
ZCOEFB = (/ -0.4124, -0.2225, 0.01178, -1.616E-04 /)
ZPOLYU = (/ 0.0981, -4.13E-03, 4.34E-5, 1.16E-08 /)
ZLIMCHAR = 0.018
ZLIMCHAR2 = 0.002
ZLIMCHAR1 = 0.1
PCHARNWA = ZCOEFU(1)*(PWIND**ZCOEFU(2))
IF (PWIND >= 7.0) THEN
    ZAA = ZCOEFA(1) + ZCOEFA(2)*PWIND    &
                    + ZCOEFA(3)*PWIND**2 &
                    + ZCOEFA(4)*PWIND**3
    ZBB = ZCOEFB(1) + ZCOEFB(2)*PWIND    &
                    + ZCOEFB(3)*PWIND**2 &
                    + ZCOEFB(4)*PWIND**3
    PCHARNWA = ZAA * (PWAGE**ZBB)
       !!  write(*,*) PCHARN,PWAGE
ENDIF
IF (PWIND >= 23.0) THEN
    ZAA = ZCOEFA2(1) + ZCOEFA2(2)*PWIND
    ZBB = ZCOEFB2(1) + ZCOEFB2(2)*PWIND
    PCHARNWA= ZAA * (PWAGE**ZBB)
    IF (PCHARNWA < ZLIMCHAR) THEN
        PCHARNWA = ZLIMCHAR
    ENDIF
ENDIF
IF (PWIND >= 25.0) THEN
    PCHARNWA = ZPOLYU(1) + ZPOLYU(2)*PWIND    &
                       + ZPOLYU(3)*PWIND**2 &
                       + ZPOLYU(4)*PWIND**3
    IF (PCHARNWA < ZLIMCHAR2) THEN
        PCHARNWA = ZLIMCHAR2
    ENDIF
ENDIF
!! Final check, to avoid too large CHAR
IF (PCHARNWA > ZLIMCHAR1) THEN
        PCHARNWA = ZLIMCHAR1
ENDIF
IF (LHOOK) CALL DR_HOOK('MODE_CHARNOCK_WA:CHARNOCK_WAGE',1,ZHOOK_HANDLE)
END FUNCTION CHARNOCK_WAGE
!
END MODULE MODE_CHARNOCK_WA
