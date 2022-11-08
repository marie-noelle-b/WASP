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
REAL                         :: ZTHR1,ZTHR2,ZTHR3,ZTHR4,ZLON1,ZLON2
REAL                         :: ZLON3,ZLON4
REAL                         :: ZCHARN1,ZCHARN2,ZCHARN3,ZCHARN4,ZCHARN5
REAL, DIMENSION(1:2)         :: ZCOEFU,ZCOEFA2,ZCOEFB2
REAL, DIMENSION(1:4)         :: ZPOLYU,ZCOEFA,ZCOEFB
REAL, DIMENSION(1:3)         :: ZPOLYU2
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
ZPOLYU2 = (/ 5.062E-3, -1.28E-04, 9.08E-7 /)
ZLIMCHAR = 0.018
ZLIMCHAR2 = 0.002
ZLIMCHAR1 = 0.1
ZTHR1 = 7.0
ZTHR2 = 23.0
ZTHR3 = 25.0
ZTHR4 = 46.0
ZLON1 = 0.5
ZLON2 = 1.0
ZLON3 = 1.0
ZLON4 = 1.0
ZCHARN1 = ZCOEFU(1)*(PWIND**ZCOEFU(2))
ZAA = ZCOEFA(1) + ZCOEFA(2)*PWIND    &
                + ZCOEFA(3)*PWIND**2 &
                + ZCOEFA(4)*PWIND**3
ZBB = ZCOEFB(1) + ZCOEFB(2)*PWIND    &
                + ZCOEFB(3)*PWIND**2 &
                + ZCOEFB(4)*PWIND**3
ZCHARN2 = ZAA * (PWAGE**ZBB)
       !!  write(*,*) PCHARN,PWAGE
ZAA = ZCOEFA2(1) + ZCOEFA2(2)*PWIND
ZBB = ZCOEFB2(1) + ZCOEFB2(2)*PWIND
ZCHARN3 = ZAA * (PWAGE**ZBB)
IF (ZCHARN3 < ZLIMCHAR) THEN
        ZCHARN3 = ZLIMCHAR
ENDIF
ZCHARN4 = ZPOLYU(1) + ZPOLYU(2)*PWIND    &
          + ZPOLYU(3)*PWIND**2 &
          + ZPOLYU(4)*PWIND**3
ZCHARN5 = ZPOLYU2(1) + ZPOLYU2(2)*PWIND    &
          + ZPOLYU2(3)*PWIND**2
PCHARNWA = (ZCHARN1*(0.5*(1.0+TANH((ZTHR1-PWIND)/ZLON1))) & !! al1 g1
         + ZCHARN2*(0.5*(1.0+TANH((PWIND-ZTHR1)/ZLON1)))) & !! + al2 g2      
                  *(0.5*(1.0+TANH((ZTHR2-PWIND)/ZLON2))) &
         + ZCHARN3*(0.5*(1.0+TANH((PWIND-ZTHR2)/ZLON2))) &
                 *(0.5*(1.0+TANH((ZTHR3-PWIND)/ZLON3))) &
        + ZCHARN4*(0.5*(1.0+TANH((PWIND-ZTHR3)/ZLON3))) &
                *(0.5*(1.0+TANH((ZTHR4-PWIND)/ZLON4))) &
       + ZCHARN5*(0.5*(1.0+TANH((PWIND-ZTHR4)/ZLON4)))
IF (PCHARNWA > ZLIMCHAR1) THEN
    PCHARNWA = ZLIMCHAR1
ENDIF
IF (LHOOK) CALL DR_HOOK('MODE_CHARNOCK_WA:CHARNOCK_WAGE',1,ZHOOK_HANDLE)
END FUNCTION CHARNOCK_WAGE
!
END MODULE MODE_CHARNOCK_WA
