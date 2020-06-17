PROGRAM OFFLINE_WASP
!
!
IMPLICIT NONE
!
!-------------------------------
CHARACTER (LEN=7), SAVE :: CSEA_FLUX='WASPWA ' ! WASP Program Offline
!
CHARACTER(LEN=50) :: outfile
!
REAL,SAVE :: XLVTT  = 2.5008E+6
REAL,SAVE :: XP00= 1.E5
REAL,SAVE :: XRD=6.0221367E+23*1.380658E-23/28.9644E-3
REAL,SAVE :: XRV=6.0221367E+23*1.380658E-23/18.0153E-3
REAL,SAVE :: XCPD=7.* (6.0221367E+23*1.380658E-23 / 28.9644E-3)  /2.
REAL,SAVE :: XCPV
REAL,SAVE :: XCL
REAL      :: XTTS   = 273.16*(1-0.5) + (273.16 - 1.8)*0.5
REAL, DIMENSION(1) :: XZ0=0.001
!
!!REAL, DIMENSION(1)  :: PZREF=5     ! height of T,q forcing                 (m)
!!REAL, DIMENSION(1)  :: PUREF=10     ! height of wind forcing                (m)
REAL, DIMENSION(1)  :: PZREF=2.
REAL, DIMENSION(1)  :: PUREF     !! a lire dans entrée
REAL, DIMENSION(1)  :: PSNOW=0     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(1)  :: PRAIN=0     ! liquid precipitation                  (kg/m2/s)
REAL, DIMENSION(1)  :: UCUR=0     ! zonal current (m/s)
REAL, DIMENSION(1)  :: VCUR=0     ! meridional current (m/s)
REAL, DIMENSION(1)  :: PSW=550.  ! SW downward flux
REAL, DIMENSION(1)  :: PLW=50.   ! LW downward flux
REAL, DIMENSION(1)  :: PCOOL=0.
!in buoys
REAL, DIMENSION(1)  :: XSST
REAL, DIMENSION(1)  :: ZQA       ! specific humidity (g/kg)
REAL, DIMENSION(1)  :: ZWIND,ZRHOA
REAL, DIMENSION(1)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(1)  :: PTWAVE     ! wave peak period
!!REAL  :: PQA       ! air humidity forcing                  (kg/m3)
!->rh relative humidity
REAL, DIMENSION(1) :: PPS,PPA   ! pressure at atmospheric model surface (Pa) and first atmospheric level (Pa)
!
!ouputs
REAL, DIMENSION(1) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(1) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
!REAL :: PSFU      ! zonal momentum flux                   (Pa)
!REAL :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(1) :: PSTRESS ! momentum flux                (Pa)
REAL, DIMENSION(1) :: ZQSAT,ZUSTAR,ZCD,ZCDN,ZCH,ZCE,ZRI,ZRESA_SEA,ZZ0H
REAL, DIMENSION(1) :: PTAU
!
!other
REAL, DIMENSION(1) :: ZEXNS,ZEXNA
!
!-------BUOY--------------
!! INTEGER,PARAMETER :: nbmax=744 ! lion
INTEGER,PARAMETER :: nbmax=2211838 ! Series Lion Met
!INTEGER,PARAMETER :: nbmax=68673 ! Series Lion GC
!INTEGER,PARAMETER :: nbmax=59 ! Cyclones
! 2002-00-01 to 2002-01-31
! 24*31
CHARACTER(LEN=6), DIMENSION(nbmax) :: DATE
REAL, DIMENSION(nbmax) :: TAIR,PMER,FFM,DDM,HUM,SST,TP,REFU,REFP
INTEGER :: jloop
!
!*********************************
!
!0. read value in buoy data
outfile='FLX.'//TRIM(CSEA_FLUX)//'.dat'
OPEN(20, FILE=TRIM(outfile))
!
OPEN(10,FILE='Input.dat')
DO jloop=1,nbmax+1
    READ(10,*) DATE(jloop-1), TAIR(jloop-1),PMER(jloop-1),&
       &FFM(jloop-1), DDM(jloop-1),&
       &HUM(jloop-1),SST(jloop-1),&
       &TP(jloop-1),REFU(jloop-1),&
       &REFP(jloop-1)
!
    PTA(1)=TAIR(jloop-1)+273.15
    PPA(1)=(PMER(jloop-1))*100.
    PPS(1)=PPA(1)
    ZWIND(1)=FFM(jloop-1)
    PTWAVE(1)=TP(jloop-1)
    PUREF(1)=REFU(jloop-1)
!!    CALL HUMR2HUMS(PTA(1),PMER(jloop-1),HUM(jloop-1),ZQA(1))
    ZQA(1)=HUM(jloop-1)
    XSST(1)=SST(jloop-1)+273.15
    ZRHOA(1)  = PPA(1) / XRD / PTA(1) / ( 1.+((XRV/XRD)-1.)*ZQA(1) )
    ZEXNS(1)     = (PPS(1)/XP00)**(XRD/XCPD)
    ZEXNA(1)     = (PPA(1)/XP00)**(XRD/XCPD)
!
!1. compute turbulent fluxes 
!
! Parameterization
!
!dbg!
      CALL WASP_FLUX(XZ0,                          &
           PTA,ZEXNA,ZRHOA,XSST,ZEXNS,ZQA,ZWIND,   &
           PZREF,PUREF,                            &
           PPS,ZQSAT,                              &
           PSFTH, PSFTQ, ZUSTAR,                   &
           ZCD, ZCDN, ZCH, ZCE, ZRI, ZRESA_SEA,    &
           PRAIN, ZZ0H, PTWAVE, PSW, PLW, PCOOL  )
 
!! Momentum fluxes
!
!2. write fluxes
      PSFTQ(1)=(XLVTT + (XCPV-XCL)*(XSST(1)-XTTS))*PSFTQ(1)
      PTAU=ZRHOA*ZUSTAR*ZUSTAR
!! Pour calcul age
!!      PTWAVE=9.81*PTWAVE/2/3.14159/ZUSTAR
!!      WRITE(20,*) jloop-1, DATE(jloop-1), PSFTH, PSFTQ, PSTRESS
      WRITE(20,*) DATE(jloop-1), ZWIND, ZUSTAR, PTAU, ZCDN, ZCD,PSFTH, PSFTQ, ZQA
ENDDO
CLOSE(10)
!---------------------------
CONTAINS
!*******************************
SUBROUTINE HUMR2HUMS(T,pl,rh,sh)
!
  REAL, INTENT(IN) :: T
  REAL, INTENT(IN) :: pl
  REAL, INTENT(IN) :: rh
  REAL, INTENT(OUT) :: sh
!
  REAL :: eilog,eilog2,eilog3,eilog4,es,ws
!---------------
! tension de vapeur saturante
  if (T-273.15<0) then ! par rapport a la glace 
    eilog = -9.09718 * ((273.16 / T) -1)
    eilog2 = -3.5654 * log10 (273.16 / T)
    eilog3 = 0.876793 * (1 - (T / 273.16 ))
    es=6.1071 * 10**(eilog+eilog2+eilog3)

  else ! par rapport à l'eau 
    eilog=-7.90298 * ((373.16 / T) - 1)
    eilog2=5.02808 * log10((373.16 / T))
    eilog3=-1.3816e-7 * ((10**(11.344*(1-(T/373.16)))) -1)
    eilog4=8.1328e-3 * ((10**( -3.49149 * ((373.16 / T) - 1) )) -1)
    es=1013.246 * 10**(eilog+eilog2+eilog3+eilog4)
  endif
! masse
  ws = 0.62197 * es / (pl - 0.378*es)
! humidité spécifique
  sh = (rh/100.0) * ws
!!  sh=sh.*1000; %%% passage en g.kg-1
!
END SUBROUTINE
!

!********************************************************
!     #########
    SUBROUTINE WASP_FLUX(PZ0SEA,PTA,PEXNA,PRHOA,PSST,PEXNS,PQA,  &
            PVMOD,PZREF,PUREF,PPS,PQSAT,PSFTH,PSFTQ,PUSTAR,PCD,PCDN,PCH,PCE,PRI,&
            PRESA,PRAIN,PZ0HSEA, PTP, PSW, PLW, PCOOL)  
!     #######################################################################
!
!
!!****  *WASP_FLUX*  
!!
!!    PURPOSE
!!    -------
!      Calculate the surface fluxes of heat, moisture, and momentum over
!      sea surface with a new bulk algorithm based on the wave-age dependency
!      of the Charnock parameter (alpha) in the surface wind range 7 -- 23 m/s
!      This algorithm mimics the Cd values obtained by Pineau-Guillou in
!      modifying the wave-wind coupling implemented in IFS (Janssen, 1989 and
!      following, see also Bidlot)
!      The purpose is to reproduce the Cd variability wrt wage age without
!      making use of the full wave spectrum (esp. the HF part) which is subject
!      to caution in the wave models. The results (Cd, u*) is fitted to
!      observations in the surface wind range 0 -- 30 m/s, as Ecume V8 or Coare
!      3.5 are.
!     
!!**  METHOD
!!    ------
!      For momentum: the iteration is made on the scale parameter u*, the
!      surface roughness length z0 and the Charnock parameter alpha (+ the
!      Obukhov length which is necessary for the stability function)
!      The method is close to the one used by Fairall et al (2003) in the Coare
!      3.0 algorithm. Especially, we use the same stability functions and take
!      advantage of the Richardson number as a proxy of the stability to obtain
!      precise first-guess and limit the number of iterations.
!      The drag coefficient Cd is computed a posteriori (diagnostic).
!      For heat and moisture: the heat (moisture) roughness length z0t (z0q) is 
!      used instead of the transfer parameter Ch (Ce). It is derived from the
!      balance between z0, Cd and Ch to obtain values of Ch corresponding to
!      observations. As heat flux datasets are not consistent with each other
!      between different experiments, we use here a constant value of 1.2 E-3
!      for Chn. This can be changed easily using the Coare 3.0 formula for z0t
!      namely z0t = min (1.15E-4, 5.5E-5/(r**0.6)); in such a case, Chn will
!      increase up to 1.5E-3 for strong winds (15 to 20 m/s) then get back to
!      1.2E-3 for cyclonic winds.
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!     MN Bouin (Meteo-France) - L. Pineau-Guillou (Ifremer)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     10/03/2017
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!USE MODD_CSTS,       ONLY : XKARMAN, XG, XSTEFAN, XRD, XRV, XPI, &
!                              XLVTT, XCL, XCPD, XCPV, XTT,XP00  
!USE MODD_SEAFLUX_n
!USE MODD_SURF_PAR,   ONLY : XUNDEF
!USE MODD_WATER_PAR
!USE MODD_SURF_ATM, ONLY : XRIMAX
!
!USE MODI_SURFACE_RI
!USE MODI_WIND_THRESHOLD
!USE MODE_COARE30_PSI
!
!USE MODE_THERMOS
!
!
!
!USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(IN)       :: PTA   ! air temperature at atm. level (K)
REAL, DIMENSION(:), INTENT(IN)       :: PQA   ! air humidity at atm. level (kg/kg)
REAL, DIMENSION(:), INTENT(IN)       :: PEXNA ! Exner function at atm. level
REAL, DIMENSION(:), INTENT(IN)       :: PRHOA ! air density at atm. level
REAL, DIMENSION(:), INTENT(INOUT)       :: PVMOD ! module of wind at atm. wind level (m/s)
REAL, DIMENSION(:), INTENT(IN)       :: PZREF ! atm. level for temp. and humidity (m)
REAL, DIMENSION(:), INTENT(IN)       :: PUREF ! atm. level for wind (m)
REAL, DIMENSION(:), INTENT(IN)       :: PSST  ! Sea Surface Temperature (K)
REAL, DIMENSION(:), INTENT(IN)       :: PEXNS ! Exner function at sea surface
REAL, DIMENSION(:), INTENT(IN)       :: PPS   ! air pressure at sea surface (Pa)
REAL, DIMENSION(:), INTENT(IN)       :: PRAIN !precipitation rate (kg/s/m2)
REAL, DIMENSION(:), INTENT(IN)       :: PSW   ! SW flux
REAL, DIMENSION(:), INTENT(IN)       :: PLW   ! LW flux
REAL, DIMENSION(:), INTENT(IN)       :: PCOOL   ! ZCOOL
!
REAL, DIMENSION(:), INTENT(INOUT)    :: PZ0SEA! roughness length over the ocean
!                                                                                 
!  surface fluxes : latent heat, sensible heat, friction fluxes
REAL, DIMENSION(:), INTENT(OUT)      :: PSFTH ! heat flux (W/m2)
REAL, DIMENSION(:), INTENT(OUT)      :: PSFTQ ! water flux (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)      :: PUSTAR! friction velocity (m/s)
!
! diagnostics
REAL, DIMENSION(:), INTENT(OUT)      :: PQSAT ! humidity at saturation
REAL, DIMENSION(:), INTENT(OUT)      :: PCD   ! heat drag coefficient
REAL, DIMENSION(:), INTENT(OUT)      :: PCDN  ! momentum drag coefficient
REAL, DIMENSION(:), INTENT(OUT)      :: PCH   ! neutral momentum drag coefficient
REAL, DIMENSION(:), INTENT(OUT)      :: PCE  !transfer coef. for latent heat flux
REAL, DIMENSION(:), INTENT(OUT)      :: PRI   ! Richardson number
REAL, DIMENSION(:), INTENT(OUT)      :: PRESA ! aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)      :: PZ0HSEA ! heat roughness length
REAL, DIMENSION(:), INTENT(OUT)      :: PTP ! wave peak period
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA))      :: ZVMOD    ! wind intensity
REAL, DIMENSION(SIZE(PTA))      :: ZPA      ! Pressure at atm. level
REAL, DIMENSION(SIZE(PTA))      :: ZQASAT   ! specific humidity at saturation  at atm. level (kg/kg)
!
REAL, DIMENSION(SIZE(PTA))      :: ZO       ! rougness length ref 
REAL, DIMENSION(SIZE(PTA))      :: ZWG      ! gustiness factor (m/s)
!
REAL, DIMENSION(SIZE(PTA))      :: ZDU,ZDT,ZDQ,ZDUWG !differences
!
REAL, DIMENSION(SIZE(PTA))      :: ZUSR        !velocity scaling parameter "ustar" (m/s) = friction velocity
REAL, DIMENSION(SIZE(PTA))      :: ZTSR        !temperature sacling parameter "tstar" (degC)
REAL, DIMENSION(SIZE(PTA))      :: ZQSR        !humidity scaling parameter "qstar" (kg/kg)
!
REAL, DIMENSION(SIZE(PTA))      :: ZU10        !vertical profils (10-m height) 
REAL, DIMENSION(SIZE(PTA))      :: ZVISA       !kinematic viscosity of dry air
REAL, DIMENSION(SIZE(PTA))      :: ZO10,ZOT10  !roughness length at 10m
REAL, DIMENSION(SIZE(PTA))      :: ZCD,ZCT,ZCC
REAL, DIMENSION(SIZE(PTA))      :: ZCD10,ZCT10 !transfer coef. at 10m
REAL, DIMENSION(SIZE(PTA))      :: ZRIBU,ZRIBCU
REAL, DIMENSION(SIZE(PTA))      :: ZETU,ZL10
!
REAL, DIMENSION(SIZE(PTA))      :: ZCHARN                      !Charnock number depends on wind module
REAL, DIMENSION(SIZE(PTA))      :: ZTWAVE,ZCWAVE      !to compute gravity waves' impact
REAL, DIMENSION(SIZE(PTA))      :: ZWAGE              !Wave age = cp/usr
REAL, DIMENSION(SIZE(PTA))      :: ZCOOL
!
REAL, DIMENSION(SIZE(PTA))      :: ZZL,ZZTL!,ZZQL    !Obukhovs stability 
                                                     !param. z/l for u,T,q
REAL, DIMENSION(SIZE(PTA))      :: ZRR
REAL, DIMENSION(SIZE(PTA))      :: ZOT,ZOQ           !rougness length ref
REAL, DIMENSION(SIZE(PTA))      :: ZPUZ,ZPTZ,ZPQZ    !PHI funct. for u,T,q 
!
REAL, DIMENSION(SIZE(PTA))      :: ZBF               !constants to compute gustiness factor
!
REAL, DIMENSION(SIZE(PTA))      :: ZTAU       !momentum flux (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZHF        !sensible heat flux (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZEF        !latent heat flux (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZWBAR      !diag for webb correction but not used here after
REAL, DIMENSION(SIZE(PTA))      :: ZTAUR      !momentum flux due to rain (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZRF        !sensible heat flux due to rain (W/m2)
REAL, DIMENSION(SIZE(PTA))      :: ZCHN,ZCEN  !neutral coef. for heat and vapor
!
REAL, DIMENSION(SIZE(PTA))      :: ZLV      !latent heat constant
!
REAL, DIMENSION(SIZE(PTA))      :: ZTAC,ZDQSDT,ZDTMP,ZDWAT,ZALFAC ! for precipitation impact
REAL, DIMENSION(SIZE(PTA))      :: ZXLR                           ! vaporisation  heat  at a given temperature
REAL, DIMENSION(SIZE(PTA))      :: ZCPLW                          ! specific heat for water at a given temperature 
!
REAL, DIMENSION(SIZE(PTA))      :: ZUSTAR2  ! square of friction velocity
!
REAL, DIMENSION(SIZE(PTA))      :: ZDIRCOSZW! orography slope cosine (=1 on water!)
REAL, DIMENSION(SIZE(PTA))      :: ZAC      ! Aerodynamical conductance
REAL, DIMENSION(SIZE(PTA))      :: ZAL, ZBIGC, ZWETC, ZDTCS, ZDQCS, ZTKT
REAL, DIMENSION(SIZE(PTA))      :: ZLAMB, ZDELS, ZQOUT, ZCOL, ZRNL, ZRNS
REAL, DIMENSION(SIZE(PTA))      :: ZALQ, ZHFI, ZEFI, ZRL
!
!
INTEGER, DIMENSION(SIZE(PTA))   :: ITERMAX             ! maximum number of iterations
!
REAL    :: ZRVSRDM1,ZRDSRV ! thermodynamic constants
REAL    :: ZBETAGUST           !gustiness factor
REAL    :: ZZBL                !atm. boundary layer depth (m)
REAL    :: ZS                  !height of rougness length ref
REAL    :: ZCH10               !transfer coef. at 10m
!
INTEGER :: J, JLOOP    !loop indice
!REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-----csts---------
REAL :: XKARMAN, XG, XRD, XRV, XLVTT, XCL, XCPD, XCPV, XTT,XP00
REAL :: XUNDEF, XPI, XRIMAX, XVZ0CM, XALBEDO
REAL :: XGAMW,XBETAW,XALPW, XBE, XVISW, XRHOW
REAL, DIMENSION(SIZE(PTA))        :: ZFOES, ZFOESA
!-----namel--------
LOGICAL :: LPWG=.FALSE. 
LOGICAL :: LPRECIP=.FALSE.
!LOGICAL :: LPWEBB=.FALSE.
!
!-------------------------------------------------------------------------------
!
 XVZ0CM=0.0
 XP00 = 1.E5
 XTT    = 273.16
 XCL    = 4.218E+3
 XLVTT  = 2.5008E+6
 XRD=6.0221367E+23*1.380658E-23/28.9644E-3
 XRV=6.0221367E+23*1.380658E-23/18.0153E-3
 XPI = 2.*ASIN(1.)
 XCPD = 7.* XRD /2.
 XCPV   = 4.* XRV
 XG = 9.80665
 XKARMAN = 0.4
 XGAMW = (4.218E+3 - 4.* XRV) / XRV
 XBETAW = (2.5008E+6/XRV) + (XGAMW * 273.16)
 XALPW = LOG(611.14) + (XBETAW /273.16) + (XGAMW *LOG(273.16))
 XUNDEF=-999.
 XRIMAX=0.20
 XBE=0.026
 XRHOW=1022.
 XVISW=0.000001
 XALBEDO=0.06
!       1.     Initializations
!              ---------------
!
!       1.1   Constants and parameters
!
!
ZRVSRDM1  = XRV/XRD-1. ! 0.607766
ZRDSRV    = XRD/XRV    ! 0.62198
ZBETAGUST = 1.2        ! value based on TOGA-COARE experiment
ZZBL      = 600.       ! Set a default value for boundary layer depth
ZS        = 10.        ! Standard heigth =10m
ZCH10     = 0.00115
!
!
!       1.2   Array initialization by undefined values
!
PSFTH (:)=XUNDEF
PSFTQ (:)=XUNDEF
PUSTAR(:)=XUNDEF
!
PCD(:) = XUNDEF
PCDN(:) = XUNDEF
PCH(:) = XUNDEF
PCE(:) =XUNDEF
PRI(:) = XUNDEF
!
PRESA(:)=XUNDEF
!
!-------------------------------------------------------------------------------
!       2. INITIAL GUESS FOR THE ITERATIVE METHOD 
!          -------------------------------------
!
!       2.1     Wind and humidity 
!
! Sea surface specific humidity 
!
ZFOES(:) = 0.98*EXP( XALPW - XBETAW/PSST(:) - XGAMW*LOG(PSST(:)) )
PQSAT(:) = XRD/XRV*ZFOES(:)/PPS(:) / (1.+(XRD/XRV-1.)*ZFOES(:)/PPS(:))      
!              
! Set a minimum value to wind 
!
ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
PVMOD(:) = ZVMOD(:)
ZCOOL(:) = PCOOL(:)
!
! Specific humidity at saturation at the atm. level 
!
ZPA(:) = XP00* (PEXNA(:)**(XCPD/XRD))
ZFOESA(:) = EXP( XALPW - XBETAW/PTA(:) - XGAMW*LOG(PTA(:)) )
ZQASAT(:) = XRD/XRV*ZFOESA(:)/ZPA(:) / (1.+(XRD/XRV-1.)*ZFOESA(:)/ZPA(:))
ZDTCS(:) = 0.3

!
!
ZO(:)  = 0.0001
ZWG(:) = 0.
ZTKT(:) = 0.001
!
IF (LPWG) ZWG(:) = 0.5
!
DO J=1,SIZE(PTA)
  !
  !      2.2       initial guess
  !    
  ZDU(J) = ZVMOD(J)   !wind speed difference with surface current(=0) (m/s)
                      !initial guess for gustiness factor
  ZLV(J) = XLVTT + (XCPV-XCL)*(PSST(J)-XTT)
  ZWETC(J) =  ZRDSRV*XLVTT*PQSAT(J)/(XRD*(PSST(J)**2))
  ZDQCS(J) = ZWETC(J)*ZDTCS(J)
  ZDT(J) = -(PTA(J)/PEXNA(J)) + (PSST(J)/PEXNS(J)) !potential temperature difference
  ZDQ(J) = PQSAT(J)-PQA(J)                         !specific humidity difference
  ZAL(J) = 2.1E-5*(PSST(J)-XTT+3.2)**0.79
  ZBIGC(J) = 16*XG*ZCPW*(ZVISW*XRHOSW)**3/(ZRVSRDM1*ZRVSRDM1*PRHOA(J)*PRHOA(J))
  ZRL(J) = PLW(J)
  ZRNS(J) = PSW(J) * (1-XALBEDO)
  
  !
  ZDUWG(J) = (ZDU(J)*ZDU(J)+ZWG(J)*ZWG(J))**0.5    !wind speed difference including gustiness ZWG
  !
  !      2.3   initialization of neutral coefficients
  !
  ZU10(J)  = ZDUWG(J)*LOG(ZS/ZO(J))/LOG(PUREF(J)/ZO(J))
  ZUSR(J)  = 0.035*ZU10(J)
  ZVISA(J) = 1.326E-5*(1.+6.542E-3*(PTA(J)-XTT)+&
             8.301E-6*(PTA(J)-XTT)**2.-4.84E-9*(PTA(J)-XTT)**3.) !Andrea (1989) CRREL Rep. 89-11
  ! 
  ZO10(J) = 0.011*ZUSR(J)*ZUSR(J)/XG+0.11*ZVISA(J)/ZUSR(J)
  ZCD(J)  = (XKARMAN/LOG(PUREF(J)/ZO10(J)))**2  !drag coefficient
  ZCD10(J)= (XKARMAN/LOG(ZS/ZO10(J)))**2
  ZCT10(J)= ZCH10/SQRT(ZCD10(J))
  ZOT10(J)= ZS/EXP(XKARMAN/ZCT10(J))
  !
  !-------------------------------------------------------------------------------
  !             Grachev and Fairall (JAM, 1997)
  ZCT(J) = XKARMAN/LOG(PZREF(J)/ZOT10(J))      !temperature transfer coefficient
  ZCC(J) = XKARMAN*ZCT(J)/ZCD(J)               !z/L vs Rib linear coef.
  !
  ZRIBCU(J) = -PUREF(J)/(ZZBL*0.004*ZBETAGUST**3) !saturation or plateau Rib
  ZRIBU(J)  = -XG*PUREF(J)*((ZDT(J)-ZDTCS(J)*ZCOOL(J))+ZRVSRDM1*PTA(J)*ZDQ(J))/&
               (PTA(J)*ZDUWG(J)**2)  
  !
  IF (ZRIBU(J)<0.) THEN
    ZETU(J) = ZCC(J)*ZRIBU(J)/(1.+ZRIBU(J)/ZRIBCU(J))    !Unstable G and F
  ELSE
    ZETU(J) = ZCC(J)*ZRIBU(J)/(1.+27./9.*ZRIBU(J)/ZCC(J))!Stable 
  ENDIF
  !
  ZL10(J) = PUREF(J)/ZETU(J) !MO length
  !
ENDDO
!
!  First guess M-O stability dependent scaling params. (u*,T*,q*) to estimate ZO and z/L (ZZL)
ZUSR(:) = ZDUWG(:)*XKARMAN/(LOG(PUREF(:)/ZO10(:))-PSIFCTU(PUREF(:)/ZL10(:)))
ZTSR(:) = -(ZDT(:)-ZCOOL(:)*ZDTCS(:))*XKARMAN/(LOG(PZREF(:)/ZOT10(:))-PSIFCTT(PZREF(:)/ZL10(:)))
ZQSR(:) = -(ZDQ(:)-ZWETC(:)*ZDTCS(:)*ZCOOL(:))*XKARMAN/(LOG(PZREF(:)/ZOT10(:))-PSIFCTT(PZREF(:)/ZL10(:)))
!
ZCHARN(:) = 0.011  
!
ZZL(:) = 0.0
ZHFI(:) = 0.0
ZEFI(:) = 0.0
!
DO J=1,SIZE(PTA)
  !
  IF (ZETU(J)>50.) THEN
    ITERMAX(J) = 1
  ELSE
    ITERMAX(J) = 3 !number of iterations
  ENDIF
  !
  !                3.  ITERATIVE LOOP TO COMPUTE USR, TSR, QSR 
  !                -------------------------------------------
  !
!!  ZTWAVE(J) = 0.729*PVMOD(J)
!!  ZTWAVE(J) = 0.5*PVMOD(J)
  ZTWAVE(J) = PTP(J)
  ZCWAVE(J) = XG*ZTWAVE(J)/(2.*XPI)
  ZWAGE(J) = ZCWAVE(J)/ZUSR(J)
!!  write(*,*) ZU10(J), ZUSR(J), ZOT10(J), ZDT(J), ZDQ(J), ZWAGE(J)
  !
ENDDO
!
   
!
DO JLOOP=1,MAXVAL(ITERMAX) ! begin of iterative loop
  !
  DO J=1,SIZE(PTA)
    !
    IF (JLOOP.GT.ITERMAX(J)) CYCLE
    !
    ZCHARN(J) = CHARNOCK_WA(ZVMOD(J), ZWAGE(J))
!!    write(*,*) ZVMOD(J), ZUSR(J), ZWAGE(J), ZCHARN(J)
    ZO(J) = ZCHARN(J)*ZUSR(J)*ZUSR(J)/XG + 0.11*ZVISA(J)/ZUSR(J) !Smith 1988
    !
    ZRR(J) = ZO(J)*ZUSR(J)/ZVISA(J)
!!    ZOQ(J) = MIN(1.15E-4 , 5.5E-5/ZRR(J)**0.6)
    ZOT(J) = PZREF(J) * EXP(-(XKARMAN**2)/(ZCH10*LOG(PUREF(J)/ZO(J))))
!!    write(*,*) PZREF(J),PUREF(J),ZO(J),ZOT(J)
    ZOQ(J) = ZOT(J)
    !
    ZZL(J) = XKARMAN * XG * PUREF(J) * &
              ( ZTSR(J)*(1.+ZRVSRDM1*PQA(J)) + ZRVSRDM1*PTA(J)*ZQSR(J) ) / &
              ( PTA(J)*ZUSR(J)*ZUSR(J)*(1.+ZRVSRDM1*PQA(J)) )  
    ZZTL(J)= ZZL(J)*PZREF(J)/PUREF(J)  ! for T 
!    ZZQL(J)=ZZL(J)*PZREF(J)/PUREF(J)  ! for Q
  ENDDO
  !
  ZPUZ(:) = PSIFCTU(ZZL(:))     
  ZPTZ(:) = PSIFCTT(ZZTL(:))
  !
  DO J=1,SIZE(PTA)
    !
    ! ZPQZ(J)=PSIFCTT(ZZQL(J))    
    ZPQZ(J) = ZPTZ(J)
!!    write(*,*) ZZL(J),ZPUZ(J), ZPTZ(J)
    !
    !             3.1 scale parameters
    !
    ZUSR(J) = ZDUWG(J)*XKARMAN/(LOG(PUREF(J)/ZO(J)) -ZPUZ(J))
    ZTSR(J) = -(ZDT(J)-ZDTCS(J)*ZCOOL(J))  *XKARMAN/(LOG(PZREF(J)/ZOT(J))-ZPTZ(J))
    ZQSR(J) = -(ZDQ(J)-ZWETC(J)*ZDTCS(J)*ZCOOL(J))  *XKARMAN/(LOG(PZREF(J)/ZOQ(J))-ZPQZ(J))
!!    write(*,*) ZUSR(J),ZTSR(J),ZQSR(J)
    !
    !             3.2 Gustiness factor (ZWG)
    !
    IF(LPWG) THEN
      ZBF(J) = -XG/PTA(J)*ZUSR(J)*(ZTSR(J)+ZRVSRDM1*PTA(J)*ZQSR(J))
      IF (ZBF(J)>0.) THEN
        ZWG(J) = ZBETAGUST*(ZBF(J)*ZZBL)**(1./3.)
      ELSE
        ZWG(J) = 0.2
      ENDIF
    ENDIF  
    ZDUWG(J) = SQRT(ZVMOD(J)**2. + ZWG(J)**2.)
    !
    !           3.3 Cool skin parameters
    !
    ZRNL(J) = 0.97*(5.67E-8*(PSST(J)-ZDTCS(J)*ZCOOL(J))**4-ZRL(J))
    ZHFI(J)  =  -PRHOA(J)*XCPD*ZUSR(J)*ZTSR(J)
    ZEFI(J)  =  -PRHOA(J)*ZLV(J)*ZUSR(J)*ZQSR(J)
    ZQOUT(J) = ZRNL(J) + ZHFI(J) + ZEFI(J)
    ZDELS(J) = ZRNS(J) * (0.065 + 11 * ZTKT(J) - 6.6E-5 / ZTKT(J) * (1 - EXP(-ZTKT(J)/8.0E-4)))
    ZCOL(J) = ZQOUT(J) - ZDELS(J)
    ZALQ(J) = ZAL(J) * ZCOL(J) + XBE * ZEFI(J) * ZCPW / ZLV(J)
    !
    IF (ZALQ(J) .GT. 0) THEN
      ZLAMB(J) = 6/(1 + (ZBIGC(J) * ZALQ(J) / ZUSR(J)**4)**0.75)**0.333  ! Eq 13 Saunders
      ZTKT(J) = ZLAMB(J) * ZVISW / (SQRT(PRHOA(J)/XRHOSW) * ZUSR(J))  ! Eq 11 sub thk
    ELSE
      ZLAMB(J) = 6.0
      ZTKT(J) = MIN(0.01,ZLAMB(J)*ZVISW/(SQRT(PRHOA(J)/XRHOSW) * ZUSR(J)))
    ENDIF
    !
    ZDTCS(J) = ZCOL(J) * ZTKT(J) / ZRVSRDM1
    ZDQCS(J) = ZWETC(J) * ZDTCS(J)
  ENDDO
  !
ENDDO
!-------------------------------------------------------------------------------
!
!            4.  COMPUTE transfer coefficients PCD, PCH, ZCE and SURFACE FLUXES
!                --------------------------------------------------------------
!
ZTAU(:) = XUNDEF
ZHF(:)  = XUNDEF
ZEF(:)  = XUNDEF
!
ZWBAR(:) = 0.
ZTAUR(:) = 0.
ZRF(:)   = 0.
!
DO J=1,SIZE(PTA)
  !
  !            4. 0  roughness length over the ocean PZOSEA
  !
  !
  !            4. 1 transfert coefficients PCD, PCH and PCE 
  !                 and neutral PCDN, ZCHN, ZCEN 
  !
  PCD(J) = (ZUSR(J)/ZDUWG(J))**2.
  PCH(J) = ZUSR(J)*ZTSR(J)/(ZDUWG(J)*(PTA(J)*PEXNS(J)/PEXNA(J)-PSST(J)+ZDTCS(J)*ZCOOL(J)))
  PCE(J) = ZUSR(J)*ZQSR(J)/(ZDUWG(J)*(PQA(J)-PQSAT(J)+ZDQCS(J)*ZCOOL(J)))
  !-ZDQCS(J)*ZCOO
  PCDN(J) = (XKARMAN/LOG(ZS/ZO(J)))**2.
  ZCHN(J) = (XKARMAN/LOG(ZS/ZO(J)))*(XKARMAN/LOG(ZS/ZOT(J)))
  ZCEN(J) = (XKARMAN/LOG(ZS/ZO(J)))*(XKARMAN/LOG(ZS/ZOQ(J)))
!
  PZ0SEA(:) = ZCHARN(J) * ZUSR(J) * ZUSR(J) / XG + XVZ0CM * PCD(J) / PCDN(J)  
  !
  ZLV(J) = XLVTT + (XCPV-XCL)*(PSST(J)-XTT)
  !
  !            4. 2 surface fluxes 
  !
  IF (ABS(PCDN(J))<=1.E-2) THEN   !!!! secure COARE3.0 CODE 
    ZTSR(J) = -ZTSR(J)
    ZQSR(J) = -ZQSR(J)
    ZTAU(J) = -PRHOA(J)*ZUSR(J)*ZUSR(J)*ZVMOD(J)/ZDUWG(J)
    ZHF(J)  =  PRHOA(J)*XCPD*ZUSR(J)*ZTSR(J)
    ZEF(J)  =  PRHOA(J)*ZLV(J)*ZUSR(J)*ZQSR(J)
!!    write (*,*) ZHF(J),ZEF(J)
    !    
    !           4.3 Contributions to surface  fluxes due to rainfall
    !
    ! SB: a priori, le facteur ZRDSRV=XRD/XRV est introduit pour
    !     adapter la formule de Clausius-Clapeyron (pour l'air
    !     sec) au cas humide.
    IF (LPRECIP) THEN
      ! 
      ! heat surface  fluxes
      !
      ZTAC(J)  = PTA(J)-XTT
      !
      ZXLR(J)  = XLVTT + (XCPV-XCL)* ZTAC(J)                            ! latent heat of rain vaporization
      ZDQSDT(J)= ZQASAT(J) * ZXLR(J) / (XRD*PTA(J)**2)                  ! Clausius-Clapeyron relation
      ZDTMP(J) = (1.0 + 3.309e-3*ZTAC(J) -1.44e-6*ZTAC(J)*ZTAC(J)) * &  !heat diffusivity
                  0.02411 / (PRHOA(J)*XCPD)
      !
      ZDWAT(J) = 2.11e-5 * (XP00/ZPA(J)) * (PTA(J)/XTT)**1.94           ! water vapour diffusivity from eq (13.3)
      !                                                                 ! of Pruppacher and Klett (1978)      
      ZALFAC(J)= 1.0 / (1.0 + &                                         ! Eq.11 in GoF95
                   ZRDSRV*ZDQSDT(J)*ZXLR(J)*ZDWAT(J)/(ZDTMP(J)*XCPD))   ! ZALFAC=wet-bulb factor (sans dim)     
      ZCPLW(J) = 4224.8482 + ZTAC(J) * &
                              ( -4.707 + ZTAC(J) * &
                                (0.08499 + ZTAC(J) * &
                                  (1.2826e-3 + ZTAC(J) * &
                                    (4.7884e-5 - 2.0027e-6* ZTAC(J))))) ! specific heat  
      !       
      ZRF(J)   = PRAIN(J) * ZCPLW(J) * ZALFAC(J) * &                    !Eq.12 in GoF95 !SIGNE?
                   (PSST(J) - PTA(J) + (PQSAT(J)-PQA(J))*ZXLR(J)/XCPD )
      !
      ! Momentum flux due to rainfall  
      !
      ZTAUR(J)=-0.85*(PRAIN(J) *ZVMOD(J)) !pp3752 in FBR96
      !
    ENDIF
    !
    !             4.4   Webb correction to latent heat flux
    ! 
    ZWBAR(J)=- (1./ZRDSRV)*ZUSR(J)*ZQSR(J) / (1.0+(1./ZRDSRV)*PQA(J)) &
               - ZUSR(J)*ZTSR(J)/PTA(J)                        ! Eq.21*rhoa in FBR96    
    !
    !             4.5   friction velocity which contains correction du to rain            
    !
    ZUSTAR2(J)= - (ZTAU(J) + ZTAUR(J)) / PRHOA(J)
    PUSTAR(J) =  SQRT(ZUSTAR2(J))
    !
    !             4.6   Total surface fluxes
    !           
    PSFTH (J) =  ZHF(J) + ZRF(J)
    PSFTQ (J) =  ZEF(J) / ZLV(J)
    ! 
  ENDIF
ENDDO                      
!-------------------------------------------------------------------------------
!
!       5.  FINAL STEP : TOTAL SURFACE FLUXES AND DERIVED DIAGNOSTICS 
!           -----------
!       5.1    Richardson number
!             
!
ZDIRCOSZW(:) = 1.
CALL SURFACE_RI(PSST,PQSAT,PEXNS,PEXNA,PTA,ZQASAT,  &
                  PZREF,PUREF,ZDIRCOSZW,ZVMOD,PRI)  
!
PRI(:) = MIN(PRI(:),XRIMAX)
!
!       5.2     Aerodynamical conductance and resistance
!             
ZAC(:) = PCH(:)*ZVMOD(:)
PRESA(:) = 1. / ZAC(:)
!
!PZ0HSEA(:) = PZ0SEA(:)
PZ0HSEA(:) = ZCHARN(:)
!IF (LHOOK) CALL DR_HOOK('COARE30_FLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WASP_FLUX
!
FUNCTION PSIFCTU(PZL)
!#######################################################################################
!
!****  *PSIFUNCTU*
!
!       PURPOSE
!       -------
!       To evaluate the stability function psi for wind speed (if KID=1) or
!       for temperature or humidity profiles (if KID.ne.1) from stability parameter
!       z/L. 
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Lik79 : Liu, W. T., K. B. Katsaros, and J. A. Businger, 1979: 
!       Bulk parameterization of air-sea exchanges of heat and water vapor including 
!       the molecular constraints at the interface. J. Atm. Sci., 36, 1722--1735.
!       DyH70 : Dyer, A. J., and B. B. Hicks, 1970: Flux-gradient relationship 
!       in the constant flux layer. Quart. J. Roy. Meteor. Soc., 96, 715--721.
!
!       AUTHOR
!       ------
!
!       MODIFICATIONS
!       -------------
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!-------------------------------------------------------------------------------
IMPLICIT NONE
!
!        0.  Declaration
!
!        0.1 declaration of arguments
!
REAL, DIMENSION(:), INTENT(IN)    :: PZL  !Obukhovs stability parameter
REAL, DIMENSION(SIZE(PZL))        :: PSIFCTU !function psi value
!        0.2 declaration of local variables
REAL, DIMENSION(SIZE(PZL)) :: ZX, ZY,ZC,ZPSIC,ZPSIK,ZF
INTEGER :: JJ
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
DO JJ=1,SIZE(PZL)
  IF(PZL(JJ)<0.) THEN
    ZX(JJ)   = (1.0 - 16. * PZL(JJ))**0.25         ! Kansas unstable
    ZPSIK(JJ)= 2.0 * LOG((1.0+ZX(JJ)       )/2.0) &
             +       LOG((1.0+ZX(JJ)*ZX(JJ))/2.0) &
             - 2.0 * atan(ZX(JJ)) &
             + 2.0 * atan(1.0)  
    !
    ZY(JJ)   = (1.0 - 10.15 * PZL(JJ))**0.3333     ! Convective
    ZPSIC(JJ)= 1.5 * LOG((ZY(JJ)*ZY(JJ)+ZY(JJ)+1.)/3.) &
             - (3.0**0.5) * atan((2.0*ZY(JJ)+1.0)/(3.0**0.5)) &
                  + 4.0        * atan(1.0)/(3.0**0.5)
      !
      ZF(JJ)   =PZL(JJ) * PZL(JJ) / (1.0+PZL(JJ)*PZL(JJ))
      !
      PSIFCTU(JJ)=(1.-ZF(JJ)) * ZPSIK(JJ) + ZF(JJ) * ZPSIC(JJ)
    ELSE
      ZC(JJ)=MIN(50.,0.35*PZL(JJ))           ! Stable
      PSIFCTU(JJ)=-((1.+1.*PZL(JJ))**1. + 0.6667*(PZL(JJ)-14.28)/EXP(ZC(JJ)) + 8.525)
  ENDIF
ENDDO
!IF (LHOOK) CALL DR_HOOK('MODE_COARE30_PSI:PSIFUNCTU',1,ZHOOK_HANDLE)

END FUNCTION PSIFCTU
!---------------------------------------------------------------------------------------
!
!#######################################################################################
FUNCTION PSIFCTT(PZL)
!#######################################################################################
!
!****  *PSIFUNCTU*
!
!       PURPOSE
!       -------
!       To evaluate the stability function psi for wind speed (if KID=1) or
!       for temperature or humidity profiles (if KID.ne.1) from stability parameter
!       z/L. 
!
!       EXTERNAL
!       --------
!
!       IMPLICIT ARGUMENTS
!       ------------------
!
!       REFERENCE
!       ---------
!       Lik79 : Liu, W. T., K. B. Katsaros, and J. A. Businger, 1979: 
!       Bulk parameterization of air-sea exchanges of heat and water vapor including 
!       the molecular constraints at the interface. J. Atm. Sci., 36, 1722--1735.
!       DyH70 : Dyer, A. J., and B. B. Hicks, 1970: Flux-gradient relationship 
!       in the constant flux layer. Quart. J. Roy. Meteor. Soc., 96, 715--721.
!
!       AUTHOR
!       ------
!
!       MODIFICATIONS
!       -------------
!-------------------------------------------------------------------------------
IMPLICIT NONE
!
!        0.  Declaration
!
!        0.1 declaration of arguments
!
REAL, DIMENSION(:), INTENT(IN)    :: PZL  !Obukhovs stability parameter
REAL, DIMENSION(SIZE(PZL))        :: PSIFCTT !function psi value
!        0.2 declaration of local variables
REAL, DIMENSION(SIZE(PZL)) :: ZX,ZY,ZC,ZPSIC,ZPSIK,ZF
INTEGER :: JJ
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
DO JJ=1,SIZE(PZL)
  IF(PZL(JJ)<0.) THEN
    ZX(JJ)   = (1. - 15. * PZL(JJ))**.5         ! Kansas unstable
    ZPSIK(JJ)= 2.0 * LOG((1.0+ZX(JJ)       )/2.0)
    !
    ZY(JJ)   = (1.0 - 34.15 * PZL(JJ))**0.3333  ! Convective
    ZPSIC(JJ)= 1.5 * LOG((ZY(JJ)*ZY(JJ)+ZY(JJ)+1.0)/3.) &
            - (3.0**0.5) * atan((2.0*ZY(JJ)+1.0)/(3.0**0.5)) &
            + 4.0        * atan(1.0)/(3.0**0.5)
    !
    ZF(JJ)   = PZL(JJ) * PZL(JJ) / (1.0+PZL(JJ)*PZL(JJ))
    !
    PSIFCTT(JJ)= (1.-ZF(JJ)) * ZPSIK(JJ) + ZF(JJ) * ZPSIC(JJ)
    ELSE
     ZC(JJ)=MIN(50.,0.35*PZL(JJ))           ! Stable
     PSIFCTT(JJ)=-((1.+2.*PZL(JJ)/3.)**1.5 + 0.6667*(PZL(JJ)-14.28)/EXP(ZC(JJ)) + 8.525)
  ENDIF
ENDDO
!IF (LHOOK) CALL DR_HOOK('MODE_COARE30_PSI:PSIFUNCTT',1,ZHOOK_HANDLE)

END FUNCTION PSIFCTT
!****************************************************************************
!****************************************************************************
 SUBROUTINE SURFACE_RI(PTG, PQS, PEXNS, PEXNA, PTA, PQA, &
 PZREF, PUREF, PDIRCOSZW, PVMOD, PRI )
 ! ######################################################################
 !
 !!**** *SURFACE_RI*
 !!
 !! PURPOSE
 !! -------
 !
 ! Computes the richardson number near the ground
 !
 !
 !!** METHOD
 !! ------
 !
 !
 !
 !
 !! EXTERNAL
 !! --------
 !!
 !!
 !! IMPLICIT ARGUMENTS
 !! ------------------
 !!
 !! MODD_CST
 !! MODD_GROUND_PAR
 !!
 !!
 !! REFERENCE
 !! ---------
 !!
 !!
 !! AUTHOR
 !! ------
 !!
 !! V. Masson * Meteo-France *
 !!
 !! MODIFICATIONS
 !! -------------
 !! Original 22/09/98
 !-------------------------------------------------------------------------------
 !
 !* 0. DECLARATIONS
 ! ------------
 !
 !
 IMPLICIT NONE
 !
 !* 0.1 declarations of arguments
 !
 !
 REAL, DIMENSION(:), INTENT(IN) :: PTG ! surface temperature
 REAL, DIMENSION(:), INTENT(IN) :: PQS ! surface specific humidity
 REAL, DIMENSION(:), INTENT(IN) :: PEXNS ! surface exner function
 REAL, DIMENSION(:), INTENT(IN) :: PTA ! temperature at the lowest level
 REAL, DIMENSION(:), INTENT(IN) :: PQA ! specific humidity
 ! at the lowest level
 REAL, DIMENSION(:), INTENT(IN) :: PEXNA ! exner function
 ! at the lowest level
 REAL, DIMENSION(:), INTENT(IN) :: PVMOD ! module of the horizontal wind
 !
 REAL, DIMENSION(:), INTENT(IN) :: PZREF ! reference height of the first
 ! atmospheric level
 REAL, DIMENSION(:), INTENT(IN) :: PUREF ! reference height of the wind
 ! ! NOTE this is different from ZZREF
 ! ! ONLY in stand-alone/forced mode,
 ! ! NOT when coupled to a model (MesoNH)
 REAL, DIMENSION(:), INTENT(IN) :: PDIRCOSZW! Cosine of the angle between
 ! ! the normal to the surface and
 ! ! the vertical
 !
 REAL, DIMENSION(:), INTENT(OUT) :: PRI ! Richardson number
 !
 !* 0.2 declarations of local variables
 !
 !
 REAL, DIMENSION(SIZE(PTG)) :: ZTHVA, ZTHVS
 REAL, DIMENSION(SIZE(PVMOD)) :: ZVMOD
 !*******
 REAL :: XRV, XRD, XG
 !USE MODI_WIND_THRESHOLD
 !-------------------------------------------------------------------------------
 !
 XRD=6.0221367E+23*1.380658E-23/28.9644E-3
 XRV=6.0221367E+23*1.380658E-23/18.0153E-3
 XG = 9.80665
 ! 1. Richardson number
 ! -----------------
 !
 ! virtual potential
 ! temperature at the
 ! first atmospheric level and
 ! at the surface
 !
 ZTHVA(:)=PTA(:)/PEXNA(:)*( 1.+(XRV/XRD-1.)*PQA(:) )
 ZTHVS(:)=PTG(:)/PEXNS(:)*( 1.+(XRV/XRD-1.)*PQS(:) )
 !
 ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
 !
 ! Richardson's number
 PRI(:) = XG * PDIRCOSZW(:) * PUREF(:) * PUREF(:) &
 * (ZTHVA(:)-ZTHVS(:)) / (0.5 * (ZTHVA(:)+ZTHVS(:)) ) &
 / (ZVMOD(:)*ZVMOD(:)) /PZREF(:)
 !-------------------------------------------------------------------------------
 !
 END SUBROUTINE SURFACE_RI
!
!
 FUNCTION CHARNOCK_WA(PWIND, PWAGE) RESULT(PCHARN)
 ! ############################################################################
 !
 !
 !!**** *CHARNOCK_WA
 !!
 !! PURPOSE
 !! -------
 !
 ! Compute the Charnock parameter using the wave age and the surface wind
 ! The formulation used for that depends on the surface wind range
 !
 !! AUTHOR
 !! ------
 !!
 !! MODIFICATIONS
 !! -------------
 !* 0. DECLARATIONS
 ! ------------
 !
 IMPLICIT NONE
 !
 !* 0.1 declarations of arguments
 !
 !
 REAL, INTENT(IN) :: PWIND ! wind
 REAL, INTENT(IN) :: PWAGE ! wave age
 !
 REAL :: PCHARN ! Charnock parameter
 REAL    :: ZLIMCHAR,ZLIMCHAR2,ZLIMCHAR1,ZAA, ZBB
 REAL, DIMENSION(1:2)         :: ZCOEFU,ZCOEFA2,ZCOEFB2
 REAL, DIMENSION(1:4)         :: ZPOLYU,ZCOEFA,ZCOEFB
 !
 !* 0.2 declarations of local variables
 !-------------------------------------------------------------------------------
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
 !
 PCHARN = ZCOEFU(1)*(PWIND**ZCOEFU(2))
 IF (PWIND >= 7.0) THEN
   ZAA = ZCOEFA(1) + ZCOEFA(2)*PWIND    &
                   + ZCOEFA(3)*PWIND**2 &
                   + ZCOEFA(4)*PWIND**3 
   ZBB = ZCOEFB(1) + ZCOEFB(2)*PWIND    &
                   + ZCOEFB(3)*PWIND**2 &
                   + ZCOEFB(4)*PWIND**3
   PCHARN = ZAA * (PWAGE**ZBB) 
!!  write(*,*) PCHARN,PWAGE
 ENDIF
 IF (PWIND >= 23.0) THEN
   ZAA = ZCOEFA2(1) + ZCOEFA2(2)*PWIND
   ZBB = ZCOEFB2(1) + ZCOEFB2(2)*PWIND
   PCHARN= ZAA * (PWAGE**ZBB)
   IF (PCHARN < ZLIMCHAR) THEN
     PCHARN = ZLIMCHAR
   ENDIF
 ENDIF
 IF (PWIND >= 25.0) THEN
   PCHARN = ZPOLYU(1) + ZPOLYU(2)*PWIND    &
                      + ZPOLYU(3)*PWIND**2 &
                      + ZPOLYU(4)*PWIND**3
   IF (PCHARN < ZLIMCHAR2) THEN
     PCHARN = ZLIMCHAR2
   ENDIF
 ENDIF
 IF (PCHARN > ZLIMCHAR1) THEN
         PCHARN = ZLIMCHAR1
 ENDIF
 END FUNCTION CHARNOCK_WA
 !****************************************************
 !
 !
 !****************************************************
 FUNCTION WIND_THRESHOLD(PWIND,PUREF) RESULT(PWIND_NEW)
 ! ############################################################################
 !
 !!**** *WIND_THRESHOLD*
 !!
 !! PURPOSE
 !! -------
 !
 ! Set a minimum value to the wind for exchange coefficient computations.
 ! This minimum value depends on the forcing height
 !
 !! AUTHOR
 !! ------
 !!
 !! V. Masson * Meteo-France *
 !!
 !! MODIFICATIONS
 !! -------------
 !! Original 09/2007
 !-------------------------------------------------------------------------------
 !
 !USE MODD_SURF_ATM, ONLY: XCISMIN, XVMODMIN, LALDTHRES
 !
 !* 0. DECLARATIONS
 ! ------------
 !
 IMPLICIT NONE
 !
 !* 0.1 declarations of arguments
 !
 !
 REAL, DIMENSION(:), INTENT(IN) :: PWIND ! wind
 REAL, DIMENSION(:), INTENT(IN) :: PUREF ! forcing level
 !
 REAL, DIMENSION(SIZE(PWIND)) :: PWIND_NEW ! modified wind
 !
 !
 !
 !* 0.2 declarations of local variables
 LOGICAL :: LALDTHRES=.FALSE.
 REAL :: XCISMIN=6.7E-5
 REAL :: XVMODMIN=0.
 !
 !-------------------------------------------------------------------------------
 !
 ! wind gradient
 !
 IF (.NOT.LALDTHRES) THEN
 !
 ! minimum value for exchange coefficients computations : 1m/s / 10m
   PWIND_NEW = MAX(PWIND , 0.1 * MIN(10.,PUREF) )
 ELSE
 ! minimum value for exchange coefficients computations : 1m/s / 10m
   PWIND_NEW = MAX( XVMODMIN, SQRT( PWIND**2 + (XCISMIN*PUREF)**2 ) )
 ENDIF
 !
 !-------------------------------------------------------------------------------
 !
 END FUNCTION WIND_THRESHOLD
!****************************************************
END PROGRAM 
