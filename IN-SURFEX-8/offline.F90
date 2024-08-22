!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
! -------------------------------------------------
PROGRAM OFFLINE
!
! -------------------------------------------------
! Driver structure
! ----------------
! 1. Initializations
! 2. Temporal loops
!   2.a Read forcing
!   2.b Interpolate forcing in time
!   2.c Run surface
!   2.d Write prognostics and diagnostics variables
!
! modifications 
! 09/2012 G. Pigeon: coherence between radiation and zenith angle because of
!                    trouble with radiation received by wall in TEB
! 03/2014 E. Martin change indices names in OMP module according to GMAP changes
! 05/2014 B. Decharme delete trip
!                     I5 format to print DAY
! 04/2013 P. Lemoigne Add XDELTA_OROG to fix the maximum difference allowed between
!                     forcing and surface file orographies if LSET_FORC_ZS=.F
! 12/2013 S.Senesi    Add call to Gelato diag files init and close
!
! 02/2016: replace DOUBLE PRECISION by REAL to handle problem for promotion of real with GMKPACK or IBM SP
!
! -------------------------------------------------
!
USE MODD_OFF_SURFEX_n
!
USE MODD_FORC_ATM,  ONLY: CSV         ,&! name of all scalar variables
                            XDIR_ALB    ,&! direct albedo for each band
                            XSCA_ALB    ,&! diffuse albedo for each band
                            XEMIS       ,&! emissivity
                            XTSRAD      ,&! radiative temperature
                            XTSUN       ,&! solar time                    (s from midnight)
                            XZS         ,&! orography                             (m)
                            XZREF       ,&! height of T,q forcing                 (m)
                            XUREF       ,&! height of wind forcing                (m)
                            XTA         ,&! air temperature forcing               (K)
                            XQA         ,&! air humidity forcing                  (kg/m3)
                            XSV         ,&! scalar variables
                            XU          ,&! zonal wind                            (m/s)
                            XV          ,&! meridian wind                         (m/s)
                            XDIR_SW     ,&! direct  solar radiation (on horizontal surf.)
                            XSCA_SW     ,&! diffuse solar radiation (on horizontal surf.)
                            XSW_BANDS   ,&! mean wavelength of each shortwave band (m)
                            XZENITH     ,&! zenithal angle       (radian from the vertical)
                            XZENITH2    ,&! zenithal angle       (radian from the vertical)
                            XAZIM       ,&! azimuthal angle      (radian from North, clockwise)
                            XLW         ,&! longwave radiation (on horizontal surf.)
                            XPS         ,&! pressure at atmospheric model surface (Pa)
                            XPA         ,&! pressure at forcing level             (Pa)
                            XRHOA       ,&! density at forcing level              (kg/m3)
                            XCO2        ,&! CO2 concentration in the air          (kg/m3)
                            XSNOW       ,&! snow precipitation                    (kg/m2/s)
                            XRAIN       ,&! liquid precipitation                  (kg/m2/s)
                            XSFTH       ,&! flux of heat                          (W/m2)
                            XSFTQ       ,&! flux of water vapor                   (kg/m2/s)
                            XSFU        ,&! zonal momentum flux                   (m/s)
                            XSFV        ,&! meridian momentum flux                (m/s)
                            XSFCO2      ,&! flux of CO2                           (kg/m2/s)
                            XSFTS       ,&! flux of scalar var.                   (kg/m2/s)
                            XPEW_A_COEF ,&! implicit coefficients
                            XPEW_B_COEF ,&! needed if HCOUPLING='I'
                            XPET_A_COEF ,&
                            XPEQ_A_COEF ,&
                            XPET_B_COEF ,&
                            XPEQ_B_COEF ,&
                            XTSURF      ,&! effective temperature                  (K)
                            XZ0         ,&! surface roughness length for momentum  (m)
                            XZ0H        ,&! surface roughness length for heat      (m)
                            XQSURF        ! specific humidity at surface           (kg/kg)
!
USE MODD_SURF_CONF,  ONLY : CPROGNAME, CSOFTWARE
USE MODD_CSTS,       ONLY : XPI, XDAY, XRV, XRD, XG
USE MODD_IO_SURF_ASC,ONLY : CFILEIN,CFILEIN_SAVE,CFILEOUT,CFILEPGD
USE MODD_SURF_PAR
USE MODD_IO_SURF_FA, ONLY : CFILEIN_FA, CFILEIN_FA_SAVE,       &
                            CFILEOUT_FA, NUNIT_FA, CDNOMC,     &
                            IVERBFA, LFANOCOMPACT, CFILEPGD_FA  
USE MODD_IO_SURF_LFI,ONLY : CFILEIN_LFI, CFILEIN_LFI_SAVE, CLUOUT_LFI, CFILEOUT_LFI, &
                            LMNH_COMPATIBLE, CFILEPGD_LFI  
USE MODD_IO_SURF_NC, ONLY : CFILEIN_NC, CFILEIN_NC_SAVE, CFILEOUT_NC, CLUOUT_NC, &
                            CFILEPGD_NC, LDEF
USE MODD_IO_SURF_OL, ONLY : XSTART, XCOUNT, XSTRIDE,           &
                              LDEFINED_NATURE, LDEFINED_SEA,    &
                              LDEFINED_WATER,  LDEFINED_TOWN,   &
                              LDEFINED_SURF_ATM, LPARTW,        &
                              XSTARTW, XCOUNTW, LTIME_WRITTEN,  &
                              NSTEP_OUTPUT  
USE MODD_WRITE_BIN,  ONLY : NWRITE
!
USE MODD_SURFEX_MPI, ONLY : NCOMM, NPROC, NRANK, NPIO, WLOG_MPI, PREP_LOG_MPI,   &
                            NINDEX, NSIZE_TASK, XTIME_NPIO_READ, XTIME_NPIO_WRITE, &
                            XTIME_COMM_READ, XTIME_COMM_WRITE, XTIME_SEA, &
                            XTIME_NATURE, XTIME_WATER, XTIME_TOWN, &
                            XTIME_INIT_SEA, XTIME_INIT_WATER, &
                            XTIME_INIT_NATURE, XTIME_INIT_TOWN, &
                            XTIME_WRITE, XTIME_CALC, XTIME_OMP_BARR, &
                            XTIME_CALC_WRITE, IDX_W, END_LOG_MPI
!
USE MODD_SURFEX_OMP, ONLY :  NINDX1SFX, NINDX2SFX, NBLOCK, NBLOCKTOT, &
                             INIT_DIM, RESET_DIM, NWORK, NWORK2, XWORK, XWORK2, XWORK3, &
                             NWORK_FULL, NWORK2_FULL, XWORK_FULL, XWORK2_FULL
!
USE MODD_COUPLING_TOPD, ONLY : NNB_TOPD, NNB_STP_RESTART, LBUDGET_TOPD, LTOPD_STEP, &
                               LCOUPL_TOPD, NTOPD_STEP, NYEAR, NMONTH, NDAY, NH, NM
USE MODD_TOPODYN, ONLY : XTOPD_STEP, NNB_TOPD_STEP, XQTOT, XQB_RUN, XQB_DR
!
USE MODD_SLOPE_EFFECT, ONLY: XZS_THREAD,XZS_XY_THREAD,XSLOPANG_THREAD,&
                             XSLOPAZI_THREAD,XSURF_TRIANGLE_THREAD
!
USE MODD_SFX_OASIS, ONLY : LOASIS, XRUNTIME
!
USE MODE_POS_SURF
!
USE MODE_CRODEBUG
!
USE MODN_IO_OFFLINE
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_TEST_NAM_VAR_SURF
USE MODI_CLOSE_NAMELIST
USE MODI_READ_ALL_NAMELISTS
USE MODI_OPEN_CLOSE_BIN_ASC_FORC
USE MODI_OPEN_FILEIN_OL
USE MODI_OL_READ_ATM_CONF
USE MODI_ABOR1_SFX
USE MODI_OL_ALLOC_ATM
USE MODI_COMPARE_OROGRAPHY
USE MODI_SUNPOS
USE MODI_INIT_INDEX_MPI
USE MODI_OL_READ_ATM
USE MODI_GET_SIZES_PARALLEL
USE MODI_IO_BUFF_CLEAN
USE MODI_INIT_SURF_ATM_n
USE MODI_INIT_SURF_LANDUSE_n
USE MODI_OL_TIME_INTERP_ATM
USE MODI_COUPLING_SURF_ATM_n
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_WRITE_SURF_ATM_n
USE MODI_WRITE_HEADER_MNH
USE MODI_FLAG_UPDATE
USE MODI_FLAG_DIAG_UPDATE
USE MODI_DIAG_SURF_ATM_n
USE MODI_WRITE_DIAG_SURF_ATM_n
USE MODI_GET_SURF_VAR_n
USE MODI_GATHER_AND_WRITE_MPI
USE MODI_CLOSE_FILEIN_OL
USE MODI_CLOSE_FILEOUT_OL
USE MODI_INIT_OUTPUT_OL_n
USE MODI_INIT_OUTPUT_NC_n
!
USE MODI_WRITE_HEADER_FA
USE MODI_ABOR1_SFX
!
USE MODI_WRITE_DISCHARGE_FILE
USE MODI_WRITE_BUDGET_COUPL_ROUT
USE MODI_PREP_RESTART_COUPL_TOPD
!
USE MODI_INIT_SLOPE_PARAM
USE MODI_SLOPE_RADIATIVE_EFFECT
!
USE MODI_SFX_OASIS_INIT
USE MODI_SFX_OASIS_READ_NAM
USE MODI_SFX_OASIS_DEF_OL
USE MODI_SFX_OASIS_RECV_OL
USE MODI_SFX_OASIS_SEND_OL
USE MODI_SFX_OASIS_END
!RJ: missing modi
USE MODI_LOCAL_SLOPE_PARAM
!
USE MODE_GLT_DIA_LU
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef AIX64
!$ USE OMP_LIB
#endif
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE 'mpif.h'
#endif
!
#ifndef AIX64
!$ INCLUDE 'omp_lib.h'
#endif
!
!*      0.    declarations of local variables
!
 CHARACTER(LEN=3), PARAMETER       :: YINIT     = 'ALL'
!
 CHARACTER(LEN=28)                 :: YLUOUT    = 'LISTING_OFFLINE             '
!
INTEGER                           :: IYEAR               ! current year (UTC)
INTEGER                           :: IMONTH              ! current month (UTC)
INTEGER                           :: IDAY                ! current day (UTC)
REAL                              :: ZTIME               ! current time since start of the run (s)
REAL                              :: ZTIMEC              ! current duration since start of the run (s)
!
INTEGER                           :: IYEAR_OUT           ! output year name
INTEGER                           :: IMONTH_OUT          ! output month name
INTEGER                           :: IDAY_OUT            ! output day name
REAL                              :: ZTIME_OUT           ! output time since start of the run (s)
!
INTEGER, DIMENSION(11)  :: IDATEF
!
 CHARACTER(LEN=28), PARAMETER      :: YATMFILE     = '                            '
 CHARACTER(LEN=6),  PARAMETER      :: YATMFILETYPE = '      '
 CHARACTER(LEN=2),  PARAMETER      :: YTEST        = 'OK'          ! must be equal to 'OK'
!
REAL, DIMENSION(:), POINTER       :: ZLAT                ! latitude                         (rad)
REAL, DIMENSION(:), POINTER       :: ZLON                ! longitude                        (rad)
REAL, DIMENSION(:), POINTER       :: ZZS_FORC            ! orography                        (m)  
REAL, DIMENSION(:), POINTER       :: ZZREF               ! Forcing level for T
REAL, DIMENSION(:), POINTER       :: ZUREF               ! Forcing level for U
!
REAL                              :: ZTSTEP              ! atmospheric time-step            (s)
!
INTEGER                           :: INI                 ! grid dimension
INTEGER                           :: JLOOP               ! loop counter
INTEGER                           :: IBANDS              ! Number of radiative bands 
INTEGER                           :: INB_STEP_ATM        ! Number of atmospheric time-steps
INTEGER                           :: INB_ATM             ! Number of Isba time-steps 
                                                         ! within a forcing time-step
INTEGER                           :: ID_FORC             ! indice of forcing in the file
INTEGER                           :: INB_LINES           ! nb of lines to read in the forcing file
INTEGER                           :: IDMAX               ! nb of lines to read in the forcing file at last 
INTEGER                           :: JFORC_STEP          ! atmospheric loop index
INTEGER                           :: JSURF_STEP          ! isba loop index
INTEGER                           :: ICOUNT              ! day counter 
REAL                              :: ZDURATION           ! duration of run                     (s)
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTA                 ! air temperature forcing               (K)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZQA                 ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWIND               ! wind speed                            (m/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSCA_SW             ! diffuse solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIR_SW             ! direct  solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLW                 ! longwave radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSNOW               ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRAIN               ! liquid precipitation                  (kg/m2/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPS                 ! pressure at forcing level             (Pa)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCO2                ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIR                ! wind direction
INTEGER                           :: ILUOUT              ! ascii output unit number
INTEGER                           :: ILUNAM              ! namelist unit number
INTEGER                           :: IRET                ! error return code
INTEGER                           :: INB 
INTEGER                           :: INW, JNW
 CHARACTER(LEN=14)                :: YTAG      
LOGICAL                           :: GFOUND              ! return logical when reading namelist
LOGICAL                           :: GSHADOWS    
REAL, DIMENSION(:),   ALLOCATABLE :: ZSW                 ! total solar radiation (on horizontal surf.)
REAL, DIMENSION(:),   ALLOCATABLE :: ZCOEF               ! coefficient for solar radiation interpolation near sunset/sunrise
!
! Flag diag :
!
INTEGER                           :: I2M, IBEQ, IDSTEQ
LOGICAL                           :: GFRAC, GDIAG_GRID, GSURF_BUDGET, GRAD_BUDGET, GCOEF,    &
                                     GSURF_VARS, GDIAG_OCEAN, GDIAG_SEAICE, GWATER_PROFILE,  &
                                     GINTERPOL_TS, GSURF_EVAP_BUDGET, GFLOOD, GPGD_ISBA,     &
                                     GCH_NO_FLUX_ISBA, GSURF_MISC_BUDGET_ISBA, GPGD_TEB,     &
                                     GSURF_MISC_BUDGET_TEB
!
! Inquiry mode arrays:
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZSEA, ZWATER, ZNATURE, ZTOWN
REAL, DIMENSION(:), ALLOCATABLE   :: ZSEA_FULL, ZWATER_FULL, ZNATURE_FULL, ZTOWN_FULL
REAL, DIMENSION(:), ALLOCATABLE   :: ZT2M, ZQ2M
REAL, DIMENSION(:), ALLOCATABLE   :: ZZ0, ZZ0H, ZQS
REAL, DIMENSION(:), ALLOCATABLE   :: ZQS_SEA, ZQS_WATER, ZQS_NATURE, ZQS_TOWN
REAL, DIMENSION(:), ALLOCATABLE   :: ZPSNG, ZPSNV
REAL, DIMENSION(:), ALLOCATABLE   :: ZZ0EFF
REAL, DIMENSION(:), ALLOCATABLE   :: ZZS
REAL, DIMENSION(:), ALLOCATABLE   :: ZZ0_FULL, ZZ0EFF_FULL, ZZS_FULL
INTEGER :: ISERIES, ISIZE
!
! MPI variables
!
 CHARACTER(LEN=100) :: YNAME
 CHARACTER(LEN=10)  :: YRANK
INTEGER :: ILEVEL, INFOMPI, J, INKPROMA, JBLOCK
INTEGER, DIMENSION(:), ALLOCATABLE :: ISIZE_OMP
REAL :: XTIME0, XTIME1, XTIME
!
! SFX - OASIS coupling variables
!
INTEGER  :: ILOCAL_COMM  ! Local communicator
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! --------------------------------------------------------------------------------------
!
!*     0.1.   MPI, OASIS, and dr_hook initializations
!
INFOMPI=1
!
#ifdef CPLOASIS
!Must be call before DRHOOK !
 CALL SFX_OASIS_INIT(CNAMELIST,ILOCAL_COMM)
#else
LOASIS   = .FALSE.
XRUNTIME = 0.0
#endif
!
#ifdef SFX_MPI
IF(.NOT.LOASIS)THEN
 CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,ILEVEL,INFOMPI)
 IF (INFOMPI /= MPI_SUCCESS) THEN 
    CALL ABOR1_SFX('OFFLINE: ERROR WHEN INITIALIZING MPI')
 ENDIF
ENDIF
#endif
!
IF (LHOOK) CALL DR_HOOK('OFFLINE',0,ZHOOK_HANDLE)
!
 CSOFTWARE='OFFLINE'
!
#ifdef SFX_MPI
IF(LOASIS)THEN
  NCOMM=ILOCAL_COMM
ELSE
  NCOMM=MPI_COMM_WORLD
ENDIF
 CALL MPI_COMM_SIZE(NCOMM,NPROC,INFOMPI)
 CALL MPI_COMM_RANK(NCOMM,NRANK,INFOMPI)
#endif
!
!RJ: init modd_surefx_omp
!$OMP PARALLEL
!$ NBLOCKTOT = OMP_GET_NUM_THREADS()
!$ NBLOCK = OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
!
 CALL PREP_LOG_MPI
!
 CALL WLOG_MPI(' ')
!
 CALL WLOG_MPI('NBLOCKTOT ',KLOG=NBLOCKTOT)
!
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!
!*      0.3.   Open ascii file for writing
!
WRITE(YRANK,FMT='(I10)') NRANK
YNAME=TRIM(YLUOUT)//ADJUSTL(YRANK)
!
 CLUOUT_LFI =  ADJUSTL(ADJUSTR(YNAME)//'.txt')
 CLUOUT_NC  =  ADJUSTL(ADJUSTR(YNAME)//'.txt')
!
 CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YNAME)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
!
IF ( NRANK==NPIO ) THEN
  !
!$OMP SINGLE
  !
!RJ: be verbose just for openmp
  IF(NBLOCKTOT==1) THEN
!$  WRITE(*,*) "CAUTION: DID YOU THINK TO SET OMP_NUM_THREADS=1?"
!$  WRITE(*,*) "PLEASE VERIFY OMP_NUM_THREADS IS INITIALIZED : TYPE ECHO $OMP_NUM_THREADS IN A TERMINAL" 
  !
!$  WRITE(ILUOUT,*) "CAUTION: DID YOU THINK TO SET OMP_NUM_THREADS=1?"
!$  WRITE(ILUOUT,*) "PLEASE VERIFY OMP_NUM_THREADS IS INITIALIZED : TYPE ECHO $OMP_NUM_THREADS IN A TERMINAL" 
  ENDIF
  !
!$OMP END SINGLE
  !
ENDIF
!
!*      0.4.   Reads namelists
!
 CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
!
 CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND,ILUOUT)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
 CALL CLOSE_NAMELIST('ASCII ',ILUNAM)
!
IF (NPROC==1) THEN 
  XIO_FRAC=1.
ELSE
  XIO_FRAC = MAX(MIN(XIO_FRAC,1.),0.)
ENDIF
!
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CSURF_FILETYPE',CSURF_FILETYPE,'ASCII ','LFI   ','FA    ','NC    ')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CTIMESERIES_FILETYPE',CTIMESERIES_FILETYPE,'NETCDF','TEXTE ','BINARY',&
                                                                            'ASCII ','LFI   ','FA    ',&
                                                                            'NONE  ','OFFLIN','NC    ')  
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CFORCING_FILETYPE',CFORCING_FILETYPE,'NETCDF','ASCII ','BINARY')
!
IF (NSCAL>59) CALL ABOR1_SFX("OFFLINE: NSCAL MUST BE LOWER THAN OR EQUAL TO 59")
!
!
IF (CTIMESERIES_FILETYPE=='NETCDF') CTIMESERIES_FILETYPE='OFFLIN'
!
IF (NRANK==NPIO) THEN
  !
  CFILEPGD = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
  CFILEIN  = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEIN_SAVE = CFILEIN
  !
  CFILEPGD_LFI = CPGDFILE
  CFILEIN_LFI  = CPREPFILE
  CFILEIN_LFI_SAVE = CFILEIN_LFI
  !
  CFILEPGD_FA = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
  CFILEIN_FA  = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEIN_FA_SAVE  = CFILEIN_FA
  !
  CFILEPGD_NC = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
  CFILEIN_NC  = ADJUSTL(ADJUSTR(CPREPFILE)//'.nc')
  CFILEIN_NC_SAVE  = CFILEIN_NC
  !
ENDIF
!
!     Allocations of Surfex Types
 CALL SURFEX_ALLOC_LIST(1)
 YSURF_CUR => YSURF_LIST(1)
!
!     Reading all namelist (also assimilation)
 CALL READ_ALL_NAMELISTS(YSURF_CUR, &
                         CSURF_FILETYPE,'ALL',.FALSE.)
!
!
!*      0.5.   Reads SFX - OASIS coupling namelists
!
 CALL SFX_OASIS_READ_NAM(CSURF_FILETYPE,XTSTEP_SURF)
!
!*      0.6   Assume FA filetype consistency 
!
 CPROGNAME = CSURF_FILETYPE
!
! --------------------------------------------------------------------------------------
!
!*      1.    Initializations
!
!       netcdf file handling
!
IF (NRANK==NPIO) THEN
  !
  XSTART            = NUNDEF
  XSTRIDE           = NUNDEF
  XCOUNT            = NUNDEF
  XSTARTW           = 0
  XCOUNTW           = 1
  LPARTW            = .TRUE.
  LDEFINED_SURF_ATM = .FALSE.
  LDEFINED_NATURE   = .FALSE.
  LDEFINED_TOWN     = .FALSE.
  LDEFINED_WATER    = .FALSE.
  LDEFINED_SEA      = .FALSE.
  !
ENDIF
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('READ NAMELISTS ',PLOG=XTIME)
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!       forcing file handling
!
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') CALL OPEN_CLOSE_BIN_ASC_FORC('CONF ',CFORCING_FILETYPE,'R')
IF (CFORCING_FILETYPE=='NETCDF') CALL OPEN_FILEIN_OL
!
!
!       splitting of the grid
!
GSHADOWS = LSHADOWS_SLOPE .OR. LSHADOWS_OTHER
 CALL INIT_INDEX_MPI(YSURF_CUR, &
                    CSURF_FILETYPE, YALG_MPI, XIO_FRAC, GSHADOWS)
!
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('TIME_NPIO_READ init_index ',PLOG=XTIME_NPIO_READ)
 CALL WLOG_MPI('TIME_COMM_READ init_index ',PLOG=XTIME_COMM_READ)
XTIME_NPIO_READ = 0.
XTIME_COMM_READ = 0.
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('INIT_INDEX_MPI ',PLOG=XTIME)
 CALL WLOG_MPI(' ')
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!       configuration of run
!
 CALL OL_READ_ATM_CONF(YSURF_CUR, &
                       CSURF_FILETYPE, CFORCING_FILETYPE,            &
                      ZDURATION, ZTSTEP, INI, IYEAR, IMONTH, IDAY,  &
                      ZTIME, ZLAT, ZLON, ZZS_FORC, ZZREF, ZUREF     )
!
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('TIME_NPIO_READ forc conf ',PLOG=XTIME_NPIO_READ)
 CALL WLOG_MPI('TIME_COMM_READ forc conf ',PLOG=XTIME_COMM_READ)
XTIME_NPIO_READ = 0.
XTIME_COMM_READ = 0.
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('OL_READ_ATM_CONF ',PLOG=XTIME)
 CALL WLOG_MPI(' ')
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!*     time steps coherence check 
!
IF ( (MOD(XTSTEP_OUTPUT,ZTSTEP)*MOD(ZTSTEP,XTSTEP_OUTPUT) /= 0) .OR. (MOD(ZTSTEP,XTSTEP_SURF) /= 0) ) THEN
   WRITE(ILUOUT,*)' FORCING  AND OUTPUT/SURFACE TIME STEP SHOULD BE MULTIPLE', &
     NINT(ZTSTEP),NINT(XTSTEP_OUTPUT),NINT(XTSTEP_SURF)    
   CALL ABOR1_SFX('OFFLINE: FORCING  AND OUTPUT/SURFACE TIME STEP SHOULD BE MULTIPLE')
ENDIF
!
IF ( ZTIME /= 0. .AND. MOD(ZTIME,XTSTEP_SURF) /= 0  ) THEN
   WRITE(ILUOUT,*)' INITIAL AND SURFACE TIME STEP SHOULD BE MULTIPLE', &
   NINT(ZTIME),NINT(XTSTEP_SURF)  
   CALL ABOR1_SFX('OFFLINE: INITIAL AND SURFACE TIME STEP SHOULD BE MULTIPLE')
ENDIF
!
IF(LOASIS.AND.ZDURATION/=XRUNTIME)THEN
   WRITE(ILUOUT,*)'Total simulated time given by Forcing field and OASIS namcouple are different'
   WRITE(ILUOUT,*)'From Forcing (s) : ',ZDURATION, 'From OASIS   (s) : ',XRUNTIME
   CALL ABOR1_SFX('OFFLINE: TOTAL SIMULATED TIME DIFFERENT BETWEEN FORCING AND OASIS')
ENDIF
!
INB_STEP_ATM  = INT(ZDURATION / ZTSTEP)
INB_ATM       = INT(ZTSTEP / XTSTEP_SURF)
NSTEP_OUTPUT  = INT(ZDURATION / XTSTEP_OUTPUT)
!
XTOPD_STEP = 0
NNB_TOPD_STEP = 0
NTOPD_STEP = 0
IF ( LCOUPL_TOPD ) THEN
  !
  XTOPD_STEP = FLOAT(NNB_TOPD)* XTSTEP_SURF
  NNB_TOPD_STEP = INT( ZDURATION / XTOPD_STEP )
  !
  IF ( NNB_STP_RESTART==0 .AND. .NOT.LRESTART ) NNB_STP_RESTART = -1
  !
  NTOPD_STEP = 1
  !
ENDIF
!
!       allocation of variables
!
IBANDS = 1
!
 CALL OL_ALLOC_ATM(INI,IBANDS,NSCAL)
!
XZS   = ZZS_FORC
XZREF = ZZREF
XUREF = ZUREF
!
!       compare orography
!
 CALL COMPARE_OROGRAPHY(YSURF_CUR, &
                         CSURF_FILETYPE, LSET_FORC_ZS, XDELTA_OROG)
!
!       miscellaneous initialization
!
ICOUNT = 0
ZTIMEC = 0.
!
ALLOCATE(ISIZE_OMP(0:NBLOCKTOT-1))
 CALL GET_SIZES_PARALLEL(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%UG, YSURF_CUR%U, &
                         NBLOCKTOT,INI,0,ISIZE_OMP, GSHADOWS)
 CALL SUNPOS(ISIZE_OMP, IYEAR, IMONTH, IDAY, ZTIME, ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
DEALLOCATE(ISIZE_OMP)
!
!number of lines read in forcing files
INB_LINES=1
IF (NB_READ_FORC.EQ.1) THEN
  INB_LINES=INB_STEP_ATM
ELSEIF (NB_READ_FORC.NE.0) THEN
  !to be sure the number of readings will be NB_READ_FORC as a maximum
  INB_LINES=CEILING(1.*(INB_STEP_ATM+1)/NB_READ_FORC)
ENDIF
!number of lines to be read effectively
IDMAX=INB_LINES+1
!effective number of readings of the forcing files
NB_READ_FORC=CEILING(1.*(INB_STEP_ATM+1)/INB_LINES)
!
!     open Gelato specific diagnostic files (if requested by 
!     Gelato wizzard user)
!
#if ! defined in_arpege
 CALL OPNDIA()
#endif
!
!       allocate local atmospheric variables
!
IF (.NOT.ALLOCATED(ZTA)) ALLOCATE(ZTA    (INI,INB_LINES+1)) 
IF (.NOT.ALLOCATED(ZQA))ALLOCATE(ZQA    (INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZWIND))ALLOCATE(ZWIND  (INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZDIR_SW))ALLOCATE(ZDIR_SW(INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZSCA_SW))ALLOCATE(ZSCA_SW(INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZLW))ALLOCATE(ZLW    (INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZSNOW))ALLOCATE(ZSNOW  (INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZRAIN))ALLOCATE(ZRAIN  (INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZPS))ALLOCATE(ZPS    (INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZCO2))ALLOCATE(ZCO2   (INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZDIR))ALLOCATE(ZDIR   (INI,INB_LINES+1))
IF (.NOT.ALLOCATED(ZCOEF))ALLOCATE(ZCOEF   (INI))
!
IF (.NOT.ALLOCATED(ZSW))ALLOCATE(ZSW    (INI))
!
!      computes initial air co2 concentration and  density
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('COMPARE_OROGRAPHY SUNPOS ',PLOG=XTIME)
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!* opens forcing files (if ASCII or BINARY)
!
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') CALL OPEN_CLOSE_BIN_ASC_FORC('OPEN ',CFORCING_FILETYPE,'R')
!
 CALL OL_READ_ATM(&
                  CSURF_FILETYPE, CFORCING_FILETYPE, 1,             &
                   ZTA,ZQA,ZWIND,ZDIR_SW,ZSCA_SW,ZLW,ZSNOW,ZRAIN,ZPS,&
                   ZCO2,ZDIR,LLIMIT_QAIR                           ) 
!
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('TIME_NPIO_READ forc ',PLOG=XTIME_NPIO_READ)
 CALL WLOG_MPI('TIME_COMM_READ forc ',PLOG=XTIME_COMM_READ)
XTIME_NPIO_READ = 0.
XTIME_COMM_READ = 0.
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('OL_READ_ATM0 ',PLOG=XTIME)
 CALL WLOG_MPI(' ')
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
XCO2(:)  = ZCO2(:,1)
XRHOA (:) = ZPS(:,1) / (XRD * ZTA(:,1) * ( 1.+((XRV/XRD)-1.)*ZQA(:,1) ) + XG * XZREF )
!                 
!       surface Initialisation     
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('CO2 RHOA ',PLOG=XTIME)
!
 CALL IO_BUFF_CLEAN
!
 CALL SURFEX_DEALLO_LIST
 CALL SURFEX_ALLOC_LIST(NBLOCKTOT)
!
 CALL GOTO_MODEL(1)
ALLOCATE(ISIZE_OMP(0:NBLOCKTOT-1))
 CALL GET_SIZES_PARALLEL(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%UG, YSURF_CUR%U, &
                         NBLOCKTOT,INI,0,ISIZE_OMP, GSHADOWS)
DO J=0,NBLOCKTOT-1
 CALL WLOG_MPI("SIZES_OMP ",KLOG=J,KLOG2=ISIZE_OMP(J))
ENDDO
!
!$OMP PARALLEL PRIVATE(INKPROMA,XTIME,XTIME0)
!
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
!$ NBLOCK = OMP_GET_THREAD_NUM()
!
IF (NBLOCK==NBLOCKTOT) THEN
  CALL INIT_DIM(ISIZE_OMP,0,INKPROMA,NINDX1SFX,NINDX2SFX)
ELSE
  CALL INIT_DIM(ISIZE_OMP,NBLOCK,INKPROMA,NINDX1SFX,NINDX2SFX)
ENDIF
!
IF (NBLOCK==0) THEN
  CALL GOTO_MODEL(NBLOCKTOT)
ELSE
  CALL GOTO_MODEL(NBLOCK)
ENDIF
!
 CALL INIT_SURF_ATM_n(YSURF_CUR, &
                            CSURF_FILETYPE, YINIT, LLAND_USE,                      &
                     INKPROMA, NSCAL, IBANDS,                               &
                     CSV,XCO2(NINDX1SFX:NINDX2SFX),XRHOA(NINDX1SFX:NINDX2SFX),          &
                     XZENITH(NINDX1SFX:NINDX2SFX),XAZIM(NINDX1SFX:NINDX2SFX),XSW_BANDS, &
                     XDIR_ALB(NINDX1SFX:NINDX2SFX,:), XSCA_ALB(NINDX1SFX:NINDX2SFX,:),  &
                     XEMIS(NINDX1SFX:NINDX2SFX), XTSRAD(NINDX1SFX:NINDX2SFX),           &
                     XTSURF(NINDX1SFX:NINDX2SFX),                           &
                     IYEAR, IMONTH, IDAY, ZTIME,                            &
                     YATMFILE, YATMFILETYPE, YTEST                          )
!
! initialization routines to compute shadows
IF (GSHADOWS) THEN
  IF (NBLOCK==0) THEN
    CALL INIT_SLOPE_PARAM(YSURF_CUR%UG, &
                          ZZS_FORC,INI,ZLAT)
  END IF 
  !$OMP BARRIER
  CALL LOCAL_SLOPE_PARAM(NINDX1SFX,NINDX2SFX)
END IF
!
 CALL RESET_DIM(INI,INKPROMA,NINDX1SFX,NINDX2SFX)
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('INIT_SURF_ATM ',PLOG=XTIME)
 CALL WLOG_MPI(' ')
!
 CALL WLOG_MPI('TIME_NPIO_READ init ',PLOG=XTIME_NPIO_READ)
 CALL WLOG_MPI('TIME_COMM_READ init ',PLOG=XTIME_COMM_READ)
 CALL WLOG_MPI(' ')
!
!$OMP END PARALLEL
!
XTIME_NPIO_READ = 0.
XTIME_COMM_READ = 0.
!
!   Land use or/and vegetation dynamic
!                  
 CALL INIT_SURF_LANDUSE_n(YSURF_CUR, &
                          CSURF_FILETYPE,YINIT,LLAND_USE,             &
                       INI, NSCAL, IBANDS,                           &
                       CSV,XCO2(NINDX1SFX:NINDX2SFX),XRHOA(NINDX1SFX:NINDX2SFX), &
                       XZENITH(NINDX1SFX:NINDX2SFX),XAZIM(NINDX1SFX:NINDX2SFX),  &
                       XSW_BANDS,XDIR_ALB(NINDX1SFX:NINDX2SFX,:),          &
                       XSCA_ALB(NINDX1SFX:NINDX2SFX,:),                    &
                       XEMIS(NINDX1SFX:NINDX2SFX),XTSRAD(NINDX1SFX:NINDX2SFX),   &
                       XTSURF(NINDX1SFX:NINDX2SFX),                           &
                       IYEAR, IMONTH, IDAY, ZTIME,                   &
                       YATMFILE, YATMFILETYPE, YTEST                 )
!
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!
 CALL INIT_CRODEBUG(YSURF_CUR%IM%I)
!
! * SURFEX - OASIS  grid, partitions and local field definitions
!
IF(LOASIS)THEN
  CALL SFX_OASIS_DEF_OL(YSURF_CUR%IM%I, &
                        YSURF_CUR%U, &
                        CSURF_FILETYPE,YALG_MPI)
ENDIF
!
IF (CTIMESERIES_FILETYPE=="OFFLIN") THEN
  CALL GOTO_MODEL(1)
  CALL INIT_OUTPUT_OL_n (YSURF_CUR)
ENDIF
! --------------------------------------------------------------------------------------
!
INW = 1
IF (CTIMESERIES_FILETYPE=="NC    ") INW = 2
!
NWRITE = 0
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
 CALL WLOG_MPI('INIT FINISHED ',PLOG=XTIME)
#ifdef SFX_MPI
XTIME0 = MPI_WTIME()
#endif
!*      2.    Temporal loops
!
XTIME_CALC(:) = 0.
XTIME_WRITE(:) = 0.
!
DO JFORC_STEP=1,INB_STEP_ATM
  !
#ifdef SFX_MPI
  XTIME1 = MPI_WTIME()
#endif
  ! read Forcing
  !
  !indice of forcing line in forcing arrays
  ID_FORC=JFORC_STEP-INT(JFORC_STEP/INB_LINES)*INB_LINES
  IF (ID_FORC==0) ID_FORC=INB_LINES
  !new forcings to read
  IF (ID_FORC==1 .AND. JFORC_STEP.NE.1) THEN
    !if last part of forcing, the last point has to be adjusted on the end of
    !files
    IF (JFORC_STEP/INB_LINES==NB_READ_FORC-1) THEN 
      IDMAX=INB_STEP_ATM-JFORC_STEP+1+1
      !for ascii and binary forcing files
      ZTA(:,IDMAX)=ZTA(:,SIZE(ZTA,2))
      ZQA(:,IDMAX)=ZQA(:,SIZE(ZTA,2))
      ZWIND(:,IDMAX)=ZWIND(:,SIZE(ZTA,2))
      ZDIR_SW(:,IDMAX)=ZDIR_SW(:,SIZE(ZTA,2))
      ZSCA_SW(:,IDMAX)=ZSCA_SW(:,SIZE(ZTA,2))
      ZLW(:,IDMAX)=ZLW(:,SIZE(ZTA,2))
      ZSNOW(:,IDMAX)=ZSNOW(:,SIZE(ZTA,2))
      ZRAIN(:,IDMAX)=ZRAIN(:,SIZE(ZTA,2))
      ZPS(:,IDMAX)=ZPS(:,SIZE(ZTA,2))
      ZCO2(:,IDMAX)=ZCO2(:,SIZE(ZTA,2))
      ZDIR(:,IDMAX)=ZDIR(:,SIZE(ZTA,2))
    ENDIF
    CALL OL_READ_ATM(&
                  CSURF_FILETYPE, CFORCING_FILETYPE, JFORC_STEP,    &
                     ZTA(:,1:IDMAX),ZQA(:,1:IDMAX),ZWIND(:,1:IDMAX), &
                     ZDIR_SW(:,1:IDMAX),ZSCA_SW(:,1:IDMAX),ZLW(:,1:IDMAX), &
                     ZSNOW(:,1:IDMAX),ZRAIN(:,1:IDMAX),ZPS(:,1:IDMAX),&
                     ZCO2(:,1:IDMAX),ZDIR(:,1:IDMAX),LLIMIT_QAIR         )
  ENDIF

#ifdef SFX_MPI
  XTIME_CALC(1) = XTIME_CALC(1) + (MPI_WTIME() - XTIME1)
  XTIME1 = MPI_WTIME()
#endif
  !
  DO JSURF_STEP=1,INB_ATM
    !
    ! time interpolation of the forcing
    !
#ifdef SFX_MPI
    XTIME1 = MPI_WTIME()
#endif
    CALL SUNPOS(ISIZE_OMP, IYEAR, IMONTH, IDAY, ZTIME, &
                ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
    CALL SUNPOS(ISIZE_OMP, IYEAR, IMONTH, IDAY, ZTIME+XTSTEP_SURF, &
                ZLON, ZLAT, XTSUN, XZENITH2, XAZIM)
#ifdef SFX_MPI
    XTIME_CALC(2) = XTIME_CALC(2) + (MPI_WTIME() - XTIME1)
    XTIME1 = MPI_WTIME()
#endif
    !interpolation between beginning and end of current forcing time step
    CALL OL_TIME_INTERP_ATM(JSURF_STEP,INB_ATM,ISIZE_OMP,            &
                            ZTA(:,ID_FORC),ZTA(:,ID_FORC+1),         &
                            ZQA(:,ID_FORC),ZQA(:,ID_FORC+1),         &
                            ZWIND(:,ID_FORC),ZWIND(:,ID_FORC+1),     &
                            ZDIR_SW(:,ID_FORC),ZDIR_SW(:,ID_FORC+1), &
                            ZSCA_SW(:,ID_FORC),ZSCA_SW(:,ID_FORC+1), &
                            ZLW(:,ID_FORC),ZLW(:,ID_FORC+1),         &
                            ZSNOW(:,ID_FORC+1),ZRAIN(:,ID_FORC+1),   &
                            ZPS(:,ID_FORC),ZPS(:,ID_FORC+1),         &
                            ZCO2(:,ID_FORC), ZCO2(:,ID_FORC+1),      &
                            ZDIR(:,ID_FORC) ,ZDIR(:,ID_FORC+1)       )
#ifdef SFX_MPI
    XTIME_CALC(3) = XTIME_CALC(3) + (MPI_WTIME() - XTIME1)
    XTIME1 = MPI_WTIME()
#endif
    !
    IF(LADAPT_SW)THEN
      !
      ! coherence between solar zenithal angle and radiation
      ! when solar beam close to horizontal -> reduction of direct radiation to
      ! the benefit of scattered radiation
      ! when pi/2 - 0.1 < ZENITH < pi/2 - 0.05 => weight of direct to scattered radiation decreases linearly with zenith 
      ! when pi/2 - 0.05 < ZENITH => all the direct radiation is converted to scattered radiation
      ! coherence between solar zenithal angle and radiation
      !
      ZCOEF(:) = (XPI/2. - XZENITH(:) - 0.05) / 0.05
      ZCOEF(:) = MAX(MIN(ZCOEF,1.),0.)
      DO JLOOP=1,SIZE(XDIR_SW,2)
        XSCA_SW(:,JLOOP) = XSCA_SW(:,JLOOP) + XDIR_SW(:,JLOOP) * (1 - ZCOEF)
        XDIR_SW(:,JLOOP) = XDIR_SW(:,JLOOP) * ZCOEF(:)
      ENDDO
      !
    ELSE
      !
      ZSW(:) = 0.
      DO JLOOP=1,SIZE(XDIR_SW,2)
        ZSW(:) = ZSW(:) + XDIR_SW(:,JLOOP) + XSCA_SW(:,JLOOP)
      END DO
      WHERE (ZSW(:)>0.)
        XZENITH  = MIN (XZENITH ,XPI/2.-0.01)
        XZENITH2 = MIN (XZENITH2,XPI/2.-0.01)
      ELSEWHERE
        XZENITH  = MAX (XZENITH ,XPI/2.)
        XZENITH2 = MAX (XZENITH2,XPI/2.)
      END WHERE
      !
    ENDIF       
    !
    ! updates time
    ZTIMEC= ZTIMEC+XTSTEP_SURF
    IF (LCOUPL_TOPD) LTOPD_STEP = ( MOD((((JFORC_STEP-1)*INB_ATM)+JSURF_STEP),NNB_TOPD) == 0 )
    !
    ! run Surface
    !
#ifdef SFX_MPI
    XTIME_CALC(4) = XTIME_CALC(4) + (MPI_WTIME() - XTIME1)
#endif
    !
    CALL IO_BUFF_CLEAN
    !
    IF(LOASIS)THEN
     ! Receive fields to other models proc by proc
     CALL SFX_OASIS_RECV_OL(YSURF_CUR%FM%F, YSURF_CUR%IM%I, YSURF_CUR%SM%S, YSURF_CUR%U, YSURF_CUR%WM%W, &
                            CSURF_FILETYPE,INI,IBANDS,ZTIMEC,       &
                            XTSTEP_SURF,ISIZE_OMP,XZENITH,XSW_BANDS,&
                            XTSRAD,XDIR_ALB,XSCA_ALB,XEMIS,XTSURF   )
    ENDIF
    !
!$OMP PARALLEL PRIVATE(INKPROMA,XTIME1)
    !
#ifdef SFX_MPI
    XTIME1 = MPI_WTIME()
#endif
    !
!$ NBLOCK = OMP_GET_THREAD_NUM()
    !
    IF (NBLOCK==NBLOCKTOT) THEN
      CALL INIT_DIM(ISIZE_OMP,0,INKPROMA,NINDX1SFX,NINDX2SFX)
    ELSE
      CALL INIT_DIM(ISIZE_OMP,NBLOCK,INKPROMA,NINDX1SFX,NINDX2SFX)
    ENDIF
    !
    IF (NBLOCK==0) THEN
      CALL GOTO_MODEL(NBLOCKTOT)
    ELSE
      CALL GOTO_MODEL(NBLOCK)
    ENDIF
    !
    IF(GSHADOWS) THEN 
      CALL SLOPE_RADIATIVE_EFFECT(XTSTEP_SURF,XZENITH(NINDX1SFX:NINDX2SFX),&
      XAZIM(NINDX1SFX:NINDX2SFX),XPS(NINDX1SFX:NINDX2SFX),XTA(NINDX1SFX:NINDX2SFX),&
      XRAIN(NINDX1SFX:NINDX2SFX),XDIR_SW(NINDX1SFX:NINDX2SFX,:),XLW(NINDX1SFX:NINDX2SFX),&
      XZS_THREAD,XZS_XY_THREAD, XSLOPANG_THREAD,XSLOPAZI_THREAD,XSURF_TRIANGLE_THREAD)
    END IF
    !
    CALL COUPLING_SURF_ATM_n(YSURF_CUR, &
                             CSURF_FILETYPE, 'E', ZTIMEC,                    &
           XTSTEP_SURF, IYEAR, IMONTH, IDAY, ZTIME, INKPROMA, NSCAL, IBANDS, &
           XTSUN(NINDX1SFX:NINDX2SFX), XZENITH(NINDX1SFX:NINDX2SFX),         &
           XZENITH2(NINDX1SFX:NINDX2SFX), XAZIM(NINDX1SFX:NINDX2SFX),        &
           XZREF(NINDX1SFX:NINDX2SFX), XUREF(NINDX1SFX:NINDX2SFX),           &
           XZS(NINDX1SFX:NINDX2SFX), XU(NINDX1SFX:NINDX2SFX),                &
           XV(NINDX1SFX:NINDX2SFX), XQA(NINDX1SFX:NINDX2SFX),                &
           XTA(NINDX1SFX:NINDX2SFX), XRHOA(NINDX1SFX:NINDX2SFX),             &
           XSV(NINDX1SFX:NINDX2SFX,:), XCO2(NINDX1SFX:NINDX2SFX), CSV,       &
           XRAIN(NINDX1SFX:NINDX2SFX),  XSNOW(NINDX1SFX:NINDX2SFX),          &
           XLW(NINDX1SFX:NINDX2SFX), XDIR_SW(NINDX1SFX:NINDX2SFX,:),          &
           XSCA_SW(NINDX1SFX:NINDX2SFX,:), XSW_BANDS, XPS(NINDX1SFX:NINDX2SFX),&
           XPA(NINDX1SFX:NINDX2SFX), XSFTQ(NINDX1SFX:NINDX2SFX),             &
           XSFTH(NINDX1SFX:NINDX2SFX), XSFTS(NINDX1SFX:NINDX2SFX,:),         &
           XSFCO2(NINDX1SFX:NINDX2SFX), XSFU(NINDX1SFX:NINDX2SFX),           &
           XSFV(NINDX1SFX:NINDX2SFX), XTSRAD(NINDX1SFX:NINDX2SFX),           &
           XDIR_ALB(NINDX1SFX:NINDX2SFX,:), XSCA_ALB(NINDX1SFX:NINDX2SFX,:), &
           XEMIS(NINDX1SFX:NINDX2SFX),                                       &
           XTSURF(NINDX1SFX:NINDX2SFX), XZ0(NINDX1SFX:NINDX2SFX),            &
           XZ0H(NINDX1SFX:NINDX2SFX), XQSURF(NINDX1SFX:NINDX2SFX),           &          
           XPEW_A_COEF(NINDX1SFX:NINDX2SFX),XPEW_B_COEF(NINDX1SFX:NINDX2SFX),&
           XPET_A_COEF(NINDX1SFX:NINDX2SFX),XPEQ_A_COEF(NINDX1SFX:NINDX2SFX),&
           XPET_B_COEF(NINDX1SFX:NINDX2SFX),XPEQ_B_COEF(NINDX1SFX:NINDX2SFX),&
           YTEST                                                             )
    !
    CALL RESET_DIM(INI,INKPROMA,NINDX1SFX,NINDX2SFX)
    !
#ifdef SFX_MPI
    XTIME_CALC(5) = XTIME_CALC(5) + (MPI_WTIME() - XTIME1)
#endif
    !
!$OMP END PARALLEL
    !
#ifdef SFX_MPI
    XTIME1 = MPI_WTIME()
#endif
    !
    IF(LOASIS)THEN
     ! Send fields to other models proc by proc
     CALL SFX_OASIS_SEND_OL(YSURF_CUR%FM%F, YSURF_CUR%IM%I, YSURF_CUR%SM%S, YSURF_CUR%U, YSURF_CUR%WM%W, YSURF_CUR%SM%DGS, &
                            CSURF_FILETYPE,INI,ZTIMEC,XTSTEP_SURF,ISIZE_OMP)
    ENDIF
    !
    ZTIME = ZTIME + XTSTEP_SURF
    CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
#ifdef SFX_MPI
    XTIME_CALC(6) = XTIME_CALC(6) + (MPI_WTIME() - XTIME1)
    !
    XTIME1 =  MPI_WTIME()
#endif
    ! ecrit Surface
    !
    IF ( LCOUPL_TOPD .AND. LTOPD_STEP ) THEN
      !
      IF (.NOT.ALLOCATED(NYEAR))  ALLOCATE(NYEAR(NNB_TOPD_STEP))
      IF (.NOT.ALLOCATED(NMONTH)) ALLOCATE(NMONTH(NNB_TOPD_STEP))
      IF (.NOT.ALLOCATED(NDAY))   ALLOCATE(NDAY(NNB_TOPD_STEP))
      IF (.NOT.ALLOCATED(NH))     ALLOCATE(NH(NNB_TOPD_STEP))
      IF (.NOT.ALLOCATED(NM))     ALLOCATE(NM(NNB_TOPD_STEP))
      !
      NYEAR (NTOPD_STEP) = IYEAR
      NMONTH(NTOPD_STEP) = IMONTH
      NDAY  (NTOPD_STEP) = IDAY
      NH    (NTOPD_STEP) = INT(ZTIME/3600.)
      NM    (NTOPD_STEP) = INT((ZTIME-NH(NTOPD_STEP)*3600.)/60.)
      !
      IF ( NM(NTOPD_STEP)==60 ) THEN
        !
        NM(NTOPD_STEP) = 0
        NH(NTOPD_STEP) = NH(NTOPD_STEP)+1
        !
      ENDIF
      !
      IF ( NH(NTOPD_STEP)==24 ) THEN
        !
        NH  (NTOPD_STEP) = 0
        NDAY(NTOPD_STEP) = NDAY(NTOPD_STEP)+1
        !
        !!AJOUT BEC 
        SELECT CASE (NMONTH(NTOPD_STEP))
          CASE(4,6,9,11)
            IF ( NDAY(NTOPD_STEP)==31 ) THEN
              NMONTH(NTOPD_STEP) = NMONTH(NTOPD_STEP)+1
              NDAY  (NTOPD_STEP) = 1
            ENDIF
          CASE(1,3,5,7:8,10)
            IF ( NDAY(NTOPD_STEP)==32 ) THEN
              NMONTH(NTOPD_STEP) = NMONTH(NTOPD_STEP)+1
              NDAY  (NTOPD_STEP) = 1
            ENDIF
          CASE(12)
            IF ( NDAY(NTOPD_STEP)==32 ) THEN
              NYEAR (NTOPD_STEP) = NYEAR(NTOPD_STEP)+1
              NMONTH(NTOPD_STEP) = 1
              NDAY  (NTOPD_STEP) = 1
            ENDIF
          CASE(2)
            IF( MOD(NYEAR(NTOPD_STEP),4)==0 .AND. MOD(NYEAR(NTOPD_STEP),100)/=0 .OR. MOD(NYEAR(NTOPD_STEP),400)==0 ) THEN 
              IF (NDAY(NTOPD_STEP)==30) THEN
                NMONTH(NTOPD_STEP) = NMONTH(NTOPD_STEP)+1
                NDAY  (NTOPD_STEP) = 1
              ENDIF
            ELSE
              IF (NDAY(NTOPD_STEP)==29) THEN
                NMONTH(NTOPD_STEP) = NMONTH(NTOPD_STEP)+1
                NDAY  (NTOPD_STEP) = 1
              ENDIF
            ENDIF
        END SELECT
        !
      ENDIF
      !
      ! * 2. Stocking date of each time step
      !
      NTOPD_STEP = NTOPD_STEP + 1 
      !     
    ENDIF
    !
    IF (MOD(ZTIMEC,XTSTEP_OUTPUT) == 0. .AND. CTIMESERIES_FILETYPE/='NONE  ') THEN
      !
      IF (NRANK==NPIO) THEN
        !
        !* name of the file
        IF (CTIMESERIES_FILETYPE=="ASCII " .OR. &
            CTIMESERIES_FILETYPE=="LFI   " .OR. &
            CTIMESERIES_FILETYPE=="FA    " .OR. &
            CTIMESERIES_FILETYPE=="NC    "    ) THEN  
          !
          !
          ZTIME_OUT  = ZTIME
          IDAY_OUT   = IDAY
          IMONTH_OUT = IMONTH
          IYEAR_OUT  = IYEAR
          !
          IF(LOUT_TIMENAME)THEN
            ! if true, change the name of output file at the end of a day
            ! (ex: 19860502_00h00 -> 19860501_24h00)                     
            IF(ZTIME==0.0)THEN
              ZTIME_OUT = 86400.
              IDAY_OUT   = IDAY-1
              IF(IDAY_OUT==0)THEN
                IMONTH_OUT = IMONTH - 1
                IF(IMONTH_OUT==0)THEN
                  IMONTH_OUT=12
                  IYEAR_OUT = IYEAR - 1
                ENDIF
                SELECT CASE (IMONTH_OUT)
                  CASE(4,6,9,11)
                    IDAY_OUT=30
                  CASE(1,3,5,7:8,10,12)
                    IDAY_OUT=31
                  CASE(2)
                    IF( ((MOD(IYEAR_OUT,4)==0).AND.(MOD(IYEAR_OUT,100)/=0)) .OR. (MOD(IYEAR_OUT,400)==0))THEN 
                      IDAY_OUT=29
                    ELSE
                      IDAY_OUT=28
                    ENDIF
                END SELECT
              ENDIF
            ENDIF
            !
          ENDIF
          !
          WRITE(YTAG,FMT='(I4.4,I2.2,I2.2,A1,I2.2,A1,I2.2)') IYEAR_OUT,IMONTH_OUT,IDAY_OUT,&
                '_',INT(ZTIME_OUT/3600.),'h',NINT(ZTIME_OUT)/60-60*INT(ZTIME_OUT/3600.)  
          CFILEOUT    = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG//'.txt')
          CFILEOUT_LFI= ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG)
          CFILEOUT_FA = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG//'.fa')
          CFILEOUT_NC = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG//'.nc')
          !
          IF (CTIMESERIES_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
            LFANOCOMPACT = LDIAG_FA_NOCOMPACT
            IDATEF(1)= IYEAR!_OUT
            IDATEF(2)= IMONTH!_OUT
            IDATEF(3)= IDAY!_OUT
            !ZTIME instead of ZTIME_OUT (FA XRD39 do not like 24h)
            IDATEF(4)= FLOOR(ZTIME/3600.)
            IDATEF(5)= FLOOR(ZTIME/60.) - IDATEF(4) * 60 
            IDATEF(6)= NINT(ZTIME) - IDATEF(4) * 3600 - IDATEF(5) * 60
            IDATEF(7:11) = 0
            IF (CSURF_FILETYPE/='FA    ') THEN
              CALL WRITE_HEADER_FA(YSURF_CUR%UG, &
                                   CSURF_FILETYPE,'ALL')
            ELSE
              CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'UNKNOWN',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
            ENDIF
            CALL FANDAR(IRET,NUNIT_FA,IDATEF)
#endif
          END IF
          !
        END IF
        !
        XSTARTW = XSTARTW + 1
        NWRITE  = NWRITE  + 1
        LTIME_WRITTEN(:)=.FALSE.
        !
      ENDIF
      !
#ifdef SFX_MPI
      XTIME_WRITE(1) = XTIME_WRITE(1) + (MPI_WTIME() - XTIME1)
#endif
      !
      LDEF = .TRUE.
      !
      IF (CTIMESERIES_FILETYPE=="NC    ") THEN
        CALL GOTO_MODEL(1)
        CALL INIT_OUTPUT_NC_n (YSURF_CUR%TM%BDD, YSURF_CUR%CHE, YSURF_CUR%CHN, YSURF_CUR%CHU, &
                               YSURF_CUR%SM%DTS, YSURF_CUR%TM%DTT, YSURF_CUR%DTZ, YSURF_CUR%IM%I, &
                               YSURF_CUR%UG, YSURF_CUR%U, YSURF_CUR%DGU)                   
      ENDIF
      !
      IDX_W = 0
      !
      DO JNW = 1,INW
        ! 
        CALL IO_BUFF_CLEAN
        !
!$OMP PARALLEL PRIVATE(INKPROMA,XTIME1)
        !
!$ NBLOCK = OMP_GET_THREAD_NUM()
        !
        IF (NBLOCK==NBLOCKTOT) THEN
          CALL INIT_DIM(ISIZE_OMP,0,INKPROMA,NINDX1SFX,NINDX2SFX)
        ELSE
          CALL INIT_DIM(ISIZE_OMP,NBLOCK,INKPROMA,NINDX1SFX,NINDX2SFX)
        ENDIF
        !
        IF (NBLOCK==0) THEN
          CALL GOTO_MODEL(NBLOCKTOT)
        ELSE
          CALL GOTO_MODEL(NBLOCK)
        ENDIF
        !
#ifdef SFX_MPI
        XTIME1 =  MPI_WTIME()
#endif
        CALL WRITE_SURF_ATM_n(YSURF_CUR, &
                              CTIMESERIES_FILETYPE,'ALL',LLAND_USE)
#ifdef SFX_MPI
        XTIME_WRITE(2) = XTIME_WRITE(2) + (MPI_WTIME() - XTIME1)
        XTIME1 =  MPI_WTIME()
#endif
        CALL DIAG_SURF_ATM_n(YSURF_CUR%IM%DGEI, YSURF_CUR%FM%DGF, YSURF_CUR%DGL, YSURF_CUR%IM%DGI, &
                             YSURF_CUR%SM%DGS, YSURF_CUR%DGU, YSURF_CUR%TM%DGT, YSURF_CUR%WM%DGW, &
                             YSURF_CUR%U, YSURF_CUR%USS, &
                             CTIMESERIES_FILETYPE)
#ifdef SFX_MPI
        XTIME_WRITE(3) = XTIME_WRITE(3) + (MPI_WTIME() - XTIME1)
        XTIME1 =  MPI_WTIME()
#endif
        CALL WRITE_DIAG_SURF_ATM_n(YSURF_CUR, &
                                   CTIMESERIES_FILETYPE,'ALL')
#ifdef SFX_MPI
        XTIME_WRITE(4) = XTIME_WRITE(4) + (MPI_WTIME() - XTIME1)
#endif
        !
        CALL RESET_DIM(INI,INKPROMA,NINDX1SFX,NINDX2SFX)
        !
!$OMP BARRIER
        !
!$OMP END PARALLEL
        !
        LDEF = .FALSE.
        !
      ENDDO
      !
      IF (LCOUPL_TOPD .AND. NTOPD_STEP > NNB_TOPD_STEP) THEN
        !
        ! Writing of file resulting of coupling with TOPMODEL or routing ****
        CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_total.txt','FORMATTED',&
                                  NYEAR,NMONTH,NDAY,NH,NM,XQTOT)
        CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_runoff.txt','FORMATTED',&
                                  NYEAR,NMONTH,NDAY,NH,NM,XQB_RUN)
        CALL WRITE_DISCHARGE_FILE(CSURF_FILETYPE,'q_drainage.txt','FORMATTED',&
                                  NYEAR,NMONTH,NDAY,NH,NM,XQB_DR)
        ! Writing of budget files 
        IF (LBUDGET_TOPD) CALL WRITE_BUDGET_COUPL_ROUT
        !
      ENDIF
      !
#ifdef SFX_MPI
      XTIME1 =  MPI_WTIME()
#endif
      !
      IF (NRANK==NPIO) THEN
        IF (CTIMESERIES_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
          CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
#endif
        END IF
        !* add informations in the file
        IF (CTIMESERIES_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH
      ENDIF
#ifdef SFX_MPI
      XTIME_WRITE(5) = XTIME_WRITE(5) + (MPI_WTIME() - XTIME1)
#endif
      !
    ENDIF
    !
  END DO
  !
  
  IF (NRANK==NPIO) THEN
    IF (LPRINT) THEN
      IF (MOD(ZTIMEC,XDAY) == 0.) THEN
!$OMP SINGLE              
        ICOUNT = ICOUNT + 1
        CALL WLOG_MPI('SFX DAY :',KLOG=ICOUNT,KLOG2=INT(ZDURATION/XDAY))
        WRITE(*,'(A10,I5,A2,I5)')'SFX  DAY :',ICOUNT,' /',INT(ZDURATION/XDAY)
!$OMP END SINGLE        
      ENDIF
    ENDIF
  ENDIF
  !
END DO
!
!$OMP PARALLEL PRIVATE(XTIME)
!
#ifdef SFX_MPI
XTIME = (MPI_WTIME() - XTIME0)
#endif
!
 CALL WLOG_MPI(' ')
 CALL WLOG_MPI('OL_READ_ATM ',PLOG=XTIME_CALC(1))
 CALL WLOG_MPI('SUNPOS ',PLOG=XTIME_CALC(2))
 CALL WLOG_MPI('OL_TIME_INTERP_ATM ',PLOG=XTIME_CALC(3))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('ZENITH ',PLOG=XTIME_CALC(4))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('COUPLING_SURF_ATM ',PLOG=XTIME_CALC(5))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('ADD_FORECAST_TO_DATE_SURF ',PLOG=XTIME_CALC(6))
 CALL WLOG_MPI('DEF_DATE ',PLOG=XTIME_WRITE(1))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('WRITE_SURF_ATM ',PLOG=XTIME_WRITE(2))
 CALL WLOG_MPI('DIAG_SURF_ATM ',PLOG=XTIME_WRITE(3))
 CALL WLOG_MPI('WRITE_DIAG_SURF_ATM ',PLOG=XTIME_WRITE(4))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('CLOSE FILES ',PLOG=XTIME_WRITE(5))
 CALL WLOG_MPI('')
 CALL WLOG_MPI('END LOOP ',PLOG=XTIME)
 CALL WLOG_MPI('')
 CALL WLOG_MPI('TIME_NPIO_WRITE ',PLOG=XTIME_NPIO_WRITE)
 CALL WLOG_MPI('TIME_COMM_WRITE ',PLOG=XTIME_COMM_WRITE)
 CALL WLOG_MPI('TIME_OMP_BARR ',PLOG=XTIME_OMP_BARR)
 CALL WLOG_MPI('TIME_CALC_WRITE ',PLOG=XTIME_CALC_WRITE)
 CALL WLOG_MPI('')
 CALL WLOG_MPI('TIME_INIT_SEA ',PLOG=XTIME_INIT_SEA)
 CALL WLOG_MPI('TIME_INIT_WATER ',PLOG=XTIME_INIT_WATER)
 CALL WLOG_MPI('TIME_INIT_NATURE ',PLOG=XTIME_INIT_NATURE)
 CALL WLOG_MPI('TIME_INIT_TOWN ',PLOG=XTIME_INIT_TOWN)
 CALL WLOG_MPI('')
 CALL WLOG_MPI('TIME_SEA ',PLOG=XTIME_SEA)
 CALL WLOG_MPI('TIME_WATER ',PLOG=XTIME_WATER)
 CALL WLOG_MPI('TIME_NATURE ',PLOG=XTIME_NATURE)
 CALL WLOG_MPI('TIME_TOWN ',PLOG=XTIME_TOWN)
!
!$OMP END PARALLEL
!
!
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') CALL OPEN_CLOSE_BIN_ASC_FORC('CLOSE',CFORCING_FILETYPE,'R')
!
IF (CFORCING_FILETYPE=='NETCDF') CALL CLOSE_FILEIN_OL
IF (CTIMESERIES_FILETYPE=='OFFLIN') CALL CLOSE_FILEOUT_OL
!
! --------------------------------------------------------------------------------------
!
!*    3.     write restart file
!            ------------------
!
IF ( LRESTART ) THEN
  !
  IF (NRANK==NPIO) THEN
    !* name of the file
    CFILEOUT    = ADJUSTL(ADJUSTR(CSURFFILE)//'.txt')
    CFILEOUT_LFI= CSURFFILE
    CFILEOUT_FA = ADJUSTL(ADJUSTR(CSURFFILE)//'.fa')
    CFILEOUT_NC = ADJUSTL(ADJUSTR(CSURFFILE)//'.nc')

    !* opens the file
    IF (CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
      LFANOCOMPACT = .TRUE.
      IDATEF(1)= IYEAR
      IDATEF(2)= IMONTH
      IDATEF(3)= IDAY
      IDATEF(4)= FLOOR(ZTIME/3600.)
      IDATEF(5)= FLOOR(ZTIME/60.) - IDATEF(4) * 60 
      IDATEF(6)= NINT(ZTIME) - IDATEF(4) * 3600 - IDATEF(5) * 60
      IDATEF(7:11) = 0    
      CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'UNKNOWN',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
      CALL FANDAR(IRET,NUNIT_FA,IDATEF)
#endif
    END IF
    !
  ENDIF
  !
  INW = 1
  IF (CSURF_FILETYPE=="NC    ") INW = 2
  !
  LDEF = .TRUE.
  !
  IF (CSURF_FILETYPE=="NC    ") THEN
    CALL GOTO_MODEL(1)
    CALL INIT_OUTPUT_NC_n (YSURF_CUR%TM%BDD, YSURF_CUR%CHE, YSURF_CUR%CHN, YSURF_CUR%CHU, &
                         YSURF_CUR%SM%DTS, YSURF_CUR%TM%DTT, YSURF_CUR%DTZ, YSURF_CUR%IM%I, &
                         YSURF_CUR%UG, YSURF_CUR%U, YSURF_CUR%DGU)
  ENDIF
  !
  DO JNW = 1,INW
    !
    CALL IO_BUFF_CLEAN
    !
!$OMP PARALLEL PRIVATE(INKPROMA,JNW)
  !
!$ NBLOCK = OMP_GET_THREAD_NUM()
    !
    IF (NBLOCK==NBLOCKTOT) THEN
      CALL INIT_DIM(ISIZE_OMP,0,INKPROMA,NINDX1SFX,NINDX2SFX)
    ELSE
      CALL INIT_DIM(ISIZE_OMP,NBLOCK,INKPROMA,NINDX1SFX,NINDX2SFX)
    ENDIF
    !
    IF (NBLOCK==0) THEN
      CALL GOTO_MODEL(NBLOCKTOT)
    ELSE
      CALL GOTO_MODEL(NBLOCK)
    ENDIF
    !  
    CALL FLAG_UPDATE(YSURF_CUR%IM%DGI, YSURF_CUR%DGU, &
                     .FALSE.,.TRUE.,.FALSE.,.FALSE.)
    !
    IF (LRESTART_2M) THEN
      I2M       = 1
      GPGD_ISBA = .TRUE.
    ELSE
      I2M       = 0
      GPGD_ISBA = .FALSE.
    ENDIF
    GFRAC                  = .TRUE.
    GDIAG_GRID             = .TRUE.
    GSURF_BUDGET           = .FALSE.
    GRAD_BUDGET            = .FALSE.
    GCOEF                  = .FALSE.
    GSURF_VARS             = .FALSE.
    IBEQ                   = 0
    IDSTEQ                 = 0
    GDIAG_OCEAN            = .FALSE.
    GDIAG_SEAICE           = .FALSE.
    GWATER_PROFILE         = .FALSE.
    GSURF_EVAP_BUDGET      = .FALSE.
    GFLOOD                 = .FALSE.
    GCH_NO_FLUX_ISBA       = .FALSE.
    GSURF_MISC_BUDGET_ISBA = .FALSE.
    GPGD_TEB               = .FALSE.
    GSURF_MISC_BUDGET_TEB  = .FALSE.  
    !
    CALL FLAG_DIAG_UPDATE(YSURF_CUR%FM%CHF, YSURF_CUR%IM%CHI, YSURF_CUR%SM%CHS, YSURF_CUR%TM%CHT, &
                          YSURF_CUR%WM%CHW, YSURF_CUR%IM%DGEI, YSURF_CUR%FM%DGF, YSURF_CUR%IM%DGI, &
                          YSURF_CUR%FM%DGMF, YSURF_CUR%IM%DGMI, YSURF_CUR%TM%DGMTO, YSURF_CUR%SM%DGO, &
                          YSURF_CUR%SM%DGS, YSURF_CUR%SM%DGSI, YSURF_CUR%DGU, YSURF_CUR%TM%DGT, &
                          YSURF_CUR%WM%DGW, YSURF_CUR%IM%I, YSURF_CUR%U, &
                          GFRAC, GDIAG_GRID, I2M, GSURF_BUDGET, GRAD_BUDGET, GCOEF,  &
                          GSURF_VARS, IBEQ, IDSTEQ, GDIAG_OCEAN, GDIAG_SEAICE,       &
                          GWATER_PROFILE,                                            &                          
                          GSURF_EVAP_BUDGET, GFLOOD,  GPGD_ISBA, GCH_NO_FLUX_ISBA,   &
                          GSURF_MISC_BUDGET_ISBA, GPGD_TEB, GSURF_MISC_BUDGET_TEB    )
    !
    !* writes into the file
    CALL WRITE_SURF_ATM_n(YSURF_CUR, &
                              CSURF_FILETYPE,'ALL',LLAND_USE)
    IF(CSURF_FILETYPE/='FA    ' .OR. LRESTART_2M) THEN
       CALL WRITE_DIAG_SURF_ATM_n(YSURF_CUR, &
                                   CSURF_FILETYPE,'ALL')
    ENDIF
    !
    CALL RESET_DIM(INI,INKPROMA,NINDX1SFX,NINDX2SFX)
    !
!$OMP END PARALLEL     
    !
    LDEF = .FALSE.
    !
  ENDDO
  !
  !* closes the file
  IF (NRANK==0 ) THEN
    IF (CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
      CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
#endif
    END IF
    !* add informations in the file
    IF (CSURF_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH
    !
  ENDIF
  !
  IF (LCOUPL_TOPD .AND. NTOPD_STEP > NNB_TOPD_STEP) CALL PREP_RESTART_COUPL_TOPD(YSURF_CUR%UG, YSURF_CUR%U, &
                                                                                 CSURF_FILETYPE,INI)
  !
END IF
!
! --------------------------------------------------------------------------------------
!
!*    4.     inquiry mode
!            ------------
!
IF ( LINQUIRE ) THEN
  !
!$OMP PARALLEL PRIVATE(ZSEA,ZWATER,ZNATURE,ZTOWN,ZT2M,ZQ2M,ZZ0,ZZ0H,ZQS_SEA, &
!$OMP ZQS_WATER,ZQS_NATURE,ZQS_TOWN,ZQS,ZPSNG,ZPSNV,ZZ0EFF,ZZS ,&
!$OMP INKPROMA,INI)
!
!$ NBLOCK = OMP_GET_THREAD_NUM()
!
  IF (NBLOCK==NBLOCKTOT) THEN
    CALL INIT_DIM(ISIZE_OMP,0,INKPROMA,NINDX1SFX,NINDX2SFX)
  ELSE
    CALL INIT_DIM(ISIZE_OMP,NBLOCK,INKPROMA,NINDX1SFX,NINDX2SFX)
  ENDIF
  INI = NINDX2SFX-NINDX1SFX+1
  !
  ALLOCATE( ZSEA       ( INI ) )
  ALLOCATE( ZWATER     ( INI ) )
  ALLOCATE( ZNATURE    ( INI ) )
  ALLOCATE( ZTOWN      ( INI ) )
  ALLOCATE( ZT2M       ( INI ) )
  ALLOCATE( ZQ2M       ( INI ) )
  ALLOCATE( ZZ0        ( INI ) )
  ALLOCATE( ZZ0H       ( INI ) )
  ALLOCATE( ZQS_SEA    ( INI ) )
  ALLOCATE( ZQS_WATER  ( INI ) )
  ALLOCATE( ZQS_NATURE ( INI ) )
  ALLOCATE( ZQS_TOWN   ( INI ) )
  ALLOCATE( ZQS        ( INI ) )
  ALLOCATE( ZPSNG      ( INI ) )
  ALLOCATE( ZPSNV      ( INI ) )
  ALLOCATE( ZZ0EFF     ( INI ) )
  ALLOCATE( ZZS        ( INI ) )
  !
  ISERIES = 0
  CALL GET_SURF_VAR_n(YSURF_CUR%FM%DGF, YSURF_CUR%IM%I, YSURF_CUR%IM%DGI, YSURF_CUR%IM%DGMI, &
                      YSURF_CUR%SM%DGS, YSURF_CUR%DGU, YSURF_CUR%TM%DGT, YSURF_CUR%WM%DGW, &
                      YSURF_CUR%FM%F, YSURF_CUR%UG, YSURF_CUR%U, YSURF_CUR%USS, &
                      CSURF_FILETYPE,INI,ISERIES,PSEA=ZSEA,PWATER=ZWATER,PNATURE=ZNATURE,PTOWN=ZTOWN, &
                        PT2M=ZT2M,PQ2M=ZQ2M,PQS=ZQS,PZ0=ZZ0,PZ0H=ZZ0H,PZ0EFF=ZZ0EFF,PQS_SEA=ZQS_SEA,  &
                        PQS_WATER=ZQS_WATER,PQS_NATURE=ZQS_NATURE,PQS_TOWN=ZQS_TOWN,                  &
                        PPSNG=ZPSNG,PPSNV=ZPSNV,PZS=ZZS                                         )  
  !
  ISIZE = SIZE(NINDEX)
  IF (NRANK==NPIO) THEN
!$OMP SINGLE
    ALLOCATE(ZSEA_FULL   (ISIZE))
    ALLOCATE(ZWATER_FULL (ISIZE))
    ALLOCATE(ZNATURE_FULL(ISIZE))
    ALLOCATE(ZTOWN_FULL  (ISIZE))
    ALLOCATE(ZZ0_FULL    (ISIZE))
    ALLOCATE(ZZ0EFF_FULL (ISIZE))
    ALLOCATE(ZZS_FULL    (ISIZE))
!$OMP END SINGLE
  ELSE
!$OMP SINGLE          
    ALLOCATE(ZSEA_FULL   (0))
    ALLOCATE(ZWATER_FULL (0))
    ALLOCATE(ZNATURE_FULL(0))
    ALLOCATE(ZTOWN_FULL  (0))
    ALLOCATE(ZZ0_FULL    (0))
    ALLOCATE(ZZ0EFF_FULL (0))
    ALLOCATE(ZZS_FULL    (0))
!$OMP END SINGLE    
  ENDIF
  CALL GATHER_AND_WRITE_MPI(ZSEA,ZSEA_FULL) 
  CALL GATHER_AND_WRITE_MPI(ZWATER,ZWATER_FULL)
  CALL GATHER_AND_WRITE_MPI(ZNATURE,ZNATURE_FULL)
  CALL GATHER_AND_WRITE_MPI(ZTOWN,ZTOWN_FULL)
  CALL GATHER_AND_WRITE_MPI(ZZ0,ZZ0_FULL)
  CALL GATHER_AND_WRITE_MPI(ZZ0EFF,ZZ0EFF_FULL)
  CALL GATHER_AND_WRITE_MPI(ZZS,ZZS_FULL)

  IF (NRANK==NPIO) THEN
!$OMP SINGLE
    WRITE(ILUOUT,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING SEA    : ',COUNT( ZSEA_FULL    (:) > 0. ),' / ',ISIZE
    WRITE(ILUOUT,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING WATER  : ',COUNT( ZWATER_FULL  (:) > 0. ),' / ',ISIZE
    WRITE(ILUOUT,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING NATURE : ',COUNT( ZNATURE_FULL (:) > 0. ),' / ',ISIZE
    WRITE(ILUOUT,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING TOWN   : ',COUNT( ZTOWN_FULL   (:) > 0. ),' / ',ISIZE
    WRITE(ILUOUT,*)'ZZ0    = ',ZZ0_FULL
    WRITE(ILUOUT,*)'ZZ0EFF = ',ZZ0EFF_FULL
    WRITE(ILUOUT,*)'ZZS = ',ZZS_FULL
    WRITE(ILUOUT,*)'MINVAL(ZZS) = ',MINVAL(ZZS_FULL),' MAXVAL(ZZS) = ',MAXVAL(ZZS_FULL)
!$OMP END SINGLE
  ENDIF
  !
  DEALLOCATE( ZSEA       )
  DEALLOCATE( ZWATER     )
  DEALLOCATE( ZNATURE    )
  DEALLOCATE( ZTOWN      )
  DEALLOCATE( ZT2M       )
  DEALLOCATE( ZQ2M       )
  DEALLOCATE( ZZ0        )
  DEALLOCATE( ZZ0H       )
  DEALLOCATE( ZQS_SEA    )
  DEALLOCATE( ZQS_WATER  )
  DEALLOCATE( ZQS_NATURE )
  DEALLOCATE( ZQS_TOWN   )
  DEALLOCATE( ZQS        )
  DEALLOCATE( ZPSNG      )
  DEALLOCATE( ZPSNV      )
  DEALLOCATE( ZZ0EFF     )
  DEALLOCATE( ZZS        )
  !
  IF (NRANK==NPIO) THEN
!$OMP SINGLE
    DEALLOCATE(ZSEA_FULL   )
    DEALLOCATE(ZWATER_FULL )
    DEALLOCATE(ZNATURE_FULL)
    DEALLOCATE(ZTOWN_FULL  )
    DEALLOCATE(ZZ0_FULL    )
    DEALLOCATE(ZZ0EFF_FULL )
    DEALLOCATE(ZZS_FULL    )
!$OMP END SINGLE
  ENDIF
  !
!$OMP END PARALLEL
  !
ENDIF   
!
! --------------------------------------------------------------------------------------
!
!    4'    Close Gelato specific diagnostic 
#if ! defined in_arpege
 CALL CLSDIA()
#endif
!
!
!*    5.     Close parallelized I/O
!            ----------------------
!
IF (NRANK==NPIO) THEN
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) '    --------------------------'
  WRITE(ILUOUT,*) '    | OFFLINE ENDS CORRECTLY |'
  WRITE(ILUOUT,*) '    --------------------------'
  WRITE(ILUOUT,*) ' '
!$OMP SINGLE
  CLOSE(ILUOUT)
!$OMP END SINGLE  
  WRITE(*,*) ' '
  WRITE(*,*) '    --------------------------'
  WRITE(*,*) '    | OFFLINE ENDS CORRECTLY |'
  WRITE(*,*) '    --------------------------'
  WRITE(*,*) ' '
ENDIF
!
 CALL SURFEX_DEALLO_LIST
!
IF (ALLOCATED(NINDEX)) DEALLOCATE(NINDEX)
IF (ALLOCATED(NSIZE_TASK)) DEALLOCATE(NSIZE_TASK)
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
 CALL END_LOG_MPI
!
IF (LHOOK) CALL DR_HOOK('OFFLINE',1,ZHOOK_HANDLE)
!
! * MPI and OASIS must be finalized after the last DR_HOOK call
!
IF(LOASIS)THEN
  CALL SFX_OASIS_END
ELSE
#ifdef SFX_MPI
 CALL MPI_FINALIZE(INFOMPI)
#endif
ENDIF
!
! --------------------------------------------------------------------------------------
!
END PROGRAM OFFLINE
