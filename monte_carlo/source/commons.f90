MODULE COMMONS

    USE ISO_C_BINDING, ONLY: C_DOUBLE
    IMPLICIT NONE
    SAVE

!   Parameters used to control the precision of floating point numbers used during the simulation.
!   Currently set to correspond to double and quadruple precision floating point numbers (DP and QP).
    INTEGER,PARAMETER          :: DP  = SELECTED_REAL_KIND(15,307)
    INTEGER,PARAMETER          :: QP  = SELECTED_REAL_KIND(33,4931)
!   Parameter used to define double precision numbers associated with the random number generator 
!   DRAND48(), which is an external C function.
    INTEGER,PARAMETER          :: CDP = C_DOUBLE
!   Integer and floating point numbers definied in the correct double precision to be used throught
!   the program.
    REAL(KIND=DP), PARAMETER   :: PI=4.0_dp*ATAN(1.0_dp),TWOPI=2.0_dp*PI,FORPI=4.0_dp*PI,HLFPI=PI/2.0_dp,INVRPI=1.D0/SQRT(PI)

    CHARACTER                  :: CHAR(80)

    ! Parameters used by readinput.f90 to read input files
    INTEGER :: LNGTH(20), LOC(20), NITEM

    INTEGER, PARAMETER :: MYUNIT = 11, MYUNIT2 = 21, VIEWUNIT = 31, INPUNIT = 15

    ! Seed for the random number generator
    INTEGER :: ISEED
    LOGICAL :: SEEDT

!=============================================================================================
!   General Monte Carlo parameters
!=============================================================================================
    INTEGER       :: ICNT, NSTEP, NEQ, BLKLNGTH, DMPFRQ, POSDMP, ISTEP, IBLOCK, EID, INDXP
    INTEGER       :: ACCPTCT, ACCPTCR, ACCPTV, NTMOVES, NRMOVES, NVMOVES, CHCK_DIV
    REAL(KIND=DP) :: BETAKB, TEMP, PRSFIX, MAXDTR, MAXDRT, MAXBOX, ACCRATC, ACCRATV, DIV_TOL, START_TIME, MODT0
    LOGICAL       :: MCT, ARATIOT, NVTT, NPTT, ISOTROPICT, SPET, TPET, VERBOSET, PESRFT, MODT, SAVEREFT, UNITVECT, CONTINUET
    LOGICAL       :: SEGREGATET
!   --------------------------------------------------------------------------------
!   Observables to be averaged over the course of standard Monte Carlo simulations
!   --------------------------------------------------------------------------------
    ! Potential energy
    REAL(KIND=DP) :: SUMPE, SUMPE2, AVPE, AVPE2, SUMPEBLK, SUMPE2BLK, AVPEBLK, AVPE2BLK, STDPESUM, STDPE2SUM, STDPE, STDPE2
    ! Volume
    REAL(KIND=DP) :: SUMVLM,SUMVLM2,AVVLM,AVVLM2,SUMVLMBLK,SUMVLM2BLK,AVVLMBLK,AVVLM2BLK,STDVLMSUM,STDVLM2SUM,STDVLM,STDVLM2
    ! Density
    REAL(KIND=DP) :: SUMRHO, AVRHO, SUMRHOBLK, AVRHOBLK, STDRHOSUM, STDRHO
    ! Pressure
    REAL(KIND=DP) :: SUMPRES, AVPRES, SUMPRESBLK, AVPRESBLK, STDPRESSUM, STDPRES
    ! Simulation cell/box
    REAL(KIND=DP) :: SUMBOX(3), AVBOX(3), SUMBOXBLK(3), AVBOXBLK(3), STDBOXSUM(3), STDBOX(3)
    ! Enthalpy
    REAL(KIND=DP) :: SUMH,SUMH2,AVH,AVH2,SUMHBLK,SUMH2BLK,AVHBLK,AVH2BLK,STDH,STDH2,SUMHV,AVHV,SUMHVBLK,AVHVBLK,STDHV,&
                     STDHSUM,STDH2SUM, STDHVSUM

!   ----------------------------------------------------------------------------------
!   Observables to be averaged over the course of free energy Monte Carlo simulations
!   ----------------------------------------------------------------------------------
    ! Largest crystalline cluster
    REAL(KIND=DP) :: SUMLRGSTXCLSTR, AVLRGSTXCLSTR, SUMLRGSTXCLSTRBLK, AVLRGSTXCLSTRBLK, STDLRGSTXCLSTRSUM, STDLRGSTXCLSTR
    ! Chemical potential / Interface pinning
    REAL(KIND=DP) :: DELMUBLK, SUMTRQ, SUMDELMU, AVTRQ, AVDELMU, SUMTRQBLK, AVTRQBLK, STDTRQSUM, STDDELMUSUM, STDTRQ, STDDELMU
    ! Observables for Einstein crystal simulations
    REAL(KIND=DP) :: SUM_EXP_U_EIN, SUM_EXP_U_EIN_BLK, AV_EXP_U_EIN_BLK, AV_EXP_U_EIN, STD_EXP_U_EIN_SUM, STD_EXP_U_EIN
    REAL(KIND=DP) :: SUM_U_EIN_TR, SUM_U_EIN_TR_BLK, AV_U_EIN_TR_BLK, AV_U_EIN_TR, STD_U_EIN_TR_SUM, STD_U_EIN_TR
    REAL(KIND=DP) :: SUM_U_EIN_OR, SUM_U_EIN_OR_BLK, AV_U_EIN_OR_BLK, AV_U_EIN_OR, STD_U_EIN_OR_SUM, STD_U_EIN_OR
    ! Observables for Schilling-Schmid simulations
    REAL(KIND=DP) :: SUM_U_SS, SUM_EXP_U_SS, AV_U_SS, AV_EXP_U_SS, STD_U_SS_SUM, STD_EXP_U_SS_SUM
    REAL(KIND=DP) :: SUM_U_SS_BLK, SUM_EXP_U_SS_BLK, AV_U_SS_BLK, AV_EXP_U_SS_BLK, STD_U_SS, STD_EXP_U_SS

    REAL(KIND=DP), ALLOCATABLE  :: SUMNUCCOUNT(:), AVNUCCOUNTBLK(:), AVNUCCOUNT(:)
    REAL(KIND=DP), ALLOCATABLE  :: STDNUCCOUNT(:), STDNUCCOUNTSUM(:)

    ! Parameters for cell list 
    REAL(KIND=DP), ALLOCATABLE  :: SCALED_R(:,:)
    LOGICAL                     :: CELLLISTT

    ! General system parameters
    INTEGER                     :: NDIM, NPART, HLFPART, NRBSITE, LID, NUCX, NUCY, NUCZ, BNRYA, NLAYERS
    REAL(KIND=DP), ALLOCATABLE  :: R(:,:), Q(:,:)
    REAL(KIND=DP)               :: PE, PEPP, VLM, BOX(3), BOXL, RHO, RCUT, RCUTSQ, VIR, VIRTEMP, PRES, BNRYR
    LOGICAL                     :: LATTICET, DENSITYT, PACKINGT, CUBICT, ORTHORHOMBICT, RIGIDT
    LOGICAL                     :: RANDQUATT, SETORTN, SCALET
    LOGICAL                     :: HEADTAILT, BINARYT, BNRYRTIOT

!=============================================================================================
!   PARAMETERS FOR CLUSTER MOVE MONTE CARLO SIMULATIONS
!=============================================================================================
!*********************************************************************************************
! Parameters related to virtual move Monte Carlo simulations
!*********************************************************************************************

!*********************************************************************************************
! Parameters related to cluster-move Monte Carlo simulations
!*********************************************************************************************
    LOGICAL                     :: CLUSTERT, CLUSTERMOVET, VLMCLUSTERMOVET, LRGCLSTRT, LRGCLSTMVT
    INTEGER                     :: CLSTRID, CLSTRSZ, ACCPTCTC, ACCPTCRC, NTCMOVES, NRCMOVES, NCLSTRS, CLSTRSITEID
    INTEGER, ALLOCATABLE        :: CLSTR(:), CLSTRADJ(:,:), PCLSTRID(:)
    REAL(KIND=DP)               :: CLSTRRATIO, MAXDTRC, MAXDRTC, LRGCLSTRRATIO
    REAL(KIND=DP), ALLOCATABLE  :: RCLSTR(:,:), CLSTRCOM(:,:), PCLSTR(:,:), CLURIJ(:,:,:)

!=============================================================================================
!   PARAMETERS RELATED TO THE CALCULATION OF ORDER PARAMETERS
!*********************************************************************************************
!   Parameters related to the calculation of BOOPs and the largest crystalline cluster. 
!*********************************************************************************************
    LOGICAL                     :: BOPT, QLT, QBARLT, GQLT, QLDOTQLT, LRGSTXCLSTRT, TWOLIMS, BOPCLSTRT
    INTEGER                     :: BOP_L, MIN_L, MAX_L, DEL_L, XSTL_NN_MIN, XSTL_NN_MAX, L_INDEX, NPBOP
    INTEGER                     :: XSTL_NN_MIN_1, XSTL_NN_MAX_1, XSTL_NN_MIN_2, XSTL_NN_MAX_2
    REAL(KIND=DP)               :: BOP_CUT, BOP_CUT_SQ, QLMIN, QLMAX, QLMIN_2, QLMAX_2
    REAL(KIND=DP), ALLOCATABLE  :: GQL(:), QL(:,:), QBARL(:,:), QLDOTQL(:,:,:), QLDIST(:,:)

!=============================================================================================
!   PARAMETERS FOR FREE ENERGY CALCULATIONS
!=============================================================================================
!*********************************************************************************************
! Parameters related to interface pinning simulations
!*********************************************************************************************
    LOGICAL                     :: NPZTT, PINT, HALFST, ORTHO_IN, COLLDT
    INTEGER                     :: PINDIM
    REAL(KIND=DP), ALLOCATABLE  :: RREF(:,:), RBREF(:,:,:)
    REAL(KIND=DP)               :: NHLF, NX, NY, NZ, PINKAP, PINA, PINKAP2, PINKAPA, PINDELQ, TRQ, DELMU
    REAL(KIND=DP)               :: ZPIN, FORCEC, FORCEO
    COMPLEX(KIND=DP)            :: RQ

!*********************************************************************************************
! Parameters related to umbrella sampling, seeding and nucleus-size pinning simulations
!*********************************************************************************************
    LOGICAL                     :: NUCSEEDT, EQUISEEDT, GETSEEDT
    INTEGER                     :: NUCSIZE, SEEDSIZE, TRGTSEED, TRJCTYLNGTH, CLSTRSIZE, EQUIID
    INTEGER, ALLOCATABLE        :: NUCCOUNT(:), TOTNUCCOUNT(:), SEEDID(:)
    REAL(KIND=DP), ALLOCATABLE  :: RSEED(:,:), QSEED(:,:)
    REAL(KIND=DP)               :: BIASK, SEEDRADIUS, NPINK, RHOX, RHOF
!*********************************************************************************************
! Parameters related to Frenkel-Ladd / Einstein crystal simulations
!*********************************************************************************************
    LOGICAL                     :: FLFET
    INTEGER                     :: FLID
    REAL(KIND=DP)               :: LAMTR, LAMOR, ITA_EIN, A_FL
    REAL(KIND=DP)               :: Q_EIN_ROT, A_EIN_T, A0
    REAL(KIND=DP)               :: U_LATT, EXP_U_EIN, U_EIN_TR, U_EIN_OR
    REAL(KIND=DP), ALLOCATABLE  :: RCOM(:), DRCOM(:)
!*********************************************************************************************
! Parameters related to Schilling-Schmid Monte Carlo simulations
!*********************************************************************************************
    LOGICAL                     :: SCHSMIT, ROT_SWITCH
    INTEGER                     :: SSID
    REAL(KIND=DP)               :: LAM_SS, RCUT_SS, RCUTSQ_SS, SSMOVERATIO(3), U_SS, EXP_U_SS

!=============================================================================================
!   PARAMETERS FOR DIFFERENT INTERACTION POTENTIALS
!=============================================================================================
!*********************************************************************************************
! Hard sphere parameters
!*********************************************************************************************
    LOGICAL       :: HST, OVERLAPT, HARDT
!*********************************************************************************************  
! Repulsive Yukawa parameters
!*********************************************************************************************
    REAL(KIND=DP) :: YUKKAP
    LOGICAL       :: YUKT
!*********************************************************************************************  
! Generalised Lennard-Jones parameters
!*********************************************************************************************
    REAL(KIND=DP) :: GLJN, TAILV, TAIL2V
    LOGICAL       :: GLJT, TAILCORT, WCAT

!*********************************************************************************************
! Kihara potential (soft repulsive spherocylinders) parameters
!*********************************************************************************************
    LOGICAL                     :: KIHARAT
    REAL(KIND=DP)               :: RLNGTH, HLFLNGTH, DCHECKSQ

!*********************************************************************************************
! General parameters for rigid bodies
!*********************************************************************************************
    ! Number of rigid body sites on the particles
    INTEGER                     :: NSITES
    ! Reference rigid body sites of the patches and current rigid body sites of the patches.
    REAL(KIND=DP), ALLOCATABLE  :: REFSITE(:,:), RBSITES(:,:,:)
!*********************************************************************************************
! Kern-Frenkel model parameters
!*********************************************************************************************
    LOGICAL                     :: KFT
    ! interaction matrix for patch-patch interactions, array of half-opening angles for the patches.
    REAL(KIND=DP), ALLOCATABLE  :: KFIJ(:,:), KFLAM2(:,:), KFDEL(:)
    ! Lambda defines the range of the square-well component of the potential
    REAL(KIND=DP)               :: KFLAMA, KFLAMB, KFLAMC, KFLAMD
    ! Define the well-depth associated with each of the patch-patch interactions
    REAL(KIND=DP)               :: KFAA, KFBB, KFCC, KFDD
    ! Define the half-opening angles of each of the patches on the particles
    REAL(KIND=DP)               :: KFDELA, KFDELB, KFDELC, KFDELD 
!------------------------------------------------------------------------------------------
!   Parameters related to quasi-2D systems interacting with a hydrophobic surface.
!------------------------------------------------------------------------------------------
    LOGICAL                     :: TDSRFT, GEN2D
    ! Define the two half-opening angles for each rectangular patch
    REAL(KIND=DP)               :: SURFZ
!------------------------------------------------------------------------------------------
!   Parameters related to systems under spherical confinement.
!------------------------------------------------------------------------------------------
    LOGICAL                     :: SPHERECNFT
    ! Define the two half-opening angles for each rectangular patch
    REAL(KIND=DP)               :: SPHERERAD, SPHERERAD2

!*********************************************************************************************
! Patchy Generalised Lennard-Jones model parameters (recycles parameters defined for the plain
! generalised Lennard-Jones potential).
! This model is taken from: The stability of a crystal with diamond structure for patchy 
! particles with tetrahedral symmetry, Eva Noya et al., J. Chem. Phys., 132, 2010.
! This model is only implemented to benchmark the interface pinning method for anisotropic
! particles. 
!*********************************************************************************************
    LOGICAL                     :: PGLJT
    ! This parameter is a measure of the width of the patches, with 2*sqrt(2)*SIGPW 
    ! being the full width at half maximum of the Gaussian.
    REAL(KIND=DP)               :: SIGPW, TSIGPW2

!*********************************************************************************************
! Continuous patchy particle model parameters
!*********************************************************************************************
    LOGICAL                     :: CPPT
    ! interaction matrix for patch-patch interactions, array of half-opening angles for the patches.
    REAL(KIND=DP), ALLOCATABLE  :: CPPIJ(:,:), CPPDEL(:), CPPMDEL(:)
    ! Lambda defines the range of the square-well component of the potential
    REAL(KIND=DP)               :: CPPLAM, CPPLAM2, CPPS, CPPINVS, PIS
    ! Define the well-depth associated with each of the patch-patch interactions
    REAL(KIND=DP)               :: CPPAA, CPPBB, CPPCC, CPPDD
    ! Define the half-opening angles of each of the patches on the particles
    REAL(KIND=DP)               :: CPPDELA, CPPDELB, CPPDELC, CPPDELD
!------------------------------------------------------------------------------------------
!   Parameters related to elongated triblock patchy particles, where the patch-patch 
!   interactions are the same as the CPP model.
!------------------------------------------------------------------------------------------
    LOGICAL                     :: ETPT, DNAT, HTPRT, POLYDT, READPOLYT
    REAL(KIND=DP)               :: ETPCUTSQ, POLYSTD
    REAL(KIND=DP), ALLOCATABLE  :: ETPDEL(:,:), ETPMDEL(:,:)

!*********************************************************************************************
! Hard Dumbbells parameters
!*********************************************************************************************
    REAL(KIND=DP)               :: LSTAR
    LOGICAL                     :: HDMBLT

!*********************************************************************************************
! Generalised Lennard-Jones Dumbbells parameters
!*********************************************************************************************
    REAL(KIND=DP)               :: SIGAA, SIGBB, SIGCC, EPSAA, EPSBB, EPSCC
    REAL(KIND=DP), ALLOCATABLE  :: SIGIJ(:,:), SIGIJ2(:,:), EPSIJ(:,:), DMBLCUTSQ(:,:)
    LOGICAL                     :: DMBLGLJT

!*********************************************************************************************
! Dipolar Discotic Liquid Crystalline particles parameters
!*********************************************************************************************
    LOGICAL                     :: GBT, EWALDT
    INTEGER                     :: NC, NCSQMAX
    REAL(KIND=DP)               :: GBK, GBKP, GBV, GBM, GBCUT, GBX, GBXP
    REAL(KIND=DP)               :: DPMU, DPMUSQ, ALPHA, ALPSQ, GUFCTR, SLFFCT
    INTEGER, ALLOCATABLE        :: NVV(:,:,:)
    REAL(KIND=DP), ALLOCATABLE  :: FCTR(:,:,:)


!*********************************************************************************************
! Stockmayer Fluid parameters
!*********************************************************************************************
    ! REAL(KIND=DP)               :: 
    LOGICAL                     :: STCKMYRT

!*********************************************************************************************
! Parameters for a multicomponent system of DNACCs
!*********************************************************************************************
     
END MODULE COMMONS 
