!     PASSION: A package for Soft Matter Simulation
!     Copyright (C) 2015 - Dwaipayan Chakrabarti
!     This file is part of PASSION.
!
SUBROUTINE KEYWORDS
    !     This subroutine reads in the input parameters, which are driven by keywords.
        USE COMMONS
    
        IMPLICIT NONE
    
        INTEGER                 :: I, IL, J1, J2, EOF, CSHFT
        REAL(KIND=DP)           :: CHAINE
        CHARACTER (LEN = 16)    :: WORD
        LOGICAL                 :: ENDT
    
        ENDT  = .FALSE.
    
        OPEN (UNIT = INPUNIT, FILE = 'input.inp', STATUS = 'UNKNOWN')
    
    !     The END argument of the READ statement has a label number after it to which the program will
    !     branch, if it encounters the end of the file while reading.
    
    200 DO WHILE (.NOT. ENDT)
    
            READ (INPUNIT, FMT = 301, END = 900) (CHAR(I), I = 1, 120)
            
    301     FORMAT (120A1)

            CALL READITEM(WORD)

            IF (WORD == 'COMMENT') THEN
                GO TO 200
    
            ELSE IF (WORD == 'ARATIO') THEN
                ARATIOT = .TRUE.
                CALL READF(2, ACCRATC)
                IF (NITEM > 2) CALL READF(3, ACCRATV)
    
            ELSE IF (WORD == 'BETA') THEN
                CALL READF(2, BETAKB)       ! defines the temperature of the system

            ELSE IF (WORD == 'BINARY') THEN           ! consider a binary system where half the particles do not interact with one another.
                BINARYT = .TRUE.
                IF (NITEM > 1) THEN
                    CALL READI(2, BINARYID)
                ELSE
                    BINARYID = 0
                ENDIF
                IF (NITEM > 2) THEN
                    BNRYRTIOT = .TRUE.
                    CALL READF(3, BNRYR)              ! Relative ratio between A and B particles
                    BNRYA = INT(BNRYR*REAL(NPART,DP))
                ELSE
                    BNRYRTIOT = .FALSE.
                ENDIF
    
            ELSE IF (WORD == 'BOP' .OR. WORD == 'CBOP') THEN          ! instruct to consider a cubic box
                BOPT = .TRUE.
                CALL READI(2, MIN_L)
                CALL READI(3, MAX_L)
                CALL READI(4, DEL_L)
                CALL READL(5, QLT)
                CALL READL(6, QBARLT)
                CALL READL(7, GQLT)
                CALL READL(8, QLDOTQLT)
                CALL READF(9, BOP_CUT)
                IF (NITEM > 9) CALL READL(10, PATCHBST)
                IF (WORD == 'CBOP') BOPCLSTRT = .TRUE.
    
            ELSE IF (WORD == 'CELL') THEN           ! use a cell-list to speed up simulations
                CELLLISTT = .TRUE.
    
            ELSE IF (WORD == 'CLUSTER') THEN        ! instruct the system to perform cluster moves
                CLUSTERMOVET = .TRUE.
                CALL READF(2, CLSTRRATIO)           ! ratio of cluster moves to single particle moves
                IF (NITEM > 2) THEN
                    CALL READI(3, CLSTRSITEID)
                ELSE
                    CLSTRSITEID = 1
                ENDIF

                IF (NITEM > 3) THEN
                    LRGCLSTRT = .TRUE.
                    CALL READF(4, LRGCLSTRRATIO)
                ENDIF

            ELSE IF (WORD == 'CONTINUE') THEN           ! use a cell-list to speed up simulations
                CONTINUET = .TRUE.

            ELSE IF (WORD == 'CONFLINK') THEN           ! print conflink file along simulation
                PRNTCNF = .TRUE.
                CALL READI(2, CNFFRQ)
    
            ELSE IF (WORD == 'CUBIC') THEN          ! instruct to consider a cubic box
                CUBICT = .TRUE.
    
            ELSE IF (WORD == 'CUTOFF') THEN         ! specifies the cutoff distance for the potential
                CALL READF(2, RCUT)
                RCUTSQ = RCUT*RCUT

            ELSE IF (WORD == 'OCYLINDER') THEN
                CALL READF(2, CYLR)                   ! defines the density of the system
                IF (NITEM > 2) CALL READL(3, GEN_CONFIGT)
                OBJCTT  = .TRUE.
                OBJCYLT = .TRUE.
                CYLR2 = CYLR*CYLR
    
            ELSE IF (WORD == 'DENSITY') THEN
                CALL READF(2, RHO)                   ! defines the density of the system
                DENSITYT = .TRUE.
            
            ELSE IF (WORD == 'DIM') THEN             ! instructs that an atomic system is under consideration
                CALL READI(2, NDIM)                  ! defines the number of dimensions for the system under consideration
    
            ELSE IF (WORD == 'DIVERGENCE') THEN      ! parameters used to check how much the energy evaluated along the 
                CALL READI(2, CHCK_DIV)              ! MC trajectory has diverged relative to the true value.
                CALL READF(3, DIV_TOL)               ! Tolerance value for assessing whether the energy has diverged

            ELSE IF (WORD == 'DNA') THEN           ! consider a system with specific interactions.
                DNAT = .TRUE.
    
            ELSE IF (WORD == 'ENSEMBLE') THEN       ! specifies the ensemble under consideration
                CALL READI(2, EID)
                IF (EID == 1) NVTT = .TRUE.          ! canonical
                IF (EID == 2) NPTT = .TRUE.          ! isothermal-isobaric isotropic
                IF (EID == 3) NPTT = .TRUE.          ! isothermal-isobaric anisotropic
                IF (EID == 4) NPTT = .TRUE.          ! isothermal-isobaric (NPzT)

                IF(NPTT) THEN
                    
                    IF(EID==2) THEN
                        ISOTROPICT = .TRUE.

                    ELSEIF(EID==3) THEN
                        ISOTROPICT = .FALSE.

                    ELSEIF(EID==4) THEN
                        NPZTT = .TRUE.
                        ISOTROPICT = .FALSE.
                        CALL READI(3, PINDIM)
                        IF (NITEM > 3) THEN
                            COLLDT = .TRUE.
                            CALL READF(4, NX)
                            CALL READF(5, NY)
                            CALL READF(6, NZ)
                            NX   = 2.0_dp*PI*NX
                            NY   = 2.0_dp*PI*NY
                            NZ   = 2.0_dp*PI*NZ
                            NHLF = SQRT(REAL(NPART,DP))
                        ENDIF

                    ENDIF
                ENDIF
    
            ELSE IF (WORD == 'FL_FE') THEN 
            !---------------------------------------------------------------------------------------------------------------- 
            ! Perform a Frenkel-Ladd (aka Einstein crystal) simulation to calculate the free-energy of a given crystal
            ! structure. 
            !---------------------------------------------------------------------------------------------------------------- 
                CALL READI(2, FLID)
                CALL READF(3, LAMTR)         
                IF (NITEM > 3) CALL READF(4, ITA_EIN) 
                FLFET    = .TRUE.
                SAVEREFT = .TRUE.

            ELSE IF (WORD == 'IPIN') THEN    
            !----------------------------------------------------------------------------------------------------------------    
            ! Perform an interface-pinning simulation, where the system is initiated with a configuration in an orthorhombic
            ! box, where particles with rz<0 are solid-like and those with rz>0 liquid-like. Two order parameters can be
            ! chosen: (i) collective density field; (ii) global Steinhardt bond-orientational order-parameter. The first is
            ! useful for structures which cannot undergo major long-wavelength displacements (i.e., close-packed structures).
            ! To use the second, The (C)BOP keyword must also be invoked. 
            !----------------------------------------------------------------------------------------------------------------
                CALL READF(2, PINKAP)             
                CALL READF(3, PINA)
                CALL READF(4, PINDELQ)
                PINT = .TRUE.
                PINKAP2 = PINKAP/2.0_dp
                PINKAPA = PINKAP*PINA

            ELSE IF (WORD == 'LATTICE' .OR. WORD == 'ALATTICE') THEN
                CALL READI(2, LID)
                
                IF (LID < 34) THEN              ! instructs to initialise from orthorhombic FCC or Diamond lattice 
                    LATTICET = .TRUE.
                    CALL READI(3,NUCX)
                    IF (NITEM > 3) THEN
                !   Orthorhmbic box with different number of unit cells in each direction
                        CALL READI(4,NUCY)
                        CALL READI(5,NUCZ)
                    ELSE
                !   Cubic box with same number of unit cells in each direction
                        NUCY = NUCX
                        NUCZ = NUCX
                    ENDIF
                ELSE
                    WRITE (*, '(A)') ' ERROR: LATTICE ID NOT RECOGNISED'
                    STOP                
                ENDIF

                IF(WORD == 'ALATTICE') ANYCELLST = .TRUE.
    
            ELSE IF (WORD == 'LRGSTXCLSTR') THEN 
            !   Calculate the largest crystalline cluster in the system
                LRGSTXCLSTRT = .TRUE.
                CALL READI(2, BOP_L)                    ! l-value for ql_dot_ql BOP to be used
                CALL READF(3, QLMIN)                    ! Minimum value for ql_dot_ql to be considered a crystalline bond
                CALL READF(4, QLMAX)                    ! Maximum value for ql_dot_ql to be considered a crystalline bond
                CALL READI(5, XSTL_NN_MIN)              ! Number of crystalline bonds a particle needs to particles needs to be crystalline
                CALL READI(6, XSTL_NN_MAX)              ! Number of crystalline bonds a particle needs to particles needs to be crystalline
                IF (NITEM > 6) THEN
                    TWOLIMS = .TRUE.                    ! Use a second set of limits to check for crystalline bonds
                    CALL READF(7, QLMIN_2)              ! Minimum value for ql_dot_ql to be considered a crystalline bond
                    CALL READF(8, QLMAX_2)              ! Maximum value for ql_dot_ql to be considered a crystalline bond
                    CALL READI(9, XSTL_NN_MIN_1)        ! Number of crystalline bonds a particle needs to particles needs to be crystalline
                    CALL READI(10, XSTL_NN_MAX_1)       ! Number of crystalline bonds a particle needs to particles needs to be crystalline
                    CALL READI(11, XSTL_NN_MIN_2)       ! Number of crystalline bonds a particle needs to particles needs to be crystalline
                    CALL READI(12, XSTL_NN_MAX_2)       ! Number of crystalline bonds a particle needs to particles needs to be crystalline
                ENDIF

            ELSE IF (WORD == 'POLYDISPERSE') THEN        ! Consider patchy particles with polydisperse patch angles
                POLYDT = .TRUE.
                CALL READF(2, POLYSTD)                   ! Standard deviation of patch angles
                IF (NITEM > 2) THEN
                    CALL READL(3, READPOLYT)
                ELSE
                    READPOLYT = .FALSE.
                ENDIF
    
            ELSE IF (WORD == 'HALFS') THEN          ! apply harmonic springs to particles with 
                CALL READF(2, FORCEC)               ! translational coupling parameter
                CALL READF(3, FORCEO)               ! orientational coupling parameter
                IF(NITEM > 3) THEN
                    CALL READF(4,ZPIN)
                ELSE
                    ZPIN = 0.0_dp
                ENDIF
                HALFST   = .TRUE.
                SAVEREFT = .TRUE.
            
            ELSE IF (WORD == 'ORTHO') THEN ! instruct to consider a orthorhombic box
                ORTHORHOMBICT = .TRUE.
                IF (NITEM > 1) THEN
                    ORTHO_IN = .TRUE.
                    CALL READF(2, BOX(1))                !BOXA BOXB BOXC define ratio of box lengths (A:B:C x:y:z)
                    CALL READF(3, BOX(2))
                    CALL READF(4, BOX(3))
                    IF (NITEM > 4) CALL READL(5, SCALET)   
                ENDIF
    
                PRINT *, BOX
    
            ELSE IF (WORD == 'MC') THEN             ! instruct to undertake Monte Carlo simulation
                MCT = .TRUE.
    
            ELSE IF (WORD == 'MOD') THEN             
                MODT = .TRUE.
    
            ELSE IF (WORD == 'NPART') THEN            ! instructs that a molcular system is under consideration
                CALL READI(2, NPART)                  ! defines the number of molecules
                HLFPART = NPART/2
    
            ELSE IF (WORD == 'NPIN' .OR. WORD == 'NPINCNT') THEN           
        !----------------------------------------------------------------------------------------------------------------  
        !   Perform an umbrella sampling simulation, where the chosen order-parameter is the size of the 
        !   largest crystalline cluster
        !----------------------------------------------------------------------------------------------------------------
                NUCSEEDT = .TRUE. 
                CALL READI(2, TRGTSEED)               ! Target crystalline seed size
                CALL READF(3, NPINK)                  ! Harmonic potential parameter kappa
                CALL READI(4, TRJCTYLNGTH)            ! Length of each "trajectory" for a nucleus-size pinning simulation
                IF(WORD == 'NPINCNT') THEN
                    NUCCNTT = .TRUE.
                    CALL READI(5, NSEEDMAX) 
                    IF (NITEM > 5) CALL READL(6, GETSEEDT)
                    IF (NITEM > 6) CALL READL(7, PATCHBST)
                ELSE
                    NUCCNTT = .FALSE.
                    IF (NITEM > 4) CALL READL(5, GETSEEDT)
                    IF (NITEM > 5) CALL READL(6, PATCHBST)
                ENDIF

            ELSE IF (WORD == 'NUCSEED') THEN    
        !----------------------------------------------------------------------------------------------------------------
        ! Insert a crystalline nucleus into the system, requires initial coordinates to be given as well. Ideally, 
        ! the simulation should start from a fluid phase into which a crystalline seed is inserted.
        !----------------------------------------------------------------------------------------------------------------
                EQUISEEDT = .TRUE. 
                SAVEREFT  = .TRUE.
                CALL READI(2, NUCSIZE)
                CALL READF(3, FORCEC)               ! translational coupling parameter
                CALL READF(4, FORCEO)               ! orientational coupling parameter
                CALL READF(5, RHOX)

            ELSE IF (WORD == 'ONEDHIST') THEN
                ONEDHISTT = .TRUE.
                CALL READI(2, NBINS)
                CALL READF(3, BIN_MIN)
                CALL READF(4, BIN_MAX)
                CALL INITIALISE_HISTOGRAM()
            
            ELSE IF (WORD == 'ONEBOND') THEN
                ONEBONDT = .TRUE.

            ELSE IF (WORD == 'PACKING') THEN
                CALL READF(2, RHO)                   ! defines the density of the system
                PACKINGT = .TRUE.
    
            ELSE IF (WORD == 'PESRF') THEN          ! calculate PES of pair potential
                PESRFT = .TRUE.

            ELSE IF (WORD == 'POLYMOVE') THEN
                POLYMOVET = .TRUE.
                CALL READF(2, PLYMVRTIO) 
                CALL READF(3, PLYRIGID)
    
            ELSE IF (WORD == 'PRESSURE') THEN
                CALL READF(2, PRSFIX)                 ! Pressure for NPT ensemble
    
            ELSE IF (WORD == 'PRINT') THEN
                CALL READI(2, DMPFRQ)                 ! defines the frequency of dumping
                IF (NITEM > 2) THEN
                    CALL READI(3, POSDMP)             ! defines the frequency of dumping positions
                ELSE
                    POSDMP = DMPFRQ
                ENDIF

            ELSE IF (WORD == 'RADIUSGYRATION') THEN
                RADIUSGYRT = .TRUE.

            ELSE IF (WORD == 'RANDQUAT') THEN
                RANDQUATT = .TRUE.

            ELSE IF (WORD == 'RNDMLAT') THEN
                CALL READI(2, LID)
                
                IF (LID == 26) THEN              ! instructs to initialise from orthorhombic FCC or Diamond lattice 
                    LATTICET = .TRUE.
                    CALL READI(3,NUCX)
                    IF (NITEM > 4) THEN
                !   Orthorhmbic box with different number of unit cells in each direction
                        CALL READI(4,NUCY)
                        CALL READI(5,NUCZ)
                        CALL READI(6,NLAYERS)
                    ELSE
                !   Cubic box with same number of unit cells in each direction
                        NUCY = NUCX
                        NUCZ = NUCX
                        CALL READI(4,NLAYERS)
                    ENDIF
                ELSE
                    WRITE (*, '(A)') ' ERROR: LATTICE ID NOT RECOGNISED'
                    STOP                
                ENDIF

            ELSE IF (WORD == 'SCHILLING') THEN 
                CALL READI(2, SSID)
                CALL READF(3, RCUT_SS) 
                CALL READF(4, LAM_SS)        
            !   Coupling parameter for the orientational harmonic field
                IF (NITEM > 4) THEN
                    CALL READF(5, ITA_EIN) 
                ELSE
                    ITA_EIN = 1.0_dp
                ENDIF
            !   Boolean to determine whether to perform orietnational switches
                IF (NITEM > 5) THEN
                    CALL READL(6, ROT_SWITCH) 
                ELSE
                    ROT_SWITCH = .FALSE.
                ENDIF
                
                SCHSMIT  = .TRUE.
                SAVEREFT = .TRUE.
    
            ELSE IF (WORD == 'STEPS') THEN
                CALL READI(2, NSTEP)                 ! Total number of MC cycles
                CALL READI(3, NEQ)                   ! Number of equilibration steps
                CALL READI(4, BLKLNGTH)              ! Length of each block for block averaging
            
            ELSE IF (WORD == 'STEPSIZE') THEN
                CALL READF(2, MAXDTR)                   ! Translational step size
                IF (NITEM > 2) CALL READF(3, MAXDRT)    ! Rotational step size
    
            ELSE IF (WORD == 'STEPVLM') THEN
                CALL READF(2, MAXBOX)                   ! Volume step size for NPT ensemble
            
            ELSE IF (WORD == 'SPHERICAL' .OR. WORD == 'OSPHERE') THEN
                OBJCTT  = .TRUE.
                OBJSPHERET = .TRUE.
                IF(WORD == 'SPHERICAL') SPHERECNFT = .TRUE.
                CALL READF(2, SPHERERAD)
                SPHERERAD2 = SPHERERAD*SPHERERAD

            ELSE IF (WORD == 'OSURF' .OR. WORD == 'OSURFW') THEN ! Add surface to the system, surface normal is along z-axis
                OBJCTT  = .TRUE.
                OBJSURFT = .TRUE.
                IF(WORD == 'OSURFW') SURFWELLST = .TRUE. ! Add surface with potential energy wells at particular locations.

            ELSE IF (WORD == 'TAILCOR') THEN
                TAILCORT = .TRUE.

            ELSE IF (WORD == 'TEMP') THEN
                CALL READF(2, TEMP)       ! defines the temperature of the system
                BETAKB = 1.0_dp / TEMP
    
            ELSE IF (WORD == 'VERBOSE') THEN       ! instructs the simulation to output verbosely
                VERBOSET = .TRUE. 
    
            ELSE IF (WORD == 'SEED') THEN       ! read in the seed for the random number generator
                SEEDT = .TRUE. 
                CALL READI(2, ISEED)

            ELSE IF (WORD == 'STRAIN') THEN       ! perform a stress-strain simulation under uniaxial load
                STRAINT = .TRUE. 
                STRESST = .TRUE.
                CALL READI(2, STRAIN_DIM)         ! Box dimension along which strain is applied
                CALL READF(3, STRAIN_DEL)         ! Size of deformation along strained axis (percentage) 
                CALL READI(4, STRAIN_CYC)         ! Frequency (in units of cycles) to apply deformation
                CALL READF(5, MAX_STRAIN)         ! Max strain to be applied to the system
            
            ELSE IF (WORD == 'STRESS') THEN
                STRESST = .TRUE.

            ELSE IF (WORD == 'UMBRELLA') THEN           
        !----------------------------------------------------------------------------------------------------------------  
        !   Perform an umbrella sampling simulation, where the chosen order-parameter is the global 
        !   Steinhardt Order Parameter and/or the density
        !----------------------------------------------------------------------------------------------------------------
                UMBRELLAT = .TRUE. 
                CALL READF(2, TRGTQ)            ! Target Ql value
                CALL READF(3, QPINK)            ! Harmonic potential parameter kappa for Ql
                IF(NITEM > 4) THEN
                    TWODUT = .TRUE.
                    CALL READF(4, TRGTR)        ! Target rho value
                    CALL READF(5, RPINK)        ! Harmonic potential parameter kappa for rho
                    CALL READI(6, TRJCTYLNGTH)  ! Length of each "trajectory" for a nucleus-size pinning simulation
                ELSE
                    TWODUT = .FALSE.
                    CALL READI(4, TRJCTYLNGTH)  ! Length of each "trajectory" for a nucleus-size pinning simulation
                ENDIF

            ELSE IF (WORD == 'UNITVEC') THEN       ! read in the orientations as unit vectors
                UNITVECT = .TRUE. 
    
    !     ======================= the following are system definitions =================================
    
            ELSE IF (WORD == 'HS') THEN
                HST = .TRUE.
                HARDT = .TRUE.

            ELSE IF (WORD == 'BINARYHS') THEN
                BHST    = .TRUE.
                HARDT   = .TRUE.
                BINARYT = .TRUE.
                CALL READF(2, BNRYR)
                CALL READF(3, SIGAA)
                CALL READF(4, SIGBB)
                RCUT = 1.0_dp

            ELSE IF (WORD == 'SQUARE') THEN
                SQSHT  = .TRUE.
                HARDT  = .TRUE.
                CALL READF(2, SQSHDEL)
                CALL READF(3, SQSHEPS)
                CALL READI(4, SQID)

            ELSE IF (WORD == 'YUKAWA') THEN
                YUKT = .TRUE.
                CALL READF(2, YUKKAP)
            
            ELSE IF (WORD == 'GLJ') THEN
                GLJT = .TRUE.
                CALL READF(2, GLJN)
                IF (NITEM > 2) CALL READL(3, WCAT)

                IF(WCAT) THEN
                    RCUT   = 2.0_dp**(1.0_dp/GLJN)
                    RCUTSQ = RCUT*RCUT
                ENDIF

            ELSE IF (WORD == 'BINARYGLJ') THEN
                BGLJT = .TRUE.
                CALL READF(2, GLJN)
                CALL READF(3, BNRYR)
                CALL READF(4, SIGAA)
                CALL READF(5, SIGBB)
                CALL READF(6, EPSAA)
                CALL READF(7, EPSBB)
            
            ELSE IF (WORD == 'KF') THEN
                KFT    = .TRUE.
                HARDT  = .TRUE.
                RIGIDT = .TRUE.
                CALL READI(2, NSITES)
    
                IF(NSITES == 1) THEN
                    CALL READF(3, KFDELA)
                    CALL READF(4, KFLAMA)
                ELSEIF(NSITES == 2) THEN
                    CALL READF(3, KFAA)
                    CALL READF(4, KFBB)
                    CALL READF(5, KFDELA)
                    CALL READF(6, KFDELB)
                    CALL READF(7, KFLAMA)
                    CALL READF(8, KFLAMB)
                ELSEIF(NSITES == 3) THEN
                    CALL READF(3,  KFAA)
                    CALL READF(4,  KFBB)
                    CALL READF(5,  KFCC)
                    CALL READF(6,  KFDELA)
                    CALL READF(7,  KFDELB)
                    CALL READF(8,  KFDELC)
                    CALL READF(9,  KFLAMA)
                    CALL READF(10, KFLAMB)
                    CALL READF(11, KFLAMC)
                ELSEIF(NSITES == 4) THEN
                    CALL READF(3,  KFAA)
                    CALL READF(4,  KFBB)
                    CALL READF(5,  KFCC)
                    CALL READF(6,  KFDD)
                    CALL READF(7,  KFDELA)
                    CALL READF(8,  KFDELB)
                    CALL READF(9,  KFDELC)
                    CALL READF(10, KFDELD)
                    CALL READF(11, KFLAMA)
                    CALL READF(12, KFLAMB)
                    CALL READF(13, KFLAMC)
                    CALL READF(14, KFLAMD)
                ELSEIF(NSITES == 6) THEN
                    CALL READF(3,  KFAA)
                    CALL READF(4,  KFBB)
                    CALL READF(5,  KFCC)
                    CALL READF(6,  KFDD)
                    CALL READF(7,  KFEE)
                    CALL READF(8,  KFFF)
                    CALL READF(9,  KFDELA)
                    CALL READF(10, KFDELB)
                    CALL READF(11, KFDELC)
                    CALL READF(12, KFDELD)
                    CALL READF(13, KFDELE)
                    CALL READF(14, KFDELF)
                    CALL READF(15, KFLAMA)
                    CALL READF(16, KFLAMB)
                    CALL READF(17, KFLAMC)
                    CALL READF(18, KFLAMD)
                    CALL READF(19, KFLAME)
                    CALL READF(20, KFLAMF)
                ELSE
                    PRINT *, "KERN-FRENKEL MODEL CURRENTLY ONLY CONFIGURED FOR PARTICLES WITH 2, 3, 4, OR 6 PATCHES"
                    STOP
                ENDIF
                
                RCUT   = MAX(KFLAMA, KFLAMB, KFLAMC, KFLAMD)
                RCUTSQ = RCUT**2
    
                IF(NSITES == 2) THEN
                    IF(KFAA == KFBB .AND. KFDELA == KFDELB .AND. KFLAMA == KFLAMB) HEADTAILT = .TRUE.
                ENDIF
    
            ELSE IF (WORD == 'PGLJ') THEN
                PGLJT  = .TRUE.
                RIGIDT = .TRUE.
                CALL READF(2, GLJN)
                CALL READF(3, SIGPW)
                CALL READI(4, NSITES)
    
            ELSE IF (WORD == 'CPP') THEN
                CPPT  = .TRUE.
                RIGIDT = .TRUE.
                CALL READI(2, NSITES)
                CALL READF(3, CPPAA)
                IF(NSITES == 2) THEN
                    CALL READF(4, CPPBB)
                    CALL READF(5, CPPDELA)
                    CALL READF(6, CPPDELB)
                    CALL READF(7, CPPLAM)
                    CALL READF(8, CPPS)
                    CALL READF(9, YUKKAP)
                ELSEIF(NSITES == 3) THEN
                    CALL READF(4, CPPBB)
                    CALL READF(5, CPPCC)
                    CALL READF(6, CPPDELA)
                    CALL READF(7, CPPDELB)
                    CALL READF(8, CPPDELC)
                    CALL READF(9, CPPLAM)
                    CALL READF(10, CPPS)
                    CALL READF(11,YUKKAP)
                ELSEIF(NSITES == 4) THEN
                    CALL READF(4, CPPBB)
                    CALL READF(5, CPPCC)
                    CALL READF(6, CPPDD)
                    CALL READF(7, CPPDELA)
                    CALL READF(8, CPPDELB)
                    CALL READF(9, CPPDELC)
                    CALL READF(10, CPPDELD)
                    CALL READF(11,CPPLAM)
                    CALL READF(12,CPPS)
                    CALL READF(13,YUKKAP)
                ELSE
                    STOP "CPP MODEL CURRENTLY ONLY CONFIGURED FOR PARTICLES WITH 2, 3 OR 4 PATCHES"
                ENDIF
    
                IF(NSITES == 2) THEN
                    IF(CPPBB == 1.0_dp .AND. CPPDELA == CPPDELB) HEADTAILT = .TRUE.
                ENDIF

            ELSE IF (WORD == 'BNRYYUKWA') THEN           ! consider a binary system where half the particles do not interact with one another.
                BYUKWAT = .TRUE.
                CALL READI(2, BNRYA)              ! Relative ratio between A and B particles
                CALL READF(3, SIGAA)
                CALL READF(4, SIGBB)
                IF(NITEM>4) CALL READF(5,RHOX)
                IF(NITEM>5) CALL READL(6,BNRYEQUIT)
                CALL DEF_BINARY_SPHERES()
            
            ELSE IF (WORD == 'KIHARA') THEN
                KIHARAT  = .TRUE.
                RIGIDT = .TRUE.
                CALL READF(2, GLJN)
                CALL READF(3, RLNGTH)
    
                HEADTAILT = .TRUE.
    
            ELSE IF (WORD == 'ETP') THEN
                ETPT  = .TRUE.
                RIGIDT = .TRUE.
                CALL READI(2, NSITES)
                CALL READF(3, GLJN)
                CALL READF(4, RLNGTH)
                CALL READF(5, CPPLAM)
                CALL READF(6, CPPS)
                CALL READF(7, CPPDELA)   
                IF(NSITES == 2) THEN
                    CALL READF(8, CPPDELB)
                    CALL READF(9, CPPBB)
                ELSEIF(NSITES > 2) THEN
                    STOP "ETP MODEL CURRENTLY ONLY CONFIGURED FOR PARTICLES WITH 1 OR 2 PATCHES"
                ENDIF
    
                IF(NSITES == 2) THEN
                    IF(CPPBB == 1.0_dp .AND. CPPDELA == CPPDELB) HEADTAILT = .TRUE.
                ENDIF

            ELSE IF (WORD == 'HTPR') THEN
                HTPRT  = .TRUE.
                HARDT  = .TRUE.
                RIGIDT = .TRUE.
                CALL READI(2, NSITES)
                CALL READF(3, RLNGTH)
                CALL READF(4, CPPLAM)
                CALL READF(5, CPPDELA)   
                IF(NSITES == 2) THEN
                    CALL READF(6, CPPDELB)
                    CALL READF(7, CPPBB)
                ELSEIF(NSITES > 2) THEN
                    STOP "ETP MODEL CURRENTLY ONLY CONFIGURED FOR PARTICLES WITH 1 OR 2 PATCHES"
                ENDIF

                CPPLAM2 = CPPLAM*CPPLAM
                
                IF(NSITES == 2) THEN
                    IF(CPPBB == 1.0_dp .AND. CPPDELA == CPPDELB) HEADTAILT = .TRUE.
                ENDIF
    
            ELSE IF (WORD == 'DMBLGLJ') THEN
                DMBLGLJT  = .TRUE.
                RIGIDT    = .TRUE.
                CALL READI(2, NSITES)
                CALL READF(3, GLJN)
                CALL READF(4, SIGBB)
                CALL READF(5, EPSBB)
                IF(NSITES == 3) THEN
                    CALL READF(6, SIGCC)
                    CALL READF(7, EPSCC)
                ELSEIF(NSITES /= 2) THEN
                    STOP "NSITES MUST BE 2 OR 3"
                ENDIF
    
                IF(NSITES == 2) THEN
                    IF(SIGBB == 1.0_dp .AND. EPSBB == 1.0_dp) HEADTAILT = .TRUE.
                ENDIF

            ELSE IF (WORD == 'HDMBL') THEN
                HDMBLT    = .TRUE.
                HARDT     = .TRUE.
                RIGIDT    = .TRUE.
                HEADTAILT = .TRUE.
                NSITES    = 2
                CALL READF(2, LSTAR)
    
            ELSE IF (WORD == 'KFREC') THEN
                KFRECT    = .TRUE.
                HARDT     = .TRUE.
                RIGIDT    = .TRUE.
                NSITES    = 6
    
                CALL READF(2, KFAA)
                CALL READF(3, KFBB)
                CALL READF(4, KFDELA1)
                CALL READF(5, KFDELA2)
                CALL READF(6, KFDELB1)
                CALL READF(7, KFDELB2)
                CALL READF(8, SKEWAB)
                CALL READF(9, KF_LAM)
                IF (NITEM > 9) THEN
                    RACEMICT = .TRUE.
                    CALL READF(10, SKEWAB2)
                ELSE
                    RACEMICT = .FALSE.
                ENDIF
                
                KF_LAM2 = KF_LAM*KF_LAM
                RCUT    = KF_LAM
                RCUTSQ  = KF_LAM2

                IF(KFAA == KFBB .AND. SKEWAB==0.0_dp) HEADTAILT = .TRUE.

            ELSE IF (WORD == 'NEMATIC_CLD') THEN
                NMTCCT    = .TRUE.
                CALL READI(2, LMIN)
                CALL READI(3, LMAX)
                CALL READI(4, LSTEP)
                CALL READF(5, K_NC)
                CALL READF(6, KAPPA)
                CALL READF(7, RCUT_YUKAWA)
                IF (NITEM > 7) THEN
                    IL = 0
                    DO I = LMIN,LMAX,LSTEP
                        IL = IL + 1
                    ENDDO
                    ALLOCATE(NC_QLM(IL))
                    IL = 7
                    DO I = LMIN,LMAX,LSTEP
                        IL = IL + 1
                        CALL READF(IL, NC_QLM(IL-7))
                    ENDDO
                ENDIF

            ELSE IF (WORD == 'OBLATESPHY') THEN
                RIGIDT      = .TRUE.
                OBLATESPHYT = .TRUE.
                CALL READF(2, DSIG)
                CALL READF(3, OBLL)
                IF (NITEM > 3) THEN
                    NSITES = 6
                    PATCHYOBLT = .TRUE.
                    CALL READF(4, POBLTHETA)
                    CALL READF(5,  KFAA)
                    CALL READF(6,  KFBB)
                    CALL READF(7,  KFCC)
                    CALL READF(8,  KFDD)
                    CALL READF(9,  KFEE)
                    CALL READF(10, KFDELA)
                    CALL READF(11, KFDELB)
                    CALL READF(12, KFDELC)
                    CALL READF(13, KFDELD)
                    CALL READF(14, KFDELE)
                    CALL READF(15, KFLAMA)
                    CALL READF(16, KFLAMB)
                    CALL READF(17, KFLAMC)
                    CALL READF(18, KFLAMD)
                    CALL READF(19, KFLAME)
                ELSE
                    NSITES = 2
                ENDIF
                HLF_DSIG = DSIG / 2.0_dp

            ELSE IF (WORD == 'POLY_CHAIN' .OR. WORD == 'ATTR_POLY_CHAIN') THEN
                POLYCHAINT = .TRUE.
                CALL READI(2, N_POLY)
                CALL READI(3, N_POLY_L)
                CALL READF(4, POLYC_L)
                CALL READF(5, POLY_KAP)
                CALL READF(6, POLY_SIG)
                IF (NITEM > 6) THEN
                    CALL READF(7, POLY_EPS)
                    CALL READF(8, POLY_DEL)
                    IF(POLY_EPS==0.0_dp .OR. POLY_DEL==0.0_dp) THEN
                        ATTR_POLY = .FALSE.
                    ELSE
                        ATTR_POLY = .TRUE.
                        POLYCUT   = POLY_DEL * POLY_SIG
                        POLYCUT2  = POLYCUT*POLYCUT
                    ENDIF
                ENDIF
                IF (NITEM > 8) THEN
                    CALL READF(9, CP_EPS)
                    CALL READF(10, CP_DEL)
                    IF(NITEM > 10) THEN
                        CALL READF(11, PATCH_RAD)
                    ELSE
                        PATCH_RAD = 0.0_dp
                    ENDIF
                ENDIF

                POLY_SIG2  = POLY_SIG*POLY_SIG
                N_POLY_TOT = N_POLY*N_POLY_L

                IF(WORD == 'ATTR_POLY_CHAIN') THEN
                    APCHAINT = .TRUE.
                    ALLOCATE(PCHAIN_EPS(0:N_POLY_L-1,0:N_POLY_L-1))
                    PCHAIN_EPS = 0.0_dp
                    OPEN (UNIT = 24, FILE ='pchain', STATUS = 'UNKNOWN')
                    DO
                        READ(24,*,IOSTAT=EOF) J1, J2, CHAINE
                        IF (EOF < 0) GOTO 1776
                        PCHAIN_EPS(J1-1,J2-1) = CHAINE
                        PCHAIN_EPS(J2-1,J1-1) = CHAINE
                    END DO
                    1776 CONTINUE   
                    POLYCUT   = 1.1_dp * POLY_SIG
                    POLYCUT2  = POLYCUT*POLYCUT
                ENDIF

            ELSE IF (WORD=='TRIANGLE'.OR.WORD=='PATCHYTRI'.OR.WORD=='ATRIANGLE'.OR.WORD=='APATCHYTRI') THEN
                TRIANGLET  = .TRUE.
                RIGIDT     = .TRUE.
                NSITES = 4
            !   Define geometry of the triangle 
                CALL READF(2, TRI_SIG)
                CALL READF(3, TRI_ALPHA)
                CALL READF(4, TRI_B)
                CALL READF(5, TRI_C)
            !   If simulating s1quare-well hard-core triangles, set range and well-depth of square-well 
                IF(WORD == 'ATRIANGLE' .OR. WORD == 'APATCHYTRI') THEN
                    SQWT = .TRUE.
                    CALL READF(6, SQWDEL)
                    CALL READF(7, SQWEPS)
                    CSHFT = 2
                ELSE
                    SQWT = .FALSE.
                    CSHFT = 0
                ENDIF
            !   If simulating patchy triangles define geometry and energetics of the patches
                IF(WORD == 'PATCHYTRI' .OR. WORD == 'APATCHYTRI') THEN
                    PATCHYTRIT = .TRUE.
                    NSITES = 10
                    IF (NITEM < 18+CSHFT) THEN
                        CALL READF(6+CSHFT,  KFAA)
                        CALL READF(7+CSHFT,  KFBB)
                        CALL READF(8+CSHFT,  KFCC)
                        CALL READF(9+CSHFT,  KFDELA)
                        CALL READF(10+CSHFT, KFDELB)
                        CALL READF(11+CSHFT, KFDELC)
                        CALL READF(12+CSHFT, KFLAMA)
                        CALL READF(13+CSHFT, KFLAMB)
                        CALL READF(14+CSHFT, KFLAMC)
                        CALL READF(15+CSHFT, TRI_THETAA)
                        IF(NITEM < 17) THEN
                            TRI_THETAB = TRI_THETAA
                            TRI_THETAC = TRI_THETAA
                        ELSE
                            CALL READF(16+CSHFT, TRI_THETAB)
                            CALL READF(17+CSHFT, TRI_THETAC)
                        ENDIF
                        TRI_PHIA = KFDELA
                        TRI_PHIB = KFDELB
                        TRI_PHIC = KFDELC
                    ELSE
                        CALL READF(6+CSHFT,  KFAA)
                        CALL READF(7+CSHFT,  KFBB)
                        CALL READF(8+CSHFT,  KFCC)
                        CALL READF(9+CSHFT,  KFDELA)   ! "Height" of the patch
                        CALL READF(10+CSHFT, TRI_PHIA) ! "Width" of the patch
                        CALL READF(11+CSHFT, KFDELB)
                        CALL READF(12+CSHFT, TRI_PHIB)
                        CALL READF(13+CSHFT, KFDELC)
                        CALL READF(14+CSHFT, TRI_PHIC)
                        CALL READF(15+CSHFT, KFLAMA)
                        CALL READF(16+CSHFT, KFLAMB)
                        CALL READF(17+CSHFT, KFLAMC)
                        CALL READF(18+CSHFT, TRI_THETAA)
                        IF(NITEM < 20+CSHFT) THEN
                            TRI_THETAB = TRI_THETAA
                            TRI_THETAC = TRI_THETAA
                        ELSE
                            CALL READF(19+CSHFT, TRI_THETAB)
                            CALL READF(20+CSHFT, TRI_THETAC)
                        ENDIF
                    ENDIF
                ENDIF

            ELSE IF (WORD == 'SPECIFICKF') THEN
                SPECIFICKFT = .TRUE.
                CALL READF(2,  KFAA)
                CALL READF(3,  KFAB)
                CALL READF(4,  KFBB)
                IF (NITEM > 4) THEN
                    CALL READF(5,  KFAC)
                    CALL READF(6,  KFBC)
                    CALL READF(7,  KFCC)
                ENDIF
                IF (NITEM > 7) THEN
                    CALL READF(8,  KFAD)
                    CALL READF(9,  KFBD)
                    CALL READF(10, KFCD)
                    CALL READF(11, KFDD)
                ENDIF
                IF (NITEM > 11) THEN
                    CALL READF(11, KFAE)
                    CALL READF(12, KFBE)
                    CALL READF(13, KFCE)
                    CALL READF(14, KFDE)
                    CALL READF(15, KFEE)
                ENDIF
                IF (NITEM > 15) THEN
                    CALL READF(16, KFAF)
                    CALL READF(17, KFBF)
                    CALL READF(18, KFCF)
                    CALL READF(19, KFDF)
                    CALL READF(20, KFEF)
                    CALL READF(21, KFFF)
                ENDIF
            ENDIF
    
        END DO
    
    900 ENDT = .TRUE.
    
        CLOSE (INPUNIT)
    
    END SUBROUTINE KEYWORDS
    
