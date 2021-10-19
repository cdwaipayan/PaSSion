!     PASSION: A package for Soft Matter Simulation
!     Copyright (C) 2015 - Dwaipayan Chakrabarti
!     This file is part of PASSION.
!
SUBROUTINE KEYWORDS
    !     This subroutine reads in the input parameters, which are driven by keywords.
        USE COMMONS
    
        IMPLICIT NONE
    
        INTEGER                 :: I
        CHARACTER (LEN = 16)    :: WORD
        LOGICAL                 :: ENDT
    
        ENDT  = .FALSE.
    
        OPEN (UNIT = INPUNIT, FILE = 'input.inp', STATUS = 'UNKNOWN')
    
    !     The END argument of the READ statement has a label number after it to which the program will
    !     branch, if it encounters the end of the file while reading.
    
    200 DO WHILE (.NOT. ENDT)
    
            READ (INPUNIT, FMT = 301, END = 900) (CHAR(I), I = 1, 80)
    301     FORMAT (80A1)
    
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
                    BNRYRTIOT = .TRUE.
                    CALL READF(2, BNRYR)              ! Relative ratio between A and B particles
                    BNRYA = INT(BNRYR*REAL(NPART,DP))
                ELSE
                    BNRYRTIOT = .FALSE.
                ENDIF
    
            ELSE IF (WORD == 'BOP' .OR. WORD == 'CBOP') THEN          ! instruct to consider a cubic box
                BOPT = .TRUE.
                CALL READI(2, MIN_L)
                CALL READI(3, MAX_L)
                CALL READI(4, DEL_L)
                CALL READI(5, QLT)
                CALL READI(6, QBARLT)
                CALL READI(7, GQLT)
                CALL READI(8, QLDOTQLT)
                CALL READF(9, BOP_CUT)
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
    
            ELSE IF (WORD == 'CUBIC') THEN          ! instruct to consider a cubic box
                CUBICT = .TRUE.
    
            ELSE IF (WORD == 'CUTOFF') THEN         ! specifies the cutoff distance for the potential
                CALL READF(2, RCUT)
                RCUTSQ = RCUT*RCUT
    
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

            ELSE IF (WORD == 'LATTICE') THEN
                CALL READI(2, LID)
                
                IF (LID < 30) THEN              ! instructs to initialise from orthorhombic FCC or Diamond lattice 
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
                    CALL READI(3, READPOLYT)
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
                    IF (NITEM > 4) CALL READF(5, SCALET)   
                ENDIF
    
                PRINT *, BOX
    
            ELSE IF (WORD == 'MC') THEN             ! instruct to undertake Monte Carlo simulation
                MCT = .TRUE.
    
            ELSE IF (WORD == 'MOD') THEN             
                MODT = .TRUE.
    
            ELSE IF (WORD == 'NPART') THEN            ! instructs that a molcular system is under consideration
                CALL READI(2, NPART)                  ! defines the number of molecules
                HLFPART = NPART/2
    
            ELSE IF (WORD == 'NPIN') THEN           
        !----------------------------------------------------------------------------------------------------------------  
        !   Perform an umbrella sampling simulation, where the chosen order-parameter is the size of the 
        !   largest crystalline cluster
        !----------------------------------------------------------------------------------------------------------------
                NUCSEEDT = .TRUE. 
                CALL READI(2, TRGTSEED)               ! Target crystalline seed size
                CALL READF(3, NPINK)                  ! Harmonic potential parameter kappa
                CALL READI(4, TRJCTYLNGTH)            ! Length of each "trajectory" for a nucleus-size pinning simulation
                IF (NITEM > 4) CALL READI(4, GETSEEDT)

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
    
            ELSE IF (WORD == 'PACKING') THEN
                CALL READF(2, RHO)                   ! defines the density of the system
                PACKINGT = .TRUE.
    
            ELSE IF (WORD == 'PESRF') THEN          ! calculate PES of pair potential
                PESRFT = .TRUE.
    
            ELSE IF (WORD == 'PRESSURE') THEN
                CALL READF(2, PRSFIX)                 ! Pressure for NPT ensemble
    
            ELSE IF (WORD == 'PRINT') THEN
                CALL READI(2, DMPFRQ)                 ! defines the frequency of dumping
                IF (NITEM > 2) THEN
                    CALL READI(3, POSDMP)             ! defines the frequency of dumping positions
                ELSE
                    POSDMP = DMPFRQ
                ENDIF

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
                    CALL READI(6, ROT_SWITCH) 
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

            ELSE IF (WORD == 'SPHERICAL') THEN
                SPHERECNFT = .TRUE.
                CALL READF(2, SPHERERAD)
                SPHERERAD2 = SPHERERAD*SPHERERAD

            ELSE IF (WORD == '2DSURF') THEN
                TDSRFT = .TRUE.
                IF (NITEM > 1) CALL READI(2, GEN2D)               
    
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
    
            ELSE IF (WORD == 'UNITVEC') THEN       ! read in the orientations as unit vectors
                UNITVECT = .TRUE. 
    
    !     ======================= the following are system definitions =================================
    
            ELSE IF (WORD == 'HS') THEN
                HST = .TRUE.
                HARDT = .TRUE.

            ELSE IF (WORD == 'YUKAWA') THEN
                YUKT = .TRUE.
                CALL READF(2, YUKKAP)
            
            ELSE IF (WORD == 'GLJ') THEN
                GLJT = .TRUE.
                CALL READF(2, GLJN)
                IF (NITEM > 2) CALL READI(3, WCAT)

                IF(WCAT) THEN
                    RCUT   = 2.0_dp**(1.0_dp/GLJN)
                    RCUTSQ = RCUT*RCUT
                ENDIF
            
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
                ELSE
                    PRINT *, "KERN-FRENKEL MODEL CURRENTLY ONLY CONFIGURED FOR PARTICLES WITH 2, 3 OR 4 PATCHES"
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
                IF(NSITES == 2) THEN
                    CALL READF(3, CPPBB)
                    CALL READF(4, CPPDELA)
                    CALL READF(5, CPPDELB)
                    CALL READF(6, CPPLAM)
                    CALL READF(7, CPPS)
                    CALL READF(8, YUKKAP)
                ELSEIF(NSITES == 3) THEN
                    CALL READF(3, CPPBB)
                    CALL READF(4, CPPCC)
                    CALL READF(5, CPPDELA)
                    CALL READF(6, CPPDELB)
                    CALL READF(7, CPPDELC)
                    CALL READF(8, CPPLAM)
                    CALL READF(9, CPPS)
                    CALL READF(10,YUKKAP)
                ELSEIF(NSITES == 4) THEN
                    CALL READF(3, CPPBB)
                    CALL READF(4, CPPCC)
                    CALL READF(5, CPPDD)
                    CALL READF(6, CPPDELA)
                    CALL READF(7, CPPDELB)
                    CALL READF(8, CPPDELC)
                    CALL READF(9, CPPDELD)
                    CALL READF(10,CPPLAM)
                    CALL READF(11,CPPS)
                    CALL READF(12,YUKKAP)
                ELSE
                    STOP "CPP MODEL CURRENTLY ONLY CONFIGURED FOR PARTICLES WITH 2, 3 OR 4 PATCHES"
                ENDIF
    
                IF(NSITES == 2) THEN
                    IF(CPPBB == 1.0_dp .AND. CPPDELA == CPPDELB) HEADTAILT = .TRUE.
                ENDIF
            
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

            ELSE IF (WORD == 'GAYBERNE') THEN
                GBT       = .TRUE.
                RIGIDT    = .TRUE.
                HEADTAILT = .TRUE.
                NSITES    = 1
                CALL READF(2, GBK)
                CALL READF(3, GBKP)
                CALL READF(4, GBV)
                CALL READF(5, GBM)
                CALL READF(6, GBCUT)
                GBX  = (GBK**2 - 1.0_dp) / (GBK**2 + 1.0_dp)
                GBXP = (GBKP**(1.0_dp/GBM) - 1.0_dp) / (GBKP**(1.0_dp/GBM) + 1.0_dp)
                IF (NITEM > 6) THEN
                    EWALDT = .TRUE.
                    CALL READF(7, DPMU)
                    DPMUSQ = DPMU*DPMU
                ENDIF
    
            ENDIF
    
        END DO
    
    900 ENDT = .TRUE.
    
        CLOSE (INPUNIT)
    
    END SUBROUTINE KEYWORDS
    