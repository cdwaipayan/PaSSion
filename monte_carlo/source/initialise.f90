!     This file contains a number of subroutines that are called once at the beginning of the
!     simulation to initialise it.

SUBROUTINE INITIALISE()

    USE COMMONS
    USE ROTATIONS_MODULE, ONLY: Q_TO_RM
    USE CELL_LIST, ONLY: INITIALIZE_LIST, MAKE_LIST
    USE CLUSTER_MOVE, ONLY: BUILD_ALL_CLUSTERS
    
    IMPLICIT NONE

    INTEGER         :: J1, J2, RC1
    REAL(KIND=DP)   :: CORP, RM(3,3), N

    INTEGER              :: SRAND48, XXX, SEED_N, SEED_J
    INTEGER, ALLOCATABLE :: SEED_ARRAY(:)
    DOUBLE PRECISION     :: SEED_U

!   *********************************************************************
    IF(.NOT. SEEDT) THEN
    ! Initialise random number generator using drand48 and seeding the 
    ! random number generator from from /dev/urandom.
    ! In Unix-like operating systems, /dev/random is a special file that 
    ! serves as pseudorandom number generators.
    ! If you read 4 bytes from this file then you will get a nicely random 
    ! and unpredictable 32-bit seed value, and you can read as many bytes 
    ! as you need. 
    ! The file collects system hardware timings (e.g. hard disk access 
    ! times and intervals between key presses) to produce entropy and is 
    ! as close as you will get to a source of true random numbers without 
    ! a radioactive source and a Geiger counter!
    !*********************************************************************
    !---------------------------------------------------------------------
    ! Generate an array of random 32-bit integers, where the size of the
    ! array is SEED_N which is defined by calling RANDOM_SEED
    !---------------------------------------------------------------------
        CALL RANDOM_SEED(SIZE=SEED_N)
        ALLOCATE(SEED_ARRAY(SEED_N))    
        OPEN(69, FILE='/dev/urandom', ACCESS='STREAM', FORM='UNFORMATTED')
            READ(69) SEED_ARRAY
        CLOSE(69)
    !---------------------------------------------------------------------
    ! Extract a random integer between 1 and SEED_N, in order to randomly
    ! select one of the random 32-bit integers generated in the previous
    ! step.
    !---------------------------------------------------------------------
        CALL RANDOM_NUMBER(SEED_U)
        SEED_J = FLOOR(SEED_N*SEED_U) + 1
        ISEED = SEED_ARRAY(SEED_J)
    !---------------------------------------------------------------------
    ! Use ISEED to seed the DRAND48 random number generator, making sure
    ! to write the seed as an output so that the data is reproducible.
    !---------------------------------------------------------------------
        IF (ISEED.LT.0) ISEED=-ISEED
        OPEN (79, FILE='seed.dat',STATUS='UNKNOWN')
            WRITE(79,*) ISEED
        CLOSE(79)
    ENDIF
!   Seed the random number generator
    XXX=SRAND48(ISEED)
!   *************************************************************************

!---------------------------------------------------------------------
!   Allocate arrays to be used during the simulation
!---------------------------------------------------------------------
    ALLOCATE( R(NDIM,NPART) )
    IF(RIGIDT) THEN
        ALLOCATE( Q(4,NPART), REFSITE(NDIM,NSITES), RBSITES(NDIM,NSITES,NPART) )
        IF(RACEMICT) ALLOCATE( REFSITE2(NDIM,NSITES) )
    ENDIF

!---------------------------------------------------------------------
!   Initialise parameters of the selected system to be simulated
!---------------------------------------------------------------------
    CALL INIT_POT()

!---------------------------------------------------------------------
!   Calculate the pair potential energy surface for the selected 
!   potential energy function.
!---------------------------------------------------------------------
    IF(PESRFT) CALL PESRF()

!   *************************************************************************
!   Set initial translational and rotational coordiates 
    IF (LATTICET) THEN
    !   Start from a predifined crystal structure
        CALL GENERATE_UNIT_CELLS()
!   ------------------------------------------------------------------------
    ELSE
    !   Read configuration from a file provided
        CALL READCONFIG()
    ENDIF
!   *************************************************************************

!   *************************************************************************
    IF (RIGIDT) THEN
    !   Check that a rotational stepsize has been given, if not set it to be 
    !   equal to the translational stepsize. 
        IF(MAXDRT == 0.0_dp) MAXDRT = MAXDTR
    !   Initialise rigid body sites
        DO J1 = 1, NPART
            RM = Q_TO_RM( Q(:,J1) )
            IF(RACEMICT) RC1 = MOD(((J1-1)-MOD((J1-1),12))/12+1,2)
            DO J2 = 1, NSITES
                IF(RACEMICT) THEN
                    IF(RC1==1) THEN
                        RBSITES(:,J2,J1) = MATMUL(RM,REFSITE(:,J2))
                    ELSE
                        RBSITES(:,J2,J1) = MATMUL(RM,REFSITE2(:,J2))
                    ENDIF
                ELSE
                    RBSITES(:,J2,J1) = MATMUL(RM,REFSITE(:,J2))
                ENDIF
            ENDDO
        ENDDO
    ENDIF

!   *************************************************************************
!   Save reference configuration for system with particles in a harmonic trap
!   (i.e., Einstein crystal, interface pinning, Schilling-Schimd, ...)
!   *************************************************************************
    IF(SAVEREFT) THEN
        ALLOCATE(RREF(NDIM,NPART))
        RREF = R
        IF(RIGIDT) THEN
            ALLOCATE(RBREF(NDIM,NSITES,NPART))
            DO J1 = 1, NPART
                DO J2 = 1, NSITES
                    RBREF(:,J2,J1) = RBSITES(:,J2,J1)
                ENDDO
            ENDDO
        ENDIF
    ENDIF

!   *************************************************************************
!    Initialise parameters related to the use of a cell list.
!   *************************************************************************
    IF(CELLLISTT) THEN
        CALL INITIALIZE_LIST( NPART, RCUT/BOX )
        ALLOCATE(SCALED_R(3,NPART))
        SCALED_R = 0.0_dp
        
        DO J1 = 1, NPART
            SCALED_R(:,J1) = R(:,J1)/BOX
        ENDDO
        
        CALL MAKE_LIST( NPART, SCALED_R )
    ENDIF

!   *************************************************************************
!    Initialise parameters related to cluster move Monte Carlo.
!   *************************************************************************
    IF(CLUSTERMOVET .OR. BOPCLSTRT) THEN
        ALLOCATE(CLSTR(NPART),RCLSTR(NDIM,NPART))
        CLSTR   = 0
        CLSTRSZ = 0
        MAXDRTC = MAXDRT
        MAXDTRC = MAXDTR
        IF(NPTT .OR. BOPCLSTRT) THEN
            ALLOCATE(CLSTRADJ(NPART,NPART), CLSTRCOM(NDIM,NPART), PCLSTR(NDIM,NPART), PCLSTRID(NPART), CLURIJ(NDIM,NPART,NPART))
            NCLSTRS  = 0      ! Number of unique clusters in the system
            CLSTRADJ = 0      ! Adjacency matrix for the system (bonding criterion is system specific)
            PCLSTRID = 0      ! Cluster ID associated with each particle.
            CLSTRCOM = 0.0_dp ! Center-of-mass for each cluster
            PCLSTR   = 0.0_dp ! Relative position of each particle to the center-of-mass of its associated cluster.
        ENDIF
    ENDIF

!   *************************************************************************
!   Calculate the initial potential energy and virial of the system
!   *************************************************************************
    IF(BOPCLSTRT) VLMCLUSTERMOVET = .TRUE.
    CALL POTENTIAL(PE)
    IF(BOPCLSTRT) VLMCLUSTERMOVET = .FALSE.
    IF(OVERLAPT) STOP "INITIAL CONFIGURATION HAS OVERLAPPING PARTICLES"
    VIR = VIRTEMP

!   *************************************************************************
!    Initialise parameters related to the calculation of bond-orientational
!    order parameters and the largest crystalline cluster.
!   *************************************************************************
    IF(BOPT) THEN

        BOP_CUT_SQ = BOP_CUT**2
        IF(MIN_L == MAX_L) THEN
            J2 = 1
            IF(LRGSTXCLSTRT) L_INDEX = J2
        ELSE
            J2 = 0
            DO J1 = MIN_L, MAX_L, DEL_L
                J2 = J2 + 1
                IF(LRGSTXCLSTRT .AND. J1 == BOP_L) L_INDEX = J2
            ENDDO
        ENDIF

        IF(BOPCLSTRT) THEN
            CALL BUILD_ALL_CLUSTERS(BOPSNT=.TRUE.)
            NPBOP = NCLSTRS
        ELSE
            NPBOP = NPART
        ENDIF

        IF(QLT)          ALLOCATE(QL(J2,NPBOP))
        IF(QBARLT)       ALLOCATE(QBARL(J2,NPBOP))
        IF(QLDOTQLT)     ALLOCATE(QLDOTQL(J2,NPBOP,NPBOP))
        IF(GQLT)         ALLOCATE(GQL(J2))
        IF(LRGSTXCLSTRT) ALLOCATE(QLDIST(NPBOP,NPBOP))

    ENDIF

!   *************************************************************************
!   Set up the variables to calculate the free energy of crystals using the 
!   the Frenkel-Ladd (Einstein Crystal) method.
!   *************************************************************************
    IF(FLFET) THEN
        LAMOR = (LAMTR * ITA_EIN)
        print *, LAMTR, LAMOR
        U_LATT = PE
        ALLOCATE(DRCOM(NDIM)); DRCOM = 0.0_dp ! shift in the center of mass
        IF(FLID == 0) THEN
        !   Calculate the free energy for the ideal Einstein crystal with fixed center of mass.
        !   We set the thermal de Broglie wavelength to equal unity (i.e., our unit of length).
            N = REAL(NPART,DP)
            A_EIN_T  = 1.5_dp*(1.0_dp-1.0_dp/N)*LOG((BETAKB*LAMTR)/PI) - (3.0_dp/(2.0_dp*N))*LOG(N)
        !   Calculate the contribution to the free energy from the translational
        !   free energy of the ideal Einstein crystal with fixed center of mass and the
        !   free energy change between an unconstrained solid and the solid with fixed center of mass.
        !   The free energy is in the units of β/N.
            A0 = A_EIN_T + LOG(RHO)/N
            PRINT *, "A0_tr:", A0
            STOP
        ELSEIF(FLID==2) THEN
            U_EIN_TR = 0.0_dp
            U_EIN_OR = 0.0_dp
        ENDIF
    ENDIF

!   *************************************************************************
!   Set up the variables to calculate the free energy of arbitrary systems
!   using the  the Schilling-Schmid method.
!   *************************************************************************
    IF(SCHSMIT) THEN
        LAMOR = (LAM_SS * ITA_EIN)
        PRINT *, LAM_SS, LAMOR
        U_LATT = PE
        IF(SSID==2) THEN
            U_SS     = -REAL(NPART,DP)
            U_EIN_OR = 0.0_dp
        ENDIF
    ENDIF

!   *************************************************************************
!   Calculate the initial total collective density of the system when 
!   performing an interface-pinning simulation.    
!   *************************************************************************
    IF(PINT) THEN
        IF(COLLDT) THEN
            CALL COLLDENSFIELD(TRQ,RQ,-1)
        ELSE
            TRQ = BOX(PINDIM)
        ENDIF
    ENDIF

!   *************************************************************************
!   Initialise the umbrella parameters when performing a nuclues-size pinning 
!   simulation.
!   *************************************************************************
    IF(NUCSEEDT) THEN
        ALLOCATE( NUCCOUNT(0:NPART), TOTNUCCOUNT(0:NPART), SUMNUCCOUNT(0:NPART) )
        ALLOCATE( AVNUCCOUNTBLK(0:NPART), AVNUCCOUNT(0:NPART), STDNUCCOUNT(0:NPART), STDNUCCOUNTSUM(0:NPART) )
        NUCCOUNT = 0
        TOTNUCCOUNT = 0
        CALL UMBRELLA_BIAS(-1)
        PRINT *, "SEEDED WITH", CLSTRSIZE
    ENDIF

!   *************************************************************************
    PRINT *, 'initial potential energy =', PE
    IF(TAILCORT) THEN
        PRINT *, 'initial pressure =', (REAL(NPART,DP)/BETAKB + VIR)/VLM + CORP(RCUT,RHO), PRSFIX
    ELSE
        PRINT *, 'initial pressure =', (REAL(NPART,DP)/BETAKB + VIR)/VLM, PRSFIX
    ENDIF
    IF(PINT) PRINT *, 'initial Q(R) =', TRQ

!   Output xyz file to visualise the initial configuration
    CALL VIEWCONFIG('initialconfig.xyz')
!   *************************************************************************
END SUBROUTINE

!     ==============================================================================================
         
SUBROUTINE READCONFIG
!   This subroutine reads in a configuration
    USE COMMONS, ONLY: DP, NPART, NDIM, R, Q, RIGIDT, BOX, DENSITYT, PACKINGT, PI, VLM, CUBICT, ORTHORHOMBICT, RHO, UNITVECT
    USE COMMONS, ONLY: CONTINUET, NPTT, SCALET
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q

    IMPLICIT NONE
    
    INTEGER         :: J1, J2
    REAL(KIND=DP)   :: UV(NDIM), SCALE

    IF(CONTINUET) THEN
        OPEN (UNIT = 33, FILE = '../finalpos.dat', STATUS = 'OLD')
        IF (RIGIDT) OPEN (UNIT = 34, FILE = '../finalortn.dat', STATUS = 'UNKNOWN')
    ELSE
        OPEN (UNIT = 33, FILE = 'initialpos.inp', STATUS = 'OLD')
        IF (RIGIDT) OPEN (UNIT = 34, FILE = 'initialortn.inp', STATUS = 'UNKNOWN')
    ENDIF
    
    IF(CONTINUET .AND. NPTT) THEN
        OPEN (UNIT = 35, FILE = '../finalbox.dat', STATUS = 'OLD')
        OPEN (UNIT = 36, FILE = '../finaldns.dat', STATUS = 'OLD')

        READ(35,*) BOX
        READ(36,*) RHO

        IF (DENSITYT) THEN
            VLM = REAL(NPART,DP)/RHO
        ELSE IF (PACKINGT) THEN
            VLM = PI*REAL(NPART,DP) / (6.0_dp * RHO)
        ENDIF

        CLOSE(UNIT = 35)
        CLOSE(UNIT = 36)

    ELSE
        !   Calculate the volume the simulation cell
        IF (DENSITYT) THEN
            VLM = REAL(NPART,DP)/RHO
        ELSE IF (PACKINGT) THEN
            VLM = PI*REAL(NPART,DP) / (6.0_dp * RHO)
        ENDIF

    !   Calculate the edge lengths of the simulation cell
        IF(CUBICT) THEN
            BOX = VLM**(1.0_dp / 3.0_dp)
        ELSEIF(ORTHORHOMBICT) THEN 
            IF(SCALET) THEN 
                SCALE = (VLM/(BOX(1)*BOX(2)*BOX(3)) )**(1.0_dp/3.0_dp)
                BOX   = BOX * SCALE
            ELSE
                VLM = BOX(1)*BOX(2)*BOX(3)
                RHO = REAL(NPART,DP)/VLM
                IF (PACKINGT) RHO = PI/6.0_dp * RHO
            ENDIF
        ENDIF
    ENDIF

    R = 0.0_dp
    IF (RIGIDT) Q = 0.0_dp

    DO J1 = 1, NPART
        
        READ(33,*) (R(J2,J1), J2 = 1, NDIM)        ! Read in the translational coordinates
        R(:,J1) = R(:,J1) - BOX*ANINT(R(:,J1)/BOX) ! Take central image
        IF (RIGIDT) THEN
            IF(UNITVECT) THEN
            ! Read in the rotational coordinates as unit vectors
                READ(34,*) UV
                Q(:,J1) = UV_TO_Q ( UV )
            ELSE
            ! Read in the rotational coordinates as quaternions
                READ(34,*) Q(:,J1)
            ENDIF
        ENDIF

    ENDDO

    IF(ORTHORHOMBICT .AND. SCALET) R = R * SCALE

    CLOSE(UNIT = 33)
    CLOSE(UNIT = 34)

END SUBROUTINE READCONFIG
         
SUBROUTINE SPHRCL_SEED()

    USE COMMONS, ONLY: DP, NPART, NDIM, NUCSIZE, SEEDSIZE, RSEED, QSEED, SEEDRADIUS, R, Q, RIGIDT, BINARYT, SEEDID, RACEMICT

    IMPLICIT NONE

    INTEGER         :: J1
    REAL(KIND=DP)   :: RXTL(NDIM,NPART), QXTL(4,NPART)

!   Reset the position of the particles to be 0.0_dp
    RXTL        = R
    R           = 0.0_dp
    IF(RIGIDT) THEN
        QXTL = Q
        Q = 0.0_dp
    ENDIF
    SEEDRADIUS  = 0.0_dp

    DO
        DO J1 = 1, NPART
            IF( NORM2(RXTL(:,J1)) < SEEDRADIUS ) SEEDSIZE = SEEDSIZE + 1
        ENDDO

        IF(SEEDSIZE < NUCSIZE) THEN
            SEEDRADIUS = SEEDRADIUS + 0.01
            SEEDSIZE = 0
        ELSE
            
            ALLOCATE(RSEED(NDIM,SEEDSIZE))
            IF(RIGIDT) ALLOCATE(QSEED(4,SEEDSIZE))
            IF(BINARYT .OR. RACEMICT) ALLOCATE(SEEDID(SEEDSIZE))
            
            SEEDSIZE = 0
            
            DO J1 = 1, NPART
                IF( NORM2( RXTL(:,J1) ) < SEEDRADIUS ) THEN
                    
                    SEEDSIZE = SEEDSIZE + 1
                    RSEED(:,SEEDSIZE) = RXTL(:,J1)
                    IF(RIGIDT) QSEED(:,SEEDSIZE) = QXTL(:,J1)
                    
                    IF(BINARYT) THEN
                        IF(MOD(J1,2)==0) THEN
                            SEEDID(SEEDSIZE) = 0
                        ELSE
                            SEEDID(SEEDSIZE) = 1
                        ENDIF
                    ELSEIF(RACEMICT) THEN
                        IF(MOD(((J1-1)-MOD((J1-1),12))/12+1,2)==0) THEN
                            SEEDID(SEEDSIZE) = 0
                        ELSE
                            SEEDID(SEEDSIZE) = 1
                        ENDIF
                    ENDIF

                ENDIF
            ENDDO
            EXIT
        ENDIF
    ENDDO

END SUBROUTINE

!     ==============================================================================================

SUBROUTINE ADD_SEED()

    USE COMMONS, ONLY: DP, CDP, NPART, NDIM, SEEDSIZE, R, Q, RSEED, QSEED, SEEDSIZE, SEEDRADIUS, RIGIDT
    USE COMMONS, ONLY: BINARYT, SEEDID, BNRYRTIOT, BNRYA, RACEMICT, BOX
    USE ROTATIONS_MODULE, ONLY: RANDOM_QUATERNION

    IMPLICIT NONE

    INTEGER                      :: J1, J2, J3, NOVRLP, DELPART(NPART), NEWNPART, NDIFF, COUNT
    INTEGER                      :: NODD, NEVN, CODD, CEVN, SHFTO, SHFTE, SHIFT, FLUIDL, FLUIDU, NFILL
    REAL(KIND=DP)                :: RDIF(NDIM), RANDR(NDIM)
    REAL(KIND=DP), ALLOCATABLE   :: NEWR(:,:), NEWQ(:,:), EXTRAR(:,:)
    REAL(KIND=CDP)               :: DRAND48
    LOGICAL                      :: COMPLETET

    INTEGER                      :: OLD_COUNT
    REAL(KIND=DP)                :: OLDR(NDIM), DIFFR, DIFSV

    NOVRLP  = 0
    DELPART = -1
    DO J1 = 1, SEEDSIZE
        DO J2 = 1, NPART
            RDIF = RSEED(:,J1) - R(:,J2)
            IF(NORM2(RDIF) < 1.0_dp .OR. NORM2(R(:,J2)) <= SEEDRADIUS) THEN
                IF(.NOT. ANY(DELPART == J2) ) THEN
                    NOVRLP = NOVRLP + 1
                    DELPART(NOVRLP) = J2
                ENDIF
            ENDIF
        ENDDO
    ENDDO

    NEWNPART = NPART + SEEDSIZE - NOVRLP
    IF(NEWNPART > NPART) THEN
        NDIFF = NEWNPART-NPART
    ELSE
        IF(NEWNPART < NPART) THEN
            NDIFF = NPART-NEWNPART
            ALLOCATE(EXTRAR(NDIM,NDIFF))
            PRINT *, "NEW NPART=", NEWNPART
            DO J1 = 1, NDIFF
                DO
                    
                    COMPLETET = .FALSE.
                    COUNT = 0
                    RANDR = 0.0_dp
                    
                    DO J2 = 1, NDIM
                        RANDR(J2) = (DRAND48()-0.5_dp)*BOX(J2)
                    ENDDO

                    IF(NORM2(RANDR) > SEEDRADIUS) THEN
                        
                        DO J2 = 1, NPART
                            RDIF = RANDR - R(:,J2)
                            RDIF = RDIF - BOX*ANINT(RDIF/BOX)
                            DIFFR= NORM2(RDIF)
                            IF(NORM2(RDIF) < 1.0_dp) THEN
                                COUNT = COUNT + 1
                                DIFSV = DIFFR
                            ENDIF
                        ENDDO
                        IF(J1 > 1) THEN
                            DO J2 = 1, J1-1
                                RDIF = RANDR - EXTRAR(:,J2)
                                RDIF = RDIF - BOX*ANINT(RDIF/BOX)
                                DIFFR= NORM2(RDIF)
                                IF(NORM2(RDIF) < 1.0_dp) THEN
                                    COUNT = COUNT + 1
                                    DIFSV = DIFFR
                                ENDIF
                            ENDDO
                        ENDIF

                        IF(COUNT == 0) THEN
                            EXTRAR(:,J1) = RANDR
                            COMPLETET    = .TRUE.
                        ELSEIF(COUNT == 1) THEN
                            DO J3 = 1,1000
                                OLD_COUNT = COUNT
                                OLDR      = RANDR
                                COUNT     = 0
                                RANDR     = RANDR+(2.0_dp*DRAND48()-1.0_dp)*DIFSV
                                RANDR     = RANDR-BOX*ANINT(RANDR/BOX)
                                DO J2 = 1, NPART
                                    RDIF = RANDR - R(:,J2)
                                    RDIF = RDIF - BOX*ANINT(RDIF/BOX)
                                    IF(NORM2(RDIF) < 1.0_dp) COUNT = COUNT + 1
                                    IF(COUNT > OLD_COUNT) EXIT
                                ENDDO
                                IF(J1 > 1) THEN
                                    DO J2 = 1, J1-1
                                        RDIF = RANDR - EXTRAR(:,J2)
                                        RDIF = RDIF - BOX*ANINT(RDIF/BOX)
                                        IF(NORM2(RDIF) < 1.0_dp) COUNT = COUNT + 1
                                        IF(COUNT > OLD_COUNT) EXIT
                                    ENDDO
                                ENDIF
                                IF(COUNT == 0) THEN
                                    EXTRAR(:,J1) = RANDR
                                    COMPLETET    = .TRUE.
                                    EXIT
                                ELSEIF(COUNT > OLD_COUNT) THEN
                                    RANDR = OLDR
                                    COUNT = OLD_COUNT
                                ENDIF
                            ENDDO
                        ENDIF

                    ENDIF
                    IF(COMPLETET) EXIT
                ENDDO
            ENDDO
            IF(.NOT. COMPLETET) STOP "Number of particles has decreased following inclusion of seed"
        ELSE
            NDIFF = 0
        ENDIF
    ENDIF
    
    ALLOCATE(NEWR(NDIM,NPART))
    IF(RIGIDT) ALLOCATE(NEWQ(4,NPART))

    IF(BINARYT) THEN
        
        NODD = SUM(SEEDID)
        NEVN = SEEDSIZE-NODD
        CODD = 0
        CEVN = 0

        DO J1 = 1, SEEDSIZE
            IF(SEEDID(J1)==0) THEN
                CEVN  = CEVN + 1
                IF(BNRYRTIOT) THEN
                    SHIFT = CEVN
                ELSE
                    SHIFT = 2*CEVN
                ENDIF
                SHFTE = SHIFT
            ELSE
                CODD  = CODD + 1
                IF(BNRYRTIOT) THEN
                    SHIFT = CODD + BNRYA
                ELSE
                    SHIFT = 2*CODD - 1
                ENDIF
                SHFTO = SHIFT
            ENDIF
            NEWR(:,SHIFT) = RSEED(:,J1)
            IF(RIGIDT) NEWQ(:,SHIFT) = QSEED(:,J1)
        ENDDO

        COUNT = 0 

        IF(BNRYRTIOT) THEN
            IF(SHFTE > BNRYA) THEN
                STOP "System size is not big enough to accomodate seed"
            ELSE
                FLUIDL = SHFTE 
                FLUIDU = SHFTO
                NFILL  = BNRYA-SHFTE
            ENDIF
        ELSE
            IF(SHFTO < SHFTE) THEN
                FLUIDL = SHFTO
                FLUIDU = SHFTE
            ELSE   
                FLUIDL = SHFTE 
                FLUIDU = SHFTO
            ENDIF
            NFILL = (FLUIDU-FLUIDL-1)/2
        ENDIF

        DO J1 = 1, NPART
            IF(COUNT == NPART-SEEDSIZE) EXIT
            IF(.NOT. ANY(DELPART(1:NOVRLP) == J1)) THEN
                COUNT = COUNT + 1
                IF(COUNT <= NFILL) THEN
                    IF(BNRYRTIOT) THEN
                        SHIFT = FLUIDL+COUNT
                    ELSE
                        SHIFT = FLUIDL+2*COUNT
                    ENDIF
                ELSE
                    SHIFT = FLUIDU+COUNT-NFILL
                ENDIF
                NEWR(:,SHIFT) = R(:,J1)
                IF(RIGIDT) NEWQ(:,SHIFT) = Q(:,J1)
            ENDIF
        ENDDO

    ELSEIF(RACEMICT) THEN
        NODD = SUM(SEEDID)
        NEVN = SEEDSIZE-NODD
        CEVN = 0
        CODD = 0

        DO J1 = 1, SEEDSIZE
            IF(SEEDID(J1)==0) THEN
                CEVN  = CEVN + 1
                SHIFT = 12+CEVN
                SHFTE = SHIFT
                IF(MOD(CEVN,12)==0) CEVN = CEVN + 12
            ELSE
                CODD  = CODD + 1
                SHIFT = MOD(((CODD-1)-MOD((CODD-1),12))/12+1,2)+CODD-1
                SHFTO = SHIFT
                IF(MOD(CODD,12)==0) CODD = CODD + 12
            ENDIF
            NEWR(:,SHIFT) = RSEED(:,J1)
            IF(RIGIDT) NEWQ(:,SHIFT) = QSEED(:,J1)
        ENDDO

        COUNT = 0 

        FLUIDL = 12-MOD(SHFTO,12)
        FLUIDU = 12-MOD(SHFTE,12)
        NFILL  = FLUIDL+FLUIDU

        DO J1 = 1, NPART
            IF(COUNT == NPART-SEEDSIZE) EXIT
            IF(.NOT. ANY(DELPART(1:NOVRLP) == J1)) THEN
                COUNT = COUNT + 1
                IF(COUNT <= NFILL) THEN
                    IF(COUNT<=FLUIDL)THEN
                        SHIFT = SHFTO+COUNT
                    ELSE
                        SHIFT = SHFTE+COUNT-FLUIDL
                    ENDIF
                ELSE
                    SHIFT = SHIFT+1
                ENDIF
                NEWR(:,SHIFT) = R(:,J1)
                IF(RIGIDT) NEWQ(:,SHIFT) = Q(:,J1)
            ENDIF
        ENDDO

    ELSE
        DO J1 = 1, SEEDSIZE
            NEWR(:,J1) = RSEED(:,J1)
            IF(RIGIDT) NEWQ(:,J1) = QSEED(:,J1)
        ENDDO
        COUNT = 0
        DO J1 = 1, NPART
            IF(COUNT == NPART-SEEDSIZE) EXIT
            IF(.NOT. ANY(DELPART(1:NOVRLP) == J1)) THEN
                COUNT = COUNT + 1
                NEWR(:,SEEDSIZE+COUNT) = R(:,J1)
                IF(RIGIDT) NEWQ(:,SEEDSIZE+COUNT) = Q(:,J1)
            ENDIF
        ENDDO
    ENDIF

    DEALLOCATE(R)
    ALLOCATE(R(NDIM,NPART))
    R     = NEWR

    IF(RIGIDT) THEN
        DEALLOCATE(Q)
        ALLOCATE(Q(4,NPART))
        Q = NEWQ
    ENDIF

    IF(ALLOCATED(EXTRAR)) THEN
        PRINT *, NPART
        DO J1 = 1,NDIFF
            R(:,J1+(NPART-NDIFF)) = EXTRAR(:,J1)
            IF(RIGIDT) THEN
                Q(:,J1+(NPART-NDIFF)) = RANDOM_QUATERNION()
            ENDIF
        ENDDO
    ENDIF
    
END SUBROUTINE