!     ==============================================================================================
SUBROUTINE POLYMER_TRANSLATION()
!   Translate and entire polymer chain as a rigid body

    USE COMMONS, ONLY: DP, CDP, NDIM, NPART, R, BETAKB, BOX, SPET, CELLLISTT, PE, POLY_CLU, &
    INDXP, OVERLAPT, CLUSTERT,  MULTICHAINET, N_POLY, N_POLY_L, RCUT, SCALED_R, MAXDTR, CLUSTERMOVET
    USE COMMONS, ONLY: CLSTR, CLSTRSZ, CLSTRID, MAXDPTR, NTPMOVES, NTMOVES, ACCPTCTCP, ACCPTCT
    USE CELL_LIST, ONLY: C_INDEX, FINALIZE_LIST, INITIALIZE_LIST, MAKE_LIST
    USE ROTATIONS_MODULE, ONLY: RANDOM_ROTATE_QUATERNION, Q_TO_RM
    
    IMPLICIT NONE
    
    INTEGER             :: J1, J2, INDXPLY, CHAINID, CLSTRO(NPART), CLSTRSZO
    REAL(KIND=DP)       :: BEAD_ENERGY, ENEW, EOLD, RO(NDIM,NPART), DELE, BLTZMN, DISP(NDIM)
    REAL(KIND=CDP)      :: DRAND48
    LOGICAL             :: REJECTT

    SPET         = .TRUE.
    REJECTT      = .FALSE.
    OVERLAPT     = .FALSE.
    MULTICHAINET = .TRUE.

    IF(CLUSTERMOVET) THEN
        IF(DRAND48()<0.5_dp) THEN
            MULTICHAINET    = .TRUE.
        ELSE
            MULTICHAINET    = .FALSE.
        ENDIF
        
        IF(DRAND48()<0.9_dp) THEN
            CLUSTERT    = .TRUE.
        ELSE
            CLUSTERT    = .FALSE.
        ENDIF
    ENDIF

    EOLD = 0.0_dp

    IF(CLUSTERT) THEN
    !   Initialise the cluster size to be 1
        CLSTRSZ = 0
    !   Array containing IDs of particles in the cluster is initialised as an array of zeros.
        CLSTR  = 0
        CLSTRID = 1
    ENDIF
!   -------------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------------
    INDXPLY  = INT(N_POLY*DRAND48())
    CHAINID  = INDXPLY*N_POLY_L + 1
    POLY_CLU = INDXPLY

    DO J1 = CHAINID, CHAINID+N_POLY_L-1
        INDXP = J1
        CALL POTENTIAL(BEAD_ENERGY)
        EOLD = EOLD + BEAD_ENERGY
    ENDDO
!   -------------------------------------------------------------------------------------
    
!   -------------------------------------------------------------------------------------
    IF(CLUSTERT) THEN
    !   Calculate potential energy associated with capsomers bonded to the polymer chain.
        DO
            IF(CLSTRSZ==0) EXIT
            INDXP = CLSTR(CLSTRID)
        !   Calculate the pair energies associated with the current index particle and add 
        !   new particles to the cluster (this is currently done separately in each of the 
        !   pair potential subroutines).
            CALL POTENTIAL(BEAD_ENERGY)
            EOLD = EOLD + BEAD_ENERGY
        !   Update counter for the next iteration
            CLSTRID = CLSTRID + 1
        !   If the current particle index is equal the last entry in the cluster
        !   list then we have no new particles in the cluster, so we exit the loop.   
            IF(INDXP == CLSTR(CLSTRSZ)) EXIT
        ENDDO

        CLSTRO   = CLSTR
        CLSTRSZO = CLSTRSZ
    ENDIF
!   -------------------------------------------------------------------------------------
!   NEED TO ADD THE ABILITY TO ADD OTHER POLYMER CHAINS AND CAPSOMERS NOT BONDED TO THE
!   ORIGINAL POLYMER CHAIN TO THE CLUSTER!!!
!   -------------------------------------------------------------------------------------

    ENEW = 0.0_dp
    DO J1 = 1, NDIM
        IF(CLUSTERT) THEN
            DISP(J1) = (2.0_dp*DRAND48()-1.0_dp)*MAXDPTR!1.0_dp!MAXDPTR
        ELSE
            DISP(J1) = (2.0_dp*DRAND48()-1.0_dp)*MAXDTR
        ENDIF
    ENDDO

    IF(CLUSTERT) THEN
    !   Initialise the cluster size to be 1
        CLSTRSZ = 0
    !   Array containing IDs of particles in the cluster is initialised as an array of zeros.
        CLSTR  = 0
        CLSTRID = 1
    ENDIF

    ! CALL VIEWCONFIG('testtraj.xyz')

!   -------------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------------
    DO J1 = CHAINID, CHAINID+N_POLY_L-1
        RO(:,J1) = R(:,J1)
        R(:,J1)  = R(:,J1) + DISP
        R(:,J1)  = R(:,J1) - BOX*ANINT(R(:,J1)/BOX)
    ENDDO
!   -------------------------------------------------------------------------------------

!   -------------------------------------------------------------------------------------
    IF(CLUSTERT) THEN
        DO  J1 = 1, CLSTRSZO
            J2 = CLSTRO(J1)
            RO(:,J2) = R(:,J2)
            R(:,J2)  = R(:,J2) + DISP
            R(:,J2)  = R(:,J2) - BOX*ANINT(R(:,J2)/BOX)
        ENDDO
    ENDIF

    IF(CELLLISTT) THEN
        CALL FINALIZE_LIST()
        CALL INITIALIZE_LIST( NPART, RCUT/BOX )
        SCALED_R = 0.0_dp
        DO J1 = 1, NPART
            SCALED_R(:,J1) = R(:,J1)/BOX
        ENDDO
        CALL MAKE_LIST( NPART, SCALED_R )
    ENDIF
!   -------------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------------

!   -------------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------------
    DO J1 = CHAINID, CHAINID+N_POLY_L-1
        INDXP = J1
        CALL POTENTIAL(BEAD_ENERGY)
        ENEW = ENEW + BEAD_ENERGY
    ENDDO
!   -------------------------------------------------------------------------------------
    IF(CLUSTERT) THEN
    !   Calculate potential energy associated with capsomers bonded to the polymer chain.
        DO
            IF(CLSTRSZ==0) EXIT
            INDXP = CLSTR(CLSTRID)
            IF(.NOT. ANY(CLSTRO(1:CLSTRSZO)==INDXP)) THEN
                REJECTT = .TRUE.
                EXIT
            ELSE
            !   Calculate the pair energies associated with the current index particle and add 
            !   new particles to the cluster (this is currently done separately in each of the 
            !   pair potential subroutines).
                CALL POTENTIAL(BEAD_ENERGY)
                ENEW = ENEW + BEAD_ENERGY
            !   Update counter for the next iteration
                CLSTRID = CLSTRID + 1
            !   If the current particle index is equal the last entry in the cluster
            !   list then we have no new particles in the cluster, so we exit the loop.   
                IF(INDXP == CLSTR(CLSTRSZ)) EXIT
            ENDIF
        ENDDO

        IF(CLSTRSZ /= CLSTRSZO) REJECTT = .TRUE.
    ENDIF
!   -------------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------------
    ! PRINT *, "IS IT OVERLAPPING", OVERLAPT
    ! CALL VIEWCONFIG('testtraj.xyz')
    ! STOP

    IF(OVERLAPT) REJECTT = .TRUE.

    IF(.NOT. REJECTT) THEN
        DELE = ENEW - EOLD
    !   Metropolis acceptence criteria
        BLTZMN = EXP( -DELE*BETAKB )
        IF ( DRAND48() .GE. BLTZMN ) REJECTT = .TRUE.
    ENDIF

    IF (REJECTT) THEN
        DO J1 = CHAINID, CHAINID+N_POLY_L-1
            R(:,J1) = RO(:,J1)
        ENDDO

        IF(CLUSTERT) THEN
            DO J1 = 1, CLSTRSZO
                R(:,CLSTRO(J1)) = RO(:,CLSTRO(J1))
            ENDDO
        ENDIF

        IF(CELLLISTT) THEN
            CALL FINALIZE_LIST()
            CALL INITIALIZE_LIST( NPART, RCUT/BOX )
            SCALED_R = 0.0_dp
            DO J1 = 1, NPART
                SCALED_R(:,J1) = R(:,J1)/BOX
            ENDDO
            CALL MAKE_LIST( NPART, SCALED_R )
        ENDIF
    ELSE
        PE  = PE + DELE
        IF(CLUSTERT) THEN
            ACCPTCTCP = ACCPTCTCP + 1
        ELSE
            ACCPTCT = ACCPTCT + 1
        ENDIF
    ENDIF

    IF(CLUSTERT) THEN
        NTPMOVES = NTPMOVES + 1
    ELSE
        NTMOVES = NTMOVES + 1
    ENDIF

    MULTICHAINET = .FALSE.
    CLUSTERT    = .FALSE.
    OVERLAPT    = .FALSE.

END SUBROUTINE

SUBROUTINE POLYMER_ROTATE()
!   Translate and entire polymer chain as a rigid body

    USE COMMONS, ONLY: DP, CDP, NDIM, NPART, R, Q, BETAKB, BOX, SPET, CELLLISTT, PE, REFSITE, POLY_CLU, &
    INDXP, OVERLAPT, CLUSTERT,  MULTICHAINET, N_POLY, N_POLY_L, PI, RCHTT, NSITES, RBSITES, HLFPI, ISTEP
    USE COMMONS, ONLY: CLSTR, CLSTRSZ, CLSTRID, RCUT, SCALED_R, N_POLY_TOT, MAXDPRT, NRPMOVES, ACCPTCRCP, CLUSTERMOVET
    USE CELL_LIST, ONLY: C_INDEX, FINALIZE_LIST, INITIALIZE_LIST, MAKE_LIST
    USE ROTATIONS_MODULE, ONLY: RANDOM_VECTOR, QUATMUL, Q_TO_RM, RANDOM_ROTATE_QUATERNION
    
    IMPLICIT NONE
    
    INTEGER        :: J1, J2, J3, J4, J5, J6, INDXPLY, CHAINID, BEADI, BEADF, CLSTRO(NPART), CLSTRSZO, CLSTRCNT, CLN(NPART)
    INTEGER        :: BI, BF
    REAL(KIND=DP)  :: BEAD_ENERGY, EOLD, ENEW, ROLD(NDIM,NPART), QOLD(4,NPART), DELE, BLTZMN, TEST_E
    REAL(KIND=DP)  :: RIJ(NDIM), RRM(3,3), AXIS(NDIM), ANGLE, RMOVEMAP(4), RM(NDIM,NDIM), COM(NDIM)
    REAL(KIND=CDP) :: DRAND48
    LOGICAL        :: REJECTT!, UPT

    SPET         = .TRUE.
    REJECTT      = .FALSE.
    OVERLAPT     = .FALSE.
    MULTICHAINET = .FALSE.

    IF(CLUSTERMOVET) THEN
        IF(DRAND48()<0.5_dp) THEN
            MULTICHAINET = .TRUE.
            CLUSTERT     = .TRUE.
        ELSE
            MULTICHAINET    = .FALSE.
            IF(DRAND48()<0.9_dp) THEN
                CLUSTERT    = .TRUE.
            ELSE
                CLUSTERT    = .FALSE.
            ENDIF
        ENDIF
    ENDIF

!   Randomly select a polymer chain for which the move will be attempted
    INDXPLY  = INT(N_POLY*DRAND48())
    POLY_CLU = INDXPLY
!   Retrive the particle ID of the first particle in the chain.
    CHAINID = INDXPLY*N_POLY_L + 1
    BEADI   = CHAINID
    BEADF   = CHAINID + N_POLY_L - 1

    IF(CLUSTERT) THEN
    !   Initialise the cluster size to be 1
        CLSTRSZ = 0
    !   Array containing IDs of particles in the cluster is initialised as an array of zeros.
        CLSTR  = 0
        CLSTRID = 1
        CLSTRCNT = 0
        CLN = 0
    ENDIF

!   Compute the energy of the beads before the move
    EOLD = 0.0_dp
    DO J1 = BEADI, BEADF
        INDXP = J1
        CLSTRCNT = CLSTRSZ
        CALL POTENTIAL(BEAD_ENERGY)
        EOLD = EOLD + BEAD_ENERGY
        CLN(J1) = CLSTRSZ - CLSTRCNT
        ROLD(:,J1) = R(:,J1)
    ENDDO

    IF(CLUSTERT) THEN
    !   Calculate potential energy associated with capsomers bonded to the polymer chain.
        DO
            IF(CLSTRSZ==0) EXIT
            INDXP = CLSTR(CLSTRID)
            ROLD(:,INDXP) = R(:,INDXP)
            IF(INDXP > N_POLY_TOT) QOLD(:,INDXP) = Q(:,INDXP)
        !   Calculate the pair energies associated with the current index particle and add 
        !   new particles to the cluster (this is currently done separately in each of the 
        !   pair potential subroutines).
            CLSTRCNT = CLSTRSZ
            CALL POTENTIAL(BEAD_ENERGY)
            EOLD = EOLD + BEAD_ENERGY
            CLN(INDXP) = CLSTRSZ - CLSTRCNT
        !   Update counter for the next iteration
            CLSTRID = CLSTRID + 1
        !   If the current particle index is equal the last entry in the cluster
        !   list then we have no new particles in the cluster, so we exit the loop.   
            IF(INDXP == CLSTR(CLSTRSZ)) EXIT
        ENDDO

        CLSTRO   = CLSTR
        CLSTRSZO = CLSTRSZ
    ENDIF

!   Rotate the beads around the pivot bead.
!   Rotational moves
    AXIS = RANDOM_VECTOR( )                         ! Choose random unit vector
    ANGLE = ( 2.0_dp*DRAND48() - 1.0_dp ) * MAXDPRT ! Uniform random angle in desired range
! Standard formula for rotation quaternion, using half angles
    RMOVEMAP(1)   = COS(0.5_dp*ANGLE)
    RMOVEMAP(2:4) = SIN(0.5_dp*ANGLE)*AXIS
!   Extract the rotation matrix associated with the the rotational move
    RRM = Q_TO_RM(RMOVEMAP)

    IF(CLUSTERT) THEN
    !   Initialise the cluster size to be 1
        CLSTRSZ = 0
    !   Array containing IDs of particles in the cluster is initialised as an array of zeros.
        CLSTR  = 0
        CLSTRID = 1
    ENDIF

!   Generate "periodic image" of the section of the chain being rotated starting from
!   the seed particle. 
    COM = R(:,BEADI)
    DO J1 = BEADI+1, BEADF
        RIJ = R(:,J1) - R(:,J1-1)
        RIJ = RIJ - BOX*ANINT(RIJ/BOX)
        R(:,J1) = RIJ + R(:,J1-1)
        COM = COM + R(:,J1)
    ENDDO
    COM = COM / REAL(N_POLY_L,DP)

    IF(CLUSTERT) THEN
        CLSTRCNT = 1
    !   Loop through beads in randomly selected polymer for the move, and 
    !   add particles directly bonded to the chain to the periodic image.
        DO J3 = BEADI, BEADF
            J4 = CLSTRCNT 
            DO  J1 = J4, J4+CLN(J3)-1
                J2 = CLSTRO(J1)
                RIJ = R(:,J2) - R(:,J3)
                RIJ = RIJ - BOX*ANINT(RIJ/BOX)
                R(:,J2) = RIJ + R(:,J3)
                CLSTRCNT = CLSTRCNT + 1
            ENDDO
        ENDDO
    !   If larger clusters were formed then add the particles bonded to the 
    !   "first nearest neighbours" of the chain to the periodic image.
    !   Being careful to add polymer chains piece by piece so as not to have issues 
    !   associated with using the minimum image convention in a periodic box.
        DO J3 = 1, CLSTRSZO
            J4 = CLSTRCNT
            J5 = CLSTRO(J3)
            DO  J1 = J4, J4+CLN(J5)-1
                J2 = CLSTRO(J1)
                IF(J2>N_POLY_TOT) THEN
                    RIJ = R(:,J2) - R(:,J5)
                    RIJ = RIJ - BOX*ANINT(RIJ/BOX)
                    R(:,J2) = RIJ + R(:,J5)
                    CLSTRCNT = CLSTRCNT + 1
                ELSE
                    IF(MOD(J2-1,N_POLY_L)==0) THEN
                        BI = J2
                        BF = J2 + N_POLY_L - 1
                        RIJ = R(:,BI) - R(:,J5)
                        RIJ = RIJ - BOX*ANINT(RIJ/BOX)
                        R(:,BI) = RIJ + R(:,J5)
                        CLSTRCNT = CLSTRCNT + 1
                        DO J6 = BI+1, BF
                            RIJ = R(:,J6) - R(:,J6-1)
                            RIJ = RIJ - BOX*ANINT(RIJ/BOX)
                            R(:,J6) = RIJ + R(:,J6-1)
                            CLSTRCNT = CLSTRCNT + 1
                        ENDDO
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
    ENDIF

!   Perform the rotation around the COM of the polymer
    DO J1 = BEADI, BEADF
        RIJ = R(:,J1) - COM
        R(:,J1) = MATMUL(RRM,RIJ) + COM
        R(:,J1) = R(:,J1) - BOX*ANINT(R(:,J1)/BOX)
    ENDDO

    IF(CLUSTERT) THEN
        DO  J1  = 1, CLSTRSZO
            J2  = CLSTRO(J1)
            RIJ = R(:,J2) - COM
            R(:,J2) = MATMUL(RRM,RIJ) + COM
            R(:,J2) = R(:,J2) - BOX*ANINT(R(:,J2)/BOX)
            IF(J2 > N_POLY_TOT) THEN
            !   Update the orientation of the particle
                Q(:,J2) = QUATMUL(RMOVEMAP, QOLD(:,J2))
            !   Update the rigid body sites of the particle being displaced
                RM   = Q_TO_RM( Q(:,J2) )
                DO J3 = 1, NSITES
                    RBSITES(:,J3,J2) = MATMUL(RM,REFSITE(:,J3))
                ENDDO
            ENDIF
        ENDDO
    ENDIF

    IF(CELLLISTT) THEN
        CALL FINALIZE_LIST()
        CALL INITIALIZE_LIST( NPART, RCUT/BOX )
        SCALED_R = 0.0_dp
        DO J1 = 1, NPART
            SCALED_R(:,J1) = R(:,J1)/BOX
        ENDDO
        CALL MAKE_LIST( NPART, SCALED_R )
    ENDIF
    
!   Compute the energy of the beads after the move
    ENEW = 0.0_dp
    DO J1 = BEADI, BEADF
        INDXP = J1
        CALL POTENTIAL(BEAD_ENERGY)
        ENEW = ENEW + BEAD_ENERGY
    ENDDO

    IF(CLUSTERT) THEN
    !   Calculate potential energy associated with capsomers bonded to the polymer chain.
        DO
            IF(CLSTRSZ==0) EXIT
            INDXP = CLSTR(CLSTRID)
            IF(.NOT. ANY(CLSTRO(1:CLSTRSZO)==INDXP)) THEN
                REJECTT = .TRUE.
                EXIT
            ELSE
            !   Calculate the pair energies associated with the current index particle and add 
            !   new particles to the cluster (this is currently done separately in each of the 
            !   pair potential subroutines).
                CALL POTENTIAL(BEAD_ENERGY)
                ENEW = ENEW + BEAD_ENERGY
            !   Update counter for the next iteration
                CLSTRID = CLSTRID + 1
            !   If the current particle index is equal the last entry in the cluster
            !   list then we have no new particles in the cluster, so we exit the loop.   
                IF(INDXP == CLSTR(CLSTRSZ)) EXIT
            ENDIF
        ENDDO
        IF(CLSTRSZ /= CLSTRSZO) REJECTT = .TRUE.
    ENDIF

    IF(OVERLAPT) THEN
        REJECTT = .TRUE.
    ELSE
        DELE = ENEW - EOLD
    !   Metropolis acceptence criteria
        BLTZMN = EXP( -DELE*BETAKB )
        IF (DRAND48() .GE. BLTZMN) REJECTT = .TRUE.
    ENDIF

    IF (REJECTT) THEN
        DO J1 = BEADI, BEADF
            R(:,J1) = ROLD(:,J1)
        ENDDO

        IF(CLUSTERT) THEN
            DO J1 = 1, CLSTRSZO
                J2 = CLSTRO(J1)
                R(:,J2) = ROLD(:,J2)
                IF(J2 > N_POLY_TOT) THEN
                    Q(:,J2) = QOLD(:,J2)
                    RM   = Q_TO_RM( Q(:,J2) )
                    DO J3 = 1, NSITES
                        RBSITES(:,J3,J2) = MATMUL(RM,REFSITE(:,J3))
                    ENDDO
                ENDIF
            ENDDO
        ENDIF

        IF(CELLLISTT) THEN
            CALL FINALIZE_LIST()
            CALL INITIALIZE_LIST( NPART, RCUT/BOX )
            SCALED_R = 0.0_dp
            DO J1 = 1, NPART
                SCALED_R(:,J1) = R(:,J1)/BOX
            ENDDO
            CALL MAKE_LIST( NPART, SCALED_R )
        ENDIF
    ELSE
        PE  = PE + DELE
        IF(CLUSTERT) ACCPTCRCP = ACCPTCRCP + 1
    ENDIF

    IF(CLUSTERT) NRPMOVES = NRPMOVES + 1
    
    SPET         = .FALSE.
    OVERLAPT     = .FALSE.
    CLUSTERT     = .FALSE.
    MULTICHAINET = .FALSE.

END SUBROUTINE

SUBROUTINE POLYMER_RATCHET()
!   Translate and entire polymer chain as a rigid body

    USE COMMONS, ONLY: DP, CDP, NDIM, NPART, R, Q, BETAKB, BOX, SPET, CELLLISTT, PE, REFSITE, OBJCTT, &
    INDXP, OVERLAPT, CLUSTERT,  MULTICHAINET, N_POLY, N_POLY_L, PI, RCHTT, NSITES, RBSITES, HLFPI, ISTEP, CLUSTERMOVET
    USE COMMONS, ONLY: CLSTR, CLSTRSZ, CLSTRID, RCUT, SCALED_R, N_POLY_TOT, MAXDPRATT, NRATCPMOVES, ACCPTCRATCP
    USE CELL_LIST, ONLY: C_INDEX, FINALIZE_LIST, INITIALIZE_LIST, MAKE_LIST
    USE ROTATIONS_MODULE, ONLY: RANDOM_VECTOR, QUATMUL, Q_TO_RM, RANDOM_ROTATE_QUATERNION
    
    IMPLICIT NONE
    
    INTEGER        :: J1, J2, J3, J4, INDXPLY, CHAINID, BEAD_SD, BEADI, BEADF, CLSTRO(NPART), CLSTRSZO, CLSTRCNT, CLN(N_POLY_TOT)
    REAL(KIND=DP)  :: BEAD_ENERGY, EOLD, ENEW, ROLD(NDIM,NPART), QOLD(4,NPART), DELE, BLTZMN
    REAL(KIND=DP)  :: RIJ(NDIM), RRM(3,3), AXIS(NDIM), ANGLE, RMOVEMAP(4), RM(NDIM,NDIM)
    REAL(KIND=CDP) :: DRAND48
    LOGICAL        :: REJECTT, UPT

    SPET         = .TRUE.
    RCHTT        = .TRUE.
    REJECTT      = .FALSE.
    OVERLAPT     = .FALSE.
    MULTICHAINET = .FALSE.
    
    IF(CLUSTERMOVET) THEN
        IF(DRAND48()<0.5_dp) THEN
            CLUSTERT    = .TRUE.
        ELSE
            CLUSTERT    = .FALSE.
        ENDIF
    ENDIF

    UPT         = .FALSE.

!   Randomly select a polymer chain for which the move will be attempted
    INDXPLY = INT(N_POLY*DRAND48())
!   Retrive the particle ID of the first particle in the chain.
    CHAINID = INDXPLY*N_POLY_L + 1
!   Randomly select a particle in the chain (excluding the first and last beads) to act as the 
!   pivot particle for the move.
    BEAD_SD = CHAINID + INT((N_POLY_L-1)*DRAND48()) + 1

!   Randomly select whether to move beads "up-stream" from the pivot or "down-stream"
    IF(DRAND48()<0.5_dp) THEN
    !   Move up-stream
        BEADI = BEAD_SD + 1
        BEADF = CHAINID + N_POLY_L - 1
        UPT = .TRUE.
    ELSE
    !   Move down-stream
        BEADI = CHAINID
        BEADF = BEAD_SD - 1
    ENDIF

    IF(CLUSTERT) THEN
    !   Initialise the cluster size to be 1
        CLSTRSZ = 0
    !   Array containing IDs of particles in the cluster is initialised as an array of zeros.
        CLSTR  = 0
        CLSTRID = 1
        CLSTRCNT = 0
        CLN = 0
    ENDIF

!   Compute the energy of the beads before the move
    EOLD = 0.0_dp
    DO J1 = BEADI, BEADF
        INDXP = J1
        CLSTRCNT = CLSTRSZ
        CALL POTENTIAL(BEAD_ENERGY)
        EOLD = EOLD + BEAD_ENERGY
        CLN(J1) = CLSTRSZ - CLSTRCNT
        ROLD(:,J1) = R(:,J1)
    ENDDO

    IF(CLUSTERT) THEN
    !   Calculate potential energy associated with capsomers bonded to the polymer chain.
        DO CLSTRID = 1, CLSTRSZ
            INDXP = CLSTR(CLSTRID)
            ROLD(:,INDXP) = R(:,INDXP)
            IF(INDXP>N_POLY_TOT) QOLD(:,INDXP) = Q(:,INDXP)
        !   Calculate the pair energies associated with the current index particle and add 
        !   new particles to the cluster (this is currently done separately in each of the 
        !   pair potential subroutines).
            CALL POTENTIAL(BEAD_ENERGY)
            EOLD = EOLD + BEAD_ENERGY
        ENDDO
        CLSTRO   = CLSTR
        CLSTRSZO = CLSTRSZ
    ENDIF

!   Rotate the beads around the pivot bead.
!   Rotational moves
    AXIS = RANDOM_VECTOR( )                        ! Choose random unit vector
    ANGLE = ( 2.0_dp*DRAND48() - 1.0_dp ) * MAXDPRATT!PI/4.0_dp!MAXDPRATT!PI     ! Uniform random angle in desired range
! Standard formula for rotation quaternion, using half angles
    RMOVEMAP(1)   = COS(0.5_dp*ANGLE)
    RMOVEMAP(2:4) = SIN(0.5_dp*ANGLE)*AXIS
!   Extract the rotation matrix associated with the the rotational move
    RRM = Q_TO_RM(RMOVEMAP)

    IF(CLUSTERT) THEN
    !   Initialise the cluster size to be 1
        CLSTRSZ = 0
    !   Array containing IDs of particles in the cluster is initialised as an array of zeros.
        CLSTR  = 0
        CLSTRID = 1
    ENDIF

!   Generate "periodic image" of the section of the chain being rotated starting from
!   the seed particle. 
    IF(UPT) THEN
        DO J1 = BEADI, BEADF
            RIJ = R(:,J1) - R(:,J1-1)
            ! RIJ = RIJ - BOX*ANINT(RIJ/BOX)
            IF(OBJCTT) THEN
                RIJ(1:2) = RIJ(1:2) - BOX(1:2)*ANINT(RIJ(1:2)/BOX(1:2))
            ELSE
                RIJ = RIJ - BOX*ANINT(RIJ/BOX)
            ENDIF
            R(:,J1) = RIJ + R(:,J1-1)
        ENDDO
    ELSE
        DO J1 = BEADF, BEADI, -1
            RIJ = R(:,J1) - R(:,J1+1)
            ! RIJ = RIJ - BOX*ANINT(RIJ/BOX)
            IF(OBJCTT) THEN
                RIJ(1:2) = RIJ(1:2) - BOX(1:2)*ANINT(RIJ(1:2)/BOX(1:2))
            ELSE
                RIJ = RIJ - BOX*ANINT(RIJ/BOX)
            ENDIF
            R(:,J1) = RIJ + R(:,J1+1)
        ENDDO
    ENDIF

    IF(CLUSTERT) THEN
        CLSTRCNT = 1
        DO J3 = BEADI, BEADF
            J4 = CLSTRCNT
            DO  J1 = J4, J4+CLN(J3)-1
                J2 = CLSTRO(J1)
                RIJ = R(:,J2) - R(:,J3)
                ! RIJ = RIJ - BOX*ANINT(RIJ/BOX)
                IF(OBJCTT) THEN
                    RIJ(1:2) = RIJ(1:2) - BOX(1:2)*ANINT(RIJ(1:2)/BOX(1:2))
                ELSE
                    RIJ = RIJ - BOX*ANINT(RIJ/BOX)
                ENDIF
                R(:,J2) = RIJ + R(:,J3)
                CLSTRCNT = CLSTRCNT + 1
            ENDDO
        ENDDO
    ENDIF

!   Perform the rotation around the pivot bead
    DO J1 = BEADI, BEADF
        RIJ = R(:,J1) - R(:,BEAD_SD)
        R(:,J1) = MATMUL(RRM,RIJ) + R(:,BEAD_SD)
        ! R(:,J1) = R(:,J1) - BOX*ANINT(R(:,J1)/BOX)
        IF(OBJCTT) THEN
            R(1:2,J1) = R(1:2,J1) - BOX(1:2)*ANINT(R(1:2,J1)/BOX(1:2))
        ELSE
            R(:,J1) = R(:,J1) - BOX*ANINT(R(:,J1)/BOX)
        ENDIF
    ENDDO

    IF(CLUSTERT) THEN
        DO  J1  = 1, CLSTRSZO
            J2  = CLSTRO(J1)
            RIJ = R(:,J2) - R(:,BEAD_SD)
            R(:,J2) = MATMUL(RRM,RIJ) + R(:,BEAD_SD)
            IF(OBJCTT) THEN
                R(1:2,J2) = R(1:2,J2) - BOX(1:2)*ANINT(R(1:2,J2)/BOX(1:2))
            ELSE
                R(:,J2) = R(:,J2) - BOX*ANINT(R(:,J2)/BOX)
            ENDIF
            IF(J2>N_POLY_TOT) THEN
            !   Update the orientation of the particle
                Q(:,J2) = QUATMUL(RMOVEMAP, QOLD(:,J2))
            !   Update the rigid body sites of the particle being displaced
                RM   = Q_TO_RM( Q(:,J2) )
                DO J3 = 1, NSITES
                    RBSITES(:,J3,J2) = MATMUL(RM,REFSITE(:,J3))
                ENDDO
            ENDIF
        ENDDO
    ENDIF

    IF(CELLLISTT) THEN
        CALL FINALIZE_LIST()
        CALL INITIALIZE_LIST( NPART, RCUT/BOX )
        SCALED_R = 0.0_dp
        DO J1 = 1, NPART
            SCALED_R(:,J1) = R(:,J1)/BOX
        ENDDO
        CALL MAKE_LIST( NPART, SCALED_R )
    ENDIF
    
!   Compute the energy of the beads after the move
    ENEW = 0.0_dp
    DO J1 = BEADI, BEADF
        INDXP = J1
        CALL POTENTIAL(BEAD_ENERGY)
        ENEW = ENEW + BEAD_ENERGY
    ENDDO

    IF(CLUSTERT) THEN
    !   Calculate potential energy associated with capsomers bonded to the polymer chain.
        DO CLSTRID = 1, CLSTRSZ
            INDXP = CLSTR(CLSTRID)
            IF(.NOT. ANY(CLSTRO(1:CLSTRSZO)==INDXP)) THEN
                REJECTT = .TRUE.
                EXIT
            ELSE
            !   Calculate the pair energies associated with the current index particle and add 
            !   new particles to the cluster (this is currently done separately in each of the 
            !   pair potential subroutines).
                CALL POTENTIAL(BEAD_ENERGY)
                ENEW = ENEW + BEAD_ENERGY
            ENDIF
        ENDDO
        IF(CLSTRSZ /= CLSTRSZO) REJECTT = .TRUE.
    ENDIF


    IF(OVERLAPT) THEN
        REJECTT = .TRUE.
    ELSE
        DELE = ENEW - EOLD
    !   Metropolis acceptence criteria
        BLTZMN = EXP( -DELE*BETAKB )
        IF (DRAND48() .GE. BLTZMN) REJECTT = .TRUE.
    ENDIF

    IF (REJECTT) THEN
        DO J1 = BEADI, BEADF
            R(:,J1) = ROLD(:,J1)
        ENDDO

        IF(CLUSTERT) THEN
            DO J1 = 1, CLSTRSZO
                J2 = CLSTRO(J1)
                R(:,J2) = ROLD(:,J2)
                IF(J2>N_POLY_TOT) THEN
                    Q(:,J2) = QOLD(:,J2)
                    RM   = Q_TO_RM( Q(:,J2) )
                    DO J3 = 1, NSITES
                        RBSITES(:,J3,J2) = MATMUL(RM,REFSITE(:,J3))
                    ENDDO
                ENDIF
            ENDDO
        ENDIF

        IF(CELLLISTT) THEN
            CALL FINALIZE_LIST()
            CALL INITIALIZE_LIST( NPART, RCUT/BOX )
            SCALED_R = 0.0_dp
            DO J1 = 1, NPART
                SCALED_R(:,J1) = R(:,J1)/BOX
            ENDDO
            CALL MAKE_LIST( NPART, SCALED_R )
        ENDIF
    ELSE
        PE  = PE + DELE
        IF(CLUSTERT) ACCPTCRATCP = ACCPTCRATCP + 1
    ENDIF

    ! IF(CLUSTERT) 
    NRATCPMOVES = NRATCPMOVES + 1
    SPET     = .FALSE.
    RCHTT    = .FALSE.
    CLUSTERT = .FALSE.
    OVERLAPT = .FALSE.

END SUBROUTINE