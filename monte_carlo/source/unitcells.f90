SUBROUTINE GENERATE_UNIT_CELLS()

    USE COMMONS, ONLY: DP, NPART, LID, VLM, RIGIDT, EQUISEEDT, RANDQUATT, SETORTN, NLAYERS
    USE COMMONS, ONLY: R, Q, BOX, NUCX, NUCY, NUCZ, ORTHO_IN, DENSITYT, PACKINGT, RHO, RHOX, RHOF, PI
    USE ROTATIONS_MODULE, ONLY: RANDOM_QUATERNION

    IMPLICIT NONE

    INTEGER                     :: J1, I, IX, IY, IZ, M, IREF, NPCELL
    REAL(KIND=DP)               :: A1, A2, A3

    IF(LID > 5) THEN
        SETORTN = .TRUE.
    ELSE
        SETORTN = .FALSE.
    ENDIF
    
    IF(EQUISEEDT) THEN
        RHOF = RHO
        RHO  = RHOX
    ENDIF

    IF(ORTHO_IN) THEN
        IF (DENSITYT) THEN
            VLM = REAL(NPART,DP)/RHO
        ELSE IF (PACKINGT) THEN
            VLM = PI*REAL(NPART,DP) / (6.0_dp * RHO)
        ENDIF
    !   Calculate the edge lengths of the simulation cell
        BOX = BOX * ( (VLM/(BOX(1)*BOX(2)*BOX(3)) )**(1.0_dp/3.0_dp) )
    ENDIF

!====================================================================================
!       "CLOSE-PACKED" STRUCTURES FORMED BY ISOTROPIC SPHERES
!====================================================================================
    IF (LID == 1) THEN
!   Start from an fcc lattice configuration
        CALL FCC(NPCELL, A1, A2, A3)

    ELSEIF (LID == 2) THEN
!   Start from a hcp lattice configuration
        NPCELL = 4
        CALL HCP(A1, A2, A3) 

    ELSEIF (LID == 3) THEN
!   Start from a dhcp lattice configuration
        NPCELL = 8
        CALL DHCP(A1, A2, A3) 

    ELSE IF (LID == 4 .OR. LID == 19) THEN
!   Start from a bcc lattice configuration
        IF (LID == 4)  NPCELL = 2
        IF (LID == 19) NPCELL = 16
        CALL BCC(A1)
        A2 = A1
        A3 = A1

    ELSE IF (LID == 5) THEN
!   Start from a simple cubic lattice configuration
        NPCELL = 1
        CALL SC(A1) 
        A2 = A1
        A3 = A1

!====================================================================================
!       OPEN STRUCTURES FORMED BY TETRAHEDRAL PATCHY PARTICLES
!====================================================================================
    ELSEIF (LID == 6) THEN
!   Start from a cubic diamond lattice configuration
        CALL DMND(NPCELL, A1, A2, A3) 

    ELSEIF (LID == 7) THEN
!   Start from a hexagonal diamond lattice configuration
        NPCELL = 8
        CALL HEX_DMND(A1, A2, A3)
        
    ELSEIF (LID == 26) THEN
!   Start from a hexagonal diamond lattice configuration
        NPCELL = NLAYERS*4
        CALL RNDM_DMND(A1, A2, A3)

!====================================================================================
!       OPEN STRUCTURES FORMED BY TRIBLOCK PATCHY PARTICLES
!====================================================================================
!------------------------------------------------------------------------------------
!   Structures related to tetrastack
!------------------------------------------------------------------------------------
    ELSEIF (LID == 8) THEN
!   Start from a cubic tetrastack lattice configuration
        CALL CUBTETRSTK(NPCELL, A1, A2, A3)

    ELSEIF (LID == 9) THEN
!   Start from a double cubic tetrastack lattice configuration
        NPCELL = 4
        CALL DBLE_CUBTETRSTK(A1)
        A2 = A1
        A3 = A1

    ELSEIF (LID == 10) THEN
    !   Start from a hexagonal tetrastack lattice configuration
        NPCELL = 16
        CALL HEXTETRSTK(A1, A2, A3) 

!------------------------------------------------------------------------------------
!   Structures related to tetrahedral diamond
!------------------------------------------------------------------------------------
    ELSEIF (LID == 11) THEN
    !   Start from a tetrahedral diamond lattice configuration
        CALL TD_DMND(NPCELL, A1, A2, A3)

    ELSEIF (LID == 12) THEN
    !   Start from a double tetrahedral diamond lattice configuration
        NPCELL = 8
        CALL DBLE_TD_DMND(A1)
        A2 = A1
        A3 = A1

    ELSEIF (LID == 13) THEN
    !   Start from a hexagonal tetrahedral diamond lattice configuration
        NPCELL = 32
        CALL HEX_TD_DMND(A1,A2,A3)

!------------------------------------------------------------------------------------
!   Structures related to octahedral simple cubic
!------------------------------------------------------------------------------------
    ELSE IF (LID == 14) THEN
    !   Start from an octahedral simple cubic lattice configuration
        NPCELL = 6
        CALL OCTAHEDRAL_SC(A1)
        A2 = A1
        A3 = A1
    ELSE IF (LID == 15) THEN
    !   Start from an octahedral BCC lattice configuration
        NPCELL = 12
        CALL OCTAHEDRAL_BCC(A1)
        A2 = A1
        A3 = A1

!------------------------------------------------------------------------------------
!   Structures related to gyroid
!------------------------------------------------------------------------------------
    ELSE IF (LID == 16 .OR. LID == 17 .OR. LID == 18) THEN
    !   Start from a gyroid or double gyroid lattice configuration
        IF(LID==18) THEN
            NPCELL = 24
        ELSE
            NPCELL = 12
        ENDIF
        CALL GYROID(A1,LID)
        A2 = A1
        A3 = A1

    ELSE IF (LID == 20) THEN
    !   Start from a 2D Kagome lattice configuration
        NPCELL = 6
        CALL KAGOME(A1,A2,A3)

    ELSE IF (LID == 21) THEN
        !   Start from a Fddd lattice configuration
            NPCELL = 24
            CALL FDDD(A1,A2,A3)
!====================================================================================
!====================================================================================
!------------------------------------------------------------------------------------
!   Structures related to clathrates
!------------------------------------------------------------------------------------

    ELSE IF (LID == 23) THEN
        NPCELL = 46
        CALL CLATHRATE_I(A1,.FALSE.)
        A2 = A1
        A3 = A1
    
    ELSE IF (LID == 24) THEN
        NPCELL = 136
        CALL CLATHRATE_II(A1,.FALSE.)
        A2 = A1
        A3 = A1
    
    ELSE IF (LID == 25) THEN
        NPCELL = 12
        CALL CLATHRATE_III(A1,.FALSE.)
        A2 = A1
        A3 = A1

    ELSE IF (LID == 27) THEN
        NPCELL = 184
        CALL CLATHRATE_I(A1,.TRUE.)
        A2 = A1
        A3 = A1

    ELSE IF (LID == 28) THEN
        NPCELL = 544
        CALL CLATHRATE_II(A1,.TRUE.)
        A2 = A1
        A3 = A1

    ELSE IF (LID == 29) THEN
        NPCELL = 48
        CALL CLATHRATE_III(A1,.TRUE.)
        A2 = A1
        A3 = A1

    ENDIF

    IF(NPCELL*NUCX*NUCY*NUCZ /= NPART .AND. (.NOT. EQUISEEDT)) THEN
        STOP "No. of unit cells does not match No. of particles"
    ENDIF

    VLM = BOX(1)*BOX(2)*BOX(3)

    !   ------------------------------------------------------------------------
    IF (RIGIDT .AND. (.NOT. SETORTN)) THEN
    !   Assign random quaternion to orientational degrees of freedom
        Q = 0.0_dp
        DO J1 = 1, NPCELL
        !   Currently hard coded initial configuration for PGLJ system with octahedral patches.
        !   NEEDS FIXING, SHOULD NOT BE LEFT LIKE THIS.
            IF(RANDQUATT) THEN
                Q(:,J1) = RANDOM_QUATERNION()
            ELSE
                Q(1,J1) = 1.0_dp
            ENDIF
        ENDDO
    ENDIF
    
    !     Construct the lattice from the unit cell
    M = 0
    DO IZ = 1, NUCZ
        DO IY = 1, NUCY
            DO IX = 1, NUCX
                DO IREF = 1, NPCELL
                    R(1,IREF+M) = R(1,IREF) + A1*REAL((IX-1),DP)
                    R(2,IREF+M) = R(2,IREF) + A2*REAL((IY-1),DP)
                    R(3,IREF+M) = R(3,IREF) + A3*REAL((IZ-1),DP)
                    IF(RIGIDT) Q(:,IREF+M) = Q(:,IREF)
                ENDDO
                M = M + NPCELL
            ENDDO
        ENDDO
    ENDDO

    !     Shift the centre of the box to the origin
    DO I = 1, NPART
        R(:,I) = R(:,I) - BOX / 2.0_dp
        R(:,I) = R(:,I) - ANINT(R(:,I)/BOX)*BOX
    ENDDO
!   ------------------------------------------------------------------------       
    IF(EQUISEEDT) THEN
        RHO = RHOF
        CALL SPHRCL_SEED()
        CALL READCONFIG()
        CALL ADD_SEED()
    ENDIF

END SUBROUTINE

SUBROUTINE FCC(NPCELL, A1, A2, A3)
!     This subroutine constructs a fcc lattice

    USE COMMONS, ONLY: DP, NDIM, R, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, ORTHORHOMBICT, RIGIDT, Q, SETORTN
    USE ROTATIONS_MODULE, ONLY: ROTATE_QUATERNION

    IMPLICIT NONE
    INTEGER         :: NPCELL
    REAL(KIND=DP)   :: UCL, A1, A2, A3, ANGLE, AXIS(NDIM)

    !     Unit cell length
    IF (DENSITYT) THEN
        UCL     = (4.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL = (4.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    IF(ORTHORHOMBICT) THEN
        NPCELL = 6

        A1     = UCL / SQRT(2.0_dp)
        A2     = SQRT(1.5_dp) * UCL
        A3     = SQRT(3.0_dp) * UCL
        BOX    = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]

        R(:,1) = 0.0_dp
        R(:,2) = [0.5_dp*A1,      0.5_dp*A2,           0.0_dp      ]
        R(:,3) = [0.5_dp*A1,      A2/6.0_dp,          A3/3.0_dp    ]
        R(:,4) = [0.0_dp,    (2.0_dp/3.0_dp)*A2,      A3/3.0_dp    ]
        R(:,5) = [0.0_dp,         A2/3.0_dp,     (2.0_dp/3.0_dp)*A3]
        R(:,6) = [0.5_dp*A1, (5.0_dp/6.0_dp)*A2, (2.0_dp/3.0_dp)*A3]
        PRINT *, "ORTHO FCC UNIT CELL:", BOX
    ELSE
        NPCELL = 4
        
        A1     = UCL
        A2     = UCL
        A3     = UCL
        BOX    = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]

        R(:,1) = 0.0_dp
        R(:,2) = [0.5_dp*UCL, 0.5_dp*UCL,   0.0_dp  ]
        R(:,3) = [0.0_dp,     0.5_dp*UCL, 0.5_dp*UCL]
        R(:,4) = [0.5_dp*UCL,   0.0_dp,   0.5_dp*UCL]
        PRINT *, "CUBIC FCC UNIT CELL:", BOX

        IF(RIGIDT) THEN
        !   Initial orientations for ordered fcc lattice with tetrahedral patchy particles.
            SETORTN = .TRUE.
            Q = 0.0_dp
            ANGLE = PI/4.0_dp
            AXIS = [1.0_dp, 0.0_dp, 0.0_dp]
        !   Particles need to be rotated by +/-45 degrees.
            Q(1,1) = 1.0_dp
            Q(:,1) = ROTATE_QUATERNION(-ANGLE,AXIS,Q(:,1))
            Q(1,2) = 1.0_dp
            Q(:,2) = ROTATE_QUATERNION(ANGLE,AXIS,Q(:,2))
            Q(1,3) = 1.0_dp
            Q(:,3) = ROTATE_QUATERNION(-ANGLE,AXIS,Q(:,3))
            Q(1,4) = 1.0_dp
            Q(:,4) = ROTATE_QUATERNION(ANGLE,AXIS,Q(:,4))
        ENDIF
    ENDIF
    
END SUBROUTINE FCC

!     ==============================================================================================

SUBROUTINE HCP(A1, A2, A3)
!     This subroutine constructs a hcp lattice

    USE COMMONS, ONLY: DP, R, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI

    IMPLICIT NONE
    REAL(KIND=DP) :: TWOR2, CELL_VLM, A1, A2, A3

    TWOR2    = 2.0_dp * SQRT(2.0_dp) 

    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 4.0_dp / RHO
    ELSE
        CELL_VLM = PI*4.0_dp / (6.0_dp * RHO)
    ENDIF

    A1 = (CELL_VLM / TWOR2)**(1.0_dp/3.0_dp)
    A2 = SQRT(3.0_dp)*A1
    A3 = SQRT(8.0_dp/3.0_dp)*A1

    R(:,1) = 0.0_dp
    R(:,2) = [0.5_dp*A1,      A2/6.0_dp,     0.5_dp*A3]
    R(:,3) = [0.5_dp*A1,      0.5_dp*A2,     0.0_dp   ]
    R(:,4) = [0.0_dp,    (2.0_dp/3.0_dp)*A2, 0.5_dp*A3]

    BOX     = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
    PRINT *, "HCP UNIT CELL:", BOX

END SUBROUTINE HCP

SUBROUTINE DHCP(A1, A2, A3P)
!     This subroutine constructs a dhcp lattice

    USE COMMONS, ONLY: DP, R, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI

    IMPLICIT NONE
    REAL(KIND=DP) :: TWOR2, CELL_VLM, A1, A2, A3, A3P

    TWOR2    = 2.0_dp * SQRT(2.0_dp) 

    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 4.0_dp / RHO
    ELSE
        CELL_VLM = PI*4.0_dp / (6.0_dp * RHO)
    ENDIF

    A1  = (CELL_VLM / TWOR2)**(1.0_dp/3.0_dp)
    A2  = SQRT(3.0_dp)*A1
    A3  = SQRT(6.0_dp)/3.0_dp*A1
    A3P = 4.0_dp*A3

    R(:,1) = 0.0_dp 
    R(:,2) = [0.5_dp*A1,      0.5_dp*A2,     0.0_dp   ]
    R(:,3) = [0.5_dp*A1,      A2/6.0_dp,     A3       ]
    R(:,4) = [0.0_dp,    (2.0_dp/3.0_dp)*A2, A3       ]
    R(:,5) = [0.0_dp,           0.0_dp,      2.0_dp*A3]
    R(:,6) = [0.5_dp*A1,      0.5_dp*A2,     2.0_dp*A3]
    R(:,7) = [0.0_dp,        A2/3.0_dp,      3.0_dp*A3]
    R(:,8) = [0.5_dp*A1, (5.0_dp/6.0_dp)*A2, 3.0_dp*A3]

    BOX     = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3P*REAL(NUCZ,DP) ]
    PRINT *, "DHCP UNIT CELL:", BOX

END SUBROUTINE DHCP

!     ==============================================================================================
            
SUBROUTINE BCC(UCL)
!     This subroutine constructs a bcc lattice

    USE COMMONS, ONLY: DP, NSITES, LID, R, Q, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, RIGIDT, SETORTN

    IMPLICIT NONE
    REAL(KIND=DP) :: UCL, UCLQ, UCLH, UCTQ

    IF(LID == 4) THEN
        !     Unit cell length
        IF (DENSITYT) THEN
            UCL     = (2.0_dp/RHO)**(1.0_dp/3.0_dp)
        ELSE
            UCL = (2.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
        ENDIF

        BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]

        R(:,1) = 0.0_dp
        R(:,2) = [0.5_dp*UCL, 0.5_dp*UCL, 0.5_dp*UCL]

        IF(RIGIDT) THEN
        !   Initial orientations for ideal cubic diamond lattice with tetrahedral patchy particles.
            Q = 0.0_dp
            SETORTN = .TRUE.
        !   Particles sitting on FCC lattice points already have the correct orientation, 
        !   so we set the initial quaternions to correspond to 0.0_dp rotation.
            Q(1,1) = 1.0_dp
        !   Particles sitting on tetrahedral vacancies, need to be rotated by 90 degrees about 
        !   the z-axis.
            Q(:,2) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
        ENDIF
    ELSE
        !   Unit cell length (remember there are 8 particles in the unit cell of diamond)
        IF (DENSITYT) THEN
            UCL = (16.0_dp/RHO)**(1.0_dp/3.0_dp)
        ELSE
            UCL = (16.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
        ENDIF

        BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
        PRINT *, "Intepenetrating cubic diamond box:", BOX
        
        UCLQ = 0.25_dp*UCL
        UCLH = 0.5_dp*UCL
        UCTQ = 0.75_dp*UCL
    !   Build the unit cell
    !   Particles sitting on FCC lattice points
        R(:,1)  = 0.0_dp
        R(:,3)  = [UCLH, UCLH, 0.0_dp]
        R(:,5)  = [UCLH, 0.0_dp, UCLH]
        R(:,7)  = [0.0_dp, UCLH, UCLH]
    !   Particles sitting on tetrahedral vacancies 
        R(:,2)  = UCLQ
        R(:,4)  = [UCTQ, UCTQ, UCLQ]
        R(:,6)  = [UCTQ, UCLQ, UCTQ]
        R(:,8)  = [UCLQ, UCTQ, UCTQ]
    !   Particles sitting on FCC lattice points
        R(:,10) = UCLH
        R(:,12) = [UCLH, UCLH, 0.0_dp] + UCLH
        R(:,14) = [UCLH, 0.0_dp, UCLH] + UCLH
        R(:,16) = [0.0_dp, UCLH, UCLH] + UCLH
    !   Particles sitting on tetrahedral vacancies 
        R(:,9 ) = UCTQ
        R(:,11) = [UCTQ, UCTQ, UCLQ] + UCLH
        R(:,13) = [UCTQ, UCLQ, UCTQ] + UCLH
        R(:,15) = [UCLQ, UCTQ, UCTQ] + UCLH

        IF(RIGIDT) THEN
            IF(NSITES == 4) THEN
            !   Initial orientations for ideal cubic diamond lattice with tetrahedral patchy particles.
                Q = 0.0_dp
            !   Particles sitting on FCC lattice points already have the correct orientation, 
            !   so we set the initial quaternions to correspond to 0.0_dp rotation.
                Q(1,1)  = 1.0_dp
                Q(1,3)  = 1.0_dp
                Q(1,5)  = 1.0_dp
                Q(1,7)  = 1.0_dp
            !   Particles sitting on tetrahedral vacancies, need to be rotated by 90 degrees about 
            !   the z-axis.
                Q(:,2)  = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,4)  = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,6)  = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,8)  = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
            !   Particles sitting on FCC lattice points already have the correct orientation, 
            !   so we set the initial quaternions to correspond to 0.0_dp rotation.
                Q(1,10) = 1.0_dp
                Q(1,12) = 1.0_dp
                Q(1,14) = 1.0_dp
                Q(1,16) = 1.0_dp
            !   Particles sitting on tetrahedral vacancies, need to be rotated by 90 degrees about 
            !   the z-axis.
                Q(:,9 ) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,11) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,13) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,15) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
            ELSE
                Q      = 0.0_dp
                Q(1,:) = 1.0_dp
            ENDIF
        ENDIF
    ENDIF

END SUBROUTINE BCC
!     ==============================================================================================

SUBROUTINE SC(UCL)
!     This subroutine constructs a simple cubic lattice

    USE COMMONS, ONLY: DP, R, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI

    IMPLICIT NONE
    REAL(KIND=DP) :: UCL

    !     Unit cell length
    IF (DENSITYT) THEN
        UCL     = (1.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL = (PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
    R(:,1) = 0.0_dp

END SUBROUTINE SC

!     ==============================================================================================

SUBROUTINE DMND(NPCELL, A1, A2, A3)
!     This subroutine constructs a cubic diamond lattice

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, PI, NSITES, ORTHORHOMBICT, ORTHO_IN, DENSITYT
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL

    IMPLICIT NONE
    INTEGER       :: J1, NPCELL
    REAL(KIND=DP) :: A, A1, A2, A3, UCL, UCLQ, UCLH, UCTQ, Q1(4),Q2(4),Q3(4),Q4(4),Q5(4),Q6(4)

!   Unit cell length (remember there are 8 particles in the unit cell)
    IF (DENSITYT) THEN
        UCL     = (8.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL = (8.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    IF(ORTHORHOMBICT) THEN

        NPCELL = 12

        IF(ORTHO_IN) THEN
            A1 = BOX(1) / REAL(NUCX,DP)
            A2 = BOX(2) / REAL(NUCY,DP)
            A3 = BOX(3) / REAL(NUCZ,DP)
        ELSE
            A  = UCL
            A1 = A / SQRT(2.0_dp)
            A2 = SQRT(1.5_dp) * A
            A3 = SQRT(3.0_dp) * A
            BOX     = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
        ENDIF

        PRINT *, "Cubic Diamond orthorhombic box:", BOX

    !   Calculate the positions of the tetrahedral centers
        R(:,1)  = [0.0_dp,         0.0_dp,                 0.0_dp     ]
        R(:,3)  = [0.5_dp*A1,    0.5_dp*A2,                0.0_dp     ]

        R(:,4)  = [0.0_dp,    (2.0_dp/3.0_dp)*A2,  (1.0_dp/12.0_dp)*A3]
        R(:,2)  = [0.5_dp*A1,    A2/6.0_dp,        (1.0_dp/12.0_dp)*A3]
        
        R(:,7)  = [0.0_dp,    (2.0_dp/3.0_dp)*A2,        A3/3.0_dp    ]
        R(:,5)  = [0.5_dp*A1,    A2/6.0_dp,              A3/3.0_dp    ]

        R(:,6)  = [0.0_dp,       A2/3.0_dp,        (5.0_dp/12.0_dp)*A3]
        R(:,8)  = [0.5_dp*A1, (5.0_dp/6.0_dp)*A2,  (5.0_dp/12.0_dp)*A3]
        
        R(:,9)  = [0.0_dp,       A2/3.0_dp,        (2.0_dp/3.0_dp)*A3 ]
        R(:,11) = [0.5_dp*A1, (5.0_dp/6.0_dp)*A2,  (2.0_dp/3.0_dp)*A3 ]
        
        R(:,10) = [0.0_dp,         0.0_dp,              0.75_dp*A3    ]
        R(:,12) = [0.5_dp*A1,    0.5_dp*A2,             0.75_dp*A3    ]


        IF(RIGIDT) THEN
            DO J1 = 1, 12

                IF(MOD(J1,2) == 0) THEN
                    Q1 = UV_TO_Q ( [0.0_dp, 0.0_dp,-1.0_dp] )
                ELSE
                    Q1 = UV_TO_Q ( [0.0_dp, 0.0_dp, 1.0_dp] )
                ENDIF
                
                Q2 = [ COS(PI/8.0_dp),0.0_dp,0.0_dp,SIN(PI/8.0_dp) ]
                
                Q3 = [ COS(ATAN(SQRT(2.0_dp))/2.0_dp),SIN(ATAN(SQRT(2.0_dp))/2.0_dp),0.0_dp,0.0_dp ]
                Q4 = QUATMUL(Q1,Q2)

                Q5 = [ COS(PI/6.0_dp),0.0_dp,0.0_dp,SIN(PI/6.0_dp) ]
                Q6 = QUATMUL(Q3,Q4)
                                
                Q(:,J1) = QUATMUL(Q5,Q6)

            ENDDO
        ENDIF

    ELSE
        NPCELL  = 8
        A1      = UCL
        A2      = UCL
        A3      = UCL
        BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
        PRINT *, "Cubic Diamond cubic box:", BOX
        
        UCLQ = 0.25_dp*UCL
        UCLH = 0.5_dp*UCL
        UCTQ = 0.75_dp*UCL
    !   Build the unit cell
    !   Particles sitting on FCC lattice points
        R(:,1) = 0.0_dp
        R(:,3) = [UCLH, UCLH, 0.0_dp]
        R(:,5) = [UCLH, 0.0_dp, UCLH]
        R(:,7) = [0.0_dp, UCLH, UCLH]
    !   Particles sitting on tetrahedral vacancies 
        R(:,2) = UCLQ
        R(:,4) = [UCTQ, UCTQ, UCLQ]
        R(:,6) = [UCTQ, UCLQ, UCTQ]
        R(:,8) = [UCLQ, UCTQ, UCTQ]

        IF(RIGIDT) THEN
            IF(NSITES == 4) THEN
            !   Initial orientations for ideal cubic diamond lattice with tetrahedral patchy particles.
                Q = 0.0_dp
            !   Particles sitting on FCC lattice points already have the correct orientation, 
            !   so we set the initial quaternions to correspond to 0.0_dp rotation.
                Q(1,1) = 1.0_dp
                Q(1,3) = 1.0_dp
                Q(1,5) = 1.0_dp
                Q(1,7) = 1.0_dp
            !   Particles sitting on tetrahedral vacancies, need to be rotated by 90 degrees about 
            !   the z-axis.
                Q(:,2) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,4) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,6) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
                Q(:,8) = [1.0_dp/SQRT(2.0_dp), 0.0_dp, 0.0_dp, 1.0_dp/SQRT(2.0_dp)]
            ELSE
                Q      = 0.0_dp
                Q(1,:) = 1.0_dp
            ENDIF
        ENDIF
    ENDIF

END SUBROUTINE DMND

!     ==============================================================================================

SUBROUTINE HEX_DMND(A1, A2, A3)
    
    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL

    IMPLICIT NONE
    INTEGER       :: J1
    REAL(KIND=DP) :: CELL_VLM, A1, A2, A3, Q1(4),Q2(4),Q3(4),Q4(4),Q5(4),Q6(4)

    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 8.0_dp / RHO
    ELSE
        CELL_VLM = PI*8.0_dp / (6.0_dp * RHO)
    ENDIF
    
    A1 = ( CELL_VLM / (2.0_dp * SQRT(2.0_dp)) )**(1.0_dp/3.0_dp)
    A2 = SQRT(3.0_dp)*A1
    A3 = SQRT(8.0_dp/3.0_dp)*A1
    BOX     = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
    PRINT *, "HEX DMND UNIT CELL:", BOX

    R(:,1) = 0.0_dp
    R(:,5) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/2.0_dp),       0.0_dp     ]

    R(:,6) = [       0.0_dp,      A2*(2.0_dp/3.0_dp), A3*(1.0_dp/8.0_dp)]
    R(:,2) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/6.0_dp), A3*(1.0_dp/8.0_dp)]

    R(:,7) = [       0.0_dp,      A2*(2.0_dp/3.0_dp), A3*(1.0_dp/2.0_dp)]
    R(:,3) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/6.0_dp), A3*(1.0_dp/2.0_dp)]

    R(:,4) = [       0.0_dp,             0.0_dp,      A3*(5.0_dp/8.0_dp)]
    R(:,8) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/2.0_dp), A3*(5.0_dp/8.0_dp)]

    IF(RIGIDT) THEN
        DO J1 = 1, 8
            IF(J1 == 1 .OR. J1 == 5) THEN
                Q1 = UV_TO_Q ( [0.0_dp, 0.0_dp, 1.0_dp] )
                Q5 = [ COS(PI/6.0_dp),0.0_dp,0.0_dp,SIN(PI/6.0_dp) ]
            ELSEIF(J1 == 2 .OR. J1 == 6) THEN
                Q1 = UV_TO_Q ( [0.0_dp, 0.0_dp, -1.0_dp] )
                Q5 = [ COS(PI/6.0_dp),0.0_dp,0.0_dp,SIN(PI/6.0_dp) ]
            ELSEIF(J1 == 3 .OR. J1 == 7) THEN
                Q1 = UV_TO_Q ( [0.0_dp, 0.0_dp, 1.0_dp] )
                Q5 = [ COS(PI/3.0_dp),0.0_dp,0.0_dp,SIN(PI/3.0_dp) ]
            ELSEIF(J1 == 4 .OR. J1 == 8) THEN
                Q1 = UV_TO_Q ( [0.0_dp, 0.0_dp, -1.0_dp] )
                Q5 = [ COS(PI/3.0_dp),0.0_dp,0.0_dp,SIN(PI/3.0_dp) ]
            ENDIF
            
            Q2 = [ COS(PI/8.0_dp),0.0_dp,0.0_dp,SIN(PI/8.0_dp) ]
            
            Q3 = [ COS(ATAN(SQRT(2.0_dp))/2.0_dp),SIN(ATAN(SQRT(2.0_dp))/2.0_dp),0.0_dp,0.0_dp ]
            Q4 = QUATMUL(Q1,Q2)

            Q6 = QUATMUL(Q3,Q4)
                            
            Q(:,J1) = QUATMUL(Q5,Q6)
        ENDDO
    ENDIF

END SUBROUTINE HEX_DMND

SUBROUTINE RNDM_DMND(A1, A2, A3)
    
    USE COMMONS, ONLY: DP, CDP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, NLAYERS
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL
    USE ROTATIONS_MODULE, ONLY: RANDOM_QUATERNION

    IMPLICIT NONE
    INTEGER        :: J1, NLID
    REAL(KIND=DP)  :: UCL, A1, A2, A3, NL, AL, SL, S1, S2
    REAL(KIND=DP)  :: Q1(4),Q2(4),Q3(4),Q4(4),Q5(4),Q6(4),Q7(4),Q8(4)
    REAL(KIND=DP)  :: QS1(4),QS2(4),QE1(4),QE2(4),QE3(4),QE4(4)
    REAL(KIND=CDP) :: DRAND48
    LOGICAL        :: UNCONNECTT

    NL = REAL(NLAYERS,DP)
    UNCONNECTT = .TRUE.

!     Unit cell length
    IF (DENSITYT) THEN
        UCL = (8.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL = (8.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF
    
    A1  = UCL / SQRT(2.0_dp)
    
    A2  = SQRT(1.5_dp)*UCL
    
    AL  = SQRT(6.0_dp)/2.0_dp*A1
    SL  = SQRT(3.0_dp)*UCL/12.0_dp
    A3  = 0.5_dp*AL*REAL(NLAYERS,DP)+REAL(NLAYERS,DP)*SL

    BOX = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
    PRINT *, "RNDM DMND UNIT CELL:", BOX

    R(:,1) = [A1*0.0_dp, A2*(2.0_dp/3.0_dp), 0.0_dp*AL]
    R(:,3) = [A1*0.5_dp, A2*(1.0_dp/6.0_dp), 0.0_dp*AL]
    R(:,2) = [A1*0.0_dp, A2*(2.0_dp/3.0_dp), 0.5_dp*AL]
    R(:,4) = [A1*0.5_dp, A2*(1.0_dp/6.0_dp), 0.5_dp*AL]

    IF(RIGIDT) THEN
        Q1 = UV_TO_Q ( [0.0_dp, 0.0_dp, 1.0_dp] )
        Q2 = UV_TO_Q ( [0.0_dp, 0.0_dp,-1.0_dp] )
        Q3 = [ COS(PI/3.0_dp),0.0_dp,0.0_dp,SIN(PI/3.0_dp) ]
        Q4 = [ COS(PI/6.0_dp),0.0_dp,0.0_dp,SIN(PI/6.0_dp) ]
        Q5 = [ COS(PI/8.0_dp),0.0_dp,0.0_dp,SIN(PI/8.0_dp) ]
        Q6 = [ COS(ATAN(SQRT(2.0_dp))/2.0_dp),SIN(ATAN(SQRT(2.0_dp))/2.0_dp),0.0_dp,0.0_dp ]

        Q7 = QUATMUL(Q1,Q5); Q8 = QUATMUL(Q6,Q7)  
        QS1 = QUATMUL(Q4,Q8) ! Staggered orientation, pointing down
        Q7 = QUATMUL(Q2,Q5); Q8 = QUATMUL(Q6,Q7)  
        QS2 = QUATMUL(Q4,Q8) ! Staggered orientation, pointing
        Q7 = QUATMUL(Q1,Q5); Q8 = QUATMUL(Q6,Q7)  
        QE1 = QUATMUL(Q4,Q8) ! Eclipsed orientation, pointing down
        Q7 = QUATMUL(Q1,Q5); Q8 = QUATMUL(Q6,Q7)  
        QE2 = QUATMUL(Q3,Q8) ! Eclipsed orientation, pointing down
        Q7 = QUATMUL(Q2,Q5); Q8 = QUATMUL(Q6,Q7)  
        QE3 = QUATMUL(Q4,Q8) ! Eclipsed orientation, pointing up
        Q7 = QUATMUL(Q2,Q5); Q8 = QUATMUL(Q6,Q7)  
        QE4 = QUATMUL(Q3,Q8) ! Eclipsed orientation, pointing up 
    ENDIF   

    DO WHILE(UNCONNECTT) 
        S1 = 2.0_dp/3.0_dp
        S2 = 1.0_dp/6.0_dp

        NLID = 0

        DO J1 = 1, NLAYERS-1
            IF(DRAND48() < 0.5_dp) THEN
        !   Add eclipsed layer
                IF(NLID==0) THEN
                    NLID=1
                    S1 = S1 + 1.0_dp/3.0_dp
                    S2 = S2 + 1.0_dp/3.0_dp
                    IF(RIGIDT) THEN
                        Q(:,2+(J1-1)*4) = QE2
                        Q(:,4+(J1-1)*4) = QE2
                        Q(:,1+J1*4)     = QE4
                        Q(:,3+J1*4)     = QE4
                    ENDIF
                ELSE
                    NLID=0
                    S1 = S1 + 2.0_dp/3.0_dp
                    S2 = S2 + 2.0_dp/3.0_dp
                    IF(RIGIDT) THEN
                        Q(:,2+(J1-1)*4) =  QE1
                        Q(:,4+(J1-1)*4) =  QE1
                        Q(:,1+J1*4)     =  QE3
                        Q(:,3+J1*4)     =  QE3
                    ENDIF
                ENDIF
            ELSE
        !   Add staggered layer
                IF(NLID==0) THEN
                    S1 = S1 + 2.0_dp/3.0_dp
                    S2 = S2 + 2.0_dp/3.0_dp
                    IF(RIGIDT) THEN
                        Q(:,2+(J1-1)*4) = QS1
                        Q(:,4+(J1-1)*4) = QS1
                        Q(:,1+J1*4)     = QS2
                        Q(:,3+J1*4)     = QS2
                    ENDIF
                ELSE
                    S1 = S1 - 2.0_dp/3.0_dp
                    S2 = S2 - 2.0_dp/3.0_dp
                    IF(RIGIDT) THEN
                        Q(:,2+(J1-1)*4) = QE2
                        Q(:,4+(J1-1)*4) = QE2
                        Q(:,1+J1*4)     = QE4
                        Q(:,3+J1*4)     = QE4
                    ENDIF
                ENDIF
            ENDIF

            IF(S1>=1.0_dp) S1 = S1-1.0_dp
            IF(S2>=1.0_dp) S2 = S2-1.0_dp
            IF(S1< 0.0_dp) S1 = S1+1.0_dp
            IF(S2< 0.0_dp) S2 = S2+1.0_dp

            R(:,1+J1*4) = [0.0_dp,    A2*S1, 0.5_dp*REAL(J1,DP)*AL+REAL(J1,DP)*SL]
            R(:,3+J1*4) = [A1*0.5_dp, A2*S2, 0.5_dp*REAL(J1,DP)*AL+REAL(J1,DP)*SL]
            R(:,2+J1*4) = [0.0_dp,    A2*S1, 0.5_dp*AL*(1.0_dp+REAL(J1,DP))+REAL(J1,DP)*SL]
            R(:,4+J1*4) = [A1*0.5_dp, A2*S2, 0.5_dp*AL*(1.0_dp+REAL(J1,DP))+REAL(J1,DP)*SL]
        ENDDO

        IF(ABS(S1-2.0_dp/3.0_dp) > 1.e-05 .AND. ABS(S2-1.0_dp/6.0_dp) > 1.e-05) THEN
            UNCONNECTT = .FALSE.
            IF( (ABS(S1-1.0_dp)<1.e-05 .OR. ABS(S1-0.0_dp)<1.e-05) .AND. ABS(S2-0.5_dp) < 1.e-05) THEN
                Q(:,1)          = QS2
                Q(:,3)          = QS2
                Q(:,2+(J1-1)*4) = QS1
                Q(:,4+(J1-1)*4) = QS1
            ELSE
                Q(:,1)          = QE4
                Q(:,3)          = QE4
                Q(:,2+(J1-1)*4) = QE2
                Q(:,4+(J1-1)*4) = QE2
            ENDIF
        ENDIF
    ENDDO

    DO J1 = 1, NLAYERS-1
        R(3,1+J1*4) = R(3,1+J1*4) - AL
        R(3,3+J1*4) = R(3,3+J1*4) - AL
        R(3,2+J1*4) = R(3,2+J1*4) - AL
        R(3,4+J1*4) = R(3,4+J1*4) - AL
    ENDDO

    R(3,1) = 0.5_dp*REAL(NLAYERS-2,DP)*AL+REAL(NLAYERS,DP)*SL
    R(3,3) = 0.5_dp*REAL(NLAYERS-2,DP)*AL+REAL(NLAYERS,DP)*SL
    R(3,2) = 0.5_dp*AL*(1.0_dp+REAL(NLAYERS-2,DP))+REAL(NLAYERS,DP)*SL
    R(3,4) = 0.5_dp*AL*(1.0_dp+REAL(NLAYERS-2,DP))+REAL(NLAYERS,DP)*SL

END SUBROUTINE RNDM_DMND

SUBROUTINE CUBTETRSTK(NPCELL, A1, A2, A3)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, NSITES, ORTHORHOMBICT
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q

    IMPLICIT NONE
    INTEGER :: J1, NPCELL
    REAL(KIND=DP) :: UCL, UCL18, UCL38, UCL58, UCL78, TANG, A1, A2, A3, CN0, CN1, CN2, CN3
    REAL(KIND=DP), ALLOCATABLE :: U(:,:)

!   Unit cell length (remember there are 16 particles in the unit cell)
    IF (DENSITYT) THEN
        UCL = (16.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL = (16.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    IF(ORTHORHOMBICT) THEN
        NPCELL = 24

        A1     = UCL / SQRT(2.0_dp)
        A2     = SQRT(1.5_dp) * UCL
        A3     = SQRT(3.0_dp) * UCL
        BOX    = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]

        R(:,1)  = [1.0_dp/2.0_dp*A1, 0.0_dp       *A2,   0.0_dp       *A3]
        R(:,2)  = [1.0_dp/4.0_dp*A1, 3.0_dp/4.0_dp*A2,   0.0_dp       *A3]
        R(:,3)  = [3.0_dp/4.0_dp*A1, 3.0_dp/4.0_dp*A2,   0.0_dp       *A3]
        R(:,4)  = [1.0_dp/2.0_dp*A1, 5.0_dp/6.0_dp*A2,   1.0_dp/6.0_dp*A3]
        R(:,5)  = [3.0_dp/4.0_dp*A1, 1.0_dp/4.0_dp*A2,   0.0_dp       *A3]
        R(:,6)  = [0.0_dp       *A1, 1.0_dp/2.0_dp*A2,   0.0_dp       *A3]
        R(:,7)  = [1.0_dp/4.0_dp*A1, 1.0_dp/4.0_dp*A2,   0.0_dp       *A3]
        R(:,8)  = [0.0_dp       *A1, 1.0_dp/3.0_dp*A2,   1.0_dp/6.0_dp*A3]
        R(:,9)  = [1.0_dp/4.0_dp*A1, 5.0_dp/12.0_dp*A2,  1.0_dp/3.0_dp*A3]
        R(:,10) = [1.0_dp/2.0_dp*A1, 2.0_dp/3.0_dp *A2,  1.0_dp/3.0_dp*A3]
        R(:,11) = [3.0_dp/4.0_dp*A1, 5.0_dp/12.0_dp*A2,  1.0_dp/3.0_dp*A3]
        R(:,12) = [1.0_dp/2.0_dp*A1, 1.0_dp/2.0_dp *A2,  1.0_dp/2.0_dp*A3]
        R(:,13) = [1.0_dp/4.0_dp*A1, 11.0_dp/12.0_dp*A2, 1.0_dp/3.0_dp*A3]
        R(:,14) = [3.0_dp/4.0_dp*A1, 11.0_dp/12.0_dp*A2, 1.0_dp/3.0_dp*A3]
        R(:,15) = [0.0_dp       *A1, 1.0_dp/6.0_dp  *A2, 1.0_dp/3.0_dp*A3]
        R(:,16) = [0.0_dp       *A1, 0.0_dp         *A2, 1.0_dp/2.0_dp*A3]
        R(:,17) = [1.0_dp/2.0_dp*A1, 1.0_dp/3.0_dp *A2,  2.0_dp/3.0_dp*A3]
        R(:,18) = [1.0_dp/4.0_dp*A1, 1.0_dp/12.0_dp*A2,  2.0_dp/3.0_dp*A3]
        R(:,19) = [3.0_dp/4.0_dp*A1, 1.0_dp/12.0_dp*A2,  2.0_dp/3.0_dp*A3]
        R(:,20) = [1.0_dp/2.0_dp*A1, 1.0_dp/6.0_dp *A2,  5.0_dp/6.0_dp*A3]
        R(:,21) = [1.0_dp/4.0_dp*A1, 7.0_dp/12.0_dp*A2,  2.0_dp/3.0_dp*A3]
        R(:,22) = [3.0_dp/4.0_dp*A1, 7.0_dp/12.0_dp*A2,  2.0_dp/3.0_dp*A3]
        R(:,23) = [0.0_dp       *A1, 5.0_dp/6.0_dp *A2,  2.0_dp/3.0_dp*A3]
        R(:,24) = [0.0_dp       *A1, 2.0_dp/3.0_dp *A2,  5.0_dp/6.0_dp*A3]

        IF(RIGIDT) THEN
            ALLOCATE(U(3,24))
    !   Initial orientations for ideal cubic tetrastack lattice with triblock patchy particles.
            IF(NSITES == 2) THEN
                
                CN0 = 1.0_dp / 3.0_dp
                CN1 = SQRT(2.0_dp) / 3.0_dp
                CN2 = 2.0_dp*SQRT(2.0_dp) / 3.0_dp 
                CN3 = SQRT(6.0_dp) / 3.0_dp

                U(:,1)  = [0.0_dp,   CN2,   -CN0 ]
                U(:,2)  = [-CN3,    -CN1,   -CN0 ]
                U(:,3)  = [ CN3,    -CN1,   -CN0 ]
                U(:,4)  = [0.0_dp, 0.0_dp, 1.0_dp]

                U(:,5)  = [-CN3,    -CN1,   -CN0 ]
                U(:,6)  = [0.0_dp,   CN2,   -CN0 ]
                U(:,7)  = [ CN3,    -CN1,   -CN0 ]
                U(:,8)  = [0.0_dp, 0.0_dp, 1.0_dp]

                U(:,9)  = [-CN3,    -CN1,   -CN0 ]
                U(:,10) = [0.0_dp,   CN2,   -CN0 ]
                U(:,11) = [ CN3,    -CN1,   -CN0 ]
                U(:,12) = [0.0_dp, 0.0_dp, 1.0_dp]

                U(:,13) = [ CN3,    -CN1,   -CN0 ]
                U(:,14) = [-CN3,    -CN1,   -CN0 ]
                U(:,15) = [0.0_dp,   CN2,   -CN0 ]
                U(:,16) = [0.0_dp, 0.0_dp, 1.0_dp]

                U(:,17) = [0.0_dp,   CN2,   -CN0 ]
                U(:,18) = [-CN3,    -CN1,   -CN0 ]
                U(:,19) = [ CN3,    -CN1,   -CN0 ]
                U(:,20) = [0.0_dp, 0.0_dp, 1.0_dp]

                U(:,21) = [ CN3,    -CN1,   -CN0 ]
                U(:,22) = [-CN3,    -CN1,   -CN0 ]
                U(:,23) = [0.0_dp,   CN2,   -CN0 ]
                U(:,24) = [0.0_dp, 0.0_dp, 1.0_dp]
        !   Convert unit vectors into quaternions, where the roation is relative to the reference 
        !   configuration of the patchy particles.
                Q    = 0.0_dp
                DO J1 = 1, 24
                    Q(:,J1) = UV_TO_Q ( U(:,J1) )
                ENDDO
            ELSE
                Q      = 0.0_dp
                Q(1,:) = 1.0_dp
            ENDIF
        ENDIF
        
    ELSE
        NPCELL  = 16
        A1      = UCL
        A2      = UCL
        A3      = UCL
        BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
        
        UCL18 = 0.125_dp*UCL
        UCL38 = 0.375_dp*UCL
        UCL58 = 0.625_dp*UCL
        UCL78 = 0.875_dp*UCL
    !   Build the unit cell
    !   Particles sitting on cubic tetrastack lattice points
        R(:,1)  = [UCL78, UCL78, UCL18]
        R(:,2)  = [UCL58, UCL58, UCL18]
        R(:,3)  = [UCL78, UCL58, UCL38]
        R(:,4)  = [UCL58, UCL78, UCL38]
        R(:,5)  = [UCL78, UCL38, UCL58]
        R(:,6)  = [UCL58, UCL18, UCL58]
        R(:,7)  = [UCL78, UCL18, UCL78]
        R(:,8)  = [UCL58, UCL38, UCL78]
        R(:,9)  = [UCL38, UCL78, UCL58]
        R(:,10) = [UCL18, UCL58, UCL58]
        R(:,11) = [UCL38, UCL58, UCL78]
        R(:,12) = [UCL18, UCL78, UCL78]
        R(:,13) = [UCL38, UCL38, UCL18]
        R(:,14) = [UCL18, UCL18, UCL18]
        R(:,15) = [UCL38, UCL18, UCL38]
        R(:,16) = [UCL18, UCL38, UCL38]

        IF(RIGIDT) THEN
            ALLOCATE(U(3,24))
    !   Initial orientations for ideal cubic tetrastack lattice with triblock patchy particles.
            TANG = 1.0_dp/SQRT(3.0_dp)
            U(:,1)  = [ TANG,  TANG, -TANG]
            U(:,2)  = [-TANG, -TANG, -TANG]
            U(:,3)  = [ TANG, -TANG,  TANG]
            U(:,4)  = [-TANG,  TANG,  TANG]

            U(:,5)  = [ TANG,  TANG, -TANG]
            U(:,6)  = [-TANG, -TANG, -TANG]
            U(:,7)  = [ TANG, -TANG,  TANG]
            U(:,8)  = [-TANG,  TANG,  TANG]

            U(:,9)  = [ TANG,  TANG, -TANG]
            U(:,10) = [-TANG, -TANG, -TANG]
            U(:,11) = [ TANG, -TANG,  TANG]
            U(:,12) = [-TANG,  TANG,  TANG]

            U(:,13) = [ TANG,  TANG, -TANG]
            U(:,14) = [-TANG, -TANG, -TANG]
            U(:,15) = [ TANG, -TANG,  TANG]
            U(:,16) = [-TANG,  TANG,  TANG]
    !   Convert unit vectors into quaternions, where the roation is relative to the reference 
    !   configuration of the patchy particles.
            Q    = 0.0_dp
            DO J1 = 1, 16
                Q(:,J1) = UV_TO_Q ( U(:,J1) )
            ENDDO
        ENDIF
    ENDIF

END SUBROUTINE CUBTETRSTK

SUBROUTINE DBLE_CUBTETRSTK(UCL)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q

    IMPLICIT NONE
    INTEGER :: J1, NPCELL
    REAL(KIND=DP) :: UCL, UCL1, UCL2, TANG, A1, A2, A3
    REAL(KIND=DP), ALLOCATABLE :: U(:,:)

    IF (DENSITYT) THEN
        UCL = (4.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL = (4.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    
    NPCELL  = 4
    A1      = UCL
    A2      = UCL
    A3      = UCL
    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
    
    UCL1 = 0.25_dp*UCL*0.66_dp
    UCL2 = 0.75_dp*UCL*0.66_dp
!   Build the unit cell
!   Particles sitting on cubic tetrastack lattice points
    R(:,1)  = [UCL2, UCL1, UCL2]
    R(:,2)  = [UCL1, UCL2, UCL2]
    R(:,3)  = [UCL1, UCL1, UCL1]
    R(:,4)  = [UCL2, UCL2, UCL1]

    IF(RIGIDT) THEN
        ALLOCATE(U(3,4))
!   Initial orientations for ideal cubic tetrastack lattice with triblock patchy particles.
        TANG = 1.0_dp/SQRT(3.0_dp)
        U(:,1)  = [ TANG, -TANG,  TANG]
        U(:,2)  = [-TANG,  TANG,  TANG]
        U(:,3)  = [-TANG, -TANG, -TANG]
        U(:,4)  = [ TANG,  TANG, -TANG]
!   Convert unit vectors into quaternions, where the roation is relative to the reference 
!   configuration of the patchy particles.
        Q    = 0.0_dp
        DO J1 = 1, 4
            Q(:,J1) = UV_TO_Q ( U(:,J1) )
        ENDDO
    ENDIF

END SUBROUTINE DBLE_CUBTETRSTK

SUBROUTINE HEXTETRSTK(A1, A2, A3)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, NSITES
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q

    IMPLICIT NONE
    INTEGER       :: J1
    REAL(KIND=DP) :: TWOR2, CELL_VLM, A1, A2, A3, U(3,16), CN0, CN1, CN2, CN3

    TWOR2    = 2.0_dp * SQRT(2.0_dp) 

    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 16.0_dp / RHO
    ELSE
        CELL_VLM = PI*16.0_dp / (6.0_dp * RHO)
    ENDIF

    A1 = (CELL_VLM / TWOR2)**(1.0_dp/3.0_dp)
    A2 = SQRT(3.0_dp)*A1
    A3 = SQRT(8.0_dp/3.0_dp)*A1

    R(:,1)  = [1.0_dp/2.0_dp*A1,     0.0_dp,           0.0_dp      ]
    R(:,2)  = [1.0_dp/4.0_dp*A1, 3.0_dp/4.0_dp*A2,     0.0_dp      ]
    R(:,3)  = [3.0_dp/4.0_dp*A1, 3.0_dp/4.0_dp*A2,     0.0_dp      ]
    R(:,4)  = [1.0_dp/2.0_dp*A1, 5.0_dp/6.0_dp*A2, 1.0_dp/4.0_dp*A3]
    R(:,5)  = [0.0_dp,           1.0_dp/2.0_dp*A2,     0.0_dp      ]
    R(:,6)  = [3.0_dp/4.0_dp*A1, 1.0_dp/4.0_dp*A2,     0.0_dp      ]
    R(:,7)  = [1.0_dp/4.0_dp*A1, 1.0_dp/4.0_dp*A2,     0.0_dp      ]
    R(:,8)  = [0.0_dp,           1.0_dp/3.0_dp*A2, 1.0_dp/4.0_dp*A3]
    R(:,9)  = [1.0_dp/2.0_dp*A1,     0.0_dp,       1.0_dp/2.0_dp*A3]
    R(:,10) = [1.0_dp/4.0_dp*A1, 1.0_dp/4.0_dp*A2, 1.0_dp/2.0_dp*A3]
    R(:,11) = [3.0_dp/4.0_dp*A1, 1.0_dp/4.0_dp*A2, 1.0_dp/2.0_dp*A3]
    R(:,12) = [1.0_dp/2.0_dp*A1, 1.0_dp/6.0_dp*A2, 3.0_dp/4.0_dp*A3]
    R(:,13) = [3.0_dp/4.0_dp*A1, 3.0_dp/4.0_dp*A2, 1.0_dp/2.0_dp*A3]
    R(:,14) = [1.0_dp/4.0_dp*A1, 3.0_dp/4.0_dp*A2, 1.0_dp/2.0_dp*A3]
    R(:,15) = [0.0_dp,           1.0_dp/2.0_dp*A2, 1.0_dp/2.0_dp*A3]
    R(:,16) = [0.0_dp,           2.0_dp/3.0_dp*A2, 3.0_dp/4.0_dp*A3]

    BOX     = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
    PRINT *, "HT UNIT CELL:", BOX

    IF(RIGIDT) THEN
!   Initial orientations for ideal hexagonal tetrastack lattice with triblock patchy particles.
        IF(NSITES == 2) THEN
            CN0 = 1.0_dp / 3.0_dp
            CN1 = SQRT(2.0_dp) / 3.0_dp
            CN2 = 2.0_dp*SQRT(2.0_dp) / 3.0_dp 
            CN3 = SQRT(6.0_dp) / 3.0_dp

            U(:,1)  = [0.0_dp,   CN2,   -CN0 ] 
            U(:,2)  = [-CN3,    -CN1,   -CN0 ]
            U(:,3)  = [ CN3,    -CN1,   -CN0 ] 
            U(:,4)  = [0.0_dp, 0.0_dp, 1.0_dp]
            U(:,5)  = [0.0_dp,   CN2,   -CN0 ] 
            U(:,6)  = [-CN3,    -CN1,   -CN0 ]
            U(:,7)  = [ CN3,    -CN1,   -CN0 ] 
            U(:,8)  = [0.0_dp, 0.0_dp, 1.0_dp]
            U(:,9)  = [0.0_dp,  -CN2,   -CN0 ] 
            U(:,10) = [-CN3,     CN1,   -CN0 ]
            U(:,11) = [ CN3,     CN1,   -CN0 ] 
            U(:,12) = [0.0_dp, 0.0_dp, 1.0_dp]
            U(:,13) = [-CN3,     CN1,   -CN0 ] 
            U(:,14) = [ CN3,     CN1,   -CN0 ]
            U(:,15) = [0.0_dp,  -CN2,   -CN0 ] 
            U(:,16) = [0.0_dp, 0.0_dp, 1.0_dp]
    !   Convert unit vectors into quaternions, where the roation is relative to the reference 
    !   configuration of the patchy particles.
            Q    = 0.0_dp
            DO J1 = 1, 16
                Q(:,J1) = UV_TO_Q ( U(:,J1) )
            ENDDO
        ELSE
            Q      = 0.0_dp
            Q(1,:) = 1.0_dp
        ENDIF
    ENDIF

END SUBROUTINE HEXTETRSTK

SUBROUTINE TD_DMND(NPCELL, A1, A2, A3)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, CPPT, CPPLAM, ORTHORHOMBICT
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q

    IMPLICIT NONE
    INTEGER                    :: J1, J2, INDX, NPCELL
    REAL(KIND=DP)              :: CELL_VLM, A, A1, A2, A3, CN0, C1, C2, C3, C4, BL, RM(3,3)
    REAL(KIND=DP), ALLOCATABLE :: DMNDR(:,:), DMNDU(:,:)

    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 32.0_dp / RHO
    ELSE
        CELL_VLM = PI*32.0_dp / (6.0_dp * RHO)
    ENDIF

    IF(ORTHORHOMBICT) THEN

        NPCELL = 48
        ALLOCATE( DMNDR(3,12), DMNDU(3,8) )

        A  = (CELL_VLM)**(1.0_dp/3.0_dp)
        A1 = A / SQRT(2.0_dp)
        A2 = SQRT(1.5_dp) * A
        A3 = SQRT(3.0_dp) * A
        BOX     = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
        PRINT *, "TD DMND UNIT CELL:", BOX

    !   Calculate the positions of the tetrahedral centers
        DMNDR(:,1)  = [0.0_dp,         0.0_dp,                 0.0_dp     ]
        DMNDR(:,2)  = [0.5_dp*A1,    0.5_dp*A2,                0.0_dp     ]
        DMNDR(:,3)  = [0.5_dp*A1,    A2/6.0_dp,              A3/3.0_dp    ]
        DMNDR(:,4)  = [0.0_dp,    (2.0_dp/3.0_dp)*A2,        A3/3.0_dp    ]
        DMNDR(:,5)  = [0.0_dp,       A2/3.0_dp,        (2.0_dp/3.0_dp)*A3 ]
        DMNDR(:,6)  = [0.5_dp*A1, (5.0_dp/6.0_dp)*A2,  (2.0_dp/3.0_dp)*A3 ]
        DMNDR(:,7)  = [0.5_dp*A1,    A2/6.0_dp,        (1.0_dp/12.0_dp)*A3]
        DMNDR(:,8)  = [0.0_dp,    (2.0_dp/3.0_dp)*A2,  (1.0_dp/12.0_dp)*A3]
        DMNDR(:,9)  = [0.0_dp,       A2/3.0_dp,        (5.0_dp/12.0_dp)*A3]
        DMNDR(:,10) = [0.5_dp*A1, (5.0_dp/6.0_dp)*A2,  (5.0_dp/12.0_dp)*A3]
        DMNDR(:,11) = [0.0_dp,         0.0_dp,              0.75_dp*A3    ]
        DMNDR(:,12) = [0.5_dp*A1,    0.5_dp*A2,             0.75_dp*A3    ]

        IF(CPPT) THEN
            BL = CPPLAM
        ELSE
            BL = 1.01_dp
        ENDIF
        C1 = SQRT(2.0_dp/3.0_dp)
        C2 = SQRT(2.0_dp)/3.0_dp
        C3 = 1.0_dp/3.0_dp
        C4 = 2.0_dp*SQRT(2.0_dp) / 3.0_dp
        DMNDU(:,1) = [ C1,      C2,     C3  ]
        DMNDU(:,2) = [-C1,      C2,     C3  ]
        DMNDU(:,3) = [0.0_dp,  -C4,     C3  ]
        DMNDU(:,4) = [0.0_dp, 0.0_dp,-1.0_dp]
        DMNDU(:,5) = [-C1,     -C2,    -C3  ]
        DMNDU(:,6) = [ C1,     -C2,    -C3  ]
        DMNDU(:,7) = [0.0_dp,   C4,    -C3  ]
        DMNDU(:,8) = [0.0_dp, 0.0_dp, 1.0_dp]

        DO J1 = 0, 11
            IF(J1 < 6) THEN
                DO J2 = 0, 3
                    INDX = J1*4 + J2 + 1
                    R(:,INDX) = DMNDR(:,J1+1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*DMNDU(:,J2+1)
                    IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( DMNDU(:,J2+1) )
                ENDDO
            ELSE
                DO J2 = 4, 7
                    INDX = J1*4 + J2 - 3
                    R(:,INDX) = DMNDR(:,J1+1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*DMNDU(:,J2+1)
                    IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( DMNDU(:,J2+1) )
                ENDDO
            ENDIF
        ENDDO

        DO J1 = 1, 48
            R(1,J1) = R(1,J1) + BL/SQRT(3.0_dp)
            R(2,J1) = R(2,J1) + BL/SQRT(3.0_dp)
            R(3,J1) = R(3,J1) + BL*(SQRT(3.0_dp) / (2.0_dp*SQRT(2.0_dp)) )
        ENDDO

    ELSE

        NPCELL = 32
        ALLOCATE( DMNDR(3,8), DMNDU(3,8) )

        A1 = (CELL_VLM)**(1.0_dp/3.0_dp)
        A2 = A1
        A3 = A1
        BOX     = [ A1*REAL(NUCX,DP), A1*REAL(NUCY,DP), A1*REAL(NUCZ,DP) ]
        PRINT *, "TD DMND UNIT CELL:", BOX

    !   Calculate the positions of the tetrahedral centers
        DMNDR(:,1) = A1*0.0_dp
        DMNDR(:,2) = A1*[0.5_dp, 0.5_dp, 0.0_dp]
        DMNDR(:,3) = A1*[0.5_dp, 0.0_dp, 0.5_dp]
        DMNDR(:,4) = A1*[0.0_dp, 0.5_dp, 0.5_dp]
        DMNDR(:,5) = A1*0.25_dp
        DMNDR(:,6) = A1*[0.75_dp, 0.75_dp, 0.25_dp]
        DMNDR(:,7) = A1*[0.75_dp, 0.25_dp, 0.75_dp]
        DMNDR(:,8) = A1*[0.25_dp, 0.75_dp, 0.75_dp]

    !   Calculate the reference position of the particles which sit of the tetrahedral
    !   centers
        CN0 = 1.0_dp / SQRT(3.0_dp)
        RM = 0.0_dp
        RM(1,2) = -1.0_dp; RM(2,1) = 1.0_dp; RM(3,3) = 1.0_dp; 
        IF(CPPT) THEN
            BL = CPPLAM
        ELSE
            BL = 1.01_dp
        ENDIF
        DMNDU(:,1) = [ CN0,  CN0,  CN0]
        DMNDU(:,2) = [-CN0, -CN0,  CN0]
        DMNDU(:,3) = [ CN0, -CN0, -CN0]
        DMNDU(:,4) = [-CN0,  CN0, -CN0]
        DMNDU(:,5) = MATMUL(RM,DMNDU(:,1))
        DMNDU(:,6) = MATMUL(RM,DMNDU(:,2))
        DMNDU(:,7) = MATMUL(RM,DMNDU(:,3))
        DMNDU(:,8) = MATMUL(RM,DMNDU(:,4))

        DO J1 = 0, 7
            IF(J1 < 4) THEN
                DO J2 = 0, 3
                    INDX = J1*4 + J2 + 1
                    R(:,INDX) = DMNDR(:,J1+1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*DMNDU(:,J2+1)
                    IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( DMNDU(:,J2+1) )
                ENDDO
            ELSE
                DO J2 = 4, 7
                    INDX = J1*4 + J2 - 3
                    R(:,INDX) = DMNDR(:,J1+1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*DMNDU(:,J2+1)
                    IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( DMNDU(:,J2+1) )
                ENDDO
            ENDIF
        ENDDO

        DO J1 = 1, 32
            R(:,J1) = R(:,J1) + BL*(SQRT(3.0_dp) / (2.0_dp*SQRT(2.0_dp)) )
        ENDDO

    ENDIF

END SUBROUTINE TD_DMND

SUBROUTINE DBLE_TD_DMND(A1)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, CPPT, CPPLAM
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q

    IMPLICIT NONE
    INTEGER       :: J1, J2, INDX
    REAL(KIND=DP) :: CELL_VLM, A1, CN0, BL, RM(3,3), BCCR(3,2), BCCU(3,8)

    IF(MOD(NUCX,2) /= 0 .OR. MOD(NUCY,2) /= 0 .OR. MOD(NUCZ,2) /= 0) THEN
        STOP "CAN ONLY HAVE AN EVEN No. OF UNIT CELLS IN EACH DIRECTION"
    ENDIF
    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 8.0_dp / RHO
    ELSE
        CELL_VLM = PI*8.0_dp / (6.0_dp * RHO)
    ENDIF

    A1 = (CELL_VLM)**(1.0_dp/3.0_dp)
    BOX     = [ A1*REAL(NUCX,DP), A1*REAL(NUCY,DP), A1*REAL(NUCZ,DP) ]
    PRINT *, "TD DMND UNIT CELL:", BOX

!   Calculate the positions of the tetrahedral centers

    BCCR(:,1) = A1*0.0_dp
    BCCR(:,2) = A1*0.5_dp

!   Calculate the reference position of the particles which sit of the tetrahedral
!   centers
    CN0 = 1.0_dp / SQRT(3.0_dp)
    RM = 0.0_dp
    RM(1,2) = -1.0_dp; RM(2,1) = 1.0_dp; RM(3,3) = 1.0_dp; 
    IF(CPPT) THEN
        BL = CPPLAM
    ELSE
        BL = 1.01_dp
    ENDIF
    BCCU(:,1) = [ CN0,  CN0,  CN0]
    BCCU(:,2) = [-CN0, -CN0,  CN0]
    BCCU(:,3) = [ CN0, -CN0, -CN0]
    BCCU(:,4) = [-CN0,  CN0, -CN0]
    BCCU(:,5) = MATMUL(RM,BCCU(:,1))
    BCCU(:,6) = MATMUL(RM,BCCU(:,2))
    BCCU(:,7) = MATMUL(RM,BCCU(:,3))
    BCCU(:,8) = MATMUL(RM,BCCU(:,4))

    DO J1 = 0, 1
        IF(J1 == 0) THEN
            DO J2 = 0, 3
                INDX = J1*4 + J2 + 1
                R(:,INDX) = BCCR(:,J1+1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*BCCU(:,J2+1)
                IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( BCCU(:,J2+1) )
            ENDDO
        ELSE
            DO J2 = 4, 7
                INDX = J1*4 + J2 - 3
                R(:,INDX) = BCCR(:,J1+1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*BCCU(:,J2+1)
                IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( BCCU(:,J2+1) )
            ENDDO
        ENDIF
    ENDDO

    DO J1 = 1, 32
        R(:,J1) = R(:,J1) + BL*(SQRT(3.0_dp) / (2.0_dp*SQRT(2.0_dp)) )
    ENDDO

END SUBROUTINE DBLE_TD_DMND

SUBROUTINE HEX_TD_DMND(A1, A2, A3)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, CPPT, CPPLAM
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q

    IMPLICIT NONE
    INTEGER       :: J1, J2, INDX
    REAL(KIND=DP) :: CELL_VLM, A1, A2, A3, CN0, CN1, CN2, CN3, BL, HDR(3,8), HDU(3,8)

    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 32.0_dp / RHO
    ELSE
        CELL_VLM = PI*32.0_dp / (6.0_dp * RHO)
    ENDIF
 
    A1 = ( CELL_VLM / (2.0_dp * SQRT(2.0_dp)) )**(1.0_dp/3.0_dp)
    A2 = SQRT(3.0_dp)*A1
    A3 = SQRT(8.0_dp/3.0_dp)*A1
    BOX     = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
    PRINT *, "HEX TD DMND UNIT CELL:", BOX

!   Calculate the positions of the tetrahedral centers

    HDR(:,1) = 0.0_dp
    HDR(:,2) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/6.0_dp), A3*(1.0_dp/8.0_dp)]
    HDR(:,3) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/6.0_dp), A3*(1.0_dp/2.0_dp)]
    HDR(:,4) = [       0.0_dp,             0.0_dp,      A3*(5.0_dp/8.0_dp)]
    HDR(:,5) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/2.0_dp),       0.0_dp     ]
    HDR(:,6) = [       0.0_dp,      A2*(2.0_dp/3.0_dp), A3*(1.0_dp/8.0_dp)]
    HDR(:,7) = [       0.0_dp,      A2*(2.0_dp/3.0_dp), A3*(1.0_dp/2.0_dp)]
    HDR(:,8) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/2.0_dp), A3*(5.0_dp/8.0_dp)]

!   Calculate the reference position of the particles which sit of the tetrahedral
!   centers
    CN0 = SQRT(2.0_dp/3.0_dp)
    CN1 = SQRT(2.0_dp) / 3.0_dp
    CN2 = 2.0_dp*SQRT(2.0_dp) / 3.0_dp
    CN3 = 1.0_dp / 3.0_dp

    IF(CPPT) THEN
        BL = CPPLAM
    ELSE
        BL = 1.01_dp
    ENDIF
    HDU(:,1) = [  CN0,     CN1,    CN3  ]
    HDU(:,2) = [ -CN0,     CN1,    CN3  ]
    HDU(:,3) = [0.0_dp,   -CN2,    CN3  ]
    HDU(:,4) = [0.0_dp,  0.0_dp, -1.0_dp]

    HDU(:,5) = [  CN0,     CN1,   -CN3  ]
    HDU(:,6) = [ -CN0,     CN1,   -CN3  ]
    HDU(:,7) = [0.0_dp,   -CN2,   -CN3  ]
    HDU(:,8) = [0.0_dp,  0.0_dp,  1.0_dp]

    DO J1 = 0, 7
        IF(J1 == 0 .OR. J1 == 4) THEN
            
            DO J2 = 0, 3
                INDX = J1*4 + J2 + 1
                R(:,INDX) = HDR(:,J1+1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*HDU(:,J2+1)
                IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( HDU(:,J2+1) )
            ENDDO

        ELSEIF(J1 == 1 .OR. J1 == 5) THEN
            
            DO J2 = 0, 3
                INDX = J1*4 + J2 + 1
                R(:,INDX) = HDR(:,J1+1) - BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*HDU(:,J2+1)
                IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( -HDU(:,J2+1) )
            ENDDO

        ELSEIF(J1 == 2 .OR. J1 == 6) THEN
            
            DO J2 = 0, 3
                INDX = J1*4 + J2 + 1
                R(:,INDX) = HDR(:,J1+1) - BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*HDU(:,J2+5)
                IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( -HDU(:,J2+5) )
            ENDDO

        ELSEIF(J1 == 3 .OR. J1 == 7) THEN
            
            DO J2 = 0, 3
                INDX = J1*4 + J2 + 1
                R(:,INDX) = HDR(:,J1+1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*HDU(:,J2+5)
                IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( HDU(:,J2+5) )
            ENDDO

        ENDIF
    ENDDO

    DO J1 = 1, 32
        R(1,J1) = R(1,J1) + BL/SQRT(3.0_dp)
        R(2,J1) = R(2,J1) + BL/SQRT(3.0_dp)
        R(3,J1) = R(3,J1) + BL*(SQRT(3.0_dp) / (2.0_dp*SQRT(2.0_dp)) )
    ENDDO

END SUBROUTINE HEX_TD_DMND

SUBROUTINE GYROID(UCL,LID)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, NSITES, RBSITES, REFSITE
    USE COMMONS, ONLY: REFSITE, REFSITE2, SKEWAB, SKEWAB2, RACEMICT, KFRECT
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL, ROTATE_QUATERNION, Q_TO_RM

    IMPLICIT NONE
    INTEGER       :: J1, J2, LID, S
    REAL(KIND=DP) :: NP, UCL, SQRTTU, ANGLE, UV(3,24), R1(3), R2(3), NM(3), RMI(3,3), SK2

    IF(LID ==18) THEN
        NP = 24.0_dp
    ELSE
        NP = 12.0_dp
    ENDIF

    IF(RACEMICT) THEN
        SK2 = SKEWAB2
    ELSE
        SK2 = SKEWAB
    ENDIF


    IF (DENSITYT) THEN
        UCL     = (NP/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL = (NP*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
    PRINT *, "Gyroid cubic box:", BOX, LID
    
    IF(LID==16 .OR. LID ==18) THEN

    !   Build the unit cell
        R(:,1 ) = [0.375_dp, 0.000_dp, 0.750_dp]*UCL
        R(:,2 ) = [0.500_dp, 0.750_dp, 0.625_dp]*UCL
        R(:,3 ) = [0.250_dp, 0.875_dp, 0.500_dp]*UCL
        R(:,4 ) = [0.000_dp, 0.750_dp, 0.375_dp]*UCL
        R(:,5 ) = [0.750_dp, 0.625_dp, 0.500_dp]*UCL
        R(:,6 ) = [0.875_dp, 0.500_dp, 0.250_dp]*UCL
        R(:,7 ) = [0.750_dp, 0.375_dp, 0.000_dp]*UCL
        R(:,8 ) = [0.625_dp, 0.500_dp, 0.750_dp]*UCL
        R(:,9 ) = [0.500_dp, 0.250_dp, 0.875_dp]*UCL
        R(:,10) = [0.000_dp, 0.250_dp, 0.125_dp]*UCL
        R(:,11) = [0.125_dp, 0.000_dp, 0.250_dp]*UCL
        R(:,12) = [0.250_dp, 0.125_dp, 0.000_dp]*UCL

        IF(RIGIDT) THEN
            SQRTTU = SQRT(2.0_dp) / 2.0_dp

            UV(:,1 ) = [ 0.0_dp,  SQRTTU,  SQRTTU]
            UV(:,2 ) = [ SQRTTU, -SQRTTU,  0.0_dp]
            UV(:,3 ) = [-SQRTTU,  0.0_dp, -SQRTTU]
            UV(:,4 ) = [ SQRTTU,  SQRTTU,  0.0_dp]
            UV(:,5 ) = [-SQRTTU,  0.0_dp,  SQRTTU]
            UV(:,6 ) = [ 0.0_dp, -SQRTTU, -SQRTTU]
            UV(:,7 ) = [ SQRTTU,  0.0_dp,  SQRTTU]
            UV(:,8 ) = [ 0.0_dp,  SQRTTU, -SQRTTU]
            UV(:,9 ) = [-SQRTTU, -SQRTTU,  0.0_dp]
            UV(:,10) = [-SQRTTU,  SQRTTU,  0.0_dp]
            UV(:,11) = [ 0.0_dp, -SQRTTU,  SQRTTU]
            UV(:,12) = [ SQRTTU,  0.0_dp, -SQRTTU]

            DO J1 = 1, 12
                Q(:,J1) = UV_TO_Q( UV(:,J1)  )
                RMI     = Q_TO_RM( Q(:,J1) )
                DO J2 = 1, NSITES
                    RBSITES(:,J2,J1) = MATMUL(RMI,REFSITE(:,J2))
                ENDDO
            ENDDO

            IF(KFRECT) THEN
                R1      = R(:,1) - R(:,2)
                R1      = R1 - UCL*ANINT( R1/UCL )
                R2      = R(:,1) - R(:,3)
                R2      = R2 - UCL*ANINT( R2/UCL )
                NM(1)   = R1(2)*R2(3) - R1(3)*R2(2)
                NM(2)   = R1(3)*R2(1) - R1(1)*R2(3)
                NM(3)   = R1(1)*R2(2) - R1(2)*R2(1)
                NM      = NM / NORM2(NM)
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(-NM,RBSITES(:,2,1)))
                Q(:,1)  = ROTATE_QUATERNION(ANGLE-SKEWAB*PI/180_dp,-UV(:,1),Q(:,1))
                ANGLE   = ACOS(DOT_PRODUCT(NM,RBSITES(:,2,2)))
                Q(:,2)  = ROTATE_QUATERNION(ANGLE,UV(:,2),Q(:,2))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,3)))
                Q(:,3)  = ROTATE_QUATERNION(ANGLE,-UV(:,3),Q(:,3))

                R1      = R(:,4) - R(:,5)
                R1      = R1 - UCL*ANINT( R1/UCL )
                R2      = R(:,4) - R(:,6)
                R2      = R2 - UCL*ANINT( R2/UCL )
                NM(1)   = R1(2)*R2(3) - R1(3)*R2(2)
                NM(2)   = R1(3)*R2(1) - R1(1)*R2(3)
                NM(3)   = R1(1)*R2(2) - R1(2)*R2(1)
                NM      = NM / NORM2(NM)
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(-NM,RBSITES(:,2,4)))
                Q(:,4)  = ROTATE_QUATERNION(ANGLE-SKEWAB*PI/180_dp,-UV(:,4),Q(:,4))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,5)))
                Q(:,5)  = ROTATE_QUATERNION(ANGLE,UV(:,5),Q(:,5))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,6)))
                Q(:,6)  = ROTATE_QUATERNION(ANGLE,-UV(:,6),Q(:,6))

                R1      = R(:,7) - R(:,8)
                R1      = R1 - UCL*ANINT( R1/UCL )
                R2      = R(:,7) - R(:,9)
                R2      = R2 - UCL*ANINT( R2/UCL )
                NM(1)   = R1(2)*R2(3) - R1(3)*R2(2)
                NM(2)   = R1(3)*R2(1) - R1(1)*R2(3)
                NM(3)   = R1(1)*R2(2) - R1(2)*R2(1)
                NM      = NM / NORM2(NM)
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(-NM,RBSITES(:,2,7)))
                Q(:,7)  = ROTATE_QUATERNION(ANGLE-SKEWAB*PI/180_dp,-UV(:,7),Q(:,7))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,8)))
                Q(:,8)  = ROTATE_QUATERNION(ANGLE,UV(:,8),Q(:,8))
                ANGLE   = ACOS(DOT_PRODUCT(NM,RBSITES(:,2,9)))
                Q(:,9)  = ROTATE_QUATERNION(ANGLE,-UV(:,9),Q(:,9))
                
                R1      = R(:,10) - R(:,11)
                R2      = R(:,10) - R(:,12)
                NM(1)   = R1(2)*R2(3) - R1(3)*R2(2)
                NM(2)   = R1(3)*R2(1) - R1(1)*R2(3)
                NM(3)   = R1(1)*R2(2) - R1(2)*R2(1)
                NM      = NM / NORM2(NM)
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(-NM,RBSITES(:,2,10)))
                Q(:,10) = ROTATE_QUATERNION(ANGLE-SKEWAB*PI/180_dp,-UV(:,10),Q(:,10))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,11)))
                Q(:,11) = ROTATE_QUATERNION(ANGLE,-UV(:,11),Q(:,11))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,12)))
                Q(:,12) = ROTATE_QUATERNION(ANGLE,-UV(:,12),Q(:,12))
            ENDIF
        ENDIF
    ENDIF

    IF(LID==17 .OR. LID ==18) THEN

        IF(LID==18) THEN
            S = 12
        ELSE
            S = 0
        ENDIF

    !   Build the Gyroid (Left-Handed)? cell
        R(:,1 +S) = [0.000_dp, 0.250_dp, 0.625_dp]*UCL
        R(:,2 +S) = [0.250_dp, 0.375_dp, 0.500_dp]*UCL
        R(:,3 +S) = [0.125_dp, 0.500_dp, 0.750_dp]*UCL
        R(:,4 +S) = [0.375_dp, 0.500_dp, 0.250_dp]*UCL
        R(:,5 +S) = [0.500_dp, 0.750_dp, 0.125_dp]*UCL
        R(:,6 +S) = [0.250_dp, 0.625_dp, 0.000_dp]*UCL
        R(:,7 +S) = [0.750_dp, 0.125_dp, 0.500_dp]*UCL
        R(:,8 +S) = [0.500_dp, 0.250_dp, 0.375_dp]*UCL
        R(:,9 +S) = [0.625_dp, 0.000_dp, 0.250_dp]*UCL
        R(:,10+S) = [0.750_dp, 0.875_dp, 0.000_dp]*UCL
        R(:,11+S) = [0.875_dp, 0.000_dp, 0.750_dp]*UCL
        R(:,12+S) = [0.000_dp, 0.750_dp, 0.875_dp]*UCL

        IF(RIGIDT) THEN
            SQRTTU = SQRT(2.0_dp) / 2.0_dp

            UV(:,1 +S) = [-SQRTTU, -SQRTTU,  0.0_dp]
            UV(:,2 +S) = [ SQRTTU,  0.0_dp, -SQRTTU]
            UV(:,3 +S) = [ 0.0_dp,  SQRTTU,  SQRTTU]
            UV(:,4 +S) = [ 0.0_dp, -SQRTTU,  SQRTTU]
            UV(:,5 +S) = [ SQRTTU,  SQRTTU,  0.0_dp]
            UV(:,6 +S) = [-SQRTTU,  0.0_dp, -SQRTTU]
            UV(:,7 +S) = [ SQRTTU,  0.0_dp,  SQRTTU]
            UV(:,8 +S) = [-SQRTTU,  SQRTTU,  0.0_dp]
            UV(:,9 +S) = [ 0.0_dp, -SQRTTU, -SQRTTU]
            UV(:,10+S) = [-SQRTTU,  0.0_dp,  SQRTTU]
            UV(:,11+S) = [ 0.0_dp,  SQRTTU, -SQRTTU]
            UV(:,12+S) = [ SQRTTU, -SQRTTU,  0.0_dp]

            DO J1 = 1, 12
                Q(:,J1+S) = UV_TO_Q( UV(:,J1+S)  )
                RMI     = Q_TO_RM( Q(:,J1+S) )
                DO J2 = 1, NSITES
                    IF(RACEMICT) THEN
                        RBSITES(:,J2,J1+S) = MATMUL(RMI,REFSITE2(:,J2))
                    ELSE
                        RBSITES(:,J2,J1+S) = MATMUL(RMI,REFSITE(:,J2))
                    ENDIF
                ENDDO
            ENDDO

            IF(KFRECT) THEN
                R1      = R(:,1+S) - R(:,2+S)
                R1      = R1 - UCL*ANINT( R1/UCL )
                R2      = R(:,1+S) - R(:,3+S)
                R2      = R2 - UCL*ANINT( R2/UCL )
                NM(1)   = R1(2)*R2(3) - R1(3)*R2(2)
                NM(2)   = R1(3)*R2(1) - R1(1)*R2(3)
                NM(3)   = R1(1)*R2(2) - R1(2)*R2(1)
                NM      = NM / NORM2(NM)
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(-NM,RBSITES(:,2,1+S)))
                Q(:,1+S)  = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,1+S),Q(:,1+S))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(-NM,RBSITES(:,2,2+S)))
                Q(:,2+S)  = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,2+S),Q(:,2+S))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,3+S)))
                Q(:,3+S)  = ROTATE_QUATERNION(ANGLE,-UV(:,3+S),Q(:,3+S))

                R1      = R(:,4+S) - R(:,5+S)
                R1      = R1 - UCL*ANINT( R1/UCL )
                R2      = R(:,4+S) - R(:,6+S)
                R2      = R2 - UCL*ANINT( R2/UCL )
                NM(1)   = R1(2)*R2(3) - R1(3)*R2(2)
                NM(2)   = R1(3)*R2(1) - R1(1)*R2(3)
                NM(3)   = R1(1)*R2(2) - R1(2)*R2(1)
                NM      = NM / NORM2(NM)
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,4+S)))
                Q(:,4+S)  = ROTATE_QUATERNION(ANGLE,UV(:,4+S),Q(:,4+S))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,5+S)))
                Q(:,5+S)  = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,5+S),Q(:,5+S))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,6+S)))
                Q(:,6+S)  = ROTATE_QUATERNION(ANGLE,UV(:,6+S),Q(:,6+S))

                R1      = R(:,7+S) - R(:,8+S)
                R1      = R1 - UCL*ANINT( R1/UCL )
                R2      = R(:,7+S) - R(:,9+S)
                R2      = R2 - UCL*ANINT( R2/UCL )
                NM(1)   = R1(2)*R2(3) - R1(3)*R2(2)
                NM(2)   = R1(3)*R2(1) - R1(1)*R2(3)
                NM(3)   = R1(1)*R2(2) - R1(2)*R2(1)
                NM      = NM / NORM2(NM)
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(-NM,RBSITES(:,2,7+S)))
                Q(:,7+S)  = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,7+S),Q(:,7+S))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,8+S)))
                Q(:,8+S)  = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,8+S),Q(:,8+S))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,9+S)))
                Q(:,9+S)  = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,9+S),Q(:,9+S))
                
                R1      = R(:,10+S) - R(:,11+S)
                R1      = R1 - UCL*ANINT( R1/UCL )
                R2      = R(:,10+S) - R(:,12+S)
                R2      = R2 - UCL*ANINT( R2/UCL )
                NM(1)   = R1(2)*R2(3) - R1(3)*R2(2)
                NM(2)   = R1(3)*R2(1) - R1(1)*R2(3)
                NM(3)   = R1(1)*R2(2) - R1(2)*R2(1)
                NM      = NM / NORM2(NM)
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(-NM,RBSITES(:,2,10+S)))
                Q(:,10+S) = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,10+S),Q(:,10+S))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,11+S)))
                Q(:,11+S) = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,11+S),Q(:,11+S))
                ANGLE   = PI/2_dp - ACOS(DOT_PRODUCT(NM,RBSITES(:,2,12+S)))
                Q(:,12+S) = ROTATE_QUATERNION(ANGLE-SK2*PI/180_dp,-UV(:,12+S),Q(:,12+S))
            ENDIF
        ENDIF
    ENDIF

END SUBROUTINE GYROID

!     ==============================================================================================

SUBROUTINE KAGOME(A1, A2, A3)
    
    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL

    IMPLICIT NONE
    REAL(KIND=DP) :: IDEAL_CELL_VLM, CELL_VLM, A1, A2, A3, CN1, CN3

    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 6.0_dp / RHO
    ELSE
        CELL_VLM = PI*6.0_dp / (6.0_dp * RHO)
    ENDIF

    IDEAL_CELL_VLM = 4.0_dp*SQRT(3.0_dp)
    
    A1  = 2.00000000000001_dp
    A2  = SQRT(3.0_dp)*A1
    A3  = 1.00000000000001_dp*(CELL_VLM/IDEAL_CELL_VLM)
    BOX = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
    PRINT *, "KAGOME UNIT CELL:", BOX

    R(:,1) = 0.0_dp
    R(:,2) = [A1*(1.0_dp/2.0_dp),        0.0_dp,      0.0_dp]
    R(:,3) = [0.0_dp,             A2*(1.0_dp/2.0_dp), 0.0_dp]
    R(:,4) = [A1*(1.0_dp/2.0_dp), A2*(1.0_dp/2.0_dp), 0.0_dp]
    R(:,5) = [A1*(1.0_dp/4.0_dp), A2*(1.0_dp/4.0_dp), 0.0_dp]
    R(:,6) = [A1*(3.0_dp/4.0_dp), A2*(3.0_dp/4.0_dp), 0.0_dp]

    IF(RIGIDT) THEN
        CN1 = SQRT(2.0_dp) / 3.0_dp
        CN3 = SQRT(6.0_dp) / 3.0_dp
        Q(:,1) = UV_TO_Q ( [-CN3,    -CN1,   0.0_dp] )
        Q(:,2) = UV_TO_Q ( [ CN3,    -CN1,   0.0_dp] )
        Q(:,3) = UV_TO_Q ( [ CN3,    -CN1,   0.0_dp] )
        Q(:,4) = UV_TO_Q ( [-CN3,    -CN1,   0.0_dp] )
        Q(:,5) = UV_TO_Q ( [0.0_dp, 1.0_dp,  0.0_dp] )
        Q(:,6) = UV_TO_Q ( [0.0_dp, 1.0_dp,  0.0_dp] )
    ENDIF

END SUBROUTINE KAGOME

!     ==============================================================================================

SUBROUTINE FDDD(A1, A2, A3)
    
    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL

    IMPLICIT NONE
    REAL(KIND=DP) :: IDEAL_CELL_VLM, CELL_VLM, SCALE, A1, A2, A3, CN1, CN2, CN3

    !     Unit cell length
    IF (DENSITYT) THEN
        CELL_VLM = 24.0_dp / RHO
    ELSE
        CELL_VLM = PI*24.0_dp / (6.0_dp * RHO)
    ENDIF
    
    A1             = 2.0_dp
    A2             = SQRT(12.0_dp)
    A3             = 2.0_dp*A2
    IDEAL_CELL_VLM = A1*A2*A3
    SCALE          = (CELL_VLM/IDEAL_CELL_VLM)**(1.0_dp/3.0_dp)
    A1             = A1*SCALE
    A2             = A2*SCALE
    A3             = A3*SCALE

    BOX = [ A1*REAL(NUCX,DP), A2*REAL(NUCY,DP), A3*REAL(NUCZ,DP) ]
    PRINT *, "FDDD UNIT CELL:", BOX

    R(:,1)  = [A1*0.125_dp, A2*0.125_dp, A3*0.125_dp]
    R(:,2)  = [A1*0.125_dp, A2*0.875_dp, A3*0.875_dp]
    R(:,3)  = [A1*0.000_dp, A2*0.000_dp, A3*0.000_dp]
    R(:,4)  = [A1*0.875_dp, A2*0.875_dp, A3*0.125_dp]
    R(:,5)  = [A1*0.750_dp, A2*0.750_dp, A3*0.250_dp]
    R(:,6)  = [A1*0.875_dp, A2*0.125_dp, A3*0.875_dp]
    R(:,7)  = [A1*0.125_dp, A2*0.625_dp, A3*0.625_dp]
    R(:,8)  = [A1*0.125_dp, A2*0.375_dp, A3*0.375_dp]
    R(:,9)  = [A1*0.000_dp, A2*0.500_dp, A3*0.500_dp]
    R(:,10) = [A1*0.875_dp, A2*0.375_dp, A3*0.625_dp]
    R(:,11) = [A1*0.750_dp, A2*0.250_dp, A3*0.750_dp]
    R(:,12) = [A1*0.875_dp, A2*0.625_dp, A3*0.375_dp]
    R(:,13) = [A1*0.625_dp, A2*0.125_dp, A3*0.625_dp]
    R(:,14) = [A1*0.625_dp, A2*0.875_dp, A3*0.375_dp]
    R(:,15) = [A1*0.500_dp, A2*0.000_dp, A3*0.500_dp]
    R(:,16) = [A1*0.375_dp, A2*0.875_dp, A3*0.625_dp]
    R(:,17) = [A1*0.250_dp, A2*0.750_dp, A3*0.750_dp]
    R(:,18) = [A1*0.375_dp, A2*0.125_dp, A3*0.375_dp]
    R(:,19) = [A1*0.625_dp, A2*0.625_dp, A3*0.125_dp]
    R(:,20) = [A1*0.625_dp, A2*0.375_dp, A3*0.875_dp]
    R(:,21) = [A1*0.500_dp, A2*0.500_dp, A3*0.000_dp]
    R(:,22) = [A1*0.375_dp, A2*0.375_dp, A3*0.125_dp]
    R(:,23) = [A1*0.250_dp, A2*0.250_dp, A3*0.250_dp]
    R(:,24) = [A1*0.375_dp, A2*0.625_dp, A3*0.875_dp]

    IF(RIGIDT) THEN
        CN1 = 0.43922692270767_dp
        CN2 = 0.74637770845439_dp
        CN3 = 0.5_dp
        Q(:,1)  = UV_TO_Q ( [ CN1,     CN2,    0.5_dp] )
        Q(:,2)  = UV_TO_Q ( [-CN1,     CN2,    0.5_dp] )
        Q(:,3)  = UV_TO_Q ( [ 0.0_dp, 0.0_dp, -1.0_dp] )
        Q(:,4)  = UV_TO_Q ( [-CN1,    -CN2,    0.5_dp] )
        Q(:,5)  = UV_TO_Q ( [ 0.0_dp, 0.0_dp, -1.0_dp] )
        Q(:,6)  = UV_TO_Q ( [ CN1,    -CN2,    0.5_dp] )
        Q(:,7)  = UV_TO_Q ( [ CN1,     CN2,    0.5_dp] )
        Q(:,8)  = UV_TO_Q ( [-CN1,     CN2,    0.5_dp] )
        Q(:,9)  = UV_TO_Q ( [ 0.0_dp, 0.0_dp, -1.0_dp] )
        Q(:,10) = UV_TO_Q ( [-CN1,    -CN2,    0.5_dp] )
        Q(:,11) = UV_TO_Q ( [ 0.0_dp, 0.0_dp, -1.0_dp] )
        Q(:,12) = UV_TO_Q ( [ CN1,    -CN2,    0.5_dp] )
        Q(:,13) = UV_TO_Q ( [ CN1,     CN2,    0.5_dp] )
        Q(:,14) = UV_TO_Q ( [-CN1,     CN2,    0.5_dp] )
        Q(:,15) = UV_TO_Q ( [ 0.0_dp, 0.0_dp, -1.0_dp] )
        Q(:,16) = UV_TO_Q ( [-CN1,    -CN2,    0.5_dp] )
        Q(:,17) = UV_TO_Q ( [ 0.0_dp, 0.0_dp, -1.0_dp] )
        Q(:,18) = UV_TO_Q ( [ CN1,    -CN2,    0.5_dp] )
        Q(:,19) = UV_TO_Q ( [ CN1,     CN2,    0.5_dp] )
        Q(:,20) = UV_TO_Q ( [-CN1,     CN2,    0.5_dp] )
        Q(:,21) = UV_TO_Q ( [ 0.0_dp, 0.0_dp, -1.0_dp] )
        Q(:,22) = UV_TO_Q ( [-CN1,    -CN2,    0.5_dp] )
        Q(:,23) = UV_TO_Q ( [ 0.0_dp, 0.0_dp, -1.0_dp] )
        Q(:,24) = UV_TO_Q ( [ CN1,    -CN2,    0.5_dp] )
    ENDIF

END SUBROUTINE FDDD

SUBROUTINE OCTAHEDRAL_SC(UCL)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL

    IMPLICIT NONE
    REAL(KIND=DP) :: UCL, C1, C2, C3, MAXRHO

    MAXRHO = 6.0_dp / (1.0_dp+SQRT(2.0_dp))**3
    IF (DENSITYT) THEN
        UCL     = (6.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL  = (6.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
    PRINT *, "Octahedral-Simple Cubic box:", BOX, MAXRHO
    
    C1 = ( 1.0_dp+SQRT(2.0_dp) ) / 2.0_dp * 1.0000000001_dp
    C2 = ( 1.0_dp+2.0_dp*SQRT(2.0_dp) ) / 2.0_dp * 1.0000000001_dp
    C3 = 0.5_dp*1.0000000001_dp
    
!   Build the unit cell
    R(:,1)  = [C3, C1, C1]
    R(:,2)  = [C1, C3, C1]
    R(:,3)  = [C1, C1, C3]
    R(:,4)  = [C1, C1, C2]
    R(:,5)  = [C2, C1, C1]
    R(:,6)  = [C1, C2, C1]

    IF(RIGIDT) THEN
        Q(:,1)  = UV_TO_Q ( [-1.0_dp,  0.0_dp, 0.0_dp] )
        Q(:,2)  = UV_TO_Q ( [ 0.0_dp, -1.0_dp, 0.0_dp] )
        Q(:,3)  = UV_TO_Q ( [ 0.0_dp,  0.0_dp,-1.0_dp] )
        Q(:,4)  = UV_TO_Q ( [ 0.0_dp,  0.0_dp, 1.0_dp] )
        Q(:,5)  = UV_TO_Q ( [ 1.0_dp,  0.0_dp, 0.0_dp] )
        Q(:,6)  = UV_TO_Q ( [ 0.0_dp,  1.0_dp, 0.0_dp] )
    ENDIF

END SUBROUTINE OCTAHEDRAL_SC

SUBROUTINE OCTAHEDRAL_BCC(UCL)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL

    IMPLICIT NONE
    REAL(KIND=DP) :: UCL, C1, C2, C3

!   Unit cell length (remember there are 24 particles in the unit cell)
    IF (DENSITYT) THEN
        UCL     = (12.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL  = (12.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
    PRINT *, "Octahedral-BCC box:", BOX
    
    C1 = ( 1.0_dp+SQRT(2.0_dp) ) / 2.0_dp * 1.0000000001_dp
    C2 = ( 1.0_dp+2.0_dp*SQRT(2.0_dp) ) / 2.0_dp * 1.0000000001_dp
    C3 = 0.5_dp*1.0000000001_dp
    
!   Build the unit cell
    R(:,1)   = [C3, C1, C1]
    R(:,2)   = [C1, C3, C1]
    R(:,3)   = [C1, C1, C3]
    R(:,4)   = [C1, C1, C2]
    R(:,5)   = [C2, C1, C1]
    R(:,6)   = [C1, C2, C1]
    R(:,7)   = [C3, C1, C1] + 0.5_dp*UCL
    R(:,8)   = [C1, C3, C1] + 0.5_dp*UCL
    R(:,9)   = [C1, C1, C3] + 0.5_dp*UCL
    R(:,10)  = [C1, C1, C2] + 0.5_dp*UCL
    R(:,11)  = [C2, C1, C1] + 0.5_dp*UCL
    R(:,12)  = [C1, C2, C1] + 0.5_dp*UCL

    IF(RIGIDT) THEN
        Q(:,1)  = UV_TO_Q (  [-1.0_dp,  0.0_dp, 0.0_dp] )
        Q(:,2)  = UV_TO_Q (  [ 0.0_dp, -1.0_dp, 0.0_dp] )
        Q(:,3)  = UV_TO_Q (  [ 0.0_dp,  0.0_dp,-1.0_dp] )
        Q(:,4)  = UV_TO_Q (  [ 0.0_dp,  0.0_dp, 1.0_dp] )
        Q(:,5)  = UV_TO_Q (  [ 1.0_dp,  0.0_dp, 0.0_dp] )
        Q(:,6)  = UV_TO_Q (  [ 0.0_dp,  1.0_dp, 0.0_dp] )
        Q(:,7)  = UV_TO_Q (  [-1.0_dp,  0.0_dp, 0.0_dp] )
        Q(:,8)  = UV_TO_Q (  [ 0.0_dp, -1.0_dp, 0.0_dp] )
        Q(:,9)  = UV_TO_Q (  [ 0.0_dp,  0.0_dp,-1.0_dp] )
        Q(:,10)  = UV_TO_Q ( [ 0.0_dp,  0.0_dp, 1.0_dp] )
        Q(:,11)  = UV_TO_Q ( [ 1.0_dp,  0.0_dp, 0.0_dp] )
        Q(:,12)  = UV_TO_Q ( [ 0.0_dp,  1.0_dp, 0.0_dp] )
    ENDIF

END SUBROUTINE OCTAHEDRAL_BCC

SUBROUTINE CLATHRATE_I(UCL,COLLMOLT)

USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, CPPT, CPPLAM
USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL, Q_TO_RM

IMPLICIT NONE

INTEGER       :: J1, J2, INDX
REAL(KIND=DP) :: UCL, RCLU(3,46), QCLU(4,46), CN0, BL, RM(3,3), RSITE(3,4), USITE(3)
LOGICAL       :: COLLMOLT

IF(COLLMOLT) THEN
    IF (DENSITYT) THEN
        UCL     = (184.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL  = (184.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
    PRINT *, "Clathrate sI box:", BOX
ELSE
    !   Unit cell length (remember there are 24 particles in the unit cell)
    IF (DENSITYT) THEN
        UCL     = (46.0_dp/RHO)**(1.0_dp/3.0_dp)
    ELSE
        UCL  = (46.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
    ENDIF

    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
    PRINT *, "Clathrate sI box:", BOX
ENDIF

!   Build the unit cell
R(:,1) = [0.25000_dp,0.50000_dp,0.00000_dp]*UCL;R(:,24)=[0.00000_dp,0.88328_dp,0.69114_dp]*UCL 
R(:,2) = [0.00000_dp,0.25000_dp,0.50000_dp]*UCL;R(:,25)=[0.30886_dp,0.00000_dp,0.11672_dp]*UCL 
R(:,3) = [0.50000_dp,0.00000_dp,0.25000_dp]*UCL;R(:,26)=[0.69114_dp,0.00000_dp,0.88328_dp]*UCL 
R(:,4) = [0.75000_dp,0.50000_dp,0.00000_dp]*UCL;R(:,27)=[0.11672_dp,0.30886_dp,0.00000_dp]*UCL 
R(:,5) = [0.00000_dp,0.75000_dp,0.50000_dp]*UCL;R(:,28)=[0.88328_dp,0.69114_dp,0.00000_dp]*UCL 
R(:,6) = [0.50000_dp,0.00000_dp,0.75000_dp]*UCL;R(:,29)=[0.00000_dp,0.11672_dp,0.69114_dp]*UCL 
R(:,7) = [0.18309_dp,0.18309_dp,0.18309_dp]*UCL;R(:,30)=[0.00000_dp,0.88328_dp,0.30886_dp]*UCL 
R(:,8) = [0.81691_dp,0.81691_dp,0.81691_dp]*UCL;R(:,31)=[0.69114_dp,0.00000_dp,0.11672_dp]*UCL 
R(:,9) = [0.18309_dp,0.81691_dp,0.81691_dp]*UCL;R(:,32)=[0.30886_dp,0.00000_dp,0.88328_dp]*UCL 
R(:,10)= [0.81691_dp,0.18309_dp,0.81691_dp]*UCL;R(:,33)=[0.11672_dp,0.69114_dp,0.00000_dp]*UCL 
R(:,11)= [0.81691_dp,0.81691_dp,0.18309_dp]*UCL;R(:,34)=[0.88328_dp,0.30886_dp,0.00000_dp]*UCL 
R(:,12)= [0.81691_dp,0.18309_dp,0.18309_dp]*UCL;R(:,35)=[0.50000_dp,0.80886_dp,0.61672_dp]*UCL 
R(:,13)= [0.18309_dp,0.81691_dp,0.18309_dp]*UCL;R(:,36)=[0.50000_dp,0.19114_dp,0.38328_dp]*UCL 
R(:,14)= [0.18309_dp,0.18309_dp,0.81691_dp]*UCL;R(:,37)=[0.50000_dp,0.19114_dp,0.61672_dp]*UCL 
R(:,15)= [0.68309_dp,0.68309_dp,0.68309_dp]*UCL;R(:,38)=[0.50000_dp,0.80886_dp,0.38328_dp]*UCL 
R(:,16)= [0.31691_dp,0.31691_dp,0.31691_dp]*UCL;R(:,39)=[0.61672_dp,0.50000_dp,0.80886_dp]*UCL 
R(:,17)= [0.68309_dp,0.31691_dp,0.31691_dp]*UCL;R(:,40)=[0.38328_dp,0.50000_dp,0.19114_dp]*UCL 
R(:,18)= [0.31691_dp,0.68309_dp,0.31691_dp]*UCL;R(:,41)=[0.61672_dp,0.50000_dp,0.19114_dp]*UCL 
R(:,19)= [0.31691_dp,0.31691_dp,0.68309_dp]*UCL;R(:,42)=[0.38328_dp,0.50000_dp,0.80886_dp]*UCL 
R(:,20)= [0.31691_dp,0.68309_dp,0.68309_dp]*UCL;R(:,43)=[0.80886_dp,0.61672_dp,0.50000_dp]*UCL 
R(:,21)= [0.68309_dp,0.31691_dp,0.68309_dp]*UCL;R(:,44)=[0.19114_dp,0.38328_dp,0.50000_dp]*UCL 
R(:,22)= [0.68309_dp,0.68309_dp,0.31691_dp]*UCL;R(:,45)=[0.19114_dp,0.61672_dp,0.50000_dp]*UCL 
R(:,23)= [0.00000_dp,0.11672_dp,0.30886_dp]*UCL;R(:,46)=[0.80886_dp,0.38328_dp,0.50000_dp]*UCL 

IF(RIGIDT) THEN
Q(:,1) = [ 0.6532815_dp,-0.2705980_dp,-0.6532815_dp, 0.2705980_dp];Q(:,24)=[ 0.3717328_dp, 0.0908920_dp,-0.2194327_dp, 0.8974423_dp] 
Q(:,2) = [ 0.6532815_dp,-0.2705981_dp,-0.2705981_dp, 0.6532814_dp];Q(:,25)=[ 0.5703172_dp,-0.1076924_dp, 0.4180172_dp, 0.6988579_dp] 
Q(:,3) = [ 0.9238795_dp, 0.0000000_dp, 0.0000000_dp, 0.3826834_dp];Q(:,26)=[ 0.6988579_dp, 0.4180172_dp,-0.1076924_dp, 0.5703171_dp] 
Q(:,4) = [ 0.6532815_dp, 0.2705981_dp, 0.2705980_dp,-0.6532815_dp];Q(:,27)=[ 0.4794251_dp,-0.1985844_dp,-0.3271251_dp, 0.7897499_dp] 
Q(:,5) = [ 0.6532815_dp, 0.2705981_dp, 0.2705980_dp, 0.6532814_dp];Q(:,28)=[ 0.8974423_dp, 0.0908920_dp,-0.3717328_dp, 0.2194328_dp] 
Q(:,6) = [ 0.3826835_dp, 0.0000000_dp, 0.0000000_dp, 0.9238795_dp];Q(:,29)=[-0.0908920_dp,-0.3717327_dp,-0.8974423_dp, 0.2194328_dp] 
Q(:,7) = [ 0.6558941_dp, 0.4358146_dp, 0.4358145_dp, 0.4358145_dp];Q(:,30)=[ 0.3271251_dp, 0.1985844_dp, 0.7897499_dp, 0.4794251_dp] 
Q(:,8) = [ 0.7719546_dp,-0.1556197_dp, 0.0000000_dp, 0.6163348_dp];Q(:,31)=[ 0.6988579_dp, 0.4180172_dp, 0.1076924_dp,-0.5703172_dp] 
Q(:,9) = [ 0.4358145_dp,-0.4358146_dp,-0.4358145_dp, 0.6558941_dp];Q(:,32)=[-0.1985844_dp,-0.7897499_dp, 0.3271251_dp, 0.4794251_dp] 
Q(:,10)= [ 0.1100398_dp,-0.1100398_dp, 0.9816689_dp, 0.1100397_dp];Q(:,33)=[ 0.3717327_dp, 0.2194328_dp, 0.8974423_dp,-0.0908920_dp] 
Q(:,11)= [ 0.9816689_dp, 0.1100397_dp, 0.1100398_dp,-0.1100398_dp];Q(:,34)=[ 0.4180172_dp,-0.1076924_dp, 0.6988579_dp, 0.5703172_dp] 
Q(:,12)= [ 0.6163349_dp, 0.0000000_dp, 0.1556197_dp, 0.7719546_dp];Q(:,35)=[ 0.3717327_dp,-0.0908920_dp, 0.8974423_dp,-0.2194328_dp] 
Q(:,13)= [ 0.7719546_dp, 0.1556197_dp, 0.0000000_dp,-0.6163348_dp];Q(:,36)=[ 0.4794251_dp,-0.7897500_dp,-0.3271251_dp,-0.1985843_dp] 
Q(:,14)= [ 0.1556197_dp,-0.6163348_dp, 0.0000000_dp, 0.7719546_dp];Q(:,37)=[-0.8974423_dp,-0.2194327_dp,-0.3717328_dp,-0.0908921_dp] 
Q(:,15)= [ 0.3257747_dp, 0.5458543_dp, 0.5458544_dp, 0.5458543_dp];Q(:,38)=[ 0.0908920_dp,-0.3717327_dp,-0.2194327_dp, 0.8974423_dp] 
Q(:,16)= [ 0.7719546_dp, 0.1556197_dp,-0.6163348_dp, 0.0000000_dp];Q(:,39)=[ 0.3717328_dp,-0.2194328_dp,-0.0908920_dp, 0.8974423_dp] 
Q(:,17)= [ 0.6558940_dp,-0.4358146_dp, 0.4358146_dp, 0.4358146_dp];Q(:,40)=[ 0.6988579_dp,-0.5703172_dp, 0.1076923_dp,-0.4180171_dp] 
Q(:,18)= [ 0.5458542_dp, 0.5458543_dp, 0.5458543_dp, 0.3257748_dp];Q(:,41)=[-0.4794251_dp,-0.1985844_dp,-0.7897499_dp, 0.3271251_dp] 
Q(:,19)= [-0.9816689_dp, 0.1100399_dp, 0.1100383_dp,-0.1100409_dp];Q(:,42)=[ 0.1076924_dp,-0.4180172_dp, 0.6988579_dp,-0.5703172_dp] 
Q(:,20)= [ 0.1556197_dp,-0.7719546_dp, 0.0000000_dp, 0.6163349_dp];Q(:,43)=[ 0.4180172_dp,-0.6988579_dp,-0.1076923_dp, 0.5703171_dp] 
Q(:,21)= [ 0.7719546_dp, 0.0000000_dp, 0.1556197_dp, 0.6163348_dp];Q(:,44)=[ 0.5703172_dp,-0.1076923_dp, 0.6988579_dp,-0.4180172_dp] 
Q(:,22)= [-0.1556197_dp,-0.7719545_dp, 0.0000000_dp, 0.6163349_dp];Q(:,45)=[ 0.1076924_dp, 0.5703172_dp, 0.4180172_dp, 0.6988579_dp] 
Q(:,23)= [ 0.2194327_dp,-0.8974423_dp,-0.3717328_dp,-0.0908920_dp];Q(:,46)=[ 0.6988579_dp, 0.4180171_dp,-0.5703172_dp,-0.1076924_dp]
ENDIF

IF(COLLMOLT) THEN
    RCLU = R(:,1:46) 
    QCLU = Q(:,1:46)
    R = 0.D0
    Q = 0.D0

    RSITE(:,1)= [  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
    RSITE(:,2)= [ -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
    RSITE(:,3)= [  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]
    RSITE(:,4)= [ -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]

    CN0 = 1.0_dp / SQRT(3.0_dp)
    RM  = 0.0_dp
    IF(CPPT) THEN
        BL = CPPLAM
    ELSE
        BL = 1.01_dp
    ENDIF

    DO J1 = 1, 46
        DO J2 = 1,4
            INDX = 4*(J1-1) + J2
            USITE = MATMUL(Q_TO_RM(QCLU(:,J1)), RSITE(:,J2))
            R(:,INDX) = RCLU(:,J1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*USITE
            IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( USITE )
        ENDDO
    ENDDO

    DO J1 = 1, 184
        R(:,J1) = R(:,J1) + BL*(SQRT(3.0_dp) / (2.0_dp*SQRT(2.0_dp)) )
    ENDDO
ENDIF

END SUBROUTINE CLATHRATE_I

SUBROUTINE CLATHRATE_II(UCL,COLLMOLT)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, CPPT, CPPLAM
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL, Q_TO_RM
    
    IMPLICIT NONE
    
    INTEGER       :: J1, J2, INDX
    REAL(KIND=DP) :: UCL, RCLU(3,136), QCLU(4,136), CN0, BL, RM(3,3), RSITE(3,4), USITE(3)
    REAL(KIND=DP) :: C0, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15
    LOGICAL       :: COLLMOLT
    
    IF(COLLMOLT) THEN
        IF (DENSITYT) THEN
            UCL     = (544.0_dp/RHO)**(1.0_dp/3.0_dp)
        ELSE
            UCL  = (544.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
        ENDIF
    
        PRINT *, "Tetrahedral Clathrate sII box:", BOX
    ELSE
        !   Unit cell length (remember there are 24 particles in the unit cell)
        IF (DENSITYT) THEN
            UCL     = (136.0_dp/RHO)**(1.0_dp/3.0_dp)
        ELSE
            UCL  = (136.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
        ENDIF

        PRINT *, "Clathrate sII box:", BOX
    ENDIF

    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]

    R(:,1)  =[0.19226_dp,0.49513_dp,0.69226_dp]*UCL;R(:,2)  =[0.80774_dp,0.99513_dp,0.80774_dp]*UCL;
    R(:,3)  =[0.94226_dp,0.05774_dp,0.74513_dp]*UCL;R(:,4)  =[0.05774_dp,0.25487_dp,0.05774_dp]*UCL
    R(:,5)  =[0.05774_dp,0.55774_dp,0.75487_dp]*UCL;R(:,6)  =[0.99513_dp,0.80774_dp,0.80774_dp]*UCL
    R(:,7)  =[0.44226_dp,0.24513_dp,0.05774_dp]*UCL;R(:,8)  =[0.50000_dp,0.00000_dp,0.00000_dp]*UCL
    R(:,9)  =[0.94226_dp,0.25487_dp,0.94226_dp]*UCL;R(:,10) =[0.55774_dp,0.24513_dp,0.94226_dp]*UCL
    R(:,11) =[0.84247_dp,0.84247_dp,0.84247_dp]*UCL;R(:,12) =[0.75487_dp,0.55774_dp,0.05774_dp]*UCL
    R(:,13) =[0.74513_dp,0.05774_dp,0.94226_dp]*UCL;R(:,14) =[0.75487_dp,0.44226_dp,0.94226_dp]*UCL
    R(:,15) =[0.90753_dp,0.40753_dp,0.90753_dp]*UCL;R(:,16) =[0.34247_dp,0.15753_dp,0.65753_dp]*UCL
    R(:,17) =[0.74513_dp,0.94226_dp,0.05774_dp]*UCL;R(:,18) =[0.30774_dp,0.00487_dp,0.69226_dp]*UCL
    R(:,19) =[0.80774_dp,0.19226_dp,0.00487_dp]*UCL;R(:,20) =[0.59247_dp,0.09247_dp,0.90753_dp]*UCL
    R(:,21) =[0.94226_dp,0.44226_dp,0.75487_dp]*UCL;R(:,22) =[0.00487_dp,0.19226_dp,0.80774_dp]*UCL
    R(:,23) =[0.69226_dp,0.30774_dp,0.00487_dp]*UCL;R(:,24) =[0.15753_dp,0.15753_dp,0.84247_dp]*UCL
    R(:,25) =[0.19226_dp,0.00487_dp,0.80774_dp]*UCL;R(:,26) =[0.49513_dp,0.19226_dp,0.69226_dp]*UCL
    R(:,27) =[0.19226_dp,0.19226_dp,0.99513_dp]*UCL;R(:,28) =[0.00487_dp,0.30774_dp,0.69226_dp]*UCL
    R(:,29) =[0.25000_dp,0.25000_dp,0.75000_dp]*UCL;R(:,30) =[0.80774_dp,0.80774_dp,0.99513_dp]*UCL
    R(:,31) =[0.59247_dp,0.90753_dp,0.09247_dp]*UCL;R(:,32) =[0.15753_dp,0.34247_dp,0.65753_dp]*UCL
    R(:,33) =[0.05774_dp,0.94226_dp,0.74513_dp]*UCL;R(:,34) =[0.40753_dp,0.09247_dp,0.09247_dp]*UCL
    R(:,35) =[0.19226_dp,0.99513_dp,0.19226_dp]*UCL;R(:,36) =[0.80774_dp,0.49513_dp,0.30774_dp]*UCL
    R(:,37) =[0.94226_dp,0.55774_dp,0.24513_dp]*UCL;R(:,38) =[0.05774_dp,0.75487_dp,0.55774_dp]*UCL
    R(:,39) =[0.05774_dp,0.05774_dp,0.25487_dp]*UCL;R(:,40) =[0.99513_dp,0.30774_dp,0.30774_dp]*UCL
    R(:,41) =[0.44226_dp,0.74513_dp,0.55774_dp]*UCL;R(:,42) =[0.50000_dp,0.50000_dp,0.50000_dp]*UCL
    R(:,43) =[0.94226_dp,0.75487_dp,0.44226_dp]*UCL;R(:,44) =[0.55774_dp,0.74513_dp,0.44226_dp]*UCL
    R(:,45) =[0.84247_dp,0.34247_dp,0.34247_dp]*UCL;R(:,46) =[0.75487_dp,0.05774_dp,0.55774_dp]*UCL
    R(:,47) =[0.74513_dp,0.55774_dp,0.44226_dp]*UCL;R(:,48) =[0.75487_dp,0.94226_dp,0.44226_dp]*UCL
    R(:,49) =[0.90753_dp,0.90753_dp,0.40753_dp]*UCL;R(:,50) =[0.34247_dp,0.65753_dp,0.15753_dp]*UCL
    R(:,51) =[0.74513_dp,0.44226_dp,0.55774_dp]*UCL;R(:,52) =[0.30774_dp,0.50487_dp,0.19226_dp]*UCL
    R(:,53) =[0.80774_dp,0.69226_dp,0.50487_dp]*UCL;R(:,54) =[0.59247_dp,0.59247_dp,0.40753_dp]*UCL
    R(:,55) =[0.94226_dp,0.94226_dp,0.25487_dp]*UCL;R(:,56) =[0.00487_dp,0.69226_dp,0.30774_dp]*UCL
    R(:,57) =[0.69226_dp,0.80774_dp,0.50487_dp]*UCL;R(:,58) =[0.15753_dp,0.65753_dp,0.34247_dp]*UCL
    R(:,59) =[0.19226_dp,0.50487_dp,0.30774_dp]*UCL;R(:,60) =[0.49513_dp,0.69226_dp,0.19226_dp]*UCL
    R(:,61) =[0.19226_dp,0.69226_dp,0.49513_dp]*UCL;R(:,62) =[0.00487_dp,0.80774_dp,0.19226_dp]*UCL
    R(:,63) =[0.25000_dp,0.75000_dp,0.25000_dp]*UCL;R(:,64) =[0.80774_dp,0.30774_dp,0.49513_dp]*UCL
    R(:,65) =[0.59247_dp,0.40753_dp,0.59247_dp]*UCL;R(:,66) =[0.15753_dp,0.84247_dp,0.15753_dp]*UCL
    R(:,67) =[0.05774_dp,0.44226_dp,0.24513_dp]*UCL;R(:,68) =[0.40753_dp,0.59247_dp,0.59247_dp]*UCL
    R(:,69) =[0.69226_dp,0.49513_dp,0.19226_dp]*UCL;R(:,70) =[0.30774_dp,0.99513_dp,0.30774_dp]*UCL
    R(:,71) =[0.44226_dp,0.05774_dp,0.24513_dp]*UCL;R(:,72) =[0.55774_dp,0.25487_dp,0.55774_dp]*UCL
    R(:,73) =[0.55774_dp,0.55774_dp,0.25487_dp]*UCL;R(:,74) =[0.49513_dp,0.80774_dp,0.30774_dp]*UCL
    R(:,75) =[0.94226_dp,0.24513_dp,0.55774_dp]*UCL;R(:,76) =[0.00000_dp,0.00000_dp,0.50000_dp]*UCL
    R(:,77) =[0.44226_dp,0.25487_dp,0.44226_dp]*UCL;R(:,78) =[0.05774_dp,0.24513_dp,0.44226_dp]*UCL
    R(:,79) =[0.34247_dp,0.84247_dp,0.34247_dp]*UCL;R(:,80) =[0.25487_dp,0.55774_dp,0.55774_dp]*UCL
    R(:,81) =[0.24513_dp,0.05774_dp,0.44226_dp]*UCL;R(:,82) =[0.25487_dp,0.44226_dp,0.44226_dp]*UCL
    R(:,83) =[0.40753_dp,0.40753_dp,0.40753_dp]*UCL;R(:,84) =[0.84247_dp,0.15753_dp,0.15753_dp]*UCL
    R(:,85) =[0.24513_dp,0.94226_dp,0.55774_dp]*UCL;R(:,86) =[0.80774_dp,0.00487_dp,0.19226_dp]*UCL
    R(:,87) =[0.30774_dp,0.19226_dp,0.50487_dp]*UCL;R(:,88) =[0.09247_dp,0.09247_dp,0.40753_dp]*UCL
    R(:,89) =[0.44226_dp,0.44226_dp,0.25487_dp]*UCL;R(:,90) =[0.50487_dp,0.19226_dp,0.30774_dp]*UCL
    R(:,91) =[0.19226_dp,0.30774_dp,0.50487_dp]*UCL;R(:,92) =[0.65753_dp,0.15753_dp,0.34247_dp]*UCL
    R(:,93) =[0.69226_dp,0.00487_dp,0.30774_dp]*UCL;R(:,94) =[0.99513_dp,0.19226_dp,0.19226_dp]*UCL
    R(:,95) =[0.69226_dp,0.19226_dp,0.49513_dp]*UCL;R(:,96) =[0.50487_dp,0.30774_dp,0.19226_dp]*UCL
    R(:,97) =[0.75000_dp,0.25000_dp,0.25000_dp]*UCL;R(:,98) =[0.30774_dp,0.80774_dp,0.49513_dp]*UCL
    R(:,99) =[0.09247_dp,0.90753_dp,0.59247_dp]*UCL;R(:,100)=[0.65753_dp,0.34247_dp,0.15753_dp]*UCL
    R(:,101)=[0.55774_dp,0.94226_dp,0.24513_dp]*UCL;R(:,102)=[0.90753_dp,0.09247_dp,0.59247_dp]*UCL
    R(:,103)=[0.69226_dp,0.99513_dp,0.69226_dp]*UCL;R(:,104)=[0.30774_dp,0.49513_dp,0.80774_dp]*UCL
    R(:,105)=[0.44226_dp,0.55774_dp,0.74513_dp]*UCL;R(:,106)=[0.55774_dp,0.75487_dp,0.05774_dp]*UCL
    R(:,107)=[0.55774_dp,0.05774_dp,0.75487_dp]*UCL;R(:,108)=[0.49513_dp,0.30774_dp,0.80774_dp]*UCL
    R(:,109)=[0.94226_dp,0.74513_dp,0.05774_dp]*UCL;R(:,110)=[0.00000_dp,0.50000_dp,0.00000_dp]*UCL
    R(:,111)=[0.44226_dp,0.75487_dp,0.94226_dp]*UCL;R(:,112)=[0.05774_dp,0.74513_dp,0.94226_dp]*UCL
    R(:,113)=[0.34247_dp,0.34247_dp,0.84247_dp]*UCL;R(:,114)=[0.25487_dp,0.05774_dp,0.05774_dp]*UCL
    R(:,115)=[0.24513_dp,0.55774_dp,0.94226_dp]*UCL;R(:,116)=[0.25487_dp,0.94226_dp,0.94226_dp]*UCL
    R(:,117)=[0.40753_dp,0.90753_dp,0.90753_dp]*UCL;R(:,118)=[0.84247_dp,0.65753_dp,0.65753_dp]*UCL
    R(:,119)=[0.24513_dp,0.44226_dp,0.05774_dp]*UCL;R(:,120)=[0.80774_dp,0.50487_dp,0.69226_dp]*UCL
    R(:,121)=[0.30774_dp,0.69226_dp,0.00487_dp]*UCL;R(:,122)=[0.09247_dp,0.59247_dp,0.90753_dp]*UCL
    R(:,123)=[0.44226_dp,0.94226_dp,0.75487_dp]*UCL;R(:,124)=[0.50487_dp,0.69226_dp,0.80774_dp]*UCL
    R(:,125)=[0.19226_dp,0.80774_dp,0.00487_dp]*UCL;R(:,126)=[0.65753_dp,0.65753_dp,0.84247_dp]*UCL
    R(:,127)=[0.69226_dp,0.50487_dp,0.80774_dp]*UCL;R(:,128)=[0.99513_dp,0.69226_dp,0.69226_dp]*UCL
    R(:,129)=[0.69226_dp,0.69226_dp,0.99513_dp]*UCL;R(:,130)=[0.50487_dp,0.80774_dp,0.69226_dp]*UCL
    R(:,131)=[0.75000_dp,0.75000_dp,0.75000_dp]*UCL;R(:,132)=[0.30774_dp,0.30774_dp,0.99513_dp]*UCL
    R(:,133)=[0.09247_dp,0.40753_dp,0.09247_dp]*UCL;R(:,134)=[0.65753_dp,0.84247_dp,0.65753_dp]*UCL
    R(:,135)=[0.55774_dp,0.44226_dp,0.74513_dp]*UCL;R(:,136)=[0.90753_dp,0.59247_dp,0.09247_dp]*UCL
    IF(RIGIDT) THEN
        C0=0.0_dp;C15=1.0_dp;C1=1.0_dp/2.0_dp;C2=1.0_dp/6.0_dp;C3=5.0_dp/6.0_dp;C4=1.0_dp/(3.0_dp*SQRT(2.0_dp))
        C5=4.0_dp/(3.0_dp*SQRT(2.0_dp));C6=SQRT(2.0_dp)/3.0_dp;C7=1.0_dp/SQRT(2.0_dp);C8=SQRT(1.0_dp/6.0_dp);C9=SQRT(2.0_dp/3.0_dp)
        C10=1.0_dp/SQRT(3.0_dp);C11=1.0_dp/3.0_dp;C12=2.0_dp/3.0_dp;C13=SQRT(3.0_dp)/6.0_dp;C14=SQRT(3.0_dp)/2.0_dp
        Q(:,1)  =[-C1 , C3 , C2 ,-C2 ];Q(:,2)  =[ C2 ,-C2 ,-C3 , C1 ];Q(:,3)  =[ C0 , C4 ,-C4 ,-C5 ];Q(:,4)  =[ C6 ,-C6 ,-C7 ,-C4 ]    
        Q(:,5)  =[ C4 , C7 ,-C6 ,-C6 ];Q(:,6)  =[ C1 , C2 ,-C3 , C2 ];Q(:,7)  =[ C7 ,-C6 , C6 ,-C4 ];Q(:,8)  =[-C7 , C7 , C0 , C0 ]    
        Q(:,9)  =[ C4 , C6 ,-C6 ,-C7 ];Q(:,10) =[ C5 ,-C4 ,-C0 ,-C4 ];Q(:,11) =[ C8 ,-C9 , C0 ,-C8 ];Q(:,12) =[-C6 , C4 ,-C7 ,-C6 ]    
        Q(:,13) =[-C7 ,-C6 , C4 , C6 ];Q(:,14) =[ C4 ,-C4 , C5 ,-C0 ];Q(:,15) =[ C0 , C10, C10, C10];Q(:,16) =[-C9 , C8 , C0 , C8 ]    
        Q(:,17) =[ C4 ,-C4 , C0 , C5 ];Q(:,18) =[ C12,-C0 , C12, C11];Q(:,19) =[ C11, C12, C12,-C0 ];Q(:,20) =[-C10, C10,-C0 , C10]     
        Q(:,21) =[ C7 , C4 ,-C6 ,-C6 ];Q(:,22) =[ C3 ,-C2 , C2 ,-C1 ];Q(:,23) =[-C2 ,-C2 , C3 ,-C1 ];Q(:,24) =[ C0 , C8 , C8 , C9 ]    
        Q(:,25) =[ C2 , C3 , C1 ,-C2 ];Q(:,26) =[-C1 , C2 ,-C2 , C3 ];Q(:,27) =[-C3 , C1 , C2 , C2 ];Q(:,28) =[-C2 ,-C3 ,-C1 ,-C2 ]    
        Q(:,29) =[-C1 , C1 ,-C1 ,-C1 ];Q(:,30) =[ C0 ,-C12,-C12, C11];Q(:,31) =[-C10, C10,-C0 ,-C10];Q(:,32) =[-C8 ,-C8 , C9 , C0 ]    
        Q(:,33) =[-C4 ,-C6 ,-C7 , C6 ];Q(:,34) =[ C13,-C14, C13,-C13];Q(:,35) =[-C2 ,-C2 , C3 , C1 ];Q(:,36) =[-C12,-C0 , C12, C11]     
        Q(:,37) =[ C6 ,-C7 , C6 , C4 ];Q(:,38) =[-C6 ,-C6 , C4 ,-C7 ];Q(:,39) =[ C4 , C0 , C5 , C4 ];Q(:,40) =[-C3 ,-C2 , C2 , C1 ]    
        Q(:,41) =[ C7 ,-C6 , C6 ,-C4 ];Q(:,42) =[-C7 ,-C7 , C0 , C0 ];Q(:,43) =[-C6 ,-C7 ,-C4 ,-C6 ];Q(:,44) =[ C5 ,-C4 ,-C0 ,-C4 ]    
        Q(:,45) =[-C8 ,-C0 , C8 , C9 ];Q(:,46) =[ C7 ,-C6 ,-C4 , C6 ];Q(:,47) =[-C6 , C7 ,-C6 , C4 ];Q(:,48) =[-C6 ,-C7 , C6 , C4 ]    
        Q(:,49) =[-C13,-C13, C14, C13];Q(:,50) =[ C8 ,-C0 ,-C9 ,-C8 ];Q(:,51) =[ C7 ,-C6 , C6 , C4 ];Q(:,52) =[ C2 ,-C1 ,-C3 ,-C2 ]    
        Q(:,53) =[ C2 , C3 ,-C2 ,-C1 ];Q(:,54) =[ C13,-C13,-C14,-C13];Q(:,55) =[-C7 ,-C6 , C4 ,-C6 ];Q(:,56) =[-C12,-C12,-C11, C0 ]    
        Q(:,57) =[-C2 , C2 , C1 , C3 ];Q(:,58) =[ C8 ,-C9 ,-C8 ,-C0 ];Q(:,59) =[ C1 ,-C3 , C2 ,-C2 ];Q(:,60) =[-C2 , C3 , C1 ,-C2 ]    
        Q(:,61) =[ C0 ,-C12,-C12,-C11];Q(:,62) =[-C11,-C0 ,-C12,-C12];Q(:,63) =[-C15, C0 , C0 , C0 ];Q(:,64) =[-C2 ,-C3 ,-C2 ,-C1 ]    
        Q(:,65) =[ C13,-C13, C14, C13];Q(:,66) =[ C8 , C0 , C9 ,-C8 ];Q(:,67) =[ C0 ,-C4 , C4 ,-C5 ];Q(:,68) =[ C13,-C13,-C13,-C14]     
        Q(:,69) =[ C1 ,-C2 , C2 , C3 ];Q(:,70) =[ C11,-C12,-C0 , C12];Q(:,71) =[-C6 , C6 ,-C4 , C7 ];Q(:,72) =[-C7 ,-C6 , C6 , C4 ]    
        Q(:,73) =[ C6 ,-C4 ,-C6 ,-C7 ];Q(:,74) =[-C1 ,-C2 , C3 ,-C2 ];Q(:,75) =[-C7 , C6 ,-C6 , C4 ];Q(:,76) =[ C0 , C0 ,-C7 ,-C7 ]    
        Q(:,77) =[ C7 ,-C4 , C6 , C6 ];Q(:,78) =[-C4 ,-C5 , C4 ,-C0 ];Q(:,79) =[ C0 ,-C8 ,-C8 , C9 ];Q(:,80) =[ C6 , C7 , C6 , C4 ]    
        Q(:,81) =[ C4 ,-C6 , C6 , C7 ];Q(:,82) =[ C5 ,-C0 ,-C4 , C4 ];Q(:,83) =[ C10, C10, C0 ,-C10];Q(:,84) =[-C9 ,-C8 ,-C8 , C0 ]    
        Q(:,85) =[-C6 ,-C4 , C7 ,-C6 ];Q(:,86) =[-C12, C0 ,-C12,-C11];Q(:,87) =[-C3 , C2 ,-C1 , C2 ];Q(:,88) =[ C13,-C14, C13, C13]     
        Q(:,89) =[ C6 ,-C6 ,-C4 , C7 ];Q(:,90) =[-C2 , C1 , C3 ,-C2 ];Q(:,91) =[-C12,-C0 ,-C11, C12];Q(:,92) =[ C8 ,-C0 ,-C9 , C8 ]    
        Q(:,93) =[-C2 ,-C2 ,-C3 ,-C1 ];Q(:,94) =[ C2 , C3 ,-C2 , C1 ];Q(:,95) =[ C0 ,-C12,-C12,-C11];Q(:,96) =[ C12,-C12,-C0 ,-C11]     
        Q(:,97) =[ C1 , C1 , C1 ,-C1 ];Q(:,98) =[ C2 , C1 ,-C2 ,-C3 ];Q(:,99) =[-C10, C0 , C10, C10];Q(:,100)=[ C8 ,-C9 , C8 ,-C0 ]    
        Q(:,101)=[ C4 ,-C0 ,-C5 ,-C4 ];Q(:,102)=[ C13, C14, C13,-C13];Q(:,103)=[-C3 ,-C1 ,-C2 ,-C2 ];Q(:,104)=[-C12,-C0 , C12, C11]     
        Q(:,105)=[-C6 ,-C4 , C6 ,-C7 ];Q(:,106)=[ C6 ,-C6 ,-C7 ,-C4 ];Q(:,107)=[-C7 , C4 ,-C6 , C6 ];Q(:,108)=[-C3 , C2 ,-C1 ,-C2 ]    
        Q(:,109)=[ C4 ,-C5 ,-C4 , C0 ];Q(:,110)=[-C7 , C0 , C7 , C0 ];Q(:,111)=[ C0 ,-C4 ,-C5 ,-C4 ];Q(:,112)=[ C4 ,-C6 , C6 ,-C7 ]    
        Q(:,113)=[ C0 ,-C9 , C8 , C8 ];Q(:,114)=[ C4 , C4 , C0 , C5 ];Q(:,115)=[ C6 , C4 , C7 ,-C6 ];Q(:,116)=[-C4 ,-C4 , C0 , C5 ]    
        Q(:,117)=[-C13,-C13, C14, C13];Q(:,118)=[-C8 ,-C9 , C8 ,-C0 ];Q(:,119)=[ C7 ,-C6 , C6 , C4 ];Q(:,120)=[ C0 ,-C12,-C11, C12]     
        Q(:,121)=[ C2 ,-C2 , C3 , C1 ];Q(:,122)=[ C13, C13,-C13,-C14];Q(:,123)=[-C6 , C7 , C6 , C4 ];Q(:,124)=[ C2 , C1 , C2 ,-C3 ]    
        Q(:,125)=[ C0 ,-C12, C12, C11];Q(:,126)=[ C0 ,-C8 ,-C8 ,-C9 ];Q(:,127)=[-C12,-C0 ,-C12, C11];Q(:,128)=[ C2 , C3 ,-C2 , C1 ]    
        Q(:,129)=[ C3 ,-C1 ,-C2 ,-C2 ];Q(:,130)=[-C2 , C1 ,-C3 , C2 ];Q(:,131)=[ C0 , C0 , C0 ,-C15];Q(:,132)=[-C2 ,-C3 ,-C2 ,-C1 ]    
        Q(:,133)=[-C10, C0 , C10, C10];Q(:,134)=[-C8 ,-C8 , C9 , C0 ];Q(:,135)=[-C4 ,-C7 ,-C6 ,-C6 ];Q(:,136)=[ C13,-C13,-C13,-C14]
    ENDIF

    IF(COLLMOLT) THEN
        RCLU = R(:,1:136) 
        QCLU = Q(:,1:136)
        R = 0.D0
        Q = 0.D0
    
        RSITE(:,1)= [  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
        RSITE(:,2)= [ -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
        RSITE(:,3)= [  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]
        RSITE(:,4)= [ -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]
    
        CN0 = 1.0_dp / SQRT(3.0_dp)
        RM  = 0.0_dp
        IF(CPPT) THEN
            BL = CPPLAM
        ELSE
            BL = 1.01_dp
        ENDIF
    
        DO J1 = 1, 136
            DO J2 = 1,4
                INDX = 4*(J1-1) + J2
                USITE = MATMUL(Q_TO_RM(QCLU(:,J1)), RSITE(:,J2))
                R(:,INDX) = RCLU(:,J1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*USITE
                IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( USITE )
            ENDDO
        ENDDO
    
        DO J1 = 1, 544
            R(:,J1) = R(:,J1) + BL*(SQRT(3.0_dp) / (2.0_dp*SQRT(2.0_dp)) )
        ENDDO
    ENDIF

END SUBROUTINE CLATHRATE_II

SUBROUTINE CLATHRATE_III(UCL,COLLMOLT)

    USE COMMONS, ONLY: DP, R, Q, RIGIDT, BOX, RHO, NUCX, NUCY, NUCZ, DENSITYT, PI, CPPT, CPPLAM
    USE ROTATIONS_MODULE, ONLY: UV_TO_Q, QUATMUL, Q_TO_RM
    
    IMPLICIT NONE
    
    INTEGER       :: J1, J2, INDX
    REAL(KIND=DP) :: UCL, RCLU(3,12), QCLU(4,12), CN0, BL, RM(3,3), RSITE(3,4), USITE(3)
    REAL(KIND=DP) :: C1, S1, C2, S2
    LOGICAL       :: COLLMOLT
    
    IF(COLLMOLT) THEN
        IF (DENSITYT) THEN
            UCL     = (48.0_dp/RHO)**(1.0_dp/3.0_dp)
        ELSE
            UCL  = (48.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
        ENDIF
        PRINT *, "Tetrahedral Clathrate sIII box:", BOX
    ELSE
        !   Unit cell length (remember there are 24 particles in the unit cell)
        IF (DENSITYT) THEN
            UCL     = (12.0_dp/RHO)**(1.0_dp/3.0_dp)
        ELSE
            UCL  = (12.0_dp*PI/(6.0_dp*RHO))**(1.0_dp/3.0_dp)
        ENDIF
        PRINT *, "Clathrate sIII box:", BOX
    ENDIF
    
    BOX     = [ UCL*REAL(NUCX,DP), UCL*REAL(NUCY,DP), UCL*REAL(NUCZ,DP) ]
    
    !   Build the unit cell
    R(:,1) = [0.50_dp, 0.25_dp, 0.00_dp]*UCL
    R(:,2) = [0.00_dp, 0.75_dp, 0.50_dp]*UCL
    R(:,3) = [0.25_dp, 0.00_dp, 0.50_dp]*UCL
    R(:,4) = [0.75_dp, 0.50_dp, 0.00_dp]*UCL
    R(:,5) = [0.00_dp, 0.50_dp, 0.25_dp]*UCL
    R(:,6) = [0.50_dp, 0.00_dp, 0.75_dp]*UCL
    R(:,7) = [0.00_dp, 0.50_dp, 0.75_dp]*UCL
    R(:,8) = [0.50_dp, 0.00_dp, 0.25_dp]*UCL
    R(:,9) = [0.50_dp, 0.75_dp, 0.00_dp]*UCL
    R(:,10)= [0.00_dp, 0.25_dp, 0.50_dp]*UCL
    R(:,11)= [0.75_dp, 0.00_dp, 0.50_dp]*UCL
    R(:,12)= [0.25_dp, 0.50_dp, 0.00_dp]*UCL
    
    IF(RIGIDT) THEN
        C1 = COS(PI/8.0_dp)
        S1 = SIN(PI/8.0_dp)
        C2 = 1.0_dp/SQRT(2.0_dp)*C1
        S2 = 1.0_dp/SQRT(2.0_dp)*S1

        Q(:,1) = [ -C1,   0.0_dp,  -S1,   0.0_dp]
        Q(:,2) = [  S1,   0.0_dp,  -C1,   0.0_dp]
        Q(:,3) = [ -S2,    -C2,    -S2,    -C2  ]
        Q(:,4) = [0.0_dp, 0.0_dp,   S1,    -C1  ]
        Q(:,5) = [ -C1,   0.0_dp, 0.0_dp,  -S1  ]
        Q(:,6) = [ -S2,    -S2,    -C2,    -C2  ]
        Q(:,7) = [0.0_dp,   C1,    -S1,   0.0_dp]
        Q(:,8) = [ -S2,     S2,    -C2,     C2  ]
        Q(:,9) = [  S2,    -C2,    -C2,     S2  ]
        Q(:,10)= [0.0_dp,   C1,   0.0_dp,   S1  ]
        Q(:,11)= [ -S2,     C2,    -S2,     C2  ]
        Q(:,12)= [  S2,    -C2,    -C2,    -S2  ]
    ENDIF

    IF(COLLMOLT) THEN
        RCLU = R(:,1:12) 
        QCLU = Q(:,1:12)
        R = 0.D0
        Q = 0.D0
    
        RSITE(:,1)= [  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
        RSITE(:,2)= [ -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
        RSITE(:,3)= [  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]
        RSITE(:,4)= [ -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]
    
        CN0 = 1.0_dp / SQRT(3.0_dp)
        RM  = 0.0_dp
        IF(CPPT) THEN
            BL = CPPLAM
        ELSE
            BL = 1.01_dp
        ENDIF
    
        DO J1 = 1, 12
            DO J2 = 1,4
                INDX = 4*(J1-1) + J2
                USITE = MATMUL(Q_TO_RM(QCLU(:,J1)), RSITE(:,J2))
                R(:,INDX) = RCLU(:,J1) + BL*( SQRT(3.0_dp)/(2.0_dp*SQRT(2.0_dp)) )*USITE
                IF(RIGIDT) Q(:,INDX) = UV_TO_Q ( USITE )
            ENDDO
        ENDDO
    
        DO J1 = 1, 48
            R(:,J1) = R(:,J1) + BL*(SQRT(3.0_dp) / (2.0_dp*SQRT(2.0_dp)) )
        ENDDO
    ENDIF

END SUBROUTINE CLATHRATE_III