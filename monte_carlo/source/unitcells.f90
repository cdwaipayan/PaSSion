SUBROUTINE GENERATE_UNIT_CELLS()

    USE COMMONS, ONLY: DP, NPART, LID, VLM, RIGIDT, EQUISEEDT, RANDQUATT, SETORTN, NLAYERS, SPHERECNFT, SPHERERAD
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

    IF(SPHERECNFT) THEN
        VLM = 4.0_dp/3.0_dp * PI * SPHERERAD**3
        BOX = SPHERERAD*2.0_dp
        RHO = REAL(NPART,DP)*PI/(6.0_dp*VLM)
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