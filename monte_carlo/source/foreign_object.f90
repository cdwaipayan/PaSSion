SUBROUTINE OBJECT_PARTICLE_ENERGY(ENERGY, J1, RI)

    USE COMMONS, ONLY: DP, NDIM, OBJCYLT, OBJSURFT, OBJSPHERET

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1
    REAL(KIND=DP), INTENT(IN)  :: RI(NDIM)
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    IF(OBJCYLT) THEN
        CALL CYLINDER(ENERGY, RI)
    ELSEIF(OBJSURFT) THEN
        CALL SURFACE(ENERGY, J1, RI)
    ELSEIF(OBJSPHERET) THEN
        CALL SPHERE(ENERGY, J1, RI)
    ENDIF

END SUBROUTINE

SUBROUTINE DEF_SURFACE_WELLS()

    USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
    USE COMMONS, ONLY: DP, NDIM, BOX, N_POLY, SURF_WELLS

    IMPLICIT NONE

    ALLOCATE( SURF_WELLS(NDIM-1,N_POLY) )

!   Divide a square plane into n_poly segments.

    IF(N_POLY==1) THEN
        SURF_WELLS(:,1) = 0.0_dp
    ELSEIF(N_POLY==4) THEN
        SURF_WELLS(:,1) = [ 0.5_dp,-0.5_dp]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,2) = [ 0.5_dp, 0.5_dp]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,3) = [-0.5_dp,-0.5_dp]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,4) = [-0.5_dp, 0.5_dp]*BOX(1:2)/2.0_dp
    ELSEIF(N_POLY==9) THEN
        SURF_WELLS(:,1) = [ 1.0_dp/3.0_dp,-1.0_dp/3.0_dp]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,2) = [ 1.0_dp/3.0_dp,    0.0_dp    ]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,3) = [ 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,4) = [    0.0_dp,    -1.0_dp/3.0_dp]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,5) = [    0.0_dp,        0.0_dp    ]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,6) = [    0.0_dp,     1.0_dp/3.0_dp]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,7) = [-1.0_dp/3.0_dp,-1.0_dp/3.0_dp]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,8) = [-1.0_dp/3.0_dp,    0.0_dp    ]*BOX(1:2)/2.0_dp
        SURF_WELLS(:,9) = [-1.0_dp/3.0_dp, 1.0_dp/3.0_dp]*BOX(1:2)/2.0_dp
    ! ELSEIF(N_POLY==16) THEN

    ELSE
        WRITE ( UNIT=ERROR_UNIT, FMT='(A)') 'INCORRECT NO. OF POLYMER CHAINS'
        STOP 'Error in def_surface'
    ENDIF


END SUBROUTINE

SUBROUTINE SURFACE(ENERGY, J1, RI)

    USE COMMONS, ONLY: DP, NDIM, OVERLAPT, BOX, KFT
    USE COMMONS, ONLY: NSITES, RBSITES, KFDEL, KFLAM, KFIJ
    USE COMMONS, ONLY: CAPST, POLYCHAINT, N_POLY_TOT, N_POLY_L, SURF_WELLS

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1
    INTEGER                    :: J3!, POLYID
    REAL(KIND=DP), INTENT(IN)  :: RI(NDIM)
    REAL(KIND=DP)              :: DIST, RHAT(NDIM), EA(NDIM), EARIJ!, HFLBOX, RXY(NDIM-1), DXY!, RIJSQ
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    DIST   = RI(3) 
    ENERGY = 0.0_dp

    IF(KFT) THEN
        IF(DIST<0.5_dp .OR. DIST > BOX(3)/2.0_dp) THEN
            OVERLAPT = .TRUE.
            RETURN
        ENDIF

        DO J3 = 1, NSITES
            IF(DIST<=(KFLAM(J3,J3)-0.5_dp)) THEN
            !   Direction of patch alpha on particle I
                EA  = RBSITES(:,J3,J1)
                EARIJ = -DOT_PRODUCT(EA,[0.0_dp,0.0_dp,1.0_dp])
            !   If normalised distance vector doesn't pass through patch alpha
            !   the conditions for bonding are not met so no need to progress.
                IF(EARIJ <= KFDEL(J3)) THEN
                    CYCLE
                ELSE
                !   The conditions for bonding are met so calculate the contribution 
                !   of the interaction between patches alpha and beta to the energy.
                    ENERGY = ENERGY - KFIJ(J3,J3)
                ENDIF
            ENDIF
        ENDDO ! Loop over each of the patches on particle i

    ELSEIF(POLYCHAINT) THEN
        STOP 'ANCHORED POLYMER CODE NEEDS UPDATING!'
    !--------------------------------------------------------------
    !   NEEDS UPDATING FOR PATCHY TRIANGLES AND PATCHY DISKS!!!!!
    !--------------------------------------------------------------
    !     IF(J1>N_POLY_TOT) THEN
    ! !   Scenario where particle I is a capsomer
    !         ENERGY = ENERGY + CAP_APX * (0.5_dp*CAP_S/DIST)**48
    !     ELSE
    ! !   Scenario where particle I is a part of a polymer chain
    !         IF(MOD((J1-1)+N_POLY_L,N_POLY_L)==0) THEN
    !     ! Check if particle I is at the start of the polymer chain, 
    !     ! if it is then it is tethered to the surface and so there is an 
    !     ! attractive interaction as well as a repulsive interaction.
    !         !   Harmonic spring along z-axis
    !             ENERGY = ENERGY + 10.0_dp*POLY_KAP*(DIST-POLYC_L*POLY_SIG)**2
    !         !   Harmonic well in xy-plane
    !             POLYID = (J1-1)/N_POLY_L + 1
    !             RXY = RI(1:2) - SURF_WELLS(:,POLYID)
    !             RXY = RXY - BOX(1:2)*ANINT(RXY/BOX(1:2))
    !             DXY = NORM2(RXY)
    !             ENERGY = ENERGY + 10.0_dp*POLY_KAP*(DXY-POLYC_L*POLY_SIG)**2
    !         ENDIF
    !     !   Repulsive interaction between the surface and the bead
    !         ENERGY = ENERGY + (0.5_dp*POLY_SIG/DIST)**32
    !     ENDIF

    ENDIF

END SUBROUTINE

SUBROUTINE SPHERE(ENERGY, J1, RI)
    USE COMMONS, ONLY: DP, NDIM, R, OVERLAPT, BOX, KFT, PATCHYTRIT
    USE COMMONS, ONLY: NSITES, RBSITES, KFDEL, KFLAM, KFIJ, SPHERERAD, SPHERECNFT
    USE COMMONS, ONLY: CAPST, POLYCHAINT, N_POLY_TOT, N_POLY_L, SURF_WELLS

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1
    INTEGER                    :: J3
    REAL(KIND=DP), INTENT(IN)  :: RI(NDIM)
    REAL(KIND=DP)              :: R2, DIST, RHAT(NDIM), EA(NDIM), EARIJ
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp

    IF(KFT) THEN
        R2   = NORM2(RI)
        IF(SPHERECNFT) THEN
            DIST = SPHERERAD-(DIST+0.5_dp)
        ELSE
            DIST = R2 - SPHERERAD
            ! PRINT *, DIST
        ENDIF

        IF (DIST < 0.5_dp) THEN
            ! ENERGY = ENERGY + 10.0_dp
            OVERLAPT = .TRUE.
            RETURN
        ENDIF

        IF(OVERLAPT) RETURN

        DO J3 = 1, NSITES
            IF(DIST>(KFLAM(J3,J3)-0.5_dp)) CYCLE
            RHAT = RI / DIST
        !   Direction of patch alpha on particle I
            EA  = RBSITES(:,J3,J1)
            ! IF(SPHERECNFT) THEN
            EARIJ = DOT_PRODUCT(EA,RHAT)
            ! ELSE
            !     EARIJ = -DOT_PRODUCT(EA,RHAT)
            ! ENDIF
        !   If normalised distance vector doesn't pass through patch alpha
        !   the conditions for bonding are not met so no need to progress.
            IF(EARIJ <= KFDEL(J3)) THEN
                CYCLE
            ELSE
            !   The conditions for bonding are met so calculate the contribution 
            !   of the interaction between patches alpha and beta to the energy.
                ENERGY = ENERGY - KFIJ(J3,J3)
            ENDIF
        ENDDO ! Loop over each of the patches on particle i
    
    ! ELSEIF(PATCHYTRIT) THEN

    ENDIF

END SUBROUTINE

SUBROUTINE CYLINDER(ENERGY, RI)

    USE COMMONS, ONLY: DP, NDIM, HARDT, OVERLAPT, GEN_CONFIGT
    USE COMMONS, ONLY: GLJT, YUKT, CYLR, CYLR2

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN)  :: RI(NDIM)
    REAL(KIND=DP)              :: RIJSQ, DIST
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    RIJSQ = DOT_PRODUCT(RI(1:2),RI(1:2)) !- 0.025_dp
    DIST  = SQRT(RIJSQ) - 0.5_dp
    RIJSQ = DIST**2

    IF (HARDT) THEN
        IF(RIJSQ < CYLR2) THEN
            IF(GEN_CONFIGT) THEN
                ENERGY = 100.0_dp
            ELSE
                OVERLAPT = .TRUE.
            ENDIF
        ENDIF
    ELSE 
!   ---------------------------------------------------
!   Pair potentials which have soft-core repulsion
!   ---------------------------------------------------    
        IF(YUKT) THEN
            CALL YUKAWA(ENERGY, RIJSQ)

        ELSEIF(GLJT) THEN
            CALL GLJ(ENERGY, RIJSQ)
            ENERGY = ENERGY * 2.0_dp!1.2055_dp
        ELSE
            PRINT *, "NO POTENTIAL SELECTED, STOPPING PROGRAM."
            STOP " Stopping in potential.f90"
        ENDIF

    ENDIF

END SUBROUTINE

SUBROUTINE CONE(ENERGY, RI)

    USE COMMONS, ONLY: DP, NDIM, HARDT, OVERLAPT, GEN_CONFIGT
    USE COMMONS, ONLY: GLJT, YUKT, CYLR, CYLR2

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN)  :: RI(NDIM)
    REAL(KIND=DP)              :: RIJSQ, DIST
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ! cone_r=cone_z*np.tan(cone_alpha)
    RIJSQ = DOT_PRODUCT(RI(1:2),RI(1:2)) !- 0.025_dp
    DIST  = SQRT(RIJSQ) - 0.5_dp
    RIJSQ = DIST**2

    IF (HARDT) THEN
        IF(RIJSQ < CYLR2) THEN
            IF(GEN_CONFIGT) THEN
                ENERGY = 100.0_dp
            ELSE
                OVERLAPT = .TRUE.
            ENDIF
        ENDIF
    ELSE 
!   ---------------------------------------------------
!   Pair potentials which have soft-core repulsion
!   ---------------------------------------------------    
        IF(YUKT) THEN
            CALL YUKAWA(ENERGY, RIJSQ)

        ELSEIF(GLJT) THEN
            CALL GLJ(ENERGY, RIJSQ)
            ENERGY = ENERGY * 1.2055_dp
        ELSE
            PRINT *, "NO POTENTIAL SELECTED, STOPPING PROGRAM."
            STOP " Stopping in potential.f90"
        ENDIF

    ENDIF

END SUBROUTINE