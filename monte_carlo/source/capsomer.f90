SUBROUTINE POLYGON_CAPS(ENERGY, RIJ, RIJSQ, J1, J2)
!   =====================================================================================================================

!   =====================================================================================================================
    USE COMMONS, ONLY: DP, NDIM, NSITES, R, RBSITES, CAP_P, CAP_RE, CAP_S, CAP_APX, CLUSTERT, VLMCLUSTERMOVET, KAPPA
    USE COMMONS, ONLY: CLSTR,CLSTRSZ,CLSTRID,CLSTRADJ,CLURIJ,CLSTRSITEID,LRGCLSTMVT,PRNTCNF,BOPCLSTRT,POLYCHAINT,SNGLCHAINET

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: J1, J2
    INTEGER             :: J3, J4

    REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM), RIJSQ
    REAL(KIND=DP) :: C1, AIJ, EI(NDIM), EJ(NDIM), EIJ(NDIM), VA, VP, DIJ

    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp
    VP     = 0.0_dp  

!   Compute Apex-Apex repulsive interaction
    VA     = CAP_APX * (CAP_S/SQRT(RIJSQ))**48
    ! DIJ = SQRT(RIJSQ)
    ! VA  = CAP_APX*EXP( -KAPPA*(DIJ-CAP_S) )/DIJ
    ! PRINT *, DIJ, CAP_S, VA
    ! STOP

!   Compute Vertex-Vertex attractive interaction
    DO J3 = 2,NSITES
        EI  = RBSITES(:,J3,J1)
        DO J4 = 2,NSITES
            EJ  = RBSITES(:,J4,J2)
            EIJ = RIJ + EI - EJ
            AIJ = NORM2(EIJ)
            IF(AIJ > CAP_P .AND. AIJ <= CAP_RE) THEN
                VP = VP - 1.0_dp
            ENDIF
            ! C1  = EXP( CAP_P*(1.0_dp-AIJ/CAP_RE) )
            ! VP  = VP + C1*(C1-2.0_dp)
        ENDDO
    ENDDO

    ENERGY = VA + VP

    IF((CLUSTERT .OR. VLMCLUSTERMOVET)) THEN
        IF(POLYCHAINT) THEN
            IF(.NOT. SNGLCHAINET) RETURN
        ENDIF
    !   If simulating triblock patchy particles then we only consider particles
    !   interacting via patch B-patch B bonds to be apart of the same cluster.
        IF(VP >= 0.0_dp) RETURN
        
    !   If performing a volume cluster move update the adjacency matrix and the pair
    !   distance matrix for the system.
        IF(VLMCLUSTERMOVET) THEN
            CLSTRADJ(J1,J2) = 1
            CLSTRADJ(J2,J1) = 1
            CLURIJ(:,J1,J2) = RIJ
            CLURIJ(:,J2,J1) = -RIJ
        ELSE
        !   Check if the particle is already in the cluster, only need to check 
        !   particle IDs that come after the current particle in the list as the
        !   previous check should take care of those that come earlier.
            IF(.NOT.( ANY(CLSTR(1:CLSTRSZ)==J2) )) THEN
                CLSTRSZ = CLSTRSZ + 1
                CLSTR(CLSTRSZ) = J2
            ENDIF
        ENDIF
    ENDIF

END SUBROUTINE POLYGON_CAPS

SUBROUTINE INTRA_POLYMER(ENERGY, J1)
!   =====================================================================================================================
!   Compute the harmonic energy associate with nearest neighbours in freely-jointed polymer chains
!   =====================================================================================================================
    USE COMMONS, ONLY: DP, NDIM, BOX, R, POLYC_L, POLY_SIG, POLY_KAP, N_POLY_L, SPET, RCHTT

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1
    REAL(KIND=DP)              :: RIJ(NDIM), DIST
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp
    
    IF(MOD(J1,N_POLY_L)/=1 .AND. SPET ) THEN
        RIJ    = R(:,J1) - R(:,J1-1)
        RIJ    = RIJ - BOX*ANINT( RIJ/BOX )
        DIST   = NORM2(RIJ)
        ENERGY = ENERGY + POLY_KAP*(DIST-POLYC_L*POLY_SIG)**2
    ENDIF

    IF(MOD(J1,N_POLY_L)/=0) THEN
        RIJ    = R(:,J1) - R(:,J1+1)
        RIJ    = RIJ - BOX*ANINT( RIJ/BOX )
        DIST   = NORM2(RIJ)
        ENERGY = ENERGY + POLY_KAP*(DIST-POLYC_L*POLY_SIG)**2
    ENDIF

END SUBROUTINE INTRA_POLYMER

SUBROUTINE INTER_POLYMER(ENERGY, RIJSQ, J1, J2)
!   =====================================================================================================================
!   Compute the soft repulsive energy associated with non-nearest neighbours in freely-jointed polymer chains
!   =====================================================================================================================
    USE COMMONS, ONLY: DP, POLY_SIG2!, N_POLY_L
    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1, J2
    REAL(KIND=DP), INTENT(IN)  :: RIJSQ
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = (POLY_SIG2/RIJSQ)**16
    
END SUBROUTINE INTER_POLYMER

SUBROUTINE CAP_POLYMER(ENERGY, RIJ, RIJSQ, J1, J2)
!   =====================================================================================================================
!   Polymer-capsomer interaction, particle I is always a polymer bead and particle J is always a capsomer
!   =====================================================================================================================
    USE COMMONS, ONLY: DP, NDIM, RBSITES, CLUSTERT, VLMCLUSTERMOVET, CP_EPS, CP_SIG1, CP_SIG2, CP_CUT2
    USE COMMONS, ONLY: CLSTR,CLSTRSZ,CLSTRID,CLSTRADJ,CLURIJ,LRGCLSTMVT,PRNTCNF,BOPCLSTRT, SNGLCHAINET

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1, J2
    REAL(KIND=DP), INTENT(IN)  :: RIJ(NDIM), RIJSQ
    REAL(KIND=DP)              :: C2(NDIM), CR, R2, VA, VR
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp

!   Capsomer has attractive site centred on particle ~COM and 
!   repulsive site centred on apex
    IF(J1 > J2) THEN
        C2 = RIJ + RBSITES(:,1,J1) 
    ELSE
        C2 = RIJ - RBSITES(:,1,J2) 
    ENDIF

    CR = DOT_PRODUCT(C2,C2)
    IF(CR<CP_CUT2) THEN
        VA = -CP_EPS
    ELSE
        VA = 0.0_dp
    ENDIF
    
    VR = (CP_SIG1/SQRT(RIJSQ))**32
    ! R2 = (CP_SIG2/NORM2(C2))**6
    ! VA = 4.0_dp*CP_EPS*(R2**2 - R2)

    ENERGY = VR + VA

    IF(CLUSTERT .OR. VLMCLUSTERMOVET) THEN
    !   If simulating triblock patchy particles then we only consider particles
    !   interacting via patch B-patch B bonds to be apart of the same cluster.
        IF( .NOT. LRGCLSTMVT ) THEN
            IF( VA >= 0.0_dp ) RETURN
        ENDIF
    !   If performing a volume cluster move update the adjacency matrix and the pair
    !   distance matrix for the system.
        IF(VLMCLUSTERMOVET) THEN
            CLSTRADJ(J1,J2) = 1
            CLSTRADJ(J2,J1) = 1
            CLURIJ(:,J1,J2) = RIJ
            CLURIJ(:,J2,J1) = -RIJ
        ELSE
        !   Check if the particle is already in the cluster, only need to check 
        !   particle IDs that come after the current particle in the list as the
        !   previous check should take care of those that come earlier.
            IF(.NOT.( ANY(CLSTR(1:CLSTRSZ)==J2) ) .AND. J2>J1) THEN
                CLSTRSZ = CLSTRSZ + 1
                CLSTR(CLSTRSZ) = J2
            ENDIF
        ENDIF
    ENDIF

    IF(PRNTCNF .OR. BOPCLSTRT) THEN
        CLSTRADJ(J1,J2) = 1
        CLSTRADJ(J2,J1) = 1
        CLURIJ(:,J1,J2) = RIJ
        CLURIJ(:,J2,J1) = -RIJ
    ENDIF

END SUBROUTINE CAP_POLYMER

!   =====================================================================================================================
!   =====================================================================================================================
!   =====================================================================================================================
!   =====================================================================================================================
SUBROUTINE DEF_POLYGON_CAPS()

        USE COMMONS, ONLY: DP, PI, NSITES, REFSITE, CAP_CP_SIG, CP_CUT2, POLYCHAINT, CAP_THETA
        USE COMMONS, ONLY: CAP_RB, CAP_H, CP_SIG1, CP_SIG2, CP_SIG12, CP_SIG22, POLY_SIG, CAP_S
    
        IMPLICIT NONE

        REAL(KIND=DP) :: PENT_ANGL, RMY(3,3), RMZ(3,3)

        PENT_ANGL = 72.0_dp*(PI/180_dp)
        CAP_THETA = (90.0_dp+CAP_THETA)*(PI/180_dp)

        RMY = 0.0_dp
        RMZ = 0.0_dp

        RMY(1,1) = COS(CAP_THETA)
        RMY(1,3) = SIN(CAP_THETA)
        RMY(2,2) = 1.0_dp
        RMY(3,1) = -SIN(CAP_THETA)
        RMY(3,3) = COS(CAP_THETA)

        RMZ(1,1) = COS(PENT_ANGL)
        RMZ(1,2) = SIN(PENT_ANGL)
        RMZ(2,1) = -SIN(PENT_ANGL)
        RMZ(2,2) = COS(PENT_ANGL)
        RMZ(3,3) = 1.0_dp

        IF(NSITES == 6) THEN
            REFSITE(:,1)= [ 0.0_dp, 0.0_dp, -CAP_H]

            REFSITE(:,2) = [ CAP_RB, 0.0_dp, 0.0_dp]
            REFSITE(:,2) = MATMUL(RMY,REFSITE(:,2))
            ! REFSITE(:,2) = REFSITE(:,2) + REFSITE(:,1)

            REFSITE(:,3)=MATMUL(RMZ,REFSITE(:,2))
            REFSITE(:,4)=MATMUL(RMZ,REFSITE(:,3))
            REFSITE(:,5)=MATMUL(RMZ,REFSITE(:,4))
            REFSITE(:,6)=MATMUL(RMZ,REFSITE(:,5))

            ! REFSITE(:,2)= [ 0.0_dp, CAP_RB, 0.0_dp]
            ! REFSITE(:,3)= [ CAP_RB*COS(PI/10.0_dp),  CAP_RB*SIN(PI/10.0_dp), 0.0_dp]
            ! REFSITE(:,4)= [ CAP_RB*COS(-3.0_dp*PI/10.0_dp), CAP_RB*SIN(-3.0_dp*PI/10.0_dp), 0.0_dp]
            ! REFSITE(:,5)= [-CAP_RB*COS(PI/10.0_dp),  CAP_RB*SIN(PI/10.0_dp), 0.0_dp]
            ! REFSITE(:,6)= [-CAP_RB*COS(-3.0_dp*PI/10.0_dp), CAP_RB*SIN(-3.0_dp*PI/10.0_dp), 0.0_dp]
        ENDIF

        IF(POLYCHAINT) THEN
            CP_SIG1   = (POLY_SIG+CAP_S)/2.0_dp    ! repulsive sigma
            CP_SIG12  = CP_SIG1**2
            CP_SIG2   = (POLY_SIG+CAP_CP_SIG)/2.0_dp  ! attractive sigma
            CP_SIG22  = CP_SIG2**2
            CP_CUT2   = (CP_SIG2*1.1_dp)**2
        ENDIF

END SUBROUTINE DEF_POLYGON_CAPS

SUBROUTINE DEF_POLYMER_CHAIN()
    
    USE COMMONS, ONLY: DP, PI, NSITES, REFSITE
    USE COMMONS, ONLY: POLYC_L

    IMPLICIT NONE

    ! INTEGER ::  M, IZ, IY, IX, ND, JY, JX

    ! ND = INT(ANINT( NSITES**(1.0_dp/3.0_dp) ))

    ! M = 0
    ! REFSITE(:,1) = 0.0_dp

    ! DO IZ = 1, ND
    !     DO IY = 1, ND
    !         IF(MOD(IZ,2)/=1) THEN
    !             JY = ND - IY + 1
    !         ELSE
    !             JY = IY
    !         ENDIF
    !         DO IX = 1, ND
    !             IF(MOD(IY,2)/=1) THEN
    !                 JX = ND - IX + 1
    !             ELSE
    !                 JX = IX
    !             ENDIF
    !             REFSITE(1,1+M) = REFSITE(1,1) + POLYC_L*REAL((JX-1),DP)
    !             REFSITE(2,1+M) = REFSITE(2,1) + POLYC_L*REAL((JY-1),DP)
    !             REFSITE(3,1+M) = REFSITE(3,1) + POLYC_L*REAL((IZ-1),DP)
    !             M = M + 1
    !         ENDDO
    !     ENDDO
    ! ENDDO

END SUBROUTINE DEF_POLYMER_CHAIN

!====================================================================================================
    
SUBROUTINE VIEW_POLYGON_CAPS()
    
    USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT
    USE ROTATIONS_MODULE, ONLY: Q_TO_RM

    IMPLICIT NONE

    INTEGER:: J1, J2
    REAL(KIND=DP) :: RWRITE(3,NPART), RM(3,3), RBCOORDS(3)

    WRITE(VIEWUNIT,*) NPART*(NSITES+1)
    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO J1 = 1, NPART
        RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX(:))*BOX(:)
    END DO

    DO J1 = 1, NPART
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)

        RM = Q_TO_RM( Q(:,J1) )
        DO J2 = 1, NSITES
            RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
            IF(J2==1) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ELSE
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ENDIF
        ENDDO
    END DO
    
END SUBROUTINE

!====================================================================================================
    
SUBROUTINE VIEW_POLYMER()
    
    USE COMMONS, ONLY: DP, NPART, R, RBSITES, NSITES, VIEWUNIT, BOX

    IMPLICIT NONE

    INTEGER:: J1, J2
    REAL(KIND=DP) :: RBCOORDS(3), RWRITE(3,NPART)

    WRITE(VIEWUNIT,*) NPART!*(NSITES)
    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO J1 = 1, NPART
        RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX(:))*BOX(:)
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
    END DO

    ! DO J1 = 1, NPART
    !     DO J2 = 1, NSITES
    !         RBCOORDS = R(:,J1) + RBSITES(:,J2,J1)
    !         WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
    !     ENDDO
    ! END DO
    
END SUBROUTINE

!====================================================================================================
    
SUBROUTINE VIEW_POLY_AND_CAPS()
    
    USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT, N_POLY_TOT
    USE ROTATIONS_MODULE, ONLY: Q_TO_RM

    IMPLICIT NONE

    INTEGER:: J1, J2
    REAL(KIND=DP) :: RWRITE(3,NPART), RM(3,3), RBCOORDS(3)

    WRITE(VIEWUNIT,*) N_POLY_TOT+(NPART-N_POLY_TOT)*(NSITES+1)
    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO J1 = 1, NPART
        RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX(:))*BOX(:)
    END DO

    DO J1 = 1, NPART
        IF(J1<= N_POLY_TOT) THEN
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'H ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
        ELSE
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            RM = Q_TO_RM( Q(:,J1) )
            DO J2 = 1, NSITES
                RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
                IF(J2==1) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSE
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ENDIF
            ENDDO
        ENDIF
    END DO
    
END SUBROUTINE