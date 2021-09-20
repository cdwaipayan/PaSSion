SUBROUTINE KF_REC(ENERGY, RIJ, RIJSQ, J1, J2)
    !==============================================================================================
    !   Subroutine to calculate the potential energy for a system of patchy particles interacting
    !   via the Kern-Frenkel potential.
    !==============================================================================================
    USE COMMONS, ONLY: DP, NDIM, RBSITES
    USE COMMONS, ONLY: KF_LAM2, KFDEL, KFIJ
    USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, CLSTRID, CLSTRSITEID, CLSTRADJ, VLMCLUSTERMOVET, CLURIJ, CLSTRSITEID, LRGCLSTMVT

    IMPLICIT NONE

    INTEGER, INTENT(IN)         :: J1, J2
    REAL(KIND=DP), INTENT(IN)   :: RIJ(NDIM), RIJSQ

    INTEGER                     :: J3, J4
    REAL(KIND=DP)               :: EA1(NDIM), EA2(NDIM), EA3(NDIM)
    REAL(KIND=DP)               :: EARIJ1(NDIM), EARIJ2(NDIM), EARIJ3
    REAL(KIND=DP)               :: RI12(NDIM), EARI12
    REAL(KIND=DP)               :: EB1(NDIM), EB2(NDIM), EB3(NDIM)
    REAL(KIND=DP)               :: EBRJI1(NDIM), EBRJI2(NDIM), EBRJI3
    REAL(KIND=DP)               :: RJ12(NDIM), EBRJ12

    REAL(KIND=DP), INTENT(OUT)  :: ENERGY

    ENERGY   = 0.0_dp

    IF (RIJSQ > KF_LAM2) RETURN

    DO J3 = 1, 2
    !   Direction of patch alpha on particle I
        EA1  = RBSITES(:,(J3-1)*3+1,J1)
    !   Direction of long-axis of rectangular patch
        EA2  = RBSITES(:,(J3-1)*3+2,J1)
    !   Direction of short-axis of the rectangular patch
        EA3  = RBSITES(:,(J3-1)*3+3,J1)

    !   Calculate the cosine of the angle between RIJ and EA3
        EARIJ3 = -DOT_PRODUCT(EA3,RIJ/NORM2(RIJ))
        IF(EARIJ3 <= -KFDEL((J3-1)*2+2) .OR. EARIJ3 >= KFDEL((J3-1)*2+2)) CYCLE !BELT CRITERIA FIRST BECAUSE IT IS STRICT AND CHEAP

    !   Project the distance vector onto the axes defining the coordinate frame of the patch
        EARIJ1 = -DOT_PRODUCT(EA1,RIJ)*EA1
        EARIJ2 = -DOT_PRODUCT(EA2,RIJ)*EA2
    !   Calculate the projection of the distance vector onto the plane defined by axes 1&2
        RI12 =  EARIJ1 + EARIJ2
    !   Calculate the cosine of the angle between the projection and the direction vector of the patch
        EARI12 = DOT_PRODUCT(RI12/NORM2(RI12), EA1)
    !   If normalised distance vector doesn't pass through patch A
    !   the conditions for bonding are not met so no need to progress.
        IF(EARI12 <= KFDEL((J3-1)*2+1)) CYCLE
        
        DO J4 = 1, 2
        !   Direction of patch beta on particle J
            EB1  = RBSITES(:,(J4-1)*3+1,J2)
        !   Direction of long-axis of rectangular patch
            EB2  = RBSITES(:,(J4-1)*3+2,J2)
        !   Direction of short-axis of the rectangular patch
            EB3  = RBSITES(:,(J4-1)*3+3,J2)

        !   Calculate the cosine of the angle between RIJ and EB3
            EBRJI3 = DOT_PRODUCT(EB3,RIJ/NORM2(RIJ))
            IF(EBRJI3 <= -KFDEL((J4-1)*2+2) .OR. EBRJI3 >= KFDEL((J4-1)*2+2)) CYCLE !BELT CRITERIA FIRST BECAUSE IT IS STRICT AND CHEAP

        !   Project the distance vector onto the axes defining the coordinate frame of the patch
            EBRJI1 = DOT_PRODUCT(EB1,RIJ)*EB1
            EBRJI2 = DOT_PRODUCT(EB2,RIJ)*EB2
        !   Calculate the projection of the distance vector onto the plane defined by axes 1&2
            RJ12 =  EBRJI1 + EBRJI2
        !   Calculate the cosine of the angle between the projection and the direction vector of the patch
            EBRJ12 = DOT_PRODUCT(RJ12/NORM2(RJ12), EB1)
        !   If normalised distance vector doesn't pass through patch B
        !   the conditions for bonding are not met so no need to progress.
            IF(EBRJ12 <= KFDEL((J4-1)*2+1)) THEN
                CYCLE
            ELSE

            !   The conditions for bonding are met so calculate the contribution 
            !   of the interaction between patches A and B to the energy.
                ENERGY = ENERGY - KFIJ(J3,J4)
            
                IF(CLUSTERT .OR. VLMCLUSTERMOVET) THEN
                !   If simulating triblock patchy particles then we only consider particles
                !   interacting via patch B-patch B bonds to be apart of the same cluster.
                    IF( .NOT. LRGCLSTMVT ) THEN
                        IF( J3 /= CLSTRSITEID .OR. J4 /= CLSTRSITEID ) CYCLE ! DO NOTHING  
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
                        IF(.NOT.( ANY(CLSTR(CLSTRID+1:CLSTRSZ)==J2) )) THEN
                            CLSTRSZ = CLSTRSZ + 1
                            CLSTR(CLSTRSZ) = J2
                        ENDIF
                    ENDIF
                ENDIF

            ENDIF
        ENDDO ! Loop over each of the patches on particle j
    ENDDO ! Loop over each of the patches on particle i

END SUBROUTINE KF_REC
    
!====================================================================================================

SUBROUTINE DEF_KF_REC()
!----------------------------------------------------------------
! Define the position of the patches on particles interacting via
! the Kern-Frenkel potential.
! Also define constants to be used during the simulations which 
! depend on the width of the patches and the well-depth of the 
! various patch-patch interactions.
! The well-depth of patch A is used to define the reduced unit
! of temperature and so KFAA=1 for all cases.
! Currently, for a given number of patches, the geometry of the
! patches is fixed.
!----------------------------------------------------------------

    USE COMMONS, ONLY: DP, PI, REFSITE, REFSITE2, RACEMICT, SKEWAB2
    USE COMMONS, ONLY: KFDELA1, KFDELA2, KFDELB1, KFDELB2, KFIJ, KFDEL, KFAA, KFBB, SKEWAB
    USE COMMONS, ONLY: TANA, TANB

    IMPLICIT NONE

    REAL(KIND=DP) :: SS, CS, RMZ(3,3)

    ALLOCATE(KFIJ(2,2), KFDEL(4))

    KFIJ(1,1) = KFAA
    KFIJ(1,2) = SQRT(KFAA*KFBB)
    KFIJ(2,1) = KFIJ(1,2)
    KFIJ(2,2) = KFBB

    KFDEL(1)  = COS(KFDELA1*PI/180_dp)
    KFDEL(2)  = COS(PI/2_dp - KFDELA2*PI/180_dp) !BELT
    KFDEL(3)  = COS(KFDELB1*PI/180_dp)
    KFDEL(4)  = COS(PI/2_dp - KFDELB2*PI/180_dp) !BELT

    TANA = TAN(KFDELA1*PI/180_dp)
    TANB = TAN(KFDELB1*PI/180_dp)
    
    SS = SIN(SKEWAB*PI/180_dp)
    CS = COS(SKEWAB*PI/180_dp)
    RMZ(1,:) = [  CS,    -SS,   0.0_dp]
    RMZ(2,:) = [  SS,     CS,   0.0_dp]
    RMZ(3,:) = [0.0_dp, 0.0_dp, 1.0_dp]

    REFSITE(:,1) = [ 0.0_dp,  0.0_dp, -1.0_dp]
    REFSITE(:,2) = [ 0.0_dp, -1.0_dp,  0.0_dp]
    REFSITE(:,3) = [-1.0_dp,  0.0_dp,  0.0_dp]
    REFSITE(:,4) = [ 0.0_dp,  0.0_dp,  1.0_dp]
    REFSITE(:,5) = MATMUL(RMZ, REFSITE(:,2))
    REFSITE(:,6) = MATMUL(RMZ, REFSITE(:,3))

    IF(RACEMICT) THEN
        SS = SIN(SKEWAB2*PI/180_dp)
        CS = COS(SKEWAB2*PI/180_dp)
        RMZ(1,:) = [  CS,    -SS,   0.0_dp]
        RMZ(2,:) = [  SS,     CS,   0.0_dp]
        RMZ(3,:) = [0.0_dp, 0.0_dp, 1.0_dp]

        REFSITE2(:,1) = [ 0.0_dp,  0.0_dp, -1.0_dp]
        REFSITE2(:,2) = [ 0.0_dp, -1.0_dp,  0.0_dp]
        REFSITE2(:,3) = [-1.0_dp,  0.0_dp,  0.0_dp]
        REFSITE2(:,4) = [ 0.0_dp,  0.0_dp,  1.0_dp]
        REFSITE2(:,5) = MATMUL(RMZ, REFSITE2(:,2))
        REFSITE2(:,6) = MATMUL(RMZ, REFSITE2(:,3))
    ENDIF

    
END SUBROUTINE

!====================================================================================================

SUBROUTINE VIEW_KF_REC()

    USE COMMONS, ONLY: DP, NPART, R, RBSITES, BOX, VIEWUNIT, TANA, TANB, RACEMICT

    IMPLICIT NONE

    INTEGER:: J1
    REAL(KIND=DP) :: RWRITE(3,NPART), RC1

    WRITE(VIEWUNIT,*) NPART*7
    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO J1 = 1, NPART
        IF(RACEMICT) RC1 = MOD(((J1-1)-MOD((J1-1),12))/12+1,2)

        RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX) * BOX
        
        IF(RACEMICT .AND. RC1 == 0) THEN
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Co ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
        ELSE
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
        ENDIF
        
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ',  RWRITE(:,J1) + 0.37_dp*RBSITES(:,1,J1)
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ',  RWRITE(:,J1) + 0.3_dp*(RBSITES(:,1,J1)+TANA*RBSITES(:,2,J1))
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ',  RWRITE(:,J1) + 0.3_dp*(RBSITES(:,1,J1)-TANA*RBSITES(:,2,J1))

        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ',  RWRITE(:,J1) + 0.37_dp*RBSITES(:,4,J1)
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ',  RWRITE(:,J1) + 0.3_dp*(RBSITES(:,4,J1)+TANB*RBSITES(:,5,J1))
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ',  RWRITE(:,J1) + 0.3_dp*(RBSITES(:,4,J1)-TANB*RBSITES(:,5,J1))
    END DO
    
END SUBROUTINE
