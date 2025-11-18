    SUBROUTINE KF(ENERGY, RIJ, RIJSQ, J1, J2)
!==============================================================================================
!   Subroutine to calculate the potential energy for a system of patchy particles interacting
!   via the Kern-Frenkel potential.
!==============================================================================================
        USE COMMONS, ONLY: DP, NDIM, RBSITES, NSITES, ONEBONDT, NPBONDS!, PNEIGHS
        USE COMMONS, ONLY: RCUTSQ, KFLAM2, KFDEL, KFIJ, DNAT, BINARYT, BINARYID, BNRYRTIOT, BNRYA, SQSHT, SQSHDEL, SQSHEPS, SQID
        USE COMMONS, ONLY: CLUSTERT,CLSTR,CLSTRSZ,CLSTRID,CLSTRADJ,VLMCLUSTERMOVET,CLURIJ,CLSTRSITEID,LRGCLSTMVT,PRNTCNF,BOPCLSTRT
    
        IMPLICIT NONE
    
        INTEGER, INTENT(IN)         :: J1, J2
        REAL(KIND=DP), INTENT(IN)   :: RIJ(NDIM), RIJSQ

        INTEGER                     :: J3, J4
        REAL(KIND=DP)               :: RHAT(NDIM), DIST
        REAL(KIND=DP)               :: EA(NDIM), EB(NDIM)
        REAL(KIND=DP)               :: EARIJ, EBRJI
    
        REAL(KIND=DP), INTENT(OUT)  :: ENERGY

        ENERGY   = 0.0_dp
    
        IF (RIJSQ <= RCUTSQ) THEN
            DIST = SQRT(RIJSQ)
            RHAT = RIJ / DIST
        ELSE
            RETURN
        ENDIF


    !   If considering a binary system then the 1st and 2nd halves of the particle list do not iteract with one another.
        IF(BINARYT) THEN
            IF(BNRYRTIOT) THEN
                IF( (J1<=BNRYA .AND. J2<=BNRYA) .OR. (J1>BNRYA .AND. J2>BNRYA) ) THEN
                    IF(SQSHT .AND. SQID==1 .AND. DIST <= SQSHDEL) ENERGY = ENERGY + SQSHEPS
                    RETURN
                ENDIF
            ELSE
                IF(BINARYID==0) THEN
                    IF( (MOD(J1,2) == 0 .AND. MOD(J2,2) == 0) .OR. (MOD(J1,2) == 1 .AND. MOD(J2,2) == 1) ) THEN
                        IF(SQSHT .AND. SQID==1 .AND. DIST <= SQSHDEL) ENERGY = ENERGY + SQSHEPS
                        RETURN
                    ENDIF
                ELSE
                    IF( (MOD(J1,2) == 0 .AND. MOD(J2,2) == 1) .OR. (MOD(J1,2) == 1 .AND. MOD(J2,2) == 0) ) THEN
                        IF(SQSHT .AND. SQID==1 .AND. DIST <= SQSHDEL) ENERGY = ENERGY + SQSHEPS
                        RETURN
                    ENDIF
                ENDIF
            ENDIF
        ENDIF

        IF(SQSHT .AND. SQID==0 .AND. DIST <= SQSHDEL) ENERGY = ENERGY + SQSHEPS

        DO J3 = 1, NSITES
        !   Direction of patch alpha on particle I
            EA  = RBSITES(:,J3,J1)
            EARIJ = -DOT_PRODUCT(EA,RHAT)
        !   If normalised distance vector doesn't pass through patch alpha
        !   the conditions for bonding are not met so no need to progress.
            IF(EARIJ <= KFDEL(J3)) CYCLE
            
            DO J4 = 1, NSITES

                IF(DNAT) THEN
                    IF(J3 /= J4) CYCLE
                ENDIF

            !   Direction of patch beta on particle J
                EB    = RBSITES(:,J4,J2)
                EBRJI = DOT_PRODUCT(EB,RHAT)
            !   Check if the normalised distance vector passes through patch beta
            !   and that the patches are close enough to interact.
                IF(EBRJI <= KFDEL(J4) .OR. RIJSQ > KFLAM2(J3,J4) ) THEN
                    CYCLE
                ELSE

                !   The conditions for bonding are met so calculate the contribution 
                !   of the interaction between patches alpha and beta to the energy.
                    ENERGY = ENERGY - KFIJ(J3,J4)
                    IF(ONEBONDT) THEN
                        NPBONDS(J3) = NPBONDS(J3) + 1
                    ENDIF

                !   If performing MC with cluster-moves, add particle J to the current 
                !   cluster (if it is not already in the cluster). 
                    IF(CLUSTERT .OR. VLMCLUSTERMOVET) THEN
                    !   If simulating triblock patchy particles then we only consider particles
                    !   interacting via patch B-patch B bonds to be apart of the same cluster.
                        IF( .NOT. LRGCLSTMVT ) THEN
                            IF( J3 /= CLSTRSITEID .OR. J4 /= CLSTRSITEID ) CYCLE ! DO NOTHING ELSE
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

                    IF(PRNTCNF .OR. BOPCLSTRT) THEN
                        IF( J3 /= CLSTRSITEID .OR. J4 /= CLSTRSITEID ) THEN
                            CLSTRADJ(J1,J2) = -1
                            CLSTRADJ(J2,J1) = -1
                        ELSE
                            CLSTRADJ(J1,J2) = 1
                            CLSTRADJ(J2,J1) = 1
                            CLURIJ(:,J1,J2) = RIJ
                            CLURIJ(:,J2,J1) = -RIJ
                        ENDIF
                    ENDIF

                ENDIF
            ENDDO ! Loop over each of the patches on particle j
        ENDDO ! Loop over each of the patches on particle i
    
    END SUBROUTINE KF

    SUBROUTINE KF_SPHERICAL_SURF(ENERGY, J1)
!==============================================================================================
!   Subroutine to calculate the potential energy for a system of patchy particles interacting
!   via the Kern-Frenkel potential.
!==============================================================================================
        USE COMMONS, ONLY: DP, NDIM, RBSITES, NSITES, R, KFDEL, KFIJ, SPHERERAD
    
        IMPLICIT NONE
    
        INTEGER, INTENT(IN)         :: J1

        INTEGER                     :: J3
        REAL(KIND=DP)               :: RI(NDIM), RHAT(NDIM), R2
        REAL(KIND=DP)               :: EA(NDIM), EARIJ
    
        REAL(KIND=DP), INTENT(OUT)  :: ENERGY

        ENERGY   = 0.0_dp

        RI = R(:,J1)
        R2 = NORM2(RI)

        IF (ABS(SPHERERAD-R2) <= 0.7_dp) THEN
            RHAT = RI / R2
        ELSE
            RETURN
        ENDIF

        DO J3 = 1, NSITES
        !   Direction of patch alpha on particle I
            EA  = RBSITES(:,J3,J1)
            EARIJ = DOT_PRODUCT(EA,RHAT)
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
    
    END SUBROUTINE KF_SPHERICAL_SURF

!====================================================================================================
    
    SUBROUTINE DEF_KF()
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
    
        USE COMMONS, ONLY: DP, PI, NSITES, REFSITE, KFIJ, KFLAM, KFLAM2, KFDEL
        USE COMMONS, ONLY: KFAA, KFBB, KFCC, KFDD, KFEE, KFFF
        USE COMMONS, ONLY: KFLAMA, KFLAMB, KFLAMC, KFLAMD, KFLAME, KFLAMF
        USE COMMONS, ONLY: KFDELA, KFDELB, KFDELC, KFDELD, KFDELE, KFDELF

        IMPLICIT NONE

        ALLOCATE( KFIJ(NSITES,NSITES), KFLAM(NSITES,NSITES), KFLAM2(NSITES,NSITES), KFDEL(NSITES) )

        IF(NSITES == 1) KFAA = 1.0_dp
        KFIJ(1,1)   = KFAA
        KFLAM2(1,1) = KFLAMA**2

        IF(NSITES == 1) THEN
        ! Janus patchy particles 
            REFSITE(:,1) = [0.0_dp, 0.0_dp,  1.0_dp]
            KFDEL(1)     = COS(KFDELA*PI/180_dp)

        ELSEIF(NSITES == 2) THEN
        ! Triblock patchy particles
            KFIJ(1,2) = SQRT(KFAA*KFBB)
            KFIJ(2,1) = KFIJ(1,2)
            KFIJ(2,2) = KFBB

            KFLAM2(1,2) = ( (KFLAMA + KFLAMB) / 2.0_dp )**2
            KFLAM2(2,1) = ( (KFLAMB + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,2) = KFLAMB**2

            REFSITE(:,1)= [0.0_dp, 0.0_dp, -1.0_dp]
            REFSITE(:,2)= [0.0_dp, 0.0_dp,  1.0_dp]

            KFDEL(1)  = COS(KFDELA*PI/180_dp)
            KFDEL(2)  = COS(KFDELB*PI/180_dp)

        ELSEIF(NSITES == 3) THEN
        ! Trivalent patchy particles with patches in a trigonal planar geometry
            REFSITE(:,1)= [ 0.0_dp,   0.0_dp,                 1.0_dp ]
            REFSITE(:,2)= [ 0.0_dp,   (SQRT(3.0_dp)/2.0_dp), -0.5_dp ]
            REFSITE(:,3)= [ 0.0_dp,  -(SQRT(3.0_dp)/2.0_dp), -0.5_dp ]

            KFDEL(1)  = COS(KFDELA*PI/180_dp)
            KFDEL(2)  = COS(KFDELB*PI/180_dp)
            KFDEL(3)  = COS(KFDELC*PI/180_dp)

            KFIJ(1,2) = SQRT(KFAA*KFBB)
            KFIJ(2,1) = KFIJ(1,2)
            KFIJ(2,2) = KFBB
            KFIJ(1,3) = SQRT(KFAA*KFCC)
            KFIJ(3,1) = KFIJ(1,3)
            KFIJ(2,3) = SQRT(KFBB*KFCC)
            KFIJ(3,2) = KFIJ(2,3)
            KFIJ(3,3) = KFCC

            KFLAM2(1,2) = ( (KFLAMA + KFLAMB) / 2.0_dp )**2
            KFLAM2(2,1) = ( (KFLAMB + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,2) = KFLAMB**2
            KFLAM2(1,3) = ( (KFLAMA + KFLAMC) / 2.0_dp )**2
            KFLAM2(3,1) = ( (KFLAMC + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,3) = ( (KFLAMB + KFLAMC) / 2.0_dp )**2
            KFLAM2(3,2) = ( (KFLAMC + KFLAMB) / 2.0_dp )**2
            KFLAM2(3,3) = KFLAMC**2
        
        ELSEIF(NSITES == 4) THEN
        ! Tetrahedral patchy particles
            REFSITE(:,1)= [  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
            REFSITE(:,2)= [ -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
            REFSITE(:,3)= [  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]
            REFSITE(:,4)= [ -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]

            KFDEL(1)  = COS(KFDELA*PI/180_dp)
            KFDEL(2)  = COS(KFDELB*PI/180_dp)
            KFDEL(3)  = COS(KFDELC*PI/180_dp)
            KFDEL(4)  = COS(KFDELD*PI/180_dp)

            KFIJ(1,2) = SQRT(KFAA*KFBB)
            KFIJ(2,1) = KFIJ(1,2)
            KFIJ(2,2) = KFBB
            KFIJ(1,3) = SQRT(KFAA*KFCC)
            KFIJ(3,1) = KFIJ(1,3)
            KFIJ(2,3) = SQRT(KFBB*KFCC)
            KFIJ(3,2) = KFIJ(2,3)
            KFIJ(3,3) = KFCC
            KFIJ(1,4) = SQRT(KFAA*KFDD)
            KFIJ(4,1) = KFIJ(1,4)
            KFIJ(2,4) = SQRT(KFBB*KFDD)
            KFIJ(4,2) = KFIJ(2,4)
            KFIJ(3,4) = SQRT(KFCC*KFDD)
            KFIJ(4,3) = KFIJ(3,4)
            KFIJ(4,4) = KFDD

            KFLAM2(1,2) = ( (KFLAMA + KFLAMB) / 2.0_dp )**2
            KFLAM2(2,1) = ( (KFLAMB + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,2) = KFLAMB**2
            KFLAM2(1,3) = ( (KFLAMA + KFLAMC) / 2.0_dp )**2
            KFLAM2(3,1) = ( (KFLAMC + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,3) = ( (KFLAMB + KFLAMC) / 2.0_dp )**2
            KFLAM2(3,2) = ( (KFLAMC + KFLAMB) / 2.0_dp )**2
            KFLAM2(3,3) = KFLAMC**2
            KFLAM2(1,4) = ( (KFLAMA + KFLAMD) / 2.0_dp )**2
            KFLAM2(4,1) = ( (KFLAMD + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,4) = ( (KFLAMB + KFLAMD) / 2.0_dp )**2
            KFLAM2(4,2) = ( (KFLAMD + KFLAMB) / 2.0_dp )**2
            KFLAM2(3,4) = ( (KFLAMC + KFLAMD) / 2.0_dp )**2
            KFLAM2(4,3) = ( (KFLAMD + KFLAMC) / 2.0_dp )**2
            KFLAM2(4,4) = KFLAMD**2

        ELSEIF(NSITES == 6) THEN
            REFSITE(:,1)= (/ 1.0_dp,  0.0_dp,  0.0_dp /)
            REFSITE(:,2)= (/ 0.0_dp,  1.0_dp,  0.0_dp /)
            REFSITE(:,3)= (/ 0.0_dp,  0.0_dp,  1.0_dp /)
            REFSITE(:,4)= (/-1.0_dp,  0.0_dp,  0.0_dp /)
            REFSITE(:,5)= (/ 0.0_dp, -1.0_dp,  0.0_dp /)
            REFSITE(:,6)= (/ 0.0_dp,  0.0_dp, -1.0_dp /)

            KFDEL(1)  = COS(KFDELA*PI/180_dp)
            KFDEL(2)  = COS(KFDELB*PI/180_dp)
            KFDEL(3)  = COS(KFDELC*PI/180_dp)
            KFDEL(4)  = COS(KFDELD*PI/180_dp)
            KFDEL(5)  = COS(KFDELE*PI/180_dp)
            KFDEL(6)  = COS(KFDELF*PI/180_dp)

            KFIJ(1,2) = SQRT(KFAA*KFBB); KFIJ(2,1) = KFIJ(1,2)
            KFIJ(2,2) = KFBB
            KFIJ(1,3) = SQRT(KFAA*KFCC); KFIJ(3,1) = KFIJ(1,3)
            KFIJ(2,3) = SQRT(KFBB*KFCC); KFIJ(3,2) = KFIJ(2,3)
            KFIJ(3,3) = KFCC
            KFIJ(1,4) = SQRT(KFAA*KFDD); KFIJ(4,1) = KFIJ(1,4)
            KFIJ(2,4) = SQRT(KFBB*KFDD); KFIJ(4,2) = KFIJ(2,4)
            KFIJ(3,4) = SQRT(KFCC*KFDD); KFIJ(4,3) = KFIJ(3,4)
            KFIJ(4,4) = KFDD
            KFIJ(1,5) = SQRT(KFAA*KFEE); KFIJ(5,1) = KFIJ(1,5)
            KFIJ(2,5) = SQRT(KFBB*KFEE); KFIJ(5,2) = KFIJ(2,5)
            KFIJ(3,5) = SQRT(KFCC*KFEE); KFIJ(5,3) = KFIJ(3,5)
            KFIJ(4,5) = SQRT(KFDD*KFEE); KFIJ(5,4) = KFIJ(4,5)
            KFIJ(5,5) = KFEE
            KFIJ(1,6) = SQRT(KFAA*KFFF); KFIJ(6,1) = KFIJ(1,6)
            KFIJ(2,6) = SQRT(KFBB*KFFF); KFIJ(6,2) = KFIJ(2,6)
            KFIJ(3,6) = SQRT(KFCC*KFFF); KFIJ(6,3) = KFIJ(3,6)
            KFIJ(4,6) = SQRT(KFDD*KFFF); KFIJ(6,4) = KFIJ(4,6)
            KFIJ(5,6) = SQRT(KFEE*KFFF); KFIJ(6,5) = KFIJ(5,6)
            KFIJ(6,6) = KFFF

            KFLAM2(1,2) = ( (KFLAMA + KFLAMB) / 2.0_dp )**2
            KFLAM2(2,1) = ( (KFLAMB + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,2) = KFLAMB**2
            KFLAM2(1,3) = ( (KFLAMA + KFLAMC) / 2.0_dp )**2
            KFLAM2(3,1) = ( (KFLAMC + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,3) = ( (KFLAMB + KFLAMC) / 2.0_dp )**2
            KFLAM2(3,2) = ( (KFLAMC + KFLAMB) / 2.0_dp )**2
            KFLAM2(3,3) = KFLAMC**2
            KFLAM2(1,4) = ( (KFLAMA + KFLAMD) / 2.0_dp )**2
            KFLAM2(4,1) = ( (KFLAMD + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,4) = ( (KFLAMB + KFLAMD) / 2.0_dp )**2
            KFLAM2(4,2) = ( (KFLAMD + KFLAMB) / 2.0_dp )**2
            KFLAM2(3,4) = ( (KFLAMC + KFLAMD) / 2.0_dp )**2
            KFLAM2(4,3) = ( (KFLAMD + KFLAMC) / 2.0_dp )**2
            KFLAM2(4,4) = KFLAMD**2
            KFLAM2(1,5) = ( (KFLAMA + KFLAME) / 2.0_dp )**2
            KFLAM2(5,1) = ( (KFLAME + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,5) = ( (KFLAMB + KFLAME) / 2.0_dp )**2
            KFLAM2(5,2) = ( (KFLAME + KFLAMB) / 2.0_dp )**2
            KFLAM2(3,5) = ( (KFLAMC + KFLAME) / 2.0_dp )**2
            KFLAM2(5,3) = ( (KFLAME + KFLAMC) / 2.0_dp )**2
            KFLAM2(4,5) = ( (KFLAMD + KFLAME) / 2.0_dp )**2
            KFLAM2(5,4) = ( (KFLAME + KFLAMD) / 2.0_dp )**2
            KFLAM2(5,5) = KFLAME**2
            KFLAM2(1,6) = ( (KFLAMA + KFLAMF) / 2.0_dp )**2
            KFLAM2(6,1) = ( (KFLAMF + KFLAMA) / 2.0_dp )**2
            KFLAM2(2,6) = ( (KFLAMB + KFLAMF) / 2.0_dp )**2
            KFLAM2(6,2) = ( (KFLAMF + KFLAMB) / 2.0_dp )**2
            KFLAM2(3,6) = ( (KFLAMC + KFLAMF) / 2.0_dp )**2
            KFLAM2(6,3) = ( (KFLAMF + KFLAMC) / 2.0_dp )**2
            KFLAM2(4,6) = ( (KFLAMD + KFLAMF) / 2.0_dp )**2
            KFLAM2(6,4) = ( (KFLAMF + KFLAMD) / 2.0_dp )**2
            KFLAM2(5,6) = ( (KFLAME + KFLAMF) / 2.0_dp )**2
            KFLAM2(6,5) = ( (KFLAMF + KFLAME) / 2.0_dp )**2
            KFLAM2(6,6) = KFLAMF**2
        ENDIF

        KFLAM = SQRT(KFLAM2)
        
    END SUBROUTINE
    
!====================================================================================================
    
    SUBROUTINE VIEW_KF()
    
        USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT, BINARYT, BNRYRTIOT, BNRYA
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
            IF(BINARYT) THEN 
                IF(BNRYRTIOT) THEN
                    IF( J1<=BNRYA ) THEN
                        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
                    ELSE
                        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Co ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
                    ENDIF
                ELSE
                    IF( MOD(J1,2) == 0 ) THEN
                        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Co ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
                    ELSE
                        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
                    ENDIF
                ENDIF
            ELSE
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            ENDIF

            RM = Q_TO_RM( Q(:,J1) )
            DO J2 = 1, NSITES
                RBCOORDS = RWRITE(:,J1) + 0.35*MATMUL(RM,REFSITE(:,J2))
                IF(J2==1)THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==2) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==3) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'K ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==4) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'B ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==5) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'P ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSE
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Ni ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ENDIF
            ENDDO
        END DO
        
    END SUBROUTINE