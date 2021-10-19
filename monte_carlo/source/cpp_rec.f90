SUBROUTINE KF_REC(ENERGY, RIJ, RIJSQ, J1, J2)
    !==============================================================================================
    !   Subroutine to calculate the potential energy for a system of patchy particles interacting
    !   via the Kern-Frenkel potential.
    !==============================================================================================
            USE COMMONS, ONLY: DP, NDIM, RBSITES, NSITES
            USE COMMONS, ONLY: KFDEL, KFMDEL, KFIJ, CPPINVS, PIS, CPPLAM, YUKKAP
            USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, CLSTRID
        
            IMPLICIT NONE
        
            INTEGER, INTENT(IN)         :: J1, J2
            REAL(KIND=DP), INTENT(IN)   :: RIJ(NDIM), RIJSQ
    
            INTEGER                     :: J3, J4
            REAL(KIND=DP)               :: EA1(NDIM), EA2(NDIM), EA3(NDIM)
            REAL(KIND=DP)               :: EARIJ1(NDIM), EARIJ2(NDIM), EARIJ3(NDIM)
            REAL(KIND=DP)               :: RI12(NDIM), RI13(NDIM), EARI12, EARI13
            REAL(KIND=DP)               :: EB1(NDIM), EB2(NDIM), EB3(NDIM)
            REAL(KIND=DP)               :: EBRJI1(NDIM), EBRJI2(NDIM), EBRJI3(NDIM)
            REAL(KIND=DP)               :: RJ12(NDIM), RJ13(NDIM), EBRJ12, EBRJ13
            REAL(KIND=DP)               :: PHIA1, PHIA2, PHIB1, PHIB2
            REAL(KIND=DP)               :: DIST, UYUK, DLAM, WPP
        
            REAL(KIND=DP), INTENT(OUT)  :: ENERGY
    
            ENERGY   = 0.0_dp
        
            DIST = SQRT(RIJSQ)
        !   Calculate repulsive contribution to the pair energy
            UYUK   = EXP( -YUKKAP*(DIST - 1.0_dp) )
            ENERGY = UYUK / DIST
        !   Check if particles are close enough for patches to interact

            DLAM = DIST - CPPLAM

            IF(DLAM <= CPPINVS) THEN
                WPP = 0.0_dp
            !   Distance dependence of the patch-patch interactions
                IF(DLAM < 0.0_dp) THEN
                    WPP   = 1.0_dp
                ELSE
                    WPP   = 0.5_dp*(1.0_dp + COS(PIS*DLAM))
                ENDIF
    
                DO J3 = 1, 2
                !   Direction of patch alpha on particle I
                    EA1  = RBSITES(:,(J3-1)*3+1,J1)
                !   Direction of long-axis of rectangular patch
                    EA2  = RBSITES(:,(J3-1)*3+2,J1)
                !   Direction of short-axis of the rectangular patch
                    EA3  = RBSITES(:,(J3-1)*3+3,J1)
                !   Project the distance vector onto the axes defining the coordinate frame of the patch
                    EARIJ1 = -DOT_PRODUCT(EA1,RIJ)*EA1
                    EARIJ2 = -DOT_PRODUCT(EA2,RIJ)*EA2
                    EARIJ3 = -DOT_PRODUCT(EA3,RIJ)*EA3
                !   Calculate the projection of the distance vector onto the planes defined by axes 1/2 and 1/3.
                    RI12 =  EARIJ1 + EARIJ2
                    RI13 =  EARIJ1 + EARIJ3
                !   Calculate the cosine of the angles between the projections and the direction vector of the patch
                    EARI12 = DOT_PRODUCT(RI12/NORM2(RI12), EA1)
                    EARI13 = DOT_PRODUCT(RI13/NORM2(RI13), EA1)
                !   If normalised distance vector doesn't pass through patch A
                !   the conditions for bonding are not met so no need to progress.
                    IF(EARI12 <= KFDEL((J3-1)*2+1) .OR. EARI13 <= KFDEL((J3-1)*2+2)) THEN
                        CYCLE
                    ELSE
                        PHIA1 = 1.0_dp - COS( KFMDEL((J3-1)*2+1) * (EARI12 - KFDEL((J3-1)*2+1)) )
                        PHIA2 = 1.0_dp - COS( KFMDEL((J3-1)*2+2) * (EARI13 - KFDEL((J3-1)*2+2)) )
                    ENDIF
                    
                    DO J4 = 1, 2
                    !   Direction of patch beta on particle J
                        EB1  = RBSITES(:,(J4-1)*3+1,J2)
                    !   Direction of long-axis of rectangular patch
                        EB2  = RBSITES(:,(J4-1)*3+2,J2)
                    !   Direction of short-axis of the rectangular patch
                        EB3  = RBSITES(:,(J4-1)*3+3,J2)
                    !   Project the distance vector onto the axes defining the coordinate frame of the patch
                        EBRJI1 = DOT_PRODUCT(EB1,RIJ)*EB1
                        EBRJI2 = DOT_PRODUCT(EB2,RIJ)*EB2
                        EBRJI3 = DOT_PRODUCT(EB3,RIJ)*EB3
                    !   Calculate the projection of the distance vector onto the planes defined by axes 1/2 and 1/3.
                        RJ12 =  EBRJI1 + EBRJI2
                        RJ13 =  EBRJI1 + EBRJI3
                    !   Calculate the cosine of the angles between the projections and the direction vector of the patch
                        EBRJ12 = DOT_PRODUCT(RJ12/NORM2(RJ12), EB1)
                        EBRJ13 = DOT_PRODUCT(RJ13/NORM2(RJ13), EB1)
                    !   If normalised distance vector doesn't pass through patch B
                    !   the conditions for bonding are not met so no need to progress.
                        IF(EBRJ12 <= KFDEL((J4-1)*2+1) .OR. EBRJ13 <= KFDEL((J4-1)*2+2)) THEN
                            CYCLE
                        ELSE
                            PHIB1 = 1.0_dp - COS( KFMDEL((J4-1)*2+1) * (EBRJ12 - KFDEL((J4-1)*2+1)) )
                            PHIB2 = 1.0_dp - COS( KFMDEL((J4-1)*2+2) * (EBRJ13 - KFDEL((J4-1)*2+2)) )
                        ENDIF
        
                    !   The conditions for bonding are met so calculate the contribution 
                    !   of the interaction between patches A and B to the energy.
                        ENERGY = ENERGY - KFIJ(J3,J4)*WPP*(PHIA1*PHIA2*PHIB1*PHIB2) / 16.0_dp
                    
                    !   If performing MC with cluster-moves, add particle J to the current 
                    !   cluster (if it is not already in the cluster). 
                        IF(CLUSTERT) THEN
                        !   If simulating triblock patchy particles then we only consider particles
                        !   interacting via patch B-patch B bonds to be apart of the same cluster.
                            IF( NSITES == 2 .AND. (J3 == 1 .OR. J4 == 1) ) CYCLE

                        !   Check if the particle is already in the cluster, only need to check 
                        !   particle IDs that come after the current particle in the list as the
                        !   previous check should take care of those that come earlier.
                            IF(.NOT.( ANY(CLSTR(CLSTRID+1:CLSTRSZ)==J2) )) THEN
                                CLSTRSZ = CLSTRSZ + 1
                                CLSTR(CLSTRSZ) = J2
                            ENDIF
                        ENDIF
        
                    ENDDO ! Loop over each of the patches on particle j
                ENDDO ! Loop over each of the patches on particle i
            ENDIF
            
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
        
            USE COMMONS, ONLY: DP, PI, NDIM, NSITES, REFSITE
            USE COMMONS, ONLY: KFDELA1, KFDELA2, KFDELB1, KFDELB2, KFIJ, KFDEL, KFMDEL, KFAA, KFBB, RBSITES, NPART, SKEWAB
            USE COMMONS, ONLY: TANA, TANB, CPPS, CPPINVS, PIS
        
            IMPLICIT NONE

            REAL(KIND=DP) :: SS, CS, RMZ(3,3)
    
            ALLOCATE(KFIJ(2,2), KFDEL(4), KFMDEL(4))
            ALLOCATE(REFSITE(NDIM,NSITES), RBSITES(NDIM,NSITES,NPART))

            CPPINVS    = 1.0_dp / CPPS
            PIS        = PI * CPPS

            KFIJ(1,1) = KFAA
            KFIJ(1,2) = SQRT(KFAA*KFBB)
            KFIJ(2,1) = KFIJ(1,2)
            KFIJ(2,2) = KFBB

            KFDEL(1)  = COS(KFDELA1*PI/180_dp)
            KFDEL(2)  = COS(KFDELA2*PI/180_dp)
            KFDEL(3)  = COS(KFDELB1*PI/180_dp)
            KFDEL(4)  = COS(KFDELB2*PI/180_dp)

            KFMDEL(1) = PI / ( 1.0_dp - KFDEL(1) )
            KFMDEL(2) = PI / ( 1.0_dp - KFDEL(2) )
            KFMDEL(3) = PI / ( 1.0_dp - KFDEL(3) )
            KFMDEL(4) = PI / ( 1.0_dp - KFDEL(4) )

            TANA = TAN(KFDELA1/2.0_dp*PI/180_dp)
            TANB = TAN(KFDELB1/2.0_dp*PI/180_dp)
            
            SS = SIN(SKEWAB*PI/180_dp)
            CS = COS(SKEWAB*PI/180_dp)
            RMZ(1,:) = [  CS,    -SS,   0.0_dp]
            RMZ(2,:) = [  SS,     CS,   0.0_dp]
            RMZ(3,:) = [0.0_dp, 0.0_dp, 1.0_dp]

            REFSITE(:,1) = [ 0.0_dp, 0.0_dp, -1.0_dp]
            REFSITE(:,2) = [ 0.0_dp,-1.0_dp,  0.0_dp]
            REFSITE(:,3) = [-1.0_dp, 0.0_dp,  0.0_dp]
            REFSITE(:,4) = [0.0_dp, 0.0_dp,  1.0_dp]
            REFSITE(:,5) = MATMUL(RMZ, REFSITE(:,2))
            REFSITE(:,6) = MATMUL(RMZ, REFSITE(:,3))
            
        END SUBROUTINE
        
    !====================================================================================================
        
        SUBROUTINE VIEW_KF_REC()
        
            USE COMMONS, ONLY: DP, NPART, R, RBSITES, BOX, VIEWUNIT, TANA, TANB
        
            IMPLICIT NONE
        
            INTEGER:: J1
            REAL(KIND=DP) :: RWRITE(3,NPART)
        
            WRITE(VIEWUNIT,*) NPART*7
            WRITE(VIEWUNIT,*)
        
            RWRITE(:,:) = R(:,:)
            
            DO J1 = 1, NPART
                RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX) * BOX
                WRITE(VIEWUNIT,'(A5,1X,3F20.10)') 'Au ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
                
                WRITE(VIEWUNIT,'(A5,1X,3F20.10)') 'N ',  RWRITE(:,J1) + 0.37_dp*RBSITES(:,1,J1)
                WRITE(VIEWUNIT,'(A5,1X,3F20.10)') 'N ',  RWRITE(:,J1) + 0.3_dp*(RBSITES(:,1,J1)+TANA*RBSITES(:,2,J1))
                WRITE(VIEWUNIT,'(A5,1X,3F20.10)') 'N ',  RWRITE(:,J1) + 0.3_dp*(RBSITES(:,1,J1)-TANA*RBSITES(:,2,J1))

                WRITE(VIEWUNIT,'(A5,1X,3F20.10)') 'O ',  RWRITE(:,J1) + 0.37_dp*RBSITES(:,4,J1)
                WRITE(VIEWUNIT,'(A5,1X,3F20.10)') 'O ',  RWRITE(:,J1) + 0.3_dp*(RBSITES(:,4,J1)+TANB*RBSITES(:,5,J1))
                WRITE(VIEWUNIT,'(A5,1X,3F20.10)') 'O ',  RWRITE(:,J1) + 0.3_dp*(RBSITES(:,4,J1)-TANB*RBSITES(:,5,J1))
            END DO
            
        END SUBROUTINE