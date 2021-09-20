    SUBROUTINE CPP(ENERGY, RIJ, RIJSQ, J1, J2)
!     ==============================================================================================
!     This subroutine is for particles interacting via a continuous Kern-Frenkel patchy model. 
!     See "Designing a Bernal Spiral from Patchy Colloids", ACS Nano, 7, 1246-1256 (2013) for the
!     patch-patch interactions.
!     ==============================================================================================
        USE COMMONS, ONLY: DP, NDIM, RBSITES, NSITES, VIRTEMP
        USE COMMONS, ONLY: YUKKAP, CPPDEL, CPPMDEL, CPPIJ, CPPLAM, CPPINVS, PIS, DNAT, BINARYT, HLFPART
        USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, CLSTRID, CLSTRADJ, VLMCLUSTERMOVET, CLURIJ, CLSTRSITEID
    
        IMPLICIT NONE

        INTEGER, INTENT(IN)       :: J1, J2
        REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM), RIJSQ

    !   Counters for the loops   
        INTEGER       :: J3, J4
    !   Parameters related to the translational coordinates of the particles         
        REAL(KIND=DP) :: RHAT(NDIM), DIST
    !   Parameters related to the orientational coordinates of the particles        
        REAL(KIND=DP) :: EA(NDIM), EB(NDIM), EARIJ, EBRJI
    !   Parameters specific to potential energy function 
        REAL(KIND=DP) :: UYUK, DLAM, WPP, PHII, PHIJ,  UPP
    !   Parameters related to calculating the virial pressure
        REAL(KIND=DP) :: DWPDR, FPIPJ, DPHIJDR, DPHIIDR

    !   PARAMETER TO BE OUTPUT FROM SUBROUTINE
        REAL(KIND=DP), INTENT(OUT) :: ENERGY
 
        ENERGY   = 0.0_dp

        DIST = SQRT(RIJSQ)
    !   Calculate repulsive contribution to the pair energy
        UYUK    = EXP( -YUKKAP*(DIST - 1.0_dp) )
        ENERGY  = ENERGY + UYUK / DIST
        VIRTEMP = VIRTEMP + UYUK*(YUKKAP + 1.0_dp/DIST)

        !   If considering a binary system then the 1st and 2nd halves of the particle list do not iteract with one another.
        IF(BINARYT) THEN
            IF( (J1 < HLFPART .AND. J2 < HLFPART) .OR. (J1 >= HLFPART .AND. J2 >= HLFPART) ) THEN
                RETURN
            ENDIF
        ENDIF

    !   Check if particles are close enough for patches to interact
        DLAM = DIST - CPPLAM

        IF(DLAM <= CPPINVS) THEN
            WPP = 0.0_dp
        !   Distance dependence of the patch-patch interactions
            IF(DLAM < 0.0_dp) THEN
                WPP   = -1.0_dp
                DWPDR = 0.0_dp
            ELSE
                WPP   = -0.5_dp*(1.0_dp + COS(PIS*DLAM))
                DWPDR = 0.5_dp*PIS*SIN(PIS*DLAM)
            ENDIF

            RHAT  = RIJ / DIST

        !   Angular dependence of the patch-patch interactions
            DO J3 = 1, NSITES
            !   Direction of patch alpha on particle I
                EA     = RBSITES(:,J3,J1)
            !   Angle between the direction of patch alpha and the inter-particle distance vector
                EARIJ  = -DOT_PRODUCT(EA,RHAT)
                
                IF(EARIJ < CPPDEL(J3)) THEN
                    CYCLE   ! There cannot be any patch-patch interaction so cycle 
                ELSE
                    PHII    = 1.0_dp - COS( CPPMDEL(J3) * (EARIJ - CPPDEL(J3)) )
                    DPHIIDR = (1.0_dp/RIJSQ)*(CPPMDEL(J3)*EARIJ)*SIN( (1.0_dp/DIST)*CPPMDEL(J3)*(EARIJ-CPPDEL(J3)) )
                ENDIF

                DO J4 = 1, NSITES

                    IF(DNAT) THEN
                        IF(J3 /= J4) CYCLE
                    ENDIF

                    UPP   = 0.0_dp
                    FPIPJ = 0.0_dp

                !   Direction of patch beta on particle J
                    EB     = RBSITES(:,J4,J2)
                !   Angle between the direction of patch beta and the inter-particle distance vector
                    EBRJI  = DOT_PRODUCT(EB,RHAT)
                    
                    IF(EBRJI < CPPDEL(J4)) THEN
                        CYCLE ! There cannot be any patch-patch interaction so cycle 
                    ELSE
                        PHIJ    = 1.0_dp - COS( CPPMDEL(J4) * (EBRJI - CPPDEL(J4)) )
                        DPHIJDR = -(1.0_dp/RIJSQ)*(CPPMDEL(J4)*EBRJI)*SIN((1.0_dp/DIST)*CPPMDEL(J4)*(EBRJI-CPPDEL(J4)))
                    ENDIF

                    UPP   = UPP + CPPIJ(J3,J4) * PHII * PHIJ
                    FPIPJ = FPIPJ + CPPIJ(J3,J4)*(PHIJ*WPP*DPHIIDR + PHII*WPP*DPHIJDR + PHII*PHIJ*DWPDR)

                    ENERGY  = ENERGY + 0.25_dp * UPP * WPP
                    VIRTEMP = VIRTEMP - 0.25_dp*DIST*FPIPJ

                !   If performing MC with cluster-moves, add particle J to the current 
                !   cluster (if it is not already in the cluster). 
                    IF(CLUSTERT .OR. VLMCLUSTERMOVET) THEN
                    !   If simulating triblock patchy particles then we only consider particles
                    !   interacting via patch B-patch B bonds to be apart of the same cluster.
                        IF( J3 /= CLSTRSITEID .OR. J4 /= CLSTRSITEID ) CYCLE
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

                ENDDO ! Loop over each of the patches on particle j
            ENDDO ! Loop over each of the patches on particle i

        ENDIF ! Check if particles are within range for patch-patch interactions
    
    END SUBROUTINE CPP

!====================================================================================================
    
    SUBROUTINE DEF_CPP()
    !----------------------------------------------------------------
    ! Define the position of the patches on particles interacting via
    ! the continuous patchy potential.
    ! Also define constants to be used during the simulations which 
    ! depend on the width of the patches and the well-depth of the 
    ! various patch-patch interactions.
    ! The well-depth of patch A is used to define the reduced unit
    ! of temperature and so CPPAA=1 for all cases.
    ! Currently, for a given number of patches, the geometry of the
    ! patches is fixed.
    !----------------------------------------------------------------
    
        USE COMMONS, ONLY: DP, NSITES, REFSITE, CPPDELA, CPPDELB, CPPDELC, CPPDELD, PI, CPPIJ, CPPDEL
        USE COMMONS, ONLY: CPPAA, CPPBB, CPPCC, CPPDD, CPPMDEL, CPPS, CPPINVS, PIS
    
        IMPLICIT NONE

        ALLOCATE( CPPIJ(NSITES,NSITES), CPPDEL(NSITES), CPPMDEL(NSITES) )

        CPPINVS    = 1.0_dp / CPPS
        PIS        = PI * CPPS

        CPPAA      = 1.0_dp
        CPPIJ(1,1) = CPPAA
        CPPIJ(1,2) = SQRT(CPPAA*CPPBB)
        CPPIJ(2,1) = CPPIJ(1,2)
        CPPIJ(2,2) = CPPBB

        ! Triblock patchy particles
        IF(NSITES == 2) THEN
            REFSITE(:,1)= [0.0_dp, 0.0_dp, -1.0_dp]
            REFSITE(:,2)= [0.0_dp, 0.0_dp,  1.0_dp]

            CPPDEL(1)  = COS(CPPDELA*PI/180.0_dp)
            CPPDEL(2)  = COS(CPPDELB*PI/180.0_dp)

            CPPMDEL(1) = PI / ( 1.0_dp - CPPDEL(1) )
            CPPMDEL(2) = PI / ( 1.0_dp - CPPDEL(2) )

        ! Trivalent patchy particles with patches in a trigonal planar geometry
        ELSEIF(NSITES == 3) THEN

            REFSITE(:,1) = [ 0.0_dp,   0.0_dp,                 1.0_dp ]
            REFSITE(:,2) = [ 0.0_dp,   (SQRT(3.0_dp)/2.0_dp), -0.5_dp ]
            REFSITE(:,3) = [ 0.0_dp,  -(SQRT(3.0_dp)/2.0_dp), -0.5_dp ]

            CPPDEL(1)  = COS(CPPDELA*PI/180.0_dp)
            CPPDEL(2)  = COS(CPPDELB*PI/180.0_dp)
            CPPDEL(3)  = COS(CPPDELC*PI/180.0_dp)

            CPPMDEL(1) = PI / ( 1.0_dp - CPPDEL(1) )
            CPPMDEL(2) = PI / ( 1.0_dp - CPPDEL(2) )
            CPPMDEL(3) = PI / ( 1.0_dp - CPPDEL(3) )

            CPPIJ(1,3) = SQRT(CPPAA*CPPCC)
            CPPIJ(3,1) = CPPIJ(1,3)
            CPPIJ(2,3) = SQRT(CPPBB*CPPCC)
            CPPIJ(3,2) = CPPIJ(2,3)
            CPPIJ(3,3) = CPPCC
        
        ! Tetrahedral patchy particles
        ELSEIF(NSITES == 4) THEN

            REFSITE(:,2) = [  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
            REFSITE(:,1) = [ -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) ]
            REFSITE(:,3) = [  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]
            REFSITE(:,4) = [ -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) ]

            CPPDEL(1)    = COS(CPPDELA*PI/180.0_dp)
            CPPDEL(2)    = COS(CPPDELB*PI/180.0_dp)
            CPPDEL(3)    = COS(CPPDELC*PI/180.0_dp)
            CPPDEL(4)    = COS(CPPDELD*PI/180.0_dp)

            CPPMDEL(1)   = PI / ( 1.0_dp - CPPDEL(1) )
            CPPMDEL(2)   = PI / ( 1.0_dp - CPPDEL(2) )
            CPPMDEL(3)   = PI / ( 1.0_dp - CPPDEL(3) )
            CPPMDEL(4)   = PI / ( 1.0_dp - CPPDEL(4) )

            CPPIJ(1,3)   = SQRT(CPPAA*CPPCC)
            CPPIJ(3,1)   = CPPIJ(1,3)
            CPPIJ(2,3)   = SQRT(CPPBB*CPPCC)
            CPPIJ(3,2)   = CPPIJ(2,3)
            CPPIJ(3,3)   = CPPCC
            CPPIJ(1,4)   = SQRT(CPPAA*CPPDD)
            CPPIJ(4,1)   = CPPIJ(1,4)
            CPPIJ(2,4)   = SQRT(CPPBB*CPPDD)
            CPPIJ(4,2)   = CPPIJ(2,4)
            CPPIJ(3,4)   = SQRT(CPPCC*CPPDD)
            CPPIJ(4,3)   = CPPIJ(3,4)
            CPPIJ(4,4)   = CPPDD

        ENDIF
        
    END SUBROUTINE
    
!====================================================================================================
    
    SUBROUTINE VIEW_CPP()
    
        USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT, BINARYT, HLFPART
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
            IF(BINARYT .AND. (J1 >= HLFPART)) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Co ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            ELSE
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            ENDIF
            
            RM = Q_TO_RM( Q(:,J1) )
            DO J2 = 1, NSITES
                RBCOORDS = RWRITE(:,J1) + 0.5_dp*MATMUL(RM,REFSITE(:,J2))
                IF(J2==1)THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==2) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==3) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'K ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSE
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Ni ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ENDIF
            ENDDO
        END DO
        
    END SUBROUTINE

!====================================================================================================
    
    SUBROUTINE PESRF_CPP()
!   =================================================================================================
!    Subroutine to calculate and output the pair potential for the CPP model, written specifically
!    for the triblock patchy particles (i.e., particles with two patches).
!    Outputs both the distance and angular dependence of the potential.
!   =================================================================================================
        USE COMMONS, ONLY: DP, NPART, NDIM, R, Q, NSITES, RBSITES, REFSITE, RCUT, PI, SPET, CELLLISTT, CPPLAM
        USE ROTATIONS_MODULE, ONLY: QUATMUL, Q_TO_RM
    
        IMPLICIT NONE
    
        INTEGER       :: J1, J2, J3, J4, J5, MAXR, MAXTHT, NPARTO
        REAL(KIND=DP) :: DELR, DELTHT, STPR, STPTHT, PE
        REAL(KIND=DP) :: RIJ(3), RIJSQ, RMOVEMAP(4), RRM(3,3), RM(3,3), QOLD(4)

        IF(ALLOCATED(R)) DEALLOCATE(R)
        IF(ALLOCATED(Q)) DEALLOCATE(Q)
        IF(ALLOCATED(RBSITES)) DEALLOCATE(RBSITES)   
        
        ALLOCATE(R(NDIM,2), Q(4,2), RBSITES(NDIM,NSITES,2) )

        SPET      = .FALSE.
        CELLLISTT = .FALSE.
        
        NPARTO = NPART
        NPART  = 2

        STPR   = 0.001_dp
        STPTHT = 0.1_dp
        
        MAXR   = INT(RCUT / STPR)
        MAXTHT = INT(90.0 / STPTHT) 

        R = 0.0_dp

        OPEN (UNIT = 81,  FILE = 'pe_ab.dat', STATUS = 'UNKNOWN')
        OPEN (UNIT = 82,  FILE = 'pe_aa.dat', STATUS = 'UNKNOWN')
        OPEN (UNIT = 83,  FILE = 'pe_bb.dat', STATUS = 'UNKNOWN')
        OPEN (UNIT = 84,  FILE = 'ang_ab.dat', STATUS = 'UNKNOWN')
        OPEN (UNIT = 85,  FILE = 'ang_aa.dat', STATUS = 'UNKNOWN')
        OPEN (UNIT = 86,  FILE = 'ang_bb.dat', STATUS = 'UNKNOWN')

        DO J3 = 1, 3
            DO J1 = 1, MAXR
                
                DELR = REAL(J1,DP) * STPR

                R(:,2) = [0.0_dp, 0.0_dp, DELR]

                IF(J3 == 1) THEN
                    Q(:,1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                    Q(:,2) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                ELSE IF(J3 == 2) THEN
                    Q(:,1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                    Q(:,2) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
                ELSE
                    Q(:,1) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
                    Q(:,2) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                ENDIF

                DO J4 = 1, NPART
                    RM = Q_TO_RM( Q(:,J4) )
                    DO J5 = 1, NSITES
                        RBSITES(:,J5,J4) = MATMUL(RM,REFSITE(:,J5))
                    ENDDO
                ENDDO
                
                RIJ = R(:,1) - R(:,2)
                RIJSQ = DOT_PRODUCT(RIJ,RIJ)
                CALL PAIR_POTENTIAL(PE, RIJ, RIJSQ, 1, 2)

                IF(PE > HUGE(DP) .OR. ISNAN(PE)) CYCLE

                IF(J3 == 1) THEN
                    WRITE (81,*) DELR, PE
                    IF(DELR == CPPLAM) WRITE (84,*) 0.0, PE
                ELSE IF(J3 == 2) THEN
                    WRITE (82,*) DELR, PE
                    IF(DELR == CPPLAM) WRITE (85,*) 0.0, PE
                ELSE
                    WRITE (83,*) DELR, PE
                    IF(DELR == CPPLAM) WRITE (86,*) 0.0, PE
                ENDIF
        
                DO J2 = 1, MAXTHT
                    
                    R(:,2) = [0.0_dp, 0.0_dp, DELR]
                    DELTHT = REAL(J2,DP) * (PI/180.0_dp) * STPTHT
                    
                    IF(J3 == 1) THEN
                        Q(:,1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                        Q(:,2) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                    ELSE IF(J3 == 2) THEN
                        Q(:,1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                        Q(:,2) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
                    ELSE
                        Q(:,1) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]
                        Q(:,2) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
                    ENDIF
                    
                    ! Standard formula for rotation quaternion, using half angles
                    RMOVEMAP(1)   = COS(0.5_dp*DELTHT)
                    RMOVEMAP(2:4) = SIN(0.5_dp*DELTHT)*[0.0_dp, 1.0_dp, 0.0_dp]

                !   Extract the rotation matrix associated with the the rotational move
                    RRM    = Q_TO_RM(RMOVEMAP)
                    R(:,2) = MATMUL(RRM, R(:,2))
                
                !   Save old orientations
                    QOLD = Q(:,2)
                !   Update the orientation of the particle
                    Q(:,2) = QUATMUL(RMOVEMAP, QOLD)

                    DO J4 = 1, NPART
                        RM = Q_TO_RM( Q(:,J4) )
                        DO J5 = 1, NSITES
                            RBSITES(:,J5,J4) = MATMUL(RM,REFSITE(:,J5))
                        ENDDO
                    ENDDO
                    
                    RIJ = R(:,1) - R(:,2)
                    RIJSQ = DOT_PRODUCT(RIJ,RIJ)
                    CALL PAIR_POTENTIAL(PE, RIJ, RIJSQ, 1, 2)

                    IF(DELR == CPPLAM) THEN
                        IF(J3 == 1) THEN
                            WRITE (84,*) DELTHT*(180.0_dp/PI), PE
                        ELSE IF(J3 == 2) THEN
                            WRITE (85,*) DELTHT*(180.0_dp/PI), PE
                        ELSE
                            WRITE (86,*) DELTHT*(180.0_dp/PI), PE
                        ENDIF
                    ENDIF

                ENDDO
            ENDDO
        ENDDO

        CLOSE (UNIT = 81, STATUS = 'KEEP')
        CLOSE (UNIT = 82, STATUS = 'KEEP')
        CLOSE (UNIT = 83, STATUS = 'KEEP')
        CLOSE (UNIT = 84, STATUS = 'KEEP')
        CLOSE (UNIT = 85, STATUS = 'KEEP')
        CLOSE (UNIT = 86, STATUS = 'KEEP')

        STOP
        
    END SUBROUTINE