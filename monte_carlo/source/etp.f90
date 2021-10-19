SUBROUTINE ETP(ENERGY, RIJ, RIJSQ, J1, J2)
    !     ==============================================================================================
    !     This subroutine is for particles interacting via a continuous Kern-Frenkel patchy model. 
    !     See "Designing a Bernal Spiral from Patchy Colloids", ACS Nano, 7, 1246-1256 (2013) for the
    !     patch-patch interactions.
    !     ==============================================================================================
        USE COMMONS, ONLY: DP, NDIM, RBSITES, NSITES, BOX, R, DNAT, POLYDT, ETPDEL, ETPMDEL
        USE COMMONS, ONLY: GLJN, DCHECKSQ, HLFLNGTH, CPPDEL, CPPMDEL, CPPIJ, CPPLAM, CPPINVS, PIS, ETPCUTSQ
        USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, CLSTRID, CLSTRADJ, VLMCLUSTERMOVET, CLURIJ, CLSTRSITEID
    
        IMPLICIT NONE
    
        INTEGER, INTENT(IN)       :: J1, J2
        REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM), RIJSQ
    
    !   Counters for the loops   
        INTEGER       :: J3, J4
    !   Parameters related to the translational coordinates of the particles      
        REAL(KIND=DP) :: RI(NDIM), RJ(NDIM), SDISTSQ
    !   Parameters related to the orientational coordinates of the particles        
        REAL(KIND=DP) :: EA(NDIM), EB(NDIM)
        REAL(KIND=DP) :: EARIJ, EBRJI
    !   Parameters specific to potential energy function 
        REAL(KIND=DP) :: DIJ(NDIM), SD2, SDN, S2DN, VLJ, EK
        REAL(KIND=DP) :: RA(NDIM), RB(NDIM), RAB(NDIM), DAB, DIST, RABN(NDIM)
        REAL(KIND=DP) :: DLAM, WPP, PHII, PHIJ,  UPP, ANGS(NSITES,2), MANGS(NSITES,2)
    
    !   PARAMETER TO BE OUTPUT FROM SUBROUTINE
        REAL(KIND=DP), INTENT(OUT) :: ENERGY

        IF(POLYDT) THEN
            ANGS(:,1)  = ETPDEL(:,J1)   ! COS_T for patches a and b on particle I
            MANGS(:,1) = ETPMDEL(:,J1)  ! PI/(1-COS_T) for particle I
            ANGS(:,2)  = ETPDEL(:,J2)   ! COS_T for patches a and b on particle J
            MANGS(:,2) = ETPMDEL(:,J2)  ! PI/(1-COS_T) for particle J
        ELSE
            ANGS(:,1)  = CPPDEL
            MANGS(:,1) = CPPMDEL
            ANGS(:,2)  = CPPDEL
            MANGS(:,2) = CPPMDEL
        ENDIF
    
        ENERGY   = 0.0_dp
    
        IF(HLFLNGTH /= 0.0_dp) THEN
    !   Calculate the shortest distance between two elongated particles
            CALL SHORTD( RIJ, RBSITES(:,1,J1), RBSITES(:,1,J2), DIJ )
            SDISTSQ = DOT_PRODUCT(DIJ,DIJ)
        ELSE
    !   Particles are spherical so no need to calculate the shortest distance 
            DIJ = RIJ
            SDISTSQ = RIJSQ
        ENDIF
    !   Calculate repulsive Kihara contribution to the pair energy
    !   This depends on the shortest distance between the spherocylinders
        IF(SDISTSQ < DCHECKSQ) THEN
            SD2     = 1.0_dp/SDISTSQ
            SDN     = SD2**(GLJN/2.0_dp)
            S2DN    = SDN**2
            VLJ     = S2DN - SDN   
            EK      = 4.0_dp*(VLJ + 0.25_dp)
            ENERGY  = ENERGY + EK 
    
            IF(ISNAN(ENERGY)) RETURN
        ELSE
    !   Check if the shortest distance between the particles is such that 2.0_dp of the
    !   rigid body sites may be close enough for patch-patch interactions to be present.
            IF(SDISTSQ > ETPCUTSQ) RETURN
        ENDIF
    
        RI = R(:,J1)
        RJ = R(:,J2)
    
    !   Angular dependence of the patch-patch interactions
        DO J3 = 1, NSITES
        !   Direction of patch alpha on particle I
            EA = RBSITES(:,J3,J1)
            RA = RI + HLFLNGTH*EA
    
            DO J4 = 1, NSITES
                IF(DNAT) THEN
                    IF(J3 /= J4) CYCLE
                ENDIF
            !   Direction of patch beta on particle J
                EB = RBSITES(:,J4,J2)
                RB = RJ + HLFLNGTH*EB
    
                RAB = RA - RB
                RAB = RAB - ANINT(RAB/BOX)*BOX
                DAB = DOT_PRODUCT(RAB,RAB)
                DIST  = SQRT(DAB)
                
                DLAM  = DIST - CPPLAM
    
                IF(DLAM <= CPPINVS) THEN
            !   Check if particles are close enough and orientated correctly for patches to interact
                    WPP   = 0.0_dp
                !   Distance dependence of the patch-patch interactions
                    IF(DLAM < 0.0_dp) THEN
                        WPP   = -1.0_dp
                    ELSE
                        WPP   = -0.5_dp*(1.0_dp + COS(PIS*DLAM))
                    ENDIF
    
                    RABN  = RAB / DIST
                    UPP   = 0.0_dp
    
                    !   Angle between the direction of patch alpha and the inter-particle distance vector
                    EARIJ  = -DOT_PRODUCT(EA,RABN)
    
                    IF(EARIJ < ANGS(J3,1)) THEN
                        EXIT   ! There cannot be any patch-patch interaction so cycle to next patch on particle I
                    ELSE
                        PHII    = 1.0_dp - COS( MANGS(J3,1) * (EARIJ - ANGS(J3,1)) )
                    ENDIF
    
                    !   Angle between the direction of patch beta and the inter-particle distance vector
                    EBRJI  = DOT_PRODUCT(EB,RABN)
    
                    IF(EBRJI < ANGS(J4,2)) THEN
                        CYCLE ! There cannot be any patch-patch interaction so cycle to next patch on particle J
                    ELSE
                        PHIJ    = 1.0_dp - COS( MANGS(J4,2) * (EBRJI - ANGS(J4,2)) )
                    ENDIF
                        
                    UPP     = UPP + CPPIJ(J3,J4) * PHII * PHIJ
                    ENERGY  = ENERGY + 0.25_dp * UPP * WPP
    
                !   If performing MC with cluster-moves, add particle J to the current 
                !   cluster (if it is not already in the cluster). 
                    IF( CLUSTERT .OR. VLMCLUSTERMOVET ) THEN
                        IF( (EARIJ >= ANGS(J3,1)) .AND. (EBRJI >= ANGS(J4,2)) ) THEN

                        !   If simulating triblock patchy particles then we only consider particles
                        !   interacting via patch B-patch B bonds to be apart of the same cluster.
                            IF( J3 /= CLSTRSITEID .OR. J4 /= CLSTRSITEID ) CYCLE ! DO NOTHING  
                            
                            IF(VLMCLUSTERMOVET) THEN
                                CLSTRADJ(J1,J2) = 1
                                CLSTRADJ(J2,J1) = 1
                                CLURIJ(:,J1,J2) = RIJ
                                CLURIJ(:,J2,J1) = -RIJ
                            ELSE
                                IF(.NOT.( ANY(CLSTR(CLSTRID+1:CLSTRSZ)==J2) )) THEN
                            !   Check if the particle is already in the cluster, only need to check 
                            !   particle IDs that come after the current particle in the list as the
                            !   previous check should take care of those that come earlier.
                                    CLSTRSZ = CLSTRSZ + 1
                                    CLSTR(CLSTRSZ) = J2
                                ENDIF
                            ENDIF

                        ENDIF
                    ENDIF
    
                ENDIF
            ENDDO ! Loop over each of the patches on particle j
        ENDDO ! Loop over each of the patches on particle i
        
    END SUBROUTINE ETP

    SUBROUTINE HARD_TPR(ENERGY, RIJ, RIJSQ, J1, J2)
    !     ==============================================================================================
    !     This subroutine is for particles interacting via a continuous Kern-Frenkel patchy model. 
    !     See "Designing a Bernal Spiral from Patchy Colloids", ACS Nano, 7, 1246-1256 (2013) for the
    !     patch-patch interactions.
    !     ==============================================================================================
        USE COMMONS, ONLY: DP, NDIM, RBSITES, NSITES, BOX, R, DNAT, OVERLAPT
        USE COMMONS, ONLY: HLFLNGTH, CPPDEL, CPPIJ, CPPLAM, CPPLAM2, POLYDT, ETPDEL
        USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, CLSTRID, CLSTRADJ, VLMCLUSTERMOVET, CLURIJ, CLSTRSITEID, LRGCLSTMVT
    
        IMPLICIT NONE
    
        INTEGER, INTENT(IN)       :: J1, J2
        REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM), RIJSQ
    
    !   Counters for the loops   
        INTEGER       :: J3, J4
    !   Parameters related to the translational coordinates of the particles      
        REAL(KIND=DP) :: RI(NDIM), RJ(NDIM), SDISTSQ
    !   Parameters related to the orientational coordinates of the particles        
        REAL(KIND=DP) :: EA(NDIM), EB(NDIM)
        REAL(KIND=DP) :: EARIJ, EBRJI
    !   Parameters specific to potential energy function 
        REAL(KIND=DP) :: RA(NDIM), RB(NDIM), RAB(NDIM), DIJ(NDIM), DAB, RABN(NDIM), ANGS(NSITES,2)
    
    !   PARAMETER TO BE OUTPUT FROM SUBROUTINE
        REAL(KIND=DP), INTENT(OUT) :: ENERGY

        IF(POLYDT) THEN
            ANGS(:,1)  = ETPDEL(:,J1)   ! COS_T for patches a and b on particle I
            ANGS(:,2)  = ETPDEL(:,J2)   ! COS_T for patches a and b on particle J
        ELSE
            ANGS(:,1)  = CPPDEL
            ANGS(:,2)  = CPPDEL
        ENDIF
    
        ENERGY   = 0.0_dp
    
        IF(HLFLNGTH /= 0.0_dp) THEN
    !   Calculate the shortest distance between two elongated particles
            CALL SHORTD( RIJ, RBSITES(:,1,J1), RBSITES(:,1,J2), DIJ )
            SDISTSQ = DOT_PRODUCT(DIJ,DIJ)
        ELSE
    !   Particles are spherical so no need to calculate the shortest distance 
            DIJ = RIJ
            SDISTSQ = RIJSQ
        ENDIF
    !   Calculate repulsive Kihara contribution to the pair energy
    !   This depends on the shortest distance between the spherocylinders
        IF(SDISTSQ < 1.0_dp) THEN
            OVERLAPT = .TRUE.
            RETURN
        ELSE
    !   Check if the shortest distance between the particles is such that two of the
    !   rigid body sites may be close enough for patch-patch interactions to be present.
            IF(SDISTSQ > CPPLAM2) RETURN
        ENDIF
    
        RI = R(:,J1)
        RJ = R(:,J2)
    
    !   Angular dependence of the patch-patch interactions
        DO J3 = 1, NSITES
        !   Direction of patch alpha on particle I
            EA = RBSITES(:,J3,J1)
            RA = RI + HLFLNGTH*EA
    
            DO J4 = 1, NSITES
                IF(DNAT) THEN
                    IF(J3 /= J4) CYCLE
                ENDIF
            !   Direction of patch beta on particle J
                EB = RBSITES(:,J4,J2)
                RB = RJ + HLFLNGTH*EB
    
                RAB = RA - RB
                RAB = RAB - ANINT(RAB/BOX)*BOX
                DAB = SQRT( DOT_PRODUCT(RAB,RAB) )
    
                IF(DAB <= CPPLAM) THEN
            !   Check if particles are close enough and orientated correctly for patches to interact
    
                    RABN  = RAB / DAB
                    
                    !   Angle between the direction of patch alpha and the inter-particle distance vector
                    EARIJ  = -DOT_PRODUCT(EA,RABN)
                    IF(EARIJ < ANGS(J3,1)) EXIT
    
                    !   Angle between the direction of patch beta and the inter-particle distance vector
                    EBRJI  = DOT_PRODUCT(EB,RABN)
                    IF(EBRJI < ANGS(J4,2)) THEN 
                        CYCLE
                    ELSE
                    
                    !   The conditions for bonding are met so calculate the contribution 
                    !   of the interaction between patches alpha and beta to the energy.
                        ENERGY = ENERGY - CPPIJ(J3,J4)
        
                    !   If performing MC with cluster-moves, add particle J to the current 
                    !   cluster (if it is not already in the cluster). 
                        IF( CLUSTERT .OR. VLMCLUSTERMOVET ) THEN
                            IF( (EARIJ >= ANGS(J3,1)) .AND. (EBRJI >= ANGS(J4,2)) ) THEN

                            !   If simulating triblock patchy particles then we only consider particles
                            !   interacting via patch B-patch B bonds to be apart of the same cluster.
                                IF( .NOT. LRGCLSTMVT ) THEN
                                    IF( J3 /= CLSTRSITEID .OR. J4 /= CLSTRSITEID ) CYCLE ! DO NOTHING  
                                ENDIF
                                
                                IF(VLMCLUSTERMOVET) THEN
                                    CLSTRADJ(J1,J2) = 1
                                    CLSTRADJ(J2,J1) = 1
                                    CLURIJ(:,J1,J2) = RIJ
                                    CLURIJ(:,J2,J1) = -RIJ
                                ELSE
                                    IF(.NOT.( ANY(CLSTR(CLSTRID+1:CLSTRSZ)==J2) )) THEN
                                !   Check if the particle is already in the cluster, only need to check 
                                !   particle IDs that come after the current particle in the list as the
                                !   previous check should take care of those that come earlier.
                                        CLSTRSZ = CLSTRSZ + 1
                                        CLSTR(CLSTRSZ) = J2
                                    ENDIF
                                ENDIF

                            ENDIF
                        ENDIF
                    ENDIF
    
                ENDIF
            ENDDO ! Loop over each of the patches on particle j
        ENDDO ! Loop over each of the patches on particle i
        
    END SUBROUTINE HARD_TPR
    
    !====================================================================================================
        
    SUBROUTINE DEF_ETP()
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
    
        USE COMMONS, ONLY: DP, NPART, NSITES, REFSITE, CPPDELA, CPPDELB, PI, CPPIJ, CPPDEL, RCUT, RCUTSQ, &
        CPPAA, CPPBB, CPPMDEL, CPPS, CPPINVS, PIS, DCHECKSQ, GLJN, ETPCUTSQ, RLNGTH, HLFLNGTH, POLYDT, POLYSTD, ETPDEL, ETPMDEL, &
        READPOLYT

        USE ROTATIONS_MODULE, ONLY: BOX_MULLER
    
        IMPLICIT NONE

        INTEGER       :: J1
        REAL(KIND=DP) :: POLYANGLE
    
        ALLOCATE( CPPIJ(NSITES,NSITES), CPPDEL(NSITES), CPPMDEL(NSITES) )
    
        CPPINVS    = 1.0_dp / CPPS
        PIS        = PI * CPPS
    
        CPPAA      = 1.0_dp
        CPPIJ(1,1) = CPPAA
    
    !   Elongated triblock patchy (ETP) particles
        REFSITE(:,1) = [0.0_dp, 0.0_dp, -1.0_dp]
        IF(POLYDT) THEN
            OPEN (UNIT = 24, FILE ='poly_angles.dat', STATUS = 'UNKNOWN')
            ALLOCATE(ETPDEL(NSITES,NPART), ETPMDEL(NSITES,NPART))
            IF(READPOLYT) THEN
                DO J1 = 1, NPART
                    READ(24,*) POLYANGLE
                    ETPDEL(1,J1)  = COS(POLYANGLE*PI/180.0_dp)
                    ETPMDEL(1,J1) = PI / ( 1.0_dp - ETPDEL(1,J1) )
                ENDDO
            ELSE
                DO J1 = 1, NPART
                    POLYANGLE = BOX_MULLER ( CPPDELA, POLYSTD )
                    ETPDEL(1,J1)  = COS(POLYANGLE*PI/180.0_dp)
                    ETPMDEL(1,J1) = PI / ( 1.0_dp - ETPDEL(1,J1) )
                    WRITE (24,*) POLYANGLE
                ENDDO
            ENDIF
        ELSE
            CPPDEL(1)    = COS(CPPDELA*PI/180.0_dp)
            CPPMDEL(1)   = PI / ( 1.0_dp - CPPDEL(1) )
        ENDIF
    
    !   The cut-off for the repulsive Kihara component of the potential 
        DCHECKSQ = 2.0_dp**(2.0_dp/GLJN)
    !   Define a cut-off for the shortest distance between the particles
        ETPCUTSQ = RCUTSQ
    !   Adjust the cut-off from the input to account for the fact that the particles
    !   have an elongated core.
        RCUT     = RCUT + RLNGTH
        RCUTSQ   = RCUT * RCUT
    !   The half-length of the ETP particles 
        HLFLNGTH = RLNGTH / 2.0_dp
    
        IF(NSITES == 2) THEN
            CPPIJ(1,2)   = SQRT(CPPAA*CPPBB)
            CPPIJ(2,1)   = CPPIJ(1,2)
            CPPIJ(2,2)   = CPPBB
    
            REFSITE(:,2) = [0.0_dp, 0.0_dp, 1.0_dp]
            IF(POLYDT) THEN
                IF(READPOLYT) THEN
                    DO J1 = 1, NPART
                        READ(24,*) POLYANGLE
                        ETPDEL(2,J1)  = COS(POLYANGLE*PI/180.0_dp)
                        ETPMDEL(2,J1) = PI / ( 1.0_dp - ETPDEL(2,J1) )
                    ENDDO
                ELSE
                    DO J1 = 1, NPART
                        POLYANGLE = BOX_MULLER ( CPPDELB, POLYSTD )
                        ETPDEL(2,J1)  = COS(POLYANGLE*PI/180.0_dp)
                        ETPMDEL(2,J1) = PI / ( 1.0_dp - ETPDEL(2,J1) )
                        WRITE (24,*) POLYANGLE
                    ENDDO
                ENDIF
                CLOSE(UNIT=24)
            ELSE
                CPPDEL(2)    = COS(CPPDELB*PI/180.0_dp)
                CPPMDEL(2)   = PI / ( 1.0_dp - CPPDEL(2) )
            ENDIF
        ENDIF
        
    END SUBROUTINE
    
    !====================================================================================================
    
    SUBROUTINE VIEW_ETP()
    
        USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT, HLFLNGTH
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
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            RM = Q_TO_RM( Q(:,J1) )
            DO J2 = 1, NSITES
                IF(HLFLNGTH /= 0.0_dp) THEN
                    RBCOORDS = RWRITE(:,J1) + HLFLNGTH*MATMUL(RM,REFSITE(:,J2))
                ELSE
                    RBCOORDS = RWRITE(:,J1) + 0.5_dp*MATMUL(RM,REFSITE(:,J2))
                ENDIF
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
    
    SUBROUTINE PESRF_ETP()
        !   =================================================================================================
        !    Subroutine to calculate and output the pair potential for the CPP model, written specifically
        !    for the triblock patchy particles (i.e., particles with two patches).
        !    Outputs both the distance and angular dependence of the potential.
        !   =================================================================================================
        USE COMMONS, ONLY: DP, NPART, NDIM, R, Q, NSITES, RBSITES, REFSITE, RCUT, PI, SPET, CELLLISTT, CPPLAM, ETPCUTSQ
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
    
                print *, J3, DELR, PE, RCUT, ETPCUTSQ
    
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
    
                    IF(DELR == 1.014_dp) THEN
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
    