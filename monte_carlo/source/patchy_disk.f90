SUBROUTINE OBLATE_SPHEROCYLINDER(ENERGY, RIJ, J1, J2)
!   =====================================================================================================================

!   =====================================================================================================================
    USE COMMONS, ONLY: DP, NDIM, NSITES, RBSITES, CLUSTERT, VLMCLUSTERMOVET, OBLL, POLYCHAINT, MULTICHAINET
    USE COMMONS, ONLY: CLSTR,CLSTRSZ,CLSTRID,CLSTRADJ,CLURIJ, OB_BONDST, OB_N_BNDS, OB_NEIGHS, SPET, MCPLYMVT
    USE COMMONS, ONLY: OVERLAPT, PATCHYOBLT, POBLC1, POBLC2, KFLAM2, KFDEL, KFIJ

    ! USE COMMONS, ONLY: NPART, Q, PI, REFSITE, SPET
    ! USE ROTATIONS_MODULE, ONLY: Q_TO_RM

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: J1, J2

    REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM)!, RIJSQ
    REAL(KIND=DP) :: DIJ

    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    INTEGER :: J3, J4, J5, J6
    REAL(KIND=DP) :: CA(NDIM), EA(NDIM), CB(NDIM), EB(NDIM), RAB(NDIM), DAB2, DAB, RABH(NDIM)
    REAL(KIND=DP) :: EARIJ, EBRJI
    ! REAL(KIND=DP) :: RSITE(NDIM), RM(3,3), RMZ(3,3)

    ENERGY = 0.0_dp
    CALL SHORTD_OBLATE(J1, J2, RIJ, RBSITES(:,1,J1), RBSITES(:,1,J2), DIJ)
    
    IF(DIJ < OBLL) THEN 
        OVERLAPT = .TRUE.
        RETURN
    ENDIF

    IF(PATCHYOBLT) THEN
        CA = POBLC1*RBSITES(:,1,J1) ! Orientation of particle I (vector normal to face of the disk)
        CB = POBLC1*RBSITES(:,1,J2) ! Orientation of particle J
        DO J3 = 2, NSITES
            J5 = J3 - 1
        !   Direction of site A on particle I, (the two reference vectors are multiplied by constants to ensure 
        !   that the resulting vector is a unit vector, saves having to normalise the vector every time).
        !   Additionally, the constants bestow a pre-chosen angle to the patch which sets the direction of the patch
        !   relative to the vector defining the direction of the site relative to the centre of the particle.
            EA = CA + POBLC2*RBSITES(:,J3,J1)
            DO J4 = 2, NSITES
                J6 = J4 - 1
            !   Centre-to-centre separation between the patchy sites
                RAB  = RIJ + RBSITES(:,J3,J1) - RBSITES(:,J4,J2)
                DAB2 = DOT_PRODUCT(RAB,RAB)
            !   Check that the patches are close enough to interact
                IF(DAB2 < KFLAM2(J5,J6)) THEN
                    DAB  = SQRT(DAB2)
                !   Make sure there is no overlap of the disks at the patchy sites
                !   (Just in case the shortest-distance algorithm fails since it hasn't been
                !   rigourously shown to always work and we expect the patchy disks to be closest 
                !   at these patchy sites anyway)
                    IF(DAB < OBLL) THEN
                        OVERLAPT = .TRUE.
                        RETURN
                    ENDIF
                    RABH = RAB / DAB

                !   Check if the centre-to-centre separation vector passes through patch A
                    EARIJ = -DOT_PRODUCT(EA,RABH) 
                    IF(EARIJ > KFDEL(J5)) THEN
                    !   Direction of site B on particle J
                        EB = CB + POBLC2*RBSITES(:,J4,J2)
                    !   Check if the centre-to-centre separation vector passes through patch B
                        EBRJI = DOT_PRODUCT(EB,RABH) 
                        IF(EBRJI > KFDEL(J6)) THEN
                        !   Update the pair energy according to the patch-patch interaction matrix
                            ENERGY = ENERGY - KFIJ(J5,J6)
                        !   If performing MC with cluster-moves, add particle J to the current 
                        !   cluster (if it is not already in the cluster). 
                            IF(CLUSTERT .OR. VLMCLUSTERMOVET) THEN
                                IF(POLYCHAINT .AND. MCPLYMVT) THEN
                                    IF(.NOT. MULTICHAINET) RETURN
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
                        !   For following the number of bonds associated with each particle
                        !   useful for Umbrella Sampling simulations.
                            IF(OB_BONDST) THEN                            
                                OB_N_BNDS(J1) = OB_N_BNDS(J1) + 1
                                OB_NEIGHS(OB_N_BNDS(J1),J1) = J2
                                IF(.NOT. SPET) THEN
                                    OB_N_BNDS(J2) = OB_N_BNDS(J2) + 1
                                    OB_NEIGHS(OB_N_BNDS(J2),J2) = J1
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF

                ENDIF
            ENDDO
        ENDDO
    ENDIF

END SUBROUTINE OBLATE_SPHEROCYLINDER

SUBROUTINE SHORTD_OBLATE(I, J, RIJ, EI, EJ, DIJ)
    !====================================================================================================
    !   Subroutine to evaluate the shortest distance between 2.0_dp spherocylinders
    !   RIJ = Vector connecting the geometrical centers of the 2.0_dp spherocylinders
    !   UI  = Unitary vector definig the orientation of particle I
    !   UJ  = Unitary vector definig the orientation of particle J
    !   DIJ = Vector giving the shortest distance between the 2.0_dp spherocylinders   
    !====================================================================================================
    USE COMMONS, ONLY: DP, NDIM, DSIG, HLF_DSIG, RBSITES, R, PI

    USE ROTATIONS_MODULE, ONLY: Q_TO_RM

    IMPLICIT NONE

    INTEGER, INTENT(IN)         :: I, J

    INTEGER                     :: J1, K, L!, COUNT
    REAL(KIND=DP)               :: EIJ, E1(NDIM), E2(NDIM), R12(NDIM), RPARL(NDIM), RPARL2, RPERP(NDIM), RPERPD
    REAL(KIND=DP)               :: U(NDIM), B(NDIM), C0, T, C1(NDIM), QP, QM, QN, XP, XM, X, C2, RP
    REAL(KIND=DP)               :: DDIFF, RT(NDIM), D1, D2, T1(NDIM), T2(NDIM), RIJ(3), EI(3), EJ(3), V1(NDIM), V2(NDIM), V3(NDIM)

    REAL(KIND=DP), INTENT(OUT)  :: DIJ

    EIJ = DOT_PRODUCT(EI,EJ)

    IF(ABS(EIJ-1.0_dp)<1.0e-7) THEN
        RPARL  = DOT_PRODUCT(RIJ,EI)*EI
        RPARL2 = DOT_PRODUCT(RPARL,RPARL)
        RPERP  = RIJ - RPARL
        RPERPD = NORM2(RPERP)

        IF(RPERPD<=DSIG) THEN
            DIJ = SQRT(RPARL2)
        ELSE
            DIJ = SQRT(RPARL2+(RPERPD-DSIG)**2)
        ENDIF
        RETURN
    ENDIF

    DO J1 = 1,2

        IF(J1==1) THEN
            E1  = EI
            E2  = EJ
            R12 = RIJ
            K   = I 
            L   = J
        ELSE
            E1  = EJ
            E2  = EI
            R12 = -RIJ
            K   = J
            L   = I
        ENDIF

        U = [E1(2)*E2(3)-E1(3)*E2(2), E1(3)*E2(1)-E1(1)*E2(3), E1(1)*E2(2)-E1(2)*E2(1)]
        U = U / NORM2(U) 
        B = [E2(2)*U(3)-E2(3)*U(2),   E2(3)*U(1)-E2(1)*U(3),   E2(1)*U(2)-E2(2)*U(1)]
        B = B / NORM2(B)

        C0 = DOT_PRODUCT(R12,E1)
        C2 = HLF_DSIG*DOT_PRODUCT(B,E1)
        T  = - C0 / C2

        C1 = HLF_DSIG * SQRT(1.0_dp - T**2) * U 
        V1 = HLF_DSIG*B
        V2 = V1*T
        V3 = R12 + V2
        
        QP = NORM2(V3 + C1)
        QM = NORM2(V3 - C1)
        QN = MIN(QP,QM)
        
        IF(ABS(T)<1.0_dp) THEN
            IF(QN<=HLF_DSIG) THEN
                DIJ = 0.0_dp
                RETURN
            ELSE
                ! NO SOLUTION
            ENDIF
        ELSE
            XP = (C0 + C2)
            XM = (C0 - C2)
            IF(ABS(XP)<ABS(XM)) THEN
                X  = XP
                RP = NORM2(R12 + V1 - XP*E1)
            ELSE
                X  = XM
                RP = NORM2(R12 - V1 - XM*E1)
            ENDIF
            
            IF(RP < HLF_DSIG) THEN
                DIJ = ABS(X)
                RETURN
            ELSE
                ! NO SOLUTION
            ENDIF
        ENDIF
    ENDDO
    
!   If no solution found yet it means that the minimum distance between 
!   the two disks is located between two points on the edges of the disks.
    DDIFF = 1.0_dp
    T1    = R(:,K)+RBSITES(:,2,K)
    
    ! COUNT = 0
    DO WHILE (DDIFF > 1.e-07)
        ! COUNT = COUNT + 1
!   Compute the shortest distance between I and J from point 1 on disk I
        RPARL  = DOT_PRODUCT(T1,E1)*E1
        RPERP  = T1 - RPARL
        RPERPD = NORM2(RPERP)

        IF(RPERPD <= HLF_DSIG) THEN
            D1 = NORM2(RPARL)
            RT = RPARL
        ELSE
            RT = (T1 - HLF_DSIG*RPERP/RPERPD)
            D1 = NORM2(RT)
        ENDIF
        
        T2 = T1 - RT - R12

!   Compute the shortest distance between I and J from point 1 on disk J
        RPARL  = DOT_PRODUCT(T2,E2)*E2
        RPERP  = T2 - RPARL
        RPERPD = NORM2(RPERP)

        IF(RPERPD <= HLF_DSIG) THEN
            D2 = NORM2(RPARL)
            RT = RPARL
        ELSE
            RT = (T2 - HLF_DSIG*RPERP/RPERPD)
            D2 = NORM2(RT)
        ENDIF

        T1 = T2 - RT + R12

        DDIFF = ABS(D1-D2)
    ENDDO

    DIJ = D2

END SUBROUTINE

! ==================================================================================================================
! ==================================================================================================================
! ==================================================================================================================
! ==================================================================================================================

SUBROUTINE INTRA_POLYMER(ENERGY, J1)
!   =====================================================================================================================
!   Compute the harmonic energy associate with nearest neighbours in freely-jointed polymer chains
!   =====================================================================================================================
    USE COMMONS, ONLY: DP, NDIM, BOX, R, POLYC_L, POLY_SIG, POLY_KAP, N_POLY_L, SPET, OBJCTT, OVERLAPT, APCHAINT, PCHAIN_EPS,POLYCUT
    USE COMMONS, ONLY: N_POLY
    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1
    INTEGER                    :: J2, BEAD_ID, BI, BJ, CHNI, BDS, BDE
    REAL(KIND=DP)              :: RIJ(NDIM), DIST
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp

    BEAD_ID = MOD(J1,N_POLY_L)
    
    IF(BEAD_ID/=1 .AND. SPET ) THEN
        RIJ    = R(:,J1) - R(:,J1-1)
        IF(OBJCTT) THEN
            RIJ(1:2) = RIJ(1:2) - BOX(1:2)*ANINT(RIJ(1:2)/BOX(1:2))
        ELSE
            RIJ = RIJ - BOX*ANINT(RIJ/BOX)
        ENDIF
        DIST   = NORM2(RIJ)
        IF(DIST<POLY_SIG) THEN
            OVERLAPT = .TRUE.
            RETURN
        ENDIF
        ENERGY = ENERGY + POLY_KAP*(DIST-POLYC_L*POLY_SIG)**2

        IF(APCHAINT) THEN
            BI = MOD(J1-1,N_POLY_L)
            IF(PCHAIN_EPS(BI,BI-1) /= 0.0_dp) THEN
                IF(DIST<POLYCUT) ENERGY = ENERGY - PCHAIN_EPS(BI,BI-1)
            ENDIF
        ENDIF

    ENDIF

    IF(BEAD_ID/=0) THEN
        RIJ    = R(:,J1) - R(:,J1+1)
        IF(OBJCTT) THEN
            RIJ(1:2) = RIJ(1:2) - BOX(1:2)*ANINT(RIJ(1:2)/BOX(1:2))
        ELSE
            RIJ = RIJ - BOX*ANINT(RIJ/BOX)
        ENDIF
        DIST   = NORM2(RIJ)
        IF(DIST<POLY_SIG) THEN
            OVERLAPT = .TRUE.
            RETURN
        ENDIF
        ENERGY = ENERGY + POLY_KAP*(DIST-POLYC_L*POLY_SIG)**2

        IF(APCHAINT) THEN
            BI = MOD(J1-1,N_POLY_L)
            IF(PCHAIN_EPS(BI,BI+1) /= 0.0_dp) THEN
                IF(DIST<POLYCUT) ENERGY = ENERGY - PCHAIN_EPS(BI,BI+1)
            ENDIF
        ENDIF

    ENDIF

    ! IF(APCHAINT) THEN
    !     BI = MOD(J1-1,N_POLY_L)
    !     CHNI = (J1-1)/N_POLY_L
    !     BDS  = CHNI*N_POLY_L + 1 
    !     BDE  = BDS + N_POLY_L - 1
    !     DO J2 = BDS, BDE
    !         IF(J2==J1) CYCLE
    !         BJ = MOD(J2-1,N_POLY_L)
    !         IF(PCHAIN_EPS(BI,BJ) /= 0.0_dp) THEN
    !             RIJ    = R(:,J1) - R(:,J2)
    !             IF(OBJCTT) THEN
    !                 RIJ(1:2) = RIJ(1:2) - BOX(1:2)*ANINT(RIJ(1:2)/BOX(1:2))
    !             ELSE
    !                 RIJ = RIJ - BOX*ANINT(RIJ/BOX)
    !             ENDIF
    !             DIST   = NORM2(RIJ)
    !             IF(DIST<POLYCUT) ENERGY = ENERGY - PCHAIN_EPS(BI,BJ)
    !         ENDIF
    !     ENDDO
    ! ENDIF

END SUBROUTINE INTRA_POLYMER

SUBROUTINE INTER_POLYMER(ENERGY, J1, J2, RIJSQ)
!   =====================================================================================================================
!   Compute the hard repulsive interaction associated with non-nearest neighbours in freely-jointed polymer chains
!   =====================================================================================================================
    USE COMMONS, ONLY: DP, POLY_SIG2, OVERLAPT, ATTR_POLY, POLYCUT2, POLY_EPS, APCHAINT, PCHAIN_EPS, RCHTT, SPET
    USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, OVERLAPT, POLY_CLU, N_POLY_L, RCHTT, MULTICHAINET, MCPLYMVT
    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1, J2
    REAL(KIND=DP), INTENT(IN)  :: RIJSQ

    INTEGER                    :: J4, POLY_ID, CHAINID, B1, B2

    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp

    IF(RIJSQ<POLY_SIG2) THEN
        OVERLAPT = .TRUE.
        RETURN
    ENDIF

    IF(ATTR_POLY .AND. RIJSQ<POLYCUT2) THEN
        ENERGY = -POLY_EPS
        IF(CLUSTERT .AND. MCPLYMVT .AND. (.NOT. RCHTT) .AND. MULTICHAINET) THEN
            POLY_ID = FLOOR( REAL((J2-1)) / REAL(N_POLY_L) )
            CHAINID = POLY_ID*N_POLY_L + 1
            IF(.NOT.(ANY(CLSTR(1:CLSTRSZ)==J2)) .AND. POLY_ID /= POLY_CLU) THEN
                DO J4 = CHAINID, CHAINID+N_POLY_L-1
                    CLSTRSZ = CLSTRSZ + 1
                    CLSTR(CLSTRSZ) = J4
                ENDDO
            ENDIF
        ENDIF
    ENDIF

    IF(APCHAINT .AND. ((J1-1)/N_POLY_L==(J2-1)/N_POLY_L)) THEN
        B1 = MOD(J1-1,N_POLY_L)
        B2 = MOD(J2-1,N_POLY_L)
        IF(PCHAIN_EPS(B1,B2) /= 0.0_dp) THEN
            IF(RIJSQ<POLYCUT2) ENERGY = ENERGY - PCHAIN_EPS(B1,B2)
        ENDIF
    ENDIF
    
END SUBROUTINE INTER_POLYMER

! ==================================================================================================================
! ==================================================================================================================
! ==================================================================================================================
! ==================================================================================================================

SUBROUTINE DISK_POLYMER(ENERGY, RIJ, J1, J2)
!   =====================================================================================================================
!   Polymer-capsomer interaction, particle I is always a polymer bead and particle J is always a capsomer
!   =====================================================================================================================
    USE COMMONS, ONLY: DP, PI, HLFPI, NDIM, R, RBSITES, CLUSTERT, CP_EPS, HLF_DSIG, CP_SIG1, CP_DEL, PATCH_RAD
    USE COMMONS, ONLY: CLSTR, CLSTRSZ, OVERLAPT, POLY_CLU, N_POLY_L, RCHTT, MULTICHAINET, MCPLYMVT

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1, J2
    INTEGER                    :: J3, J4, POLY_ID, CHAINID
    REAL(KIND=DP), INTENT(IN)  :: RIJ(NDIM)
    REAL(KIND=DP)              :: RJ(NDIM), RAB(NDIM), NI(NDIM), RN, PJ(NDIM), DJ, RP(NDIM), DRP, THETA
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp

    IF(J1 > J2) THEN
        J3  = J1             ! Index of capsomer particle
        RAB = -RIJ
    ELSE
        J3  = J2
        RAB = RIJ
    ENDIF

    RJ = R(:,J3) + RAB ! Position of polymer bead using MIC w.r.t. capsomer
    NI = RBSITES(:,1,J3)
    RN = DOT_PRODUCT(RAB,NI) 

    PJ = RAB - RN*NI
    DJ = NORM2(PJ)

    IF(DJ <= HLF_DSIG) THEN ! Check if the bead lies over the face of the disk
    !   Compute the distance between the bead and the face of the disk
        RP  = RAB - PJ
        DRP = NORM2(RP)
    !   Check if the bead and the disk are overlapping one another    
        IF(DRP<CP_SIG1) THEN
            OVERLAPT = .TRUE.
            RETURN
    !   Check if the bead is close enough to the face of the disk to interact with it, also
    !   if the patch is a ring instead of a circe make sure that the bead lies on the area covered by the ring.
        ELSEIF(DRP<=CP_SIG1+CP_DEL .AND. DJ > PATCH_RAD) THEN
            THETA = ACOS(RN/NORM2(RAB))
        !   Check the bead on the correct side of the disk to interact with the patch
            IF( THETA<HLFPI ) THEN 
                ENERGY = -CP_EPS
                IF(CLUSTERT .AND. MCPLYMVT) THEN
            !   Check if the particle is already in the cluster, only need to check 
            !   particle IDs that come after the current particle in the list as the
            !   previous check should take care of those that come earlier.
                !   If particle J is a capsomer and particle I a polymer bead then add the capsomer
                !   to the list of particles in the growing cluster.
                    IF(.NOT.( ANY(CLSTR(1:CLSTRSZ)==J2) ) .AND. J2>J1) THEN
                        CLSTRSZ = CLSTRSZ + 1
                        CLSTR(CLSTRSZ) = J2
                !   If particle J is a polymer bead then add the whole polymer chain to the cluster.
                !   However, the polymer chain is only added to the cluster if the bead belongs to a different
                !   polymer chain than the one initially selected at the start of the MC move.
                    ELSEIF(J1>J2 .AND. MULTICHAINET) THEN
                        POLY_ID = FLOOR( REAL((J2-1)) / REAL(N_POLY_L) )
                        CHAINID = POLY_ID*N_POLY_L + 1
                        IF(.NOT.(ANY(CLSTR(1:CLSTRSZ)==J2)) .AND. POLY_ID /= POLY_CLU) THEN
                            DO J4 = CHAINID, CHAINID+N_POLY_L-1
                                IF(.NOT.( ANY(CLSTR(1:CLSTRSZ)==J4) )) THEN
                                    CLSTRSZ = CLSTRSZ + 1
                                    CLSTR(CLSTRSZ) = J4
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
        ENDIF

    ELSE
    !   The bead doesn't overlap with the face of the disk, so the shortest distance 
    !   is between the edge of the disk and the center of the bead.
        PJ  = PJ/DJ * HLF_DSIG ! Find the point on the edge of the disk that is closest to the bead.
        RP  = RAB - PJ
        DRP = NORM2(RP)
    !   Check if the bead and the disk are overlapping one another
        IF(DRP<CP_SIG1) THEN
            OVERLAPT = .TRUE.
            RETURN
        ENDIF

    ENDIF

END SUBROUTINE DISK_POLYMER

SUBROUTINE DEF_OBLATE_SPHEROCYLINDER()
    
    USE COMMONS, ONLY: DP, REFSITE, HLF_DSIG, PATCHYOBLT, PI, NSITES, POBLTHETA, POBLC1, POBLC2, OBLL
    USE COMMONS, ONLY: KFIJ, KFLAM2, KFDEL, KFAA, KFBB, KFCC, KFDD, KFEE
    USE COMMONS, ONLY: KFLAMA, KFLAMB, KFLAMC, KFLAMD, KFLAME
    USE COMMONS, ONLY: KFDELA, KFDELB, KFDELC, KFDELD, KFDELE
    USE COMMONS, ONLY: POLYCHAINT, POLY_SIG, CP_SIG1, CP_SIG12

    IMPLICIT NONE

    INTEGER       :: J1
    REAL(KIND=DP) :: PENT_ANGL, RMZ(3,3)

    REFSITE(:,1)= [ 0.0_dp,  0.0_dp,  1.0_dp ]
    REFSITE(:,2)= [ HLF_DSIG, 0.0_dp, 0.0_dp]

    IF(PATCHYOBLT) THEN
        PENT_ANGL = 72.0_dp*(PI/180_dp)
        RMZ = 0.0_dp
        RMZ(1,1) = COS(PENT_ANGL)
        RMZ(1,2) = SIN(PENT_ANGL)
        RMZ(2,1) = -SIN(PENT_ANGL)
        RMZ(2,2) = COS(PENT_ANGL)
        RMZ(3,3) = 1.0_dp

        DO J1 = 3, NSITES
            REFSITE(:,J1)=MATMUL(RMZ,REFSITE(:,J1-1))
        ENDDO

        POBLC1 = SIN(POBLTHETA*PI/180_dp)
        POBLC2 = COS(POBLTHETA*PI/180_dp) / HLF_DSIG

        ALLOCATE( KFIJ(5,5), KFLAM2(5,5), KFDEL(5) )

        KFDEL(1)  = COS(KFDELA*PI/180_dp)
        KFDEL(2)  = COS(KFDELB*PI/180_dp)
        KFDEL(3)  = COS(KFDELC*PI/180_dp)
        KFDEL(4)  = COS(KFDELD*PI/180_dp)
        KFDEL(5)  = COS(KFDELE*PI/180_dp)

        KFIJ(1,1) = KFAA
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

        KFLAM2(1,1) = KFLAMA**2
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
    ENDIF

    IF(POLYCHAINT) THEN
        CP_SIG1  = (POLY_SIG+OBLL)/2.0_dp
        CP_SIG12 = CP_SIG1**2
    ENDIF

END SUBROUTINE DEF_OBLATE_SPHEROCYLINDER

SUBROUTINE VIEW_OBLATE_SPHEROCYLINDER()
        
    USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT, POLYCHAINT, N_POLY_TOT
    USE ROTATIONS_MODULE, ONLY: Q_TO_RM

    IMPLICIT NONE

    INTEGER:: J1, J2
    REAL(KIND=DP) :: RWRITE(3,NPART), RM(3,3), RBCOORDS(3)

    IF(POLYCHAINT) THEN
        WRITE(VIEWUNIT,*) N_POLY_TOT+(NPART-N_POLY_TOT)*(NSITES+1)
    ELSE
        WRITE(VIEWUNIT,*) NPART*(NSITES+1)
    ENDIF

    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO J1 = 1, NPART
        RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX(:))*BOX(:)
    END DO

    DO J1 = 1, NPART
        IF(POLYCHAINT) THEN
            IF(J1<= N_POLY_TOT) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'H ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            ELSE
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
                RM = Q_TO_RM( Q(:,J1) )
                DO J2 = 1, NSITES
                    IF(J2==1) THEN
                        RBCOORDS = RWRITE(:,J1) + 0.25_dp*MATMUL(RM,REFSITE(:,J2))
                        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ELSE
                        RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
                        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ENDIF
                ENDDO
            ENDIF
        ELSE
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            RM = Q_TO_RM( Q(:,J1) )
            DO J2 = 1, NSITES
                IF(J2==1) THEN
                    RBCOORDS = RWRITE(:,J1) + 0.25_dp*MATMUL(RM,REFSITE(:,J2))
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSE
                    RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ENDIF
            ENDDO
        ENDIF
    END DO
    
END SUBROUTINE