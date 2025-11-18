SUBROUTINE HARD_TRIANGLE(ENERGY, RIJ, J1, J2)
!   ==========================================================================================
!   ==========================================================================================

!   ==========================================================================================
!   ==========================================================================================
    USE COMMONS, ONLY: HLFPI, PI, DP, NDIM, RBSITES, OVERLAPT, TRI_SIG, TRI_SIG2, OVERLAPT
    USE COMMONS, ONLY: PATCHYTRIT, KFLAM2, TRI_DEL, KFIJ, SQWT, SQWDEL2, SQWEPS
    USE COMMONS, ONLY: CLUSTERT, VLMCLUSTERMOVET,CLSTR,CLSTRSZ,CLSTRID,CLSTRADJ,CLURIJ
    USE COMMONS, ONLY: POLYCHAINT, MULTICHAINET, MCPLYMVT

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: J1, J2
    REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM)

    INTEGER       :: J3, J4, J5, J6
    REAL(KIND=DP) :: DIJ(NDIM), D2, EDGES(2,NDIM,3), AXANGL, NV(NDIM)
    REAL(KIND=DP) :: RAB(NDIM), DAB2, DAB, RABH(NDIM)
    REAL(KIND=DP) :: EARIJ, EBRJI!, RIJ(3)

    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp
    
    CALL GET_SHORTEST_DIST_TRI_TRI(J1, J2, RIJ, DIJ, D2, EDGES)
    IF(D2<TRI_SIG2) THEN
        OVERLAPT = .TRUE.
        RETURN
    ELSEIF(SQWT .AND. D2<SQWDEL2) THEN
        ENERGY = ENERGY - SQWEPS
    ENDIF

    IF(PATCHYTRIT) THEN
        DO J3 = 5, 7
            J5 = J3 - 4
            DO J4 = 5, 7
                J6 = J4 - 4
            
            !   Centre-to-centre separation between the patchy sites
                RAB  = RIJ + RBSITES(:,J3,J1) - RBSITES(:,J4,J2)
                DAB2 = DOT_PRODUCT(RAB,RAB)
                DAB  = SQRT(DAB2)
                RABH = RAB / DAB
            
            !   Check angle between the separation vector and the edges.
                EARIJ = -DOT_PRODUCT(EDGES(1,:,J5),RABH) 
                IF(EARIJ <= -TRI_DEL(2,J5) .OR. EARIJ >= TRI_DEL(2,J5)) CYCLE
                EBRJI = DOT_PRODUCT(EDGES(2,:,J6),RABH) 
                IF(EBRJI <= -TRI_DEL(2,J6) .OR. EBRJI >= TRI_DEL(2,J6)) CYCLE
            
            !   Check if the patches close enough to interact. Must account for the fact that 
            !   the particles have cylindrical edges and so the cutoff distance will depend on
            !   the angle between the separation vector and the edge.
                IF(DAB2 > KFLAM2(J5,J6)/(1.0_dp-EARIJ**2)) CYCLE
                IF(DAB2 > KFLAM2(J5,J6)/(1.0_dp-EBRJI**2)) CYCLE

            !   Project the separation vector onto the plane containing the vector defining the 
            !   direction of the patch on each particle, and check if the patches are properly 
            !   oriented to be overlapping.
                NV = -RAB - EARIJ*DAB*EDGES(1,:,J5) ! Particle I
                IF(DOT_PRODUCT(RBSITES(:,J3+3,J1),NV/NORM2(NV)) <= TRI_DEL(1,J5)) CYCLE
                NV = RAB - EBRJI*DAB*EDGES(2,:,J6)  ! Particle J
                IF(DOT_PRODUCT(RBSITES(:,J4+3,J2),NV/NORM2(NV)) > TRI_DEL(1,J6)) THEN

                ! Update the pair energy according to the patch-patch interaction matrix
                    AXANGL = DOT_PRODUCT(EDGES(1,:,J5),EDGES(2,:,J6)) - TRI_SIG/DAB !1.0_dp
                    ENERGY = ENERGY + KFIJ(J5,J6)*AXANGL
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

                ENDIF
            ENDDO
        ENDDO
    ENDIF

    ! PRINT *, "ENERGY=",ENERGY
    ! STOP


END SUBROUTINE


SUBROUTINE HARD_TRIANGLE_SPHERE(ENERGY, RIJ, J1, J2)
!   ==========================================================================================
!   ==========================================================================================

!   ==========================================================================================
!   ==========================================================================================
    USE COMMONS, ONLY: HLFPI, PI, DP, NDIM, RBSITES, OVERLAPT, R, OVERLAPT
    USE COMMONS, ONLY: CP_SIG12, TP_CUT2, N_POLY_L, CP_EPS
    USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, MCPLYMVT, MULTICHAINET, POLY_CLU

    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J1, J2
    REAL(KIND=DP), INTENT(IN)  :: RIJ(NDIM)
    
    INTEGER                    :: J3, J4, POLY_ID, CHAINID
    REAL(KIND=DP)              :: V1(NDIM), V2(NDIM), V3(NDIM), DIJ(NDIM), DIST2
    REAL(KIND=DP)              :: RJ(NDIM), RAB(NDIM), RN, THETA
    LOGICAL                    :: FACET
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp

    IF(J1 > J2) THEN
        J3  = J1 ! Index of capsomer particle
        RAB = -RIJ
    ELSE
        J3  = J2
        RAB = RIJ
    ENDIF

    RJ = R(:,J3) + RAB              ! Position of polymer bead using MIC w.r.t. capsomer
    V1 = R(:,J3) + RBSITES(:,2,J3)  ! Position of vertex 1 of the triangle
    V2 = R(:,J3) + RBSITES(:,3,J3)  ! Position of vertex 2 of the triangle
    V3 = R(:,J3) + RBSITES(:,4,J3)  ! Position of vertex 3 of the triangle

    CALL GET_SHORTEST_DIST_TRI_POINT(RJ, V1, V2, V3, DIJ, FACET)

    DIST2 = DOT_PRODUCT(DIJ,DIJ)
    IF(DIST2 < CP_SIG12) THEN
        OVERLAPT = .TRUE.
        RETURN
    ENDIF

    IF(FACET .AND. DIST2<=TP_CUT2) THEN
        ! Check the bead on the correct side of the disk to interact with the patch
        RN = DOT_PRODUCT(RAB,RBSITES(:,1,J3)) 
        THETA = ACOS(RN/NORM2(RAB))
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
END SUBROUTINE

SUBROUTINE GET_SHORTEST_DIST_TRI_TRI(J1, J2, RIJ, MINDIJ, MIND2, EDGES)
!   ==========================================================================================
!   ==========================================================================================

!   ==========================================================================================
!   ==========================================================================================
    USE COMMONS, ONLY: NDIM, DP, PI, R, RBSITES, TRI_EL
    ! USE COMMONS, ONLY: TRI_A, TRI_B, TRI_C
    USE ROTATIONS_MODULE, ONLY: CROSS

    IMPLICIT NONE

    INTEGER, INTENT(IN)       :: J1,J2
    REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM)

    INTEGER             :: J3, J4, J5, POINT
    REAL(KIND=DP)       :: RI(NDIM), RJ(NDIM), NI(NDIM), NJ(NDIM), UI(NDIM), UJ(NDIM)
    REAL(KIND=DP)       :: VI(2,NDIM,3), VJ(2,NDIM,3), VIJ(NDIM)
    REAL(KIND=DP)       :: Z(NDIM), Y(NDIM), X(NDIM), ANG1, ANG2
    REAL(KIND=DP)       :: S(NDIM), TP(NDIM)
    REAL(KIND=DP)       :: D2, DIJ(NDIM)
    LOGICAL             :: DISJOINT

    REAL(KIND=DP), INTENT(OUT) :: MINDIJ(NDIM), MIND2, EDGES(2,NDIM,3)

    DISJOINT = .FALSE.

    RI = R(:,J1)
!   Get vectors connecting each of the vertices on triangle I 
    VI(1,:,1) = (RBSITES(:,3,J1) - RBSITES(:,2,J1)) ! Edge A-B
    VI(1,:,2) = (RBSITES(:,4,J1) - RBSITES(:,3,J1)) ! Edge B-C
    VI(1,:,3) = (RBSITES(:,2,J1) - RBSITES(:,4,J1)) ! Edge C-A
!   Get the centre of each of the vertices on triangle I 
    VI(2,:,1) = RI+(RBSITES(:,3,J1) + RBSITES(:,2,J1))/2.0 ! Edge A-B
    VI(2,:,2) = RI+(RBSITES(:,4,J1) + RBSITES(:,3,J1))/2.0 ! Edge B-C
    VI(2,:,3) = RI+(RBSITES(:,2,J1) + RBSITES(:,4,J1))/2.0 ! Edge C-A

    RJ = R(:,J1) - RIJ
!   Get vectors connecting each of the vertices on triangle J
    VJ(1,:,1) = (RBSITES(:,3,J2) - RBSITES(:,2,J2)) ! Edge A-B
    VJ(1,:,2) = (RBSITES(:,4,J2) - RBSITES(:,3,J2)) ! Edge B-C
    VJ(1,:,3) = (RBSITES(:,2,J2) - RBSITES(:,4,J2)) ! Edge C-A
!   Get the centre of each of the vertices on triangle I 
    VJ(2,:,1) = RJ+(RBSITES(:,3,J2) + RBSITES(:,2,J2))/2.0 ! Edge A-B
    VJ(2,:,2) = RJ+(RBSITES(:,4,J2) + RBSITES(:,3,J2))/2.0 ! Edge B-C
    VJ(2,:,3) = RJ+(RBSITES(:,2,J2) + RBSITES(:,4,J2))/2.0 ! Edge C-A

    DO J3 = 1, 3
        EDGES(1,:,J3) = VI(1,:,J3)/NORM2(VI(1,:,J3))
        EDGES(2,:,J3) = VJ(1,:,J3)/NORM2(VJ(1,:,J3))
    ENDDO

!   Set the initial value for the minimum squared distance to be arbitrarily high.
    MIND2 = 100.0_dp
!   Loop through each pair of edges on the two triangles and find the minimum distance
!   between each pair, keeping track of which pair were closest together.
    DO J3 = 1, 3
        UI = EDGES(1,:,J3)
        DO J4 = 1, 3
            UJ = EDGES(2,:,J4)
            VIJ = VI(2,:,J3) - VJ(2,:,J4)
            CALL SHORTD_TRI_EDGE(VIJ, UI, UJ, DIJ, NI, NJ, TRI_EL(J3), TRI_EL(J4))
            NI = VI(2,:,J3) + NI ! Position along edge I closest to edge J
            NJ = VJ(2,:,J4) - NJ ! Position along edge J closest to edge I
            D2 = DOT_PRODUCT(DIJ,DIJ) ! Square of the shortest distance between edges I and J
        !   If the current pair of edges are the closest found so far then check to see if 
        !   they satisfy the condition for being the closest points between the two triangles.
            IF(D2 <= MIND2) THEN
                MIND2  = D2
                MINDIJ = DIJ
            !   Get vector connecting point NI to the 3rd vertex on triangle I
                J5 = MOD((J3+3),3) + 3
                IF(J5 > 4) J5 = 2
                Z    = RI + RBSITES(:,J5,J1) - NI
            !   Get the angle between vector Z and the shortest distance vector. 
                ANG1 = DOT_PRODUCT(Z,DIJ) 
            !   Get vector connecting point NJ to the 3rd vertex on triangle J
                J5 = MOD((J4+3),3) + 3
                IF(J5 > 4) J5 = 2
                Z    = RJ + RBSITES(:,J5,J2) - NJ
            !   Get the angle between vector Z and the shortest distance vector.
                ANG2 = DOT_PRODUCT(Z,DIJ) 
            !   If the two vertices lie between the planes perpendicular to the shortest
            !   distance vector then points NI and NJ correspond to closest points on the 
            !   two triangles. 
                IF((ANG1>=0.0_dp) .AND. (ANG2<=0.0_dp)) RETURN

                IF(ANG1>0.0_dp) ANG1 = 0.0_dp
                IF(ANG2<0.0_dp) ANG2 = 0.0_dp
                IF((D2+ANG1+ANG2)>0.0_dp) DISJOINT = .TRUE.

            ENDIF
        ENDDO
    ENDDO

!   No edge pairs contained the closest points, therefore, either:
!   1. one of the closest points is a vertex, and the other point is interior to a face.
!   2. the triangles are overlapping.
!   3. an edge of one triangle is parallel to the other's face. If cases 1 and 2 are not true, 
!      then the closest points from the 9 edge pairs checks above can be taken as closest
!      points for the triangles.
!   4. possibly, the triangles were degenerate.  When the triangle points are nearly colinear
!      or coincident, one  of above tests might fail even though the edges tested contain the closest points.

!   -------------------------------------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------------------------------------
!   Check for case 1
!   -------------------------------------------------------------------------------------------------------------
!   Check projections of vertices of triangle J onto triangle I
!   -------------------------------------------------------------------------------------------------------------
    S = CROSS(VI(1,:,1),VI(1,:,2))!RBSITES(:,1,J2)
    
    TP(1) = DOT_PRODUCT(S, ((RI+RBSITES(:,2,J1)) - (RJ+RBSITES(:,2,J2))))
    TP(2) = DOT_PRODUCT(S, ((RI+RBSITES(:,2,J1)) - (RJ+RBSITES(:,3,J2))))
    TP(3) = DOT_PRODUCT(S, ((RI+RBSITES(:,2,J1)) - (RJ+RBSITES(:,4,J2))))

    POINT = -1

    IF( (TP(1)>0) .AND. (TP(2)>0) .AND. (TP(3)>0) ) THEN
        IF(TP(1)<TP(2)) THEN
            POINT = 2
        ELSE
            POINT = 3
        ENDIF
        IF(TP(3)<TP(POINT-1)) POINT = 4
    ELSEIF((TP(1)<0) .AND. (TP(2)<0) .AND. (TP(3)<0)) THEN
        IF(TP(1)>TP(2)) THEN
            POINT = 2
        ELSE
            POINT = 3
        ENDIF
        IF(TP(3)>TP(POINT-1)) POINT = 4
    ENDIF

    IF(POINT > 0) THEN
        DISJOINT = .TRUE.
        Y = (RJ + RBSITES(:,POINT,J2))
        
        X = Y - (RI + RBSITES(:,2,J1))
        Z = CROSS(S,VI(1,:,1))
        
        IF(DOT_PRODUCT(X,Z)>0.0_dp) THEN
            X = Y - (RI + RBSITES(:,3,J1))
            Z = CROSS(S,VI(1,:,2))
            
            IF(DOT_PRODUCT(X,Z)>0.0_dp) THEN
                X = Y - (RI + RBSITES(:,4,J1))
                Z = CROSS(S,VI(1,:,3))
                
                IF(DOT_PRODUCT(X,Z)>0.0_dp) THEN
                    NI = Y+(S*(TP(POINT-1)/DOT_PRODUCT(S,S)))
                    NJ = Y
                    MINDIJ = NI - NJ
                    MIND2 = DOT_PRODUCT(MINDIJ,MINDIJ)
                    RETURN
                ENDIF
            ENDIF
        ENDIF
    ENDIF
!   -------------------------------------------------------------------------------------------------------------
!   Check projections of vertices of triangle I onto triangle J
!   -------------------------------------------------------------------------------------------------------------
    S = CROSS(VJ(1,:,1),VJ(1,:,2))!RBSITES(:,1,J2)
    
    ! PRINT *, S
    ! PRINT *, RBSITES(:,1,J2)
    TP(1) = DOT_PRODUCT(S, ((RJ+RBSITES(:,2,J2)) - (RI+RBSITES(:,2,J1))))
    TP(2) = DOT_PRODUCT(S, ((RJ+RBSITES(:,2,J2)) - (RI+RBSITES(:,3,J1))))
    TP(3) = DOT_PRODUCT(S, ((RJ+RBSITES(:,2,J2)) - (RI+RBSITES(:,4,J1))))

    POINT = -1

    IF( (TP(1)>0) .AND. (TP(2)>0) .AND. (TP(3)>0) ) THEN
        IF(TP(1)<TP(2)) THEN
            POINT = 2
        ELSE
            POINT = 3
        ENDIF
        IF(TP(3)<TP(POINT-1)) POINT = 4
    ELSEIF((TP(1)<0) .AND. (TP(2)<0) .AND. (TP(3)<0)) THEN
        IF(TP(1)>TP(2)) THEN
            POINT = 2
        ELSE
            POINT = 3
        ENDIF
        IF(TP(3)>TP(POINT-1)) POINT = 4
    ENDIF

    IF(POINT > 0) THEN
        DISJOINT = .TRUE.
        Y = (RI + RBSITES(:,POINT,J1))
        
        X = Y - (RJ + RBSITES(:,2,J2))
        Z = CROSS(S,VJ(1,:,1))
        
        IF(DOT_PRODUCT(X,Z)>0.0_dp) THEN
            X = Y - (RJ + RBSITES(:,3,J2))
            Z = CROSS(S,VJ(1,:,2))
            
            IF(DOT_PRODUCT(X,Z)>0.0_dp) THEN
                X = Y - (RJ + RBSITES(:,4,J2))
                Z = CROSS(S,VJ(1,:,3))
                
                IF(DOT_PRODUCT(X,Z)>0.0_dp) THEN
                    NJ = Y+(S*(TP(POINT-1)/DOT_PRODUCT(S,S)))
                    NI = Y
                    MINDIJ = NI - NJ
                    MIND2 = DOT_PRODUCT(MINDIJ,MINDIJ)
                    RETURN
                ENDIF
            ENDIF
        ENDIF
    ENDIF
!   -------------------------------------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------------------------------------

    IF(.NOT. DISJOINT) THEN
        DIJ = 0.0_dp
        D2  = 0.0_dp
        RETURN
    ENDIF

!   -------------------------------------------------------------------------------------------------------------
!   -------------------------------------------------------------------------------------------------------------

END SUBROUTINE

SUBROUTINE GET_SHORTEST_DIST_TRI_POINT(P, V1, V2, V3, DIJ, FACET)
!   ==========================================================================================
!   ==========================================================================================
!   https://github.com/embree/embree/blob/master/tutorials/common/math/closest_point.h
!   ==========================================================================================
!   ==========================================================================================
    USE COMMONS, ONLY: DP

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN)   :: P(3), V1(3), V2(3), V3(3)

    REAL(KIND=DP)               :: AB(3), AC(3), AP(3), BP(3), CP(3), DENOM, W
    REAL(KIND=DP)               :: D1, D2, D3, D4, D5, D6, VC, VB, VA, V, C1, C2

    REAL(KIND=DP)               :: NT(3)

    REAL(KIND=DP), INTENT(OUT)  :: DIJ(3)
    LOGICAL, INTENT(OUT)        :: FACET

    FACET = .FALSE.

    AB = V2 - V1
    AC = V3 - V1

!   ----------------------------------------------------------------------
!   Check if the closest point is vertex 1 of the triangle
!   ----------------------------------------------------------------------
    AP = P - V1
    D1 = DOT_PRODUCT(AB,AP)
    D2 = DOT_PRODUCT(AC,AP)
    IF(D1 <= 0.0_dp .AND. D2 <= 0.0_dp) THEN
        NT = V1 ! Closest point is V1
        GO TO 200
    ENDIF

!   ----------------------------------------------------------------------
!   Check if the closest point is vertex 2 of the triangle
!   ----------------------------------------------------------------------
    BP = P - V2
    D3 = DOT_PRODUCT(AB,BP)
    D4 = DOT_PRODUCT(AC,BP)
    IF(D3 >= 0.0_dp .AND. D4 <= D3) THEN
        NT = V2 ! Closest point is V2
        GO TO 200
    ENDIF

!   ----------------------------------------------------------------------
!   Check if the closest point is vertex 2 of the triangle
!   ----------------------------------------------------------------------
    CP = P - V3
    D5 = DOT_PRODUCT(AB,CP)
    D6 = DOT_PRODUCT(AC,CP)
    IF(D6 >= 0.0_dp .AND. D5 <= D6) THEN
        NT = V3 ! Closest point is V3
        GO TO 200
    ENDIF

!   ----------------------------------------------------------------------
!   Check if the closest point is the V1-V2 edge of the triangle
!   ----------------------------------------------------------------------
    VC = D1*D4 - D3*D2
    IF(VC <= 0.0_dp .AND. D1 >= 0.0_dp .AND. D3 <= 0.0_dp) THEN
        V  = D1 / (D1 - D3)
        NT = V1 + V*AB
        GO TO 200
    ENDIF

!   ----------------------------------------------------------------------
!   Check if the closest point is the V1-V3 edge of the triangle
!   ----------------------------------------------------------------------
    VB = D5*D2 - D1*D6
    IF(VB <= 0.0_dp .AND. D2 >= 0.0_dp .AND. D6 <= 0.0_dp) THEN
        V  = D2 / (D2 - D6)
        NT = V1 + V*AC
        GO TO 200
    ENDIF

!   ----------------------------------------------------------------------
!   Check if the closest point is the V2-V3 edge of the triangle
!   ----------------------------------------------------------------------
    VA = D3*D6 - D5*D4
    C1 = D4-D3
    C2 = D5-D6
    IF(VA <= 0.0_dp .AND. C1 >= 0.0_dp .AND. C2 >= 0.0_dp) THEN
        V  = C1 / (C1 + C2)
        NT = V2 + V*(V3-V2)
        GO TO 200
    ENDIF

!   ----------------------------------------------------------------------
!   Closest point is on the face of the triangle
!   ----------------------------------------------------------------------
    FACET = .TRUE.
    DENOM = 1.0_dp / (VA + VB + VC)
    V     = VB * DENOM
    W     = VC * DENOM
    NT    = V1 + V*AB + W*AC

    200 DIJ = P - NT

END SUBROUTINE

SUBROUTINE SHORTD_TRI_EDGE(RIJ, UI, UJ, DIJ, NI, NJ, HLFL_I, HLFL_J)
    !====================================================================================================
    !   Subroutine to evaluate the shortest distance between line segmments
    !   RIJ = Vector connecting the geometrical centers of the line segmments
    !   UI  = Unitary vector definig the orientation of line segmment I
    !   UJ  = Unitary vector definig the orientation of line segmment J
    !   DIJ = Vector giving the shortest distance between the line segmments  
    !====================================================================================================
    USE COMMONS, ONLY: DP

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN)   :: RIJ(3), UI(3), UJ(3), HLFL_I, HLFL_J

    REAL(KIND=DP)               :: LAMI, LAMJ, LI, LJ
    REAL(KIND=DP)               :: RDOT, URI, URJ, UIJ, CC

    REAL(KIND=DP), INTENT(OUT)  :: DIJ(3), NI(3), NJ(3)

    DIJ  = 0.0_dp
    LAMI = 0.0_dp
    LAMJ = 0.0_dp

    RDOT = DOT_PRODUCT(RIJ, RIJ)
    URI  = DOT_PRODUCT(UI,  RIJ)
    URJ  = DOT_PRODUCT(UJ,  RIJ)
    UIJ  = DOT_PRODUCT(UI,  UJ)

    CC = (1.0_dp-UIJ**2)

    IF( CC < 1.e-06) THEN
        IF(URI /= 0.0_dp) THEN
            LAMI = SIGN(HLFL_I,URI)
            LAMJ = LAMI*UIJ - URJ
            IF(ABS(LAMJ) > HLFL_J) LAMJ = SIGN(HLFL_J,LAMJ) 
        ELSE
            LAMI = 0.0_dp
            LAMJ = 0.0_dp
        ENDIF
    ELSE
        LAMI = (URI - UIJ*URJ) / CC
        LAMJ = (UIJ*URI - URJ) / CC

        LI = ABS(LAMI) - HLFL_I
        LJ = ABS(LAMJ) - HLFL_J

        IF( (LI > 0.0_dp) .OR. (LJ > 0.0_dp) ) THEN 
            IF( LI > LJ ) THEN 
                LAMI = SIGN(HLFL_I, LAMI)
                LAMJ = LAMI*UIJ - URJ
                IF(ABS(LAMJ) > HLFL_J) LAMJ = SIGN(HLFL_J,LAMJ) 
            ELSE
                LAMJ = SIGN(HLFL_J, LAMJ)
                LAMI = LAMJ*UIJ + URI
                IF(ABS(LAMI) > HLFL_I) LAMI = SIGN(HLFL_I,LAMI) 
            ENDIF
        ENDIF
    ENDIF

    NI  = -LAMI*UI
    NJ  = LAMJ*UJ 
    DIJ = RIJ + NI + NJ

END SUBROUTINE SHORTD_TRI_EDGE

SUBROUTINE DEF_PATCHY_TRIANGLE()
!   ============================================================================================
!   ============================================================================================
!   Define the reference position of the vertices for the hard triangle particles and the
!   orientation the particle (which is taken to be normal to the face along the z-axis).
!   The vertices are defined using to the angle at the A vertex and the length of the AC and AB 
!   edges, taking the A vertex to initially be at (0,0,0), then repositioning so that the centre
!   of the triangle is at the origin.
!   For instance, if alpha=60 and AB=AC then an equilateral triangle will be defined. Instead,
!   if alpha<60 and AB=AC, an isosceles triangle will be defined where BC is shorter edge; and
!   if alpha>60 and AB=AC, an isosceles triangle will be defined where BC is longer edge.

!   We also define attractive patches along the edges of the triangle which run along the whole
!   length of the edge. 
!   ============================================================================================
!   ============================================================================================
    USE COMMONS, ONLY: DP, REFSITE, HLFPI, PI
    USE COMMONS, ONLY: TRI_ALPHA, TRI_B, TRI_C, PATCHYTRIT, TRI_THETAA, TRI_THETAB, TRI_THETAC, TRI_SIG2, TRI_EL
    USE COMMONS, ONLY: KFIJ, KFLAM2, TRI_DEL, KFAA, KFAB, KFBB, KFAC, KFBC, KFCC
    USE COMMONS, ONLY: KFAD, KFBD, KFCD, KFDD, KFAE, KFBE, KFCE, KFDE, KFEE
    USE COMMONS, ONLY: SPECIFICKFT, KFAF, KFBF, KFCF, KFDF, KFEF, KFFF
    USE COMMONS, ONLY: KFLAMA, KFLAMB, KFLAMC, RCUT, RCUTSQ, SQWT,  SQWDEL, SQWDEL2, SQWEPS
    USE COMMONS, ONLY: KFDELA, KFDELB, KFDELC, TRI_PHIA, TRI_PHIB, TRI_PHIC
    USE COMMONS, ONLY: POLYCHAINT, POLY_SIG, CP_SIG1, CP_SIG12, TRI_SIG, TP_CUT, TP_CUT2, CP_DEL

    USE ROTATIONS_MODULE, ONLY: CROSS, Q_TO_RM

    IMPLICIT NONE

    INTEGER       :: J1, VMAX
    REAL(KIND=DP) :: COM(3), CX, CY, EDGEL, AX(3), EDGE(3), RM(3,3), ANGLE, DV, DMAX
    REAL(KIND=DP) :: AXANG

    TRI_SIG2 = TRI_SIG*TRI_SIG
!   Orientation of triangle (normal to the face of the triangle)
    REFSITE(:,1)= [0.0_dp, 0.0_dp, 1.0_dp]
!   Position of vertex A 
    REFSITE(:,2)= [0.0_dp, 0.0_dp, 0.0_dp]
!   Position of vertex B 
    REFSITE(:,3)= [TRI_C,  0.0_dp, 0.0_dp]
!   Position of vertex C
    CX = TRI_B*COS(TRI_ALPHA*PI/180.0_dp)
    CY = TRI_B*SIN(TRI_ALPHA*PI/180.0_dp)
    REFSITE(:,4)= [ CX, CY, 0.0_dp]
!   Reposition the vertices so that the centre of the triangle is at the origin
    COM = 0.0_dp
    COM = (REFSITE(:,2) + REFSITE(:,3) + REFSITE(:,4))/3.0_dp
    REFSITE(:,2) = REFSITE(:,2) - COM
    REFSITE(:,3) = REFSITE(:,3) - COM
    REFSITE(:,4) = REFSITE(:,4) - COM

    DMAX = 0.0_dp
    DO J1 = 1, 3
        DV = NORM2(REFSITE(:,2))
        IF(DV > DMAX) THEN
            VMAX = J1
            DMAX = DV
        ENDIF
    ENDDO
    RCUT   = TRI_SIG + 2.0_dp*DMAX + 0.1_dp
    RCUTSQ = RCUT*RCUT

!   Find the half-lengths for each of the edges of the triangle 
    EDGEL     = NORM2(REFSITE(:,3) - REFSITE(:,2)) / 2.0_dp ! Edge A-B
    TRI_EL(1) = EDGEL
    EDGEL = NORM2(REFSITE(:,4) - REFSITE(:,3)) / 2.0_dp ! Edge B-C
    TRI_EL(2) = EDGEL
    EDGEL = NORM2(REFSITE(:,2) - REFSITE(:,4)) / 2.0_dp ! Edge C-A
    TRI_EL(3) = EDGEL

    IF(SQWT) THEN
        SQWDEL  = TRI_SIG + TRI_SIG * SQWDEL
        SQWDEL2 = SQWDEL*SQWDEL
        RCUT    = SQWDEL + 2.0_dp*DMAX + 0.1_dp
        RCUTSQ  = RCUT*RCUT
    ENDIF

    IF(PATCHYTRIT) THEN
    !   Get the centre of each of the edges on the triangle, this will be the
    !   position of the patch centres.
        REFSITE(:,5) = (REFSITE(:,3) + REFSITE(:,2))/2.0_dp ! Edge A-B
        REFSITE(:,6) = (REFSITE(:,4) + REFSITE(:,3))/2.0_dp ! Edge B-C
        REFSITE(:,7) = (REFSITE(:,2) + REFSITE(:,4))/2.0_dp ! Edge C-A

        ANGLE = -0.5_dp*TRI_THETAA*PI/180_dp
        EDGE  = REFSITE(:,3) - REFSITE(:,2)
        EDGE  = EDGE/NORM2(EDGE)
        AX    = CROSS(EDGE,REFSITE(:,1))
        AX    = AX / NORM2(AX)
        AXANG = 0.5_dp*(PI/2.0_dp-ACOS(DOT_PRODUCT((AX-REFSITE(:,5))/NORM2((AX-REFSITE(:,5))),EDGE)))
        RM    = Q_TO_RM ([COS(AXANG),SIN(AXANG)*REFSITE(1,1),SIN(AXANG)*REFSITE(2,1),SIN(AXANG)*REFSITE(3,1)])
        AX    = MATMUL(AX-REFSITE(:,5),RM)+REFSITE(:,5)
        RM    = Q_TO_RM ([COS(ANGLE),SIN(ANGLE)*EDGE(1),SIN(ANGLE)*EDGE(2),SIN(ANGLE)*EDGE(3)])
        REFSITE(:,8) = MATMUL(RM,AX-REFSITE(:,5))+REFSITE(:,5)
        REFSITE(:,8) = REFSITE(:,8) / NORM2(REFSITE(:,8))

        ANGLE = -0.5_dp*TRI_THETAB*PI/180_dp
        EDGE = REFSITE(:,4) - REFSITE(:,3)
        EDGE = EDGE/NORM2(EDGE)
        AX    = CROSS(EDGE,REFSITE(:,1))
        AX    = AX / NORM2(AX)
        AXANG = 0.5_dp*(PI/2.0_dp-ACOS(DOT_PRODUCT((AX-REFSITE(:,6))/NORM2((AX-REFSITE(:,6))),EDGE)))
        RM    = Q_TO_RM ([COS(AXANG),SIN(AXANG)*REFSITE(1,1),SIN(AXANG)*REFSITE(2,1),SIN(AXANG)*REFSITE(3,1)])
        AX    = MATMUL(AX-REFSITE(:,6),RM)+REFSITE(:,6)
        RM    = Q_TO_RM ([COS(ANGLE),SIN(ANGLE)*EDGE(1),SIN(ANGLE)*EDGE(2),SIN(ANGLE)*EDGE(3)])
        REFSITE(:,9) = MATMUL(RM,AX-REFSITE(:,6))+REFSITE(:,6)
        REFSITE(:,9) = REFSITE(:,9) / NORM2(REFSITE(:,9))

        ANGLE = -0.5_dp*TRI_THETAC*PI/180_dp
        EDGE = REFSITE(:,2) - REFSITE(:,4)
        EDGE = EDGE/NORM2(EDGE)
        AX    = CROSS(EDGE,REFSITE(:,1))
        AX    = AX / NORM2(AX)
        AXANG = 0.5_dp*(PI/2.0_dp-ACOS(DOT_PRODUCT((AX-REFSITE(:,7))/NORM2((AX-REFSITE(:,7))),EDGE)))
        RM    = Q_TO_RM ([COS(AXANG),SIN(AXANG)*REFSITE(1,1),SIN(AXANG)*REFSITE(2,1),SIN(AXANG)*REFSITE(3,1)])
        AX    = MATMUL(AX-REFSITE(:,7),RM)+REFSITE(:,7)
        RM    = Q_TO_RM ([COS(ANGLE),SIN(ANGLE)*EDGE(1),SIN(ANGLE)*EDGE(2),SIN(ANGLE)*EDGE(3)])
        REFSITE(:,10) = MATMUL(RM,AX-REFSITE(:,7))+REFSITE(:,7)
        REFSITE(:,10) = REFSITE(:,10) / NORM2(REFSITE(:,10))

        ! DO J1 = 1,10
        !     PRINT *, "C", REFSITE(:,J1)
        ! ENDDO

        ! STOP

        ALLOCATE( KFIJ(3,3), KFLAM2(3,3), TRI_DEL(2,3) )

        TRI_DEL(1,1)  = COS(KFDELA*PI/180_dp)
        TRI_DEL(1,2)  = COS(KFDELB*PI/180_dp)
        TRI_DEL(1,3)  = COS(KFDELC*PI/180_dp)
        TRI_DEL(2,1)  = COS(HLFPI-TRI_PHIA*PI/180_dp)
        TRI_DEL(2,2)  = COS(HLFPI-TRI_PHIB*PI/180_dp)
        TRI_DEL(2,3)  = COS(HLFPI-TRI_PHIC*PI/180_dp)

        KFIJ = 0.0_dp
        IF(SPECIFICKFT) THEN
            KFIJ(1,1) = KFAA; KFIJ(2,2) = KFBB; KFIJ(3,3) = KFCC
            KFIJ(1,2) = KFAB; KFIJ(2,1) = KFAB
            KFIJ(1,3) = KFAC; KFIJ(3,1) = KFAC
            KFIJ(2,3) = KFBC; KFIJ(3,2) = KFBC
        ELSE
            KFIJ(1,1) = KFAA
            KFIJ(1,2) = SQRT(KFAA*KFBB); KFIJ(2,1) = KFIJ(1,2)
            KFIJ(2,2) = KFBB
            KFIJ(1,3) = SQRT(KFAA*KFCC); KFIJ(3,1) = KFIJ(1,3)
            KFIJ(2,3) = SQRT(KFBB*KFCC); KFIJ(3,2) = KFIJ(2,3)
            KFIJ(3,3) = KFCC
        ENDIF
        
        KFLAMA = KFLAMA*TRI_SIG; KFLAMB = KFLAMB*TRI_SIG; KFLAMC = KFLAMC*TRI_SIG
        KFLAM2(1,1) = KFLAMA**2
        KFLAM2(1,2) = ( (KFLAMA + KFLAMB) / 2.0_dp )**2; KFLAM2(2,1) = ( (KFLAMB + KFLAMA) / 2.0_dp )**2
        KFLAM2(2,2) = KFLAMB**2
        KFLAM2(1,3) = ( (KFLAMA + KFLAMC) / 2.0_dp )**2; KFLAM2(3,1) = ( (KFLAMC + KFLAMA) / 2.0_dp )**2
        KFLAM2(2,3) = ( (KFLAMB + KFLAMC) / 2.0_dp )**2; KFLAM2(3,2) = ( (KFLAMC + KFLAMB) / 2.0_dp )**2
        KFLAM2(3,3) = KFLAMC**2
        

        IF(POLYCHAINT) THEN
            CP_SIG1  = (POLY_SIG+TRI_SIG)/2.0_dp
            CP_SIG12 = CP_SIG1**2
            TP_CUT   = CP_SIG1+CP_DEL
            TP_CUT2  = TP_CUT*TP_CUT
            IF(POLY_SIG > TRI_SIG) THEN
                RCUT   = CP_SIG1 + DMAX + 0.1_dp
                RCUTSQ = RCUT*RCUT
            ENDIF
        ENDIF
        
    ENDIF

END SUBROUTINE

SUBROUTINE VIEW_TRIANGLE()
        
    USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT, PATCHYTRIT, POLYCHAINT, N_POLY_TOT
    USE ROTATIONS_MODULE, ONLY: Q_TO_RM

    IMPLICIT NONE

    INTEGER:: J1, J2, CNTR
    REAL(KIND=DP) :: RWRITE(3,NPART), RM(3,3), RBCOORDS(3)

    IF(PATCHYTRIT) THEN
        CNTR = 3
    ELSE
        CNTR = 0
    ENDIF

    IF(POLYCHAINT) THEN
        WRITE(VIEWUNIT,*) N_POLY_TOT+(NPART-N_POLY_TOT)*(NSITES+1-CNTR)
    ELSE
        WRITE(VIEWUNIT,*) NPART*(NSITES+1-CNTR)
    ENDIF

    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO J1 = 1, NPART
        RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX(:))*BOX(:)
    END DO

    DO J1 = 1, NPART
        IF(POLYCHAINT .AND. J1<= N_POLY_TOT) THEN
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'H ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
        ELSE
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            RM = Q_TO_RM( Q(:,J1) )
            DO J2 = 1, NSITES-CNTR
                IF(J2==1) THEN
                    RBCOORDS = RWRITE(:,J1) + 0.25_dp*MATMUL(RM,REFSITE(:,J2))
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2<5) THEN
                    RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==5) THEN
                    RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Co ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==6) THEN
                    RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==7) THEN
                    RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Cl ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ENDIF
            ENDDO
        ENDIF
    END DO
    
END SUBROUTINE

! function tri_tri_intersect(p1, q1, r1, p2, q2, r2) result(intersects)
!     implicit none
!     real, dimension(3), intent(in) :: p1, q1, r1, p2, q2, r2
!     logical :: intersects
!     real, dimension(3) :: dir1, dir2, n1, n2, w, u, v
!     real :: det, dett, t, s, epsilon
!     integer :: i
    
!     epsilon = 1.0E-6   ! choose an appropriate epsilon value for your application
    
!     ! compute direction vectors for each triangle
!     dir1 = q1 - p1
!     dir2 = q2 - p2
    
!     ! compute the normal vectors of each triangle
!     n1 = cross_product(dir1, r1 - p1)
!     n2 = cross_product(dir2, r2 - p2)
    
!     ! compute the determinant of the system of equations
!     det = dot_product(n1, dir2)
    
!     ! check for parallel or coplanar triangles
!     if (abs(det) < epsilon) then
!       intersects = .false.
!       return
!     endif
    
!     ! compute the parameters for the intersection point
!     w = p1 - p2
!     dett = dot_product(n1, w)
    
!     ! compute the intersection point along the direction of triangle 2
!     t = dett / det
    
!     ! check if the intersection point is outside the bounds of triangle 2
!     if (t < 0.0 .or. t > 1.0) then
!       intersects = .false.
!       return
!     endif
    
!     ! compute the intersection point along the direction of triangle 1
!     u = cross_product(w, dir2)
!     s = dot_product(u, n2) / det
    
!     ! check if the intersection point is outside the bounds of triangle 1
!     if (s < 0.0 .or. s > 1.0) then
!       intersects = .false.
!       return
!     endif
    
!     ! check if the intersection point is inside the bounds of triangle 1
!     v = cross_product(dir1, w)
!     t = dot_product(v, n1) / det
    
!     ! check if the intersection point is between the face of one triangle
!     ! and the edge or vertex of the other triangle
!     if (t < 0.0 .or. t > 1.0) then
!       intersects = .false.
!       return
!     endif
    
!     intersects = .true.
    
!   contains
    
!     ! subroutine to compute the cross product of two vectors
!     function cross_product(a, b) result(c)
!       real, dimension(3), intent(in) :: a, b
!       real, dimension(3) :: c
!       c(1) = a(2)*b(3) - a(3)*b(2)
!       c(2) = a(3)*b(1) - a(1)*b(3)
!       c(3) = a(1)*b(2) - a(2)*b(1)
!     end function cross_product
    
!   end function tri_tri_intersect