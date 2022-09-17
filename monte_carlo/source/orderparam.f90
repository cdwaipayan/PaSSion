MODULE ORDERPARAM

    USE COMMONS,   ONLY: DP, PI, FORPI, NPART, NDIM, BOP_CUT, BOP_CUT_SQ, CELLLISTT, NPBOP
    USE CELL_LIST, ONLY: C_INDEX, NEIGHBOURS

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: GET_BOPS, LRGSTXCLSTR

CONTAINS

    SUBROUTINE GET_BOPS()
        
        USE COMMONS, ONLY: MIN_L, MAX_L, DEL_L, QLT, QBARLT, QLDOTQLT, GQLT, QL, QBARL, GQL, QLDOTQL, LRGSTXCLSTRT, QLDIST

        IMPLICIT NONE

        INTEGER                         :: J1, J2, NB(NPBOP)
        REAL(KIND=DP)                   :: COORDS(2,NPBOP,NPBOP), DIST(NPBOP,NPBOP)
        COMPLEX(KIND=DP), ALLOCATABLE   :: QLM(:,:)

        CALL GET_NEIGHBOURS(COORDS, DIST, NB)
        IF(LRGSTXCLSTRT) QLDIST = DIST

    !   Add loop here for multiple L values
        J2 = 0
        DO J1 = MIN_L, MAX_L, DEL_L
            
            J2 = J2 + 1
            
            IF(ALLOCATED(QLM)) THEN
                DEALLOCATE(QLM)
                ALLOCATE(QLM(-J1:J1,NPART))
            ELSE
                ALLOCATE(QLM(-J1:J1,NPART))
            ENDIF

            QLM = STEINHARDT(J1, COORDS, DIST, NB)

            IF(QLT)      QL(J2,:)        = LBOP_QL(J1, QLM)
            IF(QBARLT)   QBARL(J2,:)     = LBOP_QBARL(J1, QLM, DIST, NB)
            IF(QLDOTQLT) QLDOTQL(J2,:,:) = QL_DOT_QL(J1, QLM, DIST)
            IF(GQLT)     GQL(J2)         = GBOP_QL(J1, QLM)

        ENDDO

    END SUBROUTINE

    SUBROUTINE LRGSTXCLSTR(OUTPUTT, LRGST_X_CLSTR)
    !   ==================================================================
    !   Calculate the largest crystalline cluster. 
    !   ==================================================================

        USE COMMONS, ONLY: QBARLT, QLDOTQLT, QBARL, QLDOTQL
        USE COMMONS, ONLY: L_INDEX, QLMIN, QLMAX, QLMIN_2, QLMAX_2, TWOLIMS, QLDIST, BOP_CUT
        USE COMMONS, ONLY: XSTL_NN_MIN, XSTL_NN_MAX, XSTL_NN_MIN_1, XSTL_NN_MAX_1, XSTL_NN_MIN_2, XSTL_NN_MAX_2
    
        IMPLICIT NONE

        LOGICAL, INTENT(IN)         :: OUTPUTT

        INTEGER                     :: J1, J2, NBS(NPBOP), NBS2(NPBOP), XSTL_COUNT
        LOGICAL                     :: XBOND, XBOND2, XSTL_ADJ(NPBOP,NPBOP), CHECK1, CHECK2, CHECK3, CHECK4
        REAL(KIND=DP)               :: QLDQL

        INTEGER                     :: NCLSTRS
        INTEGER                     :: CLSTR(NPBOP), CLSTRSIZE(NPBOP), XBAR(NPBOP)

        INTEGER                     :: X_CLSTR, LRGST_X_ID

        INTEGER, INTENT(OUT)        :: LRGST_X_CLSTR

        NBS        = 0
        CLSTR      = 0
        NCLSTRS    = 0
        CLSTRSIZE  = 0
        XSTL_COUNT = 0
        XSTL_ADJ   = .FALSE.
        IF(TWOLIMS) NBS2 = 0

        IF(QBARLT .AND. QLDOTQLT) THEN
            STOP "SELECT ONLY ONE OF q_BAR_l or d_l"
        ELSE IF(QLDOTQLT)THEN
    !   ====================================================================================================
    !   Determine which particles belong to the largest crystalline cluster according to the number of 
    !   crystalline bonds they possess, which is measured using q_l(i)*q_l(j). A particle is ultimately 
    !   classified as being crystalline if it has a number of crystalline bonds which fit into some
    !   user defined limits.
    !   ====================================================================================================
    !   ----------------------------------------------------------------------------------------------------
    !   Establish how many crystalline bonds each particle has using q_l(i)*q_l(j), and build 
    !   an adjacency matrix for the crystalline particles.
    !   ----------------------------------------------------------------------------------------------------
            DO J1 = 1, NPBOP-1
                DO J2 = J1+1, NPBOP 
                !   Check if particles are close enough to form a bond.
                    IF(QLDIST(J1,J2) <= BOP_CUT) THEN
                        
                        XBOND = .FALSE.
                        IF(TWOLIMS) XBOND2 = .FALSE.
                        QLDQL = QLDOTQL(L_INDEX,J1,J2)

                    !   Check if bond-order parameter satisfies any conditions for a crystalline bond
                        IF(QLDQL > QLMIN .AND. QLDQL < QLMAX) THEN
                            XBOND = .TRUE.
                        ELSEIF(TWOLIMS) THEN
                            IF(QLDQL > QLMIN_2 .AND. QLDQL < QLMAX_2) XBOND2 = .TRUE.
                        ENDIF    

                    !   If i and j share a crystalline bond update their respective number of nearest neighbours
                    !   and update the adjacency matrix to account for their bond.
                        IF(XBOND .OR. XBOND2) THEN
                            IF(XBOND) THEN
                                NBS(J1) = NBS(J1) + 1
                                NBS(J2) = NBS(J2) + 1
                            ELSE
                                NBS2(J1) = NBS2(J1) + 1
                                NBS2(J2) = NBS2(J2) + 1
                            ENDIF
                            XSTL_ADJ(J1,J2) = .TRUE.
                            XSTL_ADJ(J2,J1) = .TRUE.
                        ENDIF

                    ENDIF
                ENDDO ! Loop over particles j
            ENDDO ! Loop over particles i

    !   ----------------------------------------------------------------------------------------------------
    !   Establish the number of "crystalline" particles, where we define a crystalline 
    !   particle as one which has at least XSTL_NN_MIN crystalline bonds and no more than 
    !   XSTL_NN_MAX crystalline bonds.
    !   If a particle is not identified as crystalline make sure to remove all bonds
    !   associated with that particle from the adjacency matrix, as it cannot be a part
    !   of a crystalline cluster.
    !   ----------------------------------------------------------------------------------------------------
            DO J1 = 1, NPBOP
                CHECK1 = NBS(J1) >= XSTL_NN_MIN .AND. NBS(J1) <= XSTL_NN_MAX
                IF( TWOLIMS .AND. (.NOT. CHECK1) ) THEN
                    CHECK2 = NBS(J1) >= XSTL_NN_MIN_1 .AND. NBS(J1) <= XSTL_NN_MAX_1
                    CHECK3 = NBS2(J1) >= XSTL_NN_MIN_2 .AND. NBS2(J1) <= XSTL_NN_MAX_2
                    CHECK4 = CHECK2 .AND. CHECK3
                ELSE
                    CHECK2 = .FALSE.
                ENDIF
                
                IF( CHECK1 .OR. CHECK4 ) THEN
                    XSTL_COUNT = XSTL_COUNT + 1
                ELSE
                    XSTL_ADJ(J1,:) = .FALSE.
                    XSTL_ADJ(:,J1) = .FALSE.
                ENDIF
            ENDDO

    !   ----------------------------------------------------------------------------------------------------
    !   Use a depth-first search algorithm to identify all clusters of crystalline particles
    !   using the adjacency matrix constructed from the previous steps. Only pass in the particles
    !   that have been labelled as crystalline so that crystalline particles not in the largest
    !   crystalline cluster can be distinguished from fluid/liquid particles.
    !   ----------------------------------------------------------------------------------------------------
            DO J1 = 1, NPBOP
                CHECK1 = NBS(J1) >= XSTL_NN_MIN .AND. NBS(J1) <= XSTL_NN_MAX
                IF( TWOLIMS .AND. (.NOT. CHECK1) ) THEN
                    CHECK2 = NBS(J1) >= XSTL_NN_MIN_1 .AND. NBS(J1) <= XSTL_NN_MAX_1
                    CHECK3 = NBS2(J1) >= XSTL_NN_MIN_2 .AND. NBS2(J1) <= XSTL_NN_MAX_2
                    CHECK4 = CHECK2 .AND. CHECK3
                ELSE
                    CHECK2 = .FALSE.
                ENDIF

                IF( CHECK1 .OR. CHECK4 ) THEN
                    NCLSTRS = NCLSTRS + 1
                    CALL DFS(J1, NCLSTRS, XSTL_ADJ, CLSTR, X_CLSTR)
                    CLSTRSIZE(NCLSTRS) = X_CLSTR
                ENDIF
            ENDDO

        ELSEIF(QBARLT) THEN
    !   ====================================================================================================
    !   Determine which particles belong to the largest crystalline cluster using q_l_bar. Two checks are 
    !   performed in order to classify a particle as being crystalline:
    !       1) If it has a q_l_bar value fits in between some user defined limits.
    !       2) If it has N neighbours also labelled as crystalline, where N_min <= N <= N_max.
    !   ====================================================================================================
    !   ----------------------------------------------------------------------------------------------------
    !   Establish which particles are crystalline according to check (1).
    !   ----------------------------------------------------------------------------------------------------
            XBAR = 0
            DO J1 = 1, NPBOP
                QLDQL = QBARL(L_INDEX,J1)
                IF(QLDQL > QLMIN .AND. QLDQL < QLMAX) THEN
                    XBAR(J1) = 1
                ELSEIF(TWOLIMS) THEN
                    IF(QLDQL > QLMIN_2 .AND. QLDQL < QLMAX_2) XBAR(J1) = 1
                ENDIF
            ENDDO
    !   ----------------------------------------------------------------------------------------------------
    !   Establish which crystalline particles from step (1) are nearest neighbours, according to a cutoff.
    !   Populate the crystalline adjaceny matrix with these neighbours.
    !   ----------------------------------------------------------------------------------------------------
            DO J1 = 1, NPBOP-1
                IF( XBAR(J1) == 0 ) CYCLE
                DO J2 = J1 + 1, NPBOP
                    IF( XBAR(J2) == 1 ) THEN
                        IF(QLDIST(J1,J2) <= BOP_CUT) THEN
                            NBS(J1) = NBS(J1) + 1
                            NBS(J2) = NBS(J2) + 1
                            XSTL_ADJ(J1,J2) = .TRUE.
                            XSTL_ADJ(J2,J1) = .TRUE.
                        ENDIF
                    ENDIF
                ENDDO
            ENDDO
    !   ----------------------------------------------------------------------------------------------------
    !   Cycle through each particle and check whether it has satisfied check (2) to remain labelled as 
    !   crystalline. If not prune the crystalline adjaceny matrix accordingly.
    !   ----------------------------------------------------------------------------------------------------
            DO J1 = 1, NPBOP
                IF(XBAR(J1) == 1) THEN
                    IF(NBS(J1) < XSTL_NN_MIN .OR. NBS(J1) > XSTL_NN_MAX) THEN
                        XBAR(J1) = 0
                        XSTL_ADJ(J1,:) = .FALSE.
                        XSTL_ADJ(:,J1) = .FALSE.
                    ELSE
                        XSTL_COUNT = XSTL_COUNT + 1
                    ENDIF
                ENDIF
            ENDDO

    !   ----------------------------------------------------------------------------------------------------
    !   Use a depth-first search algorithm to identify all clusters of crystalline particles
    !   using the adjacency matrix constructed from the previous steps. Only pass in the particles
    !   that have been labelled as crystalline so that crystalline particles not in the largest
    !   crystalline cluster can be distinguished from fluid/liquid particles.
    !   ----------------------------------------------------------------------------------------------------
            DO J1 = 1, NPBOP
                IF(XBAR(J1) == 1) THEN
                    NCLSTRS = NCLSTRS + 1
                    CALL DFS(J1, NCLSTRS, XSTL_ADJ, CLSTR, X_CLSTR)
                    CLSTRSIZE(NCLSTRS) = X_CLSTR
                ENDIF
            ENDDO

        ELSE
            STOP "MUST SELECT ONE OF q_BAR_l or d_l"
        ENDIF

    !   ====================================================================================================
    !   Identify the size and ID of the largest cluster.
    !   ====================================================================================================
        LRGST_X_CLSTR = 0
        DO J1 = 1, NCLSTRS
            IF(CLSTRSIZE(J1) > LRGST_X_CLSTR) THEN
                LRGST_X_ID = J1
                LRGST_X_CLSTR = CLSTRSIZE(J1)
            ENDIF
        ENDDO
    !   ====================================================================================================

        IF(OUTPUTT) CALL VIEW_LRGSTXCLSTR(CLSTR,LRGST_X_CLSTR,LRGST_X_ID)
        
    END SUBROUTINE LRGSTXCLSTR

    SUBROUTINE VIEW_LRGSTXCLSTR(CLSTR,LRGST_X_CLSTR,LRGST_X_ID)

!       ==================================================================
!       Visualise the largest crystalline cluster in the system and output
!       the positions of the particles in the largest crystalline cluster.
!       ==================================================================

        USE COMMONS, ONLY: R, NSITES, HLFLNGTH, Q, REFSITE, RIGIDT, ETPT, BOPCLSTRT, CLSTRCOM, PCLSTRID, BINARYT, BNRYRTIOT, BNRYA
        USE ROTATIONS_MODULE, ONLY: Q_TO_RM

        INTEGER, INTENT(IN)         :: CLSTR(NPBOP), LRGST_X_CLSTR, LRGST_X_ID

        INTEGER                     :: J1, J2
        REAL(KIND=DP)               :: RBCOORDS(3), RM(3,3)

        CHARACTER(len=2)            :: PID

        OPEN (UNIT = 173, FILE ='pos_x_clstr.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
        OPEN (UNIT = 174, FILE ='lrgst_x_clstr.xyz', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
        
        IF(BOPCLSTRT) THEN
            WRITE (174,*) NPBOP+NPART
            WRITE (174,*) "Largest Crystalline Cluster Size:",LRGST_X_CLSTR
            
            DO J1 = 1, NPBOP
            !   Particles in the largest crystalline cluster
                IF(CLSTR(J1) == LRGST_X_ID) THEN
                    WRITE (173,'(3F12.7)') CLSTRCOM(:,J1)
                    WRITE (174,'(A5,1X,3F12.7)') "O", CLSTRCOM(:,J1)
            !   Crystalline particles which are not in the largest crystalline cluster
                ELSEIF(CLSTR(J1) /= 0) THEN
                    WRITE (174,'(A5,1X,3F12.7)') "N", CLSTRCOM(:,J1)
            !   Non-crystalline particles
                ELSE
                    WRITE (174,'(A5,1X,3F12.7)') "C", CLSTRCOM(:,J1)
                ENDIF
            ENDDO
        !   Add rigid body sites to particles
            DO J1 = 1, NPART
                ! PRINT *, J1, CLSTR(PCLSTRID(J1)), PCLSTRID(J1), NPBOP
            !   Particles in the largest crystalline cluster
                IF(CLSTR(PCLSTRID(J1)) == LRGST_X_ID) THEN
                    WRITE(174,'(A5,1X,3F12.7)') "H", R(:,J1)
            !   Crystalline particles which are not in the largest crystalline cluster
                ELSEIF(CLSTR(PCLSTRID(J1)) /= 0) THEN
                    WRITE(174,'(A5,1X,3F12.7)') "He", R(:,J1)
            !   Non-crystalline particles
                ELSE
                    WRITE(174,'(A5,1X,3F12.7)') "Li", R(:,J1)
                ENDIF
            ENDDO
        ELSE
            WRITE (174,*) NPART*(NSITES+1)
            WRITE (174,*) "Largest Crystalline Cluster Size:",LRGST_X_CLSTR
            
            DO J1 = 1, NPART
            !   Particles in the largest crystalline cluster
                IF(CLSTR(J1) == LRGST_X_ID) THEN
                    WRITE (173,'(3F12.7)') R(:,J1)

                    IF(BINARYT) THEN
                        IF(BNRYRTIOT) THEN
                            IF( J1<=BNRYA ) THEN
                                WRITE (174,'(A5,1X,3F12.7)') "Kr", R(:,J1)
                            ELSE
                                WRITE (174,'(A5,1X,3F12.7)') "O", R(:,J1)
                            ENDIF
                        ELSE
                            IF( MOD(J1,2) == 0 ) THEN
                                WRITE (174,'(A5,1X,3F12.7)') "O", R(:,J1)
                            ELSE
                                WRITE (174,'(A5,1X,3F12.7)') "Kr", R(:,J1)
                            ENDIF
                        ENDIF
                    ELSE
                        WRITE (174,'(A5,1X,3F12.7)') "O", R(:,J1)
                    ENDIF
                    PID = "H"
            !   Crystalline particles which are not in the largest crystalline cluster
                ELSEIF(CLSTR(J1) /= 0) THEN
                    WRITE (174,'(A5,1X,3F12.7)') "N", R(:,J1)
                    PID = "He"
            !   Non-crystalline particles
                ELSE
                    WRITE (174,'(A5,1X,3F12.7)') "C", R(:,J1)
                    PID = "Li"
                ENDIF
            !   Add rigid body sites to particles
                IF(RIGIDT) THEN
                    DO J2 = 1, NSITES
                        RM = Q_TO_RM( Q(:,J1) )
                        IF(ETPT .AND. HLFLNGTH /= 0.0_dp) THEN
                            RBCOORDS = R(:,J1) + HLFLNGTH*MATMUL(RM,REFSITE(:,J2))
                        ELSE
                            RBCOORDS = R(:,J1) + 0.5_dp*MATMUL(RM,REFSITE(:,J2))
                        ENDIF
                        WRITE(174,'(A5,1X,3F12.7)') TRIM(ADJUSTL(PID)), RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                    ENDDO
                ENDIF
            ENDDO
        ENDIF

        CLOSE (UNIT = 173, STATUS = 'KEEP')
        CLOSE (UNIT = 174, STATUS = 'KEEP')

    END SUBROUTINE VIEW_LRGSTXCLSTR

    SUBROUTINE GET_NEIGHBOURS(COORDS, DIST, NB)

!       ==================================================================
!       Extract the nearest neighbours of each particle in the system, and
!       the required geometric parameters for the calculation of the 
!       Steinhardt Bond-Orientational Order parameters.
!       ==================================================================

        USE COMMONS, ONLY: R, BOX, BOPCLSTRT, CLSTRCOM
        USE CLUSTER_MOVE, ONLY: BUILD_ALL_CLUSTERS

        IMPLICIT NONE

        INTEGER                     :: J1, J2, J3
        INTEGER                     :: CI(3), J_LIST(NPART), JJ
        REAL(KIND=DP)               :: DIJ, RI(NDIM), RJ(NDIM), RIJ(NDIM), RJI(NDIM), RIJSQ
        REAL(KIND=DP)               :: PHI_IJ, PHI_JI, CTHETA_IJ, CTHETA_JI
        REAL(KIND=DP), INTENT(OUT)  :: DIST(NPBOP,NPBOP), COORDS(2,NPBOP,NPBOP)
        INTEGER, INTENT(OUT)        :: NB(NPBOP)

        COORDS      = 0.0_dp
        DIST(:,:)   = 1.e+10
        NB(:)       = 0

        IF(BOPCLSTRT) CALL BUILD_ALL_CLUSTERS(BOPSNT=.TRUE.)
        
        DO J1 = 1, NPBOP
        !   Position of particle I
            IF(BOPCLSTRT) THEN
                RI = CLSTRCOM(:,J1)
            ELSE
                RI  = R(:,J1)
            ENDIF

            IF(CELLLISTT .AND. (.NOT. BOPCLSTRT)) THEN
            !   Only extract the index of particles from 13 neighbouring cells
            !   for computational efficiency.
                CI = C_INDEX ( RI/BOX )
                J_LIST = NEIGHBOURS ( NPART, J1, CI, .TRUE. )
                JJ = 0
            ELSE
                JJ = J1
            ENDIF
    
            DO
                JJ = JJ + 1                     ! Next entry
    
                IF(CELLLISTT .AND. (.NOT. BOPCLSTRT)) THEN
                    J2 = J_LIST(JJ)             ! Get neighbour index
                    IF ( J2 == 0 ) EXIT         ! List exhausted
                ELSE
                    J2 = JJ
                    IF ( J2 == NPBOP+1 ) EXIT     ! List exhausted
                    IF ( J1 == J2 ) CYCLE
                ENDIF
    
            !   Position of particle J
                IF(BOPCLSTRT) THEN
                    RJ = CLSTRCOM(:,J2)
                ELSE
                    RJ    = R(:,J2)
                ENDIF

                RIJ   = RI - RJ
                RIJ   = RIJ - BOX*ANINT( RIJ/BOX )
                RIJSQ = DOT_PRODUCT(RIJ,RIJ)

                IF(RIJSQ <= BOP_CUT_SQ) THEN
                !   Compute the reverse separation vector for particle J 
                    RJI         = -RIJ
                    DIJ         = SQRT(RIJSQ)

                    DIST(J2,J1) = DIJ
                    DIST(J1,J2) = DIJ

                !   Update the number of neighbours for particles I and J                    
                    NB(J1)      = NB(J1) + 1
                    NB(J2)      = NB(J2) + 1

                !   Check for artificial negative coordiates which may alter the sign of our
                !   complex numbers later down the line.
                    DO J3 = 1, NDIM
                        IF(RIJ(J3) == -0.0_dp) RIJ(J3) = 0.0_dp
                        IF(RJI(J3) == -0.0_dp) RJI(J3) = 0.0_dp
                    ENDDO

                !   Calculate the cosine of the polar angle for RIJ and RJI
                    CTHETA_IJ = RIJ(3) / DIJ
                    CTHETA_JI = RJI(3) / DIJ
                
                !   Calculate the azimuthal angle for RIJ and RJI, making sure to add PI so 
                !   that 0 <= phi <= 2pi (as DATAN2 returns a number between -pi and pi).
                    PHI_IJ = ATAN2( RIJ(2), RIJ(1) ) + PI
                    PHI_JI = ATAN2( RJI(2), RJI(1) ) + PI

                !   Save the spherical coordinates associated with the bond between particles I and J
                    COORDS(1,J2,J1) = CTHETA_IJ
                    COORDS(1,J1,J2) = CTHETA_JI

                    COORDS(2,J2,J1) = PHI_IJ
                    COORDS(2,J1,J2) = PHI_JI
                ENDIF

            ENDDO
        ENDDO

    END SUBROUTINE

    FUNCTION LBOP_QL(L, QLM) RESULT ( QL ) 
!   ==================================================================
!   Calculate the rotationally invariant local bond order parameter
!   (aka q-L) for each particle.
!   ==================================================================

        IMPLICIT NONE

        INTEGER, INTENT(IN)             :: L
        COMPLEX(KIND=DP), INTENT(IN)    :: QLM(-L:L, NPBOP)

        INTEGER                         :: J1, M
        REAL(KIND=DP)                   :: QL(NPBOP)

        QL = 0.0_dp

        DO J1 = 1, NPBOP
                
            DO M = -L, L
                QL(J1) = QL(J1) + ABS(QLM(M,J1))**2
            ENDDO

        ENDDO

        QL = SQRT(FORPI/REAL((2*L+1),DP) * QL)

    END FUNCTION LBOP_QL

    FUNCTION LBOP_QBARL(L,QLM,DIST,NB) RESULT ( QBARL ) 
!   ==================================================================
!   Calculate the rotationally invariant local bond order parameter
!   (aka q-L) for each particle.
!   ==================================================================

        IMPLICIT NONE

        INTEGER, INTENT(IN)             :: L, NB(NPBOP)
        REAL(KIND=DP), INTENT(IN)       :: DIST(NPBOP,NPBOP)
        COMPLEX(KIND=DP), INTENT(IN)    :: QLM(-L:L, NPBOP)

        INTEGER                         :: J1, J2, M
        REAL(KIND=DP)                   :: DIJ
        COMPLEX(KIND=DP)                :: QBARLM(-L:L, NPBOP)

        REAL(KIND=DP)                   :: QBARL(NPBOP)

        QBARLM = 0.0_dp
        QBARL  = 0.0_dp

        DO J1 = 1, NPBOP
            
            DO J2 = 1, NPBOP
                
                DIJ = DIST(J2,J1)
                
                IF (J1 /= J2 .AND. DIJ > BOP_CUT) CYCLE
                
                DO M = -L, L
                    QBARLM(M, J1) = QBARLM(M, J1) + QLM(M,J2)
                ENDDO

            ENDDO
            
            QBARLM(:, J1) = QBARLM(:,J1) / REAL((1+NB(J1)),DP)
            
            DO M = -L, L
                QBARL(J1) = QBARL(J1) + ABS( QBARLM(M,J1) )**2
            ENDDO

        ENDDO

        QBARL = SQRT(FORPI/REAL((2*L+1),DP) * QBARL)

    END FUNCTION LBOP_QBARL

    FUNCTION GBOP_QL(L,QLM) RESULT ( GQL ) 
!   ==================================================================
!   Calculate the rotationally invariant global bond order parameter
!   (aka Q-L) for each particle.
!   ==================================================================

        IMPLICIT NONE
        INTEGER, INTENT(IN)             :: L
        COMPLEX(KIND=DP), INTENT(IN)    :: QLM(-L:L, NPBOP)

        INTEGER                         :: J1, M
        COMPLEX(KIND=DP)                :: GQLM(-L:L)
        REAL(KIND=DP)                   :: GQL

        GQLM = 0.0_dp
        GQL  = 0.0_dp

        DO J1 = 1, NPBOP
            DO M = -L, L
                GQLM(M) = GQLM(M) + QLM(M,J1)
            ENDDO
        ENDDO

        GQLM = GQLM / REAL(NPBOP,DP)

        DO M = -L, L
            GQL = GQL + ABS(GQLM(M))**2
        ENDDO

        GQL = SQRT(FORPI/REAL((2*L+1),DP) * GQL)

    END FUNCTION GBOP_QL

    FUNCTION QL_DOT_QL(L,QLM,DIST) RESULT ( QLDOTQL ) 
!   ==================================================================
!   Calculate which particles are to be classified as translationally
!   ordered (or solid-like) using the translational-order correlation
!   between each particle and its neighbours. 
!   ==================================================================

        IMPLICIT NONE

        INTEGER, INTENT(IN)             :: L
        COMPLEX(KIND=DP), INTENT(IN)    :: QLM(-L:L, NPBOP)
        REAL(KIND=DP), INTENT(IN)       :: DIST(NPBOP,NPBOP)

        INTEGER                         :: J1, J2
        REAL(KIND=DP)                   :: QLI_QLJ, QL_I, QL_J, DLIJ
        REAL(KIND=DP)                   :: DIJ

        REAL(KIND=DP)                   :: QLDOTQL(NPBOP,NPBOP)

        QLDOTQL = 0.0_dp

        DO J1 = 1, NPBOP-1

            QL_I = SQRT( SUM( ABS(QLM(:,J1))**2 ) )

            DO J2 = J1+1, NPBOP
                
                DIJ = DIST(J2,J1)
                IF (DIJ > BOP_CUT) CYCLE

                QL_J = SQRT( SUM( ABS(QLM(:,J2))**2 ) )
                
                QLI_QLJ = REALPART( SUM( QLM(:,J1) * CONJG(QLM(:,J2)) ) )
                DLIJ = QLI_QLJ / (QL_I*QL_J)

                QLDOTQL(J1,J2) = DLIJ
                QLDOTQL(J2,J1) = DLIJ
                
            ENDDO ! Loop over particles j
        ENDDO ! Loop over particles i

    END FUNCTION QL_DOT_QL

    FUNCTION STEINHARDT(L, COORDS, DIST, NB) RESULT ( QLM ) 
!   ==================================================================
!   Calculate the Steinhardt order parameter for each particle.
!   ==================================================================

        IMPLICIT NONE

        INTEGER, INTENT(IN)         :: L, NB(NPBOP)
        REAL(KIND=DP), INTENT(IN)   :: COORDS(2,NPBOP,NPBOP), DIST(NPBOP,NPBOP)

        INTEGER                     :: J1, J2
        INTEGER                     :: M
        REAL(KIND=DP)               :: DIJ, LM
        REAL(KIND=DP)               :: CTHETA_IJ, CTHETA_JI, PHI_IJ, PHI_JI, PLMI, PLMJ
        COMPLEX(KIND=DP)            :: YLM, YLMI, YLMJ, YLMMI, YLMMJ
        COMPLEX(KIND=DP), PARAMETER :: i = (0.0_dp, 1.0_dp)

        COMPLEX(KIND=DP)            :: QLM(-L:L, NPBOP)

        QLM = 0.0_dp

        DO J1 = 1, NPBOP-1
            DO J2 = J1+1, NPBOP
            
            !   Only consider pairs whose center-to-center distance are within the cutoff
                DIJ = DIST(J2,J1)
                IF (DIJ > BOP_CUT) CYCLE

                CTHETA_IJ = COORDS(1,J2,J1)
                CTHETA_JI = COORDS(1,J1,J2)
                
                PHI_IJ    = COORDS(2,J2,J1)
                PHI_JI    = COORDS(2,J1,J2)

                DO M = 0, L
                    PLMI = PLGNDR(L, M,  CTHETA_IJ)
                    PLMJ = PLGNDR(L, M,  CTHETA_JI)

                    LM   = FCTRL_FRAC( (L-M), (L+M) )
                    
                    YLM  = SQRT( REAL((2*L+1),DP)/FORPI * LM )

                    YLMI = YLM * PLMI * EXP(i*M*PHI_IJ)
                    YLMJ = YLM * PLMJ * EXP(i*M*PHI_JI)
                    
                    QLM(M,J1) = QLM(M,J1) + YLMI
                    QLM(M,J2) = QLM(M,J2) + YLMJ

                    IF(M > 0) THEN
                        YLMMI       = (-1)**M * CONJG(YLMI)
                        YLMMJ       = (-1)**M * CONJG(YLMJ)

                        QLM(-M,J1) = QLM(-M,J1) + YLMMI
                        QLM(-M,J2) = QLM(-M,J2) + YLMMJ
                    ENDIF
                ENDDO
                
            ENDDO ! Loop over particles j
        ENDDO ! Loop over particles i


        DO J1 = 1, NPBOP
            IF(NB(J1) > 0) THEN
                QLM(:,J1) = QLM(:,J1) / REAL(NB(J1),DP)
            ELSE
        !   If a particle has no neighbours set the value to zero
                QLM(:,J1) = 0.0_dp
            ENDIF
        ENDDO

    END FUNCTION STEINHARDT

    FUNCTION PLGNDR(L, M, X) RESULT ( PLM ) 
!   ==================================================================
!   Calculate the associated Legendre polynomial P(l,m,x).
!   Here m and l are integers satisfying 0 <= m <= l. 
!   x lies in the range -1 <= x <= 1 (it is the cosine of the polar 
!   angle between particles i and j for which the spherical harmonic 
!   is being computed for).
!   Taken from: 
!   NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING
!   ==================================================================

        IMPLICIT NONE

        INTEGER, INTENT(IN)         :: L, M
        REAL(KIND=DP), INTENT(IN)   :: X

        INTEGER                     :: I, LL
        REAL(KIND=DP)               :: FACT, PLL, PMM, PMMP1, SOMX2

        REAL(KIND=DP)               :: PLM

        IF(M < 0 .OR. M > L .OR. ABS(X) > 1.0_dp) THEN
            PRINT *, "BAD ARGUMENTS IN PLGNDR"
            STOP
        ENDIF

        PMM = 1.0_dp

        IF(M > 0) THEN
            SOMX2 = SQRT( (1.0_dp-X)*(1.0_dp+X) )
            FACT  = 1.0_dp

            DO I = 1, M
                PMM = -PMM * FACT * SOMX2
                FACT = FACT + 2.0_dp
            ENDDO
        ENDIF

        IF(L == M) THEN
            PLM = PMM
        ELSE
            PMMP1 = X * REAL((2*M+1),DP) * PMM
            IF(L == M+1) THEN
                PLM = PMMP1
            ELSE
                DO LL = M+2, L
                    PLL   = ( X*REAL((2*LL-1),DP)*PMMP1-REAL((LL+M-1),DP)*PMM ) / REAL((LL-M),DP)
                    PMM   = PMMP1
                    PMMP1 = PLL
                ENDDO
                PLM = PLL
            ENDIF
        ENDIF

    END FUNCTION PLGNDR

    FUNCTION FCTRL_FRAC(A,B) RESULT ( F ) 
!   ==================================================================
!   Calculate F = A!/B!, by first computing B!/A! and then returning
!   the its inverse. Computing the fraction directly is not only 
!   computationally more efficient than explicitly evaluating each of
!   the factorials and then calculating the fraction, but it also
!   allows us to compute F for larger values of A! and B! as we avoid
!   numerical overflow.
!   ==================================================================
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: A, B

        INTEGER             :: J1, END
        REAL(KIND=DP)       :: F

        F = 1.0_dp

        IF(A == 0 .OR. A == 1) THEN
            END = 1
        ELSE
            END = A + 1
        ENDIF
        

        DO J1 = B, END, -1
            F = F * REAL(J1,DP)
        ENDDO

        F = 1.0_dp / F
    
    END FUNCTION FCTRL_FRAC

    SUBROUTINE DFS(V, CLSTR_ID, ADJ, CLSTR, CLSTRSIZE)
!   =======================================================================
!   Use a non-recursive implementation of a depth-first search algorithm 
!   to identify connected components of a graph. 
!   Algorithm taken from https://en.wikipedia.org/wiki/Depth-first_search
!   -----------------------------------------------------------------------
!   PARAMETER DEFINITIONS:
!   * V is the starting node from which to begin traversing the graph 
!   under consideration
!   * CLSTR_ID is the index of the current cluster under consideration
!   * ADJ is the adjacency matrix of the entire graph under consideration.
!   * CLSTR is an array holding the ID of the cluster associated 
!     with each particle
!   * CLSTRSIZE is an array holding the size of each cluster
!   =======================================================================
        USE STACK 

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: V, CLSTR_ID
        INTEGER :: CLSTR(NPBOP)
        LOGICAL :: ADJ(NPBOP,NPBOP)
        INTEGER, INTENT(OUT) :: CLSTRSIZE
        INTEGER :: U, J2

        TYPE(STACK_VAR) :: S

        U = V
        CALL PUSH(S,U)

        CLSTRSIZE = 0

        DO 
            IF( EMPTY(S) ) EXIT
            U = POP(S)

            IF(CLSTR(U) == 0) THEN
                CLSTR(U) = CLSTR_ID
                CLSTRSIZE = CLSTRSIZE + 1
                DO J2 = 1, NPBOP
                    IF ( ADJ(U,J2) ) THEN
                        CALL PUSH(S,J2)
                    ENDIF
                ENDDO
            ENDIF

        ENDDO

    END SUBROUTINE DFS

END MODULE ORDERPARAM