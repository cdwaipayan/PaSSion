MODULE ORDERPARAM

    USE COMMONS,   ONLY: DP, PI, FORPI, NPART, NDIM, BOP_CUT, BOP_CUT_SQ, CELLLISTT
    USE CELL_LIST, ONLY: C_INDEX, NEIGHBOURS

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: GET_BOPS, LRGSTXCLSTR

CONTAINS

    SUBROUTINE GET_BOPS()
        
        USE COMMONS, ONLY: MIN_L, MAX_L, DEL_L, QLT, QBARLT, QLDOTQLT, GQLT, QL, QBARL, GQL, QLDOTQL

        IMPLICIT NONE

        INTEGER                         :: J1, J2, NB(NPART)
        REAL(KIND=DP)                   :: COORDS(2,NPART,NPART), DIST(NPART,NPART)
        COMPLEX(KIND=DP), ALLOCATABLE   :: QLM(:,:)

        CALL GET_NEIGHBOURS(COORDS, DIST, NB)

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

        USE COMMONS, ONLY: R, L_INDEX, QLDOTQL, QLMIN, QLMAX, XSTL_NN
    
        IMPLICIT NONE

        LOGICAL, INTENT(IN)         :: OUTPUTT

        INTEGER                     :: J1, J2, NBS(NPART), XSTL_COUNT
        LOGICAL                     :: XSTL_ADJ(NPART,NPART)

        INTEGER                     :: NCLSTRS
        INTEGER                     :: CLSTR(NPART), CLSTRSIZE(NPART)

        INTEGER                     :: X_CLSTR, LRGST_X_ID

        INTEGER, INTENT(OUT)        :: LRGST_X_CLSTR

        NBS        = 0
        CLSTR      = 0
        NCLSTRS    = 0
        XSTL_COUNT = 0
        XSTL_ADJ   = .FALSE.

    !   Establish how many crystalline bonds each particle has using q_l(i)*q_l(j), and build 
    !   an adjacency matrix for the crystalline particles.
        DO J1 = 1, NPART-1
            DO J2 = J1+1, NPART 
                IF(QLDOTQL(L_INDEX,J1,J2) > QLMIN .AND. QLDOTQL(L_INDEX,J1,J2) < QLMAX) THEN
                    NBS(J1) = NBS(J1) + 1
                    NBS(J2) = NBS(J2) + 1
                    XSTL_ADJ(J1,J2) = .TRUE.
                    XSTL_ADJ(J2,J1) = .TRUE.
                ENDIF    
            ENDDO ! Loop over particles j
        ENDDO ! Loop over particles i

    !   Establish the number of "crystalline" particles, where we define a crystalline 
    !   particle as 1.0_dp which has at least XSTL_NN crystalline bonds.
    !   If a particle is not identified as crystalline make sure to remove all bonds
    !   associated with that particle from the adjacency matrix, as it cannot be a part
    !   of a crystalline cluster.
        DO J1 = 1, NPART
            IF(NBS(J1) >= XSTL_NN) THEN
                XSTL_COUNT = XSTL_COUNT + 1
            ELSE
                XSTL_ADJ(J1,:) = .FALSE.
                XSTL_ADJ(:,J1) = .FALSE.
            ENDIF
        ENDDO

        NCLSTRS    = 0
        CLSTR      = 0
        CLSTRSIZE  = 0
    !   Use a depth-first search algorithm to identify all clusters of crystalline particles
    !   using the adjacency matrix constructed from the previous steps.
        DO J1 = 1, NPART
            IF (CLSTR(J1) == 0 .AND. NBS(J1) >= XSTL_NN) THEN
                NCLSTRS = NCLSTRS + 1
                CALL DFS(J1, NCLSTRS, XSTL_ADJ, CLSTR, X_CLSTR)
                CLSTRSIZE(NCLSTRS) = X_CLSTR
            ENDIF
        ENDDO

    !   Identify the size and ID of the largest cluster.
        LRGST_X_CLSTR = 0
        DO J1 = 1, NCLSTRS
            IF(CLSTRSIZE(J1) > LRGST_X_CLSTR) THEN
                LRGST_X_ID = J1
                LRGST_X_CLSTR = CLSTRSIZE(J1)
            ENDIF
        ENDDO

        IF(OUTPUTT) THEN

            OPEN (UNIT = 173, FILE ='pos_x_clstr.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            OPEN (UNIT = 174, FILE ='lrgst_x_clstr.xyz', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            OPEN (UNIT = 175, FILE ='nxtal.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')

            WRITE (174,*) NPART
            WRITE (174,*) 
            
            DO J1 = 1, NPART
                IF(CLSTR(J1) == LRGST_X_ID) THEN
                    WRITE (173,*) R(:,J1)
                    WRITE (174,*) "O", R(:,J1)
                ELSEIF(CLSTR(J1) /= 0) THEN
                    WRITE (174,*) "N", R(:,J1)
                ELSE
                    WRITE (174,*) "C", R(:,J1)
                ENDIF
            ENDDO

            WRITE (175,*) SUM(CLSTRSIZE)

            CLOSE (UNIT = 173, STATUS = 'KEEP')
            CLOSE (UNIT = 174, STATUS = 'KEEP')
            CLOSE (UNIT = 175, STATUS = 'KEEP')
        ENDIF
        
    END SUBROUTINE LRGSTXCLSTR

    SUBROUTINE GET_NEIGHBOURS(COORDS, DIST, NB)

        USE COMMONS, ONLY: R, BOX, NNEIGHT, NNEIGH

        IMPLICIT NONE

        INTEGER                     :: J1, J2, J3
        INTEGER                     :: CI(3), J_LIST(NPART), JJ
        
        REAL(KIND=DP)               :: DIJ, RI(NDIM), RJ(NDIM), RIJ(NDIM), RJI(NDIM), RIJSQ
        REAL(KIND=DP)               :: PHI_IJ, PHI_JI, CTHETA_IJ, CTHETA_JI
        
        INTEGER                     :: MAX_INDEX, REPLACE_INDEX, JR
        REAL(KIND=DP)               :: MINNEIGHDIST, NNEIGHS(2,NNEIGH,NPART)
        LOGICAL                     :: ADDNEIGHI, REPLACENEIGHI, ADDNEIGHJ, REPLACENEIGHJ

        REAL(KIND=DP), INTENT(OUT)  :: DIST(NPART,NPART), COORDS(2,NPART,NPART)
        INTEGER, INTENT(OUT)        :: NB(NPART)

        COORDS      = 0.0_dp
        DIST(:,:)   = 1.e+10
        NB(:)       = 0
        NNEIGHS     = 1.e+10

        DO J1 = 1, NPART
        !   Position of particle I
            RI  = R(:,J1)

            IF(CELLLISTT) THEN
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
    
                IF(CELLLISTT) THEN
                    J2 = J_LIST(JJ)             ! Get neighbour index
                    IF ( J2 == 0 ) EXIT         ! List exhausted
                ELSE
                    J2 = JJ
                    IF ( J2 == NPART+1 ) EXIT     ! List exhausted
                    IF ( J1 == J2 ) CYCLE
                ENDIF
    
            !   Position of particle J
                RJ    = R(:,J2)
                RIJ   = RI - RJ
                RIJ   = RIJ - BOX*ANINT( RIJ/BOX )
                RIJSQ = DOT_PRODUCT(RIJ,RIJ)

                IF(RIJSQ <= BOP_CUT_SQ) THEN
                    ADDNEIGHI = .FALSE.
                    ADDNEIGHJ = .FALSE.

                !   Compute the reverse separation vector for particle J 
                    RJI         = -RIJ
                    DIJ         = SQRT(RIJSQ)

                !   Particle I 
                    IF(NNEIGHT .AND. NB(J1) >= NNEIGH) THEN
                        IF(NB(J1) > NNEIGH) THEN
                            STOP "No. of nearest neighbours exceeds limit"
                        ENDIF

                        MAX_INDEX = MAXLOC(NNEIGHS(2,:,J1), DIM=1)
                        
                        IF(NNEIGHS(2,MAX_INDEX,J1) > DIJ) THEN
                            ADDNEIGHI     = .TRUE.
                            NB(J1)        = NB(J1) - 1

                            REPLACE_INDEX = INT(NNEIGHS(1,MAX_INDEX,J1))
                            DIST(REPLACE_INDEX,J1) = 1.e+10
                        ENDIF

                    ELSE
                        ADDNEIGHI = .TRUE.
                        MAX_INDEX = NB(J1) + 1
                    ENDIF

                    IF(ADDNEIGHI .EQV. .TRUE.) THEN
                        DIST(J2,J1) = DIJ
                    !   Update the number of neighbours for particles I and J                    
                        NB(J1)      = NB(J1) + 1
                    !   Check for artificial negative coordiates which may alter the sign of our
                    !   complex numbers later down the line.
                        DO J3 = 1, NDIM
                            IF(RIJ(J3) == -0.0_dp) RIJ(J3) = 0.0_dp
                        ENDDO
                    !   Calculate the cosine of the polar angle for RIJ and RJI
                        CTHETA_IJ = RIJ(3) / DIJ
                    !   Calculate the azimuthal angle for RIJ and RJI, making sure to add PI so 
                    !   that 0 <= phi <= 2pi (as DATAN2 returns a number between -pi and pi).
                        PHI_IJ = ATAN2( RIJ(2), RIJ(1) ) + PI
                    !   Save the spherical coordinates associated with the bond between particles I and J
                        COORDS(1,J2,J1) = CTHETA_IJ
                        COORDS(2,J2,J1) = PHI_IJ

                        IF(NNEIGHT) THEN
                            NNEIGHS(1,MAX_INDEX,J1) = J2
                            NNEIGHS(2,MAX_INDEX,J1) = DIJ
                        ENDIF
                    ENDIF
                    
                !   Particle J
                    IF(NNEIGHT .AND. NB(J2) >= NNEIGH) THEN
                        IF(NB(J2) > NNEIGH) THEN
                            STOP "No. of nearest neighbours exceeds limit"
                        ENDIF

                        MAX_INDEX = MAXLOC(NNEIGHS(2,:,J2), DIM=1)
                        
                        IF(NNEIGHS(2,MAX_INDEX,J2) > DIJ) THEN
                            ADDNEIGHJ     = .TRUE.
                            NB(J2)        = NB(J2) - 1

                            REPLACE_INDEX = INT(NNEIGHS(1,MAX_INDEX,J2))
                            DIST(REPLACE_INDEX,J2) = 1.e+10
                        ENDIF

                    ELSE
                        ADDNEIGHJ = .TRUE.
                        MAX_INDEX = NB(J2) + 1
                    ENDIF

                    IF(ADDNEIGHJ .EQV. .TRUE.) THEN
                        DIST(J1,J2) = DIJ
                    !   Update the number of neighbours for particles I and J    
                        NB(J2)      = NB(J2) + 1

                    !   Check for artificial negative coordiates which may alter the sign of our
                    !   complex numbers later down the line.
                        DO J3 = 1, NDIM
                            IF(RJI(J3) == -0.0_dp) RJI(J3) = 0.0_dp
                        ENDDO
                    !   Calculate the cosine of the polar angle for RIJ and RJI
                        CTHETA_JI = RJI(3) / DIJ
                    !   Calculate the azimuthal angle for RIJ and RJI, making sure to add PI so 
                    !   that 0 <= phi <= 2pi (as DATAN2 returns a number between -pi and pi).
                        PHI_JI = ATAN2( RJI(2), RJI(1) ) + PI
                    !   Save the spherical coordinates associated with the bond between particles I and J
                        COORDS(1,J1,J2) = CTHETA_JI
                        COORDS(2,J1,J2) = PHI_JI

                        IF(NNEIGHT) THEN
                            NNEIGHS(1,MAX_INDEX,J2) = J1
                            NNEIGHS(2,MAX_INDEX,J2) = DIJ
                        ENDIF
                    ENDIF

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
        COMPLEX(KIND=DP), INTENT(IN)    :: QLM(-L:L, NPART)

        INTEGER                         :: J1, M
        REAL(KIND=DP)                   :: QL(NPART)

        QL = 0.0_dp

        DO J1 = 1, NPART
                
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

        INTEGER, INTENT(IN)             :: L, NB(NPART)
        REAL(KIND=DP), INTENT(IN)       :: DIST(NPART,NPART)
        COMPLEX(KIND=DP), INTENT(IN)    :: QLM(-L:L, NPART)

        INTEGER                         :: J1, J2, M
        REAL(KIND=DP)                   :: DIJ
        COMPLEX(KIND=DP)                :: QBARLM(-L:L, NPART)

        REAL(KIND=DP)                   :: QBARL(NPART)

        QBARLM = 0.0_dp
        QBARL  = 0.0_dp

        DO J1 = 1, NPART
            
            DO J2 = 1, NPART
                
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
        COMPLEX(KIND=DP), INTENT(IN)    :: QLM(-L:L, NPART)

        INTEGER                         :: J1, M
        COMPLEX(KIND=DP)                :: GQLM(-L:L)
        REAL(KIND=DP)                   :: GQL

        GQLM = 0.0_dp
        GQL  = 0.0_dp

        DO J1 = 1, NPART
            DO M = -L, L
                GQLM(M) = GQLM(M) + QLM(M,J1)
            ENDDO
        ENDDO

        GQLM = GQLM / REAL(NPART,DP)

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
        COMPLEX(KIND=DP), INTENT(IN)    :: QLM(-L:L, NPART)
        REAL(KIND=DP), INTENT(IN)       :: DIST(NPART,NPART)

        INTEGER                         :: J1, J2
        REAL(KIND=DP)                   :: QLI_QLJ, QL_I, QL_J, DLIJ
        REAL(KIND=DP)                   :: DIJ

        REAL(KIND=DP)                   :: QLDOTQL(NPART,NPART)

        QLDOTQL = 100.0_dp

        DO J1 = 1, NPART-1

            QL_I = SQRT( SUM( ABS(QLM(:,J1))**2 ) )

            DO J2 = J1+1, NPART
                
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

        INTEGER, INTENT(IN)         :: L, NB(NPART)
        REAL(KIND=DP), INTENT(IN)   :: COORDS(2,NPART,NPART), DIST(NPART,NPART)

        INTEGER                     :: J1, J2
        INTEGER                     :: M
        REAL(KIND=DP)               :: DIJ, LM
        REAL(KIND=DP)               :: CTHETA_IJ, CTHETA_JI, PHI_IJ, PHI_JI, PLMI, PLMJ
        COMPLEX(KIND=DP)            :: YLM, YLMI, YLMJ, YLMMI, YLMMJ
        COMPLEX(KIND=DP), PARAMETER :: i = (0.0_dp, 1.0_dp)

        COMPLEX(KIND=DP)            :: QLM(-L:L, NPART)

        QLM = 0.0_dp

        DO J1 = 1, NPART-1
            DO J2 = J1+1, NPART
            
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
        
        DO J1 = 1, NPART
            IF(NB(J1) > 0) THEN
                QLM(:,J1) = QLM(:,J1) / REAL(NB(J1),DP)
            ELSE
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
        INTEGER :: CLSTR(NPART)
        LOGICAL :: ADJ(NPART,NPART)
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
                DO J2 = 1, NPART
                    IF ( ADJ(U,J2) ) THEN
                        CALL PUSH(S,J2)
                    ENDIF
                ENDDO
            ENDIF

        ENDDO

    END SUBROUTINE DFS

END MODULE ORDERPARAM


! GET TETRASTACK ENVIRONMENTS
! QLQL4 = QL_DOT_QL(4)
! QLQL6 = QL_DOT_QL(6)
! PRINT *, NPART
! PRINT *, ""
! DO J1 = 1, NPART
!     NT = 0
!     NC = 0
!     NH = 0
!     DO J2 = 1, NPART 

!         ! IF(QLQL4(J1,J2) /= 0.0_dp .AND. QLQL6(J1,J2) /= 0.0_dp) PRINT *, QLQL4(J1,J2), QLQL6(J1,J2)

!         IF(J1 == J2) CYCLE

!         IF(QLQL6(J1,J2) > 0.25 .AND. QLQL6(J1,J2) < 0.6) THEN
!             NT = NT + 1
!             IF(QLQL4(J1,J2) > -0.2 .AND. QLQL4(J1,J2) < 0.35) THEN
!                 NC = NC + 1
!             ELSEIF(QLQL4(J1,J2) > -0.65 .AND. QLQL4(J1,J2) < 0.35) THEN
!                 NH = NH + 1
!             ENDIF
!         ENDIF
        
!     ENDDO ! Loop over particles j

!     IF(NT > 5) THEN
!         IF(NC > 5) THEN
!             PRINT *, "C", R(:,J1)
!         ELSEIF(NH > 1) THEN
!             PRINT *, "O", R(:,J1)
!         ELSE
!             PRINT *, "N", R(:,J1)
!         ENDIF
!     ELSE 
!         PRINT *, "H", R(:,J1)
!     ENDIF

! ENDDO ! Loop over particles i