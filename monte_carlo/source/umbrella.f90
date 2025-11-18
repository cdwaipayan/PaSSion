!   ==============================================================================================   
    SUBROUTINE UMBRELLA_BIAS(ISTEP)
!   ==================================================================
!   Calculate the umbrella bias potetnial
!   ==================================================================

        USE COMMONS, ONLY: DP, NDIM, NPART, R, Q, RBSITES, PE, RIGIDT, BETAKB, PRES, RHO, VIR, VLM, NSITES
        USE COMMONS, ONLY: BOX, CELLLISTT, TRGTSEED, CLSTRSIZE, NPINK, SCALED_R, RCUT
        USE COMMONS, ONLY: SUMLRGSTXCLSTRBLK, NUCCOUNT, TOTNUCCOUNT, NUCSEEDT, NUCCNTT, CRITNC, NSEEDMAX
        USE COMMONS, ONLY: GQL, TRGTQ, QPINK
        USE CELL_LIST
        USE ORDERPARAM, ONLY: GET_BOPS, LRGSTXCLSTR

        IMPLICIT NONE

        INTEGER, INTENT(IN)                 :: ISTEP

        INTEGER                             :: J1, NXTLNEW
        REAL(KIND=DP)                       :: BIASNEW, DELB, QLNEW, N1, N2
        REAL(KIND=DP)                       :: DRAND48

        LOGICAL                             :: REJECTT

        INTEGER, SAVE                       :: NXTLOLD
        REAL(KIND=DP), SAVE                 :: ENEOLD, VLMOLD, VIROLD, PRESOLD, RHOOLD, BIASOLD
        REAL(KIND=DP), ALLOCATABLE, SAVE    :: ROLD(:,:), BOXO(:), QOLD(:,:), RBOLD(:,:,:), SUMBOLD(:)

        IF(.NOT. ALLOCATED(BOXO)) ALLOCATE(BOXO(NDIM))
        IF(.NOT. ALLOCATED(ROLD)) ALLOCATE(ROLD(NDIM,NPART))
        IF(.NOT. ALLOCATED(QOLD) .AND. RIGIDT) ALLOCATE(QOLD(4,NPART))
        IF(.NOT. ALLOCATED(RBOLD) .AND. RIGIDT) ALLOCATE(RBOLD(NDIM,NSITES,NPART))
        IF(.NOT. ALLOCATED(SUMBOLD)) ALLOCATE(SUMBOLD(NDIM))

        REJECTT = .FALSE.

!   ========================================================================================================
!   Compute and save system observables at the beginning of the simulation. This is done in initialise.f90. 
!   ========================================================================================================
        IF(ISTEP == -1) THEN
            
            ENEOLD  = PE
            VLMOLD  = VLM
            RHOOLD  = RHO
            PRESOLD = PRES
            VIROLD  = VIR
            BOXO    = BOX
            
            ROLD   = R
            IF(RIGIDT) THEN
                QOLD  = Q
                RBOLD = RBSITES
            ENDIF
            
            CALL GET_BOPS()

            IF(NUCSEEDT) THEN
                CALL LRGSTXCLSTR(.TRUE., NXTLOLD)
                CLSTRSIZE = NXTLOLD
                IF(NUCCNTT) THEN
                    N1 = REAL(NXTLOLD)**(2.0_dp/3.0_dp)
                    N2 = REAL(NXTLOLD)**(1.0_dp/3.0_dp)
                    BIASOLD = NPINK*N1*(N2-CRITNC)
                ELSE
                    BIASOLD = NPINK * REAL((NXTLOLD - TRGTSEED),DP)**2
                ENDIF
            ELSE
                BIASOLD = QPINK *(GQL(1) - TRGTQ)**2
            ENDIF
            
            RETURN
        ENDIF

!   ====================================================================================
!       Calculate the size of the largest crystalline cluster
!   ====================================================================================
        CALL GET_BOPS()
        IF(NUCSEEDT) THEN
            CALL LRGSTXCLSTR(.FALSE., NXTLNEW)
            IF(NUCCNTT) THEN
                IF(NXTLNEW>NSEEDMAX) THEN
                    REJECTT = .TRUE.
                ELSE
                    N1 = REAL(NXTLNEW)**(2.0_dp/3.0_dp)
                    N2 = REAL(NXTLNEW)**(1.0_dp/3.0_dp)
                    BIASNEW = NPINK*N1*(N2-CRITNC)
                ENDIF
            ELSE
            !   Calculate the new value for the harmonic bias potential, which we take to be:
            !   phi = k*(NXT_max - NXT_0)**2.
                BIASNEW = NPINK * REAL((NXTLNEW - TRGTSEED),DP)**2
            ENDIF
        ELSE
            QLNEW   = GQL(1)
            BIASNEW = QPINK *(QLNEW - TRGTQ)**2
        ENDIF
    
    !   Calculate the change in the bias potential
        IF (.NOT. REJECTT) DELB = BIASNEW - BIASOLD

!   ====================================================================================
!       Metropolis acceptence criteria
!   ====================================================================================
        IF ( REJECTT .OR. (DRAND48()>EXP(-BETAKB*DELB).AND.(NPINK>0.0_dp.OR.QPINK>0.0_dp)) ) THEN
        !   Reject the move, and reset system parameters.
            PE   = ENEOLD
            VLM  = VLMOLD
            RHO  = RHOOLD
            PRES = PRESOLD
            VIR  = VIROLD
            BOX  = BOXO

            R = ROLD
            IF(RIGIDT) THEN
                Q       = QOLD
                RBSITES = RBOLD
            ENDIF
        !   Update cell-list
            IF(CELLLISTT) THEN
                CALL FINALIZE_LIST()
                CALL INITIALIZE_LIST( NPART, RCUT/BOX )
                SCALED_R = 0.0_dp
                DO J1 = 1, NPART
                    SCALED_R(:,J1) = R(:,J1)/BOX
                ENDDO
                CALL MAKE_LIST( NPART, SCALED_R )
            ENDIF

        ELSE
        !   Accept the move and save current system parameters
            ENEOLD  = PE
            VLMOLD  = VLM
            RHOOLD  = RHO
            PRESOLD = PRES
            VIROLD  = VIR
            BOXO    = BOX

            ROLD = R
            IF(RIGIDT) THEN
                QOLD  = Q
                RBOLD = RBSITES
            ENDIF

            BIASOLD   = BIASNEW

            IF(NUCSEEDT) THEN
                NXTLOLD = NXTLNEW
            ENDIF

        ENDIF

        IF(NUCSEEDT) THEN
            CLSTRSIZE              = NXTLOLD
            SUMLRGSTXCLSTRBLK      = SUMLRGSTXCLSTRBLK + CLSTRSIZE
            NUCCOUNT(CLSTRSIZE)    = NUCCOUNT(CLSTRSIZE) + 1
            TOTNUCCOUNT(CLSTRSIZE) = TOTNUCCOUNT(CLSTRSIZE) + 1
        ENDIF

    END SUBROUTINE

!   ==============================================================================================
!    Auxillary routines related to Monte Carlo move sets. The majority of these routines are
!    related to some kind of biased sampling (i.e., Umbrella Sampling or Interface Pinning).
!   ==============================================================================================

    SUBROUTINE HARMONIC_TRAP(BIAS)
    !   ----------------------------------------------------------------------------------------------------------------
    !    Apply a harmonic trap to specific particles, and calculate the associated energy cost of displacing a particle
    !    from its ideal lattice site.
    !   ----------------------------------------------------------------------------------------------------------------
        
        USE COMMONS, ONLY: DP, NDIM, R, RBSITES, RIGIDT, NSITES, SEEDSIZE, EQUISEEDT, INDXP
        USE COMMONS, ONLY: HALFST, FORCEC, FORCEO, RREF, RBREF, ZPIN

        IMPLICIT NONE

        REAL(KIND=DP)               :: RDIFF(NDIM), THETAA, THETAB
        REAL(KIND=DP), INTENT(OUT)  :: BIAS

        BIAS = 0.0_dp

    !   If a harmonic trap is being applied to the system, calculate the energy cost for moving particle I.
        IF( (HALFST .AND. RREF(3,INDXP) < ZPIN) .OR. (EQUISEEDT .AND. INDXP <= SEEDSIZE) ) THEN
            RDIFF = R(:,INDXP) - RREF(:,INDXP)
            BIAS = BIAS + FORCEC*DOT_PRODUCT(RDIFF,RDIFF)
            
            IF(RIGIDT) THEN
                THETAA = ACOS(MIN(MAX(DOT_PRODUCT(RBREF(:,1,INDXP),RBSITES(:,1,INDXP)), -1.0_dp), 1.0_dp))

                IF(NSITES==2) THEN
                    BIAS = BIAS + FORCEO*SIN(THETAA)**2
                ELSEIF(NSITES==3) THEN
                    THETAB = ACOS(MIN(MAX(DOT_PRODUCT(RBREF(:,2,INDXP),RBSITES(:,2,INDXP)), -1.0_dp), 1.0_dp))
                    BIAS = BIAS + FORCEO*( SIN(3.0_dp*THETAA/2.0_dp)**2 + SIN(THETAB)**2 )
                ELSEIF(NSITES==4) THEN
                    THETAB = ACOS(MIN(MAX(DOT_PRODUCT(RBREF(:,2,INDXP),RBSITES(:,2,INDXP)), -1.0_dp), 1.0_dp))
                    BIAS = BIAS + FORCEO*( SIN(THETAA)**2 + SIN(THETAB)**2 )
                ENDIF
            ENDIF
        ENDIF
    END SUBROUTINE

!   ==============================================================================================
!   ==============================================================================================
    
    SUBROUTINE COLLDENSFIELD(QR,Q,INDXP)
!   ------------------------------------------------------------------------------
!   Calculate the collective density field of the system or a specific particle
!   ------------------------------------------------------------------------------
        USE COMMONS, ONLY: DP, R, BOX, NPART, PI, NX, NY, NZ, NHLF
    
        IMPLICIT NONE
    
        INTEGER, INTENT(IN)             :: INDXP

        REAL(KIND=DP), INTENT(OUT)      :: QR
        COMPLEX(KIND=DP), INTENT(OUT)   :: Q
    
        INTEGER                         :: J1
        COMPLEX(KIND=DP), PARAMETER     :: i = (0.0_dp, 1.0_dp)
        REAL(KIND=DP)                   :: K(3)
    
        Q  = 0.0_dp
        QR = 0.0_dp
    
        K = (/ NX/BOX(1), NY/BOX(2), NZ/BOX(3) /)
    
        IF(INDXP == -1) THEN
            DO J1 = 1, NPART
                Q = Q + EXP(-i*DOT_PRODUCT(K,R(:,J1)))
            ENDDO
            QR = ABS(Q)/NHLF
        ELSE
            Q = EXP(-i*DOT_PRODUCT(K,R(:,INDXP)))
        ENDIF
    END SUBROUTINE


    SUBROUTINE INITIALISE_HISTOGRAM()

        USE COMMONS,   ONLY: DP, NBINS, BIN_MIN, BIN_MAX, BINEDGES, HISTCOUNTS
        IMPLICIT NONE

        INTEGER :: J1

        ALLOCATE(BINEDGES(0:NBINS), HISTCOUNTS(0:NBINS))

        BINEDGES   = 0.0_dp
        HISTCOUNTS = 0

        DO J1 = 0, NBINS
            BINEDGES(J1) = BIN_MIN + REAL(J1,DP) * (BIN_MAX-BIN_MIN)/REAL(NBINS,DP)
        ENDDO

    END SUBROUTINE

    SUBROUTINE UPDATE_HISTOGRAM(X)

        USE COMMONS,   ONLY: DP, NBINS, BIN_MIN, BIN_MAX, HISTCOUNTS

        IMPLICIT NONE

        REAL(KIND=DP), INTENT(IN) :: X
        INTEGER                   :: H_LOC

        H_LOC = FLOOR( (X-BIN_MIN) / (BIN_MAX-BIN_MIN) * REAL(NBINS,DP) )

        IF(H_LOC>NBINS) THEN
            H_LOC = NBINS
        ELSEIF(H_LOC<0) THEN
            H_LOC = 0
        ENDIF

        HISTCOUNTS(H_LOC) = HISTCOUNTS(H_LOC) + 1

    END SUBROUTINE