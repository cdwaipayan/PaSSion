SUBROUTINE KIHARA(ENERGY, RIJ, J1, J2)
    !     ==============================================================================================
    !     This subroutine is for particles interacting via a continuous Kern-Frenkel patchy model. 
    !     See "Designing a Bernal Spiral from Patchy Colloids", ACS Nano, 7, 1246-1256 (2013) for the
    !     patch-patch interactions.
    !     ==============================================================================================
        USE COMMONS, ONLY: DP, NDIM, RBSITES, VIRTEMP
        USE COMMONS, ONLY: GLJN, DCHECKSQ
    
        IMPLICIT NONE
    
        INTEGER, INTENT(IN)       :: J1, J2
        REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM)
        
    !   Parameters related to the translational coordinates of the particles         
        REAL(KIND=DP) :: SDISTSQ
    !   Parameters specific to potential energy function 
        REAL(KIND=DP) :: DIJ(NDIM), SD2, SDN, S2DN, VLJ, EK
    !   Parameters related to calculating the virial pressure
        REAL(KIND=DP) :: WIJ
    
    !   PARAMETER TO BE OUTPUT FROM SUBROUTINE
        REAL(KIND=DP), INTENT(OUT) :: ENERGY
    
        ENERGY   = 0.0_dp

        CALL SHORTD( RIJ, RBSITES(:,1,J1), RBSITES(:,1,J2), DIJ )
        SDISTSQ = DOT_PRODUCT(DIJ,DIJ)
    !   Calculate repulsive Kihara contribution to the pair energy
    !   This depends on the shortest distance between the spherocylinders
        IF(SDISTSQ < DCHECKSQ) THEN
            SD2     = 1.0_dp/SDISTSQ
            SDN     = SD2**(GLJN/2.0_dp)
            S2DN    = SDN**2
            VLJ     = S2DN - SDN   
            WIJ     = VLJ + S2DN
            EK      = 4.0_dp*(VLJ + 0.25_dp)
            ENERGY  = ENERGY + EK 

            IF(ISNAN(ENERGY)) RETURN
            VIRTEMP = VIRTEMP + 4.0_dp*GLJN*WIJ
        ENDIF
        
    END SUBROUTINE KIHARA
    
    !====================================================================================================
        
    SUBROUTINE DEF_KIHARA()
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
    
        USE COMMONS, ONLY: DP, RCUT, RCUTSQ, DCHECKSQ, GLJN, RLNGTH, HLFLNGTH
        USE COMMONS, ONLY: NDIM, NPART, NSITES, REFSITE, RBSITES

        IMPLICIT NONE
    
        NSITES = 1

        ALLOCATE( REFSITE(NDIM,NSITES), RBSITES(NDIM,NSITES,NPART) )
    
    !   Elongated triblock patchy (ETP) particles
        REFSITE(:,1) = [0.0_dp, 0.0_dp, 1.0_dp]

    !   The cut-off for the repulsive Kihara component of the potential 
        DCHECKSQ = 2.0_dp**(2.0_dp/GLJN)
    !   Adjust the cut-off from the input to account for the fact that the particles
    !   have an elongated core.
        RCUT     = RCUT + RLNGTH
        RCUTSQ   = RCUT * RCUT
    !   The half-length of the ETP particles 
        HLFLNGTH = RLNGTH / 2.0_dp
        
    END SUBROUTINE
    
    !====================================================================================================
    
    SUBROUTINE VIEW_KIHARA()
    
        USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, BOX, VIEWUNIT, HLFLNGTH
        USE ROTATIONS_MODULE, ONLY: Q_TO_RM
    
        IMPLICIT NONE
    
        INTEGER:: J1
        REAL(KIND=DP) :: RWRITE(3,NPART), RM(3,3)
    
        WRITE(VIEWUNIT,*) NPART*5
        WRITE(VIEWUNIT,*)
    
        RWRITE(:,:) = R(:,:)
        
        DO J1 = 1, NPART
            RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX(:))*BOX(:)
        END DO
    
        DO J1 = 1, NPART
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(1,J1), RWRITE(2,J1), RWRITE(3,J1)
            RM = Q_TO_RM( Q(:,J1) )
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(:,J1) + HLFLNGTH*MATMUL(RM,REFSITE(:,1))
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(:,J1) - HLFLNGTH*MATMUL(RM,REFSITE(:,1))
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(:,J1) + 0.5_dp*HLFLNGTH*MATMUL(RM,REFSITE(:,1))
            WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Au ', RWRITE(:,J1) - 0.5_dp*HLFLNGTH*MATMUL(RM,REFSITE(:,1))
        END DO
        
    END SUBROUTINE
    
    SUBROUTINE SHORTD(RIJ, UI, UJ, DIJ)
!====================================================================================================
!   Subroutine to evaluate the shortest distance between 2.0_dp spherocylinders
!   RIJ = Vector connecting the geometrical centers of the 2.0_dp spherocylinders
!   UI  = Unitary vector definig the orientation of particle I
!   UJ  = Unitary vector definig the orientation of particle J
!   DIJ = Vector giving the shortest distance between the 2.0_dp spherocylinders   
!====================================================================================================
        USE COMMONS, ONLY: DP, HLFLNGTH
    
        IMPLICIT NONE
    
        REAL(KIND=DP), INTENT(IN)   :: RIJ(3), UI(3), UJ(3)
    
        REAL(KIND=DP)               :: LAMI, LAMJ, LI, LJ
        REAL(KIND=DP)               :: RDOT, URI, URJ, UIJ, CC
    
        REAL(KIND=DP), INTENT(OUT)  :: DIJ(3)
    
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
                LAMI = SIGN(HLFLNGTH,URI)
                LAMJ = LAMI*UIJ - URJ
                IF(ABS(LAMJ) > HLFLNGTH) LAMJ = SIGN(HLFLNGTH,LAMJ) 
            ELSE
                LAMI = 0.0_dp
                LAMJ = 0.0_dp
            ENDIF
        ELSE
            LAMI = (URI - UIJ*URJ) / CC
            LAMJ = (UIJ*URI - URJ) / CC
    
            LI = ABS(LAMI) - HLFLNGTH
            LJ = ABS(LAMJ) - HLFLNGTH
    
            IF( (LI > 0.0_dp) .OR. (LJ > 0.0_dp) ) THEN 
                IF( LI > LJ ) THEN 
                    LAMI = SIGN(HLFLNGTH, LAMI)
                    LAMJ = LAMI*UIJ - URJ
                    IF(ABS(LAMJ) > HLFLNGTH) LAMJ = SIGN(HLFLNGTH,LAMJ) 
                ELSE
                    LAMJ = SIGN(HLFLNGTH, LAMJ)
                    LAMI = LAMJ*UIJ + URI
                    IF(ABS(LAMI) > HLFLNGTH) LAMI = SIGN(HLFLNGTH,LAMI) 
                ENDIF
            ENDIF
        ENDIF
    
        DIJ = RIJ + LAMJ*UJ - LAMI*UI
    
    END SUBROUTINE SHORTD