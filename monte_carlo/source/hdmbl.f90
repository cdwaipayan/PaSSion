SUBROUTINE HDMBL(ENERGY, J1, J2)

    USE COMMONS, ONLY: DP, R, NDIM, NSITES, RBSITES, BOX, OVERLAPT

    IMPLICIT NONE

    INTEGER, INTENT(IN)       :: J1, J2

    INTEGER       :: J3, J4
    REAL(KIND=DP) :: RI(NDIM), RJ(NDIM), RABSQ
    REAL(KIND=DP) :: RA(NDIM), RB(NDIM), RAB(NDIM)

    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY   = 0.0_dp

    RI = R(:,J1)
    RJ = R(:,J2)

    DO J3 = 1, NSITES
    !   Position of site a on particle I
        RA  = RI + RBSITES(:,J3,J1)
        
        DO J4 = 1, NSITES
        !   Position of site b on particle J
            RB    = RJ + RBSITES(:,J4,J2)

            RAB   = RA - RB
            RAB   = RAB - BOX*ANINT( RAB/BOX )
            RABSQ = DOT_PRODUCT(RAB,RAB)

            IF (RABSQ <= 1.0_dp) THEN
                OVERLAPT = .TRUE.
                RETURN
            ENDIF
        ENDDO ! Loop over each of the patches on particle j
    ENDDO ! Loop over particles j

END SUBROUTINE HDMBL

!====================================================================================================
    
SUBROUTINE DEF_HDMBL()
!----------------------------------------------------------------
! 
!----------------------------------------------------------------
    
    USE COMMONS, ONLY: DP, PI, REFSITE, RCUT, RCUTSQ, LSTAR

    IMPLICIT NONE

    REFSITE(:,1)= [ 0.0_dp, 0.0_dp, -LSTAR/2.0_dp]
    REFSITE(:,2)= [ 0.0_dp, 0.0_dp,  LSTAR/2.0_dp]

    RCUT   = 1.0_dp + LSTAR
    RCUTSQ = RCUT*RCUT
    
END SUBROUTINE
    
!====================================================================================================
    
SUBROUTINE VIEW_HDMBL()

    USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT
    USE ROTATIONS_MODULE, ONLY: Q_TO_RM

    IMPLICIT NONE

    INTEGER:: J1, J2
    REAL(KIND=DP) :: RWRITE(3,NPART), RM(3,3), RBCOORDS(3)

    WRITE(VIEWUNIT,*) NPART*NSITES
    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO J1 = 1, NPART
        RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX(:))*BOX(:)
    END DO

    DO J1 = 1, NPART
        RM = Q_TO_RM( Q(:,J1) )
        DO J2 = 1, NSITES
            RBCOORDS = RWRITE(:,J1) + MATMUL(RM,REFSITE(:,J2))
            IF(J2==1)THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ELSEIF(J2==2) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ENDIF
        ENDDO
    END DO
    
END SUBROUTINE