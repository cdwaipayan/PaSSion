SUBROUTINE GLJ(ENERGY, RIJSQ)

!     This subroutine is for spherical particles interacting via generalised Lennard-Jones pair potential.
!     Reference:
!     G. A. Vliegenthart, J. F. M. Lodge and H. N. W. Lekkerkerker, Physica A 263, 378-388 (1999).

    USE COMMONS, ONLY: DP, GLJN, VIRTEMP, TAILCORT, TAILV, TAIL2V, WCAT

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN)  :: RIJSQ
    REAL(KIND=DP)              :: R2, RLJN, R2LJN, VIJ, WIJ
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    R2      = 1.0_dp/RIJSQ
    RLJN    = R2**(GLJN/2.0_dp)
    R2LJN   = RLJN*RLJN
    VIJ     = R2LJN - RLJN
    
    IF (TAILCORT) THEN
        ENERGY = 4.0_dp*(VIJ - TAILV)
        WIJ    = VIJ + R2LJN - TAIL2V
    ELSE
        IF(WCAT) THEN
            ENERGY  = 4.0_dp*(VIJ + 0.25_dp)
        ELSE
            ENERGY  = 4.0_dp*VIJ
        ENDIF
        WIJ     = VIJ + R2LJN - TAIL2V
    ENDIF

!   Pair virial function w(r) = r dv(r)/dr
    VIRTEMP = VIRTEMP + 4.0_dp*GLJN*WIJ

END SUBROUTINE GLJ

SUBROUTINE YUKAWA(ENERGY, RIJSQ)

    USE COMMONS, ONLY: DP, VIRTEMP, YUKKAP

    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN)  :: RIJSQ       
    REAL(KIND=DP)              :: DIST, UYUK
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY   = 0.0_dp

    DIST = SQRT(RIJSQ)
!   Calculate repulsive contribution to the pair energy
    UYUK    = EXP( -YUKKAP*(DIST - 1.0_dp) )
    ENERGY  = ENERGY + UYUK / DIST
    VIRTEMP = VIRTEMP + UYUK*(YUKKAP + 1.0_dp/DIST)

END SUBROUTINE YUKAWA

SUBROUTINE DEF_GLJ()

    USE COMMONS, ONLY: DP, RCUTSQ, GLJN, TAILCORT, TAILV, TAIL2V

    IMPLICIT NONE

    REAL(KIND=DP) RLJN, R2LJN

    IF(TAILCORT) THEN
        RLJN   = (1.0_dp/RCUTSQ)**(GLJN/2.0_dp)
        R2LJN  = RLJN*RLJN
        TAILV  = (R2LJN - RLJN)
        TAIL2V = (TAILV + R2LJN)
    ENDIF
    
END SUBROUTINE

!     ==============================================================================================

REAL(KIND=DP) FUNCTION CORP(RC, RHO)

    USE COMMONS, ONLY: DP, PI

    IMPLICIT NONE

    REAL(KIND=DP) :: RC, RHO, RI3

    RI3  = 1.0_dp/(RC*RC*RC)
    CORP = 16.0_dp*PI*(RHO**2)*(2.0_dp*RI3*RI3*RI3/9.0_dp - RI3/3.0_dp)

    RETURN

END FUNCTION

!     ==============================================================================================

REAL(KIND=DP) FUNCTION CORU(RC, RHO)

    USE COMMONS, ONLY: DP, PI

    IMPLICIT NONE

    REAL(KIND=DP) :: RC, RHO, RI3

    RI3  = 1.0_dp/(RC*RC*RC)
    CORU = 8.0_dp*PI*RHO*(RI3*RI3*RI3/9.0_dp - RI3/3.0_dp)

    RETURN
END FUNCTION

      !====================================================================================================

SUBROUTINE VIEW_GLJ()

    USE COMMONS, ONLY: DP, NPART, R, BOX, VIEWUNIT

    IMPLICIT NONE

    INTEGER:: I
    REAL(KIND=DP) :: RWRITE(3,NPART)

    WRITE(VIEWUNIT,*) NPART
    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO I = 1, NPART
        RWRITE(:,I) = RWRITE(:,I) - ANINT(RWRITE(:,I)/BOX(:))*BOX(:)
    END DO

    DO I = 1, NPART
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RWRITE(1,I), RWRITE(2,I), RWRITE(3,I)
    END DO
    
END SUBROUTINE