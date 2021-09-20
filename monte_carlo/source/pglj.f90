SUBROUTINE PGLJ(ENERGY, RIJ, RIJSQ, J1, J2)

!     This subroutine is for spherical particles interacting via generalised Lennard-Jones pair potential.
!     Reference:
!     G. A. Vliegenthart, J. F. M. Lodge and H. N. W. Lekkerkerker, Physica A 263, 378-388 (1999).

    USE COMMONS, ONLY: DP, NDIM, NSITES, RBSITES
    USE COMMONS, ONLY: GLJN, TSIGPW2, TAILV, TAILCORT

    IMPLICIT NONE

    INTEGER, INTENT(IN)       :: J1, J2
    REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM), RIJSQ

    INTEGER       :: J3
    REAL(KIND=DP) :: R2, RLJN, R2LJN, VIJ

    REAL(KIND=DP) :: RHAT(NDIM), THETA_IJK, THETA_JIL, TIJK_MIN, TJIL_MIN
    REAL(KIND=DP) :: EK(NDIM), EL(NDIM)

    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    VIJ      = 0.0_dp

    TIJK_MIN = 0.0_dp
    TJIL_MIN = 0.0_dp

    R2      = 1.0_dp/RIJSQ
    RLJN    = R2**(GLJN/2.0_dp)
    R2LJN   = RLJN*RLJN
    VIJ     = R2LJN - RLJN

!   Calculate the shifted potential
    IF (TAILCORT) VIJ = VIJ - TAILV
    
!   If the particles are not overlapping then calculate the angular modulation
!   due to patch-patch interactions
    IF(RIJSQ >= 1.0_dp .AND. TSIGPW2 > 0.0_dp) THEN
        RHAT     = RIJ / SQRT(RIJSQ)

    !   Find the patch on each particle which forms the narrowest angle with RHAT.
        DO J3 = 1, NSITES
        !   Direction of patch K on particle I
            EK  = RBSITES(:,J3,J1)
        !   Calculate the angle subtended by RHAT and EK, making sure to account for
        !   rounding errors in double precision floating point numbers 
            THETA_IJK = ACOS(MIN(MAX(DOT_PRODUCT(-RHAT,EK), -1.0_dp), 1.0_dp))
        !   Keep track of the smallest angle
            IF(J3 == 1 .OR. THETA_IJK < TIJK_MIN) TIJK_MIN = THETA_IJK
        
        !   Direction of patch L on particle J
            EL    = RBSITES(:,J3,J2)
        !   Calculate the angle subtended by RHAT and EL, making sure to account for
        !   rounding errors in double precision floating point numbers
            THETA_JIL = ACOS(MIN(MAX(DOT_PRODUCT(RHAT,EL), -1.0_dp), 1.0_dp))
        !   Keep track of the smallest angle
            IF(J3 == 1 .OR. THETA_JIL < TJIL_MIN) TJIL_MIN = THETA_JIL

        ENDDO
    !   Compute the product of the Gaussian functions
        VIJ = VIJ * EXP( -(TIJK_MIN**2 + TJIL_MIN**2) / TSIGPW2 )
    ENDIF

    ENERGY = 4.0_dp*VIJ

END SUBROUTINE PGLJ

SUBROUTINE DEF_PGLJ()
!----------------------------------------------------------------
! Define the position of the patches on particles interacting via
! the the patchy Lennard-Jones potential.
! We only consider a system of tetrahedral patchy particles here.
!----------------------------------------------------------------

    USE COMMONS, ONLY: DP, NSITES, REFSITE
    USE COMMONS, ONLY: SIGPW, TSIGPW2, RCUTSQ, GLJN, TAILCORT, TAILV, TAIL2V

    IMPLICIT NONE

    REAL(KIND=DP) RLJN, R2LJN

    TSIGPW2 = 2.0_dp*(SIGPW**2)

    IF(NSITES == 4) THEN
        REFSITE(:,1)= (/  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) /)
        REFSITE(:,2)= (/ -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp) /)
        REFSITE(:,3)= (/  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) /)
        REFSITE(:,4)= (/ -1.0_dp/SQRT(3.0_dp),  1.0_dp/SQRT(3.0_dp), -1.0_dp/SQRT(3.0_dp) /)
    ELSEIF(NSITES == 6) THEN
        REFSITE(:,1)= (/ 1.0_dp,  0.0_dp,  0.0_dp /)
        REFSITE(:,2)= (/ 0.0_dp,  1.0_dp,  0.0_dp /)
        REFSITE(:,3)= (/ 0.0_dp,  0.0_dp,  1.0_dp /)
        REFSITE(:,4)= (/-1.0_dp,  0.0_dp,  0.0_dp /)
        REFSITE(:,5)= (/ 0.0_dp, -1.0_dp,  0.0_dp /)
        REFSITE(:,6)= (/ 0.0_dp,  0.0_dp, -1.0_dp /)
    ELSE
        STOP "PGLJ MODEL MUST HAVE 4 OR 6 PATCHES"
    ENDIF

    !   *************************************************************************
    IF(TAILCORT) THEN
        RLJN   = (1.0_dp/RCUTSQ)**(GLJN/2.0_dp)
        R2LJN  = RLJN*RLJN
        TAILV  = (R2LJN - RLJN)
        TAIL2V = (TAILV + R2LJN)
    ENDIF
    
END SUBROUTINE

SUBROUTINE VIEW_PGLJ()

    USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, NSITES, BOX, VIEWUNIT
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
            RBCOORDS = RWRITE(:,J1) + 0.5*MATMUL(RM,REFSITE(:,J2))
            IF(J2==1)THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ELSEIF(J2==2) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ELSEIF(J2==3) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'K ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ELSEIF(J2==4) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Ni ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ELSEIF(J2==5) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Ga ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ELSE
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'Ge ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
            ENDIF
        ENDDO
    END DO
    
END SUBROUTINE