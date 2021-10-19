SUBROUTINE DIPLR_DSCS(ENERGY, RIJ, RIJSQ, J1, J2)

! This subroutine computes the potential energy for Dipolar Discotic Liquid Crystalline particles.

    USE COMMONS, ONLY: NDIM, DP, RBSITES, FORPI, TWOPI
    USE COMMONS, ONLY: GBCUT, GBK, GBV, GBM, GBX, GBXP, EWALDT, DPMUSQ, ALPHA, ALPSQ, INVRPI

    IMPLICIT NONE
    
    INTEGER, INTENT(IN)         :: J1, J2
    REAL(KIND=DP), INTENT(IN)   :: RIJ(NDIM), RIJSQ
    
    REAL(KIND=DP)              :: DIST, RHAT(NDIM), R2
    REAL(KIND=DP)              :: EI(NDIM), EJ(NDIM), ERI, ERJ, EIJ
    REAL(KIND=DP)              :: EGB, SG, EP, EP1, EP2, N1, N2, D1, D2
    REAL(KIND=DP)              :: RHOIJ, RHO6, RHO12
    REAL(KIND=DP)              :: A1, A2, B, C, EDD, THRER2
    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY = 0.0_dp
    EGB    = 0.0_dp
    EDD    = 0.0_dp

    DIST  = SQRT(RIJSQ)

    IF(DIST <= GBCUT .OR. EWALDT) THEN
        RHAT  = RIJ / DIST
        EI  = RBSITES(:,1,J1)
        EJ  = RBSITES(:,1,J2)
        EIJ = DOT_PRODUCT(EI, EJ)
        ERI = DOT_PRODUCT(RHAT, EI)
        ERJ = DOT_PRODUCT(RHAT, EJ)
    ENDIF

    IF(DIST <= GBCUT) THEN
! ------------------------------------------------------------------------------
!   Compute Gay-Berne contribution to the potential energy 
! ------------------------------------------------------------------------------
        D1    = (ERI+ERJ)**2
        D2    = (ERI-ERJ)**2
        N1    = GBX*EIJ
        N2    = GBXP*EIJ

        SG    = 1.0_dp / SQRT( (1.0_dp - GBX/2.0_dp * ( (D1/(1.0_dp+N1))+(D2/(1.0_dp-N1)) )) )
        RHOIJ = GBK / (DIST - SG + GBK)
        RHO6  = RHOIJ**6
        RHO12 = RHO6**2

        EP1   = 1.0_dp / SQRT( (1.0_dp - (GBX**2) * (EIJ**2)) )
        EP2   = 1.0_dp - GBXP/2.0_dp * ( (D1/(1.0_dp+N2))+(D2/(1.0_dp-N2)) )
        EP    = (EP1**GBV) * (EP2**GBM)

        EGB   = 4.0_dp * EP * (RHO12-RHO6)
    ENDIF

    IF(EWALDT) THEN
! ------------------------------------------------------------------------------
!   Compute Real-Space component of dipolar potential (with axial
!   dipole moment) contribution to the potential energy. 
! ------------------------------------------------------------------------------
        R2      = 1.0_dp/RIJSQ
        THRER2  = 3.0_dp*R2
        A1      = (2.0_dp*ALPHA*INVRPI* R2) * EXP(-ALPSQ*RIJSQ)
        A2      = R2*ERFC(ALPHA*DIST)/DIST
        B       = A2 + A1
        C       = THRER2*A2 + (2.0_dp*ALPSQ+THRER2)*A1
        EDD     = DPMUSQ*(EIJ*B - RIJSQ*ERI*ERJ*C)
    ENDIF

    ENERGY = EGB + EDD

END SUBROUTINE DIPLR_DSCS

! ========================================================================================================
! ========================================================================================================

SUBROUTINE DIPLR_DSCS_RECIP(J,PEK,TCOS,TSIN)

! This subroutine computes the potential energy for Dipolar Discotic Liquid Crystalline particles.

    USE COMMONS, ONLY: NPART, NDIM, DP, RBSITES, BOXL, PI, ALPHA, GUFCTR, NC, NCSQMAX, NVV, FCTR
      
    IMPLICIT NONE

    INTEGER, INTENT(IN)        :: J
    REAL(KIND=DP), INTENT(IN)  :: TCOS(NDIM,0:NC,NPART), TSIN(NDIM,0:NC,NPART)
    INTEGER                    :: NX, NY, NZ
    REAL(KIND=DP)              :: DUMMY, PC, PS, SUMC, SUMS, W
    REAL(KIND=DP)              :: VC(3), VS(3)
    REAL(KIND=DP), INTENT(OUT) :: PEK

    PEK = 0.0_dp

    W = PI*PI/(BOXL*ALPHA)**2

    DO NZ = 0, NC
        DO NY =-NC, NC
            DO NX =-NC, NC

                IF (NVV(NX,NY,NZ) == 0 .OR. NVV(NX,NY,NZ) > NCSQMAX) CYCLE
                
                SUMC = 0.0_dp; SUMS = 0.0_dp
                
                VC = (/TCOS(1,ABS(NX),J), TCOS(2,ABS(NY),J), TCOS(3,NZ,J)/)
                VS = (/TSIN(1,ABS(NX),J), TSIN(2,ABS(NY),J), TSIN(3,NZ,J)/)

                IF (NX < 0) VS(1) =-VS(1)
                IF (NY < 0) VS(2) =-VS(2)

                PC = VC(1)*VC(2)*VC(3) - VC(1)*VS(2)*VS(3) - VS(1)*VC(2)*VS(3) - VS(1)*VS(2)*VC(3)
                PS = VS(1)*VC(2)*VC(3) + VC(1)*VS(2)*VC(3) + VC(1)*VC(2)*VS(3) - VS(1)*VS(2)*VS(3)

                DUMMY = NX*RBSITES(1,1,J) + NY*RBSITES(2,1,J) + NZ*RBSITES(3,1,J)
                SUMC  = SUMC + PC*DUMMY
                SUMS  = SUMS + PS*DUMMY

                PEK = PEK + GUFCTR*FCTR(NX,NY,NZ)*(SUMC*SUMC + SUMS*SUMS)
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE DIPLR_DSCS_RECIP

! ========================================================================================================
! ========================================================================================================

SUBROUTINE DEF_GB_DSCS()
    USE COMMONS, ONLY: DP,REFSITE
    IMPLICIT NONE
    REFSITE(:,1) = [0.0_dp, 0.0_dp,  1.0_dp] 
END SUBROUTINE DEF_GB_DSCS

! ========================================================================================================
! ========================================================================================================

SUBROUTINE DEF_DIPLR()

    USE COMMONS, ONLY: DP,PI,TWOPI,DPMUSQ,BOX,GUFCTR,BOXL,RCUT,RCUTSQ,ALPHA,ALPSQ,NCSQMAX,NC,SLFFCT,FCTR,NVV

    IMPLICIT NONE

    INTEGER                    :: VN(3), NX, NY, NZ, NVI
    REAL(KIND=DP)              :: W

    BOXL    = BOX(1)
    RCUT    = BOXL / 2.0_dp
    RCUTSQ  = RCUT*RCUT
    ALPHA   = 5.0_dp / BOXL
    ALPSQ   = ALPHA*ALPHA
    NCSQMAX = INT(ALPSQ*RCUT*BOXL/PI)
    NC      = NCSQMAX
    GUFCTR  = 2.0_dp*PI*DPMUSQ/BOXL**3
    SLFFCT  = TWOPI * ALPHA**3 * DPMUSQ / (3.0_dp*SQRT(PI))

    ALLOCATE(NVV(-NC:NC,-NC:NC,0:NC),FCTR(-NC:NC,-NC:NC,0:NC))
    
    W = PI*PI/(BOXL*ALPHA)**2

    DO NZ = 0, NC
        DO NY =-NC, NC
            DO NX =-NC, NC
                VN  = (/NX, NY, NZ/)
                NVI = DOT_PRODUCT(VN,VN)
                NVV(NX,NY,NZ) = NVI
                IF (NVI == 0 .OR. NVI > NCSQMAX) CYCLE
                FCTR(NX,NY,NZ) = 2.0_dp*EXP(-W*NVI)/NVI    
                IF (NZ == 0) FCTR(NX,NY,NZ) = 0.5_dp*FCTR(NX,NY,NZ)
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE DEF_DIPLR

! ========================================================================================================
! ========================================================================================================

SUBROUTINE VIEW_DIPLR_DSCS()
    
    USE COMMONS, ONLY: DP, NPART, R, Q, REFSITE, BOX, VIEWUNIT
    USE ROTATIONS_MODULE, ONLY: Q_TO_RM

    IMPLICIT NONE

    INTEGER:: J1
    REAL(KIND=DP) :: RWRITE(3,NPART), RM(3,3), RBCOORDS1(3), RBCOORDS2(3)

    WRITE(VIEWUNIT,*) NPART*6
    WRITE(VIEWUNIT,*)

    RWRITE(:,:) = R(:,:)
    
    DO J1 = 1, NPART
        RWRITE(:,J1) = RWRITE(:,J1) - ANINT(RWRITE(:,J1)/BOX(:))*BOX(:)
    END DO

    DO J1 = 1, NPART
        RM = Q_TO_RM( Q(:,J1) )
        RBCOORDS1 = RWRITE(:,J1) + 0.05_dp*MATMUL(RM,REFSITE(:,1))
        RBCOORDS2 = RWRITE(:,J1) - 0.05_dp*MATMUL(RM,REFSITE(:,1))
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS1(1), RBCOORDS1(2), RBCOORDS1(3)
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RBCOORDS2(1), RBCOORDS2(2), RBCOORDS2(3)
        RBCOORDS1 = RWRITE(:,J1) + 0.25_dp*MATMUL(RM,[ SQRT(2.0_dp)/2.0_dp,-SQRT(2.0_dp)/2.0_dp,0.0_dp])
        RBCOORDS2 = RWRITE(:,J1) + 0.25_dp*MATMUL(RM,[-SQRT(2.0_dp)/2.0_dp, SQRT(2.0_dp)/2.0_dp,0.0_dp])
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RBCOORDS1(1), RBCOORDS1(2), RBCOORDS1(3)
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RBCOORDS2(1), RBCOORDS2(2), RBCOORDS2(3)
        RBCOORDS1 = RWRITE(:,J1) + 0.25_dp*MATMUL(RM,[ SQRT(2.0_dp)/2.0_dp, SQRT(2.0_dp)/2.0_dp,0.0_dp])
        RBCOORDS2 = RWRITE(:,J1) + 0.25_dp*MATMUL(RM,[-SQRT(2.0_dp)/2.0_dp,-SQRT(2.0_dp)/2.0_dp,0.0_dp])
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RBCOORDS1(1), RBCOORDS1(2), RBCOORDS1(3)
        WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'C ', RBCOORDS2(1), RBCOORDS2(2), RBCOORDS2(3)
    END DO
    
END SUBROUTINE