    SUBROUTINE BINARY_HS(J1, J2, RIJSQ)
        
        USE COMMONS, ONLY: DP, SIGIJ2, BNRYA, OVERLAPT

        IMPLICIT NONE
        
        INTEGER, INTENT(IN)        :: J1, J2
        REAL(KIND=DP), INTENT(IN)  :: RIJSQ
        INTEGER                    :: J3, J4

        IF(J1<=BNRYA) THEN
            J3 = 1
        ELSE
            J3 = 2
        ENDIF

        IF(J2<=BNRYA) THEN
            J4 = 1
        ELSE
            J4 = 2
        ENDIF

        IF( RIJSQ<SIGIJ2(J3,J4) ) OVERLAPT = .TRUE.

    END SUBROUTINE BINARY_HS

    SUBROUTINE BINARY_GLJ(ENERGY, J1, J2, RIJSQ)
    
        USE COMMONS, ONLY: DP, GLJN, VIRTEMP, SIGIJ2, EPSIJ, BNRYA
    
        IMPLICIT NONE
        
        INTEGER, INTENT(IN)        :: J1, J2
        REAL(KIND=DP), INTENT(IN)  :: RIJSQ
        INTEGER                    :: J3, J4
        REAL(KIND=DP)              :: DIST, R2, RLJN, R2LJN, VIJ, WIJ
        REAL(KIND=DP), INTENT(OUT) :: ENERGY

        IF(J1<=BNRYA) THEN
            J3 = 1
        ELSE
            J3 = 2
        ENDIF

        IF(J2<=BNRYA) THEN
            J4 = 1
        ELSE
            J4 = 2
        ENDIF
        
        R2      = SIGIJ2(J3,J4)/RIJSQ
        RLJN    = R2**(GLJN/2.0_dp)
        R2LJN   = RLJN*RLJN
        VIJ     = R2LJN - RLJN
        
        ENERGY  = 4.0_dp*EPSIJ(J3,J4)*VIJ
    
    !   Pair virial function w(r) = r dv(r)/dr
        WIJ     = VIJ + R2LJN
        VIRTEMP = VIRTEMP + 4.0_dp*EPSIJ(J3,J4)*GLJN*WIJ
    
    END SUBROUTINE BINARY_GLJ

    SUBROUTINE BINARY_YUKAWA(ENERGY, J1, J2, RIJSQ)
        USE COMMONS, ONLY: DP, VIRTEMP, YUKKAP, SIGIJ, BNRYA, BNRYEQUIT

        IMPLICIT NONE

        INTEGER, INTENT(IN)        :: J1, J2
        REAL(KIND=DP), INTENT(IN)  :: RIJSQ    
        INTEGER                    :: J3, J4   
        REAL(KIND=DP)              :: DIST, UYUK
        REAL(KIND=DP), INTENT(OUT) :: ENERGY

        ENERGY   = 0.0_dp

        IF(J1<=BNRYA) THEN
            J3 = 1
        ELSE
            J3 = 2
        ENDIF

        IF(J2<=BNRYA) THEN
            J4 = 1
        ELSE
            J4 = 2
        ENDIF

        DIST = SQRT(RIJSQ)
        IF(BNRYEQUIT) THEN
            IF(DIST<SIGIJ(J3,J4)) ENERGY = 10.0_dp
        ELSE
        !   Calculate repulsive contribution to the pair energy
            UYUK    = EXP( -YUKKAP*(DIST - SIGIJ(J3,J4)) )
            ENERGY  = ENERGY + UYUK / DIST
            VIRTEMP = VIRTEMP + UYUK*(YUKKAP + SIGIJ(J3,J4)/DIST)*DIST
        ENDIF

    END SUBROUTINE
    
    SUBROUTINE DEF_BINARY_SPHERES()
    
        USE COMMONS, ONLY: DP, NPART, BNRYR, BNRYA, SIGIJ, SIGIJ2, SIGAA, SIGBB, HARDT, EPSAA, EPSBB, EPSIJ, RCUT, RCUTSQ
        USE COMMONS, ONLY: BYUKWAT

        IMPLICIT NONE
        
        INTEGER :: J1, J2

        ALLOCATE(SIGIJ(2,2), SIGIJ2(2,2))

        SIGIJ(1,1) = SIGAA
        SIGIJ(1,2) = (SIGAA + SIGBB) / 2.0_dp
        SIGIJ(2,1) = SIGIJ(1,2)
        SIGIJ(2,2) = SIGBB

        IF(.NOT. BYUKWAT) THEN
            
            BNRYA = INT(BNRYR*REAL(NPART,DP))

            DO J1 = 1, 2
                DO J2 = 1, 2
                    SIGIJ2(J1,J2) = SIGIJ(J1,J2)*SIGIJ(J1,J2)
                ENDDO
            ENDDO
            
            RCUT   = MAX(SIGAA,SIGBB)*RCUT
            RCUTSQ = RCUT*RCUT

            IF(.NOT. HARDT) THEN
                ALLOCATE(EPSIJ(2,2))
                EPSIJ(1,1) = EPSAA
                EPSIJ(1,2) = SQRT(EPSAA*EPSBB)
                EPSIJ(2,1) = EPSIJ(1,2)
                EPSIJ(2,2) = EPSBB
            ENDIF
        ENDIF

    END SUBROUTINE
    
          !====================================================================================================
    
    SUBROUTINE VIEW_BINARY_SPHERES()
    
        USE COMMONS, ONLY: DP, NPART, R, BOX, VIEWUNIT, BNRYA
    
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
            IF(I<=BNRYA) THEN
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RWRITE(1,I), RWRITE(2,I), RWRITE(3,I)
            ELSE
                WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RWRITE(1,I), RWRITE(2,I), RWRITE(3,I)
            ENDIF
        END DO
        
    END SUBROUTINE