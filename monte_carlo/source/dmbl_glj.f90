SUBROUTINE DMBL_GLJ(ENERGY, J1, J2)

    USE COMMONS, ONLY: DP, R, NDIM, NSITES, GLJN, RBSITES, BOX, DNAT
    USE COMMONS, ONLY: SIGIJ2, EPSIJ, DMBLCUTSQ
    USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, CLSTRID, CLSTRADJ, VLMCLUSTERMOVET, CLURIJ, CLSTRSITEID

    IMPLICIT NONE

    INTEGER, INTENT(IN)       :: J1, J2

    INTEGER       :: J3, J4
    REAL(KIND=DP) :: RI(NDIM), RJ(NDIM), RIJ(NDIM), R2, RLJN, R2LJN, RABSQ, VIJ
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

            IF (RABSQ <= DMBLCUTSQ(J3,J4)) THEN

                R2      = SIGIJ2(J3,J4)/RABSQ
                RLJN    = R2**(GLJN/2.0_dp)
                R2LJN   = RLJN*RLJN

                IF(DNAT .AND. J3 == J4) THEN
                    VIJ = 4.0_dp*(R2LJN - RLJN + 0.25_dp)
                ELSE
                    VIJ = 4.0_dp*EPSIJ(J3,J4)*(R2LJN - RLJN)
                ENDIF
    
                ENERGY  = ENERGY + VIJ
            !   If performing MC with cluster-moves, add particle J to the current 
            !   cluster (if it is not already in the cluster). 
                IF(CLUSTERT .OR. VLMCLUSTERMOVET) THEN
                !   Only consider particles with BB bonds to be apart of the same cluster.
                    IF( J3 == CLSTRSITEID .AND. J4 == CLSTRSITEID ) THEN
                        IF(VLMCLUSTERMOVET) THEN
                            RIJ = RI - RJ
                            RIJ = RIJ - BOX*ANINT( RIJ/BOX )
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

            ENDIF
        ENDDO ! Loop over each of the patches on particle j
    ENDDO ! Loop over particles j

END SUBROUTINE DMBL_GLJ

!====================================================================================================
    
SUBROUTINE DEF_DMBL_GLJ()
    !----------------------------------------------------------------
    ! 
    !----------------------------------------------------------------
    
        USE COMMONS, ONLY: DP, PI, NSITES, REFSITE, RCUT, RCUTSQ
        USE COMMONS, ONLY: SIGAA, SIGBB, SIGCC, SIGIJ, SIGIJ2, EPSAA, EPSBB, EPSCC, EPSIJ, DMBLCUTSQ
    
        IMPLICIT NONE

        INTEGER :: J1, J2

        ALLOCATE(SIGIJ(NSITES,NSITES), SIGIJ2(NSITES,NSITES), EPSIJ(NSITES,NSITES), DMBLCUTSQ(NSITES,NSITES))

        SIGAA        = 1.0_dp
        EPSAA        = 1.0_dp
        
        SIGIJ(1,1) = SIGAA
        SIGIJ(1,2) = (SIGAA + SIGBB) / 2.0_dp
        SIGIJ(2,1) = SIGIJ(1,2)
        SIGIJ(2,2) = SIGBB

        EPSIJ(1,1) = EPSAA
        EPSIJ(1,2) = SQRT(EPSAA*EPSBB)
        EPSIJ(2,1) = EPSIJ(1,2)
        EPSIJ(2,2) = EPSBB

        REFSITE(:,1)= [ 0.0_dp, 0.0_dp, 0.0_dp     ]
        REFSITE(:,2)= [ 0.0_dp, 0.0_dp, SIGIJ(1,2) ]

        DMBLCUTSQ(1,1) = RCUT**2
        DMBLCUTSQ(1,2) = ( RCUT*SIGIJ(1,2) )**2
        DMBLCUTSQ(2,1) = DMBLCUTSQ(1,2)
        DMBLCUTSQ(2,2) = ( RCUT*SIGIJ(2,2) )**2


        IF(NSITES == 2) THEN
            RCUT   = (1.0_dp + 2.0_dp*SIGBB) + SIGBB*RCUT
            RCUTSQ = RCUT*RCUT
        ELSE
        ! Triblock particles

            SIGIJ(1,3) = (SIGAA + SIGCC) / 2.0_dp
            SIGIJ(3,1) = SIGIJ(1,3)
            SIGIJ(2,3) = (SIGBB + SIGCC) / 2.0_dp
            SIGIJ(3,2) = SIGIJ(2,3)
            SIGIJ(3,3) = SIGCC

            EPSIJ(1,3) = SQRT(EPSAA*EPSCC)
            EPSIJ(3,1) = EPSIJ(1,3)
            EPSIJ(2,3) = SQRT(EPSBB*EPSCC)
            EPSIJ(3,2) = EPSIJ(2,3)
            EPSIJ(3,3) = EPSCC

            REFSITE(:,3)= [ 0.0_dp, 0.0_dp, -SIGIJ(1,3) ]

            DMBLCUTSQ(1,3) = (RCUT*SIGIJ(1,3))**2
            DMBLCUTSQ(3,1) = DMBLCUTSQ(1,3)
            DMBLCUTSQ(2,3) = (RCUT*SIGIJ(2,3))**2
            DMBLCUTSQ(3,2) = DMBLCUTSQ(2,3)
            DMBLCUTSQ(3,3) = (RCUT*SIGIJ(3,3))**2

            RCUT   = (1.0_dp + 2.0_dp*MAX(SIGBB,SIGCC)) + MAX(SIGBB,SIGCC)*RCUT
            RCUTSQ = RCUT*RCUT

        ENDIF

        DO J1 = 1, NSITES
            DO J2 = 1, NSITES
                SIGIJ2(J1,J2) = SIGIJ(J1,J2)*SIGIJ(J1,J2)
            ENDDO
        ENDDO
        
    END SUBROUTINE
    
!====================================================================================================
    
    SUBROUTINE VIEW_DMBL_GLJ()
    
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
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'N ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==2) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'O ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ELSEIF(J2==3) THEN
                    WRITE(VIEWUNIT,'(A5,1X,3F12.7)') 'F ', RBCOORDS(1), RBCOORDS(2), RBCOORDS(3)
                ENDIF
            ENDDO
        END DO
        
    END SUBROUTINE