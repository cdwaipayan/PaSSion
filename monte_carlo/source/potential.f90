SUBROUTINE POTENTIAL(ENERGY)
    !==================================================================================
    !   Subroutine to calculate the potential energy associated either with the whole
    !   system or for a given particle whose index is given by INDXP.
    !==================================================================================
        USE COMMONS, ONLY: DP, NPART, NDIM, R, Q, BOX, SPET, INDXP, RCUTSQ, RIGIDT
        USE COMMONS, ONLY: RBSITES, REFSITE, NSITES, VIRTEMP, CELLLISTT, RACEMICT, REFSITE2
        USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRID
        
        USE CELL_LIST, ONLY: C_INDEX, NEIGHBOURS
        
        USE ROTATIONS_MODULE, ONLY: Q_TO_RM
        IMPLICIT NONE
    
    !   Counters for the loops   
        INTEGER       :: J1, J2, J1START, J1END, RC1
    !   Parameters related to the translational coordinates of the particles         
        REAL(KIND=DP) :: RI(NDIM), RJ(NDIM), RIJ(NDIM), RIJSQ
    !   Parameters related to the orientational coordinates of the particles        
        REAL(KIND=DP) :: RMI(NDIM,NDIM)
    !   Pair energy between particles J1 and J2
        REAL(KIND=DP) :: PAIR_ENERGY
    !   Parameters related to cell-list 
        INTEGER       :: CI(3), J_LIST(NPART), JJ
        
    !   PARAMETER TO BE OUTPUT FROM SUBROUTINE
        REAL(KIND=DP), INTENT(OUT) :: ENERGY
    
        ENERGY   = 0.0_dp
        VIRTEMP  = 0.0_dp
    
        IF (SPET) THEN
        !   Calculating the energy associated with a single particle Monte Carlo move.
        !   Therefore we only need to calculate the change in energy associated with 
        !   moving particle INDXP.
            J1START = INDXP 
            J1END   = INDXP
        
        ELSE
        !   Calculating the energy for the entire system
            J1START = 1
            IF(CELLLISTT) THEN
                J1END = NPART
            ELSE
                J1END = NPART-1
            ENDIF
            
            IF(RIGIDT) THEN
            !   Extract the directional unit vectors for the rigid sites on the particles
                DO J1 = 1, NPART
                    RMI     = Q_TO_RM( Q(:,J1) )
                    IF(RACEMICT) RC1 = MOD(((J1-1)-MOD((J1-1),12))/12+1,2)
                    DO J2 = 1, NSITES
                        IF(RACEMICT) THEN
                            IF(RC1==1) THEN
                                RBSITES(:,J2,J1) = MATMUL(RMI,REFSITE(:,J2))
                            ELSE
                                RBSITES(:,J2,J1) = MATMUL(RMI,REFSITE2(:,J2))
                            ENDIF
                        ELSE
                            RBSITES(:,J2,J1) = MATMUL(RMI,REFSITE(:,J2))
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
    
        ENDIF
    
        DO J1 = J1START, J1END
        !   Position of particle I
            RI  = R(:,J1)
    
            IF(CELLLISTT) THEN
    
                CI = C_INDEX ( RI/BOX )
                
            !   Put neighbours in j_list
                IF (SPET) THEN
                !   Calculating the total energy associated with particle I, therefore
                !   need to extract index of particles from ALL 26 neighbouring cells.
                    J_LIST = NEIGHBOURS ( NPART, J1, CI, .FALSE. )
                ELSE
                !   Calculating the total energy of the system, therefore we only extract
                !   the index of particles from 13 neighbouring cells for computational 
                !   efficiency.
                    J_LIST = NEIGHBOURS ( NPART, J1, CI, .TRUE. )
                ENDIF
    
                JJ = 0
    
            ELSE
                IF(SPET) THEN
                    JJ = 0
                ELSE
                    JJ = J1
                ENDIF
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
            
            !   If performing a cluster move then check if contribution to the pair energy 
            !   from the interaction between particles I and J has already been accounted for.
                IF( CLUSTERT .AND. ANY(CLSTR(1:CLSTRID-1)==J2) ) CYCLE
    
            !   Position of particle J
                RJ    = R(:,J2)
                RIJ   = RI - RJ
                RIJ   = RIJ - BOX*ANINT( RIJ/BOX )
                RIJSQ = DOT_PRODUCT(RIJ,RIJ)
            !   Check if the particle J is within the cutoff distance
                IF (RIJSQ <= RCUTSQ) THEN
                    PAIR_ENERGY = 0.0_dp
                    CALL PAIR_POTENTIAL(PAIR_ENERGY, RIJ, RIJSQ, J1, J2)
                    ENERGY = ENERGY + PAIR_ENERGY
                ENDIF ! Check if particles are within cutoff
            ENDDO ! Loop over particles j
        ENDDO ! Loop over particles i
    
        VIRTEMP = (1.0_dp/3.0_dp)*VIRTEMP
    
    END SUBROUTINE POTENTIAL
    
    SUBROUTINE PAIR_POTENTIAL(PAIR_ENERGY, RIJ, RIJSQ, J1, J2)
    
        USE COMMONS, ONLY: DP, NDIM, HARDT, OVERLAPT, MODT
        USE COMMONS, ONLY: KFT, GLJT, PGLJT, DMBLGLJT, KIHARAT, PGLJT, CPPT, ETPT, KFRECT, HDMBLT, HTPRT, YUKT
        
        IMPLICIT NONE
    
    !   Particle indices 
        INTEGER, INTENT(IN)        :: J1, J2
    !   Parameters related to the translational coordinates of the particles         
        REAL(KIND=DP), INTENT(IN)  :: RIJ(NDIM), RIJSQ
    !   Pair energy between particles J1 and J2
        REAL(KIND=DP), INTENT(OUT) :: PAIR_ENERGY
    
        PAIR_ENERGY   = 0.0_dp
    
        IF (HARDT) THEN
    !   ---------------------------------------------------
    !   Pair potentials which have hard-core repulsion
    !   ---------------------------------------------------
            IF (RIJSQ < 1.0_dp) THEN
        !   Check if hard particles are overlapping, if yes reject the move
                OVERLAPT = .TRUE.
                RETURN
            ELSE
        !   Kern-Frenkel patchy particles
                IF (KFT) THEN
                    CALL KF(PAIR_ENERGY, RIJ, RIJSQ, J1, J2)
        !   Kern-Frenkel patchy particles with rectangular patches
                ELSEIF (KFRECT) THEN
                    CALL KF_REC(PAIR_ENERGY, RIJ, RIJSQ, J1, J2)
                ELSEIF (HDMBLT) THEN
                    CALL HDMBL(PAIR_ENERGY, J1, J2)
                ELSEIF(HTPRT) THEN
                    CALL HARD_TPR(PAIR_ENERGY, RIJ, RIJSQ, J1, J2)
                ELSE
                    STOP "UNDEFINED HARD-CORE POTENTIAL"
                ENDIF
            ENDIF
    
        ELSE 
    !   ---------------------------------------------------
    !   Pair potentials which have soft-core repulsion
    !   ---------------------------------------------------    
            IF(YUKT) THEN
                CALL YUKAWA(PAIR_ENERGY, RIJSQ)

            ELSEIF(GLJT) THEN
                CALL GLJ(PAIR_ENERGY, RIJSQ)
    
            ELSEIF (KIHARAT) THEN
                CALL KIHARA(PAIR_ENERGY, RIJ, J1, J2)
    
            ELSEIF (PGLJT) THEN
                CALL PGLJ(PAIR_ENERGY, RIJ, RIJSQ, J1, J2)
    
            ELSEIF (CPPT) THEN
                CALL CPP(PAIR_ENERGY, RIJ, RIJSQ, J1, J2)
    
            ELSEIF (ETPT) THEN
                CALL ETP(PAIR_ENERGY, RIJ, RIJSQ, J1, J2)
    
            ELSEIF (DMBLGLJT) THEN
                CALL DMBL_GLJ(PAIR_ENERGY, J1, J2)
    
            ELSE
                PRINT *, "NO POTENTIAL SELECTED, STOPPING PROGRAM."
                STOP " Stopping in potential.f90"
            ENDIF
    
            IF(MODT) THEN
        !   When calculating free energies it is necessary to perform simulations
        !   at and near infinite temperature, where the system acts like an ideal gas.
        !   However, in this high temperature limit soft-core potentials are not 
        !   well defined as the energy will reach large values resulting in numerical 
        !   overflow. We can instead introduce a modified potential: Vm(r) = min{V(r),100}
        !   for temperatures above a certain limit where we expect the system to behave like
        !   an ideal gas. 
        !   See https://www.roma1.infn.it/~sciortif/didattica/SIMATOM/SIMATOM/free-energy.pdf
                PAIR_ENERGY = MIN(PAIR_ENERGY, 100.0_dp)
            ENDIF
    
        ENDIF
    
    END SUBROUTINE
    
    SUBROUTINE INIT_POT()
    
        USE COMMONS, ONLY: RCUT, RCUTSQ, GLJT, KIHARAT, KFT, PGLJT, CPPT, ETPT, DMBLGLJT, KFRECT, HDMBLT, HTPRT
        IMPLICIT NONE
    
        RCUTSQ = RCUT*RCUT
    
        IF (GLJT) THEN
            CALL DEF_GLJ()
        ELSEIF (KIHARAT) THEN
            CALL DEF_KIHARA()
        ELSEIF (KFT) THEN
            CALL DEF_KF()
        ELSEIF (PGLJT) THEN
            CALL DEF_PGLJ()
        ELSEIF (CPPT) THEN
            CALL DEF_CPP()
        ELSEIF (ETPT .OR. HTPRT) THEN
            CALL DEF_ETP()
        ELSEIF (DMBLGLJT) THEN
            CALL DEF_DMBL_GLJ()
        ELSEIF (KFRECT) THEN
            CALL DEF_KF_REC()
        ELSEIF (HDMBLT) THEN
            CALL DEF_HDMBL()
        ENDIF
    
    END SUBROUTINE INIT_POT
    
    SUBROUTINE CHECK_DIV(ISTEP)
    
        USE COMMONS, ONLY: DP, PE, SPET, DIV_TOL
        IMPLICIT NONE
    
        INTEGER         :: ISTEP
        REAL(KIND=DP)   :: TOT_ENE, ENE_DIFF
    
        SPET  = .FALSE.
        CALL POTENTIAL(TOT_ENE)
    
        ENE_DIFF = ABS(TOT_ENE - PE)
    
        IF(ENE_DIFF > DIV_TOL) THEN
            PRINT *, "ENERGY CALCULATED ALONG MC TRAJECTORY HAS DIVERGED"
            PRINT *, "ENERGY ALONG MC TRAJECTORY:", PE
            PRINT *, "TOTAL ENERGY OF THE SYSTEM:", TOT_ENE
            PRINT *, "MAGNITUDE OF DIVERGENCE:", ENE_DIFF
            PRINT *, "NUMBER OF MC CYCLES:", ISTEP
        ENDIF
    
        PE = TOT_ENE
    
    END SUBROUTINE CHECK_DIV
    
    SUBROUTINE PESRF()
    
        USE COMMONS, ONLY: CPPT, ETPT
        IMPLICIT NONE
    
        IF (CPPT) THEN
            CALL PESRF_CPP()
        ELSEIF (ETPT) THEN
            CALL PESRF_ETP()
        ELSE
            PRINT*, "NO POTENTIAL SELECTED, STOPPING PROGRAM."
            STOP " Stopping in potential.f90 (PESRF)"
        ENDIF
    
    END SUBROUTINE PESRF