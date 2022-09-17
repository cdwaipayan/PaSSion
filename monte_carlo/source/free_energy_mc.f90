MODULE FREE_ENERGY_MC
    !   ================================================================================================
    !   Calculate the free energy of crystal structures using the Frenkel-Ladd approach outlined in:
    !   Determination of phase diagrams via computer simulation: methodology and applications to water,
    !   electrolytes and proteins; Vega, C and Sanz, E and Abascal, JLF and Noya, EG;
    !   Journal of Physics: Condensed Matter; 20; 153101; 2008
    !   ------------------------------------------------------------------------------------------------
    !   Calculate the free energy of arbitrary systems using the Schilling-Schmid approach outlined in:
    !   Computing absolute free energies of disordered structures by molecular simulation; 
    !   T. Schilling and F. Schmid; J. Chem. Phys.; 131; 231102 (2009); doi: 10.1063/1.3274951
    !   ================================================================================================
    
        USE COMMONS, ONLY: DP, CDP, NDIM, NPART, NSITES, NSTEP, BLKLNGTH, PI, HLFPI, TWOPI, PE, OVERLAPT, BOX, RIGIDT, BETAKB, R, Q
        USE COMMONS, ONLY: RBSITES, REFSITE, RBREF, RREF, LAM_SS, LAMOR, RCUT_SS, RCUTSQ_SS, SSMOVERATIO, SSID, U_SS, U_EIN_OR
        USE COMMONS, ONLY: FLFET, SCHSMIT, LAMTR, FLID, LAMTR, DRCOM, U_EIN_TR, EXP_U_EIN
        USE COMMONS, ONLY: CLUSTERMOVET, CLSTRRATIO, LRGCLSTRT, LRGCLSTRRATIO, LRGCLSTMVT
        USE CELL_LIST, ONLY: C_INDEX, MOVE_IN_LIST
        IMPLICIT NONE
        
        PRIVATE
    
        PUBLIC :: FREE_ENERGY
    
    CONTAINS
    
        SUBROUTINE FREE_ENERGY()
    
            USE COMMONS, ONLY: U_LATT, EXP_U_SS, HST
            IMPLICIT NONE
    
            INTEGER        :: J1
            REAL(KIND=DP)  :: EXP_U
            LOGICAL        :: SWAPT
            REAL(KIND=CDP) :: DRAND48
    
            DO J1 = 1, NPART
                
                SWAPT = .FALSE.
                
                IF(SCHSMIT) THEN
                    IF( DRAND48() < 0.1_dp ) SWAPT = .TRUE.
                ENDIF
    
                IF(SWAPT) THEN
                    CALL SIMPLE_SWAP()
                ELSE
                    IF(FLID == 1 .OR. SSID == 1) THEN
                        CALL SIMPLE_MC_MOVE()
                    ELSE
                        
                        IF(SCHSMIT .AND. CLUSTERMOVET) THEN

                            IF(DRAND48() < CLSTRRATIO) THEN
                                IF(LRGCLSTRT) THEN
                                    IF(DRAND48() < LRGCLSTRRATIO) THEN
                                        LRGCLSTMVT = .TRUE.
                                    ELSE
                                        LRGCLSTMVT = .FALSE.
                                    ENDIF
                                ELSE
                                    LRGCLSTMVT = .FALSE.
                                ENDIF
                                CALL CLUSTER_MC_MOVE()
                            ELSE
                                
                                CALL SIMPLE_MC_MOVE()
                            ENDIF

                        ELSE
                            CALL SIMPLE_MC_MOVE()
                        ENDIF

                    ENDIF
                ENDIF
    
            ENDDO
    
            IF(FLID == 1 .OR. SSID == 1) THEN
            !   Calculate the potential energy of the interacting system
                OVERLAPT  = .FALSE.
                CALL POTENTIAL(PE)
    
                IF(OVERLAPT) THEN
            !   We can not always directly evaluate the potential energy of a hard sphere system
            !   as it is not possible to define infinity in the program.
            !   Instead we directly evaluate the value of the exponential where if U = inf, exp(-U) = 0.
                    EXP_U = 0.0_dp
                    OVERLAPT = .FALSE.
                ELSE
                    IF(HST) THEN
                        EXP_U = 1.0_dp
                    ELSE
                        EXP_U = EXP( BETAKB*(U_LATT - PE) )
                    ENDIF
                ENDIF

                IF(FLFET) THEN
                    EXP_U_EIN = EXP_U
                ELSE
                    EXP_U_SS  = EXP_U
                ENDIF

            ENDIF
    
        END SUBROUTINE
    
    !   ##############################################################################################################
    !   ##############################################################################################################
        FUNCTION LINEAR_WELL_ENERGY(RDIFF)
        !   ---------------------------------------------------------------
        !   Calculate the potential energy associated with particle I and
        !   its corresponding attractive linear well.
        !   ---------------------------------------------------------------
            IMPLICIT NONE
    
            REAL(KIND=DP), INTENT(IN)  :: RDIFF
            REAL(KIND=DP)              :: LINEAR_WELL_ENERGY
    
            LINEAR_WELL_ENERGY = MIN( RDIFF/RCUT_SS - 1.0_dp, 0.0_dp )
    
        END FUNCTION
    !   ##############################################################################################################
    !   ##############################################################################################################
        SUBROUTINE EINSTEIN_CRYSTAL_TR(PID, DELR, FIXEDCMT, ENERGY)
    !   --------------------------------------------------------------------------
    !    Calculate the translational component of the potential energy for a
    !    particle in an Einstein Crystal.
    !   --------------------------------------------------------------------------
            USE COMMONS, ONLY: R, RREF
    
            IMPLICIT NONE
    
            INTEGER, INTENT(IN)        :: PID
            REAL(KIND=DP), INTENT(IN)  :: DELR(NDIM)
            LOGICAL, INTENT(IN)        :: FIXEDCMT
    
            REAL(KIND=DP)              :: N, RI0(NDIM)
    
            REAL(KIND=DP), INTENT(OUT) :: ENERGY
    
            N   = REAL(NPART,DP)
            RI0 = R(:,PID) - RREF(:,PID) - DRCOM
    
            IF(FIXEDCMT) THEN
        !   Calculate the energy CHANGE associated with moving particle I by Δr under the constrainst
        !   of the center of mass being fixed in the Einstein crystal.
        !   This should ONLY be evoked when performing a Monte Carlo move with the center of mass fixed.
                ENERGY = ( 2.0_dp*DOT_PRODUCT(RI0,DELR) + (N-1.0_dp)/N * DOT_PRODUCT(DELR,DELR) )
            ELSE
        !   Calculate the energy associated with particle I in the Einstein crystal.
                ENERGY = DOT_PRODUCT(RI0, RI0)
            ENDIF
    
        END SUBROUTINE
    !   ##############################################################################################################
    !   ##############################################################################################################
        SUBROUTINE EINSTEIN_OR(REFSITES, NEW_RBSITES, ENERGY, REJECTT)
        !   --------------------------------------------------------------------------
        !    Calculate the potential energy associated with rotating a particle in an 
        !    Einstein Crystal with fixed center of mass.
        !   --------------------------------------------------------------------------
            USE COMMONS, ONLY: HEADTAILT
    
            IMPLICIT NONE
    
            REAL(KIND=DP), INTENT(IN)  :: NEW_RBSITES(NDIM,NSITES), REFSITES(NDIM,NSITES)
    
            INTEGER                    :: J1, JA, JB
            REAL(KIND=DP)              :: THETAA, THETAA_MIN, THETAB, THETAB_MIN
    
            REAL(KIND=DP), INTENT(OUT) :: ENERGY
            LOGICAL, INTENT(OUT)       :: REJECTT
    
            ENERGY = 0.0_dp
            
            IF(NSITES <= 2 .AND. (.NOT. HEADTAILT) ) THEN
            !   Particles with axial symmetry but do not have head-tail symmetry (C_inf-v)
                ENERGY = 1.0_dp - DOT_PRODUCT( REFSITES(:,1),NEW_RBSITES(:,1) )
                REJECTT = .FALSE.
                
            ELSE
                DO J1 = 1, NSITES
            !   Calculate the angle subtended by the position of reference site 1 and each of the new sites. 
            !   Additionally, we make sure to account for rounding errors in double precision floating point numbers 
                    THETAA = ACOS(MIN(MAX(DOT_PRODUCT(REFSITES(:,1),NEW_RBSITES(:,J1)), -1.0_dp), 1.0_dp))
                !   Keep track of which site gives the smallest angle
                    IF(J1 == 1 .OR. THETAA < THETAA_MIN) THEN
                        JA = J1
                        THETAA_MIN = THETAA
                    ENDIF
    
                    IF(NSITES > 2) THEN
            !   Calculate the angle subtended by the position of reference site 2 and each of the new sites. 
            !   Additionally, we make sure to account for rounding errors in double precision floating point numbers 
                        THETAB = ACOS(MIN(MAX(DOT_PRODUCT(REFSITES(:,2),NEW_RBSITES(:,J1)), -1.0_dp), 1.0_dp))
                    !   Keep track of the smallest angle
                        IF(J1 == 1 .OR. THETAB < THETAB_MIN) THEN
                            JB = J1
                            THETAB_MIN = THETAB
                        ENDIF
                    ENDIF
                ENDDO
                
            !   If the same rigid body site in the new configuration provides the smallest angle 
            !   between both of the reference sites, reject the move.
                IF(JA == JB) THEN
                    REJECTT = .TRUE.
                    RETURN
                ELSE
                    REJECTT = .FALSE.
                ENDIF
            
                IF(NSITES==2) THEN
            !   Particles have uniaxial symmetry and have head-tail symmetry (D_inf-h)
                    ENERGY = SIN(THETAA_MIN)**2
    
                ELSEIF(NSITES==3) THEN
            !   Particles with trigonal planar symmetry (D_3h)
                        ENERGY = ( SIN(3.0_dp*THETAA_MIN/2.0_dp)**2 + SIN(THETAB_MIN)**2 )
    
                ELSEIF(NSITES==4 .OR. NSITES==6) THEN
            !   Particles with octahedral (Oh) or tetrahedral (Td) symmetry.
                        ENERGY = ( SIN(THETAA_MIN)**2 + SIN(THETAB_MIN)**2 )
    
                ELSE
                    STOP "ERROR - NO ORIENTATIONAL FIELD SELECTED"
                ENDIF
    
            ENDIF
        
        END SUBROUTINE
    !   ##############################################################################################################
    !   ##############################################################################################################
    
    !   ##############################################################################################################
    !   ##############################################################################################################
        SUBROUTINE SIMPLE_MC_MOVE()
        !   ---------------------------------------------------------------
        !   This subroutine performs a single-particle Monte Carlo move.
        !   ---------------------------------------------------------------
            USE COMMONS, ONLY: CELLLISTT
            USE COMMONS, ONLY: INDXP, SPET, MAXDTR, MAXDRT, NTMOVES, NRMOVES, ACCPTCT, ACCPTCR
            USE ROTATIONS_MODULE, ONLY: RANDOM_ROTATE_QUATERNION, Q_TO_RM
            
            IMPLICIT NONE
    
            INTEGER             :: J1, PID, RC1
            
            REAL(KIND=DP)       :: WELL_O, WELL_N, WELL_DIFF, RI0(NDIM), DELR(NDIM), OLD_DRCOM(NDIM)
            REAL(KIND=DP)       :: DELE_EIN, OLD_EIN_OR, NEW_EIN_OR, RM(3,3), TDIFF
            REAL(KIND=DP)       :: OLD_ENE, NEW_ENE, DELE_POT
    
        !   Coordinates before the move takes places.
            REAL(KIND=DP)       :: RO(NDIM), QO(4), OLD_RBSITES(NDIM,NSITES)
        !   Coordinates after the move takes place.
            REAL(KIND=DP)       :: RN(NDIM), QN(4), NEW_RBSITES(NDIM,NSITES)
        
            !   Parameters related to the evaluation of random variables
            REAL(KIND=DP)       :: RANDMOVE
            REAL(KIND=CDP)      :: DRAND48
            LOGICAL             :: REJCT_WELL, REJCT_ORTN, REJCT_POT
    
            INTEGER             :: CI(3)
    
            PID = INT(NPART*DRAND48()) + 1
    
            REJCT_WELL = .FALSE.
            REJCT_ORTN = .FALSE.
            REJCT_POT  = .FALSE.
    
            IF(RIGIDT) THEN
        !   Randomly select whether to perform a translational or rotational move.
                RANDMOVE = DRAND48()
            ELSE
                RANDMOVE = 0.0_dp
            ENDIF 
    
        !   Calculate the energy associated with the particle before the move.
            IF(FLID == 2 .OR. SSID == 2) THEN
                INDXP     = PID
                SPET      = .TRUE.
                OVERLAPT  = .FALSE.
                CALL POTENTIAL(OLD_ENE)
                IF(FLFET) OLD_DRCOM = DRCOM
            ENDIF
    
            IF( RANDMOVE < 0.5_dp ) THEN
        !   TRANSLATIONAL MOVE
                NTMOVES = NTMOVES + 1
                DO J1 = 1, NDIM
                    DELR(J1) = (2.0_dp*DRAND48()-1.0_dp)*MAXDTR
                ENDDO
        
            !   Perform the MC move
                RO = R(:,PID)
                RN = RO + DELR
    
                IF(FLFET) THEN
                !   Calculate the change in the harmonic potential energy.
                    CALL EINSTEIN_CRYSTAL_TR(PID, DELR, .TRUE., DELE_EIN)
                    TDIFF = LAMTR*DELE_EIN
                ELSE
                !   Calculate the well energy associated with the
                !   selected particle BEFORE the MC move.
                    RI0    = RO - RREF(:,PID)
                    RI0    = RI0 - BOX*ANINT(RI0/BOX)
                    WELL_O = LINEAR_WELL_ENERGY(NORM2(RI0))
                !   Calculate the well energy associated with the
                !   selected particle AFTER the MC move.
                    RN     = RN - BOX*ANINT(RN/BOX)
                    RI0    = RN - RREF(:,PID)
                    RI0    = RI0 - BOX*ANINT(RI0/BOX)
                    WELL_N = LINEAR_WELL_ENERGY(NORM2(RI0))
                !   Calculate the difference in energy associated
                !   with the MC move.
                    WELL_DIFF = WELL_N - WELL_O
                    TDIFF = LAM_SS*WELL_DIFF
                ENDIF
        
            !   Metropolis acceptence criteria
                IF( DRAND48() < EXP(-BETAKB*TDIFF) ) THEN
                !   Accept the move so update the coordinates of the system 
                    R(:,PID) = RN
                    IF(FLFET) DRCOM = DRCOM + DELR / REAL(NPART,DP)
                !   Update the position of the particle in the cell list following the move
                    IF(CELLLISTT) THEN
                        CI  = C_INDEX ( RN/BOX ) ! NEW CELL INDEX
                        CALL MOVE_IN_LIST ( PID, CI )
                    ENDIF
                    IF(FLID == 1 .OR. SSID == 1) ACCPTCT = ACCPTCT + 1
                ELSE
                !   Reject the move and do nothing
                    REJCT_WELL = .TRUE.
                ENDIF
            ELSE
        !   ROTATIONAL MOVES
                NRMOVES    = NRMOVES + 1
            !   Save the orientational coordinates before the move
                QO          = Q(:,PID)
                OLD_RBSITES = RBSITES(:,:,PID)
            !   Calculate the energy before the move.
                CALL EINSTEIN_OR(RBREF(:,:,PID), OLD_RBSITES, OLD_EIN_OR, REJCT_ORTN)
            !   Perform a rotational move which satisfies the constrainst laid out in defining
            !   the harmonic orientational field.
                REJCT_ORTN = .TRUE.
                DO WHILE(REJCT_ORTN)
                !   Perform random perturbation in quaternion space
                    QN = RANDOM_ROTATE_QUATERNION ( MAXDRT, QO )
                !   Update the rigid body sites of the particle being displaced
                    RM   = Q_TO_RM( QN )
                    DO J1 = 1, NSITES
                        NEW_RBSITES(:,J1) = MATMUL(RM ,REFSITE(:,J1))
                    ENDDO
                !   Calculate the energy associated with the rotational move.
                    CALL EINSTEIN_OR(RBREF(:,:,PID), NEW_RBSITES, NEW_EIN_OR, REJCT_ORTN)
                ENDDO
            !   Metropolis acceptence criteria
                DELE_EIN = NEW_EIN_OR - OLD_EIN_OR
                IF( DRAND48() < EXP(-BETAKB*LAMOR*DELE_EIN) ) THEN
                !   Accept the move so update the coordinates of the system and the 
                !   shift in the center of mass.
                    Q(:,PID) = QN
                    RBSITES(:,:,PID) = NEW_RBSITES
                    IF(FLID == 1 .OR. SSID == 1) ACCPTCR = ACCPTCR + 1
                ELSE
                !   Reject the move and do nothing
                    REJCT_WELL = .TRUE.
                ENDIF
            ENDIF
    
        !----------------------------------------------------------------------------------------
        !   If calculating the free energy difference between the interacting linear well
        !   and the system of interest perform second metropolis check based on the
        !   potential energy of the system under consideration. 
        !----------------------------------------------------------------------------------------
            IF(FLID == 2 .OR. SSID == 2) THEN
                IF(.NOT. REJCT_WELL) THEN
                !   Calculate the new energy of the system
                    OVERLAPT  = .FALSE.
                    CALL POTENTIAL(NEW_ENE)
                !   Metropolis acceptence criteria
                    IF(OVERLAPT .OR. ISNAN(NEW_ENE)) THEN
                        REJCT_POT = .TRUE.
                    ELSE
                        DELE_POT = NEW_ENE - OLD_ENE
                        IF(DRAND48() < EXP(-BETAKB*DELE_POT)) THEN
                            REJCT_POT = .FALSE.
                            IF( RANDMOVE < 0.5_dp ) THEN
                            !   Keep track of the energy associated with the harmonic springs or linear well.
                                IF(FLFET) THEN
                                    U_EIN_TR = U_EIN_TR + DELE_EIN
                                ELSE
                                    U_SS = U_SS + WELL_DIFF
                                ENDIF
                            ELSE
                            !   Keep track of the energy associated with the orientational Einstein crystal.
                                U_EIN_OR = U_EIN_OR + DELE_EIN
                            ENDIF
                        !   Keep track of the energy associated with the interacting crystal.
                            PE    = PE + DELE_POT
                            IF(RANDMOVE < 0.5_dp) THEN
                                ACCPTCT = ACCPTCT + 1
                            ELSE 
                                ACCPTCR = ACCPTCR + 1
                            ENDIF
                        ELSE
                            REJCT_POT = .TRUE.
                        ENDIF
                    ENDIF
    
                    IF (REJCT_POT) THEN
                        IF(RANDMOVE < 0.5_dp) THEN
                        !   Restore previous position of the particle
                            R(:,PID) = RO
                            IF(FLFET) DRCOM = OLD_DRCOM
                        !   Update the position of the particle in the cell list following the move
                            IF(CELLLISTT) THEN
                                CI  = C_INDEX ( RO/BOX ) ! NEW CELL INDEX
                                CALL MOVE_IN_LIST ( PID, CI )
                            ENDIF
                        ELSE 
                        !   Restore previous orientation of the particle
                            Q(:,PID) = QO
                            RBSITES(:,:,PID) = OLD_RBSITES
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
            SPET = .FALSE.
        END SUBROUTINE
    !   ##############################################################################################################
    !   ##############################################################################################################
    
    !   ##############################################################################################################
    !   ##############################################################################################################
        SUBROUTINE CLUSTER_MC_MOVE()
        !   ---------------------------------------------------------------
        !   This subroutine performs a single-particle Monte Carlo move.
        !   ---------------------------------------------------------------
            USE COMMONS, ONLY: CELLLISTT, INDXP, SPET, CLUSTERT, CLSTR, CLSTRSZ, ACCPTCTC, ACCPTCRC
            USE ROTATIONS_MODULE, ONLY: RANDOM_ROTATE_QUATERNION, Q_TO_RM
            USE CLUSTER_MOVE
            
            IMPLICIT NONE
    
            INTEGER             :: J1, J2, J3, INDCLSTR, CLSTRO(NPART), CLSTRSZO, CI(3)
            REAL(KIND=DP)       :: WELL_O, WELL_N, OEIN_O, OEIN_N, WELL_DIFF, RI0(NDIM)
            REAL(KIND=DP)       :: DELE_EIN, EIN_OR, TOT_DIFF
            REAL(KIND=DP)       :: ENRGN, ENRGO, DELE, WS
            REAL(KIND=DP)       :: ROLD(NDIM,NPART), QOLD(4,NPART)
    
            REAL(KIND=DP), ALLOCATABLE :: ORBSITES(:,:,:)
    
        !   Parameters related to the evaluation of random variables
            REAL(KIND=DP)       :: RANDMOVE
            REAL(KIND=CDP)      :: DRAND48
            LOGICAL             :: REJCT_WELL, REJCT_ORTN, REJCT_POT
    
            INDXP = INT(NPART*DRAND48()) + 1
    
            SPET       = .TRUE.
            CLUSTERT   = .TRUE.
            OVERLAPT   = .FALSE.
            REJCT_WELL = .FALSE.
            REJCT_ORTN = .FALSE.
            REJCT_POT  = .FALSE.
    
            WELL_O = 0.0_dp
            WELL_N = 0.0_dp
            OEIN_O = 0.0_dp 
            OEIN_N = 0.0_dp
    
            IF(RIGIDT) THEN
        !   Randomly select whether to perform a translational or rotational move.
                RANDMOVE = DRAND48()
            ELSE
                RANDMOVE = 0.0_dp
            ENDIF 
    
            INDCLSTR = INDXP
        !   Calculate the potential energy of the cluster with particle I as the seed before the move.
            CALL CLUSTERENERGY(RANDMOVE, ENRGO, WS)
        !   Save the size of the cluster and the IDs of particles in the cluster before the move.
            CLSTRO   = CLSTR
            CLSTRSZO = CLSTRSZ
        !   Save the old rigid body sites of the particles
            IF(RANDMOVE >= 0.5_dp) THEN
                ALLOCATE(ORBSITES(NDIM,NSITES,CLSTRSZ))
                DO J1 = 1, CLSTRSZ
                    J2 = CLSTR(J1)
                    ORBSITES(:,:,J1) = RBSITES(:,:,J2)
                ENDDO
            ENDIF
        !   Move the cluster which we have just constructed. 
        !   We also return the previous positions and orientations of the particles in the cluster.
            CALL CLUSTERMOVE(RANDMOVE, ROLD, QOLD)
        !   Reset the current index to the seed particle
            INDXP = INDCLSTR
        !   Calculate the potential energy of the cluster with particle I as the seed before the move.
            CALL CLUSTERENERGY(RANDMOVE, ENRGN, WS)
    
            IF(OVERLAPT .OR. ISNAN(ENRGN) .OR. (CLSTRSZO /= CLSTRSZ)) THEN
                REJCT_POT = .TRUE.
            ELSE
                DELE = ENRGN - ENRGO
            !   Metropolis acceptence criteria
                IF (DRAND48() .GE. EXP(-BETAKB*DELE)) REJCT_POT = .TRUE.
            ENDIF
    
            IF(.NOT. REJCT_POT) THEN
                DO J1 = 1, CLSTRSZ
                    J2 = CLSTRO(J1)
                    J3 = CLSTR(J1)
                !   LINEAR WELL ENERGY
                !   Calculate the energy of cluster in their respective linear wells before the move.
                    RI0    = ROLD(:,J2) - RREF(:,J2)
                    RI0    = RI0 - BOX*ANINT(RI0/BOX)
                    WELL_O = WELL_O + LINEAR_WELL_ENERGY(NORM2(RI0))
                !   Calculate the energy of cluster in their respective linear wells before the move.
                    RI0    = R(:,J3) - RREF(:,J3)
                    RI0    = RI0 - BOX*ANINT(RI0/BOX)
                    WELL_N = WELL_N + LINEAR_WELL_ENERGY(NORM2(RI0))
                !   ORIENTATIONAL HARMONIC FIELD ENERGY
                    IF(RANDMOVE >= 0.5_dp) THEN
                    !   Calculate the orientational harmonic energy of the particles before the move.
                        CALL EINSTEIN_OR(RBREF(:,:,J2), ORBSITES(:,:,J1), EIN_OR, REJCT_ORTN)
                        OEIN_O = OEIN_O + EIN_OR
                    !   Calculate the orientational harmonic energy of the particles after the move.
                        CALL EINSTEIN_OR(RBREF(:,:,J3), RBSITES(:,:,J3), EIN_OR, REJCT_ORTN)
                        IF(REJCT_ORTN) EXIT
                        OEIN_N = OEIN_N + EIN_OR
                    ENDIF
                ENDDO
            !   Metropolis acceptence criteria
                IF(.NOT. REJCT_ORTN) THEN
                    WELL_DIFF = WELL_N - WELL_O
                    TOT_DIFF  = LAM_SS*WELL_DIFF
                    IF(RANDMOVE >= 0.5_dp) THEN
                        DELE_EIN = OEIN_N - OEIN_O
                        TOT_DIFF = TOT_DIFF + LAMOR*DELE_EIN
                    ENDIF
                    IF( DRAND48() .GE. EXP( -BETAKB * TOT_DIFF ) ) REJCT_WELL = .TRUE.
                ENDIF
            ENDIF
    
            IF(REJCT_POT .OR. REJCT_WELL .OR. REJCT_ORTN) THEN
                DO J1 = 1, CLSTRSZO
                    J2 = CLSTRO(J1)
                    R(:,J2) = ROLD(:,J2)
                    IF(RANDMOVE >= 0.5_dp) THEN
                        Q(:,J2) = QOLD(:,J2)
                        RBSITES(:,:,J2) = ORBSITES(:,:,J1)
                    ENDIF
                !   Restore the position of the particle in the cell list
                    IF(CELLLISTT) THEN
                        CI  = C_INDEX ( R(:,J2)/BOX ) ! NEW CELL INDEX
                        CALL MOVE_IN_LIST ( J2, CI )
                    ENDIF
                ENDDO
            ELSE
                IF(RANDMOVE < 0.5_dp) THEN
                    ACCPTCTC = ACCPTCTC + 1
                ELSE
                    ACCPTCRC = ACCPTCRC + 1
                ENDIF
            !   Keep track of the energy associated with the interacting crystal.
                PE   = PE + DELE
                U_SS = U_SS + WELL_DIFF
                IF(RANDMOVE >= 0.5_dp) U_EIN_OR = U_EIN_OR + DELE_EIN
            ENDIF
    
            IF(RANDMOVE >= 0.5_dp) DEALLOCATE(ORBSITES)
            SPET     = .FALSE.
            CLUSTERT = .FALSE.
    
        END SUBROUTINE
    !   ##############################################################################################################
    !   ##############################################################################################################
    
    !   ##############################################################################################################
    !   ##############################################################################################################
        SUBROUTINE SIMPLE_SWAP()

            USE COMMONS, ONLY: ROT_SWITCH
        
            IMPLICIT NONE
    
            INTEGER             :: PID1, PID2
        !   Translational parameters
            REAL(KIND=DP)       :: WELL1_O, WELL2_O, WELL1_N, WELL2_N, WELL_DIFF
            REAL(KIND=DP)       :: R1(NDIM), R2(NDIM), RRF1(NDIM), RRF2(NDIM)
        !   Rotational parameters
            REAL(KIND=DP)       :: R1RR1(NDIM), R2RR2(NDIM), R1RR2(NDIM), R2RR1(NDIM)
            REAL(KIND=DP)       :: RBS1(NDIM,NSITES), RBR1(NDIM,NSITES), RBS2(NDIM,NSITES), RBR2(NDIM,NSITES)
            REAL(KIND=DP)       :: EIN1_O, EIN2_O, EIN1_N, EIN2_N, EIN_DIFF
        !   Parameters related to the evaluation of random variables
            REAL(KIND=DP)       :: RANDMOVE, DIFFTOL
            REAL(KIND=CDP)      :: DRAND48
            LOGICAL             :: REJCT_ORTN, TMOVE, RMOVE
    
    
        !   Randomly select whether to perform a translational and/or rotational swap.
            IF(RIGIDT .AND. ROT_SWITCH) THEN
                RANDMOVE = DRAND48()
            ELSE
                RANDMOVE = 0.0_dp
            ENDIF 
    
            IF( RANDMOVE < (1.0_dp/3.0_dp) ) THEN
            !   Translational swap
                TMOVE = .TRUE. 
                RMOVE = .FALSE.
            ELSEIF( RANDMOVE < (2.0_dp/3.0_dp) ) THEN
            !   Rotational swap
                TMOVE = .FALSE. 
                RMOVE = .TRUE.
            ELSE
            !   Translational + rotational swap
                TMOVE = .TRUE. 
                RMOVE = .TRUE.
            ENDIF
    
            WELL_DIFF = 0.0_dp
            EIN_DIFF  = 0.0_dp
    
        !   Randomly select particle I
            PID1 = INT(NPART*DRAND48()) + 1
            PID2 = PID1
        !   Randomly select particle J, such that I /= J
            DO WHILE(PID1 == PID2)
                PID2 = INT(NPART*DRAND48()) + 1
            ENDDO
    
        !   TRANSLATIONAL SWAP
            IF( TMOVE ) THEN
                R1      = R(:,PID1)
                RRF1    = RREF(:,PID1)
                R2      = R(:,PID2)
                RRF2    = RREF(:,PID2)
    
                R1RR1   = R1 - RRF1
                R1RR1   = R1RR1 - BOX*ANINT(R1RR1/BOX)
                WELL1_O = LINEAR_WELL_ENERGY(NORM2(R1RR1))
    
                R2RR2   = R2 - RRF2
                R2RR2   = R2RR2 - BOX*ANINT(R2RR2/BOX)
                WELL2_O = LINEAR_WELL_ENERGY(NORM2(R2RR2))
    
                R1RR2   = R1 - RRF2
                R1RR2   = R1RR2 - BOX*ANINT(R1RR2/BOX)
                WELL1_N = LINEAR_WELL_ENERGY(NORM2(R1RR2))
    
                R2RR1   = R2 - RRF1
                R2RR1   = R2RR1 - BOX*ANINT(R2RR1/BOX)
                WELL2_N = LINEAR_WELL_ENERGY(NORM2(R2RR1))
    
                WELL_DIFF = (WELL1_N+WELL2_N)-(WELL1_O+WELL2_O)
            ENDIF
    
        !   ROTATIONAL SWAP
            IF( RMOVE ) THEN
                RBS1 = RBSITES(:,:,PID1)
                RBR1 = RBREF(:,:,PID1)
    
                RBS2 = RBSITES(:,:,PID2)
                RBR2 = RBREF(:,:,PID2)
            
                CALL EINSTEIN_OR(RBR1, RBS1, EIN1_O, REJCT_ORTN)
                CALL EINSTEIN_OR(RBR2, RBS2, EIN2_O, REJCT_ORTN)
    
                CALL EINSTEIN_OR(RBR2, RBS1, EIN1_N, REJCT_ORTN)
                IF(REJCT_ORTN .EQV. .TRUE.) RETURN
                CALL EINSTEIN_OR(RBR1, RBS2, EIN2_N, REJCT_ORTN)
                IF(REJCT_ORTN .EQV. .TRUE.) RETURN
    
                EIN_DIFF = (EIN1_N+EIN2_N)-(EIN1_O+EIN2_O)
            ENDIF
    
            IF( RANDMOVE < (1.0_dp/3.0_dp) ) THEN
                DIFFTOL = LAM_SS*WELL_DIFF
            ELSEIF( RANDMOVE < (2.0_dp/3.0_dp) ) THEN
                DIFFTOL = LAMOR*EIN_DIFF
            ELSE
                DIFFTOL = LAM_SS*WELL_DIFF + LAMOR*EIN_DIFF
            ENDIF
    
        !   Metropolis acceptence criteria
            IF( DRAND48() < EXP(-BETAKB*DIFFTOL) ) THEN
                IF(TMOVE) THEN
                    RREF(:,PID1) = RRF2
                    RREF(:,PID2) = RRF1
                    IF(SSID == 2) U_SS = U_SS + WELL_DIFF
                ENDIF
    
                IF(RMOVE) THEN
                    RBREF(:,:,PID1) = RBR2
                    RBREF(:,:,PID2) = RBR1
                    IF(SSID == 2) U_EIN_OR = U_EIN_OR + EIN_DIFF
                ENDIF
            ENDIF
    
        END SUBROUTINE
    !   ##############################################################################################################
    !   ##############################################################################################################
    END MODULE
    