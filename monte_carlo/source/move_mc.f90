

    SUBROUTINE MV_MC()

    !     This subroutine performs a single-particle Monte Carlo move.

        USE COMMONS, ONLY: DP, CDP, NDIM, NPART, R, Q, BETAKB, BOX, MAXDTR, MAXDRT, SPET, RACEMICT, REFSITE2, &
        NTMOVES, NRMOVES, ACCPTCT, ACCPTCR, VIRTEMP, COLLDT, PINT, PINKAP2, PINKAPA, TRQ, CELLLISTT, RQ, RIGIDT, PE, &
        INDXP, VIR, OVERLAPT, NSITES, RBSITES, REFSITE, CLUSTERT, HALFST, EQUISEEDT, CLUSTERMOVET, CLSTRRATIO
        USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, ACCPTCTC, ACCPTCRC, LRGCLSTRT, LRGCLSTRRATIO, LRGCLSTMVT
        USE CELL_LIST, ONLY: C_INDEX, MOVE_IN_LIST
        USE ROTATIONS_MODULE, ONLY: RANDOM_ROTATE_QUATERNION, Q_TO_RM
        USE CLUSTER_MOVE
        
        IMPLICIT NONE
        
        INTEGER             :: J1, RC1
        REAL(KIND=DP)       :: RO(NDIM)
        REAL(KIND=DP)       :: WSN, WSO, ENRGN, ENRGO, DELE, BLTZMN, VIROLD, BIASO, BIASN
        REAL(KIND=DP)       :: RANDMOVE, TRQN, DELQR, DELQR2, UQR, TRQOJ, TRQNJ
        REAL(KIND=DP)       :: QO(4), RM(3,3)

        REAL(KIND=CDP)      :: DRAND48

        LOGICAL             :: REJECTT

        COMPLEX(KIND=DP)    :: QJO, QJN, RQN
        INTEGER             :: CI(3)

        INTEGER             :: J2, J3, INDCLSTR, CLSTRO(NPART), CLSTRSZO
        REAL(KIND=DP)       :: ROLD(NDIM,NPART), QOLD(4,NPART)

        SPET     = .TRUE.
        REJECTT  = .FALSE.
        OVERLAPT = .FALSE.
        CLUSTERT = .FALSE.

        ENRGN   = 0.0_dp ! New energy
        ENRGO   = 0.0_dp ! Old energy
        TRQNJ   = 0.0_dp 
        TRQOJ   = 0.0_dp
        QJO     = 0.0_dp
        QJN     = 0.0_dp
        VIROLD  = VIR

        WSO = 0.0_dp
        WSN = 0.0_dp

        IF(CLUSTERMOVET) THEN
            IF(DRAND48() < CLSTRRATIO) CLUSTERT = .TRUE.
            IF(LRGCLSTRT .AND. CLUSTERT) THEN
                IF(DRAND48() < LRGCLSTRRATIO) THEN
                    LRGCLSTMVT = .TRUE.
                ELSE
                    LRGCLSTMVT = .FALSE.
                ENDIF
            ELSE
                LRGCLSTMVT = .FALSE.
            ENDIF
        ENDIF

        INDXP   = INT(NPART*DRAND48()) + 1

        IF(RIGIDT) THEN
            RANDMOVE = DRAND48()
        ELSE
            RANDMOVE = 0.0_dp
        ENDIF 

    !   #####################################################################################################################
    !   CLUSTER MOVE
        IF(CLUSTERT) THEN
            INDCLSTR = INDXP
        !   Calculate the potential energy of the cluster with particle I as the seed before the move.
            CALL CLUSTERENERGY(RANDMOVE, ENRGO, WSO, QJO, TRQOJ, BIASO)
        !   Save the size of the cluster and the IDs of particles in the cluster before the move.
            CLSTRO   = CLSTR
            CLSTRSZO = CLSTRSZ
        !   Move the cluster which we have just constructed. 
        !   We also return the previous positions and orientations of the particles in the cluster.
            CALL CLUSTERMOVE(RANDMOVE, ROLD, QOLD)
        !   Reset the current index to the seed particle
            INDXP = INDCLSTR
        !   Calculate the potential energy of the cluster with particle I as the seed before the move.
            CALL CLUSTERENERGY(RANDMOVE, ENRGN, WSN, QJN, TRQNJ, BIASN)
    !   ------------------------------------------------------------------------------------------------------------
    !   SINGLE-PARTICLE MOVE
    !   ------------------------------------------------------------------------------------------------------------
        ELSE
        !   Calculate energy associated with particle I before the move
            CALL POTENTIAL(ENRGO)
            WSO  = VIRTEMP

        !   Calculate bias potentials
            IF(PINT .AND. COLLDT) THEN
            !   If performing an interface-pinning simulation calculate the collective
            !   density field for the particle.
                CALL COLLDENSFIELD(TRQOJ,QJO,INDXP)
            ELSEIF( HALFST .OR. EQUISEEDT ) THEN
            !   If performing an equilibration of the interface between two phases
            !   (bulk solid or solid seed + bulk liquid) calculate the harmonic potential
            !   energy associated with "trapped" particles. 
                CALL HARMONIC_TRAP(BIASO)
            ENDIF

        !   Perform random move
            IF( RANDMOVE < 0.5_dp ) THEN
            !   Translational moves
                DO J1 = 1, NDIM
                    RO(J1)      = R(J1,INDXP)
                    R(J1,INDXP) = R(J1,INDXP) + (2.0_dp*DRAND48()-1.0_dp)*MAXDTR
            !   Pick up the central image
                    R(J1,INDXP) = R(J1,INDXP) - BOX(J1)*ANINT(R(J1,INDXP)/BOX(J1))
                ENDDO
            !   Update the position of the particle in the cell list following the move
                IF(CELLLISTT) THEN
                    CI  = C_INDEX ( R(:,INDXP)/BOX ) ! NEW CELL INDEX
                    CALL MOVE_IN_LIST ( INDXP, CI )
                ENDIF
                NTMOVES = NTMOVES + 1
            ELSE
            !   Rotational moves
            !   Perform random perturbation in quaternion space
                QO         = Q(:,INDXP)
                Q(:,INDXP) = RANDOM_ROTATE_QUATERNION ( MAXDRT, QO )
            !   Update the rigid body sites of the particle being displaced
                RM   = Q_TO_RM( Q(:,INDXP) )
                IF(RACEMICT) RC1 = MOD(((INDXP-1)-MOD((INDXP-1),12))/12+1,2)
                DO J1 = 1, NSITES
                    IF(RACEMICT) THEN
                        IF(RC1==1) THEN
                            RBSITES(:,J1,INDXP) = MATMUL(RM ,REFSITE(:,J1))
                        ELSE
                            RBSITES(:,J1,INDXP) = MATMUL(RM ,REFSITE2(:,J1))
                        ENDIF
                    ELSE
                        RBSITES(:,J1,INDXP) = MATMUL(RM ,REFSITE(:,J1))
                    ENDIF
                ENDDO
                NRMOVES    = NRMOVES + 1
            ENDIF
        !   Calculate energy associated with particle I after the move
            CALL POTENTIAL(ENRGN)

            IF(PINT .AND. COLLDT) THEN
                CALL COLLDENSFIELD(TRQNJ,QJN,INDXP)
            ELSEIF( HALFST .OR. EQUISEEDT ) THEN
                CALL HARMONIC_TRAP(BIASN)
            ENDIF

        ENDIF
    !   #####################################################################################################################

        IF(OVERLAPT .OR. ISNAN(ENRGN)) THEN
            REJECTT = .TRUE.
        ELSE
            IF( CLUSTERT .AND. CLSTRSZO /= CLSTRSZ) THEN
                REJECTT = .TRUE.
            ELSE
                DELE = ENRGN - ENRGO
                IF(.NOT. CLUSTERT) WSN  = VIRTEMP
                
                IF(PINT .AND. COLLDT) THEN
                    RQN    = RQ + (QJN - QJO)
                    TRQN   = ABS(RQN)/SQRT(REAL(NPART,DP))
                    DELQR  = TRQN - TRQ
                    DELQR2 = TRQN**2 - TRQ**2
                    UQR    = PINKAP2*DELQR2 - PINKAPA*DELQR
                    DELE   = DELE + UQR
                ELSEIF( HALFST .OR. EQUISEEDT ) THEN
                    UQR    = BIASN - BIASO
                    DELE   = DELE + UQR
                ENDIF

            !   Metropolis acceptence criteria
                BLTZMN = EXP( -DELE*BETAKB )
                IF (DRAND48() .GE. BLTZMN) THEN
                    REJECTT = .TRUE.
                ENDIF
            ENDIF
        ENDIF

        IF (REJECTT) THEN
            IF(RANDMOVE < 0.5_dp) THEN
                IF(CLUSTERT) THEN
                    DO J1 = 1, CLSTRSZO
                        J2 = CLSTRO(J1)
                        R(:,J2) = ROLD(:,J2)
                    !   Restore the position of the particle in the cell list
                        IF(CELLLISTT) THEN
                            CI  = C_INDEX ( R(:,J2)/BOX ) ! NEW CELL INDEX
                            CALL MOVE_IN_LIST ( J2, CI )
                        ENDIF
                    ENDDO
                ELSE
                    R(:,INDXP) = RO
                !   Restore the position of the particle in the cell list
                    IF(CELLLISTT) THEN
                        CI  = C_INDEX ( R(:,INDXP)/BOX ) ! NEW CELL INDEX
                        CALL MOVE_IN_LIST ( INDXP, CI )
                    ENDIF
                ENDIF
            ELSE 
                IF(CLUSTERT) THEN
                    DO J1 = 1, CLSTRSZO
                        J2 = CLSTRO(J1)
                        R(:,J2) = ROLD(:,J2)
                        Q(:,J2) = QOLD(:,J2)
                    !   Update rigid body sites to previous positions
                        RM    = Q_TO_RM( Q(:,J2) )
                        IF(RACEMICT) RC1 = MOD(((J2-1)-MOD((J2-1),12))/12+1,2)
                        DO J3 = 1, NSITES
                            IF(RACEMICT) THEN
                                IF(RC1==1) THEN
                                    RBSITES(:,J3,J2) = MATMUL(RM,REFSITE(:,J3))
                                ELSE
                                    RBSITES(:,J3,J2) = MATMUL(RM,REFSITE2(:,J3))
                                ENDIF
                            ELSE
                                RBSITES(:,J3,J2) = MATMUL(RM,REFSITE(:,J3))
                            ENDIF
                        ENDDO
                    !   Restore the position of the particle in the cell list
                        IF(CELLLISTT) THEN
                            CI  = C_INDEX ( R(:,J2)/BOX ) ! NEW CELL INDEX
                            CALL MOVE_IN_LIST ( J2, CI )
                        ENDIF
                    ENDDO
                ELSE
                !   Restore previous orientation of the particle
                    Q(:,INDXP) = QO
                !   Update rigid body sites to previous positions
                    RM    = Q_TO_RM( Q(:,INDXP) )
                    IF(RACEMICT) RC1 = MOD(((INDXP-1)-MOD((INDXP-1),12))/12+1,2)
                    DO J1 = 1, NSITES
                        IF(RACEMICT) THEN
                            IF(RC1==1) THEN
                                RBSITES(:,J1,INDXP) = MATMUL(RM ,REFSITE(:,J1))
                            ELSE
                                RBSITES(:,J1,INDXP) = MATMUL(RM ,REFSITE2(:,J1))
                            ENDIF
                        ELSE
                            RBSITES(:,J1,INDXP) = MATMUL(RM ,REFSITE(:,J1))
                        ENDIF
                    ENDDO
                ENDIF
            ENDIF
            VIR = VIROLD
        ELSE
            
            IF(RANDMOVE < 0.5_dp) THEN
                IF(CLUSTERT) THEN
                    ACCPTCTC = ACCPTCTC + 1
                ELSE
                    ACCPTCT = ACCPTCT + 1
                ENDIF
            ELSE
                IF(CLUSTERT) THEN
                    ACCPTCRC = ACCPTCRC + 1
                ELSE
                    ACCPTCR = ACCPTCR + 1
                ENDIF
            ENDIF
            
            PE  = PE + DELE
            VIR = VIROLD + (WSN - WSO)
        
            IF(PINT .AND. COLLDT) THEN
                TRQ = TRQN
                RQ  = RQN
                PE  = PE - UQR
            ELSEIF( HALFST .OR. EQUISEEDT ) THEN
                PE  = PE - UQR
            ENDIF
            
        ENDIF

        SPET     = .FALSE.
        CLUSTERT = .FALSE.

    END SUBROUTINE MV_MC

    SUBROUTINE VMOVE()
!   ================================================================================================
!   This routines performs a volume scaling move for Monte Carlo simulations in the NPT ensemble.
!   The simulation cell to use this routine must either be cubic or orthorhombic.
!   ================================================================================================
        USE COMMONS, ONLY: DP, CDP, BOX, PRSFIX, VLM, R, BETAKB, NPART, NVMOVES, ACCPTV, PE
        USE COMMONS, ONLY: VIR, VIRTEMP, MAXBOX, ISOTROPICT, NPZTT, NDIM, CELLLISTT, SPET, OVERLAPT, SCALED_R, RCUT
        USE COMMONS, ONLY: NCLSTRS, CLSTRADJ, CLUSTERMOVET, VLMCLUSTERMOVET
        USE COMMONS, ONLY: PINT, PINKAP2, PINA, PINDIM, TRQ, COLLDT
        USE CELL_LIST
        USE CLUSTER_MOVE, ONLY: VLM_CLUSTERMOVE
     
        IMPLICIT NONE

        INTEGER          :: J1, J2, BDIM, OLDADJ(NPART,NPART)
        REAL(KIND=DP)    :: VLMNEW, VLMOLD, LNBL
        REAL(KIND=DP)    :: PEO, PEN, VIRO, VIRN, DELB, ARG, DELPE
        REAL(KIND=DP)    :: BOXO(NDIM), ROLD(NDIM, NPART), RANDDIM, BOXDIFF
        REAL(KIND=CDP)   :: DRAND48
        REAL(KIND=DP)    :: TRQO, TRQN, UPINN, UPINO
        COMPLEX(KIND=DP) :: RQO, RQN
        LOGICAL          :: REJECTT
        
        NVMOVES = NVMOVES + 1

        SPET = .FALSE.; REJECTT = .FALSE.; OVERLAPT = .FALSE.
        PE = 0.0_dp; PEO = 0.0_dp; PEN = 0.0_dp
        VIRO = 0.0_dp; VIRN = 0.0_dp

        IF(CLUSTERMOVET) THEN
            VLMCLUSTERMOVET = .TRUE.
            CLSTRADJ = 0      ! Adjacency matrix for the system (bonding criterion is system specific)
        ENDIF

    !   Calculate the energy before the volume move.
        CALL POTENTIAL(PEO)

    !   Save variables in-case the move is rejected.
        VIRO   = VIRTEMP
        VLMOLD = VLM
        ROLD   = R
     
    !   PERFORM VOLUME MOVE AND SCALE THE POSITION OF THE PARTICLES ACCORDINGLY.
        IF(ISOTROPICT) THEN
    !   Random Walk in log(V), all box lengths remain equal
            BOXO   = BOX
            DELB   = (DRAND48() - 0.5_dp) * MAXBOX
            VLMNEW = EXP( LOG(VLMOLD) + DELB )
            BOX    = VLMNEW**(1.0_dp/3.0_dp)
        
            IF(CLUSTERMOVET) THEN
                CALL VLM_CLUSTERMOVE(.TRUE., BOXDIFF=BOX/BOXO)
            ELSE
                DO J1 = 1, NPART
                    R(:,J1) = R(:,J1) * BOX / BOXO
                ENDDO
            ENDIF
        ELSE
    !   Random Walk in log(V), box lengths are allowed to fluctuate independently
    !   Steps in volume are taken by randomly perturbing one of the box lengths.
            BOXO = BOX

            IF(NPZTT) THEN
        !   In the NPzT ensemble, so the cell lengths along the X and Y axes remain fixed.
                BDIM = PINDIM
            ELSE
                RANDDIM = DRAND48()
                IF(RANDDIM < 1.0_dp/3.0_dp) THEN
                    BDIM = 1
                ELSEIF(RANDDIM < 2.0_dp/3.0_dp) THEN
                    BDIM = 2
                ELSE
                    BDIM = 3
                ENDIF
            ENDIF

            LNBL     = (DRAND48() - 0.5_dp)*MAXBOX
            BOX(BDIM) = BOX(BDIM)*EXP(LNBL)
            BOXDIFF  = BOX(BDIM) / BOXO(BDIM)
            VLMNEW   = VLMOLD*BOXDIFF
            
            IF(CLUSTERMOVET) THEN
                CALL VLM_CLUSTERMOVE(.FALSE., BDIM=BDIM, BOXD=BOXDIFF)
            ELSE
                DO J1 = 1, NPART
                    R(BDIM,J1) = R(BDIM,J1) * BOXDIFF
                ENDDO
            ENDIF

        ENDIF 
     
        IF(CELLLISTT) THEN
            CALL FINALIZE_LIST()
            CALL INITIALIZE_LIST( NPART, RCUT/BOX )
            SCALED_R = 0.0_dp
            DO J1 = 1, NPART
                SCALED_R(:,J1) = R(:,J1)/BOX
            ENDDO
            CALL MAKE_LIST( NPART, SCALED_R )
        ENDIF

        IF(CLUSTERMOVET) THEN
            OLDADJ     = CLSTRADJ
            CLSTRADJ = 0
        ENDIF
     
    !   Calculate new energy
        CALL POTENTIAL(PEN)
        VIRN = VIRTEMP

    !   Check if the cluster move has altered the bonding within the system. If it has reject the move to
    !   ensure that detailed balanced is satisfied. Only need to check the upper triangular portion as an
    !   adjacency matrix is symmetric.
        IF(CLUSTERMOVET) THEN
            DO J1 = 1, NPART-1
                DO J2 = J1+1, NPART
                    IF(OLDADJ(J1,J2) /= CLSTRADJ(J1,J2)) REJECTT = .TRUE.
                ENDDO
            ENDDO
        ENDIF
        
    !   Check if the move should be rejected. For hard-core potentials, if two particles overlap the move is
    !   rejected. For soft-core potentials if two particles overlap enough the energy can reach a value which the 
    !   computer is unable to handle and so we also check for NaNs (this is only an issue for steeply repulsive potentials).
        IF(.NOT. REJECTT) THEN
            IF (OVERLAPT .OR. ISNAN(PEN)) THEN ! Rejected
                REJECTT = .TRUE.
            ELSE
        !   Calculate Boltzmann factor for isothermal-isobaric ensemble, if performing a cluster move the number of particles
        !   is replaced by the number of clusters.
                
                DELPE = PEN - PEO

                IF(CLUSTERMOVET) THEN
                    ARG = DELPE + ( PRSFIX * (VLMNEW - VLMOLD) )-( REAL((NCLSTRS+1),DP)*LOG(VLMNEW/VLMOLD)/BETAKB )
                ELSE
                    ARG = DELPE + ( PRSFIX * (VLMNEW - VLMOLD) )-( REAL((NPART+1),DP)*LOG(VLMNEW/VLMOLD)/BETAKB )
                ENDIF
                IF ( DRAND48() > EXP(-ARG*BETAKB) ) REJECTT = .TRUE.
            ENDIF
        ENDIF
        
        IF(REJECTT) THEN !   Reset parameters
            IF(ISOTROPICT) THEN
                BOX = BOXO
            ELSE
                BOX(BDIM) = BOXO(BDIM)
            ENDIF
     
            PE  = PEO
            VIR = VIRO
            VLM = VLMOLD
            R   = ROLD

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
           PE     = PEN     ! Volume change accepted
           VIR    = VIRN
           VLM    = VLMNEW
           ACCPTV = ACCPTV + 1
        ENDIF

        IF(CLUSTERMOVET) VLMCLUSTERMOVET = .FALSE.
     
    END SUBROUTINE VMOVE

    SUBROUTINE ADJUST()
!   ==============================================================================================
!   Subroutine to adjust the step sizes for the various Monte Carlo move sets in order to maintain
!   a preset acceptance ratio (usually takes a value of 30%-50%).
!   ==============================================================================================
        USE COMMONS, ONLY: DP, NPART, DMPFRQ, MAXDTR, MAXDRT, MAXBOX, ACCRATC, ACCRATV, NPTT, &
                           RIGIDT, PI, NVMOVES, NTMOVES, NRMOVES, ACCPTCT, ACCPTCR, ACCPTV, &
                           MAXDTRC, MAXDRTC, NTCMOVES, NRCMOVES, ACCPTCTC, ACCPTCRC, CLUSTERMOVET

        IMPLICIT NONE

        REAL(KIND=DP) :: RACCPCT, RACCPCR, RACCPV, RACCPCTC, RACCPCRC

        OPEN (UNIT=61, FILE='run.dat', STATUS='UNKNOWN', ACCESS='APPEND')
        IF(RIGIDT) THEN
            RACCPCT = REAL(ACCPTCT,DP) / REAL(NTMOVES,DP)
            RACCPCR = REAL(ACCPTCR,DP) / REAL(NRMOVES,DP)
        ELSE
            RACCPCT = REAL(ACCPTCT,DP) / REAL((DMPFRQ*NPART),DP)
        ENDIF

        IF(NPTT) RACCPV = REAL(ACCPTV,DP) / REAL(NVMOVES,DP)

        IF (RACCPCT < ACCRATC) THEN
            MAXDTR = MAXDTR*0.975_dp !Translational step
        ELSE
            IF(MAXDTR*1.025_dp <= 2.0_dp) MAXDTR = MAXDTR*1.025_dp
        ENDIF

        IF (RACCPCR < ACCRATC) THEN
            MAXDRT = MAXDRT*0.975_dp !Rotational step
        ELSE
            IF(MAXDRT*1.025 <= PI) MAXDRT = MAXDRT*1.025_dp
        ENDIF

        IF(CLUSTERMOVET) THEN
            RACCPCTC = REAL(ACCPTCTC,DP) / REAL(NTCMOVES,DP)
            RACCPCRC = REAL(ACCPTCRC,DP) / REAL(NRCMOVES,DP)
            
            IF (RACCPCTC < ACCRATC) THEN
                MAXDTRC = MAXDTRC*0.975_dp !Translational step
            ELSE
                IF(MAXDTRC*1.025 <= 2.0_dp) MAXDTRC = MAXDTRC*1.025
            ENDIF
    
            IF (RACCPCRC < ACCRATC) THEN
                MAXDRTC = MAXDRTC*0.975_dp !Rotational step
            ELSE
                IF(MAXDRTC*1.025_dp <= PI) MAXDRTC = MAXDRTC*1.025_dp
            ENDIF

        ENDIF

        IF (NPTT) THEN
            
            IF (RACCPV < ACCRATV) THEN
                MAXBOX = MAXBOX*0.98_dp
            ELSE
                MAXBOX = MAXBOX*1.02_dp
            ENDIF
            
            IF(RIGIDT) THEN
                WRITE (61, '(A3,1X,4F12.7)') "SPM", RACCPCT, RACCPCR, MAXDTR, MAXDRT
                WRITE (61, '(A3,1X,2F12.7)') "VLM", RACCPV, MAXBOX
            ELSE
                WRITE (61, '(A3,1X,2F12.7)') "SPM", RACCPCT, MAXDTR
                WRITE (61, '(A3,1X,2F12.7)') "VLM", RACCPV, MAXBOX
            ENDIF
        ELSE
            IF(RIGIDT) THEN
                WRITE (61, '(A3,1X,4F12.7)') "SPM", RACCPCT, RACCPCR, MAXDTR, MAXDRT
            ELSE
                WRITE (61, '(A3,1X,2F12.7)') "SPM", RACCPCT, MAXDTR
            ENDIF
        ENDIF
        IF(CLUSTERMOVET) WRITE (61, '(A3,1X,4F12.7)') "CLM", RACCPCTC, RACCPCRC, MAXDTRC, MAXDRTC

        CLOSE(UNIT=61, STATUS='KEEP')

        NTMOVES  = 0
        NRMOVES  = 0
        NVMOVES  = 0
        ACCPTCT  = 0
        ACCPTCR  = 0
        ACCPTV   = 0
        ACCPTCTC = 0
        ACCPTCRC = 0
        NTCMOVES = 0
        NRCMOVES = 0
        
    END SUBROUTINE ADJUST