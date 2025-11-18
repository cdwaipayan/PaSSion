

    SUBROUTINE MV_MC()

    !     This subroutine performs a single-particle and cluster Monte Carlo move.

        USE COMMONS, ONLY: DP, CDP, NDIM, NPART, R, Q, BETAKB, BOX, MAXDTR, MAXDRT, SPET, RACEMICT, REFSITE2, ONEBONDT, NPBONDS, &
        NTMOVES, NRMOVES, ACCPTCT, ACCPTCR, COLLDT, PINT, PINKAP2, PINKAPA, TRQ, CELLLISTT, RQ, RIGIDT, PE, N_POLY_TOT, N_POLY, &
        INDXP, OVERLAPT, NSITES, RBSITES, REFSITE, CLUSTERT, HALFST, EQUISEEDT, CLUSTERMOVET, CLSTRRATIO, POLYCHAINT, POLY_SIG2
        USE COMMONS, ONLY: CLUSTERT, CLSTR, CLSTRSZ, ACCPTCTC, ACCPTCRC, LRGCLSTRT, LRGCLSTRRATIO, LRGCLSTMVT, N_POLY_L
        USE COMMONS, ONLY: VIR, VIRTEMP, STRESST, STRESSTEMP, VIR_TENS, BYUKWAT, BNRYA, OBLATESPHYT, TRIANGLET
        USE COMMONS, ONLY: POLYMOVET, PLYMVRTIO, PLYRIGID, MCPLYMVT, OB_BONDST, OB_N_BNDS, OB_NEIGHS,OBJSURFT,SPHERECNFT,SPHERERAD
        USE CELL_LIST, ONLY: C_INDEX, MOVE_IN_LIST
        USE ROTATIONS_MODULE, ONLY: RANDOM_ROTATE_QUATERNION, Q_TO_RM
        USE CLUSTER_MOVE
        
        IMPLICIT NONE
        
        INTEGER             :: J1, RC1, INDXPLY, CHAINID
        REAL(KIND=DP)       :: RO(NDIM)
        REAL(KIND=DP)       :: WSN, WSO, ENRGN, ENRGO, DELE, BLTZMN, VIROLD, BIASO, BIASN
        REAL(KIND=DP)       :: STRESSOLD(NDIM,NDIM), STRSO(NDIM,NDIM), STRSN(NDIM,NDIM)
        REAL(KIND=DP)       :: RANDMOVE, TRQN, DELQR, DELQR2, UQR, TRQOJ, TRQNJ, PBNDO, PBNDN
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
        MCPLYMVT = .FALSE.

        ENRGN   = 0.0_dp ! New energy
        ENRGO   = 0.0_dp ! Old energy
        TRQNJ   = 0.0_dp 
        TRQOJ   = 0.0_dp
        QJO     = 0.0_dp
        QJN     = 0.0_dp
        VIROLD  = VIR
        IF(STRESST) STRESSOLD = VIR_TENS

        WSO = 0.0_dp
        WSN = 0.0_dp
        IF(STRESST) THEN
            STRSO = 0.0_dp
            STRSN = 0.0_dp
        ENDIF

        IF(POLYMOVET) THEN
            IF(DRAND48() < PLYMVRTIO) THEN ! 33_dp
                MCPLYMVT = .TRUE.
                IF(.NOT. OBJSURFT) THEN
                    IF(DRAND48() < PLYRIGID) THEN
                        IF(DRAND48()<0.5_dp) THEN
                            CALL POLYMER_TRANSLATION()
                        ELSE
                            CALL POLYMER_ROTATE()
                        ENDIF
                        RETURN
                    ENDIF
                ENDIF
                CALL POLYMER_RATCHET()
                RETURN 
            ENDIF
        ENDIF

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

        IF(POLYMOVET .AND. (OBLATESPHYT .OR. TRIANGLET)) THEN
            IF(DRAND48()<0.5_dp) THEN
            !   Randomly select a polymer chain for which the move will be attempted
                INDXPLY = INT(N_POLY*DRAND48())
            !   Retrive the particle ID of the first particle in the chain.
                CHAINID = INDXPLY*N_POLY_L + 1
            !   Randomly select a particle in the chain (excluding the first and last beads) to act as the 
            !   pivot particle for the move.
                INDXP = CHAINID + INT((N_POLY_L-1)*DRAND48()) + 1
            ELSE
            !   Randomly select a capsomer to move
                INDXP = INT((NPART-N_POLY_TOT)*DRAND48()) + 1 + N_POLY_TOT
            ENDIF
        ELSE
            INDXP   = INT(NPART*DRAND48()) + 1
        ENDIF

        IF(RIGIDT) THEN
            IF(POLYCHAINT .AND. INDXP <= N_POLY_TOT) THEN
                RANDMOVE = 0.0_dp
            ELSEIF(BYUKWAT .AND. INDXP>BNRYA) THEN
                RANDMOVE = 0.0_dp
            ELSE
                RANDMOVE = DRAND48()
            ENDIF
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
            IF(STRESST) STRSO = STRESSTEMP

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

            IF(ONEBONDT) THEN
                PBNDO = 0.0_dp
                DO J3 = 1, NSITES
                    IF(NPBONDS(J3)>1) THEN
                        PBNDO = PBNDO + REAL(NPBONDS(J3),DP) - 1.0_dp
                    ENDIF
                ENDDO
            ENDIF 

        !   Perform random move
            IF( RANDMOVE < 0.5_dp ) THEN
            !   Translational moves
                DO J1 = 1, NDIM
                    RO(J1)      = R(J1,INDXP)
                    IF(POLYCHAINT .AND. INDXP <= N_POLY_TOT) THEN
                        R(J1,INDXP) = R(J1,INDXP) + (2.0_dp*DRAND48()-1.0_dp)*MAXDTR*POLY_SIG2
                    ELSE
                        R(J1,INDXP) = R(J1,INDXP) + (2.0_dp*DRAND48()-1.0_dp)*MAXDTR
                    ENDIF
            !   Pick up the central image
                    IF((OBJSURFT .AND. J1 /= NDIM).OR.(.NOT. OBJSURFT).OR.(.NOT. SPHERECNFT)) THEN
                        R(J1,INDXP) = R(J1,INDXP) - BOX(J1)*ANINT(R(J1,INDXP)/BOX(J1))
                    ENDIF
                ENDDO

                IF(SPHERECNFT) THEN
                    IF(NORM2(R(:,INDXP))+0.5_dp > SPHERERAD) REJECTT = .TRUE.
                ENDIF

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

            IF(ONEBONDT) THEN
                PBNDN = 0.0_dp
                DO J3 = 1, NSITES
                    IF(NPBONDS(J3)>1) THEN
                        PBNDN = PBNDN + REAL(NPBONDS(J3),DP) - 1.0_dp
                    ENDIF
                ENDDO
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
                IF(.NOT. CLUSTERT) THEN
                    WSN  = VIRTEMP
                    IF(STRESST) STRSN = STRESSTEMP
                ENDIF
                
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

                IF(ONEBONDT) DELE = DELE + (PBNDN-PBNDO)

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
            IF(STRESST) VIR_TENS = STRESSOLD

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
            IF(STRESST) VIR_TENS = STRESSOLD + (STRSN-STRSO)
        
            IF(PINT .AND. COLLDT) THEN
                TRQ = TRQN
                RQ  = RQN
                PE  = PE - UQR
            ELSEIF( HALFST .OR. EQUISEEDT ) THEN
                PE  = PE - UQR
            ENDIF

            IF(ONEBONDT) PE = PE - (PBNDN-PBNDO)
            
        ENDIF

        SPET     = .FALSE.
        CLUSTERT = .FALSE.
        MCPLYMVT = .FALSE.
        OVERLAPT = .FALSE.

    END SUBROUTINE MV_MC

    SUBROUTINE VMOVE()
!   ================================================================================================
!   This routines performs a volume scaling move for Monte Carlo simulations in the NPT ensemble.
!   The simulation cell to use this routine must either be cubic or orthorhombic.
!   ================================================================================================
        USE COMMONS, ONLY: DP, CDP, BOX, PRSFIX, RHO, VLM, R, BETAKB, NPART, NVMOVES, ACCPTV, PE, PRNTCNF
        USE COMMONS, ONLY: VIR, VIRTEMP, MAXBOX, ISOTROPICT, NPZTT, NDIM, CELLLISTT, SPET, OVERLAPT, SCALED_R, RCUT
        USE COMMONS, ONLY: NCLSTRS, CLSTRADJ, CLUSTERMOVET, VLMCLUSTERMOVET, UMBRELLAT, TRGTR, RPINK, PINDIM
        USE COMMONS, ONLY: STRESST, STRAIN_DIM, STRESSTEMP, STRAINT, VIR_TENS
        USE COMMONS, ONLY: PI, SPHERECNFT, SPHERERAD, SPHERERAD2, OBJSURFT
        ! USE COMMONS, ONLY: PINT, PINKAP2, PINA, TRQ, COLLDT
        USE CELL_LIST
        USE CLUSTER_MOVE, ONLY: VLM_CLUSTERMOVE
     
        IMPLICIT NONE

        INTEGER            :: J1, J2, BDIM, OLDADJ(NPART,NPART)
        REAL(KIND=DP)      :: VLMNEW, VLMOLD, RHONEW, RHOOLD, LNBL
        REAL(KIND=DP)      :: PEO, PEN, VIRO, VIRN, DELB, ARG, DELPE
        REAL(KIND=DP)      :: STRESSO(NDIM,NDIM), STRESSN(NDIM,NDIM)
        REAL(KIND=DP)      :: BOXO(NDIM), ROLD(NDIM, NPART), RANDDIM, BOXDIFF
        REAL(KIND=DP)      :: RADO, RSOLD(NDIM, NPART)
        REAL(KIND=CDP)     :: DRAND48
        ! REAL(KIND=DP)    :: TRQO, TRQN, UPINN, UPINO
        REAL(KIND=DP)      :: BIASOLD, BIASNEW
        ! COMPLEX(KIND=DP) :: RQO, RQN
        LOGICAL            :: REJECTT
        
        NVMOVES = NVMOVES + 1

        SPET = .FALSE.; REJECTT = .FALSE.; OVERLAPT = .FALSE.
        PE = 0.0_dp; PEO = 0.0_dp; PEN = 0.0_dp
        VIRO = 0.0_dp; VIRN = 0.0_dp
        IF(STRESST) THEN
            STRESSO = 0.0_dp; STRESSN = 0.0_dp
        ENDIF

        IF(CLUSTERMOVET) THEN
            VLMCLUSTERMOVET = .TRUE.
            CLSTRADJ = 0      ! Adjacency matrix for the system (bonding criterion is system specific)
        ENDIF

    !   Calculate the energy before the volume move.
        CALL POTENTIAL(PEO)

    !   Save variables in-case the move is rejected.
        VIRO   = VIRTEMP
        IF(STRESST) STRESSO = STRESSTEMP
        VLMOLD = VLM
        ROLD   = R
        IF(UMBRELLAT) RHOOLD = RHO

        IF(SPHERECNFT) THEN
            RADO = SPHERERAD
            BOXO = BOX
            DO J1 = 1, NPART
                RSOLD(1,J1) = NORM2(R(:,J1))
                RSOLD(2,J1) = ATAN2( R(2,J1),R(1,J1) )
                RSOLD(3,J1) = ACOS( R(3,J1)/RSOLD(1,J1) )
            ENDDO
            DELB       = (DRAND48() - 0.5_dp) * MAXBOX
            VLMNEW     = EXP( LOG(VLMOLD) + DELB )
            SPHERERAD  = (3.0_dp*VLMNEW/(4.0_dp*PI))**(1.0_dp/3.0_dp)
            RSOLD(1,:) = RSOLD(1,:)*SPHERERAD/RADO
            DO J1 = 1, NPART
                R(1,J1) = RSOLD(1,J1)*COS(RSOLD(2,J1))*SIN(RSOLD(3,J1))
                R(2,J1) = RSOLD(1,J1)*SIN(RSOLD(2,J1))*SIN(RSOLD(3,J1))
                R(3,J1) = RSOLD(1,J1)*COS(RSOLD(3,J1))
            ENDDO
            BOX = SPHERERAD*2.0_dp

        ELSE
     
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
                ELSEIF(OBJSURFT) THEN
                    IF(RANDDIM < 0.5_dp) THEN
                        BDIM = 1
                    ELSE
                        BDIM = 2
                    ENDIF
                ELSEIF(STRAINT) THEN
            !   When performing a stress-strain simulation allow the box lengths perpendicular 
            !   to the axis along which stress is being appplied to fluctuate in order to keep 
            !   the pressure of the system constant.
                    BDIM = STRAIN_DIM
                    DO WHILE(BDIM==STRAIN_DIM)
                        RANDDIM = DRAND48()
                        IF(RANDDIM < 1.0_dp/3.0_dp) THEN
                            BDIM = 1
                        ELSEIF(RANDDIM < 2.0_dp/3.0_dp) THEN
                            BDIM = 2
                        ELSE
                            BDIM = 3
                        ENDIF
                    ENDDO
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
        IF(STRESST) STRESSN = STRESSTEMP

    !   Check if the cluster move has altered the bonding within the system. If it has reject the move to
    !   ensure that detailed balanced is satisfied. Only need to check the upper triangular portion as an
    !   adjacency matrix is symmetric.
        IF(CLUSTERMOVET) THEN
            DO J1 = 1, NPART-1
                DO J2 = J1+1, NPART
                    IF(PRNTCNF) THEN
                        IF(OLDADJ(J1,J2) /= CLSTRADJ(J1,J2)) THEN
                            IF(OLDADJ(J1,J2)==1 .OR. CLSTRADJ(J1,J2)==1) THEN
                                REJECTT = .TRUE.
                            ENDIF
                        ENDIF
                    ELSE
                        IF(OLDADJ(J1,J2) /= CLSTRADJ(J1,J2)) REJECTT = .TRUE.
                    ENDIF
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
        !   Calculate difference in energy between the old and new configuration
                DELPE = PEN - PEO
        !   If performing umbrella sampling with density as an order-parameter, calculate the bias potentials
                IF(UMBRELLAT .AND. RPINK>0.0_dp) THEN
                    BIASOLD = RPINK *(RHOOLD - TRGTR)**2
                    RHONEW  = REAL(NPART,DP) / VLMNEW
                    BIASNEW = RPINK *(RHONEW - TRGTR)**2
                ENDIF
        !   Calculate Boltzmann factor for isothermal-isobaric ensemble, if performing a cluster move the number of particles
        !   is replaced by the number of clusters.
                IF(CLUSTERMOVET) THEN
                    ARG = DELPE + ( PRSFIX * (VLMNEW - VLMOLD) )-( REAL((NCLSTRS+1),DP)*LOG(VLMNEW/VLMOLD)/BETAKB )
                ELSE
                    ARG = DELPE + ( PRSFIX * (VLMNEW - VLMOLD) )-( REAL((NPART+1),DP)*LOG(VLMNEW/VLMOLD)/BETAKB )
                ENDIF
                IF(UMBRELLAT .AND. RPINK>0.0_dp) ARG = ARG + (BIASNEW - BIASOLD)
                IF ( DRAND48() > EXP(-ARG*BETAKB) ) REJECTT = .TRUE.
            ENDIF
        ENDIF
        
        IF(REJECTT) THEN !   Reset parameters
            IF(ISOTROPICT) THEN
                BOX = BOXO
            ELSE
                BOX(BDIM) = BOXO(BDIM)
            ENDIF

            IF(SPHERECNFT) SPHERERAD = RADO
     
            PE  = PEO
            VIR = VIRO
            IF(STRESST) VIR_TENS = STRESSO
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
           IF(STRESST) VIR_TENS = STRESSN
           VLM    = VLMNEW
           IF(SPHERECNFT) SPHERERAD2 = SPHERERAD*SPHERERAD
           ACCPTV = ACCPTV + 1
        ENDIF

        IF(CLUSTERMOVET) VLMCLUSTERMOVET = .FALSE.
     
    END SUBROUTINE VMOVE

    SUBROUTINE APPLY_STRAIN()
        
        USE COMMONS, ONLY: DP, BOX, R, NPART, VLM
        USE COMMONS, ONLY: CELLLISTT, SCALED_R, RCUT, STRAINT, NPTT
        USE COMMONS, ONLY: STRAIN_DIM, TOT_STRAIN, MAX_STRAIN, STRAIN_DEL
        USE CELL_LIST
     
        IMPLICIT NONE

        INTEGER            :: J1
        REAL(KIND=DP)      :: BOXO
        LOGICAL            :: APPLY_STRN

        APPLY_STRN = .FALSE.

        IF(MAX_STRAIN>0.0_dp) THEN
            IF(TOT_STRAIN<MAX_STRAIN) APPLY_STRN = .TRUE.
        ELSE
            IF(TOT_STRAIN>MAX_STRAIN) APPLY_STRN = .TRUE.
        ENDIF

        IF(APPLY_STRN) THEN
        !   Update box length along axis where strain is being applied
            BOXO = BOX(STRAIN_DIM)
            BOX(STRAIN_DIM) = BOX(STRAIN_DIM) + BOX(STRAIN_DIM)*STRAIN_DEL
        !   Update the volume of the system following application of strain
            VLM = VLM*BOX(STRAIN_DIM)/BOXO
        !   Update the positions of the particles accordingly
            DO J1 = 1, NPART
                R(STRAIN_DIM,J1) = R(STRAIN_DIM,J1) * BOX(STRAIN_DIM)/BOXO
            ENDDO
        !   Update cell-list (if being used)
            IF(CELLLISTT) THEN
                CALL FINALIZE_LIST()
                CALL INITIALIZE_LIST( NPART, RCUT/BOX )
                SCALED_R = 0.0_dp
                DO J1 = 1, NPART
                    SCALED_R(:,J1) = R(:,J1)/BOX
                ENDDO
                CALL MAKE_LIST( NPART, SCALED_R )
            ENDIF
        !   Update the count of the total strain applied to the system
            TOT_STRAIN = TOT_STRAIN + STRAIN_DEL
            PRINT *, TOT_STRAIN, BOX
        ! ELSE
        !     STRAINT = .FALSE.
        !     NPTT    = .FALSE.
        ENDIF

    END SUBROUTINE




    SUBROUTINE ADJUST()
!   ==============================================================================================
!   Subroutine to adjust the step sizes for the various Monte Carlo move sets in order to maintain
!   a preset acceptance ratio (usually takes a value of 30%-50%).
!   ==============================================================================================
        USE COMMONS, ONLY: DP, NPART, DMPFRQ, MAXDTR, MAXDRT, MAXBOX, ACCRATC, ACCRATV, NPTT, POLYCHAINT, MAXDPTR, MAXDPRT, &
                            RIGIDT, PI, NVMOVES, NTMOVES, NRMOVES, ACCPTCT, ACCPTCR, ACCPTV, NTPMOVES, NRPMOVES, &
                            MAXDTRC, MAXDRTC, NTCMOVES, NRCMOVES, ACCPTCTC, ACCPTCRC, CLUSTERMOVET, ACCPTCTCP, ACCPTCRCP, &
                            OBLATESPHYT, MAXDPRATT, NRATCPMOVES, ACCPTCRATCP

        IMPLICIT NONE

        REAL(KIND=DP) :: RACCPCT, RACCPCR, RACCPV, RACCPCTC, RACCPCRC, RACCPCRATC, RACCPCTP, RACCPCRP

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
        
        IF(POLYCHAINT) THEN
            RACCPCTP   = REAL(ACCPTCTCP,DP) / REAL(NTPMOVES,DP)
            RACCPCRP   = REAL(ACCPTCRCP,DP) / REAL(NRPMOVES,DP)
            RACCPCRATC = REAL(ACCPTCRATCP,DP) / REAL(NRATCPMOVES,DP)
            
            IF (RACCPCTP < ACCRATC) THEN
                MAXDPTR = MAXDPTR*0.975_dp !Translational step
            ELSE
                IF(MAXDPTR*1.025 <= 2.0_dp) MAXDPTR = MAXDPTR*1.025
            ENDIF
    
            IF (RACCPCRP < ACCRATC) THEN
                MAXDPRT = MAXDPRT*0.975_dp !Rotational step
            ELSE
                IF(MAXDPRT*1.025_dp <= PI) MAXDPRT = MAXDPRT*1.025_dp
            ENDIF

            IF (RACCPCRATC < ACCRATC) THEN
                MAXDPRATT = MAXDPRATT*0.975_dp !Rotational step
            ELSE
                IF(MAXDPRATT*1.025_dp <= PI) MAXDPRATT = MAXDPRATT*1.025_dp
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
        IF((POLYCHAINT)) WRITE (61, '(A3,1X,6F12.7)') "PLY", RACCPCTP, RACCPCRP, RACCPCRATC, MAXDPTR, MAXDPRT, MAXDPRATT

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
        ACCPTCTCP= 0 
        ACCPTCRCP= 0
        NTPMOVES = 0
        NRPMOVES = 0
        ACCPTCRATCP = 0
        NRATCPMOVES = 0
        
    END SUBROUTINE ADJUST