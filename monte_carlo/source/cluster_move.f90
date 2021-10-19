MODULE CLUSTER_MOVE
    IMPLICIT NONE
CONTAINS    
!     ==============================================================================================
    SUBROUTINE CLUSTERENERGY(RANDMOVE, CLSTR_ENERGY, WS, QJ, TRQJ, BIAS)
    !   This subroutine performs a cluster-move Monte Carlo cycle.
    !   To satisfy detailed balance, if the move results in a change in the cluster, the move is rejected.

        USE COMMONS, ONLY: DP, NDIM, R, BOX, VIRTEMP, INDXP, OVERLAPT, PINT, HALFST, EQUISEEDT, COLLDT
        USE COMMONS, ONLY: CLSTR, CLSTRSZ, CLSTRID, RCLSTR

        USE CELL_LIST, ONLY: C_INDEX, MOVE_IN_LIST

        USE ROTATIONS_MODULE, ONLY: RANDOM_ROTATE_QUATERNION, Q_TO_RM
        
        IMPLICIT NONE

        REAL(KIND=DP), INTENT(IN)               :: RANDMOVE
        
        INTEGER             :: J1, CLSTRCNT
        REAL(KIND=DP)       :: RIJ(NDIM), ECLSTR, BIASCLSTR, TRQCLSTR
        COMPLEX(KIND=DP)    :: QJCLSTR

        REAL(KIND=DP), INTENT(OUT)              :: CLSTR_ENERGY 
        REAL(KIND=DP), INTENT(OUT), OPTIONAL    :: WS, TRQJ, BIAS
        COMPLEX(KIND=DP), INTENT(OUT), OPTIONAL :: QJ

    !   Initialise the energy of the cluster to be 0
        CLSTR_ENERGY = 0.0_dp
        WS           = 0.0_dp
    !   Initialise the cluster size to be 1.0_dp
        CLSTRSZ = 1
    !   Counter to keep track of where we are in the list of particles in the cluster.
        CLSTRID = 1
    !   Array containing IDs of particles in the cluster is initialised
    !   as an array of zeros.
        CLSTR  = 0
    !   Set the 1st entry of the cluster ID array to be the seed particle.
        CLSTR(CLSTRID)  = INDXP
        RCLSTR(:,INDXP) = R(:,INDXP)

    !   Initialise variables used for interface-pinning simulations
        IF(PINT .AND. COLLDT) THEN
        !   Value of OP for particles in the cluster before/after the move
            TRQJ = 0.0_dp 
        ! Value of complex OP for particles in the cluster before/after the move
            QJ   = 0.0_dp 
        ELSEIF( HALFST .OR. EQUISEEDT ) THEN
            BIAS = 0.0_dp
        ENDIF

    !   Calculate potential energy associated with each particle in the cluster.
    !   The cluster is built iteratively, with the list extended each time the 
    !   the potential energy is calculated.
        DO
        !   Keep track of the size of the cluster at the beginning of the iteration so we can figure out by 
        !   how much the cluster has groen from iteration to iteration.
            CLSTRCNT = CLSTRSZ
        !   Calculate the pair energies associated with the current index particle and add 
        !   new particles to the cluster (this is currently done separately in each of the 
        !   pair potential subroutines).
            CALL POTENTIAL(ECLSTR)
            CLSTR_ENERGY = CLSTR_ENERGY + ECLSTR
            WS           = WS + VIRTEMP
        !   If there are particle overlaps there is no need to proceed further as we reject the move.
        !   Check for NaN which may arise due to overlapping particles interacting via a 
        !   repulsive Yukawa/WCA/Kihara potential with a steep gradient.
            IF( OVERLAPT .OR. ISNAN(CLSTR_ENERGY)) RETURN
        
            IF(PINT .AND. COLLDT) THEN
        !   If performing an interface-pinning simulation calculate the collective
        !   density field for the particles in the cluster. 
                CALL COLLDENSFIELD(TRQCLSTR,QJCLSTR,INDXP)
                TRQJ = TRQJ + TRQCLSTR
                QJ   = QJ   + QJCLSTR
            ELSEIF( HALFST .OR. EQUISEEDT ) THEN
        !   If performing an equilibration of the interface between two phases
        !   (bulk solid or solid seed + bulk liquid) calculate the harmonic potential
        !   energy associated with "trapped" particles. 
                CALL HARMONIC_TRAP(BIASCLSTR)
                BIAS = BIAS + BIASCLSTR
            ENDIF

        !   Save the positions of the particles added to the cluster relative to the current seed particle.
        !   This ensures that a single network will be moved when performing a rotational move.
            IF(RANDMOVE >= 0.5_dp) THEN
                DO J1 = CLSTRCNT+1, CLSTRSZ
                    RIJ = RCLSTR(:,CLSTR(CLSTRID)) - R(:,CLSTR(J1))
                    RIJ = RIJ - BOX*ANINT(RIJ/BOX)
                    RCLSTR(:,CLSTR(J1)) = RCLSTR(:,CLSTR(CLSTRID)) - RIJ
                ENDDO   
            ENDIF
            
        !   Update counter for the next iteration
            CLSTRID = CLSTRID + 1
        !   If the current particle index is equal the last entry in the cluster
        !   list then we have no new particles in the cluster, so we exit the loop.   
            IF(INDXP == CLSTR(CLSTRSZ)) THEN
                EXIT
            ELSE
                INDXP = CLSTR(CLSTRID)
            ENDIF
        ENDDO

    END SUBROUTINE

    SUBROUTINE CLUSTERMOVE (RANDMOVE, ROLD, QOLD)

        USE COMMONS, ONLY: DP, CDP, NPART, NDIM, MAXDTRC, MAXDRTC, NTCMOVES, NRCMOVES
        USE COMMONS, ONLY: R, Q, BOX, CELLLISTT, NSITES, RBSITES, REFSITE, CLSTR, CLSTRSZ, RCLSTR
        USE CELL_LIST, ONLY: C_INDEX, MOVE_IN_LIST
        USE ROTATIONS_MODULE, ONLY: RANDOM_VECTOR, QUATMUL, Q_TO_RM, RANDOM_ROTATE_QUATERNION
        
        IMPLICIT NONE

        REAL(KIND=DP), INTENT(IN)    :: RANDMOVE
        REAL(KIND=DP), INTENT(OUT)   :: ROLD(NDIM,NPART), QOLD(4,NPART)

        INTEGER       :: J1, J2, J3, CI(3), SEED
        REAL(KIND=DP) :: TMOVEMAP(NDIM), RMOVEMAP(4)
        REAL(KIND=DP) :: RM(3,3), RRM(3,3), AXIS(NDIM), ANGLE, RCOM(NDIM)
        REAL(KIND=CDP):: DRAND48

        ROLD = 0.0_dp
        QOLD = 0.0_dp

        IF(RANDMOVE < 0.5_dp) THEN
        !   Translational moves
            DO J1 = 1, NDIM
                TMOVEMAP(J1) = (2.0_dp*DRAND48()-1.0_dp)*MAXDTRC
            ENDDO

            DO J1 = 1, CLSTRSZ
                J2 = CLSTR(J1)
                IF(J2 == 0) EXIT

            !   Save old positions
                ROLD(:,J2) = R(:,J2)
            !   Translate particles in the cluster
                R(:,J2) = R(:,J2) + TMOVEMAP
            !   Pick up the central image
                R(:,J2) = R(:,J2) - BOX*ANINT(R(:,J2)/BOX)

            !   Update the position of the particle in the cell list following the move
                IF(CELLLISTT) THEN
                    CI  = C_INDEX ( R(:,J2)/BOX ) ! NEW CELL INDEX
                    CALL MOVE_IN_LIST ( J2, CI )
                ENDIF
            ENDDO

            NTCMOVES = NTCMOVES + 1

        ELSE
        !   Rotational moves
            AXIS = RANDOM_VECTOR( )                             ! Choose random unit vector
            ANGLE = ( 2.0_dp*DRAND48() - 1.0_dp ) * MAXDRTC     ! Uniform random angle in desired range

            ! Standard formula for rotation quaternion, using half angles
            RMOVEMAP(1)   = COS(0.5_dp*ANGLE)
            RMOVEMAP(2:4) = SIN(0.5_dp*ANGLE)*AXIS

        !   Extract the rotation matrix associated with the the rotational move
            RRM = Q_TO_RM(RMOVEMAP)
            
        !   Calculate the center of mass of the cluster, using the seed particle as 
        !   the reference position.
            SEED = CLSTR(1)
            RCOM = RCLSTR(:,SEED)
        !   Save old position of the seed particle
            ROLD(:,SEED) = R(:,SEED)
            DO J1 = 2, CLSTRSZ
                J2 = CLSTR(J1)
                IF(J2 == 0) EXIT ! We have exhausted our list of particles
            !   Save old positions
                ROLD(:,J2) = R(:,J2)
            !   Extract the position of each particle in the cluster
                R(:,J2) = RCLSTR(:,J2)
            !   Update the center of mass sum
                RCOM = RCOM + R(:,J2)
            ENDDO

        !   Normalise the sum to get the center of mass
            RCOM = RCOM / REAL(CLSTRSZ,DP)

            DO J1 = 1, CLSTRSZ
                J2 = CLSTR(J1)
                IF(J2 == 0) EXIT
            !   Rotate particle around the center of mass of the cluster
                R(:,J2) = R(:,J2) - RCOM
                R(:,J2) = MATMUL(RRM,R(:,J2)) + RCOM
            !   Pick up the central image
                R(:,J2) = R(:,J2) - BOX*ANINT(R(:,J2)/BOX)
            
            !   Save old orientations
                QOLD(:,J2) = Q(:,J2)
            !   Update the orientation of the particle
                Q(:,J2) = QUATMUL(RMOVEMAP, QOLD(:,J2))
            !   Update the rigid body sites of the particle being displaced
                RM   = Q_TO_RM( Q(:,J2) )
                DO J3 = 1, NSITES
                    RBSITES(:,J3,J2) = MATMUL(RM,REFSITE(:,J3))
                ENDDO
            !   Update the position of the particle in the cell list following the move
                IF(CELLLISTT) THEN
                    CI  = C_INDEX ( R(:,J2)/BOX ) ! NEW CELL INDEX
                    CALL MOVE_IN_LIST ( J2, CI )
                ENDIF
            ENDDO

            NRCMOVES    = NRCMOVES + 1
        ENDIF

    END SUBROUTINE

    SUBROUTINE BUILD_ALL_CLUSTERS(BOPSNT)

        USE COMMONS, ONLY: DP, R, NDIM, NPART, BOX, NCLSTRS, CLSTRADJ, CLSTRCOM, PCLSTR, PCLSTRID, CLURIJ, NPTT
        IMPLICIT NONE

        INTEGER                         :: J1, J2, CLSTRSZ(NPART), ICLU, PREF(NPART)
        REAL(KIND=DP)                   :: RIJ(NDIM)
        LOGICAL, OPTIONAL, VALUE        :: BOPSNT

        IF(.NOT. PRESENT(BOPSNT)) BOPSNT = .FALSE.

        NCLSTRS  = 0      ! Number of unique clusters in the system
        PCLSTRID = 0      ! Cluster ID associated with each particle.
        PREF     = 0      ! ID of the seed particle associated with each cluster 
        CLSTRSZ  = 0      ! Size of each of the clusters
        CLSTRCOM = 0.0_dp ! Center-of-mass for each cluster
        PCLSTR   = 0.0_dp ! Relative position of each particle to the center-of-mass of its associated cluster.

        DO J1 = 1, NPART
        !   If the particle is yet to be added to a cluster, create a new cluster with the curent particle as the seed. 
        !   Increment the total number of clusters by one. Initialise the sum for the center-of-mass for this cluster 
        !   to the position of the seed particle.
            IF( PCLSTRID(J1) == 0 ) THEN
                NCLSTRS = NCLSTRS + 1
                CLSTRCOM(:,NCLSTRS) = R(:,J1)
                CLSTRSZ(NCLSTRS) = 1
                PCLSTRID(J1)  = NCLSTRS
                PREF(NCLSTRS) = J1
            ENDIF
        !   Convenient variable to use as the ID of the current cluster of interest.
            ICLU = PCLSTRID(J1)
        !   Check to which other particles the current particle is bonded to, making sure not to re-add particles
        !   which are already in the cluster. Update the sum for the center-of-mass for the cluster, using the seed particle
        !   as the anchor point (this avoids erroneous COMs for clusters which traverse boundaries).
            DO J2 = 1, NPART
                IF(J1==J2) CYCLE
                IF( CLSTRADJ(J1,J2) == 1 .AND. ICLU /= PCLSTRID(J2) ) THEN
                    CLSTRSZ(ICLU) = CLSTRSZ(ICLU) + 1
                    PCLSTRID(J2) = ICLU
                !   If building clusters for the computation of Steinhardt order parameter must also compute the
                !   distance between the seed particle of the cluster and the remaining particles.
                    IF(BOPSNT) THEN
                        RIJ = R(:,J1) - R(:,J2)
                        CLURIJ(:,J1,J2) = RIJ - BOX*ANINT( RIJ/BOX )
                    ENDIF
                !   Update the sum for the COM of the cluster
                    IF(CLSTRADJ(PREF(ICLU),J2) == 0) CLURIJ(:,PREF(ICLU),J2) = CLURIJ(:,J1,J2) + CLURIJ(:,PREF(ICLU),J1)
                    CLSTRCOM(:,ICLU) = CLSTRCOM(:,ICLU) + R(:,PREF(ICLU)) - CLURIJ(:,PREF(ICLU),J2)
                ENDIF
            ENDDO
        ENDDO
         
        DO J1 = 1, NCLSTRS
        !   Calculate the COM for the clusters and the relative positions of the particles to the COM.
            IF(CLSTRSZ(J1) > 1) THEN
                CLSTRCOM(:,J1) = CLSTRCOM(:,J1) / REAL(CLSTRSZ(J1),DP)
                IF(NPTT) THEN
            !   Compute the vector connecting each particle in the cluster to the COM. This will be used later to place the 
            !   particles in the correct position following a volume scaling move.
                    DO J2 = 1, NPART
                        IF(PCLSTRID(J2) == J1) THEN
                            PCLSTR(:,J2) = R(:,PREF(J1)) - CLURIJ(:,PREF(J1),J2) - CLSTRCOM(:,J1)
                        ENDIF
                    ENDDO
                ENDIF
            !   If building clusters for the computation of Steinhardt order parameter fold the COM back into the simulation cell.
                IF(BOPSNT) CLSTRCOM(:,J1) = CLSTRCOM(:,J1) - BOX*ANINT( CLSTRCOM(:,J1)/BOX )
            ENDIF
        ENDDO

    END SUBROUTINE

    SUBROUTINE VLM_CLUSTERMOVE(ISOTROPICT, BOXDIFF, BDIM, BOXD)

        USE COMMONS, ONLY: DP, BOX, R, NDIM, NPART
        USE COMMONS, ONLY: NCLSTRS, CLSTRCOM, PCLSTR, PCLSTRID
     
        IMPLICIT NONE

        INTEGER                              :: J1
        INTEGER, INTENT(IN), OPTIONAL        :: BDIM
        REAL(KIND=DP), INTENT(IN), OPTIONAL  :: BOXDIFF(NDIM), BOXD
        LOGICAL, INTENT(IN)                  :: ISOTROPICT

        CALL BUILD_ALL_CLUSTERS()

    !   PERFORM VOLUME MOVE AND SCALE THE POSITION OF THE PARTICLES ACCORDINGLY.
        DO J1 = 1, NCLSTRS
            IF(ISOTROPICT) THEN
                CLSTRCOM(:,J1) = CLSTRCOM(:,J1) * BOXDIFF
            ELSE
                CLSTRCOM(BDIM,J1) = CLSTRCOM(BDIM,J1) * BOXD
            ENDIF
        ENDDO
    !   Update the positions of the particles 
        DO J1 = 1, NPART
            R(:,J1) = PCLSTR(:,J1) + CLSTRCOM(:,PCLSTRID(J1))
            R(:,J1) = R(:,J1) - BOX*ANINT(R(:,J1)/BOX)
        ENDDO

    END SUBROUTINE

END MODULE