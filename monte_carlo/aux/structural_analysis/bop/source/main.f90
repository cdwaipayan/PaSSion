program orderparam
    USE COMMONS
    USE CELL_LIST
    USE ORDERPARAM
    USE ROTATIONS_MODULE

    IMPLICIT NONE

    INTEGER         :: J1, J2, J3, LRGST_X_CLSTR, DUMMY, NT, NC, NCL, NH, FC, FCL, FH, FI, SKIP
    REAL(KIND=DP)   :: RM(3,3), AVE_C, AVE_H, AVE_I, AVE_CL
    REAL(KIND=DP), ALLOCATABLE :: QLQL3(:,:), QLQL4(:,:), QLQL6(:,:)

!   Paramters to be set before compiling and running the program
    NDIM    = 3                         ! Spatial dimensions of the system (i.e., 2 or 3)
    NPART   = 1000                       ! Number of particles in the system.
    RIGIDT  = .FALSE.                    ! True if simulating rigid bodies.
    NRBSITE = 2                         ! Number of rigid body sites per particle
    RHO     = 0.5_dp        ! Number density / Packing fraction of the system
    RCUT    = 1.2_dp          ! Cut-off for nearest neighbours
    RCUTSQ  = RCUT*RCUT

    DENSITYT = .TRUE.   ! True is using number density, false if using packing fraction
    CUBICT   = .TRUE.   ! True if simulating in a cubic box
    ORTHOT   = .FALSE.
    TETRAT   = .FALSE.  ! True is using tetrahedral centres

    NNEIGHT  = .FALSE.  ! True is using fixed number of nearest neighbours
    NNEIGH   = 8        ! Number of nearest neighbours to consider

    NDUMP = 1000        ! Number of snapshots to consider  
    SKIP  = 100         ! Number of snapshots to skip when printing assigned trajectory.    

!   Calculate the volume the simulation cell
    IF (DENSITYT) THEN
        VLM = REAL(NPART,DP)/RHO
    ELSE
        VLM = PI*REAL(NPART,DP) / (6.0_dp * RHO)
    ENDIF

!   Calculate the edge lengths of the simulation cell
    IF(CUBICT) THEN
        BOX = VLM**(1.0_dp / 3.0_dp)
    ELSEIF(ORTHOT) THEN
        BOX(1) = 10.620448472056D0
        BOX(2) = 12.263437568512D0
        BOX(3) = 86.715598653524D0
        VLM    = BOX(1)*BOX(2)*BOX(3) 
        ! STOP "NOT YET SET UP FOR ORTHORHOMBIC BOXES"    
        ! BOX = BOX * ( (VLM/(BOX(1)*BOX(2)*BOX(3)) )**(1.0_dp/3.0_dp) )
    ENDIF

    ALLOCATE( R(NDIM,NPART) )
    IF(RIGIDT) THEN
        ALLOCATE( Q(4,NPART), RBSITES(NDIM,NRBSITE,NPART) )
    ENDIF

    ALLOCATE(QLQL3(NPART,NPART), QLQL4(NPART,NPART), QLQL6(NPART,NPART))

!   *************************************************************************
!    Initialise parameters related to the use of a cell list.
!   *************************************************************************
    CELLLISTT = .TRUE.
    IF(CELLLISTT) THEN
        CALL INITIALIZE_LIST( NPART, RCUT/BOX )
        ALLOCATE(SCALED_R(NDIM,NPART))
        SCALED_R = 0.0_dp
    ENDIF
!   *************************************************************************
!    Initialise parameters related to the calculation of bond-orientational
!    order parameters and the largest crystalline cluster.
!   *************************************************************************
    MIN_L      = 3
    MAX_L      = 3
    DEL_L      = 1
    QLT        = .FALSE.
    QBARLT     = .FALSE.
    GQLT       = .TRUE.
    QLDOTQLT   = .TRUE.
    BOP_CUT    = RCUT
    BOP_CUT_SQ = BOP_CUT**2

!   Parameters related to the calculation of the largest crystalline cluster in the system
    LRGSTXCLSTRT = .FALSE.
    BOP_L   = 3             !6          ! l-value for ql_dot_ql BOP to be used
    QLMIN   = -1.1_dp       !0.2_dp     ! Minimum value for ql_dot_ql to be considered a crystalline bond
    QLMAX   = -0.85_dp      !0.65_dp     ! Maximum value for ql_dot_ql to be considered a crystalline bond
    XSTL_NN = 3             !6  

    IF(MIN_L == MAX_L) THEN
        J2 = 1
        IF(LRGSTXCLSTRT) L_INDEX = J2
    ELSE
        J2 = 0
        DO J1 = MIN_L, MAX_L, DEL_L
            J2 = J2 + 1
            IF(LRGSTXCLSTRT .AND. J1 == BOP_L) L_INDEX = J2
        ENDDO
    ENDIF

    IF(QLT)      ALLOCATE(QL(J2,NPART))
    IF(QBARLT)   ALLOCATE(QBARL(J2,NPART))
    IF(QLDOTQLT) ALLOCATE(QLDOTQL(J2,NPART,NPART))
    IF(GQLT)     ALLOCATE(GQL(J2))


    ! OPEN (UNIT = 1717, FILE ='tetrastack.xyz', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
    OPEN (UNIT = 1717, FILE ='diamond.xyz', STATUS = 'UNKNOWN', ACCESS = 'APPEND')

    IF(TETRAT) THEN
        OPEN (UNIT = 33, FILE = '../tet.dat', STATUS = 'OLD')
    ELSE
        OPEN (UNIT = 33, FILE = '../pos.dat', STATUS = 'OLD')
        IF(RIGIDT) OPEN (UNIT = 34, FILE = '../ortn.dat', STATUS = 'OLD')
    ENDIF

    AVE_C = 0.0_dp
    AVE_H = 0.0_dp
    AVE_I = 0.0_dp

    ISTEP = 0 

    DO J1 = 1, NDUMP
        
        R = 0.0_dp
        IF(RIGIDT) Q = 0.0_dp

        IF(TETRAT) READ(33,*) DUMMY
        DO J2 = 1, NPART
            READ(33,*) R(:,J2) ! Read in the translational coordinates
            R(:,J2) = R(:,J2) - BOX*ANINT(R(:,J2)/BOX) ! Take central image
            IF(CELLLISTT) SCALED_R(:,J2) = R(:,J2)/BOX

            IF(RIGIDT) THEN
                READ(34,*) Q(:,J2)
                RM = Q_TO_RM( Q(:,J2) )
                DO J3 = 1, NRBSITE
                    RBSITES(:,J3,J2) = MATMUL(RM,[0.0_dp, 0.0_dp, 1.0_dp])
                ENDDO
            ENDIF

        ENDDO

        IF(CELLLISTT) THEN
            IF(J1 == 1) THEN
                CALL MAKE_LIST( NPART, SCALED_R )
            ELSE
                CALL FINALIZE_LIST()
                CALL INITIALIZE_LIST( NPART, RCUT/BOX )
                CALL MAKE_LIST( NPART, SCALED_R )
            ENDIF
        ENDIF

        !   Calculate Bond-orientational order parameters (BOOP) based on the Steinhardt order parameters
        CALL GET_BOPS()
    !   Local BOOP q_l based just on the local environment of the particle    
        IF(QLT) THEN
            OPEN (UNIT = 23, FILE ='ql.out', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            DO J2 = 1, NPART
                WRITE (23,*) QL(:,J2)
            ENDDO
            CLOSE (UNIT = 23, STATUS = 'KEEP')
        ENDIF
    !   Local BOOP qbar_l which takes into account the local environment of the particle and its neearest neighbours.
        IF(QBARLT) THEN
            OPEN (UNIT = 24, FILE ='ql_bar.out', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            DO J2 = 1, NPART
                WRITE (24,*) QBARL(:,J2)
            ENDDO
            CLOSE (UNIT = 24, STATUS = 'KEEP')
        ENDIF
    !   Real component of the complex dot product based on q_l
        IF(QLDOTQLT) THEN
            OPEN (UNIT = 25, FILE ='ql_dot_ql.out', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            DO J2 = 1, NPART-1
                DO J3 = J2 + 1, NPART
                    IF(ANY(QLDOTQL(:,J2,J3) /= 100.0_dp)) WRITE (25,*) QLDOTQL(:,J2,J3)
                ENDDO
            ENDDO
            CLOSE (UNIT = 25, STATUS = 'KEEP')
        ENDIF
    !   Global BOOP Q_l
        IF(GQLT) THEN
            OPEN (UNIT = 26, FILE ='global_Ql.out', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            WRITE (26,*) GQL
            CLOSE (UNIT = 26, STATUS = 'KEEP')
        ENDIF
    !   Calculate the size of the largest crystalline cluster based on the complex dot product of q_l
        IF(LRGSTXCLSTRT) THEN
            OPEN (UNIT = 27, FILE ='lrgst_x_clstr.out', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            CALL LRGSTXCLSTR(.TRUE., LRGST_X_CLSTR)
            WRITE (27,*) LRGST_X_CLSTR
            CLOSE (UNIT = 27, STATUS = 'KEEP')
        ENDIF

!   ######################################################################################################
!       ASSIGN TETRASTACK ENVIRONMENTS
!   ######################################################################################################
        
    !     QLQL4 = QLDOTQL(1,:,:)
    !     QLQL6 = QLDOTQL(2,:,:)

    !     FC = 0
    !     FH = 0
    !     FI = 0
        
    !     IF(MOD(J1,SKIP)==0) THEN
    !         WRITE (1717,*) 3*NPART
    !         WRITE (1717,*) ""
    !     ENDIF

    !     DO J2 = 1, NPART
    !         NT = 0
    !         NC = 0
    !         NH = 0
    !         DO J3 = 1, NPART 

    !             IF(J2 == J3) CYCLE

    !             IF(QLQL6(J2,J3) > 0.2 .AND. QLQL6(J2,J3) < 0.65) THEN
    !                 NT = NT + 1
    !                 IF(QLQL4(J2,J3) > -0.2 .AND. QLQL4(J2,J3) < 0.35) THEN
    !                     NC = NC + 1
    !                 ELSEIF(QLQL4(J2,J3) > -0.65 .AND. QLQL4(J2,J3) < 0.35) THEN
    !                     NH = NH + 1
    !                 ENDIF
    !             ENDIF

    !         ENDDO ! Loop over particles j

    !         IF(NT > 5) THEN
    !             IF(NC > 5) THEN
    !                 FC = FC + 1
    !                 IF(MOD(J1,SKIP)==0) THEN
    !                     WRITE (1717,*) "C", R(:,J2)
    !                     WRITE (1717,*) "C", R(:,J2) + 0.5*RBSITES(:,1,J2)
    !                     WRITE (1717,*) "C", R(:,J2) - 0.5*RBSITES(:,2,J2)
    !                 ENDIF
    !             ELSEIF(NH > 1) THEN
    !                 FH = FH + 1
    !                 IF(MOD(J1,SKIP)==0) THEN
    !                     WRITE (1717,*) "O", R(:,J2)
    !                     WRITE (1717,*) "O", R(:,J2) + 0.5*RBSITES(:,1,J2)
    !                     WRITE (1717,*) "O", R(:,J2) - 0.5*RBSITES(:,2,J2)
    !                 ENDIF
    !             ELSE
    !                 FI = FI + 1
    !                 IF(MOD(J1,SKIP)==0) THEN
    !                     WRITE (1717,*) "N", R(:,J2)
    !                     WRITE (1717,*) "N", R(:,J2) + 0.5*RBSITES(:,1,J2)
    !                     WRITE (1717,*) "N", R(:,J2) - 0.5*RBSITES(:,2,J2)
    !                 ENDIF
    !             ENDIF
    !        ELSE 
    !             IF(MOD(J1,SKIP)==0) THEN
    !                 WRITE (1717,*) "H", R(:,J2)
    !                 WRITE (1717,*) "H", R(:,J2) + 0.5*RBSITES(:,1,J2)
    !                 WRITE (1717,*) "H", R(:,J2) - 0.5*RBSITES(:,2,J2)
    !             ENDIF
    !        ENDIF

    !     ENDDO
        
    !     AVE_C = AVE_C + REAL(FC)/REAL(NPART)
    !     AVE_H = AVE_H + REAL(FH)/REAL(NPART)
    !     AVE_I = AVE_I + REAL(FI)/REAL(NPART)
    !     PRINT *, REAL(FC)/REAL(NPART),REAL(FH)/REAL(NPART),REAL(FI)/REAL(NPART) 

!   ######################################################################################################
!       ASSIGN DIAMOND ENVIRONMENTS
!   ######################################################################################################

        QLQL3 = QLDOTQL(1,:,:)

        FC  = 0
        FH  = 0
        FI  = 0
        FCL = 0
        
        IF(MOD(J1,SKIP)==0) THEN
            WRITE (1717,*) NPART
            WRITE (1717,*) ""
        ENDIF

        DO J2 = 1, NPART
            NT = 0
            NC = 0
            NH = 0

            DO J3 = 1, NPART 

                IF(J2 == J3) CYCLE

                IF(QLQL3(J2,J3) > -1.1 .AND. QLQL3(J2,J3) < -0.85) THEN
                    NC = NC + 1
                ELSEIF(QLQL3(J2,J3) > -0.3 .AND. QLQL3(J2,J3) < 0.1) THEN
                    NH = NH + 1
                ENDIF

            ENDDO ! Loop over particles j

            IF(NC+NH >= 3 .AND. NC+NH <= 4) THEN
                IF(NC == 4 .AND. NH == 0) THEN
                    FC = FC + 1
                    IF(MOD(J1,SKIP)==0) WRITE (1717,*) "C", R(:,J2)
                ELSEIF(NC == 3 .AND. NH == 1) THEN
                    FH = FH + 1
                    IF(MOD(J1,SKIP)==0) WRITE (1717,*) "O", R(:,J2)
                ELSEIF(NC == 0 .AND. NH == 4) THEN
                    FCL = FCL + 1
                    IF(MOD(J1,SKIP)==0) WRITE (1717,*) "Cl", R(:,J2)
                ELSE
                    FI = FI + 1
                    IF(MOD(J1,SKIP)==0) WRITE (1717,*) "N", R(:,J2)
                ENDIF
           ELSE 
                IF(MOD(J1,SKIP)==0) WRITE (1717,*) "H", R(:,J2)
           ENDIF

        ENDDO
        
        AVE_C  = AVE_C  + REAL(FC)/REAL(NPART)
        AVE_H  = AVE_H  + REAL(FH)/REAL(NPART)
        AVE_CL = AVE_CL + REAL(FCL)/REAL(NPART)
        AVE_I  = AVE_I  + REAL(FI)/REAL(NPART)
        PRINT *, REAL(FC)/REAL(NPART),REAL(FH)/REAL(NPART),REAL(FCL)/REAL(NPART),REAL(FI)/REAL(NPART) 

    ENDDO

    PRINT *, "AVERAGE(C,H,CL,I): ", AVE_C/REAL(NDUMP), AVE_H/REAL(NDUMP), AVE_CL/REAL(NDUMP), AVE_I/REAL(NDUMP)

    ! OPEN (UNIT = 2006, FILE ='test.xyz', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
    ! WRITE (2006,*) NPART
    ! WRITE (2006,*) 
    ! DO J2 = 1, NPART
    !     WRITE (2006,*) "C", R(:,J2) ! Read in the translational coordinates
    ! ENDDO

end program
