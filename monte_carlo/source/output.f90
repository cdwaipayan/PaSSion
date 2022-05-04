SUBROUTINE OUTPUT()

    USE COMMONS
    USE ORDERPARAM
    IMPLICIT NONE

    INTEGER            :: J1, J2, NSTART, LRGST_X_CLSTR
    character(len=100) :: fmt1, snbop

    IF (.NOT. SEGREGATET) THEN
        OPEN (UNIT = 3, FILE ='data.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
        !   Current Monte Carlo cycle/step 
            WRITE (3,*) "STEP: ", ISTEP
        !   Instantaneous energy
            WRITE (3,*) "ENERGY: ", PE
        !   Instantaneous pressure (only valid for continuous models)       
            WRITE (3,*) "PRESSURE: ", PRES
            IF(NPTT) THEN
            !   Instantaneous density      
                WRITE (3,*) "DENSITY: ", RHO
            !   Instantaneous volume  
                WRITE (3,*) "VOLUME: ", VLM
            !   Collective density field
                IF(NPZTT .OR. PINT) WRITE (3,*) "COLL_DENS: ", TRQ
            ENDIF
        CLOSE (UNIT = 3, STATUS = 'KEEP')

    ELSE !PRINT SEPERATE FILES
        OPEN (UNIT = 54, FILE ='energy.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
        !   Current Monte Carlo cycle/step & instantaneous energy
            WRITE (54,*) ISTEP, PE
        CLOSE (UNIT = 54, STATUS = 'KEEP')

        OPEN (UNIT = 55, FILE ='pressure.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
        !   Instantaneous pressure (only valid for continuous models)       
            WRITE (55,*) ISTEP, PRES
        CLOSE (UNIT = 55, STATUS = 'KEEP')

        IF(NPTT) THEN
            OPEN (UNIT = 56, FILE ='density.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            !   Instantaneous density      
                WRITE (56,*) ISTEP, RHO
            CLOSE (UNIT = 56, STATUS = 'KEEP')

            OPEN (UNIT = 57, FILE ='volume.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            !   Instantaneous volume  
                WRITE (57,*) ISTEP, VLM
            CLOSE (UNIT = 57, STATUS = 'KEEP')

            IF(NPZTT .OR. PINT) THEN
                OPEN (UNIT = 58, FILE ='coll_dens.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
                !   Collective density field
                    WRITE (58,*) TRQ
                CLOSE (UNIT = 58, STATUS = 'KEEP')
            ENDIF
        ENDIF
    ENDIF

    IF(BOPT) THEN
        WRITE(fmt1,'(I6)') (MAX_L-MIN_L)/DEL_L+1
        snbop = '('//trim(adjustl(fmt1))//'F12.7)'
    !   Calculate Bond-orientational order parameters (BOOP) based on the Steinhardt order parameters
        CALL GET_BOPS()
    !   Local BOOP q_l based just on the local environment of the particle    
        IF(QLT) THEN
            OPEN (UNIT = 23, FILE ='ql.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            DO J1 = 1, NPBOP
                IF(EQUISEEDT) THEN
                    IF(J1 > SEEDSIZE) CYCLE
                ENDIF
                WRITE (23,snbop) QL(:,J1)
            ENDDO
            CLOSE (UNIT = 23, STATUS = 'KEEP')
        ENDIF
    !   Local BOOP qbar_l which takes into account the local environment of the particle and its neearest neighbours.
        IF(QBARLT) THEN
            OPEN (UNIT = 24, FILE ='ql_bar.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            DO J1 = 1, NPBOP
                IF(EQUISEEDT) THEN
                    IF(J1 > SEEDSIZE) CYCLE
                ENDIF
                WRITE (24,snbop) QBARL(:,J1)
            ENDDO
            CLOSE (UNIT = 24, STATUS = 'KEEP')
        ENDIF
    !   Real component of the complex dot product based on q_l
        IF(QLDOTQLT) THEN
            OPEN (UNIT = 25, FILE ='ql_dot_ql.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            DO J1 = 1, NPBOP-1
                IF(EQUISEEDT) THEN
                    IF(J1 > SEEDSIZE) CYCLE
                ENDIF
                DO J2 = J1 + 1, NPBOP
                    IF(ANY(QLDOTQL(:,J1,J2) /= 0.0_dp)) WRITE (25,snbop) QLDOTQL(:,J1,J2)
                ENDDO
            ENDDO
            CLOSE (UNIT = 25, STATUS = 'KEEP')
        ENDIF
    !   Global BOOP Q_l
        IF(GQLT) THEN
            OPEN (UNIT = 26, FILE ='global_Ql.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            WRITE (26,snbop) GQL
            CLOSE (UNIT = 26, STATUS = 'KEEP')
        ENDIF
    !   Calculate the size of the largest crystalline cluster based on the complex dot product of q_l
        IF(LRGSTXCLSTRT) THEN
            OPEN (UNIT = 27, FILE ='lrgst_x_clstr.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            CALL LRGSTXCLSTR(.TRUE., LRGST_X_CLSTR)
            
            IF(NUCSEEDT) THEN
                WRITE (27,'(F12.7,I4)') SUMLRGSTXCLSTR/REAL(ICNT,DP), LRGST_X_CLSTR
            ELSE
                WRITE (27,'(I4)') LRGST_X_CLSTR
            ENDIF
            CLOSE (UNIT = 27, STATUS = 'KEEP')

            IF(GETSEEDT) THEN
                IF(LRGST_X_CLSTR == TRGTSEED) THEN
                    OPEN (UNIT = 28, FILE ='seedpos.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
                    DO J1 = 1, NPART
                        WRITE (28,*) R(:,J1)
                    ENDDO
                    IF (RIGIDT) THEN
                        OPEN (UNIT = 29,  FILE = 'seedortn.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
                        DO J1 = 1, NPART
                            WRITE (29,*)  Q(:,J1)
                        ENDDO
                        CLOSE (UNIT = 29, STATUS = 'KEEP')
                    ENDIF
                    STOP
                ENDIF
            ENDIF

        ENDIF

    ENDIF

    NSTART = NEQ

    IF (ISTEP > NSTART .AND. (.NOT. FLFET) .AND. (.NOT. SCHSMIT) .AND. MOD(ISTEP,POSDMP)==0) THEN

        OPEN (UNIT = 7,  FILE = 'pos.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
        DO J1 = 1, NPART
            WRITE (7,'(3F12.7)') R(:,J1)
        ENDDO
        CLOSE (UNIT = 7, STATUS = 'KEEP')

        IF (RIGIDT) THEN
            OPEN (UNIT = 8,  FILE = 'ortn.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            DO J1 = 1, NPART
                WRITE (8,'(4F12.7)')  Q(:,J1)
            ENDDO
            CLOSE (UNIT = 8, STATUS = 'KEEP')
        ENDIF

        IF(NPTT) THEN
            OPEN (UNIT = 19, FILE = 'box.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            WRITE (19,'(3F12.7)') BOX
            CLOSE (UNIT = 19, STATUS = 'KEEP')
        ENDIF

    ENDIF

    IF (VERBOSET) THEN
        WRITE(*,'(I8,A1,I8,A5,F12.5,A8,F6.2,A7,F9.3,A7,F5.3,A12,F8.3)') ISTEP, "/", NSTEP, " PE: ", PE, "  BETA: ", BETAKB, &
            & "  VOL: ", VLM, "  RHO: ", RHO, "  AVE PRES: ", AVPRES
        CALL VIEWCONFIG('traj.xyz')
    ENDIF !Verbose

    IF (ISTEP == NSTEP) THEN

    !   Output final translational coordinates of the particles
        OPEN(UNIT = 37,  FILE = 'finalpos.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
        DO J1 = 1, NPART
            WRITE (37,*) R(:,J1)
        ENDDO
        CLOSE (UNIT = 37, STATUS = 'KEEP')

    !   If the particles are rigid bodies, output final rotational coordinates of the particles
        IF (RIGIDT) THEN
            OPEN (UNIT = 38,  FILE = 'finalortn.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            DO J1 = 1, NPART
                WRITE (38,*) Q(:,J1)
            ENDDO
        ENDIF
        CLOSE (UNIT = 38, STATUS = 'KEEP')

        IF(NPTT) THEN
            OPEN (UNIT = 39, FILE = 'finalbox.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            OPEN (UNIT = 40, FILE = 'finaldns.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            WRITE (39,*) BOX
            WRITE (40,*) RHO
            CLOSE (UNIT = 39, STATUS = 'KEEP')
            CLOSE (UNIT = 40, STATUS = 'KEEP')
        ENDIF

    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    !   Output values of sums to file in order to continue calculating averages in
    !   further simulations.
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
        OPEN(UNIT = 41,  FILE = 'finalave.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')

        WRITE (41,*) IBLOCK
        WRITE (41,*) SUMPE, STDPESUM
        WRITE (41,*) SUMPE2, STDPE2SUM
        WRITE (41,*) SUMPRES, STDPRESSUM
        IF(NPTT) THEN
            WRITE (41,*) SUMVLM, STDVLMSUM
            WRITE (41,*) SUMVLM2, STDVLM2SUM
            WRITE (41,*) SUMH, STDHSUM
            WRITE (41,*) SUMH2, STDH2SUM
            WRITE (41,*) SUMHV, STDRHOSUM
            WRITE (41,*) SUMRHO, STDHSUM
            WRITE (41,*) SUMBOX, STDBOXSUM
        ENDIF

    !----------------------------------------------------------------------------------
    !   Chemical potential and order parameter for interface-pinning simulations
    !----------------------------------------------------------------------------------
        IF(PINT) THEN
            WRITE (41,*) SUMTRQ, STDTRQSUM 
            WRITE (41,*) SUMDELMU, STDDELMUSUM
        ENDIF

    !----------------------------------------------------------------------------------
    !   Count of largest crystalline cluster for umbrella sampling and nucleus-size
    !   pinning simulations.
    !----------------------------------------------------------------------------------
        IF (NUCSEEDT) THEN
            WRITE (41,*) SUMLRGSTXCLSTR, STDLRGSTXCLSTRSUM
            WRITE (41,*) SUMNUCCOUNT, STDNUCCOUNTSUM
            WRITE (41,*) NUCCOUNT, TOTNUCCOUNT
        ENDIF

    !----------------------------------------------------------------------------------
    !   Free energies from Einstein-crystal integration simulations
    !----------------------------------------------------------------------------------
        IF(FLFET) THEN
            IF(FLID==1) THEN
                WRITE (41,*) SUM_EXP_U_EIN, STD_EXP_U_EIN_SUM
            ELSE
                WRITE (41,*) SUM_U_EIN_TR, STD_U_EIN_TR_SUM
                WRITE (41,*) SUM_U_EIN_OR, STD_U_EIN_OR_SUM
            ENDIF
        ENDIF
    
    !----------------------------------------------------------------------------------
    !   Free energies from Schilling-Schmid integration simulations
    !----------------------------------------------------------------------------------
        IF(SCHSMIT) THEN
            IF(SSID==1) THEN
                WRITE (41,*) SUM_EXP_U_SS, STD_EXP_U_SS_SUM
            ELSE
                WRITE (41,*) SUM_U_SS, STD_U_SS_SUM
            ENDIF
        ENDIF

    !   Output xyz file to visualise the final configuration
        CALL VIEWCONFIG('finalconfig.xyz')

    ENDIF

END SUBROUTINE OUTPUT

SUBROUTINE SUMMARY_START()
    
    USE COMMONS

    IMPLICIT NONE

    character(len=100) :: fmt1, fmt2, fmt3
    character(len=100) :: string

    INTEGER             :: J1

    OPEN (UNIT = 1771, FILE = 'output.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
    WRITE (1771,*) "***********************************************************************************"
    WRITE(fmt1,'(I6)') NPART
    string = 'Simulating '//trim(adjustl(fmt1))//' particles'

    IF(LATTICET) THEN
        IF (LID == 1)  WRITE (1771,*) "Starting from a Face-Centered Cubic (fcc) crystal ..."
        IF (LID == 2)  WRITE (1771,*) "Starting from a Hexagonal Close-Packed (hcp) crystal ..."
        IF (LID == 3)  WRITE (1771,*) "Starting from a Double Hexagonal Close-Packed (dhcp) crystal ..."
        IF (LID == 4)  WRITE (1771,*) "Starting from a Body-Centered Cubic (bcc) crystal ..."
        IF (LID == 5)  WRITE (1771,*) "Starting from a Simple Cubic crystal ..."
        IF (LID == 6)  WRITE (1771,*) "Starting from a Cubic Diamond crystal ..."
        IF (LID == 7)  WRITE (1771,*) "Starting from a Hexagonal Diamond crystal ..."
        IF (LID == 8)  WRITE (1771,*) "Starting from a Cubic Tetrastack crystal ..."
        IF (LID == 9)  WRITE (1771,*) "Starting from a Double Cubic Tetrastack (aka fcc) crystal ..."
        IF (LID == 10) WRITE (1771,*) "Starting from a Hexagonal Tetrastack crystal ..."
        IF (LID == 11) WRITE (1771,*) "Starting from a Cubic Tetrahedral Diamond crystal ..."
        IF (LID == 12) WRITE (1771,*) "Starting from a Double Tetrahedral Diamond crystal ..."
        IF (LID == 13) WRITE (1771,*) "Starting from a Hexagonal Tetrahedral Diamond crystal ..."
        IF (LID == 14) WRITE (1771,*) "Starting from a Gyroid crystal ..."
        IF (LID == 15) WRITE (1771,*) "Starting from a Double Gyroid crystal ..."
    ENDIF

    IF(FLFET) THEN
        IF(FLID == 1) THEN
            WRITE (1771,*) "Performing Frenkel-Ladd free energy calculation (DA1)"
        ELSEIF(FLID == 2) THEN
            WRITE (1771,*) "Performing Frenkel-Ladd free energy calculation (DA2)"
        ENDIF
    ELSEIF(PINT) THEN
        WRITE(fmt1,'(F20.4)') PINKAP
        WRITE(fmt2,'(F20.4)') PINA
        string = 'Starting an interface pinning simulation with k=' &
        //trim(adjustl(fmt1))//' and a='//trim(adjustl(fmt2))
        WRITE (1771,*) string
    ELSEIF(NUCSEEDT) THEN
        WRITE(fmt1,'(I3)') TRGTSEED
        WRITE(fmt2,'(F20.4)') NPINK
        WRITE(fmt3,'(I3)') TRJCTYLNGTH
        string = 'Umbrella sampling simulation with window centre at ' &
        //trim(adjustl(fmt1))//', kappa='//trim(adjustl(fmt2))//' and trajectory length of ' &
        //trim(adjustl(fmt3))
        WRITE (1771,*) string
    ELSEIF(EQUISEEDT) THEN
        IF (LID == 1)  WRITE (1771,*) "Seeding fluid with a Face-Centered Cubic (fcc) seed and applying a harmonic trap ..."
        IF (LID == 2)  WRITE (1771,*) "Seeding fluid with a Hexagonal Close-Packed (hcp) seed and applying a harmonic trap ..."
        IF (LID == 3)  WRITE (1771,*) "Seeding fluid with a Double HCP seed and applying a harmonic trap ..."
        IF (LID == 4)  WRITE (1771,*) "Seeding fluid with a Body-Centered Cubic (bcc) seed and applying a harmonic trap ..."
        IF (LID == 5)  WRITE (1771,*) "Seeding fluid with a Simple Cubic seed and applying a harmonic trap ..."
        IF (LID == 6)  WRITE (1771,*) "Seeding fluid with a Cubic Diamond seed and applying a harmonic trap ..."
        IF (LID == 7)  WRITE (1771,*) "Seeding fluid with a Hexagonal Diamond seed and applying a harmonic trap ..."
        IF (LID == 8)  WRITE (1771,*) "Seeding fluid with a Cubic Tetrastack seed and applying a harmonic trap ..."
        IF (LID == 9)  WRITE (1771,*) "Seeding fluid with a Double Cubic Tetrastack (aka fcc) seed and applying a harmonic trap ..."
        IF (LID == 10) WRITE (1771,*) "Seeding fluid with a Hexagonal Tetrastack seed and applying a harmonic trap ..."
        IF (LID == 11) WRITE (1771,*) "Seeding fluid with a Cubic Tetrahedral Diamond seed and applying a harmonic trap ..."
        IF (LID == 12) WRITE (1771,*) "Seeding fluid with a Double Tetrahedral Diamond seed and applying a harmonic trap ..."
        IF (LID == 13) WRITE (1771,*) "Seeding fluid with a Hexagonal Tetrahedral Diamond seed and applying a harmonic trap ..."
        IF (LID == 14) WRITE (1771,*) "Seeding fluid with a Gyroid seed and applying a harmonic trap ..."
        IF (LID == 15) WRITE (1771,*) "Seeding fluid with a Double Gyroid seed and applying a harmonic trap ..."
    ELSEIF(NPTT) THEN
        WRITE (1771,*) "Starting a Monte-Carlo simulation in the NPT ensemble ..."
        WRITE(fmt1,'(F20.12)') PRSFIX
        string = 'Pressure set to '//trim(adjustl(fmt1))
        write(1771,*) string
    ELSEIF(NPZTT) THEN
        WRITE (1771,*) "Starting a Monte-Carlo simulation in the NPzT ensemble ..."
        WRITE(fmt1,'(F20.12)') PRSFIX
        string = 'Pressure set to '//trim(adjustl(fmt1))
        write(1771,*) string
    ELSE
        WRITE (1771,*) "Starting a Monte-Carlo simulation in the NVT ensemble ..."
    ENDIF

    WRITE(fmt1,'(F20.12)') 1.0_dp/BETAKB
    WRITE(fmt2,'(F20.12)') BETAKB
    string = 'Temperature set to '//trim(adjustl(fmt1))//' (beta='//trim(adjustl(fmt2))//')'
    write(1771,*) string

    WRITE(fmt1,'(I9)') NSTEP
    WRITE(fmt2,'(I9)') BLKLNGTH
    WRITE(fmt3,'(I9)') NEQ
    string = 'Performing '//trim(adjustl(fmt1))//' cycles with blocks of width '//trim(adjustl(fmt2))//' and ' &
    & //trim(adjustl(fmt3))//' equilibration steps'
    write(1771,*) string

    WRITE(fmt1,'(I9)') CHCK_DIV
    WRITE(fmt2,'(F20.12)') DIV_TOL
    string = 'Energy divergence checked every '//trim(adjustl(fmt1))//' MC cycles with a tolerance of '//trim(adjustl(fmt2))
    write(1771,*) string

    IF(CELLLISTT) write(1771,*) "Using a cell-list for efficiency"
    IF(CLUSTERMOVET) write(1771,*) "Using rigid-body cluster moves"

    WRITE(fmt1,'(F20.12)') RCUT
    string = 'Using a radial cut-off of: '//trim(adjustl(fmt1))
    write(1771,*) string

    IF (HST) THEN
        write(1771,*) 'Simulating hard-sphere particles'

    ELSEIF (GLJT) THEN
        write(1771,*) 'Simulating GLJ particles'
        WRITE(fmt1,'(F20.12)') GLJN
        string = '  + Using an exponent of: '//trim(adjustl(fmt1))
        write(1771,*) string

    ELSEIF (KIHARAT) THEN
        write(1771,*) 'Simulating kihara particles'
        WRITE(fmt1,'(F20.12)') GLJN
        WRITE(fmt2,'(F20.12)') RLNGTH
        string = '  + Using an exponent of: '//trim(adjustl(fmt1))//' with L* = '//trim(adjustl(fmt2))
        write(1771,*) string

    ELSEIF (KFT) THEN
        write(1771,*) 'Simulating Kern-Frenkel particles'
        write(1771,*) string
        DO J1 = 1, NSITES
            IF(J1 == 1) THEN
                WRITE(fmt2,'(F20.12)') KFDELA
                string = '  + Patch A has a width of: '//trim(adjustl(fmt2))
                write(1771,*) string
                WRITE(fmt2,'(F20.12)') KFAA
                string = '  + Patch A has a well-depth of: '//trim(adjustl(fmt2))
                write(1771,*) string
                WRITE(fmt2,'(F20.12)') KFLAMA
                string = '  + Patch A has a range of: '//trim(adjustl(fmt2))
                write(1771,*) string
            ELSEIF(J1 == 2) THEN
                WRITE(fmt2,'(F20.12)') KFDELB
                string = '  + Patch B has a width of: '//trim(adjustl(fmt2))
                write(1771,*) string
                WRITE(fmt2,'(F20.12)') KFBB
                string = '  + Patch B has a well-depth of: '//trim(adjustl(fmt2))
                write(1771,*) string
                WRITE(fmt2,'(F20.12)') KFLAMB
                string = '  + Patch B has a range of: '//trim(adjustl(fmt2))
                write(1771,*) string
            ELSEIF(J1 == 3) THEN
                WRITE(fmt2,'(F20.12)') KFDELC
                string = '  + Patch C has a width of: '//trim(adjustl(fmt2))
                write(1771,*) string
                WRITE(fmt2,'(F20.12)') KFCC
                string = '  + Patch C has a well-depth of: '//trim(adjustl(fmt2))
                write(1771,*) string
                WRITE(fmt2,'(F20.12)') KFLAMC
                string = '  + Patch C has a range of: '//trim(adjustl(fmt2))
                write(1771,*) string
            ELSEIF(J1 == 4) THEN
                WRITE(fmt2,'(F20.12)') KFDELD
                string = '  + Patch D has a width of: '//trim(adjustl(fmt2))
                write(1771,*) string
                WRITE(fmt2,'(F20.12)') KFDD
                string = '  + Patch D has a well-depth of: '//trim(adjustl(fmt2))
                write(1771,*) string
                WRITE(fmt2,'(F20.12)') KFLAMD
                string = '  + Patch D has a range of: '//trim(adjustl(fmt2))
                write(1771,*) string
            ENDIF
        ENDDO

    ELSEIF (GYROIDALT) THEN
        write(1771,*) 'Simulating "gyroidal" Kern-Frenkel particles'
        write(1771,*) string
        !PAIR A
        WRITE(fmt2,'(F20.12)') KFDELA
        string = '  + Patch pair A have widths of: '//trim(adjustl(fmt2))
        write(1771,*) string
        WRITE(fmt2,'(F20.12)') GYBETA1
        string = '  + Patch pair A have a beta of: '//trim(adjustl(fmt2))
        write(1771,*) string
        WRITE(fmt2,'(F20.12)') KFAA
        string = '  + Patch pair A have well-depths of: '//trim(adjustl(fmt2))
        write(1771,*) string
        WRITE(fmt2,'(F20.12)') KFLAMA
        string = '  + Patch pair A have ranges of: '//trim(adjustl(fmt2))
        write(1771,*) string
        !PAIR B
        WRITE(fmt2,'(F20.12)') KFDELC
        string = '  + Patch pair B have widths of: '//trim(adjustl(fmt2))
        write(1771,*) string
        WRITE(fmt2,'(F20.12)') GYBETA2
        string = '  + Patch pair B have a beta of: '//trim(adjustl(fmt2))
        write(1771,*) string
        WRITE(fmt2,'(F20.12)') KFCC
        string = '  + Patch pair B have well-depths of: '//trim(adjustl(fmt2))
        write(1771,*) string
        WRITE(fmt2,'(F20.12)') KFLAMC
        string = '  + Patch pair B have ranges of: '//trim(adjustl(fmt2))
        write(1771,*) string

    ELSEIF (PGLJT) THEN
        write(1771,*) 'Simulating patchy GLJ particles'
    
    ELSEIF (CPPT .OR. ETPT) THEN

        IF (CPPT) THEN
            write(1771,*) 'Simulating continuous Kern-Frenkel particles'
            WRITE(fmt1,'(F20.12)') YUKKAP
            string = '  + Using a kappa of: '//trim(adjustl(fmt1))
            write(1771,*) string
        ELSEIF(ETPT) THEN
            write(1771,*) 'Simulating ETP particles'
            WRITE(fmt1,'(F20.12)') RLNGTH
            string = '  + Using an L* of: '//trim(adjustl(fmt1))
            write(1771,*) string
            WRITE(fmt1,'(F20.12)') SQRT(ETPCUTSQ)
            string = '  + Using a shortest distance cut-off of: '//trim(adjustl(fmt1))
            write(1771,*) string
        ENDIF

        WRITE(fmt1,'(F20.12)') CPPS
        string = '  + Using an S value of: '//trim(adjustl(fmt1))
        write(1771,*) string

        WRITE(fmt1,'(F20.12)') CPPLAM
        string = '  + Using a lambda value of: '//trim(adjustl(fmt1))
        write(1771,*) string

        DO J1 = 1, NSITES
            IF(J1 == 1) THEN

                WRITE(fmt2,'(F20.12)') CPPAA
                string = '  + Patch A has a well-depth of: '//trim(adjustl(fmt2))
                write(1771,*) string
                
                WRITE(fmt2,'(F20.12)') CPPDELA
                string = '  + Patch A has a width of: '//trim(adjustl(fmt2))
                write(1771,*) string

            ELSEIF(J1 == 2) THEN

                WRITE(fmt2,'(F20.12)') CPPBB
                string = '  + Patch B has a well-depth of: '//trim(adjustl(fmt2))
                write(1771,*) string

                WRITE(fmt2,'(F20.12)') CPPDELB
                string = '  + Patch B has a width of: '//trim(adjustl(fmt2))
                write(1771,*) string

            ELSEIF(J1 == 3) THEN

                WRITE(fmt2,'(F20.12)') CPPCC
                string = '  + Patch C has a well-depth of: '//trim(adjustl(fmt2))
                write(1771,*) string

                WRITE(fmt2,'(F20.12)') CPPDELC
                string = '  + Patch C has a width of: '//trim(adjustl(fmt2))
                write(1771,*) string

            ELSEIF(J1 == 4) THEN

                WRITE(fmt2,'(F20.12)') CPPDD
                string = '  + Patch D has a well-depth of: '//trim(adjustl(fmt2))
                write(1771,*) string

                WRITE(fmt2,'(F20.12)') CPPDELD
                string = '  + Patch D has a width of: '//trim(adjustl(fmt2))
                write(1771,*) string

            ENDIF
        ENDDO
    ENDIF

    WRITE(fmt1,'(F20.12)') PE
    string = 'Initial potential energy: '//trim(adjustl(fmt1))
    write(1771,*) string

    WRITE(fmt1,'(F20.12)') VLM
    string = 'Initial volume: '//trim(adjustl(fmt1))
    write(1771,*) string

    WRITE(fmt1,'(F20.12)') RHO
    IF(PACKINGT) THEN
        string = 'Initial packing fraction: '//trim(adjustl(fmt1))
    ELSE
        string = 'Initial density: '//trim(adjustl(fmt1))
    ENDIF
    write(1771,*) string

    WRITE(fmt1,'(F20.12)') BOX(1)
    WRITE(fmt2,'(F20.12)') BOX(2)
    WRITE(fmt3,'(F20.12)') BOX(3)
    string = 'Initial box lengths: '//trim(adjustl(fmt1))//' '//trim(adjustl(fmt2))//' '//trim(adjustl(fmt3))
    write(1771,*) string

    IF(LRGSTXCLSTRT) THEN
        WRITE(fmt1,'(I9)') BOP_L
        string = 'Following the size of the largest crystalline cluster &
        &with q'//trim(adjustl(fmt1))//'(i)*q'//trim(adjustl(fmt1))//'(j)'
        write(1771,*) string
    ENDIF

    FLUSH(1771)

END SUBROUTINE

SUBROUTINE SUMMARY_END()
    
    USE COMMONS

    IMPLICIT NONE

    character(len=100) :: fmt1, fmt2, fmt3
    character(len=100) :: string
    REAL(KIND=DP)      :: FINISH_TIME

    WRITE(1771,*) ""

    IF(PINT) THEN
        WRITE (1771,*) "Finishing the interface pinning simulation ..."
        WRITE(fmt1,'(F20.12)') AVDELMU
        string = 'Final chemical potential difference (solid-liquid): '//trim(adjustl(fmt1))
        write(1771,*) string
    ELSEIF(NUCSEEDT) THEN
        WRITE (1771,*) "Finishing the nucleus-size pinning simulation ..."
    ELSE
        WRITE (1771,*) "Finishing the Monte-Carlo simulation ..."
    ENDIF

    CELLLISTT = .FALSE.
    CALL POTENTIAL(PE)
    WRITE(fmt1,'(F20.12)') PE
    string = 'Final total potential energy:                        '//trim(adjustl(fmt1))
    write(1771,*) string

    WRITE(fmt1,'(F20.12)') AVPE
    string = 'Final average potential energy: '//trim(adjustl(fmt1))
    write(1771,*) string

    WRITE(fmt1,'(F20.12)') AVPRES
    string = 'Final average pressure: '//trim(adjustl(fmt1))
    write(1771,*) string

    WRITE(fmt1,'(F20.12)') VLM
    string = 'Final volume: '//trim(adjustl(fmt1))
    write(1771,*) string

    WRITE(fmt1,'(F20.12)') RHO
    IF(PACKINGT) THEN
        string = 'Final packing fraction: '//trim(adjustl(fmt1))
    ELSE
        string = 'Final density: '//trim(adjustl(fmt1))
    ENDIF
    write(1771,*) string

    WRITE(fmt1,'(F20.12)') BOX(1)
    WRITE(fmt2,'(F20.12)') BOX(2)
    WRITE(fmt3,'(F20.12)') BOX(3)
    string = 'Final box lengths: '//trim(adjustl(fmt1))//' '//trim(adjustl(fmt2))//' '//trim(adjustl(fmt3))
    write(1771,*) string
    
    CALL CPU_TIME(FINISH_TIME)
    WRITE(fmt1,'(F20.12)') ((FINISH_TIME - START_TIME)/60.0) / 24.0
    string = 'Total run time was '//trim(adjustl(fmt1))//' hours'
    write(1771,*) string

    WRITE (1771,*) "**********************************************************************************"

END SUBROUTINE
