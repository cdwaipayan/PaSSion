PROGRAM FULL_MONTE

!===================================================================================|
!-----------------------------------------------------------------------------------|
! Copyright (c) 2019 Andreas N. Neophytou                                           |
!                                                                                   |
!   THIS IS A MONTE CARLO PACKAGE DEVELOPED DURING THE COURSE OF MY PHD THAT IS     |
!   PRIMARILY AIMED AT INVESTIGATING THE SELF-ASSEMBLY OF SOFT MATTER SYSTEMS.      |
!                                                                                   |    
!   This package has the following functionalities:                                 |
!       1/ Standard Monte Carlo simulations in the NVT and NPT ensembles.           |
!       2/ Cluster Move and Virtual Move Monte Carlo                                |
!       3/ Crystal structure prediction using floppy box Monte Carlo                |
!       4/ Solid-liquid Gibbs free energy difference using interface pinning.       |
!       5/ Nucleation free energy barriers using nucleus-size pinning.              |
!       6/ Free energy calculations using Hamiltonian thermodynamic integration     |
!       7/ Free energy of solids using Frenkel-Ladd and Einstein Molecule methods   |
!       8/ Free energy difference between crystal structures using Lattice-switch   |
!          Monte Carlo                                                              |
!-----------------------------------------------------------------------------------|
!===================================================================================|
      
    USE COMMONS

    IMPLICIT NONE

    SEEDT    = .FALSE.
    CHCK_DIV = 50000
    DIV_TOL  = 1.e-10
    CALL CPU_TIME(START_TIME)

!---------------------------------------------------------------------
! Read input file
!---------------------------------------------------------------------
    CALL KEYWORDS()
!---------------------------------------------------------------------
! Initialise the system (i.e., set initial coordinates and conditions)
!---------------------------------------------------------------------
    CALL INITIALISE()
!---------------------------------------------------------------------
! Write a summary of the initial conditions to output file
!---------------------------------------------------------------------
    CALL SUMMARY_START()
!---------------------------------------------------------------------
! MAIN LOOP OF THE PROGRAM - Perform NSTEP MC cycles 
!---------------------------------------------------------------------
    ISTEP = 0
    CALL ACCUMULATORS()
    DO WHILE (ISTEP < NSTEP) 
        ISTEP = ISTEP + 1
        ICNT  = ICNT + 1
          
        CALL MOVE()

        IF (MOD(ISTEP, CHCK_DIV) == 0) CALL CHECK_DIV(ISTEP)

        IF (NUCSEEDT) THEN
            IF(MOD(ISTEP, TRJCTYLNGTH) == 0 ) THEN
                CALL UMBRELLA_BIAS(ISTEP)
                CALL ACCUMULATORS()
            ENDIF
        ELSE
            CALL ACCUMULATORS()
        ENDIF

        ! Writing to output
        IF (MOD(ISTEP,DMPFRQ) == 0 .OR. ISTEP==NSTEP) THEN
        !   Adjust step-sizes to obtain a preset acceptance criterion.
        !   Step-sizes are only adjusted during equilibration to ensure detailed balance
        !   is satisfied when calculating thermodynamic properties.
            IF(ISTEP <= NEQ .AND. ARATIOT) CALL ADJUST()
            CALL OUTPUT()
        ENDIF
    ENDDO

!---------------------------------------------------------------------
! Finalise the simulation
!---------------------------------------------------------------------
    CALL SUMMARY_END()

END PROGRAM FULL_MONTE
