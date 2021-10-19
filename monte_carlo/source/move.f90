
SUBROUTINE MOVE()
!     This subroutine moves the trajectory based on the simulation technique chosen.
    USE COMMONS, ONLY: DP, CDP, NPART, NPTT, FLFET, SCHSMIT
    USE FREE_ENERGY_MC

    IMPLICIT NONE

    INTEGER        :: I
    REAL(KIND=CDP) :: DRAND48

    IF(FLFET .OR. SCHSMIT) THEN
!   Calculate the free energy of crystal structures using the Einstein crystal or Schilling-Schmid approach.
        CALL FREE_ENERGY()
    ELSE
!   Perform a standard Monte Carlo simulation.
        DO I = 1, NPART
            
            IF (NPTT) THEN
                IF(DRAND48() < (1.0_dp/REAL(NPART,DP))) THEN
            !   Perform a volume scaing move on average once every N steps for simulations
            !   run under constant pressure.
                    CALL VMOVE
                    CYCLE
                ENDIF
            ENDIF

            CALL MV_MC

        ENDDO
    ENDIF

END SUBROUTINE MOVE
