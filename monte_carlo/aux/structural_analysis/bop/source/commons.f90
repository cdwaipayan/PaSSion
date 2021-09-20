MODULE COMMONS

    USE ISO_C_BINDING, ONLY: C_DOUBLE
    IMPLICIT NONE
    SAVE

!   Parameters used to control the precision of floating point numbers used during the simulation.
!   Currently set to correspond to double and quadruple precision floating point numbers (DP and QP).
    INTEGER,PARAMETER          :: DP  = SELECTED_REAL_KIND(15,307)
    INTEGER,PARAMETER          :: QP  = SELECTED_REAL_KIND(33,4931)
!   Integer and floating point numbers definied in the correct double precision to be used throught
!   the program.
    REAL(KIND=DP), PARAMETER   :: PI = 4.0_dp*ATAN(1.0_dp), TWOPI = 2.0_dp*PI, FORPI = 4.0_dp*PI, HLFPI = PI/2.0_dp

    CHARACTER                  :: CHAR(80)

    ! Parameters used by readinput.f90 to read input files
    INTEGER :: LNGTH(20), LOC(20), NITEM

    INTEGER, PARAMETER :: MYUNIT = 11, MYUNIT2 = 21, VIEWUNIT = 31, INPUNIT = 15

!=============================================================================================
!   General Monte Carlo parameters
!=============================================================================================

    ! Parameters for cell list 
    LOGICAL                     :: CELLLISTT

    ! General system parameters
    INTEGER                     :: NDIM, NPART, NRBSITE, ISTEP, NDUMP, IDUMP, NNEIGH
    REAL(KIND=DP), ALLOCATABLE  :: R(:,:), SCALED_R(:,:), Q(:,:), RBSITES(:,:,:), REFSITE(:,:)
    REAL(KIND=DP)               :: VLM, BOX(3), RHO, RCUT, RCUTSQ
    LOGICAL                     :: DENSITYT, CUBICT, ORTHOT, RIGIDT, TETRAT, NNEIGHT

!------------------------------------------------------------------------------------------
!   Parameters related to the calculation of bond-orientational order parameters and 
!   the largest crystalline cluster for nucleation free-energy barrier calculations.
!------------------------------------------------------------------------------------------
    LOGICAL                     :: BOPT, QLT, QBARLT, GQLT, QLDOTQLT, LRGSTXCLSTRT
    INTEGER                     :: BOP_L, MIN_L, MAX_L, DEL_L, XSTL_NN, L_INDEX
    REAL(KIND=DP)               :: BOP_CUT, BOP_CUT_SQ, QLMIN, QLMAX
    REAL(KIND=DP), ALLOCATABLE  :: GQL(:), QL(:,:), QBARL(:,:), QLDOTQL(:,:,:)
    
END MODULE COMMONS 