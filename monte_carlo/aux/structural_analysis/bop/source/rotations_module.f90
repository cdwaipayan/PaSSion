MODULE ROTATIONS_MODULE

  !------------------------------------------------------------------------------------------------!
  ! Adapted from https://github.com/Allen-Tildesley/examples/blob/master/maths_module.f90
  !------------------------------------------------------------------------------------------------!

    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY : ERROR_UNIT, IOSTAT_END, IOSTAT_EOR
    USE COMMONS, ONLY: DP, PI, TWOPI

    IMPLICIT NONE
    PRIVATE

!   Private data
    REAL(KIND=DP), PARAMETER :: TOL = 1.e+06

!   Public routines which operate on quaternions
    PUBLIC :: ROTATE_QUATERNION, QUATMUL, Q_TO_RM
!   Public routines to convert between quaternions and unit vectors
    PUBLIC :: Q_TO_UV, UV_TO_Q

CONTAINS

    FUNCTION ROTATE_QUATERNION ( ANGLE, AXIS, OLD ) RESULT ( E )
        IMPLICIT NONE
        REAL(KIND=DP), DIMENSION(0:3)             :: E     ! Returns a quaternion rotated by
        REAL(KIND=DP),                 INTENT(in) :: ANGLE ! specified rotation angle (in radians) about
        REAL(KIND=DP), DIMENSION(3),   INTENT(in) :: AXIS  ! specified rotation axis relative to
        REAL(KIND=DP), DIMENSION(0:3), INTENT(in) :: OLD   ! old quaternion

!   Note that the axis vector should be normalized and we test for this
!   In general, the old quaternion need not be normalized, and the same goes for the result
!   although in our applications we only ever use unit quaternions (to represent orientations)

        REAL(KIND=DP)                 :: NORM
        REAL(KIND=DP), DIMENSION(0:3) :: ROT

        NORM = SUM ( AXIS**2 ) ! Axis squared length
        
        IF ( ABS ( NORM - 1.0_dp ) > TOL ) THEN
            WRITE ( UNIT=ERROR_UNIT, fmt='(a,2es20.8)' ) 'axis normalization error', norm, tol
            STOP 'Error in rotate_quaternion'
        END IF

        ! Standard formula for rotation quaternion, using half angles
        ROT(0)   = COS(0.5_dp*ANGLE)
        ROT(1:3) = SIN(0.5_dp*ANGLE)*AXIS

        E = QUATMUL ( ROT, OLD ) ! Apply rotation to old quaternion

    END FUNCTION ROTATE_QUATERNION

    FUNCTION QUATMUL ( A, B ) RESULT ( C )
    ! Multiply 2.0_dp quaternions, A and B, together to give quaternion C
    ! See https://www.mathworks.com/help/aeroblks/quaternionmultiplication.html
        IMPLICIT NONE
        REAL(KIND=DP), DIMENSION(0:3)             :: C    ! Returns quaternion product of
        REAL(KIND=DP), DIMENSION(0:3), INTENT(in) :: A, B ! 2.0_dp supplied quaternions

        C(0) = A(0)*B(0) - A(1)*B(1) - A(2)*B(2) - A(3)*B(3)
        C(1) = A(1)*B(0) + A(0)*B(1) - A(3)*B(2) + A(2)*B(3)
        C(2) = A(2)*B(0) + A(3)*B(1) + A(0)*B(2) - A(1)*B(3)
        C(3) = A(3)*B(0) - A(2)*B(1) + A(1)*B(2) + A(0)*B(3)

    END FUNCTION QUATMUL

    FUNCTION Q_TO_RM ( Q ) RESULT ( RM )
        IMPLICIT NONE
        REAL(KIND=DP), DIMENSION(3,3)             :: RM ! Returns a 3x3 rotation matrix calculated from
        REAL(KIND=DP), DIMENSION(0:3), INTENT(in) :: Q ! supplied quaternion

        ! The rows of the rotation matrix correspond to unit vectors of the molecule in the space-fixed frame
        ! The third row  a(:,3) is "the" axis of the molecule, for uniaxial molecules
        ! Use a to convert space-fixed to body-fixed axes thus: db = matmul(a,ds)
        ! Use transpose of a to convert body-fixed to space-fixed axes thus: ds = matmul(db,a)

        ! The supplied quaternion should be normalized and we check for this

        REAL(KIND=DP) :: NORM

        NORM = SUM ( Q**2 ) ! Quaternion squared length
        
        IF ( ABS ( NORM - 1.0_dp ) > TOL ) THEN
            WRITE ( UNIT=ERROR_UNIT, fmt='(a,2es20.8)' ) 'quaternion normalization error', norm, tol
            STOP 'Error in q_to_a'
        END IF

        ! Write out column by column, for clarity

        RM(:,1) = [ Q(0)**2+Q(1)**2-Q(2)**2-Q(3)**2,   2.0_dp*( Q(1)*Q(2)+Q(0)*Q(3) ),     2.0_dp*( Q(1)*Q(3)-Q(0)*Q(2) )     ] ! 1st Column
        RM(:,2) = [  2.0_dp*( Q(1)*Q(2)-Q(0)*Q(3) ),    Q(0)**2-Q(1)**2+Q(2)**2-Q(3)**2,   2.0_dp*( Q(2)*Q(3)+Q(0)*Q(1) )     ] ! 2nd Column
        RM(:,3) = [  2.0_dp*( Q(1)*Q(3)+Q(0)*Q(2) ),     2.0_dp*( Q(2)*Q(3)-Q(0)*Q(1) ),    Q(0)**2-Q(1)**2-Q(2)**2+Q(3)**2   ] ! 3rd Column
    END FUNCTION Q_TO_RM

    !====================================================================================================

    FUNCTION Q_TO_UV ( Q ) RESULT ( UV )
        IMPLICIT NONE
        REAL(KIND=DP), DIMENSION(3)               :: UV ! Returns a unit vector calculated from
        REAL(KIND=DP), DIMENSION(0:3), INTENT(in) :: Q ! supplied quaternion

        ! We assume that the reference orientation of the particles is defined by a unit vector pointing along the 
        ! z-axis (i.e., UV(ref) = [0,0,1]). 
        ! Then the unit vector associated with a quaternion is simply the 3rd column of the rotation matrix defined 
        ! in the above FUNCTION Q_TO_RM.

        ! The supplied quaternion should be normalized and we check for this

        REAL(KIND=DP) :: NORM

        NORM = SUM ( Q**2 ) ! Quaternion squared length
        
        IF ( ABS ( NORM - 1.0_dp ) > TOL ) THEN
            WRITE ( UNIT=ERROR_UNIT, fmt='(a,2es20.8)' ) 'quaternion normalization error', norm, tol
            STOP 'Error in q_to_a'
        END IF

        ! Write out row by row, for clarity

        UV = [  2.0_dp*( Q(1)*Q(3)+Q(0)*Q(2) ), 2.0_dp*( Q(2)*Q(3)-Q(0)*Q(1) ), Q(0)**2-Q(1)**2-Q(2)**2+Q(3)**2   ] ! 3rd Column
    END FUNCTION Q_TO_UV

    FUNCTION UV_TO_Q ( UV ) RESULT ( Q )
    ! Calculates the quaternion required to convert one unit vector into another. Algorithm taken from:
    ! https://stackoverflow.com/questions/1171849/finding-quaternion-representing-the-rotation-from-one-vector-to-another
    ! We are looking to calculate the quaternion required to rotate a reference unit vector (which defines the 
    ! reference orientation of our particles) into a new unit vector which is the input to this function (i.e., UV).
        
        USE COMMONS, ONLY: REFSITE
        
        IMPLICIT NONE
        REAL(KIND=DP), DIMENSION(3), INTENT(in) :: UV 
        REAL(KIND=DP), DIMENSION(0:3)           :: Q

        REAL(KIND=DP)                           :: UREF(3), K, ANGLE

        Q    = 0.0_dp
        
        UREF = REFSITE(:,1)
        K     = SQRT( DOT_PRODUCT(UREF,UREF) * DOT_PRODUCT(UV,UV) )
        ANGLE = -UV(3)

        IF(ANGLE / K == -1.0_dp) THEN
            Q = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
        ELSE
            Q(0)   = ANGLE + K
            Q(1:3) = [ UV(2), -UV(1), 0.0_dp ]
        ENDIF

        Q = Q / NORM2(Q)

    END FUNCTION UV_TO_Q

END MODULE ROTATIONS_MODULE