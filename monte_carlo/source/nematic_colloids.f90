SUBROUTINE NEMATIC_COLLOIDS(ENERGY, RIJ, RIJSQ)

    USE COMMONS, ONLY: DP, NDIM
    USE COMMONS, ONLY: K4PI, KAPPA, AY, NC_QLM, LMIN, LMAX, LSTEP

    IMPLICIT NONE

    INTEGER       :: LI, MI, L, M

    REAL(KIND=DP), INTENT(IN) :: RIJ(NDIM), RIJSQ

    REAL(KIND=DP) :: DIST, RHAT(NDIM)
    REAL(KIND=DP) :: UY, COS_T, PLM(1,20), UEL, FCTRL, XLM

    REAL(KIND=DP), INTENT(OUT) :: ENERGY

    ENERGY   = 0.0_dp

    UEL = 0.0_dp

    DIST = SQRT(RIJSQ)
    RHAT = RIJ/DIST

!   Calculate Yukawa repulsion
    UY    = AY*EXP( -KAPPA*(DIST-1.0_DP) )/DIST

!   Calculate multipolar approximation for nematic colloids
!   Nematic director is taken to be [0,0,1]
    COS_T = RHAT(NDIM)

    CALL P_POLYNOMIAL_VALUE ( 1, LMAX+LMAX, [COS_T], PLM )

    LI = 0
    DO L = LMIN, LMAX, LSTEP
        MI = 0
        LI = LI + 1

        DO M = LMIN, LMAX, LSTEP
            MI = MI + 1

            XLM = K4PI * NC_QLM(LI) * NC_QLM(MI) * (-1)**M * FCTRL(L+M)
            UEL = UEL + XLM * PLM(1,L+M+1) / DIST**(L+M+1)   
        ENDDO
    ENDDO

    ENERGY = UY + UEL

END SUBROUTINE NEMATIC_COLLOIDS

SUBROUTINE DEF_NEMATIC_COLLOIDS()
    !----------------------------------------------------------------
    ! 
    !----------------------------------------------------------------
    
        USE COMMONS, ONLY: DP, PI
        USE COMMONS, ONLY: K_NC, K4PI, AY, NC_QLM
    
        IMPLICIT NONE

        K4PI = 4.0_dp * PI * K_NC
        AY   = 1.0_dp
        IF(.NOT. ALLOCATED(NC_QLM)) THEN
            ALLOCATE(NC_QLM(3))
            NC_QLM = [-0.017*0.5**3, -0.092*0.5**5, 0.003*0.5**7]
        ENDIF

        PRINT *, NC_QLM, K4PI

END SUBROUTINE DEF_NEMATIC_COLLOIDS

SUBROUTINE P_POLYNOMIAL_VALUE ( M, N, X, V )

    !*****************************************************************************80
    !
    !! P_POLYNOMIAL_VALUE EVALUATES THE LEGENDRE POLYNOMIALS P(N,X).
    !
    !  DISCUSSION:
    !
    !    P(N,1) = 1.
    !    P(N,-1) = (-1)^N.
    !    | P(N,X) | <= 1 IN [-1,1].
    !
    !    THE N ZEROES OF P(N,X) ARE THE ABSCISSAS USED FOR GAUSS-LEGENDRE
    !    QUADRATURE OF THE INTEGRAL OF A FUNCTION F(X) WITH WEIGHT FUNCTION 1
    !    OVER THE INTERVAL [-1,1].
    !
    !    THE LEGENDRE POLYNOMIALS ARE ORTHOGONAL UNDER THE INNER PRODUCT DEFINED
    !    AS INTEGRATION FROM -1 TO 1:
    !
    !      INTEGRAL ( -1 <= X <= 1 ) P(I,X) * P(J,X) DX 
    !        = 0 IF I =/= J
    !        = 2 / ( 2*I+1 ) IF I = J.
    !
    !    EXCEPT FOR P(0,X), THE INTEGRAL OF P(I,X) FROM -1 TO 1 IS 0.
    !
    !    A FUNCTION F(X) DEFINED ON [-1,1] MAY BE APPROXIMATED BY THE SERIES
    !      C0*P(0,X) + C1*P(1,X) + ... + CN*P(N,X)
    !    WHERE
    !      C(I) = (2*I+1)/(2) * INTEGRAL ( -1 <= X <= 1 ) F(X) P(I,X) DX.
    !
    !    THE FORMULA IS:
    !
    !      P(N,X) = (1/2^N) * SUM ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
    !
    !  DIFFERENTIAL EQUATION:
    !
    !    (1-X*X) * P(N,X)'' - 2 * X * P(N,X)' + N * (N+1) = 0
    !
    !  FIRST TERMS:
    !
    !    P( 0,X) =      1
    !    P( 1,X) =      1 X
    !    P( 2,X) = (    3 X^2 -       1)/2
    !    P( 3,X) = (    5 X^3 -     3 X)/2
    !    P( 4,X) = (   35 X^4 -    30 X^2 +     3)/8
    !    P( 5,X) = (   63 X^5 -    70 X^3 +    15 X)/8
    !    P( 6,X) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
    !    P( 7,X) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
    !    P( 8,X) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
    !    P( 9,X) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
    !    P(10,X) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
    !
    !  RECURSION:
    !
    !    P(0,X) = 1
    !    P(1,X) = X
    !    P(N,X) = ( (2*N-1)*X*P(N-1,X)-(N-1)*P(N-2,X) ) / N
    !
    !    P'(0,X) = 0
    !    P'(1,X) = 1
    !    P'(N,X) = ( (2*N-1)*(P(N-1,X)+X*P'(N-1,X)-(N-1)*P'(N-2,X) ) / N
    !
    !  LICENSING:
    !
    !    THIS CODE IS DISTRIBUTED UNDER THE GNU LGPL LICENSE. 
    !
    !  MODIFIED:
    !
    !    10 MARCH 2012
    !
    !  AUTHOR:
    !
    !    JOHN BURKARDT
    !
    !  REFERENCE:
    !
    !    MILTON ABRAMOWITZ, IRENE STEGUN,
    !    HANDBOOK OF MATHEMATICAL FUNCTIONS,
    !    NATIONAL BUREAU OF STANDARDS, 1964,
    !    ISBN: 0-486-61272-4,
    !    LC: QA47.A34.
    !
    !    DANIEL ZWILLINGER, EDITOR,
    !    CRC STANDARD MATHEMATICAL TABLES AND FORMULAE,
    !    30TH EDITION,
    !    CRC PRESS, 1996.
    !
    !  PARAMETERS:
    !
    !    INPUT, INTEGER ( KIND = 4 ) M, THE NUMBER OF EVALUATION POINTS.
    !
    !    INPUT, INTEGER ( KIND = 4 ) N, THE HIGHEST ORDER POLYNOMIAL TO EVALUATE.
    !    NOTE THAT POLYNOMIALS 0 THROUGH N WILL BE EVALUATED.
    !
    !    INPUT, REAL ( KIND = RK ) X(M), THE EVALUATION POINTS.
    !
    !    OUTPUT, REAL ( KIND = RK ) V(M,0:N), THE VALUES OF THE LEGENDRE POLYNOMIALS 
    !    OF ORDER 0 THROUGH N AT THE POINTS X.
    !
      IMPLICIT NONE
    
      INTEGER, PARAMETER :: RK = KIND ( 1.0D+00 )
    
      INTEGER ( KIND = 4 ) M
      INTEGER ( KIND = 4 ) N
    
      INTEGER ( KIND = 4 ) I
      REAL ( KIND = RK ) V(M,0:N)
      REAL ( KIND = RK ) X(M)
    
      IF ( N < 0 ) THEN
        RETURN
      END IF
    
      V(1:M,0) = 1.0D+00
    
      IF ( N < 1 ) THEN
        RETURN
      END IF
    
      V(1:M,1) = X(1:M)
     
      DO I = 2, N
     
        V(1:M,I) = ( REAL ( 2 * I - 1, KIND = RK ) * X(1:M) * V(1:M,I-1)   &
                   - REAL (     I - 1, KIND = RK ) *          V(1:M,I-2) ) &
                   / REAL (     I,     KIND = RK )
     
      END DO
     
      RETURN
END SUBROUTINE

FUNCTION FCTRL ( N )

    USE COMMONS, ONLY: DP

    IMPLICIT NONE

    REAL ( KIND = DP ) FCTRL
    INTEGER ( KIND = 4 ) I
    INTEGER ( KIND = 4 ) N

    FCTRL = 1.0D+00

    DO I = 1, N
        FCTRL = FCTRL * REAL ( I, KIND = DP )
    END DO

    RETURN
END