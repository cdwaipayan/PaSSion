MODULE EWALD

    USE COMMONS, ONLY: DP, NPART, NDIM, R, Q, BOX, SPET, INDXP, RCUTSQ, RIGIDT
    USE COMMONS, ONLY: RBSITES, REFSITE, NSITES, VIRTEMP, CELLLISTT

    IMPLICIT NONE

    INTEGER, PRIVATE             :: I, J, A, B, C, IA, IB, IC, IK, JL, I1, J1
    REAL(KIND=DP), PRIVATE    :: DIST, RBOX_SQ, NBOX_SQ
    REAL(KIND=DP), PRIVATE    :: G, E_1, E_2, KDOTR, K, K_SQ
    REAL(KIND=DP), PRIVATE    :: QSUM1, QSUM2

    REAL(KIND=DP), DIMENSION(3), PRIVATE  :: RIJ, RHAT, D, KVEC
    REAL(KIND=DP), PARAMETER, PRIVATE     :: PI = 4.D0*DATAN(1.D0), TWOPI = 2.D0*PI, TWOPI_SQ = TWOPI**2
    REAL(KIND=DP), PARAMETER, PRIVATE     :: Y = 0.5772156649D0 !Euler's constant

    PUBLIC              :: BRUTE_EXP, EW_REAL, EW_RECIP, EW_SELF

CONTAINS

!=====================================================================================
!=====================================================================================
!     recursive subroutine to calculate the upper incomplete gamma function.
!     i_gamma is the returned value
!     n must be an integer or half integer
!     x must be a positive value
!=====================================================================================
!=====================================================================================
    RECURSIVE FUNCTION INC_GAMMA(N, X) RESULT(I_GAMMA)

        IMPLICIT NONE

        REAL(KIND=DP)              :: I_GAMMA
        REAL(KIND=DP), INTENT(in)  :: N, X
        REAL(KIND=DP)              :: R_GAMMA, SUM, COUNTER, I_FACT

        I_GAMMA = 0.D0
        R_GAMMA = 0.D0
!
!           gamma(0,x) equals the exponential integral (Ei(x))
!           Ei(x) = y + ln(x) + infinitesum(x**i/i*i!)
!           becomes infinite when x >= 80
!
        IF(N==0.D0) THEN
            SUM = 0.D0; COUNTER = -1.D0
            I = 1; I_FACT = 1.D0

            DO WHILE(SUM - COUNTER > 1.D-15)
                I_FACT = I_FACT*DBLE(I)
                COUNTER = SUM
                SUM = SUM + ((-1.D0)**I)*(X**I)/(DBLE(I)*I_FACT)
                I = I+1
            ENDDO

            I_GAMMA = -Y - DLOG(X) - SUM

!           gamma(0.5,x) = sqrt(pi)*erfc(sqrt(x))
        ELSE IF(N==0.5D0) THEN
            I_GAMMA = DSQRT(PI)*( DERFC( DSQRT(X) ) )

!           gamma(1,x) = exp(-x)
        ELSE IF(N==1.D0) THEN
            I_GAMMA = DEXP(-X)

!--------------------------------------------------------------------
!     Recursive part of the function is computed here
!--------------------------------------------------------------------
!           if n is positive use decreasing recurrence relationship
!           gamma(n,x) = (n-1)*gamma(n-1,x) + x**(n-1)*exp(-x)
        ELSE IF (N > 0.D0) THEN
            R_GAMMA = INC_GAMMA(N-1.D0, X)
            I_GAMMA = (N - 1.D0)*R_GAMMA + (X**(N-1))*DEXP(-X)

!           if n is negative use increasing recurrence relationship
!           gamma(n,x) = {gamma(n+1,x) - [x**(n)*exp(-x)]}/n
        ELSE
            R_GAMMA = INC_GAMMA(N+1.D0, X)
            I_GAMMA = (R_GAMMA - (X**N)*DEXP(-X)) / N
        ENDIF
!--------------------------------------------------------------------
!--------------------------------------------------------------------
        RETURN

    END FUNCTION


    SUBROUTINE EW_REAL(N, NBOX, ALPHA, RC, BIJ, E_REAL)
!=============================================================================================
!     Calculates real component of general ewald sum and forces using the equation:
!     
!     INPUTS:
!     - N:      Exponent of distance vector (r^-n)
!     - NBOX:   Array containing number of images to be considered in each direction
!     - ALPHA:  Ewald damping parameter
!     - RC:     Real-space cut-off value
!
!     OUTPUTS:
!     - E_REAL: Real space contribution to energy
!=============================================================================================

        IMPLICIT NONE

        INTEGER, INTENT(in)           :: NBOX(3), N
        REAL(KIND=DP), INTENT(in)     :: ALPHA, RC, BIJ
        REAL(KIND=DP), INTENT(out)    :: E_REAL  
        REAL(KIND=DP)                 :: G1, GA, DBN
        REAL(KIND=DP), DIMENSION(3)   :: NVEC

        E_REAL = 0.D0

        DBN = DBLE(N)

        NBOX_SQ = INT(RC)**2
        IA = NBOX(1); IB = NBOX(2); IC = NBOX(3)

    !--------------------------------------------------------------------
    !  in order to calculate energy and forces in the same loop the outer
    !  loop is over the images of the unit cell to be considered.
    !--------------------------------------------------------------------

        DO A = -IA, IA
            DO B = -IB, IB
                DO C = -IC, IC
                    
                    RBOX_SQ = NVEC(1)**2 + NVEC(2)**2 + NVEC(3)**2
                    IF ( RBOX_SQ > NBOX_SQ ) CYCLE ! skip if outside maximum sphere of boxes           
                    
                    IF ( I==J .AND. RBOX_SQ==0.0_dp) CYCLE ! ignore |RIJ + D| = 0 terms							                           
                    !Calculate Rj = RJ + Rj0*RMJ
                    
                    RIJ = R(:,I1) - R(:,J1) - D   
                    DIST = NORM2(RIJ)

                    IF(DIST.GT.RC) CYCLE

                    GA = (ALPHA*DIST)**2

                    G1 = INC_GAMMA(DBN/2.D0, GA)

                    E_REAL = E_REAL + ( BIJ * G1 ) / DIST**N 

                END DO 
            END DO
        END DO

        E_REAL = E_REAL / (2.D0 * DGAMMA(DBN/2.D0) ) 

        RETURN
    END SUBROUTINE


    SUBROUTINE EW_RECIP(N, NBOX, R, DR1, DR2, DR3, ALPHA, P, V, KC, GTEST, E_RECIP, F_RECIP, F_CELL)

        IMPLICIT NONE

        INTEGER, INTENT(in)                                        :: NBOX(3), N
        REAL(KIND=DP), INTENT(in)                               :: ALPHA, KC, V
        REAL(KIND=DP), DIMENSION(3,NRBSITES*NMOL), INTENT(in)   :: R 
        REAL(KIND=DP), DIMENSION(6), INTENT(in)                 :: P   
        LOGICAL, INTENT(in)                                        :: GTEST
        REAL(KIND=DP), DIMENSION(3,NRBSITES*NMOL), INTENT(in)   :: DR1, DR2, DR3 

        REAL(KIND=DP), INTENT(out)                              :: E_RECIP
        REAL(KIND=DP), DIMENSION(3,NATOMS), INTENT(out)         :: F_RECIP
        REAL(KIND=DP), DIMENSION(6), INTENT(out)                :: F_CELL
        
        
        REAL(KIND=DP)                                    :: KR, RECI, G1, G2, AI, AA, A2, F1, DBN
        REAL(KIND=DP), DIMENSION(3)                      :: F2, F3, F4, F5, F6, KAP
        REAL(KIND=DP), DIMENSION(3)                      :: DUKL1, DUKL2, DUKL3, DUKA1, DUKA2, DUKA3
        REAL(KIND=DP), DIMENSION(3,3)                    :: L1, U1  
        REAL(KIND=DP), DIMENSION(3,3)                    :: DUDL1, DUDL2, DUDL3, DUDA1, DUDA2, DUDA3
        REAL(KIND=DP)                                    :: C1, C2, C3, S1, S2, S3, W, KN3

        DBN = DBLE(N)

        G       = 0.D0
        G1      = (3.D0-DBN)/2.D0

        AA      = ALPHA**2
        A2      = 2.D0*ALPHA   
        AI      = 1.D0/(2.D0*AA)
        
        E_RECIP = 0.D0
        F_RECIP = 0.D0
        F_CELL  = 0.D0

        NBOX_SQ = KC**2
        IA = NBOX(1); IB = NBOX(2); IC = NBOX(3)

        CALL CELL_ARRAYS(P, L1, U1)
        CALL ICELL_DER(P, DUDL1, DUDL2, DUDL3, DUDA1, DUDA2, DUDA3, W, C1, C2, C3, S1, S2, S3)

        DO A = -IA, IA
            
            KAP(1) = DBLE(A)
            KVEC(1) = TWOPI*(KAP(1)*U1(1,1))
            
            DO B = -IB, IB
                
                KAP(2) = DBLE(B)
                KVEC(2) = TWOPI*(KAP(1)*U1(2,1) + KAP(2)*U1(2,2))
                
                DO C = -IC, IC
                    
                    KAP(3) = DBLE(C)
                    KVEC(3) = TWOPI*(KAP(1)*U1(3,1) + KAP(2)*U1(3,2) + KAP(3)*U1(3,3))
					
                    K_SQ = KVEC(1)**2 + KVEC(2)**2 + KVEC(3)**2
                    
                    IF (K_SQ==0 .OR. K_SQ.GT.((KC**2)+2)) CYCLE

                    K = DSQRT(K_SQ)
                    KN3 = K**(N - 3)

                    F4 = 0.D0

                    E_1 = 0.D0
                    E_2 = 0.D0

                    G2 = K_SQ / (4.D0*AA)
                    G = INC_GAMMA( G1, G2)

                    F5 = (( (DBN-3.D0)*KVEC ) / K_SQ) * G
                    F2 = DEXP(-G2)*KVEC*AI*(K/A2)**(1-N)
            
                    DO I = 1, NMOL
                        DO IK = 1, NRBSITES

                            IF(QB(IK)==0.D0) CYCLE

                            I1 = ((I-1)*NRBSITES) + IK
                        
                            F1 = 0.D0

                            KDOTR = DOT_PRODUCT(KVEC,R(:,I1))
                            
                            E_1 = E_1 + QB(IK)*DCOS(KDOTR)
                            E_2 = E_2 + QB(IK)*DSIN(KDOTR)

                            IF(GTEST) THEN
                                
                                DO J = 1, NMOL        
                                    DO JL = 1, NRBSITES     
                                        
                                        IF(QB(JL)==0.D0) CYCLE

                                        J1 = ((J-1)*NRBSITES) + JL  
                                        
                                        RIJ = R(:,I1) - R(:,J1)
                                        KR = DOT_PRODUCT(KVEC,RIJ)

                                        F1 = F1 + QB(IK)*QB(JL)*DSIN(KR)
                                        F3 = DSIN(KR)*RIJ*G
                                        F4 = F4 + QB(IK)*QB(JL)*KN3*(DCOS(KR)*(F5-F2) - F3)
                                    ENDDO
                                ENDDO
                                
                                F6 = KVEC * KN3 * G * F1
                                
                                F_RECIP(:,I) = F_RECIP(:,I) + F6
                                IF(RIGIDT) THEN
                                    F_RECIP(1,I+NMOL) = F_RECIP(1,I+NMOL) + DOT_PRODUCT(F6, DR1(:,I1))
                                    F_RECIP(2,I+NMOL) = F_RECIP(2,I+NMOL) + DOT_PRODUCT(F6, DR2(:,I1))
                                    F_RECIP(3,I+NMOL) = F_RECIP(3,I+NMOL) + DOT_PRODUCT(F6, DR3(:,I1))
                                ENDIF
                            ENDIF
                        ENDDO
                    ENDDO

                    IF(GTEST) THEN

                        DUKL1 = MATMUL(DUDL1, KAP)
                        DUKL2 = MATMUL(DUDL2, KAP)
                        DUKL3 = MATMUL(DUDL3, KAP)
                        DUKA1 = MATMUL(DUDA1, KAP)
                        DUKA2 = MATMUL(DUDA2, KAP)
                        DUKA3 = MATMUL(DUDA3, KAP)

                        F_CELL(1) = F_CELL(1) + DOT_PRODUCT( F4, DUKL1 )
                        F_CELL(2) = F_CELL(2) + DOT_PRODUCT( F4, DUKL2 )
                        F_CELL(3) = F_CELL(3) + DOT_PRODUCT( F4, DUKL3 )
                        F_CELL(4) = F_CELL(4) + DOT_PRODUCT( F4, DUKA1 )
                        F_CELL(5) = F_CELL(5) + DOT_PRODUCT( F4, DUKA2 )
                        F_CELL(6) = F_CELL(6) + DOT_PRODUCT( F4, DUKA3 )
                    ENDIF

                    E_RECIP = E_RECIP + K**(N - 3) * G * ( (E_1**2) + (E_2**2) )
                END DO
            END DO
        END DO

        RECI = ( (PI**1.5D0) / ( 2.D0**(N-2) * DGAMMA(DBN/2.D0) * V ))
        
        E_RECIP =  RECI * E_RECIP
        F_RECIP = -(RECI*2.D0) * F_RECIP
        F_CELL  =  TWOPI * RECI * F_CELL

        F_CELL(1) = F_CELL(1) - (E_RECIP / P(1))
        F_CELL(2) = F_CELL(2) - (E_RECIP / P(2))
        F_CELL(3) = F_CELL(3) - (E_RECIP / P(3))

        F_CELL(4) = F_CELL(4) - (E_RECIP / W) * (S1 * (C1 - C2*C3))
        F_CELL(5) = F_CELL(5) - (E_RECIP / W) * (S2 * (C2 - C1*C3))
        F_CELL(6) = F_CELL(6) - (E_RECIP / W) * (S3 * (C3 - C1*C2))

        RETURN

    END SUBROUTINE

!=====================================================================================
!=====================================================================================
!     recursive subroutine to calculate the upper incomplete gamma function.
!     i_gamma is the returned value
!     n must be an integer or half integer
!     x must be a positive value
!=====================================================================================
!=====================================================================================
    SUBROUTINE EW_SELF(N, P, ALPHA, V, R, E_SELF, C_SELF, GTEST) 

        IMPLICIT NONE

        INTEGER, INTENT(in)                                        :: N
        REAL(KIND=DP), INTENT(in)                               :: ALPHA, V
        REAL(KIND=DP), DIMENSION(6), INTENT(in)                 :: P
        REAL(KIND=DP), DIMENSION(3,NRBSITES*NMOL), INTENT(in)   :: R
        LOGICAL, INTENT(in)                                        :: GTEST

        REAL(KIND=DP), INTENT(out)                              :: E_SELF
        REAL(KIND=DP), DIMENSION(6), INTENT(out)                :: C_SELF

        REAL(KIND=DP)                                :: C1, C2, C3, S1, S2, S3, W, MSELF, GAMA, M1, AD, DBN

        E_SELF = 0.D0
        C_SELF = 0.D0
        QSUM1  = 0.D0
        QSUM2  = 0.D0

        DBN = DBLE(N)
	
        DO I = 1, NMOL
            DO IK = 1, NRBSITES

                QSUM1 = QSUM1 + QB(IK)*QB(IK)

                DO J = 1, NMOL
                    DO JL = 1, NRBSITES
                        QSUM2 = QSUM2 + QB(IK)*QB(JL)
                    ENDDO
                ENDDO
            
            ENDDO
        ENDDO

        MSELF = 0.D0
        IF(RIGIDT) THEN
            GAMA = DGAMMA(DBN/2.D0)
            
            DO I = 1, NMOL
                DO IK = 1, NRBSITES
                    I1 = ((I-1)*NRBSITES) + IK       
                    DO JL = 1, NRBSITES
                        J1 = ((I-1)*NRBSITES) + JL
                        
                        IF(IK == JL) CYCLE
                        RIJ = R(:,I1) - R(:,J1)   
                        DIST = NORM2(RIJ)
                        AD = ALPHA*DIST
                        M1 = INC_GAMMA(DBN/2.D0, AD**2)
                        MSELF = MSELF + QB(IK)*QB(JL)*(1.D0 - M1/GAMA)/(2.D0*DIST**N)
                    ENDDO
                ENDDO
            ENDDO

        ENDIF
		
        E_SELF = ( (PI**1.5D0)  * (ALPHA**(N-3)) )  / ( (DBN-3.D0) * GAMMA(DBN/2.D0) * V ) * QSUM2


        IF(GTEST) THEN
            C_SELF = E_SELF  
        
            C1 = DCOS(P(4)); C2 = DCOS(P(5)); C3 = DCOS(P(6)); 
            S1 = DSIN(P(4)); S2 = DSIN(P(5)); S3 = DSIN(P(6)); 

            W = 1.D0 - (C1**2 + C2**2 + C3**2) + 2.D0*(C1*C2*C3)  
        
            C_SELF(1) = -( C_SELF(1) / P(1) )
            C_SELF(2) = -( C_SELF(2) / P(2) )
            C_SELF(3) = -( C_SELF(3) / P(3) )
            C_SELF(4) = -( C_SELF(4) / W ) * ( S1 * (C1 - C2*C3) )
            C_SELF(5) = -( C_SELF(5) / W ) * ( S2 * (C2 - C1*C3) )
            C_SELF(6) = -( C_SELF(6) / W ) * ( S3 * (C3 - C1*C2) )
        ENDIF

        E_SELF = E_SELF - ( ALPHA**N / (DBN*DGAMMA(DBN/2.D0)) ) * QSUM1

        IF(RIGIDT) THEN
            E_SELF = E_SELF - MSELF
        ENDIF

        RETURN

    END SUBROUTINE

END MODULE