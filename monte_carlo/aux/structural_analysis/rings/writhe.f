C FILE: WRITHE.F
        SUBROUTINE CALC_WRITHE(R,N,M,WR)
C
C     INCREMENT THE FIRST ROW AND DECREMENT THE FIRST COLUMN OF A
C
        INTEGER          :: N,M,I,J
        DOUBLE PRECISION :: R(N,M),R1(M),R2(M),R3(M),R4(M),R1234(M)
        DOUBLE PRECISION :: N1(M),N2(M),N3(M),N4(M),R12(M),R34(M)
        DOUBLE PRECISION :: A,B,C,D,X,Y,Z
        DOUBLE PRECISION :: WR
Cf2py   intent(in) R
Cf2py   integer intent(hide),depend(R) :: n=shape(R,0), m=shape(R,1)
Cf2py   intent(out) WR
            WR = 0.D0
            DO I = 1,N
                
                IF(I==N-1) THEN
                    IN = N
                ELSE
                    IN = MOD(I+1,N)!+1
                ENDIF

                DO J = 1,N

                    IF(J==N-1) THEN
                        JN = N
                    ELSE
                        JN = MOD(J+1,N)!+1
                    ENDIF

                    IF(I==J .OR. I==JN .OR. IN==J .OR. IN==JN ) CYCLE

                    R1 = R(I,:)  - R(JN,:)
                    R2 = R(I,:)  - R(J,:)
                    R3 = R(IN,:) - R(JN,:)
                    R4 = R(IN,:) - R(J,:)

                    X = R2(2)*R1(3)-R2(3)*R1(2)
                    Y = R2(3)*R1(1)-R2(1)*R1(3)
                    Z = R2(1)*R1(2)-R2(2)*R1(1)
                    N1=[X,Y,Z]

                    X = R1(2)*R3(3)-R1(3)*R3(2)
                    Y = R1(3)*R3(1)-R1(1)*R3(3)
                    Z = R1(1)*R3(2)-R1(2)*R3(1)
                    N2 = [X,Y,Z]

                    X = R3(2)*R4(3)-R3(3)*R4(2)
                    Y = R3(3)*R4(1)-R3(1)*R4(3)
                    Z = R3(1)*R4(2)-R3(2)*R4(1)
                    N3 = [X,Y,Z]

                    X = R4(2)*R2(3)-R4(3)*R2(2)
                    Y = R4(3)*R2(1)-R4(1)*R2(3)
                    Z = R4(1)*R2(2)-R4(2)*R2(1)
                    N4 = [X,Y,Z]

                    N1 = N1 / NORM2(N1)
                    N2 = N2 / NORM2(N2)
                    N3 = N3 / NORM2(N3)
                    N4 = N4 / NORM2(N4)

                    A = DASIN( DOT_PRODUCT(N1,N2) )
                    B = DASIN( DOT_PRODUCT(N2,N3) )
                    C = DASIN( DOT_PRODUCT(N3,N4) )
                    D = DASIN( DOT_PRODUCT(N4,N1) )

                    R12 = R(I,:) - R(IN,:)
                    R34 = R(J,:) - R(JN,:)

                    X = R12(2)*R34(3)-R12(3)*R34(2)
                    Y = R12(3)*R34(1)-R12(1)*R34(3)
                    Z = R12(1)*R34(2)-R12(2)*R34(1)
                    R1234 = [X,Y,Z]

                    S = DSIGN(DABS(A+B+C+D),DOT_PRODUCT(R1234,R2))

                    WR = WR + S

                ENDDO
            ENDDO

            WR = WR / (4.D0*4.D0*DATAN(1.D0))
        END

        SUBROUTINE CALC_LINKING_NUMBER(R1,R2,NI,MI,NJ,MJ,LK)
C
C     INCREMENT THE FIRST ROW AND DECREMENT THE FIRST COLUMN OF A
C
        INTEGER          :: NI,MI,NJ,MJ,I,J
        DOUBLE PRECISION :: R1(NI,MI),R2(NJ,MJ)
        DOUBLE PRECISION :: A(MI),B(MI),C(MI),D(MI)
        DOUBLE PRECISION :: AN,BN,CN,DN,X,Y,Z,BC(MI),DA(MI)
        DOUBLE PRECISION :: N1,D1,N2,D2
        DOUBLE PRECISION :: LK
Cf2py   intent(in) R1, R2
Cf2py   integer intent(hide),depend(R1) :: ni=shape(R1,0), mi=shape(R1,1)
Cf2py   integer intent(hide),depend(R2) :: nj=shape(R2,0), mj=shape(R2,1)
Cf2py   intent(out) LK
            LK = 0.D0
            DO I = 1,NI
                
                IF(I==NI-1) THEN
                    IN = NI
                ELSE
                    IN = MOD(I+1,NI)!+1
                ENDIF

                DO J = 1,NJ

                    IF(J==NJ-1) THEN
                        JN = NJ
                    ELSE
                        JN = MOD(J+1,NJ)!+1
                    ENDIF

                    A = R2(J,:)  - R1(I,:)
                    B = R2(J,:)  - R1(IN,:)
                    C = R2(JN,:) - R1(IN,:)
                    D = R2(JN,:) - R1(I,:)

                    AN = NORM2(A)
                    BN = NORM2(B)
                    CN = NORM2(C)
                    DN = NORM2(D)
                    
                    X =B(2)*C(3)-B(3)*C(2)
                    Y =B(3)*C(1)-B(1)*C(3)
                    Z =B(1)*C(2)-B(2)*C(1)
                    BC=[X,Y,Z]
                    N1=DOT_PRODUCT(A,BC)
                    D1=AN*BN*CN+DOT_PRODUCT(A,B)*CN
                    D1=D1+DOT_PRODUCT(C,A)*BN+DOT_PRODUCT(B,C)*AN

                    X =D(2)*A(3)-D(3)*A(2)
                    Y =D(3)*A(1)-D(1)*A(3)
                    Z =D(1)*A(2)-D(2)*A(1)
                    DA=[X,Y,Z]
                    N2=DOT_PRODUCT(C,DA)
                    D2=CN*DN*AN+DOT_PRODUCT(C,D)*AN
                    D2=D2+DOT_PRODUCT(A,C)*DN+DOT_PRODUCT(D,A)*CN
                    
                    LK = LK + ( DATAN2(N1,D1)+DATAN2(N2,D2) )

                ENDDO
            ENDDO

            LK = LK / (4.D0*DATAN(1.D0))
        END

        SUBROUTINE PERC_CHECK(R,BOXL,DCUT,N,M,MAXD,RP,COM)
C
C     INCREMENT THE FIRST ROW AND DECREMENT THE FIRST COLUMN OF A
C
        INTEGER          :: N,M,I,NI,IN
        DOUBLE PRECISION :: DCUT,BOXL,R(N,M),RIJ(M),RD
        DOUBLE PRECISION :: MAXD,RP(N,M),COM(M)
Cf2py   intent(in) R, BOXL,DCUT
Cf2py   integer intent(hide),depend(R) :: n=shape(R,0), m=shape(R,1)
Cf2py   intent(out) MAXD, RP,COM
            MAXD    = 0.D0
            RP(1,:) = R(1,:)  
            DO I=2,N
                RIJ = RP(1,:) - R(I,:)
                RIJ = RIJ - BOXL*ANINT(RIJ/BOXL)
                RP(I,:) = RP(1,:) - RIJ
            ENDDO

            DO I=1,N
                IF(I==N-1) THEN
                    IN = NI
                ELSE
                    IN = MOD(I+1,N)
                ENDIF
                RIJ = RP(I,:) - RP(IN,:)
                RD  = DOT_PRODUCT(RIJ,RIJ)
                IF(RD > MAXD) MAXD = RD
            ENDDO

            COM = 0.D0 
            IF(MAXD <= DCUT) THEN
                DO I=1,N
                    COM = COM + RP(I,:)
                ENDDO
                COM = COM / DBLE(N)
            ENDIF

        END

        SUBROUTINE RING_IMAGE(C1,R2,C2,BOXL,N,M,RPJ)
C
C     INCREMENT THE FIRST ROW AND DECREMENT THE FIRST COLUMN OF A
C
        INTEGER          :: N,M,I
        DOUBLE PRECISION :: BOXL,C1(M),R2(N,M),C2(M),RIJ(M)
        DOUBLE PRECISION :: RPS(N,M),RPJ(N,M),CJ(M)
Cf2py   intent(in) C1,R2,C2,BOXL
Cf2py   integer intent(hide),depend(R2) :: n=shape(R2,0), m=shape(R2,1)
Cf2py   intent(out) RPJ
            
            CJ = C2

            DO I=1,N
                RPS(I,:) = CJ - R2(I,:)
            ENDDO

            RIJ = C1 - CJ
            RIJ = RIJ - BOXL*ANINT(RIJ/BOXL)
            CJ  = C1 - RIJ

            DO I=1,N
                RPJ(I,:) = CJ - RPS(I,:)
            ENDDO

        END
C END OF FILE WRITHE.F