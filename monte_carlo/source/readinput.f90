SUBROUTINE READITEM(WORD)

    USE COMMONS, ONLY: LNGTH, LOC, NITEM, CHAR

    IMPLICIT NONE

    INTEGER :: I, ITEM, L, LAST
    LOGICAL :: BLANK, TCOMMA
    CHARACTER :: SPACE, COMMA, SQUOTE, DQUOTE
    CHARACTER (LEN = 16) :: WORD
    DATA SPACE /' '/, COMMA /','/, SQUOTE /''''/, DQUOTE /'"'/

    LAST = 80
!     Find last non-blank character
30  IF (CHAR(LAST) == SPACE) THEN
        LAST = LAST - 1
        IF (LAST > 1) GOTO 30
    ENDIF
!     Analyse input
    ITEM  = 1
    NITEM = 1
    L     = 1
    LNGTH(ITEM) = 0
    LOC(ITEM) = L
    TCOMMA = .FALSE.
    BLANK  = .FALSE.

    DO WHILE (L <= LAST)

        IF (BLANK) THEN
            BLANK = .FALSE.
            LOC(ITEM) = L
        ENDIF
            TCOMMA = (CHAR(L) == COMMA)
        IF (TCOMMA) THEN
            WRITE(*, *) 'COMMA found in the input, but should not be used' 
            STOP
        ENDIF
        BLANK = (CHAR(L) == SPACE)
        IF (BLANK) THEN
            IF (L == 1) THEN
                WRITE (*, *) 'cannot have space at the beginning of an input line'
                STOP
            ENDIF
            IF (ITEM == 1) THEN
                WORD = CHAR(1)
!      TRIM(string) returns the string with the trailing blanks removed
                IF (L > 2) THEN
                    DO I = 2, L-1
                        WORD = TRIM(WORD)//CHAR(I)             ! the operator // concatenates two or more strings
                    ENDDO
                ENDIF
            ENDIF
            IF (CHAR(L-1) == SPACE) THEN
                WRITE(*, *) 'more than one space has been used between items' 
                STOP
            ENDIF 
!     Item found
            ITEM = ITEM + 1
            NITEM = NITEM + 1
!     Looking for next item
            LNGTH(ITEM) = 0
        ELSE
            LNGTH(ITEM) = LNGTH(ITEM) + 1
        ENDIF
        L = L + 1
    ENDDO

    IF (NITEM == 1) THEN
        WORD = CHAR(1)
        IF (L > 2) THEN
            DO I = 2, L-1
                WORD = TRIM(WORD)//CHAR(I) 
            ENDDO
        ENDIF
    ENDIF

END SUBROUTINE READITEM

!     ----------------------------------------------------------------------------------------------

SUBROUTINE READI(ITEM, INTGR)

    USE COMMONS, ONLY: LNGTH, LOC, CHAR

    IMPLICIT NONE

    INTEGER :: INTGR, ITEM, J
    CHARACTER (LEN=16) :: SUBSTR


    SUBSTR = CHAR(LOC(ITEM))
    IF (LNGTH(ITEM) /=1) THEN
        DO J = LOC(ITEM) + 1, LOC(ITEM) + LNGTH(ITEM) - 1 
            SUBSTR = TRIM(SUBSTR)//CHAR(J)        
        ENDDO
    ENDIF
    READ (SUBSTR, *) INTGR

END SUBROUTINE READI

!     ----------------------------------------------------------------------------------------------

SUBROUTINE READF(ITEM, DPNMBR)

    USE COMMONS, ONLY: LNGTH, LOC, CHAR, DP

    IMPLICIT NONE

    INTEGER :: ITEM, J
    CHARACTER (LEN=16) :: SUBSTR
    REAL(KIND=DP) :: DPNMBR

    SUBSTR = CHAR(LOC(ITEM))
    IF (LNGTH(ITEM) /=1) THEN
        DO J = LOC(ITEM) + 1, LOC(ITEM) + LNGTH(ITEM) - 1 
            SUBSTR = TRIM(SUBSTR)//CHAR(J)        
        ENDDO
    ENDIF
    READ (SUBSTR, *) DPNMBR

END SUBROUTINE READF
