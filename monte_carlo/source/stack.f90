MODULE STACK
 
    PUBLIC
 
    ! DEFINE THE DATA-STRUCTURE TO HOLD THE DATA
    TYPE STACK_VAR
        INTEGER, ALLOCATABLE :: DATA(:)
        INTEGER              :: SIZE = 0
    END TYPE STACK_VAR
    
    ! SET THE SIZE OF ALLOCATED MEMORY BLOCKS
    INTEGER, PARAMETER, PRIVATE :: BLOCK_SIZE = 10
 
CONTAINS
 
    ! PUSH ----------------------------------------------------------------------
    SUBROUTINE PUSH(S, E)
        TYPE(STACK_VAR), INTENT(INOUT) :: S
        INTEGER, INTENT(IN)            :: E
        INTEGER, ALLOCATABLE :: WK(:)
        IF (.NOT. ALLOCATED(S%DATA)) THEN
        ! ALLOCATE SPACE IF NOT YET DONE
        ALLOCATE(S%DATA(BLOCK_SIZE))
    
        ELSEIF (S%SIZE == SIZE(S%DATA)) THEN
        ! GROW THE ALLOCATED SPACE
        ALLOCATE(WK(SIZE(S%DATA)+BLOCK_SIZE))
        WK(1:S%SIZE) = S%DATA
        CALL MOVE_ALLOC(WK,S%DATA)
    
        END IF
    
        ! STORE THE DATA IN THE STACK
        S%SIZE = S%SIZE + 1
        S%DATA(S%SIZE) = E
    END SUBROUTINE PUSH
 
    ! POP -----------------------------------------------------------------------
    FUNCTION POP(S)
        INTEGER :: POP
        TYPE(STACK_VAR), INTENT(INOUT) :: S
        IF (S%SIZE == 0 .OR. .NOT. ALLOCATED(S%DATA)) THEN
        POP = 0
        RETURN
        END IF
        POP = S%DATA(S%SIZE)
        S%SIZE = S%SIZE - 1
    END FUNCTION POP
 
    ! PEEK ----------------------------------------------------------------------
    INTEGER FUNCTION PEEK(S)
        TYPE(STACK_VAR), INTENT(INOUT) :: S
        IF (S%SIZE == 0 .OR. .NOT. ALLOCATED(S%DATA)) THEN
            PEEK = 0
            RETURN
        END IF
        PEEK = S%DATA(S%SIZE)
    END FUNCTION PEEK
    
    ! EMPTY ---------------------------------------------------------------------
    LOGICAL FUNCTION EMPTY(S)
        TYPE(STACK_VAR), INTENT(INOUT) :: S
        EMPTY = (S%SIZE == 0 .OR. .NOT. ALLOCATED(S%DATA))
    END FUNCTION EMPTY
 
END MODULE STACK