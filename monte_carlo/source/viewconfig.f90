SUBROUTINE VIEWCONFIG(FILENAME)

    USE COMMONS, ONLY: VIEWUNIT, HST, YUKT, GLJT, KIHARAT, KFT, PGLJT, CPPT, ETPT, DMBLGLJT, KFRECT, HDMBLT, HTPRT

    IMPLICIT NONE

    CHARACTER (LEN=*) :: FILENAME

    OPEN(VIEWUNIT,FILE=FILENAME,STATUS='UNKNOWN',ACCESS='APPEND')

    IF(HST .OR. YUKT .OR. GLJT) THEN
        CALL VIEW_GLJ()
    ELSEIF(KIHARAT) THEN
        CALL VIEW_KIHARA()
    ELSEIF(KFT) THEN
        CALL VIEW_KF()
    ELSEIF(PGLJT) THEN
        CALL VIEW_PGLJ()
    ELSEIF(CPPT) THEN
        CALL VIEW_CPP()
    ELSEIF(ETPT .OR. HTPRT) THEN
        CALL VIEW_ETP()
    ELSEIF(DMBLGLJT) THEN
        CALL VIEW_DMBL_GLJ()
    ELSEIF(KFRECT) THEN
        CALL VIEW_KF_REC()
    ELSEIF(HDMBLT) THEN
        CALL VIEW_HDMBL()
    ENDIF

    CLOSE(UNIT=VIEWUNIT)

END SUBROUTINE VIEWCONFIG