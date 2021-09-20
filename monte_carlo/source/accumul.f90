
SUBROUTINE ACCUMULATORS ()
!     This subroutine updates various accumulators for average calculations.
    USE COMMONS
    USE ORDERPARAM, ONLY: GET_BOPS
    
    IMPLICIT NONE
    INTEGER       :: J1
    REAL(KIND=DP) :: CORP, INVBLK, INVBLK2, BENP, NP, BLK, DE, B2, H, DH, DV, DVH, NV
    REAL(KIND=DP) :: CV, CP, KT, AL, STDCV, STDCP, STDKT, STDAL
    LOGICAL       :: FILE_EXISTS

    IF (ISTEP == 0) THEN
    !----------------------------------------------------------------------------------   
    !   Observables for standard NVT or NPT Monte Carlo simulations
    !----------------------------------------------------------------------------------
        ICNT    = 0   
        INQUIRE(FILE="../finalave.dat", EXIST=FILE_EXISTS)
        
        IF(CONTINUET .AND. FILE_EXISTS) THEN
            OPEN (UNIT = 41, FILE = '../finalave.dat', STATUS = 'OLD')
            READ (41,*) IBLOCK
            READ (41,*) SUMPE, STDPESUM
            READ (41,*) SUMPE2, STDPE2SUM
            READ (41,*) SUMPRES, STDPRESSUM
            IF(NPTT) THEN
                READ (41,*) SUMVLM, STDVLMSUM
                READ (41,*) SUMVLM2, STDVLM2SUM
                READ (41,*) SUMH, STDHSUM
                READ (41,*) SUMH2, STDH2SUM
                READ (41,*) SUMHV, STDRHOSUM
                READ (41,*) SUMRHO, STDHSUM
                READ (41,*) SUMBOX, STDBOXSUM
            ENDIF
            IF(PINT) THEN
                READ (41,*) SUMTRQ, STDTRQSUM 
                READ (41,*) SUMDELMU, STDDELMUSUM
            ENDIF
            IF (NUCSEEDT) THEN
                READ (41,*) SUMLRGSTXCLSTR, STDLRGSTXCLSTRSUM
                READ (41,*) SUMNUCCOUNT, STDNUCCOUNTSUM
                READ (41,*) NUCCOUNT, TOTNUCCOUNT
            ENDIF
            IF(FLFET) THEN
                IF(FLID==1) THEN
                    READ (41,*) SUM_EXP_U_EIN, STD_EXP_U_EIN_SUM
                ELSE
                    READ (41,*) SUM_U_EIN_TR, STD_U_EIN_TR_SUM
                    READ (41,*) SUM_U_EIN_OR, STD_U_EIN_OR_SUM
                ENDIF
            ENDIF
            IF(SCHSMIT) THEN
                IF(SSID==1) THEN
                    READ (41,*) SUM_EXP_U_SS, STD_EXP_U_SS_SUM
                ELSE
                    READ (41,*) SUM_U_SS, STD_U_SS_SUM
                ENDIF
            ENDIF
            CLOSE(UNIT = 41)
        ELSE

            IBLOCK  = 0
            SUMPE   = 0.0_dp; STDPESUM   = 0.0_dp
            SUMPE2  = 0.0_dp; STDPE2SUM  = 0.0_dp
            SUMPRES = 0.0_dp; STDPRESSUM = 0.0_dp
            SUMVLM  = 0.0_dp; STDVLMSUM  = 0.0_dp
            SUMVLM2 = 0.0_dp; STDVLM2SUM = 0.0_dp
            SUMH    = 0.0_dp; STDHSUM    = 0.0_dp
            SUMH2   = 0.0_dp; STDH2SUM   = 0.0_dp
            SUMHV   = 0.0_dp; STDHVSUM   = 0.0_dp
            SUMRHO  = 0.0_dp; STDRHOSUM  = 0.0_dp
            SUMBOX  = 0.0_dp; STDBOXSUM  = 0.0_dp

        !----------------------------------------------------------------------------------
        !   Chemical potential and order parameter for interface-pinning simulations
        !----------------------------------------------------------------------------------
            SUMTRQ   = 0.0_dp; STDTRQSUM   = 0.0_dp
            SUMDELMU = 0.0_dp; STDDELMUSUM = 0.0_dp

        !----------------------------------------------------------------------------------
        !   Count of largest crystalline cluster for umbrella sampling and nucleus-size
        !   pinning simulations.
        !----------------------------------------------------------------------------------
            SUMLRGSTXCLSTR = 0.0_dp; STDLRGSTXCLSTRSUM = 0.0_dp
            IF (NUCSEEDT) THEN
                NUCCOUNT    = 0 
                SUMNUCCOUNT = 0; STDNUCCOUNTSUM    = 0.0_dp
                TOTNUCCOUNT = 0
            ENDIF

        !----------------------------------------------------------------------------------
        !   Free energies from Einstein-crystal integration simulations
        !----------------------------------------------------------------------------------
            SUM_EXP_U_EIN  = 0.0_dp; STD_EXP_U_EIN_SUM  = 0.0_dp
            SUM_U_EIN_TR   = 0.0_dp; STD_U_EIN_TR_SUM   = 0.0_dp
            SUM_U_EIN_OR   = 0.0_dp; STD_U_EIN_OR_SUM   = 0.0_dp
        
        !----------------------------------------------------------------------------------
        !   Free energies from Schilling-Schmid integration simulations
        !----------------------------------------------------------------------------------
            SUM_U_SS     = 0.0_dp; STD_U_SS_SUM     = 0.0_dp 
            SUM_EXP_U_SS = 0.0_dp; STD_EXP_U_SS_SUM = 0.0_dp 
        ENDIF

    !   Initialise counters for block averages
        SUMPEBLK          = 0.0_dp
        SUMPE2BLK         = 0.0_dp

        SUMPRESBLK        = 0.0_dp
        SUMRHOBLK         = 0.0_dp
        SUMVLMBLK         = 0.0_dp
        SUMVLM2BLK        = 0.0_dp
        SUMHBLK           = 0.0_dp
        SUMH2BLK          = 0.0_dp
        SUMHVBLK          = 0.0_dp
        SUMBOXBLK         = 0.0_dp
        
        SUMTRQBLK         = 0.0_dp
        SUMLRGSTXCLSTRBLK = 0.0_dp
        
        SUM_EXP_U_EIN_BLK = 0.0_dp
        SUM_U_EIN_TR_BLK  = 0.0_dp
        SUM_U_EIN_OR_BLK  = 0.0_dp

        SUM_U_SS_BLK     = 0.0_dp
        SUM_EXP_U_SS_BLK = 0.0_dp

        RETURN
        
    ENDIF

!----------------------------------------------------------------------------------
!   Update the running sum for each of the system observables for the current block
!----------------------------------------------------------------------------------
!   Potential energy per particle
    PEPP        = PE/REAL(NPART,DP)
    SUMPEBLK    = SUMPEBLK  + PE
    SUMPE2BLK   = SUMPE2BLK + PE*PE
!   Volume
    SUMVLMBLK   = SUMVLMBLK + VLM
    SUMVLM2BLK  = SUMVLM2BLK + VLM*VLM
!   Enthalpy
    H           = (PE+PRSFIX*VLM)
    SUMHBLK     = SUMHBLK  + H
    SUMH2BLK    = SUMH2BLK + H**2
    SUMHVBLK    = SUMHVBLK + H*VLM
!   Density / packing fraction
    IF (DENSITYT) THEN
        RHO   = REAL(NPART,DP)/VLM
    ELSE IF (PACKINGT) THEN
        RHO   = REAL(NPART,DP)*PI/(6.0_dp*VLM)
    ENDIF
    SUMRHOBLK   = SUMRHOBLK + RHO
!   Virial pressure
    PRES     = (REAL(NPART,DP)/BETAKB + VIR)/VLM
    IF (TAILCORT) THEN
        IF (PACKINGT) STOP ' PACKINGT cannot be TRUE'
        PRES = PRES + CORP(RCUT,RHO)
    ENDIF
    SUMPRESBLK  = SUMPRESBLK + PRES
!   Box lengths
    IF(NPTT) SUMBOXBLK = SUMBOXBLK + BOX
!   Translation order-parameter for interface pinning
    IF (NPZTT) THEN
        IF(.NOT. PINT) THEN
            IF(COLLDT) THEN
                CALL COLLDENSFIELD(TRQ,RQ,-1)
            ENDIF
        ENDIF
        SUMTRQBLK  = SUMTRQBLK + TRQ
    ENDIF

    IF (FLFET) THEN
        SUM_EXP_U_EIN_BLK = SUM_EXP_U_EIN_BLK + EXP_U_EIN
        SUM_U_EIN_TR_BLK  = SUM_U_EIN_TR_BLK + U_EIN_TR
        SUM_U_EIN_OR_BLK  = SUM_U_EIN_OR_BLK + ITA_EIN*U_EIN_OR
    ENDIF

    IF (SCHSMIT) THEN
        SUM_EXP_U_SS_BLK  = SUM_EXP_U_SS_BLK + EXP_U_SS
        SUM_U_SS_BLK      = SUM_U_SS_BLK + U_SS
        SUM_U_EIN_OR_BLK  = SUM_U_EIN_OR_BLK + ITA_EIN*U_EIN_OR
    ENDIF

!   Check if we have reached the end of a block
    IF ( MOD(ISTEP,BLKLNGTH) == 0 ) THEN
    !   Update the index of the current block
        IBLOCK = IBLOCK + 1

        IF( NUCSEEDT ) THEN
            BLK = REAL(BLKLNGTH/TRJCTYLNGTH,DP)
        ELSE
            BLK = REAL(BLKLNGTH,DP)
        ENDIF
    !--------------------------------------------------------------------------------
    !   Calculate the average values for the system observables for the current block
    !--------------------------------------------------------------------------------
        AVPEBLK    = SUMPEBLK/BLK
        AVPE2BLK   = SUMPE2BLK/BLK
        AVPRESBLK  = SUMPRESBLK/BLK
        AVRHOBLK   = SUMRHOBLK/BLK
        AVVLMBLK   = SUMVLMBLK/BLK 
        AVVLM2BLK  = SUMVLM2BLK/BLK 
        AVHBLK     = SUMHBLK/BLK 
        AVH2BLK    = SUMH2BLK/BLK 
        AVHVBLK    = SUMHVBLK/BLK 
        IF (NPTT) AVBOXBLK = SUMBOXBLK/BLK
        IF (NPZTT .OR. PINT) THEN
            AVTRQBLK  = SUMTRQBLK/BLK
        !   Calculate the chemical potential difference for the current block
        !   based on the average value for the collective density field.
            IF(PINT) DELMUBLK = -(PINKAP*PINDELQ)/REAL(NPART,DP) * (AVTRQBLK - PINA)
        ENDIF

        IF (NUCSEEDT) THEN
            AVLRGSTXCLSTRBLK = SUMLRGSTXCLSTRBLK/BLK
            AVNUCCOUNTBLK    = REAL(NUCCOUNT,DP)/BLK
        ENDIF    
        
        IF (FLFET) THEN
            AV_EXP_U_EIN_BLK = SUM_EXP_U_EIN_BLK/BLK 
            AV_U_EIN_TR_BLK  = SUM_U_EIN_TR_BLK/BLK 
            AV_U_EIN_OR_BLK  = SUM_U_EIN_OR_BLK/BLK 
        ENDIF

        IF (SCHSMIT) THEN
            AV_EXP_U_SS_BLK  = SUM_EXP_U_SS_BLK/BLK 
            AV_U_SS_BLK      = SUM_U_SS_BLK/BLK 
            AV_U_EIN_OR_BLK  = SUM_U_EIN_OR_BLK/BLK 
        ENDIF

    !------------------------------------------------------------------
    !   Add the computed averages to the running sum of average values. 
    !------------------------------------------------------------------
        SUMPE   = SUMPE + AVPEBLK
        SUMPE2  = SUMPE2 + AVPE2BLK
        SUMPRES = SUMPRES + AVPRESBLK
        SUMRHO  = SUMRHO + AVRHOBLK
        SUMVLM  = SUMVLM + AVVLMBLK
        SUMVLM2 = SUMVLM2 + AVVLM2BLK
        SUMH    = SUMH + AVHBLK
        SUMH2   = SUMH2 + AVH2BLK
        SUMHV   = SUMHV + AVHVBLK
        IF(NPTT) SUMBOX = SUMBOX + AVBOXBLK
        IF (NPZTT .OR. PINT) THEN
            SUMTRQ  = SUMTRQ + AVTRQBLK
            IF(PINT) SUMDELMU = SUMDELMU + DELMUBLK
        ENDIF
        
        IF (NUCSEEDT) THEN
            SUMLRGSTXCLSTR = SUMLRGSTXCLSTR + AVLRGSTXCLSTRBLK
            SUMNUCCOUNT    = SUMNUCCOUNT + AVNUCCOUNTBLK
        ENDIF

        IF (FLFET) THEN
            SUM_EXP_U_EIN = SUM_EXP_U_EIN + AV_EXP_U_EIN_BLK
            SUM_U_EIN_TR  = SUM_U_EIN_TR + AV_U_EIN_TR_BLK
            SUM_U_EIN_OR  = SUM_U_EIN_OR + AV_U_EIN_OR_BLK
        ENDIF

        IF (SCHSMIT) THEN
            SUM_EXP_U_SS  = SUM_EXP_U_SS + AV_EXP_U_SS_BLK
            SUM_U_SS      = SUM_U_SS + AV_U_SS_BLK
            SUM_U_EIN_OR  = SUM_U_EIN_OR + AV_U_EIN_OR_BLK
        ENDIF

    !------------------------------------------------------------------
    !   Calculate the block average value for the observables
    !------------------------------------------------------------------
    !   Average potential energy
        AVPE   = SUMPE/REAL(IBLOCK,DP)
        AVPE2  = SUMPE2/REAL(IBLOCK,DP)
    !   Average virial pressure
        AVPRES = SUMPRES/REAL(IBLOCK,DP)
    !   Average density
        AVRHO  = SUMRHO/REAL(IBLOCK,DP)
    !   Average volume
        AVVLM  = SUMVLM/REAL(IBLOCK,DP) 
        AVVLM2 = SUMVLM2/REAL(IBLOCK,DP) 
        AVH    = SUMH/REAL(IBLOCK,DP) 
        AVH2   = SUMH2/REAL(IBLOCK,DP) 
        AVHV   = SUMHV/REAL(IBLOCK,DP) 
    !   Average box lengths
        IF (NPTT) AVBOX = SUMBOX/REAL(IBLOCK,DP)
        IF (NPZTT .OR. PINT) THEN
        !   Average collective density field
            AVTRQ  = SUMTRQ/REAL(IBLOCK,DP)
        !   Average chemical potential difference
            IF(PINT) AVDELMU = SUMDELMU/REAL(IBLOCK,DP)
        ENDIF

    !   Average size of largest crystalline cluster
        IF (NUCSEEDT) THEN
            AVLRGSTXCLSTR = SUMLRGSTXCLSTR/REAL(IBLOCK,DP)
            AVNUCCOUNT = SUMNUCCOUNT/REAL(IBLOCK,DP)
        ENDIF

        IF (FLFET) THEN
            AV_EXP_U_EIN = SUM_EXP_U_EIN/REAL(IBLOCK,DP)
            AV_U_EIN_TR  = SUM_U_EIN_TR/REAL(IBLOCK,DP)
            AV_U_EIN_OR  = SUM_U_EIN_OR/REAL(IBLOCK,DP)
        ENDIF

        IF (SCHSMIT) THEN
            AV_EXP_U_SS  = SUM_EXP_U_SS/REAL(IBLOCK,DP)
            AV_U_SS      = SUM_U_SS/REAL(IBLOCK,DP)
            AV_U_EIN_OR  = SUM_U_EIN_OR/REAL(IBLOCK,DP)
        ENDIF

    !-----------------------------------------------------------------------------
    !   Calculate the standard deviation for the average value of the observables
    !-----------------------------------------------------------------------------
        STDPESUM   = STDPESUM   + (AVPEBLK   - AVPE  )**2
        STDPE2SUM  = STDPE2SUM  + (AVPE2BLK  - AVPE2 )**2
        STDPRESSUM = STDPRESSUM + (AVPRESBLK - AVPRES)**2
        STDRHOSUM  = STDRHOSUM  + (AVRHOBLK  - AVRHO )**2
        STDVLMSUM  = STDVLMSUM  + (AVVLMBLK  - AVVLM )**2
        STDVLM2SUM = STDVLM2SUM + (AVVLM2BLK - AVVLM2)**2
        STDHSUM    = STDHSUM    + (AVHBLK    - AVH   )**2
        STDH2SUM   = STDH2SUM   + (AVH2BLK   - AVH2  )**2
        STDHVSUM   = STDHVSUM   + (AVHVBLK   - AVHV  )**2
        IF (NPTT) STDBOXSUM = STDBOXSUM + (AVBOXBLK - AVBOX)**2
        IF (NPZTT .OR. PINT) THEN
            STDTRQSUM = STDTRQSUM + (AVTRQBLK - AVTRQ)**2
            IF(PINT) STDDELMUSUM = STDDELMUSUM + (DELMUBLK - AVDELMU)**2
        ENDIF
        
        IF (NUCSEEDT) THEN
            STDLRGSTXCLSTRSUM = STDLRGSTXCLSTRSUM + (AVLRGSTXCLSTRBLK - AVLRGSTXCLSTR)**2
            STDNUCCOUNTSUM    = STDNUCCOUNTSUM + (AVNUCCOUNTBLK - AVNUCCOUNT)**2
        ENDIF

        IF (FLFET) THEN
            STD_EXP_U_EIN_SUM = STD_EXP_U_EIN_SUM + (AV_EXP_U_EIN_BLK - AV_EXP_U_EIN)**2
            STD_U_EIN_TR_SUM  = STD_U_EIN_TR_SUM + (AV_U_EIN_TR_BLK - AV_U_EIN_TR)**2
            STD_U_EIN_OR_SUM  = STD_U_EIN_OR_SUM + (AV_U_EIN_OR_BLK - AV_U_EIN_OR)**2
        ENDIF

        IF (SCHSMIT) THEN
            STD_EXP_U_SS_SUM  = STD_EXP_U_SS_SUM + ( AV_EXP_U_SS_BLK - AV_EXP_U_SS )**2
            STD_U_SS_SUM      = STD_U_SS_SUM + ( AV_U_SS_BLK - AV_U_SS )**2
            STD_U_EIN_OR_SUM  = STD_U_EIN_OR_SUM + (AV_U_EIN_OR_BLK - AV_U_EIN_OR)**2
        ENDIF

        IF(IBLOCK > 0) THEN
            INVBLK  = 1.0_dp / SQRT(REAL(IBLOCK,DP))
            INVBLK2 = 1.0_dp / REAL((IBLOCK - 1),DP)
            STDPE   = INVBLK*SQRT(INVBLK2*STDPESUM)
            STDPE2  = INVBLK*SQRT(INVBLK2*STDPE2SUM)
            STDPRES = INVBLK*SQRT(INVBLK2*STDPRESSUM)
            STDRHO  = INVBLK*SQRT(INVBLK2*STDRHOSUM)
            STDVLM  = INVBLK*SQRT(INVBLK2*STDVLMSUM)
            STDVLM2 = INVBLK*SQRT(INVBLK2*STDVLM2SUM)
            STDH    = INVBLK*SQRT(INVBLK2*STDHSUM)
            STDH2   = INVBLK*SQRT(INVBLK2*STDH2SUM)
            STDHV   = INVBLK*SQRT(INVBLK2*STDHVSUM)
            IF (NPTT) STDBOX = INVBLK*SQRT(INVBLK2*STDBOXSUM)
            IF (NPZTT .OR. PINT) THEN
                STDTRQ = INVBLK*SQRT(INVBLK2*STDTRQSUM)
                IF(PINT) STDDELMU = INVBLK*SQRT(INVBLK2*STDDELMUSUM)
            ENDIF
            
            IF (NUCSEEDT) THEN
                STDLRGSTXCLSTR = INVBLK*SQRT(INVBLK2*STDLRGSTXCLSTRSUM)
                STDNUCCOUNT = INVBLK*SQRT(INVBLK2*STDNUCCOUNTSUM)
            ENDIF

            IF (FLFET) THEN
                STD_EXP_U_EIN = INVBLK*SQRT(INVBLK2*STD_EXP_U_EIN_SUM)
                STD_U_EIN_TR  = INVBLK*SQRT(INVBLK2*STD_U_EIN_TR_SUM)
                STD_U_EIN_OR  = INVBLK*SQRT(INVBLK2*STD_U_EIN_OR_SUM)
            ENDIF

            IF (SCHSMIT) THEN
                STD_EXP_U_SS = INVBLK*SQRT(INVBLK2*STD_EXP_U_SS_SUM)
                STD_U_SS     = INVBLK*SQRT(INVBLK2*STD_U_SS_SUM)
                STD_U_EIN_OR = INVBLK*SQRT(INVBLK2*STD_U_EIN_OR_SUM)
            ENDIF

        ENDIF
        ! CV, CP, KT, AL, STDCV, STDCP, STDKT, STDAL
        NP = REAL(NPART,DP)
        OPEN (UNIT = 3, FILE ='ave_data.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
        WRITE (3,*) "BLOCK: ", IBLOCK, ISTEP
    !   Instantaneous energy
        WRITE (3,*) "ENERGY: ", AVPE/NP, STDPE/NP
        IF(NVTT) THEN
            B2 = BETAKB**2
            DE = AVPE2-AVPE**2
            CV = B2*DE/NP
            WRITE (3,*) "HEAT CAPACITY: ", CV!, B2*SQRT( STDPE2**2 + (2.0_dp*AVPE*STDPE)**2 )/NP, AVPE2, STDPE2
        ENDIF
    !   Average and instantaneous pressure (only valid for continuous models)       
        WRITE (3,*) "PRESSURE: ", AVPRES, STDPRES
        IF(NPTT) THEN
        !   Average and instantaneous density      
            WRITE (3,*) "DENSITY: ", AVRHO, STDRHO
        !   Average and instantaneous volume  
            WRITE (3,*) "VOLUME: ", AVVLM, STDVLM
            B2  = BETAKB**2
            NV  = NP*AVVLM
            DH  = AVH2-AVH**2
            DV  = AVVLM2-AVVLM**2
            DVH = AVHV-AVVLM*AVH

            CP     = B2*DH/NP
            STDCP  = (B2/NP)*SQRT( STDH2**2 + (2.D0*AVH)**2 * STDH**2 )
            KT     = BETAKB*DV/NV
            STDKT  = (BETAKB/NP) * SQRT( ((AVVLM2+AVVLM**2)/AVVLM**2)**2*STDVLM**2 + (STDVLM2**2)/AVVLM**2 )
            AL     = B2*DVH/NV
            ! STDAL  = 

            WRITE (3,*) "HEAT CAPACITY: ",     CP, STDCP
            WRITE (3,*) "COMPRESSIBILITY: ",   KT, STDKT
            WRITE (3,*) "THERMAL EXPANSION: ", AL, STDAL
        !   Average and instantaneous volume  
            IF(NPZTT .OR. PINT) THEN
                WRITE (3,*) "COLL_DENS: ", AVTRQ, STDTRQ
                IF(PINT) WRITE (3,*) "DEL_MU: ", AVDELMU, STDDELMU
            ENDIF
        !   Average box lengths
            OPEN (UNIT = 5, FILE = 'avebox.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            WRITE (5,*) "BLOCK: ", IBLOCK, ISTEP
            WRITE (5,*) AVBOX
            WRITE (5,*) STDBOX
            CLOSE (UNIT = 5, STATUS = 'KEEP')
        ENDIF

        IF (NUCSEEDT) THEN
            WRITE (3,*) "LRGST_CLSTR: ", AVLRGSTXCLSTR, STDLRGSTXCLSTR
            OPEN (UNIT = 4, FILE ='umbrella.dat', STATUS = 'UNKNOWN', ACCESS = 'APPEND')
            WRITE (4,*) "STEP: ", ISTEP
            DO J1 = 0, NPART
                IF(TOTNUCCOUNT(J1) > 0) WRITE (4,*) J1, TOTNUCCOUNT(J1), AVNUCCOUNT(J1), STDNUCCOUNT(J1)
            ENDDO
            CLOSE (UNIT = 4, STATUS = 'KEEP')
        ENDIF
        
        IF (FLFET) THEN
            
            IF(FLID==1) THEN
                WRITE (3,*) "dA1: ", ( BETAKB*U_LATT-LOG(AV_EXP_U_EIN) )/NP, (STD_EXP_U_EIN/AV_EXP_U_EIN)/NP, &
                AV_EXP_U_EIN, STD_EXP_U_EIN

            ELSEIF(FLID==2) THEN
                BENP = BETAKB / NP
                WRITE (3,*) "EIN_TR: ", BENP*AV_U_EIN_TR, BENP*STD_U_EIN_TR, AV_U_EIN_TR, STD_U_EIN_TR
                IF(RIGIDT) THEN
                    WRITE (3,*) "EIN_OR: ", BENP*AV_U_EIN_OR, BENP*STD_U_EIN_OR, AV_U_EIN_OR, STD_U_EIN_OR
                ENDIF
            ENDIF
        ENDIF

        IF (SCHSMIT) THEN
            NP = REAL(NPART,DP)
            
            IF(SSID==1) THEN
                WRITE (3,*) "dA1: ", ( BETAKB*U_LATT-LOG(AV_EXP_U_SS) )/NP, (STD_EXP_U_SS/AV_EXP_U_SS)/NP, &
                AV_EXP_U_SS, STD_EXP_U_SS

            ELSEIF(SSID==2) THEN
                BENP = BETAKB / NP
                WRITE (3,*) "WELL: ", BENP*AV_U_SS, BENP*STD_U_SS, AV_U_SS, STD_U_SS
                IF(RIGIDT) THEN
                    WRITE (3,*) "EIN_OR: ", BENP*AV_U_EIN_OR, BENP*STD_U_EIN_OR, AV_U_EIN_OR, STD_U_EIN_OR
                ENDIF
            ENDIF
        ENDIF

        CLOSE (UNIT = 3, STATUS = 'KEEP')

    !   Reset parameters to calculate the averages for the next block 
        SUMPEBLK          = 0.0_dp
        SUMPE2BLK         = 0.0_dp

        SUMPRESBLK        = 0.0_dp
        SUMRHOBLK         = 0.0_dp
        SUMVLMBLK         = 0.0_dp
        SUMVLM2BLK        = 0.0_dp
        SUMHBLK           = 0.0_dp
        SUMH2BLK          = 0.0_dp
        SUMHVBLK          = 0.0_dp
        SUMBOXBLK         = 0.0_dp
        
        SUMTRQBLK         = 0.0_dp
        SUMLRGSTXCLSTRBLK = 0.0_dp

        IF (NUCSEEDT) NUCCOUNT = 0

        SUM_EXP_U_EIN_BLK = 0.0_dp
        SUM_U_EIN_TR_BLK  = 0.0_dp
        SUM_U_EIN_OR_BLK  = 0.0_dp

        SUM_EXP_U_SS_BLK = 0.0_dp
        SUM_U_SS_BLK     = 0.0_dp

    ENDIF

!   ==============================================================================================
!   Following an initial equilibration period reset all averages and start computation again.
!   ==============================================================================================
    IF ( ISTEP == NEQ ) THEN
    !----------------------------------------------------------------------------------   
    !   Observables for standard NVT or NPT Monte Carlo simulations
    !----------------------------------------------------------------------------------
        ICNT    = 0;      IBLOCK     = 0
        SUMPE   = 0.0_dp; STDPESUM   = 0.0_dp
        SUMPE2  = 0.0_dp; STDPE2SUM  = 0.0_dp
        SUMPRES = 0.0_dp; STDPRESSUM = 0.0_dp
        SUMVLM  = 0.0_dp; STDVLMSUM  = 0.0_dp
        SUMVLM2 = 0.0_dp; STDVLM2SUM = 0.0_dp
        SUMH    = 0.0_dp; STDHSUM    = 0.0_dp
        SUMH2   = 0.0_dp; STDH2SUM   = 0.0_dp
        SUMHV   = 0.0_dp; STDHVSUM   = 0.0_dp
        SUMRHO  = 0.0_dp; STDRHOSUM  = 0.0_dp
        SUMBOX  = 0.0_dp; STDBOXSUM  = 0.0_dp

    !----------------------------------------------------------------------------------
    !   Chemical potential and order parameter for interface-pinning simulations
    !----------------------------------------------------------------------------------
        SUMTRQ         = 0.0_dp; STDTRQSUM         = 0.0_dp
        SUMDELMU       = 0.0_dp; STDDELMUSUM       = 0.0_dp

    !----------------------------------------------------------------------------------
    !   Count of largest crystalline cluster for umbrella sampling and nucleus-size
    !   pinning simulations.
    !----------------------------------------------------------------------------------
        SUMLRGSTXCLSTR = 0.0_dp; STDLRGSTXCLSTRSUM = 0.0_dp
        IF (NUCSEEDT) THEN
            NUCCOUNT    = 0 
            SUMNUCCOUNT = 0; STDNUCCOUNTSUM    = 0.0_dp
            TOTNUCCOUNT = 0
        ENDIF

    !----------------------------------------------------------------------------------
    !   Free energies from Einstein-crystal integration simulations
    !----------------------------------------------------------------------------------
        SUM_EXP_U_EIN  = 0.0_dp; STD_EXP_U_EIN_SUM  = 0.0_dp
        SUM_U_EIN_TR   = 0.0_dp; STD_U_EIN_TR_SUM   = 0.0_dp
        SUM_U_EIN_OR   = 0.0_dp; STD_U_EIN_OR_SUM   = 0.0_dp

    !----------------------------------------------------------------------------------
    !   Free energies from Schilling-Schmid integration simulations
    !----------------------------------------------------------------------------------
        SUM_U_SS     = 0.0_dp; STD_U_SS_SUM     = 0.0_dp 
        SUM_EXP_U_SS = 0.0_dp; STD_EXP_U_SS_SUM = 0.0_dp 
        
    ENDIF

END SUBROUTINE ACCUMULATORS
