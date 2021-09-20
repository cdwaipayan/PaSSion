PROGRAM RADDIS

    IMPLICIT NONE
    
    INTEGER, PARAMETER :: GINT = 150, NMAX = 256, NFRAMES = 1000
    DOUBLE PRECISION, PARAMETER :: PI = 4.D0*DATAN(1.D0)!, PHI = 0.2D0
    INTEGER :: I, J, A, B, IGNT, GR, NPART
    DOUBLE PRECISION :: DIS(3), DISSQ, GDIS, GRCUT, GRCUTSQ, BOXL, PHI, DUMMY1, DUMMY2
    DOUBLE PRECISION :: GRCNT(GINT), GRTEMP(GINT), R(NMAX,3), VLM, DGR, G, DA, Y

!   As the number of effective particles in the box is variable so is the
!   packing fraction. Therefore must input the edge length of the simulation to
!   allow the packing fraction to be calculated for each frame.    
    
    OPEN(UNIT = 14, FILE = 'tet.dat', STATUS = 'UNKNOWN')
    ! OPEN(UNIT = 15, FILE = '../box.dat', STATUS = 'UNKNOWN')
    BOXL = 14.167906088869D0!17.099759466767D0
    DO A = 1, NFRAMES
!   Read in number of tetrahedra identified for frame A
        READ(14,*) NPART
!   Read in the positions of the tetrahedra        
        DO B = 1, NPART
            READ(14,*) R(B,1), R(B,2), R(B,3)
        END DO
!   Calculate the packing fraction associated with frame A
        ! READ(15,*) BOXL, DUMMY1, DUMMY2
        GRCUT    = BOXL/2.D0
        DGR      = GRCUT/DBLE(GINT)
        GRCNT(:) = 0.D0
    
        GRCUTSQ  = GRCUT*GRCUT
        PHI      = DBLE(NPART)*PI/((BOXL**3)*6.D0)

        GRTEMP = 0.D0
    
        DO I = 1, (NPART-1)
            DO J = (I+1), NPART
                DIS(:) = R(I,:) - R(J,:)
                DIS(:) = DIS(:) - BOXL*ANINT(DIS(:)/BOXL)
                DISSQ  = DOT_PRODUCT(DIS(:),DIS(:))
                IF (DISSQ < GRCUTSQ) THEN
                    GDIS      = SQRT(DISSQ)
                    GR        = ANINT(GDIS/DGR)
                    GRTEMP(GR) = GRTEMP(GR) + 2.D0
                END IF
            END DO
        END DO
    
        DO I = 1, GINT
            DA        = ( (DBLE(I+1)**3) - (DBLE(I)**3) ) * (DGR**3)
            DA        = DA*PI*(4.D0/3.D0)*PHI
            GRCNT(I)  = GRCNT(I) + GRTEMP(I)/(2.D0*(DA*NPART))
        END DO
    
    END DO
    
    OPEN(UNIT = 15, FILE = 'rdf.dat', STATUS = 'UNKNOWN')
    
    DO I = 1, GINT
        G  = DGR*(DBLE(I)+0.5D0)
        Y  = GRCNT(I)!/DBLE(NFRAMES)
        WRITE(15,*) G, Y
    END DO
    
END