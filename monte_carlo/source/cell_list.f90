! link_list_module.f90
! Link clist handling routines for MC or MD simulation
MODULE CELL_LIST

  !------------------------------------------------------------------------------------------------!
  ! This software was written in 2016/17                                                           !
  ! by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        !
  ! and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             !
  ! to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     !
  ! published by Oxford University Press ("the publishers").                                       !
  !                                                                                                !
  ! LICENCE                                                                                        !
  ! Creative Commons CC0 Public Domain Dedication.                                                 !
  ! To the extent possible under law, the authors have dedicated all copyright and related         !
  ! and neighboring rights to this software to the PUBLIC domain worldwide.                        !
  ! This software is distributed without any warranty.                                             !
  ! You should have received a copy of the CC0 Public Domain Dedication along with this software.  !
  ! If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               !
  !                                                                                                !
  ! DISCLAIMER                                                                                     !
  ! The authors and publishers make no warranties about the software, and disclaim liability       !
  ! for all uses of the software, to the fullest extent permitted by applicable law.               !
  ! The authors and publishers do not recommend use of this software for any purpose.              !
  ! It is made freely available, solely to clarify points made in the text. When using or citing   !
  ! the software, you should not imply endorsement by the authors or publishers.                   !
  !------------------------------------------------------------------------------------------------!

    USE, INTRINSIC :: iso_fortran_env, ONLY : output_unit, error_unit
    USE COMMONS, ONLY: DP

    IMPLICIT NONE
    PRIVATE

    ! Public routines
    PUBLIC :: initialize_list, finalize_list, make_list
    PUBLIC :: move_in_list, create_in_list, destroy_in_list, c_index, neighbours

    ! Public (protected) data
    INTEGER,                                PROTECTED, SAVE, PUBLIC :: scd(3)     ! dimensions of head array
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, PROTECTED, SAVE, PUBLIC :: head   ! head(0:scd-1,0:scd-1,0:scd-1)
    INTEGER, DIMENSION(:),     ALLOCATABLE, PROTECTED, SAVE, PUBLIC :: clist   ! clist(n)
    INTEGER, DIMENSION(:,:),   ALLOCATABLE, PROTECTED, SAVE, PUBLIC :: c      ! c(3,n) 3D cell index of each atom

CONTAINS

    SUBROUTINE INITIALIZE_LIST ( N, R_CUT_BOX ) ! Routine to allocate clist arrays
        USE COMMONS, ONLY: ISTEP
        IMPLICIT NONE
        INTEGER, INTENT(in)        :: N         ! Number of particles
        REAL(KIND=DP),  INTENT(in) :: R_CUT_BOX(3) ! rcut/box, assume never changes

        SCD = FLOOR ( 1.0_dp / R_CUT_BOX ) ! Number of cells in each dimension

        IF(ISTEP == 0) PRINT *, SCD, N, R_CUT_BOX

        ALLOCATE ( CLIST(N), C(3,N) )
        ALLOCATE ( HEAD(0:SCD(1)-1, 0:SCD(2)-1, 0:SCD(3)-1) )

    END SUBROUTINE INITIALIZE_LIST

    SUBROUTINE FINALIZE_LIST ! Routine to deallocate clist arrays
        IMPLICIT NONE

        DEALLOCATE ( clist, c )
        DEALLOCATE ( head )

    END SUBROUTINE FINALIZE_LIST

    SUBROUTINE MAKE_LIST ( N, R ) ! Routine to make clist
        IMPLICIT NONE
        INTEGER,                 INTENT(in)         :: N ! Number of atoms
        REAL(KIND=DP),  DIMENSION(3,n), INTENT(in)  :: R ! Atom coordinates

        INTEGER :: I
        
        HEAD(:,:,:) = 0

        DO I = 1, N ! Loop over all atoms
            C(:,I) = C_INDEX ( R(:,I) )       ! Index function allocating atom i to cell
            CALL CREATE_IN_LIST ( I, C(:,I) ) ! This does the work of adding atom i to clist
        END DO ! End loop over all atoms

    END SUBROUTINE MAKE_LIST

    FUNCTION C_INDEX ( RI ) RESULT ( CI )
        USE COMMONS, ONLY: FLFET, TDSRFT, SPHERECNFT
        IMPLICIT NONE
        INTEGER, DIMENSION(3)                       :: CI ! Returns 3D cell index in range 0..scd-1, calculated from
        REAL(KIND=DP),  DIMENSION(3), INTENT(in)    :: RI ! position in box = 1 units
        REAL(KIND=DP),  DIMENSION(3)                :: RO

        RO = RI
    !   Check that RI is within the bounds of the simulation box.
        IF ( ANY ( ABS(RO) > 0.501 ) .AND. (.NOT. TDSRFT) .AND. (.NOT. SPHERECNFT)) THEN ! SHOULD NEVER HAPPEN
            IF(FLFET) THEN
        !   This condition should ONLY EVER fail for a simulation with a fixed center
        !   of mass as, in this particular case, the particles are never returned
        !   to the simulation cell when they leave.
                RO = RO - ANINT(RO)
            ELSE
                WRITE ( UNIT=ERROR_UNIT, FMT='(A,3F25.10)') 'ATOM NOT IN MAIN BOX', RO
                STOP 'Error in c_index'
            ENDIF
        END IF
        
        CI(:) = FLOOR ( ( RO(:) + 0.5_DP ) * REAL(SCD,DP) ) ! The index formula
        ! Guard against small chance of roundoff error
        WHERE ( CI(:) < 0    ) CI(:) = 0
        WHERE ( CI(:) > SCD(:)-1 ) CI(:) = SCD(:)-1

    END FUNCTION c_index

    SUBROUTINE create_in_list ( i, ci ) ! Routine to create atom i in cell ci
        IMPLICIT NONE
        INTEGER,               INTENT(in) :: i  ! Index of atom
        INTEGER, DIMENSION(3), INTENT(in) :: ci ! 3D index of cell in which i lies
        
        clist(i)                 = head(ci(1),ci(2),ci(3)) ! Transfer old head to clist
        head(ci(1),ci(2),ci(3)) = i                       ! Atom i becomes new head for this clist
        c(:,i)                  = ci(:)                   ! Store 3D index in array

    END SUBROUTINE create_in_list

    SUBROUTINE destroy_in_list ( i, ci ) ! Routine to destroy atom i in cell ci
        IMPLICIT NONE
        INTEGER,               INTENT(in) :: i  ! Index of atom
        INTEGER, DIMENSION(3), INTENT(in) :: ci ! 3D index of cell in which i lies

        INTEGER :: this, next
        
        this = head(ci(1),ci(2),ci(3)) ! Locate head of clist corresponding to cell

        IF ( this == i ) THEN ! Atom i is the head atom in that cell

        head(ci(1),ci(2),ci(3)) = clist(i) ! Simply point head at next atom, we're done

        ELSE ! Atom i lies further down the clist

        DO ! Loop traversing link-clist
            next = clist(this) ! Look ahead to the next entry

            IF ( next == i ) THEN ! Found our atom, just link over it
                clist(this) = clist(i)
                EXIT ! Leave the loop

            ELSE IF ( next == 0 ) THEN ! This should never happen
                WRITE ( unit=error_unit, fmt='(a,4i15)') 'Could not find particle in its cell', i, ci
                STOP 'Error in destroy_in_list'

            ELSE ! Move on to the next
                this = next ! Keep this index for next iteration
            END IF

        END DO ! End loop traversing link-clist

        END IF

    END SUBROUTINE destroy_in_list

    SUBROUTINE move_in_list ( i, ci ) ! Routine to move atom i from current cell to ci
        IMPLICIT NONE
        INTEGER,               INTENT(in) :: i
        INTEGER, DIMENSION(3), INTENT(in) :: ci
        
        IF ( ALL ( ci(:) == c(:,i) ) ) RETURN ! No need to do anything

        CALL destroy_in_list ( i, c(:,i) ) ! Remove atom i from old cell
        CALL create_in_list  ( i, ci(:)  ) ! Add atom i to new cell

    END SUBROUTINE move_in_list

    FUNCTION neighbours ( n, i, ci, half ) RESULT ( j_list )
        IMPLICIT NONE
        INTEGER,               INTENT(in) :: n      ! Number of atoms
        INTEGER,               INTENT(in) :: i      ! Atom whose neighbours are required
        INTEGER, DIMENSION(3), INTENT(in) :: ci     ! Cell of atom of interest
        LOGICAL,               INTENT(in) :: half   ! Determining the range of neighbours searched
        INTEGER, DIMENSION(n)             :: j_list ! Resulting clist of indices

        ! This routine uses the link-clist cell structure to fill out the array j_list
        ! with possible neighbours of atom i, padding with zeroes
        ! If half==.false., cell ci and all 26 surrounding cells are searched.
        ! If half==.true., cell ci, and just 13 of the neighbour cells, are searched
        ! and moreover, in ci, we only look down-clist making use of clist(i)
        ! There is a subtlety: using clist(i) assumes that our interest is in the cells that
        ! are neighbours of c(:,i), i.e. that ci(:)==c(:,i), and we check for this explicitly.
        ! In other words, we assume that atom i has not moved since clist was constructed.
        ! If half==.false., particle i might be in a very different position, and ci might be
        ! very different to c(:,i) but in that case we make no use of clist(i), in normal use

        ! We have a cubic cell lattice
        ! Set up vectors to each cell in the 3x3x3 neighbourhood of a given cell
        ! To work properly, these are listed with inversion symmetry about (0,0,0)
        INTEGER,                      PARAMETER :: nk = 13 
        INTEGER, DIMENSION(3,-nk:nk), PARAMETER :: d = RESHAPE( [ &
            &   -1,-1,-1,    0,-1,-1,    1,-1,-1, &
            &   -1, 1,-1,    0, 1,-1,    1, 1,-1, &
            &   -1, 0,-1,    1, 0,-1,    0, 0,-1, &
            &    0,-1, 0,    1,-1, 0,   -1,-1, 0, &
            &   -1, 0, 0,    0, 0, 0,    1, 0, 0, &
            &    1, 1, 0,   -1, 1, 0,    0, 1, 0, &
            &    0, 0, 1,   -1, 0, 1,    1, 0, 1, &
            &   -1,-1, 1,    0,-1, 1,    1,-1, 1, &
            &   -1, 1, 1,    0, 1, 1,    1, 1, 1    ], [ 3, 2*nk+1 ] )

        INTEGER               :: k1, k2, k, j, nj, err
        INTEGER, DIMENSION(3) :: cj

        IF ( half ) THEN ! Check half neighbour cells and j downlist from i in current cell
            k1 = 0
            k2 = nk
            IF ( ANY ( ci(:) /= c(:,i) ) ) THEN ! should never happen
                WRITE ( unit=error_unit, fmt='(a,7i15)' ) 'Cell mismatch ', i, ci(:), c(:,i)
                STOP 'Error in get_neighbours'
            END IF
        ELSE ! Check every atom other than i in all cells
            k1 = -nk
            k2 =  nk
        END IF

        j_list = 0 ! Initialize with zero values everywhere
        nj     = 0 ! Next position in clist to be filled

        DO k = k1, k2 ! Begin loop over neighbouring cells

            cj(:) = ci(:) + d(:,k)       ! Neighbour cell index
            cj(:) = MODULO ( cj(:), scd ) ! Periodic boundary correction

            IF ( k == 0 .AND. half ) THEN
                j = clist(i) ! Check down-clist from i in i-cell
            ELSE
                j = head(cj(1),cj(2),cj(3)) ! Check entire j-cell
            END IF

            DO ! Begin loop over j atoms in clist

                IF ( j == 0 ) EXIT  ! Exhausted clist

                IF ( j /= i ) THEN ! Skip self

                    nj = nj + 1 ! Increment count of j atoms

                    IF ( nj >= n ) THEN ! Check more than n-1 neighbours (should never happen)
                        WRITE ( unit=error_unit, fmt='(a,2i15)' ) 'Neighbour error for j_list', nj, n, j
                        DO ERR = 1, n
                            WRITE ( unit=error_unit, fmt='(a,2i15)' ) ERR, j_list(err)
                        ENDDO
                        STOP 'Impossible error in get_neighbours'
                    END IF ! End check more than n-1 neighbours

                    j_list(nj) = j       ! Store new j atom
                END IF

                j = clist(j) ! Next atom in j cell

            END DO ! End loop over j atoms in clist

        END DO ! End loop over neighbouring cells 

    END FUNCTION neighbours

END MODULE CELL_LIST