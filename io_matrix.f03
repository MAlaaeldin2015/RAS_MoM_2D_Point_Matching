! =============================================================================
! Module: io_matrix
! Purpose: Output utilities for writing simulation results to text files.
!          All routines open the target file, write row-by-row, and close it.
!          Files are opened in REPLACE mode â existing files are overwritten.
!
! Author : Mohamed A. Moharram Hassan
! License: Apache 2.0 â see LICENSE
! =============================================================================
module io_matrix
    implicit none

    contains

    ! -------------------------------------------------------------------------
    ! Subroutine: write_matrix
    ! Purpose   : Writes a complex NÃM matrix to a formatted text file, one row
    !             per line. Each element is written as "(real, imag)" pairs in
    !             Fortran list-directed format.
    !
    ! Inputs: B         â NÃM complex matrix to write
    !         file_name â output filename (character string)
    ! -------------------------------------------------------------------------
    subroutine write_matrix(B, file_name)
        complex*16, dimension(:,:), intent(in) :: B
        complex*16, allocatable, dimension(:,:) :: A
        character, intent(in) :: file_name*(*)
        integer :: fd1, ii, NN, MM, kk

        NN = size(B,1)
        MM = size(B,2)
        allocate(A(NN,MM))
        A = B
        fd1 = 20
        open(fd1, FILE=file_name, FORM='FORMATTED', STATUS="REPLACE", RECL=SIZEOF(A))
        do ii = 1, NN
            write(fd1,*) (A(ii,kk), kk=1,MM)
        enddo
        close(fd1)
        deallocate(A)
    end subroutine write_matrix

    ! -------------------------------------------------------------------------
    ! Subroutine: write_rcs
    ! Purpose   : Writes a two-column real array to a text file. Intended for
    !             RCS (Radar Cross-Section) data: column 1 = angle or frequency,
    !             column 2 = RCS value.
    !
    ! Inputs: rcs_data  â NÃ2 real array
    !         file_name â output filename
    ! -------------------------------------------------------------------------
    subroutine write_rcs(rcs_data, file_name)
        real*8, dimension(:,:), intent(in) :: rcs_data
        character, intent(in) :: file_name*(*)
        integer :: fd1, ii, NN, kk

        fd1 = 60
        NN  = size(rcs_data,1)
        open(fd1, FILE=file_name, FORM='FORMATTED', STATUS="REPLACE", RECL=SIZEOF(rcs_data))
        do ii = 1, NN
            write(fd1,*) (rcs_data(ii,kk), kk=1,2)
        enddo
        close(fd1)
    end subroutine write_rcs

    ! -------------------------------------------------------------------------
    ! Subroutine: write_fields_axis
    ! Purpose   : Writes a two-column real array to a text file. Intended for
    !             field values along an axis: column 1 = position, column 2 = value.
    !
    ! Inputs: fields    â NÃ2 real array
    !         file_name â output filename
    ! -------------------------------------------------------------------------
    subroutine write_fields_axis(fields, file_name)
        real*8, dimension(:,:), intent(in) :: fields
        character, intent(in) :: file_name*(*)
        integer :: fd1, ii, NN, kk

        NN  = size(fields,1)
        fd1 = 100
        open(fd1, FILE=file_name, STATUS="REPLACE")
        do ii = 1, NN
            write(fd1,*) (fields(ii,kk), kk=1,2)
        enddo
        close(fd1)
    end subroutine write_fields_axis

    ! -------------------------------------------------------------------------
    ! Subroutine: write_matrix_MATLAB
    ! Purpose   : Splits a complex NÃM matrix into its real and imaginary parts
    !             and writes each to a separate file with suffixes "_real.dat"
    !             and "_imag.dat". The resulting pair can be loaded directly
    !             into MATLAB with load().
    !
    ! Inputs: V         â NÃM complex matrix
    !         file_name â base filename (suffix appended automatically)
    ! -------------------------------------------------------------------------
    subroutine write_matrix_MATLAB(V, file_name)
        complex*16, dimension(:,:), intent(in) :: V
        real(8), allocatable, dimension(:,:) :: Vreal, Vimag
        character, intent(in) :: file_name*(*)
        character*30 :: f1, f2
        integer :: fd1, fd2, ii, NN, MM, kk

        NN = size(V,1)
        MM = size(V,2)
        allocate(Vreal(NN,MM))
        allocate(Vimag(NN,MM))

        Vreal = dreal(V)
        Vimag = dimag(V)

        f1 = file_name // '_real.dat' // ''
        f2 = file_name // '_imag.dat' // ''

        fd1 = 100
        fd2 = 110
        open(fd1, FILE=f1, FORM='FORMATTED', STATUS="REPLACE", RECL=SIZEOF(Vreal))
        open(fd2, FILE=f2, FORM='FORMATTED', STATUS="REPLACE", RECL=SIZEOF(Vimag))
        do ii = 1, NN
            write(fd1,*) (Vreal(ii,kk), kk=1,MM)
            write(fd2,*) (Vimag(ii,kk), kk=1,MM)
        enddo
        close(fd1)
        close(fd2)
        deallocate(Vreal)
        deallocate(Vimag)
    end subroutine write_matrix_MATLAB

    ! -------------------------------------------------------------------------
    ! Subroutine: write_vector_MATLAB
    ! Purpose   : Splits a complex N-vector into real and imaginary parts and
    !             writes each to a separate file (suffixes "_real.dat" and
    !             "_imag.dat"), one value per line. Useful for exporting source
    !             current vectors to MATLAB.
    !
    ! Inputs: V         â N-element complex vector
    !         file_name â base filename (suffix appended automatically)
    ! -------------------------------------------------------------------------
    subroutine write_vector_MATLAB(V, file_name)
        complex*16, dimension(:), intent(in) :: V
        real(8), allocatable, dimension(:) :: Vreal, Vimag
        character, intent(in) :: file_name*(*)
        character*20 :: f1, f2
        integer :: fd1, fd2, ii, NN

        NN = size(V)
        allocate(Vreal(NN))
        allocate(Vimag(NN))

        Vreal = dreal(V)
        Vimag = dimag(V)

        f1 = file_name // '_real.dat' // ''
        f2 = file_name // '_imag.dat' // ''

        fd1 = 49
        fd2 = 50
        open(fd1, FILE=f1, FORM='FORMATTED', STATUS="REPLACE", RECL=SIZEOF(Vreal))
        open(fd2, FILE=f2, FORM='FORMATTED', STATUS="REPLACE", RECL=SIZEOF(Vimag))
        do ii = 1, NN
            write(fd1,*) Vreal(ii)
            write(fd2,*) Vimag(ii)
        enddo
        close(fd1)
        close(fd2)
        deallocate(Vreal)
        deallocate(Vimag)
    end subroutine write_vector_MATLAB

end module io_matrix
