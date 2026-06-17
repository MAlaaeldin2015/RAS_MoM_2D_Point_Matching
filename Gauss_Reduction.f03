! =============================================================================
! Module: Gauss_Reduction
! Purpose: Dense complex matrix utilities â inversion, elimination, determinant.
!          Used as a fallback solver and for small auxiliary systems within the
!          RAS/MoM solver. For the main MoM matrix, the LAPACK-based solvers in
!          Matrices_Storage_Handling are preferred.
!
! Author : Mohamed A. Moharram Hassan
! License: Apache 2.0 â see LICENSE
! =============================================================================
module Gauss_Reduction
    implicit none

    ! Pivot threshold: diagonal elements with |value| < THRESH are treated as zero
    real(8), parameter :: MAX_REAL      = 1.d10
    integer, parameter :: MAX_ITERATIONS = 200
    real(8), parameter :: THRESH        = 1.d-8

    contains

    ! -------------------------------------------------------------------------
    ! Function: Mat_Inv
    ! Purpose : Computes the inverse of a complex square matrix A using
    !           Gauss-Jordan elimination with partial pivoting.
    !           Augments A with the identity matrix [A | I], then row-reduces
    !           to [I | A^{-1}].
    !
    ! Input  : A  â NÃN complex matrix (non-singular)
    ! Returns: Y  â NÃN complex matrix (inverse of A)
    !
    ! Stops with an error if A is detected as singular (|pivot| < THRESH after
    ! all row-swap attempts are exhausted).
    ! -------------------------------------------------------------------------
    function Mat_Inv(A) result(Y)
    complex*16, dimension(:,:), intent (in) :: A
    complex*16, dimension(:,:), allocatable :: D,Y
    integer :: m,k,dim1,dim2
    logical :: FLAG

    dim1 = size(A,1)
    dim2 = 2 * dim1
    allocate (D(dim1,dim2))
    allocate(Y(dim1,dim1))
    D = 0.d0
    D(1:dim1, 1:dim1) = A(1:dim1, 1:dim1)
    do m = 1, dim1
        D(m, dim1+m) = 1.d0
    enddo

    flag = .true.
    do m = 1, dim1
         if (abs(D(m,m)) < THRESH) then
                FLAG = .false.
                do k = (m+1), dim1
                    if (abs(D(k,m)) > THRESH) then
                        D(m,:) = D(m,:) + D(k,:)
                        FLAG = .true.
                        exit
                    endif
                enddo
            endif
            if (FLAG .eqv. .false.) then
                write(*,*) 'Matrix A is Singular'
                STOP
            endif

        D(m, 1:dim2) = D(m, 1:dim2) / D(m,m)
        do k = 1, dim1
            if (k /= m) then
                D(k, 1:dim2) = D(k, 1:dim2) - D(k,m) * D(m, 1:dim2)
            endif
        enddo
    enddo

    Y(1:dim1, 1:dim1) = D(1:dim1, dim1+1:dim2)
    D(1:dim1, 1:dim1) = matmul(Y, A)

    deallocate (D)
    end function Mat_Inv

    ! -------------------------------------------------------------------------
    ! Function: Gauss_Elimination
    ! Purpose : Solves the complex linear system  AÂ·x = b  using Gaussian
    !           elimination with partial pivoting, followed by back-substitution.
    !
    ! Inputs : A  â NÃN complex coefficient matrix
    !          b  â N-vector right-hand side
    !          N  â system dimension
    ! Returns: x  â N-vector solution
    !
    ! Returns x = 0 and prints a warning if A is detected as singular.
    ! Note: the augmented matrix D is stack-allocated with a fixed upper bound;
    !       N must not exceed the hard-coded array bounds.
    ! -------------------------------------------------------------------------
    function Gauss_Elimination(A, b, N) result(x)
        integer, intent(in) :: N
        complex*16, dimension(:,:), intent(in) :: A
        complex*16, dimension(:),   intent(in) :: b
        complex*16, dimension(size(b)) :: x
        logical :: FLAG
        complex*16, dimension(2*N, 3*N+1) :: D
        integer :: mm, m, k, dim1, dim2

        dim1 = N
        dim2 = 2*N + 1
        D = 0.0
        D(1:dim1, 1:dim1) = A(1:dim1, 1:dim1)
        D(1:dim1, dim2)   = b(1:dim1)
        do m = 1, dim1
            D(m, dim1+m) = 1.d0
        enddo
        x = 0.d0

        FLAG = .true.
        do mm = 1, dim1
            if (abs(D(mm,mm)) < THRESH) then
                FLAG = .false.
                do k = (mm+1), dim1
                    if (abs(D(k,mm)) > THRESH) then
                        D(mm,:) = D(mm,:) + D(k,:)
                        FLAG = .true.
                        exit
                    endif
                enddo
            endif
            if (FLAG .eqv. .false.) then
                write(*,*) 'Matrix A is Singular'
                return
            endif
            D(mm, 1:dim2) = D(mm, 1:dim2) / D(mm,mm)
            do k = mm, dim1
                if (k /= mm) then
                    D(k, 1:dim2) = D(k, 1:dim2) - D(k,mm) * D(mm, 1:dim2)
                endif
            enddo
        enddo

        do m = dim1, 1, -1
           x(m) = D(m, dim2)
           do k = dim1, (m+1), -1
               x(m) = x(m) - D(m,k) * x(k)
           enddo
        enddo

    end function Gauss_Elimination

    ! -------------------------------------------------------------------------
    ! Function: Det
    ! Purpose : Computes the determinant of a complex square matrix using
    !           upper-triangular reduction (Gaussian elimination without
    !           back-substitution). Row swaps are tracked via sign flips.
    !
    ! Input  : mat â NÃN complex matrix
    ! Returns: Det â complex scalar determinant
    !
    ! Prints a warning and returns 0 if the matrix is found to be singular.
    ! -------------------------------------------------------------------------
    complex*16 FUNCTION Det(mat)
        IMPLICIT NONE
        complex*16, DIMENSION(:,:), intent(in) :: mat
        complex*16, DIMENSION(size(mat,1), size(mat,1)) :: matrix
        INTEGER :: n
        complex*16 :: m, temp
        INTEGER :: i, j, k, l
        LOGICAL :: DetExists = .TRUE.

        matrix = mat
        n = size(matrix,1)
        l = 1
        DO k = 1, n-1
            IF (matrix(k,k) == 0) THEN
                DetExists = .FALSE.
                DO i = k+1, n
                    IF (matrix(i,k) /= 0) THEN
                        DO j = 1, n
                            temp = matrix(i,j)
                            matrix(i,j) = matrix(k,j)
                            matrix(k,j) = temp
                        END DO
                        DetExists = .TRUE.
                        l = -l
                        EXIT
                    ENDIF
                END DO
                IF (DetExists .EQV. .FALSE.) THEN
                    Det = 0.d0
                    write(*,*) 'Warning: Determinant Does not exist'
                    return
                END IF
            ENDIF
            DO j = k+1, n
                m = matrix(j,k) / matrix(k,k)
                DO i = k+1, n
                    matrix(j,i) = matrix(j,i) - m * matrix(k,i)
                END DO
            END DO
        END DO

        Det = 1.d0
        DO i = 1, n
            Det = Det * matrix(i,i)
        END DO

    END FUNCTION Det

end module Gauss_Reduction
