! ============================================================================
! Module: Matrices_Storage_Handling
! Purpose: Linear algebra back-end for the RAS/MoM solver.
!          Provides:
!
!   1. Matrix_Storage_Linked_List (derived type)
!      Singly-linked list caching factored system matrices (LDL* or SVD) and
!      Taylor-series solution vectors Y(0:N_order) keyed by group index, so
!      already-factored sub-systems are not re-computed on later iterations.
!      Type-bound procedures:
!        put_entry  -- store a new entry keyed by group index
!        get_entry  -- retrieve a previously stored entry
!        destroy    -- recursively deallocate the entire list
!
!   2. LDL* (modified Cholesky) factorisation and solve
!      get_LU_Cholesky      -- factorise Hermitian matrix (packed storage)
!      forward_subst        -- forward substitution for LDL* (vector RHS)
!      forward_subst_matrix -- forward substitution (full-matrix form)
!      backward_subst       -- backward substitution via conjugate-transpose L
!      solve_LU_Cholesky    -- one-shot factorise + solve
!      matrix_mult          -- Hermitian matrix-vector product (packed)
!
!   3. LAPACK-based solvers (wrapping zhesv / zhetrs2 / zgesvd)
!      solve_Matrix_SVD                  -- SVD pseudoinverse, cond-truncated
!      solve_Matrix_preconditioned       -- diagonal-preconditioned zhesv
!      solve_Matrix_partitioning         -- block-partitioned zhesv
!      solve_Matrix_SVD_preconditioned   -- preconditioning + SVD combined
!
!   4. Utility
!      mat2vec -- maps (i,j) Hermitian index to packed 1-D vector index
!
! Storage convention: Hermitian matrices are stored in packed lower-triangular
!   form as 1D complex vectors Y of length N*(N+1)/2.
!   Element (i,j) with i>=j is at index Y( mat2vec(N,i,j) ).
!
! Author : Mohamed A. Moharram Hassan
! License: Apache 2.0 -- see LICENSE
! ============================================================================
module Matrices_Storage_Handling
    use omp_lib
    use sim_par
    implicit none

        type,public :: Matrix_Storage_Linked_List !! A linked list to store the computed and solved LU matrices for the multiple usage by the software
            class(Matrix_Storage_Linked_List), pointer :: next =>null()
            integer,pointer :: i_key =>null()
            complex*16,allocatable,dimension(:) :: Dm,Lm
            real(8),allocatable,dimension(:) :: Y_mat_norm
            integer,pointer :: Nm,N_Lm
            complex*16,allocatable,dimension(:,:) :: Y_mat,Yinv
            complex*16,allocatable,dimension(:,:) :: Y
            integer, allocatable,dimension(:) ::IWORK

        contains
            procedure :: put_entry => set_list_entry
            procedure :: get_entry => get_list_entry
            procedure :: destroy => destroy_list


        end type Matrix_Storage_Linked_List


    contains

    ! -----------------------------------------------------------------------
    ! Subroutine: destroy_list (type-bound as 'destroy')
    ! Purpose   : Recursively walks the linked list and deallocates every
    !             node's arrays (Dm, Lm, Y_mat, Yinv, IWORK, Y, i_key, Nm).
    !             Must be called before the root node goes out of scope.
    ! -----------------------------------------------------------------------
    recursive subroutine destroy_list(this)
        class(Matrix_Storage_Linked_List) :: this

        if(associated(this%next)) then
            call this%next%destroy()
            deallocate(this%next)
        endif
        if(associated(this%i_key)) then
            deallocate(this%i_key,this%Nm,this%N_Lm,this%Y)
            if(allocated(this%Dm)) then
                deallocate(this%Dm,this%Lm)
            endif
            if(allocated(this%Y_mat)) then
                deallocate(this%Y_mat,this%Yinv)
            endif

            if(allocated(this%IWORK)) then
                deallocate(this%IWORK,this%Y_mat_norm)
            endif
        endif

    end subroutine destroy_list

    RECURSIVE subroutine get_list_entry(this,key,Nm,Dm,Lm,Y,Y_mat,Yinv,N_ORDER,IWORK,Y_mat_norm)
        class(Matrix_Storage_Linked_List) :: this
        integer,intent(in) :: key !! the key is related to iteration number
        complex*16,allocatable,intent(inout) :: Dm(:),Lm(:)
        real(8),allocatable,intent(inout) :: Y_mat_norm(:)
        integer,allocatable,intent(inout) :: IWORK(:)
        complex*16,allocatable,intent(inout) :: Y_mat(:,:),Yinv(:,:)
        complex*16,allocatable,intent(inout) :: Y(:,:)
        integer :: Nm,N_ORDER



        if(associated(this%i_key)) then
!            write(*,*) 'Required Linked List key',key,'Current Key',this%i_key
            if(this%i_key == key) then !! the item was found
                Nm = this%Nm
                if(allocated(this%Dm)) then
                    allocate(Dm(Nm),Lm(this%N_Lm))
                    Dm = this%Dm
                    Lm = this%Lm
                endif
                if(allocated(this%Y_mat)) then
                    allocate(Y_mat(Nm,Nm),Yinv(Nm,Nm))
                    Y_mat = this%Y_mat
                    Yinv = this%Yinv
                endif
                if(allocated(this%IWORK)) then
                    allocate(IWORK(Nm),Y_mat_norm(Nm))
                    IWORK = this%IWORK
                    Y_mat_norm = this%Y_mat_norm
                endif
                allocate(Y(this%N_Lm,0:N_ORDER))
                Y = this%Y
            else !! not found
            !! look through the next pointer
!                write(*,*) associated(this%next)
                if(associated(this%next)) then
                    call this%next%get_entry(key,Nm,Dm,Lm,Y,Y_mat,Yinv,N_ORDER,IWORK,Y_mat_norm)
                else
                    write(*,*) 'Next Pointer Was Not Found!'
                    write(*,*) 'The item is not found, the consistency of the algorithm in questionalble (get_list_entry) routine'
                    stop
                endif
            endif

        else
            write(*,*) 'The item is not found, the consistency of the algorithm in questionalble (get_list_entry) routine'
            stop
        endif

    end subroutine get_list_entry

    RECURSIVE subroutine set_list_entry(this,key,Nm,Dm,Lm,Y,Y_mat,Yinv,N_ORDER,IWORK,Y_mat_norm)
        class(Matrix_Storage_Linked_List) :: this
        integer :: key !! the key is related to iteration number
        integer :: N_ORDER !! the order of considered Taylor Series Coefficients
        complex*16,allocatable,intent(in) :: Dm(:),Lm(:)
        real(8),allocatable,intent(inout) :: Y_mat_norm(:)
        integer,allocatable,intent(inout) :: IWORK(:)
        complex*16,allocatable,intent(in) :: Y_mat(:,:),Yinv(:,:)
        complex*16,allocatable,intent(in) :: Y(:,:)
        integer :: Nm
        integer :: N_Lm

        N_Lm = (Nm*(Nm+1))/2

        if(associated(this%i_key)) then
            if(this%i_key == key) then !! logical error as the keys should be different at each entry
                write(*,*) 'In consistency in the input keys in (set_list_entry) routine'
                stop
            else
                if(.not. associated(this%next)) then
                    allocate(this%next)
                endif
                call this%next%put_entry(key,Nm,Dm,Lm,Y,Y_mat,Yinv,N_ORDER,IWORK,Y_mat_norm)

!                write(*,*) key,associated(this%next)

            endif
        else
            allocate(this%i_key)
            this%i_key = key
            allocate(this%Nm,this%N_Lm)
            this%Nm = Nm
            this%N_Lm = N_Lm

            if(allocated(Dm)) then
                allocate(this%Dm(Nm),this%Lm(N_Lm))
                this%Dm = Dm
                this%Lm = Lm
            endif

            if(allocated(Y_mat)) then
                allocate(this%Y_mat(Nm,Nm),this%Yinv(Nm,Nm))
                this%Y_mat = Y_mat
                this%Yinv = Yinv
            endif
            if(allocated(IWORK)) then
                allocate(this%IWORK(Nm),this%Y_mat_norm(Nm))
                this%Y_mat_norm = Y_mat_norm
                this%IWORK = IWORK
            endif

            allocate(this%Y(N_Lm,0:N_ORDER))
            this%Y = Y
        endif


    end subroutine set_list_entry


    ! -----------------------------------------------------------------------
    ! Subroutine: forward_subst
    ! Purpose   : Forward substitution step of the LDL* solve.
    !             Computes y = (L D)^{-1} b using the packed lower-triangular
    !             factor Lm and diagonal Dm stored from get_LU_Cholesky.
    ! Inputs : Lm(N*(N+1)/2), Dm(N) -- LDL* factors; bm(N) -- RHS
    ! Output : y_out(N) -- intermediate solution vector
    ! -----------------------------------------------------------------------
    subroutine forward_subst(Lm,Dm,bm,Nm,y_out)
        complex*16,allocatable,intent(inout) :: y_out(:)
        complex*16,allocatable,intent(inout)::Lm(:),Dm(:),bm(:)
        integer :: Nm,ii,kk
        complex*16 :: summ

        !        Nm = size(bm,1)
        allocate(y_out(Nm))
        do ii =1,Nm
            summ = bm(ii)
            do kk = 1,ii-1
!                summ = summ - Dm(kk)*y_out(kk)*Lm((kk-1)*Nm+ii-(kk*(kk-1))/2)
                summ = summ - Dm(kk)*y_out(kk)*Lm((kk-1)*Nm+ii-rshift(kk*(kk-1),1))
            enddo
            y_out(ii) = summ/Dm(ii)
        enddo


    end subroutine forward_subst

    ! -----------------------------------------------------------------------
    ! Subroutine: forward_subst_matrix
    ! Purpose   : In-place forward substitution using a full (dense) lower-
    !             triangular matrix A. Overwrites bm with the solution.
    !             Used for the full-matrix solver path.
    ! -----------------------------------------------------------------------
    subroutine forward_subst_matrix(A,bm,Nm)
        complex*16,allocatable,dimension(:) :: y_out
        complex*16,allocatable,intent(inout)::A(:,:),bm(:)
        integer :: Nm,ii,kk
        complex*16 :: summ

        !        Nm = size(bm,1)
        allocate(y_out(Nm))
        do ii =1,Nm
            summ = bm(ii)
            do kk = 1,ii-1
!                summ = summ - Dm(kk)*y_out(kk)*Lm((kk-1)*Nm+ii-(kk*(kk-1))/2)
!                summ = summ - Dm(kk)*y_out(kk)*Lm((kk-1)*Nm+ii-rshift(kk*(kk-1),1))
                summ = summ - y_out(kk)*A(ii,kk)
            enddo
            y_out(ii) = summ/A(ii,ii)
        enddo

        bm = y_out
        deallocate(y_out)
    end subroutine forward_subst_matrix

    ! -----------------------------------------------------------------------
    ! Subroutine: backward_subst
    ! Purpose   : Backward substitution step of the LDL* solve.
    !             Computes x = L^{-H} y using the conjugate-transpose of the
    !             packed lower-triangular factor Lm.
    ! Inputs : Lm(N*(N+1)/2) -- LDL* factor; y_out(N) -- from forward_subst
    ! Output : x(N) -- final solution vector
    ! -----------------------------------------------------------------------
    subroutine backward_subst(Lm,y_out,Nm,x)
        complex*16, allocatable,intent(inout) :: Lm(:),y_out(:)
        complex*16,dimension(:),intent(out) :: x
        integer :: ii,k,Nm,ind_ii_dep
        complex*16 :: summ

        !        Nm = size(y,1)
        !        allocate(x(Nm))

        do ii = Nm,1,-1
            summ = y_out(ii)
            ind_ii_dep = (ii-1)*Nm -rshift(ii*(ii-1),1)
            do k= Nm,ii+1,-1
                summ = summ - conjg(Lm(ind_ii_dep+k))*x(k)
            enddo
            x(ii) = summ
        enddo

    end subroutine backward_subst

    ! -----------------------------------------------------------------------
    ! Function: mat2vec
    ! Purpose : Maps 2D lower-triangular index (i,j) with i>=j to the
    !           corresponding 1D packed-storage index.
    !           indx = (j-1)*N + i - j*(j-1)/2
    ! -----------------------------------------------------------------------
    function mat2vec(Nm,i,j) result(indx)
        integer :: Nm,i,j,indx

        indx = (j-1)*Nm+i-rshift(j*(j-1),1)
    end function mat2vec

    ! -----------------------------------------------------------------------
    ! Function: solve_LU_Cholesky
    ! Purpose : Full LDL* solve: Y*x = b.
    !   If use_old=.false. : calls get_LU_Cholesky to factorise Y first,
    !                        then deallocates Y.
    !   If use_old=.true.  : reuses pre-existing Lm, Dm factors.
    !   Calls forward_subst then backward_subst.
    ! Inputs : Y (packed Hermitian), b_vector, Nm (dimension), use_old
    ! Output : x_out(Nm) -- solution
    ! -----------------------------------------------------------------------
    function solve_LU_Cholesky(Y,b_vector,Nm,Lm,Dm,use_old) result(x_out)
        integer,intent(in) :: Nm
        complex*16,dimension(Nm) :: x_out
        complex*16,allocatable,intent(inout)::b_vector(:)
        complex*16,allocatable,intent(inout) :: Lm(:),Dm(:)
        complex*16,allocatable,dimension(:) :: y_out
        complex*16,allocatable,intent(inout) :: Y(:)
        logical :: use_old

        !        allocate(b_vector(Nm))
        !
        !        b_vector = matmul(conjg(transpose(Z_matrix)),Ex_vector)


        if(.not. use_old) then
            call get_LU_Cholesky(Y,Lm,Dm,Nm)
            deallocate(Y)
!            write(*,*) 'evaluating LU decomposition'
        endif
!        write(*,*) 'Lm',size(Lm,1),'Dm',size(Dm,1)
        if(.not. allocated(Lm)) then
            write(*,*) 'Error ... Lm,Dm are not allocated in function solve_LU_Cholesky'
            STOP
        endif

        call forward_subst(Lm,Dm,b_vector,Nm,y_out)
        deallocate(b_vector)
        call backward_subst(Lm,y_out,Nm,x_out)

        deallocate(y_out)
    end function solve_LU_Cholesky

    ! -----------------------------------------------------------------------
    ! Subroutine: get_LU_Cholesky
    ! Purpose   : Computes the LDL* (modified Cholesky) factorisation of a
    !             Hermitian positive-definite matrix Y stored in packed form.
    !             Produces diagonal D and unit lower-triangular L (packed).
    !             Uses a parallelised inner loop (OpenMP) for the diagonal
    !             update. Inner subroutine calculate_Lm_raw performs the
    !             off-diagonal column update.
    ! Input : Y(N*(N+1)/2) -- packed lower-triangular Hermitian matrix
    ! Output: Lm(N*(N+1)/2), Dm(N)
    ! -----------------------------------------------------------------------
    subroutine get_LU_Cholesky(Y,Lm,Dm,Nm)

        integer :: Nm,ii,jj,kk,ind_jj_dep
!        integer :: ind_kk_dep
        complex*16,allocatable,intent(inout) :: Y(:)
        complex*16,allocatable,intent(inout) :: Lm(:),Dm(:)
!        complex*16 :: summ
        complex*16,dimension(0:Nm) :: summ_array
!        complex*16,allocatable,dimension(:,:) :: summ_mat

        !        Nm = (sqrt(1.+8.*size(Y,1))-1.)/2.
        allocate(Dm(Nm),Lm(size(Y,1)))
!        allocate(summ_array(0:Nm))
!        allocate(summ_mat(Nm-1,Nm-1))

!        summ_mat = 0.d0

        Dm(1) = Y(1)
        do ii=2,Nm
            Lm((ii-1)*Nm+ii-rshift(ii*(ii-1),1)) = 1.d0

!            do jj = 1,ii-1
!                ind_jj_dep = (jj-1)*Nm+ii-rshift(jj*(jj-1),1)
!                summ = Y(ind_jj_dep)
!                do kk = 1,jj-1
!                    ind_kk_dep = (kk-1)*Nm-rshift(kk*(kk-1),1)
!                    summ = summ - Dm(kk)*Lm(ind_kk_dep+ii)*conjg(Lm(ind_kk_dep+jj))
!                enddo
!                Lm(ind_jj_dep) = summ/Dm(jj);
!            enddo
!            !$OMP PARALLEL do private(jj,ind_jj_dep)
            do jj = 1,ii-1
                ind_jj_dep = (jj-1)*Nm+ii-rshift(jj*(jj-1),1)

                call calculate_Lm_raw(ii,jj,ind_jj_dep)
            enddo
!            !$OMP END PARALLEL do
    !            summ_array = 0.d0
    !!            !$OMP PARALLEL do private(jj,ind_jj_dep,kk,ind_kk_dep,summ_array)
    !            do jj = 1,ii-1
    !
    !                ind_jj_dep = (jj-1)*Nm+ii-rshift(jj*(jj-1),1)
    !!                summ_array = 0.d0
    !!                summ = Y(ind_jj_dep)
    !                summ_array(1) = Y(ind_jj_dep)
    !                !$OMP PARALLEL do private(kk,ind_kk_dep)
    !                do kk = 1,jj-1
    !!                    write(*,*) jj,kk,omp_get_num_threads()
    !                    ind_kk_dep = (kk-1)*Nm-rshift(kk*(kk-1),1)
    !!                    summ_mat(jj,kk+1) = -Dm(kk)*Lm(ind_kk_dep+ii)*conjg(Lm(ind_kk_dep+jj))
    !                    summ_array(kk+1) = - Dm(kk)*Lm(ind_kk_dep+ii)*conjg(Lm(ind_kk_dep+jj))
    !                enddo
    !                !$OMP END PARALLEL do
    !!                write(*,*)
    !
    !                Lm(ind_jj_dep) = sum(summ_array(1:jj))/Dm(jj);
    !
    !!                !$OMP END CRITICAL
    !
    !
    !!                Lm(ind_jj_dep) = sum(summ_mat(jj,1:jj))/Dm(jj);
    !
    !
    !            enddo
!            !$OMP END PARALLEL do

!            summ = Y((ii-1)*Nm+ii-rshift(ii*(ii-1),1))
!            do kk = 1,ii-1
!                summ = summ - Dm(kk)*abs(Lm((kk-1)*Nm+ii-rshift(kk*(kk-1),1)))**2.d0
!            enddo
!            Dm(ii) = summ
            summ_array = 0.d0
            summ_array(0) = Y((ii-1)*Nm+ii-rshift(ii*(ii-1),1))
            !$OMP PARALLEL do private(kk)
            do kk = 1,ii-1
                summ_array(kk) = - Dm(kk)*abs(Lm((kk-1)*Nm+ii-rshift(kk*(kk-1),1)))**2.d0
            enddo
            !$OMP END PARALLEL do
            Dm(ii) = sum(summ_array)

        enddo
!        deallocate(summ_array)
        contains
        subroutine calculate_Lm_raw(ii,jj,ind_jj_dep)
            integer :: ind_jj_dep
            complex*16 :: summ
            integer :: kk,jj,ind_kk_dep,ii
            summ = Y(ind_jj_dep)
            do kk = 1,jj-1
                ind_kk_dep = (kk-1)*Nm-rshift(kk*(kk-1),1)
                summ = summ - Dm(kk)*Lm(ind_kk_dep+ii)*conjg(Lm(ind_kk_dep+jj))
            enddo
            Lm(ind_jj_dep) = summ/Dm(jj);
        end subroutine calculate_Lm_raw
    end subroutine get_LU_Cholesky

    ! -----------------------------------------------------------------------
    ! Function: matrix_mult
    ! Purpose : Hermitian matrix-vector product y = A*x for a matrix stored
    !           in packed lower-triangular form. Exploits symmetry: the
    !           conjugate of sub-diagonal elements fills the upper part.
    ! -----------------------------------------------------------------------
    function matrix_mult(A,x,Nm) result(y)
        complex*16,dimension(:),intent(in) :: A,x
        complex*16,dimension(Nm) :: y
        integer :: Nm,i,j
        complex*16 :: A_val

        y = 0.d0
        do i=1,Nm
            A_val = A(mat2vec(Nm,i,i))
            y(i) = y(i) + A_val*x(i)
            do j = 1,(i-1)
                A_val = A(mat2vec(Nm,i,j))
                y(i) = y(i) + A_val*x(j)
                y(j) = y(j) + conjg(A_val)*x(i)
            enddo
        enddo
!        do i=1,Nm
!            do j = (i+1),Nm
!                A_val = conjg(A(mat2vec(Nm,j,i)))
!                y(i) = y(i) + A_val*x(j)
!            enddo
!        enddo

    end function matrix_mult

    ! -----------------------------------------------------------------------
    ! Function: solve_Matrix_SVD
    ! Purpose : Solves Y*x = b via SVD-based pseudoinverse.
    !   If use_old=.false.: unpacks Y to full Y_mat, calls LAPACK ZGESVD,
    !     truncates singular values below S_max/MAX_ALLOWED_CONDITION_NUMBER,
    !     and stores the resulting Yinv = V * S^{-1} * U^H.
    !   If use_old=.true.: reuses stored Yinv.
    ! Use this solver for ill-conditioned systems (Matrix_Solution_Method=3).
    ! -----------------------------------------------------------------------
        function solve_Matrix_SVD(Y,b_vector,Nm,Y_mat,Yinv,use_old) result(x_out)
        integer,intent(in) :: Nm
        complex*16,dimension(Nm) :: x_out
        complex*16,allocatable,intent(inout)::b_vector(:)
        complex*16,allocatable,intent(inout) :: Yinv(:,:),Y_mat(:,:)
        complex*16,allocatable,intent(inout) :: Y(:)
        integer :: ii,jj,pp,ok
        logical :: use_old
        integer,parameter :: LWORK_factor = 66
        integer :: LWORK
        complex*16,dimension(Nm,Nm) :: U,VT,U_modified
        real(8),dimension(Nm) :: S
        integer,dimension(LWORK_factor*Nm) :: IWORK
        complex*16,dimension(LWORK_factor*Nm) :: WORK,RWORK
        real(8) :: S_min_allowed


        if(.not. use_old) then
            LWORK = LWORK_factor*Nm
            !! obtain the full matrix Y_mat
            allocate(Y_mat(Nm,Nm),Yinv(Nm,Nm))

            pp = 1
            do jj = 1,Nm
                do ii = jj,Nm
    !                 write(*,*) Nm,ii,jj
                    Y_mat(ii,jj) = Y(pp)
                    Y_mat(jj,ii) = conjg(Y(pp))
                    pp = pp+1
                enddo
            enddo

    !        write(*,*) 'done',size(Y),pp,use_old
            deallocate(Y)
!            call get_LU_Cholesky(Y,Lm,Dm,Nm)
!            call ZGESDD( 'A', Nm, Nm, Y_mat, Nm, S, U, Nm, VT, Nm, WORK, -1, RWORK, IWORK, ok )
!            write(*,*) WORK(1),LWORK
!            write(*,*) size(Y_mat),size(WORK),size(RWORK)
            call ZGESVD( 'A','A', Nm, Nm, Y_mat, Nm, S, U, Nm, VT, Nm, WORK, LWORK, RWORK, IWORK, ok )
!            write(*,*) 'evaluating LU decomposition'

            VT = transpose(conjg(VT))
            U = transpose(conjg(U))

            S_min_allowed = S(1)/MAX_ALLOWED_CONDITION_NUMBER
            U_modified = 0.d0
            do jj = 1,Nm
                if(S(jj) > S_min_allowed) then
                    U_modified(jj,:) = U(jj,:)/S(jj)
                else
                    write(*,*) '==================== INFORMATION ============================================='
                    write(*,'(A,E10.2,I5,I5)') ' Matrix Truncated (SVD-Based Inversion), MAX ALLOWED CONDITION NUMBER: ',&
                        MAX_ALLOWED_CONDITION_NUMBER
                    write(*,*) 'Matix Dimension: ',Nm,', Truncated Dimension',jj
                    write(*,*) '=============================================================================='
                    exit
                endif
            enddo
            Yinv = matmul(VT,U_modified)
            LWORK = LWORK_factor*Nm

        endif



        x_out = matmul(Yinv,b_vector)
!        do jj = 1,Nm
!            write(*,*) S(jj)
!        enddo
!        stop
!        write(*,*) 'Lm',size(Lm,1),'Dm',size(Dm,1)
!        if(.not. allocated(Lm)) then
!            write(*,*) 'Error ... Lm,Dm are not allocated in function solve_LU_Cholesky'
!            STOP
!        endif
!
!        call forward_subst(Lm,Dm,b_vector,Nm,y_out)
!        deallocate(b_vector)
!        call backward_subst(Lm,y_out,Nm,x_out)
!        write(*,*) 'Inversion SVD done'

    end function solve_Matrix_SVD

    ! -----------------------------------------------------------------------
    ! Function: solve_Matrix_preconditioned
    ! Purpose : Solves Y*x = b using LAPACK ZHESV with diagonal
    !           preconditioning. Rows/columns whose normalised diagonal falls
    !           below Y_mat_norm_limit are zeroed out (removed DOFs).
    !           For N_eq > 1, normalisation is done per equation group to
    !           handle mixed-scale multi-physics systems.
    !   If use_old=.true.: re-solves using the stored factored matrix Yinv
    !     via ZHETRS2, avoiding re-factorisation.
    ! Use this solver for dielectric / IBC problems (Matrix_Solution_Method=4).
    ! -----------------------------------------------------------------------
    function solve_Matrix_preconditioned(Y,b_vector,Nm,Y_mat,Yinv,use_old,y_mat_norm,IWORK,N_eq) result(x_out)
        integer,intent(in) :: Nm
        complex*16,dimension(Nm) :: x_out
        complex*16,allocatable,intent(inout)::b_vector(:)
        complex*16,allocatable,intent(inout) :: Yinv(:,:),Y_mat(:,:)
        complex*16,allocatable,intent(inout) :: Y(:)
        complex*16,allocatable,dimension(:)::b_vector_int
        integer :: ii,jj,pp,ok,num_truncated
        logical :: use_old
        integer,parameter :: LWORK_factor = 66
        integer :: LWORK,N_eq,N_calc,nn
        real(8),dimension(Nm),intent(inout) :: Y_mat_norm
        integer,dimension(Nm),intent(inout) :: IWORK
        complex*16,dimension(LWORK_factor*Nm) :: WORK
        real(8) :: Y_mat_max
        real(8),dimension(3*N_eq) :: Y_mat_norm_eq
        real(8),dimension(3) :: Y_mat_norm_eq_global

        allocate(b_vector_int(Nm))
        b_vector_int = b_vector
        if(.not. use_old) then
            LWORK = LWORK_factor*Nm
            !! obtain the full matrix Y_mat
            allocate(Y_mat(Nm,Nm),Yinv(Nm,Nm))

            pp = 1
            do jj = 1,Nm
                do ii = jj,Nm
    !                 write(*,*) Nm,ii,jj
                    Y_mat(ii,jj) = Y(pp)
                    Y_mat(jj,ii) = conjg(Y(pp))
                    pp = pp+1
                enddo
            enddo

    !        write(*,*) 'done',size(Y),pp,use_old
            deallocate(Y)
            Y_mat_max = maxval(maxval(abs(Y_mat),2))
            do jj = 1,Nm
                Y_mat_norm(jj) = abs(Y_mat(jj,jj))
            enddo
            if(N_eq == 1) then
                Y_mat_norm = Y_mat_norm/maxval(Y_mat_norm)
            else
                N_calc = int(dble(Nm)/dble(3*N_eq))

                do jj = 1,3*N_eq
                    Y_mat_norm_eq(jj) = maxval(Y_mat_norm(((jj-1)*N_calc+1):(jj*N_calc)))
                enddo
!                write(*,*) Y_mat_norm_eq
                do jj = 1,N_eq

                    do nn = 1,3
!                        write(*,*) (nn-1)*N_eq+1
                        Y_mat_norm_eq_global(nn) = Y_mat_norm_eq(jj -1 + (nn-1)*N_eq+1)
                    enddo
!                    write(*,*) Y_mat_norm_eq_global
                    do nn = 1,3
                        Y_mat_norm_eq(jj -1 + (nn-1)*N_eq+1) = maxval(Y_mat_norm_eq_global)
                    enddo
!                    Y_mat_norm_eq(jj:N_eq:(3*N_eq)) = maxval(Y_mat_norm_eq_global)
                enddo
!                 write(*,*) Y_mat_norm_eq
!                 stop
                do jj = 1,3*N_eq
                    Y_mat_norm(((jj-1)*N_calc+1):(jj*N_calc)) = Y_mat_norm(((jj-1)*N_calc+1):(jj*N_calc))/Y_mat_norm_eq(jj)
                enddo
!                write(*,*) N_calc
!                write(*,*) Y_mat_norm
!                stop
            endif
            num_truncated = 0
            do ii = 1,Nm
                if(Y_mat_norm(ii) <= Y_mat_norm_limit) then
                   Y_mat(:,ii) = 0;
                   Y_mat(ii,:) = 0;
                   Y_mat(ii,ii) = Y_mat_max;
                   b_vector_int(ii) = 0;
                   num_truncated = num_truncated + 1
               endif
            enddo

            write(*,*) '==================== INFORMATION ============================================='
            write(*,*) 'Matix Dimension: ',Nm,', Truncated Elements',num_truncated
            write(*,*) '=============================================================================='

!            call ZGESV( Nm, 1,Y_mat, Nm, IWORK, b_vector, Nm, ok )
!            b_vector_int1 = b_vector_int
            call zhesv ('L', Nm, 1, Y_mat, Nm, IWORK, b_vector_int, Nm, WORK, LWORK, ok)
!            call ZGESVD( 'A','A', Nm, Nm, Y_mat, Nm, S, U, Nm, VT, Nm, WORK, LWORK, RWORK, IWORK, ok )



            Yinv = Y_mat

            x_out = b_vector_int
!            write(*,*) x_out(1:20)
!            CALL zhetrs2( 'L',Nm,1,Yinv,Nm,IWORK,b_vector_int1,Nm,work,ok )
!
!            write(*,*) b_vector_int1(1:20)
!            stop
        else
!            Y_mat_max = minval(minval(abs(Y_mat),2))
!            do jj = 1,Nm
!                Y_mat_norm(jj) = norm2(abs(Y_mat(:,jj)))
!            enddo
!            Y_mat_norm = Y_mat_norm/maxval(Y_mat_norm)
            do ii = 1,Nm
                if(Y_mat_norm(ii) <= Y_mat_norm_limit) then
                   b_vector_int(ii) = 0;
               endif
            enddo

!            write(*,*) 'Resolving',allocated(Yinv)
            CALL zhetrs2( 'L',Nm,1,Yinv,Nm,IWORK,b_vector_int,Nm,work,ok )
!            CALL zhetrs( 'L',Nm,1,Yinv,Nm,IWORK,b_vector_int,Nm,ok )
!            call forward_subst_matrix(Yinv,b_vector_int,Nm)
            x_out = b_vector_int
!            write(*,*) 'Resolving done'
        endif





        deallocate(b_vector_int)
    end function solve_Matrix_preconditioned

    ! -----------------------------------------------------------------------
    ! Function: solve_Matrix_partitioning
    ! Purpose : Block-partitioned ZHESV solver for multi-equation systems.
    !           Splits the N x N system into N_eq blocks of sizes N_1 and N_2
    !           per block-row, solves each block independently, and iterates.
    !           Used when N_eq > 1 (coupled E/H boundary conditions).
    ! -----------------------------------------------------------------------
    function solve_Matrix_partitioning(Y,b_vector,Nm,Y_mat,Yinv,use_old,y_mat_norm,IWORK,N_eq,N_1,N_2) result(x_out)
        integer,intent(in) :: Nm
        complex*16,dimension(Nm) :: x_out
        complex*16,allocatable,intent(inout)::b_vector(:)
        complex*16,allocatable,intent(inout) :: Yinv(:,:),Y_mat(:,:)
        complex*16,allocatable,intent(inout) :: Y(:)
        complex*16,allocatable,dimension(:) :: b_vector_update
        integer :: N_eq,ii,jj,pp,ok
        logical :: use_old
        integer,parameter :: LWORK_factor = 66
        integer :: LWORK,N_1,N_2,tt
        real(8),dimension(Nm),intent(inout) :: Y_mat_norm
        integer,dimension(Nm),intent(inout) :: IWORK
        integer,dimension(:),allocatable :: IWORK_1,IWORK_2
        complex*16,dimension(LWORK_factor*Nm) :: WORK
        complex*16,allocatable,dimension(:,:) :: Y_mat_1,Y_mat_2
        complex*16,allocatable,dimension(:) :: x_1,x_2
        real(8) :: norm_old !,norm_new

        allocate(Y_mat_1(N_1,N_1),Y_mat_2(N_2,N_2))
        allocate(IWORK_1(N_1),IWORK_2(N_2))
        allocate(x_1(N_1),x_2(N_2),b_vector_update(Nm))

        if(.not. use_old) then
            LWORK = LWORK_factor*Nm
            !! obtain the full matrix Y_mat
            allocate(Y_mat(Nm,Nm),Yinv(Nm,Nm))

            pp = 1
            do jj = 1,Nm
                do ii = jj,Nm
                    Y_mat(ii,jj) = Y(pp)
                    Y_mat(jj,ii) = conjg(Y(pp))
                    pp = pp+1
                enddo
            enddo

            deallocate(Y)

            if(N_eq > 1) then
                Yinv = Y_mat
!                write(*,*) size(Y_mat,1),size(Y_mat,2)
                do jj = 1,3

!                    write(*,*) ((jj-1)*N_1+1+(jj-1)*N_2),(jj*N_1+(jj-1)*N_2),(jj*N_1+1+(jj-1)*N_2),(jj*N_1+jj*N_2)

                    Y_mat_1 = Y_mat(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2),&
                        ((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2))

!                    write(*,*) 'done_1'

                    Y_mat(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2),&
                        ((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2)) = 0.d0

!                    write(*,*) 'done_2'

                    Y_mat_2 = Y_mat((jj*N_1+1+(jj-1)*N_2):(jj*N_1+jj*N_2),&
                        (jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2))

!                    write(*,*) 'done_3'

                    Y_mat((jj*N_1+1+(jj-1)*N_2):(jj*N_1+jj*N_2),&
                        (jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2)) = 0.d0

!                    write(*,*) 'done_4'

                    x_1 = b_vector(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2))

!                    write(*,*) 'done_5'

                    call zhesv ('L', N_1, 1, Y_mat_1, N_1, IWORK_1,x_1, N_1, WORK, LWORK, ok)

!                    write(*,*) ok,'done_6'

                    IWORK(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2)) = IWORK_1


                    Yinv(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2),&
                        ((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2)) = Y_mat_1

                    x_out(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2)) = x_1

!                    write(*,*) x_1
                    x_2 =  b_vector((jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2))
                    call zhesv ('L', N_2, 1, Y_mat_2, N_2, IWORK_2,x_2, N_2, WORK, LWORK, ok)
                    IWORK((jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2)) = IWORK_2

                    Yinv((jj*N_1+1+(jj-1)*N_2):(jj*N_1+jj*N_2),&
                        (jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2)) = Y_mat_2

                    x_out((jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2)) = x_2


                enddo
                tt = 1

!                write(*,*) Y_mat(60:80,60:80)
!                stop
                do while(tt < 100)

!                    write(*,*) b_vector(1:5)
                    b_vector_update = -matmul(Y_mat,x_out)
!                    write(*,*) b_vector_update(1:5)
                    b_vector_update = b_vector_update + b_vector
!                    write(*,*) b_vector_update(1:5)
                    norm_old = norm2(abs(b_vector_update))

                    write(*,*) norm_old
                    write(*,*) 'Not a reliable solutions method, aborted'
                    stop
                    do jj = 1,3


                        Y_mat_1 = Yinv(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2),&
                            ((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2))


                        Y_mat_2 = Yinv((jj*N_1+1+(jj-1)*N_2):(jj*N_1+jj*N_2),&
                            (jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2))


                        x_1 = b_vector_update(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2))
                        IWORK_1 = IWORK(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2))
!                        call zhesv ('L', N_1, 1, Y_mat_1, N_1, IWORK_1,x_1, N_1, WORK, LWORK, ok)
                        CALL zhetrs2( 'L',N_1, 1, Y_mat_1, N_1, IWORK_1,x_1, N_1,work,ok )

                        x_out(((jj-1)*N_1+1+(jj-1)*N_2):(jj*N_1+(jj-1)*N_2)) = x_1

                        x_2 =  b_vector_update((jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2))
                        IWORK_2 = IWORK((jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2))
!                        call zhesv ('L', N_2, 1, Y_mat_2, N_2, IWORK_2,x_2, N_2, WORK, LWORK, ok)
                        CALL zhetrs2( 'L',N_2, 1, Y_mat_2, N_2, IWORK_2,x_2, N_2,work,ok )


                        x_out((jj*N_1+1+(jj-1)*N_2):(jj*N_1+(jj)*N_2)) = x_2


                    enddo
                    tt = tt + 1
                enddo



            else


                call zhesv ('L', Nm, 1, Y_mat, Nm, IWORK, b_vector, Nm, WORK, LWORK, ok)
                Yinv = Y_mat
                x_out = b_vector
            endif
        else


            if(N_eq > 1) then

            else

                CALL zhetrs2( 'L',Nm,1,Yinv,Nm,IWORK,b_vector,Nm,work,ok )
                x_out = b_vector
            endif
        endif


        deallocate(Y_mat_1,Y_mat_2)
        deallocate(IWORK_1,IWORK_2,x_1,x_2,b_vector_update)


    end function solve_Matrix_partitioning


    ! -----------------------------------------------------------------------
    ! Function: solve_Matrix_SVD_preconditioned
    ! Purpose : Combined diagonal preconditioning + SVD pseudoinverse.
    !           Zeros out weak DOFs (Y_mat_norm_limit), then calls ZGESVD
    !           with singular-value truncation (MAX_ALLOWED_CONDITION_NUMBER).
    !           Most robust option; use when both conditioning and near-zero
    !           DOFs are problematic.
    ! -----------------------------------------------------------------------
   function solve_Matrix_SVD_preconditioned(Y,b_vector,Nm,Y_mat,Yinv,use_old,y_mat_norm,IWORK,N_eq) result(x_out)
        integer,intent(in) :: Nm
        complex*16,dimension(Nm) :: x_out
        complex*16,allocatable,intent(inout)::b_vector(:)
        complex*16,allocatable,intent(inout) :: Yinv(:,:),Y_mat(:,:)
        complex*16,allocatable,intent(inout) :: Y(:)
        integer :: ii,jj,pp,ok
        logical :: use_old
        integer,parameter :: LWORK_factor = 66
        integer :: LWORK,N_eq,num_truncated,N_calc,nn
        complex*16,dimension(Nm,Nm) :: U,VT,U_modified
        real(8),dimension(Nm) :: S
        real(8),dimension(Nm),intent(inout) :: Y_mat_norm
        integer,dimension(Nm),intent(inout) :: IWORK
        complex*16,dimension(LWORK_factor*Nm) :: WORK,RWORK
        real(8) :: S_min_allowed,Y_mat_max
        real(8),dimension(3*N_eq) :: Y_mat_norm_eq
        real(8),dimension(3) :: Y_mat_norm_eq_global


        if(.not. use_old) then
            LWORK = LWORK_factor*Nm
            !! obtain the full matrix Y_mat
            allocate(Y_mat(Nm,Nm),Yinv(Nm,Nm))

            pp = 1
            do jj = 1,Nm
                do ii = jj,Nm
    !                 write(*,*) Nm,ii,jj
                    Y_mat(ii,jj) = Y(pp)
                    Y_mat(jj,ii) = conjg(Y(pp))
                    pp = pp+1
                enddo
            enddo

    !        write(*,*) 'done',size(Y),pp,use_old
            deallocate(Y)
            Y_mat_max = maxval(maxval(abs(Y_mat),2))

            do jj = 1,Nm
                Y_mat_norm(jj) = abs(Y_mat(jj,jj))
            enddo
            if(N_eq == 1) then
                Y_mat_norm = Y_mat_norm/maxval(abs(Y_mat_norm))
            else
                N_calc = int(dble(Nm)/dble(3*N_eq))

                do jj = 1,3*N_eq
                    Y_mat_norm_eq(jj) = maxval(abs(Y_mat_norm(((jj-1)*N_calc+1):(jj*N_calc))))
                enddo
!                write(*,*) Y_mat_norm_eq
                do jj = 1,N_eq

                    do nn = 1,3
!                        write(*,*) (nn-1)*N_eq+1
                        Y_mat_norm_eq_global(nn) = Y_mat_norm_eq(jj -1 + (nn-1)*N_eq+1)
                    enddo
!                    write(*,*) Y_mat_norm_eq_global
                    do nn = 1,3
                        Y_mat_norm_eq(jj -1 + (nn-1)*N_eq+1) = maxval(Y_mat_norm_eq_global)
                    enddo
!                    Y_mat_norm_eq(jj:N_eq:(3*N_eq)) = maxval(Y_mat_norm_eq_global)
                enddo
!                 write(*,*) Y_mat_norm_eq
!                 stop
                do jj = 1,3*N_eq
                    Y_mat_norm(((jj-1)*N_calc+1):(jj*N_calc)) = Y_mat_norm(((jj-1)*N_calc+1):(jj*N_calc))/Y_mat_norm_eq(jj)
                enddo
!                write(*,*) N_calc
!                write(*,*) Y_mat_norm
!                stop
            endif
            num_truncated = 0
            do ii = 1,Nm
                if(abs(Y_mat_norm(ii)) <= Y_mat_norm_limit) then
                   Y_mat(:,ii) = 0;
                   Y_mat(ii,:) = 0;
                   Y_mat(ii,ii) = Y_mat_max;
                   b_vector(ii) = 0;
                   num_truncated = num_truncated + 1
               endif
            enddo

            write(*,*) '==================== INFORMATION ============================================='
            write(*,*) 'Matix Dimension: ',Nm,', Truncated Elements',num_truncated
            write(*,*) '=============================================================================='


!            do jj = 1,Nm
!                Y_mat_norm(jj) = norm2(abs(Y_mat(:,jj)))
!            enddo
!            Y_mat_norm = Y_mat_norm/maxval(Y_mat_norm)
!            do ii = 1,Nm
!                if(Y_mat_norm(ii) <= Y_mat_norm_limit) then
!                   Y_mat(:,ii) = 0;
!                   Y_mat(ii,:) = 0;
!                   Y_mat(ii,ii) = Y_mat_max;
!                   b_vector(ii) = 0;
!               endif
!            enddo
!            call get_LU_Cholesky(Y,Lm,Dm,Nm)
!            call ZGESDD( 'A', Nm, Nm, Y_mat, Nm, S, U, Nm, VT, Nm, WORK, -1, RWORK, IWORK, ok )
!            write(*,*) WORK(1),LWORK
!            write(*,*) size(Y_mat),size(WORK),size(RWORK)
            call ZGESVD( 'A','A', Nm, Nm, Y_mat, Nm, S, U, Nm, VT, Nm, WORK, LWORK, RWORK, IWORK, ok )
!            write(*,*) 'evaluating LU decomposition'

            VT = transpose(conjg(VT))
            U = transpose(conjg(U))

            S_min_allowed = S(1)/MAX_ALLOWED_CONDITION_NUMBER
            U_modified = 0.d0
            do jj = 1,Nm
                if(S(jj) > S_min_allowed) then
                    U_modified(jj,:) = U(jj,:)/S(jj)
                else
                    write(*,*) '==================== INFORMATION ============================================='
                    write(*,'(A,E10.2,I5,I5)') ' Matrix Truncated (SVD-Based Inversion), MAX ALLOWED CONDITION NUMBER: ',&
                        MAX_ALLOWED_CONDITION_NUMBER
                    write(*,*) 'Matix Dimension: ',Nm,', Truncated Dimension',jj
                    write(*,*) '=============================================================================='
                    exit
                endif
            enddo
            Yinv = matmul(VT,U_modified)
            LWORK = LWORK_factor*Nm

        else

            do ii = 1,Nm
                if(abs(Y_mat_norm(ii)) <= Y_mat_norm_limit) then
                   b_vector(ii) = 0;
               endif
            enddo
        endif



        x_out = matmul(Yinv,b_vector)


    end function solve_Matrix_SVD_preconditioned

end module Matrices_Storage_Handling
