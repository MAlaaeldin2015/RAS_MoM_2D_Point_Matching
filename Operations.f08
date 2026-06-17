! ============================================================================
! Module: Operations
! Purpose: Mathematical and electromagnetic kernel library.
!          Provides:
!   - Overloaded vector arithmetic for real (type vector) and complex
!     (type vector_c) 3-component vectors: +, -, *, dot, cross, absolute,
!     norm, vec constructors.
!   - Hankel and Bessel function evaluation:
!       besselh2_0/1   -- H_0^(2)(x), H_1^(2)(x) for real argument
!       besselh2_01    -- both orders in one call (real or complex argument)
!       hankarray      -- array of H_n^(2)(x) and H_n^(2)'(x) via recurrence
!       besselarray    -- array of J_n(z) and J_n'(z) for complex argument
!       Hankel_derivatives  -- d^n/dx^n [x*H_p(a*x)], used in AWE kernels
!       Hankel_derivatives2 -- d^n/dx^n [H_p(x)], used in AWE excitation
!       CJY01          -- J_0, J_1, Y_0, Y_1 for complex argument (series)
!   - Chebyshev derivative matrix: chebyshevT_D
!   - Adaptive Simpson quadrature: module Simpson_Quad (Integrate, Quadstep,
!     Integrate_self)
!   - Utilities: factorial, num2str, ceil, dist, init_random_seed,
!     check_approx, vector2array
!
! Author : Mohamed A. Moharram Hassan
! License: Apache 2.0 -- see LICENSE
! ============================================================================
module Operations
    use constants
    implicit none

    integer, parameter:: fnum=7
    type multival
        complex*16, dimension(fnum):: f
    end type multival

    type vector
        real(8),dimension(3) :: v
    end type vector

    type vector_c
    complex*16,dimension(3) :: v
    end type vector_c



    interface operator (+)
        module procedure plus_r,plus_c,plus_multival,plus_vc,plus_vvc
    end interface
    interface operator (-)
        module procedure minus,minus_multival,minus_c
    end interface
    interface operator (*)
        module procedure mult, mult_c, mult_rc, mult_cr,mult_multival,mult_vc,mult_rca
    end interface

    interface operator(==)
        module procedure check_equal_vector_r
    end interface


    interface absolute
        module procedure absolute_r,absolute_c
    end interface absolute

    interface cross_prod
        module procedure cross_prod_c, cross_prod_r,cross_prod_v,cross_prod_vv
    end interface

    interface dot
        module procedure dot_r,dot_c,dot_rc,dot_rca,dot_rra
    end interface dot

    interface vec
        module procedure vec_r,vec_c,vec_ra,vec_ca
    end interface vec

    interface norm
        module procedure norm_ar,norm_ac
    end interface norm

    interface hankarray
        module procedure hankarray_r,hankarray_c
    end interface hankarray

    interface besselh2_01
        module procedure besselh2_01_r,besselh2_01_c
    end interface besselh2_01

    interface Hankel_derivatives2
        module procedure Hankel_derivatives2_r,Hankel_derivatives2_c
    end interface Hankel_derivatives2

    interface Hankel_derivatives
        module procedure Hankel_derivatives_r,Hankel_derivatives_c
    end interface Hankel_derivatives

    interface besselarray
        module procedure besselarray_c
    end interface besselarray

    contains

    ! -----------------------------------------------------------------------
    ! Function: factorial
    ! Purpose : Returns n! as a double-precision real. Used in Taylor-series
    !           and binomial coefficient calculations.
    ! -----------------------------------------------------------------------
    function factorial(n) result(y)
        integer :: n,i
        real(8) :: y

        y = 1.d0
        do i = 2,n
            y = y*dble(i)
        enddo

    end function factorial



    ! -----------------------------------------------------------------------
    ! Function: num2str
    ! Purpose : Converts a non-negative integer i to a character string of
    !           the given 'length'. Supports 1, 2, or 3-digit integers.
    !           Used for constructing output filenames.
    ! -----------------------------------------------------------------------
    function num2str(i,length) result(str)
        integer :: i,length
        character(length) ::str

        if(i < 10) then
            write(str,"(i1)") i

        elseif(i < 100) then
!            write(*,"(A5,i2)") 'here ',i
            write(str,"(i2)") i

        else
            write(str,"(i3)") i

        endif
!        str = trim(str)
!        write(*,*) i,str
    end function num2str


    ! -----------------------------------------------------------------------
    ! Subroutine: chebyshevT_D
    ! Purpose   : Computes the Chebyshev derivative matrix T(0:N, 0:N)
    !             where T(k,n) = d^k/dx^k [T_n(x)], T_n being the n-th
    !             Chebyshev polynomial of the first kind. Used when
    !             converting the AWE Taylor solution to a Chebyshev
    !             expansion for improved wideband accuracy.
    ! Inputs : N -- maximum derivative order; x -- expansion point
    ! Output : T(0:N, 0:N) -- derivative matrix (allocated on exit)
    ! -----------------------------------------------------------------------
    subroutine chebyshevT_D(N,x,T)
        real(8) :: x
        integer :: N
        real(8),allocatable,dimension(:,:) :: T
        real(8),allocatable,dimension(:) :: F,D,F_updated,D_init
        integer :: kk,ii,nn

        allocate(T(0:N,0:N))
        allocate(F(0:N),D_init(0:N),D(0:N),F_updated(0:N))
        T = 0.d0
        T(0,0) = 1;


        F = 0.d0
        F(0:1) = [1.d0,-(1.d0-x)];
        T(0,1) = x;
        T(1,1) = 1.d0;



        D_init = 0.d0
        D_init(1) = -1.d0/(1.d0-x)

        do kk = 1,(N-1)
           D_init(kk+1) = D_init(kk)*(kk+1)/(kk);
        enddo
        D = 0.d0
        do nn = 2,N
            F_updated = F
            do kk = 0,(nn-1)
                F_updated(kk) = F_updated(kk)*(nn+kk-1)/(nn-kk);
            enddo

            !-2*F_updated(nn)*(1-x)*(nn-kk)*(nn+kk)/((2*kk+2)*(2*kk+1));
            kk = kk - 1;
            F_updated(nn) = -2.d0*F_updated(nn-1)*(1-x)*dble((nn-kk)*(nn+kk))/dble((2*kk+2)*(2*kk+1));

            F = F_updated;
        !        write(*,*) F
            T(0,nn) = nn*sum(F);
            D(0:nn) = D_init(0:nn);
            do ii = 1,nn
    !            write(*,*) D
                T(ii,nn) = nn*sum(F*D);


                do kk = ii,nn
                   D(kk) = -D(kk)*dble(kk-ii)/(1.d0-x);
                enddo
            enddo
            D = 0.d0

        enddo
    !     write(*,*)

        deallocate(F,D,F_updated,D_init)
    end subroutine chebyshevT_D








    function check_equal_vector_r(vec1,vec2) result(ch)
        type(vector),intent(in) :: vec1,vec2
        logical :: ch
        if((vec1%v(1) == vec2%v(1)).and.(vec1%v(2) == vec2%v(2)).and. (vec1%v(3) == vec2%v(3))  ) then
            ch = .true.
        else
            ch = .false.
        endif
    end function check_equal_vector_r


    function norm_ar(a) result(y)
    real(8),dimension(:) :: a
    integer :: m
    real(8) :: y
    y = 0.d0
        do m=1,size(a,1)
            y=y+abs(a(m))
        enddo
    end function norm_ar

    function norm_ac(a) result(y)
    complex*16,dimension(:) :: a
    integer :: m
    real(8) :: y
    y = 0.d0
        do m=1,size(a,1)
            y=y+abs(a(m))
        enddo
    end function norm_ac

    type (multival) function subset(mvf1,i1,i2)
        type (multival), intent(in):: mvf1
        integer,intent(in) :: i1,i2
        integer:: n,c
        c = 1
        do n=i1, i2
            subset%f(c)=mvf1%f(n)
            c =c+1
        enddo
    end function subset

    type(multival) function plus_multival(mvf1,mvf2)
        type (multival), intent(in):: mvf1,mvf2
        integer:: n
        do n=1, fnum
            plus_multival%f(n)=mvf1%f(n)+mvf2%f(n)
        enddo
    end function plus_multival

    type(multival) function minus_multival(mvf1,mvf2)
        type (multival), intent(in):: mvf1,mvf2
        integer:: n
        do n=1, fnum
            minus_multival%f(n)=mvf1%f(n)-mvf2%f(n)
        enddo
    end function minus_multival

    type(multival) function mult_multival(s,mvf)
        type (multival), intent(in):: mvf
        real(8), intent(in):: s
        integer:: n
        do n=1, fnum
            mult_multival%f(n)=s*mvf%f(n)
        enddo
    end function mult_multival


    real(8) function dot_r(v1,v2)
        type(vector) :: v1,v2
        integer :: i
        dot_r = 0.d0
        do i=1,3
            dot_r = dot_r + v1%v(i)*v2%v(i)
        enddo
    end function dot_r

    complex*16 function dot_c(v1,v2)
        type(vector_c) :: v1,v2
        integer :: i
        dot_c = 0.d0
        do i=1,3
            dot_c = dot_c + v1%v(i)*v2%v(i)
        enddo
    end function dot_c

    complex*16 function dot_rc(v1,v2)
        type(vector_c) :: v2
        type(vector) :: v1
        integer :: i
        dot_rc = 0.d0
        do i=1,3
            dot_rc = dot_rc + v1%v(i)*v2%v(i)
        enddo
    end function dot_rc


    function dot_rca(v1,v2,N) result(dot_res)
        integer :: N
        type(vector_c),intent(in) :: v2(:)
        type(vector),intent(in) :: v1
        integer :: i
        complex*16,dimension(N) :: dot_res

        do i = 1,N
            dot_res(i) = v1%v(1)*v2(i)%v(1) + v1%v(2)*v2(i)%v(2)  + v1%v(3)*v2(i)%v(3)
        enddo
!        do i=1,3
!            dot_rc = dot_rc + v1%v(i)*v2%v(i)
!        enddo
    end function dot_rca


    function dot_rra(v1,v2,N) result(dot_res)
        integer :: N
        type(vector),intent(in) :: v2(:)
        type(vector),intent(in) :: v1
        integer :: i
        real(8),dimension(N) :: dot_res

        do i = 1,N
            dot_res(i) = v1%v(1)*v2(i)%v(1) + v1%v(2)*v2(i)%v(2)  + v1%v(3)*v2(i)%v(3)
        enddo
!        do i=1,3
!            dot_rc = dot_rc + v1%v(i)*v2%v(i)
!        enddo
    end function dot_rra

    type(vector) function vec_r(x,y,z)
        real(8) :: x,y,z
        vec_r%v(1) = x
        vec_r%v(2) = y
        vec_r%v(3) = z

    end function vec_r


    function vec_ra(x,y,z,N) result(vres)
        integer :: N
        real(8),dimension(N) :: x,y,z
        type(vector),dimension(N) :: vres
        vres(:)%v(1) = x
        vres(:)%v(2) = y
        vres(:)%v(3) = z

    end function vec_ra

    function vec_ca(x,y,z,N) result(vres)
        integer :: N
        complex*16,intent(in) :: x(:),y(:),z(:)
        type(vector_c),dimension(N) :: vres
        vres(:)%v(1) = x
        vres(:)%v(2) = y
        vres(:)%v(3) = z

    end function vec_ca

    type(vector_c) function vec_c(x,y,z)
        complex*16 :: x,y,z
        vec_c%v(1) = x
        vec_c%v(2) = y
        vec_c%v(3) = z

    end function vec_c

    function absolute_r(mvf) result(absolt)
    type(vector), intent(in) ::mvf
    real(8) :: absolt
    absolt = sqrt(mvf%v(1)**2.d0 + mvf%v(2)**2.d0 + mvf%v(3)**2.d0)
    end function absolute_r

    function absolute_c(mvf) result(absolt)
    type(vector_c), intent(in) ::mvf
    real(8) :: absolt
    absolt = sqrt(abs(mvf%v(1))**2.d0 + abs(mvf%v(2))**2.d0 + abs(mvf%v(3))**2.d0)
    end function absolute_c

        type (vector), intent(in):: mvf1,mvf2

        plus_r%v=mvf1%v+mvf2%v

    end function plus_r

    function plus_vc(mvf1,mvf2) result(res)
        type (vector_c),dimension(:,:), intent(in):: mvf1,mvf2
        type(vector_c),dimension(size(mvf1,1),size(mvf1,2)) :: res
        integer :: i,nn

!        do i = 1,size(mvf1,1)
!        res(1)=mvf1%v+mvf2%v
!        enddo
        do i = 1,size(mvf1,1)
            do nn = 1,size(mvf1,2)
                res(i,nn)%v = mvf1(i,nn)%v+mvf2(i,nn)%v
            enddo
        enddo

    end function plus_vc

    function plus_vvc(mvf1,mvf2) result(res)
        type (vector_c),dimension(:), intent(in):: mvf1,mvf2
        type(vector_c),dimension(size(mvf1,1)) :: res
        integer :: i
        do i = 1,size(mvf1,1)

            res(i)%v = mvf1(i)%v+mvf2(i)%v

        enddo

    end function plus_vvc


        type (vector_c), intent(in):: mvf1,mvf2

        plus_c%v=mvf1%v+mvf2%v

    end function plus_c

    type (vector) function minus(mvf1,mvf2)
        type (vector), intent(in):: mvf1,mvf2

        minus%v=mvf1%v-mvf2%v

    end function minus

    type (vector_c) function minus_c(mvf1,mvf2)
        type (vector_c), intent(in):: mvf1,mvf2

        minus_c%v=mvf1%v-mvf2%v

    end function minus_c

    type(vector) function mult(s,mvf)
        type (vector), intent(in):: mvf
        real(8), intent(in):: s

        mult%v= s * mvf%v

    end function mult

    type(vector_c) function mult_c(s,mvf)
        type (vector_c), intent(in):: mvf
        complex*16, intent(in):: s

        mult_c%v= s * mvf%v

    end function mult_c

    type(vector_c) function mult_cr(s,mvf)
        type (vector), intent(in):: mvf
        complex*16, intent(in):: s

            mult_cr%v= (s * mvf%v)

    end function mult_cr

    type(vector_c) function mult_rc(s,mvf)
        type (vector_c), intent(in):: mvf
        real(8), intent(in):: s

        mult_rc%v= (s * mvf%v)

    end function mult_rc

    type(vector_c) function mult_vc(mvf1,mvf)
        type (vector_c), intent(in):: mvf,mvf1

        mult_vc%v= (mvf1%v * mvf%v)

    end function mult_vc

    function mult_rca(s,v_cca) result(r_cca)
        type(vector_c),dimension(:),intent(in) :: v_cca
        real(8),intent(in) :: s
        type(vector_c),dimension(size(v_cca,1)) :: r_cca

        r_cca(1:size(v_cca,1))%v(1) = s*v_cca(1:size(v_cca,1))%v(1)
        r_cca(1:size(v_cca,1))%v(2) = s*v_cca(1:size(v_cca,1))%v(2)
        r_cca(1:size(v_cca,1))%v(3) = s*v_cca(1:size(v_cca,1))%v(3)

    !        r_cca =

    end function mult_rca

    function cross_prod_c(A_v,B_v) result(C_v)
    complex*16,dimension(3) :: B_v,C_v
    real(8),dimension(3) :: A_v

    C_v(1) = A_v(2)*B_v(3) - A_v(3)*B_v(2) !! AyBz - AzBy
    C_v(2) = A_v(3)*B_v(1) - A_v(1)*B_v(3) !! AzBx - AxBz
    C_v(3) = A_v(1)*B_v(2) - A_v(2)*B_v(1) !! AxBy - AyBx

    end function cross_prod_c

    function cross_prod_r(A_v,B_v) result(C_v)
    real(8),dimension(3) :: B_v,C_v
    real(8),dimension(3) :: A_v

    C_v(1) = A_v(2)*B_v(3) - A_v(3)*B_v(2) !! AyBz - AzBy
    C_v(2) = A_v(3)*B_v(1) - A_v(1)*B_v(3) !! AzBx - AxBz
    C_v(3) = A_v(1)*B_v(2) - A_v(2)*B_v(1) !! AxBy - AyBx

    end function cross_prod_r

    function cross_prod_v(A_v,B_v) result(C_v)
    type(vector) :: B_v,C_v
    type(vector) :: A_v

    C_v%v(1) = A_v%v(2)*B_v%v(3) - A_v%v(3)*B_v%v(2) !! AyBz - AzBy
    C_v%v(2) = A_v%v(3)*B_v%v(1) - A_v%v(1)*B_v%v(3) !! AzBx - AxBz
    C_v%v(3) = A_v%v(1)*B_v%v(2) - A_v%v(2)*B_v%v(1) !! AxBy - AyBx

    end function cross_prod_v

    function cross_prod_vv(A_v,B_v) result(C_v)
    type(vector_c) :: B_v,C_v
    type(vector) :: A_v

    C_v%v(1) = A_v%v(2)*B_v%v(3) - A_v%v(3)*B_v%v(2) !! AyBz - AzBy
    C_v%v(2) = A_v%v(3)*B_v%v(1) - A_v%v(1)*B_v%v(3) !! AzBx - AxBz
    C_v%v(3) = A_v%v(1)*B_v%v(2) - A_v%v(2)*B_v%v(1) !! AxBy - AyBx

    end function cross_prod_vv

    function vector2array(mvf) result(a)
        type (vector), intent(in):: mvf
        real(8),dimension(3) :: a

        a = mvf%v

    end function vector2array

    function ceil(real_arg) result(int_arg)
        real(8),intent(in) :: real_arg
        integer :: int_arg
        integer :: temp_int
        real(8) :: temp_real

        temp_int = int(real_arg)
        temp_real = real_arg - dble(temp_int)
        if(temp_real > 0.0) then
            int_arg = temp_int + 1
        else
            int_arg = temp_int
        endif
    end function ceil

    function get_R(v2,v1) result(R) !v2 is distination, v1 is the source
        type(vector) :: v1,v2,R
        R = v2 - v1
    end function get_R

    function besselh2_1(x) result(y)
        ! calculates the hankel function of the zeros order of the second kind
        real(8) :: x
        complex*16 :: y
        y = bessel_j1(x) - cj*bessel_y1(x)
    end function besselh2_1

    function besselh1_1(x) result(y)
        ! calculates the hankel function of the zeros order of the second kind
        real(8) :: x
        complex*16 :: y
        y = bessel_j1(x) + cj*bessel_y1(x)
    end function besselh1_1

    function besselh1_0(x) result(y)
        ! calculates the hankel function of the zeros order of the second kind
        real(8) :: x
        complex*16 :: y
        y = bessel_j0(x) + cj*bessel_y0(x)
    end function besselh1_0

    function besselh2_0(x) result(y)
        ! calculates the hankel function of the zeros order of the second kind
        real(8) :: x
        complex*16 :: y
        y = bessel_j0(x) - cj*bessel_y0(x)
    end function besselh2_0
    function dist(a,b)
    real(8) :: dist
    real(8),dimension(3) :: a,b

    dist = sqrt((a(1)-b(1))**2.d0 + (a(2)-b(2))**2.d0 + (a(3)-b(3))**2.d0)

end function dist

SUBROUTINE init_random_seed()

        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
END SUBROUTINE


function check_approx(vec1,vec2) result(ch)
    type(vector),intent(in) :: vec1,vec2
    logical :: ch
    real(8) :: val

    val = dot(vec1,vec2)/(absolute(vec1)*absolute(vec2))

    if(abs(val) > 0.999d0) then
        ch = .true.
    else
        ch = .false.
    endif

end function check_approx



    ! -----------------------------------------------------------------------
    ! Subroutine: hankarray_r  (real argument)
    ! Purpose   : Computes arrays H_n^(2)(x) and d/dx[H_n^(2)(x)] for
    !             n = 0..n_terms using the upward recurrence relation
    !             H_{n+1}(x) = (2n/x)*H_n(x) - H_{n-1}(x)
    !             and the derivative recurrence H_n' = (H_{n-1}-H_{n+1})/2.
    ! Inputs : n_terms (max order), x (real argument)
    ! Output : h2(0:n_terms), h2p(0:n_terms)
    ! -----------------------------------------------------------------------
    subroutine hankarray_r(n_terms,x,h2,h2p)
        integer :: n_terms
        real(8) :: x
        integer :: i
        complex*16,dimension(0:n_terms) :: h2,h2p
!        real(8),dimension(0:(n_terms+1)) :: jx,yx
!        complex*16 :: hn,hnm1,hnp1
        complex*16 :: hn


        call besselh2_01(x,h2(0),h2(1))

        h2p(0) = -h2(1)

        do i = 1,(n_terms-1)
            h2(i+1) = 2.d0*dble(i)/x*h2(i) - h2(i-1)
            h2p(i) = (h2(i-1)-h2(i+1))/2.d0
        enddo
        i = n_terms
        hn = 2.d0*dble(i)/x*h2(i) - h2(i-1)
        h2p(i) = (h2(i-1)-hn)/2.d0
!        do i=0,(n_terms+1)
!            jx(i) = bessel_jn(i,x)
!            yx(i) = bessel_yn(i,x)
!        enddo
!
!        h2(0) = jx(0)-cj*yx(0)
!        hnm1 = h2(0)
!        h2(1) = jx(1)-cj*yx(1)
!        hn = h2(1)
!        h2p(0) = -h2(1)
!        do i=1,n_terms
!            h2(i) = hn
!            hnp1 = jx(i+1)-cj*yx(i+1)
!            h2p(i) = 0.5d0*(hnm1-hnp1)
!            hnm1 = hn
!            hn = hnp1
!
!        enddo

    end subroutine hankarray_r

    ! -----------------------------------------------------------------------
    ! Subroutine: besselarray_c  (complex argument)
    ! Purpose   : Computes J_n(z) and J_n'(z) for n = 0..n_terms using the
    !             same recurrence as hankarray but for Bessel J (needed for
    !             dielectric scatterer interior fields with complex k).
    ! Inputs : n_terms, z (complex)
    ! Output : J(0:n_terms), Jp(0:n_terms)
    ! -----------------------------------------------------------------------
    subroutine besselarray_c(n_terms,x,J,Jp)
        integer :: n_terms
        complex*16 :: x
        integer :: i
        complex*16,dimension(0:n_terms) :: J,Jp
        complex*16 :: Jn,Y0,Y1

        call CJY01(x,J(0),J(1),Y0,Y1)

        Jp(0) = -J(1)

        do i = 1,(n_terms-1)
            J(i+1) = 2.d0*dble(i)/x*J(i) - J(i-1)
            Jp(i) = (J(i-1)-J(i+1))/2.d0
        enddo
        i = n_terms
        Jn = 2.d0*dble(i)/x*J(i) - J(i-1)
        Jp(i) = (J(i-1)-Jn)/2.d0

    end subroutine besselarray_c

    subroutine hankarray_c(n_terms,x,h2,h2p)
        integer :: n_terms
        complex*16 :: x
        integer :: i
        complex*16,dimension(0:n_terms) :: h2,h2p
        complex*16 :: hn


        call besselh2_01(x,h2(0),h2(1))

        h2p(0) = -h2(1)

        do i = 1,(n_terms-1)
            h2(i+1) = 2.d0*dble(i)/x*h2(i) - h2(i-1)
            h2p(i) = (h2(i-1)-h2(i+1))/2.d0
        enddo
        i = n_terms
        hn = 2.d0*dble(i)/x*h2(i) - h2(i-1)
        h2p(i) = (h2(i-1)-hn)/2.d0



!        do i=0,(n_terms+1)
!            jx(i) = bessel_jn(i,x)
!            yx(i) = bessel_yn(i,x)
!        enddo
!
!        h2(0) = jx(0)-cj*yx(0)
!        hnm1 = h2(0)
!        h2(1) = jx(1)-cj*yx(1)
!        hn = h2(1)
!        h2p(0) = -h2(1)
!        do i=1,n_terms
!            h2(i) = hn
!            hnp1 = jx(i+1)-cj*yx(i+1)
!            h2p(i) = 0.5d0*(hnm1-hnp1)
!            hnm1 = hn
!            hn = hnp1
!
!        enddo

    end subroutine hankarray_c


    subroutine besselh2_01_c(z,h2_0,h2_1)
        complex*16 :: z,h2_0,h2_1
        complex*16 :: J0,J1,Y0,Y1

        call CJY01(z,J0,J1,Y0,Y1)

        h2_0 = J0 - cj*Y0
        h2_1 = J1 - cj*Y1

    end subroutine besselh2_01_c

    subroutine besselh2_01_r(z,h2_0,h2_1)
        real(8) ::z
        complex*16 :: h2_0,h2_1

        h2_0 = besselh2_0(z)
        h2_1 = besselh2_1(z)

    end subroutine besselh2_01_r


    ! -----------------------------------------------------------------------
    ! Subroutine: Hankel_derivatives2_r  (real argument)
    ! Author    : Mohamed A. Moharram  -- Feb. 25, 2013
    ! Purpose   : Computes d^n/dx^n [H_0^(2)(x)] and d^n/dx^n [H_1^(2)(x)]
    !             for n = 0..N using the chain-rule recurrence. These are
    !             the AWE frequency derivatives of the 2D Green's function
    !             when the argument depends linearly on frequency (kr).
    ! Inputs : x (real), N (max derivative order)
    ! Output : H0(0:N), H1(0:N)
    ! -----------------------------------------------------------------------
subroutine Hankel_derivatives2_r(x,N,H0,H1)
    !! this subroutine computes the derivatives of the hankel function that are required by MoM-AWE procedure
    !! these derivatives are for the functions y=H0(x) and y=H1(x) with respect to "x"
    !! Author :: Mohamed A. Moharram
    !! Date :: Feb. 25, 2013
        real(8),intent(in) :: x !! a is a constant and x is the argument to be differentiated w.r.t it
!        class(*),intent(in) :: x !! a is a constant and x is the argument to be differentiated w.r.t it
        integer,intent(in) :: N !! the maximum order of the required derivatives
        complex*16,dimension(0:N),intent(inout) :: H0,H1 !! these vectors are supposed to be allocated before invoking the routine
        !! dimensions of these arrays should be from 0 -> N
        complex*16,dimension(0:(N+1)) :: y0_n
        integer :: i,nn
        real(8) :: k_i


        y0_n(0) = besselh2_0(x)
        y0_n(1) = -besselh2_1(x)

        do nn = 2,(N+1)
            k_i = 1.d0/(x)
            y0_n(nn) = -y0_n(nn-2)
            do i= (nn-2),0,-1
                y0_n(nn) = y0_n(nn) - k_i*y0_n(i+1)
                k_i = -k_i*dble(i)/x
            enddo
        enddo


        H0 = y0_n(0:N)
        H1 = -y0_n(1:(N+1))


    end subroutine Hankel_derivatives2_r

    ! -----------------------------------------------------------------------
    ! Subroutine: Hankel_derivatives2_c  (complex argument)
    ! Author    : Mohamed A. Moharram  -- Apr. 1, 2013
    ! Purpose   : Complex-argument version of Hankel_derivatives2_r.
    !             Used for dielectric or IBC scatterers where the wavenumber
    !             k = ak * k0 is complex (lossy or double-negative media).
    ! -----------------------------------------------------------------------
    subroutine Hankel_derivatives2_c(x,N,H0,H1)
    !! this subroutine computes the derivatives of the hankel function that are required by MoM-AWE procedure
    !! these derivatives are for the functions y=H0(x) and y=H1(x) with respect to "x"
    !! Author :: Mohamed A. Moharram
    !! Date :: Apr. 1, 2013
        complex*16,intent(in) :: x !! a is a constant and x is the argument to be differentiated w.r.t it
!        class(*),intent(in) :: x !! a is a constant and x is the argument to be differentiated w.r.t it
        integer,intent(in) :: N !! the maximum order of the required derivatives
        complex*16,dimension(0:N),intent(inout) :: H0,H1 !! these vectors are supposed to be allocated before invoking the routine
        !! dimensions of these arrays should be from 0 -> N
        complex*16,dimension(0:(N+1)) :: y0_n
        integer :: i,nn
        complex*16 :: k_i

        call besselh2_01(x,y0_n(0),y0_n(1))
!        y0_n(0) = besselh2_0(x)
!        y0_n(1) = -besselh2_1(x)
        y0_n(1) = -y0_n(1)

        do nn = 2,(N+1)
            k_i = 1.d0/(x)
            y0_n(nn) = -y0_n(nn-2)
            do i= (nn-2),0,-1
                y0_n(nn) = y0_n(nn) - k_i*y0_n(i+1)
                k_i = -k_i*dble(i)/x
            enddo
        enddo


        H0 = y0_n(0:N)
        H1 = -y0_n(1:(N+1))


    end subroutine Hankel_derivatives2_c



    ! -----------------------------------------------------------------------
    ! Subroutine: Hankel_derivatives_r  (real argument)
    ! Author    : Mohamed A. Moharram  -- Jan. 15, 2013
    ! Purpose   : Computes the frequency derivatives of the 2D kernel
    !             functions x*H_0(a*x), x*H_1(a*x), x*H_2(a*x), and
    !             H_1(a*x) up to order N with respect to x = k*rho.
    !             These are the fundamental building blocks of the AWE
    !             matrix/excitation Taylor expansion.
    ! Inputs : a (wavenumber ratio), x (kr), N (max order)
    ! Output : xH0, xH1, xH2, H1  -- each of length 0:N
    ! -----------------------------------------------------------------------
    subroutine Hankel_derivatives_r(a,x,N,xH0,xH1,xH2,H1)
    !! this subroutine computes the derivatives of the hankel function that are required by the eval_kernel_* function
    !! these derivatives are for the functions y=x*H0(a*x),y=x*H1(a*x), y=x*H2(a*x), y=H1(a*x) with respect to "x"
    !! Author :: Mohamed A. Moharram
    !! Date :: Jan. 15, 2013
        real(8),intent(in) :: a,x !! a is a constant and x is the argument to be differentiated w.r.t it
        integer,intent(in) :: N !! the maximum order of the required derivatives
        complex*16,dimension(0:N),intent(inout) :: xH0,xH1,xH2,H1 !! these vectors are supposed to be allocated before invoking the routine
        !! dimensions of these arrays should be from 0 -> N
        complex*16,dimension(0:(N+1)) :: y0_n
        complex*16,dimension(0:N) :: H0,H2
        integer :: i,nn
        real(8) :: ax,k_i,psi_i
        real(8),dimension(0:N) :: a_power

        ax = a*x

        y0_n(0) = besselh2_0(ax)
        y0_n(1) = -besselh2_1(ax)

        do nn = 2,(N+1)
            k_i = 1.d0/(ax)
            y0_n(nn) = -y0_n(nn-2)
            do i= (nn-2),0,-1
                y0_n(nn) = y0_n(nn) - k_i*y0_n(i+1)
                k_i = -k_i*dble(i)/ax
            enddo
        enddo

        a_power(0) = 1.d0
        do i=1,N
            a_power(i) = a_power(i-1)*a
        enddo

        H0 = y0_n(0:N)*a_power
        H1 = -y0_n(1:(N+1))*a_power

        do nn = 0,N
            H2(nn) = -H0(nn)
            psi_i = 2.d0/(ax)
            do i=nn,0,-1
                H2(nn) = H2(nn) + psi_i*H1(i)
                psi_i = -psi_i*dble(i)/x
            enddo
        enddo

        xH0(0) = x*H0(0)
        xH1(0) = x*H1(0)
        xH2(0) = x*H2(0)
        do nn = 1,N
            xH0(nn) = x*H0(nn) + nn * H0(nn-1)
            xH1(nn) = x*H1(nn) + nn * H1(nn-1)
            xH2(nn) = x*H2(nn) + nn * H2(nn-1)
        enddo

!        do i = 0,N
!            write(*,*) H1(i),xH0(i),xH1(i),xH2(i)
!        enddo

    end subroutine Hankel_derivatives_r

    ! -----------------------------------------------------------------------
    ! Subroutine: Hankel_derivatives_c  (complex wavenumber)
    ! Author    : Mohamed A. Moharram  -- Apr. 1, 2013
    ! Purpose   : Complex-argument version of Hankel_derivatives_r.
    !             Used for dielectric/IBC scatterers with complex k.
    ! -----------------------------------------------------------------------
    subroutine Hankel_derivatives_c(a,x,N,xH0,xH1,xH2,H1)
    !! this subroutine computes the derivatives of the hankel function that are required by the eval_kernel_* function
    !! these derivatives are for the functions y=x*H0(a*x),y=x*H1(a*x), y=x*H2(a*x), y=H1(a*x) with respect to "x"
    !! Author :: Mohamed A. Moharram
    !! Date :: Apr. 1, 2013
        real(8),intent(in) :: x !! a is a constant and x is the argument to be differentiated w.r.t it
        complex*16,intent(in) :: a
        integer,intent(in) :: N !! the maximum order of the required derivatives
        complex*16,dimension(0:N),intent(inout) :: xH0,xH1,xH2,H1 !! these vectors are supposed to be allocated before invoking the routine
        !! dimensions of these arrays should be from 0 -> N
        complex*16,dimension(0:(N+1)) :: y0_n
        complex*16,dimension(0:N) :: H0,H2
        integer :: i,nn
        complex*16 :: ax,k_i,psi_i
        complex*16,dimension(0:N) :: a_power
        real(8) :: nnn

        ax = a*x

        call besselh2_01(ax,y0_n(0),y0_n(1))
        y0_n(1) = -y0_n(1)

!        y0_n(0) = besselh2_0(ax)
!        y0_n(1) = -besselh2_1(ax)

        do nn = 2,(N+1)
            k_i = 1.d0/(ax)
            y0_n(nn) = -y0_n(nn-2)
            do i= (nn-2),0,-1
                y0_n(nn) = y0_n(nn) - k_i*y0_n(i+1)
                k_i = -k_i*dble(i)/ax
            enddo
        enddo

        a_power(0) = 1.d0
        do i=1,N
            a_power(i) = a_power(i-1)*a
        enddo

        H0 = y0_n(0:N)*a_power
        H1 = -y0_n(1:(N+1))*a_power

        do nn = 0,N
            H2(nn) = -H0(nn)
            psi_i = 2.d0/(ax)
            do i=nn,0,-1
                H2(nn) = H2(nn) + psi_i*H1(i)
                psi_i = -psi_i*dble(i)/x
            enddo
        enddo

        xH0(0) = x*H0(0)
        xH1(0) = x*H1(0)
        xH2(0) = x*H2(0)
        do nn = 1,N
            nnn = dble(nn)
            xH0(nn) = x*H0(nn) + nnn * H0(nn-1)
            xH1(nn) = x*H1(nn) + nnn * H1(nn-1)
            xH2(nn) = x*H2(nn) + nnn * H2(nn-1)
        enddo

!        do i = 0,N
!            write(*,*) H1(i),xH0(i),xH1(i),xH2(i)
!        enddo

    end subroutine Hankel_derivatives_c


    ! -----------------------------------------------------------------------
    ! Subroutine: CJY01
    ! Purpose   : Computes Bessel and Neumann functions J_0(z), J_1(z),
    !             Y_0(z), Y_1(z) for complex argument z. Uses Taylor series
    !             for |z|<=12 and asymptotic expansion for |z|>12.
    !             Handles negative real-part of z by applying Hankel
    !             reflection formulas. Used by besselh2_01_c and
    !             all complex-argument Hankel/Bessel routines.
    !             Algorithm adapted from standard special-function libraries.
    ! Inputs : Z (complex)
    ! Output : CBJ0, CBJ1, CBY0, CBY1 (complex)
    ! -----------------------------------------------------------------------
SUBROUTINE CJY01(Z,CBJ0,CBJ1,CBY0,CBY1)
implicit none
!C
!C       =======================================================
!C       Purpose: Compute Bessel functions J0(z), J1(z), Y0(z),
!C                Y1(z), and their derivatives for a complex
!C                argument
!C       Input :  z --- Complex argument
!C       Output:  CBJ0 --- J0(z)
!C                CDJ0 --- J0'(z)
!C                CBJ1 --- J1(z)
!C                CDJ1 --- J1'(z)
!C                CBY0 --- Y0(z)
!C                CDY0 --- Y0'(z)
!C                CBY1 --- Y1(z)
!C                CDY1 --- Y1'(z)
!C       =======================================================
!C

        COMPLEX*16 :: Z
        complex*16 :: CBJ0,CBJ1,CBY0,CBY1
        real(8),DIMENSION(12) :: A,B,A1,B1
        real(8) :: EL,RP2,A0,W0,W1
        complex*16 :: Z2,Z1,CI,CR,CS,CP,CP0,CT1,CQ0,CU,CT2,CP1,CQ1
        integer :: k,K0
!        PI=3.141592653589793D0
        EL=0.5772156649015329D0
        RP2=2.0D0/PI
        CI=(0.0D0,1.0D0)
        A0=ABS(Z)
        Z2=Z*Z
        Z1=Z
        IF (A0.EQ.0.0D0) THEN
           CBJ0=(1.0D0,0.0D0)
           CBJ1=(0.0D0,0.0D0)
!           CDJ0=(0.0D0,0.0D0)
!           CDJ1=(0.5D0,0.0D0)
           CBY0=-(1.0D300,0.0D0)
           CBY1=-(1.0D300,0.0D0)
!           CDY0=(1.0D300,0.0D0)
!           CDY1=(1.0D300,0.0D0)
           RETURN
        ENDIF
        IF (REAL(Z).LT.0.0) Z1=-Z
        IF (A0.LE.12.0) THEN
           CBJ0=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO K=1,40
              CR=-0.25D0*CR*Z2/(K*K)
              CBJ0=CBJ0+CR
              IF (CDABS(CR).LT.CDABS(CBJ0)*1.0D-15) GO TO 15
           enddo
15         CBJ1=(1.0D0,0.0D0)
           CR=(1.0D0,0.0D0)
           DO K=1,40
              CR=-0.25D0*CR*Z2/(K*(K+1.0D0))
              CBJ1=CBJ1+CR
              IF (CDABS(CR).LT.CDABS(CBJ1)*1.0D-15) GO TO 25
           enddo
25         CBJ1=0.5D0*Z1*CBJ1
           W0=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(0.0D0,0.0D0)
           DO K=1,40
              W0=W0+1.0D0/K
              CR=-0.25D0*CR/(K*K)*Z2
              CP=CR*W0
              CS=CS+CP
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 35
           enddo
35         CBY0=RP2*(CDLOG(Z1/2.0D0)+EL)*CBJ0-RP2*CS
           W1=0.0D0
           CR=(1.0D0,0.0D0)
           CS=(1.0D0,0.0D0)
           DO K=1,40
              W1=W1+1.0D0/K
              CR=-0.25D0*CR/(K*(K+1))*Z2
              CP=CR*(2.0D0*W1+1.0D0/(K+1.0D0))
              CS=CS+CP
              IF (CDABS(CP).LT.CDABS(CS)*1.0D-15) GO TO 45
            enddo
45         CBY1=RP2*((CDLOG(Z1/2.0D0)+EL)*CBJ1-1.0D0/Z1-.25D0*Z1*CS)
        ELSE
           DATA A/-.703125D-01,.112152099609375D+00,&
                 -.5725014209747314D+00,.6074042001273483D+01,&
                 -.1100171402692467D+03,.3038090510922384D+04,&
                 -.1188384262567832D+06,.6252951493434797D+07,&
                 -.4259392165047669D+09,.3646840080706556D+11,&
                 -.3833534661393944D+13,.4854014686852901D+15/
           DATA B/ .732421875D-01,-.2271080017089844D+00,&
                  .1727727502584457D+01,-.2438052969955606D+02,&
                  .5513358961220206D+03,-.1825775547429318D+05,&
                  .8328593040162893D+06,-.5006958953198893D+08,&
                  .3836255180230433D+10,-.3649010818849833D+12,&
                  .4218971570284096D+14,-.5827244631566907D+16/
           DATA A1/.1171875D+00,-.144195556640625D+00,&
                  .6765925884246826D+00,-.6883914268109947D+01,&
                  .1215978918765359D+03,-.3302272294480852D+04,&
                  .1276412726461746D+06,-.6656367718817688D+07,&
                  .4502786003050393D+09,-.3833857520742790D+11,&
                  .4011838599133198D+13,-.5060568503314727D+15/
           DATA B1/-.1025390625D+00,.2775764465332031D+00,&
                  -.1993531733751297D+01,.2724882731126854D+02,&
                  -.6038440767050702D+03,.1971837591223663D+05,&
                  -.8902978767070678D+06,.5310411010968522D+08,&
                  -.4043620325107754D+10,.3827011346598605D+12,&
                  -.4406481417852278D+14,.6065091351222699D+16/
           K0=12
           IF (A0.GE.35.0) K0=10
           IF (A0.GE.50.0) K0=8
           CT1=Z1-.25D0*PI
           CP0=(1.0D0,0.0D0)
           DO K=1,K0
                CP0=CP0+A(K)*Z1**(-2*K)
           enddo
           CQ0=-0.125D0/Z1
           DO K=1,K0
            CQ0=CQ0+B(K)*Z1**(-2*K-1)
           enddo
           CU=CDSQRT(RP2/Z1)
           CBJ0=CU*(CP0*CDCOS(CT1)-CQ0*CDSIN(CT1))
           CBY0=CU*(CP0*CDSIN(CT1)+CQ0*CDCOS(CT1))
           CT2=Z1-.75D0*PI
           CP1=(1.0D0,0.0D0)
           DO K=1,K0
            CP1=CP1+A1(K)*Z1**(-2*K)
           enddo
           CQ1=0.375D0/Z1
           DO K=1,K0
            CQ1=CQ1+B1(K)*Z1**(-2*K-1)
           enddo
           CBJ1=CU*(CP1*CDCOS(CT2)-CQ1*CDSIN(CT2))
           CBY1=CU*(CP1*CDSIN(CT2)+CQ1*CDCOS(CT2))
        ENDIF
        IF (REAL(Z).LT.0.0) THEN
           IF (DIMAG(Z).LT.0.0) CBY0=CBY0-2.0D0*CI*CBJ0
           IF (DIMAG(Z).GE.0.0) CBY0=CBY0+2.0D0*CI*CBJ0
           IF (DIMAG(Z).LT.0.0) CBY1=-(CBY1-2.0D0*CI*CBJ1)
           IF (DIMAG(Z).GE.0.0) CBY1=-(CBY1+2.0D0*CI*CBJ1)
           CBJ1=-CBJ1
        ENDIF
!        CDJ0=-CBJ1
!        CDJ1=CBJ0-1.0D0/Z*CBJ1
!        CDY0=-CBY1
!        CDY1=CBY0-1.0D0/Z*CBY1

END SUBROUTINE CJY01


end module Operations


module Simpson_Quad
!*************************************************
! Uses Adaptive Simpson Quadrature method
!   to perform an n-tuple numerical integral
!*************************************************
    use Operations
    implicit none

    private Quadstep

    real(8), parameter :: MAX_REAL=1.d10
    contains

    recursive function Integrate (Func, nvi, lim, nvar, var, tol) result(Integrate_R)
    ! *******************************************************************************
    ! Integrates "Func" on "nvi" variables of integration
    ! lim: array holding the lower and upper limits of the variables of integration
    !       # lim(0) should be set to any arbitrary value
    !       # the n-th variable has the lower limit in lim(2n-1),
    !                                   upper limit in lim(2n)
    ! nvar: the number of passed parameters
    ! var : array holding the values of the passed parameters
    ! tol : required tolerance
    ! *******************************************************************************
    interface
        type (multival) function Func(n,x)
            use Operations
            integer, intent(in) :: n
            real(8), intent(in), dimension(n) :: x
        end function Func
    end interface
    integer, intent(in) :: nvi,nvar
    real(8), intent(in), dimension(0:2*nvi) :: lim
    real(8), intent(in), dimension(nvar) :: var
    real(8), intent(in) :: tol
    type (multival):: Integrate_R
    integer, parameter :: POINTS = 7
    real(8) :: h, hmin, x1, x2
    real(8), dimension(POINTS) :: xp
    type (multival), dimension(POINTS):: yp
    type (multival), dimension(3):: q
    integer :: n, nlims

    if (nvi==0) then

        Integrate_R=Func(nvar,(/var/))
        return
    endif

    x2=lim(2*nvi); x1=lim(2*nvi-1)
    nlims=2*(nvi-1)
    h=0.13579d0*(x2-x1)

    ! Seven point evaluation
    xp=(/x1,x1+h,x1+2.d0*h,(x1+x2)/2.d0,x2-2.d0*h,x2-h,x2/)
    do n=1,POINTS
        yp(n)=Integrate(Func, nvi-1, lim(0:nlims), nvar+1,(/xp(n),var/), tol )
    enddo
!    write(*,*) xp
    ! Avoid end-point singularities
    if (abs(yp(1)%f(1))>MAX_REAL .or. isnan(abs(yp(1)%f(1)) ) ) then
!        write(*,*) 'avoid end point singularity'
!        write(*,*) xp(1)
        x1 = x1+EPSILON(tol)*(x2-x1)
!        write(*,*) xp
        yp(1) = Integrate(Func, nvi-1, lim(0:nlims), nvar+1, (/x1,var/), tol)

    endif

    if (abs(yp(POINTS)%f(1))>MAX_REAL .or. isnan(abs(yp(POINTS)%f(1)) ))  then
!         write(*,*) 'avoid end point singularity'
         x2 = x2-EPSILON(tol)*(x2-x1)
!         write(*,*) 'xp(7) =',xp(7)
        yp(POINTS) = Integrate(Func, nvi-1, lim(0:nlims), nvar+1, (/x2,var/), tol)
!        if(abs(yp(POINTS)%f(1))>MAX_REAL .or.isnan(abs(yp(POINTS)%f(1)) )) then
!            write(*,*) 'Still HAVING SING'
!            read(*,*)
!        endif
    endif

!    write(*,*) 'yp(6) = ',yp(6)%f(1)
    ! Evaluate each interval separately
!    hmin=max(EPSILON(hmin)*abs(xp(7)-xp(1)),1.d-9)
    hmin=EPSILON(hmin)*abs(xp(7)-xp(1))
    q(1)=Quadstep(Func,nvi-1,lim(0:nlims),xp(1),xp(3),yp(1),yp(2),yp(3),tol,hmin,nvar,var )
    q(2)=Quadstep(Func,nvi-1,lim(0:nlims),xp(3),xp(5),yp(3),yp(4),yp(5),tol,hmin,nvar,var )
    q(3)=Quadstep(Func,nvi-1,lim(0:nlims),xp(5),xp(7),yp(5),yp(6),yp(7),tol,hmin,nvar,var )

    ! Final result
    Integrate_R=q(1)+q(2)+q(3)
!    write(*,*) Integrate_R%f(1)
!    write(*,*) Integrate_R%f(2)
!    write(*,*) Integrate_R%f(3)
!    write(*,*) Integrate_R%f(4)
    return
    end function Integrate

    recursive function Quadstep(Func,nvi,lim,a,b,fa,fc,fb,tol,hmin,nvar,var ) result(q)
    ! *******************************************************************************
    ! Applies Simpson's Formula to evaluate the integral of f(x) over [a,b] given
    !   f(a), f(c=(b+a)/2), f(b)
    ! f(x) is the integral of "Func" on "nvi" variables of integration
    ! lim: array holding the lower and upper limits of the variables of integration
    !       # lim(0) should be set to any arbitrary value
    !       # the n-th variable has the lower limit in lim(2n-1),
    !                                   upper limit in lim(2n)
    ! nvar: the number of passed parameters
    ! var : array holding the values of the passed parameters
    ! tol : required tolerance
    ! *******************************************************************************
    integer, intent(in):: nvi
    real(8), dimension(0:nvi), intent(in):: lim
    real(8), intent(in) :: a, b, hmin, tol
    type (multival), intent(in) :: fa, fc, fb
    integer, intent(in) :: nvar
    real(8), intent(in), dimension(nvar) :: var
    interface
         type (multival) function Func(n,x)
            use Operations
            integer, intent(in) :: n
            real(8), intent(in), dimension(n) :: x
        end function Func
    end interface
    type (multival) :: q, q1, q2, qac, qcb, fd, fe
    real(8) :: h, c

    h=b-a
    c=(a+b)/2.d0

    ! If "zero"-width interval, end evaluation
    if (abs(h)<hmin .or. c==a .or. c==b) then
        q=h*fc
        return
    endif

    ! Evaluate f(x) at d=(a+c)/2, e=(c+b)/2
    fd=Integrate(Func,nvi,lim,nvar+1, (/(a+c)*0.5d0,var/), tol )
!    if(isnan(abs(fd%f(1)))) then
!        write(*,*) 'fd is nan'
!    endif
    fe=Integrate(Func,nvi,lim,nvar+1, (/(c+b)*0.5d0,var/), tol )
!    if(isnan(abs(fe%f(1)))) then
!        write(*,*) 'fe is nan'
!    endif

!    if(isnan(abs(fc%f(1)))) then
!        write(*,*) 'fc is nan'
!    endif


    ! Apply Simpson's Formulae
!    q1=(h/6.d0)*(fa+4.d0*fc+fb)
    q1=(h*0.166666666666667d0)*(fa+4.d0*fc+fb)

!    q2=(h/12.d0)*(fa+fb+4.d0*fd+2.d0*fc+4.d0*fe)
    q2=(h*0.083333333333333d0)*(fa+fb+4.d0*fd+2.d0*fc+4.d0*fe)

!    q=(1.d0/15.d0)*(q2-q1)+q2
    q=(0.066666666666667d0)*(q2-q1)+q2

    ! Check for the error
    if (abs((q2%f(1)-q%f(1)))<=TOL) then !!!!!!!!!!! tolerance changed
!        write(*,*) 'error free'
        return
    else
!        write(*,*) 'still error'
!        write(*,*) 'nvi',nvi
!       write(*,*) 'h',h
!       write(*,*) 'fa = ',fa%f(1)
!       write(*,*) 'fd = ',fd%f(1)
!       write(*,*) 'fc = ',fc%f(1)
!       write(*,*) 'fe = ',fe%f(1)
!       write(*,*) 'fb = ',fb%f(1)
        ! Subdivide into two more intervals (a->c and c->b)
        qac=Quadstep(Func,nvi,lim,a,c,fa,fd,fc,tol,hmin,nvar,var )
        qcb=Quadstep(Func,nvi,lim,c,b,fc,fe,fb,tol,hmin,nvar,var )
        q=qac+qcb
    endif
    end function Quadstep


     recursive function Integrate_self (Func, nvi, lim, nvar, var, tol) result(Integrate_R)
    ! *******************************************************************************
    ! Integrates "Func" on "nvi" variables of integration
    ! lim: array holding the lower and upper limits of the variables of integration
    !       # lim(0) should be set to any arbitrary value
    !       # the n-th variable has the lower limit in lim(2n-1),
    !                                   upper limit in lim(2n)
    ! nvar: the number of passed parameters
    ! var : array holding the values of the passed parameters
    ! tol : required tolerance
    ! *******************************************************************************
    interface
        type (multival) function Func(n,x)
            use Operations
            integer, intent(in) :: n
            real(8), intent(in), dimension(n) :: x
        end function Func
    end interface
    integer, intent(in) :: nvi,nvar
    real(8), intent(in), dimension(0:2*nvi) :: lim
    real(8), intent(in), dimension(nvar) :: var
    real(8), intent(in) :: tol
    type (multival):: Integrate_R,y1,y2
    integer, parameter :: POINTS = 7
    real(8) :: h, hmin, x1, x2
    real(8), dimension(POINTS) :: xp
    type (multival), dimension(POINTS):: yp
    type (multival), dimension(3):: q
    integer :: n, nlims

    if (nvi==0) then

        Integrate_R=Func(nvar,(/var/))
        return
    endif

    x2=lim(2*nvi); x1=lim(2*nvi-1)
    nlims=2*(nvi-1)
    h=0.13579d0*(x2-x1)

    ! Seven point evaluation
    xp=(/x1,x1+h,x1+2.d0*h,(x1+x2)/2.d0,x2-2.d0*h,x2-h,x2/)
    do n=1,POINTS
        y1=Integrate(Func, nvi-1, (/lim(0:(nlims-1)) , xp(n)/), nvar+1,(/xp(n),var/), tol )
        y2=Integrate(Func, nvi-1, (/lim(0:(nlims-2)),xp(n),lim(nlims)/), nvar+1,(/xp(n),var/), tol )
        yp(n) = y1+y2
    enddo
!    write(*,*) xp
    ! Avoid end-point singularities
    if (abs(yp(1)%f(1))>MAX_REAL .or. isnan(abs(yp(1)%f(1)) ) ) then
!        write(*,*) 'avoid end point singularity'
!        write(*,*) xp(1)
        x1 = x1+EPSILON(tol)*(x2-x1)
!        write(*,*) xp

        yp(1) = Integrate(Func, nvi-1, lim(0:nlims), nvar+1, (/x1,var/), tol)

    endif

    if (abs(yp(POINTS)%f(1))>MAX_REAL .or. isnan(abs(yp(POINTS)%f(1)) ))  then
!         write(*,*) 'avoid end point singularity'
         x2 = x2-EPSILON(tol)*(x2-x1)
!         write(*,*) 'xp(7) =',xp(7)
        yp(POINTS) = Integrate(Func, nvi-1, lim(0:nlims), nvar+1, (/x2,var/), tol)
!        if(abs(yp(POINTS)%f(1))>MAX_REAL .or.isnan(abs(yp(POINTS)%f(1)) )) then
!            write(*,*) 'Still HAVING SING'
!            read(*,*)
!        endif
    endif

!    write(*,*) 'yp(6) = ',yp(6)%f(1)
    ! Evaluate each interval separately
!    hmin=max(EPSILON(hmin)*abs(xp(7)-xp(1)),1.d-9)
    hmin=EPSILON(hmin)*abs(xp(7)-xp(1))
    q(1)=Quadstep(Func,nvi-1,lim(0:nlims),xp(1),xp(3),yp(1),yp(2),yp(3),tol,hmin,nvar,var )
    q(2)=Quadstep(Func,nvi-1,lim(0:nlims),xp(3),xp(5),yp(3),yp(4),yp(5),tol,hmin,nvar,var )
    q(3)=Quadstep(Func,nvi-1,lim(0:nlims),xp(5),xp(7),yp(5),yp(6),yp(7),tol,hmin,nvar,var )

    ! Final result
    Integrate_R=q(1)+q(2)+q(3)
!    write(*,*) Integrate_R%f(1)
!    write(*,*) Integrate_R%f(2)
!    write(*,*) Integrate_R%f(3)
!    write(*,*) Integrate_R%f(4)
    return
    end function Integrate_self


end module Simpson_Quad
