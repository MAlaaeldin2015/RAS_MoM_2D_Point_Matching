! ============================================================================
! File   : RAS_MoM_2D_Point_Matching.f08
! Purpose: Main solver file. Contains three modules:
!
!   discretizer_3D  -- geometry I/O and 2D scatterer segmentation
!   scatterer_mod   -- the Scatterer derived type plus all RAS and MoM
!                      solver engines, field evaluation, and post-processing
!   RSM_3D_sph      -- top-level driver (run_main_program) that reads input,
!                      initialises scatterers, runs solvers, and writes output
!
! ---- discretizer_3D key routines ----------------------------------------
!   read_parameters          : reads Parameters.dat into sim_par globals
!   read_materials           : reads materials.dat into materials struct
!   segment_superquadratic   : discretises superquadric boundary
!   segment / get_equation   : recursive adaptive boundary discretisation
!   get_intersection         : detects segment crossing for containment tests
!
! ---- scatterer_mod key routines -----------------------------------------
!   solve_scattering              : single-scatterer iterative RAS loop
!   solve_scattering_multiscatterer     : multi-scatterer iterative IFB loop
!   solve_scattering_multiscatterer_once: multi-scatterer merged single solve
!   solve_scattering_multiscatterer_2   : multi-scatterer full IFB iteration
!   eval_matrices_RAS            : builds RAS Gram matrix and excitation
!   get_LSM_matrices             : builds MoM system matrix (packed Hermitian)
!   get_LSM_excitation           : builds MoM right-hand side vector
!   eval_kernel_Electric_2D      : 2D electric line-source field kernel
!   eval_kernel_Magnetic_2D      : 2D magnetic line-source field kernel
!   eval_near_field_Electric/Magnetic_2D: scattered near-field at a point
!   eval_bouncing_field          : computes field from one scatterer on another
!   set_excitation               : projects incident field onto testing points
!   set_Et                       : assembles total excitation (incident + bounce)
!   eval_error                   : normalised boundary residual
!   initialize_parameters        : allocates all Scatterer arrays
!   initialize_stage             : resets solver state for a new excitation
!   add_sources_rect/contour/outside_random: source placement routines
!   discretize_scatterer_file/superquad/rect: testing-point placement
!   eval_surface_current_MoM_once: full MoM matrix solve
!   eval_current_MoM             : evaluates MoM surface currents
!   set_MoM_BC_matrices          : assembles MoM impedance/admittance blocks
!   set_MoM_AWE_matrices         : AWE Taylor expansion of MoM matrices
!   eval_Pade_coefficients       : computes Pade [L/M] coefficients from Taylor
!   eval_MoM_Pade_coefficients   : MoM-specific Pade expansion
!   Eval_chebyshev_expansion_coeff_generic: Chebyshev fit of Taylor series
!   eval_far_field_electric/magnetic: far-field contributions per source
!   compare_far_field_scatterers : writes bistatic RCS comparison file
!   Export_RAS_AWE/MoM_AWE_Chebyshev_Monostatic_RCS: monostatic RCS sweep
!   get_solution_bandwidth       : adaptive AWE bandwidth estimator
!   post_processing_colormap_imaging: near-field colour map
!
! Author : Mohamed A. Moharram Hassan
! License: Apache 2.0 -- see LICENSE
! ============================================================================


module discretizer_3D
    use sim_par
    use Operations
    use constants
    implicit none

contains



    ! -----------------------------------------------------------------------
    ! Subroutine: read_materials
    ! Purpose   : Reads 'materials.dat' and populates the global 'materials'
    !             struct (type material_lib) with dielectric and IBC entries.
    !             Format: see materials.dat header comments.
    ! -----------------------------------------------------------------------
    subroutine read_materials()
        integer :: fd,ii
        real(8):: zz_r,zz_i,zt_r,zt_i,tz_r,tz_i,tt_r,tt_i


        fd=49
        open(fd, FILE="materials.dat", STATUS = "OLD")
        read(fd,*) materials%N_diel,materials%N_ibc
        allocate(materials%eta_zz(materials%N_ibc),materials%eta_zt(materials%N_ibc),&
        materials%eta_tz(materials%N_ibc),materials%eta_tt(materials%N_ibc),materials%ID_ibc(materials%N_ibc))
        allocate(materials%ID_diel(materials%N_diel),materials%eps_r_p(materials%N_diel),&
        materials%eps_r_pp(materials%N_diel),materials%mu_r_pp(materials%N_diel),materials%mu_r_p(materials%N_diel))
        read(fd,*)
        read(fd,*)
        do ii = 1,materials%N_diel
            read(fd,*) materials%ID_diel(ii),materials%eps_r_p(ii),materials%eps_r_pp(ii),materials%mu_r_p(ii),materials%mu_r_pp(ii)
        !            write(*,*) materials%eps_r_p(ii)
        enddo
        read(fd,*)
        read(fd,*)
        do ii = 1,materials%N_ibc
            read(fd,*) materials%ID_ibc(ii),zz_r,zz_i,zt_r,zt_i,tz_r,tz_i,tt_r,tt_i
            materials%eta_zz(ii) = zz_r + cj*zz_i
            materials%eta_zt(ii) = zt_r + cj*zt_i
            materials%eta_tz(ii) = tz_r + cj*tz_i
            materials%eta_tt(ii) = tt_r + cj*tt_i
        !            write(*,*) materials%eta_zz(ii),materials%eta_zt(ii),materials%eta_tz(ii),materials%eta_tt(ii)
        enddo
        close(fd)

    end subroutine read_materials

    ! -----------------------------------------------------------------------
    ! Subroutine: read_Line_Source_Model
    ! Purpose   : Reads 'Line_Sources_Model.dat' when Source_Model=1.
    !             Populates Line_sources(:) array (type, orientation,
    !             amplitude, position).
    ! -----------------------------------------------------------------------
    subroutine read_Line_Source_Model()
        integer :: fd,ii
        real(8) :: A_r,A_i
        fd=49
        open(fd, FILE="Line_Sources_Model.dat", STATUS = "OLD")
        read(fd,*) N_Line_sources

        allocate(Line_sources(N_line_sources))
        do ii = 1,N_Line_sources
            read(fd,*) Line_sources(ii)%Source_Type,Line_sources(ii)%Source_orientation,A_r,A_i,&
            Line_sources(ii)%x_s,Line_sources(ii)%y_s
            Line_sources(ii)%Amp = A_r+cj*A_i
!            write(*,*)Line_sources(ii)%Source_Type,Line_sources(ii)%Amp,Line_sources(ii)%x_s,Line_sources(ii)%y_s
            Line_Sources(ii)%is_normalized = .false.
        enddo
        close(fd)
!        stop
    end subroutine read_Line_Source_Model

    ! -----------------------------------------------------------------------
    ! Subroutine: read_parameters
    ! Purpose   : Master input reader. Reads Parameters.dat line-by-line
    !             into sim_par globals, allocates sweep arrays, and computes
    !             derived wave quantities (k0, k1, az, ar, lambda, theta_i).
    !             Then calls read_materials() and (if needed)
    !             read_Line_Source_Model().
    ! -----------------------------------------------------------------------
    subroutine read_parameters()
        !! This routine read the parameters file "Parameters.dat"; it should be in this format
        !! R_sph
        !! maximum Allowed Number of Iterations (max_iteration), N_group
        !! bound
        !! theta_i in degrees
        !! phi_i in degrees
        !! alpha_i in degrees
        !! Samples_per_Wavelength  !! the number of desired samples per wavelength
        !! er_1 !! relative permittivity of the the outside medium
        !! mur_1 !! relative permeability of the the outside medium
        !! frequency !! the operating frequency of the problem

        integer :: fd
        integer :: theta_i_degree,phi_i_degree,alpha_i_degree
        integer :: i,j


        fd=50
        open(fd, FILE="Parameters.dat", STATUS = "OLD")
        read(fd,*) Source_Model,theta_i_degree ,phi_i_degree,alpha_i_degree
        read(fd,*) number_of_scatterers,RAS_solution_method,MoM_solution_method


        allocate(physical_parameters(number_of_scatterers,5))
        allocate(scatterers_input_method(number_of_scatterers))
        allocate(scatterer_input_file_names(number_of_scatterers),material_id_read(number_of_scatterers))

        read(fd,*) (scatterers_input_method(i),i=1,number_of_scatterers)
        !        read(fd,*) superquad_center(1),superquad_center(2), a_superquad,b_superquad,g_superquad
        do i = 1,number_of_scatterers
            if(scatterers_input_method(i) == 1) then
                read(fd,*) material_id_read(i), (physical_parameters(i,j),j=1,5)
            elseif(scatterers_input_method(i) == 2) then
                read(fd,*) material_id_read(i), scatterer_input_file_names(i)
            !                write(*,*) scatterer_input_file_names(i)
            !                stop
            else
                write(*,*) 'ERROR: Bad argument in scatterers input method: check your input values'
                stop
            endif
        enddo

        !        write(*,*) 'a   ,b   ,g               = ',a_superquad,b_superquad,g_superquad
        write(*,*) 'theta_i,phi_i,alpha_i       = ',theta_i_degree,phi_i_degree,alpha_i_degree
        !        write(*,*) 'Materials Parameters      = ', er_1,mur_1,er_2,mur_2
        read(fd,*) frequency
        write(*,*) 'Frequency                   = ',frequency

        read(fd,*) simulation_mode,plot_current,samples_N,samples_S,N_iterations
        !        write(*,*) 'Simulation Mode =',simulation_mode,'Plot Current =',plot_current
        if(simulation_mode == 0 .or. simulation_mode == 5) then
            allocate(N_group_read(1),SPW_read(1))
            samples_N = 1
            samples_S = 1
            if(simulation_mode == 0) then
                N_iterations = 1
            endif
        else

            if(simulation_mode == 1) then
                allocate(N_group_read(1),SPW_read(1))
                samples_N = 1
                samples_S = 1

            elseif(simulation_mode == 2) then
                allocate(N_group_read(samples_N),SPW_read(1))

                samples_S = 1
            elseif(simulation_mode == 3) then
                allocate(SPW_read(samples_S),N_group_read(1))

                samples_N = 1
            else
                allocate(N_group_read(samples_N),SPW_read(samples_S))
            endif

        endif

        read(fd,*) (N_group_read(i),i=1,samples_N)
        read(fd,*) max_iteration,max_bouncing_iteration
        !        write(*,*) 'Max. Number of Iterations =',max_iteration
        write(*,*) 'Group Number of Sources     =',(N_group_read(i),i=1,samples_N)
        read(fd,*) (SPW_read(i),i=1,samples_S)
        read(fd,*) MoM_Samples_per_wavelength
        read(fd,*) non_uniform_type_read,N_o_read,A_segmentation_read
        write(*,*) 'Samples per Wavelength      =',(SPW_read(i),i=1,samples_S)
        write(*,*) 'Samples per Wavelength (MoM)=',MoM_Samples_per_wavelength
        read(fd,*) bound,outside_bound,corner_sep,corner_sep2
        read(fd,*) samples_per_wavelength_contour,samples_per_wavelength_outside
        read(fd,*) TOL, TOL_segmentation
        read(fd,*) MoM_activation_flag,R_probes_in,R_calc_in
        read(fd,*) FormulationType
        read(fd,*) Source_Placement
        read(fd,*) Contour_Sources_Type_read
        read(fd,*) Matrix_Solution_Method,MAX_ALLOWED_CONDITION_NUMBER,Y_mat_norm_limit
        read(fd,*) Wideband_type,N_taylor,Pade_L,fl_r,fh_r,n_points_freq,BW_limit
        write(*,*)
        close(fd)


        call read_materials()
        if(Source_Model == 1) then
            call read_Line_Source_Model()
        endif



        k0 = tpi*frequency/3d8



        k1 = k0
        ak = 1.d0

        lambda = tpi/k1
        theta_i = theta_i_degree/dpr
        phi_i = phi_i_degree/dpr
        alpha_i = alpha_i_degree/dpr


        eta1 = eta0

        az = ak*cos(theta_i)
        ar = sqrt(1.d0 - az**2.d0)

        if(MoM_activation_flag >= 0 .and. MoM_activation_flag <= 3) then
            Wideband_type = 0
            N_taylor = 0 !! discard whatever input from the file
        endif

        if(Wideband_type == 0) then !! a wide band solution not is required

            N_taylor = 0 !! discard whatever input from the file

        endif



    !        write(*,*) ar,sqrt(k1**2.d0-(az_l*k0)**2.d0)/k0

    end subroutine read_parameters

    subroutine deallocate_materials()
        if(allocated(materials%ID_ibc)) then
            deallocate(materials%eta_zz,materials%eta_zt,materials%eta_tz,materials%eta_tt,materials%ID_ibc)
        endif
        if(allocated(materials%ID_diel)) then
            deallocate(materials%ID_diel,materials%eps_r_p,materials%mu_r_p,materials%eps_r_pp,materials%mu_r_pp)
        endif
    end subroutine deallocate_materials

    ! -----------------------------------------------------------------------
    ! Subroutine: segment_superquadratic
    ! Purpose   : Generates a set of 2D boundary points describing the
    !             superquadric curve  (x/a)^g + (y/b)^g = 1 using the
    !             recursive 'segment' subroutine with chord-error tolerance
    !             TOL_s * lambda. Points are ordered counter-clockwise.
    !             g=2 gives an ellipse; larger g approaches a rectangle.
    ! Inputs : a_s, b_s (semi-axes), g_s (shape), TOL_s (tolerance/lambda)
    ! Output : points(:,2) -- boundary coordinate array (allocated on exit)
    ! -----------------------------------------------------------------------
    subroutine segment_superquadratic(a_s,b_s,g_s,points,TOL_s)
        real(8),allocatable,dimension(:,:) :: points_cylinder_quarter,points1,points2,points_temp
        real(8),allocatable,intent(inout) :: points(:,:)
        real(8) :: a_s,b_s,g_s,TOL_s
        integer :: Sp,ii,S1,S2
        real(8) :: r_mid,x_mid,y_mid,tolerence



        r_mid =sqrt(a_s**2.d0+b_s**2.d0)/(2.d0)**(1.d0/g_s)
        x_mid = r_mid*a_s/sqrt(a_s**2.d0+b_s**2.d0)
        y_mid = r_mid*b_s/sqrt(a_s**2.d0+b_s**2.d0)


        tolerence = TOL_s*lambda


        call segment(0.d0,x_mid,b_s,y_mid,a_s,b_s,g_s,tolerence,points1)
        call segment(y_mid,0.d0,x_mid,a_s,b_s,a_s,g_s,tolerence,points2)

        S1 = size(points1,1)
        S2 = size(points2,1)

        allocate(points_cylinder_quarter(S1+S2-1,2))

        points_cylinder_quarter(1:S1,:) = points1
        points_cylinder_quarter((S1+1):(S1+S2-1),1) = points2(2:S2,2)
        points_cylinder_quarter((S1+1):(S1+S2-1),2) = points2(2:S2,1)

        Sp = size(points_cylinder_quarter,1)

        allocate(points(4*Sp-4,2))
        points = 0.d0

        points(1:Sp,:) = points_cylinder_quarter
        do ii=1,Sp
            if((ii > 1) .and. (ii < Sp)) then
                points(Sp+ii-1,:) = (/points_cylinder_quarter(Sp-ii+1,1),-points_cylinder_quarter(Sp-ii+1,2)/)
                !                write(*,*) Sp+ii-1,points(Sp+ii-1,:)
                points(3*Sp-2+ii-1,:) = (/-points_cylinder_quarter(Sp-ii+1,1),points_cylinder_quarter(Sp-ii+1,2)/)
            endif

            points(2*Sp-2+ii,:) = (/-points_cylinder_quarter(ii,1),-points_cylinder_quarter(ii,2)/)
        !            write(*,*) 2*Sp-1+ii,points(2*Sp-1+ii,:)
        enddo
        deallocate(points_cylinder_quarter,points1,points2)
        !! invert points to be anti-clockwise rotation
        allocate(points_temp(size(points,1),2))
        points_temp = points
        Sp = size(points,1)
        points(1,:) = points_temp(1,:)
        do ii = 1,(Sp-1)
            points(ii+1,:) = points_temp(Sp-ii+1,:)

        enddo
        deallocate(points_temp)
        !! writing points to file
    !        OPEN(17, FILE='superquad_points.dat')
    !        write(17,*) a_s,b_s,g_s
    !        do ii = 1,Sp
    !            write(17,*) (/points(ii,:), 0.d0/)
    !        enddo
    !        close(17)

    !        write(*,*) size(points,1),'points :)'
    !        do ii=1,size(points,1)
    !            write(*,*) points(ii,:)
    !        enddo
    end subroutine segment_superquadratic

    recursive subroutine segment(x1,x2,y1,y2,a,b,g,tolerence,points)
        !! a function to generate points to descibe the superquadratic curve defined by a,b, and g. It uses recursive iteratios
        real(8),intent(in) :: x1,x2,y2,y1,a,b,g,tolerence
        real(8),allocatable,dimension(:,:) :: points1,points3
        real(8),allocatable,intent(inout)::points(:,:)
        real(8) :: x,y,y_mid_curve
        integer :: S1,S3

        x = (x1+x2)/2.d0
        y = (y1+y2)/2.d0

        y_mid_curve = b *(1.d0-(x/a)**g)**(1.d0/g)

        if(abs(y-y_mid_curve) < tolerence) then
            allocate(points(2,2))
            points(1,:) = (/x1,y1/)
            points(2,:) = (/x2,y2/)

        else
            call segment(x1,x,y1,y_mid_curve,a,b,g,tolerence,points1)
            call segment(x,x2,y_mid_curve,y2,a,b,g,tolerence,points3)


            S1 = size(points1,1)
            S3 = size(points3,1)

            allocate(points(S1+S3-1,2))

            points(1:S1,:) = points1
            points((S1+1):(S1+S3-1),:) = points3(2:S3,:)


            deallocate(points1,points3)
        endif

    end subroutine segment

    subroutine get_equation(x,y,xb,xe,yb,ye,a,b,c,d)
        real(8),allocatable,intent(inout) :: x(:),y(:)
        real(8),allocatable,intent(inout) :: xb(:),xe(:),yb(:),ye(:),a(:),b(:),c(:),d(:)
        integer :: N

        N = size(x,1)
        allocate(xb(N),xe(N),yb(N),ye(N),a(N),b(N),c(N),d(N))
        xb = x
        xe(1:(N-1)) = x(2:N)
        xe(N) = x(1)

        yb = y
        ye(1:(N-1)) = y(2:N)
        ye(N) = y(1)

        a = (ye-yb)/(xe-xb)
        b = yb - xb*a

        c = (xe-xb)/(ye-yb)
        d = xb - yb*c
    end subroutine get_equation

    subroutine get_intersection(xb_s,xe_s,yb_s,ye_s,a_s,b_s,c_s,d_s,&
    xb_q,xe_q,yb_q,ye_q,a_q,b_q,c_q,d_q,int_pt)
        real(8),allocatable,intent(inout) ::int_pt(:,:)
        real(8),allocatable,intent(inout) :: xb_s(:),xe_s(:),yb_s(:),ye_s(:),a_s(:),b_s(:),c_s(:),d_s(:)
        real(8),allocatable,intent(inout) :: xb_q(:),xe_q(:),yb_q(:),ye_q(:),a_q(:),b_q(:),c_q(:),d_q(:)
        integer :: i,j,N_s,N_q,N_allocated
        real(8),allocatable,dimension(:,:) :: int_pt_temp
        real(8) :: x_min,x_max,y_max,y_min,x_x,y_x
        real(8) :: x2_min,x2_max,y2_max,y2_min
        real(8),parameter :: MAX_LIMIT=100

        N_s = size(xb_s,1)
        N_q = size(xb_q,1)
        do i = 1,N_s
            x_min = xb_s(i)
            if(x_min > xe_s(i)) then
                x_min = xe_s(i)
                x_max = xb_s(i)
            else
                x_max = xe_s(i)
            endif
            y_min = yb_s(i)
            if(y_min > ye_s(i)) then
                y_min = ye_s(i)
                y_max = yb_s(i)
            else
                y_max = ye_s(i)
            endif

            do j =1,N_q
                x2_min = xb_q(j)
                if(x2_min > xe_q(j)) then
                    x2_min = xe_q(j)
                    x2_max = xb_q(j)
                else
                    x2_max = xe_q(j)
                endif
                y2_min = yb_q(j)
                if(y2_min > ye_q(j)) then
                    y2_min = ye_q(j)
                    y2_max = yb_q(j)
                else
                    y2_max = ye_q(j)
                endif
                x_x = -(b_s(i)-b_q(j))/(a_s(i)-a_q(j))
                if(isnan(x_x) .or. x_x > MAX_LIMIT) then
                    !!! check the y
                    y_x = -(d_s(i)-d_q(j))/(c_s(i)-c_q(j))
                    if((y_x>=y_min) .and. (y_x <=y_max)) then
                        if((y_x>=y2_min) .and. (y_x <=y2_max)) then
                            x_x = c_s(i)*y_x+d_s(i)
                            !! add point sequence
                            if(.not. allocated(int_pt)) then
                                allocate(int_pt(1,4))
                                int_pt(1,:) = (/x_x,y_x,dble(i),dble(j)/)
                            else
                                N_allocated = size(int_pt,1)
                                allocate(int_pt_temp(N_allocated,4))
                                int_pt_temp = int_pt
                                deallocate(int_pt)
                                allocate(int_pt(N_allocated+1,4))
                                int_pt(1:N_allocated,:) = int_pt_temp
                                int_pt(1+N_allocated,:) = (/x_x,y_x,dble(i),dble(j)/)

                            endif

                        endif
                    endif
                elseif((x_x>=x_min) .and. (x_x <=x_max)) then
                    if((x_x>=x2_min) .and. (x_x <=x2_max)) then
                        y_x = a_s(i)*x_x+b_s(i)
                        !! add point sequence
                        if(.not. allocated(int_pt)) then
                            allocate(int_pt(1,2))
                            int_pt(1,:) = (/x_x,y_x/)
                        else
                            N_allocated = size(int_pt,1)
                            allocate(int_pt_temp(N_allocated,2))
                            int_pt_temp = int_pt
                            deallocate(int_pt)
                            allocate(int_pt(N_allocated+1,2))
                            int_pt(1:N_allocated,:) = int_pt_temp
                            int_pt(1+N_allocated,:) = (/x_x,y_x/)

                        endif
                    endif
                endif

            enddo

        enddo
        if(allocated(int_pt_temp)) then
            deallocate(int_pt_temp)
        endif
    end subroutine get_intersection

end module discretizer_3D


module scatterer_mod
    use constants
    use sim_par
    use Operations
    use discretizer_3D
    use Gauss_Reduction
    use Simpson_Quad
    use io_matrix
    use Matrices_Storage_Handling

    implicit none

    type, public :: Scatterer
        !! material parameters
        !        real(8) :: k_local,eta_local,lambda_local
        real(8) :: lambda_local
        complex*16 ::k_local,eta_local
        complex*16 :: ar_local !! the constat to compute kr_local = ar_local*k0
        complex*16 :: ak_local
        !        real(8) :: ar_local !! the constat to compute kr_local = ar_local*k0
        !        real(8) :: ak_local
        integer :: region_ID !! a number to indicate which region is this scatterer
        integer :: Problem_type !! indicate which Boundary condition type is this scatterer
        !! 1- PEC
        !! 2- PMC
        !! 3- IBC
        !! 4- Dielectric

        !! physical parameters
        real(8),dimension(3) :: center !! the center of the spherical scatterer
        real(8) :: radius !! the radius of the scatterer
        real(8) :: a_superquad,b_superquad,g_superquad
        real(8),allocatable,dimension(:,:) :: segments_points !! segments of the scatterer
        !! Solution Parameters
        integer :: N_curr !! the current number of sources introduced inside the scatterer
        integer :: N_group !! Added number of sources (Group size)
        integer :: N_max !! maximum allowed number of line sources
        integer :: M  ! number of testing points
        integer :: N_New_Iterations !! the allowed number of new Iterations
        integer :: N_sub_groups !! number of subgroups of sources that are intorduced one by one as new sources positions
        integer :: group_index !! the current group index from 1 -> N_sub_groups
        real(8),allocatable, dimension(:) :: delta_n !! the separation, on the primeter of the 2D scatterer, between testing points
        integer :: itr_counter = 1
        !! Excitation Parameters
        type(vector_c), allocatable,dimension(:,:) :: Ei,Hi !! incident electric and magnetic fields at the surface of the condutor at the testing points
        complex*16,allocatable,dimension(:,:) :: E,H !! excitation vectors
        complex*16,allocatable,dimension(:,:) :: E_current,H_current !! excitation vectors
        complex*16,allocatable,dimension(:,:):: Eu,Ev,Hu,Hv

        type(Matrix_Storage_Linked_List),pointer :: Matrix_Storage
        !! Bouncing Fields (used for multiscatterer problems). They are added to the incident fields
        complex*16,allocatable,dimension(:,:) :: E_bounce_current,H_bounce_current
        complex*16,allocatable,dimension(:,:) :: E_bounce,H_bounce

        real(8),allocatable,dimension(:,:) :: source_pos !! the introduced sources positions
        integer,allocatable,dimension(:) :: allowed_orientation !! possible orientations 0-all, 1- (Z) only
        integer,allocatable,dimension(:) :: active_region  !! indicates which region to contribute to
        type(vector_c),allocatable,dimension(:,:) :: I_sources !!currents on the sources (Electric or magnetic)
        integer,allocatable,dimension(:) :: I_stat !! a flag array to determine which currents are for Electric and which are for Magnetic
        !! 1- Electric sources, 2- Magnetic sources, 0- Undetermined yet

        real(8),allocatable,dimension(:,:) :: contour_points !! contour points allocated to compensate for flat surfaces
        integer,allocatable,dimension(:) :: contour_points_type,contour_points_orientation,contour_point_active_region
        integer :: N_con !! the number of contour points
        integer :: Contour_Cnt !! a counter for the used contour points

        real(8),allocatable,dimension(:,:) :: Inside_sources_dielectric !! Contour points allocated for the purpose of representing the dielctric fields inside the scatterer
        integer :: N_inside_sources !! the number of sources for the inside region problem
        integer :: Contour_out_Cnt !! a counter for the used contour points used for the outer region solution

        !! the total currents
        type(vector_c),allocatable,dimension(:,:) :: I_sources_total_current !!currents on the sources (Electric or magnetic)
        !! collective LU decomposition matrices (vectors) and their lengths
        type(vector_c),allocatable,dimension(:,:) :: I_sources_total
        integer, allocatable, dimension(:) :: Nm_collect,n_sources_array,n_sources_outside_array
        logical,allocatable,dimension(:) :: is_sub_group_used
        complex*16, allocatable,dimension(:) :: Lm_collect,Dm_collect
        complex*16,allocatable,dimension(:,:) :: Y_collect
        integer :: N !! the total used number of sources


        !! Testing Points
        real(8),allocatable,dimension(:,:) :: testing_pt ! testing points on the surface of the conductor (x,y)
        integer,allocatable,dimension(:,:) :: testing_pt_status
        integer,dimension(3) :: region_status !! this to distinguish between different situations of the scatterer
        !! 0 X X -> is a stand alone scatterer
        !! 1 i X -> is a contained scatterer by region i
        type(vector),allocatable,dimension(:) :: testing_pt_MoM
        !! vectors depending on the geometry (set by disctritiziation routines)
        type(vector),allocatable,dimension(:) :: norm_v,tang_u,tang_v !! normal and tangential vectors u,v to the conductor evaluated at the testing points

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !! Evaluated Field on the surface points
        complex*16,allocatable,dimension(:,:) :: Es_current,EH_excit
        complex*16,allocatable,dimension(:,:) :: Hs_current !! the scattered magnetic field at the surface of the conductor
        complex*16,allocatable,dimension(:,:) :: Hs,Es,Hd,Ed
        complex*16,allocatable,dimension(:,:) :: Hd_current,Ed_current !! The fields inside the dielectric material (evaluated only when Problem Type = Dielectric only)

        !! to monitor the divergence of the solution of the problem, the error should be calculated first then an update for the boundary conditions is to be made
        !! So, a new variables for the current Es, Hs, Hd, Ed are to be introduced upon evaluating the scattered and dielectric fields. They are added to the final
        !! Es, ... when the error is proved to be convergent after passing the limit of 1%
        complex*16,allocatable,dimension(:,:) :: Hd_c,Ed_c,Es_c,Hs_c

        real(8),allocatable,dimension(:) :: normalize_E,normalize_H !! normalization factor for the error criterion sum(abs(Ei_vector)^2),sum(abs(Hi_vector)^2)
        logical :: converging
        complex*16,allocatable,dimension(:) :: Lm_x,Dm_x

        !! source addition rules for rect style segmentation
        real(8),allocatable,dimension(:,:) :: Str_Rules !! adding points for sources rules (based on structure)
        integer :: N_rules !! number of rules

        integer :: Contour_Sources_Type !! Type of the contour source
        !! 1- Electric, 2- Magnetic, 3- Hybrid (Randomly Chosen), 4- Auto (Determined by the program for best convergence).

        !! Impedance Boundary Conditions matrices
        complex*16,allocatable,dimension(:) :: eta_vv,eta_uu,eta_uv,eta_vu
        real(8),allocatable,dimension(:) :: norm_eta_v,norm_eta_u

        !! MoM matrix inverted
        complex*16,allocatable,dimension(:,:) :: MoM_matrix_inverted
        complex*16,allocatable,dimension(:,:,:) :: MoM_ZY
        complex*16,allocatable,dimension(:,:):: Eu_MoM,Ev_MoM,Hu_MoM,Hv_MoM
        complex*16,allocatable,dimension(:):: combined_bouncing_fields
        !! MoM surface currents
        complex*16,allocatable,dimension(:,:) :: I_u,I_v,M_u,M_v
        complex*16,allocatable,dimension(:,:) :: I_u_bounce,I_v_bounce,M_u_bounce,M_v_bounce

        type(vector_c),allocatable,dimension(:,:) :: Pade_a,Pade_b,Chebyshev_c,Chebyshev_d
        complex*16,allocatable,dimension(:,:) :: MoM_Pade_a_I_u,MoM_Pade_b_I_u
        complex*16,allocatable,dimension(:,:) :: MoM_Pade_a_I_v,MoM_Pade_b_I_v
        complex*16,allocatable,dimension(:,:) :: MoM_Pade_a_M_u,MoM_Pade_b_M_u
        complex*16,allocatable,dimension(:,:) :: MoM_Pade_a_M_v,MoM_Pade_b_M_v
        complex*16,allocatable,dimension(:,:) :: MoM_Chebyshev_d_I_u,MoM_Chebyshev_c_I_u
        complex*16,allocatable,dimension(:,:) :: MoM_Chebyshev_d_M_u,MoM_Chebyshev_c_M_u
        complex*16,allocatable,dimension(:,:) :: MoM_Chebyshev_d_I_v,MoM_Chebyshev_c_I_v
        complex*16,allocatable,dimension(:,:) :: MoM_Chebyshev_d_M_v,MoM_Chebyshev_c_M_v
        integer :: Pade_L,Pade_M

        !! parameters required to add random points outside the scatterer
        real(8) :: y_bound_max,y_bound_min,x_bound_max,x_bound_min

        complex*16 :: eta_zz,eta_zt,eta_tz,eta_tt

    end type Scatterer

    interface eval_near_field_Magnetic_2D
        module procedure eval_near_field_Magnetic_2D_c,eval_near_field_Magnetic_2D_r
    end interface eval_near_field_Magnetic_2D

    interface eval_near_field_Electric_2D
        module procedure eval_near_field_Electric_2D_c,eval_near_field_Electric_2D_r
    end interface eval_near_field_Electric_2D


    interface set_MoM_AWE_matrices
        module procedure set_MoM_AWE_matrices_r,set_MoM_AWE_matrices_c
    end interface set_MoM_AWE_matrices

    interface set_MoM_matrices
        module procedure set_MoM_matrices_r,set_MoM_matrices_c
    end interface set_MoM_matrices

contains


    subroutine deallocate_matrices_storage_linked_list(this)
        class(Scatterer) :: this
        integer :: f_counter

        call this%Matrix_Storage%destroy()
    end subroutine deallocate_matrices_storage_linked_list

    subroutine remove_redundant_testing_points(these)
        type(scatterer),allocatable,intent(inout) :: these(:)
        integer :: ii,jj,mii,mjj

        do ii=1,number_of_scatterers
            do jj = ii+1,number_of_scatterers
                do mii = 1,these(ii)%M
                    do mjj = 1,these(jj)%M
                        if((these(ii)%testing_pt(mii,1) == these(jj)%testing_pt(mjj,1)) .and. &
                        (these(ii)%testing_pt(mii,2) == these(jj)%testing_pt(mjj,2)) ) then
                            these(jj)%testing_pt_status(mjj,1) = -1 !! this indicates a cancelled point
                            these(ii)%testing_pt_status(mii,3) = jj
                            write(*,*) 'redundant',jj,mjj
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine remove_redundant_testing_points

    subroutine segment_scatterers(these)
        type(scatterer),allocatable,intent(inout) :: these(:)
        integer:: ss,ii,jj,kk
        real(8) :: value_max,xi_min,xi_max,xj_min,xj_max,x_s,y_s,y_max,y_min

        do ii = 1,number_of_scatterers
            if(scatterers_input_method(ii) == 1) then
                these(ii)%radius = these(ii)%a_superquad
                call segment_superquadratic(these(ii)%a_superquad,these(ii)%b_superquad,these(ii)%g_superquad,&
                these(ii)%segments_points,TOL_segmentation)


                do ss = 1,size(these(ii)%segments_points,1)

                    these(ii)%segments_points(ss,:) = these(ii)%segments_points(ss,:) + these(ii)%Center(1:2)
                enddo
                if(these(ii)%a_superquad > these(ii)%b_superquad) then
                    value_max = these(ii)%a_superquad
                else
                    value_max = these(ii)%b_superquad
                endif
                these(ii)%x_bound_max = these(ii)%Center(1) + 3.d0*value_max
                these(ii)%x_bound_min = these(ii)%Center(1) - 3.d0*value_max

                these(ii)%y_bound_max = these(ii)%Center(2) + 3.d0*value_max
                these(ii)%y_bound_min = these(ii)%Center(2) - 3.d0*value_max

            elseif(scatterers_input_method(ii) == 2) then
                call read_points_file(these(ii),these(ii)%segments_points)
            else
                write(*,*) 'BAD input for the scatterers_input_method ',ii,', it should be either 1 or 2 '
                stop
            endif
            call  define_rules(these(ii),these(ii)%segments_points)
        enddo
        !!1- CHECK INTERSECTION and remove common regions and do the same as touching regions
        !!2- CHECK TOUCHING and indicate outside region for each segment
        !!3- CHECK CONTAINMENT and define the exterior region for each scatterer

        !! check containment

        do ii = 1,number_of_scatterers
            these(ii)%region_status = (/0,0,ii/)
            do jj = ii+1,number_of_scatterers


                xi_max = these(ii)%Str_Rules(these(ii)%N_rules,2)
                xi_min = these(ii)%Str_Rules(1,1)
                xj_max = these(jj)%Str_Rules(these(jj)%N_rules,2)
                xj_min = these(jj)%Str_Rules(1,1)
                if((xi_min <= xj_min) .or. (xj_max <= xi_max)) then
                    !! this cannot be if they are coinside
                    cycle
                endif
                do ss = 1,size(these(ii)%segments_points,1)
                    x_s = these(ii)%segments_points(ss,1)
                    y_s = these(ii)%segments_points(ss,2)
                    do kk = 1,these(jj)%N_rules
                        if(x_s < these(jj)%Str_Rules(kk,2)) then
                            !                        i_rule = kk
                            exit
                        endif
                    enddo
                    y_max = these(jj)%Str_Rules(kk,3)*x_s+these(jj)%Str_Rules(kk,4)
                    y_min = these(jj)%Str_Rules(kk,5)*x_s+these(jj)%Str_Rules(kk,6)
                    if((y_s <= y_min) .or. (y_max <= y_s)) then
                        exit
                    endif
                enddo
                !            write(*,*) 'inside containment check',ss,size(these(ii)%segments_points,1)
                if((ss-1) == size(these(ii)%segments_points,1)) then
                    !! ii is contained in jj
                    if(these(ii)%region_status(1) == 0) then
                        these(ii)%region_status = (/1,jj,ii/)
                    !                    write(*,*) these(ii)%region_status
                    elseif(these(ii)%region_status(1) == 1) then
                        write(*,*) 'WARNING: make sure that coinside regions are defined from the inside out!'
                    endif
                endif


            enddo
        enddo

    end subroutine segment_scatterers

    ! -----------------------------------------------------------------------
    ! Subroutine: set_IBC_impedance_matrices
    ! Purpose   : Initialises the four components of the 2x2 IBC surface
    !             impedance tensor (eta_vv, eta_uu, eta_uv, eta_vu) at each
    !             testing point from the material-level tensor entries
    !             (eta_zz, eta_zt, eta_tz, eta_tt) stored in the Scatterer.
    !
    ! IBC formulation reference:
    !   A. A. Kishk and P.-S. Kildal, "Electromagnetic scattering from two
    !   dimensional anisotropic impedance objects under oblique plane wave
    !   incidence," Applied Computational Electromagnetics Society Journal,
    !   vol. 10, no. 3, pp. 81-92, 1995.
    !   Note: typographical errors in the original paper have been corrected
    !   in this implementation.
    ! -----------------------------------------------------------------------
    subroutine set_IBC_impedance_matrices(this)
        type(Scatterer) :: this
        !        integer :: ii

        allocate(this%eta_vv(this%M),this%eta_uu(this%M),this%eta_uv(this%M),this%eta_vu(this%M))

        this%eta_vv = this%eta_zz
        this%eta_vu = this%eta_zt
        this%eta_uv = this%eta_tz
        this%eta_uu = this%eta_tt

        allocate(this%norm_eta_v(this%M),this%norm_eta_u(this%M))
        this%norm_eta_v = abs(this%eta_vv)**2.d0 + abs(this%eta_vu)**2.d0
        this%norm_eta_u = abs(this%eta_uv)**2.d0 + abs(this%eta_uu)**2.d0

        !        do ii = 10,30
        !            write(*,*) this%norm_eta_v(ii),this%norm_eta_u(ii)
        !        enddo

        if(this%norm_eta_v(1) < 0.5d0) then
            this%norm_eta_v = 1.d0
        else
            this%norm_eta_v = 1.d0/sqrt(this%norm_eta_v)
        endif
        if(this%norm_eta_u(1) < 0.5d0) then
            this%norm_eta_u = 1.d0
        else
            this%norm_eta_u = 1.d0/sqrt(this%norm_eta_u)
        endif


    !        this%norm_eta_u = 1.d0
    !        this%norm_eta_v = 1.d0
    !        write(*,*) 'norm_eta_u =',sum(this%norm_eta_u)/this%M,'norm_eta_v =',sum(this%norm_eta_v)/this%M
    !        stop
    end subroutine set_IBC_impedance_matrices

    subroutine increment_group_index(this)
        type(Scatterer) :: this

        this%group_index = this%group_index + 1
        if(this%group_index > this%N_sub_groups) then
            this%group_index = 1
        endif

    end subroutine increment_group_index

    subroutine initialize_counters(this,N_order)
        type(Scatterer) :: this
        integer :: i,N_order

        this%N_curr = 0
        do i = 0,N_order
            this%I_sources(:,i) = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
        enddo
    end subroutine initialize_counters

    ! -----------------------------------------------------------------------
    ! Subroutine: solve_scattering_multiscatterer
    ! Purpose   : Simultaneous iterative RAS for multiple scatterers.
    !             Each iteration: add sources to every scatterer, solve each,
    !             evaluate inter-scatterer bouncing fields (IFB step), then
    !             update excitation and check convergence.
    !             Corresponds to RAS_solution_method = 1.
    ! -----------------------------------------------------------------------
    subroutine solve_scattering_multiscatterer(these,N_order,Err_bound,error_val,cnt)
        type(Scatterer),allocatable,intent(inout) :: these(:)
        real(8) :: Err_bound,error_val
        integer :: N_order
        real(8),allocatable,dimension(:) :: error
        real(8),allocatable,dimension(:) :: error_prev
        integer :: cnt,i,j,N_calc
        logical :: use_old
        integer,allocatable,dimension(:) :: n_added,n_sources,n_outside,n_added_outside
        logical :: bouncing_flag,converging_flag

        allocate(n_added(number_of_scatterers),n_sources(number_of_scatterers),&
        n_outside(number_of_scatterers),n_added_outside(number_of_scatterers))
        allocate(error(number_of_scatterers),error_prev(number_of_scatterers))


        error = 1.d0
        cnt = 0
        n_added = 0
        n_outside = 0
        error_val = 1.d0

        do i=1,number_of_scatterers
            n_sources(i) = these(i)%N_group
        enddo
        bouncing_flag = .false.

        do i =1,number_of_scatterers
            call add_all_sources_scatterer(these(i))
            these(i)%N_curr = 0
        enddo
        write(*,*) '========================================'
        write(*,*) '=============Iterative RAS=============='
        converging_flag = .true.
        do while(error_val > Err_bound)
            cnt = cnt+1
            do i = 1,number_of_scatterers

                call increment_group_index(these(i))
                if(these(i)%is_sub_group_used(these(i)%group_index)) then
                    use_old = .true.
                !                    write(*,*) 'Using old group',these(i)%group_index,' for scatterer',these(i)%region_ID-1
                else
                    use_old = .false.
                endif
                N_calc = these(i)%n_sources_array(these(i)%group_index) + &
                these(i)%n_sources_outside_array(these(i)%group_index)
                !                write(*,*) N_calc,these(i)%N_curr
                these(i)%N_curr = these(i)%N_curr + N_calc


                call eval_matrices_RAS(these(i),N_order,N_calc,use_old,these)

                call update_scattered_fields(these(i),N_order,N_calc,these)
            enddo

            !! Bouncing is to be placed here
            do i=1,number_of_scatterers
                do j=1,number_of_scatterers
                    if(i == j) then
                        cycle
                    endif
                    N_calc = these(j)%n_sources_array(these(j)%group_index) + &
                    these(j)%n_sources_outside_array(these(j)%group_index)
                    !                    write(*,*) 'bouncing'
                    call eval_bouncing_field_specified(these(j),N_order,these(i),N_calc,these)
                enddo
            enddo



            do i =1,number_of_scatterers
                call set_Et(these(i),N_order,.true.,these)
                error_prev(i) = error(i)
                !                write(*,*) 'N_calc = ',n_sources(i)+n_added(i)+n_outside(i)
                N_calc = these(i)%n_sources_array(these(i)%group_index) + &
                these(i)%n_sources_outside_array(these(i)%group_index)
                error(i) = eval_error(these(i),N_calc,error_prev(i),.false.)
                !                write(*,*) 'Scatterer #',i, 'error = ',error(i)


                if(error(i) == -2.d0) then
                    write(*,*) 'Best Accuracy Reached without reaching the specified bound at iteration',cnt,'!'
                    error(i) = error_prev(i)
                    converging_flag = .false.
                endif

                if((cnt > 10) .and. (error(i) >= 1.d0)) then
                    write(*,*) 'Solving problem Diverged!'
                    return
                endif

                if(these(i)%group_index == these(i)%N_sub_groups) then
                    call initialize_counters(these(i),N_taylor)


                endif

                !            call cpu_time(end_t)
                !            write(*,*) 'Iteration execution time =',end_t-start_t,'seconds'
                if(these(i)%N_curr > these(i)%N) then
                    these(i)%N = these(i)%N_curr
                endif
            enddo
            error_val = maxval(error)

            write(*,*) 'Iteration #',cnt,' Max. Error =',error_val

            if( .not. converging_flag) then
                exit !! breaks the loop
            endif

            if(cnt > max_iteration) then
                write(*,*) 'Maximum Iterations Reached'
                exit
            endif
        enddo
        if(cnt < max_iteration .and. converging_flag) then
            write(*,*) 'Accuracy Reached in', cnt, 'Iterations'
        endif
        deallocate(n_added,n_sources,n_outside,n_added_outside,error,error_prev)

        do i = 1,number_of_scatterers
            call collect_scattering_fields(these(i))
        enddo

    end subroutine solve_scattering_multiscatterer

    ! -----------------------------------------------------------------------
    ! Subroutine: solve_scattering_multiscatterer_once
    ! Purpose   : Merges all scatterers into a single large Scatterer object
    !             (Scat_big) and solves the combined system once. Sources,
    !             testing points, and excitations from all objects are
    !             concatenated in block order [Scat1 | Scat2 | ... | ScatN].
    !             Results are then redistributed back to the individual
    !             scatterer objects.
    !             Corresponds to RAS_solution_method = 2.
    ! -----------------------------------------------------------------------
    subroutine solve_scattering_multiscatterer_once(these,N_order,Err_bound,error_val,cnt,Scat_big)
        !! this subroutine merges all the scatterers into one big problem then solves it once
        type(Scatterer),allocatable,intent(inout) :: these(:)
        type(Scatterer),intent(inout) :: Scat_big !! the big scatterer to solve
        real(8) :: Err_bound,error_val
        integer :: cnt,i,ind,N_all_sources,N_calc,ii,ex_beg,N_order
        integer,allocatable,dimension(:) :: n_beg_pos


        !! merging technique [Scatterer1, Scatterer2, .... ScattererN]
        Scat_big%M = 0
        do i =1,number_of_scatterers
            Scat_big%M = Scat_big%M + these(i)%M
            call add_all_sources_scatterer(these(i))
        enddo
        allocate(Scat_big%testing_pt(Scat_big%M,3),Scat_big%tang_u(Scat_big%M),&
        Scat_big%tang_v(Scat_big%M),&
        Scat_big%norm_v(Scat_big%M),Scat_big%testing_pt_MoM(Scat_big%M))
        allocate(Scat_big%testing_pt_status(Scat_big%M,3))
        allocate(Scat_big%norm_eta_v(Scat_big%M),Scat_big%norm_eta_u(Scat_big%M))
        allocate(Scat_big%eta_uu(Scat_big%M),Scat_big%eta_vv(Scat_big%M),Scat_big%eta_uv(Scat_big%M),Scat_big%eta_vu(Scat_big%M))
        ind = 1

        Scat_big%eta_uu = 0.d0
        Scat_big%eta_vv = 0.d0
        Scat_big%eta_uv = 0.d0
        Scat_big%eta_vu = 0.d0
        N_all_sources = 0
        Scat_big%N_sub_groups = minval(these(:)%N_sub_groups)
        do i =1,number_of_scatterers
            N_all_sources = N_all_sources + size(these(i)%Source_pos,1)
            Scat_big%testing_pt(ind:(ind+these(i)%M),:) = these(i)%testing_pt
            Scat_big%tang_u(ind:(ind+these(i)%M-1)) = these(i)%tang_u
            Scat_big%tang_v(ind:(ind+these(i)%M-1)) = these(i)%tang_v
            Scat_big%norm_v(ind:(ind+these(i)%M-1)) = these(i)%norm_v
            Scat_big%testing_pt_status(ind:(ind+these(i)%M-1),:) = these(i)%testing_pt_status
            Scat_big%norm_eta_v(ind:(ind+these(i)%M-1)) = these(i)%norm_eta_v
            Scat_big%norm_eta_u(ind:(ind+these(i)%M-1)) = these(i)%norm_eta_u
            if(these(i)%Problem_Type == 3) then !!IBC
                Scat_big%eta_uu(ind:(ind+these(i)%M-1))=these(i)%eta_uu
                Scat_big%eta_vu(ind:(ind+these(i)%M-1))=these(i)%eta_vu
                Scat_big%eta_uv(ind:(ind+these(i)%M-1))=these(i)%eta_uv
                Scat_big%eta_vv(ind:(ind+these(i)%M-1))=these(i)%eta_vv
            endif
            if(these(i)%N_sub_groups < Scat_big%N_sub_groups) then !! find the minimum
                Scat_big%N_sub_groups =these(i)%N_sub_groups
            endif
            ind = ind + these(i)%M
        enddo

        !    do i=1,Scat_big%M
        !        write(*,*) Scat_big%testing_pt_status(i,:)
        !    enddo
        Scat_big%Problem_Type = these(1)%Problem_Type
        Scat_big%k_local = these(1)%k_local
        Scat_big%eta_local = these(1)%eta_local
        Scat_big%N_max = N_all_sources
        call allocate_arrays(Scat_big,.false.,these)
        ind = 1


        do ii =0,1
            do i=1,number_of_scatterers

                ex_beg = ii*these(i)%M+1
                Scat_big%E(ind:(ind+these(i)%M-1),:) = these(i)%E(ex_beg:(ex_beg+these(i)%M-1),:)
                Scat_big%H(ind:(ind+these(i)%M-1),:) = these(i)%H(ex_beg:(ex_beg+these(i)%M-1),:)



                ind = ind + these(i)%M
            !        write(*,*) these(i)%N_sub_groups
            enddo
        enddo

        !    write(*,*) 'H=', sum(abs(Scat_big%H))/sqrt(Scat_big%normalize_H)
        !    write(*,*) 'E=', sum(abs(Scat_big%E))/sqrt(Scat_big%normalize_E)
        allocate(Scat_big%normalize_E(0:N_order),Scat_big%normalize_H(0:N_order))

        Scat_big%normalize_E = 0.d0
        Scat_big%normalize_H = 0.d0
        do i = 1,number_of_scatterers
            Scat_big%normalize_E =Scat_big%normalize_E +these(i)%normalize_E
            !        write(*,*) Scat_big%normalize_E,these(i)%normalize_E
            Scat_big%normalize_H =Scat_big%normalize_H +these(i)%normalize_H
        enddo

        call initialize_stage(Scat_big,N_order,Scat_big%E,Scat_big%H,.true.,these)



        !! for the sources that are allready added

        allocate(Scat_big%is_sub_group_used(Scat_big%N_sub_groups))

        allocate(Scat_big%n_sources_array(Scat_big%N_sub_groups),Scat_big%n_sources_outside_array(Scat_big%N_sub_groups))
        Scat_big%n_sources_outside_array = 0
        Scat_big%n_sources_array = 0
        Scat_big%Source_pos = 0.d0
        Scat_big%I_stat = 0
        Scat_big%active_region = 0

        !    write(*,*) Scat_big%n_sources_array
        do i=1,number_of_scatterers
            Scat_big%n_sources_array = Scat_big%n_sources_array + these(i)%n_sources_array(1:Scat_big%N_sub_groups)
            Scat_big%n_sources_outside_array = Scat_big%n_sources_outside_array + &
            these(i)%n_sources_outside_array(1:Scat_big%N_sub_groups)
        enddo
        !    write(*,*) Scat_big%n_sources_array
        !    write(*,*) Scat_big%n_sources_outside_array
        allocate(n_beg_pos(number_of_scatterers))
        n_beg_pos = 0
        Scat_big%is_sub_group_used = .false.
        ind = 1

        do ii = 1,Scat_big%N_sub_groups
            do i=1,number_of_scatterers

                N_calc = these(i)%n_sources_array(ii) + these(i)%n_sources_outside_array(ii)

                Scat_big%Source_pos(ind:(ind-1+N_calc),:) = &
                these(i)%Source_pos((n_beg_pos(i)+1):(n_beg_pos(i)+N_calc),:)

                Scat_big%active_region(ind:(ind-1+N_calc)) = &
                these(i)%active_region((n_beg_pos(i)+1):(n_beg_pos(i)+N_calc))

                Scat_big%I_stat(ind:(ind-1+N_calc)) = &
                these(i)%I_stat((n_beg_pos(i)+1):(n_beg_pos(i)+N_calc))

                Scat_big%allowed_orientation(ind:(ind-1+N_calc)) = &
                these(i)%allowed_orientation((n_beg_pos(i)+1):(n_beg_pos(i)+N_calc))

                ind = ind + N_calc
                n_beg_pos(i) = n_beg_pos(i) + N_calc

            enddo
        enddo




        !    OPEN(25, FILE='all_sources.dat')
        !    do i=1,N_all_sources
        !        write(25,*) Scat_big%Source_pos(i,:),Scat_big%I_stat(i),Scat_big%In_out_flag(i)
        !    enddo
        !    CLOSE(25)




        !    write(*,*) 'NOTE: Scat_big method only solves problems with similar boundary conditions'

        call solve_scattering(Scat_big,N_order,Err_bound,error_val,cnt,.false.,these)

        call collect_scattering_fields(Scat_big)

        ind = 1

        do ii = 0,1
            do i=1,number_of_scatterers
                ex_beg = ii*these(i)%M + 1

                these(i)%Es(ex_beg:(ex_beg+these(i)%M-1),:) = Scat_big%Es(ind:(ind-1+these(i)%M),:)

                these(i)%Hs(ex_beg:(ex_beg+these(i)%M-1),:) = Scat_big%Hs(ind:(ind-1+these(i)%M),:)
                ind = ind + these(i)%M


            enddo
        enddo

        n_beg_pos = 0
        ind = 1

        do ii = 1,Scat_big%N_sub_groups
            do i=1,number_of_scatterers

                N_calc = these(i)%n_sources_array(ii) + these(i)%n_sources_outside_array(ii)

                these(i)%I_sources_total((n_beg_pos(i)+1):(n_beg_pos(i)+N_calc),:) =&
                Scat_big%I_sources_total(ind:(ind-1+N_calc),:)

                these(i)%N =these(i)%N + N_calc

                ind = ind + N_calc
                n_beg_pos(i) = n_beg_pos(i) + N_calc

            enddo
        enddo
        deallocate(n_beg_pos)
    !    stop
    !    write(*,*) sum(abs(these(1)%Es))




    end subroutine solve_scattering_multiscatterer_once

    ! -----------------------------------------------------------------------
    ! Subroutine: solve_scattering_multiscatterer_2
    ! Purpose   : Two-stage IFB (Iterative Farfield Bouncing) multi-scatterer
    !             solver. Stage 1: solve each scatterer independently under
    !             the incident wave. Stage 2: iterate IFB loops -- at each
    !             step each scatterer's scattered field is used as the
    !             excitation for its neighbours until bouncing error < Err_bound.
    !             Implements the method from Moharram & Kishk (APSURSI 2013).
    ! -----------------------------------------------------------------------
    subroutine solve_scattering_multiscatterer_2(these,N_order,Err_bound,error_val,cnt)
        !! this subroutine solves the multiscatterer problem using the IFB technique introduced in
        !! Volakis paper 2002. It solves first with the incident wave only. then it solves each scatterer
        !! for each bouncing level updating the dipole moments
        !! (Same as the MoM)
        type(Scatterer),allocatable,intent(inout) :: these(:)
        real(8) :: Err_bound,error_val,error_val1
        real(8),allocatable,dimension(:) :: error,error_bouncing
        real(8),allocatable,dimension(:) :: error_prev
        integer :: cnt,i,j,itr_bouncing,N_order
        integer,allocatable,dimension(:) :: itr_counter

        allocate(itr_counter(number_of_scatterers))
        allocate(error(number_of_scatterers),error_prev(number_of_scatterers))
        !! solve scatterers with incident wave only
        do i =1,number_of_scatterers
            call add_all_sources_scatterer(these(i))

            call solve_scattering(these(i),N_order,Err_bound,error(i),itr_counter(i),.false.,these)

            call collect_scattering_fields(these(i))
        enddo
        if(number_of_scatterers == 1) then
            write(*,*) 'RAS Error =',error(1), '  in  ',itr_counter(1), 'Iterations'
        endif


        if(number_of_scatterers > 1) then !! bouncing is required
            allocate(error_bouncing(number_of_scatterers))
            write(*,*) '=============================================='
            write(*,*) '==========RAS-IFB Iterations start============'
            !! Evaluate the bouncing fields for the first time
            do i=1,number_of_scatterers
                !! nullify the bouncing fields
                these(i)%E_bounce_current = 0.d0
                these(i)%H_bounce_current = 0.d0
                do j=1,number_of_scatterers
                    if(i == j) then
                        cycle
                    endif
                    call eval_bouncing_field(these(j),N_taylor,these(i),these)
                enddo
                !                error_bouncing(i) = sum(abs(these(i)%E_bounce)**2.d0)/these(i)%normalize_E
                error_bouncing(i) = sum(abs(these(i)%E_bounce_current)**2.d0)/sum(abs(these(i)%Es)**2.d0)
            enddo
            error_val = maxval(error_bouncing)

            itr_bouncing = 1
            do while(error_val > Err_bound)
                if(itr_bouncing > max_bouncing_iteration) then
                    write(*,*) 'WARNING: RAS_IFB did not converge to the desired error bound'
                    write(*,*) 'Please, consider increasing MAX_BOUNCING_LIMIT variable from parameters.dat file'
                    exit
                endif

                !! modify moments to account for the bouncing
                do i =1,number_of_scatterers
                    !! invert signs
                    !                    these(i)%E_bounce = these(i)%E_bounce
                    !                    these(i)%H_bounce = these(i)%H_bounce
                    call initialize_stage(these(i),N_order,these(i)%E_bounce_current,these(i)%H_bounce_current,.false.,these)
                    these(i)%itr_counter = 1
                    these(i)%converging = .false.
                    call solve_scattering(these(i),N_order,Err_bound,error(i),itr_counter(i),.true.,these)
                    call collect_scattering_fields(these(i))
                enddo
                error_val1 = maxval(error)
                !                write(*,*) 'Error after solution',error_val1
                !                write(*,*) 'Es MAG=',sum(abs(these(1)%Es))/sqrt(these(1)%normalize_E)
                !                write(*,*) abs(these(1)%Es(1:5))

                !! evaluate bouncing fields
                do i=1,number_of_scatterers
                    !! nullify the bouncing fields
                    these(i)%E_bounce_current = 0.d0
                    these(i)%H_bounce_current = 0.d0
                    do j=1,number_of_scatterers
                        if(i == j) then
                            cycle
                        endif
                        call eval_bouncing_field(these(j),N_order,these(i),these)
                    enddo
                    !                    call set_Et(these(i))
                    !                error_bouncing(i) = sum(abs(these(i)%E_bounce)**2.d0)/these(i)%normalize_E
                    error_bouncing(i) = sum(abs(these(i)%E_bounce_current)**2.d0)/sum(abs(these(i)%Es)**2.d0)
                enddo
                error_val = maxval(error_bouncing)
                write(*,*) 'RAS-IFB Iteration#',itr_bouncing,' Error',error_val

                itr_bouncing =itr_bouncing +1
            enddo
            if(itr_bouncing <= max_bouncing_iteration) then
                write(*,*) 'RAS-IFB converged in (',itr_bouncing-1,' ) Iterations.'

            endif
            !            if(itr_bouncing >= max_bouncing_iteration) then
            !
            !            endif
            deallocate(error_bouncing)
        else
            error_val = error(1)
            cnt = itr_counter(1)
        endif
        deallocate(error,error_prev)
        deallocate(itr_counter)
    end subroutine solve_scattering_multiscatterer_2

    subroutine collect_scattering_fields(this)
        type(Scatterer) :: this





        !        write(*,*) abs(this%Es_current(1:5))
        this%Es = this%Es + this%Es_current
        this%Hs = this%Hs + this%Hs_current

        this%E_bounce = this%E_bounce + this%E_bounce_current
        this%H_bounce = this%H_bounce + this%H_bounce_current

        !        write(*,*) abs(this%Es(1:5))
        if(this%Problem_Type == 4) then !! only dielectric case
            this%Ed = this%Ed + this%Ed_current
            this%Hd = this%Hd + this%Hd_current
        endif
        this%I_sources_total = this%I_sources_total + this%I_sources_total_current


    !        write(*,*) 'Currents Addition value =',sum(abs(this%I_sources_total_current))
    end subroutine collect_scattering_fields

    subroutine add_all_sources_scatterer(this)
        type(Scatterer) :: this
        integer :: n_sub_g,n_sources,n_outside
        integer :: n_added_ext,n_added_int,n_added_outside,Scr_Type


        allocate(this%n_sources_array(this%N_sub_groups),this%n_sources_outside_array(this%N_sub_groups))
        allocate(this%is_sub_group_used(this%N_sub_groups))



        n_added_ext = 0
        n_added_int = 0
        Scr_Type = 1
        do n_sub_g = 1,this%N_sub_groups
            this%is_sub_group_used(n_sub_g) = .false.
            n_sources = this%N_group

            if(Source_Placement == 1 .or. Source_Placement == 3) then
                if(this%Problem_Type == 4) then
                    call add_sources_contour(this,n_sources,n_added_ext,n_added_int,.false.)
                else
                    !                n_sources = 2
                    call add_sources_contour(this,n_sources,n_added_ext,n_added_int,.true.)
                endif
            endif
!            write(*,*) 'n_added_ext,n_added_int',n_added_ext,n_added_int
            !            write(*,*) 'n_added=',n_added,this%N_con
            !            if(This%Problem_Type == 1) then
            !                if(this%Contour_Cnt >= this%N_con) then
            !                    this%Contour_Cnt = 1
            !
            !                endif
            !            endif
            !            write(*,*) 'N_group',this%N_group,'N_added',n_added

            n_sources = this%N_group - n_added_ext
!            n_sources = this%N_group

            call add_sources_rect(this,n_sources)

            !            this%n_sources_array(n_sub_g) = n_sources
            this%n_sources_array(n_sub_g) = this%N_group !+ n_added_ext !!1 remove n_added_ext
            if(this%Problem_Type == 4) then !! for Dielectric Cases only

                !                n_outside = this%N_inside_sources
                !                call add_sources_outside(this,n_outside,n_added_outside)
                n_outside = this%N_group  - n_added_int
                call add_sources_outside_random(this,n_outside,n_added_outside)
                this%Contour_out_Cnt = 1
                this%n_sources_outside_array(n_sub_g) = n_added_outside + n_added_int
            else
                this%n_sources_outside_array(n_sub_g) = 0
            endif
        enddo


    end subroutine add_all_sources_scatterer



    ! -----------------------------------------------------------------------
    ! Subroutine: solve_scattering
    ! Purpose   : Core iterative RAS solver for a single scatterer.
    !             Each iteration adds N_group sources (one sub-group),
    !             calls eval_matrices_RAS to build/update the Gram matrix,
    !             solves for new source currents, evaluates the scattered
    !             fields, and computes the normalised boundary error.
    !             Converges when error < Err_bound or max_iteration is hit.
    !   disable_prim_source=.true.: uses only secondary (bounce) excitation
    !             (used in IFB multi-scatterer bouncing iterations).
    ! -----------------------------------------------------------------------
    subroutine solve_scattering(this,N_order,Err_bound,error,cnt,disable_prim_source,Scatterers_pointer)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)
        real(8) :: Err_bound,error,error_prev
        integer :: cnt,N_calc,N_order
        logical :: use_old,disable_prim_source
        !        real :: start_t,end_t
        !! disable_prim_source logical variable indicates whether to use the primary source (plane wave)
        !! or just to iterate to satisfy the boundary conditions (may be due to other secondary sources)
        error = 1.d0
        cnt = 0


        !        n_sources = this%N_group
        if(disable_prim_source) then
            !            this%E = 0.d0
            !            this%H = 0.d0
            error = 2.d0

        endif
        this%N_curr = 0
        do while(error > Err_bound)
            !            call cpu_time(start_t)
            cnt = cnt+1

            call increment_group_index(this)
            if(this%is_sub_group_used(this%group_index)) then
                use_old = .true.
            !                write(*,*) 'Using old group',this%group_index,' for scatterer',this%region_ID-1
            else
                use_old = .false.
            endif
            N_calc = this%n_sources_array(this%group_index) + &
            this%n_sources_outside_array(this%group_index)

            this%N_curr = this%N_curr + N_calc

            !            write(*,*) sum(abs(this%EH_excit))/sqrt(this%normalize_E)

            !            write(*,*) 'Iteration#',cnt,'this%N_curr=',this%N_curr,'N_sources =',n_sources,'group_index =',&
            !                    this%group_index,'N_outside',n_outside

            call eval_matrices_RAS(this,N_order,N_calc,use_old,Scatterers_pointer)

            call update_scattered_fields(this,N_order,N_calc,Scatterers_pointer)

            call set_Et(this,N_order,.false.,Scatterers_pointer)

            error_prev = error
            error = eval_error(this,N_calc,error_prev,.false.)

!                         write(*,*) 'error = ',error, N_calc,this%N_curr


            if(error == -2.d0) then
                write(*,*) 'Best Accuracy Reached without reaching the specified bound at iteration',cnt,'!'
                error = error_prev
                return
            endif

            if((cnt > 10) .and. (error >= 1.d0)) then
                write(*,*) 'Solving problem Diverged!'
                !                stop
                error = -2.d0
                return
            endif



            if(this%group_index == this%N_sub_groups) then
                call initialize_counters(this,N_taylor)

            endif
            !            call cpu_time(end_t)
            !            write(*,*) 'Iteration execution time =',end_t-start_t,'seconds'
            if(this%N_curr > this%N) then
                this%N = this%N_curr
            endif

            if(cnt > max_iteration) then
                write(*,*) 'Maximum Iterations Reached'
                exit
            endif

        enddo
    !        if(.not. disable_prim_source) then
    !            if(cnt < max_iteration) then
    !                write(*,*) 'Accuracy Reached in', cnt, 'Iterations'
    !            endif
    !            write(*,*) 'Source Scatterer',this%region_ID-1,'N sources =',this%N
    !        endif


    end subroutine solve_scattering


    ! -----------------------------------------------------------------------
    ! Subroutine: get_LSM_excitation
    ! Purpose   : Builds the RAS right-hand side (excitation) vector b for
    !             the N_calc most recently added sources.
    !             For each testing point kk and each source j:
    !               - evaluates the electric/magnetic kernel (Electric or
    !                 Magnetic depending on I_stat)
    !               - projects onto tangential unit vectors tang_u, tang_v
    !               - applies boundary condition weighting (PEC/PMC/IBC/Diel)
    !               - accumulates into b using the AWE binomial coefficients
    !                 n_C_i_matrix for Taylor-order N_order.
    !             b has dimensions (3*N_calc, 0:N_order).
    ! -----------------------------------------------------------------------
    subroutine get_LSM_excitation(this,N_order,N_calc,b,Scatterers_Lib)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        integer :: N_calc,index_mat
        integer :: N_order !! number of orders to be considered for the Taylor expansion
        complex*16,allocatable,intent(inout) :: b(:,:)
        integer :: j,kk,jj,Nm,nn,ii
        type(vector_c),allocatable,dimension(:) :: E_Sx,H_Sx,E_Sy,H_Sy,E_Sz,H_Sz

        complex*16,allocatable,dimension(:,:) :: EH_Sx_u,EH_Sy_u,EH_Sz_u
        complex*16,allocatable,dimension(:,:) :: EH_Sx_v,EH_Sy_v,EH_Sz_v
        complex*16,allocatable,dimension(:,:) :: E_Sx_u,E_Sy_u,E_Sz_u
        complex*16,allocatable,dimension(:,:) :: E_Sx_v,E_Sy_v,E_Sz_v
        complex*16,allocatable,dimension(:,:) :: H_Sx_u,H_Sy_u,H_Sz_u
        complex*16,allocatable,dimension(:,:) :: H_Sx_v,H_Sy_v,H_Sz_v
        complex*16,allocatable,dimension(:,:) :: conjg_EH_Sx_u,conjg_EH_Sy_u,conjg_EH_Sz_u
        complex*16,allocatable,dimension(:,:) :: conjg_EH_Sx_v,conjg_EH_Sy_v,conjg_EH_Sz_v
        complex*16,allocatable,dimension(:,:) :: conjg_H_Sx_u,conjg_H_Sy_u,conjg_H_Sz_u
        complex*16,allocatable,dimension(:,:) :: conjg_H_Sx_v,conjg_H_Sy_v,conjg_H_Sz_v

        complex*16,allocatable,dimension(:) :: EH1,EH2,EH3,EH4
        complex*16,allocatable,dimension(:) :: value_to_add
        integer :: current_BC,inner_region_id,outer_region_id
        complex*16 :: ar_in,ak_in,ar_ex,ak_ex,eta_ex,eta_in

        Nm = 3*N_calc
        allocate(b(Nm,0:N_order))
        allocate(EH_Sx_u(N_calc,0:N_order),EH_Sy_u(N_calc,0:N_order),EH_Sz_u(N_calc,0:N_order))
        allocate(EH_Sx_v(N_calc,0:N_order),EH_Sy_v(N_calc,0:N_order),EH_Sz_v(N_calc,0:N_order))
        allocate(E_Sx_u(N_calc,0:N_order),E_Sy_u(N_calc,0:N_order),E_Sz_u(N_calc,0:N_order))
        allocate(E_Sx_v(N_calc,0:N_order),E_Sy_v(N_calc,0:N_order),E_Sz_v(N_calc,0:N_order))
        allocate(H_Sx_u(N_calc,0:N_order),H_Sy_u(N_calc,0:N_order),H_Sz_u(N_calc,0:N_order))
        allocate(H_Sx_v(N_calc,0:N_order),H_Sy_v(N_calc,0:N_order),H_Sz_v(N_calc,0:N_order))
        allocate(conjg_EH_Sx_u(N_calc,0:N_order),conjg_EH_Sy_u(N_calc,0:N_order),conjg_EH_Sz_u(N_calc,0:N_order))
        allocate(conjg_EH_Sx_v(N_calc,0:N_order),conjg_EH_Sy_v(N_calc,0:N_order),conjg_EH_Sz_v(N_calc,0:N_order))
        allocate(conjg_H_Sx_u(N_calc,0:N_order),conjg_H_Sy_u(N_calc,0:N_order),conjg_H_Sz_u(N_calc,0:N_order))
        allocate(conjg_H_Sx_v(N_calc,0:N_order),conjg_H_Sy_v(N_calc,0:N_order),conjg_H_Sz_v(N_calc,0:N_order))
        allocate(EH1(0:N_order),EH2(0:N_order),EH4(0:N_order),EH3(0:N_order))
        allocate(E_Sx(0:N_order),H_Sx(0:N_order),E_Sy(0:N_order),H_Sy(0:N_order),E_Sz(0:N_order),H_Sz(0:N_order))
        allocate(value_to_add(0:N_order))
        b = 0.d0

        do kk=1,this%M
            current_BC = this%testing_pt_status(kk,1)
            if(current_BC == -1) then
                cycle
            endif
            inner_region_id = this%testing_pt_status(kk,2)
            outer_region_id = this%testing_pt_status(kk,3)
            if(outer_region_id == 0) then
                ar_ex = ar
                ak_ex = ak
                eta_ex = eta1
            else
                ar_ex = Scatterers_Lib(outer_region_id)%ar_local
                ak_ex = Scatterers_Lib(outer_region_id)%ak_local
                eta_ex = Scatterers_Lib(outer_region_id)%eta_local
            endif
            ar_in = Scatterers_Lib(inner_region_id)%ar_local
            ak_in = Scatterers_Lib(inner_region_id)%ak_local
            eta_in = Scatterers_Lib(inner_region_id)%eta_local

            EH1 = this%EH_excit(kk,:)
            EH2 = this%EH_excit(kk+this%M,:)
            EH3 = this%EH_excit(kk+2*this%M,:)
            EH4 = this%EH_excit(kk+3*this%M,:)
            do j =1,N_calc  !! evaluate source contributions once
                jj = j + (this%N_curr-N_calc)
                if(this%active_region(jj)==outer_region_id) then !! inside sources
                    if(this%I_stat(jj) == 1) then
                        call eval_kernel_Electric_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,eta_ex,az,ar_ex,ak_ex,kk,jj)
                    else
                        call eval_kernel_Magnetic_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,eta_ex,az,ar_ex,ak_ex,kk,jj)
                    endif
                    E_Sx_u(j,:) = dot(this%tang_u(kk),E_Sx,N_order+1)/eta_ex
                    E_Sx_v(j,:) = dot(this%tang_v(kk),E_Sx,N_order+1)/eta_ex
                    E_Sy_u(j,:) = dot(this%tang_u(kk),E_Sy,N_order+1)/eta_ex
                    E_Sy_v(j,:) = dot(this%tang_v(kk),E_Sy,N_order+1)/eta_ex
                    E_Sz_u(j,:) = dot(this%tang_u(kk),E_Sz,N_order+1)/eta_ex
                    E_Sz_v(j,:) = dot(this%tang_v(kk),E_Sz,N_order+1)/eta_ex

                    H_Sx_u(j,:) = dot(this%tang_u(kk),H_Sx,N_order+1)
                    H_Sx_v(j,:) = dot(this%tang_v(kk),H_Sx,N_order+1)
                    H_Sy_u(j,:) = dot(this%tang_u(kk),H_Sy,N_order+1)
                    H_Sy_v(j,:) = dot(this%tang_v(kk),H_Sy,N_order+1)
                    H_Sz_u(j,:) = dot(this%tang_u(kk),H_Sz,N_order+1)
                    H_Sz_v(j,:) = dot(this%tang_v(kk),H_Sz,N_order+1)
                elseif(this%active_region(jj) == inner_region_id) then !! outside sources producing the dielectric fields
                    if(this%I_stat(jj) == 1) then
                        call eval_kernel_Electric_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        eta_in,az,ar_in,ak_in,kk,jj)
                    else
                        call eval_kernel_Magnetic_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        eta_in,az,ar_in,ak_in,kk,jj)
                    endif
                    E_Sx_u(j,:) = -dot(this%tang_u(kk),E_Sx,N_order+1)/eta_ex
                    E_Sx_v(j,:) = -dot(this%tang_v(kk),E_Sx,N_order+1)/eta_ex
                    E_Sy_u(j,:) = -dot(this%tang_u(kk),E_Sy,N_order+1)/eta_ex
                    E_Sy_v(j,:) = -dot(this%tang_v(kk),E_Sy,N_order+1)/eta_ex
                    E_Sz_u(j,:) = -dot(this%tang_u(kk),E_Sz,N_order+1)/eta_ex
                    E_Sz_v(j,:) = -dot(this%tang_v(kk),E_Sz,N_order+1)/eta_ex

                    H_Sx_u(j,:) = -dot(this%tang_u(kk),H_Sx,N_order+1)
                    H_Sx_v(j,:) = -dot(this%tang_v(kk),H_Sx,N_order+1)
                    H_Sy_u(j,:) = -dot(this%tang_u(kk),H_Sy,N_order+1)
                    H_Sy_v(j,:) = -dot(this%tang_v(kk),H_Sy,N_order+1)
                    H_Sz_u(j,:) = -dot(this%tang_u(kk),H_Sz,N_order+1)
                    H_Sz_v(j,:) = -dot(this%tang_v(kk),H_Sz,N_order+1)
                else !! neither outer or inner regions [EVERYTHING should be zero in this case]
                    E_Sx_u(j,:) = 0.d0
                    E_Sx_v(j,:) = 0.d0
                    E_Sy_u(j,:) = 0.d0
                    E_Sy_v(j,:) = 0.d0
                    E_Sz_u(j,:) = 0.d0
                    E_Sz_v(j,:) = 0.d0

                    H_Sx_u(j,:) = 0.d0
                    H_Sx_v(j,:) = 0.d0
                    H_Sy_u(j,:) = 0.d0
                    H_Sy_v(j,:) = 0.d0
                    H_Sz_u(j,:) = 0.d0
                    H_Sz_v(j,:) = 0.d0
                endif

            enddo

            if(current_BC == 1) then
                EH_Sx_u = E_Sx_u
                EH_Sx_v = E_Sx_v
                EH_Sy_u = E_Sy_u
                EH_Sy_v = E_Sy_v
                EH_Sz_u = E_Sz_u
                EH_Sz_v = E_Sz_v

                conjg_EH_Sx_u = conjg(EH_Sx_u)
                conjg_EH_Sx_v = conjg(EH_Sx_v)
                conjg_EH_Sy_u = conjg(EH_Sy_u)
                conjg_EH_Sy_v = conjg(EH_Sy_v)
                conjg_EH_Sz_u = conjg(EH_Sz_u)
                conjg_EH_Sz_v = conjg(EH_Sz_v)

                conjg_H_Sx_u = 0.d0
                conjg_H_Sx_v = 0.d0
                conjg_H_Sy_u = 0.d0
                conjg_H_Sy_v = 0.d0
                conjg_H_Sz_u = 0.d0
                conjg_H_Sz_v = 0.d0
                H_Sx_u = 0.d0
                H_Sx_v = 0.d0
                H_Sy_u = 0.d0
                H_Sy_v = 0.d0
                H_Sz_u = 0.d0
                H_Sz_v = 0.d0
            elseif(current_BC == 2) then
                EH_Sx_u = 0.d0
                EH_Sx_v = 0.d0
                EH_Sy_u = 0.d0
                EH_Sy_v = 0.d0
                EH_Sz_u = 0.d0
                EH_Sz_v = 0.d0

                conjg_EH_Sx_u = 0.d0
                conjg_EH_Sx_v = 0.d0
                conjg_EH_Sy_u = 0.d0
                conjg_EH_Sy_v = 0.d0
                conjg_EH_Sz_u = 0.d0
                conjg_EH_Sz_v = 0.d0

                conjg_H_Sx_u = conjg(H_Sx_u)
                conjg_H_Sx_v = conjg(H_Sx_v)
                conjg_H_Sy_u = conjg(H_Sy_u)
                conjg_H_Sy_v = conjg(H_Sy_v)
                conjg_H_Sz_u = conjg(H_Sz_u)
                conjg_H_Sz_v = conjg(H_Sz_v)
            elseif(current_BC == 3) then
                EH_Sx_u = this%norm_eta_u(kk)*(E_Sx_u + (this%eta_uu(kk)*H_Sx_v - this%eta_uv(kk)*H_Sx_u))
                EH_Sx_v = this%norm_eta_v(kk)*(E_Sx_v + (this%eta_vu(kk)*H_Sx_v - this%eta_vv(kk)*H_Sx_u))
                EH_Sy_u = this%norm_eta_u(kk)*(E_Sy_u + (this%eta_uu(kk)*H_Sy_v - this%eta_uv(kk)*H_Sy_u))
                EH_Sy_v = this%norm_eta_v(kk)*(E_Sy_v + (this%eta_vu(kk)*H_Sy_v - this%eta_vv(kk)*H_Sy_u))
                EH_Sz_u = this%norm_eta_u(kk)*(E_Sz_u + (this%eta_uu(kk)*H_Sz_v - this%eta_uv(kk)*H_Sz_u))
                EH_Sz_v = this%norm_eta_v(kk)*(E_Sz_v + (this%eta_vu(kk)*H_Sz_v - this%eta_vv(kk)*H_Sz_u))

                !                ZY_Sx_u(j,:) = this%norm_eta_u*(E_Sx_u(j,:) + (this%eta_uu(kk)*H_Sx_v(j,:) - this%eta_uv(kk)*H_Sx_u(j,:)))
                !                ZY_Sy_u(j,:) = this%norm_eta_u*(E_Sy_u(j,:) + (this%eta_uu(kk)*H_Sy_v(j,:) - this%eta_uv(kk)*H_Sy_u(j,:)))
                !                ZY_Sz_u(j,:) = this%norm_eta_u*(E_Sz_u(j,:) + (this%eta_uu(kk)*H_Sz_v(j,:) - this%eta_uv(kk)*H_Sz_u(j,:)))
                !
                !                ZY_Sx_v(j,:) = this%norm_eta_v*(E_Sx_v(j,:) + (this%eta_vu(kk)*H_Sx_v(j,:) - this%eta_vv(kk)*H_Sx_u(j,:)))
                !                ZY_Sy_v(j,:) = this%norm_eta_v*(E_Sy_v(j,:) + (this%eta_vu(kk)*H_Sy_v(j,:) - this%eta_vv(kk)*H_Sy_u(j,:)))
                !                ZY_Sz_v(j,:) = this%norm_eta_v*(E_Sz_v(j,:) + (this%eta_vu(kk)*H_Sz_v(j,:) - this%eta_vv(kk)*H_Sz_u(j,:)))


                conjg_EH_Sx_u = conjg(EH_Sx_u)
                conjg_EH_Sx_v = conjg(EH_Sx_v)
                conjg_EH_Sy_u = conjg(EH_Sy_u)
                conjg_EH_Sy_v = conjg(EH_Sy_v)
                conjg_EH_Sz_u = conjg(EH_Sz_u)
                conjg_EH_Sz_v = conjg(EH_Sz_v)

                conjg_H_Sx_u = 0.d0
                conjg_H_Sx_v = 0.d0
                conjg_H_Sy_u = 0.d0
                conjg_H_Sy_v = 0.d0
                conjg_H_Sz_u = 0.d0
                conjg_H_Sz_v = 0.d0
                H_Sx_u = 0.d0
                H_Sx_v = 0.d0
                H_Sy_u = 0.d0
                H_Sy_v = 0.d0
                H_Sz_u = 0.d0
                H_Sz_v = 0.d0
            else
                EH_Sx_u = E_Sx_u
                EH_Sx_v = E_Sx_v
                EH_Sy_u = E_Sy_u
                EH_Sy_v = E_Sy_v
                EH_Sz_u = E_Sz_u
                EH_Sz_v = E_Sz_v

                conjg_EH_Sx_u = conjg(EH_Sx_u)
                conjg_EH_Sx_v = conjg(EH_Sx_v)
                conjg_EH_Sy_u = conjg(EH_Sy_u)
                conjg_EH_Sy_v = conjg(EH_Sy_v)
                conjg_EH_Sz_u = conjg(EH_Sz_u)
                conjg_EH_Sz_v = conjg(EH_Sz_v)

                conjg_H_Sx_u = conjg(H_Sx_u)
                conjg_H_Sx_v = conjg(H_Sx_v)
                conjg_H_Sy_u = conjg(H_Sy_u)
                conjg_H_Sy_v = conjg(H_Sy_v)
                conjg_H_Sz_u = conjg(H_Sz_u)
                conjg_H_Sz_v = conjg(H_Sz_v)
            endif




            do j=1,N_calc
                value_to_add = 0.d0
                do nn = 0,N_order
                    do ii = 0,nn
                        value_to_add(nn) = value_to_add(nn) - dble(n_C_i_matrix(nn,ii))*(conjg_EH_Sx_u(j,nn-ii)*EH1(ii) +&
                        conjg_EH_Sx_v(j,nn-ii)*EH2(ii) + conjg_H_Sx_u(j,nn-ii)*EH3(ii) +&
                        conjg_H_Sx_v(j,nn-ii)*EH4(ii))
                    enddo
                enddo
                b(j,:)= b(j,:)+value_to_add
            enddo
            index_mat = 1 + N_calc
            do j=1,N_calc
                value_to_add = 0.d0
                do nn = 0,N_order
                    do ii = 0,nn
                        value_to_add(nn) = value_to_add(nn) - dble(n_C_i_matrix(nn,ii))*(conjg_EH_Sy_u(j,nn-ii)*EH1(ii) +&
                        conjg_EH_Sy_v(j,nn-ii)*EH2(ii) + conjg_H_Sy_u(j,nn-ii)*EH3(ii) +&
                        conjg_H_Sy_v(j,nn-ii)*EH4(ii))
                    enddo
                enddo
                b(index_mat,:) = b(index_mat,:) + value_to_add
                index_mat = index_mat +1
            enddo
            do j=1,N_calc
                value_to_add = 0.d0
                do nn = 0,N_order
                    do ii = 0,nn
                        value_to_add(nn) = value_to_add(nn) - dble(n_C_i_matrix(nn,ii))*(conjg_EH_Sz_u(j,nn-ii)*EH1(ii) +&
                        conjg_EH_Sz_v(j,nn-ii)*EH2(ii) + conjg_H_Sz_u(j,nn-ii)*EH3(ii) +&
                        conjg_H_Sz_v(j,nn-ii)*EH4(ii))
                    enddo
                enddo
                b(index_mat,:) = b(index_mat,:) + value_to_add

                index_mat = index_mat +1
            enddo
        enddo


        deallocate(E_Sx,H_Sx,E_Sy,H_Sy,E_Sz,H_Sz)
        deallocate(EH_Sx_u,EH_Sy_u,EH_Sz_u,EH_Sx_v,EH_Sy_v,EH_Sz_v)
        deallocate(E_Sx_u,E_Sy_u,E_Sz_u,E_Sx_v,E_Sy_v,E_Sz_v)
        deallocate(H_Sx_u,H_Sy_u,H_Sz_u,H_Sx_v,H_Sy_v,H_Sz_v)
        deallocate(conjg_EH_Sx_u,conjg_EH_Sy_u,conjg_EH_Sz_u,conjg_EH_Sx_v,conjg_EH_Sy_v,conjg_EH_Sz_v)
        deallocate(conjg_H_Sx_u,conjg_H_Sy_u,conjg_H_Sz_u,conjg_H_Sx_v,conjg_H_Sy_v,conjg_H_Sz_v)
        deallocate(EH1,EH2,EH3,EH4)
        deallocate(value_to_add)
    end subroutine get_LSM_excitation



    ! -----------------------------------------------------------------------
    ! Subroutine: get_LSM_matrices
    ! Purpose   : Builds the RAS Gram matrix Y_vec and excitation vector b
    !             for the N_calc most recently added sources.
    !             Y_vec is stored in packed Hermitian form (lower triangle)
    !             of size N_out = (3*N_calc)*(3*N_calc+1)/2.
    !             For each (testing point, source) pair, computes the kernel
    !             and accumulates into Y and b using AWE binomial weights.
    !             Boundary condition type (PEC/PMC/IBC/Dielectric) is applied
    !             per testing point via testing_pt_status.
    ! -----------------------------------------------------------------------
    subroutine get_LSM_matrices(this,N_order,N_calc,Y_vec,b,N_out,Scatterers_Lib)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        integer :: N_calc,N_out,index_mat
        integer :: N_order !! number of orders to be considered for the Taylor expansion
        complex*16,allocatable,intent(inout) :: Y_vec(:,:),b(:,:)
        integer :: i,j,kk,jj,Nm,nn,ii
        type(vector_c),allocatable,dimension(:) :: E_Sx,H_Sx,E_Sy,H_Sy,E_Sz,H_Sz

        complex*16,allocatable,dimension(:,:) :: EH_Sx_u,EH_Sy_u,EH_Sz_u
        complex*16,allocatable,dimension(:,:) :: EH_Sx_v,EH_Sy_v,EH_Sz_v
        complex*16,allocatable,dimension(:,:) :: E_Sx_u,E_Sy_u,E_Sz_u
        complex*16,allocatable,dimension(:,:) :: E_Sx_v,E_Sy_v,E_Sz_v
        complex*16,allocatable,dimension(:,:) :: H_Sx_u,H_Sy_u,H_Sz_u
        complex*16,allocatable,dimension(:,:) :: H_Sx_v,H_Sy_v,H_Sz_v
        complex*16,allocatable,dimension(:,:) :: conjg_EH_Sx_u,conjg_EH_Sy_u,conjg_EH_Sz_u
        complex*16,allocatable,dimension(:,:) :: conjg_EH_Sx_v,conjg_EH_Sy_v,conjg_EH_Sz_v
        complex*16,allocatable,dimension(:,:) :: conjg_H_Sx_u,conjg_H_Sy_u,conjg_H_Sz_u
        complex*16,allocatable,dimension(:,:) :: conjg_H_Sx_v,conjg_H_Sy_v,conjg_H_Sz_v

        complex*16,allocatable,dimension(:) :: EH1,EH2,EH3,EH4
        complex*16,allocatable,dimension(:) :: value_to_add
        integer :: current_BC,inner_region_id,outer_region_id
        complex*16 :: ar_in,ak_in,ar_ex,ak_ex,eta_ex,eta_in


        Nm = 3*N_calc
        N_out = (Nm*(Nm+1))/2
        allocate(Y_vec(N_out,0:N_order),b(Nm,0:N_order))
        allocate(EH_Sx_u(N_calc,0:N_order),EH_Sy_u(N_calc,0:N_order),EH_Sz_u(N_calc,0:N_order))
        allocate(EH_Sx_v(N_calc,0:N_order),EH_Sy_v(N_calc,0:N_order),EH_Sz_v(N_calc,0:N_order))
        allocate(E_Sx_u(N_calc,0:N_order),E_Sy_u(N_calc,0:N_order),E_Sz_u(N_calc,0:N_order))
        allocate(E_Sx_v(N_calc,0:N_order),E_Sy_v(N_calc,0:N_order),E_Sz_v(N_calc,0:N_order))
        allocate(H_Sx_u(N_calc,0:N_order),H_Sy_u(N_calc,0:N_order),H_Sz_u(N_calc,0:N_order))
        allocate(H_Sx_v(N_calc,0:N_order),H_Sy_v(N_calc,0:N_order),H_Sz_v(N_calc,0:N_order))
        allocate(conjg_EH_Sx_u(N_calc,0:N_order),conjg_EH_Sy_u(N_calc,0:N_order),conjg_EH_Sz_u(N_calc,0:N_order))
        allocate(conjg_EH_Sx_v(N_calc,0:N_order),conjg_EH_Sy_v(N_calc,0:N_order),conjg_EH_Sz_v(N_calc,0:N_order))
        allocate(conjg_H_Sx_u(N_calc,0:N_order),conjg_H_Sy_u(N_calc,0:N_order),conjg_H_Sz_u(N_calc,0:N_order))
        allocate(conjg_H_Sx_v(N_calc,0:N_order),conjg_H_Sy_v(N_calc,0:N_order),conjg_H_Sz_v(N_calc,0:N_order))
        allocate(EH1(0:N_order),EH2(0:N_order),EH4(0:N_order),EH3(0:N_order))
        allocate(E_Sx(0:N_order),H_Sx(0:N_order),E_Sy(0:N_order),H_Sy(0:N_order),E_Sz(0:N_order),H_Sz(0:N_order))
        allocate(value_to_add(0:N_order))
        b = 0.d0
        Y_vec = 0.d0

        do kk=1,this%M
            current_BC = this%testing_pt_status(kk,1)
            if(current_BC == -1) then
                cycle
                !!know that this point has been cancelled
            endif
            inner_region_id = this%testing_pt_status(kk,2)
            outer_region_id = this%testing_pt_status(kk,3)
            if(outer_region_id == 0) then
                ar_ex = ar
                ak_ex = ak
                eta_ex = eta1
            else
                ar_ex = Scatterers_Lib(outer_region_id)%ar_local
                ak_ex = Scatterers_Lib(outer_region_id)%ak_local
                eta_ex =Scatterers_Lib(outer_region_id)%eta_local
            endif
            ar_in = Scatterers_Lib(inner_region_id)%ar_local
            ak_in = Scatterers_Lib(inner_region_id)%ak_local
            eta_in = Scatterers_Lib(inner_region_id)%eta_local
            !            write(*,* ) kk,current_BC,inner_region_id,outer_region_id,ar_in,ak_in,ar_ex,ak_ex
            EH1 = this%EH_excit(kk,:)
            EH2 = this%EH_excit(kk+this%M,:)
            EH3 = this%EH_excit(kk+2*this%M,:)
            EH4 = this%EH_excit(kk+3*this%M,:)
            !            if(kk == 2) then
            !            write(*,*) 'Ex',kk,abs(EH1),abs(EH2),abs(EH3),abs(EH4)
            !            endif
            do j =1,N_calc  !! evaluate source contributions once
                jj = j + (this%N_curr-N_calc)
                if(this%active_region(jj)==outer_region_id) then !! inside sources
                    if(this%I_stat(jj) == 1) then
                        call eval_kernel_Electric_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,eta_ex,az,ar_ex,ak_ex,kk,jj)
                    else
                        call eval_kernel_Magnetic_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,eta_ex,az,ar_ex,ak_ex,kk,jj)
                    endif
                    E_Sx_u(j,:) = dot(this%tang_u(kk),E_Sx,N_order+1)/eta_ex
                    E_Sx_v(j,:) = dot(this%tang_v(kk),E_Sx,N_order+1)/eta_ex
                    E_Sy_u(j,:) = dot(this%tang_u(kk),E_Sy,N_order+1)/eta_ex
                    E_Sy_v(j,:) = dot(this%tang_v(kk),E_Sy,N_order+1)/eta_ex
                    E_Sz_u(j,:) = dot(this%tang_u(kk),E_Sz,N_order+1)/eta_ex
                    E_Sz_v(j,:) = dot(this%tang_v(kk),E_Sz,N_order+1)/eta_ex

                    H_Sx_u(j,:) = dot(this%tang_u(kk),H_Sx,N_order+1)
                    H_Sx_v(j,:) = dot(this%tang_v(kk),H_Sx,N_order+1)
                    H_Sy_u(j,:) = dot(this%tang_u(kk),H_Sy,N_order+1)
                    H_Sy_v(j,:) = dot(this%tang_v(kk),H_Sy,N_order+1)
                    H_Sz_u(j,:) = dot(this%tang_u(kk),H_Sz,N_order+1)
                    H_Sz_v(j,:) = dot(this%tang_v(kk),H_Sz,N_order+1)
                elseif(this%active_region(jj) == inner_region_id) then !! outside sources producing the dielectric fields
                    if(this%I_stat(jj) == 1) then
                        call eval_kernel_Electric_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        eta_in,az,ar_in,ak_in,kk,jj)
                    else
                        call eval_kernel_Magnetic_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        eta_in,az,ar_in,ak_in,kk,jj)
                    endif
                    E_Sx_u(j,:) = -dot(this%tang_u(kk),E_Sx,N_order+1)/eta_ex
                    E_Sx_v(j,:) = -dot(this%tang_v(kk),E_Sx,N_order+1)/eta_ex
                    E_Sy_u(j,:) = -dot(this%tang_u(kk),E_Sy,N_order+1)/eta_ex
                    E_Sy_v(j,:) = -dot(this%tang_v(kk),E_Sy,N_order+1)/eta_ex
                    E_Sz_u(j,:) = -dot(this%tang_u(kk),E_Sz,N_order+1)/eta_ex
                    E_Sz_v(j,:) = -dot(this%tang_v(kk),E_Sz,N_order+1)/eta_ex

                    H_Sx_u(j,:) = -dot(this%tang_u(kk),H_Sx,N_order+1)
                    H_Sx_v(j,:) = -dot(this%tang_v(kk),H_Sx,N_order+1)
                    H_Sy_u(j,:) = -dot(this%tang_u(kk),H_Sy,N_order+1)
                    H_Sy_v(j,:) = -dot(this%tang_v(kk),H_Sy,N_order+1)
                    H_Sz_u(j,:) = -dot(this%tang_u(kk),H_Sz,N_order+1)
                    H_Sz_v(j,:) = -dot(this%tang_v(kk),H_Sz,N_order+1)
                else !! neither outer or inner regions [EVERYTHING should be zero in this case]
                    E_Sx_u(j,:) = 0.d0
                    E_Sx_v(j,:) = 0.d0
                    E_Sy_u(j,:) = 0.d0
                    E_Sy_v(j,:) = 0.d0
                    E_Sz_u(j,:) = 0.d0
                    E_Sz_v(j,:) = 0.d0

                    H_Sx_u(j,:) = 0.d0
                    H_Sx_v(j,:) = 0.d0
                    H_Sy_u(j,:) = 0.d0
                    H_Sy_v(j,:) = 0.d0
                    H_Sz_u(j,:) = 0.d0
                    H_Sz_v(j,:) = 0.d0
                endif

            enddo

            if(current_BC == 1) then
                EH_Sx_u = E_Sx_u
                EH_Sx_v = E_Sx_v
                EH_Sy_u = E_Sy_u
                EH_Sy_v = E_Sy_v
                EH_Sz_u = E_Sz_u
                EH_Sz_v = E_Sz_v

                conjg_EH_Sx_u = conjg(EH_Sx_u)
                conjg_EH_Sx_v = conjg(EH_Sx_v)
                conjg_EH_Sy_u = conjg(EH_Sy_u)
                conjg_EH_Sy_v = conjg(EH_Sy_v)
                conjg_EH_Sz_u = conjg(EH_Sz_u)
                conjg_EH_Sz_v = conjg(EH_Sz_v)

                conjg_H_Sx_u = 0.d0
                conjg_H_Sx_v = 0.d0
                conjg_H_Sy_u = 0.d0
                conjg_H_Sy_v = 0.d0
                conjg_H_Sz_u = 0.d0
                conjg_H_Sz_v = 0.d0
                H_Sx_u = 0.d0
                H_Sx_v = 0.d0
                H_Sy_u = 0.d0
                H_Sy_v = 0.d0
                H_Sz_u = 0.d0
                H_Sz_v = 0.d0
            elseif(current_BC == 2) then
                EH_Sx_u = 0.d0
                EH_Sx_v = 0.d0
                EH_Sy_u = 0.d0
                EH_Sy_v = 0.d0
                EH_Sz_u = 0.d0
                EH_Sz_v = 0.d0

                conjg_EH_Sx_u = 0.d0
                conjg_EH_Sx_v = 0.d0
                conjg_EH_Sy_u = 0.d0
                conjg_EH_Sy_v = 0.d0
                conjg_EH_Sz_u = 0.d0
                conjg_EH_Sz_v = 0.d0

                conjg_H_Sx_u = conjg(H_Sx_u)
                conjg_H_Sx_v = conjg(H_Sx_v)
                conjg_H_Sy_u = conjg(H_Sy_u)
                conjg_H_Sy_v = conjg(H_Sy_v)
                conjg_H_Sz_u = conjg(H_Sz_u)
                conjg_H_Sz_v = conjg(H_Sz_v)
            elseif(current_BC == 3) then
                EH_Sx_u = this%norm_eta_u(kk)*(E_Sx_u + (this%eta_uu(kk)*H_Sx_v - this%eta_uv(kk)*H_Sx_u))
                EH_Sx_v = this%norm_eta_v(kk)*(E_Sx_v + (this%eta_vu(kk)*H_Sx_v - this%eta_vv(kk)*H_Sx_u))
                EH_Sy_u = this%norm_eta_u(kk)*(E_Sy_u + (this%eta_uu(kk)*H_Sy_v - this%eta_uv(kk)*H_Sy_u))
                EH_Sy_v = this%norm_eta_v(kk)*(E_Sy_v + (this%eta_vu(kk)*H_Sy_v - this%eta_vv(kk)*H_Sy_u))
                EH_Sz_u = this%norm_eta_u(kk)*(E_Sz_u + (this%eta_uu(kk)*H_Sz_v - this%eta_uv(kk)*H_Sz_u))
                EH_Sz_v = this%norm_eta_v(kk)*(E_Sz_v + (this%eta_vu(kk)*H_Sz_v - this%eta_vv(kk)*H_Sz_u))

                !                ZY_Sx_u(j,:) = this%norm_eta_u*(E_Sx_u(j,:) + (this%eta_uu(kk)*H_Sx_v(j,:) - this%eta_uv(kk)*H_Sx_u(j,:)))
                !                ZY_Sy_u(j,:) = this%norm_eta_u*(E_Sy_u(j,:) + (this%eta_uu(kk)*H_Sy_v(j,:) - this%eta_uv(kk)*H_Sy_u(j,:)))
                !                ZY_Sz_u(j,:) = this%norm_eta_u*(E_Sz_u(j,:) + (this%eta_uu(kk)*H_Sz_v(j,:) - this%eta_uv(kk)*H_Sz_u(j,:)))
                !
                !                ZY_Sx_v(j,:) = this%norm_eta_v*(E_Sx_v(j,:) + (this%eta_vu(kk)*H_Sx_v(j,:) - this%eta_vv(kk)*H_Sx_u(j,:)))
                !                ZY_Sy_v(j,:) = this%norm_eta_v*(E_Sy_v(j,:) + (this%eta_vu(kk)*H_Sy_v(j,:) - this%eta_vv(kk)*H_Sy_u(j,:)))
                !                ZY_Sz_v(j,:) = this%norm_eta_v*(E_Sz_v(j,:) + (this%eta_vu(kk)*H_Sz_v(j,:) - this%eta_vv(kk)*H_Sz_u(j,:)))


                conjg_EH_Sx_u = conjg(EH_Sx_u)
                conjg_EH_Sx_v = conjg(EH_Sx_v)
                conjg_EH_Sy_u = conjg(EH_Sy_u)
                conjg_EH_Sy_v = conjg(EH_Sy_v)
                conjg_EH_Sz_u = conjg(EH_Sz_u)
                conjg_EH_Sz_v = conjg(EH_Sz_v)

                conjg_H_Sx_u = 0.d0
                conjg_H_Sx_v = 0.d0
                conjg_H_Sy_u = 0.d0
                conjg_H_Sy_v = 0.d0
                conjg_H_Sz_u = 0.d0
                conjg_H_Sz_v = 0.d0
                H_Sx_u = 0.d0
                H_Sx_v = 0.d0
                H_Sy_u = 0.d0
                H_Sy_v = 0.d0
                H_Sz_u = 0.d0
                H_Sz_v = 0.d0
            else
                EH_Sx_u = E_Sx_u
                EH_Sx_v = E_Sx_v
                EH_Sy_u = E_Sy_u
                EH_Sy_v = E_Sy_v
                EH_Sz_u = E_Sz_u
                EH_Sz_v = E_Sz_v

                conjg_EH_Sx_u = conjg(EH_Sx_u)
                conjg_EH_Sx_v = conjg(EH_Sx_v)
                conjg_EH_Sy_u = conjg(EH_Sy_u)
                conjg_EH_Sy_v = conjg(EH_Sy_v)
                conjg_EH_Sz_u = conjg(EH_Sz_u)
                conjg_EH_Sz_v = conjg(EH_Sz_v)

                conjg_H_Sx_u = conjg(H_Sx_u)
                conjg_H_Sx_v = conjg(H_Sx_v)
                conjg_H_Sy_u = conjg(H_Sy_u)
                conjg_H_Sy_v = conjg(H_Sy_v)
                conjg_H_Sz_u = conjg(H_Sz_u)
                conjg_H_Sz_v = conjg(H_Sz_v)
            !                conjg_H_Sx_u = 0.d0
            !                conjg_H_Sx_v = 0.d0
            !                conjg_H_Sy_u = 0.d0
            !                conjg_H_Sy_v = 0.d0
            !                conjg_H_Sz_u = 0.d0
            !                conjg_H_Sz_v = 0.d0
            !                H_Sx_u = 0.d0
            !                H_Sx_v = 0.d0
            !                H_Sy_u = 0.d0
            !                H_Sy_v = 0.d0
            !                H_Sz_u = 0.d0
            !                H_Sz_v = 0.d0
            endif




            do j=1,N_calc
                value_to_add = 0.d0
                do nn = 0,N_order
                    do ii = 0,nn
                        value_to_add(nn) = value_to_add(nn) - dble(n_C_i_matrix(nn,ii))*(conjg_EH_Sx_u(j,nn-ii)*EH1(ii) +&
                        conjg_EH_Sx_v(j,nn-ii)*EH2(ii) + conjg_H_Sx_u(j,nn-ii)*EH3(ii) +&
                        conjg_H_Sx_v(j,nn-ii)*EH4(ii))
                    enddo
                enddo
                b(j,:)= b(j,:)+value_to_add
            enddo
            index_mat = 1 + N_calc
            do j=1,N_calc
                value_to_add = 0.d0
                do nn = 0,N_order
                    do ii = 0,nn
                        value_to_add(nn) = value_to_add(nn) - dble(n_C_i_matrix(nn,ii))*(conjg_EH_Sy_u(j,nn-ii)*EH1(ii) +&
                        conjg_EH_Sy_v(j,nn-ii)*EH2(ii) + conjg_H_Sy_u(j,nn-ii)*EH3(ii) +&
                        conjg_H_Sy_v(j,nn-ii)*EH4(ii))
                    enddo
                enddo
                b(index_mat,:) = b(index_mat,:) + value_to_add
                index_mat = index_mat +1
            enddo
            do j=1,N_calc
                value_to_add = 0.d0
                do nn = 0,N_order
                    do ii = 0,nn
                        value_to_add(nn) = value_to_add(nn) - dble(n_C_i_matrix(nn,ii))*(conjg_EH_Sz_u(j,nn-ii)*EH1(ii) +&
                        conjg_EH_Sz_v(j,nn-ii)*EH2(ii) + conjg_H_Sz_u(j,nn-ii)*EH3(ii) +&
                        conjg_H_Sz_v(j,nn-ii)*EH4(ii))
                    enddo
                enddo
                b(index_mat,:) = b(index_mat,:) + value_to_add

                index_mat = index_mat +1
            enddo

            do j=1,N_calc
                index_mat = mat2vec(Nm,j,j)
                do i=j,N_calc
                    value_to_add = 0.d0
                    do nn = 0,N_order
                        do ii = 0,nn
                            value_to_add(nn) = value_to_add(nn) + dble(n_C_i_matrix(nn,ii))*(&
                            conjg_EH_Sx_u(i,nn-ii)*EH_Sx_u(j,ii) + conjg_EH_Sx_v(i,nn-ii)*EH_Sx_v(j,ii) +&
                            conjg_H_Sx_u(i,nn-ii)*H_Sx_u(j,ii) + conjg_H_Sx_v(i,nn-ii)*H_Sx_v(j,ii))
                        enddo
                    enddo
                    Y_vec(index_mat,:)= Y_vec(index_mat,:) + value_to_add

                    index_mat = index_mat +1
                enddo
                do i = 1,N_calc
                    value_to_add = 0.d0
                    do nn = 0,N_order
                        do ii = 0,nn
                            value_to_add(nn) = value_to_add(nn) + dble(n_C_i_matrix(nn,ii))*(&
                            conjg_EH_Sy_u(i,nn-ii)*EH_Sx_u(j,ii) + conjg_EH_Sy_v(i,nn-ii)*EH_Sx_v(j,ii) +&
                            conjg_H_Sy_u(i,nn-ii)*H_Sx_u(j,ii) + conjg_H_Sy_v(i,nn-ii)*H_Sx_v(j,ii))
                        enddo
                    enddo
                    Y_vec(index_mat,:)= Y_vec(index_mat,:) + value_to_add
                    index_mat = index_mat + 1
                enddo

                do i = 1,N_calc
                    value_to_add = 0.d0
                    do nn = 0,N_order
                        do ii = 0,nn
                            value_to_add(nn) = value_to_add(nn) + dble(n_C_i_matrix(nn,ii))*(&
                            conjg_EH_Sz_u(i,nn-ii)*EH_Sx_u(j,ii) + conjg_EH_Sz_v(i,nn-ii)*EH_Sx_v(j,ii) +&
                            conjg_H_Sz_u(i,nn-ii)*H_Sx_u(j,ii) + conjg_H_Sz_v(i,nn-ii)*H_Sx_v(j,ii))
                        enddo
                    enddo
                    Y_vec(index_mat,:)= Y_vec(index_mat,:) + value_to_add
                    index_mat = index_mat + 1
                enddo


                index_mat = mat2vec(Nm,j+N_calc,j+N_calc)
                do i=j,N_calc
                    value_to_add = 0.d0
                    do nn = 0,N_order
                        do ii = 0,nn
                            value_to_add(nn) = value_to_add(nn) + dble(n_C_i_matrix(nn,ii))*(&
                            conjg_EH_Sy_u(i,nn-ii)*EH_Sy_u(j,ii) + conjg_EH_Sy_v(i,nn-ii)*EH_Sy_v(j,ii) +&
                            conjg_H_Sy_u(i,nn-ii)*H_Sy_u(j,ii) + conjg_H_Sy_v(i,nn-ii)*H_Sy_v(j,ii))
                        enddo
                    enddo
                    Y_vec(index_mat,:)= Y_vec(index_mat,:) + value_to_add
                    index_mat = index_mat +1
                enddo

                do i = 1,N_calc
                    value_to_add = 0.d0
                    do nn = 0,N_order
                        do ii = 0,nn
                            value_to_add(nn) = value_to_add(nn) + dble(n_C_i_matrix(nn,ii))*(&
                            conjg_EH_Sz_u(i,nn-ii)*EH_Sy_u(j,ii) + conjg_EH_Sz_v(i,nn-ii)*EH_Sy_v(j,ii) +&
                            conjg_H_Sz_u(i,nn-ii)*H_Sy_u(j,ii) + conjg_H_Sz_v(i,nn-ii)*H_Sy_v(j,ii))
                        enddo
                    enddo
                    Y_vec(index_mat,:)= Y_vec(index_mat,:) + value_to_add
                    index_mat = index_mat + 1
                enddo


                index_mat = mat2vec(Nm,j+2*N_calc,j+2*N_calc)
                do i=j,N_calc
                    value_to_add = 0.d0
                    do nn = 0,N_order
                        do ii = 0,nn
                            value_to_add(nn) = value_to_add(nn) + dble(n_C_i_matrix(nn,ii))*(&
                            conjg_EH_Sz_u(i,nn-ii)*EH_Sz_u(j,ii) + conjg_EH_Sz_v(i,nn-ii)*EH_Sz_v(j,ii) +&
                            conjg_H_Sz_u(i,nn-ii)*H_Sz_u(j,ii) + conjg_H_Sz_v(i,nn-ii)*H_Sz_v(j,ii))
                        enddo
                    enddo
                    Y_vec(index_mat,:)= Y_vec(index_mat,:) + value_to_add
                    index_mat = index_mat +1
                enddo
            enddo
        enddo

        deallocate(E_Sx,H_Sx,E_Sy,H_Sy,E_Sz,H_Sz)
        deallocate(EH_Sx_u,EH_Sy_u,EH_Sz_u,EH_Sx_v,EH_Sy_v,EH_Sz_v)
        deallocate(E_Sx_u,E_Sy_u,E_Sz_u,E_Sx_v,E_Sy_v,E_Sz_v)
        deallocate(H_Sx_u,H_Sy_u,H_Sz_u,H_Sx_v,H_Sy_v,H_Sz_v)
        deallocate(conjg_EH_Sx_u,conjg_EH_Sy_u,conjg_EH_Sz_u,conjg_EH_Sx_v,conjg_EH_Sy_v,conjg_EH_Sz_v)
        deallocate(conjg_H_Sx_u,conjg_H_Sy_u,conjg_H_Sz_u,conjg_H_Sx_v,conjg_H_Sy_v,conjg_H_Sz_v)
        deallocate(EH1,EH2,EH3,EH4)
        deallocate(value_to_add)
    end subroutine get_LSM_matrices








    subroutine eval_near_field_Magnetic_2D_r(N,E_Mx,E_My,E_Mz,H_Mx,H_My,H_Mz,ak_l,az_l,ar_l,eta_l,x,y,source_n,pos)
        !        type(Scatterer) :: this
        !! a function that evaluates the kernel of radiation
        real(8),allocatable,intent(inout) :: pos(:,:)
        type(vector_c),allocatable,intent(inout) :: E_Mx(:),E_My(:),E_Mz(:),H_Mx(:),H_My(:),H_Mz(:)
        integer :: source_n
        real(8) :: R,x,y,xp,yp,xd,yd,az_l
        integer :: N
        complex*16,dimension(0:N) :: k0H0,k0H1,k0H2,H1
        real(8) :: ak_l,arR,ar_l,eta_l
        complex*16 :: const1
        complex*16 :: const2
        complex*16,dimension(0:N) :: const3,const4,CE_Mz,CH_Mz,const5,const6
        complex*16,dimension(0:N) :: c_zero
        c_zero = (0.d0,0.d0)

        xp = pos(source_n,1)
        yp = pos(source_n,2)

        xd = x-xp
        yd = y-yp
        R = sqrt(xd**2.d0+yd**2.d0)


        arR = ar_l*R


        call Hankel_derivatives(arR,k0,N,k0H0,k0H1,k0H2,H1)


        !        CE_Mz = kr*eta_l/(4.d0*cj)*h1/R

        CE_Mz = ar_l*eta_l/((0.d0,4.d0))*k0H1/R
        E_Mz = vec(CE_Mz*yd,-CE_Mz*xd,c_zero,N+1)



        !        CH_Mz = cj*kz*kr/(4*k_l*R)*h1
        CH_Mz = -az_l*ar_l/((0.d0,4.d0)*ak_l*R)*k0H1
        H_Mz = vec(CH_Mz*xd,CH_Mz*yd,&
        -ar_l*ar_l/((4.d0,0.d0)*ak_l)*k0H0,N+1)
        !        -kr**2.d0/(4*k_l)*h0)


        const3 = az_l*eta_l/(4.d0,0.d0)*k0H0
        const4 = ar_l*eta_l/((0.d0,4.d0)*R)*k0H1

        E_Mx = vec(c_zero,&
        -const3,&
        -const4*yd,N+1)

        E_My = vec(const3,&
        c_zero,&
        const4*xd,N+1)


        !        E_My = vec(kz*eta_l/4.d0*h0,&
        !        c_zero,&
        !        kr*eta_l/(4.d0*cj*R)*h1*xd)


        !        CE_Iy = kr**2.d0*eta_l/(4.d0*k_l)

        !        H_Mx = vec(-k_l/4.d0*h0-CE_Iy/eta_l*(h2*(xd/R)**2.d0 - h1/kRR),&
        !        -CE_Iy/eta_l*xd*yd*h2/R**2.d0 ,&
        !        cj*kr*kz/(4.d0*k_l)*h1/R*xd)

!        CE_Iy = ar_l*ar_l/((4.d0,0.d0)*ak_l)
        const1 = ar_l*ar_l/((4.d0,0.d0)*ak_l)
        const2 = ar_l*az_l/((0.d0,-4.d0)*ak_l)
        const5 = -const1*xd*yd*k0H2/R**2.d0
        const4 = -ak_l/(4.d0,0.d0)*k0H0
        const3 = const2*k0H1/R
        const6 = H1/arR
        H_Mx = vec(const4-const1*(k0H2*(xd/R)**2.d0 - const6),&
        const5 ,&
        const3*xd,N+1)

        H_My = vec(const5 ,&
        const4-const1*(k0H2*(yd/R)**2.d0 - const6),&
        const3*yd ,N+1)

    end subroutine eval_near_field_Magnetic_2D_r

    subroutine eval_near_field_Electric_2D_r(N,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,ak_l,az_l,ar_l,eta_l,x,y,source_n,pos)
        !        type(Scatterer) :: this
        !! a function that evaluates the kernel of radiation
        type(vector_c),allocatable,intent(inout)  :: E_Iz(:),E_Ix(:),E_Iy(:),H_Iz(:),H_Ix(:),H_Iy(:)
        real(8),allocatable,intent(inout) :: pos(:,:)
        real(8) :: eta_l !! the wave-number and characteristic impedance inside the medium
        integer :: source_n
        real(8) :: R,x,y,xp,yp,xd,yd,az_l

        integer :: N
        complex*16,dimension(0:N) :: k0H0,k0H1,k0H2,h1
        real(8) :: ak_l,arR,ar_l
        complex*16,dimension(0:N) :: const1,const2,const3,const4,const5,const6,c_zero,CE_Iy,CE_Iz,CH_Iz
        c_zero = (0.d0,0.d0)


        xp = pos(source_n,1)
        yp = pos(source_n,2)
        xd = x-xp
        yd = y-yp

        R = sqrt(xd**2.d0+yd**2.d0)

        arR = ar_l*R

        !        kRR = kr*R
        !
        !        h0 = besselh2_0(kRR)
        !        h1 = besselh2_1(kRR)
        !        h2 = 2.d0/kRR*h1-h0
        !
        !        k0H0 = k0*h0
        !        k0H1 = k0*h1
        !        k0H2 = k0*h2
        call Hankel_derivatives(arR,k0,N,k0H0,k0H1,k0H2,H1)

        CE_Iy = ar_l*ar_l*eta_l/((4.d0,0.d0)*ak_l)
        const3 = -ak_l*eta_l/(4.d0,0.d0)*k0H0
        const4 = -CE_Iy*xd*yd*k0H2/R**2.d0
        const5 = -ar_l*az_l*eta_l/((0.d0,4.d0)*ak_l*R)*k0H1
        const6 = H1/arR

        E_Iy = vec(const4 ,&
        const3-CE_Iy*(k0H2*(yd/R)**2.d0 - const6),&
        const5*yd ,N+1)

        E_Ix = vec(const3-CE_Iy*(k0h2*(xd/R)**2.d0 - const6),&
        const4 ,&
        const5*xd,N+1)



        CE_Iz = -az_l*ar_l*k0H1*eta_l/((0.d0,4.d0)*ak_l*R)
        E_Iz = vec(CE_Iz*xd,CE_Iz*yd,&
        -ar_l*ar_l*eta_l/((4.d0,0.d0)*ak_l)*k0H0,N+1)

        CH_Iz = ar_l/(0.d0,4.d0)*k0H1/R
        H_Iz = vec(-CH_Iz*yd,CH_Iz*xd,c_zero,N+1)

        const1 = az_l/(4.d0,0.d0)*k0H0
        const2 = ar_l/((0.d0,4.d0)*R)*k0H1

        H_Ix = vec(c_zero,&
        const1,&
        const2*yd,N+1)

        H_Iy = vec(-const1, &
        c_zero,&
        -const2*xd,N+1)
    end subroutine eval_near_field_Electric_2D_r

    subroutine eval_near_field_Magnetic_2D_c(N,E_Mx,E_My,E_Mz,H_Mx,H_My,H_Mz,ak_l,az_l,ar_l,eta_l,x,y,source_n,pos)
        !        type(Scatterer) :: this
        !! a function that evaluates the kernel of radiation
        real(8),allocatable,intent(inout) :: pos(:,:)
        type(vector_c),allocatable,intent(inout) :: E_Mx(:),E_My(:),E_Mz(:),H_Mx(:),H_My(:),H_Mz(:)
        integer :: source_n
        real(8) :: R,x,y,xp,yp,xd,yd,az_l
        integer :: N
        complex*16,dimension(0:N) :: k0H0,k0H1,k0H2,H1
        complex*16 :: ak_l,arR,ar_l,eta_l
        complex*16 :: const1
        complex*16 :: const2
        complex*16,dimension(0:N) :: const3,const4,CE_Mz,CH_Mz,const5,const6
        complex*16,dimension(0:N) :: c_zero
        c_zero = (0.d0,0.d0)

        xp = pos(source_n,1)
        yp = pos(source_n,2)

        xd = x-xp
        yd = y-yp
        R = sqrt(xd**2.d0+yd**2.d0)


        arR = ar_l*R


        call Hankel_derivatives(arR,k0,N,k0H0,k0H1,k0H2,H1)


        !        CE_Mz = kr*eta_l/(4.d0*cj)*h1/R

        CE_Mz = ar_l*eta_l/((0.d0,4.d0))*k0H1/R
        E_Mz = vec(CE_Mz*yd,-CE_Mz*xd,c_zero,N+1)



        !        CH_Mz = cj*kz*kr/(4*k_l*R)*h1
        CH_Mz = -az_l*ar_l/((0.d0,4.d0)*ak_l*R)*k0H1
        H_Mz = vec(CH_Mz*xd,CH_Mz*yd,&
        -ar_l*ar_l/((4.d0,0.d0)*ak_l)*k0H0,N+1)
        !        -kr**2.d0/(4*k_l)*h0)


        const3 = az_l*eta_l/(4.d0,0.d0)*k0H0
        const4 = ar_l*eta_l/((0.d0,4.d0)*R)*k0H1

        E_Mx = vec(c_zero,&
        -const3,&
        -const4*yd,N+1)

        E_My = vec(const3,&
        c_zero,&
        const4*xd,N+1)


        !        E_My = vec(kz*eta_l/4.d0*h0,&
        !        c_zero,&
        !        kr*eta_l/(4.d0*cj*R)*h1*xd)


        !        CE_Iy = kr**2.d0*eta_l/(4.d0*k_l)

        !        H_Mx = vec(-k_l/4.d0*h0-CE_Iy/eta_l*(h2*(xd/R)**2.d0 - h1/kRR),&
        !        -CE_Iy/eta_l*xd*yd*h2/R**2.d0 ,&
        !        cj*kr*kz/(4.d0*k_l)*h1/R*xd)

!        CE_Iy = ar_l*ar_l/((4.d0,0.d0)*ak_l)
        const1 = ar_l*ar_l/((4.d0,0.d0)*ak_l)
        const2 = ar_l*az_l/((0.d0,-4.d0)*ak_l)
        const5 = -const1*xd*yd*k0H2/R**2.d0
        const4 = -ak_l/(4.d0,0.d0)*k0H0
        const3 = const2*k0H1/R
        const6 = H1/arR
        H_Mx = vec(const4-const1*(k0H2*(xd/R)**2.d0 - const6),&
        const5 ,&
        const3*xd,N+1)

        H_My = vec(const5 ,&
        const4-const1*(k0H2*(yd/R)**2.d0 - const6),&
        const3*yd ,N+1)

    end subroutine eval_near_field_Magnetic_2D_c

    subroutine eval_near_field_Electric_2D_c(N,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,ak_l,az_l,ar_l,eta_l,x,y,source_n,pos)
        !        type(Scatterer) :: this
        !! a function that evaluates the kernel of radiation
        type(vector_c),allocatable,intent(inout)  :: E_Iz(:),E_Ix(:),E_Iy(:),H_Iz(:),H_Ix(:),H_Iy(:)
        real(8),allocatable,intent(inout) :: pos(:,:)
        complex*16 :: eta_l !! the wave-number and characteristic impedance inside the medium
        integer :: source_n
        real(8) :: R,x,y,xp,yp,xd,yd,az_l

        integer :: N
        complex*16,dimension(0:N) :: k0H0,k0H1,k0H2,h1
        complex*16 :: ak_l,arR,ar_l
        complex*16,dimension(0:N) :: const1,const2,const3,const4,const5,const6,c_zero,CE_Iy,CE_Iz,CH_Iz
        c_zero = (0.d0,0.d0)


        xp = pos(source_n,1)
        yp = pos(source_n,2)
        xd = x-xp
        yd = y-yp

        R = sqrt(xd**2.d0+yd**2.d0)

        arR = ar_l*R

        !        kRR = kr*R
        !
        !        h0 = besselh2_0(kRR)
        !        h1 = besselh2_1(kRR)
        !        h2 = 2.d0/kRR*h1-h0
        !
        !        k0H0 = k0*h0
        !        k0H1 = k0*h1
        !        k0H2 = k0*h2
        call Hankel_derivatives(arR,k0,N,k0H0,k0H1,k0H2,H1)

        CE_Iy = ar_l*ar_l*eta_l/((4.d0,0.d0)*ak_l)
        const3 = -ak_l*eta_l/(4.d0,0.d0)*k0H0
        const4 = -CE_Iy*xd*yd*k0H2/R**2.d0
        const5 = -ar_l*az_l*eta_l/((0.d0,4.d0)*ak_l*R)*k0H1
        const6 = H1/arR

        E_Iy = vec(const4 ,&
        const3-CE_Iy*(k0H2*(yd/R)**2.d0 - const6),&
        const5*yd ,N+1)

        E_Ix = vec(const3-CE_Iy*(k0h2*(xd/R)**2.d0 - const6),&
        const4 ,&
        const5*xd,N+1)



        CE_Iz = -az_l*ar_l*k0H1*eta_l/((0.d0,4.d0)*ak_l*R)
        E_Iz = vec(CE_Iz*xd,CE_Iz*yd,&
        -ar_l*ar_l*eta_l/((4.d0,0.d0)*ak_l)*k0H0,N+1)

        CH_Iz = ar_l/(0.d0,4.d0)*k0H1/R
        H_Iz = vec(-CH_Iz*yd,CH_Iz*xd,c_zero,N+1)

        const1 = az_l/(4.d0,0.d0)*k0H0
        const2 = ar_l/((0.d0,4.d0)*R)*k0H1

        H_Ix = vec(c_zero,&
        const1,&
        const2*yd,N+1)

        H_Iy = vec(-const1, &
        c_zero,&
        -const2*xd,N+1)
    end subroutine eval_near_field_Electric_2D_c

    ! -----------------------------------------------------------------------
    ! Subroutine: eval_bouncing_field
    ! Purpose   : Evaluates the scattered field from scatterer 'this' at the
    !             testing points of scatterer 'Scatterer_2' and accumulates
    !             it into Scatterer_2%E_bounce_current / H_bounce_current.
    !             Called during IFB iterations to propagate inter-scatterer
    !             interactions. All N sources currently active in 'this' are
    !             used.
    ! -----------------------------------------------------------------------
    subroutine eval_bouncing_field(this,N_order,Scatterer_2,Scatterers_Lib)
        !! evaluates the bouncing fields on Scatterer(Scatterer_2) due to the total number of sources inside Scatterer(this)
        type(Scatterer),intent(inout) :: this,Scatterer_2
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        integer :: j,q,N_order,nn
        type(vector_c),allocatable,dimension(:) ::E_tot,H_tot
        type(vector_c),allocatable,dimension(:) :: E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz
        real(8) :: x,y
        integer :: current_BC,inner_region_id,outer_region_id
        complex*16 :: ar_ex,ak_ex,ar_in,ak_in,eta_ex,eta_in

        allocate(E_Sx(0:N_order),H_Sx(0:N_order),E_Sy(0:N_order),H_Sy(0:N_order),E_Sz(0:N_order),H_Sz(0:N_order))
        allocate(E_tot(0:N_order),H_tot(0:N_order))

        do j=1,Scatterer_2%M !! loop over the second scatterers testing points
            current_BC = Scatterer_2%testing_pt_status(j,1)
            if(current_BC == -1) then
                cycle !! do not calculate scattered field for this cancelled point
            endif
            inner_region_id = Scatterer_2%testing_pt_status(j,2)
            outer_region_id = Scatterer_2%testing_pt_status(j,3)

            if(outer_region_id == 0) then
                ar_ex = ar
                ak_ex = ak
                eta_ex = eta1
            else
                ar_ex = Scatterers_Lib(outer_region_id)%ar_local
                ak_ex = Scatterers_Lib(outer_region_id)%ak_local
                eta_ex = Scatterers_Lib(outer_region_id)%eta_local
            endif
            ar_in = Scatterers_Lib(inner_region_id)%ar_local
            ak_in = Scatterers_Lib(inner_region_id)%ak_local
            eta_in = Scatterers_Lib(inner_region_id)%eta_local


            x = Scatterer_2%testing_pt(j,1)
            y = Scatterer_2%testing_pt(j,2)
            do nn = 0,N_order
                E_tot(nn) = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
                H_tot(nn) = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
            enddo

            do q = 1,this%N
                if(this%active_region(q) == outer_region_id) then !! the sources are inside
                    if(this%I_stat(q) == 1) then
                        !                    call eval_kernel_Electric_2D(this,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,k1,eta1,j,qq)
                        call eval_near_field_Electric_2D(N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        ak_ex,az,ar_ex,eta_ex,x,y,q,this%Source_pos)
                    else
                        !                    call eval_kernel_Magnetic_2D(this,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,k1,eta1,j,qq)
                        call eval_near_field_Magnetic_2D(N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        ak_ex,az,ar_ex,eta_ex,x,y,q,this%Source_pos)
                    endif
                    !                    write(*,*) (1/eta1)*E_Sx
                    !                    E_tot = E_tot + (this%I_sources_total_current(q,1)*E_Sx + &
                    !                    this%I_sources_total_current(q,2)*E_Sy +this%I_sources_total_current(q,3)*E_Sz)
                    !                    H_tot = H_tot + (this%I_sources_total_current(q,1)*H_Sx + &
                    !                    this%I_sources_total_current(q,2)*H_Sy +this%I_sources_total_current(q,3)*H_Sz)

                    do nn = 0,N_order
                        !                    write(*,*) this%I_sources_total_current(q,nn)%v(1)
                        E_tot(nn) = E_tot(nn) + (this%I_sources_total_current(q,nn)%v(1)*E_Sx(nn) + &
                        this%I_sources_total_current(q,nn)%v(2)*E_Sy(nn) +&
                        this%I_sources_total_current(q,nn)%v(3)*E_Sz(nn))
                        H_tot(nn) = H_tot(nn) + (this%I_sources_total_current(q,nn)%v(1)*H_Sx(nn) + &
                        this%I_sources_total_current(q,nn)%v(2)*H_Sy(nn) +this%I_sources_total_current(q,nn)%v(3)*H_Sz(nn))
                    enddo
                elseif(this%active_region(q) == inner_region_id) then
                    write(*,*) 'Warning: multi-region scattering is not implemented in routine eval_bouncing_field'
                endif
            enddo
            Scatterer_2%E_bounce_current(j,:) = Scatterer_2%E_bounce_current(j,:) + &
            dot(Scatterer_2%tang_u(j),E_tot,N_order+1)
            Scatterer_2%E_bounce_current(j+Scatterer_2%M,:) = Scatterer_2%E_bounce_current(j+Scatterer_2%M,:) +&
            dot(Scatterer_2%tang_v(j),E_tot,N_order+1)
            Scatterer_2%H_bounce_current(j,:) = Scatterer_2%H_bounce_current(j,:) + &
            dot(Scatterer_2%tang_u(j),H_tot,N_order+1)
            Scatterer_2%H_bounce_current(j+Scatterer_2%M,:) = Scatterer_2%H_bounce_current(j+Scatterer_2%M,:) + &
            dot(Scatterer_2%tang_v(j),H_tot,N_order+1)
        enddo
        !        write(*,*) 'Source Scatterer',this%region_ID-1,'N sources =',this%N
        !        write(*,*) 'E_bounce MAG=',sum(abs(Scatterer_2%E_bounce))/sqrt(Scatterer_2%normalize_E)
        deallocate(E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz)
        deallocate(E_tot,H_tot)
    end subroutine eval_bouncing_field


    subroutine eval_bouncing_field_specified(this,N_order,Scatterer_2,N_calc,Scatterers_Lib)
        !! evaluates the bouncing fields on Scatterer(Scatterer_2) due to the total number of sources inside Scatterer(this)
        !! due to the last N_calc sources (excluding the sources for the internal problem solution)
        type(Scatterer) :: this,Scatterer_2
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        integer :: j,q,N_calc,qq,N_order,nn
        type(vector_c),allocatable,dimension(:) :: E_tot,H_tot
        type(vector_c),allocatable,dimension(:) :: E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz
        real(8) :: x,y
        integer :: current_BC,inner_region_id,outer_region_id
        complex*16 :: ar_ex,ak_ex,ar_in,ak_in,eta_ex,eta_in

        allocate(E_Sx(0:N_order),H_Sx(0:N_order),E_Sy(0:N_order),H_Sy(0:N_order),E_Sz(0:N_order),H_Sz(0:N_order))
        allocate(E_tot(0:N_order),H_tot(0:N_order))

        do j=1,Scatterer_2%M !! loop over the second scatterers testing points
            current_BC = Scatterer_2%testing_pt_status(j,1)
            if(current_BC == -1) then
                cycle !! do not calculate scattered field for this cancelled point
            endif
            inner_region_id = Scatterer_2%testing_pt_status(j,2)
            outer_region_id = Scatterer_2%testing_pt_status(j,3)

            if(outer_region_id == 0) then
                ar_ex = ar
                ak_ex = ak
                eta_ex = eta1
            else
                ar_ex = Scatterers_Lib(outer_region_id)%ar_local
                ak_ex = Scatterers_Lib(outer_region_id)%ak_local
                eta_ex = Scatterers_Lib(outer_region_id)%eta_local
            endif
            ar_in = Scatterers_Lib(inner_region_id)%ar_local
            ak_in = Scatterers_Lib(inner_region_id)%ak_local
            eta_in = Scatterers_Lib(inner_region_id)%eta_local
            x = Scatterer_2%testing_pt(j,1)
            y = Scatterer_2%testing_pt(j,2)
            do nn = 0,N_order
                E_tot(nn) = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
                H_tot(nn) = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
            enddo
            do q = 1,N_calc
                qq = q + (this%N_curr- N_calc)
                if(this%active_region(qq) == outer_region_id) then !! the sources are inside
                    if(this%I_stat(qq) == 1) then
                        call eval_near_field_Electric_2D(N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        ak_ex,az,ar_ex,eta_ex,x,y,q,this%Source_pos)
                    else
                        call eval_near_field_Magnetic_2D(N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        ak_ex,az,ar_ex,eta_ex,x,y,q,this%Source_pos)
                    endif
                    !                    write(*,*) this%I_sources(qq,1)
                    !                    E_tot = E_tot + (this%I_sources(qq,1)*E_Sx + this%I_sources(qq,2)*E_Sy +this%I_sources(qq,3)*E_Sz)
                    !                    H_tot = H_tot + (this%I_sources(qq,1)*H_Sx + this%I_sources(qq,2)*H_Sy +this%I_sources(qq,3)*H_Sz)
                    do nn = 0,N_order
                        !                        write(*,*) this%I_sources(qq,nn)%v(1)
                        E_tot(nn) = E_tot(nn) + (this%I_sources(qq,nn)%v(1)*E_Sx(nn) + &
                        this%I_sources(qq,nn)%v(2)*E_Sy(nn) +this%I_sources(qq,nn)%v(3)*E_Sz(nn))
                        H_tot(nn) = H_tot(nn) + (this%I_sources(qq,nn)%v(1)*H_Sx(nn) + &
                        this%I_sources(qq,nn)%v(2)*H_Sy(nn) +this%I_sources(qq,nn)%v(3)*H_Sz(nn))
                    !                        write(*,*) E_tot(nn),H_tot(nn)
                    enddo
                elseif(this%active_region(qq) == inner_region_id) then !! the sources are outside (activated only in the dieletric case)


                endif
            enddo
            Scatterer_2%E_bounce_current(j,:) = Scatterer_2%E_bounce_current(j,:) + dot(Scatterer_2%tang_u(j),E_tot,N_order+1)
            Scatterer_2%E_bounce_current(j+Scatterer_2%M,:) = Scatterer_2%E_bounce_current(j+Scatterer_2%M,:) +&
            dot(Scatterer_2%tang_v(j),E_tot,N_order+1)
            Scatterer_2%H_bounce_current(j,:) = Scatterer_2%H_bounce_current(j,:) + dot(Scatterer_2%tang_u(j),H_tot,N_order+1)
            Scatterer_2%H_bounce_current(j+Scatterer_2%M,:) = Scatterer_2%H_bounce_current(j+Scatterer_2%M,:) + &
            dot(Scatterer_2%tang_v(j),H_tot,N_order+1)

        !            write(*,*) Scatterer_2%H_bounce_current(j,:)
        enddo
        !        write(*,*) 'bouncing magnitude',sum(abs(Scatterer_2%E_bounce))/sqrt(Scatterer_2%normalize_E)
        deallocate(E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz)
        deallocate(E_tot,H_tot)
    end subroutine eval_bouncing_field_specified

    subroutine update_scattered_fields(this,N_order,N_calc,Scatterers_Lib)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        integer :: N_calc,q,qq,j,ii
        integer :: N_order,nn
        type(vector_c),allocatable,dimension(:) :: E_tot,H_tot,Ed_tot,Hd_tot
        type(vector_c),allocatable,dimension(:) :: E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz
        integer :: current_BC,inner_region_id,outer_region_id
        complex*16 :: ar_ex,ak_ex,ar_in,ak_in,eta_ex,eta_in

        allocate(E_Sx(0:N_order),E_Sy(0:N_order),E_Sz(0:N_order),H_Sx(0:N_order),H_Sy(0:N_order),H_Sz(0:N_order))
        allocate(E_tot(0:N_order),H_tot(0:N_order),Ed_tot(0:N_order),Hd_tot(0:N_order))
        this%Es_c = 0.d0
        this%Hs_c = 0.d0
        this%Ed_c = 0.d0
        this%Hd_c = 0.d0



        do j=1,this%M
            current_BC = this%testing_pt_status(j,1)
            if(current_BC == -1) then
                cycle !! do not calculate scattered field for this cancelled point
            endif
            do nn = 0,N_order
                H_tot(nn)%v = 0.d0
                Ed_tot(nn)%v = 0.d0
                Hd_tot(nn)%v = 0.d0
                E_tot(nn)%v = 0.d0
            enddo

            inner_region_id = this%testing_pt_status(j,2)
            outer_region_id = this%testing_pt_status(j,3)

            if(outer_region_id == 0) then
                ar_ex = ar
                ak_ex = ak
                eta_ex = eta1
            else
                ar_ex = Scatterers_Lib(outer_region_id)%ar_local
                ak_ex = Scatterers_Lib(outer_region_id)%ak_local
                eta_ex = Scatterers_Lib(outer_region_id)%eta_local
            endif
            ar_in = Scatterers_Lib(inner_region_id)%ar_local
            ak_in = Scatterers_Lib(inner_region_id)%ak_local
            eta_in = Scatterers_Lib(inner_region_id)%eta_local

            do q = 1,N_calc
                qq = q + (this%N_curr- N_calc)

                if(this%active_region(qq) == outer_region_id) then !! the sources are inside
                    !                    write(*,*) 'background'
                    if(this%I_stat(qq) == 1) then
                        call eval_kernel_Electric_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,eta_ex,az,ar_ex,ak_ex,j,qq)
                    else
                        call eval_kernel_Magnetic_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,eta_ex,az,ar_ex,ak_ex,j,qq)
                    endif

                    do nn = 0,N_order
                        do ii = 0,nn
                            E_tot(nn) = E_tot(nn) + dble(n_P_i_matrix(nn,ii))*(this%I_sources(qq,ii)%v(1)*E_Sx(nn-ii)+&
                            this%I_sources(qq,ii)%v(2)*E_Sy(nn-ii)+this%I_sources(qq,ii)%v(3)*E_Sz(nn-ii))
                            H_tot(nn) = H_tot(nn) + dble(n_P_i_matrix(nn,ii))*(this%I_sources(qq,ii)%v(1)*H_Sx(nn-ii)+&
                            this%I_sources(qq,ii)%v(2)*H_Sy(nn-ii)+this%I_sources(qq,ii)%v(3)*H_Sz(nn-ii))
                        enddo
                    enddo
                elseif(this%active_region(qq) == inner_region_id) then !! the sources are outside (activated only in the dieletric case)
                    !                    write(*,*) 'dielectric'
                    if(this%I_stat(qq) == 1) then
                        call eval_kernel_Electric_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        eta_in,az,ar_in,ak_in,j,qq)
                    else
                        call eval_kernel_Magnetic_2D(this,this,N_order,E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz,&
                        eta_in,az,ar_in,ak_in,j,qq)
                    endif

                    do nn = 0,N_order
                        do ii = 0,nn
                            Ed_tot(nn) = Ed_tot(nn) + dble(n_P_i_matrix(nn,ii))*(this%I_sources(qq,ii)%v(1)*E_Sx(nn-ii)+&
                            this%I_sources(qq,ii)%v(2)*E_Sy(nn-ii)+this%I_sources(qq,ii)%v(3)*E_Sz(nn-ii))
                            Hd_tot(nn) = Hd_tot(nn) + dble(n_P_i_matrix(nn,ii))*(this%I_sources(qq,ii)%v(1)*H_Sx(nn-ii)+&
                            this%I_sources(qq,ii)%v(2)*H_Sy(nn-ii)+this%I_sources(qq,ii)%v(3)*H_Sz(nn-ii))
                        enddo
                    enddo
                else
                    cycle
                endif

            enddo
            this%Es_c(j,:)   =  dot(this%tang_u(j),E_tot,N_order+1)
            this%Es_c(j+this%M,:) =  dot(this%tang_v(j),E_tot,N_order+1)
            this%Hs_c(j,:)   =  dot(this%tang_u(j),H_tot,N_order+1)
            this%Hs_c(j+this%M,:) =  dot(this%tang_v(j),H_tot,N_order+1)

            !            if(this%Problem_type == 4) then
            !            write(*,*) 'Es_tot',E_tot,'Ed_tot',Ed_tot
            this%Ed_c(j,:) =   dot(this%tang_u(j),Ed_tot,N_order+1)
            this%Ed_c(j+this%M,:) =  dot(this%tang_v(j),Ed_tot,N_order+1)
            this%Hd_c(j,:)   =  dot(this%tang_u(j),Hd_tot,N_order+1)
            this%Hd_c(j+this%M,:) = dot(this%tang_v(j),Hd_tot,N_order+1)
        !              endif

        enddo
        deallocate(E_Sx,E_Sy,E_Sz,H_Sx,H_Sy,H_Sz)
        deallocate(E_tot,H_tot,Ed_tot,Hd_tot)
    end subroutine update_scattered_fields


    ! -----------------------------------------------------------------------
    ! Subroutine: eval_matrices_RAS
    ! Purpose   : Builds or extends the RAS linear system for the current
    !             group of N_calc new sources.
    !   If use_old=.false.: assembles new Y and b from scratch for the current
    !     group, factorises, and caches the result in Matrix_Storage.
    !   If use_old=.true. : retrieves the cached factorisation and reuses it.
    !   Selects the matrix solver based on Matrix_Solution_Method:
    !     1 = solve_LU_Cholesky, 2 = solve_Matrix_preconditioned,
    !     3 = solve_Matrix_SVD, 4 = solve_Matrix_SVD_preconditioned.
    !   Stores solution currents in I_sources(:, 0:N_order).
    ! -----------------------------------------------------------------------
    subroutine eval_matrices_RAS(this,N_order,N_calc,use_old,Scatterers_pointer)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)
        integer :: N_order,nn,ii
        integer :: N_calc,p,q,N_out,Nm,Lm_start_index,Dm_start_index,Nm1,N_out1
        complex*16,allocatable,dimension(:) ::Dm,Lm,Dm_temp,Lm_temp,y_n,Y_0_vec,b_0_vec
        complex*16,allocatable,dimension(:,:) :: Y_collect_temp
        complex*16,allocatable,dimension(:,:) :: I_M_n,Y,b,Yinv,Y_mat
        integer,allocatable,dimension(:) :: Nm_collect_temp,rm_list
        real(8),dimension(:),allocatable :: Y_mat_norm
        integer,dimension(:),allocatable :: IWORK
        logical :: use_old
        !        real :: start_t,finish_t


        Nm = 3*N_calc

        !        call cpu_time(start_t)
        if(.not. use_old) then
            !            write(*,*) 'before LSM' ,sum(abs(this%EH_excit(1:this%M,0)))
            call get_LSM_matrices(this,N_order,N_calc,Y,b,N_out,Scatterers_pointer)

            if(Matrix_Solution_Method == 3 .or. Matrix_Solution_Method == 4 .or. Matrix_Solution_Method == 5) then
                allocate(IWORK(Nm),Y_mat_norm(Nm))
            endif
        else
            call get_LSM_excitation(this,N_order,N_calc,b,Scatterers_pointer)
            !            if(this%Problem_Type == 1) then !! the conventional electric field matching for PEC
            !                call get_LSM_excitation_vector_PEC(this,N_order,N_calc,b)
            !            elseif(this%Problem_Type == 2) then !! PMC
            !                call get_LSM_excitation_vector_PMC(this,N_order,N_calc,b)
            !            elseif(this%Problem_Type == 3) then !! IBC
            !                        call get_LSM_excitation_vector_IBC(this,N_order,N_calc,b)
            !            else !! Dielectric
            !        !                       call get_LSM_matrices_Dielectric(this,N_order,N_calc,Y,b,N_out)
            !                call get_LSM_excitation_vector_Dielectric(this,N_order,N_calc,b)
            !            endif
            !            write(*,*) b(1:10)
            !            stop

!            Lm_start_index = 1
!            Dm_start_index = 1
!            do p=1,this%group_index-1
!                Nm1 = this%Nm_collect(p)
!                N_out1 = (Nm1*(Nm1+1))/2
!                Lm_start_index = Lm_start_index + N_out1
!                Dm_start_index = Dm_start_index + Nm1
!            enddo
!            !            N_out = (Nm*(Nm+1))/2
!            Nm = this%Nm_collect(this%group_index)
!            N_out = (Nm*(Nm+1))/2
!
!            !            write(*,*) Nm,N_out
!            allocate(Lm(N_out),Dm(Nm),Y(N_out,0:N_order))
!
!
!            !            write(*,*) 'Dm_start_index',Dm_start_index,'  Lm_start_index',Lm_start_index
!            Lm = this%Lm_collect(Lm_start_index : (Lm_start_index+N_out-1))
!            Dm = this%Dm_collect(Dm_start_index : (Dm_start_index+Nm-1))
!            Y = this%Y_collect(Lm_start_index : (Lm_start_index+N_out-1),:)
            call this%Matrix_Storage%get_entry(this%group_index,Nm,Dm,Lm,Y,Y_mat,Yinv,&
                N_Taylor,IWORK,Y_mat_norm)

        endif
        !        write(*,*) b(1:10,0)
        !        write(*,*) Y(1:10,0)
        !! modification in the excitation vector, if required
        call modify_excitation(this,b,N_order,N_calc,rm_list,Nm)
        !! modify excitation vectors and Y matrix(vector) to remove some of the unknows from the requirements



        if(.not. use_old) then
            !! here, the Y matrices are modified removing certain rows and columns
            call modify_Y_matrix(N_order,Y,3*N_calc,Nm,rm_list)

            N_out = (Nm*(Nm+1))/2
        !        else
        !            if(allocated(rm_list)) then
        !                deallocate(rm_list)
        !            endif
        endif

        allocate(I_M_n(Nm,0:N_order))

        !        write(*,*) 'Nm',Nm,'N_out',N_out
        !        write(*,*) b(1:10,0)
        allocate(Y_0_vec(size(Y,1)),b_0_vec(size(b,1)))
        Y_0_vec = Y(:,0)
        b_0_vec = b(:,0)
        if(Matrix_Solution_Method == 1) then
            write(*,*) 'MATRIX SOLUTION: LU-Cholesky'
            I_M_n(:,0) = solve_LU_Cholesky(Y_0_vec,b_0_vec,Nm,Lm,Dm,use_old)
        elseif(Matrix_Solution_Method == 2) then
            write(*,*) 'MATRIX SOLUTION: SVD with Matrix Truncation'
            I_M_n(:,0) = solve_Matrix_SVD(Y_0_vec,b_0_vec,Nm,Y_mat,Yinv,use_old)

        else
            write(*,*) 'MATRIX SOLUTION: ERROR in choosing the matrix solution method'
            stop
        endif


        !        write(*,*) sum(abs(I_M_n(:,0)))

        if(N_order > 0) then
            allocate(y_n(Nm))


            do nn = 1,N_order
                if(.not. allocated(y_n)) then
                    !                    write(*,*) 'allocating y_m'
                    allocate(y_n(Nm))
                endif

                if(.not. allocated(Y_0_vec)) then
                    !                    write(*,*) 'allocating Y_0_vec'
                    allocate(Y_0_vec(size(Y,1)))
                endif

                Y_0_vec = Y(:,0)

                y_n = b(:,nn)
                do ii = 1,nn
                    y_n = y_n - dble(n_X_i_matrix(nn,ii))*(matrix_mult(Y(:,ii),I_M_n(:,nn-ii),Nm))
                enddo

                !                write(*,*) abs(sum(y_n(1:30)))
                if(Matrix_Solution_Method == 1) then
                    I_M_n(:,nn) = solve_LU_Cholesky(Y_0_vec,y_n,Nm,Lm,Dm,.true.)/(dble(n_X_i_matrix(nn,1)))
                elseif(Matrix_Solution_Method == 2) then
                    I_M_n(:,nn) = solve_Matrix_SVD(Y_0_vec,y_n,Nm,Y_mat,Yinv,.true.)/(dble(n_X_i_matrix(nn,1)))

                endif


            !                write(*,*) abs(I_M_n(15,nn))
            enddo
            if(allocated(y_n)) then
                !                write(*,*) 'deallocating y_m'
                deallocate(y_n)
            endif

            if(allocated(Y_0_vec)) then
                !                write(*,*) 'deallocating Y_0_vec'
                deallocate(Y_0_vec)
            endif

        endif

        !        write(*,*) sum(abs(Dm))

        if(.not. use_old) then
            this%is_sub_group_used(this%group_index) = .true.
!            if(allocated(this%Lm_collect)) then
!                allocate(Dm_temp(size(this%Dm_collect,1)),Lm_temp(size(this%Lm_collect,1)),&
!                Nm_collect_temp(size(this%Nm_collect,1)))
!                allocate(Y_collect_temp(size(this%Lm_collect,1),0:N_order))
!                Lm_temp = this%Lm_collect
!                Dm_temp = this%Dm_collect
!                Nm_collect_temp = this%Nm_collect
!                Y_collect_temp = this%Y_collect
!                deallocate(this%Y_collect)
!                deallocate(this%Lm_collect,this%Dm_collect,this%Nm_collect)
!                allocate(this%Dm_collect(size(Dm_temp,1)+Nm),this%Lm_collect(size(Lm_temp,1)+N_out),&
!                this%Nm_collect(size(Nm_collect_temp,1)+1) )
!                allocate(this%Y_collect(size(Lm_temp,1)+N_out,0:N_order))
!                this%Dm_collect = (/Dm_temp,Dm/)
!                this%Lm_collect = (/Lm_temp,Lm/)
!                this%Nm_collect = (/Nm_collect_temp,Nm/)
!                this%Y_collect(1:size(Lm_temp,1),:) = Y_collect_temp
!                this%Y_collect((size(Lm_temp,1)+1):(size(Lm_temp,1)+N_out),:) = Y
!                deallocate(Dm_temp,Lm_temp,Nm_collect_temp,Y_collect_temp)
!
!            else
!                allocate(this%Dm_collect(Nm),this%Lm_collect(N_out),this%Nm_collect(1) )
!                allocate(this%Y_collect(N_out,0:N_order))
!                this%Dm_collect = Dm
!                this%Lm_collect = Lm
!                this%Nm_collect = Nm
!                this%Y_collect = Y
!
!            !
!            endif
            call this%Matrix_Storage%put_entry(this%group_index,Nm,Dm,Lm,Y,Y_mat,&
                Yinv,N_order,IWORK,Y_mat_norm)

        endif
        call regain_current_length(N_order,I_M_n,rm_list)

        !        call cpu_time(finish_t)
        !        write(*,*) 'Matrix Evaluation and Inversion time =',finish_t-start_t,'seconds'
        p = 1
        do q=(this%N_curr-N_calc+1),this%N_curr
            this%I_sources(q,:) = vec(I_M_n(p,:),I_M_n(p+N_calc,:),I_M_n(p+2*N_calc,:),N_order+1)
            !            this%I_sources_total(q,:) = this%I_sources_total(q,:) + this%I_sources(q,:)

            !            write(*,*) abs(I_M_n(p,0)),abs(I_M_n(p+N_calc,0)),abs(I_M_n(p+2*N_calc,0))
            !                        write(*,*) q,abs(this%I_sources(q,:))
            p = p+1
        enddo
        !        write(*,*) '==============='
        !        stop
        !        write(*,*) 'average_currents',sum(abs(this%I_sources((this%N_curr-N_calc+1):this%N_curr,1)))
        deallocate(I_M_n)
        if(allocated(Lm)) then
            deallocate(Lm,Dm)
        endif
        deallocate(Y,b)
        if(allocated(Y_mat)) then
            deallocate(Y_mat,Yinv)
        endif

        if(Matrix_Solution_Method == 3.or. Matrix_Solution_Method == 4 .or. Matrix_Solution_Method == 5) then
            deallocate(IWORK,Y_mat_norm)
        endif
    end subroutine eval_matrices_RAS

    subroutine regain_current_length(N_order,I_m_n,rm_list)
        complex*16,allocatable,intent(inout) :: I_m_n(:,:)
        integer,allocatable,intent(inout) :: rm_list(:)
        complex*16,allocatable,dimension(:,:) :: I_temp
        integer :: ii,a_cnt,N_old,N_add,N_new,rm_cnt,N_order

        if(.not.  allocated(rm_list)) then
            return
        endif

        N_old = size(I_m_n,1)
        N_add = size(rm_list,1)
        N_new = N_old+N_add

        allocate(I_temp(N_new,0:N_order))

        !        do ii = 1,N_old
        !            write(*,*) I_m_n(ii,0)
        !        enddo

        a_cnt = 1
        rm_cnt = 1
        do ii=1,N_new
            if(rm_cnt <= N_add) then
                if(ii == rm_list(rm_cnt)) then
                    I_temp(ii,:) = 0.d0
                    rm_cnt = rm_cnt + 1
                    cycle
                endif
            endif
            I_temp(ii,:) = I_m_n(a_cnt,:)
            a_cnt = a_cnt + 1
        enddo
        deallocate(I_m_n)
        allocate(I_m_n(N_new,0:N_order))
        I_m_n = I_temp

        !        write(*,*) '=================='
        !        do ii = 1,N_new
        !            write(*,*) I_m_n(ii,0)
        !        enddo
        !        write(*,*) '=================='
        deallocate(rm_list,I_temp)
    end subroutine regain_current_length

    subroutine modify_Y_matrix(N_order,Y,Nm_old,Nm_new,remove_list)
        complex*16,allocatable,intent(inout) :: Y(:,:)
        complex*16,allocatable,dimension(:,:) :: Y_temp
        integer,allocatable,intent(inout) :: remove_list(:)
        integer :: ii,jj,N_remove_list,Nm_old,Nm_new,N_out,N_order
        integer :: rm_cnt,int_temp,N_out_old,y_cnt
        integer,allocatable,dimension(:) :: rm_list


        if(.not. allocated(remove_list)) then
            return
        endif

        N_remove_list = size(remove_list,1)
        N_out =  (Nm_new*(Nm_new+1))/2
        N_out_old =  (Nm_old*(Nm_old+1))/2
        allocate(Y_temp(N_out,0:N_order))
        allocate(rm_list(N_remove_list*Nm_old))

        !        write(*,*) '2*N_remove_list*Nm_old',2*N_remove_list*Nm_old

        rm_cnt = 1
        rm_list = 0
        do ii=1,N_remove_list
            do jj = 1,remove_list(ii)
                rm_list(rm_cnt) = mat2vec(Nm_old,remove_list(ii),jj)
                rm_cnt = rm_cnt + 1
            enddo
            !            write(*,*) remove_list(ii)
            do jj = (remove_list(ii)+1),Nm_old
                rm_list(rm_cnt) = mat2vec(Nm_old,jj,remove_list(ii))
                rm_cnt = rm_cnt + 1
            enddo
        enddo


        !! sorting rm_list array
        do ii =1,N_remove_list*Nm_old

            do jj = ii+1,N_remove_list*Nm_old
                if(rm_list(ii) == rm_list(jj)) then
                    rm_list(jj) = 0
                endif
                if(rm_list(ii) > rm_list(jj)) then !! swap condition
                    int_temp = rm_list(ii)
                    rm_list(ii) = rm_list(jj)
                    rm_list(jj) = int_temp
                endif
            enddo
        enddo
        !        write(*,*) 'N_out',N_out
        !        do ii=1,N_remove_list*Nm_old
        !            write(*,*) rm_list(ii)
        !        enddo
        !        stop
        rm_cnt = 1
        y_cnt = 1
        do ii = 1,N_out_old
            do while(rm_list(rm_cnt) == 0)
                rm_cnt = rm_cnt+1
                if(rm_cnt > size(rm_list,1)) then
                    exit
                endif
            enddo
            if(rm_list(rm_cnt) == ii) then
                rm_cnt = rm_cnt + 1
                cycle
            endif
            Y_temp(y_cnt,:) = Y(ii,:)
            y_cnt = y_cnt+1
        enddo
        deallocate(Y)
        allocate(Y(N_out,0:N_order))
        Y =Y_temp
        deallocate(Y_temp)
        deallocate(rm_list)
    end subroutine modify_Y_matrix

    subroutine modify_excitation(this,b,N_order,N_calc,remove_list,N_new)
        type(Scatterer) :: this
        complex*16,allocatable,intent(inout) :: b(:,:)
        complex*16,allocatable,dimension(:,:) :: b_temp
        integer,allocatable,intent(inout) :: remove_list(:)
        integer :: ii,N_remove,N_in,N_order,N_calc,b_cnt,jj,int_temp
        integer,intent(inout) :: N_new
        logical :: found_flag

        N_in = size(b,1)
        N_remove = 0
        do ii = (this%N_curr-N_calc+1),this%N_curr
            if(this%allowed_orientation(ii) == 1) then
                N_remove = N_remove + 2
            endif
        enddo
        if(N_remove == 0) then
            return
        endif
        allocate(remove_list(N_remove))
        b_cnt = 0
        do ii = (this%N_curr-N_calc+1),this%N_curr
            if(this%allowed_orientation(ii) == 1) then
                b_cnt = b_cnt + 1
                remove_list(b_cnt) = ii - this%N_curr + N_calc !! mark the x component
                b_cnt = b_cnt + 1
                remove_list(b_cnt) = remove_list(b_cnt-1) + N_calc !! mark the y component
            endif
        enddo
        !! sorting the remove_list vector
        do ii = 1,N_remove
            do jj = ii+1,N_remove
                if(remove_list(ii) > remove_list(jj)) then
                    int_temp = remove_list(ii)
                    remove_list(ii) = remove_list(jj)
                    remove_list(jj) = int_temp
                endif
            enddo
        enddo
        !        do ii = 1,N_remove
        !            write(*,*) remove_list(ii)
        !        enddo


        N_new = N_in - N_remove
        allocate(b_temp(N_new,0:N_order))
        b_cnt = 0
        found_flag = .false.
        do ii = 1,N_in
            do jj = 1,N_remove
                if(remove_list(jj) == ii) then

                    found_flag = .true.
                    exit
                endif
            enddo
            if(found_flag) then
                found_flag = .false.
                cycle
            endif

            b_cnt = b_cnt+1
            b_temp(b_cnt,:) = b(ii,:)
        enddo
        deallocate(b)
        allocate(b(N_new,0:N_order))
        b = b_temp
        deallocate(b_temp)
    end subroutine modify_excitation

    ! -----------------------------------------------------------------------
    ! Subroutine: eval_kernel_Magnetic_2D
    ! Purpose   : Evaluates the 2D electromagnetic field (E and H) radiated
    !             by a unit-amplitude MAGNETIC line source at position
    !             source_pos(source_n, :) evaluated at testing point test_j.
    !             Handles x-, y-, and z-oriented magnetic sources.
    !             When N_order > 0, also computes Taylor frequency-derivative
    !             arrays of order 0..N using Hankel_derivatives.
    !             thisD is the scatterer owning the source (may differ from
    !             'this' in multi-scatterer problems).
    ! Output: E_Mx/My/Mz, H_Mx/My/Mz -- each vector(0:N) of type vector_c
    ! -----------------------------------------------------------------------
    subroutine eval_kernel_Magnetic_2D(this,thisD,N,E_Mx,E_My,E_Mz,H_Mx,H_My,H_Mz,eta_l,az_l,ar_l,ak_l,test_j,source_n)
        type(Scatterer) :: this,thisD
        !! a function that evaluates the kernel of radiation
        type(vector_c),allocatable,intent(inout) :: E_Mx(:),E_My(:),E_Mz(:),H_Mx(:),H_My(:),H_Mz(:)
        real(8) :: az_l
        complex*16 :: ar_l,eta_l
        integer :: test_j,source_n
        real(8) :: R,x,y,xp,yp,xd,yd
        !        real(8) :: kr,kz,kRR
        integer :: N
        !        complex*16 :: h0,h2

        complex*16,dimension(0:N) :: k0H0,k0H1,k0H2,H1
        complex*16 :: ak_l,arR

        complex*16 :: const1
        complex*16 :: const2!!,CE_Iy
        complex*16,dimension(0:N) :: const3,const4,CE_Mz,CH_Mz,const5,const6
        complex*16,dimension(0:N) :: c_zero
        c_zero = (0.d0,0.d0)

        xp = this%source_pos(source_n,1)
        yp = this%source_pos(source_n,2)
        x = thisD%testing_pt(test_j,1)
        y = thisD%testing_pt(test_j,2)
        xd = x-xp
        yd = y-yp
        R = sqrt(xd**2.d0+yd**2.d0)


        !        kz = k0*cos(theta_i)
        !        kr = sqrt(k_l**2.d0 - kz**2.d0)

        !        kz = az_l*k0
        !        kr = ar_l*k0

        arR = ar_l*R

        !        kRR = kr*R
        !
        !        h0 = besselh2_0(kRR)
        !        h1(0) = besselh2_1(kRR)
        !        h2 = 2.d0/kRR*h1(0)-h0
        !
        !        k0H0 = k0*h0
        !        k0H1 = k0*h1
        !        k0H2 = k0*h2

        call Hankel_derivatives(arR,k0,N,k0H0,k0H1,k0H2,H1)


        !        CE_Mz = kr*eta_l/(4.d0*cj)*h1/R

        CE_Mz = ar_l*eta_l/((0.d0,4.d0))*k0H1/R
        E_Mz = vec(CE_Mz*yd,-CE_Mz*xd,c_zero,N+1)



        !        CH_Mz = cj*kz*kr/(4*k_l*R)*h1
        CH_Mz = -az_l*ar_l/((0.d0,4.d0)*ak_l*R)*k0H1
        H_Mz = vec(CH_Mz*xd,CH_Mz*yd,&
        -ar_l*ar_l/((4.d0,0.d0)*ak_l)*k0H0,N+1)
        !        -kr**2.d0/(4*k_l)*h0)


        const3 = az_l*eta_l/(4.d0,0.d0)*k0H0
        const4 = ar_l*eta_l/((0.d0,4.d0)*R)*k0H1

        E_Mx = vec(c_zero,&
        -const3,&
        -const4*yd,N+1)

        E_My = vec(const3,&
        c_zero,&
        const4*xd,N+1)


        !        E_My = vec(kz*eta_l/4.d0*h0,&
        !        c_zero,&
        !        kr*eta_l/(4.d0*cj*R)*h1*xd)


        !        CE_Iy = kr**2.d0*eta_l/(4.d0*k_l)

        !        H_Mx = vec(-k_l/4.d0*h0-CE_Iy/eta_l*(h2*(xd/R)**2.d0 - h1/kRR),&
        !        -CE_Iy/eta_l*xd*yd*h2/R**2.d0 ,&
        !        cj*kr*kz/(4.d0*k_l)*h1/R*xd)

!        CE_Iy = ar_l*ar_l/((4.d0,0.d0)*ak_l)
        const1 = ar_l*ar_l/((4.d0,0.d0)*ak_l)
        const2 = ar_l*az_l/((0.d0,-4.d0)*ak_l)
        const5 = -const1*xd*yd*k0H2/R**2.d0
        const4 = -ak_l/(4.d0,0.d0)*k0H0
        const3 = const2*k0H1/R
        const6 = H1/arR
        H_Mx = vec(const4-const1*(k0H2*(xd/R)**2.d0 - const6),&
        const5 ,&
        const3*xd,N+1)

        H_My = vec(const5 ,&
        const4-const1*(k0H2*(yd/R)**2.d0 - const6),&
        const3*yd ,N+1)

    !        CE_Iy = kr**2.d0*eta_l/(4.d0*k_l)
    !        H_My = vec(-CE_Iy/eta_l*xd*yd*h2/R**2.d0 ,&
    !        -k_l/4.d0*h0-CE_Iy/eta_l*(h2*(yd/R)**2.d0 - h1/kRR),&
    !        cj*kr*kz/(4.d0*k_l)*h1/R*yd )



    end subroutine eval_kernel_Magnetic_2D


    ! -----------------------------------------------------------------------
    ! Subroutine: eval_kernel_Electric_2D
    ! Purpose   : Evaluates the 2D electromagnetic field radiated by a unit-
    !             amplitude ELECTRIC line source at source_pos(source_n, :)
    !             evaluated at testing point test_j.
    !             Mirror of eval_kernel_Magnetic_2D but for electric sources.
    !             When N_order > 0, returns Taylor derivative arrays 0..N.
    ! Output: E_Ix/Iy/Iz, H_Ix/Iy/Iz -- each vector(0:N) of type vector_c
    ! -----------------------------------------------------------------------
    subroutine eval_kernel_Electric_2D(this,thisD,N,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_l,az_l,ar_l,ak_l,test_j,source_n)
        type(Scatterer) :: this,thisD
        !! a function that evaluates the kernel of radiation
        type(vector_c),allocatable,intent(inout)  :: E_Iz(:),E_Ix(:),E_Iy(:),H_Iz(:),H_Ix(:),H_Iy(:)
        complex*16 :: ar_l,eta_l !! the wave-number and characteristic impedance inside the medium
        real(8) :: az_l
        integer :: N
        integer :: test_j,source_n
        real(8) :: R,x,y,xp,yp,xd,yd
        !        real(8) :: kRR,kr,kz,
        !        complex*16 :: h0,h2

        complex*16,dimension(0:N) :: k0H0,k0H1,k0H2,h1
        complex*16 :: ak_l,arR
        complex*16,dimension(0:N) :: const1,const2,const3,const4,const5,const6,c_zero,CE_Iy,CE_Iz,CH_Iz

        c_zero = (0.d0,0.d0)
        !        complex*16 :: const2,const3,const4

        xp = this%source_pos(source_n,1)
        yp = this%source_pos(source_n,2)
        x = thisD%testing_pt(test_j,1)
        y = thisD%testing_pt(test_j,2)
        xd = x-xp
        yd = y - yp
        R = sqrt(xd**2.d0+yd**2.d0)

        !        kz = k0*cos(theta_i)
        !        kr = sqrt(k_l**2.d0 - kz**2.d0)

        !        kz = az_l*k0
        !        kr = ar_l*k0


        arR = ar_l*R

        !        kRR = kr*R
        !
        !        h0 = besselh2_0(kRR)
        !        h1 = besselh2_1(kRR)
        !        h2 = 2.d0/kRR*h1-h0
        !
        !        k0H0 = k0*h0
        !        k0H1 = k0*h1
        !        k0H2 = k0*h2
        call Hankel_derivatives(arR,k0,N,k0H0,k0H1,k0H2,H1)

        CE_Iy = ar_l*ar_l*eta_l/((4.d0,0.d0)*ak_l)
        const3 = -ak_l*eta_l/(4.d0,0.d0)*k0H0
        const4 = -CE_Iy*xd*yd*k0H2/R**2.d0
        const5 = -ar_l*az_l*eta_l/((0.d0,4.d0)*ak_l*R)*k0H1
        const6 = H1/arR

        E_Iy = vec(const4 ,&
        const3-CE_Iy*(k0H2*(yd/R)**2.d0 - const6),&
        const5*yd ,N+1)

        E_Ix = vec(const3-CE_Iy*(k0h2*(xd/R)**2.d0 - const6),&
        const4 ,&
        const5*xd,N+1)



        CE_Iz = -az_l*ar_l*k0H1*eta_l/((0.d0,4.d0)*ak_l*R)
        E_Iz = vec(CE_Iz*xd,CE_Iz*yd,&
        -ar_l*ar_l*eta_l/((4.d0,0.d0)*ak_l)*k0H0,N+1)

        CH_Iz = ar_l/(0.d0,4.d0)*k0H1/R
        H_Iz = vec(-CH_Iz*yd,CH_Iz*xd,c_zero,N+1)

        const1 = az_l/(4.d0,0.d0)*k0H0
        const2 = ar_l/((0.d0,4.d0)*R)*k0H1

        H_Ix = vec(c_zero,&
        const1,&
        const2*yd,N+1)

        H_Iy = vec(-const1, &
        c_zero,&
        -const2*xd,N+1)

    !CE_Iy = kr**2.d0*eta_l/(4.d0*k_l)
    !
    !        E_Iy = vec(-CE_Iy*xd*yd*h2/R**2.d0 ,&
    !        -k_l*eta_l/4.d0*h0-CE_Iy*(h2*(yd/R)**2.d0 - h1/kRR),&
    !        cj*kr*kz*eta_l/(4.d0*k_l)*h1/R*yd )
    !
    !        E_Ix = vec(-k_l*eta_l/4.d0*h0-CE_Iy*(h2*(xd/R)**2.d0 - h1/kRR),&
    !        -CE_Iy*xd*yd*h2/R**2.d0 ,&
    !        cj*kr*kz*eta_l/(4.d0*k_l)*h1/R*xd)
    !
    !        H_Iy = vec(-kz/4.d0*h0, &
    !        c_zero,&
    !        -kr/(4.d0*cj*R)*h1*xd)
    !
    !        CE_Iz = cj*kz*kr*h1*eta_l/(4.d0*k_l*R)
    !        E_Iz = vec(CE_Iz*xd,CE_Iz*yd,&
    !        -kr**2.d0*eta_l/(4.d0*k_l)*h0)
    !
    !        CH_Iz = kr/(4.d0*cj)*h1/R
    !        H_Iz = vec(-CH_Iz*yd,CH_Iz*xd,c_zero)
    !
    !
    !        H_Ix = vec(c_zero,&
    !        kz/4.d0*h0,&
    !        kr/(4.d0*cj*R)*h1*yd)

    end subroutine eval_kernel_Electric_2D


    subroutine deallocate_initialized_arrays(this)
        type(Scatterer) :: this
        deallocate(this%source_pos,this%E,this%H,this%Es_current,this%EH_excit,this%Hs_current,&
        this%E_bounce_current,this%H_bounce_current)
        if(allocated(this%allowed_orientation)) then
            deallocate(this%allowed_orientation)
        endif

        !        if(allocated(this%segments_points)) then
        !            deallocate(this%segments_points)
        !        endif


        if(allocated(this%I_u_bounce))then
            deallocate(this%I_u_bounce,this%I_v_bounce,this%M_u_bounce,this%M_v_bounce)

        endif
        if(allocated(this%combined_bouncing_fields))then
            deallocate(this%combined_bouncing_fields)
        endif
        if(allocated(this%normalize_H))then
            deallocate(this%normalize_E,this%normalize_H)
        endif

        if(allocated(this%contour_point_active_region)) then
            deallocate(this%contour_point_active_region)
        endif


        if(allocated(this%testing_pt)) then
            deallocate(this%testing_pt)
        endif
        if(allocated(this%norm_v)) then
            deallocate(this%norm_v)
        endif

        if(allocated(this%tang_u)) then
            deallocate(this%tang_u)
        endif

        if(allocated(this%delta_n)) then
            deallocate(this%delta_n)
        endif

        if(allocated(this%tang_v)) then
            deallocate(this%tang_v)
        endif


        deallocate(this%E_current,this%H_current)
        deallocate(this%Es,this%Hs)
        deallocate(this%I_sources, this%I_stat)
        deallocate(this%I_sources_total)
        deallocate(this%I_sources_total_current)
        if(allocated(this%Lm_collect)) then
            deallocate(this%Nm_collect,this%Lm_collect,this%Dm_collect)
        endif
        if(allocated(this%Y_collect)) then
            deallocate(this%Y_collect)
        endif
        if(allocated(this%is_sub_group_used)) then
            deallocate(this%is_sub_group_used)
        endif

        if(allocated(this%norm_eta_v)) then
            deallocate(this%norm_eta_v,this%norm_eta_u)
        endif

        if(allocated(this%n_sources_array)) then
            deallocate(this%n_sources_array)
        endif
        if(allocated(this%n_sources_outside_array)) then
            deallocate(this%n_sources_outside_array)
        endif


        call this%Matrix_Storage%destroy()

        deallocate(this%Matrix_Storage)

        !        if(allocated(this%Str_Rules)) then
        !            deallocate(this%Str_Rules)
        !        endif

        if(allocated(this%eta_uu)) then
            deallocate(this%eta_vv,this%eta_uu,this%eta_vu,this%eta_uv)
        endif
        if(allocated(this%Ed_current)) then !! Dielectric Case
            deallocate(this%Ed_current,this%Hd_current)
            deallocate(this%Ed,this%Hd)

        endif
        if(allocated(this%Inside_sources_dielectric)) then
            deallocate(this%Inside_sources_dielectric)
        endif
        deallocate(this%active_region)
        if(allocated(this%Hd_c)) then
            deallocate(this%Hd_c,this%Ed_c)
        endif
        if(allocated(this%Hs_c)) then
            deallocate(this%Es_c,this%Hs_c)
        endif

        if(allocated(this%MoM_matrix_inverted)) then
            deallocate(this%MoM_matrix_inverted)
        endif
        if(allocated(this%testing_pt_MoM)) then
            deallocate(this%testing_pt_MoM)
        endif
        if(allocated(this%M_v)) then
            deallocate(this%I_u,this%I_v,this%M_u,this%M_v)
        endif
        if(allocated(this%Eu_MoM)) then
            deallocate(this%Eu_MoM,this%Ev_MoM,this%Hu_MoM,this%Hv_MoM)
        endif
        deallocate(this%E_bounce,this%H_bounce)

        if(allocated(this%Pade_a)) then
            deallocate(this%Pade_a,this%Pade_b)
        endif

        if(allocated(this%MoM_Pade_a_I_u)) then
            deallocate(this%MoM_Pade_a_I_u,this%MoM_Pade_b_I_u)
            deallocate(this%MoM_Pade_a_M_u,this%MoM_Pade_b_M_u)
            deallocate(this%MoM_Pade_a_I_v,this%MoM_Pade_b_I_v)
            deallocate(this%MoM_Pade_a_M_v,this%MoM_Pade_b_M_v)
        endif

        if(allocated(this%Chebyshev_c)) then
            deallocate(this%Chebyshev_c,this%Chebyshev_d)
        endif

        if(allocated(this%MoM_Chebyshev_c_I_u)) then
            deallocate(this%MoM_Chebyshev_c_I_u,this%MoM_Chebyshev_d_I_u)
            deallocate(this%MoM_Chebyshev_c_M_u,this%MoM_Chebyshev_d_M_u)
            deallocate(this%MoM_Chebyshev_c_I_v,this%MoM_Chebyshev_d_I_v)
            deallocate(this%MoM_Chebyshev_c_M_v,this%MoM_Chebyshev_d_M_v)
        endif

        if(allocated(this%contour_points_type)) then
            deallocate(this%contour_points_type,this%contour_points)
        endif
        if(allocated(this%MoM_ZY)) then
            deallocate(this%MoM_ZY)
        endif

        if(allocated(this%contour_points_orientation)) then
            deallocate(this%contour_points_orientation)
        endif

        if(allocated(this%testing_pt_status)) then
            deallocate(this%testing_pt_status)
        endif
    end subroutine deallocate_initialized_arrays

    subroutine allocate_arrays(this,eval_inc_flag,Scatterers_pointer)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)
        logical :: eval_inc_flag



        allocate(this%Matrix_Storage)


        allocate(this%E(2*this%M,0:N_taylor),this%H(2*this%M,0:N_taylor),&
        this%E_bounce_current(2*this%M,0:N_taylor),this%H_bounce_current(2*this%M,0:N_taylor))
        allocate(this%E_bounce(2*this%M,0:N_taylor),this%H_bounce(2*this%M,0:N_taylor))
        allocate(this%E_current(2*this%M,0:N_taylor),this%H_current(2*this%M,0:N_taylor))
        allocate(this%Es_current(2*this%M,0:N_taylor),this%Hs_current(2*this%M,0:N_taylor))
        allocate(this%Es(2*this%M,0:N_taylor),this%Hs(2*this%M,0:N_taylor))
        allocate(this%source_pos(this%N_max,3))
        allocate(this%allowed_orientation(this%N_max))
        allocate(this%I_sources(this%N_max,0:N_taylor),this%I_sources_total_current(this%N_max,0:N_taylor),this%I_stat(this%N_max))
        allocate(this%I_sources_total(this%N_max,0:N_taylor))
        allocate(this%active_region(this%N_max))
        allocate(this%Hd_c(2*this%M,0:N_taylor),this%Ed_c(2*this%M,0:N_taylor),&
        this%Es_c(2*this%M,0:N_taylor),this%Hs_c(2*this%M,0:N_taylor))
        !        if(this%Problem_Type == 4) then !! Dielectric Case
        allocate(this%Ed_current(2*this%M,0:N_taylor),this%Hd_current(2*this%M,0:N_taylor))
        allocate(this%Ed(2*this%M,0:N_taylor),this%Hd(2*this%M,0:N_taylor))
            !! both fields should be used as boundary conditions
        !        else
        !            allocate(this%EH_excit(2*this%M,0:N_taylor))
        !        endif
        allocate(this%EH_excit(4*this%M,0:N_taylor))
        if(eval_inc_flag) then
            allocate(this%Eu(this%M,0:N_taylor),this%Ev(this%M,0:N_taylor),this%Hu(this%M,0:N_taylor),this%Hv(this%M,0:N_taylor))
            call  set_excitation_vectors(this,this%Eu,this%Ev,this%Hu,this%Hv)
            this%E_current(1:this%M,:) = this%Eu
            this%E_current((this%M+1):(2*this%M),:) = this%Ev
            this%H_current(1:this%M,:) = this%Hu
            this%H_current((this%M+1):(2*this%M),:) = this%Hv

            !        this%E_current = (/this%Eu,this%Ev/)
            !        this%H_current = (/this%Hu,this%Hv/)
            this%E = this%E_current
            this%H = this%H_current
!            write(*,*) this%H
            deallocate(this%Eu,this%Ev,this%Hu,this%Hv)
        endif



        call initialize_parameters(this,N_taylor,Scatterers_pointer)
    end subroutine allocate_arrays

    subroutine set_excitation_vectors(this,Eu,Ev,Hu,Hv)
        type(Scatterer) :: this

        integer :: j
        complex*16,allocatable,intent(inout):: Eu(:,:),Ev(:,:),Hu(:,:),Hv(:,:)
!        type(vector) :: z_vec

        if(allocated(this%normalize_E)) then
            deallocate(this%normalize_E,this%normalize_H)
        endif
        allocate(this%normalize_E(0:N_taylor),this%normalize_H(0:N_taylor))

        Eu = 0.d0
        Ev = 0.d0
        Hu = 0.d0
        Hv = 0.d0

        this%normalize_E = 0.d0
        this%normalize_H = 0.d0
!        z_vec = vec(0.d0,0.d0,1.d0)
        do j =1,this%M
            if(this%testing_pt_status(j,3) /= 0) then !! the incident electric field doesn't exsist at this region
                cycle
            endif
            Eu(j,:) = dot(this%tang_u(j),this%Ei(j,:),N_taylor+1)
            Ev(j,:) = dot(this%tang_v(j),this%Ei(j,:),N_taylor+1)
            Hu(j,:) = dot(this%tang_u(j),this%Hi(j,:),N_taylor+1)
            Hv(j,:) = dot(this%tang_v(j),this%Hi(j,:),N_taylor+1)
            !            write(*,*) Ht(j)
            this%normalize_E = this%normalize_E + abs(Eu(j,:)/eta1)**2.d0 + abs(Ev(j,:)/eta1)**2.d0
            this%normalize_H = this%normalize_H + abs(Hu(j,:))**2.d0 + abs(Hv(j,:))**2.d0
        enddo
        deallocate(this%Ei,this%Hi)
    end subroutine set_excitation_vectors

    ! -----------------------------------------------------------------------
    ! Subroutine: set_excitation
    ! Purpose   : Projects the incident electric (Ei) and magnetic (Hi)
    !             field vectors onto the tangential directions tang_u, tang_v
    !             at each testing point and stores the result in
    !             EH_excit(1:4*M, 0:N_taylor). Also computes and stores the
    !             normalisation factors normalize_E and normalize_H used in
    !             the error criterion. For N_taylor > 0, Taylor-series
    !             derivatives of the incident field are also accumulated.
    ! -----------------------------------------------------------------------
    subroutine set_excitation(this,points,Ei,Hi)
        type(Scatterer) :: this
        real(8),allocatable,intent(inout) :: points(:,:)
        integer:: j,nn,ii
        type(vector) :: Ep,En,rho_hat,test_pt
        complex*16 :: exponential,exp_arg_const
        type(vector_c),allocatable,intent(inout) :: Ei(:,:),Hi(:,:)
        type(vector_c),allocatable,dimension(:)  :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,E_cont,H_cont
        real(8),allocatable,dimension(:,:) :: scr_pos
        real(8) :: Excitation_norm,Ei_max,Ei_curr

        if(Wideband_type == 0) then !! a wide band solution not is required

            N_taylor = 0 !! discard whatever input from the file

        endif



        allocate(Ei(this%M,0:N_taylor),Hi(this%M,0:N_taylor))

        Ep = vec(-cos(theta_i)*cos(phi_i),-cos(theta_i)*sin(phi_i),sin(theta_i))
        En = vec(sin(phi_i),-cos(phi_i),0.d0)
        rho_hat = vec(sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i))

        do j=1,this%M
            do nn = 1,N_taylor
                Ei(j,nn) = vec(c_zero,c_zero,c_zero)
                Hi(j,nn) = vec(c_zero,c_zero,c_zero)
            enddo
        enddo
        if(Source_Model == 0) then
            do j=1,this%M
                if(this%testing_pt_status(j,3) /= 0) then !! the incident electric field doesn't exsist at this region
                    cycle
                endif
                test_pt = vec(points(j,1),points(j,2),points(j,3))
                exp_arg_const = cj*ak*dot(rho_hat,test_pt)
                exponential = exp(exp_arg_const*k0)
                !                write(*,*) dot(rho_hat,testing_pt_MoM(j))
                Ei(j,0) = eta1*exponential*(cos(alpha_i)*Ep+sin(alpha_i)*En)
                !                write(*,*) Ei(j)%v
                Hi(j,0) = exponential*(sin(alpha_i)*Ep-cos(alpha_i)*En)

                do nn = 1,N_taylor
                    Ei(j,nn) = exp_arg_const*Ei(j,nn-1)
                    Hi(j,nn) = exp_arg_const*Hi(j,nn-1)
                enddo

            enddo
        elseif(Source_Model == 1) then
            allocate(E_Ix(0:N_taylor),E_Iy(0:N_taylor),E_Iz(0:N_taylor),&
            H_Ix(0:N_taylor),H_Iy(0:N_taylor),H_Iz(0:N_taylor),E_cont(0:N_taylor),H_cont(0:N_taylor))
            allocate(scr_pos(1,2))
            do j=1,this%M
                if(this%testing_pt_status(j,3) /= 0) then !! the incident electric field doesn't exsist at this region
!                    write(*,*) 'ba-cycle',j
                    cycle
                endif
                do nn = 0,N_taylor
                    E_cont(nn) = vec(c_zero,c_zero,c_zero)
                enddo
                H_cont = E_cont
                do ii = 1,N_Line_sources
                    scr_pos(1,:) = (/Line_sources(ii)%x_s,Line_sources(ii)%y_s /)
                    if(Line_sources(ii)%Source_Type == 2) then !! magnetic Line source
                        call eval_near_field_Magnetic_2D(N_taylor,E_Ix,E_Iy,E_Iz,H_Ix,&
                        H_Iy,H_Iz,ak,az,ar,eta1,points(j,1),points(j,2),1,scr_pos)
                    else !! Electric line source
                        call eval_near_field_Electric_2D(N_taylor,E_Ix,E_Iy,E_Iz,H_Ix,&
                        H_Iy,H_Iz,ak,az,ar,eta1,points(j,1),points(j,2),1,scr_pos)
                    endif
                    if(Line_sources(ii)%Source_Orientation == 1) then !! z-oriented
                        do nn = 0,N_taylor
                            E_cont(nn) = E_cont(nn)+Line_sources(ii)%Amp*E_Iz(nn)
                            H_cont(nn) = H_cont(nn)+Line_sources(ii)%Amp*H_Iz(nn)
                        enddo
                    elseif(Line_sources(ii)%Source_Orientation == 2) then !! x-oriented
                        do nn = 0,N_taylor
                            E_cont(nn) = E_cont(nn)+Line_sources(ii)%Amp*E_Ix(nn)
                            H_cont(nn) = H_cont(nn)+Line_sources(ii)%Amp*H_Ix(nn)
                        enddo
                    else !! y-oriented
                        do nn = 0,N_taylor
                            E_cont(nn) = E_cont(nn)+Line_sources(ii)%Amp*E_Iy(nn)
                            H_cont(nn) = H_cont(nn)+Line_sources(ii)%Amp*H_Iy(nn)
                        enddo
                    endif
                enddo

                Ei(j,:) = E_cont
                Hi(j,:) = H_cont

            enddo
            deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,scr_pos,E_cont,H_cont)
            if(.not. Line_Sources(1)%is_normalized) then
!            Ei_max = absolute(Ei(1,0))
!            do j = 2,this%M
!                Ei_curr = absolute(Ei(j,0))
!                if(Ei_curr > Ei_max) then
!                    Ei_max = Ei_curr
!                endif
!            enddo
!            if(Ei_max == 0.d0) then
!                return
!            endif
!            Excitation_norm = eta1/Ei_max
            Ei_max = 0.d0
            do j=1,this%M
                Ei_curr = absolute(Ei(j,0))/eta1
                Ei_max = Ei_max + Ei_curr
            enddo
            if(Ei_max < 1.d-10) then
                return
            endif
            Excitation_norm = dble(this%M/4)/Ei_max
!            write(*,*) Excitation_norm
            do j=1,this%M
                do nn = 0,N_taylor
                    Ei(j,nn) = Excitation_norm*Ei(j,nn)
                    Hi(j,nn) = Excitation_norm*Hi(j,nn)

                enddo
            enddo
            do ii = 1,N_Line_Sources
                Line_Sources(ii)%Amp = Line_Sources(ii)%Amp*Excitation_norm
                Line_Sources(ii)%is_normalized = .true.
            enddo
            endif

        else
            write(*,*) 'ERROR: unidentified source model; it should be either 0 or 1'
            stop
        endif
    end subroutine set_excitation

    subroutine set_excitation_v(this,points_v,Ei,Hi)
        type(Scatterer) :: this
        type(vector),allocatable,intent(inout) :: points_v(:)
        integer:: j,nn,ii
        type(vector) :: Ep,En,rho_hat,test_pt
        complex*16 :: exponential,exp_arg_const
        type(vector_c),allocatable,intent(inout) :: Ei(:,:),Hi(:,:)
        type(vector_c),allocatable,dimension(:)  :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,E_cont,H_cont
        real(8),allocatable,dimension(:,:) :: scr_pos
!        real(8) :: Excitation_norm,Ei_max,Ei_curr

        allocate(Ei(this%M,0:N_taylor),Hi(this%M,0:N_taylor))

        Ep = vec(-cos(theta_i)*cos(phi_i),-cos(theta_i)*sin(phi_i),sin(theta_i))
        En = vec(sin(phi_i),-cos(phi_i),0.d0)
        rho_hat = vec(sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i))

        do j=1,this%M
            do nn = 1,N_taylor
                Ei(j,nn) = vec(c_zero,c_zero,c_zero)
                Hi(j,nn) = vec(c_zero,c_zero,c_zero)
            enddo
        enddo

        if(Source_Model == 0) then
            do j=1,this%M
                if(this%testing_pt_status(j,3) /= 0) then !! the incident electric field doesn't exsist at this region
                    cycle
                endif
                test_pt = points_v(j)
                exp_arg_const = cj*ak*dot(rho_hat,test_pt)
                exponential = exp(exp_arg_const*k0)
                !                write(*,*) dot(rho_hat,testing_pt_MoM(j))
                Ei(j,0) = eta1*exponential*(cos(alpha_i)*Ep+sin(alpha_i)*En)
!                                write(*,*) Ei(j,0)%v
                Hi(j,0) = exponential*(sin(alpha_i)*Ep-cos(alpha_i)*En)


                do nn = 1,N_taylor
                    Ei(j,nn) = exp_arg_const*Ei(j,nn-1)
                    Hi(j,nn) = exp_arg_const*Hi(j,nn-1)
                enddo
            enddo
!            stop
        elseif(Source_Model == 1) then
!            write(*,*) 'here', N_taylor

            allocate(E_Ix(0:N_taylor),E_Iy(0:N_taylor),E_Iz(0:N_taylor),&
            H_Ix(0:N_taylor),H_Iy(0:N_taylor),H_Iz(0:N_taylor),E_cont(0:N_taylor),H_cont(0:N_taylor))
            allocate(scr_pos(1,2))
            do j=1,this%M
                if(this%testing_pt_status(j,3) /= 0) then !! the incident electric field doesn't exsist at this region
!                    write(*,*) 'ba-cycle',j
                    cycle
                endif
                do nn = 0,N_taylor
                    E_cont(nn) = vec(c_zero,c_zero,c_zero)
                enddo
                H_cont = E_cont
!                write(*,*) 'N_Line_sources',N_Line_sources
                do ii = 1,N_Line_sources
                    scr_pos(1,:) = (/Line_sources(ii)%x_s,Line_sources(ii)%y_s /)
!                    write(*,*) scr_pos(1,:),Line_sources(ii)%Amp
                    if(Line_sources(ii)%Source_Type == 2) then !! magnetic Line source
                        call eval_near_field_Magnetic_2D(N_taylor,E_Ix,E_Iy,E_Iz,H_Ix,&
                        H_Iy,H_Iz,ak,az,ar,eta1,points_v(j)%v(1),points_v(j)%v(2),1,scr_pos)
                    else !! Electric line source
                        call eval_near_field_Electric_2D(N_taylor,E_Ix,E_Iy,E_Iz,H_Ix,&
                        H_Iy,H_Iz,ak,az,ar,eta1,points_v(j)%v(1),points_v(j)%v(2),1,scr_pos)
                    endif
                    if(Line_sources(ii)%Source_Orientation == 1) then !! z-oriented
                        do nn = 0,N_taylor
                            E_cont(nn) = E_cont(nn)+Line_sources(ii)%Amp*E_Iz(nn)
                            H_cont(nn) = H_cont(nn)+Line_sources(ii)%Amp*H_Iz(nn)
                        enddo
                    elseif(Line_sources(ii)%Source_Orientation == 2) then !! x-oriented
                        do nn = 0,N_taylor
                            E_cont(nn) = E_cont(nn)+Line_sources(ii)%Amp*E_Ix(nn)
                            H_cont(nn) = H_cont(nn)+Line_sources(ii)%Amp*H_Ix(nn)
                        enddo
                    else !! y-oriented
                        do nn = 0,N_taylor
                            E_cont(nn) = E_cont(nn)+Line_sources(ii)%Amp*E_Iy(nn)
                            H_cont(nn) = H_cont(nn)+Line_sources(ii)%Amp*H_Iy(nn)
                        enddo
                    endif
                enddo

                Ei(j,:) = E_cont
                Hi(j,:) = H_cont

!                write(*,*) Line_sources(ii)%Amp
!                write(*,*) E_cont

            enddo
            deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,scr_pos,E_cont,H_cont)
!           Ei_max = absolute(Ei(1,0))
!            do j = 2,this%M
!                Ei_curr = absolute(Ei(j,0))
!                if(Ei_curr > Ei_max) then
!                    Ei_max = Ei_curr
!                endif
!            enddo
!
!            if(Ei_max == 0.d0) then
!                Ei_max = eta1
!            endif
!            Excitation_norm = eta1/Ei_max
!            write(*,*) Excitation_norm
!            do j=1,this%M
!                do nn = 0,N_taylor
!                    Ei(j,nn) = Excitation_norm*Ei(j,nn)
!                    Hi(j,nn) = Excitation_norm*Hi(j,nn)
!
!                enddo
!            enddo
        else
            write(*,*) 'ERROR: unidentified source model; it should be either 0 or 1'
            stop
        endif
    end subroutine set_excitation_v

    ! -----------------------------------------------------------------------
    ! Subroutine: set_Et
    ! Purpose   : Assembles the total excitation EH_excit = incident +
    !             scattered (from other sources already accepted into the
    !             solution) + bouncing (from other scatterers if
    !             include_bouncing=.true.).
    !             This is the residual that the next group of sources must
    !             satisfy; called after each iteration before eval_error.
    ! -----------------------------------------------------------------------
    subroutine set_Et(this,N_order,include_bouncing,Scatterers_pointer)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)
        integer :: N_order,jj
        logical :: include_bouncing
        real(8) :: bounce_flag
        complex*16 :: eta_ex
        !        complex*16,allocatable,dimension(:,:) :: H_u,H_v
        complex*16,dimension(0:N_order) :: H_u,H_v

        if(include_bouncing) then
            bounce_flag = 1.d0
        else
            bounce_flag = 0.d0
        endif


        do jj = 1,this%M
            if(this%testing_pt_status(jj,1) /= -1) then
                if(this%testing_pt_status(jj,3) == 0) then
                    eta_ex = eta1
                else
                    eta_ex = Scatterers_pointer(this%testing_pt_status(jj,3))%eta_local
                endif

                this%EH_excit(jj,:) = (this%E_current(jj,:)+ bounce_flag*this%E_bounce_current(jj,:) + &
                this%Es_current(jj,:) + this%Es_c(jj,:)- this%Ed_current(jj,:) - this%Ed_c(jj,:))/eta_ex


                this%EH_excit(jj+this%M,:) = (this%E_current(jj+this%M,:)+ bounce_flag*this%E_bounce_current(jj+this%M,:) + &
                this%Es_current(jj+this%M,:) + this%Es_c(jj+this%M,:)-this%Ed_current(jj+this%M,:) - this%Ed_c(jj+this%M,:))/eta_ex
            else !! the point is cancelled
                this%EH_excit(jj,:) = 0.d0
                this%EH_excit(jj+this%M,:) = 0.d0
            endif
        enddo


        this%EH_excit((2*this%M+1):(4*this%M),:) = (this%H_current+ bounce_flag*this%H_bounce_current + &
        this%Hs_current - this%Hd_current + this%Hs_c - this%Hd_c)

        do jj = 1,this%M
            if(this%testing_pt_status(jj,1) /= 3) then !! only IBC can do this loop
                cycle
            endif

            H_u = this%EH_excit(2*this%M+jj,:)
            H_v = this%EH_excit(3*this%M+jj,:)

            !            H_u = (this%H_current(1:this%M,:)+ bounce_flag*this%H_bounce_current(1:this%M,:)  + &
            !            this%Hs_current(1:this%M,:) + this%Hs_c(1:this%M,:))
            !            H_v = (this%H_current((this%M+1):(2*this%M),:)+ bounce_flag*this%H_bounce_current((this%M+1):(2*this%M),:) + &
            !            this%Hs_current((this%M+1):(2*this%M),:) + this%Hs_c((this%M+1):(2*this%M),:))
            !            do nn = 0,N_order
            !                this%EH_excit(1:this%M,nn) = this%norm_eta_u*(this%EH_excit(1:this%M,nn) +  &
            !                this%eta_uu*H_v(:,nn) - this%eta_uv*H_u(:,nn)) !! E_u
            !                this%EH_excit((this%M+1):2*this%M,nn) = this%norm_eta_v*(this%EH_excit((this%M+1):2*this%M,nn) + &
            !                this%eta_vu*H_v(:,nn) - this%eta_vv*H_u(:,nn)) !! E_v
            !            enddo
            this%EH_excit(jj,:) = this%norm_eta_u(jj)*(this%EH_excit(jj,:) +  &
            this%eta_uu(jj)*H_v - this%eta_uv(jj)*H_u) !! E_u
            this%EH_excit(this%M+jj,:) = this%norm_eta_v(jj)*(this%EH_excit(this%M+jj,:) + &
            this%eta_vu(jj)*H_v - this%eta_vv(jj)*H_u) !! E_v

        enddo
    !        write(*,*)  'E_current ',sum(abs(this%E_current))
    !        write(*,*)  'Es_current',sum(abs(this%Es_current))
    !        write(*,*)  'Es_c      ',sum(abs(this%Es_c))
    !        write(*,*)  'Ed_current',sum(abs(this%Ed_current))
    !        write(*,*)  'Ed_c      ',sum(abs(this%Ed_c))
    !        write(*,*)  'H_current ',sum(abs(this%H_current))
    !        write(*,*)  'Hs_current',sum(abs(this%Hs_current))
    !        write(*,*)  'Hs_c      ',sum(abs(this%Hs_c))
    !        write(*,*)  'Hd_current',sum(abs(this%Hd_current))
    !        write(*,*)  'Hd_c      ',sum(abs(this%Hd_c))
    !        write(*,*)  'EH        ',sum(abs(this%EH_excit(1:2*this%M,0)))
    !        write(*,*) sum(abs(this%EH_excit(1:2*this%M,0)))
    !        if(this%Problem_Type == 1) then !! PEC
    !            this%EH_excit = (this%E_current+ bounce_flag*this%E_bounce_current + this%Es_current + this%Es_c)/eta1
    !
    !        elseif(this%Problem_Type == 2) then !! PMC
    !            this%EH_excit = this%H_current+ bounce_flag*this%H_bounce_current + this%Hs_current + this%Hs_c
    !        elseif(this%Problem_Type == 3) then !! IBC
    !            this%EH_excit = (this%E_current+ bounce_flag*this%E_bounce_current + this%Es_current + this%Es_c)/eta1
    !            this%EH_excit(1:this%M,:) = this%norm_eta_u*this%EH_excit(1:this%M,:)
    !            this%EH_excit((this%M+1):(2*this%M),:) = this%norm_eta_v*this%EH_excit((this%M+1):(2*this%M),:)
    !            !! then add the magnetic field part in case of having IBC
    !            allocate(H_u(this%M,0:N_taylor),H_v(this%M,0:N_taylor))
    !            H_u = (this%H_current(1:this%M,:)+ bounce_flag*this%H_bounce_current(1:this%M,:)  + &
    !            this%Hs_current(1:this%M,:) + this%Hs_c(1:this%M,:))
    !            H_v = (this%H_current((this%M+1):(2*this%M),:)+ bounce_flag*this%H_bounce_current((this%M+1):(2*this%M),:) + &
    !            this%Hs_current((this%M+1):(2*this%M),:) &
    !            + this%Hs_c((this%M+1):(2*this%M),:))
    !            do ii = 1,this%M
    !                this%EH_excit(ii,:) = this%EH_excit(ii,:) +  this%norm_eta_u*(this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:)) !! E_u
    !                this%EH_excit(this%M+ii,:) = this%EH_excit(this%M+ii,:) + &
    !                this%norm_eta_v*(this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:)) !! E_v
    !        !
    !        !                this%EH_excit(ii,:) = this%EH_excit(ii,:) -  (this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:)) !! E_u
    !        !                this%EH_excit(this%M+ii,:) = this%EH_excit(this%M+ii,:) - (this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:)) !! E_v
    !            enddo
    !            deallocate(H_u,H_v)
    !        else !! Dielectric
    !            !            write(*,*) sum(abs(this%Es))/sqrt(this%normalize_E),sum(abs(this%Ed))/sqrt(this%normalize_E),&
    !            !            sum(abs(this%Hs))/sqrt(this%normalize_H),sum(abs(this%Hd))/sqrt(this%normalize_H)
    !            this%EH_excit(1:(2*this%M),:) = (this%E_current+ bounce_flag*this%E_bounce_current  + this%Es_current -&
    !             this%Ed_current + this%Es_c - this%Ed_c)/eta1
    !            this%EH_excit((2*this%M+1):(4*this%M),:) = (this%H_current+ bounce_flag*this%H_bounce_current + &
    !            this%Hs_current - this%Hd_current + this%Hs_c - this%Hd_c)
    !        endif


    end subroutine set_Et

    ! -----------------------------------------------------------------------
    ! Function: eval_error
    ! Purpose : Computes the normalised boundary-condition residual:
    !             error = sqrt( sum|EH_excit|^2 / normalize_E )
    !           Returns -2.0 if the error has stagnated (no improvement over
    !           the previous iteration), signalling best-achievable accuracy.
    !   error_only=.true.: evaluates without updating the scattered fields
    !             (used during post-processing checks).
    ! -----------------------------------------------------------------------
    function eval_error(this,N_calc,error_prev,error_only) result(error)
        type(Scatterer) :: this
        real(8) :: error,error_prev
        integer :: N_calc
        integer :: j,p,q
        logical :: error_only
        real(8),dimension(0:N_taylor) :: error_local,error_u,error_v



        error_u = 0.d0
        error_v = 0.d0
        do j=1,this%M

            if(this%testing_pt_status(j,1) == 2) then !! PMC boundary at this specific testing point
                error_u = error_u + abs(this%EH_excit(j+2*this%M,:))**2.d0
                error_v = error_v + abs(this%EH_excit(j+3*this%M,:))**2.d0
            else
                error_u = error_u + abs(this%EH_excit(j,:))**2.d0
                error_v = error_v + abs(this%EH_excit(j+this%M,:))**2.d0
            endif



        enddo


        if(this%Problem_Type == 2) then !! PMC Boundary
            error_local = (error_u+error_v)/this%normalize_H(0)
        else
            error_local = (error_u+error_v)/this%normalize_E(0)
        endif
        !        write(*,*) 'Error Array = ',error_local
        error = error_local(0)
        !        if(.not. error_only)then
        !        write(*,*) 'error =',error
!        if(This%Problem_Type == 1) then
!            write(*,*) 'error =',error,'        Error Proposed =',sum(abs(this%EH_excit(1:2*this%M,0))**2.d0)/(2*this%M)
!        elseif(This%Problem_Type == 2) then
!            write(*,*) 'error =',error,'        Error Proposed =',&
!            sum(abs(this%EH_excit((1+2*this%M):4*this%M,0))**2.d0)/(2*this%M)
!        elseif(This%Problem_Type == 3) then
!            write(*,*) 'error =',error,'        Error Proposed =',sum(abs(this%EH_excit(1:2*this%M,0))**2.d0)/(2*this%M)
!        elseif(This%Problem_Type == 4) then
!            write(*,*) 'error =',error,'        Error Proposed =',sum(abs(this%EH_excit(:,0))**2.d0)/size(this%EH_excit,1)
!        endif

        !        endif
        !        if(error > 1.d0) then
        !            error = 1.d0
        !        endif
        if(.not. error_only) then
            if(.not. this%converging) then !! the error is still above the limit TOL_internal
                p = 1
                do q=(this%N_curr-N_calc+1),this%N_curr
                    this%I_sources_total_current(q,:) = this%I_sources_total_current(q,:) + this%I_sources(q,:)
                    p = p+1
                enddo
                this%Es_current = this%Es_current + this%Es_c
                this%Hs_current = this%Hs_current + this%Hs_c

                this%Ed_current = this%Ed_current + this%Ed_c
                this%Hd_current = this%Hd_current + this%Hd_c

                if(RAS_Solution_Method == 1) then
                    if(this%itr_counter >= this%N_sub_groups/2) then
                        if(error < 1.d0) then
                            this%converging = .true.
                        else
                            error = -2.d0
                        endif
                    endif
                elseif(RAS_Solution_Method == 2) then
                    if(this%itr_counter >= this%N_sub_groups) then
                        if(error < 1.d0) then
                            this%converging = .true.
                        else
                            error = -2.d0
                        endif
                    endif
                else
                    if(this%itr_counter >= 2*this%N_sub_groups) then
                        if(error < 1.d0) then
                            this%converging = .true.
                        else
                            error = -2.d0
                        endif
                    endif
                endif
            else
                if((error - error_prev) > -5*TOL ) then !! problem still converging
!                if((error - error_prev) > -error_prev ) then !! problem still converging
                    p = 1
                    !                write(*,*) 'modifying currents'
                    do q=(this%N_curr-N_calc+1),this%N_curr
                        this%I_sources_total_current(q,:) = this%I_sources_total_current(q,:) + this%I_sources(q,:)
                        p = p+1
                    enddo
                    this%Es_current = this%Es_current + this%Es_c
                    this%Hs_current = this%Hs_current + this%Hs_c
                    !                if(this%Problem_Type == 4) then !! dielectric case
                    this%Ed_current = this%Ed_current + this%Ed_c
                    this%Hd_current = this%Hd_current + this%Hd_c
                !                endif
                else
                    write(*,*) 'started to diverge' !,error - error_prev,-5*TOL
                    error = -2.d0

                endif
            endif
        endif
        !        write(*,*) this%itr_counter
        this%itr_counter = this%itr_counter + 1
    end function eval_error


    ! -----------------------------------------------------------------------
    ! Subroutine: initialize_parameters
    ! Purpose   : Allocates and initialises all Scatterer member arrays based
    !             on the testing-point count M, source count N_max, and
    !             AWE order N_order. Also discretises the boundary, places
    !             contour points, sets material parameters (k_local, eta_local,
    !             ar_local, ak_local), and initialises IBC impedance matrices.
    !             Must be called once before any solve routine.
    ! -----------------------------------------------------------------------
    subroutine initialize_parameters(this,N_order,Scatterers_pointer)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)
        complex*16,dimension(0:N_order) :: H_u,H_v
        integer :: ii,N_order,jj
        complex*16 :: eta_ex

        this%Contour_Cnt = 1
        this%Contour_out_Cnt = 1
        this%N_curr = 0
        this%N = 0
        this%Es = 0.d0
        this%Es_current = 0.d0
        this%Hs = 0.d0
        this%Hs_current = 0.d0
        this%Hs_c = 0.d0
        this%Es_c = 0.d0
        this%E_bounce_current = 0.d0
        this%H_bounce_current = 0.d0
        this%E_bounce = 0.d0
        this%H_bounce = 0.d0
        this%Hd_c = 0.d0
        this%Ed_c = 0.d0
        do ii = 0,N_taylor
            this%I_sources(:,ii) = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
        enddo
        this%I_sources_total_current = this%I_sources
        this%I_sources_total = this%I_sources
        this%I_stat = 0
        this%converging = .false.
        this%group_index = 0
        this%active_region = 0



        do jj = 1,this%M
            if(this%testing_pt_status(jj,3) == 0) then
                eta_ex = eta1
            else
                eta_ex = Scatterers_pointer(this%testing_pt_status(jj,3))%eta_local
            endif

            this%EH_excit(jj,:) = this%E(jj,:)/eta_ex
            this%EH_excit(jj+this%M,:) = this%E(jj+this%M,:)/eta_ex
        enddo



        !        this%EH_excit = this%E/eta1
        !        do nn = 0,N_order
        !            this%EH_excit(1:this%M,nn)= this%norm_eta_u*this%E(1:this%M,nn)/eta1
        !            this%EH_excit((1+this%M):(2*this%M),nn)= this%norm_eta_v*this%E((1+this%M):(2*this%M),nn)/eta1
        !        enddo
        this%EH_excit((2*this%M+1):(4*this%M),:) = this%H





        do jj = 1,this%M
            if(this%testing_pt_status(jj,1) /= 3) then !! only IBC can do this loop
                cycle
            endif

            H_u = this%EH_excit(2*this%M+jj,:)
            H_v = this%EH_excit(3*this%M+jj,:)

            this%EH_excit(jj,:) = this%norm_eta_u(jj)*(this%EH_excit(jj,:) +  &
            this%eta_uu(jj)*H_v - this%eta_uv(jj)*H_u) !! E_u
            this%EH_excit(this%M+jj,:) = this%norm_eta_v(jj)*(this%EH_excit(this%M+jj,:) + &
            this%eta_vu(jj)*H_v - this%eta_vv(jj)*H_u) !! E_v

        enddo



        !        if(this%Problem_Type  == 3) then !! IBC
        !            !! then add the magnetic field part in case of having IBC
        !            allocate(H_u(this%M,0:N_taylor),H_v(this%M,0:N_taylor))
        !            H_u = this%H(1:this%M,:)
        !            H_v = this%H((this%M+1):(2*this%M),:)
        !            do ii = 1,this%M
        !
        !                this%EH_excit(ii,:) = this%norm_eta_u(ii)*(this%EH_excit(ii,:) + &
        !                this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:) )
        !                this%EH_excit(this%M+ii,:) = this%norm_eta_v(ii)*(this%EH_excit(this%M+ii,:) + &
        !                this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:) )
        !
        !
        !            enddo
        !
        !
        !            deallocate(H_u,H_v)
        !        endif

        !        if(this%Problem_Type == 1) then !! PEC
        !            this%EH_excit=this%E/eta1
        !        elseif(this%Problem_Type == 2) then !! PMC
        !            this%EH_excit=this%H
        !        elseif(this%Problem_Type  == 3) then !! IBC
        !            this%EH_excit(1:this%M,:)= this%norm_eta_u*this%E(1:this%M,:)/eta1
        !            this%EH_excit((1+this%M):(2*this%M),:)= this%norm_eta_v*this%E((1+this%M):(2*this%M),:)/eta1
        !            !! then add the magnetic field part in case of having IBC
        !            allocate(H_u(this%M,0:N_taylor),H_v(this%M,0:N_taylor))
        !            H_u = this%H(1:this%M,:)
        !            H_v = this%H((this%M+1):(2*this%M),:)
        !            do ii = 1,this%M
        !                this%EH_excit(ii,:) = this%EH_excit(ii,:) + &
        !                this%norm_eta_u*(this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:) )
        !                this%EH_excit(this%M+ii,:) = this%EH_excit(this%M+ii,:) + &
        !                this%norm_eta_v*(this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:) )
        !
        !
        !            !                this%EH_excit(ii,:) = this%EH_excit(ii,:) - &
        !            !                (this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:) )
        !            !                this%EH_excit(this%M+ii,:) = this%EH_excit(this%M+ii,:) - &
        !            !                (this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:) )
        !            enddo
        !            deallocate(H_u,H_v)
        !        elseif(this%Problem_Type == 4) then !! Dielectric Case
        !            this%EH_excit(1:(2*this%M),:)=this%E/eta1
        !            this%EH_excit((2*this%M+1):(4*this%M),:) = this%H
        !
        !        endif
        this%Ed_current = 0.d0
        this%Hd_current = 0.d0
    end subroutine initialize_parameters

    subroutine initialize_stage(this,N_order,E_inc,H_inc,is_first_time,Scatterers_pointer)
        logical :: is_first_time
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)
        type(Scatterer) :: this
        complex*16,allocatable,intent(inout) :: E_inc(:,:),H_inc(:,:)
        complex*16,dimension(0:N_order) :: H_u,H_v
        integer :: ii,N_order,jj
        complex*16 :: eta_ex
        this%N_curr = 0
        this%Es_current = 0.d0
        this%Hs_current = 0.d0
        this%Hs_c = 0.d0
        this%Es_c = 0.d0
        this%Hd_c = 0.d0
        this%Ed_c = 0.d0
        if(is_first_time) then
            this%N = 0
            this%Hs = 0.d0
            this%Es = 0.d0
            this%E_bounce_current = 0.d0
            this%H_bounce_current = 0.d0
            this%E_bounce = 0.d0
            this%H_bounce = 0.d0
            this%converging = .false.
            do ii = 0,N_taylor
                this%I_sources(:,ii) = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
            enddo
            this%I_sources_total_current = this%I_sources
            this%I_sources_total = this%I_sources
        endif
        do ii = 0,N_taylor
            this%I_sources(:,ii) = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
        enddo
        this%I_sources_total_current = this%I_sources
        this%group_index = 0
        this%E_current = E_inc
        this%H_current = H_inc

        do jj = 1,this%M
            if(this%testing_pt_status(jj,3) == 0) then
                eta_ex = eta1
            else
                eta_ex = Scatterers_pointer(this%testing_pt_status(jj,3))%eta_local
            endif

            this%EH_excit(jj,:) = E_inc(jj,:)/eta_ex
            this%EH_excit(jj+this%M,:) = E_inc(jj+this%M,:)/eta_ex
        enddo



        !        this%EH_excit = this%E/eta1
        !        do nn = 0,N_order
        !            this%EH_excit(1:this%M,nn)= this%norm_eta_u*this%E(1:this%M,nn)/eta1
        !            this%EH_excit((1+this%M):(2*this%M),nn)= this%norm_eta_v*this%E((1+this%M):(2*this%M),nn)/eta1
        !        enddo
        this%EH_excit((2*this%M+1):(4*this%M),:) = H_inc





        do jj = 1,this%M
            if(this%testing_pt_status(jj,1) /= 3) then !! only IBC can do this loop
                cycle
            endif

            H_u = this%EH_excit(2*this%M+jj,:)
            H_v = this%EH_excit(3*this%M+jj,:)

            this%EH_excit(jj,:) = this%norm_eta_u(jj)*(this%EH_excit(jj,:) +  &
            this%eta_uu(jj)*H_v - this%eta_uv(jj)*H_u) !! E_u
            this%EH_excit(this%M+jj,:) = this%norm_eta_v(jj)*(this%EH_excit(this%M+jj,:) + &
            this%eta_vu(jj)*H_v - this%eta_vv(jj)*H_u) !! E_v

        enddo


        !        this%EH_excit = E_inc/eta1
        !    !        do nn = 0,N_order
        !    !        this%EH_excit(1:this%M,nn)= this%norm_eta_u*E_inc(1:this%M,nn)/eta1
        !    !        this%EH_excit((1+this%M):(2*this%M),nn)= this%norm_eta_v*E_inc((1+this%M):(2*this%M),nn)/eta1
        !    !        enddo
        !    !        write(*,*) this%norm_eta_v
        !        this%EH_excit((2*this%M+1):(4*this%M),:) = H_inc
        !        if(This%Problem_Type == 3) then
        !            allocate(H_u(this%M,0:N_taylor),H_v(this%M,0:N_taylor))
        !            H_u = H_inc(1:this%M,:)
        !            H_v = H_inc((this%M+1):(2*this%M),:)
        !            do ii = 1,this%M
        !                this%EH_excit(ii,:) = this%norm_eta_u(ii)*(this%EH_excit(ii,:) + &
        !                this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:) )
        !                this%EH_excit(this%M+ii,:) = this%norm_eta_v(ii)*(this%EH_excit(this%M+ii,:) + &
        !                this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:) )
        !            enddo
        !            deallocate(H_u,H_v)
        !        endif

        !        if(this%Problem_Type == 1) then !! PEC
        !            this%EH_excit=this%E_current/eta1
        !        elseif(this%Problem_Type == 2) then !! PMC
        !            this%EH_excit=this%H_current
        !        elseif(this%Problem_Type  == 3) then !! IBC
        !            this%EH_excit(1:this%M,:)= this%norm_eta_u*this%E_current(1:this%M,:)/eta1
        !            this%EH_excit((1+this%M):(2*this%M),:)= this%norm_eta_v*this%E_current((1+this%M):(2*this%M),:)/eta1
        !            !! then add the magnetic field part in case of having IBC
        !            allocate(H_u(this%M,0:N_taylor),H_v(this%M,0:N_taylor))
        !            H_u = this%H_current(1:this%M,:)
        !            H_v = this%H_current((this%M+1):(2*this%M),:)
        !            do ii = 1,this%M
        !        !                this%EH_excit(ii,:) = this%EH_excit(ii,:) + &
        !        !                (this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:) )
        !        !                this%EH_excit(this%M+ii,:) = this%EH_excit(this%M+ii,:) + &
        !        !                (this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:) )
        !
        !                this%EH_excit(ii,:) = this%EH_excit(ii,:) + &
        !                this%norm_eta_u*(this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:) )
        !                this%EH_excit(this%M+ii,:) = this%EH_excit(this%M+ii,:) + &
        !                this%norm_eta_v*(this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:) )
        !
        !        !                this%EH_excit(ii,:) = this%EH_excit(ii,:) - &
        !        !                (this%eta_uu(ii)*H_v(ii,:) - this%eta_uv(ii)*H_u(ii,:) )
        !        !                this%EH_excit(this%M+ii,:) = this%EH_excit(this%M+ii,:) - &
        !        !                (this%eta_vu(ii)*H_v(ii,:) - this%eta_vv(ii)*H_u(ii,:) )
        !
        !            enddo
        !            deallocate(H_u,H_v)
        !        elseif(this%Problem_Type == 4) then !! Dielectric Case
        !            this%EH_excit(1:(2*this%M),:)=this%E_current/eta1
        !            this%EH_excit((2*this%M+1):(4*this%M),:) = this%H_current
        !
        !        endif
        this%Ed_current = 0.d0
        this%Hd_current = 0.d0
    !        !! re-evaluate the normalization values
    !        this%normalize_E = sum(abs(this%E_current)**2.d0)
    !        this%normalize_H = sum(abs(this%H_current)**2.d0)

    end subroutine initialize_stage

    ! -----------------------------------------------------------------------
    ! Subroutine: add_sources_rect
    ! Purpose   : Adds n_sources random interior sources to the Scatterer
    !             using a uniform random placement within the bounding box,
    !             constrained by the 'bound' parameter (fraction of lambda).
    !             Source orientations and types (electric/magnetic) are
    !             assigned based on Source_Placement and Contour_Sources_Type.
    ! -----------------------------------------------------------------------
    subroutine add_sources_rect(this,n_sources)
        ! add 'n_sources' to the current sources
        ! it can work if there is no sources, N should be set to 0 before usage for the first time
        type(Scatterer) :: this
        integer :: n_sources
        integer :: i,N_old,p
        real(8) :: x_s,y_s
        integer :: seed
        integer,dimension(3) :: time
        real(8),dimension(n_sources,3) :: one
        real(8),dimension(n_sources) :: two
        real(8) :: x_range,y_range,x_min,x_max,y_max,y_min,x_mid,y_mid
        integer :: i_rule,kk
        real(8) :: rand_bound_shift = 0.0d0

        if(n_sources < 1) then
            return
        endif

        N_old = this%N_curr
        this%N_curr = this%N_curr + n_sources
        !        write(*,*) 'current N =',N
        if(N_old == 0) then !! for the first time only
            call itime(time)
            seed = product(time)
            CALL init_random_seed()
            call srand(seed)
        endif
        !    write(*,*) 'N_old,N = ',N_old,N
        CALL RANDOM_NUMBER(one)
        CALL RANDOM_NUMBER(two)
        x_max = this%Str_Rules(this%N_rules,2)
        x_min = this%Str_Rules(1,1)
        x_mid = (x_max + x_min)/2.d0
        x_range = x_max - x_min
        !        write(*,*) x_range,x_min
        p =0
        do i=N_old+1,this%N_curr
            p = p+1

            !            r_s = (ibound)*one(p,1)
            !            th_s = 2.d0*pi*one(p,2)


            x_s = (x_range-2.d0*lambda*(1.d0-bound+rand_bound_shift))*(one(p,1)-0.5d0)+x_mid
            !! find i_rule
            do kk = 1,this%N_rules
                if(x_s < this%Str_Rules(kk,2)) then
                    i_rule = kk
                    exit
                endif
            enddo
            y_max = this%Str_Rules(kk,3)*x_s+this%Str_Rules(kk,4)
            y_min = this%Str_Rules(kk,5)*x_s+this%Str_Rules(kk,6)
            y_mid = (y_max + y_min)/2.d0
            y_range = y_max - y_min

            y_s = (y_range-2.d0*(1.d0-bound +rand_bound_shift)*lambda)*(one(p,2)-0.5d0)+y_mid


            this%source_pos(i,1) = x_s
            this%source_pos(i,2) = y_s
            this%source_pos(i,3) = 0.d0

            this%active_region(i) = this%region_status(2)

            !            write(*,*) this%source_pos(i,:)
            this%allowed_orientation(i) = 0
            if(two(p) < 0.5d0) then
                this%I_stat(i) = 1 !! electric sources
            else
                this%I_stat(i) = 2 !! magnetic sources
            endif
            write(12,*) this%source_pos(i,:),this%I_stat(i),this%active_region(i)
        enddo
    end subroutine add_sources_rect

    subroutine add_sources_outside_random(this,n_sources,n_added)
        type(Scatterer) :: this
        integer :: kk,n_sources,n_added
        integer :: i,N_old,p
        real(8) :: x_s,y_s
        integer :: seed
        integer,dimension(3) :: time
        real(8),dimension(n_sources,3) :: one
        real(8),dimension(n_sources) :: two
        real(8) :: x_range,y_range,y_max,y_min,x_min_str,x_max_str
        integer :: i_rule
        real(8) :: dist_max,dist_min

        if(n_sources < 1) then
            return
        endif

        N_old = this%N_curr
        this%N_curr = this%N_curr + n_sources
        n_added = n_sources
        !        write(*,*) 'current N =',N
        if(N_old == 0) then !! for the first time only
            call itime(time)
            seed = product(time)
            CALL init_random_seed()
            call srand(seed)
        endif
        !    write(*,*) 'N_old,N = ',N_old,N
        CALL RANDOM_NUMBER(one)
        CALL RANDOM_NUMBER(two)
        x_max_str = this%Str_Rules(this%N_rules,2)
        x_min_str = this%Str_Rules(1,1)
        x_range = this%x_bound_max - this%x_bound_min

        p =0
        do i=N_old+1,this%N_curr
            p = p+1

            x_s = x_range*one(p,1)+this%x_bound_min

            if((x_s > (x_min_str)) .and. (x_s < x_max_str)) then !! x_s is within the structure range
                do kk = 1,this%N_rules
                    if(x_s < this%Str_Rules(kk,2)) then
                        i_rule = kk
                        exit
                    endif
                enddo
                y_max = this%Str_Rules(kk,3)*x_s+this%Str_Rules(kk,4)
                y_min = this%Str_Rules(kk,5)*x_s+this%Str_Rules(kk,6)

                y_range = this%y_bound_max - this%y_bound_min
                y_s = y_range*one(p,2)+this%y_bound_min

                if((y_s > (y_min)) .and. (y_s < (y_max))) then

                    dist_max = y_max - y_s
                    dist_min = y_s - y_min

                    if(dist_max > dist_min) then !! closer to min point
                        y_range = y_min - this%y_bound_min
                        y_s = y_range*one(p,2)+this%y_bound_min - outside_bound*lambda
                    else !! closer to max point
                        y_range = this%y_bound_max - y_max
                        y_s = y_range*one(p,2)+y_max + outside_bound*lambda
                    endif

                endif


            else !! outside the structure's range
                if(abs(x_min_str - x_s) < outside_bound*lambda) then
                    x_s = x_s - outside_bound*lambda
                elseif(abs(x_s - x_max_str) < outside_bound*lambda) then
                    x_s = x_s + outside_bound*lambda
                endif

                y_range = this%y_bound_max - this%y_bound_min
                y_s = y_range*one(p,2)+this%y_bound_min


            endif
            this%source_pos(i,1) = x_s
            this%source_pos(i,2) = y_s
            this%source_pos(i,3) = 0.d0

            this%active_region(i) = this%Region_ID
            this%allowed_orientation(i) = 0
            !            write(*,*) this%source_pos(i,:)

            if(two(p) < 0.5d0) then
                this%I_stat(i) = 1 !! electric sources
            else
                this%I_stat(i) = 2 !! magnetic sources
            endif
            write(12,*) this%source_pos(i,:),this%I_stat(i),this%active_region(i)
        enddo
    end subroutine add_sources_outside_random


    ! -----------------------------------------------------------------------
    ! Subroutine: add_sources_contour
    ! Purpose   : Adds MAS-style sources on an interior contour parallel to
    !             the scatterer boundary (offset inward by 'bound'). Used
    !             when Source_Placement=1 or 3. For corners (is_corner=.true.)
    !             additional corner-point sources are inserted.
    ! -----------------------------------------------------------------------
    subroutine add_sources_contour(this,n_sources,n_added_ext,n_added_int,is_corner)
        type(Scatterer) :: this
        integer :: N_old,n_sources,n_added_ext,n_added_int,mm
        real(8),dimension(n_sources) :: two
        integer :: seed
        integer,dimension(3) :: time
        logical :: is_corner
        !        logical :: add_opposit_also = .true. !! a flag o add also the opposite type of the contour sources in another iteration to insure generality

        N_old = this%N_curr
        !        if((this%Contour_Cnt > this%N_con)) then
        !            if(.not. add_opposit_also) then
        !            n_added = 0
        !            return
        !            endif
        !        endif
        n_added_int = 0
        n_added_ext = 0
        if((this%Contour_Cnt > this%N_con)) then
            return
        endif


        !        if((this%Contour_Cnt + n_sources -1) > this%N_con) then
        !            n_added = this%N_con - this%Contour_Cnt + 1
        !            this%N_curr = this%N_curr + n_added
        !        else
        !            this%N_curr = this%N_curr + n_sources
        !            n_added = n_sources
        !        endif

        !        if((this%Contour_Cnt + n_sources -1) > 2*this%N_con) then

        mm = N_old+1

        if(N_old == 0) then !! for the first time only
            call itime(time)
            seed = product(time)
            CALL init_random_seed()
            call srand(seed)
        endif
        !    write(*,*) 'N_old,N = ',N_old,N
        CALL RANDOM_NUMBER(two)

        !        do ii = 1,this%N_con
        !            write(*,*) this%contour_points(ii,:)
        !        enddo
        !        write(*,*)

!        if((this%Contour_Cnt + n_sources -1) > this%N_con) then
!            !            n_added = 2*this%N_con - this%Contour_Cnt + 1
!            n_added_ext = this%N_con - this%Contour_Cnt + 1
!            this%N_curr = this%N_curr + n_added_ext
!        else
!            this%N_curr = this%N_curr + n_sources
!            n_added_ext = n_sources
!        endif



        do while(this%Contour_Cnt <= this%N_con)

            this%source_pos(mm,:) = this%contour_points(this%Contour_Cnt,:)



            this%I_stat(mm) = this%contour_points_type(this%Contour_Cnt)




            this%active_region(mm) = this%contour_point_active_region(this%Contour_Cnt)
            write(12,*) this%source_pos(mm,:),this%I_stat(mm),this%active_region(mm)
!            write(*,*) this%source_pos(mm,:),this%I_stat(mm),this%active_region(mm)



            if(is_corner) then
                this%allowed_orientation(mm) = this%contour_points_orientation(this%Contour_Cnt)
            else
                this%allowed_orientation(mm) = 0
            endif
            if(this%active_region(mm) == this%region_status(2)) then
                n_added_ext = n_added_ext+1
            else
                n_added_int = n_added_int+1
            endif

            this%N_curr = this%N_curr + 1
            this%Contour_Cnt = this%Contour_Cnt + 1
            mm = mm + 1
            if((n_added_ext + n_added_int) > n_sources ) then
                exit !! breaks the loop
            endif
        enddo


!        do kk=1,N_added_ext
!            this%source_pos(mm,:) = this%contour_points(this%Contour_Cnt,:)
!
!            this%I_stat(mm) = this%contour_points_type(this%Contour_Cnt)
!
!
!            this%active_region(mm) = this%contour_point_active_region(this%Contour_Cnt)
!            write(12,*) this%source_pos(mm,:),this%I_stat(mm),this%active_region(mm)
!
!            if(is_corner) then
!                this%allowed_orientation(mm) = this%contour_points_orientation(this%Contour_Cnt)
!            else
!                this%allowed_orientation(mm) = 0
!            endif
!            this%Contour_Cnt = this%Contour_Cnt + 1
!            mm = mm + 1
!
!        enddo


    end subroutine add_sources_contour



    subroutine add_sources_outside(this,n_sources,n_added)
        type(Scatterer) :: this
        integer :: N_old,kk,n_sources,n_added,mm
        !        real(8),dimension(n_sources) :: two
        integer :: seed
        integer,dimension(3) :: time

        if(n_sources < 1) then
            return
        endif

        N_old = this%N_curr
        this%N_curr = this%N_curr + n_sources

        n_added = n_sources

        mm = N_old+1

        if(N_old == 0) then !! for the first time only
            call itime(time)
            seed = product(time)
            CALL init_random_seed()
            call srand(seed)
        endif
        !        !    write(*,*) 'N_old,N = ',N_old,N
        !        CALL RANDOM_NUMBER(two)


        do kk=1,N_added
            this%source_pos(mm,:) = this%Inside_sources_dielectric(this%Contour_out_Cnt,:)
            !            write(12,*) this%source_pos(mm,:)
            this%Contour_out_Cnt = this%Contour_out_Cnt + 1
            !! Type Selection
            if(this%Contour_Sources_Type == 1) then !! Electric Sources
                this%I_stat(mm) = 1
            else!! Magnetic Sources
                this%I_stat(mm) = 2

            endif
            this%active_region(mm) = this%Region_ID
            !            write(*,*)  this%source_pos(mm,:)
            write(12,*) this%source_pos(mm,:),this%I_stat(mm),this%active_region(mm)
            mm = mm + 1
        enddo


    end subroutine add_sources_outside

    subroutine set_Z_Y_matrices(Z,Y,Z_f,Y_f,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt)
        complex*16,allocatable,intent(inout) :: Zzz(:,:,:),Zzt(:,:,:),Ztz(:,:,:),Ztt(:,:,:),&
        Yzt(:,:,:),Ytz(:,:,:),Ytt(:,:,:),Z(:,:,:),Y(:,:,:),Y_f(:,:,:),Z_f(:,:,:)
        integer :: n_rows,n_cols

        n_rows = size(Zzz,1)
        n_cols = size(Zzz,2)

        Y(1:n_rows,1:n_cols,:) = 0.d0
        Y(1:n_rows,(n_cols+1):2*n_cols,:) = Y(1:n_rows,(n_cols+1):2*n_cols,:) + Yzt
        Y((n_rows+1):2*n_rows,1:n_cols,:) = Y((n_rows+1):2*n_rows,1:n_cols,:)+ Ytz
        Y((n_rows+1):2*n_rows,(n_cols+1):2*n_cols,:) = Y((n_rows+1):2*n_rows,(n_cols+1):2*n_cols,:) + Ytt

        Y_f(1:n_rows,1:n_cols,:) = Y_f(1:n_rows,1:n_cols,:)+Ytz
        Y_f(1:n_rows,(n_cols+1):2*n_cols,:) = Y_f(1:n_rows,(n_cols+1):2*n_cols,:)+Ytt
        Y_f((n_rows+1):2*n_rows,1:n_cols,:) = 0.d0
        Y_f((n_rows+1):2*n_rows,(n_cols+1):2*n_cols,:) = Y_f((n_rows+1):2*n_rows,(n_cols+1):2*n_cols,:)-Yzt

        Z(1:n_rows,1:n_cols,:) = Z(1:n_rows,1:n_cols,:)+Zzz
        Z(1:n_rows,(n_cols+1):2*n_cols,:) = Z(1:n_rows,(n_cols+1):2*n_cols,:)+Zzt
        Z((n_rows+1):2*n_rows,1:n_cols,:) = Z((n_rows+1):2*n_rows,1:n_cols,:)+Ztz
        Z((n_rows+1):2*n_rows,(n_cols+1):2*n_cols,:) = Z((n_rows+1):2*n_rows,(n_cols+1):2*n_cols,:)+Ztt


        Z_f(1:n_rows,1:n_cols,:) = Z_f(1:n_rows,1:n_cols,:)+Ztz
        Z_f(1:n_rows,(n_cols+1):2*n_cols,:) = Z_f(1:n_rows,(n_cols+1):2*n_cols,:)+Ztt
        Z_f((n_rows+1):2*n_rows,1:n_cols,:) = Z_f((n_rows+1):2*n_rows,1:n_cols,:) -Zzz
        Z_f((n_rows+1):2*n_rows,(n_cols+1):2*n_cols,:) = Z_f((n_rows+1):2*n_rows,(n_cols+1):2*n_cols,:) -Zzt


    end subroutine set_Z_Y_matrices

    function eval_surface_current_MoM_once(Scat_big,N_order,these) result(errors)
        type(Scatterer) :: Scat_big
        type(Scatterer),allocatable,intent(inout) :: these(:)
        !        type(Scatterer),allocatable,dimension(:) :: these
        complex*16,allocatable,dimension(:,:,:) :: Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt
        complex*16,allocatable,dimension(:,:,:) :: Z,Y_f,Y,Z_f,ZY,ZY_l
        complex*16,allocatable,dimension(:,:) :: EH_l,EH,E,H,I_j
        integer :: N_order,ii,jj,nn
        real(8),allocatable,dimension(:) :: errors
        complex*16 :: ar_ex,ar_in,ak_ex,ak_in,eta_ex,eta_in,eta_used,eta_med,eta_norm,eta_norm_s
        integer :: BC_type_t,inner_region_id_t,outer_region_id_t
        integer :: BC_type_s,inner_region_id_s,outer_region_id_s
        integer :: n_unknowns,ind_beg,j_cnt,divisor,divisor_t,divisor_s
        integer,allocatable,dimension(:)::n_unknowns_array
        type surface_impedance
            complex*16,allocatable,dimension(:,:) :: eta_s,eta_f
        end type surface_impedance
        type(surface_impedance),allocatable,dimension(:) :: anis_imp

        allocate(anis_imp(number_of_scatterers))
        !        write(*,*) 'WARNING: MoM solution once does not support impedance boundary condition'
        allocate(errors(number_of_scatterers))
        if(allocated(Scat_big%Lm_collect)) then
            deallocate(Scat_big%Nm_collect,Scat_big%Lm_collect,Scat_big%Dm_collect)
        endif
        if(allocated(Scat_big%Y_collect)) then
            deallocate(Scat_big%Y_collect)
        endif
        if(allocated(Scat_big%Hd_c)) then
            deallocate(Scat_big%Hd_c,Scat_big%Ed_c)
        endif
        if(allocated(Scat_big%Hd_c)) then
            deallocate(Scat_big%Es_c,Scat_big%Hs_c)
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Deallocation and allocataion of discretization parameters
        do ii = 1,number_of_scatterers
            deallocate(these(ii)%tang_u,these(ii)%tang_v,these(ii)%norm_v,these(ii)%delta_n)
            deallocate(these(ii)%testing_pt,these(ii)%testing_pt_MoM)
            deallocate(these(ii)%testing_pt_status)

            if(allocated(these(ii)%contour_points)) then
            deallocate(these(ii)%contour_points,these(ii)%contour_points_type)
            endif

            if(allocated(these(ii)%contour_point_active_region)) then
            deallocate(these(ii)%contour_point_active_region)
            endif

            if(allocated(these(ii)%contour_points_orientation)) then
                deallocate(these(ii)%contour_points_orientation)
            endif

            if(allocated(these(ii)%Inside_sources_dielectric)) then
                deallocate(these(ii)%Inside_sources_dielectric)
            endif



            if(scatterers_input_method(ii) == 1) then
                call discretize_superquad(these(ii),MoM_Samples_per_wavelength,these)
            elseif(scatterers_input_method(ii) == 2) then
                call discretize_scatterer_file(these(ii),MoM_Samples_per_wavelength,these)
            endif

            if(allocated(these(ii)%eta_vv)) then
            deallocate(these(ii)%eta_vv,these(ii)%eta_uu,these(ii)%eta_uv,these(ii)%eta_vu)
            endif

            deallocate(these(ii)%norm_eta_u,these(ii)%norm_eta_v)
            if(these(ii)%Problem_Type == 3) then !! IBC
                call set_IBC_impedance_matrices(these(ii))
            else
                allocate(these(ii)%norm_eta_u(these(ii)%M),these(ii)%norm_eta_v(these(ii)%M))
                these(ii)%norm_eta_u = 1.d0
                these(ii)%norm_eta_v = 1.d0
            endif
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        !        do ii = 1,number_of_scatterers
        !            deallocate(these(ii)%tang_u,these(ii)%tang_v,these(ii)%norm_v)
        !            deallocate(these(ii)%testing_pt_MoM,these(ii)%testing_pt,these(ii)%delta_n)
        !            call discretize_rect(these(ii),these(ii)%segments_points,MoM_Samples_per_wavelength,&
        !            these,non_uniform_type_read,N_o_read,A_segmentation_read)
        !        enddo

        !        call assign_scatterers(these_in,these)

        n_unknowns = 0

        allocate(n_unknowns_array(number_of_scatterers))
        n_unknowns_array = 0
        do ii = 1,number_of_scatterers

            do jj=1,these(ii)%M
                if(these(ii)%testing_pt_status(jj,1) == -1) then
                    cycle !! cancelled testing point
                endif
                if(these(ii)%Problem_Type == 4) then
                    n_unknowns = n_unknowns + 4
                    n_unknowns_array(ii) = n_unknowns_array(ii)+4
                else
                    n_unknowns = n_unknowns + 2
                    n_unknowns_array(ii) = n_unknowns_array(ii)+2
                endif
            enddo

        !            if(these(ii)%Problem_Type == 3) then !! IBC boundary
        !                call eval_IBC_impedance_matrices(these(ii),eta_s,eta_f)
        !            endif
        enddo

        allocate(ZY(n_unknowns,n_unknowns,0:N_order),EH(n_unknowns,0:N_order),I_j(n_unknowns,0:N_order))

        do ii = 1,number_of_scatterers
            BC_type_t = these(ii)%Problem_Type

            if(BC_type_t /= 3) then !! only IBC is allowed to be allocated
                cycle
            endif
            call eval_IBC_impedance_matrices(these(ii),anis_imp(ii)%eta_s,anis_imp(ii)%eta_f)

        enddo


        do ii=1,number_of_scatterers
            if(these(ii)%Problem_Type == 4) then
                divisor = 4
                divisor_t = 2
                allocate(E(n_unknowns_array(ii)/2,0:N_order),H(n_unknowns_array(ii)/2,0:N_order))
            else
                divisor = 2
                divisor_t = 1
                allocate(E(n_unknowns_array(ii),0:N_order),H(n_unknowns_array(ii),0:N_order))
            endif

            call set_excitation_v(these(ii),these(ii)%testing_pt_MoM,these(ii)%Ei,these(ii)%Hi)
            allocate(these(ii)%Eu_MoM(these(ii)%M,0:N_order),these(ii)%Ev_MoM(these(ii)%M,0:N_order),&
            these(ii)%Hu_MoM(these(ii)%M,0:N_order),these(ii)%Hv_MoM(these(ii)%M,0:N_order))
            call set_excitation_vectors(these(ii),these(ii)%Eu_MoM,these(ii)%Ev_MoM,these(ii)%Hu_MoM,these(ii)%Hv_MoM)
            allocate(these(ii)%I_v(these(ii)%M,0:N_order),these(ii)%I_u(these(ii)%M,0:N_order),&
            these(ii)%M_v(these(ii)%M,0:N_order),these(ii)%M_u(these(ii)%M,0:N_order))

            these(ii)%I_v = 0.d0
            these(ii)%I_u = 0.d0
            these(ii)%M_v = 0.d0
            these(ii)%M_u = 0.d0

            BC_type_t = these(ii)%Problem_Type
            inner_region_id_t = these(ii)%region_status(3)
            outer_region_id_t = these(ii)%region_status(2)
            if(outer_region_id_t == 0) then
                eta_ex = eta1
            else
                eta_ex = these(outer_region_id_t)%eta_local
            endif



            allocate(EH_l(n_unknowns_array(ii),0:N_order))

            j_cnt = 1

            do jj = 1,these(ii)%M
                if(these(ii)%testing_pt_status(jj,1) == -1) then
                    cycle
                endif

                E(j_cnt,:) = these(ii)%Ev_MoM(jj,:)
                E(j_cnt+n_unknowns_array(ii)/divisor,:) = these(ii)%Eu_MoM(jj,:)
                H(j_cnt,:) = these(ii)%Hu_MoM(jj,:)
                H(j_cnt+n_unknowns_array(ii)/divisor,:) = -these(ii)%Hv_MoM(jj,:)

                j_cnt = j_cnt + 1
            enddo
!            write(*,*) sum(abs(E(:,0)))
            if(BC_type_t == 1) then !! PEC
                if(FormulationType == 1)then !! EFIE
                    EH_l = E
                elseif(FormulationType == 2)then !! MFIE
                    EH_l = H
                else !! CFIE
                    EH_l = E + eta_ex*H
                endif
            elseif(BC_type_t == 2) then !! PMC
                if(FormulationType == 1) then !! EFIE
                    EH_l = -E
                elseif(FormulationType == 2) then !! MFIE
                    EH_l = H
                else !! CFIE
                    EH_l = H - E/eta_ex
                endif
            elseif(BC_type_t == 3) then !! IBC
                if(FormulationType == 1) then !! EFIE
                    EH_l = E
                elseif(FormulationType == 2) then !! MFIE
                    EH_l = H
                elseif(FormulationType == 3) then !! CFIE
                    EH_l = E + eta_ex*H
                else
                    EH_l = E - eta_ex*matmul(anis_imp(ii)%eta_s,H)
                endif
            else
                EH_l(1:2*these(ii)%M,:) = E/eta_ex
                EH_l((1+2*these(ii)%M):4*these(ii)%M,:) = H
            endif
            if(ii ==1) then
                EH(1:n_unknowns_array(ii),:) = EH_l
            else
                EH((1+sum(n_unknowns_array(1:(ii-1)))):sum(n_unknowns_array(1:ii)),:) = EH_l
            endif
            !            write(*,*) ii,sum(abs(EH_l))

            do jj = 1,number_of_scatterers
                if(these(jj)%Problem_Type == 4) then
                    divisor_s = 2
                else
                    divisor_s = 1
                endif
                allocate(ZY_l(n_unknowns_array(ii),n_unknowns_array(jj),0:N_order))
                allocate(Z(n_unknowns_array(ii)/divisor_t,n_unknowns_array(jj)/divisor_s,0:N_order),&
                Y_f(n_unknowns_array(ii)/divisor_t,n_unknowns_array(jj)/divisor_s,0:N_order))
                allocate(Zzz(n_unknowns_array(ii)/(2*divisor_t),n_unknowns_array(jj)/(2*divisor_s),0:N_order),&
                Zzt(n_unknowns_array(ii)/(2*divisor_t),n_unknowns_array(jj)/(2*divisor_s),0:N_order),&
                Ztz(n_unknowns_array(ii)/(2*divisor_t),n_unknowns_array(jj)/(2*divisor_s),0:N_order),&
                Ztt(n_unknowns_array(ii)/(2*divisor_t),n_unknowns_array(jj)/(2*divisor_s),0:N_order),&
                Yzt(n_unknowns_array(ii)/(2*divisor_t),n_unknowns_array(jj)/(2*divisor_s),0:N_order),&
                Ytz(n_unknowns_array(ii)/(2*divisor_t),n_unknowns_array(jj)/(2*divisor_s),0:N_order),&
                Ytt(n_unknowns_array(ii)/(2*divisor_t),n_unknowns_array(jj)/(2*divisor_s),0:N_order))
                allocate(Y(n_unknowns_array(ii)/divisor_t,n_unknowns_array(jj)/divisor_s,0:N_order),&
                Z_f(n_unknowns_array(ii)/divisor_t,n_unknowns_array(jj)/divisor_s,0:N_order))
                Z = 0.d0
                Y = 0.d0
                Z_f = 0.d0
                Y_f = 0.d0
                Zzz = 0.d0
                Ztz = 0.d0
                Zzt= 0.d0
                Ztt= 0.d0
                Yzt= 0.d0
                Ytz= 0.d0
                Ytt= 0.d0

                BC_type_s = these(jj)%Problem_Type
                inner_region_id_s = these(jj)%region_status(3)
                outer_region_id_s = these(jj)%region_status(2)

                ar_in = these(inner_region_id_s)%ar_local
                ak_in = these(inner_region_id_s)%ak_local
                eta_in = these(inner_region_id_s)%eta_local
                if(outer_region_id_s == 0) then
                    ar_ex = ar
                    ak_ex = ak
                    eta_ex = eta1
                else
                    ar_ex = these(outer_region_id_s)%ar_local
                    ak_ex = these(outer_region_id_s)%ak_local
                    eta_ex = these(outer_region_id_s)%eta_local
                endif

                if(ii == jj) then

                    call set_MoM_matrices_separate(these(ii),Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az,.true.,these)
                    if(Wideband_Type > 0) then
                        call set_MoM_AWE_matrices(these(ii),N_order,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ar_ex,az,eta_ex,.true.)
                    endif

                    eta_used = eta_ex
                    eta_med = eta_ex
                    eta_norm = eta_ex
                    eta_norm_s = eta_norm
                    call set_Z_Y_matrices(Z,Y,Z_f,Y_f,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt)
!                    write(*,*) 'determinants =',abs(det(Zzz(:,:,0))),abs(det(Ztt(:,:,0))),abs(det(Ztz(:,:,0))),abs(det(Zzt(:,:,0)))
!                    write(*,*) 'determinant =',det(Z_f(:,:,0)),det(Z(:,:,0)),det(Y_f(:,:,0))
                    if(BC_type_s == 4) then
                        call set_MoM_matrices_separate(these(ii),Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az,.false.,these)
                        if(Wideband_Type > 0) then
                            call set_MoM_AWE_matrices(these(ii),N_order,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,&
                            ar_in,az,eta_in,.false.)
                        endif

                        Z_f = ((eta_in/eta_ex)*(eta_in/eta_ex))*Z_f

                        call set_Z_Y_matrices(Z,Y,Z_f,Y_f,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt)
                        Z_f = Z_f/(eta_in*eta_in)
!                        write(*,*) 'determinant =',det(Z_f(:,:,0)),det(Z(:,:,0)),det(Y_f(:,:,0))
                    endif

                elseif((inner_region_id_s == outer_region_id_t) .or. (inner_region_id_t == outer_region_id_s) .or.&
                (outer_region_id_t == outer_region_id_s)) then
                    !! in this case --- matrices are to be evaluated
                    call set_MoM_bouncing_matrices(these(jj),these(ii),Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az,these)
                    call set_Z_Y_matrices(Z,Y,Z_f,Y_f,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt)
                    if((outer_region_id_s == inner_region_id_t)) then

                        if(outer_region_id_t == 0) then
                            eta_used = eta1
                        else
                            eta_used = these(outer_region_id_t)%eta_local
                        endif
                        eta_med = eta_ex
                        eta_norm = eta_used
                        if(outer_region_id_s == 0) then
                            eta_norm_s = eta1
                        else
                            eta_norm_s = these(outer_region_id_s)%eta_local
                        endif

                        Z = -Z
                        Y = -Y
                        Y_f = -Y_f
                        Z_f = -Z_f
                    !                        write(*,*) ii,jj,'outer_region_id_s == inner_region_id_t'
                    elseif((outer_region_id_t == outer_region_id_s)) then

                        eta_used = eta_ex
                        eta_med = eta_ex
                        eta_norm = eta_ex
                        if(outer_region_id_s == 0) then
                            eta_norm_s = eta1
                        else
                            eta_norm_s = these(outer_region_id_s)%eta_local
                        endif
                    !                        write(*,*) ii,jj,'outer_region_id_t == outer_region_id_s'
                    elseif(inner_region_id_s == outer_region_id_t) then

                        eta_used = eta_ex
                        eta_med = eta_in
                        if(outer_region_id_t == 0) then
                            eta_norm = eta1
                        else
                            eta_norm = these(outer_region_id_t)%eta_local
                        endif
                        if(outer_region_id_s == 0) then
                            eta_norm_s = eta1
                        else
                            eta_norm_s = these(outer_region_id_s)%eta_local
                        endif
                        !                        write(*,*) ii,jj,'inner_region_id_s == outer_region_id_t'
                        Z = -Z
                        Y = -Y
                        Y_f = -Y_f
                        Z_f = -Z_f
                    !                        write(*,*) ii,jj,'inner_region_id_s == outer_region_id_t'
                    endif
                endif
                deallocate(Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt)
                !!!!!!!
                !                write(*,*) ii,jj,'eta_used',eta_used
                if(BC_type_t == 1) then !! PEC
                    if(BC_type_s == 1) then
                        if(FormulationType == 1)then !! EFIE
                            ZY_l = Z
                        elseif(FormulationType == 2)then !! MFIE
                            ZY_l = Y_f
                        else !! CFIE
                            ZY_l = Z + eta_used*Y_f
                        endif
                    elseif(BC_type_s == 2) then
                        if(FormulationType == 1)then !! EFIE
                            ZY_l = -Y
                        elseif(FormulationType == 2)then !! MFIE
                            ZY_l = Z_f/eta_used**2.d0
                        else !! CFIE
                            ZY_l = -Y + Z_f/eta_used
                        endif
                    elseif(BC_type_s == 3) then
                        if(FormulationType == 1) then !! EFIE
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Z(:,:,nn) + matmul(Y(:,:,nn),anis_imp(jj)%eta_f)*eta_used
                            enddo
                        elseif(FormulationType == 2) then !! MFIE
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Y_f(:,:,nn) - matmul(Z_f(:,:,nn),anis_imp(jj)%eta_f)/eta_used
                            enddo
                        else !! CFIE
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Z(:,:,nn) + matmul(Y(:,:,nn),anis_imp(jj)%eta_f)*eta_used + &
                                eta_used*Y_f(:,:,nn) - matmul(Z_f(:,:,nn),anis_imp(jj)%eta_f)
                            enddo
                        endif
                    else
                        !!!!!!Dielectric
                        !                        write(*,*) 'eta_used/eta1 ',eta_used/eta1,'eta_med/eta1 ',eta_med/eta1
                        if(FormulationType == 1) then !! EFIE
                            ind_beg = 1
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = Z
                            ind_beg = ind_beg+2*these(jj)%M
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = -eta_used*Y
                        elseif(FormulationType == 2) then !! MFIE
                            ind_beg = 1
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = Y_f
                            ind_beg = ind_beg+2*these(jj)%M
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = eta_used*(Z_f/eta_med**2.d0)
                        else
                            ind_beg = 1
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = Z + eta_med*Y_f
                            ind_beg = ind_beg+2*these(jj)%M
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = -eta_used*Y + (eta_used/eta_med)*Z_f
                        endif
                    endif
                elseif(BC_type_t == 2) then !! PMC
                    if(BC_type_s == 1) then
                        if(FormulationType == 1) then !! EFIE
                            ZY_l = -Z
                        elseif(FormulationType == 2) then !! MFIE
                            ZY_l = Y_f
                        else !! CFIE
                            ZY_l = Y_f  - Z/eta_used
                        endif

                    elseif(BC_type_s == 2) then
                        if(FormulationType == 1) then !! EFIE
                            ZY_l = Y
                        elseif(FormulationType == 2) then !! MFIE
                            ZY_l = Z_f/eta_used**2.d0
                        else !! CFIE
                            ZY_l = Z_f/eta_used**2.d0 + Y/eta_ex
                        endif
                    elseif(BC_type_s == 3) then
                        if(FormulationType == 1) then !! EFIE
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = -Z(:,:,nn) - matmul(Y(:,:,nn),anis_imp(jj)%eta_f)*eta_used
                            enddo
                        elseif(FormulationType == 2) then !! MFIE
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Y_f(:,:,nn) - matmul(Z_f(:,:,nn),anis_imp(jj)%eta_f)/eta_used
                            enddo
                        else !! CFIE
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Y_f(:,:,nn) - matmul(Z_f(:,:,nn),anis_imp(jj)%eta_f)/eta_used -&
                                Z(:,:,nn)/eta_used - matmul(Y(:,:,nn),anis_imp(jj)%eta_f)
                            enddo
                        endif
                    else
                        !!!!!!Dielectric
                        !                        write(*,*) 'eta_used/eta1 ',eta_used/eta1,'eta_med/eta1 ',eta_med/eta1
                        if(FormulationType == 1) then !! EFIE
                            ind_beg = 1
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = -Z
                            ind_beg = ind_beg+2*these(jj)%M
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = eta_used*Y
                        elseif(FormulationType == 2) then !! MFIE
                            ind_beg = 1
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = Y_f
                            ind_beg = ind_beg+2*these(jj)%M
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = eta_used*Z_f/eta_med**2.d0
                        else
                            ind_beg = 1
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = -Z/eta_med + Y_f
                            ind_beg = ind_beg+2*these(jj)%M
                            ZY_l(:,(ind_beg):(ind_beg-1+2*these(jj)%M),:) = (eta_used/eta_med)*Y + eta_used*Z_f/eta_med**2.d0
                        endif
                    endif
                elseif(BC_type_t == 3) then !! IBC
                    if(BC_type_s == 1) then !! source PEC
                        if(FormulationType == 1) then
                            ZY_l = Z
                        elseif(FormulationType == 2) then
                            ZY_l = Y_f
                        elseif(FormulationType == 3) then
                            ZY_l = Z+eta_used*Y_f
                        else
                            do nn =0,N_order
                                ZY_l(:,:,nn) = Z(:,:,nn) - eta_used*matmul(anis_imp(ii)%eta_s,Y_f(:,:,nn))
                            enddo
                        endif

                    elseif(BC_type_s == 2) then !! source PMC
                        if(FormulationType == 1) then
                            ZY_l = -Y
                        elseif(FormulationType == 2) then
                            ZY_l = Z_f/eta_used**2.d0
                        elseif(FormulationType == 3) then
                            ZY_l = -Y + Z_f/eta_used
                        else
                            do nn =0,N_order
                                ZY_l(:,:,nn) = -Y(:,:,nn) - matmul(anis_imp(ii)%eta_s,Z_f(:,:,nn))/eta_used
                            enddo
                        endif
                    elseif(BC_type_s == 3) then !! source IBC
                        if(FormulationType == 1) then
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Z(:,:,nn) + eta_used*matmul(Y(:,:,nn),anis_imp(jj)%eta_f)
                            enddo
                        elseif(FormulationType == 2) then
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Y_f(:,:,nn) - matmul(Z_f(:,:,nn),anis_imp(jj)%eta_f)/eta_used
                            enddo
                        elseif(FormulationType == 3) then
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Z(:,:,nn) + matmul(Y(:,:,nn),anis_imp(jj)%eta_f)*eta_used + &
                                eta_used*Y_f(:,:,nn) - matmul(Z_f(:,:,nn),anis_imp(jj)%eta_f)
                            enddo
                        else
                            do nn = 0,N_order
                                ZY_l(:,:,nn) = Z(:,:,nn) + matmul(Y(:,:,nn),anis_imp(jj)%eta_f)*eta_used - &
                                eta_used*matmul(anis_imp(ii)%eta_s,(Y_f(:,:,nn) - matmul(Z_f(:,:,nn),anis_imp(jj)%eta_f)/eta_used))
                            enddo
                        endif



                    else !! dielectric boundary source
                        !                        write(*,*) 'eta_used/eta1 ',eta_used/eta1,'eta_med/eta1 ',eta_med/eta1
                        if(FormulationType == 1) then
                            ZY_l(1:2*these(ii)%M,1:2*these(jj)%M,:) = Z
                            ZY_l(1:2*these(ii)%M,(1+2*these(jj)%M):4*these(jj)%M,:) = -eta_used*Y
                        elseif(FormulationType == 2) then
                            ZY_l(1:2*these(ii)%M,1:2*these(jj)%M,:) = Y_f
                            ZY_l(1:2*these(ii)%M,(1+2*these(jj)%M):4*these(jj)%M,:) = eta_used*Z_f/eta_med**2.d0
                        elseif(FormulationType == 3) then
                            ZY_l(1:2*these(ii)%M,1:2*these(jj)%M,:) = Z + eta_med*Y_f
                            ZY_l(1:2*these(ii)%M,(1+2*these(jj)%M):4*these(jj)%M,:) = -eta_used*Y + (eta_used/eta_med)*Z_f
                        else
                            do nn = 0,N_order
                                ZY_l(1:2*these(ii)%M,1:2*these(jj)%M,nn) = Z(:,:,nn) - &
                                eta_med*matmul(anis_imp(ii)%eta_s,Y_f(:,:,nn))
                                ZY_l(1:2*these(ii)%M,(1+2*these(jj)%M):4*these(jj)%M,nn) = &
                                -eta_used*Y(:,:,nn) - eta_used/eta_med*matmul(anis_imp(ii)%eta_s,Z_f(:,:,nn))
                            enddo
                        endif


                    endif
                else !! Dielectric boundary
                    !                    write(*,*) 'eta_norm_s/eta1 ',eta_norm_s/eta1,'eta_med/eta1 ',eta_med/eta1,'eta_norm/eta1 ',eta_norm/eta1
                    if(BC_type_s == 4) then
                        ZY_l(1:2*these(ii)%M,1:2*these(jj)%M,:) = Z/eta_norm
                        ZY_l(1:2*these(ii)%M,(1+2*these(jj)%M):4*these(jj)%M,:) = -eta_norm_s/eta_norm*Y !
                        ZY_l((1+2*these(ii)%M):4*these(ii)%M,1:2*these(jj)%M,:) = Y_f
                        if(ii == jj) then
                            ZY_l((1+2*these(ii)%M):4*these(ii)%M,(1+2*these(jj)%M):4*these(jj)%M,:) = eta_norm*Z_f
                        else
                            ZY_l((1+2*these(ii)%M):4*these(ii)%M,(1+2*these(jj)%M):4*these(jj)%M,:) = eta_norm_s*Z_f/eta_med**2.d0
                        endif
                    elseif(BC_type_s == 1) then
                        ZY_l(1:2*these(ii)%M,:,:) = Z/eta_used
                        ZY_l((1+2*these(ii)%M):4*these(ii)%M,:,:) = Y_f
                    elseif(BC_type_s == 2) then
                        ZY_l(1:2*these(ii)%M,:,:) = -Y/eta_used
                        ZY_l((1+2*these(ii)%M):4*these(ii)%M,:,:) = Z_f/eta_med**2.d0
                    elseif(BC_type_s == 3) then
                        do nn = 0,N_order
                            ZY_l(1:2*these(ii)%M,:,nn) = Z(:,:,nn)/eta_used + eta_med/eta_used*matmul(Y(:,:,nn),anis_imp(jj)%eta_f)
                            ZY_l((1+2*these(ii)%M):4*these(ii)%M,:,nn) = Y_f(:,:,nn) - &
                            matmul(Z_f(:,:,nn),anis_imp(jj)%eta_f)/eta_med
                        enddo
                    endif





                endif

                if((ii == jj)) then
                    if(ii == 1) then
                        ZY(1:n_unknowns_array(1),1:n_unknowns_array(1),:) = ZY_l
                    else
                        ZY((1+sum(n_unknowns_array(1:(ii-1)))):sum(n_unknowns_array(1:ii)),&
                        (1+sum(n_unknowns_array(1:(ii-1)))):sum(n_unknowns_array(1:ii)) ,:) = ZY_l
                    endif
                else
                    if(ii == 1) then
                        ZY(1:n_unknowns_array(1),&
                        (1+sum(n_unknowns_array(1:(jj-1)))):sum(n_unknowns_array(1:jj)) ,:) = ZY_l
                    elseif(jj == 1) then
                        ZY((1+sum(n_unknowns_array(1:(ii-1)))):sum(n_unknowns_array(1:ii)),&
                        1:n_unknowns_array(1) ,:) = ZY_l
                    else
                        ZY((1+sum(n_unknowns_array(1:(ii-1)))):sum(n_unknowns_array(1:ii)),&
                        (1+sum(n_unknowns_array(1:(jj-1)))):sum(n_unknowns_array(1:jj)) ,:) = ZY_l
                    endif

                endif

                deallocate(Z,Y_f,Y,Z_f,ZY_l)
            enddo
            deallocate(E,H,EH_l)
        enddo
!        call write_matrix(ZY(:,:,0),'ZY.dat')
!        call write_matrix_MATLAB(ZY(:,:,0),'ZY_MATLAB')
!        stop
!        write(*,*) 'determinant =',det(ZY(:,:,0))
        allocate(these(1)%MoM_matrix_inverted(n_unknowns,n_unknowns))


!        write(*,*) 'determinant =',det(ZY(:,:,0))
        these(1)%MoM_matrix_inverted = Mat_Inv(ZY(:,:,0))
        I_j(:,0) = matmul(these(1)%MoM_matrix_inverted,EH(:,0))

        do ii = 1,number_of_scatterers
            if(ii == 1) then
                ind_beg = 1
            else
                ind_beg = sum(n_unknowns_array(1:(ii-1)))+1
            endif
            if(these(ii)%Problem_Type == 1) then
                these(ii)%I_v = I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
                ind_beg = ind_beg+these(ii)%M
                these(ii)%I_u = I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
            elseif(these(ii)%Problem_Type == 2) then
                these(ii)%M_v = I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
                ind_beg = ind_beg+these(ii)%M
                these(ii)%M_u = I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
            elseif(these(ii)%Problem_Type == 3) then
                these(ii)%I_v = I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
                ind_beg = ind_beg+these(ii)%M
                these(ii)%I_u = I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
                do jj = 1,these(ii)%M

                    outer_region_id_t = these(ii)%region_status(2)
                    if(outer_region_id_t == 0) then
                        eta_ex = eta1
                    else
                        eta_ex = these(outer_region_id_t)%eta_local
                    endif

                    do nn = 0,N_order
                        these(ii)%M_u(jj,nn) = eta_ex*(these(ii)%eta_vv(jj)*these(ii)%I_v(jj,nn) +&
                        these(ii)%eta_vu(jj)*these(ii)%I_u(jj,nn)) !! M_t
                        these(ii)%M_v(jj,nn) = -eta_ex*(these(ii)%eta_uu(jj)*these(ii)%I_u(jj,nn) +&
                        these(ii)%eta_uv(jj)*these(ii)%I_v(jj,nn)) !! M_z
                    enddo
                enddo
            else
                outer_region_id_t = these(ii)%region_status(2)
                if(outer_region_id_t == 0) then
                    eta_ex = eta1
                else
                    eta_ex = these(outer_region_id_t)%eta_local
                endif

                these(ii)%I_v = I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
                ind_beg = ind_beg+these(ii)%M
                these(ii)%I_u = I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
                ind_beg = ind_beg+these(ii)%M
                these(ii)%M_v = eta_ex*I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
                ind_beg = ind_beg+these(ii)%M
                these(ii)%M_u = eta_ex*I_j(ind_beg:(ind_beg+these(ii)%M-1),:)
            endif
        enddo

        if(MoM_activation_flag == 1) then !! it is required to Compare farfields and currents between MoM and RAS
!            call compare_far_field_scatterers(these,180)
            call plot_far_field_MoM(these,180)
            !! plot currents
            do ii = 1,number_of_scatterers
!                errors(ii) = error_plot_current_comparison_with_MoM(these(ii))
                errors(ii) = -1.d0
                call plot_current_MoM(these(ii))
                if(allocated(anis_imp(ii)%eta_s)) then
                    deallocate(anis_imp(ii)%eta_s,anis_imp(ii)%eta_f)
                endif
            enddo

        endif
        deallocate(EH,ZY)


        deallocate(anis_imp)
    end function eval_surface_current_MoM_once

    function eval_surface_current_error_multiscatterer_MoM(these,N_order) result(errors)
        type(Scatterer),allocatable,intent(inout) :: these(:)
        integer :: N_order
        real(8),allocatable,dimension(:) :: errors,error_bouncing
        complex*16,allocatable,dimension(:,:,:) :: Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt
        complex*16,allocatable,dimension(:,:) :: Z_l,Y_l
        type :: Mutual_Matrix
            integer :: Source_ID,Dist_ID
            integer :: m_rows,n_cols
            complex*16,allocatable,dimension(:,:) :: ZY
            complex*16,allocatable,dimension(:) :: fields_vector,current_vector
        end type
        type(Mutual_matrix),allocatable,dimension(:,:) :: Mut_matrices
        integer :: i,j,itr,ii
        logical :: bouncing_flag
        real(8) :: max_error_value
!        complex*16:: temp_number

        do i=1,number_of_scatterers
            if(allocated(these(i)%Lm_collect)) then
                deallocate(these(i)%Nm_collect,these(i)%Lm_collect,these(i)%Dm_collect)
            endif
            if(allocated(these(i)%Y_collect)) then
                deallocate(these(i)%Y_collect)
            endif
            if(allocated(these(i)%Hd_c)) then
                deallocate(these(i)%Hd_c,these(i)%Ed_c)
            endif
            if(allocated(these(i)%Hd_c)) then
                deallocate(these(i)%Es_c,these(i)%Hs_c)
            endif
        enddo


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!! Deallocation and allocataion of discretization parameters
        do ii = 1,number_of_scatterers
            deallocate(these(ii)%tang_u,these(ii)%tang_v,these(ii)%norm_v,these(ii)%delta_n)
            deallocate(these(ii)%testing_pt,these(ii)%testing_pt_MoM)
            deallocate(these(ii)%testing_pt_status)

            if(allocated(these(ii)%contour_points)) then
            deallocate(these(ii)%contour_points,these(ii)%contour_points_type)
            endif

            if(allocated(these(ii)%contour_point_active_region)) then
            deallocate(these(ii)%contour_point_active_region)
            endif

            if(allocated(these(ii)%contour_points_orientation)) then
                deallocate(these(ii)%contour_points_orientation)
            endif

            if(allocated(these(ii)%Inside_sources_dielectric)) then
                deallocate(these(ii)%Inside_sources_dielectric)
            endif



            if(scatterers_input_method(ii) == 1) then
                call discretize_superquad(these(ii),MoM_Samples_per_wavelength,these)
            elseif(scatterers_input_method(ii) == 2) then
                call discretize_scatterer_file(these(ii),MoM_Samples_per_wavelength,these)
            endif

            if(allocated(these(ii)%eta_vv)) then
            deallocate(these(ii)%eta_vv,these(ii)%eta_uu,these(ii)%eta_uv,these(ii)%eta_vu)
            endif

            deallocate(these(ii)%norm_eta_u,these(ii)%norm_eta_v)
            if(these(ii)%Problem_Type == 3) then !! IBC
                call set_IBC_impedance_matrices(these(ii))
            else
                allocate(these(ii)%norm_eta_u(these(ii)%M),these(ii)%norm_eta_v(these(ii)%M))
                these(ii)%norm_eta_u = 1.d0
                these(ii)%norm_eta_v = 1.d0
            endif
        enddo
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        allocate(Mut_matrices(number_of_scatterers,number_of_scatterers))
        allocate(errors(number_of_scatterers),error_bouncing(number_of_scatterers))
        errors = 0.d0

        !! evaluate the mutual coupling matrices that can be multiplied to the vector [I_u,I_v,M_u,M_v]^t
        !! to produce [Eu_MoM, Ev_MoM, Hu_MoM, Hv_MoM]^t
        !! [ Ztt   Ztz     -Ytt        -Ytz    ]
        !! [ Zzt   Zzz     -Ytz        -Yzz    ]
        !! [ Ytt   Ytz   Ztt/eta^2   Ztz/eta^2 ]
        !! [ Yzt   Yzz   Zzt/eta^2   Zzz/eta^2 ]
        do i=1,number_of_scatterers !! testing
            do j=1,number_of_scatterers !! sources
                if(i == j) then
                    cycle
                endif
                !! allocating matrices
                allocate(Zzz(these(i)%M,these(j)%M,0:0),Zzt(these(i)%M,these(j)%M,0:0),Ztz(these(i)%M,these(j)%M,0:0),&
                Ztt(these(i)%M,these(j)%M,0:0),Yzt(these(i)%M,these(j)%M,0:0),&
                Ytz(these(i)%M,these(j)%M,0:0),Ytt(these(i)%M,these(j)%M,0:0))
                allocate(Z_l(2*these(i)%M,2*these(j)%M),Y_l(2*these(i)%M,2*these(j)%M))
                allocate(Mut_matrices(i,j)%ZY(4*these(i)%M,4*these(j)%M))
                allocate(Mut_matrices(i,j)%fields_vector(4*these(i)%M))
                allocate(Mut_matrices(i,j)%current_vector(4*these(j)%M))
                Mut_matrices(i,j)%m_rows = 4*these(i)%M
                Mut_matrices(i,j)%n_cols = 4*these(j)%M
                Mut_matrices(i,j)%Source_ID = these(j)%region_ID
                Mut_matrices(i,j)%Dist_ID = these(i)%region_ID
                call set_MoM_bouncing_matrices(these(j),these(i),Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az,these)


                Z_l(1:these(i)%M,1:these(j)%M) = Ztt(:,:,0)


                Z_l(1:these(i)%M,(these(j)%M+1):(2*these(j)%M)) = Ztz(:,:,0)
                Z_l((these(i)%M+1):(2*these(i)%M),1:these(j)%M) = Zzt(:,:,0)
                Z_l((these(i)%M+1):(2*these(i)%M),(these(j)%M+1):(2*these(j)%M)) = Zzz(:,:,0)

                Y_l(1:these(i)%M,1:these(j)%M) = Ytt(:,:,0)
                Y_l(1:these(i)%M,(these(j)%M+1):(2*these(j)%M)) = Ytz(:,:,0)
                Y_l((these(i)%M+1):(2*these(i)%M),1:these(j)%M) = Yzt(:,:,0)
                Y_l((these(i)%M+1):(2*these(i)%M),(these(j)%M+1):(2*these(j)%M)) = 0.d0

!                                write(*,*) 'determinants of Z and Y',det(Z_l),det(Y_l)
                Mut_matrices(i,j)%ZY(1:(2*these(i)%M),1:(2*these(j)%M)) = Z_l
                Mut_matrices(i,j)%ZY((2*these(i)%M+1):(4*these(i)%M),(2*these(j)%M+1):(4*these(j)%M)) = Z_l/eta1**2
                Mut_matrices(i,j)%ZY(1:(2*these(i)%M),(2*these(j)%M+1):(4*these(j)%M)) = -Y_l
                Mut_matrices(i,j)%ZY((2*these(i)%M+1):(4*these(i)%M),1:(2*these(j)%M)) = Y_l


                !! deallocating matrices
                deallocate(Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt)
                deallocate(Z_l,Y_l)
            enddo
        enddo

        do i=1,number_of_scatterers
            allocate(these(i)%combined_bouncing_fields(4*these(i)%M))
        enddo


        !! Performing Iterative Bouncing Fields
        bouncing_flag = .false.
        do itr = 0,max_bouncing_iteration
            if(bouncing_flag) then
                do i = 1,number_of_scatterers
                    these(i)%Eu_MoM = 0.d0
                    these(i)%Ev_MoM = 0.d0
                    these(i)%Hu_MoM = 0.d0
                    these(i)%Hv_MoM = 0.d0
                enddo
            endif
            do j = 1,number_of_scatterers
                if(.not. bouncing_flag) then
                    !! first iteration only (scattering from individual scatterers cosidering incident fields only)
                    allocate(these(j)%I_u(these(j)%M,0:N_order),these(j)%I_v(these(j)%M,0:N_order),&
                    these(j)%M_v(these(j)%M,0:N_order),these(j)%M_u(these(j)%M,0:N_order))
                    allocate(these(j)%I_u_bounce(these(j)%M,0:0),these(j)%I_v_bounce(these(j)%M,0:0),&
                    these(j)%M_u_bounce(these(j)%M,0:0),these(j)%M_v_bounce(these(j)%M,0:0))



                    call eval_current_MoM(these(j),N_order,these(j)%I_u,these(j)%I_v,these(j)%M_u,these(j)%M_v,&
                    bouncing_flag,these)

!                    write(*,*) 'here'

                    these(j)%I_u_bounce(:,0) = these(j)%I_u(:,0)
                    these(j)%I_v_bounce(:,0) = these(j)%I_v(:,0)
                    these(j)%M_u_bounce(:,0) = these(j)%M_u(:,0)
                    these(j)%M_v_bounce(:,0) = these(j)%M_v(:,0)

!                    write(*,*) these(j)%I_v(1,:)
                else

                    do i = 1,number_of_scatterers
                        if(i == j) then
                            cycle
                        !                        write(*,*) ' da 3abeet'
                        endif
                        Mut_matrices(i,j)%current_vector(1:(these(j)%M)) =  these(j)%I_u_bounce(:,0)
                        Mut_matrices(i,j)%current_vector((these(j)%M+1):(2*these(j)%M)) =  these(j)%I_v_bounce(:,0)
                        Mut_matrices(i,j)%current_vector((2*these(j)%M+1):(3*these(j)%M)) =  these(j)%M_u_bounce(:,0)
                        Mut_matrices(i,j)%current_vector((3*these(j)%M+1):(4*these(j)%M)) =  these(j)%M_v_bounce(:,0)
                        Mut_matrices(i,j)%fields_vector = matmul(Mut_matrices(i,j)%ZY,Mut_matrices(i,j)%current_vector)
!                                        write(*,*) 'bouncing between source',j,'and object',i,'has magnitude of',&
!                                             sum(abs(Mut_matrices(i,j)%fields_vector(1:100)))
                    enddo
                endif
            enddo
            if(number_of_scatterers == 1) then !! terminate the IBF loop if there is only one scatterer
                exit
            endif
            if(.not. bouncing_flag) then
                !            write(*,*) 'flipping at itr',itr
                bouncing_flag = .true.
                write(*,*) '=========================================='
                write(*,*) '==============MoM IFB starts=============='
                cycle
            endif

            !! now the new incident fields (due to bouncing) can be evaulated and thus the new bouncing currents
            do i = 1,number_of_scatterers
                these(i)%combined_bouncing_fields = 0.d0
                do j=1,number_of_scatterers
                    if(j == i) then
                        cycle
                    endif
                    these(i)%combined_bouncing_fields = these(i)%combined_bouncing_fields + Mut_matrices(i,j)%fields_vector
                enddo
                these(i)%Eu_MoM(:,0) = -these(i)%combined_bouncing_fields(1:(these(i)%M))
                these(i)%Ev_MoM(:,0) = -these(i)%combined_bouncing_fields((these(i)%M+1):(2*these(i)%M))
                these(i)%Hu_MoM(:,0) = -these(i)%combined_bouncing_fields((2*these(i)%M+1):(3*these(i)%M))
                these(i)%Hv_MoM(:,0) = -these(i)%combined_bouncing_fields((3*these(i)%M+1):(4*these(i)%M))
                !            write(*,*) 'fields due to bouncing are calculated',i
                call eval_current_MoM(these(i),N_order,these(i)%I_u_bounce,these(i)%I_v_bounce,these(i)%M_u_bounce,&
                these(i)%M_v_bounce,bouncing_flag,these)
                !            write(*,*) 'bouncing currents are calculated',i
                !! evaluating the error
                if(these(i)%Problem_Type == 1) then !! PEC
                    error_bouncing(i) = sum(abs((/these(i)%I_u_bounce,these(i)%I_v_bounce/))**2.d0)/&
                    sum(abs((/these(i)%I_u,these(i)%I_v/))**2.d0)
                elseif(these(i)%Problem_Type == 2) then !! PMC
                    error_bouncing(i) = sum(abs((/these(i)%M_u_bounce,these(i)%M_v_bounce/))**2.d0)/&
                    sum(abs((/these(i)%M_u,these(i)%M_v/))**2.d0)
                elseif(these(i)%Problem_Type == 3) then !! IBC
                    error_bouncing(i) = sum(abs((/these(i)%I_u_bounce,these(i)%I_v_bounce/))**2.d0)/&
                    sum(abs((/these(i)%I_u,these(i)%I_v/))**2.d0)
                else !! Dielectric
                    error_bouncing(i) = sum(abs((/these(i)%I_u_bounce,these(i)%I_v_bounce/))**2.d0)/&
                    sum(abs((/these(i)%I_u,these(i)%I_v/))**2.d0)
                endif
                !! update the total surface currents
                these(i)%I_u = these(i)%I_u + these(i)%I_u_bounce
                these(i)%I_v = these(i)%I_v + these(i)%I_v_bounce
                these(i)%M_u = these(i)%M_u + these(i)%M_u_bounce
                these(i)%M_v = these(i)%M_v + these(i)%M_v_bounce
            !            write(*,*) 'iteration #',itr,' Scatterer #',i,' Error Bouncing =',error_bouncing(i)
            enddo

            !! set exit conditions
            max_error_value = maxval(error_bouncing)
            write(*,*) 'iteration #',itr,' Max. Error Bouncing =',max_error_value
            if(max_error_value <= TOL) then
                write(*,*) 'MoM Accuracy reached after',itr,' IFB iterations'
                exit
            endif
        enddo
        if(itr >= max_bouncing_iteration) then
            write(*,*) 'WARNING: MoM IFB did not reach the desired accuracy tolerance;'
            write(*,*) '         Please, increase the max_bouncing_iteration limit'
        endif
        !! deallocating mutual matrices
        do i = 1,number_of_scatterers
            do j=1,number_of_scatterers
                if(i == j) then
                    cycle
                endif
                deallocate(Mut_matrices(i,j)%ZY)
                deallocate(Mut_matrices(i,j)%fields_vector,Mut_matrices(i,j)%current_vector)
            enddo
            deallocate(these(i)%combined_bouncing_fields)
            deallocate(these(i)%I_u_bounce,these(i)%I_v_bounce)
            deallocate(these(i)%M_u_bounce,these(i)%M_v_bounce)
        enddo
        deallocate(Mut_matrices)
        deallocate(error_bouncing)

        if(MoM_activation_flag == 1) then !! it is required to Compare farfields and currents between MoM and RAS
!            call compare_far_field_scatterers(these,180)
            call plot_far_field_MoM(these,180)
            !! plot currents
            do i = 1,number_of_scatterers
                errors(i) = error_plot_current_comparison_with_MoM(these(i))
            enddo

        endif

    end function eval_surface_current_error_multiscatterer_MoM

    function error_plot_current_comparison_with_MoM(this) result(error)
        type(Scatterer) :: this
        complex*16,allocatable,dimension(:) :: Jv,Ju,H_total,E_total,Mu,Mv
        integer ::  fd1 = 50
        real(8),allocatable,dimension(:) :: phase_u,phase_v,Jv_phase,Ju_phase,&
        M_phase_u,M_phase_v,Mv_phase,Mu_phase
        real(8) :: t,error,normalize
        integer :: i
        character*30 :: file_name


        if(plot_current /= 0) then
            file_name = 'currents_plotting_RID_'//num2str(this%region_ID,1)//'.dat'
            OPEN(fd1, FILE=trim(file_name))
        endif
        !    if(.not. allocated(this%I_u)) then
        !        allocate(this%I_u(this%M,0:N_order),this%I_v(this%M,0:N_order),&
        !        this%M_u(this%M,0:N_order),this%M_v(this%M,0:N_order))
        !        this%I_u = 0.d0
        !        this%I_v = 0.d0
        !        this%M_u = 0.d0
        !        this%M_v = 0.d0
        !    endif

        allocate(phase_u(this%M),phase_v(this%M),Jv(this%M),Ju(this%M),M_phase_u(this%M),&
        M_phase_v(this%M),Mv(this%M),Mu(this%M))
        allocate(Jv_phase(this%M),Ju_phase(this%M),Mv_phase(this%M),Mu_phase(this%M))
        !! then check the availability of RSM currents for this scatterer

        phase_u = atan2(aimag(this%I_u(:,0)),real(this%I_u(:,0)))*180.d0/pi
        phase_v = atan2(aimag(this%I_v(:,0)),real(this%I_v(:,0)))*180.d0/pi
        M_phase_u = atan2(aimag(this%M_u(:,0)),real(this%M_u(:,0)))*180.d0/pi
        M_phase_v = atan2(aimag(this%M_v(:,0)),real(this%M_v(:,0)))*180.d0/pi

        allocate(H_total(2*this%M),E_total(2*this%M))

!        write(*,*) this%Hs
!        write(*,*) this%H
!        write(*,*)

        E_total = this%Es(:,0) + this%E(:,0) + this%E_bounce(:,0)
        H_total = this%Hs(:,0) + this%H(:,0) + this%H_bounce(:,0)

        Mu = E_total((this%M+1):(2*this%M))
        Mv = -E_total(1:this%M)

        Ju = -H_total((this%M+1):(2*this%M))
        Jv = H_total(1:this%M)
        Jv_phase = 180.d0/pi*atan2(aimag(Jv),real(Jv))
        Ju_phase = 180.d0/pi*atan2(aimag(Ju),real(Ju))
        Mv_phase = 180.d0/pi*atan2(aimag(Mv),real(Mv))
        Mu_phase = 180.d0/pi*atan2(aimag(Mu),real(Mu))
        deallocate(H_total,E_total)

        t = -this%delta_n(1)/2.d0
        error = 0.d0
        normalize = 0.d0
        do i=1,this%M
            t = t+this%delta_n(i)

            normalize = normalize + abs(this%I_u(i,0))**2.d0 + abs(this%I_v(i,0))**2.d0 + &
            abs(this%M_u(i,0))**2.d0 + abs(this%M_v(i,0))**2.d0
            error = error + abs(Ju(i)-this%I_u(i,0))**2.d0 + abs(Jv(i)-this%I_v(i,0))**2.d0 + &
            abs(Mu(i)-this%M_u(i,0))**2.d0 + abs(Mv(i)-this%M_v(i,0))**2.d0

            if(plot_current /= 0) then
                write(fd1,*) t/lambda,abs(Jv(i)),Jv_phase(i),abs(this%I_v(i,0)),phase_v(i),&
                abs(Ju(i)),Ju_phase(i),abs(this%I_u(i,0)),phase_u(i),&
                abs(Mv(i)),Mv_phase(i),abs(this%M_v(i,0)),M_phase_v(i),&
                abs(Mu(i)),Mu_phase(i),abs(this%M_u(i,0)),M_phase_u(i)
            endif
        enddo
        error = error/normalize

        if(plot_current /= 0) then
            close(fd1)
        endif
        deallocate(phase_u,phase_v,Jv_phase,Ju_phase,Mv_phase,Mu_phase)
    end function error_plot_current_comparison_with_MoM

    subroutine plot_current_RAS(this)
        type(Scatterer) :: this
        complex*16,allocatable,dimension(:) :: Jv,Ju,H_total,E_total,Mu,Mv
        integer ::  fd1 = 50
        real(8),allocatable,dimension(:) :: Jv_phase,Ju_phase,Mv_phase,Mu_phase
        real(8) :: t
        integer :: i
        character*30 :: file_name


        if(plot_current /= 0) then
            file_name = 'currents_plotting_RID_'//&
                num2str(this%region_ID,1)//'_RAS.dat'
            OPEN(fd1, FILE=trim(file_name))
        endif


        allocate(Jv(this%M),Ju(this%M),Mv(this%M),Mu(this%M))
        allocate(Jv_phase(this%M),Ju_phase(this%M),Mv_phase(this%M),Mu_phase(this%M))
        allocate(H_total(2*this%M),E_total(2*this%M))

        E_total = this%Es(:,0) + this%E(:,0) + this%E_bounce(:,0)
        H_total = this%Hs(:,0) + this%H(:,0) + this%H_bounce(:,0)

        Mu = E_total((this%M+1):(2*this%M))
        Mv = -E_total(1:this%M)

        Ju = -H_total((this%M+1):(2*this%M))
        Jv = H_total(1:this%M)
        Jv_phase = 180.d0/pi*atan2(aimag(Jv),real(Jv))
        Ju_phase = 180.d0/pi*atan2(aimag(Ju),real(Ju))
        Mv_phase = 180.d0/pi*atan2(aimag(Mv),real(Mv))
        Mu_phase = 180.d0/pi*atan2(aimag(Mu),real(Mu))
        deallocate(H_total,E_total)

        t = -this%delta_n(1)/2.d0
        do i=1,this%M
            t = t+this%delta_n(i)

            if(plot_current /= 0) then
                write(fd1,*) t/lambda,abs(Jv(i)),Jv_phase(i),&
                abs(Ju(i)),Ju_phase(i),&
                abs(Mv(i)),Mv_phase(i),&
                abs(Mu(i)),Mu_phase(i)
            endif
        enddo

        if(plot_current /= 0) then
            close(fd1)
        endif
        deallocate(Jv_phase,Ju_phase,Mv_phase,Mu_phase)
    end subroutine plot_current_RAS


    subroutine plot_current_MoM(this)
        type(Scatterer) :: this
        integer ::  fd1 = 50
        real(8),allocatable,dimension(:) :: phase_u,phase_v,M_phase_u,M_phase_v
        real(8) :: t
        integer :: i
        character*30 :: file_name




        if(plot_current /= 0) then
            file_name = 'currents_plotting_RID_'//&
                num2str(this%region_ID,1)//'_MoM.dat'
            OPEN(fd1, FILE=trim(file_name))
        endif



        allocate(phase_u(this%M),phase_v(this%M),M_phase_u(this%M),&
        M_phase_v(this%M))
        !! then check the availability of RSM currents for this scatterer

        phase_u = atan2(aimag(this%I_u(:,0)),real(this%I_u(:,0)))*180.d0/pi
        phase_v = atan2(aimag(this%I_v(:,0)),real(this%I_v(:,0)))*180.d0/pi
        M_phase_u = atan2(aimag(this%M_u(:,0)),real(this%M_u(:,0)))*180.d0/pi
        M_phase_v = atan2(aimag(this%M_v(:,0)),real(this%M_v(:,0)))*180.d0/pi

!        write(*,*) 'here'

        t = -this%delta_n(1)/2.d0
        do i=1,this%M
            t = t+this%delta_n(i)

            if(plot_current /= 0) then
                write(fd1,*) t/lambda,abs(this%I_v(i,0)),phase_v(i),&
                abs(this%I_u(i,0)),phase_u(i),&
                abs(this%M_v(i,0)),M_phase_v(i),&
                abs(this%M_u(i,0)),M_phase_u(i)
            endif
        enddo

        if(plot_current /= 0) then
            close(fd1)
        endif
        deallocate(phase_u,phase_v,M_phase_u,M_phase_v)
    end subroutine plot_current_MoM



    function eval_surface_current_error_MoM(this,N_order,Scatterers_Lib) result(error)
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        type(Scatterer) :: this
        integer :: N_order
        complex*16,allocatable,dimension(:) :: Jv,Ju,H_total,E_total,Mu,Mv
        integer ::  fd1 = 50
        real(8),allocatable,dimension(:) :: phase_u,phase_v,Jv_phase,Ju_phase,&
        M_phase_u,M_phase_v,Mv_phase,Mu_phase
        real(8) :: t,error,normalize
        integer :: i

        if(plot_current /= 0) then
            OPEN(fd1, FILE='currents_plotting.dat')
        endif
        allocate(this%I_u(this%M,0:N_order),this%I_v(this%M,0:N_order),&
        this%M_u(this%M,0:N_order),this%M_v(this%M,0:N_order))


        allocate(phase_u(this%M),phase_v(this%M),Jv(this%M),Ju(this%M),M_phase_u(this%M),&
        M_phase_v(this%M),Mv(this%M),Mu(this%M))
        allocate(Jv_phase(this%M),Ju_phase(this%M),Mv_phase(this%M),Mu_phase(this%M))



        call eval_current_MoM(this,N_order,this%I_u,this%I_v,this%M_u,this%M_v,.false.,Scatterers_Lib)


        phase_u = atan2(aimag(this%I_u(:,0)),real(this%I_u(:,0)))*180.d0/pi
        phase_v = atan2(aimag(this%I_v(:,0)),real(this%I_v(:,0)))*180.d0/pi
        M_phase_u = atan2(aimag(this%M_u(:,0)),real(this%M_u(:,0)))*180.d0/pi
        M_phase_v = atan2(aimag(this%M_v(:,0)),real(this%M_v(:,0)))*180.d0/pi


        allocate(H_total(2*this%M),E_total(2*this%M))

        E_total = this%Es(:,0) + this%E(:,0)
        H_total = this%Hs(:,0) + this%H(:,0)

        Mu = E_total((this%M+1):(2*this%M))
        Mv = -E_total(1:this%M)

        Ju = -H_total((this%M+1):(2*this%M))
        Jv = H_total(1:this%M)
        Jv_phase = 180.d0/pi*atan2(aimag(Jv),real(Jv))
        Ju_phase = 180.d0/pi*atan2(aimag(Ju),real(Ju))
        Mv_phase = 180.d0/pi*atan2(aimag(Mv),real(Mv))
        Mu_phase = 180.d0/pi*atan2(aimag(Mu),real(Mu))
        deallocate(H_total,E_total)

        t = -this%delta_n(1)/2.d0
        error = 0.d0
        normalize = 0.d0
        do i=1,this%M
            t = t+this%delta_n(i)

            normalize = normalize + abs(this%I_u(i,0))**2.d0 + abs(this%I_v(i,0))**2.d0 + &
            abs(this%M_u(i,0))**2.d0 + abs(this%M_v(i,0))**2.d0
            error = error + abs(Ju(i)-this%I_u(i,0))**2.d0 + abs(Jv(i)-this%I_v(i,0))**2.d0 + &
            abs(Mu(i)-this%M_u(i,0))**2.d0 + abs(Mv(i)-this%M_v(i,0))**2.d0

            if(plot_current /= 0) then
                write(fd1,*) t,abs(Jv(i)),Jv_phase(i),abs(this%I_v(i,0)),phase_v(i),&
                abs(Ju(i)),Ju_phase(i),abs(this%I_u(i,0)),phase_u(i),&
                abs(Mv(i)),Mv_phase(i),abs(this%M_v(i,0)),M_phase_v(i),&
                abs(Mu(i)),Mu_phase(i),abs(this%M_u(i,0)),M_phase_u(i)
            endif
        enddo
        error = error/normalize

        if(plot_current /= 0) then
            close(fd1)
        endif
        deallocate(phase_u,phase_v,Jv_phase,Ju_phase,Mv_phase,Mu_phase)

    end function eval_surface_current_error_MoM

    function eval_surface_current_Dielectric(this) result(error)
        type(Scatterer) :: this
        real(8) :: error,normalize
        complex*16,dimension(this%M) :: Jz,Jt,Mz,Mt
        integer ::  fd1 = 50
        complex*16,dimension(2*this%M) :: H_total,E_total
        integer :: j,kk
        real(8),dimension(this%M) :: Jz_phase,Jt_phase,Mz_phase,Mt_phase
        real(8) :: t,kz
        complex*16 :: kr1,kr2
        complex*16,allocatable,dimension(:)::h2_1,h2_2,h2p_1,h2p_2,besj_1,besjp_1,besj_2,besjp_2,Denum,Denum_TE
        integer :: n_terms
        complex*16 :: Ja_z_TM,Ja_t_TE,Ja_z_TE,Ja_z,Ja_t !! analytical surface current
        complex*16 :: Ma_z_TM,Ma_t_TM,Ma_z,Ma_t !! analytical surface magnetic current
        complex*16 :: sum_Jz_TM,sum_Jz_TE,sum_Jt_TE,cjj,const_Ja_z_TM,const_Ja_z_TE,const_Ja_t_TE,const_Ma_z_TE,sum_Mz_TE
        complex*16 :: const_Ma_z_TM,const_Ma_t_TM,sum_Mz_TM,sum_Mt_TM,Ma_z_TE
        real(8) :: x_s,y_s,phi,Ja_z_phase,Ja_t_phase,Ma_z_phase,Ma_t_phase,R_cylind
        character*30 :: file_name

        if(plot_current /= 0) then
            file_name = 'currents_plotting'//num2str(this%region_ID,1)//'.dat'
            OPEN(fd1, FILE=trim(file_name))
        endif

        !------------------------------------------------------
        !! Initialization of the exact current distribution

        R_cylind = this%a_superquad
        kz = k0*cos(theta_i)
        kr1 = sqrt(k1**2.d0 - kz**2.d0)
        kr2 = sqrt(this%k_local - kz)*sqrt(this%k_local + kz)

        const_Ja_z_TM = (-2.d0,0.d0)*eta1*k1*this%k_local/(pi*R_cylind*kr2*kr1**2.d0)*sin(theta_i)*cos(alpha_i)
        const_Ma_z_TM =(4.d0,0.d0)*eta1*this%eta_local/(pi*R_cylind**2.d0)*kz/(kr1**3.d0)*cos(alpha_i) !! oblique
        const_Ma_t_TM = (0.d0,-2.d0)*eta1*this%eta_local/(pi*R_cylind)*k1/(kr1**2.d0)*sin(theta_i)*cos(alpha_i)

        const_Ja_z_TE = (-4.d0,0.d0)*sin(alpha_i)*sin(theta_i)*eta1*k1*kz/(pi*kr2*kr2*kr1**2.d0*R_cylind**2.d0) !! oblique
        const_Ja_t_TE = (0.d0,2.d0)*eta1/(pi*R_cylind)*k1/kr1**2.d0*sin(theta_i)*sin(alpha_i)
        const_Ma_z_TE = (-2.d0,0.d0)*eta1*this%eta_local/(pi*R_cylind)*k1*this%k_local/(kr1**2.d0*kr2*kr2)*sin(theta_i)*sin(alpha_i)

        n_terms = int(2*k1*R_cylind)+1
        allocate(this%I_u(this%M,0:0),this%I_v(this%M,0:0),this%M_u(this%M,0:0),this%M_v(this%M,0:0))
        allocate(h2_1(0:n_terms),h2p_1(0:n_terms),besj_1(0:n_terms),besjp_1(0:n_terms),Denum(0:n_terms))
        allocate(h2_2(0:n_terms),h2p_2(0:n_terms),besj_2(0:n_terms),besjp_2(0:n_terms),Denum_TE(0:n_terms))

        call hankarray(n_terms,kr1*R_cylind,h2_1,h2p_1)
        besj_1 = real(h2_1)
        besjp_1 = real(h2p_1)
        call hankarray(n_terms,kr2*R_cylind,h2_2,h2p_2)
        call besselarray(n_terms,kr2*R_cylind,besj_2,besjp_2)
        !    besj_2 = real(h2_2)
        !    besjp_2 = real(h2p_2)
        Denum = this%eta_local*k1/kr1*h2p_1*besj_2-eta1*this%k_local/kr2*h2_1*besjp_2
        Denum_TE = eta1*k1/kr1*h2p_1*besj_2-this%eta_local*this%k_local/kr2*h2_1*besjp_2
        !-------------------------------------------------------
        !! computed current
        H_total = this%Hs(:,0) + this%H(:,0)
        E_total = this%Es(:,0) + this%E(:,0)
        Jt = -H_total((this%M+1):(2*this%M))
        Jz = H_total(1:this%M)
        Mt = E_total((this%M+1):(2*this%M))
        Mz = -E_total(1:this%M)
        Jz_phase = 180.d0/pi*atan2(aimag(Jz),real(Jz))
        Jt_phase = 180.d0/pi*atan2(aimag(Jt),real(Jt))
        Mz_phase = 180.d0/pi*atan2(aimag(Mz),real(Mz))
        Mt_phase = 180.d0/pi*atan2(aimag(Mt),real(Mt))
        t=-this%delta_n(1)
        error = 0.d0
        normalize = 0.d0
        do j =1,this%M

            x_s = this%testing_pt(j,1)
            y_s = this%testing_pt(j,2)
            phi = atan2(y_s,x_s)

            sum_Jz_TM = 0.d0
            sum_Jz_TE = 0.d0
            sum_Jt_TE = 0.d0
            sum_Mz_TM = 0.d0
            sum_Mt_TM = 0.d0
            sum_Mz_TE = 0.d0

            cjj = cj
            do kk=1,n_terms
                sum_Jz_TM = sum_Jz_TM + cjj*besjp_2(kk)/Denum(kk)*cos(kk*(phi-phi_i))
                sum_Mt_TM = sum_Mt_TM + cjj*besj_2(kk)/Denum(kk)*cos(kk*(phi-phi_i))
                sum_Mz_TM = sum_Mz_TM + kk*cjj*besj_2(kk)/Denum(kk)*sin(kk*(phi-phi_i))

                sum_Jz_TE = sum_Jz_TE + kk*cjj*besj_2(kk)/Denum_TE(kk)*sin(kk*(phi-phi_i))
                sum_Jt_TE = sum_Jt_TE + cjj*besj_2(kk)/Denum_TE(kk)*cos(kk*(phi-phi_i))
                sum_Mz_TE = sum_Mz_TE + cjj*besjp_2(kk)/Denum_TE(kk)*cos(kk*(phi-phi_i))
                cjj = cjj*cj
            enddo
            Ja_z_TM = const_Ja_z_TM*(besjp_2(0)/Denum(0) + (2.d0,0.d0)*sum_Jz_TM )
            Ma_t_TM = const_Ma_t_TM*(besj_2(0)/Denum(0) + (2.d0,0.d0)*sum_Mt_TM )
            Ma_z_TM = const_Ma_z_TM*sum_Mz_TM

            Ja_t_TE = const_Ja_t_TE*(besj_2(0)/Denum_TE(0) +(2.d0,0.d0)*sum_Jt_TE)
            Ja_z_TE = const_Ja_z_TE*sum_Jz_TE
            Ma_z_TE = kr2*const_Ma_z_TE*( besjp_2(0)/Denum_TE(0) + (2.d0,0.d0)*sum_Mz_TE)

            Ja_z = Ja_z_TM + Ja_z_TE
            Ja_t = Ja_t_TE

            Ma_z = Ma_z_TM + Ma_z_TE
            Ma_t = Ma_t_TM

            Ja_z_phase = 180.d0/pi*atan2(aimag(Ja_z),real(Ja_z))
            Ja_t_phase = 180.d0/pi*atan2(aimag(Ja_t),real(Ja_t))
            Ma_z_phase = 180.d0/pi*atan2(aimag(Ma_z),real(Ma_z))
            Ma_t_phase = 180.d0/pi*atan2(aimag(Ma_t),real(Ma_t))

            t = t + this%delta_n(1)

            normalize = normalize + abs(Ja_z)**2.d0 + abs(Ja_t)**2.d0 + abs(Ma_z)**2.d0 + abs(Ma_t)**2.d0
            error = error + abs(Jt(j)-Ja_t)**2.d0 + abs(Jz(j)-Ja_z)**2.d0 + &
            abs(Mt(j)-Ma_t)**2.d0 + abs(Mz(j)-Ma_z)**2.d0

            this%I_u(j,0) = Ja_t
            this%I_v(j,0) = Ja_z
            this%M_u(j,0) = Ma_t
            this%M_v(j,0) = Ma_z

            if(plot_current /= 0) then
                write(fd1,*) t/lambda,abs(Jz(j)),Jz_phase(j),abs(Ja_z),Ja_z_phase,&
                abs(Jt(j)),Jt_phase(j),abs(Ja_t),Ja_t_phase,&
                abs(Mz(j)),Mz_phase(j),abs(Ma_z),Ma_z_phase,&
                abs(Mt(j)),Mt_phase(j),abs(Ma_t),Ma_t_phase
            endif

        enddo
        close(fd1)
        deallocate(h2_1,h2p_1,besj_1,besjp_1,Denum)
        deallocate(h2_2,h2p_2,besj_2,besjp_2)
    end function eval_surface_current_Dielectric



    subroutine eval_current_MoM(this,N_order,I_u,I_v,M_u,M_v,bouncing_flag,Scatterers_Lib)

        logical :: bouncing_flag
        !! false for the first time calculation for bouncing fields or a single scatterer problem
        !! true if bouncing field is computed
        !! -> here no need to recompute the matrices or allocating them. Just used the pre-stored matrix inversions and use it to computed
        !!    the currents.
        !! -> for computing the bouncing currents Eu_MoM, Ev_MoM ... shoud be computed before statrting this function
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        integer :: N_order,nn,ii
        complex*16,allocatable,intent(inout) :: I_u(:,:),I_v(:,:),M_u(:,:),M_v(:,:)
        complex*16,allocatable,dimension(:,:) :: I_j,E,H,EH
        complex*16,allocatable,dimension(:) :: b_n
        complex*16,allocatable,dimension(:,:,:) :: Z,Y_f,Y,Z_f,Z_d,Y_f_d,Y_d,Z_f_d
        complex*16,allocatable,dimension(:,:,:) :: Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt
        complex*16,allocatable,dimension(:,:,:) :: Zzz_d,Zzt_d,Ztz_d,Ztt_d,Yzt_d,Ytz_d,Ytt_d
        complex*16,allocatable,dimension(:,:) ::eta_s,eta_f

        allocate(E(2*this%M,0:N_order),H(2*this%M,0:N_order),EH(2*this%M,0:N_order))

        if(.not. bouncing_flag) then
            allocate(Z(2*this%M,2*this%M,0:N_order),Y_f(2*this%M,2*this%M,0:N_order))
            allocate(this%Eu_MoM(this%M,0:N_order),this%Ev_MoM(this%M,0:N_order),&
            this%Hu_MoM(this%M,0:N_order),this%Hv_MoM(this%M,0:N_order))
            allocate(Zzz(this%M,this%M,0:N_order),Zzt(this%M,this%M,0:N_order),&
            Ztz(this%M,this%M,0:N_order),Ztt(this%M,this%M,0:N_order),&
            Yzt(this%M,this%M,0:N_order),Ytz(this%M,this%M,0:N_order),Ytt(this%M,this%M,0:N_order))



            allocate(Y(2*this%M,2*this%M,0:N_order),Z_f(2*this%M,2*this%M,0:N_order))
            if(this%Problem_Type == 4) then
                allocate(Y_d(2*this%M,2*this%M,0:N_order),Z_f_d(2*this%M,2*this%M,0:N_order))
                allocate(Y_f_d(2*this%M,2*this%M,0:N_order),Z_d(2*this%M,2*this%M,0:N_order))
                allocate(Zzz_d(this%M,this%M,0:N_order),Zzt_d(this%M,this%M,0:0),&
                Ztz_d(this%M,this%M,0:N_order),Ztt_d(this%M,this%M,0:N_order),&
                Yzt_d(this%M,this%M,0:N_order),Ytz_d(this%M,this%M,0:N_order),Ytt_d(this%M,this%M,0:N_order))
            endif

            if(Wideband_Type > 0) then
                if(this%Problem_Type == 4) then
                    allocate(this%MoM_ZY(4*this%M,4*this%M,0:N_order))
                else
                    allocate(this%MoM_ZY(2*this%M,2*this%M,0:N_order))
                endif
            else
                if(this%Problem_Type == 4) then
                    allocate(this%MoM_ZY(4*this%M,4*this%M,0:0))
                else
                    allocate(this%MoM_ZY(2*this%M,2*this%M,0:0))
                endif
            endif





            !        allocate(testing_pt_MoM(this%M))

            !        do i=1,this%M
            !            testing_pt_MoM(i) = vec(this%testing_pt(i,1),this%testing_pt(i,2),0.d0)
            !        enddo

            !        write(*,*) testing_pt_MoM(10),this%testing_pt(10,:)
            call set_excitation_v(this,this%testing_pt_MoM,this%Ei,this%Hi)

!            write(*,*) this%Ei

            call set_excitation_vectors(this,this%Eu_MoM,this%Ev_MoM,this%Hu_MoM,this%Hv_MoM)

!            write(*,*) 'here',this%M
!            stop


            call set_MoM_matrices_separate(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az,.true.,Scatterers_Lib)

!            write(*,*) 'here',this%M

            if(this%Problem_Type == 4) then !! Dielectric Boundary
                !            write(*,*) 'Dielectric matrix filling'
                call set_MoM_matrices_separate(this,Zzz_d,Zzt_d,Ztz_d,Ztt_d,Yzt_d,Ytz_d,Ytt_d,az,&
                .false.,Scatterers_Lib)

            endif

!            write(*,*) 'here',this%M

            if(Wideband_Type > 0) then
                call set_MoM_AWE_matrices(this,N_order,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ar,az,eta1,.true.)
                if(this%Problem_Type == 4) then !! Dielectric
                    call set_MoM_AWE_matrices(this,N_order,Zzz_d,Zzt_d,Ztz_d,Ztt_d,Yzt_d,Ytz_d,Ytt_d,&
                    this%ar_local,az,this%eta_local,.false.)
                endif
            endif

!        write(*,*) 'Zzz'
!        do nn =0,N_order
!            write(*,*) det(Zzz(:,:,0)),det(Ztt(:,:,0)),det(Ztz(:,:,0))
!        enddo
!        stop
            !
            !        write(*,*) 'Ztt'
            !        do nn =0,N_order
            !            write(*,*) Ztt(1,1:3,nn)
            !        enddo





            Y(1:this%M,1:this%M,:) = 0.d0
            Y(1:this%M,(this%M+1):2*this%M,:) = Yzt
            Y((this%M+1):2*this%M,1:this%M,:) = Ytz
            Y((this%M+1):2*this%M,(this%M+1):2*this%M,:) = Ytt

            Y_f(1:this%M,1:this%M,:) = Ytz
            Y_f(1:this%M,(this%M+1):2*this%M,:) = Ytt
            Y_f((this%M+1):2*this%M,1:this%M,:) = 0.d0
            Y_f((this%M+1):2*this%M,(this%M+1):2*this%M,:) = -Yzt

            deallocate(Yzt,Ytz,Ytt)
            Z(1:this%M,1:this%M,:) = Zzz
            Z(1:this%M,(this%M+1):2*this%M,:) = Zzt
            Z((this%M+1):2*this%M,1:this%M,:) = Ztz
            Z((this%M+1):2*this%M,(this%M+1):2*this%M,:) = Ztt


!             write(*,*) 'Z'
!            do nn =0,N_order
!                write(*,*) det(Z(:,:,0))
!            enddo
!            stop
!            call write_matrix(Z(:,:,0),'Z.dat')
!            stop

            Z_f(1:this%M,1:this%M,:) = Ztz
            Z_f(1:this%M,(this%M+1):2*this%M,:) = Ztt
            Z_f((this%M+1):2*this%M,1:this%M,:) = -Zzz
            Z_f((this%M+1):2*this%M,(this%M+1):2*this%M,:) = -Zzt

            deallocate(Zzz,Zzt,Ztz,Ztt)
            if(this%Problem_Type == 4) then
                Z_d(1:this%M,1:this%M,:) = Zzz_d
                Z_d(1:this%M,(this%M+1):2*this%M,:) =Zzt_d
                Z_d((this%M+1):2*this%M,1:this%M,:) =Ztz_d
                Z_d((this%M+1):2*this%M,(this%M+1):2*this%M,:) =Ztt_d

                Z_f_d(1:this%M,1:this%M,:) =  Ztz_d
                Z_f_d(1:this%M,(this%M+1):2*this%M,:) = Ztt_d
                Z_f_d((this%M+1):2*this%M,1:this%M,:) = -Zzz_d
                Z_f_d((this%M+1):2*this%M,(this%M+1):2*this%M,:) =-Zzt_d

                deallocate(Zzz_d,Zzt_d,Ztz_d,Ztt_d)
                Y_d(1:this%M,1:this%M,:) = 0.d0
                Y_d(1:this%M,(this%M+1):2*this%M,:) = Yzt_d
                Y_d((this%M+1):2*this%M,1:this%M,:) = Ytz_d
                Y_d((this%M+1):2*this%M,(this%M+1):2*this%M,:) = Ytt_d

                Y_f_d(1:this%M,1:this%M,:) = Ytz_d
                Y_f_d(1:this%M,(this%M+1):2*this%M,:) = Ytt_d
                Y_f_d((this%M+1):2*this%M,1:this%M,:) = 0.d0
                Y_f_d((this%M+1):2*this%M,(this%M+1):2*this%M,:) = -Yzt_d

                deallocate(Yzt_d,Ytz_d,Ytt_d)
            endif





        endif
        !! MoM Problem formulation
        do nn =0,N_order
            E(:,nn) = (/this%Ev_MoM(:,nn),this%Eu_MoM(:,nn)/)
            H(:,nn) = (/this%Hu_MoM(:,nn),-this%Hv_MoM(:,nn)/)
        enddo

        !! MoM Problem formulation
        if(this%Problem_Type == 3) then !! IBC
            call eval_IBC_impedance_matrices(this,eta_s,eta_f)
            if(FormulationType == 2) then !! IBCH

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: IBC, Formulation Type: IBCH'
                    do nn = 0,N_order
                        this%MoM_ZY(:,:,nn) = Y_f(:,:,nn) - matmul(Z_f(:,:,nn),eta_f)/eta1
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = H(:,:) ! E!-H
            elseif(FormulationType == 3) then !! IBCC

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: IBC, Formulation Type: IBCC'
                    do nn = 0,N_order
                        this%MoM_ZY(:,:,nn) = Z(:,:,nn)+eta1*matmul(Y(:,:,nn),eta_f) + eta1*Y_f(:,:,nn) - matmul(Z_f(:,:,nn),eta_f)
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = E(:,:) + eta1*H(:,:)
            elseif(FormulationType == 4) then !! IBC

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: IBC, Formulation Type: IBC'
                    do nn = 0,N_order
                        this%MoM_ZY(:,:,nn) = Z(:,:,nn)+eta1*matmul(Y(:,:,nn),eta_f) - &
                        eta1*matmul(eta_s,(Y_f(:,:,nn) - matmul(Z_f(:,:,nn),eta_f)/eta1))
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                do nn = 0,N_order
                    EH(:,nn) = E(:,nn) - eta1*matmul(eta_s,H(:,nn))
                enddo
            else !! IBCE

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: IBC, Formulation Type: IBCE'
                    do nn = 0,N_order
                        this%MoM_ZY(:,:,nn) = Z(:,:,nn) + eta1*matmul(Y(:,:,nn),eta_f)
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = E(:,:)! E!-H
            endif
            deallocate(eta_s,eta_f)
        elseif(this%Problem_Type == 1) then !! PEC
            if(FormulationType == 1) then !! EFIE

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: PEC, EFIE'
                    do nn = 0,N_order
                        this%MoM_ZY(:,:,nn) = Z(:,:,nn)
                    !                write(*,*) this%MoM_ZY(1,1:3,nn)
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = E(:,:)
            elseif(FormulationType == 2) then !! MFIE

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: PEC, MFIE'
                    do nn = 0,N_order
                        this%MoM_ZY(:,:,nn) = Y_f(:,:,nn)
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = H(:,:)
            else !! CFIE

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: PEC, CFIE'
                    do nn = 0,N_order
                        this%MoM_ZY(:,:,nn) = Z(:,:,nn) + eta1*Y_f(:,:,nn)
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = E(:,:) + eta1*H(:,:)
            endif

        elseif(this%Problem_Type == 2) then !! PMC
            if(FormulationType == 1) then !! EFIE

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: PMC, EFIE'
                    do nn =0,N_order
                        this%MoM_ZY(:,:,nn) = Y(:,:,nn)
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = -E(:,:)
            elseif(FormulationType == 2) then !! MFIE

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: PMC, MFIE'
                    do nn=0,N_order
                        this%MoM_ZY(:,:,nn) = Z_f(:,:,nn)/eta1**2.d0
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = H(:,:)
            else !! CFIE

                if(.not. bouncing_flag) then
                    write(*,*) 'Scatterer #',this%region_ID,'  B.C.: PMC, CFIE'
                    do nn=0,N_order
                        this%MoM_ZY(:,:,nn) = Z_f(:,:,nn)/eta1**2.d0 + Y(:,:,nn)/eta1
                    enddo
                    deallocate(Z,Y_f,Y,Z_f)

                endif
                EH(:,:) = H(:,:) - E(:,:)/eta1
            endif

        else !! Dielectric
            !            call write_matrix(Z_f,'Y_mat.dat')
            !            call write_matrix(Z_f_d,'Y_d_mat.dat')
            if(.not. bouncing_flag) then
                write(*,*) 'Scatterer #',this%region_ID,'  B.C.: Dielectric'
            endif
        endif


        if(this%Problem_Type == 4) then !! Only Dielectric
            deallocate(EH)
            allocate(EH(4*this%M,0:N_order))
            if(.not. bouncing_flag) then
                this%MoM_ZY(1:(2*this%M),1:(2*this%M),:) = (Z + Z_d)/eta1
                this%MoM_ZY(1:(2*this%M),(2*this%M+1):(4*this%M),:) = -(Y + Y_d )
                this%MoM_ZY((2*this%M+1):(4*this%M),1:(2*this%M),:) = (Y_f + Y_f_d)
                this%MoM_ZY((2*this%M+1):(4*this%M),(2*this%M+1):(4*this%M),:) = eta1*(Z_f/eta1**2.d0 + &
                Z_f_d/this%eta_local**2.d0)
                deallocate(Z,Y_f,Y,Z_f)
                deallocate(Z_d,Y_d,Z_f_d,Y_f_d)
                allocate(this%MoM_matrix_inverted(4*this%M,4*this%M))
                this%MoM_matrix_inverted = Mat_Inv(this%MoM_ZY(:,:,0))
            endif
            EH(1:(2*this%M),:) = E(:,:)/eta1
            EH((2*this%M+1):(4*this%M),:) = H(:,:)
            allocate(I_j(4*this%M,0:N_order))
            !!I_j = Gauss_Elimination(Z_big,Ext_big,4*this%M)
            I_j(:,0) = matmul(this%MoM_matrix_inverted,EH(:,0))

            if(Wideband_Type > 0) then
                allocate(b_n(4*this%M))
                do nn = 1,N_order
                    b_n = EH(:,nn)
                    do ii = 1,nn
                        b_n = b_n - dble(n_X_i_matrix(nn,ii))*matmul(This%MoM_ZY(:,:,ii),I_j(:,nn-ii))
                    enddo
                    I_j(:,nn) = matmul(this%MoM_matrix_inverted,b_n)/(dble(n_X_i_matrix(nn,1)))
                enddo
                deallocate(b_n)
            endif

            I_v = I_j(1:this%M,:)
            I_u = I_j((this%M+1):2*this%M,:)
            M_v = eta1*I_j((2*this%M+1):(3*this%M),:)
            M_u = eta1*I_j((3*this%M+1):(4*this%M),:)


        elseif(this%Problem_Type == 3) then !! Only IBC
            allocate(I_j(2*this%M,0:N_order))
            !            I_j = Gauss_Elimination(ZY,EH,2*this%M)
            if(.not. bouncing_flag) then
                allocate(this%MoM_matrix_inverted(2*this%M,2*this%M))
                this%MoM_matrix_inverted = Mat_Inv(this%MoM_ZY(:,:,0))
            endif
            I_j(:,0) = matmul(this%MoM_matrix_inverted,EH(:,0))
            if(Wideband_Type > 0) then
                allocate(b_n(2*this%M))
                do nn = 1,N_order
                    b_n = EH(:,nn)
                    do ii = 1,nn
                        b_n = b_n - dble(n_X_i_matrix(nn,ii))*matmul(This%MoM_ZY(:,:,ii),I_j(:,nn-ii))
                    enddo
                    I_j(:,nn) = matmul(this%MoM_matrix_inverted,b_n)/(dble(n_X_i_matrix(nn,1)))
                enddo
                deallocate(b_n)
            endif
            I_v = I_j(1:this%M,:)
            I_u = I_j((this%M+1):2*this%M,:)
            do nn = 0,N_order
                M_u(:,nn) = eta1*(this%eta_vv*I_v(:,nn) + this%eta_vu*I_u(:,nn)) !! M_t
                M_v(:,nn) = -eta1*(this%eta_uu*I_u(:,nn) + this%eta_uv*I_v(:,nn)) !! M_z
            enddo
        !        write(*,*) abs(I_v(1,:))
        !        write(*,*) abs(I_u(1,:))
        elseif(this%Problem_Type == 1) then !! PEC
            allocate(I_j(2*this%M,0:N_order))
            !            I_j = Gauss_Elimination(ZY,EH,2*this%M)
            if(.not. bouncing_flag) then
                allocate(this%MoM_matrix_inverted(2*this%M,2*this%M))
                this%MoM_matrix_inverted = Mat_Inv(this%MoM_ZY(:,:,0))
            endif
            I_j(:,0) = matmul(this%MoM_matrix_inverted,EH(:,0))
            if(Wideband_Type > 0) then
                allocate(b_n(2*this%M))
                do nn = 1,N_order
                    b_n = EH(:,nn)
                    do ii = 1,nn
                        b_n = b_n - dble(n_X_i_matrix(nn,ii))*matmul(This%MoM_ZY(:,:,ii),I_j(:,nn-ii))
                    enddo
                    I_j(:,nn) = matmul(this%MoM_matrix_inverted,b_n)/(dble(n_X_i_matrix(nn,1)))
                enddo
                deallocate(b_n)
            endif

            I_v = I_j(1:this%M,:)
            I_u = I_j((this%M+1):2*this%M,:)
!        do nn = 0,N_order
!            write(*,*) nn,sum(abs(I_v(:,nn)))
!            write(*,*) nn,sum(abs(I_u(:,nn)))
!        enddo
             M_u = 0.d0
            M_v = 0.d0
        elseif(this%Problem_Type == 2) then !! PMC
            allocate(I_j(2*this%M,0:N_order))
            !            I_j = Gauss_Elimination(ZY,EH,2*this%M)
            if(.not. bouncing_flag) then
                allocate(this%MoM_matrix_inverted(2*this%M,2*this%M))
                this%MoM_matrix_inverted = Mat_Inv(this%MoM_ZY(:,:,0))
            endif
            I_j(:,0) = matmul(this%MoM_matrix_inverted,EH(:,0))
            if(Wideband_Type > 0) then
                allocate(b_n(2*this%M))
                do nn = 1,N_order
                    b_n = EH(:,nn)
                    do ii = 1,nn
                        b_n = b_n - dble(n_X_i_matrix(nn,ii))*matmul(This%MoM_ZY(:,:,ii),I_j(:,nn-ii))
                    enddo
                    I_j(:,nn) = matmul(this%MoM_matrix_inverted,b_n)/(dble(n_X_i_matrix(nn,1)))
                enddo
                deallocate(b_n)
            endif
            I_v = 0.d0
            I_u = 0.d0
            M_u = I_j((this%M+1):2*this%M,:)
            M_v = I_j(1:this%M,:)
        !        write(*,*) abs(M_v(1,:))
        !        write(*,*) abs(M_u(1,:))
        endif



        deallocate(I_j,E,H,EH)
    !    if(.not. bouncing_flag) then
    !    !        deallocate(Z,Y_f,Y,Z_f)
    !    !        deallocate(Z_d,Y_d,Z_f_d,Y_f_d)
    !        deallocate(ZY)
    !    endif



    end subroutine eval_current_MoM

    subroutine eval_IBC_impedance_matrices(this,eta_s,eta_f)
        type(Scatterer) :: this
        complex*16,allocatable,intent(inout) :: eta_s(:,:),eta_f(:,:)
        integer :: ii

        allocate(eta_s(2*this%M,2*this%M),eta_f(2*this%M,2*this%M))
        eta_s = 0.d0
        eta_f = 0.d0

        do ii = 1,this%M
            eta_s(ii,ii) = this%eta_vv(ii)
            eta_s(ii,ii+this%M) = this%eta_vu(ii)
            eta_s(ii+this%M,ii) = this%eta_uv(ii)
            eta_s(ii+this%M,ii+this%M) = this%eta_uu(ii)

            eta_f(ii,ii) = this%eta_uv(ii)
            eta_f(ii,ii+this%M) = this%eta_uu(ii)
            eta_f(ii+this%M,ii) = -this%eta_vv(ii)
            eta_f(ii+this%M,ii+this%M) = -this%eta_vu(ii)


        !            eta_f(ii,ii) = this%eta_vv(ii)
        !            eta_f(ii,ii+this%M) = this%eta_vu(ii)
        !            eta_f(ii+this%M,ii) = this%eta_uv(ii)
        !            eta_f(ii+this%M,ii+this%M) = this%eta_uu(ii)
        !
        !            eta_s(ii,ii) = this%eta_uv(ii)
        !            eta_s(ii,ii+this%M) = this%eta_uu(ii)
        !            eta_s(ii+this%M,ii) = -this%eta_vv(ii)
        !            eta_s(ii+this%M,ii+this%M) = -this%eta_vu(ii)

        enddo

    end subroutine eval_IBC_impedance_matrices

    subroutine set_MoM_bouncing_matrices(this,Scat_dist,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az_l,Scatterers_Lib)
        type(Scatterer),intent(inout) :: this,Scat_dist
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        complex*16,allocatable,intent(inout) :: Zzz(:,:,:),Zzt(:,:,:),Ztz(:,:,:),Ztt(:,:,:),&
        Yzt(:,:,:),Ytz(:,:,:),Ytt(:,:,:)
        complex*16 :: k_med,eta_med
        real(8) :: TOL1,az_l
        real(8),dimension(0:2) :: lim
        real(8),dimension(12) :: var
        complex*16 :: k_z,k_r
        type(multival) :: Integ
        integer :: i,j,t_outer_region,t_inner_region,s_outer_region,s_inner_region
        integer :: i_cnt,j_cnt,t_status,s_status

        TOL1 = 1.d-6
        !    k_med = k0*ak_l
        !    eta_med = eta1

        k_z = az_l*k0
        !    k_r = sqrt(k_med**2.d0 - k_z**2.d0)
        !    k_r = ar_l*k0
        !        write(*,*) k1,k_r,k_z
        Zzt = 0.d0
        Zzz = 0.d0
        Ztz = 0.d0
        Ztt = 0.d0
        Yzt = 0.d0
        Ytz = 0.d0
        Ytt = 0.d0

        i_cnt = 1

        do i=1,Scat_dist%M
            t_status =Scat_dist%testing_pt_status(i,1)
            if(t_status == -1) then
                cycle
            endif
            t_outer_region = Scat_dist%testing_pt_status(i,3)
            t_inner_region = Scat_dist%testing_pt_status(i,2)
            j_cnt = 1
            do j= 1,this%M
                s_status =this%testing_pt_status(j,1)
                if(s_status == -1) then
                    cycle
                endif
                s_outer_region = this%testing_pt_status(j,3)
                s_inner_region = this%testing_pt_status(j,2)

                if(s_outer_region == t_outer_region .or. s_outer_region == &
                t_inner_region .or. s_inner_region == t_outer_region) then
                    if(s_outer_region == t_outer_region) then
                        if(t_outer_region == 0) then
                            k_r = ar*k0
                            eta_med = eta1
                            k_med = ak*k0
                        else
                            k_r = k0*Scatterers_Lib(t_outer_region)%ar_local
                            k_med = k0*Scatterers_Lib(t_outer_region)%ak_local
                            eta_med = Scatterers_Lib(t_outer_region)%eta_local
                        endif
                    elseif(s_outer_region == t_inner_region) then

                        k_r = k0*Scatterers_Lib(s_outer_region)%ar_local
                        k_med = k0*Scatterers_Lib(s_outer_region)%ak_local
                        eta_med = Scatterers_Lib(s_outer_region)%eta_local

                    elseif(s_inner_region == t_outer_region) then
                        if(t_outer_region == 0) then
                            k_r = ar*k0
                            k_med = ak*k0
                            eta_med = eta1
                        else
                            k_r = k0*Scatterers_Lib(t_outer_region)%ar_local
                            k_med = k0*Scatterers_Lib(t_outer_region)%ak_local
                            eta_med = Scatterers_Lib(t_outer_region)%eta_local
                        endif
                    endif
                    lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                    var = (/dble(i),dble(j),real(k_r),aimag(k_r),Scat_dist%tang_u(i)%v(1),&
                    Scat_dist%tang_u(i)%v(2),this%tang_u(j)%v(1),&
                    this%tang_u(j)%v(2),&
                    Scat_dist%testing_pt_MoM(i)%v(1),Scat_dist%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                    this%testing_pt_MoM(j)%v(2)/)
                    Integ = Integrate(Fext2,1,lim,12,var,TOL1)

                    !                Zzz(i,j) = delta_n(j)*h0
                    !                Zzz(j,i) = delta_n(i)*h0
                    Zzz(i_cnt,j_cnt,0) = k_r**2.d0*eta_med/(4.d0*k_med)*Integ%f(1)

                    !                Zzz(j,i) = Zzz(i,j)

                    !                Ztz(i,j) = delta_n(j)*h1/R*dot(tang_v(i),R_v)
                    !                Ztz(j,i) = -delta_n(i)*h1/R*dot(tang_v(j),R_v)
                    Ztz(i_cnt,j_cnt,0) = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Integ%f(2)

                    !                Ztz(j,i) = -Integ%f(3)

                    !                Zzt(i,j) = delta_n(j)*h1/R*dot(tang_v(j),R_v)
                    !                Zzt(j,i) = -delta_n(i)*h1/R*dot(tang_v(i),R_v)
                    Zzt(i_cnt,j_cnt,0) = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Integ%f(3)
                    !                Zzt(j,i) = -Integ%f(2)


                    Ztt(i_cnt,j_cnt,0) = eta_med/(4.d0*k_med)*( (k_med**2.d0*Integ%f(1)*dot(Scat_dist%tang_u(i),this%tang_u(j)) - &
                    k_r*Integ%f(6)*dot(Scat_dist%norm_v(i),this%norm_v(j))) +&
                    Integ%f(7))


                    Ytz(i_cnt,j_cnt,0) = cj*k_r/4.d0*Integ%f(4)

                    Yzt(i_cnt,j_cnt,0) = -cj*k_r/4.d0*Integ%f(5)
                    Ytt(i_cnt,j_cnt,0) = k_z/4.d0*dot(this%norm_v(j),Scat_dist%tang_u(i))*Integ%f(1)
                !            write(*,*) i,j,Zzz(i,j),Ztt(i,j)
                else
                    Zzz(i_cnt,j_cnt,0) = 0.d0
                    Ztz(i_cnt,j_cnt,0) = 0.d0
                    Zzt(i_cnt,j_cnt,0) = 0.d0
                    Ztt(i_cnt,j_cnt,0) = 0.d0
                    Ytz(i_cnt,j_cnt,0) = 0.d0
                    Yzt(i_cnt,j_cnt,0) = 0.d0
                    Ytt(i_cnt,j_cnt,0) = 0.d0
                endif

                j_cnt = j_cnt + 1
            enddo
            i_cnt = i_cnt + 1
        enddo
    !    Zzz = k_r**2.d0*eta_med/(4.d0*k_med)*Zzz
    !    Ztz = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Ztz
    !    Zzt = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Zzt
    !    Ytt = k_z/4.d0*Ytt
    end subroutine set_MoM_bouncing_matrices

    ! -----------------------------------------------------------------------
    ! Subroutine: set_MoM_BC_matrices
    ! Purpose   : Fills the 2D MoM impedance/admittance sub-matrices
    !             Zzz, Zzt, Ztz, Ztt (electric) and Yzt, Ytz, Ytt (magnetic)
    !             for a single-frequency MoM solve.
    !             out_flag selects interior (.true.) or exterior (.false.)
    !             Green's function. Uses set_MoM_matrices_r or _c depending
    !             on whether ar_l is real or complex.
    !
    ! IBC formulation reference:
    !   A. A. Kishk and P.-S. Kildal, "Electromagnetic scattering from two
    !   dimensional anisotropic impedance objects under oblique plane wave
    !   incidence," Applied Computational Electromagnetics Society Journal,
    !   vol. 10, no. 3, pp. 81-92, 1995.
    !   Note: typographical errors in the original paper have been corrected
    !   in this implementation.
    ! -----------------------------------------------------------------------
    subroutine set_MoM_BC_matrices(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ar_l,az_l,ak0_l,eta_med,out_flag)
        type(Scatterer) :: this
        real(8) :: eta_med,ar_l,az_l,sign_out_in,ak0_l
        complex*16,allocatable,intent(inout) :: Zzz(:,:),Zzt(:,:),Ztz(:,:),Ztt(:,:),Yzt(:,:),Ytz(:,:),Ytt(:,:)
        logical :: out_flag
        integer :: i,j,i_trap
        complex*16,allocatable,dimension(:) :: k0H0,k0H1,k0H2,H1,H0
        real(8) :: R,delta_trap
        type(vector) :: R_v,src_trap
        integer :: N_trap = 7
        real(8),allocatable,dimension(:) :: x_trap,w_trap

        allocate(x_trap(N_trap+1),w_trap(N_trap+1))
        w_trap = 2.d0
        w_trap(1) = 1.d0
        w_trap(N_trap+1) = 1.d0

        allocate(k0H0(0:0),k0H1(0:0),k0H2(0:0),H1(0:0),H0(0:0))
        if(out_flag) then !! the sources are intended to radiate in the outer region
            sign_out_in = 1.d0
        else !! the sources are intended to radiate in the inner region
            sign_out_in = -1.d0
        endif
        Zzz = 0.d0
        Ztz = 0.d0
        Zzt = 0.d0
        Ztt = 0.d0
        Ytz = 0.d0
        Yzt = 0.d0
        Ytt = 0.d0


        do i = 1,this%M

            Zzz(i,i) = this%delta_n(i)*ak0_l*k0*(1.d0-2.d0*cj/pi*log(this%delta_n(i)*&
            ar_l*ak0_l*k0*exp(gamma_const)/(4.d0*exp(1.d0))))

            Ztt(i,i) = Zzz(i,i) - 2.d0*ar_l*besselh2_1(this%delta_n(i)*ar_l*ak0_l*k0/2.d0)

            Ytz(i,i) = sign_out_in*0.5d0
            Yzt(i,i) = sign_out_in*(-0.5d0)

            do j=(i+1),this%M
                delta_trap = this%delta_n(j)/dble(N_trap)
                do i_trap = 1,(N_trap+1)
                    src_trap = this%testing_pt_MoM(j) - (this%delta_n(j)/2.d0-dble(i_trap-1)*delta_trap)*this%tang_u(j)
                    R_v = this%testing_pt_MoM(i) - src_trap
                    R = absolute(R_v) !! magnitude of the vector R_v
                    R_v = (1.d0/R)*R_v !! making that vector as a unit vector
                    call Hankel_derivatives(ar_l*R,ak0_l*k0,0,k0H0,k0H1,k0H2,H1)

                    Zzz(i,j) = Zzz(i,j) + w_trap(i_trap)*k0H0(0)

                    Ztz(i,j) = Ztz(i,j) + w_trap(i_trap)*k0H1(0)*dot(this%tang_u(i),R_v)

                    Zzt(i,j) = Zzt(i,j) + w_trap(i_trap)*k0H1(0)*dot(this%tang_u(j),R_v)

                    Ztt(i,j) = Ztt(i,j) + w_trap(i_trap)*((k0H0(0) - ar_l*H1(0)/R)*&
                    dot(this%tang_u(i),this%tang_u(j)) +ar_l**2.d0*k0H2(0)*dot(this%tang_u(i),R_v)*&
                    dot(this%tang_u(j),R_v))

                    Ytz(i,j) = Ytz(i,j) + w_trap(i_trap)*k0H1(0)*dot(this%norm_v(i),R_v)

                    Yzt(i,j) = Yzt(i,j) + w_trap(i_trap)*k0H1(0)*dot(this%norm_v(j),R_v)

                    Ytt(i,j) = Ytt(i,j) + w_trap(i_trap)*k0H0(0)*dot(this%tang_u(i),this%norm_v(j))
                enddo

                Zzz(i,j) = Zzz(i,j)*this%delta_n(j)/dble(2*N_trap)
                Ztz(i,j) = Ztz(i,j)*this%delta_n(j)/dble(2*N_trap)
                Zzt(i,j) = Zzt(i,j)*this%delta_n(j)/dble(2*N_trap)
                Ztt(i,j) = Ztt(i,j)*this%delta_n(j)/dble(2*N_trap)
                Ytz(i,j) = cj*ar_l/4.d0 * Ytz(i,j)*this%delta_n(j)/dble(2*N_trap)
                Yzt(i,j) = -cj*ar_l/4.d0 * Yzt(i,j)*this%delta_n(j)/dble(2*N_trap)
                Ytt(i,j) = Ytt(i,j)*this%delta_n(j)/dble(2*N_trap)




            !            R_v = vec(this%testing_pt_MoM(i)%v(1) - this%testing_pt_MoM(j)%v(1),&
            !            this%testing_pt_MoM(i)%v(2) - this%testing_pt_MoM(j)%v(2),0.d0)
            !            R = absolute(R_v) !! magnitude of the vector R_v
            !            R_v = (1.d0/R)*R_v !! making that vector as a unit vector
            !            call Hankel_derivatives(ar_l*R,k0*ak0_l,0,k0H0,k0H1,k0H2,H1)
            !
            !            Zzz(i,j) = k0H0(0)*this%delta_n(j)
            !                Zzz(j,i) = k0H0(0)*this%delta_n(i)
            !
            !            Ztz(i,j) = k0H1(0)*this%delta_n(j)*dot(this%tang_u(i),R_v)
            !                Ztz(j,i) = k0H1(0)*this%delta_n(i)*dot(this%tang_u(j),R_v)
            !
            !            Zzt(i,j) = k0H1(0)*this%delta_n(j)*dot(this%tang_u(j),R_v)
            !                Zzt(j,i) = k0H1(0)*this%delta_n(i)*dot(this%tang_u(i),R_v)
            !
            !            Ztt(i,j) = ((k0H0(0) - ar_l*H1(0)/R)*dot(this%tang_u(i),this%tang_u(j)) +&
            !                        ar_l**2.d0*k0H2(0)*dot(this%tang_u(i),R_v)*dot(this%tang_u(j),R_v))*this%delta_n(j)
            !
            !                Ztt(j,i) = ((k0H0(0) - ar_l*H1(0)/R)*dot(this%tang_u(i),this%tang_u(j)) +&
            !                            ar_l**2.d0*k0H2(0)*dot(this%tang_u(i),R_v)*dot(this%tang_u(j),R_v))*this%delta_n(i)
            !
            !
            !            Ytz(i,j) = cj*ar_l/4.d0 * k0H1(0)*this%delta_n(j)*dot(this%norm_v(i),R_v)
            !                Ytz(j,i) = cj*ar_l/4.d0 * k0H1(0)*this%delta_n(i)*dot(this%norm_v(j),R_v)
            !
            !            Yzt(i,j) = -cj*ar_l/4.d0 *k0H1(0)*this%delta_n(j)*dot(this%norm_v(j),R_v)
            !                Yzt(j,i) = -cj*ar_l/4.d0 *k0H1(0)*this%delta_n(i)*dot(this%norm_v(i),R_v)
            !
            !            Ytt(i,j) = k0H0(0)*this%delta_n(j)*dot(this%tang_u(i),this%norm_v(j))
            !                Ytt(j,i) = k0H0(0)*this%delta_n(i)*dot(this%tang_u(j),this%norm_v(i))
            enddo

            do j=1,(i-1)
                delta_trap = this%delta_n(j)/dble(N_trap)
                do i_trap = 1,(N_trap+1)
                    src_trap = this%testing_pt_MoM(j) - (this%delta_n(j)/2.d0-dble(i_trap-1)*delta_trap)*this%tang_u(j)
                    R_v = this%testing_pt_MoM(i) - src_trap
                    R = absolute(R_v) !! magnitude of the vector R_v
                    R_v = (1.d0/R)*R_v !! making that vector as a unit vector
                    call Hankel_derivatives(ar_l*R,ak0_l*k0,0,k0H0,k0H1,k0H2,H1)

                    Zzz(i,j) = Zzz(i,j) + w_trap(i_trap)*k0H0(0)

                    Ztz(i,j) = Ztz(i,j) + w_trap(i_trap)*k0H1(0)*dot(this%tang_u(i),R_v)

                    Zzt(i,j) = Zzt(i,j) + w_trap(i_trap)*k0H1(0)*dot(this%tang_u(j),R_v)

                    Ztt(i,j) = Ztt(i,j) + w_trap(i_trap)*((k0H0(0) - ar_l*H1(0)/R)*&
                    dot(this%tang_u(i),this%tang_u(j)) +ar_l**2.d0*k0H2(0)*dot(this%tang_u(i),R_v)*&
                    dot(this%tang_u(j),R_v))

                    Ytz(i,j) = Ytz(i,j) + w_trap(i_trap)*k0H1(0)*dot(this%norm_v(i),R_v)

                    Yzt(i,j) = Yzt(i,j) + w_trap(i_trap)*k0H1(0)*dot(this%norm_v(j),R_v)

                    Ytt(i,j) = Ytt(i,j) + w_trap(i_trap)*k0H0(0)*dot(this%tang_u(i),this%norm_v(j))
                enddo

                Zzz(i,j) = Zzz(i,j)*this%delta_n(j)/dble(2*N_trap)
                Ztz(i,j) = Ztz(i,j)*this%delta_n(j)/dble(2*N_trap)
                Zzt(i,j) = Zzt(i,j)*this%delta_n(j)/dble(2*N_trap)
                Ztt(i,j) = Ztt(i,j)*this%delta_n(j)/dble(2*N_trap)
                Ytz(i,j) = cj*ar_l/4.d0 * Ytz(i,j)*this%delta_n(j)/dble(2*N_trap)
                Yzt(i,j) = -cj*ar_l/4.d0 * Yzt(i,j)*this%delta_n(j)/dble(2*N_trap)
                Ytt(i,j) = Ytt(i,j)*this%delta_n(j)/dble(2*N_trap)

            enddo

        enddo
        Zzz = Zzz*ar_l**2.d0*eta_med/(4.d0)
        Ztz = -cj*ar_l*az_l*eta_med/(4.d0) * Ztz
        Zzt = -cj*ar_l*az_l*eta_med/(4.d0) * Zzt
        Ztt = eta_med/4.d0 * Ztt

        !    Ytz = Ytz
        !    Yzt =  Yzt
        Ytt = az_l/4.d0*Ytt

        deallocate(k0H0,k0H1,k0H2,H1,H0)
        deallocate(x_trap,w_trap)
    end subroutine set_MoM_BC_matrices

    subroutine set_MoM_AWE_matrices_c(this,N_order,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ar_l,az_l,eta_med,out_flag)
        type(Scatterer) :: this
        integer :: N_order
        real(8) :: az_l,sign_out_in
        complex*16 :: eta_med,ar_l,arg,const_Ztt,const_Zzz
        complex*16,allocatable,intent(inout) :: Zzz(:,:,:),Zzt(:,:,:),Ztz(:,:,:),Ztt(:,:,:),Yzt(:,:,:),Ytz(:,:,:),Ytt(:,:,:)
        logical :: out_flag
        integer :: i,j,nn,i_trap,ii
        complex*16,allocatable,dimension(:) :: k0H0,k0H1,k0H2,H1,H0
        real(8) :: R,delta_trap
        type(vector) :: R_v,src_trap
        integer :: N_trap = 7
        real(8),allocatable,dimension(:) :: x_trap,w_trap

        allocate(x_trap(N_trap+1),w_trap(N_trap+1))
        w_trap = 2.d0
        w_trap(1) = 1.d0
        w_trap(N_trap+1) = 1.d0
        allocate(k0H0(0:N_order),k0H1(0:N_order),k0H2(0:N_order),H1(0:N_order),H0(0:N_order))
        if(out_flag) then !! the sources are intended to radiate in the outer region
            sign_out_in = 1.d0
        else !! the sources are intended to radiate in the inner region
            sign_out_in = -1.d0
        endif


        do i = 1,this%M
            arg = ar_l*k0*this%delta_n(i)/2.d0
            call Hankel_derivatives2(arg,N_order,H0,H1)
            const_Zzz = this%delta_n(i)
            !        const_Ztt =  ar_l**2.d0*this%delta_n(i)
            const_Ztt = 2.d0*ar_l/k0

            do nn = 1,N_order
                Zzz(i,i,nn) = const_Zzz*H0(nn-1)
                !            Ztt(i,i,nn) = Zzz(i,i,nn)-const_Ztt*H1(nn)

                Ztt(i,i,nn) = 0.d0
                do ii=0,nn
                    Ztt(i,i,nn) = Ztt(i,i,nn) + (-1.d0)**dble(nn-ii)*dble(n_P_i_matrix(nn,ii))*&
                    (arg**dble(ii)*H1(ii) - dble(ii)*arg**dble(ii-1)*H0(ii) )
                enddo

                Ztt(i,i,nn) = Zzz(i,i,nn) - const_Ztt*Ztt(i,i,nn)
                const_Zzz = const_Zzz * ar_l* this%delta_n(i)/2.d0
                !            const_Ztt = const_Ztt * ar_l* this%delta_n(i)/2.d0
                const_Ztt = const_Ztt/k0
            enddo
            Ytz(i,i,1:N_order) = 0.0d0
            Yzt(i,i,1:N_order) = 0.0d0



            do j=(i+1),this%M
                delta_trap = this%delta_n(j)/dble(N_trap)
                do i_trap = 1,(N_trap+1)
                    src_trap = this%testing_pt_MoM(j) - (this%delta_n(j)/2.d0-dble(i_trap-1)*delta_trap)*this%tang_u(j)
                    R_v = this%testing_pt_MoM(i) - src_trap
                    R = absolute(R_v) !! magnitude of the vector R_v
                    R_v = (1.d0/R)*R_v !! making that vector as a unit vector
                    call Hankel_derivatives(ar_l*R,k0,N_order,k0H0,k0H1,k0H2,H1)

                    Zzz(i,j,1:N_order) = Zzz(i,j,1:N_order) + w_trap(i_trap)*k0H0(1:N_order)

                    Ztz(i,j,1:N_order) = Ztz(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%tang_u(i),R_v)

                    Zzt(i,j,1:N_order) = Zzt(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%tang_u(j),R_v)

                    Ztt(i,j,1:N_order) = Ztt(i,j,1:N_order) + w_trap(i_trap)*((k0H0(1:N_order) - ar_l*H1(1:N_order)/R)*&
                    dot(this%tang_u(i),this%tang_u(j)) +ar_l**2.d0*k0H2(1:N_order)*dot(this%tang_u(i),R_v)*&
                    dot(this%tang_u(j),R_v))

                    Ytz(i,j,1:N_order) = Ytz(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%norm_v(i),R_v)

                    Yzt(i,j,1:N_order) = Yzt(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%norm_v(j),R_v)

                    Ytt(i,j,1:N_order) = Ytt(i,j,1:N_order) + w_trap(i_trap)*k0H0(1:N_order)*dot(this%tang_u(i),this%norm_v(j))
                enddo

                Zzz(i,j,1:N_order) = Zzz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ztz(i,j,1:N_order) = Ztz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Zzt(i,j,1:N_order) = Zzt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ztt(i,j,1:N_order) = Ztt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ytz(i,j,1:N_order) = cj*ar_l/4.d0 * Ytz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Yzt(i,j,1:N_order) = -cj*ar_l/4.d0 * Yzt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ytt(i,j,1:N_order) = Ytt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)


            enddo
            do j=1,(i-1)
                delta_trap = this%delta_n(j)/dble(N_trap)
                do i_trap = 1,(N_trap+1)
                    src_trap = this%testing_pt_MoM(j) - (this%delta_n(j)/2.d0-dble(i_trap-1)*delta_trap)*this%tang_u(j)
                    R_v = this%testing_pt_MoM(i) - src_trap
                    R = absolute(R_v) !! magnitude of the vector R_v
                    R_v = (1.d0/R)*R_v !! making that vector as a unit vector
                    call Hankel_derivatives(ar_l*R,k0,N_order,k0H0,k0H1,k0H2,H1)

                    Zzz(i,j,1:N_order) = Zzz(i,j,1:N_order) + w_trap(i_trap)*k0H0(1:N_order)

                    Ztz(i,j,1:N_order) = Ztz(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%tang_u(i),R_v)

                    Zzt(i,j,1:N_order) = Zzt(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%tang_u(j),R_v)

                    Ztt(i,j,1:N_order) = Ztt(i,j,1:N_order) + w_trap(i_trap)*((k0H0(1:N_order) - ar_l*H1(1:N_order)/R)*&
                    dot(this%tang_u(i),this%tang_u(j)) +ar_l**2.d0*k0H2(1:N_order)*dot(this%tang_u(i),R_v)*&
                    dot(this%tang_u(j),R_v))

                    Ytz(i,j,1:N_order) = Ytz(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%norm_v(i),R_v)

                    Yzt(i,j,1:N_order) = Yzt(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%norm_v(j),R_v)

                    Ytt(i,j,1:N_order) = Ytt(i,j,1:N_order) + w_trap(i_trap)*k0H0(1:N_order)*dot(this%tang_u(i),this%norm_v(j))
                enddo

                Zzz(i,j,1:N_order) = Zzz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ztz(i,j,1:N_order) = Ztz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Zzt(i,j,1:N_order) = Zzt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ztt(i,j,1:N_order) = Ztt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ytz(i,j,1:N_order) = cj*ar_l/4.d0 * Ytz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Yzt(i,j,1:N_order) = -cj*ar_l/4.d0 * Yzt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ytt(i,j,1:N_order) = Ytt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)

            enddo

        enddo
        Zzz(:,:,1:N_order) = Zzz(:,:,1:N_order)*ar_l**2.d0*eta_med/4.d0
        Ztz(:,:,1:N_order) = -cj*ar_l*az_l*eta_med/4.d0 * Ztz(:,:,1:N_order)
        Zzt(:,:,1:N_order) = -cj*ar_l*az_l*eta_med/4.d0 * Zzt(:,:,1:N_order)
        Ztt(:,:,1:N_order) = eta_med/4.d0 * Ztt(:,:,1:N_order)

        Ytt(:,:,1:N_order) = az_l/4.d0*Ytt(:,:,1:N_order)

        deallocate(k0H0,k0H1,k0H2,H1,H0)
        deallocate(x_trap,w_trap)
    end subroutine set_MoM_AWE_matrices_c

    ! -----------------------------------------------------------------------
    ! Subroutine: set_MoM_AWE_matrices_r  (real wavenumber)
    ! Purpose   : Builds N_order+1 Taylor-series coefficients of the MoM
    !             sub-matrices Zzz..Ytt with respect to the normalised
    !             frequency variable kappa = k/k0. Each matrix has dimensions
    !             (M, M, 0:N_order). Uses Hankel_derivatives_r for the
    !             frequency derivatives of the 2D kernel integrands.
    ! -----------------------------------------------------------------------
    subroutine set_MoM_AWE_matrices_r(this,N_order,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ar_l,az_l,eta_med,out_flag)
        type(Scatterer) :: this
        integer :: N_order
        real(8) :: eta_med,ar_l,az_l,sign_out_in
        complex*16,allocatable,intent(inout) :: Zzz(:,:,:),Zzt(:,:,:),Ztz(:,:,:),Ztt(:,:,:),Yzt(:,:,:),Ytz(:,:,:),Ytt(:,:,:)
        logical :: out_flag
        integer :: i,j,nn,i_trap,ii
        complex*16,allocatable,dimension(:) :: k0H0,k0H1,k0H2,H1,H0
        real(8) :: const_Zzz,const_Ztt,R,delta_trap,arg
        type(vector) :: R_v,src_trap
        integer :: N_trap = 7
        real(8),allocatable,dimension(:) :: x_trap,w_trap

        allocate(x_trap(N_trap+1),w_trap(N_trap+1))
        w_trap = 2.d0
        w_trap(1) = 1.d0
        w_trap(N_trap+1) = 1.d0
        allocate(k0H0(0:N_order),k0H1(0:N_order),k0H2(0:N_order),H1(0:N_order),H0(0:N_order))
        if(out_flag) then !! the sources are intended to radiate in the outer region
            sign_out_in = 1.d0
        else !! the sources are intended to radiate in the inner region
            sign_out_in = -1.d0
        endif


        do i = 1,this%M
            arg = ar_l*k0*this%delta_n(i)/2.d0
            call Hankel_derivatives2(arg,N_order,H0,H1)
            const_Zzz = this%delta_n(i)
            !        const_Ztt =  ar_l**2.d0*this%delta_n(i)
            const_Ztt = 2.d0*ar_l/k0

            do nn = 1,N_order
                Zzz(i,i,nn) = const_Zzz*H0(nn-1)
                !            Ztt(i,i,nn) = Zzz(i,i,nn)-const_Ztt*H1(nn)

                Ztt(i,i,nn) = 0.d0
                do ii=0,nn
                    Ztt(i,i,nn) = Ztt(i,i,nn) + (-1.d0)**dble(nn-ii)*dble(n_P_i_matrix(nn,ii))*&
                    (arg**dble(ii)*H1(ii) - dble(ii)*arg**dble(ii-1)*H0(ii) )
                enddo

                Ztt(i,i,nn) = Zzz(i,i,nn) - const_Ztt*Ztt(i,i,nn)
                const_Zzz = const_Zzz * ar_l* this%delta_n(i)/2.d0
                !            const_Ztt = const_Ztt * ar_l* this%delta_n(i)/2.d0
                const_Ztt = const_Ztt/k0
            enddo
            Ytz(i,i,1:N_order) = 0.0d0
            Yzt(i,i,1:N_order) = 0.0d0


            do j=(i+1),this%M
                delta_trap = this%delta_n(j)/dble(N_trap)
                do i_trap = 1,(N_trap+1)
                    src_trap = this%testing_pt_MoM(j) - (this%delta_n(j)/2.d0-dble(i_trap-1)*delta_trap)*this%tang_u(j)
                    R_v = this%testing_pt_MoM(i) - src_trap
                    R = absolute(R_v) !! magnitude of the vector R_v
                    R_v = (1.d0/R)*R_v !! making that vector as a unit vector
                    call Hankel_derivatives(ar_l*R,k0,N_order,k0H0,k0H1,k0H2,H1)

                    Zzz(i,j,1:N_order) = Zzz(i,j,1:N_order) + w_trap(i_trap)*k0H0(1:N_order)

                    Ztz(i,j,1:N_order) = Ztz(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%tang_u(i),R_v)

                    Zzt(i,j,1:N_order) = Zzt(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%tang_u(j),R_v)

                    Ztt(i,j,1:N_order) = Ztt(i,j,1:N_order) + w_trap(i_trap)*((k0H0(1:N_order) - ar_l*H1(1:N_order)/R)*&
                    dot(this%tang_u(i),this%tang_u(j)) +ar_l**2.d0*k0H2(1:N_order)*dot(this%tang_u(i),R_v)*&
                    dot(this%tang_u(j),R_v))

                    Ytz(i,j,1:N_order) = Ytz(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%norm_v(i),R_v)

                    Yzt(i,j,1:N_order) = Yzt(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%norm_v(j),R_v)

                    Ytt(i,j,1:N_order) = Ytt(i,j,1:N_order) + w_trap(i_trap)*k0H0(1:N_order)*dot(this%tang_u(i),this%norm_v(j))
                enddo

                Zzz(i,j,1:N_order) = Zzz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ztz(i,j,1:N_order) = Ztz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Zzt(i,j,1:N_order) = Zzt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ztt(i,j,1:N_order) = Ztt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ytz(i,j,1:N_order) = cj*ar_l/4.d0 * Ytz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Yzt(i,j,1:N_order) = -cj*ar_l/4.d0 * Yzt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ytt(i,j,1:N_order) = Ytt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)

            !            R_v = vec(this%testing_pt_MoM(i)%v(1) - this%testing_pt_MoM(j)%v(1),&
            !            this%testing_pt_MoM(i)%v(2) - this%testing_pt_MoM(j)%v(2),0.d0)
            !            R = absolute(R_v) !! magnitude of the vector R_v
            !            R_v = (1.d0/R)*R_v !! making that vector as a unit vector
            !            call Hankel_derivatives(ar_l*R,k0,N_order,k0H0,k0H1,k0H2,H1)
            !
            !            Zzz(i,j,1:N_order) = k0H0(1:N_order)*this%delta_n(j)
            !    !                Zzz(j,i,1:N_order) = k0H0(1:N_order)*this%delta_n(i)

            !            Ztz(i,j,1:N_order) = k0H1(1:N_order)*this%delta_n(j)*dot(this%tang_u(i),R_v)
            !    !                Ztz(j,i,1:N_order) = k0H1(1:N_order)*this%delta_n(i)*dot(this%tang_u(j),R_v)

            !            Zzt(i,j,1:N_order) = k0H1(1:N_order)*this%delta_n(j)*dot(this%tang_u(j),R_v)
            !    !                Zzt(j,i,1:N_order) = k0H1(1:N_order)*this%delta_n(i)*dot(this%tang_u(i),R_v)

            !            Ztt(i,j,1:N_order) = ((k0H0(1:N_order) - ar_l*H1(1:N_order)/R)*dot(this%tang_u(i),this%tang_u(j)) +&
            !                        ar_l**2.d0*k0H2(1:N_order)*dot(this%tang_u(i),R_v)*dot(this%tang_u(j),R_v))*this%delta_n(j)
            !
            !    !                Ztt(j,i,1:N_order) = ((k0H0(1:N_order) - ar_l*H1(1:N_order)/R)*dot(this%tang_u(i),this%tang_u(j)) +&
            !    !                            ar_l**2.d0*k0H2(1:N_order)*dot(this%tang_u(i),R_v)*dot(this%tang_u(j),R_v))*this%delta_n(i)

            !            Ytz(i,j,1:N_order) = k0H1(1:N_order)*this%delta_n(j)*dot(this%norm_v(i),R_v)
            !    !                Ytz(j,i,1:N_order) = k0H1(1:N_order)*this%delta_n(i)*dot(this%norm_v(j),R_v)
            !
            !            Yzt(i,j,1:N_order) = k0H1(1:N_order)*this%delta_n(j)*dot(this%norm_v(j),R_v)
            !    !                Yzt(j,i,1:N_order) = k0H1(1:N_order)*this%delta_n(i)*dot(this%norm_v(i),R_v)
            !
            !            Ytt(i,j,1:N_order) = k0H0(1:N_order)*this%delta_n(j)*dot(this%tang_u(i),this%norm_v(j))
            !    !                Ytt(j,i,1:N_order) = k0H0(1:N_order)*this%delta_n(i)*dot(this%tang_u(j),this%norm_v(i))



            enddo
            do j=1,(i-1)
                delta_trap = this%delta_n(j)/dble(N_trap)
                do i_trap = 1,(N_trap+1)
                    src_trap = this%testing_pt_MoM(j) - (this%delta_n(j)/2.d0-dble(i_trap-1)*delta_trap)*this%tang_u(j)
                    R_v = this%testing_pt_MoM(i) - src_trap
                    R = absolute(R_v) !! magnitude of the vector R_v
                    R_v = (1.d0/R)*R_v !! making that vector as a unit vector
                    call Hankel_derivatives(ar_l*R,k0,N_order,k0H0,k0H1,k0H2,H1)

                    Zzz(i,j,1:N_order) = Zzz(i,j,1:N_order) + w_trap(i_trap)*k0H0(1:N_order)

                    Ztz(i,j,1:N_order) = Ztz(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%tang_u(i),R_v)

                    Zzt(i,j,1:N_order) = Zzt(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%tang_u(j),R_v)

                    Ztt(i,j,1:N_order) = Ztt(i,j,1:N_order) + w_trap(i_trap)*((k0H0(1:N_order) - ar_l*H1(1:N_order)/R)*&
                    dot(this%tang_u(i),this%tang_u(j)) +ar_l**2.d0*k0H2(1:N_order)*dot(this%tang_u(i),R_v)*&
                    dot(this%tang_u(j),R_v))

                    Ytz(i,j,1:N_order) = Ytz(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%norm_v(i),R_v)

                    Yzt(i,j,1:N_order) = Yzt(i,j,1:N_order) + w_trap(i_trap)*k0H1(1:N_order)*dot(this%norm_v(j),R_v)

                    Ytt(i,j,1:N_order) = Ytt(i,j,1:N_order) + w_trap(i_trap)*k0H0(1:N_order)*dot(this%tang_u(i),this%norm_v(j))
                enddo

                Zzz(i,j,1:N_order) = Zzz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ztz(i,j,1:N_order) = Ztz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Zzt(i,j,1:N_order) = Zzt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ztt(i,j,1:N_order) = Ztt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ytz(i,j,1:N_order) = cj*ar_l/4.d0 * Ytz(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Yzt(i,j,1:N_order) = -cj*ar_l/4.d0 * Yzt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)
                Ytt(i,j,1:N_order) = Ytt(i,j,1:N_order)*this%delta_n(j)/dble(2*N_trap)

            enddo

        enddo
        Zzz(:,:,1:N_order) = Zzz(:,:,1:N_order)*ar_l**2.d0*eta_med/4.d0
        Ztz(:,:,1:N_order) = -cj*ar_l*az_l*eta_med/4.d0 * Ztz(:,:,1:N_order)
        Zzt(:,:,1:N_order) = -cj*ar_l*az_l*eta_med/4.d0 * Zzt(:,:,1:N_order)
        Ztt(:,:,1:N_order) = eta_med/4.d0 * Ztt(:,:,1:N_order)

        !    Ytz(:,:,1:N_order) = cj*ar_l/4.d0 * Ytz(:,:,1:N_order)
        !    Yzt(:,:,1:N_order) = -cj*ar_l/4.d0 * Yzt(:,:,1:N_order)
        Ytt(:,:,1:N_order) = az_l/4.d0*Ytt(:,:,1:N_order)

        deallocate(k0H0,k0H1,k0H2,H1,H0)
        deallocate(x_trap,w_trap)
    end subroutine set_MoM_AWE_matrices_r

    subroutine set_MoM_matrices_c(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az_l,ar_l,ak_l,eta_med,out_flag)
        type(Scatterer) :: this
        complex*16 :: k_med,eta_med
        complex*16,allocatable,intent(inout) :: Zzz(:,:,:),Zzt(:,:,:),Ztz(:,:,:),Ztt(:,:,:),Yzt(:,:,:),Ytz(:,:,:),Ytt(:,:,:)
        logical :: out_flag
        real(8) :: sign_out_in
        integer :: i,j
        real(8) :: TOL1
        real(8),dimension(0:2) :: lim
        real(8),dimension(12) :: var
        complex*16 :: k_r,ar_l,ak_l,h0_x,h1_x
        real(8) :: k_z,az_l
        type(multival) :: Integ
        TOL1 = 1.d-6


        if(out_flag) then !! the sources are intended to radiate in the outer region
            sign_out_in = 1.d0
        else !! the sources are intended to radiate in the inner region
            sign_out_in = -1.d0
        endif


        k_z = az_l*k0
        k_r = ar_l*k0
        k_med = ak_l*k0
        !        write(*,*) k1,k_r,k_z
        Zzt = 0.d0
        Zzz = 0.d0
        Ztz = 0.d0
        Ztt = 0.d0
        Yzt = 0.d0
        Ytz = 0.d0
        Ytt = 0.d0

        do i=1,this%M

            Zzz(i,i,0) = this%delta_n(i)*(1.d0-2.d0*cj/pi*(log(this%delta_n(i)*exp(gamma_const)*k_r/(4.d0*exp(1.d0)) )-0.0d0) )
            call besselh2_01(k_r*this%delta_n(i)/2.d0,h0_x,h1_x)
            Ztt(i,i,0) = -k_r*eta_med/(2.d0*k_med)*h1_x + &
            k_med*eta_med*this%delta_n(i)/4.d0*(1.d0-2.d0*cj/pi*( log(this%delta_n(i)*&
            exp(gamma_const)*k_r/(4.d0*exp(1.d0)))-0.0d0))
            Ytz(i,i,0) = sign_out_in*0.5d0
            Yzt(i,i,0) = sign_out_in*(-0.5d0)


            do j=(i+1),this%M
                !                R_v = get_R(testing_pt_MoM(i),testing_pt_MoM(j))
                !                R = absolute(R_v)
                !                h0 = besselh2_0(kr*R)
                !                h1 = besselh2_1(kr*R)
                !                h2 = 2.d0/(kr*R)*h1-h0
                lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                var = (/dble(i),dble(j),real(k_r),aimag(k_r),this%tang_u(i)%v(1),this%tang_u(i)%v(2),this%tang_u(j)%v(1),&
                this%tang_u(j)%v(2),this%testing_pt_MoM(i)%v(1),this%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                this%testing_pt_MoM(j)%v(2)/)
                Integ = Integrate(Fext2,1,lim,12,var,TOL1)

                !                Zzz(i,j) = delta_n(j)*h0
                !                Zzz(j,i) = delta_n(i)*h0
                Zzz(i,j,0) = Integ%f(1)
                !                Zzz(j,i) = Zzz(i,j)

                !                Ztz(i,j) = delta_n(j)*h1/R*dot(tang_v(i),R_v)
                !                Ztz(j,i) = -delta_n(i)*h1/R*dot(tang_v(j),R_v)
                Ztz(i,j,0) = Integ%f(2)
                !                Ztz(j,i) = -Integ%f(3)

                !                Zzt(i,j) = delta_n(j)*h1/R*dot(tang_v(j),R_v)
                !                Zzt(j,i) = -delta_n(i)*h1/R*dot(tang_v(i),R_v)
                Zzt(i,j,0) = Integ%f(3)
                !                Zzt(j,i) = -Integ%f(2)


                Ztt(i,j,0) = eta_med/(4.d0*k_med)*( (k_med**2.d0*Integ%f(1)*dot(this%tang_u(i),this%tang_u(j)) - &
                k_r*Integ%f(6)*dot(this%norm_v(i),this%norm_v(j))) +&
                Integ%f(7))


                Ytz(i,j,0) = cj*k_r/4.d0*Integ%f(4)

                Yzt(i,j,0) = -cj*k_r/4.d0*Integ%f(5)
                Ytt(i,j,0) = dot(this%norm_v(j),this%tang_u(i))*Integ%f(1)


            enddo
            do j=1,(i-1)
                lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                var = (/dble(i),dble(j),real(k_r),aimag(k_r),this%tang_u(i)%v(1),this%tang_u(i)%v(2),this%tang_u(j)%v(1),&
                this%tang_u(j)%v(2),this%testing_pt_MoM(i)%v(1),this%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                this%testing_pt_MoM(j)%v(2)/)
                Integ = Integrate(Fext2,1,lim,12,var,TOL1)

                Zzz(i,j,0) = Integ%f(1)


                Ztz(i,j,0) = Integ%f(2)

                Zzt(i,j,0) = Integ%f(3)

                Ztt(i,j,0) = eta_med/(4.d0*k_med)*( (k_med**2.d0*Integ%f(1)*dot(this%tang_u(i),this%tang_u(j)) - &
                k_r*Integ%f(6)*dot(this%norm_v(i),this%norm_v(j))) +&
                Integ%f(7))


                Ytz(i,j,0) = cj*k_r/4.d0*Integ%f(4)


                Yzt(i,j,0) = -cj*k_r/4.d0*Integ%f(5)


                Ytt(i,j,0) = dot(this%norm_v(j),this%tang_u(i))*Integ%f(1)
            enddo
        enddo
        Zzz = k_r**2.d0*eta_med/(4.d0*k_med)*Zzz
        Ztz = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Ztz
        Zzt = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Zzt
        Ytt = k_z/4.d0*Ytt
    end subroutine set_MoM_matrices_c


    subroutine set_MoM_matrices_r(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az_l,ar_l,ak_l,eta_med,out_flag)
        type(Scatterer) :: this
        real(8) :: k_med,eta_med
        complex*16,allocatable,intent(inout) :: Zzz(:,:,:),Zzt(:,:,:),Ztz(:,:,:),Ztt(:,:,:),Yzt(:,:,:),Ytz(:,:,:),Ytt(:,:,:)
        logical :: out_flag
        real(8) :: sign_out_in
        integer :: i,j
        real(8) :: TOL1
        real(8),dimension(0:2) :: lim
        real(8),dimension(11) :: var
        real(8) :: k_z,k_r,ar_l,ak_l,az_l
        type(multival) :: Integ
        TOL1 = 1.d-6


        if(out_flag) then !! the sources are intended to radiate in the outer region
            sign_out_in = 1.d0
        else !! the sources are intended to radiate in the inner region
            sign_out_in = -1.d0
        endif


        k_z = az_l*k0
        k_r = ar_l*k0
        k_med = ak_l*k0
        !        write(*,*) k1,k_r,k_z
        Zzt = 0.d0
        Zzz = 0.d0
        Ztz = 0.d0
        Ztt = 0.d0
        Yzt = 0.d0
        Ytz = 0.d0
        Ytt = 0.d0

        do i=1,this%M

            Zzz(i,i,0) = this%delta_n(i)*(1.d0-2.d0*cj/pi*(log(this%delta_n(i)*exp(gamma_const)*k_r/(4.d0*exp(1.d0)) )-0.0d0) )
            Ztt(i,i,0) = -k_r*eta_med/(2.d0*k_med)*besselh2_1(k_r*this%delta_n(i)/2.d0) + &
            k_med*eta_med*this%delta_n(i)/4.d0*(1.d0-2.d0*cj/pi*( log(this%delta_n(i)*&
            exp(gamma_const)*k_r/(4.d0*exp(1.d0)))-0.0d0))
            Ytz(i,i,0) = sign_out_in*0.5d0
            Yzt(i,i,0) = sign_out_in*(-0.5d0)


            do j=(i+1),this%M
                !                R_v = get_R(testing_pt_MoM(i),testing_pt_MoM(j))
                !                R = absolute(R_v)
                !                h0 = besselh2_0(kr*R)
                !                h1 = besselh2_1(kr*R)
                !                h2 = 2.d0/(kr*R)*h1-h0
                lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                var = (/dble(i),dble(j),k_r,this%tang_u(i)%v(1),this%tang_u(i)%v(2),this%tang_u(j)%v(1),this%tang_u(j)%v(2),&
                this%testing_pt_MoM(i)%v(1),this%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                this%testing_pt_MoM(j)%v(2)/)
                Integ = Integrate(Fext1,1,lim,11,var,TOL1)

                !                Zzz(i,j) = delta_n(j)*h0
                !                Zzz(j,i) = delta_n(i)*h0
                Zzz(i,j,0) = Integ%f(1)
                !                Zzz(j,i) = Zzz(i,j)

                !                Ztz(i,j) = delta_n(j)*h1/R*dot(tang_v(i),R_v)
                !                Ztz(j,i) = -delta_n(i)*h1/R*dot(tang_v(j),R_v)
                Ztz(i,j,0) = Integ%f(2)
                !                Ztz(j,i) = -Integ%f(3)

                !                Zzt(i,j) = delta_n(j)*h1/R*dot(tang_v(j),R_v)
                !                Zzt(j,i) = -delta_n(i)*h1/R*dot(tang_v(i),R_v)
                Zzt(i,j,0) = Integ%f(3)
                !                Zzt(j,i) = -Integ%f(2)


                Ztt(i,j,0) = eta_med/(4.d0*k_med)*( (k_med**2.d0*Integ%f(1)*dot(this%tang_u(i),this%tang_u(j)) - &
                k_r*Integ%f(6)*dot(this%norm_v(i),this%norm_v(j))) +&
                Integ%f(7))


                Ytz(i,j,0) = cj*k_r/4.d0*Integ%f(4)

                Yzt(i,j,0) = -cj*k_r/4.d0*Integ%f(5)
                Ytt(i,j,0) = dot(this%norm_v(j),this%tang_u(i))*Integ%f(1)


            enddo
            do j=1,(i-1)
                lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                var = (/dble(i),dble(j),k_r,this%tang_u(i)%v(1),this%tang_u(i)%v(2),this%tang_u(j)%v(1),this%tang_u(j)%v(2),&
                this%testing_pt_MoM(i)%v(1),this%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                this%testing_pt_MoM(j)%v(2)/)
                Integ = Integrate(Fext1,1,lim,11,var,TOL1)

                Zzz(i,j,0) = Integ%f(1)


                Ztz(i,j,0) = Integ%f(2)

                Zzt(i,j,0) = Integ%f(3)

                Ztt(i,j,0) = eta_med/(4.d0*k_med)*( (k_med**2.d0*Integ%f(1)*dot(this%tang_u(i),this%tang_u(j)) - &
                k_r*Integ%f(6)*dot(this%norm_v(i),this%norm_v(j))) +&
                Integ%f(7))


                Ytz(i,j,0) = cj*k_r/4.d0*Integ%f(4)


                Yzt(i,j,0) = -cj*k_r/4.d0*Integ%f(5)


                Ytt(i,j,0) = dot(this%norm_v(j),this%tang_u(i))*Integ%f(1)
            enddo
        enddo
        Zzz = k_r**2.d0*eta_med/(4.d0*k_med)*Zzz
        Ztz = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Ztz
        Zzt = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Zzt
        Ytt = k_z/4.d0*Ytt
    end subroutine set_MoM_matrices_r


    subroutine set_MoM_matrices_separate(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az_l,out_flag,Scatterers_Lib)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        complex*16 :: k_med,eta_med
        complex*16,allocatable,intent(inout) :: Zzz(:,:,:),Zzt(:,:,:),Ztz(:,:,:),Ztt(:,:,:),Yzt(:,:,:),Ytz(:,:,:),Ytt(:,:,:)
        logical :: out_flag
        real(8) :: sign_out_in
        integer :: i,j,t_region_id,region_index,s_region_id,t_status,s_status,i_cnt,j_cnt
        real(8) :: TOL1
        real(8),dimension(0:2) :: lim
        real(8),dimension(12) :: var
        real(8) :: k_z,az_l
        complex*16 :: h0_x,h1_x,k_r
        type(multival) :: Integ
        TOL1 = 1.d-6


        if(out_flag) then !! the sources are intended to radiate in the outer region
            sign_out_in = 1.d0
            region_index = 3
        else !! the sources are intended to radiate in the inner region
            sign_out_in = -1.d0
            region_index = 2
        endif


        k_z = az_l*k0
        !    k_r = ar_l*k0
        !    k_med = ak_l*k0
        !        write(*,*) k1,k_r,k_z
        Zzt = 0.d0
        Zzz = 0.d0
        Ztz = 0.d0
        Ztt = 0.d0
        Yzt = 0.d0
        Ytz = 0.d0
        Ytt = 0.d0

        i_cnt = 1
        do i=1,this%M
            t_status = this%testing_pt_status(i,1)
            if(t_status == -1) then
                cycle
            endif
            t_region_id = this%testing_pt_status(i,region_index)
            if(t_region_id == 0) then
                k_r = ar*k0
                k_med = ak*k0
                eta_med = eta1
            else
                k_r = k0*Scatterers_Lib(t_region_id)%ar_local
                k_med = k0*Scatterers_Lib(t_region_id)%ak_local
                eta_med = Scatterers_Lib(t_region_id)%eta_local

            endif
!            write(*,*) k_r,eta_med,&
!            log(this%delta_n(i)*exp(gamma_const)*k_r/(4.d0*exp(1.d0)))

            Zzz(i_cnt,i_cnt,0) = this%delta_n(i)*((1.d0,0.d0)-((0.d0,2.d0)/pi)*(log(this%delta_n(i)*&
            exp(gamma_const)*k_r/((4.d0,0.d0)*exp((1.d0,0.d0))) )) )
            call besselh2_01(k_r*this%delta_n(i)/(2.d0,0.d0),h0_x,h1_x)
            Ztt(i_cnt,i_cnt,0) = -k_r*eta_med/((2.d0,0.d0)*k_med)*h1_x + &
            k_med*eta_med*this%delta_n(i)/(4.d0,0.d0)*((1.d0,0.d0)-((0.d0,2.d0)/pi)*( log(this%delta_n(i)*&
            exp(gamma_const)*k_r/((4.d0,0.d0)*exp((1.d0,0.d0))))))
            Ytz(i_cnt,i_cnt,0) = sign_out_in*0.5d0
            Yzt(i_cnt,i_cnt,0) = sign_out_in*(-0.5d0)

            j_cnt = i_cnt+1
            do j=(i+1),this%M
                s_status = this%testing_pt_status(j,1)
                if(s_status == -1) then
                    cycle
                endif
                s_region_id = this%testing_pt_status(j,region_index)
                if(s_region_id == t_region_id) then
                    lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                    var = (/dble(i),dble(j),real(k_r),aimag(k_r),this%tang_u(i)%v(1),this%tang_u(i)%v(2),this%tang_u(j)%v(1),&
                    this%tang_u(j)%v(2),&
                    this%testing_pt_MoM(i)%v(1),this%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                    this%testing_pt_MoM(j)%v(2)/)
                    Integ = Integrate(Fext2,1,lim,12,var,TOL1)

                    !                Zzz(i,j) = delta_n(j)*h0
                    !                Zzz(j,i) = delta_n(i)*h0
                    Zzz(i_cnt,j_cnt,0) = Integ%f(1)
                    !                Zzz(j,i) = Zzz(i,j)

                    !                Ztz(i,j) = delta_n(j)*h1/R*dot(tang_v(i),R_v)
                    !                Ztz(j,i) = -delta_n(i)*h1/R*dot(tang_v(j),R_v)
                    Ztz(i_cnt,j_cnt,0) = Integ%f(2)
                    !                Ztz(j,i) = -Integ%f(3)

                    !                Zzt(i,j) = delta_n(j)*h1/R*dot(tang_v(j),R_v)
                    !                Zzt(j,i) = -delta_n(i)*h1/R*dot(tang_v(i),R_v)
                    Zzt(i_cnt,j_cnt,0) = Integ%f(3)
                    !                Zzt(j,i) = -Integ%f(2)


                    Ztt(i_cnt,j_cnt,0) = eta_med/((4.d0,0.d0)*k_med)*( (k_med*k_med*&
                    Integ%f(1)*dot(this%tang_u(i),this%tang_u(j)) - &
                    k_r*Integ%f(6)*dot(this%norm_v(i),this%norm_v(j))) +&
                    Integ%f(7))


                    Ytz(i_cnt,j_cnt,0) = (0.d0,0.25d0)*k_r*Integ%f(4)

                    Yzt(i_cnt,j_cnt,0) = (0.d0,-0.25d0)*k_r*Integ%f(5)
                    Ytt(i_cnt,j_cnt,0) = dot(this%norm_v(j),this%tang_u(i))*Integ%f(1)
                endif
                j_cnt = j_cnt + 1
            enddo
            j_cnt = 1
            do j=1,(i-1)
                s_status = this%testing_pt_status(j,1)
                if(s_status == -1) then
                    cycle
                endif
                s_region_id = this%testing_pt_status(j,region_index)
                if(s_region_id == t_region_id) then



                    lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                    var = (/dble(i),dble(j),real(k_r),aimag(k_r),this%tang_u(i)%v(1),this%tang_u(i)%v(2),this%tang_u(j)%v(1),&
                    this%tang_u(j)%v(2),&
                    this%testing_pt_MoM(i)%v(1),this%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                    this%testing_pt_MoM(j)%v(2)/)
                    Integ = Integrate(Fext2,1,lim,12,var,TOL1)

                    Zzz(i_cnt,j_cnt,0) = Integ%f(1)


                    Ztz(i_cnt,j_cnt,0) = Integ%f(2)

                    Zzt(i_cnt,j_cnt,0) = Integ%f(3)

                    Ztt(i_cnt,j_cnt,0) = eta_med/((4.d0,0.d0)*k_med)*( (k_med*k_med*Integ%f(1)*&
                    dot(this%tang_u(i),this%tang_u(j)) - &
                    k_r*Integ%f(6)*dot(this%norm_v(i),this%norm_v(j))) +&
                    Integ%f(7))


                    Ytz(i_cnt,j_cnt,0) = (0.d0,0.25d0)*k_r*Integ%f(4)


                    Yzt(i_cnt,j_cnt,0) = (0.d0,-0.25d0)*k_r*Integ%f(5)


                    Ytt(i_cnt,j_cnt,0) = dot(this%norm_v(j),this%tang_u(i))*Integ%f(1)
                endif
                j_cnt = j_cnt + 1
            enddo
            Zzz(i_cnt,:,0) = (0.25d0,0.d0)*k_r*k_r*eta_med/(k_med)*Zzz(i_cnt,:,0)
            Ztz(i_cnt,:,0) = k_r*k_z*eta_med/((0.d0,4.d0)*k_med)*Ztz(i_cnt,:,0)
            Zzt(i_cnt,:,0) = k_r*k_z*eta_med/((0.d0,4.d0)*k_med)*Zzt(i_cnt,:,0)
            Ytt(i_cnt,:,0) = (0.25d0,0.d0)*k_z*Ytt(i_cnt,:,0)
            i_cnt = i_cnt + 1
        enddo

    end subroutine set_MoM_matrices_separate


    subroutine set_MoM_matrices_separate_r(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,az_l,out_flag,Scatterers_Lib)
        type(Scatterer) :: this
        type(Scatterer),allocatable,intent(inout) :: Scatterers_Lib(:)
        real(8) :: k_med,eta_med
        complex*16,allocatable,intent(inout) :: Zzz(:,:,:),Zzt(:,:,:),Ztz(:,:,:),Ztt(:,:,:),Yzt(:,:,:),Ytz(:,:,:),Ytt(:,:,:)
        logical :: out_flag
        real(8) :: sign_out_in
        integer :: i,j,t_region_id,region_index,s_region_id,t_status,s_status,i_cnt,j_cnt
        real(8) :: TOL1
        real(8),dimension(0:2) :: lim
        real(8),dimension(11) :: var
        real(8) :: k_z,az_l,k_r
        complex*16 :: h0_x,h1_x
        type(multival) :: Integ
        TOL1 = 1.d-6


        if(out_flag) then !! the sources are intended to radiate in the outer region
            sign_out_in = 1.d0
            region_index = 3
        else !! the sources are intended to radiate in the inner region
            sign_out_in = -1.d0
            region_index = 2
        endif


        k_z = az_l*k0
        !    k_r = ar_l*k0
        !    k_med = ak_l*k0
        !        write(*,*) k1,k_r,k_z
        Zzt = 0.d0
        Zzz = 0.d0
        Ztz = 0.d0
        Ztt = 0.d0
        Yzt = 0.d0
        Ytz = 0.d0
        Ytt = 0.d0

        i_cnt = 1
        do i=1,this%M
            t_status = this%testing_pt_status(i,1)
            if(t_status == -1) then
                cycle
            endif
            t_region_id = this%testing_pt_status(i,region_index)
            if(t_region_id == 0) then
                k_r = ar*k0
                k_med = ak*k0
                eta_med = eta1
            else
                k_r = k0*real(Scatterers_Lib(t_region_id)%ar_local)
                k_med = k0*real(Scatterers_Lib(t_region_id)%ak_local)
                eta_med = real(Scatterers_Lib(t_region_id)%eta_local)
            endif


            Zzz(i_cnt,i_cnt,0) = this%delta_n(i)*(1.d0-2.d0*cj/pi*(log(this%delta_n(i)*&
            exp(gamma_const)*k_r/(4.d0*exp(1.d0)) )-0.0d0) )
            call besselh2_01(k_r*this%delta_n(i)/2.d0,h0_x,h1_x)
            Ztt(i_cnt,i_cnt,0) = -k_r*eta_med/(2.d0*k_med)*h1_x + &
            k_med*eta_med*this%delta_n(i)/4.d0*(1.d0-2.d0*cj/pi*( log(this%delta_n(i)*&
            exp(gamma_const)*k_r/(4.d0*exp(1.d0)))-0.0d0))
            Ytz(i_cnt,i_cnt,0) = sign_out_in*0.5d0
            Yzt(i_cnt,i_cnt,0) = sign_out_in*(-0.5d0)

            j_cnt = i_cnt+1
            do j=(i+1),this%M
                s_status = this%testing_pt_status(j,1)
                if(s_status == -1) then
                    cycle
                endif
                s_region_id = this%testing_pt_status(j,region_index)
                if(s_region_id == t_region_id) then
                    lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                    var = (/dble(i),dble(j),k_r,this%tang_u(i)%v(1),this%tang_u(i)%v(2),this%tang_u(j)%v(1),&
                    this%tang_u(j)%v(2),&
                    this%testing_pt_MoM(i)%v(1),this%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                    this%testing_pt_MoM(j)%v(2)/)
                    Integ = Integrate(Fext1,1,lim,11,var,TOL1)

                    !                Zzz(i,j) = delta_n(j)*h0
                    !                Zzz(j,i) = delta_n(i)*h0
                    Zzz(i_cnt,j_cnt,0) = Integ%f(1)
                    !                Zzz(j,i) = Zzz(i,j)

                    !                Ztz(i,j) = delta_n(j)*h1/R*dot(tang_v(i),R_v)
                    !                Ztz(j,i) = -delta_n(i)*h1/R*dot(tang_v(j),R_v)
                    Ztz(i_cnt,j_cnt,0) = Integ%f(2)
                    !                Ztz(j,i) = -Integ%f(3)

                    !                Zzt(i,j) = delta_n(j)*h1/R*dot(tang_v(j),R_v)
                    !                Zzt(j,i) = -delta_n(i)*h1/R*dot(tang_v(i),R_v)
                    Zzt(i_cnt,j_cnt,0) = Integ%f(3)
                    !                Zzt(j,i) = -Integ%f(2)


                    Ztt(i_cnt,j_cnt,0) = eta_med/(4.d0*k_med)*( (k_med**2.d0*Integ%f(1)*dot(this%tang_u(i),this%tang_u(j)) - &
                    k_r*Integ%f(6)*dot(this%norm_v(i),this%norm_v(j))) +&
                    Integ%f(7))


                    Ytz(i_cnt,j_cnt,0) = cj*k_r/4.d0*Integ%f(4)

                    Yzt(i_cnt,j_cnt,0) = -cj*k_r/4.d0*Integ%f(5)
                    Ytt(i_cnt,j_cnt,0) = dot(this%norm_v(j),this%tang_u(i))*Integ%f(1)
                endif
                j_cnt = j_cnt + 1
            enddo
            j_cnt = 1
            do j=1,(i-1)
                s_status = this%testing_pt_status(j,1)
                if(s_status == -1) then
                    cycle
                endif
                s_region_id = this%testing_pt_status(j,region_index)
                if(s_region_id == t_region_id) then



                    lim = (/0.d0,-this%delta_n(j)/2.d0,this%delta_n(j)/2.d0/)
                    var = (/dble(i),dble(j),k_r,this%tang_u(i)%v(1),this%tang_u(i)%v(2),this%tang_u(j)%v(1),&
                    this%tang_u(j)%v(2),&
                    this%testing_pt_MoM(i)%v(1),this%testing_pt_MoM(i)%v(2),this%testing_pt_MoM(j)%v(1),&
                    this%testing_pt_MoM(j)%v(2)/)
                    Integ = Integrate(Fext1,1,lim,11,var,TOL1)

                    Zzz(i_cnt,j_cnt,0) = Integ%f(1)


                    Ztz(i_cnt,j_cnt,0) = Integ%f(2)

                    Zzt(i_cnt,j_cnt,0) = Integ%f(3)

                    Ztt(i_cnt,j_cnt,0) = eta_med/(4.d0*k_med)*( (k_med**2.d0*Integ%f(1)*dot(this%tang_u(i),this%tang_u(j)) - &
                    k_r*Integ%f(6)*dot(this%norm_v(i),this%norm_v(j))) +&
                    Integ%f(7))


                    Ytz(i_cnt,j_cnt,0) = cj*k_r/4.d0*Integ%f(4)


                    Yzt(i_cnt,j_cnt,0) = -cj*k_r/4.d0*Integ%f(5)


                    Ytt(i_cnt,j_cnt,0) = dot(this%norm_v(j),this%tang_u(i))*Integ%f(1)
                endif
                j_cnt = j_cnt + 1
            enddo
            i_cnt = i_cnt + 1
        enddo
        Zzz = k_r**2.d0*eta_med/(4.d0*k_med)*Zzz
        Ztz = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Ztz
        Zzt = -cj*k_r*k_z*eta_med/(4.d0*k_med)*Zzt
        Ytt = k_z/4.d0*Ytt
    end subroutine set_MoM_matrices_separate_r


    function eval_surface_current_error(this) result(error)
        type(Scatterer) :: this
        complex*16,allocatable,dimension(:) :: Jz,Jt
        integer ::  fd1 = 50
        complex*16,allocatable,dimension(:) :: H_total
        integer :: j,kk
        real(8),allocatable,dimension(:) :: Jz_phase,Jt_phase
        real(8) :: t,kr,kz,phi_dev
        complex*16,allocatable,dimension(:)::h2,h2p
        integer :: n_terms
        real(8) :: error,normalize
        character*30 :: file_name

        complex*16 :: Ja_z_TM,Ja_t_TE,Ja_z_TE,Ja_z,Ja_t !! analytical surface current
        complex*16 :: sum_z_TM,sum_z_TE,sum_t_TE,cjj,const_Ja_z_TM,const_Ja_z_TE,const_Ja_t_TE
        real(8) :: x_s,y_s,phi,Ja_z_phase,Ja_t_phase

        if(plot_current /= 0) then
            file_name = 'currents_plotting'//num2str(this%region_ID,1)//'.dat'
            OPEN(fd1, FILE=trim(file_name))
        endif

        error = 0.d0
        normalize = 0.d0

        this%radius = this%a_superquad

        kz = k0*cos(theta_i)
        kr = sqrt(k1**2.d0 - kz**2.d0)

        allocate(Jz(this%M),Jt(this%M),H_total(2*this%M),Jz_phase(this%M),Jt_phase(this%M))
        !------------------------------------------------------
        !! Initialization of the exact current distribution
        const_Ja_z_TM = 2.d0*cos(alpha_i)/(pi*kr*this%radius)
        const_Ja_z_TE = -4.d0*sin(alpha_i)*cos(theta_i)/(pi*kr**2.d0*this%radius**2.d0)
        const_Ja_t_TE = cj*2.d0*sin(alpha_i)/(pi*k1*this%radius)

        n_terms = int(2*k1*this%radius)+1
        allocate(h2(0:n_terms),h2p(0:n_terms))

        call hankarray(n_terms,kr*this%radius,h2,h2p)


        !-------------------------------------------------------
        !! computed current
        H_total = this%Hs(:,0) + this%H(:,0)
        Jt = -H_total((this%M+1):(2*this%M))
        Jz = H_total(1:this%M)
        Jz_phase = 180.d0/pi*atan2(aimag(Jz),real(Jz))
        Jt_phase = 180.d0/pi*atan2(aimag(Jt),real(Jt))
        phi_dev = tpi/dble(this%M)
        allocate(this%I_u(this%M,0:0),this%I_v(this%M,0:0),this%M_u(this%M,0:0),this%M_v(this%M,0:0))
        t=-this%delta_n(1)
        do j =1,this%M

            x_s = this%testing_pt(j,1)
            y_s = this%testing_pt(j,2)
            phi = atan2(y_s,x_s)

            sum_z_TM = 0.d0
            sum_z_TE = 0.d0
            sum_t_TE = 0.d0

            cjj = cj
            do kk=1,n_terms
                sum_z_TM = sum_z_TM + cjj/h2(kk)*cos(kk*(phi-phi_i))
                sum_z_TE = sum_z_TE + kk*cjj/h2p(kk)*sin(kk*(phi-phi_i))
                sum_t_TE = sum_t_TE + cjj/h2p(kk)*cos(kk*(phi-phi_i))
                cjj = cjj*cj
            enddo
            Ja_z_TM = const_Ja_z_TM*(1.d0/h2(0) + 2.d0*sum_z_TM )
            Ja_t_TE = const_Ja_t_TE*(1.d0/h2p(0) + 2.d0*sum_t_TE)
            Ja_z_TE = const_Ja_z_TE*sum_z_TE

            Ja_z = Ja_z_TM + Ja_z_TE
            Ja_t = Ja_t_TE

            normalize = normalize + abs(Ja_z)**2.d0 + abs(Ja_t)**2.d0
            error = error + abs(Jz(j)-Ja_z)**2.d0 + abs(Jt(j)-Ja_t)**2.d0

            Ja_z_phase = 180.d0/pi*atan2(aimag(Ja_z),real(Ja_z))
            Ja_t_phase = 180.d0/pi*atan2(aimag(Ja_t),real(Ja_t))

            this%I_u(j,0) = Ja_t
            this%I_v(j,0) = Ja_z
            this%M_u(j,0) = 0.d0
            this%M_v(j,0) = 0.d0

            t = t + this%delta_n(j)
            if(plot_current /= 0) then
                write(fd1,*) t/lambda,abs(Jz(j)),Jz_phase(j),abs(Ja_z),Ja_z_phase,&
                abs(Jt(j)),Jt_phase(j),abs(Ja_t),Ja_t_phase,&
                0,0,0,0,&
                0,0,0,0


            endif
        enddo




        error = error/ normalize
        if(plot_current /= 0) then
            close(fd1)
        endif
        deallocate(Jz,Jt,H_total,Jz_phase,Jt_phase,h2,h2p)
    end function eval_surface_current_error


    function eval_surface_current_error_IBC(this) result(error)
        type(Scatterer) :: this
        complex*16,allocatable,dimension(:) :: Jz,Jt,Mz,Mt
        integer ::  fd1 = 50
        complex*16,allocatable,dimension(:) :: H_total,E_total
        integer :: j,kk
        real(8),allocatable,dimension(:) :: Jz_phase,Jt_phase, Mz_phase,Mt_phase
        real(8) :: t,kr,kz,phi_dev
        complex*16,allocatable,dimension(:)::h2,h2p
        integer :: n_terms
        real(8) :: error,normalize
        character*30 :: file_name

        complex*16 :: Ja_z_TM,Ja_t_TE,Ja_z,Ja_t !! analytical surface current
        complex*16 :: sum_z_TM,sum_t_TE,cjj,const_Ja_z_TM,const_Ja_t_TE,eta_zz,eta_pp
        real(8) :: x_s,y_s,phi,Ja_z_phase,Ja_t_phase,Ma_z_phase,Ma_t_phase

        if(plot_current /= 0) then
            file_name = 'currents_plotting'//num2str(this%region_ID,1)//'.dat'
            OPEN(fd1, FILE=trim(file_name))
        endif

        error = 0.d0
        normalize = 0.d0

        this%radius = this%a_superquad
        eta_zz = this%eta_vv(1)
        eta_pp = this%eta_uu(1)
        kz = k0*cos(theta_i)
        kr = sqrt(k1**2.d0 - kz**2.d0)
        allocate(this%I_u(this%M,0:0),this%I_v(this%M,0:0),this%M_u(this%M,0:0),this%M_v(this%M,0:0))
        allocate(Jz(this%M),Jt(this%M),Mz(this%M),Mt(this%M),H_total(2*this%M),&
        E_total(2*this%M),Jz_phase(this%M),Jt_phase(this%M),Mz_phase(this%M),Mt_phase(this%M))
        !------------------------------------------------------
        !! Initialization of the exact current distribution
        const_Ja_z_TM = 2.d0*cos(alpha_i)/(pi*kr*this%radius)
        const_Ja_t_TE = 2.d0*cj*sin(alpha_i)/(pi*kr*this%radius)

        n_terms = int(2*k1*this%radius)+1
        allocate(h2(0:n_terms),h2p(0:n_terms))

        call hankarray(n_terms,kr*this%radius,h2,h2p)


        !-------------------------------------------------------
        !! computed current
        H_total = this%Hs(:,0) + this%H(:,0)
        Jt = -H_total((this%M+1):(2*this%M))
        Jz = H_total(1:this%M)

        E_total = this%Es(:,0) + this%E(:,0)
        Mt = E_total((this%M+1):(2*this%M))
        Mz = -E_total(1:this%M)

        Jz_phase = 180.d0/pi*atan2(aimag(Jz),real(Jz))
        Jt_phase = 180.d0/pi*atan2(aimag(Jt),real(Jt))

        Mz_phase = 180.d0/pi*atan2(aimag(Mz),real(Mz))
        Mt_phase = 180.d0/pi*atan2(aimag(Mt),real(Mt))

        phi_dev = tpi/dble(this%M)

        t=-this%delta_n(1)
        do j =1,this%M

            x_s = this%testing_pt(j,1)
            y_s = this%testing_pt(j,2)
            phi = atan2(y_s,x_s)

            sum_z_TM = 0.d0
            sum_t_TE = 0.d0

            cjj = cj
            do kk=1,n_terms
                sum_z_TM = sum_z_TM + cjj/(h2(kk) + cj*eta_zz*h2p(kk) )*cos(kk*(phi-phi_i))
                sum_t_TE = sum_t_TE + cjj/(h2p(kk) - cj*eta_pp*h2(kk))*cos(kk*(phi-phi_i))
                cjj = cjj*cj
            enddo
            Ja_z_TM = const_Ja_z_TM*(1.d0/(h2(0)+ cj*eta_zz*h2p(0)  )  + 2.d0*sum_z_TM )
            Ja_t_TE = const_Ja_t_TE*(1.d0/(h2p(0) - cj*eta_pp*h2(0)) + 2.d0*sum_t_TE)


            Ja_z = Ja_z_TM
            Ja_t = Ja_t_TE

            normalize = normalize + abs(Ja_z)**2.d0 + abs(Ja_t)**2.d0
            error = error + abs(Jz(j)-Ja_z)**2.d0 + abs(Jt(j)-Ja_t)**2.d0

            Ja_z_phase = 180.d0/pi*atan2(aimag(Ja_z),real(Ja_z))
            Ja_t_phase = 180.d0/pi*atan2(aimag(Ja_t),real(Ja_t))

            this%I_u(j,0) = Ja_t
            this%I_v(j,0) = Ja_z
            this%M_u(j,0) = eta1*(this%eta_vv(j)*this%I_v(j,0) + this%eta_vu(j)*this%I_u(j,0))
            this%M_v(j,0) = -eta1*(this%eta_uu(j)*this%I_u(j,0) + this%eta_uv(j)*this%I_v(j,0))

            Ma_z_phase = 180.d0/pi*atan2(aimag(this%M_v(j,0)),real(this%M_v(j,0)))
            Ma_t_phase = 180.d0/pi*atan2(aimag(this%M_u(j,0)),real(this%M_u(j,0)))

            t = t + this%delta_n(j)
            if(plot_current /= 0) then
                write(fd1,*) t/lambda,abs(Jz(j)),Jz_phase(j),abs(Ja_z),Ja_z_phase,&
                abs(Jt(j)),Jt_phase(j),abs(Ja_t),Ja_t_phase,&
                abs(Mz(j)),Mz_phase(j),abs(this%M_v(j,0)),Ma_z_phase,&
                abs(Mt(j)),Mt_phase(j),abs(this%M_u(j,0)),Ma_t_phase


            endif
        enddo
        error = error/ normalize
        if(plot_current /= 0) then
            close(fd1)
        endif
        deallocate(Jz,Jt,H_total,Mz,Mt,E_total,Jz_phase,Jt_phase,Mz_phase,Mt_phase,h2,h2p)
    end function eval_surface_current_error_IBC


    function eval_surface_current_error_PMC(this) result(error)
        type(Scatterer) :: this
        complex*16,allocatable,dimension(:) :: Mz,Mt
        integer ::  fd1 = 50
        complex*16,allocatable,dimension(:) :: E_total
        integer :: j,kk
        real(8),allocatable,dimension(:) :: Mz_phase,Mt_phase
        real(8) :: t,kr,kz,phi_dev
        complex*16,allocatable,dimension(:)::h2,h2p
        integer :: n_terms
        real(8) :: error,normalize
        character*30 :: file_name

        complex*16 :: Ma_t_TM,Ma_z_TE,Ma_z,Ma_t !! analytical surface current
        complex*16 :: sum_t_TM,sum_z_TE,cjj,const_Ma_t_TM,const_Ma_z_TE,const_Ma_z_TM,sum_z_TM,Ma_z_TM
        real(8) :: x_s,y_s,phi,Ma_z_phase,Ma_t_phase

        if(plot_current /= 0) then
            file_name = 'currents_plotting'//num2str(this%region_ID,1)//'.dat'
            OPEN(fd1, FILE=trim(file_name))
        endif

        error = 0.d0
        normalize = 0.d0

        this%radius = this%a_superquad

        kz = k0*cos(theta_i)
        kr = sqrt(k1**2.d0 - kz**2.d0)
        allocate(this%I_u(this%M,0:0),this%I_v(this%M,0:0),this%M_u(this%M,0:0),this%M_v(this%M,0:0))
        allocate(Mz(this%M),Mt(this%M),E_total(2*this%M),Mz_phase(this%M),Mt_phase(this%M))
        !------------------------------------------------------
        !! Initialization of the exact current distribution
        const_Ma_t_TM = -2.d0*cj*eta1*cos(alpha_i)/(pi*k1*this%radius)
        const_Ma_z_TE = 2.d0*sin(alpha_i)*eta1/(pi*kr*this%radius)
        const_Ma_z_TM = 4.d0*cos(alpha_i)*cos(theta_i)*eta1/(pi*kr**2.d0*this%radius**2.d0)


        n_terms = int(2*k1*this%radius)+1
        allocate(h2(0:n_terms),h2p(0:n_terms))

        call hankarray(n_terms,kr*this%radius,h2,h2p)


        !-------------------------------------------------------
        !! computed current
        E_total = this%Es(:,0) + this%E(:,0)
        Mt = E_total((this%M+1):(2*this%M))
        Mz = -E_total(1:this%M)
        Mz_phase = 180.d0/pi*atan2(aimag(Mz),real(Mz))
        Mt_phase = 180.d0/pi*atan2(aimag(Mt),real(Mt))
        phi_dev = tpi/dble(this%M)

        t=-this%delta_n(1)
        do j =1,this%M

            x_s = this%testing_pt(j,1)
            y_s = this%testing_pt(j,2)
            phi = atan2(y_s,x_s)

            sum_z_TE = 0.d0
            sum_t_TM = 0.d0
            sum_z_TM = 0.d0

            cjj = cj
            do kk=1,n_terms
                sum_z_TE = sum_z_TE + cjj/h2(kk)*cos(kk*(phi-phi_i))
                sum_z_TM = sum_z_TM + kk*cjj/h2p(kk)*sin(kk*(phi-phi_i))
                sum_t_TM = sum_t_TM + cjj/h2p(kk)*cos(kk*(phi-phi_i))
                cjj = cjj*cj
            enddo
            Ma_t_TM = const_Ma_t_TM*(1.d0/h2p(0) + 2.d0*sum_t_TM )
            Ma_z_TE = const_Ma_z_TE*(1.d0/h2(0) + 2.d0*sum_z_TE)
            Ma_z_TM = const_Ma_z_TM*sum_z_TM

            Ma_z = Ma_z_TE + Ma_z_TM
            Ma_t = Ma_t_TM

            normalize = normalize + abs(Ma_z)**2.d0 + abs(Ma_t)**2.d0
            error = error + abs(Mz(j)-Ma_z)**2.d0 + abs(Mt(j)-Ma_t)**2.d0

            Ma_z_phase = 180.d0/pi*atan2(aimag(Ma_z),real(Ma_z))
            Ma_t_phase = 180.d0/pi*atan2(aimag(Ma_t),real(Ma_t))

            this%I_u(j,0) = 0.d0
            this%I_v(j,0) = 0.d0
            this%M_u(j,0) = Ma_t
            this%M_v(j,0) = Ma_z

            t = t + this%delta_n(j)
            if(plot_current /= 0) then
                write(fd1,*) t/lambda,0,0,0,0,&
                0,0,0,0,&
                abs(Mz(j)),Mz_phase(j),abs(Ma_z),Ma_z_phase,&
                abs(Mt(j)),Mt_phase(j),abs(Ma_t),Ma_t_phase


            endif
        enddo
        error = error/ normalize
        if(plot_current /= 0) then
            close(fd1)
        endif
        deallocate(Mz,Mt,E_total,Mz_phase,Mt_phase,h2,h2p)
    end function eval_surface_current_error_PMC


    subroutine read_points_file(this,file_points)
        type(scatterer) :: this
        real(8),allocatable,intent(inout) :: file_points(:,:)
        integer :: j,n_points
        integer :: fd = 81
        real(8) :: x_range,y_range,x_max,x_min,y_max,y_min

        OPEN(fd, FILE=trim(scatterer_input_file_names(this%region_ID)))
        read(fd,*) n_points
        read(fd,*)    !! x y of the points

        allocate(file_points(n_points,2))

        do j = 1,n_points
            read(fd,*) file_points(j,1),file_points(j,2)
        enddo

        close(fd)

        x_max = maxval(file_points(:,1))
        x_min = minval(file_points(:,1))

        y_max = maxval(file_points(:,2))
        y_min = minval(file_points(:,2))

        y_range = y_max - y_min
        x_range = x_max - x_min

        if(x_range > y_range) then
            y_range = x_range
        endif

        if(y_range < (lambda/2.d0)) then
            y_range = lambda/2.d0
        endif

        this%y_bound_max = y_max + y_range
        this%y_bound_min = y_min - y_range
        this%x_bound_max = x_max + y_range
        this%x_bound_min = x_min - y_range

    end subroutine read_points_file

    ! -----------------------------------------------------------------------
    ! Subroutine: discretize_scatterer_file
    ! Purpose   : Discretises an arbitrary scatterer defined by a points file.
    !             Distributes SPW testing points per wavelength uniformly
    !             along the boundary arc length. Computes normal and tangential
    !             vectors at each testing point using finite differences.
    ! -----------------------------------------------------------------------
    subroutine discretize_scatterer_file(this,SPW,Scatterers_pointer)
        type(scatterer) :: this
        integer :: SPW
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)
        real(8),allocatable,dimension(:,:) :: file_points1
        type(vector),allocatable,dimension(:) :: t_vector,n_vector
        type(vector) :: n_vector_temp
        integer :: ii,n_points



        call discretize_rect(this,this%segments_points,SPW,Scatterers_pointer,&
        non_uniform_type_read,N_o_read,A_segmentation_read)
        !    call  define_rules(this,this%segments_points)

        allocate(this%testing_pt_status(this%M,3))
        this%testing_pt_status(:,1) = this%Problem_Type
        this%testing_pt_status(:,2) = this%region_ID
        if(this%region_status(1) == 1 ) then
            this%testing_pt_status(:,3) = this%region_status(2)
        else
            this%testing_pt_status(:,3) = 0
        endif


        if(this%Problem_Type == 4 ) then !! only dielectric problems

            n_points = size(this%segments_points,1)
            allocate(file_points1(n_points,3))
            !        do ii = 1,n_points
            !    !            file_points1(ii) = (file_points(ii) - center_estimate)*(1.d0 + 0.25d0) + center_estimate
            !            file_points1(ii,1) = file_points(ii,1) +  (0.25d0)*lambda*sign(1.d0,file_points(ii,1))
            !            file_points1(ii,2) = file_points(ii,2) +  (0.25d0)*lambda*sign(1.d0,file_points(ii,2))
            !        enddo

            allocate(t_vector(n_points),n_vector(n_points))

            do ii = 1,(n_points-1)
                t_vector(ii) = vec(this%segments_points(ii+1,1) - this%segments_points(ii,1),this%segments_points(ii+1,2) -&
                this%segments_points(ii,2),0.d0)
                t_vector(ii) = (1.d0/absolute(t_vector(ii)))*t_vector(ii)
            !        write(*,*) t_vector(ii)
            enddo

            ii = n_points
            t_vector(ii) = vec(this%segments_points(1,1) - this%segments_points(ii,1),this%segments_points(1,2) -&
            this%segments_points(ii,2),0.d0)

            do ii = 1,n_points
                n_vector(ii) = cross_prod(t_vector(ii),z_vec)
            enddo

            do ii = 1,(n_points-1)

                n_vector_temp = n_vector(ii) + n_vector(ii+1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                file_points1(ii,:) = (outside_bound*lambda)*vector2array(n_vector_temp) +&
                (/this%segments_points(ii+1,1),this%segments_points(ii+1,2),0.d0/)

            enddo
            ii = n_points
            n_vector_temp = n_vector(ii) + n_vector(1)
            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            file_points1(ii,:) = (outside_bound*lambda)*vector2array(n_vector_temp) +&
            (/this%segments_points(1,1),this%segments_points(1,2),0.d0/)

            deallocate(t_vector,n_vector)

            !        do ii = 1,size(file_points1,1)
            !            write(*,*)  file_points1(ii,:)
            !        enddo
            !        stop
            !        call set_contour_points_larger(this,0.75d0)
            call  discretize_contour(this,file_points1,Samples_per_wavelength_outside)

            this%N_inside_sources = this%N_con
            allocate(this%Inside_sources_dielectric(size(this%contour_points,1),size(this%contour_points,2)))
            this%Inside_sources_dielectric = this%contour_points

            !        write(*,*) 'Outside Sources Number = ',this%N_inside_sources


            deallocate(this%contour_points)
            deallocate(file_points1)
            this%N_con = 0

            call allocate_corner_points(this,this%segments_points)

        elseif((this%Problem_Type == 1) .or. (this%Problem_Type == 2)) then
            this%N_inside_sources = 0
            call allocate_corner_points(this,this%segments_points)
        else
            this%N_inside_sources = 0
            this%N_con = 0
        endif



        if((Source_Placement == 1) .or.(Source_Placement == 3)) then

            n_points = size(this%segments_points,1)
            allocate(file_points1(n_points,2))


            do ii = 1,n_points
                !            file_points1(ii) = (file_points(ii) - center_estimate)*(1.d0 + 0.25d0) + center_estimate
                file_points1(ii,1) = this%segments_points(ii,1) - (1.d0-bound)*lambda*sign(1.d0,this%segments_points(ii,1))
                file_points1(ii,2) = this%segments_points(ii,2) - (1.d0-bound)*lambda*sign(1.d0,this%segments_points(ii,2))
            !            write(*,*) sign(1.d0,file_points(ii,1)),sign(1.d0,file_points(ii,2))
            enddo

            !        do ii = 1,size(file_points1,1)
            !            write(*,*)  file_points1(ii,:)
            !        enddo

            call  discretize_contour(this,file_points1,Samples_per_wavelength_contour)
            do ii = 1,size(this%contour_points,1)
                write(*,*)  this%contour_points(ii,:)
            enddo
            write(*,*) 'Contour Sources Number = ',this%N_con
            deallocate(file_points1)

        else

        !        this%N_con = 0
        !        call allocate_corner_points(this,file_points)
        endif

    end subroutine discretize_scatterer_file

    subroutine allocate_corner_points(this,file_points)
        type(Scatterer) :: this
        real(8),allocatable,intent(inout) :: file_points(:,:)
        type(vector),allocatable,dimension(:) :: t_vector,n_vector
        type(vector) :: n_vector_temp
        integer :: ii,n_points,spear_type,rear_type
        real(8) :: local_sep,local_sep2

        spear_type = 1
        rear_type = 1



        n_points = size(file_points,1)

!        write(*,*) 'Our Problem Type now is', this%Problem_Type

        allocate(t_vector(n_points),n_vector(n_points))

        do ii = 1,(n_points-1)
            t_vector(ii) = vec(file_points(ii+1,1) - file_points(ii,1),file_points(ii+1,2) - file_points(ii,2),0.d0)
            t_vector(ii) = (1.d0/absolute(t_vector(ii)))*t_vector(ii)
        !        write(*,*) t_vector(ii)
        enddo

        ii = n_points
        t_vector(ii) = vec(file_points(1,1) - file_points(ii,1),file_points(1,2) - file_points(ii,2),0.d0)

        do ii = 1,n_points
            n_vector(ii) = cross_prod(t_vector(ii),z_vec)
        enddo

        if(this%Problem_Type == 4) then

            local_sep = 0.15d0

            this%N_con = 4*n_points
            allocate(this%contour_points(this%N_con,3),this%contour_points_type(this%N_con))
            allocate(this%contour_point_active_region(this%N_con))

            do ii = 1,(n_points-1)

                n_vector_temp = n_vector(ii) + n_vector(ii+1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
                this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
                (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
                this%contour_points_type(ii) = 1
                this%contour_point_active_region(ii) = this%region_status(2)
            enddo
            ii = n_points
            n_vector_temp = n_vector(ii) + n_vector(1)
            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
            (/file_points(1,1),file_points(1,2),0.d0/)
            this%contour_points_type(ii) = 1
            this%contour_point_active_region(ii) = this%region_status(2)

            this%contour_points((n_points+1):(2*n_points),:) = this%contour_points(1:n_points,:)
            this%contour_points_type((n_points+1):(2*n_points)) = 2
            this%contour_point_active_region((n_points+1):(2*n_points)) = this%contour_point_active_region(1:n_points)
            !! for the interior equivalence
            do ii = 1,(n_points-1)

                n_vector_temp = n_vector(ii) + n_vector(ii+1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
                this%contour_points(ii+2*n_points,:) = (local_sep*this%lambda_local)*vector2array(n_vector_temp) +&
                (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
                this%contour_points_type(ii+2*n_points) = 1
                this%contour_point_active_region(ii+2*n_points) = this%region_status(3)
            enddo
            ii = n_points
            n_vector_temp = n_vector(ii) + n_vector(1)
            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            this%contour_points(ii+2*n_points,:) = (local_sep*this%lambda_local)*vector2array(n_vector_temp) +&
            (/file_points(1,1),file_points(1,2),0.d0/)
            this%contour_points_type(ii+2*n_points) = 1
            this%contour_point_active_region(ii+2*n_points) = this%region_status(3)

            this%contour_points((3*n_points+1):(4*n_points),:) = this%contour_points((2*n_points+1):(3*n_points),:)
            this%contour_points_type((3*n_points+1):(4*n_points)) = 2
            this%contour_point_active_region((3*n_points+1):(4*n_points)) = &
            this%contour_point_active_region((2*n_points+1):(3*n_points))



!            this%N_con = 0
        elseif(this%Problem_Type == 1) then

            !        this%N_con = n_points
            !        allocate(this%contour_points(n_points,3),this%contour_points_type(n_points),&
            !        this%contour_points_orientation(this%N_con))

            !        do ii = 1,(n_points-1)
            !
            !            n_vector_temp = n_vector(ii) + n_vector(ii+1)
            !            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            !        !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
            !            this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
            !            (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
            !            this%contour_points_type(ii) = 1
            !            this%contour_points_orientation(ii) = 0
            !        enddo
            !        ii = n_points
            !        n_vector_temp = n_vector(ii) + n_vector(1)
            !        n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            !        this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
            !        (/file_points(1,1),file_points(1,2),0.d0/)
            !        this%contour_points_type(ii) = 1
            !        this%contour_points_orientation(ii) = 0
            local_sep = corner_sep
            local_sep2 = corner_sep2


            if(local_sep2 > 0.d0) then
            this%N_con = 2*n_points
            allocate(this%contour_points(this%N_con,3),this%contour_points_type(this%N_con),&
            this%contour_points_orientation(this%N_con))
            allocate(this%contour_point_active_region(this%N_con))

            do ii = 1,(n_points-1)

                n_vector_temp = n_vector(ii) + n_vector(ii+1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
                this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
                (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
                this%contour_points_type(ii) = 1
                this%contour_points_orientation(ii) = 0

            enddo
            ii = n_points
            n_vector_temp = n_vector(ii) + n_vector(1)
            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
            (/file_points(1,1),file_points(1,2),0.d0/)
            this%contour_points_type(ii) = 1

            this%contour_points_orientation(ii) = 0


            do ii = 1,(n_points-1)

                n_vector_temp = n_vector(ii) + n_vector(ii+1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
                this%contour_points(ii+n_points,:) = (-local_sep2*lambda)*vector2array(n_vector_temp) +&
                (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
                this%contour_points_type(ii+n_points) = 1

                this%contour_points_orientation(ii+n_points) = 1
            enddo
            ii = n_points
            n_vector_temp = n_vector(ii) + n_vector(1)
            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            this%contour_points(ii+n_points,:) = (-local_sep2*lambda)*vector2array(n_vector_temp) +&
            (/file_points(1,1),file_points(1,2),0.d0/)
            this%contour_points_type(ii+n_points) = 1

            this%contour_points_orientation(ii+n_points) = 1

            this%contour_point_active_region(:) = this%region_status(2)

            else

            this%N_con = n_points
            allocate(this%contour_points(this%N_con,3),this%contour_points_type(this%N_con),&
            this%contour_points_orientation(this%N_con))
            allocate(this%contour_point_active_region(this%N_con))

            do ii = 1,(n_points-1)

                n_vector_temp = n_vector(ii) + n_vector(ii+1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
                this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
                (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
                this%contour_points_type(ii) = 1
                this%contour_points_orientation(ii) = 0

            enddo
            ii = n_points
            n_vector_temp = n_vector(ii) + n_vector(1)
            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
            (/file_points(1,1),file_points(1,2),0.d0/)
            this%contour_points_type(ii) = 1

            this%contour_points_orientation(ii) = 0

            this%contour_point_active_region(:) = this%region_status(2)

            endif


        elseif(this%Problem_Type == 2) then

!            write(*,*) 'ben7ot corner points aho'
            local_sep = corner_sep
            local_sep2 = corner_sep2

            if(local_sep2 > 0.d0) then
            this%N_con = 2*n_points
            allocate(this%contour_points(this%N_con,3),this%contour_points_type(this%N_con),&
            this%contour_points_orientation(this%N_con))
            allocate(this%contour_point_active_region(this%N_con))
            do ii = 1,(n_points-1)

                n_vector_temp = n_vector(ii) + n_vector(ii+1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
                this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
                (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
                this%contour_points_type(ii) = 2
                this%contour_points_orientation(ii) = 0

            enddo
            ii = n_points
            n_vector_temp = n_vector(ii) + n_vector(1)
            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
            (/file_points(1,1),file_points(1,2),0.d0/)
            this%contour_points_type(ii) = 2

            this%contour_points_orientation(ii) = 0


            do ii = 1,(n_points-1)

                n_vector_temp = n_vector(ii) + n_vector(ii+1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
                this%contour_points(ii+n_points,:) = (-local_sep2*lambda)*vector2array(n_vector_temp) +&
                (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
                this%contour_points_type(ii+n_points) = 2

                this%contour_points_orientation(ii+n_points) = 1
            enddo
            ii = n_points
            n_vector_temp = n_vector(ii) + n_vector(1)
            n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
            this%contour_points(ii+n_points,:) = (-local_sep2*lambda)*vector2array(n_vector_temp) +&
            (/file_points(1,1),file_points(1,2),0.d0/)
            this%contour_points_type(ii+n_points) = 2

            this%contour_points_orientation(ii+n_points) = 1

            this%contour_point_active_region(:) = this%region_status(2)
            else
                this%N_con = n_points
                allocate(this%contour_points(this%N_con,3),this%contour_points_type(this%N_con),&
                this%contour_points_orientation(this%N_con))
                allocate(this%contour_point_active_region(this%N_con))
                do ii = 1,(n_points-1)

                    n_vector_temp = n_vector(ii) + n_vector(ii+1)
                    n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                    !        write(*,*) t_vector(ii),local_sep*lambda,(-local_sep*lambda)*vector2array(n_vector_temp)
                    this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
                    (/file_points(ii+1,1),file_points(ii+1,2),0.d0/)
                    this%contour_points_type(ii) = 2
                    this%contour_points_orientation(ii) = 0

                enddo
                ii = n_points
                n_vector_temp = n_vector(ii) + n_vector(1)
                n_vector_temp = (1.d0/absolute(n_vector_temp))*n_vector_temp
                this%contour_points(ii,:) = (-local_sep*lambda)*vector2array(n_vector_temp) +&
                (/file_points(1,1),file_points(1,2),0.d0/)
                this%contour_points_type(ii) = 2

                this%contour_points_orientation(ii) = 0

                this%contour_point_active_region(:) = this%region_status(2)
            endif

        else
            this%N_con = 0
        endif



        !    do ii = 1,this%N_con
        !        write(*,*) this%contour_points(ii,:)
        !    enddo

        deallocate(t_vector,n_vector)
    end subroutine allocate_corner_points



    subroutine discretize_superquad(this,SPW,Scatterers_pointer)
        type(scatterer) :: this
        integer :: SPW
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)


        call discretize_rect(this,this%segments_points,SPW,Scatterers_pointer,&
        non_uniform_type_read,N_o_read,A_segmentation_read)

        allocate(this%testing_pt_status(this%M,3))
        this%testing_pt_status(:,1) = this%Problem_Type
        this%testing_pt_status(:,2) = this%region_ID
        if(this%region_status(1) == 1 ) then
            this%testing_pt_status(:,3) = this%region_status(2)
        else
            this%testing_pt_status(:,3) = 0
        endif

        if(this%Problem_Type == 4) then !! Discretization of a larger countour for the inside problem solution
            if(Source_Placement == 1 .or. Source_Placement == 3) then
                call set_contour_points_larger(this,0.75d0)
                this%N_inside_sources = this%N_con
                allocate(this%Inside_sources_dielectric(size(this%contour_points,1),size(this%contour_points,2)))
                this%Inside_sources_dielectric = this%contour_points

                !        write(*,*) 'Outside Sources Number = ',this%N_inside_sources

                deallocate(this%contour_points)
            endif
            this%N_con = 0
        !            write(*,*) this%N_inside_sources
        !! this is just a trick to use a previously existing function to generate the contour points for the inside region problem solution
        else
            this%N_inside_sources = 0
        endif

        if((Source_Placement == 1) .or.(Source_Placement == 3)) then
            call set_contour_points(this)
            write(*,*) 'Contour Sources Number = ',this%N_con
        else
            this%N_con = 0
        endif

    !        do ii = 1,this%M-1
    !            write(*,*) ii,acos(min(dot(this%tang_u(ii),this%tang_u(ii+1)),1.d0))*180/pi
    !        enddo



    end subroutine discretize_superquad

    subroutine set_contour_points_larger(this,bound_o)
        type(Scatterer) :: this
        real(8) :: a_c,b_c,g_c,bound_o
        real(8),allocatable,dimension(:,:) :: contour
        integer :: ii

        a_c = this%a_superquad + (1.d0-bound_o)*lambda
        b_c = this%b_superquad + (1.d0-bound_o)*lambda
        g_c = this%g_superquad

        !        call segment_superquadratic(a_c,b_c,g_c,contour,TOL_segmentation*10.d0)
        call segment_superquadratic(a_c,b_c,g_c,contour,5.d-1*lambda)

        do ii = 1,size(contour,1)

            contour(ii,:) = contour(ii,:) + this%Center(1:2)
        enddo

        call  discretize_contour(this,contour,Samples_per_wavelength_outside)
        deallocate(contour)
    end subroutine set_contour_points_larger

    subroutine set_contour_points(this)
        type(Scatterer) :: this
        real(8) :: a_c,b_c,g_c
        real(8),allocatable,dimension(:,:) :: contour
        integer :: ii

        a_c = this%a_superquad - (1.d0-bound)*lambda
        b_c = this%b_superquad - (1.d0-bound)*lambda
        g_c = this%g_superquad

        !        call segment_superquadratic(a_c,b_c,g_c,contour,TOL_segmentation*10.d0)
        call segment_superquadratic(a_c,b_c,g_c,contour,9.d-2)
        do ii = 1,size(contour,1)

            contour(ii,:) = contour(ii,:) + this%Center(1:2)
        enddo
        call  discretize_contour(this,contour,Samples_per_wavelength_contour)
        deallocate(contour)
        allocate(this%contour_points_type(this%N_con),this%contour_point_active_region(this%N_con),&
        this%contour_points_orientation(this%N_con))


        if(this%Problem_Type == 1) then !! PEC
            this%contour_points_type = 1
        elseif(this%Problem_Type == 2) then !! PMC
            this%contour_points_type = 2
        else
            this%contour_points_type = 1
        endif
        this%contour_point_active_region = 0
        this%contour_points_orientation = 0
    end subroutine set_contour_points

    subroutine discretize_contour(this,contour,SPW)
        type(Scatterer) :: this
        real(8),allocatable,intent(inout)::contour(:,:)
        integer :: SPW

        integer :: ii,kk,mm
        real(8) :: sep,L
        type (vector) :: t_vector,temp_v,temp_v1,temp_v2
        real(8), allocatable, dimension(:) :: segment_length
        integer, allocatable, dimension(:) :: number_of_elements
        integer :: number_points,M_local


        ! initialization of the discretizer vaiables
        number_points = size(contour,1)

        L = lambda/dble(SPW)
        !        write(*,*) 'Samples_per_wavelength =',Samples_per_wavelength,'L=',L
        allocate(segment_length(number_points))
        allocate(number_of_elements(number_points))
        t_vector%v(3) = 0.d0
        do ii = 1,(number_points-1)
            t_vector%v(1) =contour(ii+1,1) - contour(ii,1)
            t_vector%v(2) =contour(ii+1,2) - contour(ii,2)
            !            write(*,*) absolute(t_vector)
            segment_length(ii) = absolute(t_vector)
            number_of_elements(ii)  = ceil(segment_length(ii)/L)
        !            write(*,*) 'number_of_elements(ii)=',number_of_elements(ii),'segment length=',segment_length(ii)
        enddo
        t_vector%v(1) =contour(1,1) - contour(number_points,1)
        t_vector%v(2) =contour(1,2) - contour(number_points,2)

        segment_length(number_points) = absolute(t_vector)
        number_of_elements(number_points)  = ceil(segment_length(number_points)/L)

        M_local = sum(number_of_elements)

        allocate(this%contour_points(M_local,3))
        this%N_con = M_local
        !        write(*,*) '=================================================================='
        !        write(*,*) 'number of contour points, N_con =',this%N_con
        !        write(*,*) '=================================================================='

        mm = 1
        do ii = 1,number_points
            if(ii<number_points) then
                t_vector%v(1) =(contour(ii+1,1) - contour(ii,1))/segment_length(ii)
                t_vector%v(2) =(contour(ii+1,2) - contour(ii,2))/segment_length(ii)
            else
                t_vector%v(1) =(contour(1,1) - contour(ii,1))/segment_length(ii)
                t_vector%v(2) =(contour(1,2) - contour(ii,2))/segment_length(ii)
            endif
            t_vector%v(3) = 0.d0
            sep = segment_length(ii)/number_of_elements(ii)
            !            write(*,*) 'sep =',sep,t_vector

            temp_v1 = vec(contour(ii,1),contour(ii,2),0.d0)
            !            write(*,*) points(ii,:)

            do kk = 1,number_of_elements(ii)
                temp_v = temp_v1 + ((dble(kk)-1.d0)*sep*t_vector)
                temp_v2 = temp_v1 + ((dble(kk)-0.5d0)*sep*t_vector)
                this%contour_points(mm,1) = temp_v%v(1)
                this%contour_points(mm,2) = temp_v%v(2)
                this%contour_points(mm,3) = 0.d0

                !                write(*,*) t_vector%v(1:2),n_vector%v(1:2)
                !                write(14,*) testing_pt(mm,:),norm_v(mm)
                !                write(*,*) mm,this%contour_points(mm,:)
                mm = mm + 1
            enddo
        enddo

        deallocate(segment_length)
        deallocate(number_of_elements)

    end subroutine discretize_contour

    subroutine get_flat_compensate_contour(this)
        type(Scatterer) :: this
        integer :: jj,pp,rep_cnt,N_flat,mm
        type(vector),allocatable,dimension(:) :: norm_array,temp_n_array,tang_array,temp_t_array
        real(8),allocatable,dimension(:,:) :: beg_end_flat_points,beg_end_flat_points_temp
        integer,allocatable,dimension(:) :: N_contours
        real(8) :: L,sep
        real(8),dimension(3) :: point_temp

        pp = 1
        rep_cnt = 0
        do jj = 2,this%M
            !            if(this%norm_v(jj) == this%norm_v(pp)) then
            if(check_approx(this%norm_v(jj),this%norm_v(pp))) then
                !                write(*,*) 'found equal :D'
                rep_cnt = rep_cnt + 1
            else
                if(rep_cnt > Samples_per_wavelength) then
                    if(allocated(beg_end_flat_points)) then
                        N_flat = size(norm_array,1)
                        allocate(beg_end_flat_points_temp(N_flat,6))
                        allocate(temp_n_array(N_flat),temp_t_array(N_flat))
                        beg_end_flat_points_temp = beg_end_flat_points
                        temp_n_array = norm_array
                        temp_t_array = tang_array
                        deallocate(norm_array,tang_array)
                        deallocate(beg_end_flat_points)
                        allocate(beg_end_flat_points(N_flat+1,6))
                        allocate(norm_array(N_flat+1),tang_array(N_flat+1))

                        beg_end_flat_points(1:N_flat,:) = beg_end_flat_points_temp
                        beg_end_flat_points((1+N_flat),:) = (/this%testing_pt(pp,:),this%testing_pt(jj-1,:)/)
                        norm_array(1:N_flat) = temp_n_array
                        !                        write(*,*) pp,rep_cnt,(2*pp+rep_cnt)/2
                        norm_array(N_flat+1) = this%norm_v((2*pp+rep_cnt)/2)
                        tang_array(1:N_flat) = temp_t_array
                        tang_array(N_flat+1) = this%tang_u((2*pp+rep_cnt)/2)
                        deallocate(beg_end_flat_points_temp)
                        deallocate(temp_t_array,temp_n_array)
                    else !! just the first time to allocate arrays
                        allocate(beg_end_flat_points(1,6))
                        allocate(norm_array(1),tang_array(1))

                        beg_end_flat_points(1,1:3) = this%testing_pt(pp,:)
                        beg_end_flat_points(1,4:6) = this%testing_pt(jj-1,:)
                        !                        write(*,*) pp,rep_cnt,(2*pp+rep_cnt)/2
                        norm_array(1) = this%norm_v((2*pp+rep_cnt)/2)
                        tang_array(1) = this%tang_u((2*pp+rep_cnt)/2)
                    endif
                endif

                pp = jj
                rep_cnt = 0
            endif
        enddo
        if(rep_cnt > Samples_per_wavelength) then !! considering the last point case
            if(allocated(beg_end_flat_points)) then
                N_flat = size(norm_array,1)
                allocate(beg_end_flat_points_temp(N_flat,6))
                allocate(temp_n_array(N_flat),temp_t_array(N_flat))
                beg_end_flat_points_temp = beg_end_flat_points
                temp_n_array = norm_array
                temp_t_array = tang_array
                deallocate(norm_array,tang_array)
                deallocate(beg_end_flat_points)
                allocate(beg_end_flat_points(N_flat+1,6))
                allocate(norm_array(N_flat+1),tang_array(N_flat+1))

                beg_end_flat_points(1:N_flat,:) = beg_end_flat_points_temp
                beg_end_flat_points((1+N_flat),:) = (/this%testing_pt(pp,:),this%testing_pt(jj-1,:)/)
                norm_array(1:N_flat) = temp_n_array
                norm_array(N_flat+1) = this%norm_v(pp)
                tang_array(1:N_flat) = temp_t_array
                tang_array(N_flat+1) = this%tang_u(pp)
                deallocate(beg_end_flat_points_temp)
                deallocate(temp_t_array,temp_n_array)
            else !! just the first time to allocate arrays
                allocate(beg_end_flat_points(1,6))
                allocate(norm_array(1),tang_array(1))

                beg_end_flat_points(1,1:3) = this%testing_pt(pp,:)
                beg_end_flat_points(1,4:6) = this%testing_pt(jj-1,:)
                norm_array(1) = this%norm_v(pp)
                tang_array(1) = this%tang_u(pp)
            endif
        endif
        N_flat = size(tang_array,1)
        !         do jj = 1,N_flat
        !            write(*,*) beg_end_flat_points(jj,:)
        !        enddo
        !        write(*,*) '******after changing*******'
        if(check_approx(norm_array(N_flat),norm_array(1)) .and. (beg_end_flat_points(1,1) == beg_end_flat_points(N_flat,4))) then
            ! the first and last points can be combined because the first segment started at its middle
            allocate(beg_end_flat_points_temp(N_flat,6))
            allocate(temp_n_array(N_flat),temp_t_array(N_flat))
            beg_end_flat_points_temp = beg_end_flat_points
            temp_n_array = norm_array
            temp_t_array = tang_array
            deallocate(norm_array,tang_array)
            deallocate(beg_end_flat_points)

            allocate(beg_end_flat_points(N_flat-1,6))
            allocate(norm_array(N_flat-1),tang_array(N_flat-1))

            beg_end_flat_points(1:(N_flat-1),:) = beg_end_flat_points_temp(1:(N_flat-1),:)
            beg_end_flat_points(1,1:2) = beg_end_flat_points_temp(N_flat,1:2)
            norm_array = temp_n_array(1:(N_flat-1))

            tang_array = temp_t_array(1:(N_flat-1))

            deallocate(beg_end_flat_points_temp)
            deallocate(temp_t_array,temp_n_array)
            N_flat = N_flat-1
        endif

        !        do jj = 1,N_flat
        !            write(*,*) beg_end_flat_points(jj,:)
        !        enddo
        allocate(N_contours(N_flat))
        L = lambda/samples_per_wavelength_contour
        do jj=1,N_flat
            N_contours(jj) = ceil(dist(beg_end_flat_points(jj,1:3),beg_end_flat_points(jj,4:6))/L)+1
        enddo
        !        write(*,*) N_contours
        this%N_con = sum(N_contours)
        allocate(this%contour_points(this%N_con,3))
        mm = 1
        do jj = 1,N_flat !! beware that this method can lead to wrong answers in case of having sharp corners
            sep = dist(beg_end_flat_points(jj,1:3),beg_end_flat_points(jj,4:6))/dble(N_contours(jj)-1)
            point_temp = beg_end_flat_points(jj,1:3) - 0.35d0*lambda*norm_array(jj)%v
            !            write(*,*) N_contours(jj),beg_end_flat_points(jj,:)
            do pp = 1,N_contours(jj)
                this%contour_points(mm,:) = point_temp + dble(pp-1)*sep*tang_array(jj)%v
                mm = mm+1
            !                write(*,*) this%contour_points(mm,1:2)
            enddo
        enddo

        deallocate(norm_array,tang_array,beg_end_flat_points)
        deallocate(N_contours)
    end subroutine get_flat_compensate_contour







    subroutine discretize_rect(this,points,SPW,Scatterers_pointer,non_uniform_type,N_o,A)
        type(Scatterer) :: this
        integer :: non_uniform_type,N_o,A
        type(Scatterer),allocatable,intent(inout) :: Scatterers_pointer(:)
        integer :: ii,kk,mm,SPW
        real(8) :: L,sep_sum
        real(8),allocatable,dimension(:) :: nn
        type (vector) :: t_vector,n_vector,z_vector,temp_v,temp_v1,temp_v2
        real(8), allocatable, dimension(:) :: segment_length
        integer, allocatable, dimension(:) :: number_of_elements
        real(8), dimension(:,:) :: points
        real(8) :: lambda_local
        integer :: number_points,M_local,ext_region_id,N_lim
        real(8),dimension(1:2) :: lambda_array
        type segment
            integer :: N,N_added,N_t
            real(8),allocatable,dimension(:) :: dn,dn_m
            real(8) :: Length,E
        end type segment
        type(segment),allocatable,dimension(:) :: segments


        ! initialization of the discretizer vaiables
        number_points = size(points,1)

        ext_region_id = this%region_status(2)
        !    write(*,*) Scatterers_pointer(1)%lambda_local
        !    write(*,*) Scatterers_pointer(2)%lambda_local
        !    write(*,*) Scatterers_pointer(ext_region_id)%lambda_local


        if(ext_region_id == 0) then
            lambda_array(2) = lambda
        else
            lambda_array(2) = Scatterers_pointer(ext_region_id)%lambda_local
        !        write(*,*) 'aho',lambda_array(2),Scatterers_pointer(ext_region_id)%lambda_local
        endif

        if(this%Problem_Type == 4) then !! Dielectric case
            lambda_array(1) = this%lambda_local
        else
            lambda_array(1) = lambda
        endif
        lambda_local = minval(lambda_array)

        !    lambda_local = lambda
        !    lambda_local = sqrt(lambda*minval(lambda_array))

        !    write(*,*) lambda_local,'lambda_array',lambda_array,lambda,ext_region_id
        !    stop

        L = lambda_local/dble(SPW)
        !        write(*,*) 'Samples_per_wavelength =',Samples_per_wavelength,'L=',L
        allocate(segment_length(number_points))
        allocate(number_of_elements(number_points))
        allocate(segments(number_points))
        t_vector%v(3) = 0.d0
        do ii = 1,(number_points-1)
            t_vector%v(1) =points(ii+1,1) - points(ii,1)
            t_vector%v(2) =points(ii+1,2) - points(ii,2)
            !            write(*,*) absolute(t_vector)
            segment_length(ii) = absolute(t_vector)
            number_of_elements(ii)  = ceil(segment_length(ii)/L)

            segments(ii)%Length = segment_length(ii)
            segments(ii)%N = number_of_elements(ii)

            allocate(segments(ii)%dn(1:segments(ii)%N),nn(1:segments(ii)%N))
            do kk = 1,segments(ii)%N
                nn(kk) = dble(kk)
            enddo

            if(non_uniform_type == 0) then
                segments(ii)%dn = segment_length(ii)/segments(ii)%N
            elseif(non_uniform_type == 1) then
                segments(ii)%dn = lambda_local/(SPW*(1.d0 + (dble(A)*nn -&
                dble(A)*dble(segments(ii)%N-1)/2.d0 )**dble(N_o)/dble(segments(ii)%N-1)**dble(N_o)))
            elseif(non_uniform_type == 2) then !! logarithmic
                segments(ii)%dn = lambda_local/(SPW*(1.d0 + exp(log(abs((nn/dble(A)) -&
                (dble(segments(ii)%N)/(dble(A)*2.d0))))/(log(dble(segments(ii)%N)/dble(A)))) ))
            else
                write(*,*) 'ERROR in segmentation: choose a good value for non_uniform_type parameter'
                stop
            endif
            deallocate(nn)
        !            write(*,*) 'number_of_elements(ii)=',number_of_elements(ii),'segment length=',segment_length(ii)
        enddo
        t_vector%v(1) =points(1,1) - points(number_points,1)
        t_vector%v(2) =points(1,2) - points(number_points,2)
        ii = number_points
        segment_length(number_points) = absolute(t_vector)
        number_of_elements(number_points)  = ceil(segment_length(number_points)/L)

        segments(ii)%Length = segment_length(ii)
        segments(ii)%N = number_of_elements(ii)

        allocate(segments(ii)%dn(1:segments(ii)%N),nn(1:segments(ii)%N))
        do kk = 1,segments(ii)%N
            nn(kk) = dble(kk)
        enddo

        if(non_uniform_type == 0) then
            segments(ii)%dn = segment_length(ii)/segments(ii)%N
        elseif(non_uniform_type == 1) then
            segments(ii)%dn = lambda_local/(SPW*(1.d0 + (dble(A)*nn -&
            dble(A)*dble(segments(ii)%N-1)/2.d0 )**dble(N_o)/dble(segments(ii)%N-1)**dble(N_o)))
        elseif(non_uniform_type == 2) then !! logarithmic
            segments(ii)%dn = lambda_local/(SPW*(1.d0 + exp(log(abs((nn/dble(A)) -&
            (dble(segments(ii)%N)/(dble(A)*2.d0))))/(log(dble(segments(ii)%N)/dble(A)))) ))
        endif
        deallocate(nn)

        M_local = 0
        do ii = 1,number_points
            sep_sum = sum(segments(ii)%dn)

            if(sep_sum == segments(ii)%Length) then
                segments(ii)%N_added = 0
                segments(ii)%N_t = segments(ii)%N
                allocate(segments(ii)%dn_m(1:segments(ii)%N_t))
                segments(ii)%dn_m = segments(ii)%dn
                deallocate(segments(ii)%dn)
            else
                segments(ii)%E = segments(ii)%Length - sep_sum
                if(segments(ii)%E < 0.d0) then
                    segments(ii)%N_added = 0
                    segments(ii)%N_t = segments(ii)%N
                    allocate(segments(ii)%dn_m(1:segments(ii)%N_t))
                    segments(ii)%dn_m = segments(ii)%dn + segments(ii)%E/dble(segments(ii)%N)

                else
                    segments(ii)%N_added = ceil(SPW*segments(ii)%E/lambda_local)
                    segments(ii)%N_t = segments(ii)%N + segments(ii)%N_added
                    allocate(segments(ii)%dn_m(1:segments(ii)%N_t))
                    N_lim = ceil(dble(segments(ii)%N)/2.d0)
                    segments(ii)%dn_m(1:N_lim) = segments(ii)%dn(1:N_lim)
                    segments(ii)%dn_m((segments(ii)%N_t-segments(ii)%N + N_lim +1):segments(ii)%N_t) =&
                    segments(ii)%dn((N_lim+1):segments(ii)%N)
                    segments(ii)%dn_m((N_lim+1):(segments(ii)%N_t-segments(ii)%N + N_lim)) = segments(ii)%E/segments(ii)%N_added
                endif


                deallocate(segments(ii)%dn)

            endif
            !        write(*,*) segments(ii)%N,segments(ii)%N_t,segments(ii)%N_added
            M_local = M_local + segments(ii)%N_t
        enddo
        !    stop
        !    do ii = 1,number_points
        !        do kk =1,segments(ii)%N_t
        !            write(*,*) segments(ii)%dn_m(kk)
        !        enddo
        !        write(*,*) '=================='
        !        write(*,*) sum(segments(ii)%dn_m),segments(ii)%Length
        !    enddo
        !
        !    write(*,*) 'M_local',M_local,sum(number_of_elements)
        !    stop

        !    M_local = sum(number_of_elements)

        z_vector%v(1) = 0.d0
        z_vector%v(2) = 0.d0
        z_vector%v(3) = 1.d0
        allocate(this%testing_pt(M_local,3),this%tang_u(M_local),this%tang_v(M_local),this%norm_v(M_local))
        allocate(this%delta_n(M_local))
        allocate(this%testing_pt_MoM(M_local))
        this%M = M_local
        write(*,*) '=================================================================='
        write(*,*) 'number of testing points, M =',this%M
        write(*,*) '=================================================================='

        mm = 1
        do ii = 1,number_points
            if(ii<number_points) then
                t_vector%v(1) =(points(ii+1,1) - points(ii,1))/segment_length(ii)
                t_vector%v(2) =(points(ii+1,2) - points(ii,2))/segment_length(ii)
            else
                t_vector%v(1) =(points(1,1) - points(ii,1))/segment_length(ii)
                t_vector%v(2) =(points(1,2) - points(ii,2))/segment_length(ii)
            endif
            t_vector%v(3) = 0.d0
            n_vector = cross_prod(t_vector,z_vector)
            !        allocate(sep(1:number_of_elements(ii)),nn(1:number_of_elements(ii)))
            !        N = number_of_elements(ii)
            !        do kk = 1,N
            !            nn(kk) = dble(kk)
            !        enddo
            !
            !        !            write(*,*) 'sep =',sep,t_vector
            !        if(non_uniform_type == 0) then
            !            sep = segment_length(ii)/number_of_elements(ii)
            !        elseif(non_uniform_type == 1) then
            !            sep = lambda/(SPW*(1.d0 + (dble(A)*nn - dble(A)*dble(N-1)/2.d0 )**dble(N_o)/dble(N-1)**dble(N_o)))
            !        elseif(non_uniform_type == 2) then !! logarithmic
            !            sep = lambda/(SPW*(1.d0 + exp(log(abs((nn/dble(A)) -&
            !             (dble(N)/(dble(A)*2.d0))))/(log(dble(N)/dble(A)))) ))
            !        else
            !            write(*,*) 'ERROR in segmentation: choose a good value for non_uniform_type parameter'
            !            stop
            !        endif
            !        sep_sum = sum(sep)
            !        sep = sep + (segment_length(ii) - sep_sum)/N

            !        do kk=0,(N-1)
            !        write(*,*) segment_length(ii)/number_of_elements(ii),sep(kk)
            !        enddo
            !        write(*,*) '=================>',segment_length(ii),sum(sep)

            temp_v1 = vec(points(ii,1),points(ii,2),0.d0)
            !            write(*,*) points(ii,:)

            do kk = 0,segments(ii)%N_t-1
                !            temp_v = temp_v1 + ((dble(kk)-1.0d0)*sep*t_vector)


                !            temp_v = temp_v1 + ((dble(kk)+0.5d0)*sep(kk)*t_vector)
                !            temp_v2 = temp_v1 + ((dble(kk)+0.5d0)*sep(kk)*t_vector)
                !            if(kk == 0) then
                !                temp_v = temp_v1+0.5d0*sep(1)*t_vector
                !                temp_v2 = temp_v1 +0.5d0*sep(1)*t_vector
                !            else
                !                temp_v = temp_v1 + (sum(sep(1:kk)) + 0.5d0*sep(kk+1))*t_vector
                !                temp_v2 = temp_v1 + (sum(sep(1:kk)) + 0.5d0*sep(kk+1))*t_vector
                !            endif

                if(kk == 0) then
                    temp_v = temp_v1+0.5d0*segments(ii)%dn_m(1)*t_vector
                    temp_v2 = temp_v1 +0.5d0*segments(ii)%dn_m(1)*t_vector
                else
                    temp_v = temp_v1 + (sum(segments(ii)%dn_m(1:kk)) + 0.5d0*segments(ii)%dn_m(kk+1))*t_vector
                    temp_v2 = temp_v1 + (sum(segments(ii)%dn_m(1:kk)) + 0.5d0*segments(ii)%dn_m(kk+1))*t_vector
                endif

                this%testing_pt(mm,1) = temp_v%v(1)
                this%testing_pt(mm,2) = temp_v%v(2)
                this%testing_pt(mm,3) = 0.d0
                this%delta_n(mm) = segments(ii)%dn_m(kk+1)
                this%testing_pt_MoM(mm) = vec(temp_v2%v(1),temp_v2%v(2),0.d0)

                this%norm_v(mm) = n_vector
                this%tang_u(mm) = t_vector
                this%tang_v(mm) = z_vector
                !                write(*,*) t_vector%v(1:2),n_vector%v(1:2)
                !                write(14,*) testing_pt(mm,:),norm_v(mm)
!                              write(*,*) mm,this%testing_pt(mm,:),this%norm_v(mm)
                mm = mm + 1
            enddo
            deallocate(segments(ii)%dn_m)
        enddo
!        write(*,*) 100,this%testing_pt(100,:),this%norm_v(mm)
        deallocate(segments)
    end subroutine discretize_rect



    subroutine define_rules(this,points)
        type(Scatterer) :: this
        real(8),dimension(:,:) :: points
        integer :: ii,N_actual,kk
        real(8),allocatable,dimension(:,:) :: segments_equ,segments_equ_actual
        real(8) :: x1,y1,x2,y2
        integer , dimension(2) :: curr_equ
        real(8) :: x2_1,x2_2,y2_1,y2_2,x_beg,x_max
        real(8),allocatable,dimension(:,:) :: Str_Rules_temp !! adding points for sources rules (based on structure) temp
        integer :: number_points

        number_points = size(points,1)

        allocate(segments_equ(number_points,4),segments_equ_actual(number_points,4))
        do ii = 1,(number_points-1)
            x1 = points(ii,1)
            x2 = points(ii+1,1)
            y1 = points(ii,2)
            y2 = points(ii+1,2)
            if(x1 == x2) then !! neglect this segment
                segments_equ(ii,1) = 0.d0
                segments_equ(ii,2) = 0.d0
                segments_equ(ii,3) = 0.d0
                segments_equ(ii,4) = 0.d0
            else
                if(x2 < x1) then ! swap
                    call swap_points(x1,y1,x2,y2)
                endif
                segments_equ(ii,1) = x1 ! x_min
                segments_equ(ii,2) = x2 ! x_max
                segments_equ(ii,3) = (y2-y1)/(x2-x1)
                segments_equ(ii,4) = y1-x1*(y2-y1)/(x2-x1)
            endif
        enddo
        ii = number_points
        x1 = points(ii,1)
        x2 = points(1,1)
        y1 = points(ii,2)
        y2 = points(1,2)
        if(x1 == x2) then !! neglect this segment
            segments_equ(ii,1) = 0.d0
            segments_equ(ii,2) = 0.d0
            segments_equ(ii,3) = 0.d0
            segments_equ(ii,4) = 0.d0
        else
            if(x2 < x1) then ! swap
                call swap_points(x1,y1,x2,y2)
            endif
            segments_equ(ii,1) = x1 ! x_min
            segments_equ(ii,2) = x2 ! x_max
            segments_equ(ii,3) = (y2-y1)/(x2-x1)
            segments_equ(ii,4) = y1-x1*(y2-y1)/(x2-x1)
        endif
        kk = 1
        do ii = 1,number_points
            x1 = segments_equ(ii,1)
            x2 = segments_equ(ii,2)
            if(x1 == 0.d0 .and. x2 == 0.d0) then
            !                write(*,*) 'Ha...continue!'
            continue
        else
            segments_equ_actual(kk,:) = segments_equ(ii,:)
            kk = kk + 1
        endif

    enddo
    N_actual = kk - 1

    !        do ii = 1,N_actual
    !            write(*,*) segments_equ_actual(ii,:)
    !        enddo
    allocate(Str_Rules_temp(N_actual-1,6))
    !        write(*,*) segments_equ_actual(1:N_actual,1)
    x1 = minval(segments_equ_actual(1:N_actual,1))
    x_max = maxval(segments_equ_actual(1:N_actual,2))
    kk = 1
    do ii = 1,N_actual
        if(kk > 2) then
            exit
        endif
        x2 = segments_equ_actual(ii,1)
        !            write(*,*) x1,x2
        if(x2 == x1) then
            curr_equ(kk) = ii
            kk = kk +1
        endif

    enddo
    !        write(*,*)  curr_equ
    x_beg = segments_equ_actual(curr_equ(1),1)
    kk = 1
    do ii = 1,(N_actual-1)
        !! find the nearest x2
        x2_1 = segments_equ_actual(curr_equ(1),2)
        x2_2 = segments_equ_actual(curr_equ(2),2)
        x2 = (x_beg+min(x2_1,x2_2))/2.d0 !! get the intermediate point
        !! evaluate the y for each segment
        y2_1 = segments_equ_actual(curr_equ(1),3)*x2 + segments_equ_actual(curr_equ(1),4)
        y2_2 = segments_equ_actual(curr_equ(2),3)*x2 + segments_equ_actual(curr_equ(2),4)
        if(y2_1 > y2_2) then ! set the maximum and minimum segment equations
            Str_Rules_temp(ii,3:4) = segments_equ_actual(curr_equ(1),3:4)
            Str_Rules_temp(ii,5:6) = segments_equ_actual(curr_equ(2),3:4)
        else
            Str_Rules_temp(ii,3:4) = segments_equ_actual(curr_equ(2),3:4)
            Str_Rules_temp(ii,5:6) = segments_equ_actual(curr_equ(1),3:4)
        endif
        Str_Rules_temp(ii,1) = x_beg
        if(x2_2 == x2_1) then !! the two segments are the same from the bottom and top
            Str_Rules_temp(ii,2) = segments_equ_actual(curr_equ(1),2)
            curr_equ(1) = curr_equ(1) - 1
            if(curr_equ(1) < 1) then
                curr_equ(1) = N_actual
            endif
            curr_equ(2) = curr_equ(2) + 1
            if(curr_equ(2) > N_actual) then
                curr_equ(2) = 1
            endif
        elseif(x2_1 < x2_2) then !! segment 1 is shorter than segment 2
            Str_Rules_temp(ii,2) = segments_equ_actual(curr_equ(1),2)
            curr_equ(1) = curr_equ(1) - 1
            if(curr_equ(1) < 1) then
                curr_equ(1) = N_actual
            endif

        else !! segment 2 is shorter
            Str_Rules_temp(ii,2) = segments_equ_actual(curr_equ(2),2)
            curr_equ(2) = curr_equ(2) + 1
            if(curr_equ(2) > N_actual) then
                curr_equ(2) = 1
            endif
        endif
        x_beg = Str_Rules_temp(ii,2)
        if(x_beg == x_max) then
            exit
        endif
        kk = kk + 1
    !            write(*,*) Str_Rules_temp(ii,:)
    enddo
    !    write(*,*) ' '
    this%N_rules = kk
    allocate(this%Str_Rules(this%N_rules,6))
    do ii=1,this%N_rules
        this%Str_Rules(ii,:)=Str_Rules_temp(ii,:)
    !            write(*,*) Str_Rules(ii,:)
    enddo
end subroutine define_rules

subroutine swap_points(x1,y1,x2,y2)
    real(8),intent(inout) :: x1,y1,x2,y2
    real(8) :: xt,yt

    xt = x2
    yt = y2
    x2 = x1
    y2 = y1
    x1 = xt
    y1 = yt
end subroutine swap_points

type (multival) function Fext1(n,x)
    !*****************************************
    ! Kernel of integration
    !*****************************************
    ! Function Arguments
    integer, intent(in) :: n                            ! number of variables and parameters
    real(8), intent(in), dimension(n) :: x     ! variables and parameters
    !! variables
    !! x(1) t'
    !! x(2) i
    !! x(3) j
    !! x(4) k_r
    !! x(5) tang_u(i)%v(1)
    !! x(6) tang_u(i)%v(2)
    !! x(7) tang_u(j)%v(1)
    !! x(8) tang_u(j)%v(2)
    !! x(9) x_t -> testing point x
    !! x(10) y_t -> testing point y
    !! x(11) xp (middle of the segment)
    !! x(12) yp (middle of the segment)
    real(8) :: xp,yp,R,x_t,y_t,k_r
    type(vector) :: R_v,v_t,v_tp,v_n,v_np
    complex*16 :: h0,h1,h2,h1_R
    integer :: i,j

    i = int(x(2))
    j = int(x(3))

    k_r = x(4)
    !    x_t = testing_pt_MoM(i)%v(1)
    !    y_t = testing_pt_MoM(i)%v(2)
    x_t = x(9)
    y_t = x(10)


    v_t = vec(x(5),x(6),0.d0)
    v_tp = vec(x(7),x(8),0.d0)

    !    write(*,*) z_vec

    v_n = cross_prod(v_t,z_vec)
    v_np = cross_prod(v_tp,z_vec)

    !    write(*,*) v_n,norm_v_MoM(i)

    !    xp = testing_pt_MoM(j)%v(1) + x(1)*tang_v_MoM(j)%v(1)
    xp = x(11) + x(1)*x(7)
    !    yp = testing_pt_MoM(j)%v(2) + x(1)*tang_v_MoM(j)%v(2)
    yp = x(12) + x(1)*x(8)

    R_v = vec(x_t-xp,y_t-yp,0.d0)
    R = absolute(R_v)
    call besselh2_01(k_r*R,h0,h1)
    !    h0 = besselh2_0(k_r*R)
    !    h1 = besselh2_1(k_r*R)
    h2 = 2.d0/(k_r*R)*h1-h0
    h1_R = h1/R

    Fext1%f(1) = h0
    !    Fext1%f(2) = h1_R*dot(tang_v_MoM(i),R_v)
    Fext1%f(2) = h1_R*dot(v_t,R_v)
    !    Fext1%f(3) = h1_R*dot(tang_v_MoM(j),R_v)
    Fext1%f(3) = h1_R*dot(v_tp,R_v)
    !    Fext1%f(4) = h1_R*dot(norm_v_MoM(i),R_v)
    Fext1%f(4) = h1_R*dot(v_n,R_v)
    !    Fext1%f(5) = h1_R*dot(norm_v_MoM(j),R_v)
    Fext1%f(5) = h1_R*dot(v_np,R_v)
    Fext1%f(6) = h1_R
    !    Fext1%f(7) = h2/R**2.d0*dot(tang_v_MoM(i),R_v)*dot(tang_v_MoM(j),R_v)
    !    Fext1%f(7) = (-h0*k_r**2.d0+2.d0*k_r*h1/R)*1.d0/R**2.d0*dot(tang_v_MoM(i),R_v)*dot(tang_v_MoM(j),R_v)
    Fext1%f(7) = (-h0*k_r**2.d0+2.d0*k_r*h1_R)*1.d0/R**2.d0*dot(v_t,R_v)*dot(v_tp,R_v)
!    write(*,*) i,j,R,Fext1%f(1)

end function Fext1


type (multival) function Fext2(n,x)
    !*****************************************
    ! Kernel of integration
    !*****************************************
    ! Function Arguments
    integer, intent(in) :: n                            ! number of variables and parameters
    real(8), intent(in), dimension(n) :: x     ! variables and parameters
    !! variables
    !! x(1) t'
    !! x(2) i
    !! x(3) j
    !! x(4) k_r_r
    !! x(5) k_r_i
    !! x(6) tang_u(i)%v(1)
    !! x(7) tang_u(i)%v(2)
    !! x(8) tang_u(j)%v(1)
    !! x(9) tang_u(j)%v(2)
    !! x(10) x_t -> testing point x
    !! x(11) y_t -> testing point y
    !! x(12) xp (middle of the segment)
    !! x(13) yp (middle of the segment)
    real(8) :: xp,yp,R,x_t,y_t
    complex*16 :: k_r
    type(vector) :: R_v,v_t,v_tp,v_n,v_np
    complex*16 :: h0,h1,h2,h1_R
    integer :: i,j

    i = int(x(2))
    j = int(x(3))

    k_r = x(4) + cj*x(5)
    !    x_t = testing_pt_MoM(i)%v(1)
    !    y_t = testing_pt_MoM(i)%v(2)
    x_t = x(10)
    y_t = x(11)


    v_t = vec(x(6),x(7),0.d0)
    v_tp = vec(x(8),x(9),0.d0)

    !    write(*,*) z_vec

    v_n = cross_prod(v_t,z_vec)
    v_np = cross_prod(v_tp,z_vec)

    !    write(*,*) v_n,norm_v_MoM(i)

    !    xp = testing_pt_MoM(j)%v(1) + x(1)*tang_v_MoM(j)%v(1)
    xp = x(12) + x(1)*x(8)
    !    yp = testing_pt_MoM(j)%v(2) + x(1)*tang_v_MoM(j)%v(2)
    yp = x(13) + x(1)*x(9)

    R_v = vec(x_t-xp,y_t-yp,0.d0)
    R = absolute(R_v)
    call besselh2_01(k_r*R,h0,h1)
    !    h0 = besselh2_0(k_r*R)
    !    h1 = besselh2_1(k_r*R)
    h2 = (2.d0,0.d0)/(k_r*R)*h1-h0
    h1_R = h1/R

    Fext2%f(1) = h0
    !    Fext1%f(2) = h1_R*dot(tang_v_MoM(i),R_v)
    Fext2%f(2) = h1_R*dot(v_t,R_v)
    !    Fext1%f(3) = h1_R*dot(tang_v_MoM(j),R_v)
    Fext2%f(3) = h1_R*dot(v_tp,R_v)
    !    Fext1%f(4) = h1_R*dot(norm_v_MoM(i),R_v)
    Fext2%f(4) = h1_R*dot(v_n,R_v)
    !    Fext1%f(5) = h1_R*dot(norm_v_MoM(j),R_v)
    Fext2%f(5) = h1_R*dot(v_np,R_v)
    Fext2%f(6) = h1_R
    !    Fext1%f(7) = h2/R**2.d0*dot(tang_v_MoM(i),R_v)*dot(tang_v_MoM(j),R_v)
    !    Fext1%f(7) = (-h0*k_r**2.d0+2.d0*k_r*h1/R)*1.d0/R**2.d0*dot(tang_v_MoM(i),R_v)*dot(tang_v_MoM(j),R_v)
    Fext2%f(7) = (-h0*k_r*k_r+(2.d0,0.d0)*k_r*h1_R)/R**2.d0*dot(v_t,R_v)*dot(v_tp,R_v)
!    write(*,*) i,j,R,Fext1%f(1)

end function Fext2


subroutine post_processing_colormap_imaging(this,R_probs,R_calc,SPW,delta_R)

    type(Scatterer) :: this
    real(8) :: R_probs,cosphi,sinphi,delta_R,xy_plot,R_calc
    integer :: SPW
    type(vector_c),allocatable,dimension(:) :: E_tot,H_tot
    type(vector_c) :: E_J,E_M
    type(vector_c),allocatable,dimension(:) :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz
    complex*16,allocatable,dimension(:) :: Mx,My,Mz,Jx,Jy,Jz
    real(8),allocatable,dimension(:,:) :: pos
    type(vector_c) :: E_cross_n,n_cross_H
    type(vector) :: norm_vector,co_pol,cross_pol,Ep,En
    complex*16 :: E_J_co,E_J_cross,E_M_co,E_M_cross
    integer :: M,mm,N_R,pp,kk
    real(8) :: x,y,r_pos,sign_factor,r_plot,delta,delta_x,delta_y
    integer :: fd = 55
    integer :: fd2 = 56
    integer :: fd3 = 57
    integer :: Samples_per_wavelength1

    allocate(E_Ix(0:0),E_Iy(0:0),E_Iz(0:0),H_Ix(0:0),H_Iy(0:0),H_Iz(0:0))

    Ep = vec(-cos(theta_i)*cos(phi_i),-cos(theta_i)*sin(phi_i),sin(theta_i))
    En = vec(sin(phi_i),-cos(phi_i),0.d0)
    co_pol = (cos(alpha_i)*Ep+sin(alpha_i)*En)
    cross_pol = (sin(alpha_i)*Ep-cos(alpha_i)*En)

    OPEN(fd, FILE='field_color_map.dat')
    OPEN(fd2,FILE='field_co_cross_color_map.dat')
    OPEN(fd3,FILE='field_exact_co_cross_color_map.dat')

    Samples_per_wavelength1 = 10
    call get_post_processing_boundary_condition(this,R_probs,pos,E_tot,H_tot,Samples_per_wavelength1,.false.)
    M = size(E_tot,1)
    delta = lambda/dble(Samples_per_wavelength1)
    allocate(Mx(M),My(M),Mz(M),Jx(M),Jy(M),Jz(M))
    do mm = 1,M
        cosphi = pos(mm,1)/R_probs
        sinphi = pos(mm,2)/R_probs
        delta_x = delta*abs(sinphi)
        delta_y = delta*abs(cosphi)

        norm_vector = vec(-cosphi,-sinphi,0.d0)
        n_cross_H = delta*cross_prod(norm_vector,H_tot(mm))
        E_cross_n = (-1.d0)*delta*cross_prod(norm_vector,E_tot(mm))

        Mx(mm) = E_cross_n%v(1)
        My(mm) = E_cross_n%v(2)
        Mz(mm) = E_cross_n%v(3)
        Jx(mm) = n_cross_H%v(1)
        Jy(mm) = n_cross_H%v(2)
        Jz(mm) = n_cross_H%v(3)
    enddo
    deallocate(E_tot,H_tot)
    SPW = Samples_per_wavelength
    !        xy_plot = (R_probs-delta_R)*sqrt(2.d0)
    xy_plot = R_calc*sqrt(2.d0)
    N_R = ceil(xy_plot/lambda*SPW)


    do mm= 1,N_R+1

        x = -xy_plot/2.d0 + dble(mm - 1)/N_R*xy_plot

        do pp = 1,N_R+1

            y = -xy_plot/2.d0 + dble(pp - 1)/N_R*xy_plot
            r_plot = sqrt(x**2.d0 + y**2.d0)
            E_J = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
            E_M = E_J
            if(abs(r_plot-R_probs) > delta_R) then
                do kk = 1,M
                    call eval_near_field_Electric_2D(0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,ak,az,ar,eta1,x,y,kk,pos)
                    E_J = E_J + Jx(kk)*E_Ix(0) + Jy(kk)*E_Iy(0) + Jz(kk)*E_Iz(0)
                    E_M = E_M - Mx(kk)*H_Ix(0) - My(kk)*H_Iy(0) - Mz(kk)*H_Iz(0)   !! Using duality

                enddo
            endif
            if(r_plot > R_probs) then !! outside region
                sign_factor = -1.d0
            else !! inside the probes region
                sign_factor = 1.d0
            endif
            E_J_co = sign_factor*dot(co_pol,E_J)
            E_J_cross = sign_factor*dot(cross_pol,E_J)
            E_M_co = sign_factor*dot(co_pol,E_M)
            E_M_cross = sign_factor*dot(cross_pol,E_M)

            !! now the output part
            !                write(*,*) x,y
            write(fd,*) x,y,real(E_J%v(1)),aimag(E_J%v(1)),real(E_J%v(2)),aimag(E_J%v(2)),&
            real(E_J%v(3)),aimag(E_J%v(3)),real(E_M%v(1)),aimag(E_M%v(1)),&
            real(E_M%v(2)),aimag(E_M%v(2)),real(E_M%v(3)),aimag(E_M%v(3))

            write(fd2,*) x,y,real(E_J_co),aimag(E_J_co),real(E_J_cross),aimag(E_J_cross),&
            real(E_M_co),aimag(E_M_co),real(E_M_cross),aimag(E_M_cross)


            E_J = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
            r_pos = ((x-this%center(1))/dble(this%a_superquad))**dble(this%g_superquad)+&
            ((y-this%center(2))/dble(this%b_superquad))**dble(this%g_superquad)
            if(r_pos > 1.d0) then
                !                if(r_plot > R_probs+delta_R) then
                !                    write(*,*) x,y,r_pos
                do kk = 1,this%N
                    if(this%active_region(kk)==0) then
                        if(this%I_stat(kk) == 1) then
                            call eval_near_field_Electric_2D(0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,ak,az,ar,eta1,x,y,kk,this%source_pos)
                        else
                            call eval_near_field_Magnetic_2D(0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,ak,az,ar,eta1,x,y,kk,this%source_pos)
                        endif
                        E_J = E_J + this%I_sources_total(kk,1)*E_Ix(0) + this%I_sources_total(kk,2)*E_Iy(0) +&
                        this%I_sources_total(kk,3)*E_Iz(0)
                    endif
                enddo
            endif
            E_J_co = dot(co_pol,E_J)
            E_J_cross = dot(cross_pol,E_J)
            write(fd3,*) x,y,real(E_J_co),aimag(E_J_co),real(E_J_cross),aimag(E_J_cross)
        enddo
    enddo
    close(fd)
    close(fd2)
    close(fd3)
    deallocate(pos,Jx,Jy,Jz,Mx,My,Mz)
    deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz)
end subroutine post_processing_colormap_imaging



subroutine get_post_processing_boundary_condition(this,R_probs,pos,E_tot,H_tot,Samples_per_wavelength1,include_incidence)
    type(Scatterer) :: this
    real(8) :: R_probs,phi,x,y,prim,const1
    logical :: include_incidence
    real(8),allocatable,intent(inout) :: pos(:,:)
    integer :: pp,M,kk,Samples_per_wavelength1
    type(vector_c),allocatable,intent(inout) :: E_tot(:),H_tot(:)
    type(vector) :: Ep,En,rho_hat,test_pt
    complex*16 :: exponential
    type(vector_c),allocatable,dimension(:)  :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz

    allocate(E_Ix(0:0),E_Iy(0:0),E_Iz(0:0),H_Ix(0:0),H_Iy(0:0),H_Iz(0:0))

    Ep = vec(-cos(theta_i)*cos(phi_i),-cos(theta_i)*sin(phi_i),sin(theta_i))
    En = vec(sin(phi_i),-cos(phi_i),0.d0)
    rho_hat = vec(sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i))

    prim = tpi*R_probs
    M = ceil(prim/lambda*dble(Samples_per_wavelength1))
    allocate(H_tot(M),E_tot(M),pos(M,2))
    const1 = 360.d0/dble(M)/dpr
    E_tot = vec((0.d0,0.d0),(0.d0,0.d0),(0.d0,0.d0))
    H_tot = E_tot
    do pp = 1,M
        phi = const1*dble((pp-1))
        x = R_probs*cos(phi)
        y = R_probs*sin(phi)
        pos(pp,:) = (/x,y/)

        if(include_incidence) then
            test_pt = vec(x,y,0.d0)
            exponential = exp(cj*k1*dot(rho_hat,test_pt))
            E_tot(pp) = eta1*exponential*(cos(alpha_i)*Ep+sin(alpha_i)*En)

            H_tot(pp) = exponential*(sin(alpha_i)*Ep-cos(alpha_i)*En)
        endif

        do kk = 1,this%N
            if(this%active_region(kk)==0) then
                if(this%I_stat(kk) == 1) then
                    call eval_near_field_Electric_2D(0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,ak,az,ar,eta1,x,y,kk,this%source_pos)
                else
                    call eval_near_field_Magnetic_2D(0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,ak,az,ar,eta1,x,y,kk,this%source_pos)
                endif
                E_tot(pp) = E_tot(pp) + this%I_sources_total(kk,0)%v(1)*E_Ix(0) + this%I_sources_total(kk,0)%v(2)*E_Iy(0) +&
                this%I_sources_total(kk,0)%v(1)*E_Iz(0)
                H_tot(pp) = H_tot(pp) + this%I_sources_total(kk,0)%v(1)*H_Ix(0) + this%I_sources_total(kk,0)%v(2)*H_Iy(0) +&
                this%I_sources_total(kk,0)%v(1)*H_Iz(0)

            endif
        enddo

    enddo
    deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz)
end subroutine get_post_processing_boundary_condition

    ! -----------------------------------------------------------------------
    ! Subroutine: Export_RAS_AWE_Chebyshev_Monostatic_RCS
    ! Purpose   : Sweeps the monostatic RCS over nl+nh+1 frequencies in
    !             [fl_f0*f0, fh_f0*f0] using the Chebyshev AWE approximant
    !             for the RAS currents. Writes results to
    !             'Monostatic_RCS_RAS_AWE.dat' (angle vs RCS [dBsm]).
    ! -----------------------------------------------------------------------
subroutine Export_RAS_AWE_Chebyshev_Monostatic_RCS(this,N_order,fl_f0,fh_f0,nl,nh)
    type(Scatterer) :: this
    integer :: N_order
    real(8) :: fl_f0,fh_f0
    real(8) :: f_loc
    integer :: nl,nh,ii,mm,N_total,nn,scr_cnt,cmp_cnt,jj,T_order
    real(8) :: step_l,step_h
    real(8) :: ak0,ar0,az0,const_freq,Kappa
    real(8) , allocatable,dimension(:) :: freq_array,phi_points
    real(8) , allocatable,dimension(:,:) :: T
    complex*16, allocatable, dimension(:) :: E_z_AWE,E_t_AWE,E_z_Taylor,E_t_Taylor
    complex*16, allocatable, dimension(:) :: E_z_Chebyshev,E_t_Chebyshev
    type(vector_c),allocatable,dimension(:) :: I_sources_AWE,I_sources_Taylor,I_sources_Chebyshev
    complex*16 :: numerator,denumerator

    open(25,FILE='Monostatic_RCS_RAS_AWE.dat')
    write(25,*) 'Freq (GHz)','                  Ez                Et'

    N_total = nh+nl+1
    allocate(freq_array(nh+nl+1))
    allocate(phi_points(1),E_z_Taylor(1),E_t_Taylor(1),E_z_AWE(1),E_t_AWE(1))
    allocate(E_z_Chebyshev(1),E_t_Chebyshev(1))
    allocate(I_sources_AWE(this%N),I_sources_Taylor(this%N),I_sources_Chebyshev(this%N))
    step_l = frequency*(1.d0-fl_f0)/dble(nl)
    step_h = frequency*(fh_f0-1.d0)/dble(nh)

    if(Number_of_scatterers > 1) then
        write(*,*) '-----------------------------------------------------'
        write(*,*) 'WARNING, solution bandwidth routine is only implemented for single scatterer problems'
        write(*,*) '-----------------------------------------------------'
    endif

    mm = 1
    do ii = 1,nl+1
        freq_array(mm) = frequency*fl_f0 + dble(ii-1)*step_l
        mm = mm +1
    enddo

    do ii = 1,nh
        freq_array(mm) = frequency + dble(ii)*step_h
        mm = mm +1
    enddo

    T_order = max(this%Pade_L,this%Pade_M)

    phi_points = phi_i
!    write(*,*) 'phi_points =',phi_points
    do mm = 1,(nl+nh+1)
        f_loc = freq_array(mm)

        ak0 = ak*f_loc/Frequency
        ar0 = ar*f_loc/Frequency
        az0 = az*f_loc/Frequency

        I_sources_Taylor = this%I_sources_total(1:this%N,0)
        const_freq = (ak0-ak)*k0
        do nn = 1,N_order
            I_sources_Taylor = I_sources_Taylor + const_freq*this%I_sources_total(1:this%N,nn)
            const_freq = const_freq*(ak0-ak)*k0
        enddo
        call eval_far_field_scatterers_RAS_AWE(this,phi_points,ak0,ar0,az0,E_z_Taylor,E_t_Taylor,I_sources_Taylor)

        do scr_cnt = 1,this%N
            do cmp_cnt = 1,3
                numerator = 0.d0
                denumerator = 1.d0
                do ii = 0,this%Pade_L
                    numerator = numerator + this%Pade_a(scr_cnt,ii)%v(cmp_cnt)*((ak0-ak)*k0)**dble(ii)
                enddo
                do jj = 1,this%Pade_M
                    denumerator = denumerator + this%Pade_b(scr_cnt,jj)%v(cmp_cnt)*((ak0-ak)*k0)**dble(jj)
                enddo
                I_sources_AWE(scr_cnt)%v(cmp_cnt) = numerator/denumerator
            enddo
        enddo


        call eval_far_field_scatterers_RAS_AWE(this,phi_points,ak0,ar0,az0,E_z_AWE,E_t_AWE,I_sources_AWE)


!        Kappa  = Get_Kappa(ak0,fh_r,fl_r)
!
!
!
!
!        call ChebyshevT_D(T_order,Kappa,T)
!
!        do scr_cnt = 1,this%N
!            do cmp_cnt = 1,3
!                numerator = 0.d0
!                denumerator = 1.d0
!                do ii = 0,this%Pade_L
!                    numerator = numerator + this%Chebyshev_c(scr_cnt,ii)%v(cmp_cnt)*T(0,ii)
!                enddo
!                do jj = 1,this%Pade_M
!                    denumerator = denumerator + this%Chebyshev_d(scr_cnt,jj)%v(cmp_cnt)*T(0,jj)
!                enddo
!                I_sources_Chebyshev(scr_cnt)%v(cmp_cnt) = numerator/denumerator
!            enddo
!        enddo
!        deallocate(T)

!        call eval_far_field_scatterers_RAS_AWE(this,phi_points,ak0,ar0,az0,&
!            E_z_Chebyshev,E_t_Chebyshev,I_sources_Chebyshev)

!        write(25,*) f_loc,abs(E_z_Taylor),abs(E_t_Taylor),abs(E_z_AWE),abs(E_t_AWE),&
!        abs(E_z_Chebyshev),abs(E_t_Chebyshev)

        write(25,*) f_loc,abs(E_z_Taylor),abs(E_t_Taylor),abs(E_z_AWE),abs(E_t_AWE)
    enddo

    close(25)
    deallocate(I_sources_AWE,E_z_AWE,E_t_AWE,phi_points,E_z_Taylor,E_t_Taylor)
    deallocate(E_z_Chebyshev,E_t_Chebyshev)
    deallocate(freq_array,I_sources_Taylor,I_sources_Chebyshev)
end subroutine Export_RAS_AWE_Chebyshev_Monostatic_RCS

subroutine Evaluate_RAS_AWE_bandwidth(this,N_order,fl_f0,fh_f0,nl,nh,BW_RAS_Taylor,BW_RAS_AWE)
    type(Scatterer) :: this
    integer :: N_order
    real(8) :: fl_f0,fh_f0
    real(8) :: f_loc
    integer :: nl,nh,ii,mm
    real(8) :: step_l,step_h, BW_RAS_Taylor,BW_RAS_AWE
    real(8) , allocatable,dimension(:) :: freq_array
    real(8) , allocatable,dimension(:) :: error_array_taylor,error_array_pade
    complex*16,allocatable,dimension(:) :: Ev,Hv,Eu,Hu
    real(8) :: norm_E,norm_H

!    open(25,FILE='solution_spectrum_RAS_AWE.dat')
!    write(25,*) 'Freq (GHz)','                  Error Taylor                Error AWE'
    allocate(freq_array(nh+nl+1),error_array_taylor(nl+nh+1),error_array_pade(nl+nh+1))
    step_l = frequency*(1.d0-fl_f0)/dble(nl)
    step_h = frequency*(fh_f0-1.d0)/dble(nh)

    if(Number_of_scatterers > 1) then
        write(*,*) '-----------------------------------------------------'
        write(*,*) 'WARNING, solution bandwidth routine is only implemented for single scatterer problems'
        write(*,*) '-----------------------------------------------------'
    endif

    mm = 1
    do ii = 1,nl+1
        freq_array(mm) = frequency*fl_f0 + dble(ii-1)*step_l
        mm = mm +1
    enddo

    do ii = 1,nh
        freq_array(mm) = frequency + dble(ii)*step_h
        mm = mm +1
    enddo
!    write(*,*) N_order
    call eval_Pade_coefficients(this,N_order,Pade_L)

    allocate(Ev(this%M),Eu(this%M),Hv(this%M),Hu(this%M))
    do mm = 1,(nl+nh+1)
        f_loc = freq_array(mm)

        call eval_excitation_frequency_for_bandwidth_test(this,f_loc,Ev,Eu,Hv,Hu,norm_E,norm_H)

        call eval_error_frequency_RAS_Taylor(this,f_loc,N_order,&
            error_array_taylor(mm),Ev,Eu,Hv,Hu,norm_E,norm_H)

        call eval_error_frequency_RAS_AWE(this,f_loc,Eu,Ev,Hu,Hv,&
        norm_E,norm_H,error_array_pade(mm))
!        write(25,*) f_loc,error_array_taylor(mm),error_array_pade(mm)
    enddo



    BW_RAS_Taylor = Evaluate_Bandwidth(freq_array,error_array_taylor)
    BW_RAS_AWE = Evaluate_Bandwidth(freq_array,error_array_pade)

!    write(*,*) 'RAS-Taylor Bandwidth = ', BW_RAS_Taylor*100
!    write(*,*) 'RAS-AWE Bandwidth = ', BW_RAS_AWE*100

!    close(25)
    deallocate(Ev,Hv,Eu,Hu)
    deallocate(freq_array)
end subroutine Evaluate_RAS_AWE_bandwidth


subroutine Evaluate_MoM_AWE_Chebyshev_bandwidth(this,N_order,fl_f0,fh_f0,nl,nh,&
BW_MoM_Taylor,BW_MoM_AWE,BW_MoM_Chebyshev)
    type(Scatterer) :: this
    integer :: N_order
    real(8) :: fl_f0,fh_f0
    real(8) :: f_loc
    integer :: nl,nh,ii,mm
    real(8) :: step_l,step_h, BW_MoM_Taylor,BW_MoM_AWE,BW_MoM_Chebyshev
    real(8) , allocatable,dimension(:) :: freq_array,error_array_chebyshev
    real(8) , allocatable,dimension(:) :: error_array_taylor,error_array_pade
    complex*16,allocatable,dimension(:) :: Ev,Hv,Eu,Hu
    real(8) :: norm_E,norm_H
    integer :: fd2 = 94

!    open(25,FILE='solution_spectrum_RAS_AWE.dat')
!    write(25,*) 'Freq (GHz)','                  Error Taylor                Error AWE'
    allocate(freq_array(nh+nl+1),error_array_taylor(nl+nh+1),error_array_pade(nl+nh+1))
    allocate(error_array_chebyshev(nl+nh+1))
    step_l = frequency*(1.d0-fl_f0)/dble(nl)
    step_h = frequency*(fh_f0-1.d0)/dble(nh)

    if(Number_of_scatterers > 1) then
        write(*,*) '-----------------------------------------------------'
        write(*,*) 'WARNING, solution bandwidth routine is only implemented for single scatterer problems'
        write(*,*) '-----------------------------------------------------'
    endif

    mm = 1
    do ii = 1,nl+1
        freq_array(mm) = frequency*fl_f0 + dble(ii-1)*step_l
        mm = mm +1
    enddo

    do ii = 1,nh
        freq_array(mm) = frequency + dble(ii)*step_h
        mm = mm +1
    enddo
!    write(*,*) 'evaluating MoM pade coefficients'
!    call eval_Pade_coefficients(this,N_order,Pade_L)
!    write(*,*) N_order,Pade_L
    call eval_MoM_Pade_coefficients(this,N_order,Pade_L)
    call Eval_MoM_chebyshev_expansion_coefficeints(this,N_order,Pade_L,fl_f0,fh_f0)

    open(fd2, FILE = 'MoM_Formulations_BC_error_Wideband_Frequencies.dat')

    allocate(Ev(this%M),Eu(this%M),Hv(this%M),Hu(this%M))
    do mm = 1,(nl+nh+1)
        f_loc = freq_array(mm)

        call eval_excitation_frequency_for_bandwidth_test(this,f_loc,Ev,Eu,Hv,Hu,norm_E,norm_H)

        call eval_error_frequency_MoM_AWE_Chebyshev(this,f_loc,Eu,Ev,Hu,Hv,&
            norm_E,norm_H,error_array_taylor(mm),error_array_pade(mm),error_array_chebyshev(mm))


        write(fd2,*) f_loc/1d9,error_array_taylor(mm),error_array_pade(mm),error_array_chebyshev(mm)
!        write(*,*) f_loc/1d9,error_array_taylor(mm),error_array_pade(mm),error_array_chebyshev(mm)
    enddo

    close(fd2)

    BW_MoM_Taylor = Evaluate_Bandwidth(freq_array,error_array_taylor)
    BW_MoM_AWE = Evaluate_Bandwidth(freq_array,error_array_pade)
    BW_MoM_Chebyshev = Evaluate_Bandwidth(freq_array,error_array_chebyshev)
    write(*,*) '=============================================='
    write(*,*) 'Wideband Assessement: MoM '
    write(*,*) 'Taylor Bandwidth: ',BW_MoM_Taylor,' Percent'
    write(*,*) 'Pade Bandwidth  : ',BW_MoM_AWE,' Percent'
    write(*,*) '=============================================='
!    write(*,*) 'RAS-Taylor Bandwidth = ', BW_MoM_Taylor
!    write(*,*) 'RAS-AWE Bandwidth = ', BW_MoM_AWE

!    close(25)
    deallocate(Ev,Hv,Eu,Hu)
    deallocate(freq_array)
    deallocate(error_array_chebyshev,error_array_pade,error_array_taylor)
end subroutine Evaluate_MoM_AWE_Chebyshev_bandwidth


function Evaluate_Bandwidth(freq_array,error_array) result(BW)
    real(8),dimension(:) :: freq_array,error_array
    real(8):: BW
    real(8) :: freq_start1,freq_end1
    integer :: L_array,ii

    L_array = size(error_array,1)

    do ii = 1,L_array
        if(error_array(ii) <= BW_limit) then
            freq_start1 = freq_array(ii)
            exit
        endif
    enddo

    do ii = L_array,1,-1
        if(error_array(ii) <= BW_limit) then
            freq_end1 = freq_array(ii)
            exit
        endif
    enddo

    BW = (freq_end1 - freq_start1)/frequency*100d0 !(freq_end1 + freq_start1)*2.d0

end function Evaluate_Bandwidth

subroutine get_solution_bandwidth(this,N_order,fl_f0,fh_f0,nl,nh)
    type(Scatterer) ::this
    integer :: N_order
    real(8) :: fl_f0,fh_f0
    real(8) :: f_loc
    integer :: nl,nh,ii,mm
    real(8) :: step_l,step_h
    real(8) , allocatable,dimension(:) :: freq_array,MoM_error_array, err_t
    real(8) , allocatable,dimension(:,:) :: error_array
    integer :: N_f_Calculations
    real(8) :: BW,freq_start1,freq_end1

    open(19,FILE='solution_spectrum.dat')
    write(19,*) 'Freq (GHz)','                  Error'
    allocate(freq_array(nh+nl+1),error_array(nl+nh+1,0:N_order),MoM_error_array(nl+nh+1),err_t(0:N_order))
    step_l = frequency*(1.d0-fl_f0)/dble(nl)
    step_h = frequency*(fh_f0-1.d0)/dble(nh)

    if(Number_of_scatterers > 1) then
        write(*,*) '-----------------------------------------------------'
        write(*,*) 'WARNING, solution bandwidth routine is only implemented for single scatterer problems'
        write(*,*) '-----------------------------------------------------'
    endif

    mm = 1
    do ii = 1,nl+1
        freq_array(mm) = frequency*fl_f0 + dble(ii-1)*step_l
        mm = mm +1
    enddo

    do ii = 1,nh
        freq_array(mm) = frequency + dble(ii)*step_h
        mm = mm +1
    enddo



    do mm = 1,(nl+nh+1)
        f_loc = freq_array(mm)

        call eval_error_frequency(this,f_loc,N_order, err_t,MoM_error_array(mm))
        error_array(mm,:) = err_t
        !        write(19,*) f_loc/1d9,(error_array(mm,ii),ii=0,N_order),MoM_error_array(mm)
        write(19,*) f_loc/1d9,error_array(mm,N_order),MoM_error_array(mm)
    !        write(*,*) f_loc/1d9,error_array(mm)
    enddo
    close(19)

    N_f_Calculations = nl+nh+1
    freq_start1 = 0.d0
    freq_end1 = 0.d0
    do mm = 1,N_f_Calculations
        if(error_array(mm,N_order) <= BW_limit) then
            freq_start1 = freq_array(mm)
            exit
        endif
    enddo

    do mm = N_f_Calculations,1,-1
        if(error_array(mm,N_order) <= BW_limit) then
            freq_end1 = freq_array(mm)
            exit
        endif
    enddo

    BW = (freq_end1 - freq_start1)/frequency*100.d0
    write(*,*) '=============================================='
    write(*,*) 'Wideband Assessement: Taylor Expansion (RAS) '
    write(*,*) 'Bandwidth: ',BW,' Percent'
    write(*,*) '=============================================='



    deallocate(freq_array,error_array,MoM_error_array, err_t)
end subroutine get_solution_bandwidth
!

subroutine eval_excitation_frequency_for_bandwidth_test(this,freq1,Ev,Eu,Hv,Hu,norm_E,norm_H)
    type(Scatterer) :: this
    type(vector) :: Ep,En,rho_hat,test_pt
    integer :: j
    complex*16 :: exponential
    complex*16,allocatable,intent(inout) :: Eu(:),Ev(:),Hv(:),Hu(:)
    type(vector_c),allocatable,dimension(:) :: E_tot,H_tot
    real(8) :: norm_E, norm_H
    real(8) :: ak0,ar0,az0,az1,freq1


    ak0 = ak*freq1/frequency
    ar0 = ar*freq1/frequency
    az0 = az*freq1/frequency

    allocate(E_tot(this%M),H_tot(this%M))

    Ep = vec(-cos(theta_i)*cos(phi_i),-cos(theta_i)*sin(phi_i),sin(theta_i))
    En = vec(sin(phi_i),-cos(phi_i),0.d0)
    rho_hat = vec(sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i))
    norm_E = 0.d0
    norm_H =0.d0

    do j=1,this%M
        test_pt = vec(this%testing_pt(j,1),this%testing_pt(j,2),0.d0)
        exponential = exp(cj*ak0*k0*dot(rho_hat,test_pt))
        !                write(*,*) dot(rho_hat,testing_pt_MoM(j))
        E_tot(j) = eta1*exponential*(cos(alpha_i)*Ep+sin(alpha_i)*En)
        !                write(*,*) Ei(j)%v
        H_tot(j) = exponential*(sin(alpha_i)*Ep-cos(alpha_i)*En)


        Ev(j) = dot(this%tang_v(j),E_tot(j))
        Eu(j) = dot(this%tang_u(j),E_tot(j))
        Hv(j) = dot(this%tang_v(j),H_tot(j))
        Hu(j) = dot(this%tang_u(j),H_tot(j))
        norm_E = norm_E + abs(Eu(j))**2.d0 + abs(Ev(j))**2.d0
        norm_H = norm_H + abs(Hv(j))**2.d0 + abs(Hu(j))**2.d0
    enddo

    deallocate(E_tot,H_tot)
end subroutine eval_excitation_frequency_for_bandwidth_test

subroutine eval_error_frequency_RAS_Taylor(this,freq1,N_order,error_t,Ev,Eu,Hv,Hu,norm_e_loc,norm_h_loc)
    type(Scatterer) :: this
    integer,intent(in) :: N_order
    real(8) :: error_t,norm_e_loc,norm_h_loc,error_v,error_u

    type(vector_c),allocatable,dimension(:) :: E_tot,H_tot
    complex*16,allocatable,intent(inout) :: Eu(:),Ev(:),Hv(:),Hu(:)
    complex*16,allocatable,dimension(:) :: H_v,H_u
    integer:: j,kk,nn

    real(8) :: ak0,ar0,az0,az1,freq1,const_freq
    complex*16 :: ak1,ar1

    complex*16,allocatable,dimension(:) :: BC_total
    type(vector_c),allocatable,dimension(:)  :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz
    type(vector_c),allocatable,dimension(:) :: I_current
    complex*16 :: eta_ex,ar_ex,ak_ex

    ak0 = ak*freq1/frequency
    ar0 = ar*freq1/frequency
    az0 = az*freq1/frequency

    ak1 = this%ak_local*freq1/frequency
    ar1 = this%ar_local*freq1/frequency
    az1 = az0

    allocate(E_tot(this%M),H_tot(this%M))
    allocate(BC_total(2*this%M))
    allocate(I_current(this%N))
    allocate(E_Ix(0:0),E_Iy(0:0),E_Iz(0:0),H_Ix(0:0),H_Iy(0:0),H_Iz(0:0))

    !! computation of the currents vector
    I_current = this%I_sources_total(1:this%N,0)
    const_freq = (ak0-ak)*k0

    do nn = 1,N_order
        I_current = I_current + const_freq*this%I_sources_total(1:this%N,nn)
        const_freq = const_freq*(ak0-ak)*k0
    enddo


    do j = 1,this%M
        E_tot(j)%v(:) = 0.d0
        H_tot(j)%v(:) = 0.d0
        do kk = 1,this%N
            if(this%active_region(kk)==0) then !! inside sources
                eta_ex = eta1
                ar_ex = ar0
                ak_ex = ak0
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                endif

                E_tot(j) = E_tot(j) + I_current(kk)%v(1)*E_Ix(0) + I_current(kk)%v(2)*E_Iy(0) +&
                    I_current(kk)%v(3)*E_Iz(0)
                H_tot(j) = H_tot(j) + I_current(kk)%v(1)*H_Ix(0) + I_current(kk)%v(2)*H_Iy(0) +&
                    I_current(kk)%v(3)*H_Iz(0)

            else !! outside sources
                eta_ex = this%eta_local
                ar_ex = ar1
                ak_ex = ak1
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                endif


                E_tot(j) = E_tot(j) - I_current(kk)%v(1)*E_Ix(0) - I_current(kk)%v(2)*E_Iy(0) -&
                    I_current(kk)%v(3)*E_Iz(0)
                H_tot(j) = H_tot(j) - I_current(kk)%v(1)*H_Ix(0) - I_current(kk)%v(2)*H_Iy(0) -&
                    I_current(kk)%v(3)*H_Iz(0)

            endif

        enddo

    enddo


    !! Evaluate the new boundary condition
    if(this%Problem_Type == 1) then !! PEC
        BC_total(1:this%M) = Eu
        BC_total((1+this%M):(2*this%M)) = Ev
        do j =1,this%M
            BC_total(j) = BC_total(j) + dot(this%tang_u(j),E_tot(j))
            BC_total(j+this%M) = BC_total(j+this%M) + dot(this%tang_v(j),E_tot(j))
        enddo



    elseif(this%Problem_Type == 2) then !! PMC
        BC_total(1:this%M) = Hu
        BC_total((1+this%M):(2*this%M)) = Hv
        do j =1,this%M
            BC_total(j) = BC_total(j)+dot(this%tang_u(j),H_tot(j))
            BC_total(j+this%M) = BC_total(j+this%M)+dot(this%tang_v(j),H_tot(j))
        enddo



    elseif(this%Problem_Type == 3) then !! IBC

        BC_total(1:this%M) = Eu +  eta1*(this%eta_uu*Hv - this%eta_uv*Hu)
        BC_total((1+this%M):(2*this%M)) = Ev +  eta1*(this%eta_vu*Hv - this%eta_vv*Hu)

        allocate(H_u(this%M),H_v(this%M))


        do j = 1,this%M
            H_u(j) = dot(this%tang_u(j),H_tot(j))
            H_v(j) = dot(this%tang_v(j),H_tot(j))
            BC_total(j) = BC_total(j)+dot(this%tang_u(j),E_tot(j))
            BC_total(this%M+j) = BC_total(this%M+j)+dot(this%tang_v(j),E_tot(j))
        enddo
        H_u = eta1*H_u
        H_v = eta1*H_v
        do j = 1,this%M
            BC_total(j) = BC_total(j) +  this%eta_uu(j)*H_v(j) - this%eta_uv(j)*H_u(j) !! E_u
            BC_total(this%M+j) = BC_total(this%M+j) +  this%eta_vu(j)*H_v(j) - this%eta_vv(j)*H_u(j) !! E_v
        enddo

        deallocate(H_u,H_v)


    elseif(this%Problem_Type == 4) then !! Dielectric
        BC_total(1:this%M) = Eu
        BC_total((1+this%M):(2*this%M)) = Ev
        do j =1,this%M
            BC_total(j) = BC_total(j) + dot(this%tang_u(j),E_tot(j))
            BC_total(j+this%M) = BC_total(j+this%M) + dot(this%tang_v(j),E_tot(j))
        enddo


    else
        write(*,*) 'ERROR: wrong problem type, it should be in the range of [1,4]'
        STOP
    endif

    error_u = 0.d0
    error_v = 0.d0
    do j=1,this%M
        error_u = error_u + abs(BC_total(j))**2.d0
        error_v = error_v + abs(BC_total(j+this%M))**2.d0

    enddo
    if(this%Problem_Type == 2) then !! PMC
        error_t = (error_v + error_u)/norm_H_loc
    else
        error_t = (error_v + error_u)/norm_E_loc
    endif



    deallocate(BC_total)
    deallocate(I_current)
    deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz)
    deallocate(E_tot,H_tot)
end subroutine eval_error_frequency_RAS_Taylor

    ! -----------------------------------------------------------------------
    ! Subroutine: Export_MoM_AWE_Chebyshev_Monostatic_RCS
    ! Purpose   : Same as Export_RAS_AWE_Chebyshev_Monostatic_RCS but uses
    !             the MoM/Chebyshev AWE surface current expansion.
    !             Writes to 'Monostatic_RCS_MoM_AWE_Chebyshev.dat'.
    ! -----------------------------------------------------------------------
subroutine Export_MoM_AWE_Chebyshev_Monostatic_RCS(this,N_order,fl_f0,fh_f0,nl,nh)
    complex*16,allocatable,dimension(:) :: MoM_I_z,MoM_I_t,MoM_M_z,MoM_M_t
    complex*16,allocatable,dimension(:) :: MoM_I_z_AWE,MoM_I_t_AWE,MoM_M_z_AWE,MoM_M_t_AWE
!    complex*16,allocatable,dimension(:) :: MoM_I_z_Chebyshev,MoM_I_t_Chebyshev
!    complex*16,allocatable,dimension(:) :: MoM_M_z_Chebyshev,MoM_M_t_Chebyshev
    type(Scatterer) :: this
    integer :: ii,jj,nn,N_order,scr_cnt,i,T_order
    real(8),allocatable,dimension(:,:) :: T
    complex*16 :: numerator,denumerator,eta_ex,ar_ex,ak_ex
    real(8) :: Kappa,step_l,step_h,const_freq
    real(8) :: f_loc,fl_f0,fh_f0,ak0,ar0,az0,ak1,ar1,az1
    integer :: nl,nh,mm
    real(8) , allocatable,dimension(:) :: freq_array,phi_points
    complex*16, allocatable, dimension(:) :: E_z_AWE,E_t_AWE,E_z_Taylor,E_t_Taylor
!    complex*16, allocatable, dimension(:) :: E_z_Chebyshev,E_t_Chebyshev

    allocate(MoM_I_z(this%M),MoM_I_t(this%M),MoM_M_z(this%M),MoM_M_t(this%M))
    allocate(MoM_I_z_AWE(this%M),MoM_I_t_AWE(this%M),MoM_M_z_AWE(this%M),MoM_M_t_AWE(this%M))
!    allocate(MoM_I_z_Chebyshev(this%M),MoM_I_t_Chebyshev(this%M),&
!    MoM_M_z_Chebyshev(this%M),MoM_M_t_Chebyshev(this%M))
    allocate(freq_array(nh+nl+1))
    allocate(phi_points(1),E_z_Taylor(1),E_t_Taylor(1),E_z_AWE(1),E_t_AWE(1))
!    allocate(E_z_Chebyshev(1),E_t_Chebyshev(1))

    open(25,FILE='Monostatic_RCS_MoM_AWE_Chebyshev.dat')
    write(25,*) 'Freq (GHz)','                  Ez                Et'

    step_l = frequency*(1.d0-fl_f0)/dble(nl)
    step_h = frequency*(fh_f0-1.d0)/dble(nh)

    if(Number_of_scatterers > 1) then
        write(*,*) '-----------------------------------------------------'
        write(*,*) 'WARNING, solution bandwidth routine is only implemented for single scatterer problems'
        write(*,*) '-----------------------------------------------------'
    endif

    mm = 1
    do ii = 1,nl+1
        freq_array(mm) = frequency*fl_f0 + dble(ii-1)*step_l
        mm = mm +1
    enddo

    do ii = 1,nh
        freq_array(mm) = frequency + dble(ii)*step_h
        mm = mm +1
    enddo

    T_order = max(this%Pade_L,this%Pade_M)

    phi_points = phi_i
!    write(*,*) 'phi_points =',phi_points
    do mm = 1,(nl+nh+1)
        f_loc = freq_array(mm)

        ak0 = ak*f_loc/frequency
        ar0 = ar
        az0 = az

        ak1 = this%ak_local*f_loc/frequency
        ar1 = this%ar_local
        az1 = az0

        do scr_cnt = 1,this%M

            numerator = 0.d0
            denumerator = 1.d0
            do ii = 0,this%Pade_L
                numerator = numerator + this%MoM_Pade_a_I_u(scr_cnt,ii)*((ak0-ak)*k0)**dble(ii)
            enddo
            do jj = 1,this%Pade_M
                denumerator = denumerator + this%MoM_Pade_b_I_u(scr_cnt,jj)*((ak0-ak)*k0)**dble(jj)
            enddo
            MoM_I_t_AWE(scr_cnt) = numerator/denumerator

            numerator = 0.d0
            denumerator = 1.d0
            do ii = 0,this%Pade_L
                numerator = numerator + this%MoM_Pade_a_M_u(scr_cnt,ii)*((ak0-ak)*k0)**dble(ii)
            enddo
            do jj = 1,this%Pade_M
                denumerator = denumerator + this%MoM_Pade_b_M_u(scr_cnt,jj)*((ak0-ak)*k0)**dble(jj)
            enddo
            MoM_M_t_AWE(scr_cnt) = numerator/denumerator

            numerator = 0.d0
            denumerator = 1.d0
            do ii = 0,this%Pade_L
                numerator = numerator + this%MoM_Pade_a_M_v(scr_cnt,ii)*((ak0-ak)*k0)**dble(ii)
            enddo
            do jj = 1,this%Pade_M
                denumerator = denumerator + this%MoM_Pade_b_M_v(scr_cnt,jj)*((ak0-ak)*k0)**dble(jj)
            enddo
            MoM_M_z_AWE(scr_cnt) = numerator/denumerator

            numerator = 0.d0
            denumerator = 1.d0
            do ii = 0,this%Pade_L
                numerator = numerator + this%MoM_Pade_a_I_v(scr_cnt,ii)*((ak0-ak)*k0)**dble(ii)
            enddo
            do jj = 1,this%Pade_M
                denumerator = denumerator + this%MoM_Pade_b_I_v(scr_cnt,jj)*((ak0-ak)*k0)**dble(jj)
            enddo
            MoM_I_z_AWE(scr_cnt) = numerator/denumerator

        enddo
!        Kappa  = Get_Kappa(ak0,fh_r,fl_r)
!
!!        write(*,*) Kappa
!
!        call ChebyshevT_D(T_order,Kappa,T)
!
!         do scr_cnt = 1,this%M
!            numerator = 0.d0
!            denumerator = 1.d0
!            do ii = 0,this%Pade_L
!                numerator = numerator + this%MoM_Chebyshev_c_I_u(scr_cnt,ii)*T(0,ii)
!            enddo
!            do jj = 1,this%Pade_M
!                denumerator = denumerator + this%MoM_Chebyshev_d_I_u(scr_cnt,jj)*T(0,jj)
!            enddo
!            MoM_I_t_Chebyshev(scr_cnt) = numerator/denumerator
!
!         enddo
!
!         do scr_cnt = 1,this%M
!            numerator = 0.d0
!            denumerator = 1.d0
!            do ii = 0,this%Pade_L
!                numerator = numerator + this%MoM_Chebyshev_c_I_v(scr_cnt,ii)*T(0,ii)
!            enddo
!            do jj = 1,this%Pade_M
!                denumerator = denumerator + this%MoM_Chebyshev_d_I_v(scr_cnt,jj)*T(0,jj)
!            enddo
!            MoM_I_z_Chebyshev(scr_cnt) = numerator/denumerator
!
!         enddo
!
!         do scr_cnt = 1,this%M
!            numerator = 0.d0
!            denumerator = 1.d0
!            do ii = 0,this%Pade_L
!                numerator = numerator + this%MoM_Chebyshev_c_M_v(scr_cnt,ii)*T(0,ii)
!            enddo
!            do jj = 1,this%Pade_M
!                denumerator = denumerator + this%MoM_Chebyshev_d_M_v(scr_cnt,jj)*T(0,jj)
!            enddo
!            MoM_M_z_Chebyshev(scr_cnt) = numerator/denumerator
!
!         enddo
!
!         do scr_cnt = 1,this%M
!            numerator = 0.d0
!            denumerator = 1.d0
!            do ii = 0,this%Pade_L
!                numerator = numerator + this%MoM_Chebyshev_c_M_u(scr_cnt,ii)*T(0,ii)
!            enddo
!            do jj = 1,this%Pade_M
!                denumerator = denumerator + this%MoM_Chebyshev_d_M_u(scr_cnt,jj)*T(0,jj)
!            enddo
!            MoM_M_t_Chebyshev(scr_cnt) = numerator/denumerator
!
!         enddo

!        deallocate(T)

        MoM_I_z = this%I_v(:,0)
        MoM_I_t = this%I_u(:,0)
        MoM_M_z = this%M_v(:,0)
        MoM_M_t = this%M_u(:,0)
        const_freq = (ak0-ak)*k0

        do nn = 1,N_order
            MoM_I_z = MoM_I_z + const_freq*this%I_v(:,nn)
            MoM_I_t = MoM_I_t + const_freq*this%I_u(:,nn)
            MoM_M_z = MoM_M_z + const_freq*this%M_v(:,nn)
            MoM_M_t = MoM_M_t + const_freq*this%M_u(:,nn)

            const_freq = const_freq*(ak0-ak)*k0
        enddo

        call eval_far_fields_scatterer_MoM_Wideband(this,phi_points,ak0,ar0,az0,&
        E_z_Taylor,E_t_Taylor,MoM_I_z,MoM_I_t,MoM_M_z,MoM_M_t)

        call eval_far_fields_scatterer_MoM_Wideband(this,phi_points,ak0,ar0,az0,&
        E_z_AWE,E_t_AWE,MoM_I_z_AWE,MoM_I_t_AWE,MoM_M_z_AWE,MoM_M_t_AWE)

!        call eval_far_fields_scatterer_MoM_Wideband(this,phi_points,ak0,ar0,az0,&
!        E_z_Chebyshev,E_t_Chebyshev,MoM_I_z_Chebyshev,MoM_I_t_Chebyshev,&
!        MoM_M_z_Chebyshev,MoM_M_t_Chebyshev)


!        write(25,*) f_loc,abs(E_z_Taylor),abs(E_t_Taylor),abs(E_z_AWE),abs(E_t_AWE),&
!        abs(E_z_Chebyshev),abs(E_t_Chebyshev)

        write(25,*) f_loc,abs(E_z_Taylor),abs(E_t_Taylor),abs(E_z_AWE),abs(E_t_AWE)
!        write(*,*) f_loc,abs(E_z_Taylor),abs(E_t_Taylor),abs(E_z_AWE),abs(E_t_AWE)

    enddo


    close(25)
    deallocate(MoM_I_z_AWE,MoM_I_t_AWE,MoM_M_z_AWE,MoM_M_t_AWE)
    deallocate(MoM_I_z,MoM_I_t,MoM_M_z,MoM_M_t)
!    deallocate(MoM_I_z_Chebyshev,MoM_I_t_Chebyshev,MoM_M_z_Chebyshev,MoM_M_t_Chebyshev)
    deallocate(freq_array,phi_points)
    deallocate(E_z_AWE,E_t_AWE,E_z_Taylor,E_t_Taylor)!,E_z_Chebyshev,E_t_Chebyshev)
end subroutine Export_MoM_AWE_Chebyshev_Monostatic_RCS

subroutine eval_error_frequency_MoM_AWE_Chebyshev(this,freq1,Eu,Ev,Hu,Hv,norm_E_loc,&
    norm_H_loc,MoM_error_t,MoM_error_AWE,MoM_error_Chebyshev)
    type(Scatterer) :: this
    real(8) :: MoM_error_t,MoM_error_Chebyshev
    complex*16,allocatable,dimension(:) :: H_u,H_v
    real(8) :: ak0,ar0,az0,az1,freq1,t_l
    complex*16 :: ak1,ar1
    real(8) :: norm_E_loc, norm_H_loc,const_freq,MoM_error_AWE
    integer :: ii,jj,nn,N_order,scr_cnt,i,T_order
    real(8),allocatable,dimension(:,:) :: T
    complex*16 :: numerator,denumerator,eta_ex,ar_ex,ak_ex
    complex*16,allocatable,intent(inout) :: Eu(:),Ev(:),Hu(:),Hv(:)
    complex*16,allocatable,dimension(:) :: MoM_I_z,MoM_I_t,MoM_M_z,MoM_M_t
    complex*16,allocatable,dimension(:) :: MoM_I_z_AWE,MoM_I_t_AWE,MoM_M_z_AWE,MoM_M_t_AWE
    complex*16,allocatable,dimension(:) :: MoM_I_z_Chebyshev,MoM_I_t_Chebyshev
    complex*16,allocatable,dimension(:) :: MoM_M_z_Chebyshev,MoM_M_t_Chebyshev
    complex*16,allocatable,dimension(:,:) :: Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt
    complex*16,allocatable,dimension(:)  :: MoM_Ez,MoM_Et,MoM_Hz,MoM_Ht
    complex*16,allocatable,dimension(:)  :: MoMs_Ez,MoMs_Et,MoMs_Hz,MoMs_Ht
    complex*16,allocatable,dimension(:)  :: MoMs_AWE_Ez,MoMs_AWE_Et,MoMs_AWE_Hz,MoMs_AWE_Ht
    complex*16,allocatable,dimension(:)  :: MoMs_Chebyshev_Ez,MoMs_Chebyshev_Et
    complex*16,allocatable,dimension(:)  :: MoMs_Chebyshev_Hz,MoMs_Chebyshev_Ht
    character*100 :: file_name
    integer :: fd1=91
    real(8) :: Kappa






    N_order = this%Pade_L+This%Pade_M
!    write(*,*) N_order

    ak0 = ak*freq1/frequency
!    ar0 = ar*freq1/frequency
!    az0 = az*freq1/frequency
    ar0 = ar
    az0 = az

    ak1 = this%ak_local*freq1/frequency
!    ar1 = this%ar_local*freq1/frequency
!    az1 = az0
    ar1 = this%ar_local
    az1 = az0

!    write(*,*) 'evaluating error'

    allocate(MoM_Ez(this%M),MoM_Et(this%M),MoM_Hz(this%M),MoM_Ht(this%M))
    allocate(MoMs_Ez(this%M),MoMs_Et(this%M),MoMs_Hz(this%M),MoMs_Ht(this%M))
    allocate(MoMs_AWE_Ez(this%M),MoMs_AWE_Et(this%M),MoMs_AWE_Hz(this%M),MoMs_AWE_Ht(this%M))
    allocate(MoMs_Chebyshev_Ez(this%M),MoMs_Chebyshev_Et(this%M),MoMs_Chebyshev_Hz(this%M),MoMs_Chebyshev_Ht(this%M))
    allocate(MoM_I_z(this%M),MoM_I_t(this%M),MoM_M_z(this%M),MoM_M_t(this%M))
    allocate(MoM_I_z_AWE(this%M),MoM_I_t_AWE(this%M),MoM_M_z_AWE(this%M),MoM_M_t_AWE(this%M))

    allocate(MoM_I_z_Chebyshev(this%M),MoM_I_t_Chebyshev(this%M),&
    MoM_M_z_Chebyshev(this%M),MoM_M_t_Chebyshev(this%M))

!    write(*,*)

    do scr_cnt = 1,this%M

        numerator = 0.d0
        denumerator = 1.d0
        do ii = 0,this%Pade_L
            numerator = numerator + this%MoM_Pade_a_I_u(scr_cnt,ii)*((ak0-ak)*k0)**dble(ii)
        enddo
        do jj = 1,this%Pade_M
            denumerator = denumerator + this%MoM_Pade_b_I_u(scr_cnt,jj)*((ak0-ak)*k0)**dble(jj)
        enddo
        MoM_I_t_AWE(scr_cnt) = numerator/denumerator

        numerator = 0.d0
        denumerator = 1.d0
        do ii = 0,this%Pade_L
            numerator = numerator + this%MoM_Pade_a_M_u(scr_cnt,ii)*((ak0-ak)*k0)**dble(ii)
        enddo
        do jj = 1,this%Pade_M
            denumerator = denumerator + this%MoM_Pade_b_M_u(scr_cnt,jj)*((ak0-ak)*k0)**dble(jj)
        enddo
        MoM_M_t_AWE(scr_cnt) = numerator/denumerator

        numerator = 0.d0
        denumerator = 1.d0
        do ii = 0,this%Pade_L
            numerator = numerator + this%MoM_Pade_a_M_v(scr_cnt,ii)*((ak0-ak)*k0)**dble(ii)
        enddo
        do jj = 1,this%Pade_M
            denumerator = denumerator + this%MoM_Pade_b_M_v(scr_cnt,jj)*((ak0-ak)*k0)**dble(jj)
        enddo
        MoM_M_z_AWE(scr_cnt) = numerator/denumerator

        numerator = 0.d0
        denumerator = 1.d0
        do ii = 0,this%Pade_L
            numerator = numerator + this%MoM_Pade_a_I_v(scr_cnt,ii)*((ak0-ak)*k0)**dble(ii)
        enddo
        do jj = 1,this%Pade_M
            denumerator = denumerator + this%MoM_Pade_b_I_v(scr_cnt,jj)*((ak0-ak)*k0)**dble(jj)
        enddo
        MoM_I_z_AWE(scr_cnt) = numerator/denumerator

    enddo


    Kappa  = Get_Kappa(ak0,fh_r,fl_r)




    T_order = max(this%Pade_L,this%Pade_M)

    call ChebyshevT_D(T_order,Kappa,T)

!     write(*,*) Kappa
!    write(*,*) T(0,:)

     do scr_cnt = 1,this%M
        numerator = 0.d0
        denumerator = 1.d0
        do ii = 0,this%Pade_L
            numerator = numerator + this%MoM_Chebyshev_c_I_u(scr_cnt,ii)*T(0,ii)
        enddo
        do jj = 1,this%Pade_M
            denumerator = denumerator + this%MoM_Chebyshev_d_I_u(scr_cnt,jj)*T(0,jj)
        enddo
        MoM_I_t_Chebyshev(scr_cnt) = numerator/denumerator

     enddo

     do scr_cnt = 1,this%M
        numerator = 0.d0
        denumerator = 1.d0
        do ii = 0,this%Pade_L
            numerator = numerator + this%MoM_Chebyshev_c_I_v(scr_cnt,ii)*T(0,ii)
        enddo
        do jj = 1,this%Pade_M
            denumerator = denumerator + this%MoM_Chebyshev_d_I_v(scr_cnt,jj)*T(0,jj)
        enddo
        MoM_I_z_Chebyshev(scr_cnt) = numerator/denumerator

     enddo

     do scr_cnt = 1,this%M
        numerator = 0.d0
        denumerator = 1.d0
        do ii = 0,this%Pade_L
            numerator = numerator + this%MoM_Chebyshev_c_M_v(scr_cnt,ii)*T(0,ii)
        enddo
        do jj = 1,this%Pade_M
            denumerator = denumerator + this%MoM_Chebyshev_d_M_v(scr_cnt,jj)*T(0,jj)
        enddo
        MoM_M_z_Chebyshev(scr_cnt) = numerator/denumerator

     enddo

     do scr_cnt = 1,this%M
        numerator = 0.d0
        denumerator = 1.d0
        do ii = 0,this%Pade_L
            numerator = numerator + this%MoM_Chebyshev_c_M_u(scr_cnt,ii)*T(0,ii)
        enddo
        do jj = 1,this%Pade_M
            denumerator = denumerator + this%MoM_Chebyshev_d_M_u(scr_cnt,jj)*T(0,jj)
        enddo
        MoM_M_t_Chebyshev(scr_cnt) = numerator/denumerator

     enddo

    deallocate(T)

    MoM_Ez = Ev
    MoM_Et = Eu

    MoM_Hz = Hv
    MoM_Ht = Hu


    MoM_I_z = this%I_v(:,0)
    MoM_I_t = this%I_u(:,0)
    MoM_M_z = this%M_v(:,0)
    MoM_M_t = this%M_u(:,0)
    const_freq = (ak0-ak)*k0

    do nn = 1,N_order
        MoM_I_z = MoM_I_z + const_freq*this%I_v(:,nn)
        !        write(*,*) 'MoM_I_z',sum(MoM_I_z),sum(this%I_v(:,nn))
        !        write(*,*) MoM_I_z(10)
        MoM_I_t = MoM_I_t + const_freq*this%I_u(:,nn)
        MoM_M_z = MoM_M_z + const_freq*this%M_v(:,nn)
        MoM_M_t = MoM_M_t + const_freq*this%M_u(:,nn)

        const_freq = const_freq*(ak0-ak)*k0
    !        write(*,*) const_freq
    enddo

     if(plot_current /= 0) then
    file_name = 'currents_plotting_Taylor_AWE_Frequency_'//num2str(int(freq1/1d9),2)//'.dat'
!    write(*,*) int(freq1/1d9),file_name
    OPEN(fd1, FILE=trim(file_name))
    endif
    t_l = -this%delta_n(1)/2.d0
    do i=1,this%M
        t_l = t_l+this%delta_n(i)

        if(plot_current /= 0) then
            write(fd1,*) t_l,abs(MoM_I_z(i)),180.d0/pi*atan2(aimag(MoM_I_z(i)),real(MoM_I_z(i))),&
            abs(MoM_I_t(i)),180.d0/pi*atan2(aimag(MoM_I_t(i)),real(MoM_I_t(i))),&
            abs(MoM_M_z(i)),180.d0/pi*atan2(aimag(MoM_M_z(i)),real(MoM_M_z(i))),&
            abs(MoM_M_t(i)),180.d0/pi*atan2(aimag(MoM_M_t(i)),real(MoM_M_t(i))),&
            abs(MoM_I_z_AWE(i)),180.d0/pi*atan2(aimag(MoM_I_z_AWE(i)),real(MoM_I_z_AWE(i))),&
            abs(MoM_I_t_AWE(i)),180.d0/pi*atan2(aimag(MoM_I_t_AWE(i)),real(MoM_I_t_AWE(i))),&
            abs(MoM_M_z_AWE(i)),180.d0/pi*atan2(aimag(MoM_M_z_AWE(i)),real(MoM_M_z_AWE(i))),&
            abs(MoM_M_t_AWE(i)),180.d0/pi*atan2(aimag(MoM_M_t_AWE(i)),real(MoM_M_t_AWE(i))),&
            abs(MoM_I_z_Chebyshev(i)),180.d0/pi*atan2(aimag(MoM_I_z_Chebyshev(i)),real(MoM_I_z_Chebyshev(i))),&
            abs(MoM_I_t_Chebyshev(i)),180.d0/pi*atan2(aimag(MoM_I_t_Chebyshev(i)),real(MoM_I_t_Chebyshev(i))),&
            abs(MoM_M_z_Chebyshev(i)),180.d0/pi*atan2(aimag(MoM_M_z_Chebyshev(i)),real(MoM_M_z_Chebyshev(i))),&
            abs(MoM_M_t_Chebyshev(i)),180.d0/pi*atan2(aimag(MoM_M_t_Chebyshev(i)),real(MoM_M_t_Chebyshev(i)))
        endif
    enddo
    if(plot_current /= 0) then
    close(fd1)
    endif

!    write(*,*) norm(MoM_I_z_AWE),norm(MoM_I_t_AWE),norm(MoM_M_z_AWE),norm(MoM_M_t_AWE)
!    write(*,*) norm(MoM_I_z),norm(MoM_I_t),norm(MoM_M_z),norm(MoM_M_t)

!    write(*,*) norm(MoM_I_z - MoM_I_z_AWE),norm(MoM_I_t - MoM_I_t_AWE),&
!        norm(MoM_M_z - MoM_M_z_AWE),norm(MoM_M_t - MoM_M_t_AWE)
!    write(*,*) MoM_I_t(10),MoM_I_t_AWE(10)
!    write(*,*) MoM_I_t(10:12) - MoM_I_t_AWE(10:12)

    allocate(Zzz(this%M,this%M),Zzt(this%M,this%M),Ztz(this%M,this%M),&
    Ztt(this%M,this%M),Yzt(this%M,this%M),Ytz(this%M,this%M),Ytt(this%M,this%M))

    !    !    call set_MoM_matrices(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ak0*k0,eta1,.true.)
    call set_MoM_BC_matrices(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ar0,az0,ak0,eta1,.true.)
    !    write(*,*) freq1/1d9
    !    write(*,*) Ztt(1,1:3)
    !    !    write(*,*) MoM_I_z(1:10)
    MoMs_Ez =  matmul(Zzz,MoM_I_z) + matmul(Zzt,MoM_I_t) - matmul(Yzt,MoM_M_t)
    MoMs_Et =  matmul(Ztz,MoM_I_z) + matmul(Ztt,MoM_I_t) - &
    matmul(Ytz,MoM_M_z) - matmul(Ytt,MoM_M_t)
    MoMs_Hz =  matmul(Yzt,MoM_I_t) + 1.d0/eta1**2.d0*(matmul(Zzz,MoM_M_z) + matmul(Zzt,MoM_M_t))
    MoMs_Ht =  matmul(Ytz,MoM_I_z) + matmul(Ytt,MoM_I_t) + &
    1.d0/eta1**2.d0*(matmul(Ztz,MoM_M_z)+matmul(Ztt,MoM_M_t))


    MoMs_AWE_Ez =  matmul(Zzz,MoM_I_z_AWE) + matmul(Zzt,MoM_I_t_AWE) - matmul(Yzt,MoM_M_t_AWE)
    MoMs_AWE_Et =  matmul(Ztz,MoM_I_z_AWE) + matmul(Ztt,MoM_I_t_AWE) - &
    matmul(Ytz,MoM_M_z_AWE) - matmul(Ytt,MoM_M_t_AWE)
    MoMs_AWE_Hz =  matmul(Yzt,MoM_I_t_AWE) + 1.d0/eta1**2.d0*(matmul(Zzz,MoM_M_z_AWE) + matmul(Zzt,MoM_M_t_AWE))
    MoMs_AWE_Ht =  matmul(Ytz,MoM_I_z_AWE) + matmul(Ytt,MoM_I_t_AWE) + &
    1.d0/eta1**2.d0*(matmul(Ztz,MoM_M_z_AWE)+matmul(Ztt,MoM_M_t_AWE))

    MoMs_Chebyshev_Ez =  matmul(Zzz,MoM_I_z_Chebyshev) + matmul(Zzt,MoM_I_t_Chebyshev) - &
        matmul(Yzt,MoM_M_t_Chebyshev)
    MoMs_Chebyshev_Et =  matmul(Ztz,MoM_I_z_Chebyshev) + matmul(Ztt,MoM_I_t_Chebyshev) - &
        matmul(Ytz,MoM_M_z_Chebyshev) - matmul(Ytt,MoM_M_t_Chebyshev)
    MoMs_Chebyshev_Hz =  matmul(Yzt,MoM_I_t_Chebyshev) + 1.d0/eta1**2.d0*(matmul(Zzz,MoM_M_z_Chebyshev) +&
        matmul(Zzt,MoM_M_t_Chebyshev))
    MoMs_Chebyshev_Ht =  matmul(Ytz,MoM_I_z_Chebyshev) + matmul(Ytt,MoM_I_t_Chebyshev) + &
        1.d0/eta1**2.d0*(matmul(Ztz,MoM_M_z_Chebyshev)+matmul(Ztt,MoM_M_t_Chebyshev))

    deallocate(Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(this%Problem_Type == 1) then !! PEC
        MoM_error_t = (sum(abs(MoM_Ez-MoMs_Ez)**2.d0)+sum(abs(MoM_Et - MoMs_Et)**2.d0))/norm_E_loc

        MoM_error_AWE = (sum(abs(MoM_Ez-MoMs_AWE_Ez)**2.d0)+sum(abs(MoM_Et - MoMs_AWE_Et)**2.d0))/norm_E_loc

        MoM_error_Chebyshev = (sum(abs(MoM_Ez-MoMs_Chebyshev_Ez)**2.d0)+&
            sum(abs(MoM_Et - MoMs_Chebyshev_Et)**2.d0))/norm_E_loc
    elseif(this%Problem_Type == 2) then !! PMC
        MoM_error_t = (sum(abs(MoM_Hz - MoMs_Hz)**2.d0)+sum(abs(MoM_Ht- MoMs_Ht)**2.d0))/norm_H_loc

        MoM_error_AWE = (sum(abs(MoM_Hz - MoMs_AWE_Hz)**2.d0)+sum(abs(MoM_Ht- MoMs_AWE_Ht)**2.d0))/norm_H_loc

        MoM_error_Chebyshev = (sum(abs(MoM_Hz - MoMs_Chebyshev_Hz)**2.d0)+&
            sum(abs(MoM_Ht- MoMs_Chebyshev_Ht)**2.d0))/norm_H_loc

    elseif(this%Problem_Type == 3) then !! IBC
        write(*,*) 'MoM AWE is not implemented for IBC'
        MoM_Et = MoM_Et + eta1*(this%eta_uu(1)*MoM_Hz - this%eta_uv(1)*MoM_Ht)
        MoM_Ez = MoM_Ez + eta1*(this%eta_vu(1)*MoM_Hz - this%eta_vv(1)*MoM_Ht)

        MoM_error_t = (sum(abs(MoM_Ez)**2.d0)+sum(abs(MoM_Et)**2.d0))/norm_E_loc

        MoM_error_AWE = 1e10
        MoM_error_Chebyshev = 1e10
    elseif(this%Problem_Type == 4) then !! Dielectric
        write(*,*) 'MoM AWE is not implemented for Dielectric Materials'
        MoM_error_t = (sum(abs(MoM_Ez)**2.d0)+sum(abs(MoM_Et)**2.d0))/norm_E_loc
        MoM_error_AWE = 1e10
        MoM_error_Chebyshev = 1e10
    else
        write(*,*) 'ERROR: wrong problem type, it should be in the range of [1,4]'
        STOP
    endif





    deallocate(MoMs_AWE_Ez,MoMs_AWE_Et,MoMs_AWE_Hz,MoMs_AWE_Ht)
    deallocate(MoMs_Chebyshev_Ez,MoMs_Chebyshev_Et,MoMs_Chebyshev_Hz,MoMs_Chebyshev_Ht)
    deallocate(MoM_Ez,MoM_Et,MoM_Hz,MoM_Ht,MoMs_Ez,MoMs_Et,MoMs_Hz,MoMs_Ht)
    deallocate(MoM_I_z_AWE,MoM_I_t_AWE,MoM_M_z_AWE,MoM_M_t_AWE)
    deallocate(MoM_I_z,MoM_I_t,MoM_M_z,MoM_M_t)
    deallocate(MoM_I_z_Chebyshev,MoM_I_t_Chebyshev,MoM_M_z_Chebyshev,MoM_M_t_Chebyshev)
end subroutine eval_error_frequency_MoM_AWE_Chebyshev

subroutine eval_error_frequency(this,freq1,N_order,error_t,MoM_error_t)
    type(Scatterer) :: this
    integer,intent(in) :: N_order
    real(8),allocatable,intent(inout) :: error_t(:)
    real(8) :: MoM_error_t
    complex*16,allocatable,dimension(:) :: H_u,H_v
    type(vector_c),allocatable,dimension(:,:) :: E_tot,H_tot
    complex*16,allocatable,dimension(:,:) :: BC_total
    integer:: j,kk,nn
    type(vector) :: Ep,En,rho_hat,test_pt
    real(8) :: ak0,ar0,az0,az1,freq1
    complex*16 :: ak1,ar1
    complex*16 :: exponential
    type(vector_c),allocatable,dimension(:)  :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz
    complex*16,allocatable,dimension(:)  :: MoM_Ez,MoM_Et,MoM_Hz,MoM_Ht
    complex*16,allocatable,dimension(:)  :: MoMs_Ez,MoMs_Et,MoMs_Hz,MoMs_Ht
    type(vector_c),allocatable,dimension(:,:) :: I_current
    real(8) :: error_u,error_v,norm_E_loc, norm_H_loc,const_freq
    complex*16,allocatable,dimension(:) :: MoM_I_z,MoM_I_t,MoM_M_z,MoM_M_t
    complex*16,allocatable,dimension(:,:) :: Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt
    complex*16 :: eta_ex,ar_ex,ak_ex

    ak0 = ak*freq1/frequency
    ar0 = ar*freq1/frequency
    az0 = az*freq1/frequency

    ak1 = this%ak_local*freq1/frequency
    ar1 = this%ar_local*freq1/frequency
    az1 = az0
    allocate(MoM_Ez(this%M),MoM_Et(this%M),MoM_Hz(this%M),MoM_Ht(this%M))
    allocate(MoMs_Ez(this%M),MoMs_Et(this%M),MoMs_Hz(this%M),MoMs_Ht(this%M))
    allocate(MoM_I_z(this%M),MoM_I_t(this%M),MoM_M_z(this%M),MoM_M_t(this%M))
    allocate(BC_total(2*this%M,0:N_order),E_tot(this%M,0:N_order),H_tot(this%M,0:N_order))
    allocate(I_current(this%N,0:N_order))
    allocate(E_Ix(0:0),E_Iy(0:0),E_Iz(0:0),H_Ix(0:0),H_Iy(0:0),H_Iz(0:0))
    Ep = vec(-cos(theta_i)*cos(phi_i),-cos(theta_i)*sin(phi_i),sin(theta_i))
    En = vec(sin(phi_i),-cos(phi_i),0.d0)
    rho_hat = vec(sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i))
    norm_E_loc = 0.d0
    norm_H_loc =0.d0
    do j=1,this%M
        test_pt = vec(this%testing_pt(j,1),this%testing_pt(j,2),0.d0)
        exponential = exp(cj*ak0*k0*dot(rho_hat,test_pt))
        !                write(*,*) dot(rho_hat,testing_pt_MoM(j))
        E_tot(j,:) = eta1*exponential*(cos(alpha_i)*Ep+sin(alpha_i)*En)
        !                write(*,*) Ei(j)%v
        H_tot(j,:) = exponential*(sin(alpha_i)*Ep-cos(alpha_i)*En)
        MoM_Ez(j) = dot(this%tang_v(j),E_tot(j,0))
        MoM_Et(j) = dot(this%tang_u(j),E_tot(j,0))
        MoM_Hz(j) = dot(this%tang_v(j),H_tot(j,0))
        MoM_Ht(j) = dot(this%tang_u(j),H_tot(j,0))
        norm_E_loc = norm_E_loc + abs(MoM_Et(j))**2.d0 + abs(MoM_Ez(j))**2.d0
        norm_H_loc = norm_H_loc + abs(MoM_Hz(j))**2.d0 + abs(MoM_Ht(j))**2.d0
    enddo
    !! computation of the currents vector
    I_current = vec(c_zero,c_zero,c_zero)
    I_current(1:this%N,0) = this%I_sources_total(1:this%N,0)
    MoM_I_z = this%I_v(:,0)
    MoM_I_t = this%I_u(:,0)
    MoM_M_z = this%M_v(:,0)
    MoM_M_t = this%M_u(:,0)
    const_freq = (ak0-ak)*k0
    !    write(*,*) 'MoM_I_z',sum(MoM_I_z),sum(this%I_v(:,0))
    do nn = 1,N_order
        I_current(1:this%N,nn) = I_current(1:this%N,nn-1) + const_freq*this%I_sources_total(1:this%N,nn)
        MoM_I_z = MoM_I_z + const_freq*this%I_v(:,nn)
        !        write(*,*) 'MoM_I_z',sum(MoM_I_z),sum(this%I_v(:,nn))
        !        write(*,*) MoM_I_z(10)
        MoM_I_t = MoM_I_t + const_freq*this%I_u(:,nn)
        MoM_M_z = MoM_M_z + const_freq*this%M_v(:,nn)
        MoM_M_t = MoM_M_t + const_freq*this%M_u(:,nn)

        const_freq = const_freq*(ak0-ak)*k0
    !        write(*,*) const_freq
    enddo

!    do nn = 1,this%N
!        write(*,*) I_current(nn,1),I_current(nn,2)
!    enddo
!    stop

    allocate(Zzz(this%M,this%M),Zzt(this%M,this%M),Ztz(this%M,this%M),&
    Ztt(this%M,this%M),Yzt(this%M,this%M),Ytz(this%M,this%M),Ytt(this%M,this%M))

    !    !    call set_MoM_matrices(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ak0*k0,eta1,.true.)
    call set_MoM_BC_matrices(this,Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt,ar0,az0,ak0,eta1,.true.)
    !    write(*,*) freq1/1d9
    !    write(*,*) Ztt(1,1:3)
    !    !    write(*,*) MoM_I_z(1:10)
!    MoMs_Ez =  matmul(Zzz,MoM_I_z) + matmul(Zzt,MoM_I_t) - matmul(Yzt,MoM_M_t)
!    MoMs_Et =  matmul(Ztz,MoM_I_z) + matmul(Ztt,MoM_I_t) - &
!    matmul(Ytz,MoM_M_z) - matmul(Ytt,MoM_M_t)
!    MoMs_Hz =  matmul(Yzt,MoM_I_t) + 1.d0/eta1**2.d0*(matmul(Zzz,MoM_M_z) + matmul(Zzt,MoM_M_t))
!    MoMs_Ht =  matmul(Ytz,MoM_I_z) + matmul(Ytt,MoM_I_t) + &
!    1.d0/eta1**2.d0*(matmul(Ztz,MoM_M_z)+matmul(Ztt,MoM_M_t))


    MoMs_Ez =  matmul(Zzz,MoM_I_z) + matmul(Zzt,MoM_I_t)
    MoMs_Et =  matmul(Ztz,MoM_I_z) + matmul(Ztt,MoM_I_t)
    MoMs_Hz =  matmul(Yzt,MoM_I_t)
    MoMs_Ht =  matmul(Ytz,MoM_I_z) + matmul(Ytt,MoM_I_t)

    deallocate(Zzz,Zzt,Ztz,Ztt,Yzt,Ytz,Ytt)
    do j = 1,this%M
        do kk = 1,this%N
            if(this%active_region(kk)==0) then !! inside sources
                eta_ex = eta1
                ar_ex = ar0
                ak_ex = ak0
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                endif
                do nn = 0,N_taylor
                    E_tot(j,nn) = E_tot(j,nn) + I_current(kk,nn)%v(1)*E_Ix(0) + I_current(kk,nn)%v(2)*E_Iy(0) +&
                    I_current(kk,nn)%v(3)*E_Iz(0)
                    H_tot(j,nn) = H_tot(j,nn) + I_current(kk,nn)%v(1)*H_Ix(0) + I_current(kk,nn)%v(2)*H_Iy(0) +&
                    I_current(kk,nn)%v(3)*H_Iz(0)
                enddo
            else !! outside sources
                eta_ex = this%eta_local
                ar_ex = ar1
                ak_ex = ak1
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                endif
                !                E_tot(j) = E_tot(j) - this%I_sources_total(kk,1)*E_Ix(0) - this%I_sources_total(kk,2)*E_Iy(0) -&
                !                this%I_sources_total(kk,3)*E_Iz(0)
                !                H_tot(j) = H_tot(j) - this%I_sources_total(kk,1)*H_Ix(0) - this%I_sources_total(kk,2)*H_Iy(0) -&
                !                this%I_sources_total(kk,3)*H_Iz(0)
                !                E_tot(j) = E_tot(j) - I_current(kk)%v(1)*E_Ix(0) - I_current(kk)%v(2)*E_Iy(0) -&
                !                I_current(kk)%v(3)*E_Iz(0)
                !                H_tot(j) = H_tot(j) - I_current(kk)%v(1)*H_Ix(0) - I_current(kk)%v(2)*H_Iy(0) -&
                !                I_current(kk)%v(3)*H_Iz(0)
                do nn = 0,N_order
                    E_tot(j,nn) = E_tot(j,nn) - I_current(kk,nn)%v(1)*E_Ix(0) - I_current(kk,nn)%v(2)*E_Iy(0) -&
                    I_current(kk,nn)%v(3)*E_Iz(0)
                    H_tot(j,nn) = H_tot(j,nn) - I_current(kk,nn)%v(1)*H_Ix(0) - I_current(kk,nn)%v(2)*H_Iy(0) -&
                    I_current(kk,nn)%v(3)*H_Iz(0)
                enddo
            endif

        enddo

    enddo

    if(this%Problem_Type == 1) then !! PEC
        !        do nn = 0,N_order
        do j =1,this%M
            BC_total(j,:) = dot(this%tang_u(j),E_tot(j,:),N_order+1)
            BC_total(j+this%M,:) = dot(this%tang_v(j),E_tot(j,:),N_order+1)
        enddo

        !        write(*,*) abs(MoM_Et(1:5))
        !        write(*,*) abs(MoMs_Et(1:5))
        !        write(*,*) abs(MoM_Et(1:5) - MoMs_Et(1:5) )
        !        write(*,*) '=================================='
        !        do j=1,this%M
        !            MoMs_Hz(j) = dot(this%tang_v(j),H_tot(j,N_order))
        !
        !            MoMs_Ht(j) = dot(this%tang_u(j),H_tot(j,N_order))
        !        enddo

        MoM_error_t = (sum(abs(MoM_Ez-MoMs_Ez)**2.d0)+sum(abs(MoM_Et - MoMs_Et)**2.d0))/norm_E_loc
    !        write(*,*) 'norm_H',norm_H_loc
    !        write(*,*) 'sum I_z',sum(abs(MoM_I_z))
    !        write(*,*) 'sum Ht',sum(abs(MoM_Ht))
    !        MoM_error_t = (sum(abs(MoMs_Hz+MoM_I_t)**2.d0)+sum(abs(MoMs_Ht - MoM_I_z)**2.d0))/norm_H_loc
    !        write(*,*) MoM_error_t,sum(abs(MoM_Ez)**2.d0),sum(abs(MoM_Et)**2.d0),norm_E_loc
    !        enddo

    elseif(this%Problem_Type == 2) then !! PMC
        do nn = 0,N_order
            do j =1,this%M
                BC_total(j,nn) = dot(this%tang_u(j),H_tot(j,nn))
                BC_total(j+this%M,nn) = dot(this%tang_v(j),H_tot(j,nn))
            enddo
        enddo

        MoM_error_t = (sum(abs(MoM_Hz - MoMs_Hz)**2.d0)+sum(abs(MoM_Ht- MoMs_Ht)**2.d0))/norm_H_loc

    elseif(this%Problem_Type == 3) then !! IBC
        allocate(H_u(this%M),H_v(this%M))

        do nn = 0,N_order
            do j = 1,this%M
                H_u(j) = dot(this%tang_u(j),H_tot(j,nn))
                H_v(j) = dot(this%tang_v(j),H_tot(j,nn))
                BC_total(j,nn) = dot(this%tang_u(j),E_tot(j,nn))
                BC_total(this%M+j,nn) = dot(this%tang_v(j),E_tot(j,nn))
            enddo
            H_u = eta1*H_u
            H_v = eta1*H_v
            do j = 1,this%M
                BC_total(j,nn) = BC_total(j,nn) +  this%eta_uu(j)*H_v(j) - this%eta_uv(j)*H_u(j) !! E_u
                BC_total(this%M+j,nn) = BC_total(this%M+j,nn) +  this%eta_vu(j)*H_v(j) - this%eta_vv(j)*H_u(j) !! E_v
            enddo
        enddo
        deallocate(H_u,H_v)

        MoM_Et = MoM_Et + eta1*(this%eta_uu(1)*MoM_Hz - this%eta_uv(1)*MoM_Ht)
        MoM_Ez = MoM_Ez + eta1*(this%eta_vu(1)*MoM_Hz - this%eta_vv(1)*MoM_Ht)

        MoM_error_t = (sum(abs(MoM_Ez)**2.d0)+sum(abs(MoM_Et)**2.d0))/norm_E_loc

    elseif(this%Problem_Type == 4) then !! Dielectric
        do nn=0,N_order
            do j =1,this%M
                BC_total(j,nn) = dot(this%tang_u(j),E_tot(j,nn))
                BC_total(j+this%M,nn) = dot(this%tang_v(j),E_tot(j,nn))
            enddo
        enddo
        MoM_error_t = (sum(abs(MoM_Ez)**2.d0)+sum(abs(MoM_Et)**2.d0))/norm_E_loc
    else
        write(*,*) 'ERROR: wrong problem type, it should be in the range of [1,4]'
        STOP
    endif

    do nn = 0,N_order
        error_u = 0.d0
        error_v = 0.d0
        do j=1,this%M
            error_u = error_u + abs(BC_total(j,nn))**2.d0
            error_v = error_v + abs(BC_total(j+this%M,nn))**2.d0

        enddo
        if(this%Problem_Type == 2) then !! PMC
            error_t(nn) = (error_v + error_u)/norm_H_loc
        else
            error_t(nn) = (error_v + error_u)/norm_E_loc
        endif
    enddo


    !        write(*,*) 'norm_E=',this%normalize_E,'norm_E_loc =',norm_E_loc

    deallocate(BC_total,E_tot,H_tot)
    deallocate(I_current)
    deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz)
    deallocate(MoM_I_z,MoM_I_t,MoM_M_z,MoM_M_t)
    deallocate(MoM_Ez,MoM_Et,MoM_Hz,MoM_Ht,MoMs_Ez,MoMs_Et,MoMs_Hz,MoMs_Ht)
end subroutine eval_error_frequency

subroutine eval_MoM_Pade_coefficients(this,Q,L)
    type(Scatterer) :: this
    integer :: Q !! it should be the evaluated taylor series coefficients max. order
    integer :: L,M
    complex*16,allocatable,dimension(:,:) :: Pade_matrix
    complex*16,allocatable,dimension(:) :: Pade_excit,m_vec,a_vec
    integer :: ii,jj,scr_cnt,cmp_cnt
    real(8) :: norm_val

    if(L > Q) then
        write(*,*) 'ERROR: Pade argument L should not be larger than Q'
        return
    endif
    M = Q-L

!    write(*,*) Q,L
    this%Pade_L = L
    this%Pade_M = M

    allocate(this%MoM_Pade_a_I_u(this%M,0:L),this%MoM_Pade_b_I_u(this%M,0:M))
    allocate(this%MoM_Pade_a_I_v(this%M,0:L),this%MoM_Pade_b_I_v(this%M,0:M))
    allocate(this%MoM_Pade_a_M_u(this%M,0:L),this%MoM_Pade_b_M_u(this%M,0:M))
    allocate(this%MoM_Pade_a_M_v(this%M,0:L),this%MoM_Pade_b_M_v(this%M,0:M))

    allocate(Pade_matrix(1:M,1:M),Pade_excit(1:M))
    allocate(m_vec(0:Q),a_vec(0:L))

    this%MoM_Pade_a_I_u = 0.d0
    this%MoM_Pade_b_I_u = 0.d0
    this%MoM_Pade_a_M_u = 0.d0
    this%MoM_Pade_b_M_u = 0.d0

    this%MoM_Pade_a_I_v = 0.d0
    this%MoM_Pade_b_I_v = 0.d0
    this%MoM_Pade_a_M_v = 0.d0
    this%MoM_Pade_b_M_v = 0.d0

    this%MoM_Pade_b_I_u(:,0) = 1.d0
    this%MoM_Pade_b_I_v(:,0) = 1.d0
    this%MoM_Pade_b_M_u(:,0) = 1.d0
    this%MoM_Pade_b_M_v(:,0) = 1.d0

!    write(*,*) size(this%I_u(scr_cnt,:))

    if(this%Problem_Type == 1 .or. this%Problem_Type == 3 .or. this%Problem_Type == 4 ) then
!        write(*,*) 'I am here 1'

        do scr_cnt = 1,this%M
            m_vec = this%I_u(scr_cnt,:)
            norm_val = dble(Q)/sum(abs(m_vec(L:(Q-1))))
!            write(*,*) 'I am here 2'
            do ii = 1,M
                Pade_excit(ii) = -m_vec(L+ii)
                do jj = 1,M
                    Pade_matrix(ii,jj) = norm_val*m_vec(L+ii-jj)
                enddo

            enddo
!            write(*,*) 'I am here 3'
            this%MoM_Pade_b_I_u(scr_cnt,0) = 1.d0
            this%MoM_Pade_b_I_u(scr_cnt,1:M) = norm_val*Gauss_Elimination(Pade_matrix,Pade_excit,M)

            a_vec = 0.d0
            do ii = 0,L
                do jj = 0,ii
                    a_vec(ii) = a_vec(ii) + this%MoM_Pade_b_I_u(scr_cnt,jj)*m_vec(ii-jj)
                enddo
            enddo
            this%MoM_Pade_a_I_u(scr_cnt,:) = a_vec
!            write(*,*) 'I am here 4'
        enddo

!        write(*,*) 'I am here 5'
        do scr_cnt = 1,this%M
            m_vec = this%I_v(scr_cnt,:)
            norm_val = dble(Q)/sum(abs(m_vec(L:(Q-1))))

            do ii = 1,M
                Pade_excit(ii) = -m_vec(L+ii)
                do jj = 1,M
                    Pade_matrix(ii,jj) = norm_val*m_vec(L+ii-jj)
                enddo

            enddo
            this%MoM_Pade_b_I_v(scr_cnt,0) = 1.d0
            this%MoM_Pade_b_I_v(scr_cnt,1:M) = norm_val*Gauss_Elimination(Pade_matrix,Pade_excit,M)

            a_vec = 0.d0
            do ii = 0,L
                do jj = 0,ii
                    a_vec(ii) = a_vec(ii) + this%MoM_Pade_b_I_v(scr_cnt,jj)*m_vec(ii-jj)
                enddo
            enddo
            this%MoM_Pade_a_I_v(scr_cnt,:) = a_vec
        enddo

    endif


    if(this%Problem_Type == 2 .or. this%Problem_Type == 3 .or. this%Problem_Type == 4 ) then
        do scr_cnt = 1,this%M
            m_vec = this%M_u(scr_cnt,:)
            norm_val = dble(Q)/sum(abs(m_vec(L:(Q-1))))

            do ii = 1,M
                Pade_excit(ii) = -m_vec(L+ii)
                do jj = 1,M
                    Pade_matrix(ii,jj) = norm_val*m_vec(L+ii-jj)
                enddo

            enddo
            this%MoM_Pade_b_M_u(scr_cnt,0) = 1.d0
            this%MoM_Pade_b_M_u(scr_cnt,1:M) = norm_val*Gauss_Elimination(Pade_matrix,Pade_excit,M)

            a_vec = 0.d0
            do ii = 0,L
                do jj = 0,ii
                    a_vec(ii) = a_vec(ii) + this%MoM_Pade_b_M_u(scr_cnt,jj)*m_vec(ii-jj)
                enddo
            enddo
            this%MoM_Pade_a_M_u(scr_cnt,:) = a_vec
        enddo

        do scr_cnt = 1,this%M
            m_vec = this%M_v(scr_cnt,:)
            norm_val = dble(Q)/sum(abs(m_vec(L:(Q-1))))

            do ii = 1,M
                Pade_excit(ii) = -m_vec(L+ii)
                do jj = 1,M
                    Pade_matrix(ii,jj) = norm_val*m_vec(L+ii-jj)
                enddo

            enddo
            this%MoM_Pade_b_M_v(scr_cnt,0) = 1.d0
            this%MoM_Pade_b_M_v(scr_cnt,1:M) = norm_val*Gauss_Elimination(Pade_matrix,Pade_excit,M)

            a_vec = 0.d0
            do ii = 0,L
                do jj = 0,ii
                    a_vec(ii) = a_vec(ii) + this%MoM_Pade_b_M_v(scr_cnt,jj)*m_vec(ii-jj)
                enddo
            enddo
            this%MoM_Pade_a_M_v(scr_cnt,:) = a_vec
        enddo
    endif

    deallocate(Pade_matrix,Pade_excit,m_vec,a_vec)
end subroutine eval_MoM_Pade_coefficients

    ! -----------------------------------------------------------------------
    ! Subroutine: eval_Pade_coefficients
    ! Purpose   : Computes the Pade [L/M] rational-function approximant
    !             coefficients (Pade_a, Pade_b) for the RAS source currents
    !             from the stored Taylor series I_sources_total(:, 0:Q).
    !             Solves the linear Pade system using Gauss elimination.
    !             L poles, M = Q-L zeros.
    ! -----------------------------------------------------------------------
subroutine eval_Pade_coefficients(this,Q,L)
    type(Scatterer) :: this
    integer :: Q !! it should be the evaluated taylor series coefficients max. order
    integer :: L,M
    complex*16,allocatable,dimension(:,:) :: Pade_matrix
    complex*16,allocatable,dimension(:) :: Pade_excit,m_vec,a_vec
    integer :: ii,jj,scr_cnt,cmp_cnt
    real(8) :: norm_val
    !    complex*16 :: matrix_det

!    write(*,*) Q,L

    if(L > Q) then
        write(*,*) 'ERROR: Pade argument L should not be larger than Q'
        return
    endif
    M = Q-L

    allocate(this%Pade_a(this%N,0:L),this%Pade_b(this%N,0:M))
    allocate(Pade_matrix(1:M,1:M),Pade_excit(1:M))
    allocate(m_vec(0:Q),a_vec(0:L))


!    write(*,*) this%I_sources_total(10,:)

    do scr_cnt = 1,this%N
        do cmp_cnt = 1,3
            m_vec = this%I_sources_total(scr_cnt,:)%v(cmp_cnt)
            !            norm_val = 1.d0/maxval(abs(m_vec(2:(Q-1))))
            norm_val = dble(Q)/sum(abs(m_vec(L:(Q-1))))
            !            if(sum(abs(m_vec((L+1):Q))) > 1d-7 ) then
            do ii = 1,M
                Pade_excit(ii) = -m_vec(L+ii)
                do jj = 1,M
                    Pade_matrix(ii,jj) = norm_val*m_vec(L+ii-jj)
                enddo
            !                matrix_det = det(Pade_matrix)
            !                if(abs(matrix_det) < 1d-8) then
            !                    write(*,*) abs(matrix_det)
            !                endif
            enddo
            this%Pade_b(scr_cnt,0)%v(cmp_cnt) = 1.d0
            this%Pade_b(scr_cnt,1:M)%v(cmp_cnt) = norm_val*Gauss_Elimination(Pade_matrix,Pade_excit,M)
            !            else
            !            write(*,*) norm_val,abs(m_vec)
            !            this%Pade_b(scr_cnt,0)%v(cmp_cnt) = 1.d0
            !            this%Pade_b(scr_cnt,1:M)%v(cmp_cnt) = 0.d0
            !            endif
            a_vec = 0.d0
            do ii = 0,L
                do jj = 0,ii
                    a_vec(ii) = a_vec(ii) + this%Pade_b(scr_cnt,jj)%v(cmp_cnt)*m_vec(ii-jj)
                enddo
            enddo
            this%Pade_a(scr_cnt,:)%v(cmp_cnt) = a_vec
        enddo
    enddo

!    write(*,*) this%Pade_a(scr_cnt,:)

    this%Pade_L = L
    this%Pade_M = M

    deallocate(Pade_matrix,Pade_excit,m_vec,a_vec)
end subroutine eval_Pade_coefficients

function Get_Kappa(ak_l,akh_r,akl_r) result(Kappa)
    real(8) :: ak_l,akh_r,akl_r,Kappa

    !    Kappa = (2.d0 - (fh_f0+fl_f0) )/(fh_f0-fl_f0)
    !    Kappa = 0.d0
!    !!!!! Fomulation 1 --> Not working
!    Kappa = (ak_l - 1.d0)*k0
    !!!!! Fomulation 2
    Kappa = (2.d0*ak_l - (akh_r+akl_r)) /(akh_r-akl_r) + 2.d0
!    !!!!! Fomulation 3
!    Kappa = (ak_l - 1.d0)/(akh_r-akl_r)
!    !!!!! Fomulation 4
!    Kappa = (ak_l - 1.d0)/(akh_r-1.d0)
!    !!!!! Fomulation 5
!    Kappa = -1.d0*ak_l/1.d0
!    !!!!! Fomulation 6
!    Kappa = (2.d0 - (akh_r+akl_r)) /(akh_r-akl_r)
!    !!!!! Fomulation 7
!    Kappa = -1.d0*ak_l*(akh_r-akl_r)/(akh_r+akl_r)*k0 + 1.d0

    if(Kappa == 1.d0) then
        Kappa = 0.9999d0
    endif

end function Get_Kappa

function Get_dKappa_dk(ak_l,akh_r,akl_r) result(dKappa_dk)
    real(8) :: ak_l,akh_r,akl_r,dKappa_dk
!    dKappa_dk = 2.d0/(k0*(fh_f0-fl_f0))

!    !!!!! Fomulation 1 --> Not working
!    dKappa_dk = 1.d0
    !!!!! Fomulation 2
    dKappa_dk = 2.d0/(k0*(akh_r-akl_r))
!    !!!!! Fomulation 3
!    dKappa_dk = 1.d0/(k0*(akh_r-akl_r))
!    !!!!! Fomulation 4
!    dKappa_dk = 1.d0/(k0*(akh_r-1.d0))
!    !!!!! Fomulation 5
!    dKappa_dk = -1.d0/(k0)
!    !!!!! Fomulation 6
!    dKappa_dk = 1.d0/(k0*(akh_r-akl_r))
!    !!!!! Fomulation 7
!    dKappa_dk = -1.d0*(akh_r-akl_r)/(akh_r+akl_r)
end function Get_dKappa_dk


    ! -----------------------------------------------------------------------
    ! Subroutine: Eval_chebyshev_expansion_coeff_generic
    ! Purpose   : Converts a Taylor polynomial of order Q into a Chebyshev
    !             series expansion of order L on the normalised frequency
    !             interval [-1, 1]. Uses the Chebyshev derivative matrix T
    !             from chebyshevT_D. Outputs coefficient vectors c_vec and
    !             d_vec for the numerator and denominator of the rational
    !             Chebyshev approximant. Used for AWE wideband reconstruction.
    ! -----------------------------------------------------------------------
subroutine Eval_chebyshev_expansion_coeff_generic(Q,L,T,x_vec,c_vec,d_vec)
    integer :: Q,L,M
    real(8),allocatable,intent(inout) :: T(:,:)
    complex*16,allocatable,intent(inout) :: x_vec(:)
    integer :: mm,qq,ii,q_ind,T_order_max
    complex*16,allocatable,dimension(:,:) :: Ch_matrix,Ch_matrix_M
    complex*16,allocatable,dimension(:) :: CH_excit,c_vec_excit
    complex*16,allocatable,intent(inout) :: c_vec(:),d_vec(:)
    real(8) :: norm_val

    M = Q - L
    allocate(CH_excit(M),Ch_matrix(M,M),c_vec_excit(0:L))
    allocate(Ch_matrix_M(L+1,M))

    d_vec(0) = 1.d0

    T_order_max = max(L,M)

    if(M > 0) then
        CH_excit = x_vec((Q-M+1):Q)
        norm_val = dble(Q)/sum(abs(CH_excit))
!        norm_val = 1.d0
    !    write(*,*) norm_val,T(T_order_max,T_order_max)
    !    norm_val = T(Q,Q)
        do mm = 1,M
            q_ind = 1
            do qq = (Q-M+1),Q
                Ch_matrix(q_ind,mm) = 0.d0
                do ii = 0,qq
                    Ch_matrix(q_ind,mm) = Ch_matrix(q_ind,mm) -  &
                    norm_val*dble(n_C_i_matrix(qq,ii))*x_vec(qq-ii)*T(ii,mm)
                enddo
    !            Ch_matrix(q_ind,mm) =   Ch_matrix(q_ind,mm)*T(0,q_ind)

                q_ind = q_ind+1
            enddo
    !        CH_excit(mm) = CH_excit(mm)*T(0,mm)
        enddo

        d_vec(1:M) = norm_val*Gauss_Elimination(Ch_matrix,CH_excit,M)


        write(*,*) det(Ch_matrix)
        call write_matrix_MATLAB(Ch_matrix,'Ch_matrix')
        call  write_vector_MATLAB(norm_val*CH_excit,'CH_excit')


        do mm = 1,M
            do qq = 0,L
                Ch_matrix_M(qq+1,mm) = 0.d0
                do ii = 0,qq
                    Ch_matrix_M(qq+1,mm) = Ch_matrix_M(qq+1,mm) -  &
                    dble(n_C_i_matrix(qq,ii))*x_vec(qq-ii)*T(ii,mm)
                enddo
            enddo
        enddo
        c_vec_excit(0:L) = x_vec(0:L) - matmul(Ch_matrix_M,d_vec(1:M))

    else
        c_vec_excit = x_vec
        c_vec(L) = c_vec_excit(L)/T(L,L)
        do ii = L-1,0,-1
            c_vec(ii) = c_vec_excit(ii)
            do mm = L,ii+1,-1
                c_vec(ii) = c_vec(ii) - c_vec(mm)*T(ii,mm)
            enddo
            c_vec(ii) = c_vec(ii)/T(ii,ii)
        enddo
    endif


    deallocate(CH_excit,Ch_matrix,c_vec_excit,Ch_matrix_M)
end subroutine Eval_chebyshev_expansion_coeff_generic

subroutine Eval_MoM_chebyshev_expansion_coefficeints(this,Q,L,fl_f0,fh_f0)
    type(Scatterer) :: this
    real(8) :: fl_f0,fh_f0
    integer :: Q !! it should be the evaluated taylor series coefficients max. order
    integer :: M,L
    complex*16,allocatable,dimension(:,:) :: Ch_matrix,ch_matrix_inverted
    real(8),allocatable,dimension(:,:) :: T
    real(8),allocatable,dimension(:) ::  norm_array
    complex*16,allocatable,dimension(:) :: CH_excit,c_vec,x_vec,c_vec_temp,d_vec_temp
    integer :: ii,nn,scr_cnt,cmp_cnt,qq,mm
    real(8) :: d_Kappa_by_dk,Kappa,d_Kappa_by_dk_current

    !    complex*16 :: matrix_det

    M = Q-L
    if(L > Q) then
        CALL PERROR('Bad Choice for L and Q in The Chebyshev Wideband Representation of the Solution')
        stop
    endif


    allocate(norm_array(0:Q))
    allocate(this%MoM_Chebyshev_c_I_u(1:this%M,0:L),this%MoM_Chebyshev_d_I_u(1:this%M,0:M))
    allocate(this%MoM_Chebyshev_c_M_u(1:this%M,0:L),this%MoM_Chebyshev_d_M_u(1:this%M,0:M))
    allocate(this%MoM_Chebyshev_c_I_v(1:this%M,0:L),this%MoM_Chebyshev_d_I_v(1:this%M,0:M))
    allocate(this%MoM_Chebyshev_c_M_v(1:this%M,0:L),this%MoM_Chebyshev_d_M_v(1:this%M,0:M))

    this%MoM_Chebyshev_c_I_u = 0.d0
    this%MoM_Chebyshev_d_I_u = 0.d0
    this%MoM_Chebyshev_c_M_u = 0.d0
    this%MoM_Chebyshev_d_M_u = 0.d0

    this%MoM_Chebyshev_c_I_v = 0.d0
    this%MoM_Chebyshev_d_I_v = 0.d0
    this%MoM_Chebyshev_c_M_v = 0.d0
    this%MoM_Chebyshev_d_M_v = 0.d0

    this%MoM_Chebyshev_d_I_u(:,0) = 1.d0
    this%MoM_Chebyshev_d_I_v(:,0) = 1.d0
    this%MoM_Chebyshev_d_M_u(:,0) = 1.d0
    this%MoM_Chebyshev_d_M_v(:,0) = 1.d0

    allocate(c_vec(1:(1+Q)),CH_excit(Q+1),Ch_matrix(1+Q,1+Q),Ch_matrix_inverted(1+Q,1+Q))
    allocate(x_vec(0:Q))

    allocate(c_vec_temp(0:L),d_vec_temp(0:M))





    Kappa = get_Kappa(1.d0,fh_f0,fl_f0)


    d_Kappa_by_dk = Get_dKappa_dk(1.d0,fh_f0,fl_f0)

    Ch_matrix(1,1) = 1.d0
    Ch_matrix(2:(Q+1),1) = 0.d0

!    write(*,*) Kappa
    call chebyshevT_D(Q,Kappa,T)
!    do ii=0,Q
!        write(*,*) T(ii,:)
!    enddo
!    write(*,*)
    norm_array(0) = 1.d0
    d_Kappa_by_dk_current = d_Kappa_by_dk
    do qq=1,Q
        T(qq,:) = T(qq,:)*d_Kappa_by_dk_current
        norm_array(qq) = 1.d0/d_Kappa_by_dk_current
        d_Kappa_by_dk_current = d_Kappa_by_dk_current*d_Kappa_by_dk
    enddo

!    do ii=0,Q
!        write(*,*) T(ii,:)
!    enddo

    Ch_matrix = 0.d0

    Ch_matrix(1:(L+1),1:(L+1)) = T(0:L,0:L)
    do qq=1,L
        Ch_matrix(qq+1,1:(L+1)) = Ch_matrix(qq+1,1:(L+1))*norm_array(qq)
    enddo

    if(this%Problem_Type == 1 .or. this%Problem_Type == 3 .or. this%Problem_Type == 4 ) then

        do scr_cnt = 1,this%M
            CH_excit = this%I_u(scr_cnt,:)
            x_vec = CH_excit
            do mm = 1,M
                do qq = 0,Q
                    Ch_matrix(qq+1,L+1+mm) = 0.d0
                    do ii = 0,qq
                        Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm) -  &
                        dble(n_C_i_matrix(qq,ii))*x_vec(qq-ii)*T(ii,mm)
                    enddo
                    Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm)*norm_array(qq)
                enddo
            enddo
            CH_excit = CH_excit*norm_array
            c_vec = Gauss_Elimination(Ch_matrix,CH_excit,Q+1)

            this%MoM_Chebyshev_c_I_u(scr_cnt,0:L) = c_vec(1:(L+1))
            this%MoM_Chebyshev_d_I_u(scr_cnt,0) = (1.d0,0.d0)
            this%MoM_Chebyshev_d_I_u(scr_cnt,1:M) = c_vec((L+2):(Q+1))

!            call Eval_chebyshev_expansion_coeff_generic(Q,L,T,x_vec,c_vec_temp,d_vec_temp)
!
!            write(*,*) this%MoM_Chebyshev_c_I_u(scr_cnt,:)
!            write(*,*) c_vec_temp
!
!            stop
        enddo

        do scr_cnt = 1,this%M
            CH_excit = this%I_v(scr_cnt,:)
            x_vec = CH_excit
            do mm = 1,M
                do qq = 0,Q
                    Ch_matrix(qq+1,L+1+mm) = 0.d0
                    do ii = 0,qq
                        Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm) -  &
                        dble(n_C_i_matrix(qq,ii))*x_vec(qq-ii)*T(ii,mm)
                    enddo
                    Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm)*norm_array(qq)
                enddo
            enddo
            CH_excit = CH_excit*norm_array
            c_vec = Gauss_Elimination(Ch_matrix,CH_excit,Q+1)

            this%MoM_Chebyshev_c_I_v(scr_cnt,0:L) = c_vec(1:(L+1))
            this%MoM_Chebyshev_d_I_v(scr_cnt,0) = (1.d0,0.d0)
            this%MoM_Chebyshev_d_I_v(scr_cnt,1:M) = c_vec((L+2):(Q+1))
        enddo


    endif

    if(this%Problem_Type == 2 .or. this%Problem_Type == 3  .or. this%Problem_Type == 4) then
        do scr_cnt = 1,this%M
            CH_excit = this%M_u(scr_cnt,:)
            x_vec = CH_excit
            do mm = 1,M
                do qq = 0,Q
                    Ch_matrix(qq+1,L+1+mm) = 0.d0
                    do ii = 0,qq
                        Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm) -  &
                        dble(n_C_i_matrix(qq,ii))*x_vec(qq-ii)*T(ii,mm)
                    enddo
                    Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm)*norm_array(qq)
                enddo
            enddo
            CH_excit = CH_excit*norm_array
            c_vec = Gauss_Elimination(Ch_matrix,CH_excit,Q+1)

            this%MoM_Chebyshev_c_M_u(scr_cnt,0:L) = c_vec(1:(L+1))
            this%MoM_Chebyshev_d_M_u(scr_cnt,0) = (1.d0,0.d0)
            this%MoM_Chebyshev_d_M_u(scr_cnt,1:M) = c_vec((L+2):(Q+1))
        enddo

        do scr_cnt = 1,this%M
            CH_excit = this%M_v(scr_cnt,:)
            x_vec = CH_excit
            do mm = 1,M
                do qq = 0,Q
                    Ch_matrix(qq+1,L+1+mm) = 0.d0
                    do ii = 0,qq
                        Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm) -  &
                        dble(n_C_i_matrix(qq,ii))*x_vec(qq-ii)*T(ii,mm)
                    enddo
                    Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm)*norm_array(qq)
                enddo
            enddo
            CH_excit = CH_excit*norm_array
            c_vec = Gauss_Elimination(Ch_matrix,CH_excit,Q+1)

            this%MoM_Chebyshev_c_M_v(scr_cnt,0:L) = c_vec(1:(L+1))
            this%MoM_Chebyshev_d_M_v(scr_cnt,0) = (1.d0,0.d0)
            this%MoM_Chebyshev_d_M_v(scr_cnt,1:M) = c_vec((L+2):(Q+1))
        enddo

    endif

    deallocate(T,x_vec,norm_array)
    deallocate(CH_excit,c_vec,Ch_matrix_inverted,c_vec_temp,d_vec_temp)
end subroutine Eval_MoM_chebyshev_expansion_coefficeints


subroutine Eval_chebyshev_expansion_coefficeints(this,Q,L,fl_f0,fh_f0)
    type(Scatterer) :: this
    real(8) :: fl_f0,fh_f0
    integer :: Q !! it should be the evaluated taylor series coefficients max. order
    integer :: M,L
    complex*16,allocatable,dimension(:,:) :: Ch_matrix,ch_matrix_inverted
    real(8),allocatable,dimension(:,:) :: T
    real(8),allocatable,dimension(:) ::  norm_array
    complex*16,allocatable,dimension(:) :: CH_excit,c_vec,x_vec
    integer :: ii,nn,scr_cnt,cmp_cnt,qq,mm
    real(8) :: d_Kappa_by_dk,Kappa,d_Kappa_by_dk_current

    !    complex*16 :: matrix_det

    M = Q-L
    if(L > Q) then
        CALL PERROR('Bad Choice for L and Q in The Chebyshev Wideband Representation of the Solution')
        stop
    endif


    allocate(norm_array(0:Q))
    allocate(this%Chebyshev_c(1:this%N,0:L),this%Chebyshev_d(1:this%N,0:M))
    allocate(c_vec(1:(1+Q)),CH_excit(Q+1),Ch_matrix(1+Q,1+Q),Ch_matrix_inverted(1+Q,1+Q))
    allocate(x_vec(0:Q))



!    write(*,*) Q,L,M
!    return

    Kappa = get_Kappa(1.d0,fh_f0,fl_f0)

    d_Kappa_by_dk = Get_dKappa_dk(1.d0,fh_f0,fl_f0)

    Ch_matrix(1,1) = 1.d0
    Ch_matrix(2:(Q+1),1) = 0.d0

!    write(*,*) Kappa
    call chebyshevT_D(Q,Kappa,T)
!    do ii=0,Q
!        write(*,*) T(ii,:)
!    enddo
!    write(*,*)
    norm_array(0) = 1.d0
    d_Kappa_by_dk_current = d_Kappa_by_dk
    do qq=1,Q
        T(qq,:) = T(qq,:)*d_Kappa_by_dk_current
        norm_array(qq) = 1.d0/d_Kappa_by_dk_current
        d_Kappa_by_dk_current = d_Kappa_by_dk_current*d_Kappa_by_dk
    enddo

!    do ii=0,Q
!        write(*,*) T(ii,:)
!    enddo

    Ch_matrix = 0.d0

    Ch_matrix(1:(L+1),1:(L+1)) = T(0:L,0:L)
    do qq=1,L
        Ch_matrix(qq+1,1:(L+1)) = Ch_matrix(qq+1,1:(L+1))*norm_array(qq)
    enddo

    do scr_cnt = 1,this%N
        do cmp_cnt = 1,3
            CH_excit = this%I_sources_total(scr_cnt,:)%v(cmp_cnt)
            x_vec = this%I_sources_total(scr_cnt,:)%v(cmp_cnt)

!            dble(n_C_i_matrix(nn,ii))
            do mm = 1,M
                do qq = 0,Q
                    Ch_matrix(qq+1,L+1+mm) = 0.d0
                    do ii = 0,qq
                        Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm) -  &
                        dble(n_C_i_matrix(qq,ii))*x_vec(qq-ii)*T(ii,mm)
                    enddo
                    Ch_matrix(qq+1,L+1+mm) = Ch_matrix(qq+1,L+1+mm)*norm_array(qq)
                enddo
            enddo
            CH_excit = CH_excit*norm_array

            c_vec = Gauss_Elimination(Ch_matrix,CH_excit,Q+1)

            this%Chebyshev_c(scr_cnt,0:L)%v(cmp_cnt) = c_vec(1:(L+1))
            this%Chebyshev_d(scr_cnt,0)%v(cmp_cnt)  = (1.d0,0.d0)
            this%Chebyshev_d(scr_cnt,1:M)%v(cmp_cnt)  = c_vec((L+2):(Q+1))

!            write(*,*) det(Ch_matrix)
!            do ii= 1,Q+1
!                write(*,*) Ch_matrix(ii,:)
!            enddo
!            stop
        enddo
    enddo

!    do scr_cnt = 1,this%N
!        write(*,*) this%Chebyshev_c(scr_cnt,:)%v(1)
!
!    enddo
!    stop
!    write(*,*) this%Chebyshev_c(1,:)%v(1)
!    write(*,*) this%Chebyshev_d(1,:)%v(1)


    deallocate(T,x_vec,norm_array)
    deallocate(CH_excit,c_vec,Ch_matrix_inverted)
end subroutine Eval_chebyshev_expansion_coefficeints

subroutine plot_solution_bandwidth_Chebyshev(this,N_order,fl_f0,fh_f0,nl,nh)
    type(Scatterer) :: this
    integer :: N_order
    real(8) :: fl_f0,fh_f0
    real(8) :: f_loc
    integer :: nl,nh,ii,mm
    real(8) :: step_l,step_h
    real(8) , allocatable,dimension(:) :: freq_array
    real(8) , allocatable,dimension(:) :: error_array

    open(19,FILE='solution_spectrum_Chebyshev.dat')
    write(19,*) 'Freq (GHz)','                  Error'
    allocate(freq_array(nh+nl+1),error_array(nl+nh+1))
    step_l = frequency*(1.d0-fl_f0)/dble(nl)
    step_h = frequency*(fh_f0-1.d0)/dble(nh)

    if(Number_of_scatterers > 1) then
        write(*,*) '-----------------------------------------------------'
        write(*,*) 'WARNING, solution bandwidth routine is only implemented for single scatterer problems'
        write(*,*) '-----------------------------------------------------'
    endif

    mm = 1
    do ii = 1,nl+1
        freq_array(mm) = frequency*fl_f0 + dble(ii-1)*step_l
        mm = mm +1
    enddo

    do ii = 1,nh
        freq_array(mm) = frequency + dble(ii)*step_h
        mm = mm +1
    enddo

!    write(*,*) freq_array

    do mm = 1,(nl+nh+1)
        f_loc = freq_array(mm)
!        write(*,*)  f_loc
        error_array(mm) = eval_error_frequency_Chebyshev(this,f_loc)
        write(19,*) f_loc/1d9,error_array(mm)
!            write(*,*) f_loc/1d9,error_array(mm)
    enddo
    close(19)
    deallocate(freq_array,error_array)
end subroutine plot_solution_bandwidth_Chebyshev


subroutine plot_solution_bandwidth_Pade(this,N_order,fl_f0,fh_f0,nl,nh)
    type(Scatterer) :: this
    integer :: N_order
    real(8) :: fl_f0,fh_f0
    real(8) :: f_loc
    integer :: nl,nh,ii,mm
    real(8) :: step_l,step_h
    real(8) , allocatable,dimension(:) :: freq_array
    real(8) , allocatable,dimension(:) :: error_array

    integer :: N_f_Calculations
    real(8) :: BW,freq_start1,freq_end1

    open(19,FILE='solution_spectrum_Pade.dat')
    write(19,*) 'Freq (GHz)','                  Error'
    allocate(freq_array(nh+nl+1),error_array(nl+nh+1))
    step_l = frequency*(1.d0-fl_f0)/dble(nl)
    step_h = frequency*(fh_f0-1.d0)/dble(nh)

    if(Number_of_scatterers > 1) then
        write(*,*) '-----------------------------------------------------'
        write(*,*) 'WARNING, solution bandwidth routine is only implemented for single scatterer problems'
        write(*,*) '-----------------------------------------------------'
    endif

    mm = 1
    do ii = 1,nl+1
        freq_array(mm) = frequency*fl_f0 + dble(ii-1)*step_l
        mm = mm +1
    enddo

    do ii = 1,nh
        freq_array(mm) = frequency + dble(ii)*step_h
        mm = mm +1
    enddo


!    write(*,*) freq_array
!    call eval_MoM_Pade_coefficients(this,N_Taylor,Pade_L)
    call eval_Pade_coefficients(this,N_order,Pade_L)

    do mm = 1,(nl+nh+1)
        f_loc = freq_array(mm)

        call eval_error_frequency_Pade(this,f_loc,error_array(mm))
        write(19,*) f_loc/1d9,error_array(mm)
!            write(*,*) f_loc/1d9,error_array(mm)
    enddo
    close(19)

    N_f_Calculations = nl+nh+1
    freq_start1 = 0.d0
    freq_end1 = 0.d0
    do mm = 1,N_f_Calculations
        if(error_array(mm) <= BW_limit) then
            freq_start1 = freq_array(mm)
            exit
        endif
    enddo

    do mm = N_f_Calculations,1,-1
        if(error_array(mm) <= BW_limit) then
            freq_end1 = freq_array(mm)
            exit
        endif
    enddo

    BW = (freq_end1 - freq_start1)/frequency*100.d0
    write(*,*) '=============================================='
    write(*,*) 'Wideband Assessement: Pade Approximation (RAS) '
    write(*,*) 'Bandwidth: ',BW,' Percent'
    write(*,*) '=============================================='


    deallocate(freq_array,error_array)
end subroutine plot_solution_bandwidth_Pade




subroutine eval_error_frequency_RAS_AWE(this,freq1,Eu,Ev,Hu,Hv,norm_E_loc,norm_H_loc,error_t)
    type(Scatterer) :: this
    real(8) :: error_t
    complex*16,allocatable,dimension(:) :: H_u,H_v
    type(vector_c),allocatable,dimension(:) :: E_tot,H_tot
    complex*16,allocatable,dimension(:) :: BC_total
    integer:: j,kk
    complex*16,allocatable,intent(inout) :: Eu(:),Ev(:),Hu(:),Hv(:)

    real(8) :: ak0,ar0,az0,az1,freq1
    complex*16 :: ak1,ar1
    type(vector_c),allocatable,dimension(:)  :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz
    type(vector_c),allocatable,dimension(:) :: I_current
    real(8) :: error_u,error_v,norm_E_loc, norm_H_loc
    complex*16 :: numerator,denumerator,eta_ex,ar_ex,ak_ex
    integer :: ii,jj,scr_cnt,cmp_cnt

    ak0 = ak*freq1/frequency
    ar0 = ar*freq1/frequency
    az0 = az*freq1/frequency

    ak1 = this%ak_local*freq1/frequency
    ar1 = this%ar_local*freq1/frequency
    az1 = az0

    allocate(BC_total(2*this%M),E_tot(this%M),H_tot(this%M))
    allocate(I_current(this%N))
    allocate(E_Ix(0:0),E_Iy(0:0),E_Iz(0:0),H_Ix(0:0),H_Iy(0:0),H_Iz(0:0))



    do scr_cnt = 1,this%N
        do cmp_cnt = 1,3
            numerator = 0.d0
            denumerator = 1.d0
            do ii = 0,this%Pade_L
                numerator = numerator + this%Pade_a(scr_cnt,ii)%v(cmp_cnt)*((ak0-ak)*k0)**dble(ii)
            enddo
            do jj = 1,this%Pade_M
                denumerator = denumerator + this%Pade_b(scr_cnt,jj)%v(cmp_cnt)*((ak0-ak)*k0)**dble(jj)
            enddo
            I_current(scr_cnt)%v(cmp_cnt) = numerator/denumerator
        enddo
    enddo


    do j = 1,this%M
        E_tot(j)%v(:) = 0.d0
        H_tot(j)%v(:) = 0.d0
        do kk = 1,this%N
            if(this%active_region(kk)==0) then !! inside sources
                eta_ex = eta1
                ar_ex = ar0
                ak_ex = ak0
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                endif

!                write(*,*) E_Ix

                E_tot(j) = E_tot(j) + I_current(kk)%v(1)*E_Ix(0) + I_current(kk)%v(2)*E_Iy(0) +&
                I_current(kk)%v(3)*E_Iz(0)
                H_tot(j) = H_tot(j) + I_current(kk)%v(1)*H_Ix(0) + I_current(kk)%v(2)*H_Iy(0) +&
                I_current(kk)%v(3)*H_Iz(0)

            else !! outside sources
                eta_ex = this%eta_local
                ar_ex = ar1
                ak_ex = ak1
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                endif

                E_tot(j) = E_tot(j) - I_current(kk)%v(1)*E_Ix(0) - I_current(kk)%v(2)*E_Iy(0) -&
                I_current(kk)%v(3)*E_Iz(0)
                H_tot(j) = H_tot(j) - I_current(kk)%v(1)*H_Ix(0) - I_current(kk)%v(2)*H_Iy(0) -&
                I_current(kk)%v(3)*H_Iz(0)

            endif

        enddo

    enddo


    !! Evaluate the new boundary condition
    if(this%Problem_Type == 1) then !! PEC
        BC_total(1:this%M) = Eu
        BC_total((1+this%M):(2*this%M)) = Ev
        do j =1,this%M
            BC_total(j) = BC_total(j) + dot(this%tang_u(j),E_tot(j))
            BC_total(j+this%M) = BC_total(j+this%M) + dot(this%tang_v(j),E_tot(j))
        enddo



    elseif(this%Problem_Type == 2) then !! PMC
        BC_total(1:this%M) = Hu
        BC_total((1+this%M):(2*this%M)) = Hv
        do j =1,this%M
            BC_total(j) = BC_total(j)+dot(this%tang_u(j),H_tot(j))
            BC_total(j+this%M) = BC_total(j+this%M)+dot(this%tang_v(j),H_tot(j))
        enddo



    elseif(this%Problem_Type == 3) then !! IBC

        BC_total(1:this%M) = Eu +  eta1*(this%eta_uu*Hv - this%eta_uv*Hu)
        BC_total((1+this%M):(2*this%M)) = Ev +  eta1*(this%eta_vu*Hv - this%eta_vv*Hu)

        allocate(H_u(this%M),H_v(this%M))


        do j = 1,this%M
            H_u(j) = dot(this%tang_u(j),H_tot(j))
            H_v(j) = dot(this%tang_v(j),H_tot(j))
            BC_total(j) = BC_total(j)+dot(this%tang_u(j),E_tot(j))
            BC_total(this%M+j) = BC_total(this%M+j)+dot(this%tang_v(j),E_tot(j))
        enddo
        H_u = eta1*H_u
        H_v = eta1*H_v
        do j = 1,this%M
            BC_total(j) = BC_total(j) +  this%eta_uu(j)*H_v(j) - this%eta_uv(j)*H_u(j) !! E_u
            BC_total(this%M+j) = BC_total(this%M+j) +  this%eta_vu(j)*H_v(j) - this%eta_vv(j)*H_u(j) !! E_v
        enddo

        deallocate(H_u,H_v)


    elseif(this%Problem_Type == 4) then !! Dielectric
        BC_total(1:this%M) = Eu
        BC_total((1+this%M):(2*this%M)) = Ev
        do j =1,this%M
            BC_total(j) = BC_total(j) + dot(this%tang_u(j),E_tot(j))
            BC_total(j+this%M) = BC_total(j+this%M) + dot(this%tang_v(j),E_tot(j))
        enddo


    else
        write(*,*) 'ERROR: wrong problem type, it should be in the range of [1,4]'
        STOP
    endif

    error_u = 0.d0
    error_v = 0.d0
    do j=1,this%M
        error_u = error_u + abs(BC_total(j))**2.d0
        error_v = error_v + abs(BC_total(j+this%M))**2.d0

    enddo
    if(this%Problem_Type == 2) then !! PMC
        error_t = (error_v + error_u)/norm_H_loc
    else
        error_t = (error_v + error_u)/norm_E_loc
    endif


    deallocate(BC_total,E_tot,H_tot)
    deallocate(I_current)
    deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz)


end subroutine eval_error_frequency_RAS_AWE


function eval_error_frequency_Chebyshev(this,freq1) result(error_t)
    type(Scatterer) :: this

    real(8) :: error_t
    complex*16,allocatable,dimension(:) :: H_u,H_v
    type(vector_c),allocatable,dimension(:) :: E_tot,H_tot
    complex*16,allocatable,dimension(:) :: BC_total
    integer:: j,kk
    type(vector) :: Ep,En,rho_hat,test_pt
    real(8) :: ak0,ar0,az0,az1,freq1
    complex*16 :: exponential,ak1,ar1
    type(vector_c),allocatable,dimension(:)  :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz
    type(vector_c),allocatable,dimension(:) :: I_current
    real(8),allocatable,dimension(:,:) :: T
    real(8) :: error_u,error_v,norm_E_loc, norm_H_loc,Kappa
    complex*16 :: numerator,denumerator,eta_ex,ar_ex,ak_ex
    integer :: ii,jj,scr_cnt,cmp_cnt,T_order

    ak0 = ak*freq1/frequency
    ar0 = ar*freq1/frequency
    az0 = az*freq1/frequency

    ak1 = this%ak_local*freq1/frequency
    ar1 = this%ar_local*freq1/frequency
    az1 = az0

    allocate(BC_total(2*this%M),E_tot(this%M),H_tot(this%M))
    allocate(I_current(this%N))
    allocate(E_Ix(0:0),E_Iy(0:0),E_Iz(0:0),H_Ix(0:0),H_Iy(0:0),H_Iz(0:0))
    Ep = vec(-cos(theta_i)*cos(phi_i),-cos(theta_i)*sin(phi_i),sin(theta_i))
    En = vec(sin(phi_i),-cos(phi_i),0.d0)
    rho_hat = vec(sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i))
    norm_E_loc = 0.d0
    norm_H_loc =0.d0
    do j=1,this%M
        test_pt = vec(this%testing_pt(j,1),this%testing_pt(j,2),0.d0)
        exponential = exp(cj*ak0*k0*dot(rho_hat,test_pt))
        !                write(*,*) dot(rho_hat,testing_pt_MoM(j))
        E_tot(j) = eta1*exponential*(cos(alpha_i)*Ep+sin(alpha_i)*En)
        !                write(*,*) Ei(j)%v
        H_tot(j) = exponential*(sin(alpha_i)*Ep-cos(alpha_i)*En)
        norm_E_loc = norm_E_loc + abs(dot(this%tang_u(j),E_tot(j))/eta1)**2.d0 + abs(dot(this%tang_v(j),E_tot(j))/eta1)**2.d0
        norm_H_loc = norm_H_loc + abs(dot(this%tang_u(j),H_tot(j)))**2.d0 + abs(dot(this%tang_v(j),H_tot(j)))**2.d0
    enddo

!    write(*,*) norm_E_loc,norm_H_loc, this%normalize_E, this%normalize_H

        !! computation of the currents vector
!        write(*,*) this%N,this%N_curr
!        stop
!    Kappa = (2.d0*ak0 - (fh_r+fl_r))/(fh_r-fl_r)

!    Kappa = (ak0 - 1.d0)*k0

!    Kappa = (ak0 - 1.d0)/(fh_r-fl_r)
    Kappa  = Get_Kappa(ak0,fh_r,fl_r)


!    write(*,*) Kappa

!    write(*,*) this%Pade_L,this%Pade_M,size(this%Chebyshev_c)

    T_order = max(this%Pade_L,this%Pade_M)

    call ChebyshevT_D(T_order,Kappa,T)

    do scr_cnt = 1,this%N
        do cmp_cnt = 1,3
            numerator = 0.d0
            denumerator = 1.d0
            do ii = 0,this%Pade_L
                numerator = numerator + this%Chebyshev_c(scr_cnt,ii)%v(cmp_cnt)*T(0,ii)
            enddo
            do jj = 1,this%Pade_M
                denumerator = denumerator + this%Chebyshev_d(scr_cnt,jj)%v(cmp_cnt)*T(0,jj)
            enddo
            I_current(scr_cnt)%v(cmp_cnt) = numerator/denumerator
        enddo
    enddo

!    write(*,*) I_current(1:3)%v(1)
!    do scr_cnt = 1,this%N
!        write(*,*) I_current(scr_cnt)%v(:)
!    enddo
!    stop
    !        I_current(1:this%N) = this%I_sources_total(1:this%N,0)
    !
    !        I_current(1:this%N) = I_current(1:this%N)
!    write(*,*) I_current(3)%v(:)
!    stop
!    write(*,*) E_tot(3)

    do j = 1,this%M
        do kk = 1,this%N
            if(this%active_region(kk)==0) then !! inside sources
                eta_ex = eta1
                ar_ex = ar0
                ak_ex = ak0
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                endif

!                write(*,*) E_Ix

                E_tot(j) = E_tot(j) + I_current(kk)%v(1)*E_Ix(0) + I_current(kk)%v(2)*E_Iy(0) +&
                I_current(kk)%v(3)*E_Iz(0)
                H_tot(j) = H_tot(j) + I_current(kk)%v(1)*H_Ix(0) + I_current(kk)%v(2)*H_Iy(0) +&
                I_current(kk)%v(3)*H_Iz(0)

            else !! outside sources
                eta_ex = this%eta_local
                ar_ex = ar1
                ak_ex = ak1
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                endif

                E_tot(j) = E_tot(j) - I_current(kk)%v(1)*E_Ix(0) - I_current(kk)%v(2)*E_Iy(0) -&
                I_current(kk)%v(3)*E_Iz(0)
                H_tot(j) = H_tot(j) - I_current(kk)%v(1)*H_Ix(0) - I_current(kk)%v(2)*H_Iy(0) -&
                I_current(kk)%v(3)*H_Iz(0)

            endif

        enddo

    enddo

!    stop
!    write(*,*) E_tot(3)

    if(this%Problem_Type == 1) then !! PEC

!        write(*,*) abs(E_tot(1:5)%v(1)),this%tang_u(10)
        do j =1,this%M
            BC_total(j) = dot(this%tang_u(j),E_tot(j))/eta1
            BC_total(j+this%M) = dot(this%tang_v(j),E_tot(j))/eta1
        enddo


    elseif(this%Problem_Type == 2) then !! PMC

        do j =1,this%M
            BC_total(j) = dot(this%tang_u(j),H_tot(j))
            BC_total(j+this%M) = dot(this%tang_v(j),H_tot(j))
        enddo


    elseif(this%Problem_Type == 3) then !! IBC
        allocate(H_u(this%M),H_v(this%M))


        do j = 1,this%M
            H_u(j) = dot(this%tang_u(j),H_tot(j))
            H_v(j) = dot(this%tang_v(j),H_tot(j))
            BC_total(j) = dot(this%tang_u(j),E_tot(j))/eta1
            BC_total(this%M+j) = dot(this%tang_v(j),E_tot(j))/eta1
        enddo
        H_u = H_u
        H_v = H_v
        do j = 1,this%M
            BC_total(j) = BC_total(j) +  this%eta_uu(j)*H_v(j) - this%eta_uv(j)*H_u(j) !! E_u
            BC_total(this%M+j) = BC_total(this%M+j) +  this%eta_vu(j)*H_v(j) - this%eta_vv(j)*H_u(j) !! E_v
        enddo



    elseif(this%Problem_Type == 4) then !! Dielectric

        do j =1,this%M
            BC_total(j) = dot(this%tang_u(j),E_tot(j))/eta1
            BC_total(j+this%M) = dot(this%tang_v(j),E_tot(j))/eta1
        enddo


    else
        write(*,*) 'ERROR: wrong problem type, it should be in the range of [1,4]'
        STOP
    endif


    error_u = 0.d0
    error_v = 0.d0
    do j=1,this%M
        error_u = error_u + abs(BC_total(j))**2.d0
        error_v = error_v + abs(BC_total(j+this%M))**2.d0

    enddo

!    write(*,*) error_u,error_v

    if(this%Problem_Type == 2) then !! PMC
        error_t = (error_v + error_u)/norm_H_loc
    else
        error_t = (error_v + error_u)/norm_E_loc
    endif



!            write(*,*) 'norm_E=',this%normalize_E(0),'norm_E_loc =',norm_E_loc

    deallocate(BC_total,E_tot,H_tot)
    deallocate(I_current)
    deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz)


end function eval_error_frequency_Chebyshev



subroutine eval_error_frequency_Pade(this,freq1,error_t)
    type(Scatterer) :: this

    real(8) :: error_t
    complex*16,allocatable,dimension(:) :: H_u,H_v
    type(vector_c),allocatable,dimension(:) :: E_tot,H_tot
    complex*16,allocatable,dimension(:) :: BC_total
    integer:: j,kk
    type(vector) :: Ep,En,rho_hat,test_pt
    real(8) :: ak0,ar0,az0,az1,freq1
    complex*16 :: exponential,ak1,ar1
    type(vector_c),allocatable,dimension(:)  :: E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz
    type(vector_c),allocatable,dimension(:) :: I_current
    real(8) :: error_u,error_v,norm_E_loc, norm_H_loc
    complex*16 :: numerator,denumerator,eta_ex,ar_ex,ak_ex
    integer :: ii,jj,scr_cnt,cmp_cnt


    ak0 = ak*freq1/frequency
    ar0 = ar*freq1/frequency
    az0 = az*freq1/frequency

    ak1 = this%ak_local*freq1/frequency
    ar1 = this%ar_local*freq1/frequency
    az1 = az0

    allocate(BC_total(2*this%M),E_tot(this%M),H_tot(this%M))
    allocate(I_current(this%N))
    allocate(E_Ix(0:0),E_Iy(0:0),E_Iz(0:0),H_Ix(0:0),H_Iy(0:0),H_Iz(0:0))
    Ep = vec(-cos(theta_i)*cos(phi_i),-cos(theta_i)*sin(phi_i),sin(theta_i))
    En = vec(sin(phi_i),-cos(phi_i),0.d0)
    rho_hat = vec(sin(theta_i)*cos(phi_i),sin(theta_i)*sin(phi_i),cos(theta_i))
    norm_E_loc = 0.d0
    norm_H_loc =0.d0
    do j=1,this%M
        test_pt = vec(this%testing_pt(j,1),this%testing_pt(j,2),0.d0)
        exponential = exp(cj*ak0*k0*dot(rho_hat,test_pt))
        !                write(*,*) dot(rho_hat,testing_pt_MoM(j))
        E_tot(j) = eta1*exponential*(cos(alpha_i)*Ep+sin(alpha_i)*En)
        !                write(*,*) Ei(j)%v
        H_tot(j) = exponential*(sin(alpha_i)*Ep-cos(alpha_i)*En)
        norm_E_loc = norm_E_loc + abs(dot(this%tang_u(j),E_tot(j))/eta1)**2.d0 + abs(dot(this%tang_v(j),E_tot(j))/eta1)**2.d0
        norm_H_loc = norm_H_loc + abs(dot(this%tang_u(j),H_tot(j)))**2.d0 + abs(dot(this%tang_v(j),H_tot(j)))**2.d0
    enddo

!    write(*,*) norm_E_loc,norm_H_loc, this%normalize_E, this%normalize_H

        !! computation of the currents vector
!        write(*,*) this%N,this%N_curr
!        stop

    do scr_cnt = 1,this%N
        do cmp_cnt = 1,3
            numerator = 0.d0
            denumerator = 1.d0
            do ii = 0,this%Pade_L
                numerator = numerator + this%Pade_a(scr_cnt,ii)%v(cmp_cnt)*((ak0-ak)*k0)**dble(ii)
            enddo
            do jj = 1,this%Pade_M
                denumerator = denumerator + this%Pade_b(scr_cnt,jj)%v(cmp_cnt)*((ak0-ak)*k0)**dble(jj)
            enddo
            I_current(scr_cnt)%v(cmp_cnt) = numerator/denumerator
        enddo
    enddo





    do j = 1,this%M
        do kk = 1,this%N
            if(this%active_region(kk)==0) then !! inside sources
                eta_ex = eta1
                ar_ex = ar0
                ak_ex = ak0
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,eta_ex,az0,ar_ex,ak_ex,j,kk)
                endif

!                write(*,*) E_Ix

                E_tot(j) = E_tot(j) + I_current(kk)%v(1)*E_Ix(0) + I_current(kk)%v(2)*E_Iy(0) +&
                I_current(kk)%v(3)*E_Iz(0)
                H_tot(j) = H_tot(j) + I_current(kk)%v(1)*H_Ix(0) + I_current(kk)%v(2)*H_Iy(0) +&
                I_current(kk)%v(3)*H_Iz(0)

            else !! outside sources
                eta_ex = this%eta_local
                ar_ex = ar1
                ak_ex = ak1
                if(this%I_stat(kk) == 1) then
                    call eval_kernel_Electric_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                else
                    call eval_kernel_Magnetic_2D(this,this,0,E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz,&
                    eta_ex,az1,ar_ex,ak_ex,j,kk)
                endif

                E_tot(j) = E_tot(j) - I_current(kk)%v(1)*E_Ix(0) - I_current(kk)%v(2)*E_Iy(0) -&
                I_current(kk)%v(3)*E_Iz(0)
                H_tot(j) = H_tot(j) - I_current(kk)%v(1)*H_Ix(0) - I_current(kk)%v(2)*H_Iy(0) -&
                I_current(kk)%v(3)*H_Iz(0)

            endif

        enddo

    enddo

!    stop
!    write(*,*) E_tot(3)

    if(this%Problem_Type == 1) then !! PEC

!        write(*,*) abs(E_tot(1:5)%v(1)),this%tang_u(10)
        do j =1,this%M
            BC_total(j) = dot(this%tang_u(j),E_tot(j))/eta1
            BC_total(j+this%M) = dot(this%tang_v(j),E_tot(j))/eta1
        enddo


    elseif(this%Problem_Type == 2) then !! PMC

        do j =1,this%M
            BC_total(j) = dot(this%tang_u(j),H_tot(j))
            BC_total(j+this%M) = dot(this%tang_v(j),H_tot(j))
        enddo


    elseif(this%Problem_Type == 3) then !! IBC
        allocate(H_u(this%M),H_v(this%M))


        do j = 1,this%M
            H_u(j) = dot(this%tang_u(j),H_tot(j))
            H_v(j) = dot(this%tang_v(j),H_tot(j))
            BC_total(j) = dot(this%tang_u(j),E_tot(j))/eta1
            BC_total(this%M+j) = dot(this%tang_v(j),E_tot(j))/eta1
        enddo
        H_u = H_u
        H_v = H_v
        do j = 1,this%M
            BC_total(j) = BC_total(j) +  this%eta_uu(j)*H_v(j) - this%eta_uv(j)*H_u(j) !! E_u
            BC_total(this%M+j) = BC_total(this%M+j) +  this%eta_vu(j)*H_v(j) - this%eta_vv(j)*H_u(j) !! E_v
        enddo



    elseif(this%Problem_Type == 4) then !! Dielectric

        do j =1,this%M
            BC_total(j) = dot(this%tang_u(j),E_tot(j))/eta1
            BC_total(j+this%M) = dot(this%tang_v(j),E_tot(j))/eta1
        enddo


    else
        write(*,*) 'ERROR: wrong problem type, it should be in the range of [1,4]'
        STOP
    endif


    error_u = 0.d0
    error_v = 0.d0
    do j=1,this%M
        error_u = error_u + abs(BC_total(j))**2.d0
        error_v = error_v + abs(BC_total(j+this%M))**2.d0

    enddo

!    write(*,*) error_u,error_v

    if(this%Problem_Type == 2) then !! PMC
        error_t = (error_v + error_u)/norm_H_loc
    else
        error_t = (error_v + error_u)/norm_E_loc
    endif




    deallocate(BC_total,E_tot,H_tot)
    deallocate(I_current)
    deallocate(E_Ix,E_Iy,E_Iz,H_Ix,H_Iy,H_Iz)


end subroutine eval_error_frequency_Pade

subroutine eval_far_field_electric(I_vec,xp,yp,ak_l,ar_l,az_l,phi_points,E_z,E_t,H_z,H_t)
    type(vector_c) :: I_vec
    real(8) :: xp,yp,ak_l,ar_l,az_l
    real(8),allocatable,intent(inout) :: phi_points(:)
    complex*16,dimension(size(phi_points,1)) :: E_z,E_t,H_z,H_t
    complex*16,dimension(size(phi_points,1)) :: exponential

    !    exponential = ak_l*k0*eta1*exp(cj*ar_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))*sqrt(cj/(8*pi*ar_l*k0))
    !    exponential = ak_l*k0*eta1*exp(cj*ar_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))/4.d0
    !    exponential = eta1/(4.d0*ak_l*k0)*exp(cj*ar_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))
    exponential = exp(cj*ar_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))
    !    E_z = I_vec%v(3)*exponential
    !    E_t = exponential*(-I_vec%v(1)*sin(phi_points) + I_vec%v(2)*cos(phi_points) )
    E_z = -exponential/ak_l*ar_l*(ar_l*I_vec%v(3) + az_l*(I_vec%v(1)*cos(phi_points) + I_vec%v(2)*sin(phi_points)))
    E_t = exponential*ak_l*(I_vec%v(1)*sin(phi_points) - I_vec%v(2)*cos(phi_points))

    H_z = -exponential*ar_l*(-I_vec%v(1)*sin(phi_points) + I_vec%v(2)*cos(phi_points))
    H_t = exponential*(ar_l*I_vec%v(3) + az_l*(I_vec%v(1)*cos(phi_points) + I_vec%v(2)*sin(phi_points) ))

end subroutine eval_far_field_electric

subroutine eval_far_field_magnetic(I_vec,xp,yp,ak_l,ar_l,az_l,phi_points,E_z,E_t,H_z,H_t)
    type(vector_c) :: I_vec
    real(8) :: xp,yp,ak_l,ar_l,az_l
    real(8),allocatable,intent(inout) :: phi_points(:)
    complex*16,dimension(size(phi_points,1)) :: E_z,E_t,H_z,H_t
    complex*16,dimension(size(phi_points,1)) :: exponential

    !    exponential = ak_l*k0*eta1*exp(cj*ar_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))*sqrt(cj/(8*pi*ar_l*k0))
    !    exponential = ak_l*k0*eta1*exp(cj*ar_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))/4.d0
    exponential = exp(cj*ar_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))

    E_z = exponential*ar_l*(-I_vec%v(1)*sin(phi_points) + I_vec%v(2)*cos(phi_points))
    E_t = -exponential*(ar_l*I_vec%v(3) + az_l*(I_vec%v(1)*cos(phi_points) + I_vec%v(2)*sin(phi_points) ))

    H_z = -exponential/ak_l*ar_l*(ar_l*I_vec%v(3) + az_l*(I_vec%v(1)*cos(phi_points) + I_vec%v(2)*sin(phi_points)))
    H_t = exponential*ak_l*(I_vec%v(1)*sin(phi_points) - I_vec%v(2)*cos(phi_points))
!    E_z = exponential*(I_vec%v(1)*sin(phi_points) - I_vec%v(2)*cos(phi_points) )
!    E_t = I_vec%v(3)*exponential
end subroutine eval_far_field_magnetic


subroutine eval_far_field_scatterers_RAS_AWE(this,phi_points,ak_l,ar_l,az_l,E_z,E_t,I_sources_AWE)
    type(Scatterer) :: this
    real(8),allocatable,intent(inout) :: phi_points(:)
    real(8) :: xp,yp,ak_l,ar_l,az_l
    complex*16,allocatable,intent(inout) :: E_z(:),E_t(:)
    complex*16,allocatable,dimension(:) :: E_z_temp,E_t_temp,H_z_temp,H_t_temp,H_z,H_t
    type(vector_c) :: I_vec
    type(vector_c),allocatable,intent(inout) :: I_sources_AWE(:)
    integer :: ii,nn

    allocate(H_z(size(E_z,1)),H_t(size(E_z,1)))
    E_z = 0.d0
    E_t = 0.d0
    H_z = 0.d0
    H_t = 0.d0

    allocate(E_z_temp(size(E_z,1)),E_t_temp(size(E_z,1)),H_z_temp(size(E_z,1)),H_t_temp(size(E_z,1)))
    !    write(*,*) 'gowwa el RAS'

    do nn = 1,this%N
        !            write(*,*) Scat(ii)%active_region(nn)
        if(this%active_region(nn)==0) then
            I_vec = I_sources_AWE(nn)
            !                write(*,*) abs(I_vec%v)
            xp = this%Source_pos(nn,1)
            yp = this%Source_pos(nn,2)
            if(this%I_stat(nn) == 1) then !! electric source
                call eval_far_field_electric(I_vec,xp,yp,ak_l,ar_l,az_l,phi_points,E_z_temp,E_t_temp,H_z_temp,H_t_temp)
            else
                call eval_far_field_magnetic(I_vec,xp,yp,ak_l,ar_l,az_l,phi_points,E_z_temp,E_t_temp,H_z_temp,H_t_temp)
            endif
            E_z = E_z + E_z_temp
            E_t = E_t + E_t_temp
            H_z = H_z + H_z_temp
            H_t = H_t + H_t_temp
        endif

    enddo


    H_z = 0.25d0*k0*H_z
    H_t = 0.25d0*k0*H_t

    E_z = 0.25d0*eta1*k0*E_z
    E_t = 0.25d0*eta1*k0*E_t

!    E_z = (-eta1*H_t+ E_z)/2.d0
!    E_t = (eta1*H_z+E_t)/2.d0

    deallocate(E_z_temp,E_t_temp,H_z_temp,H_t_temp,H_z,H_t)
end subroutine eval_far_field_scatterers_RAS_AWE



subroutine eval_far_field_scatterers_RAS(Scat,phi_points,ak_l,ar_l,az_l,E_z,E_t)
    type(Scatterer),allocatable,intent(inout) :: Scat(:)
    real(8),allocatable,intent(inout) :: phi_points(:)
    real(8) :: xp,yp,ak_l,ar_l,az_l
    complex*16,allocatable,intent(inout) :: E_z(:),E_t(:)
    complex*16,allocatable,dimension(:) :: E_z_temp,E_t_temp,H_z_temp,H_t_temp,H_z,H_t
    type(vector_c) :: I_vec
    integer :: ii,nn

    allocate(H_z(size(E_z,1)),H_t(size(E_z,1)))
    E_z = 0.d0
    E_t = 0.d0
    H_z = 0.d0
    H_t = 0.d0

    allocate(E_z_temp(size(E_z,1)),E_t_temp(size(E_z,1)),H_z_temp(size(E_z,1)),H_t_temp(size(E_z,1)))
    !    write(*,*) 'gowwa el RAS'
    do ii = 1,number_of_scatterers
        do nn = 1,Scat(ii)%N
            !            write(*,*) Scat(ii)%active_region(nn)
            if(Scat(ii)%active_region(nn)==0) then
                I_vec = Scat(ii)%I_sources_total(nn,0)
                !                write(*,*) abs(I_vec%v)
                xp = Scat(ii)%Source_pos(nn,1)
                yp = Scat(ii)%Source_pos(nn,2)
                if(Scat(ii)%I_stat(nn) == 1) then !! electric source
                    call eval_far_field_electric(I_vec,xp,yp,ak_l,ar_l,az_l,phi_points,E_z_temp,E_t_temp,H_z_temp,H_t_temp)
                else
                    call eval_far_field_magnetic(I_vec,xp,yp,ak_l,ar_l,az_l,phi_points,E_z_temp,E_t_temp,H_z_temp,H_t_temp)
                endif
                E_z = E_z + E_z_temp
                E_t = E_t + E_t_temp
                H_z = H_z + H_z_temp
                H_t = H_t + H_t_temp
            endif

        enddo
    enddo

    if(Source_Model == 1) then  !! line source excitation
    !! include far-fields of the sources
        do ii=1,N_Line_sources
            if(Line_Sources(ii)%Source_Orientation == 1) then !1 z_oriented
                I_vec = vec((0.d0,0.d0),(0.d0,0.d0),Line_Sources(ii)%Amp)
            elseif(Line_Sources(ii)%Source_Orientation == 2)then !! x_oriented
                I_vec = vec(Line_Sources(ii)%Amp,(0.d0,0.d0),(0.d0,0.d0))
            else
                I_vec = vec((0.d0,0.d0),Line_Sources(ii)%Amp,(0.d0,0.d0))
            endif
            xp = Line_Sources(ii)%x_s
            yp = Line_Sources(ii)%y_s
            if(Line_Sources(ii)%Source_Type == 1) then !! electric source
                call eval_far_field_electric(I_vec,xp,yp,ak_l,ar_l,az_l,phi_points,E_z_temp,E_t_temp,H_z_temp,H_t_temp)
            else
                call eval_far_field_magnetic(I_vec,xp,yp,ak_l,ar_l,az_l,phi_points,E_z_temp,E_t_temp,H_z_temp,H_t_temp)
            endif
            E_z = E_z + E_z_temp
            E_t = E_t + E_t_temp
            H_z = H_z + H_z_temp
            H_t = H_t + H_t_temp
        enddo
    endif


    H_z = 0.25d0*k0*H_z
    H_t = 0.25d0*k0*H_t

    E_z = 0.25d0*eta1*k0*E_z
    E_t = 0.25d0*eta1*k0*E_t

!    E_z = (-eta1*H_t+ E_z)/2.d0
!    E_t = (eta1*H_z+E_t)/2.d0

    deallocate(E_z_temp,E_t_temp,H_z_temp,H_t_temp,H_z,H_t)
end subroutine eval_far_field_scatterers_RAS


subroutine eval_far_fields_scatterer_MoM_Wideband(Scat,phi_points,ak_l,ar_l,az_l,E_z,E_t,I_v,I_u,M_v,M_u)
    type(Scatterer) :: Scat
    real(8),allocatable,intent(inout) :: phi_points(:)
    complex*16,allocatable,intent(inout) :: I_v(:),I_u(:),M_v(:),M_u(:)
    real(8) :: ak_l,ar_l,xp,yp,az_l
    complex*16,allocatable,intent(inout) :: E_z(:),E_t(:)
    complex*16,dimension(size(phi_points,1)) :: exponential
    type(vector),dimension(size(phi_points,1)) :: phi_hat,rho_hat
    real(8),dimension(size(phi_points,1))  :: r_zero
    integer :: ss,ii,M

    M = size(phi_points,1)
    r_zero = 0.d0
    phi_hat = vec(-sin(phi_points),cos(phi_points),r_zero,M)
    rho_hat = vec(cos(phi_points),sin(phi_points),r_zero,M)

    E_z = 0.d0
    E_t = E_z



    do ii = 1,Scat%M
        if(Scat%testing_pt_status(ii,1) == -1) then
            cycle !! cancelled point
        endif
        if(Scat%testing_pt_status(ii,3) /= 0) then
            cycle !! the exterior point is not the free-space region
        endif
        xp = Scat%testing_pt_MoM(ii)%v(1)
        yp = Scat%testing_pt_MoM(ii)%v(2)

        exponential = Scat%delta_n(ii)*exp(cj*ar_l*ak_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))

        E_z = E_z - exponential*ar_l*ak_l*k0/4.d0*(eta1*(ar_l*I_v(ii) + &
        az_l*I_u(ii)*dot(Scat%tang_u(ii),rho_hat,M)) - M_u(ii)*&
        dot(Scat%tang_u(ii),phi_hat,M))

        E_t = E_t - exponential*ak_l*k0/4.d0*(eta1*I_u(ii)*dot(Scat%tang_u(ii),phi_hat,M) +&
        (ar_l*M_v(ii) - az_l*M_u(ii)*dot(Scat%norm_v(ii),phi_hat,M)))
    enddo


    do ii=1,N_Line_Sources
        xp = Line_Sources(ii)%x_s
        yp = Line_Sources(ii)%y_s
        exponential = exp(cj*ar_l*ak_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))
        if(Line_Sources(ii)%Source_Type == 1) then !! Electric
            if(Line_Sources(ii)%Source_Orientation == 1) then !! z-oriented
                E_z = E_z - exponential*ak_l*k0*(eta1/4.d0*ar_l*ar_l*Line_Sources(ii)%Amp)
            elseif(Line_Sources(ii)%Source_Orientation == 2) then !! x-oriented
                E_z = E_z - exponential*ak_l*ar_l*k0*(eta1/4.d0*(az_l*Line_Sources(ii)%Amp*cos(phi_points)) )

                E_t = E_t - exponential*(eta1*ak_l*k0/4.d0*Line_Sources(ii)%Amp*(-sin(phi_points)) )
            else !! y-oriented
                E_z = E_z - exponential*ak_l*ar_l*k0*(eta1/4.d0*(az_l*Line_Sources(ii)%Amp*sin(phi_points)) )

                E_t = E_t - exponential*(eta1*ak_l*k0/4.d0*Line_Sources(ii)%Amp*(cos(phi_points)) )
            endif

        else !! Magnetic
            exponential = exponential*eta1
            if(Line_Sources(ii)%Source_Orientation == 1) then !! z-oriented
!                write(*,*) 'MoM-Magnetic Line Source'
                E_t = E_t - exponential*0.25d0*k0*ar_l*Line_Sources(ii)%Amp
            elseif(Line_Sources(ii)%Source_Orientation == 2) then !! x-oriented
                E_z = E_z - exponential*ak_l*k0*0.25d0*Line_Sources(ii)%Amp*(sin(phi_points))

                E_t = E_t - exponential*0.25d0*k0*az_l*ak_l*Line_Sources(ii)%Amp*(cos(phi_points))
            else !! y-oriented
                E_z = E_z - exponential*ak_l*k0*0.25d0*Line_Sources(ii)%Amp*(-cos(phi_points))

                E_t = E_t - exponential*0.25d0*k0*az_l*ak_l*Line_Sources(ii)%Amp*(sin(phi_points))

!                E_z = E_z - exponential*ak_l*k0*(eta1/(4.d0*ak_l)*(ar_l*Scat(ss)%I_v(ii,0) + &
!                az_l*Scat(ss)%I_u(ii,0)*dot(Scat(ss)%tang_u(ii),rho_hat,M)) - 0.25d0*Scat(ss)%M_u(ii,0)*&
!                dot(Scat(ss)%tang_u(ii),phi_hat,M))
!
!                E_t = E_t - exponential*(eta1*ak_l*k0/4.d0*Scat(ss)%I_u(ii,0)*dot(Scat(ss)%tang_u(ii),phi_hat,M) +&
!                0.25d0*k0*(ar_l*Scat(ss)%M_v(ii,0) - az_l*Scat(ss)%M_u(ii,0)*dot(Scat(ss)%norm_v(ii),phi_hat,M)))
            endif
        endif
    enddo


end subroutine eval_far_fields_scatterer_MoM_Wideband


subroutine eval_far_fields_scatterer_MoM(Scat,phi_points,ak_l,ar_l,az_l,E_z,E_t)
    type(Scatterer),allocatable,intent(inout) :: Scat(:)
    real(8),allocatable,intent(inout) :: phi_points(:)
    real(8) :: ak_l,ar_l,xp,yp,az_l
    complex*16,allocatable,intent(inout) :: E_z(:),E_t(:)
    complex*16,dimension(size(phi_points,1)) :: exponential
    type(vector),dimension(size(phi_points,1)) :: phi_hat,rho_hat
    real(8),dimension(size(phi_points,1))  :: r_zero
    integer :: ss,ii,M

    M = size(phi_points,1)
    r_zero = 0.d0
    phi_hat = vec(-sin(phi_points),cos(phi_points),r_zero,M)
    rho_hat = vec(cos(phi_points),sin(phi_points),r_zero,M)

    E_z = 0.d0
    E_t = E_z
    do ss = 1,number_of_scatterers
        do ii = 1,Scat(ss)%M
            if(Scat(ss)%testing_pt_status(ii,1) == -1) then
                cycle !! cancelled point
            endif
            if(Scat(ss)%testing_pt_status(ii,3) /= 0) then
                cycle !! the exterior point is not the free-space region
            endif
            xp = Scat(ss)%testing_pt_MoM(ii)%v(1)
            yp = Scat(ss)%testing_pt_MoM(ii)%v(2)
            !            exponential = Scat(ss)%delta_n(ii)*ak_l*k0*sqrt(cj/(8.d0*pi*ar_l*k0))*&
            !            exponential = Scat(ss)%delta_n(ii)*ak_l*k0/4.d0*&
            !            exp(cj*ar_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))
            !
            !            E_z = E_z + exponential*(Scat(ss)%I_v(ii)*eta1 - dot(Scat(ss)%tang_u(ii),phi_hat,M)*Scat(ss)%M_u(ii))
            !            E_t = E_t + exponential*(Scat(ss)%M_v(ii) + dot(Scat(ss)%tang_u(ii),phi_hat,M)*eta1*Scat(ss)%I_u(ii))
            exponential = Scat(ss)%delta_n(ii)*exp(cj*ar_l*ak_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))

            E_z = E_z - exponential*ar_l*ak_l*k0/4.d0*(eta1*(ar_l*Scat(ss)%I_v(ii,0) + &
            az_l*Scat(ss)%I_u(ii,0)*dot(Scat(ss)%tang_u(ii),rho_hat,M)) - Scat(ss)%M_u(ii,0)*&
            dot(Scat(ss)%tang_u(ii),phi_hat,M))

            E_t = E_t - exponential*ak_l*k0/4.d0*(eta1*Scat(ss)%I_u(ii,0)*dot(Scat(ss)%tang_u(ii),phi_hat,M) +&
            (ar_l*Scat(ss)%M_v(ii,0) - az_l*Scat(ss)%M_u(ii,0)*dot(Scat(ss)%norm_v(ii),phi_hat,M)))
        enddo
    enddo

    do ii=1,N_Line_Sources
        xp = Line_Sources(ii)%x_s
        yp = Line_Sources(ii)%y_s
        exponential = exp(cj*ar_l*ak_l*k0*(xp*cos(phi_points) + yp*sin(phi_points)))
        if(Line_Sources(ii)%Source_Type == 1) then !! Electric
            if(Line_Sources(ii)%Source_Orientation == 1) then !! z-oriented
                E_z = E_z - exponential*ak_l*k0*(eta1/4.d0*ar_l*ar_l*Line_Sources(ii)%Amp)
            elseif(Line_Sources(ii)%Source_Orientation == 2) then !! x-oriented
                E_z = E_z - exponential*ak_l*ar_l*k0*(eta1/4.d0*(az_l*Line_Sources(ii)%Amp*cos(phi_points)) )

                E_t = E_t - exponential*(eta1*ak_l*k0/4.d0*Line_Sources(ii)%Amp*(-sin(phi_points)) )
            else !! y-oriented
                E_z = E_z - exponential*ak_l*ar_l*k0*(eta1/4.d0*(az_l*Line_Sources(ii)%Amp*sin(phi_points)) )

                E_t = E_t - exponential*(eta1*ak_l*k0/4.d0*Line_Sources(ii)%Amp*(cos(phi_points)) )
            endif

        else !! Magnetic
            exponential = exponential*eta1
            if(Line_Sources(ii)%Source_Orientation == 1) then !! z-oriented
!                write(*,*) 'MoM-Magnetic Line Source'
                E_t = E_t - exponential*0.25d0*k0*ar_l*Line_Sources(ii)%Amp
            elseif(Line_Sources(ii)%Source_Orientation == 2) then !! x-oriented
                E_z = E_z - exponential*ak_l*k0*0.25d0*Line_Sources(ii)%Amp*(sin(phi_points))

                E_t = E_t - exponential*0.25d0*k0*az_l*ak_l*Line_Sources(ii)%Amp*(cos(phi_points))
            else !! y-oriented
                E_z = E_z - exponential*ak_l*k0*0.25d0*Line_Sources(ii)%Amp*(-cos(phi_points))

                E_t = E_t - exponential*0.25d0*k0*az_l*ak_l*Line_Sources(ii)%Amp*(sin(phi_points))

!                E_z = E_z - exponential*ak_l*k0*(eta1/(4.d0*ak_l)*(ar_l*Scat(ss)%I_v(ii,0) + &
!                az_l*Scat(ss)%I_u(ii,0)*dot(Scat(ss)%tang_u(ii),rho_hat,M)) - 0.25d0*Scat(ss)%M_u(ii,0)*&
!                dot(Scat(ss)%tang_u(ii),phi_hat,M))
!
!                E_t = E_t - exponential*(eta1*ak_l*k0/4.d0*Scat(ss)%I_u(ii,0)*dot(Scat(ss)%tang_u(ii),phi_hat,M) +&
!                0.25d0*k0*(ar_l*Scat(ss)%M_v(ii,0) - az_l*Scat(ss)%M_u(ii,0)*dot(Scat(ss)%norm_v(ii),phi_hat,M)))
            endif
        endif
    enddo


end subroutine eval_far_fields_scatterer_MoM

subroutine eval_monostatic_RCS_comparison(Scat)
    type(Scatterer),allocatable,intent(inout) :: Scat(:)
    real(8),allocatable,dimension(:) :: phi_points
    complex*16,allocatable,dimension(:)  ::  E_z_RAS,E_t_RAS,E_z_MoM,E_t_MoM
    integer :: fd1 = 86

    allocate(E_z_RAS(1),E_t_RAS(1),phi_points(1),E_z_MoM(1),E_t_MoM(1))

    phi_points = phi_i

    call eval_far_fields_scatterer_MoM(Scat,phi_points,ak,ar,az,E_z_MoM,E_t_MoM)
    call eval_far_field_scatterers_RAS(Scat,phi_points,ak,ar,az,E_z_RAS,E_t_RAS)

    OPEN(fd1, FILE='Monostatic_RCS_Single_Freq_Comparison_RAS_MOM.dat')
        write(fd1,*) 'Frequency,        abs(E_z_RAS),       abs(E_t_RAS),        abs(E_z_MoM),       abs(E_t_MoM)'
        write(fd1,*) Frequency,abs(E_z_RAS),abs(E_t_RAS),abs(E_z_MoM),abs(E_t_MoM)
    close(fd1)
    deallocate(E_z_RAS,E_t_RAS,phi_points,E_z_MoM,E_t_MoM)
end subroutine eval_monostatic_RCS_comparison

    ! -----------------------------------------------------------------------
    ! Subroutine: compare_far_field_scatterers
    ! Purpose   : Evaluates and writes the bistatic far-field RCS pattern
    !             for M_phi azimuth angles from both RAS and MoM solutions.
    !             Outputs to 'far_field_comparison_RAS.dat' and
    !             'far_field_comparison_MoM.dat'.
    ! -----------------------------------------------------------------------
subroutine compare_far_field_scatterers(Scat,M_phi)
    type(Scatterer),allocatable,intent(inout) :: Scat(:)
    integer :: M_phi,mm
    real(8),allocatable,dimension(:) :: phi_points
    complex*16,allocatable,dimension(:) :: E_z_RAS,E_t_RAS,E_z_MoM,E_t_MoM
    real(8),allocatable,dimension(:) :: phase_E_z_RAS,phase_E_t_RAS
    real(8),allocatable,dimension(:) :: phase_E_z_MoM,phase_E_t_MoM
    real(8) :: delta_phi
    integer :: fd1 = 86

    allocate(phi_points(M_phi),E_z_RAS(M_phi),E_t_RAS(M_phi),E_z_MoM(M_phi),E_t_MoM(M_phi))
    allocate(phase_E_z_RAS(M_phi),phase_E_t_RAS(M_phi),phase_E_z_MoM(M_phi),phase_E_t_MoM(M_phi))
    phi_points = 0.d0
    delta_phi = tpi/(M_phi-1)
    do mm = 1,(M_phi-1)
        phi_points(mm+1) = mm*delta_phi
    enddo
    !    write(*,*) 'ana hena ho barsem el far fields'
    call eval_far_field_scatterers_RAS(Scat,phi_points,ak,ar,az,E_z_RAS,E_t_RAS)
    call eval_far_fields_scatterer_MoM(Scat,phi_points,ak,ar,az,E_z_MoM,E_t_MoM)

    phase_E_z_RAS = atan2(aimag(E_z_RAS),real(E_z_RAS))*180.d0/pi
    phase_E_t_RAS = atan2(aimag(E_t_RAS),real(E_t_RAS))*180.d0/pi
    phase_E_z_MoM = atan2(aimag(E_z_MoM),real(E_z_MoM))*180.d0/pi
    phase_E_t_MoM = atan2(aimag(E_t_MoM),real(E_t_MoM))*180.d0/pi

    OPEN(fd1, FILE='far_field_comparison.dat')
    do mm = 1,M_phi
        write(fd1,*) phi_points(mm)*180.d0/pi,abs(E_z_RAS(mm)),abs(E_z_MoM(mm)),abs(E_t_RAS(mm)),abs(E_t_MoM(mm)),&
        phase_E_z_RAS(mm),phase_E_z_MoM(mm),phase_E_t_RAS(mm),phase_E_t_MoM(mm)
    enddo

    close(fd1)
    deallocate(phi_points,E_z_RAS,E_t_RAS,E_z_MoM,E_t_MoM)
    deallocate(phase_E_z_RAS,phase_E_t_RAS,phase_E_z_MoM,phase_E_t_MoM)
end subroutine compare_far_field_scatterers


subroutine plot_far_field_MoM(Scat,M_phi)
    type(Scatterer),allocatable,intent(inout) :: Scat(:)
    integer :: M_phi,mm
    real(8),allocatable,dimension(:) :: phi_points
    complex*16,allocatable,dimension(:) :: E_z_MoM,E_t_MoM
    real(8),allocatable,dimension(:) :: phase_E_z_MoM,phase_E_t_MoM
    real(8) :: delta_phi
    integer :: fd1 = 86

    allocate(phi_points(M_phi),E_z_MoM(M_phi),E_t_MoM(M_phi))
    allocate(phase_E_z_MoM(M_phi),phase_E_t_MoM(M_phi))
    phi_points = 0.d0
    delta_phi = tpi/(M_phi-1)
    do mm = 1,(M_phi-1)
        phi_points(mm+1) = mm*delta_phi
    enddo
    call eval_far_fields_scatterer_MoM(Scat,phi_points,ak,ar,az,E_z_MoM,E_t_MoM)

    phase_E_z_MoM = atan2(aimag(E_z_MoM),real(E_z_MoM))*180.d0/pi
    phase_E_t_MoM = atan2(aimag(E_t_MoM),real(E_t_MoM))*180.d0/pi

    OPEN(fd1, FILE='far_field_comparison_MoM.dat')
    do mm = 1,M_phi
        write(fd1,*) phi_points(mm)*180.d0/pi,abs(E_z_MoM(mm)),abs(E_t_MoM(mm)),&
        phase_E_z_MoM(mm),phase_E_t_MoM(mm)
    enddo

    close(fd1)
    deallocate(phi_points,E_z_MoM,E_t_MoM)
    deallocate(phase_E_z_MoM,phase_E_t_MoM)
end subroutine plot_far_field_MoM

subroutine plot_far_field_RAS(Scat,M_phi)
    type(Scatterer),allocatable,intent(inout) :: Scat(:)
    integer :: M_phi,mm
    real(8),allocatable,dimension(:) :: phi_points
    complex*16,allocatable,dimension(:) :: E_z_RAS,E_t_RAS
    real(8),allocatable,dimension(:) :: phase_E_z_RAS,phase_E_t_RAS
    real(8) :: delta_phi
    integer :: fd1 = 86

    allocate(phi_points(M_phi),E_z_RAS(M_phi),E_t_RAS(M_phi))
    allocate(phase_E_z_RAS(M_phi),phase_E_t_RAS(M_phi))
    phi_points = 0.d0
    delta_phi = tpi/(M_phi-1)
    do mm = 1,(M_phi-1)
        phi_points(mm+1) = mm*delta_phi
    enddo
    !    write(*,*) 'ana hena ho barsem el far fields'
    call eval_far_field_scatterers_RAS(Scat,phi_points,ak,ar,az,E_z_RAS,E_t_RAS)

    phase_E_z_RAS = atan2(aimag(E_z_RAS),real(E_z_RAS))*180.d0/pi
    phase_E_t_RAS = atan2(aimag(E_t_RAS),real(E_t_RAS))*180.d0/pi

    OPEN(fd1, FILE='far_field_comparison_RAS.dat')
    do mm = 1,M_phi
        write(fd1,*) phi_points(mm)*180.d0/pi,abs(E_z_RAS(mm)),abs(E_t_RAS(mm)),&
        phase_E_z_RAS(mm),phase_E_t_RAS(mm)
    enddo

    close(fd1)
    deallocate(phi_points,E_z_RAS,E_t_RAS)
    deallocate(phase_E_z_RAS,phase_E_t_RAS)
end subroutine plot_far_field_RAS

end module scatterer_mod







module RSM_3D_sph
    use sim_par
    use Operations
    !    use Gauss_Reduction
    use constants
    use discretizer_3D
    use scatterer_mod
    implicit none

contains
    subroutine set_scatterers_IDs(Scats)
        type(Scatterer),allocatable,intent(inout) :: Scats(:)
        integer :: i,n
        complex*16:: eps_r,mu_r

        allocate(Scats(number_of_scatterers))



        write(*,*) '======================================================================='
        write(*,*) '===================Scatterers Phisical Parameters======================'
        do i=1,number_of_scatterers
            Scats(i)%region_ID = i
            if(material_ID_read(i) == 0) then !! free-space (background material)
                Scats(i)%Problem_Type = 4

                Scats(i)%ak_local = 1.d0
                Scats(i)%k_local = k0
                Scats(i)%ar_local = sqrt(1.d0 - az**2.d0)
                Scats(i)%lambda_local = tpi/real(Scats(i)%k_local)
                Scats(i)%eta_local = eta0
            elseif(material_ID_read(i) == 1) then !! PEC
                Scats(i)%Problem_Type = 1

                Scats(i)%ak_local = 0.d0
                Scats(i)%k_local = 0.d0
                Scats(i)%ar_local = 0.d0
                Scats(i)%lambda_local = 0.d0
                Scats(i)%eta_local = 0.d0
            elseif(material_ID_read(i) == -1) then !! PMC
                Scats(i)%Problem_Type = 2

                Scats(i)%ak_local = 0.d0
                Scats(i)%k_local = 0.d0
                Scats(i)%ar_local = 0.d0
                Scats(i)%lambda_local = 0.d0
                Scats(i)%eta_local = eta0*100.d0

            elseif(material_ID_read(i) > 1) then !! Dielectric
                Scats(i)%Problem_Type = 4
                do n=1,materials%N_diel
                    if(materials%ID_diel(n) == material_ID_read(i)) then
                        exit
                    endif
                enddo
                eps_r = materials%eps_r_p(n) + cj*materials%eps_r_pp(n)
                mu_r = materials%mu_r_p(n) + cj*materials%mu_r_pp(n)
                Scats(i)%ak_local = sqrt(eps_r)*sqrt(mu_r)
                Scats(i)%k_local = k0*Scats(i)%ak_local
                !                Scats(i)%ar_local = sqrt(materials%eps_r_p(n)*materials%mu_r_p(n) - az**2.d0)
                Scats(i)%ar_local = sqrt(Scats(i)%ak_local - az)*sqrt(Scats(i)%ak_local + az)
                Scats(i)%lambda_local = tpi/abs(real(Scats(i)%k_local))
                Scats(i)%eta_local = eta0*sqrt(mu_r)/sqrt(eps_r)
            !                write(*,*) materials%mu_r_p(n),materials%eps_r_p(n)
            !                write(*,*) Scats(i)%k_local,Scats(i)%eta_local,Scats(i)%ar_local,Scats(i)%lambda_local
            !                stop
            elseif(material_ID_read(i) < -1) then !! IBC
                Scats(i)%Problem_Type = 3
                do n=1,materials%N_ibc
                    if(materials%ID_ibc(n) == material_ID_read(i)) then
                        exit
                    endif
                enddo
                Scats(i)%eta_zz = materials%eta_zz(n)
                Scats(i)%eta_zt = materials%eta_zt(n)
                Scats(i)%eta_tz = materials%eta_tz(n)
                Scats(i)%eta_tt = materials%eta_tt(n)

                Scats(i)%ak_local = ak
                Scats(i)%k_local = k0
                Scats(i)%ar_local = ar
                Scats(i)%lambda_local = lambda
                Scats(i)%eta_local = eta1


            endif
            !            Scats(i)%Problem_Type = Problem_Type_array(i)
            !            Scats(i)%k_local = k_array(i+1)
            !            Scats(i)%ak_local = Scats(i)%k_local/k0
            !            Scats(i)%ar_local = sqrt(er_array(i+1)*mur_array(i+1) - az**2.d0)
            !            write(*,*) Scats(i)%ar_local, Scats(i)%k_local/k0
            !            Scats(i)%lambda_local = tpi/Scats(i)%k_local
            !            Scats(i)%eta_local = eta_array(i+1)
            !            write(*,*) 'Scatterer #',i
            if(Scats(i)%Problem_Type == 1) then
                write(*,*) 'Scatterer #',i,'-> B.C. Type: PEC'
            elseif(Scats(i)%Problem_Type == 2) then
                write(*,*) 'Scatterer #',i,'-> B.C. Type: PMC'
            elseif(Scats(i)%Problem_Type == 3) then
                write(*,*) 'Scatterer #',i,'-> B.C. Type: IBC'
            elseif(Scats(i)%Problem_Type == 4) then
                write(*,*) 'Scatterer #',i,'-> B.C. Type: Dielectric'
            else
                write(*,*) 'ERROR: problem type not supported'
                Stop
            endif


            if(scatterers_input_method(i) == 1) then
                Scats(i)%Center = (/physical_parameters(i,1),physical_parameters(i,2),0.d0/)
                Scats(i)%a_superquad = physical_parameters(i,3)
                Scats(i)%b_superquad = physical_parameters(i,4)
                Scats(i)%g_superquad = physical_parameters(i,5)
                write(*,*) 'Center (',Scats(i)%Center,')'
                write(*,*) 'Superquadratic Par a,b,g =',Scats(i)%a_superquad,Scats(i)%b_superquad,Scats(i)%g_superquad
            else
                write(*,*) 'Scatterer #',i,' is read from file: ',scatterer_input_file_names(i)
                Scats(i)%Center = (/0.d0,0.d0,0.d0/)
                Scats(i)%a_superquad = 0.d0
                Scats(i)%b_superquad = 0.d0
                Scats(i)%g_superquad = 2


            endif




            !            write(*,*) 'Medium Parameters k,eta',Scats(i)%k_local,Scats(i)%eta_local
            write(*,*) '_______________________________________________________________________'

        enddo

        call segment_scatterers(Scats)

        if(Contour_Sources_Type_read == 4) then !! Auto set for the contour sources type for best convergence
            do i =1,number_of_scatterers
                if(Scats(i)%Problem_type == 1) then !! PEC Boundary
                    Scats(i)%Contour_Sources_Type = 1
                elseif(Scats(i)%Problem_Type == 2) then !! PMC Boundary
                    Scats(i)%Contour_Sources_Type = 2
                elseif(Scats(i)%Problem_Type == 3) then !! IBC Boundary
                    Scats(i)%Contour_Sources_Type = 3
                else !! Dielectric
                    Scats(i)%Contour_Sources_Type = 1
                endif
            enddo
        else
            do i =1,number_of_scatterers
                Scats(i)%Contour_Sources_Type = Contour_Sources_Type_read
            enddo
        endif

        if(Wideband_type > 0) then
            allocate(n_C_i_matrix(0:N_taylor,0:N_taylor),n_P_i_matrix(0:N_taylor,0:N_taylor))
            allocate(n_X_i_matrix(0:N_taylor,0:N_taylor))
            n_C_i_matrix = 0
            n_P_i_matrix = 0
            n_X_i_matrix = 0
            do n =0,N_taylor
                n_C_i_matrix(n,0) = 1
                n_P_i_matrix(n,0) = 1

                do i = 1,n
                    n_C_i_matrix(n,i) = (n_C_i_matrix(n,i-1)*(n-i+1))/i
                    n_P_i_matrix(n,i) = n_P_i_matrix(n,i-1)*(n-i+1)
                enddo
                n_X_i_matrix(n,n) = 1
                do i = n-1,1,-1
                    n_X_i_matrix(n,i) = (i+1)*n_X_i_matrix(n,i+1)
                enddo
            !                write(*,*) n,'n_P_i', n_P_i_matrix(n,0:n)
            !                 write(*,*) n,'n_C_i', n_C_i_matrix(n,0:n)
            enddo
        else
            N_taylor = 0
            allocate(n_C_i_matrix(0:0,0:0),n_P_i_matrix(0:0,0:0),n_X_i_matrix(0:0,0:0))
            n_C_i_matrix = 1
            n_P_i_matrix = 1
            n_X_i_matrix = 1
        endif

    end subroutine set_scatterers_IDs

    ! -----------------------------------------------------------------------
    ! Subroutine: run_main_program  (module RSM_3D_sph, entry point)
    ! Purpose   : Top-level driver. Sequence:
    !   1. read_parameters() -- load config from Parameters.dat
    !   2. Allocate Scatterer array; call set_scatterers_IDs
    !   3. segment_scatterers -- discretise all boundaries
    !   4. Loop over SPW/N_group parameter sweeps:
    !      a. initialize_parameters per scatterer
    !      b. Evaluate incident field (plane wave or line sources)
    !      c. solve_scattering_multiscatterer* (RAS)
    !      d. eval_surface_current_MoM_once (MoM, if enabled)
    !      e. Post-processing: far field, RCS, wideband, colour map
    !   5. Write output files and deallocate
    ! -----------------------------------------------------------------------
    subroutine run_main_program()
        real :: start_t, finish_t,start_t_MoM, finish_t_MoM
        real(8) :: length,R_probs
        integer :: cnt,i
        integer :: cnt_N,cnt_S
        integer :: cnt_Itr,N_Itr_Total
        type(Scatterer),allocatable,dimension(:) :: Scat
        type(Scatterer) :: Scat_big
        real(8),allocatable,dimension(:)::error_curr,error_array
        real(8),allocatable,dimension(:) ::error
        integer,allocatable,dimension(:) :: itr_counter
        real(8) :: BW_RAS_Taylor, BW_RAS_AWE
        real(8) :: BW_MoM_Taylor,BW_MoM_AWE,BW_MoM_Chebyshev
        real :: finish_t_RAS_AWE,start_t_RAS_AWE


        call read_parameters()
        call set_scatterers_IDs(Scat)

        allocate(error_array(number_of_scatterers+1))

        z_vec = vec(0.d0,0.d0,1.d0)
        if(simulation_mode == 0) then
            OPEN(11, FILE='sources.dat')
        endif
        !        OPEN(12, FILE='UNKNOWNS.dat')
        OPEN(13, FILE='error.dat')
        write(13,*) 'Radius of cylinder           ','Curve Length                     ','SPW     ',&
        '#Testing pts    ','#Group           ','Error               ',&
        'Error in Current         ','Time Consumed (Seconds)','     #Iterations'
        open(12,FILE='source_positions.dat')
        OPEN(14, FILE='error_observation.dat')
        if(max_iteration > 0) then
            OPEN(15, FILE='Convergence_rate.dat')
            write(15,*) '       Loop Counters   ','Convergence Number'
        endif

        if(MoM_activation_flag == 5) then
            OPEN(16,FILE='BW_Potential_RAS_2D_AWE.dat')
            write(16,*) '       #Run             SPW             #Testing pts',&
            '          #Group           BW-RAS-Taylor               BW-RAS-AWE',&
            '         Consumed AWE Calculation Time             Center Frequency Error'
        endif



        N_Itr_Total = N_Iterations*samples_S*samples_N
        write(14,*) N_Itr_Total
        allocate( error_curr(N_Itr_Total))
        allocate(error(0:N_Itr_Total))
        allocate(itr_counter(number_of_scatterers))


        cnt_loops = 1
        error = 1.d0
        error_curr = 1.d0

        do cnt_Itr = 1,N_Iterations
            do cnt_S = 1,samples_S

                do cnt_N = 1,samples_N
                    call cpu_time(start_t)

                    Samples_per_wavelength = SPW_read(cnt_S)

                    !                        Scatterer1%N_group = N_group_read(cnt_N)
                    !                        Scatterer1%N_max = 2*Scatterer1%N_group*(max_iteration+1)

                    !! superquadratic 2D
                    do i=1,number_of_scatterers

                        if(scatterers_input_method(i) == 1) then
                            call discretize_superquad(Scat(i),Samples_per_wavelength,Scat)
                        elseif(scatterers_input_method(i) == 2) then
                            call discretize_scatterer_file(Scat(i),Samples_per_wavelength,Scat)
                        endif


                        if(Scat(i)%Problem_Type == 3) then !! IBC
                            call set_IBC_impedance_matrices(Scat(i))
                        else
                            allocate(Scat(i)%norm_eta_u(Scat(i)%M),Scat(i)%norm_eta_v(Scat(i)%M))
                            Scat(i)%norm_eta_u = 1.d0
                            Scat(i)%norm_eta_v = 1.d0
                        endif


                        if(Source_Placement == 1) then !! MAS
                            Scat(i)%N_group = Scat(1)%N_con !! Use all the available contour points for the inside problem
                            max_iteration = 0 !! Solution is made with only one Iteeration
                            Scat(i)%N_sub_groups =  1
                        elseif(Source_Placement == 2) then !! RAS
                            if(N_group_read(cnt_N) < 0) then
                                Scat(i)%N_group = 20*Scat(i)%M/(abs(N_group_read(cnt_N))*Samples_per_wavelength)
                            else
                                Scat(i)%N_group = N_group_read(cnt_N)
                            endif
                            Scat(i)%N_sub_groups =  maxval((/ceil(dble(Scat(i)%M*2)/dble(3*Scat(i)%N_group)),5/))
                        else
                            if(N_group_read(cnt_N) < 0) then
                                Scat(i)%N_group = 20*Scat(i)%M/(abs(N_group_read(cnt_N))   *Samples_per_wavelength)
                            else
                                Scat(i)%N_group = N_group_read(cnt_N)
                            endif
                            Scat(i)%N_sub_groups =  maxval((/ceil(dble(Scat(i)%M*2)/dble(3*Scat(i)%N_group)),5/))
                        endif


                        !

                        Scat(i)%N_max = 2*(Scat(i)%N_sub_groups+1)*(Scat(i)%N_group+Scat(i)%N_inside_sources)
                        if(RAS_solution_method /= 0) then
                            write(*,*) 'N_group=',Scat(i)%N_group,'Ourside Sources',Scat(i)%N_inside_sources,&
                            'number of groups=',Scat(i)%N_sub_groups,'N_max',Scat(i)%N_max
                        endif
                        call set_excitation(Scat(i),Scat(i)%testing_pt,Scat(i)%Ei,Scat(i)%Hi)
                        call allocate_arrays(Scat(i),.true.,Scat)
                    enddo

                    !                    call remove_redundant_testing_points(Scat)
                    !                    stop


                    if(RAS_solution_method == 1) then
                        call solve_scattering_multiscatterer_once(Scat,N_taylor,Tol,error(cnt_loops),itr_counter(1),Scat_big)
                    elseif(RAS_solution_method == 2) then
                        call solve_scattering_multiscatterer_2(Scat,N_taylor,Tol,error(cnt_loops),itr_counter(1))
                    elseif(RAS_solution_method == 3) then
                        call solve_scattering_multiscatterer(Scat,N_taylor,Tol,error(cnt_loops),itr_counter(1))
                    else
                        !! nothing
                    endif

                    if(RAS_solution_method /= 0) then
                        if(plot_current /= 0) then
                            do i = 1,number_of_scatterers
                                call plot_current_RAS(Scat(i))
                            enddo
                            call plot_far_field_RAS(Scat,180)
                        endif

                        call cpu_time(finish_t)
                        write(*,*) 'Elapsed Time = ', finish_t - start_t, 'seconds'

                        write(15,*) cnt_loops,'     ',cnt
                    endif
                    !                        Area = 2.d0*tpi*R_cylinder(cnt_R)**2.d0
                    length = tpi*Scat(1)%a_superquad
                    !                        error_curr(cnt_loops) = eval_error_current_Electric_or_Magnetic(Scatterer1,trim(current_compare_file))

                    !                    do i = 1,number_of_scatterers
                    !                        write(*,*) 'Scatterer',i,' took',Scat(i)%itr_counter,' iterations'
                    !                    enddo

                    if(MoM_activation_flag == 1) then
                        if(error(cnt_loops) == -2.d0) then
                            write(*,*) 'No need to evaluate MoM'
                        else
                            call cpu_time(start_t_MoM)
                            if(MoM_solution_method == 1) then

                                error_array = eval_surface_current_MoM_once(Scat_big,N_taylor,Scat)
                            else
                                error_array = eval_surface_current_error_multiscatterer_MoM(Scat,N_taylor)
                            endif
                            error_curr(cnt_loops) = error_array(1)

                            call cpu_time(finish_t_MoM)
                        endif
                        write(*,*) 'MoM Execution time =', finish_t_MoM - start_t_MoM, ' seconds'
                        call eval_monostatic_RCS_comparison(Scat)
                    elseif(MoM_activation_flag == 2) then !! analytic solution
                        if(number_of_scatterers > 1) then
                            write(*,*) 'WARNING: Analytic solutions are not available for multi-scattering problems'
                        else
                            if(Scat(1)%Problem_Type == 1) then        !! PEC
                                error_curr(cnt_loops) =  eval_surface_current_error(Scat(1))
                            elseif(Scat(1)%Problem_Type == 2) then !! PMC
                                error_curr(cnt_loops) =  eval_surface_current_error_PMC(Scat(1))
                            elseif(Scat(1)%Problem_Type == 3) then !! IBC
                                error_curr(cnt_loops) =  eval_surface_current_error_IBC(Scat(1))
                            elseif(Scat(1)%Problem_Type == 4) then    !! Dielectric
                                error_curr(cnt_loops) =  eval_surface_current_Dielectric(Scat(1))
                            endif
                            call compare_far_field_scatterers(Scat,721)

                            deallocate(Scat(1)%I_u,Scat(1)%I_v,Scat(1)%M_u,Scat(1)%M_v)
                        endif
                    elseif(MoM_activation_flag == 3) then !! color map requist
                        R_probs = R_probes_in
                        call post_processing_colormap_imaging(Scat(1),R_probs,R_calc_in,5,0.04*lambda)

                    elseif(MoM_activation_flag == 4) then !! bandwidth measurments
                        if(error(cnt_loops) == -2.d0) then
                            write(*,*) 'No need to evaluate MoM'
                        else
                            write(*,*) 'Beware that the SPW in RAS may be different than what is in MoM'
                            call cpu_time(start_t_MoM)
                            error_array = eval_surface_current_error_multiscatterer_MoM(Scat,N_taylor)
                            error_curr(cnt_loops) = error_array(1)
                            call cpu_time(finish_t_MoM)
                            write(*,*) 'MoM-AWE Execution time =', finish_t_MoM - start_t_MoM, ' seconds'
                            if(Wideband_type > 0) then
                                write(*,*) '================================================================='
                                call cpu_time(start_t)
                                call get_solution_bandwidth(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                n_points_freq-n_points_freq/2)
                                call cpu_time(finish_t)
                                write(*,*) 'Taylor Series Solution Bandwidth Time      = ', finish_t - start_t, 'seconds'
!                                call cpu_time(start_t)
                                call cpu_time(start_t)

                                call plot_solution_bandwidth_Pade(Scat(1),N_taylor,fl_r,&
                                fh_r,n_points_freq/2,n_points_freq-n_points_freq/2)
                                call cpu_time(finish_t)
                                write(*,*) 'AWE Solution Bandwidth Time      = ', finish_t - start_t, 'seconds'


!                                call cpu_time(finish_t)
!                                write(*,*) 'Pade Approximation Solution Bandwidth Time = ', finish_t - start_t, 'seconds'
!                                write(*,*) '================================================================='
!                                call cpu_time(start_t)
!                                call Eval_chebyshev_expansion_coefficeints(Scat(1),N_taylor,Pade_L,fl_r,fh_r)
!                                write(*,*) '================================================================='
!                                call plot_solution_bandwidth_Chebyshev(Scat(1),N_taylor,fl_r,&
!                                fh_r,n_points_freq/2,n_points_freq-n_points_freq/2)
!                                call cpu_time(finish_t)
!                                write(*,*) 'Chebyshev Solution Bandwidth Time      = ', finish_t - start_t, 'seconds'
!                                write(*,*) '================================================================='
                                call Export_RAS_AWE_Chebyshev_Monostatic_RCS(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                    n_points_freq-n_points_freq/2)
!                                Export_MoM_AWE_Chebyshev_Monostatic_RCS(this,N_order,fl_f0,fh_f0,nl,nh)
!                                call Export_MoM_AWE_Chebyshev_Monostatic_RCS(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
!                                    n_points_freq-n_points_freq/2)
                                call cpu_time(start_t_RAS_AWE)
                                call Evaluate_MoM_AWE_Chebyshev_bandwidth(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                        n_points_freq-n_points_freq/2,BW_MoM_Taylor,BW_MoM_AWE,BW_MoM_Chebyshev)
                                call cpu_time(finish_t_RAS_AWE)
                                write(*,*) 'MoM-AWE Time = ',finish_t_RAS_AWE - start_t_RAS_AWE
        !                        write(16,*) cnt_loops,'    ',Samples_per_Wavelength,'    ',Scat(1)%M,'    ',&
        !                        Scat(1)%N_group,'    ','    ',BW_RAS_Taylor,&
        !                        '    ',BW_RAS_AWE,finish_t_RAS_AWE - start_t_RAS_AWE,'              ',error(cnt_loops)
                                call Export_MoM_AWE_Chebyshev_Monostatic_RCS(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                            n_points_freq-n_points_freq/2)

                            endif



                        endif
                    elseif(MoM_activation_flag == 5) then !! bandwidth measurements for RAS only

!                        write(16,*) '#Run             SPW             #Testing pts          #Group           BW-RAS-Taylor               BW-RAS-AWE         Consumed AWE Calculation Time'
                        call cpu_time(start_t_RAS_AWE)
                        call Evaluate_RAS_AWE_bandwidth(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                n_points_freq-n_points_freq/2,BW_RAS_Taylor,BW_RAS_AWE)
                        call cpu_time(finish_t_RAS_AWE)
                        write(*,*) 'RAS-AWE Time = ',finish_t_RAS_AWE - start_t_RAS_AWE
                        write(16,*) cnt_loops,'    ',Samples_per_Wavelength,'    ',Scat(1)%M,'    ',&
                        Scat(1)%N_group,'    ','    ',BW_RAS_Taylor,&
                        '    ',BW_RAS_AWE,finish_t_RAS_AWE - start_t_RAS_AWE,'              ',error(cnt_loops)


                    elseif(MoM_activation_flag == 6) then !! bandwidth measurements using different Methods AWE, Chebyshev
                        if(error(cnt_loops) == -2.d0) then
                            write(*,*) 'No need to evaluate MoM'
                        else
                            write(*,*) 'Beware that the SPW in RAS may be different than what is in MoM'
                            call cpu_time(start_t_MoM)
                            error_array = eval_surface_current_error_multiscatterer_MoM(Scat,N_taylor)
                            error_curr(cnt_loops) = error_array(1)
                            call cpu_time(finish_t_MoM)
                            write(*,*) 'MoM-AWE Execution time =', finish_t_MoM - start_t_MoM, ' seconds'
                            if(Wideband_type > 0) then
                                write(*,*) '================================================================='
                                call cpu_time(start_t)
                                call get_solution_bandwidth(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                n_points_freq-n_points_freq/2)
                                call cpu_time(finish_t)
                                write(*,*) 'Taylor Series Solution Bandwidth Time      = ', finish_t - start_t, 'seconds'
!                                call cpu_time(start_t)
                                call plot_solution_bandwidth_Pade(Scat(1),N_taylor,fl_r,&
                                fh_r,n_points_freq/2,n_points_freq-n_points_freq/2)
!                                call cpu_time(finish_t)
!                                write(*,*) 'Pade Approximation Solution Bandwidth Time = ', finish_t - start_t, 'seconds'
                                write(*,*) '================================================================='


                                call cpu_time(start_t)
                                call Eval_chebyshev_expansion_coefficeints(Scat(1),N_taylor,Pade_L,fl_r,fh_r)
                                call plot_solution_bandwidth_Chebyshev(Scat(1),N_taylor,fl_r,&
                                fh_r,n_points_freq/2,n_points_freq-n_points_freq/2)
                                call cpu_time(finish_t)
                                write(*,*) 'Chebyshev Solution Bandwidth Time      = ', finish_t - start_t, 'seconds'

                                call Export_RAS_AWE_Chebyshev_Monostatic_RCS(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                    n_points_freq-n_points_freq/2)
                                write(*,*) '================================================================='
                            endif



                        endif

                    elseif(MoM_activation_flag == 7) then !! bandwidth measurements for MoM only
                        write(*,*) 'Beware that the SPW in RAS may be different than what is in MoM'
                            call cpu_time(start_t_MoM)
                            error_array = eval_surface_current_error_multiscatterer_MoM(Scat,N_taylor)
                            error_curr(cnt_loops) = error_array(1)
                            call cpu_time(finish_t_MoM)
                            write(*,*) 'MoM-AWE Execution time =', finish_t_MoM - start_t_MoM, ' seconds'
                        write(*,*) 'WARNING: There is no need for multiple iterations',&
                        ' for this routine because it is deterministic'
!                        write(16,*) '#Run             SPW             #Testing pts          #Group           BW-RAS-Taylor               BW-RAS-AWE         Consumed AWE Calculation Time'
                        call cpu_time(start_t_RAS_AWE)
                        call Evaluate_MoM_AWE_Chebyshev_bandwidth(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                n_points_freq-n_points_freq/2,BW_MoM_Taylor,BW_MoM_AWE,BW_MoM_Chebyshev)
                        call cpu_time(finish_t_RAS_AWE)
                        write(*,*) 'MoM-AWE Time = ',finish_t_RAS_AWE - start_t_RAS_AWE
!                        write(16,*) cnt_loops,'    ',Samples_per_Wavelength,'    ',Scat(1)%M,'    ',&
!                        Scat(1)%N_group,'    ','    ',BW_RAS_Taylor,&
!                        '    ',BW_RAS_AWE,finish_t_RAS_AWE - start_t_RAS_AWE,'              ',error(cnt_loops)
                        call Export_MoM_AWE_Chebyshev_Monostatic_RCS(Scat(1),N_taylor,fl_r,fh_r,n_points_freq/2,&
                                    n_points_freq-n_points_freq/2)

                    endif

                    write(13,*) Scat(1)%a_superquad,length,Samples_per_Wavelength,Scat(1)%M,&
                    Scat(1)%N_group,'    ',error(cnt_loops),&
                    error_curr(cnt_loops),finish_t - start_t,'      ',itr_counter
                    do i = 1,number_of_scatterers
                        call deallocate_initialized_arrays(Scat(i))
                    enddo
                    if(RAS_solution_method == 1) then
                        call deallocate_initialized_arrays(Scat_big)
                    endif
                    cnt_loops = cnt_loops+1
                enddo

            enddo
        enddo

        if(simulation_mode == 0) then
            close(11)
        endif
        close(13)
        close(12)
        close(14)
        close(15)

        if(MoM_activation_flag == 5) then
            CLOSE(16)
        endif

        deallocate(error_array,error_curr,itr_counter)
        deallocate(Scat,error)
        deallocate(scatterers_input_method,scatterer_input_file_names)
        call deallocate_materials()
        if(allocated(n_C_i_matrix)) then
            deallocate(n_C_i_matrix,n_P_i_matrix)
        endif
        if(allocated(n_X_i_matrix)) then
            deallocate(n_X_i_matrix)
        endif
        if(allocated(Line_sources)) then
            deallocate(Line_sources)
        endif
    end subroutine run_main_program


end module RSM_3D_sph


program RSM_3D_cond_shp_Elec_or_Mag
    use RSM_3D_sph
    implicit none

    call run_main_program()

!    write(*,*) 'Press Any Key to Exit'
!    read(*,*)

!    write(*,*) log((-1.d0,0.d0))
!    complex*16,allocatable,dimension(:) :: xH0,xH1,xH2,H1,H0
!    integer :: N=5
!    integer :: nn
!    complex*16 :: h0,h1

!    call besselh2_01(0.1d0,h0,h1)
!    write(*,*) h0,h1
!    call besselh2_01((-1.0d0,-0.0d0),h0,h1)
!    write(*,*) h0,h1
!    allocate(xH0(0:N),xH1(0:N),xH2(0:N),H1(0:N),H0(0:N))
    !    call Hankel_derivatives(1.5d0,1.5d0,N,xH0,xH1,xH2,H1)
!    CALL Hankel_derivatives2(0.314d0,N,H0,H1)
!    do nn = 0,N
!    write(*,*) H0(nn),H1(nn)
!    enddo
!    deallocate(xH0,xH1,xH2,H1,H0)
end program RSM_3D_cond_shp_Elec_or_Mag
