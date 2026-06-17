! ============================================================================
! Module: sim_par
! Purpose: Global simulation parameter store. Acts as shared state for the
!          entire solver. All modules USE it to access problem configuration.
!          Variables are populated by read_parameters() at program start and
!          remain constant thereafter, except for n_C_i_matrix / n_P_i_matrix /
!          n_X_i_matrix which are computed once during AWE initialisation.
!
! Key variable groups:
!   Wave parameters  : k0, k1, az, ar, ak, eta1, lambda, frequency
!   Iteration limits : max_iteration, max_bouncing_iteration
!   Solver flags     : RAS_solution_method, MoM_solution_method,
!                      Matrix_Solution_Method, MoM_activation_flag
!   Geometry bounds  : bound, outside_bound, corner_sep, corner_sep2
!   AWE / wideband   : Wideband_type, N_taylor, Pade_L, fl_r, fh_r,
!                      n_points_freq, n_C_i_matrix
!   Material library : materials (type material_lib)
!   Source model     : Source_Model, Line_sources(:)
!
! Author : Mohamed A. Moharram Hassan
! License: Apache 2.0 -- see LICENSE
! ============================================================================
module sim_par
    use Operations
    implicit none

    integer :: number_of_scatterers
    integer, allocatable, dimension(:) :: scatterers_input_method
    character*30,allocatable,dimension(:) :: scatterer_input_file_names
    integer :: Source_Model !! 0-Plane Wave , 1- Line Sources
    integer :: N_Line_Sources
    type Line_source
        integer :: Source_Type
        integer :: Source_orientation
        real(8) ::x_s,y_s
        complex*16 :: Amp
        logical :: is_normalized
    end type Line_source
    type(Line_source),allocatable,dimension(:) :: Line_sources

    real(8) :: TOL !! acceptable error in the satisfaction of the boundary condition
    real(8) :: TOL_segmentation !! Tolerance in the segmentation of the structure
    real(8) :: k0 ! free-space wavelength
    real(8) :: k1 !! the wave-numbers of the background medium
    real(8) :: az !! the constant the is used to computed kz = az_l*k0
    real(8) :: ar !! the constant to compute kr = ar*k0 for the background material
    real(8) :: ak !! the constant to compute k from k0 k = ak*k0, ak = sqrt(mu_r*eps_r)


    real(8) :: MAX_ALLOWED_CONDITION_NUMBER
    real(8) :: Y_mat_norm_limit
    integer :: Matrix_Solution_Method

    real(8) :: BW_limit = 0.1d0

    real(8) :: eta1 !! the intrinsic impedance of the background medium
    integer :: Samples_per_wavelength,MoM_Samples_per_wavelength
    real(8) :: frequency,lambda
    integer :: excitation_type
    real(8) :: theta_i,phi_i,alpha_i
    real(8),allocatable,dimension(:,:) :: physical_parameters
    integer :: samples_per_wavelength_contour
    integer :: samples_per_wavelength_outside
    integer :: non_uniform_type_read,N_o_read,A_segmentation_read

    !! Parameters for RSM
    integer :: max_iteration !! maximum number of iterations
    real(8) :: bound,outside_bound,corner_sep,corner_sep2 ! the maximum allowable distance from the conductor a source can exist in
    real(8),parameter :: delta_pos = 1.d-1 !! distance from corner sources to corners
    real(8),parameter :: TOL_internal = 1.d-2 !! the internal tolerance of satisfaction of boundary conditions
    real(8) :: R_cyl
    !!!!!!!!!!!!

    !! Simulation mode parameters
    integer,allocatable,dimension(:) :: N_group_read !! the simulaing N_group array
    integer,allocatable,dimension(:) :: SPW_read !!Samples per wavelength array
    integer :: simulation_mode !! 0 - just a single run for one simulation parameters set
    !! 1 - sweep R_cylinder - obsolute  (deleted)
    !! 2 - sweep N_group
    !! 3 - sweep Samples_per_Wavelength
    !! 4 - all parameters
    !! 5 - Iterations (Using all given values for parameters)
    !! 6 - sweep Electric_sources_ratio
    integer :: plot_current !! if plot_current ~= 0, the currents are to be plotted at the end
    integer :: samples_R,samples_N,samples_S !! the number of sampling points for parameter sweep for R_cylinder, N_group, Samples_per_Wavelength
    integer :: cnt_loops,N_iterations


    !! Formulation type of the problem for IBC
    integer :: FormulationType !! 1- IBCE, 2- IBCH, 3- IBCC, 4- IBC as appeared in Dr. Kihsk's paper
    integer :: Source_Placement !! Source Placement Procedure
    !! 1- MAS (Contour), 2- RAS (All Random), 3- MAS + RAS (Contour and Random)
    integer :: MoM_activation_flag !! a flag to set the required type of verification
    !! 0- none, 1- MoM, 2- Analytic

    real(8) :: R_probes_in !! the position of the probes
    real(8) :: R_calc_in !! the radius to calculate within the scattered fields
    real(8),dimension(2) :: superquad_center

    integer :: Contour_Sources_Type_read

    type(vector) :: z_vec
    complex*16 :: c_zero = (0.d0,0.d0)
    !! parameters used by bouncing field method
    integer :: max_bouncing_iteration

    !! Wide band impelemtation
    integer :: Wideband_type
    integer :: N_taylor,Pade_L
    real(8) :: fl_r,fh_r
    integer :: n_points_freq
    integer,allocatable,dimension(:,:) :: n_C_i_matrix
    integer,allocatable,dimension(:,:) :: n_P_i_matrix,n_X_i_matrix

    type material_lib
        integer :: N_diel,N_ibc
        integer,allocatable,dimension(:) :: ID_diel,ID_ibc
        complex*16,allocatable,dimension(:) :: eta_zz,eta_zt,eta_tz,eta_tt
        real(8),allocatable,dimension(:) ::eps_r_p,mu_r_p,eps_r_pp,mu_r_pp
    end type material_lib

    type(material_lib) :: materials
    integer,allocatable,dimension(:) :: material_id_read

    integer :: RAS_solution_method,MoM_solution_method

!    real(8) :: BW_limit

end module sim_par
