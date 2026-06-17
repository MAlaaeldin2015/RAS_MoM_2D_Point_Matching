! =============================================================================
! Module: constants
! Purpose: Defines universal physical and mathematical constants used throughout
!          the RAS/MoM solver. All values are double-precision (real(8)).
!
! Author : Mohamed A. Moharram Hassan
! License: Apache 2.0 â see LICENSE
! =============================================================================
module constants
implicit none

  ! Imaginary unit  j = sqrt(-1)
  COMPLEX*16, PARAMETER :: cj = (0.d0, 1.d0)

  ! Mathematical constants
  REAL(8), PARAMETER :: pi  = 3.141592653589793238462643d0  ! Ï
  REAL(8), PARAMETER :: tpi = 2.d0 * pi                     ! 2Ï
  REAL(8), PARAMETER :: dpr = 180.d0 / pi                   ! degrees-per-radian

  ! Electromagnetic constants (SI)
  REAL(8), PARAMETER :: mu0      = 4d-7 * pi                        ! Permeability of free space [H/m]
  REAL(8), PARAMETER :: cspeed   = 2.99792458d8                     ! Speed of light [m/s]
  REAL(8), PARAMETER :: epsilon0 = 1.d0 / (mu0 * cspeed**2)        ! Permittivity of free space [F/m]
  REAL(8), PARAMETER :: eta0     = 1.d0 / (cspeed * epsilon0)       ! Intrinsic impedance of free space [Î©] â 377 Î©

  ! EulerâMascheroni constant (used in Bessel/Hankel series expansions)
  REAL(8), PARAMETER :: gamma_const = &
      0.57721566490153286060651209008240243104215933593992d0

end module constants
