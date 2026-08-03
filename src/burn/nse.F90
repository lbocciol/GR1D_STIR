! nse.F90
!
! Nuclear Statistical Equilibrium (NSE) composition solver for the burn network.
!
! PURPOSE
!   At high temperature the composition is set by NSE rather than by the reaction
!   network: it is a unique function of (rho, T, Ye).  GR1D assumes NSE above T_eos_high
!   and integrates the network below it; when a zone cools back below T_eos_high the
!   network needs a starting composition.  This module provides that seed so the
!   network has a valid handoff (see the burn loop in Step.F90).
!
! METHOD  (Seitenzahl, Townsley, Peng & Truran 2009, ADNDT 95, 96; Coulomb fit from
!          Chabrier & Potekhin 1998 / Calder et al. 2007, ApJ 656, 313.)
!   Detailed balance for  (A_i,Z_i) <-> Z_i p + N_i n  gives each nucleus its Saha
!   abundance in terms of the proton and neutron kinetic chemical potentials
!   (their Eq. 8).  In molar abundances Y_i = n_i/(rho N_A), with the dimensionless
!   kinetic potentials  u_p = mu_p^kin/kT,  u_n = mu_n^kin/kT,  this reads
!
!       ln Y_i = lnpre_i + Z_i u_p + N_i u_n ,
!       lnpre_i = ln g_i + 3/2 ln A_i + ln(2 pi m_u kT/h^2) - ln(rho N_A)
!                 + Q_i/kT + Z_i f_p - f_i .
!
!   Here Q_i = (Z_i m_p + N_i m_n - m_i)c^2 is the binding energy, g_i = 2J_i+1 the
!   ground-state spin factor, and (Z_i f_p - f_i) the Chabrier-Potekhin Coulomb
!   correction (f = mu^coul/kT, evaluated at Gamma_i = Z_i^{5/3} Gamma_e).  The two
!   unknowns (u_p,u_n) are fixed by baryon and charge conservation,
!
!       sum_i A_i Y_i = 1 ,     sum_i Z_i Y_i = ye ,
!
!   by a standard 2-D Newton-Raphson, warm-started by a 1-D bisection (u_p = u_n).
!
! SCOPE / FAILURE
!   This targets the NSE regime, T >~ 5 GK.  Below that NSE does not hold and the
!   solver is not designed to converge.  It does NOT force a result: if Newton fails
!   to converge it returns ierr /= 0.  In particular, a network whose nuclei cannot
!   represent the requested Ye (e.g. any Ye /= 0.5 for the all-N=Z alpha chain, whose
!   charge constraint is not independent of baryon number) makes the Jacobian
!   singular and the solve fails -- which is the correct response to that physics,
!   not something to paper over.
!
! CONVENTIONS (chosen to feed straight into burn_newton)
!   rho [g/cm^3], T [K], ye = sum_i Z_i Y_i, Y_i = X_i/A_i [mol/g], same ordering as
!   pynet's aion/zion.  We keep the A_i baryon basis (sum A_i Y_i = 1, X_i = A_i Y_i)
!   rather than Seitenzahl's exact-mass refinement, because the network defines baryon
!   number via A_i; the ~1e-4 difference is far below the Coulomb effect.

module nse

  use composition, only: nspec, aion, zion, nuclei_binding_energy  ! full set: network plus any appended n,p
  implicit none
  private

  public :: nse_init, nse_solve, nse_solve_core

  ! --- physical constants (CGS) -------------------------------------------------
  real(8), parameter :: pi          = 3.14159265358979324d0
  real(8), parameter :: m_u         = 1.66053906660d-24   ! atomic mass unit [g]
  real(8), parameter :: k_boltz     = 1.380649d-16        ! Boltzmann   [erg/K]
  real(8), parameter :: h_planck    = 6.62607015d-27      ! Planck      [erg s]
  real(8), parameter :: n_avogadro  = 6.02214076d23       ! [1/mol]
  real(8), parameter :: mev_to_erg  = 1.602176634d-6      ! [erg/MeV]
  real(8), parameter :: e_charge    = 4.80320425d-10      ! elementary charge [esu]

  ! Binding energies Q_i [MeV] come straight from the composition module (FLASH
  ! `bion`; free p and n have Q = 0).

  ! --- solver controls ----------------------------------------------------------
  real(8), parameter :: exp_cap  = 700.0d0   ! clamp exp() argument (avoid overflow)
  real(8), parameter :: nse_tol  = 1.0d-10   ! Newton convergence tolerance
  real(8), parameter :: step_max = 2.0d0     ! cap on |Newton step| in (u_p,u_n)
  integer, parameter :: max_newton = 200     ! Newton iteration cap (else fail)
  integer, parameter :: max_bisect = 200     ! 1-D bisection iteration cap

  ! --- module data for the pynet network, filled once in nse_init ---------------
  real(8), save, allocatable :: binding_mev(:) ! Q_i for the pynet species [MeV]
  real(8), save, allocatable :: gspin(:)       ! g_i = 2J_i+1 for the pynet species
  logical, save :: initialized = .false.

contains

  ! ---------------------------------------------------------------------------
  ! One-time setup for the pynet network.  Called once at startup (start.F90).
  subroutine nse_init()
    integer :: i
    if (initialized) return
    allocate(binding_mev(nspec), gspin(nspec))
    binding_mev = nuclei_binding_energy
    ! Ground-state spin factors g_i = 2J_i+1: J=0 -> 1 for every nucleus in the
    ! present alpha network; a free nucleon (A=1) is J=1/2 -> 2.
    do i = 1, nspec
      if (nint(aion(i)) == 1) then
        gspin(i) = 2.0d0
      else
        gspin(i) = 1.0d0
      end if
    end do
    initialized = .true.
  end subroutine nse_init

  ! ---------------------------------------------------------------------------
  ! Solve NSE for the pynet network at (rho, T, ye).  Thin wrapper around
  ! nse_solve_core with the cached pynet species data.
  !   out: Y(nspec) [mol/g];  ierr (0 ok, 1 not initialized, else solver failure)
  subroutine nse_solve(rho, T, ye, Y, ierr)
    real(8), intent(in)  :: rho, T, ye
    real(8), intent(out) :: Y(nspec)
    integer, intent(out) :: ierr
    if (.not. initialized) then
      Y = 0.0d0; ierr = 1; return
    end if
    call nse_solve_core(aion, zion, binding_mev, gspin, rho, T, ye, Y, ierr)
  end subroutine nse_solve

  ! ---------------------------------------------------------------------------
  ! General, network-agnostic NSE solver.  Operates on caller-supplied species
  ! data, so it works for any network.
  !   in : aion_(:), zion_(:)   mass and proton numbers
  !        binding_(:) [MeV]    binding energies Q_i
  !        gspin_(:)            spin factors g_i = 2J_i+1
  !        rho [g/cc], T [K], ye
  !   out: Y(:) [mol/g] (sum A_i Y_i = 1, sum Z_i Y_i = ye); ierr (0 = ok)
  subroutine nse_solve_core(aion_, zion_, binding_, gspin_, rho, T, ye, Y, ierr)
    real(8), intent(in)  :: aion_(:), zion_(:), binding_(:), gspin_(:)
    real(8), intent(in)  :: rho, T, ye
    real(8), intent(out) :: Y(:)
    integer, intent(out) :: ierr

    real(8) :: nion_(size(aion_)), lnpre(size(aion_))
    real(8) :: eta, u_p, u_n

    Y     = 0.0d0
    nion_ = aion_ - zion_

    ! (rho,T,ye)-dependent part of ln(Y_i), including Coulomb corrections.
    call saha_lnpre(aion_, zion_, gspin_, binding_, rho, T, ye, lnpre)

    ! Warm start from the 1-D solution (u_p = u_n): satisfies baryon conservation.
    call solve_1d(aion_, lnpre, eta, ierr)
    if (ierr /= 0) return
    u_p = eta
    u_n = eta

    ! Refine to also satisfy charge conservation.  For a network that cannot reach
    ! the requested ye this does not converge and returns ierr /= 0.
    call solve_2d(aion_, zion_, nion_, lnpre, ye, u_p, u_n, ierr)
    if (ierr /= 0) return

    call saha_Y(lnpre, zion_, nion_, u_p, u_n, Y)
    Y = Y / sum(aion_ * Y)          ! make baryon conservation exact
  end subroutine nse_solve_core

  ! ===========================================================================
  ! Internal helpers
  ! ===========================================================================

  ! Chabrier & Potekhin (1998) fit to the Coulomb free energy per ion of a
  ! one-component plasma, f(Gamma) = mu^coul/kT (Seitenzahl Eq. 14 / Calder Eq. A1).
  pure function coulomb_f(gamma) result(f)
    real(8), intent(in) :: gamma
    real(8) :: f, a3
    real(8), parameter :: a1 = -0.9052d0, a2 = 0.6322d0
    if (gamma <= 0.0d0) then
      f = 0.0d0
      return
    end if
    a3 = -sqrt(3.0d0)/2.0d0 - a1/sqrt(a2)
    f =  a1 * ( sqrt(gamma*(a2 + gamma))                                  &
                - a2 * log( sqrt(gamma/a2) + sqrt(1.0d0 + gamma/a2) ) )   &
       + 2.0d0*a3 * ( sqrt(gamma) - atan(sqrt(gamma)) )
  end function coulomb_f

  ! Build lnpre_i, the (rho,T,ye)-dependent part of ln(Y_i) (see the header).
  ! Coulomb: Gamma_e = e^2 (4 pi n_e/3)^{1/3}/kT with n_e = rho N_A ye; the proton
  ! term f_p = f(Gamma_e) is the reference, f_i = f(Z_i^{5/3} Gamma_e) the nucleus.
  subroutine saha_lnpre(aion_, zion_, gspin_, binding_, rho, T, ye, lnpre)
    real(8), intent(in)  :: aion_(:), zion_(:), gspin_(:), binding_(:)
    real(8), intent(in)  :: rho, T, ye
    real(8), intent(out) :: lnpre(:)
    real(8) :: ln_theta, ln_rhoNA, kT_erg, n_e, a_e, gamma_e, f_p
    integer :: i

    ln_theta = 1.5d0 * log( 2.0d0*pi*m_u*k_boltz*T / h_planck**2 )
    ln_rhoNA = log( rho * n_avogadro )
    kT_erg   = k_boltz * T

    n_e     = rho * n_avogadro * ye
    a_e     = ( 3.0d0 / (4.0d0*pi*n_e) )**(1.0d0/3.0d0)
    gamma_e = e_charge**2 / (a_e * kT_erg)
    f_p     = coulomb_f(gamma_e)

    do i = 1, size(aion_)
      lnpre(i) = log(gspin_(i)) + 1.5d0*log(aion_(i)) + ln_theta - ln_rhoNA  &
                 + binding_(i)*mev_to_erg/kT_erg                             &
                 + zion_(i)*f_p - coulomb_f( zion_(i)**(5.0d0/3.0d0) * gamma_e )
    end do
  end subroutine saha_lnpre

  ! Y_i = exp( lnpre_i + Z_i u_p + N_i u_n ), argument clamped against overflow.
  subroutine saha_Y(lnpre, zion_, nion_, u_p, u_n, Y)
    real(8), intent(in)  :: lnpre(:), zion_(:), nion_(:), u_p, u_n
    real(8), intent(out) :: Y(:)
    Y = exp( min(lnpre + zion_*u_p + nion_*u_n, exp_cap) )
  end subroutine saha_Y

  ! Baryon residual of the 1-D problem (u_p = u_n = eta): g(eta) = sum A_i Y_i - 1
  ! with Y_i = exp(lnpre_i + A_i eta).
  pure function g_baryon(aion_, lnpre, eta) result(g)
    real(8), intent(in) :: aion_(:), lnpre(:), eta
    real(8) :: g
    g = sum( aion_ * exp(min(lnpre + aion_*eta, exp_cap)) ) - 1.0d0
  end function g_baryon

  ! ---------------------------------------------------------------------------
  ! 1-D solve with u_p = u_n = eta.  g_baryon is monotonically increasing in eta,
  ! so plain bisection finds its unique root.  Exact when the charge constraint is
  ! not an independent degree of freedom, and the warm start for the 2-D Newton.
  subroutine solve_1d(aion_, lnpre, eta, ierr)
    real(8), intent(in)  :: aion_(:), lnpre(:)
    real(8), intent(out) :: eta
    integer, intent(out) :: ierr
    real(8) :: lo, hi
    integer :: it

    lo = -400.0d0
    hi =  100.0d0
    if (g_baryon(aion_, lnpre, lo) > 0.0d0 .or. &
        g_baryon(aion_, lnpre, hi) < 0.0d0) then
      ierr = 2; eta = 0.0d0; return        ! root not bracketed -> fail
    end if

    do it = 1, max_bisect
      eta = 0.5d0*(lo + hi)
      if (hi - lo < 1.0d-13) exit
      if (g_baryon(aion_, lnpre, eta) < 0.0d0) then
        lo = eta
      else
        hi = eta
      end if
    end do
    ierr = 0
  end subroutine solve_1d

  ! ---------------------------------------------------------------------------
  ! 2-D Newton-Raphson for (u_p,u_n) on baryon and charge conservation.  The only
  ! safeguard is a cap on the step length: near an iron-peak state one nucleus
  ! dominates and the Jacobian is ill-conditioned, so the raw Newton step can
  ! overshoot wildly.  Capping keeps the Newton direction while preventing the
  ! overshoot; near the solution steps are small and convergence is quadratic.
  ! If it still does not converge within max_newton (e.g. a network that cannot
  ! reach the requested ye, where the Jacobian is singular) it reports failure.
  subroutine solve_2d(aion_, zion_, nion_, lnpre, ye, u_p, u_n, ierr)
    real(8), intent(in)    :: aion_(:), zion_(:), nion_(:), lnpre(:), ye
    real(8), intent(inout) :: u_p, u_n
    integer, intent(out)   :: ierr

    real(8) :: Y(size(aion_))
    real(8) :: g1, g2, j11, j12, j21, j22, det, dup, dun, scale
    integer :: it

    do it = 1, max_newton
      call saha_Y(lnpre, zion_, nion_, u_p, u_n, Y)
      g1 = sum(aion_*Y) - 1.0d0
      g2 = sum(zion_*Y) - ye
      if (abs(g1) < nse_tol .and. abs(g2) < nse_tol) then
        ierr = 0; return
      end if

      ! analytic 2x2 Jacobian  d(g1,g2)/d(u_p,u_n)  (dY_i/du_p = Z_i Y_i, etc.)
      j11 = sum(aion_*zion_*Y);  j12 = sum(aion_*nion_*Y)
      j21 = sum(zion_*zion_*Y);  j22 = sum(zion_*nion_*Y)
      det = j11*j22 - j12*j21

      ! Newton step  delta = -J^{-1} (g1,g2), length-capped for robustness
      dup = -( j22*g1 - j12*g2) / det
      dun = -(-j21*g1 + j11*g2) / det
      scale = max(abs(dup), abs(dun))
      if (scale > step_max) then
        dup = dup * (step_max/scale)
        dun = dun * (step_max/scale)
      end if
      u_p = u_p + dup
      u_n = u_n + dun
    end do
    ierr = 3                               ! did not converge -> fail
  end subroutine solve_2d

end module nse
