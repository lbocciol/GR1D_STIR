! burn.F90
!
! Units:  rho [g/cc], T [K], X [mass fractions], edot [erg/g/s], e_released [erg/g]
! Convention: pynet RHS takes mass fractions X and returns dY/dt (Y = X/A).
!             We convert X->Y on entry and Y->X on exit.
! Energy: integrated via trapezoidal rule using ener_gener_f output.
! Integrator: single backward-Euler step with Newton-Raphson; LAPACK dgesv solves
!             the linear system at each Newton iteration.
! Ye: this network has only strong/EM rates so Ye is conserved; we assert it.

module burn

  use pynet, only: nspec, aion, zion, &
                   rhs_f, jac_f, ener_gener_f, network_init

  implicit none
  private

  public :: burn_init, burn_rhs, burn_state, burn_newton

  integer,  parameter :: max_newton = 50
  real(8),  parameter :: newton_tol = 1.0d-10
  real(8),  parameter :: ye_tol     = 1.0d-8

contains

  ! ---------------------------------------------------------------------------
  subroutine burn_init()
    call network_init()
  end subroutine burn_init

  ! ---------------------------------------------------------------------------
  ! Single RHS + energy evaluation (no time integration).
  !   in : rho, T, X(nspec) mass fractions
  !   out: dYdt(nspec) [mol/g/s], edot [erg/g/s]
  subroutine burn_rhs(rho, T, X, dYdt, edot)
    real(8), intent(in)  :: rho, T, X(nspec)
    real(8), intent(out) :: dYdt(nspec), edot
    real(8) :: enuc, enu_weak

    enu_weak = 0.0d0
    dYdt     = 0.0d0
    call rhs_f(rho, T, X, dYdt, enu_weak)
    call ener_gener_f(dYdt, enuc)
    edot = enuc + enu_weak
  end subroutine burn_rhs

  ! ---------------------------------------------------------------------------
  ! Advance composition over hydro substep dt.
  !   in : rho, T, X_in(nspec), dt [s]
  !   out: X_out(nspec), e_released [erg/g], ierr (0=ok), Ye_out (optional)
  subroutine burn_state(rho, T, X_in, dt, X_out, e_released, ierr, Ye_out)
    real(8), intent(in)            :: rho, T, X_in(nspec), dt
    real(8), intent(out)           :: X_out(nspec), e_released
    integer, intent(out)           :: ierr
    real(8), intent(out), optional :: Ye_out

    real(8) :: Y(nspec), Ye_in, Ye_out_loc
    logical :: ok

    ierr = 0

    ! X -> Y
    Y = X_in / aion
    Ye_in = sum(zion * Y)

    call burn_newton(rho, T, Y, dt, e_released, ok)
    if (.not. ok) ierr = 1

    ! floor negatives
    where (Y < 0.0d0) Y = 0.0d0

    ! Y -> X and renormalize if needed
    X_out = Y * aion
    if (abs(sum(X_out) - 1.0d0) > 1.0d-10) X_out = X_out / sum(X_out)

    ! Ye conservation check
    Ye_out_loc = sum(zion * (X_out / aion))
    if (abs(Ye_out_loc - Ye_in) > ye_tol) &
      write(*,'(a,es10.3)') "burn WARNING: Ye drift =", Ye_out_loc - Ye_in

    if (present(Ye_out)) Ye_out = Ye_out_loc

  end subroutine burn_state

  ! ---------------------------------------------------------------------------
  ! Backward-Euler step with Newton-Raphson.
  ! Solves  Y^{n+1} - Y^n - dt*f(Y^{n+1}) = 0
  ! using   [I - dt*J] * delta = -(Y^k - Y^n - dt*f(Y^k))
  ! Energy via trapezoidal: e = 0.5*(edot_in + edot_out)*dt
  subroutine burn_newton(rho, T, Y, dt, e_step, converged)
    real(8), intent(in)    :: rho, T, dt
    real(8), intent(inout) :: Y(nspec)     ! in: Y^n; out: Y^{n+1}
    real(8), intent(out)   :: e_step
    logical, intent(out)   :: converged

    real(8) :: Y_n(nspec), Y_k(nspec), X_k(nspec)
    real(8) :: f_n(nspec), f_k(nspec), delta(nspec)
    real(8) :: J(nspec,nspec), A(nspec,nspec)
    real(8) :: enu_weak, enuc_n, enuc_k
    real(8) :: res
    integer :: iter, info, ipiv(nspec), i

    converged = .false.
    Y_n = Y
    Y_k = Y

    ! edot at start
    X_k = Y_n * aion
    f_n = 0.0d0; enu_weak = 0.0d0
    call rhs_f(rho, T, X_k, f_n, enu_weak)
    call ener_gener_f(f_n, enuc_n)
    enuc_n = enuc_n + enu_weak

    do iter = 1, max_newton
      X_k = Y_k * aion
      f_k = 0.0d0; enu_weak = 0.0d0
      call rhs_f(rho, T, X_k, f_k, enu_weak)

      J = 0.0d0
      call jac_f(rho, T, X_k, J)

      ! A = I - dt*J
      A = -dt * J
      do i = 1, nspec
        A(i,i) = A(i,i) + 1.0d0
      end do

      ! rhs = -(Y_k - Y_n - dt*f_k)
      delta = -(Y_k - Y_n - dt * f_k)

      call dgesv(nspec, 1, A, nspec, ipiv, delta, nspec, info)
      if (info /= 0) return

      Y_k = Y_k + delta

      res = maxval(abs(delta) / (abs(Y_k) + 1.0d-30))
      if (res < newton_tol) then
        converged = .true.
        exit
      end if
    end do

    if (.not. converged) return

    ! edot at end
    X_k = Y_k * aion
    f_k = 0.0d0; enu_weak = 0.0d0
    call rhs_f(rho, T, X_k, f_k, enu_weak)
    call ener_gener_f(f_k, enuc_k)
    enuc_k = enuc_k + enu_weak

    e_step = 0.5d0 * (enuc_n + enuc_k) * dt
    Y = Y_k

  end subroutine burn_newton

end module burn
