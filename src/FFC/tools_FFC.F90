subroutine FindAlpha(f, tol, a_root)
  implicit none
  real*8, intent(in) :: f    ! The value of f to solve for
  real*8, intent(in) :: tol  ! Convergence tolerance
  real*8, intent(out) :: a_root  ! The root value

  real*8 :: a_low, a_high, a_guess, a_new, f_guess, df_guess
  real*8 :: a_mid, diff_low, diff_high, diff_mid, f_cut
  integer :: iter, max_iter

  if (f > 1.1d0 .or. f < -1.1d0) then
    write(*,*) 'f out of bounds', f
    stop
  endif

  f_cut = MAX(-0.99d0, MIN(0.99d0, f))
  a_guess = 15.0d0 * f_cut / (5.0d0 - 3.0d0 * f_cut**2 + f_cut**3 - 3.0d0 * f_cut**4)

  a_low = -105.0d0
  a_high = 110.0d0
  max_iter = 100

  ! Evaluate bounds
  diff_low = MinerboFluxFactor(a_low, tol) - f_cut
  diff_high = MinerboFluxFactor(a_high, tol) - f_cut

  a_new = MERGE(a_guess, 0.5d0 * (a_low + a_high), a_guess > a_low .and. a_guess < a_high)

  do iter = 1, max_iter
    f_guess = MinerboFluxFactor(a_new, tol) - f_cut
    df_guess = MinerboFluxFactorDerivative(a_new)

    if (abs(f_guess) < tol) then
      a_root = a_new
      return
    endif

    if (abs(df_guess) > 1.0d-10) then
      a_new = a_new - f_guess / df_guess
    else
      exit
    endif

    if (a_new < a_low .or. a_new > a_high) exit
  end do

  ! Fallback: Bisection
  do iter = 1, max_iter
    a_mid = 0.5d0 * (a_low + a_high)
    diff_mid = MinerboFluxFactor(a_mid, tol) - f_cut

    if (abs(diff_mid) < tol) then
      a_root = a_mid
      return
    endif

    if (diff_mid * diff_low < 0.0d0) then
      a_high = a_mid
      diff_high = diff_mid
    else
      a_low = a_mid
      diff_low = diff_mid
    endif
  end do

  print *, 'FindAlpha did not converge.'
  stop

contains

  function MinerboFluxFactor(x, tol) result(y)
    real*8, intent(in) :: x, tol
    real*8 :: y, tanh_x

    if (abs(x) < tol) then
      y = x / 3.0d0
    else
      tanh_x = tanh(x)
      y = 1.0d0 / tanh_x - 1.0d0 / x
    end if
  end function MinerboFluxFactor

  function MinerboFluxFactorDerivative(x) result(dy)
    real*8, intent(in) :: x
    real*8 :: dy

    dy = -1.0d0 / sinh(x)**2 + 1.0d0 / (x * x)
  end function MinerboFluxFactorDerivative

end subroutine FindAlpha

subroutine FindEta(alpha, e, eta)
  real*8, intent(in) :: alpha, e
  real*8, intent(out) :: eta
  real*8 :: sinh_alpha

  if (e <= 0.0d0 .or. alpha == 0.0d0) then
    eta = 1.0d99  ! avoid log(0)
  else
    if (abs(alpha) < 1.0d-4) then
      sinh_alpha = alpha + alpha**3 / 6.0d0
    else
      sinh_alpha = sinh(alpha)
    endif
    eta = -log(e * alpha / sinh_alpha)
  endif

end subroutine FindEta

