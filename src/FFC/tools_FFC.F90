subroutine FindAlpha(f, tol, a_root)
  implicit none
  real*8, intent(in) :: f    ! The value of f to solve for
  real*8, intent(in) :: tol  ! Convergence tolerance
  real*8, intent(out) :: a_root  ! The root value

  real*8 :: a_low, a_high ! Bounds
  real*8 :: a_guess, a_new, f_guess, df_guess
  real*8 :: a_mid, diff_low, diff_high, diff_mid
  integer :: max_iter, iter
  real*8 :: f_cut
  logical :: use_bisection

  if ( f > 1.1d0 .or. f < -1.1d0 ) then
    write(*,*) 'f out of bounds', f
    stop
  endif

  f_cut = MIN(f, 0.99d0)
  f_cut = MAX(f_cut, -0.99d0)

  ! Use polynomial approximation for initial guess
  a_guess = 15.0d0 * f_cut / (5.0d0 - 3.0d0 * f_cut**2 + f_cut**3 - 3.0d0 * f_cut**4)

  ! Set initial bounds
  a_low  = -105.0d0
  a_high = 110.0d0
  max_iter = 100
  use_bisection = .false.

  ! Evaluate function at bounds
  diff_low = MinerboFluxFactor(a_low, tol) - f_cut
  diff_high = MinerboFluxFactor(a_high, tol) - f_cut

  ! Check if a_guess is within bounds and use it
  if (a_guess > a_low .and. a_guess < a_high) then
    a_new = a_guess
  else
    a_new = 0.5d0 * (a_low + a_high) ! Fall back to midpoint
  endif

  ! Newton-Raphson iterations
  do iter = 1, max_iter
    f_guess = MinerboFluxFactor(a_new, tol) - f_cut
    df_guess = MinerboFluxFactorDerivative(a_new) ! Compute derivative

    ! Check for convergence
    if (abs(f_guess) < tol) then
      a_root = a_new
      return
    endif

    ! Compute Newton-Raphson step
    if (abs(df_guess) > 1.0d-10) then
      a_new = a_new - f_guess / df_guess
    else
      use_bisection = .true.
    endif

    ! Check if NR step goes out of bounds
    if (a_new < a_low .or. a_new > a_high) then
      use_bisection = .true.
    endif

    if (use_bisection) exit ! Switch to bisection if NR fails
  end do

  ! If NR fails, use bisection
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

  ! Function to calculate coth(x)
  function MinerboFluxFactor(x, tol) result(y)
    real*8, intent(in) :: x, tol
    real*8 :: y, tanh_x

    if ( abs(x) < tol ) then
      y = 1.0d0 / tanh(tol) - 1.0d0 / tol  ! Approximate coth(x) near zero
    else
      tanh_x = tanh(x)
      y = 1.0d0 / tanh_x - 1.0d0 / x
    end if

  end function MinerboFluxFactor

  function MinerboFluxFactorDerivative(x) result(dy)
    real*8, intent(in) :: x
    real*8 :: dy, sinh_x

    sinh_x = sinh(x)
    dy = -1.0d0 / (sinh_x * sinh_x) + 1.0d0 / (x * x)

  end function MinerboFluxFactorDerivative

end subroutine FindAlpha

! This is for Maxwell Boltzmann approximation
subroutine FindEta(alpha, e, eta)
  real*8, intent(in) :: alpha, e !here e is the occupancy number (i.e. from 0 o 1)
  real*8, intent(out) :: eta

  eta = -log(e * alpha / sinh(alpha))
  
end subroutine FindEta
