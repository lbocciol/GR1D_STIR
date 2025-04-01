subroutine FindAlpha(f, tol, a_root)

  implicit none
  real*8, intent(in) :: f    ! The value of f to solve for
  real*8, intent(in) :: tol  ! Convergence tolerance
  real*8, intent(out) :: a_root  ! The root value

  real*8 :: a_low, a_high ! Initial bounds
  real*8 :: a_mid, diff_low, diff_high, diff_mid
  integer :: max_iter, iter
  real*8 :: f_cut

  if ( f > 1.1d0 .or. f < -1.1d0 ) then
    write(*,*) 'f out of bounds', f
    stop
  endif

  f_cut = MIN(f, 0.99d0)
  f_cut = MAX(f_cut, -0.99d0)

  ! Just use simple polynomial fit:
  ! a_root = 15.0d0*abs(f_cut) / (5.0d0 + 3.0d0*abs(f_cut)**2 + abs(f_cut)**3 - 3.0d0*abs(f_cut)**4)
  ! return

  ! Determine bounds. f is a monotonic function of alpha
  ! Never end up on zero where tanh = 0
  a_low  = -105.0d0
  a_high = 110.0d0

  ! Define the maximum iterations for the bisection method
  max_iter = 100

  ! Evaluate the function at the bounds
  diff_low = MinerboFluxFactor(a_low, tol) - f_cut
  diff_high = MinerboFluxFactor(a_high, tol) - f_cut

  ! Check if the function values at the bounds have opposite signs
  if (diff_low * diff_high > 0.0d0) then
     print *, 'No root found in the given interval.'
     stop
  end if

  ! Bisection method loop
  do iter = 1, max_iter
    a_mid = 0.5d0 * (a_low + a_high)
    diff_mid = MinerboFluxFactor(a_mid, tol) - f_cut

    ! Check if we've found the root
    if (abs(diff_mid) < tol) then
      a_root = a_mid
      return
    end if

    ! Update the bounds based on the sign of the function
    if (diff_mid * diff_low < 0.0d0) then
      a_high = a_mid
      diff_high = diff_mid
    else
      a_low = a_mid
      diff_low = diff_mid
    end if
  end do

  ! If maximum iterations are reached without finding the root
  print *, 'FindAlpha Bisection method did not converge.'
  write(*,*) abs(diff_mid - f_cut), f, f_cut, a_high, a_low, diff_mid
  stop

contains

  ! Function to calculate coth(x)
  function MinerboFluxFactor(x, tol) result(y)
    real*8, intent(in) :: x, tol
    real*8 :: y

    if ( abs(x) < tol) then
      y = 1.0d0 / tanh(tol) - 1.0d0/tol  ! Handle things below tolerance where tanh blows up
    else
      y = 1.0d0 / tanh(x) - 1.0d0/x
    end if

  end function MinerboFluxFactor

end subroutine FindAlpha

! This is for Maxwell Boltzmann approximation
subroutine FindEta(alpha, e, eta)
  real*8, intent(in) :: alpha, e !here e is the occupancy number (i.e. from 0 o 1)
  real*8, intent(out) :: eta

  eta = -log(e * alpha / sinh(alpha))
  
end subroutine FindEta