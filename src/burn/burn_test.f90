program burn_test

  use pynet, only: nspec, spec_names
  use burn,  only: burn_init, burn_rhs, burn_state

  implicit none

  real(8) :: rho, T, dt
  real(8) :: X(nspec), X_out(nspec), dYdt(nspec)
  real(8) :: edot, e_step, e_total, Ye
  integer :: i, step, ierr

  rho = 1.0d8   ! g/cc
  T   = 5.0d8   ! K
  dt  = 0.1d0   ! s

  X      = 0.0d0
  X(1)   = 0.6d0   ! He4
  X(2)   = 0.2d0   ! C12
  X(3)   = 0.2d0   ! O16

  call burn_init()

  print *, "rho =", rho, "  T =", T
  print *, ""

  ! single RHS eval
  call burn_rhs(rho, T, X, dYdt, edot)
  print *, "edot [erg/g/s] =", edot
  do i = 1, nspec
    if (abs(dYdt(i)) > 0.0d0) &
      print '(a,a6,a,es12.4)', "  dYdt(", spec_names(i), ") =", dYdt(i)
  end do
  print *, ""

  ! time integration
  e_total = 0.0d0
  do step = 1, 5
    call burn_state(rho, T, X, dt, X_out, e_step, ierr, Ye)
    if (ierr /= 0) then
      print *, "burn_state failed at step", step
      stop 1
    end if
    e_total = e_total + e_step
    X = X_out
    print '(a,i2,a,es10.3,a,es12.4,a,f8.6)', &
      "step", step, "  e_step=", e_step, "  e_total=", e_total, "  Ye=", Ye
  end do

  print *, ""
  print *, "final X:"
  do i = 1, nspec
    if (X(i) > 1.0d-12) &
      print '(a,a6,a,es12.4)', "  ", spec_names(i), " =", X(i)
  end do
  print *, "sum(X) =", sum(X)

end program burn_test
