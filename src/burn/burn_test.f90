program burn_test

  use pynet, only: nspec, spec_names, aion, zion
  use burn,  only: burn_init, burn_rhs, burn_state, burn_newton
  use composition, only: composition_init
  use nse,   only: nse_init, nse_solve, nse_solve_core

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

  ! =========================================================================
  ! NSE solver tests
  ! =========================================================================
  print *, ""
  print *, "========================================================"
  print *, " NSE solver tests"
  print *, "========================================================"

  ! bare network composition (no appended n,p): these checks probe the
  ! all-N=Z alpha chain, including the expected ye/=0.5 failure below
  call composition_init(.false.)
  call nse_init()

  ! (1) conservation + detailed balance at a representative handoff state
  call nse_check(1.0d8, 6.0d9, 0.5d0)

  ! (2) high-T limit: photodisintegration -> light nuclei (he4) dominate
  call nse_check(1.0d7, 1.2d10, 0.5d0)

  ! (3) high-rho limit: iron peak (Ni56-group) dominates
  call nse_check(1.0d9, 6.0d9, 0.5d0)

  ! (4) the all-alpha network has no charge degree of freedom (every species is
  !     N=Z), so any Ye /= 0.5 is unreachable and the solver MUST fail rather than
  !     return a wrong-Ye composition.
  call nse_solve(1.0d8, 6.0d9, 0.48d0, X_out, ierr)
  print *, ""
  print *, " unreachable-Ye check: nse_solve(ye=0.48) on the all-alpha network"
  if (ierr /= 0) then
    print '(a,i0,a)', "   OK: solver failed as expected (ierr=", ierr, ")"
  else
    print *, "  UNEXPECTED: solver returned success for an unreachable Ye"
  end if

  ! (5) any-Ye / general-network path, exercised with a synthetic species set
  call nse_synthetic()

contains

  ! Run nse_solve at (rho,T,ye) and report the physics checks:
  !   - baryon  conservation  sum A_i Y_i = 1
  !   - charge  conservation  sum Z_i Y_i = ye (pinned to 0.5 here)
  !   - detailed balance: NSE composition fed to rhs_f gives dY/dt ~ 0
  !   - dominant species
  subroutine nse_check(rho_in, T_in, ye_in)
    real(8), intent(in) :: rho_in, T_in, ye_in
    real(8) :: Ynse(nspec), Xnse(nspec), dY(nspec), edoT_eos_high
    real(8) :: Yrel(nspec), e_rel, dXmax
    real(8) :: mass_err, ye_nse, ydot_max, ydot_rel
    integer :: jerr, k, kmax, irelax
    logical :: ok

    print *, ""
    print '(a,es9.2,a,es9.2,a,f6.3)', " rho=", rho_in, "  T=", T_in, "  ye_in=", ye_in
    call nse_solve(rho_in, T_in, ye_in, Ynse, jerr)
    if (jerr /= 0) then
      print *, "  nse_solve FAILED, ierr =", jerr
      return
    end if

    Xnse     = Ynse * aion
    mass_err = abs(sum(aion*Ynse) - 1.0d0)
    ye_nse   = sum(zion*Ynse)

    ! detailed balance: dY/dt should ~vanish at the NSE composition
    call burn_rhs(rho_in, T_in, Xnse, dY, edoT_eos_high)
    ydot_max = maxval(abs(dY))
    ydot_rel = ydot_max / (maxval(abs(Ynse)) + 1.0d-30)

    print '(a,es10.3)', "   |sum A_i Y_i - 1|   =", mass_err
    print '(a,f10.6)',  "   Ye = sum Z_i Y_i    =", ye_nse
    print '(a,es10.3)', "   max|dY/dt| [mol/g/s]=", ydot_max
    print '(a,es10.3)', "   max|dY/dt|/max(Y)/s =", ydot_rel
    print '(a,es10.3)', "   edot @ NSE [erg/g/s]=", edoT_eos_high

    ! Cross-check: relax the Saha NSE composition through the network's own
    ! (backward-Euler) integrator to its true fixed point.  A correct NSE seed
    ! barely moves; the residual drift is the expected ideal-Saha vs. screened-
    ! network offset, which the network erases on handoff anyway.
    Yrel = Ynse
    do irelax = 1, 20
      call burn_newton(rho_in, T_in, Yrel, 1.0d-3, e_rel, ok)
      if (.not. ok) exit
    end do
    dXmax = maxval(abs(Yrel*aion - Xnse))
    print '(a,es10.3)', "   max|dX| after relax =", dXmax

    ! dominant species
    kmax = 1
    do k = 2, nspec
      if (Xnse(k) > Xnse(kmax)) kmax = k
    end do
    print '(a,a6,a,f8.5)', "   dominant: ", spec_names(kmax), "  X =", Xnse(kmax)
    print *, "   composition (X > 1e-3):"
    do k = 1, nspec
      if (Xnse(k) > 1.0d-3) &
        print '(a,a6,a,f9.5)', "     ", spec_names(k), " =", Xnse(k)
    end do
  end subroutine nse_check

  ! ===========================================================================
  ! Exercise the general (any-Ye) 2-D solver with a synthetic network that DOES
  ! have a charge degree of freedom: free n, free p, he4, and Si/Fe-peak nuclei
  ! spanning a range of Z/A.  Checks the qualitative regimes described in the
  ! papers (Seitenzahl et al. 2009, 2008):
  !   - exact baryon & charge conservation at arbitrary ye
  !   - ye = 0.5, high T : free-nucleon + he4 plasma, n <-> p symmetric
  !   - ye = 0.5, low  T : iron peak (ni56) dominates
  !   - ye < 0.5, low  T : neutron-rich Fe-peak nuclei favoured over ni56
  !   - ye > 0.5, low  T : a large free-PROTON fraction appears (proton-rich NSE)
  ! ===========================================================================
  subroutine nse_synthetic()
    integer, parameter :: ns = 9, nc = 4
    character(len=5) :: nm(ns)
    real(8) :: a(ns), z(ns), bind(ns), g(ns)
    real(8) :: Yn(ns), Xn(ns), dm(ns)
    real(8), parameter :: dmn = 8.0713181d0, dmh = 7.2889706d0
    character(len=56) :: label(nc)
    real(8) :: crho(nc), cT(nc), cye(nc), mass_err, ye_err
    integer :: c, k, jerr

    !            n      p      he4    si28    fe52    fe54    fe56    ni56    ni58
    nm   = [character(len=5):: "n","p","he4","si28","fe52","fe54","fe56","ni56","ni58"]
    a    = [ 1.d0,  1.d0,  4.d0,  28.d0,  52.d0,  54.d0,  56.d0,  56.d0,  58.d0]
    z    = [ 0.d0,  1.d0,  2.d0,  14.d0,  26.d0,  26.d0,  26.d0,  28.d0,  28.d0]
    ! atomic mass excesses [MeV] (AME) -> binding via Q = Z*dmh + N*dmn - dm
    dm   = [ 8.07132d0, 7.28897d0, 2.42492d0, -21.49279d0, -48.32982d0, &
            -56.25060d0, -60.60655d0, -53.90380d0, -60.22754d0]
    bind = z*dmh + (a - z)*dmn - dm
    g    = 1.0d0
    g(1) = 2.0d0;  g(2) = 2.0d0          ! free nucleons: 2J+1 = 2

    ! all temperatures >= 5 GK (the NSE regime this solver targets)
    label = [character(len=56) :: &
      "ye=0.5  high T : free nucleons + he4, n<->p symmetric", &
      "ye=0.5         : iron peak (ni56)",                     &
      "ye=0.46        : neutron-rich Fe-peak favoured",        &
      "ye=0.54        : proton-rich -> free protons"]
    crho = [ 1.0d7, 1.0d9, 1.0d9, 1.0d9 ]
    cT   = [ 9.0d9, 6.0d9, 6.0d9, 6.0d9 ]
    cye  = [ 0.5d0, 0.5d0, 0.46d0, 0.54d0 ]

    print *, ""
    print *, "========================================================"
    print *, " NSE general-network (any-Ye) tests  [synthetic species]"
    print *, "========================================================"

    do c = 1, nc
      call nse_solve_core(a, z, bind, g, crho(c), cT(c), cye(c), Yn, jerr)
      print *, ""
      print '(a)', " "//trim(label(c))
      print '(a,es9.2,a,es9.2,a,f5.2)', "   rho=",crho(c),"  T=",cT(c),"  ye=",cye(c)
      if (jerr /= 0) then
        print *, "   nse_solve_core FAILED, ierr =", jerr
        cycle
      end if
      Xn = Yn * a
      mass_err = abs(sum(a*Yn) - 1.0d0)
      ye_err   = abs(sum(z*Yn) - cye(c))
      print '(a,es10.3,a,es10.3)', "   |sum A_i Y_i - 1|=", mass_err, &
                                   "   |Ye - ye|=", ye_err
      do k = 1, ns
        if (Xn(k) > 1.0d-3) print '(a,a5,a,f9.5)', "     ", nm(k), " =", Xn(k)
      end do
    end do

  end subroutine nse_synthetic

end program burn_test
