program helm_test
  ! Composite (nuc_eos + Helmholtz) temperature-inversion test.
  !
  ! This mirrors the production blend in eos.F90 (nuc_helm_eos_short /
  ! blend_weight) as a standalone round-trip check:
  !
  !   given two transition temperatures T_lo (= T_eos_low) and T_hi (= T_eos_high),
  !       T <= T_lo            -> pure Helmholtz
  !       T >= T_hi            -> pure nuc_eos
  !       T_lo < T < T_hi      -> linear blend of BOTH backends
  !
  ! For a grid of (rho, T_true, Ye) we:
  !   1. forward: e_true = e_blend(rho, T_true, Ye)
  !   2. inverse: feed e_true back and Newton on the piecewise blended energy
  !      from a perturbed T guess (regime decided purely by temperature)
  !   3. check the recovered T against T_true.
  !
  ! The Helmholtz and nuc_eos energies use different zero points, so the helm
  ! branch carries an offset chosen to match nuc_eos at the lower seam T_lo
  ! (the standalone analog of eos.F90's energy-offset table).  Forward and
  ! inverse use the SAME offset, so the round trip is self-consistent.
  use wlHelmholtzEOS, only: ReadHelmTable, HelmEOS, HelmholtzStateType, &
                            eos_input_rt, eos_input_re
  use eosmodule,      only: eos_rhomin, eos_rhomax, eos_tempmin, eos_tempmax, &
                            eos_yemin, eos_yemax
  implicit none

  external :: nuc_eos_short, readtable

  ! MeV <-> K (GR1D_module convention)
  real(8), parameter :: mev2k = 1.1604447522806d10

  ! ---- the two blend temperatures (Kelvin) -------------------------------
  ! Change these to test other windows; defaults match GR1D's T_eos_low/T_eos_high.
  real(8), parameter :: T_lo = 3.5d9    ! T_eos_low : below this -> pure Helmholtz
  real(8), parameter :: T_hi = 5.8d9    ! T_eos_high    : above this -> pure nuc_eos

  ! ---- grid: crosses both seams (pure-helm, blend, pure-nuc) -------------
  integer, parameter :: nT   = 41       ! T_true, 1e9 .. 6e10 K (log10)
  integer, parameter :: nrho = 3        ! rho,    1e8 .. 1e12 g/cc (log10)
  integer, parameter :: nye  = 3        ! Ye,     0.40 .. 0.50 (linear)

  real(8), parameter :: logT_lo   =  9.0d0, logT_hi   = 10.778d0  ! ~6e10 K
  real(8), parameter :: logrho_lo =  8.0d0, logrho_hi = 12.0d0
  real(8), parameter :: ye_lo     =  0.40d0, ye_hi    =  0.50d0

  real(8), parameter :: abar = 28.0d0   ! fixed; zbar = Ye*abar so ye = Ye
  real(8), parameter :: tol  = 1.0d-5   ! pass/fail threshold on |dT|/T

  integer :: iT, irho, iye, npts, nfail, nshown, nskip
  integer :: nH, nB, nN
  real(8) :: T_true, rho, ye, e_offset, e_true, T_rec, errT
  real(8) :: maxErrT, maxErrH, maxErrB, maxErrN
  real(8) :: maxErrT_rho, maxErrT_T, maxErrT_ye
  real(8) :: e_n_lo, e_h_lo0, e_off0
  character(len=1) :: region_f, region_i
  character(len=256) :: nuc_table

  ! nuc_eos table: command-line arg 1, else the symlink at the repo root
  if (command_argument_count() >= 1) then
     call get_command_argument(1, nuc_table)
  else
     nuc_table = '../../Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5'
  end if

  call ReadHelmTable('helm_table.dat')
  call readtable(trim(nuc_table))

  npts = 0; nfail = 0; nshown = 0; nskip = 0
  nH = 0; nB = 0; nN = 0
  maxErrT = 0.0d0; maxErrH = 0.0d0; maxErrB = 0.0d0; maxErrN = 0.0d0
  maxErrT_rho = 0.0d0; maxErrT_T = 0.0d0; maxErrT_ye = 0.0d0

  write(*,'(a)') "=== nuc_eos + Helmholtz blended T-inversion test ==="
  write(*,'(a,es10.3,a,es10.3,a)') "blend window: T_lo=", T_lo, " K   T_hi=", T_hi, " K"
  write(*,'(a,es10.3,a,es10.3)')   "nuc_eos rho range : ", eos_rhomin, " .. ", eos_rhomax
  write(*,'(a,es10.3,a,es10.3)')   "nuc_eos T   range : ", eos_tempmin*mev2k, " .. ", eos_tempmax*mev2k
  write(*,'(a,es10.3,a,es10.3)')   "nuc_eos Ye  range : ", eos_yemin, " .. ", eos_yemax
  write(*,'(a,i0,a,i0,a,i0,a,i0)') "grid: nrho=", nrho, " nT=", nT, " nYe=", nye, &
                                   " total=", nrho*nT*nye

  do iye = 1, nye
     ye = ye_lo + (ye_hi - ye_lo) * dble(iye-1) / dble(nye-1)

     do irho = 1, nrho
        rho = 10.0d0 ** ( logrho_lo + (logrho_hi - logrho_lo) &
                                      * dble(irho-1) / dble(nrho-1) )

        ! --- offset: match Helmholtz to nuc_eos at the lower seam T_lo ---
        call nuc_e (rho, T_lo/mev2k, ye, e_n_lo)
        call helm_e(rho, T_lo, ye, 0.0d0, e_h_lo0)
        e_offset = e_n_lo - e_h_lo0

        do iT = 1, nT
           T_true = 10.0d0 ** ( logT_lo + (logT_hi - logT_lo) &
                                          * dble(iT-1) / dble(nT-1) )

           ! Guard the rho,e->T Helmholtz STOP on negative internal energy:
           ! in the pure-Helmholtz regime the inverse subtracts the offset, so
           ! skip points whose raw Helmholtz energy is non-positive.
           if (T_true <= T_lo) then
              call helm_e(rho, T_true, ye, 0.0d0, e_off0)
              if (e_off0 <= 0.0d0) then
                 nskip = nskip + 1
                 cycle
              end if
           end if

           ! ---- forward: (rho,T) -> e_blend ----
           call forward(rho, T_true, ye, e_offset, e_true, region_f)

           ! ---- inverse: e_blend -> T, from a deliberately wrong T guess ----
           call invert(rho, e_true, ye, e_offset, 1.3d0*T_true, T_rec, region_i)

           errT = abs(T_rec - T_true) / T_true
           npts = npts + 1

           select case (region_i)
           case ('H'); nH = nH + 1; if (errT > maxErrH) maxErrH = errT
           case ('B'); nB = nB + 1; if (errT > maxErrB) maxErrB = errT
           case ('N'); nN = nN + 1; if (errT > maxErrN) maxErrN = errT
           end select

           if (errT > maxErrT) then
              maxErrT = errT
              maxErrT_rho = rho; maxErrT_T = T_true; maxErrT_ye = ye
           end if

           if (errT > tol) then
              nfail = nfail + 1
              if (nshown < 20) then
                 nshown = nshown + 1
                 write(*,'(a,es10.3,a,es10.3,a,f5.3,a,a,a,a,a,es10.3,a,es10.3)') &
                      "  FAIL rho=", rho, " T=", T_true, " Ye=", ye, &
                      " fwd=", region_f, " inv=", region_i, &
                      " | T_rec=", T_rec, " errT=", errT
              end if
           end if
        end do
     end do
  end do

  write(*,'(a)') "--------------------------------------------------"
  write(*,'(a,i0)')   "points tested      : ", npts
  write(*,'(a,i0,a)') "points skipped     : ", nskip, "  (raw helm e<=0 below T_lo)"
  write(*,'(a,i0)')   "points failing tol : ", nfail
  write(*,'(a,i0,a,es10.3)') "  pure Helmholtz   : ", nH, "  max |dT|/T=", maxErrH
  write(*,'(a,i0,a,es10.3)') "  blended          : ", nB, "  max |dT|/T=", maxErrB
  write(*,'(a,i0,a,es10.3)') "  pure nuc_eos     : ", nN, "  max |dT|/T=", maxErrN
  write(*,'(a,es12.4)') "tolerance (|dT|/T) : ", tol
  write(*,'(a,es12.4,a,es10.3,a,es10.3,a,f5.3)') &
       "max |dT|/T         : ", maxErrT, &
       "  at rho=", maxErrT_rho, " T=", maxErrT_T, " Ye=", maxErrT_ye

  if (nfail == 0 .and. maxErrT <= tol) then
     write(*,'(a)') "ALL CHECKS PASSED"
  else
     write(*,'(a,i0,a)') "FAIL: ", nfail, " point(s) exceeded the tolerance"
     stop 1
  end if

contains

  ! Helmholtz specific energy in rho-T mode [erg/g], with energy offset.
  subroutine helm_e(rho, T_K, ye, e_offset, e)
    real(8), intent(in)  :: rho, T_K, ye, e_offset
    real(8), intent(out) :: e
    type(HelmholtzStateType) :: st
    st%rho      = rho
    st%T        = T_K
    st%abar     = abar
    st%zbar     = ye * abar
    st%ye       = ye
    st%e_offset = e_offset
    call HelmEOS(eos_input_rt, st)
    e = st%e
  end subroutine helm_e

  ! nuc_eos specific energy in rho-T mode [erg/g]; T in MeV.
  subroutine nuc_e(rho, T_mev, ye, e)
    real(8), intent(in)  :: rho, T_mev, ye
    real(8), intent(out) :: e
    real(8) :: t, p, ent, cs2, dedt, dpde, dpdr, munu
    integer :: keyerr
    t = T_mev
    call nuc_eos_short(rho, t, ye, e, p, ent, cs2, dedt, dpde, dpdr, munu, &
                       1, keyerr, 1.0d-9)
  end subroutine nuc_e

  ! linear blend weight w(T_K) in [0,1] (Perego et al. 2015)
  function blend_w(T_K) result(w)
    real(8), intent(in) :: T_K
    real(8) :: w
    if (T_K >= T_hi) then
       w = 1.0d0
    else if (T_K <= T_lo) then
       w = 0.0d0
    else
       w = (T_K - T_lo) / (T_hi - T_lo)
    end if
  end function blend_w

  ! blended specific energy at a single T (both backends at the same T)
  subroutine blend_e(rho, T_K, ye, e_offset, e)
    real(8), intent(in)  :: rho, T_K, ye, e_offset
    real(8), intent(out) :: e
    real(8) :: w, e_h, e_n
    w = blend_w(T_K)
    if (w <= 0.0d0) then
       call helm_e(rho, T_K, ye, e_offset, e)
    else if (w >= 1.0d0) then
       call nuc_e(rho, T_K/mev2k, ye, e)
    else
       call helm_e(rho, T_K, ye, e_offset, e_h)
       call nuc_e (rho, T_K/mev2k, ye, e_n)
       e = w*e_n + (1.0d0-w)*e_h
    end if
  end subroutine blend_e

  ! forward: (rho,T) -> e, reporting which regime was used
  subroutine forward(rho, T_K, ye, e_offset, e, region)
    real(8), intent(in)  :: rho, T_K, ye, e_offset
    real(8), intent(out) :: e
    character(len=1), intent(out) :: region
    real(8) :: w
    w = blend_w(T_K)
    if (w <= 0.0d0) then
       region = 'H'
    else if (w >= 1.0d0) then
       region = 'N'
    else
       region = 'B'
    end if
    call blend_e(rho, T_K, ye, e_offset, e)
  end subroutine forward

  ! inverse: e -> T, mirroring eos.F90's nuc_helm_eos_short keytemp=0 path.
  ! The regime is decided FIRST, from the two SEAM energies
  !     e_lo = e_helm(rho,T_lo)   e_hi = e_nuc(rho,T_hi)
  ! (each backend's energy is monotone in T inside its own regime), and the
  ! inversion is then done WITHIN that regime: Helmholtz inverts itself below
  ! T_lo (eos_input_re), nuc_eos inverts itself above T_hi (findtemp), and only
  ! the window needs a root find -- which is bracketed by construction.
  subroutine invert(rho, e_t, ye, e_offset, T_g, T_K, region)
    real(8), intent(in)  :: rho, e_t, ye, e_offset, T_g
    real(8), intent(out) :: T_K
    character(len=1), intent(out) :: region
    real(8) :: e_lo, e_hi
    logical :: cold_ok, hot_ok

    call helm_e(rho, T_lo,        ye, e_offset, e_lo)
    call nuc_e (rho, T_hi/mev2k,  ye,           e_hi)

    cold_ok = (e_t <= e_lo)
    hot_ok  = (e_t >= e_hi)

    if (cold_ok .and. hot_ok) then
       ! inverted window (e_lo > e_hi): choose by the guess, as production does
       if (T_g <= T_lo) then
          region = 'H'
       else if (T_g >= T_hi) then
          region = 'N'
       else
          region = 'B'
       end if
    else if (cold_ok) then
       region = 'H'
    else if (hot_ok) then
       region = 'N'
    else
       region = 'B'      ! e_lo < e_t < e_hi: sign change guaranteed
    end if

    select case (region)
    case ('H')
       call helm_invert(rho, e_t, ye, e_offset, min(max(T_g,1.0d4),T_lo), T_K)
    case ('N')
       call nuc_invert (rho, e_t, ye, max(T_g,T_hi), T_K)
    case default
       call window_invert(rho, e_t, ye, e_offset, T_g, e_lo, e_hi, T_K)
    end select
  end subroutine invert

  ! rho,e -> T by Helmholtz' own Newton (offset stripped once, convergence on
  ! the step size in T -- no residual on an offset-carrying energy)
  subroutine helm_invert(rho, e_t, ye, e_offset, T_guess, T_K)
    real(8), intent(in)  :: rho, e_t, ye, e_offset, T_guess
    real(8), intent(out) :: T_K
    type(HelmholtzStateType) :: st
    st%rho      = rho
    st%T        = T_guess
    st%abar     = abar
    st%zbar     = ye * abar
    st%ye       = ye
    st%e_offset = e_offset
    st%e        = e_t
    call HelmEOS(eos_input_re, st)
    T_K = st%T
  end subroutine helm_invert

  ! rho,e -> T by nuc_eos' own findtemp
  subroutine nuc_invert(rho, e_t, ye, T_guess, T_K)
    real(8), intent(in)  :: rho, e_t, ye, T_guess
    real(8), intent(out) :: T_K
    real(8) :: t, e, p, ent, cs2, dedt, dpde, dpdr, munu
    integer :: keyerr
    t = T_guess / mev2k
    e = e_t
    call nuc_eos_short(rho, t, ye, e, p, ent, cs2, dedt, dpde, dpdr, munu, &
                       0, keyerr, 1.0d-9)
    T_K = t * mev2k
  end subroutine nuc_invert

  ! rho,e -> T inside the blend window: rtsafe-style Newton with bisection
  ! fallback on the bracket [T_lo,T_hi].  The relative-T exit and the ULP floor
  ! on the residual are the two criteria the production solver relies on: the
  ! blended energy carries the composition offset, so the residual alone can be
  ! quantised above any sensible tolerance and never converge.
  subroutine window_invert(rho, e_t, ye, e_offset, T_g, e_lo, e_hi, T_K)
    real(8), intent(in)  :: rho, e_t, ye, e_offset, T_g, e_lo, e_hi
    real(8), intent(out) :: T_K
    real(8) :: a, b, fa, fb, f, t, t_new, e, ep, em, dfdt, h
    integer :: it
    integer, parameter :: maxit = 100

    a = T_lo ; fa = e_lo - e_t
    b = T_hi ; fb = e_hi - e_t
    t = min(max(T_g, a), b)

    do it = 1, maxit
       call blend_e(rho, t, ye, e_offset, e)
       f = e - e_t
       if (abs(f) <= max(1.0d-10*abs(e_t), 8.0d0*spacing(abs(e)))) exit

       if (sign(1.0d0,f) == sign(1.0d0,fa)) then
          a = t ; fa = f
       else
          b = t ; fb = f
       end if
       if (b - a <= 1.0d-10*t) exit

       h = 1.0d-4 * t
       call blend_e(rho, t+h, ye, e_offset, ep)
       call blend_e(rho, t-h, ye, e_offset, em)
       dfdt = (ep - em) / (2.0d0*h)
       if (dfdt == 0.0d0) then
          t_new = 0.5d0*(a + b)
       else
          t_new = t - f/dfdt
          if (t_new <= a .or. t_new >= b) t_new = 0.5d0*(a + b)
       end if
       t = t_new
    end do
    T_K = t
  end subroutine window_invert

end program helm_test
