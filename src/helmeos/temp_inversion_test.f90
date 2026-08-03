program temp_inversion_test
  ! Sweep test for the Helmholtz rho,e -> T inversion (eos_input_re).
  !
  ! For a grid of (rho, T_true, Ye) we:
  !   1. evaluate rho-T mode to get e_true,
  !   2. feed e_true back in rho-e mode from a deliberately wrong T guess,
  !   3. check the recovered T against T_true and round-trip the energy.
  !
  ! The Helmholtz Newton solver (FullHelmEOS) silently returns the last
  ! temperature if it fails to converge, so this is the only way to catch
  ! inversion failures across the (rho,T,Ye) space.
  use wlHelmholtzEOS, only: ReadHelmTable, FullHelmEOS, HelmEOS, HelmholtzStateType, &
                            eos_input_rt, eos_input_re
  implicit none

  type(HelmholtzStateType) :: st

  ! grid (broad table coverage, ~1 decade inside the table edges)
  integer, parameter :: nT   = 40      ! T points,   1e4 .. 1e12 K  (log10)
  integer, parameter :: nrho = 40      ! rho points, 1e-4 .. 1e13 g/cc (log10)
  integer, parameter :: nye  = 5       ! Ye points,  0.1 .. 0.5 (linear)

  real(8), parameter :: logT_lo   =  4.0d0, logT_hi   = 12.0d0
  real(8), parameter :: logrho_lo = -4.0d0, logrho_hi = 13.0d0
  real(8), parameter :: ye_lo     =  0.1d0, ye_hi     =  0.5d0

  real(8), parameter :: abar = 12.0d0  ! fixed; zbar = Ye*abar so ye = Ye
  real(8), parameter :: tol  = 1.0d-6  ! pass/fail threshold on |dT|/T

  integer :: iT, irho, iye, npts, nfail, nshown, nskip
  real(8) :: T_true, rho, ye, e_true, T_rec, e_rec
  real(8) :: errT, errE
  real(8) :: maxErrT, maxErrE
  real(8) :: maxErrT_rho, maxErrT_T, maxErrT_ye
  real(8) :: maxErrE_rho, maxErrE_T, maxErrE_ye

  call ReadHelmTable('helm_table.dat')

  npts    = 0
  nfail   = 0
  nshown  = 0
  nskip   = 0
  maxErrT = 0.0d0
  maxErrE = 0.0d0
  maxErrT_rho = 0.0d0; maxErrT_T = 0.0d0; maxErrT_ye = 0.0d0
  maxErrE_rho = 0.0d0; maxErrE_T = 0.0d0; maxErrE_ye = 0.0d0

  write(*,'(a)') "=== Helmholtz rho,e -> T inversion sweep ==="
  write(*,'(a,i0,a,i0,a,i0,a,i0)') &
       "grid: nrho=", nrho, " nT=", nT, " nYe=", nye, " total=", nrho*nT*nye

  do iye = 1, nye
     ye = ye_lo + (ye_hi - ye_lo) * dble(iye-1) / dble(nye-1)

     do irho = 1, nrho
        rho = 10.0d0 ** ( logrho_lo + (logrho_hi - logrho_lo) &
                                      * dble(irho-1) / dble(nrho-1) )

        do iT = 1, nT
           T_true = 10.0d0 ** ( logT_lo + (logT_hi - logT_lo) &
                                          * dble(iT-1) / dble(nT-1) )

           ! ---- forward: rho,T -> e ----
           st%rho      = rho
           st%T        = T_true
           st%abar     = abar
           st%zbar     = ye * abar
           st%ye       = ye
           st%e_offset = 0.0d0
           call HelmEOS(eos_input_rt, st)
           e_true = st%e

           ! The Helmholtz internal energy is negative at low T / low rho
           ! (its zero-point); eos_input_re is undefined there (it STOPs on
           ! non-positive energy).  Such points are outside the inverse API's
           ! domain -- tally and skip rather than feed them in.
           if (e_true <= 0.0d0) then
              nskip = nskip + 1
              cycle
           end if

           ! ---- inverse: rho,e -> T (from a wrong initial guess) ----
           st%T = T_true * 0.3d0
           st%e = e_true
           call HelmEOS(eos_input_re, st)
           T_rec = st%T

           ! ---- energy round-trip at the recovered T ----
           st%T        = T_rec
           st%e_offset = 0.0d0
           call HelmEOS(eos_input_rt, st)
           e_rec = st%e

           errT = abs(T_rec - T_true) / T_true
           if (e_true /= 0.0d0) then
              errE = abs(e_rec - e_true) / abs(e_true)
           else
              errE = abs(e_rec - e_true)
           end if

           npts = npts + 1

           if (errT > maxErrT) then
              maxErrT = errT
              maxErrT_rho = rho; maxErrT_T = T_true; maxErrT_ye = ye
           end if
           if (errE > maxErrE) then
              maxErrE = errE
              maxErrE_rho = rho; maxErrE_T = T_true; maxErrE_ye = ye
           end if

           if (errT > tol) then
              nfail = nfail + 1
              if (nshown < 20) then
                 nshown = nshown + 1
                 write(*,'(a,es10.3,a,es10.3,a,f5.3,a,es10.3,a,es10.3,a,es10.3)') &
                      "  FAIL rho=", rho, " T=", T_true, " Ye=", ye, &
                      " | T_rec=", T_rec, " errT=", errT, " errE=", errE
              end if
           end if
        end do
     end do
  end do

  write(*,'(a)') "--------------------------------------------------"
  write(*,'(a,i0)')        "points tested      : ", npts
  write(*,'(a,i0,a)')      "points skipped     : ", nskip, "  (e<=0, outside eos_input_re domain)"
  write(*,'(a,i0)')        "points failing tol : ", nfail
  write(*,'(a,es12.4)')    "tolerance (|dT|/T) : ", tol
  write(*,'(a,es12.4,a,es10.3,a,es10.3,a,f5.3)') &
       "max |dT|/T         : ", maxErrT, &
       "  at rho=", maxErrT_rho, " T=", maxErrT_T, " Ye=", maxErrT_ye
  write(*,'(a,es12.4,a,es10.3,a,es10.3,a,f5.3)') &
       "max |de|/e         : ", maxErrE, &
       "  at rho=", maxErrE_rho, " T=", maxErrE_T, " Ye=", maxErrE_ye

  if (nfail == 0 .and. maxErrT <= tol) then
     write(*,'(a)') "ALL CHECKS PASSED"
  else
     write(*,'(a,i0,a)') "FAIL: ", nfail, " point(s) exceeded the tolerance"
     stop 1
  end if

end program temp_inversion_test
