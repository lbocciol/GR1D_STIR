!-*-f90-* 
!purpose is to update q_M1(:,:,:,3), the Eddington tensor, and to
!computer higher order and other second moments to store in
!q_M1_extra(:,:,:,:).  We do this for the zone center and the
!reconstructed values

subroutine M1_closure

  use timers
  use GR1D_module
  implicit none

  ! --- local variables ---
  real*8 :: tol, err
  real*8 :: oneM1flux, oneM1eddy_guess
  real*8 :: alp2, invalp, invalp2, X2, invX, invX2, W2, v2, oneW, onev, onealp, oneX
  real*8 :: invX4, WvX, W2X2, WoneX, twoWX2v
  real*8 :: udown(2)
  real*8 :: littlehupdown(2,2)
  real*8 :: localJ, H2, Hup1, Hup2, Kthin22, Kthick22, Kup22
  real*8 :: ff2, ff3, ff4, chi, oldprrguess
  real*8 :: Lthin, Lthick, Luprrr, Ldownfupfr
  real*8 :: Wuprrr, Wdownfupfr

  ! Section 3: precomputed constant parts of the tensor contractions.
  ! The 2x2 sums split into a part that depends only on T(1,1),T(1,2),T(2,1)
  ! (fixed for a given flux) and a part proportional to T(2,2) = eddy*invX4
  ! (the only term that changes each iteration).
  ! This eliminates the ii/jj inner loops from the while loop entirely.
  !
  ! Naming convention: _c = constant part (flux-dependent, pre-while)
  !                    _e = eddy coefficient (multiplied by eddy each iteration)
  real*8 :: T11, T12           ! Tupmunu fixed elements
  real*8 :: localJ_c, localJ_e ! localJ = localJ_c + localJ_e * eddy
  real*8 :: Hup1_c,  Hup1_e   ! Hup1   = Hup1_c  + Hup1_e  * eddy
  real*8 :: Hup2_c,  Hup2_e   ! Hup2   = Hup2_c  + Hup2_e  * eddy
  ! H2 similarly: H2 = H2_cc + H2_ce*eddy + H2_ee*eddy^2
  ! but since Hup1/Hup2 are linear in eddy, H2 is quadratic — keep as scalar
  real*8 :: u1, u2             ! udown(1), udown(2) unpacked for clarity
  real*8 :: h11, h12, h21, h22 ! littlehupdown elements unpacked

  integer :: h, k, i, j
  integer :: count
  ! Section 3: per-(k,i,j,h) count stored in a full array so
  ! the max is taken correctly across all j, not just the last j.
  integer, allocatable :: count_khij(:,:,:,:)  ! (k,i,j,h)
  real*8  :: justshyofone
  real(8) :: t1, t2

  ! -------------------------------------------------------
  ! Local contiguous pack arrays (j=1..ngroups, h=1..3)
  ! Layout: (j,h) — j is fast index, contiguous in memory.
  ! Copy-in loads all 3 h-states in one j-pass.
  ! -------------------------------------------------------
  real*8, allocatable :: loc_flux (:,:)  ! (ngroups,3): F/E per h-state
  real*8, allocatable :: loc_eddy (:,:)  ! (ngroups,3): converged P_rr/E
  real*8, allocatable :: loc_chi  (:,:)  ! (ngroups,3): converged chi
  real*8, allocatable :: loc_E1   (:)    ! (ngroups): E center, write-back
  real*8, allocatable :: loc_E2   (:)    ! (ngroups): E minus,  write-back
  real*8, allocatable :: loc_E3   (:)    ! (ngroups): E plus,   write-back
  real*8, allocatable :: loc_Hup2 (:)    ! (ngroups): Hup(2) at conv, h=2
  real*8, allocatable :: loc_Kup22(:)    ! (ngroups): Kup(2,2) at conv, h=2
  real*8, allocatable :: loc_localJ(:)   ! (ngroups): localJ at conv, h=2

  ! -------------------------------------------------------
  ! Count tracking:
  !   count_khij(k,i,j,h) — exact iteration count per (k,i,j,h).
  !   closure_counts(k,i,j) — max over h, for post-loop reporting.
  ! Both local, printed as warnings only, no OMP races
  ! (each (k,i,j,h) written by exactly one thread).
  ! -------------------------------------------------------
  integer, allocatable :: closure_counts(:,:,:)    ! (k,i,j)   max over h
  integer :: max_count_thresh

  max_count_thresh = 20        ! warn if any zone needs >20 iterations

  tol           = 1.0d-8
  justshyofone  = 0.9999999999d0

  if (M1_testcase_number.eq.3) then
     q_M1m(:,:,:,3,1)      = justshyofone
     q_M1_extram(:,:,:,1,1)= justshyofone
     q_M1(:,:,:,3)         = justshyofone
     q_M1_extra(:,:,:,1)   = 0.0d0
     q_M1_extra(:,:,:,4)   = justshyofone
     q_M1_extra(:,:,:,2)   = justshyofone*q_M1(:,:,:,2)
     q_M1_extra(:,:,:,3)   = 0.0d0
     q_M1p(:,:,:,3,1)      = justshyofone
     q_M1_extrap(:,:,:,1,1)= justshyofone
     return
  endif

  ! --- allocate local arrays ---
  allocate(loc_flux  (number_groups, 3))
  allocate(loc_eddy  (number_groups, 3))
  allocate(loc_chi   (number_groups, 3))
  allocate(loc_E1    (number_groups))
  allocate(loc_E2    (number_groups))
  allocate(loc_E3    (number_groups))
  allocate(loc_Hup2  (number_groups))
  allocate(loc_Kup22 (number_groups))
  allocate(loc_localJ(number_groups))
  allocate(count_khij    (M1_imaxradii, number_species_to_evolve, number_groups, 3))
  allocate(closure_counts(M1_imaxradii, number_species_to_evolve, number_groups))
  count_khij     = 0
  closure_counts = 0

  CALL GetThisTime(t1)

  ! -------------------------------------------------------
  ! Section 3: precompute metric-dependent coefficients
  ! for the tensor contraction split.
  !
  ! The 2x2 sums localJ, Hup1, Hup2 are linear in T(2,2):
  !
  !   localJ = localJ_c + localJ_e * eddy
  !   Hup1   = Hup1_c   + Hup1_e  * eddy
  !   Hup2   = Hup2_c   + Hup2_e  * eddy
  !
  ! where the _c terms depend on T(1,1), T(1,2)=T(2,1) (flux-dependent,
  ! fixed for a given j before the while loop), and _e terms are the
  ! coefficients of T(2,2) = eddy*invX4, which are pure metric quantities
  ! (depend only on k — computed once here).
  !
  ! This lets us completely eliminate the ii/jj inner loops from the
  ! while loop, replacing them with 3 scalar FMAs.
  ! -------------------------------------------------------

  !$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,8) &
  !$OMP PRIVATE(alp2,onealp,invalp,invalp2,X2,oneX,invX,invX2,invX4, &
  !$OMP         WvX,W2X2,WoneX,twoWX2v, &
  !$OMP         W2,oneW,v2,onev,udown,littlehupdown,h,j, &
  !$OMP         u1,u2,h11,h12,h21,h22, &
  !$OMP         T11,T12,localJ_c,localJ_e,Hup1_c,Hup1_e,Hup2_c,Hup2_e, &
  !$OMP         oneM1flux,oneM1eddy_guess, &
  !$OMP         err,count,localJ,Hup1,Hup2,H2,Kup22,Kthin22,Kthick22, &
  !$OMP         ff2,ff3,ff4,chi,oldprrguess, &
  !$OMP         Lthin,Lthick,Luprrr,Ldownfupfr,Wuprrr,Wdownfupfr, &
  !$OMP         loc_flux,loc_eddy,loc_chi,loc_E1,loc_E2,loc_E3, &
  !$OMP         loc_Hup2,loc_Kup22,loc_localJ)
  do k = ghosts1+1, M1_imaxradii
    do i = 1, number_species_to_evolve

      ! -------------------------------------------------------
      ! Metric factors (depend only on k)
      ! -------------------------------------------------------
      if (GR) then
        alp2    = alp(k)*alp(k)
        onealp  = alp(k)
        invalp  = 1.0d0/alp(k)
        invalp2 = 1.0d0/alp2
        X2      = X(k)*X(k)
        oneX    = X(k)
        invX    = 1.0d0/X(k)
        invX2   = 1.0d0/X2
      else
        alp2    = 1.0d0
        onealp  = 1.0d0
        invalp  = 1.0d0
        invalp2 = 1.0d0
        X2      = 1.0d0
        oneX    = 1.0d0
        invX    = 1.0d0
        invX2   = 1.0d0
      endif

      if (v_order.eq.-1) then
        if (GR) then
          W2   = W(k)**2
          oneW = W(k)
          v2   = v(k)**2
          onev = v(k)
        else
          if (v1(k).gt.0.5d0) stop "very relativistic and you are assuming newtonian in transport"
          W2   = 1.0d0/(1.0d0-v1(k)**2)
          oneW = sqrt(W2)
          v2   = v1(k)**2
          onev = v1(k)
        endif
      else if (v_order.eq.0) then
        W2   = 1.0d0
        oneW = 1.0d0
        v2   = 0.0d0
        onev = 0.0d0
      else
        stop "add in vorder"
      endif

      invX4   = invX2 * invX2
      WvX     = oneW * onev * oneX
      WoneX   = oneW * oneX
      W2X2    = W2 * X2
      twoWX2v = 2.0d0 * oneW * X2 * onev * oneX

      udown(1) = -oneW * onealp
      udown(2) =  oneW * onev * oneX
      littlehupdown(1,1) = -v2*W2
      littlehupdown(2,1) =  udown(2)*invX2*udown(1)
      littlehupdown(1,2) = -udown(1)*invalp2*udown(2)
      littlehupdown(2,2) =  W2

      ! Unpack for readability and to help the compiler avoid
      ! indirect addressing inside the while loop
      u1  = udown(1);          u2  = udown(2)
      h11 = littlehupdown(1,1); h12 = littlehupdown(1,2)
      h21 = littlehupdown(2,1); h22 = littlehupdown(2,2)

      ! -------------------------------------------------------
      ! Precompute the eddy-coefficient (_e) parts of the
      ! contractions.  These depend only on metric (k-only).
      ! T(2,2) = eddy * invX4, so:
      !
      !   localJ_e = u2*u2*invX4          (coefficient of eddy in localJ)
      !   Hup1_e   = -u2*h11_col2*invX4  where col2 means jj=2 terms
      !            = -(u1*h(1,2) + u2*h(1,2))*invX4 ... expand fully:
      !
      ! Full expansion of the jj=2 column (T(1,2)=T(2,1)=T12 is flux-dep,
      ! only T(2,2)=eddy*invX4 is the varying part):
      !
      !   localJ_e = u2 * u2 * invX4
      !   Hup1_e   = -(u1*h12 + u2*h12) * invX4   [ii=1,jj=2 + ii=2,jj=2]
      !            wait — need to be careful. The sum is over ii AND jj.
      !            Only jj=2 AND ii=2 contributes to T(2,2):
      !
      !   localJ += u(ii)*u(jj)*T(ii,jj)
      !   T(2,2) term: ii=2,jj=2 → u2*u2*T22
      !   so localJ_e = u2*u2
      !
      !   Hup1 += -u(ii)*h(1,jj)*T(ii,jj)
      !   T(2,2) term: ii=2,jj=2 → -u2*h12*T22
      !   so Hup1_e = -u2*h12
      !
      !   Hup2 += -u(ii)*h(2,jj)*T(ii,jj)
      !   T(2,2) term: ii=2,jj=2 → -u2*h22*T22
      !   so Hup2_e = -u2*h22
      !
      ! All multiplied by invX4 (since T22 = eddy*invX4):
      ! -------------------------------------------------------
      localJ_e = u2 * u2 * invX4
      Hup1_e   = -u2 * h12 * invX4
      Hup2_e   = -u2 * h22 * invX4

      ! -------------------------------------------------------
      ! COPY-IN: one j-pass, all 3 h-states
      ! -------------------------------------------------------
      do j = 1, number_groups
        loc_flux(j,1) = q_M1m(k,i,j,2,1) / q_M1m(k,i,j,1,1)
        loc_eddy(j,1) = q_M1(k,i,j,3)
        loc_E2(j)     = q_M1m(k,i,j,1,1)

        loc_flux(j,2) = q_M1(k,i,j,2) / q_M1(k,i,j,1)
        loc_eddy(j,2) = q_M1(k,i,j,3)
        loc_E1(j)     = q_M1(k,i,j,1)

        loc_flux(j,3) = q_M1p(k,i,j,2,1) / q_M1p(k,i,j,1,1)
        loc_eddy(j,3) = q_M1(k,i,j,3)
        loc_E3(j)     = q_M1p(k,i,j,1,1)
      enddo

      ! -------------------------------------------------------
      ! COMPUTE: h outer, j inner.
      ! For each (h,j): converge the while loop using the
      ! pre-split tensor contractions — no ii/jj loops inside.
      ! -------------------------------------------------------
      do h = 1, 3
        do j = 1, number_groups

          oneM1flux       = loc_flux(j,h)
          oneM1eddy_guess = loc_eddy(j,h)
          err   = 1.0d0
          count = 0

          ! Fixed Tupmunu elements (flux-dependent, don't change in while loop)
          T11 = invalp2                        ! T(1,1) = 1/alp^2  (oneM1en=1)
          T12 = oneM1flux * invX2 * invalp    ! T(1,2) = T(2,1)

          ! Constant parts of contractions (T11 and T12=T21 terms only)
          ! Expand the ii,jj sums for T(1,1) and T(1,2)+T(2,1):
          !   localJ_c = u1*u1*T11 + (u1*u2 + u2*u1)*T12
          !   Hup1_c   = -(u1*h11 + u2*h21)*T11 - (u1*h12 + u2*h22)*T12
          !              wait: -(u(ii)*h(1,jj)*T(ii,jj)) summed:
          !              ii=1,jj=1: -u1*h11*T11
          !              ii=2,jj=1: -u2*h21*T11
          !              ii=1,jj=2: -u1*h12*T12
          !              ii=2,jj=2: -u2*h12*T12  ← this is the _e part handled above
          !              so _c gets: -(u1*h11+u2*h21)*T11 - u1*h12*T12
          !                          (only ii=2,jj=2 excluded → that's the _e part)
          !   Hup2_c: similarly
          !              ii=1,jj=1: -u1*h21*T11
          !              ii=2,jj=1: -u2*h21*T11  ← h(2,1)
          !              ii=1,jj=2: -u1*h22*T12
          !              ii=2,jj=2: excluded (_e part)
          !   Full expansion:
          !     Hup2_c = -(u1*h21 + u2*h21)*T11 - u1*h22*T12
          !
          ! Note T(2,1)=T(1,2)=T12 by symmetry.
          ! Full expansion (ii=2,jj=2 term excluded — that's the _e part):
          !
          ! localJ: u(1)*u(1)*T11 + u(1)*u(2)*T12 + u(2)*u(1)*T12
          !       = u1^2*T11 + 2*u1*u2*T12
          !
          ! Hup1: -(u1*h(1,1)*T11 + u2*h(1,1)*T11   <- jj=1 column
          !         + u1*h(1,2)*T12)                 <- jj=2,ii=1 (ii=2 excluded)
          !      = -(u1*h11 + u2*h11)*T11 ... wait: h(1,jj), so:
          !   ii=1,jj=1: -u1*h11*T11
          !   ii=2,jj=1: -u2*h21*T11   <- h(1,jj=1) = h11? No: h(row=1,col=jj)
          !   Careful: Hup(m) = -sum_ii sum_jj u(ii)*h(m,jj)*T(ii,jj)
          !   so h(m,jj) means row=m, col=jj:
          !   Hup1: m=1  → h(1,jj):  h11(jj=1), h12(jj=2)
          !   ii=1,jj=1: -u1*h11*T11
          !   ii=2,jj=1: -u2*h11*T11   ← h(1,1)=h11 for both ii
          !   ii=1,jj=2: -u1*h12*T12
          !   ii=2,jj=2: excluded (_e)
          !   Hup1_c = -(u1+u2)*h11*T11 - u1*h12*T12
          !
          !   Hup2: m=2  → h(2,jj):  h21(jj=1), h22(jj=2)
          !   ii=1,jj=1: -u1*h21*T11
          !   ii=2,jj=1: -u2*h21*T11
          !   ii=1,jj=2: -u1*h22*T12
          !   ii=2,jj=2: excluded (_e)
          !   Hup2_c = -(u1+u2)*h21*T11 - u1*h22*T12
          localJ_c = u1*u1*T11 + 2.0d0*u1*u2*T12
          Hup1_c   = -(u1+u2)*h11*T11 - u1*h12*T12
          Hup2_c   = -(u1+u2)*h21*T11 - u1*h22*T12

          ! -------------------------------------------------------
          ! WHILE LOOP: now just scalar arithmetic, no inner loops.
          ! localJ, Hup1, Hup2 computed as linear functions of eddy.
          ! H2 is quadratic in eddy (Hup1, Hup2 are linear) — kept scalar.
          ! GCC can vectorize over j by peeling this into a j-loop version;
          ! for now the while loop is scalar but the data is cache-hot.
          ! -------------------------------------------------------
          do while (err > tol)
            count = count + 1

            ! Linear update — replaces the entire ii/jj double loops
            localJ = localJ_c + localJ_e * oneM1eddy_guess
            Hup1   = Hup1_c   + Hup1_e  * oneM1eddy_guess
            Hup2   = Hup2_c   + Hup2_e  * oneM1eddy_guess

            ! H2 = -alp2 * h(1,jj)*Hup1*Hup1_or_2 + X2 * h(2,jj)*Hup2*...
            ! Expanded (ii=1 terms get -alp2, ii=2 terms get +X2):
            !   ii=1: -alp2*(h11*Hup1*Hup1 + h12*Hup1*Hup2)
            !   ii=2: +X2* (h21*Hup2*Hup1 + h22*Hup2*Hup2)
            H2 = -alp2*(h11*Hup1*Hup1 + h12*Hup1*Hup2) &
                 +X2  *(h21*Hup2*Hup1 + h22*Hup2*Hup2)

            ! ff2: flux factor squared (unchanged formula, now uses scalar Hup1/Hup2)
            ff2 = ((onealp*WvX*Hup1)**2 + (WoneX*Hup2)**2 - &
                   2.0d0*onealp*W2*oneX*onev*Hup2*Hup1) / localJ**2 * 0.999999999999d0

            ! NaN guard: only needed if localJ=0, which is a degenerate zone.
            ! Moved the check here with early-exit to avoid branch every iter.
            if (localJ .eq. 0.0d0) then
              ff2 = 0.0d0
            else if (ff2 .ne. ff2) then
              write(*,*) oneM1flux, nt, k, i, j
              stop "ff2 is NaNing (localJ nonzero)"
            endif

            ! Closure relation
            if (M1closure .eq. 'ME') then
              ff3 = ff2*sqrt(ff2)
              ff4 = ff2*ff2
              chi = 1.0d0/3.0d0 + (3.0d0*ff2 - ff3 + 3.0d0*ff4)*0.1333333333333333333d0
            else if (M1closure .eq. 'LP') then
              if (ff2 .gt. 1.0d0) ff2 = 1.0d0
              chi = (3.0d0+4.0d0*ff2)/(5.0d0+2.0d0*sqrt(4.0d0-3.0d0*ff2))
            else
              stop "define closure"
            endif

            ! Kthin(2,2): thin-limit pressure tensor component
            ! Branch on Hup1 first (optically thin free-streaming direction)
            if (Hup1 .eq. 0.0d0) then
              Kthin22 = localJ * invX2
            else if (Hup2 .eq. 0.0d0 .or. H2 .eq. 0.0d0) then
              Kthin22 = 0.0d0
            else
              Kthin22 = Hup2*Hup2*localJ / H2
            endif

            ! Kthin NaN check — H2=0 path above already handles the
            ! degenerate case, so NaN here is a genuine error
            if (Kthin22 .ne. Kthin22) then
              write(*,*) k, i, j, h
              write(*,*) Kthin22, "kthin22", oneM1flux, H2, Hup1, Hup2, oneM1eddy_guess
              stop "Kthin22 is NaN"
            endif

            ! Kthick(2,2): thick-limit (isotropic) pressure; h22=W2
            Kthick22 = localJ * onethird * h22 * invX2

            ! Interpolated pressure tensor
            Kup22 = (3.0d0*chi-1.0d0)*0.5d0*Kthin22 + &
                    3.0d0*(1.0d0-chi)*0.5d0*Kthick22

            ! Update eddy guess: T(2,2)*X^4, simplified (oneM1en=1)
            oldprrguess     = oneM1eddy_guess
            oneM1eddy_guess = (W2X2*v2*localJ + twoWX2v*Hup2 + X2*X2*Kup22) * invX4 * X2*X2

            if (oneM1eddy_guess .ne. oneM1eddy_guess) then
              write(*,*) "eddy is NaN", localJ, Hup2, Kup22, Kthin22, Kthick22, chi
              stop
            endif

            err = abs(oneM1eddy_guess - oldprrguess) / oneM1eddy_guess

            ! --- fallback paths (logic preserved exactly) ---
            if (count .gt. 900) then
              if (oneM1flux .lt. -0.99d0) then
                oneM1flux = 0.7d0
                ! flux reset changes T12 and therefore _c terms — recompute
                T12      = oneM1flux * invX2 * invalp
                localJ_c = u1*u1*T11 + 2.0d0*u1*u2*T12
                Hup1_c   = -(u1*h11 + u2*h21)*T11 - u1*h12*T12
                Hup2_c   = -(u1*h21 + u2*h21)*T11 - u1*h22*T12
                write(*,*) "closure is failing and flux is -1, reset to 0.7d0"
              endif
            endif

            if (count .gt. 1000) then
              if (err .lt. 1.0d-5) then
                write(*,*) "accepting lower error", err, k, i, j, h, onev, oneM1eddy_guess, Hup2
                err = 1.0d-9
              else
                if (chi .gt. 0.9d0) then
                  write(*,*) "forcing p_rr to X^2", k, chi
                  write(*,*) k, i, j, h, oneM1flux, oneM1flux/oneX, q_M1(k,i,j,3), X2
                  oneM1eddy_guess = X2
                  err = 1.0d-9
                else
                  write(*,*) "Error at", count, ":", err, k, i, j, h, chi
                  if (count .gt. 1020) then
                    write(*,*) "If this is happening close to BH formation (or late stages when", &
                               " shock has receded a lot, consider using higher order explicit flux", &
                               " calculation"
                    stop "Closure is not converging after 1000 iterations"
                  endif
                endif
              endif
            endif

          enddo  ! while

          ! Store converged values into local arrays
          loc_eddy(j,h) = oneM1eddy_guess
          loc_chi (j,h) = chi

          ! h=2 extras needed for write-back
          if (h .eq. 2) then
            loc_Hup2  (j) = Hup2
            loc_Kup22 (j) = Kup22
            loc_localJ(j) = localJ
          endif

          ! Exact count for this (k,i,j,h) — stored directly, no last-j bug
          count_khij(k,i,j,h) = count

        enddo  ! j
      enddo    ! h

      ! -------------------------------------------------------
      ! COPY-OUT: one j-pass, all results
      ! -------------------------------------------------------
      do j = 1, number_groups

        ! h=1 (minus)
        q_M1m(k,i,j,3,1)       = loc_eddy(j,1)
        q_M1_extram(k,i,j,1,1) = loc_chi (j,1)

        ! h=2 (center) — extra quantities
        chi    = loc_chi  (j,2)
        Hup2   = loc_Hup2 (j)
        Kup22  = loc_Kup22(j)
        localJ = loc_localJ(j)

        if (Hup2 .eq. 0.0d0) then
          Lthin = 0.0d0
        else
          Lthin = Hup2 * invX2
        endif
        Lthick     = 0.6d0  * Hup2 * h22 * invX2     ! h22 = littlehupdown(2,2) = W2
        Luprrr     = (3.0d0*chi-1.0d0)*0.5d0*Lthin + 3.0d0*(1.0d0-chi)*0.5d0*Lthick
        Ldownfupfr = 3.0d0*(1.0d0-chi)*0.5d0*0.2d0*Hup2

        Wuprrr     = Luprrr + 3.0d0*oneW*onev*invX*Kup22 + &
                     3.0d0*W2*v2*invX2*Hup2 + W2*oneW*v2*onev*invX2*invX*localJ
        Wdownfupfr = oneW*onev*invX*3.0d0*(1.0d0-chi)*0.5d0*localJ*onethird + Ldownfupfr

        q_M1(k,i,j,3)       = loc_eddy(j,2)
        q_M1_extra(k,i,j,1) = 3.0d0*(1.0d0-chi)*0.5d0*localJ*onethird
        q_M1_extra(k,i,j,4) = chi
        q_M1_extra(k,i,j,2) = Wuprrr     * loc_E1(j)
        q_M1_extra(k,i,j,3) = Wdownfupfr * loc_E1(j)

        ! h=3 (plus)
        q_M1p(k,i,j,3,1)       = loc_eddy(j,3)
        q_M1_extrap(k,i,j,1,1) = loc_chi (j,3)

        ! Max count over h-states for this zone (correct: uses full array)
        closure_counts(k,i,j) = max(count_khij(k,i,j,1), &
                                    count_khij(k,i,j,2), &
                                    count_khij(k,i,j,3))
      enddo  ! j copy-out

    enddo  ! i
  enddo    ! k
  !$OMP END PARALLEL DO

  CALL GetThisTime(t2)
  timer_M1_clo = timer_M1_clo + (t2 - t1)

  ! -------------------------------------------------------
  ! Post-loop: print convergence warnings for slow zones.
  ! Done serially after OMP region so no race conditions.
  ! -------------------------------------------------------
  do k = ghosts1+1, M1_imaxradii
    do i = 1, number_species_to_evolve
      do j = 1, number_groups
        if (closure_counts(k,i,j) .gt. max_count_thresh) then
          write(*,'(A,3(A,I4),A,I6)') &
            "[M1_closure] slow convergence: ", &
            " k=", k, " i=", i, " j=", j, &
            " max_count=", closure_counts(k,i,j)
        endif
      enddo
    enddo
  enddo

  deallocate(loc_flux, loc_eddy, loc_chi)
  deallocate(loc_E1, loc_E2, loc_E3)
  deallocate(loc_Hup2, loc_Kup22, loc_localJ)
  deallocate(count_khij, closure_counts)

end subroutine M1_closure
