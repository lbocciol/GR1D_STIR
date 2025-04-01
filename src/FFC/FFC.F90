MODULE FFC

  USE GR1D_module, ONLY: &
    M1_moment_to_distro, q_M1, q_M1_fluid, &
    ghosts1, M1_imaxradii, number_groups, &
    number_species, Lmax, &
    GPQ_Lmax_weights, GPQ_Lmax_roots, &
    clite, pi4, pi, M1_moment_to_distro_inverse, &
    twothirds, fermi0, energy_gf, length_gf, &
    nulib_energy_gf, hbarc_mevcm, Keep_FFC_P_fixed, &
    FFC_Pn_value

  USE nulibtable, ONLY: &
    nulibtable_ewidths, nulibtable_energies

  IMPLICIT NONE
  PRIVATE 

  PUBLIC :: ApplyFastFlavorConversions

CONTAINS

  SUBROUTINE ApplyFastFlavorConversions(closure, method, dts)
    CHARACTER(LEN=64), INTENT(IN) :: closure, method
    REAL*8, INTENT(IN) :: dts

    IF (TRIM(ADJUSTL(method)) == 'Box3D') THEN
      CALL ApplyFFC_Box3D(closure, dts)
    ELSE
      WRITE(*,*) 'FFC Method not supported'
    ENDIF

  END SUBROUTINE ApplyFastFlavorConversions

  SUBROUTINE ApplyFFC_Box3D(closure, dts)
    CHARACTER(LEN=64), INTENT(IN) :: closure
    REAL*8, INTENT(IN) :: dts
    INTEGER :: i, j, k, l
    REAL*8  :: I_plus, I_minus
    REAL*8 :: alpha, eta 
    REAL*8 :: mu, mom0, mom1, new_f, growth_rate, modulation_term
    REAL*8 :: IntegrationFactor(number_groups)
    REAL*8 :: G_of_n(Lmax), P_of_n(Lmax)
    LOGICAL :: Gamma_minus(Lmax), Gamma_plus(Lmax)
    LOGICAL :: AlwaysSurvive = .false. ! This basically shuts off FFC (use for testing)
    LOGICAL :: debug_print = .false.
    REAL*8 :: extra_p(number_species,number_groups)
    REAL*8 :: ang_distr(number_species,Lmax,number_groups)
    REAL*8 :: ang_distr_asy(number_species,Lmax,number_groups)

    ! Precompute conversion factors in cgs units!
    IntegrationFactor(:) = SQRT(2.0d0) * fermi0 * (clite / hbarc_mevcm) &
         / (8.0d0 * pi**3) * nulibtable_ewidths(:) * &
         nulibtable_energies(:)**2 / nulib_energy_gf**2

    !$OMP PARALLEL DO PRIVATE(k, ang_distr, ang_distr_asy, I_plus, I_minus, Gamma_plus, Gamma_minus, &
    !$OMP G_of_n, P_of_n, growth_rate, new_f, mom0, mom1, i, l, j, mu, alpha, eta)
    !DO k = 140, 141
    DO k = ghosts1, M1_imaxradii
      I_plus = 0.0d0
      I_minus = 0.0d0
      Gamma_minus(:) = .false.
      Gamma_plus(:) = .false.
      G_of_n(:) = 0.0d0
      ang_distr(:,:,:) = 0.0d0
      ang_distr_asy(:,:,:) = 0.0d0

      DO i = 1, number_species
        DO j = 1, number_groups
          if (debug_print) WRITE(*,*) 'Before FFC', k, i, j, q_M1(k, i, j, 1),q_M1(k, i, j, 2)
        ENDDO
      ENDDO

      ! Compute G_of_n values
      DO l = 1, Lmax
        mu = GPQ_Lmax_roots(l)
        DO j = 1, number_groups

          ! Nue 
          CALL ComputeDistribution(k, j, 1, mu, closure, 1.0d0, alpha, eta, extra_p(1,j))
          ang_distr(1, l, j) = 1.0d0/EXP(eta - alpha * mu)
          if (debug_print) WRITE(*,*) 'angdis', l, j, eta, alpha, mu, ang_distr(1, l, j)

          ! Anue
          CALL ComputeDistribution(k, j, 2, mu, closure, 1.0d0, alpha, eta, extra_p(2,j))
          ang_distr(2, l, j) = 1.0d0/EXP(eta - alpha * mu)

          ! Nux
          CALL ComputeDistribution(k, j, 3, mu, closure, 4.0d0, alpha, eta, extra_p(3,j))
          ang_distr(3, l, j) = 1.0d0/EXP(eta - alpha * mu)

          ! Compute G_of_n. Notice that G_x = 0 since we have 3 species
          G_of_n(l) = G_of_n(l) + IntegrationFactor(j) * &
              (ang_distr(1, l, j) - ang_distr(2, l, j))

          if (debug_print) write(*,*) l, j, G_of_n(l)
        ENDDO
      ENDDO

      ! Compute I_plus and I_minus
      DO l = 1, Lmax
        IF (G_of_n(l) >= 0.0d0) THEN
          I_plus = I_plus + G_of_n(l) * GPQ_Lmax_weights(l)
          Gamma_plus(l) = .true.
        ELSE
          I_minus = I_minus - G_of_n(l) * GPQ_Lmax_weights(l)
          Gamma_minus(l) = .true.
        ENDIF
      ENDDO
      
      if (debug_print) write(*,*) I_plus, I_minus
      CALL SurvivalProbability(I_minus, I_plus, Gamma_minus, P_of_n)
      if (Keep_FFC_P_fixed) P_of_n(:) = FFC_Pn_value
      
      if (debug_print) WRITE(*,*) 'Probability'
      if (debug_print) WRITE(*,*) P_of_n(:)

      ! Compute asymptotic angular distributions
      DO j = 1, number_groups
        DO l = 1, Lmax
          ang_distr_asy(1, l, j) = P_of_n(l) * ang_distr(1, l, j) + (1.0d0 - P_of_n(l)) * ang_distr(3, l, j)
          ang_distr_asy(2, l, j) = P_of_n(l) * ang_distr(2, l, j) + (1.0d0 - P_of_n(l)) * ang_distr(3, l, j)
          ! Times four at the end because you have four heavy lepton neutrinos
          ang_distr_asy(3, l, j) = 0.25d0 * ((1.0d0 - P_of_n(l)) * ang_distr(1, l, j) + &
                                             (1.0d0 + P_of_n(l)) * ang_distr(3, l, j) + &
                                             (1.0d0 - P_of_n(l)) * ang_distr(2, l, j) + &
                                             (1.0d0 + P_of_n(l)) * ang_distr(3, l, j)) * 4.0d0
        ENDDO
      ENDDO

      growth_rate = SQRT(I_plus * I_minus)
      modulation_term = (1.0d0 - EXP(-growth_rate * dts))

      if (Keep_FFC_P_fixed) modulation_term = 0.0d0 
      
      if (debug_print) WRITE(*,*) 'sigma, t', growth_rate, dts
      DO i = 1, number_species
        DO j = 1, number_groups
          mom0 = 0.0d0
          mom1 = 0.0d0
          DO l = 1, Lmax
            new_f = ang_distr(i, l, j) - modulation_term * &
                    (ang_distr(i, l, j) - ang_distr_asy(i, l, j))
            new_f = ang_distr(i, l, j)

            if (abs(ang_distr(i, l, j) - new_f)/abs(new_f) > 1.0d-10) THEN
              if (debug_print) write(*,*) ang_distr(i, l, j), new_f, growth_rate
            endif

            mom0 = mom0 + new_f * GPQ_Lmax_weights(l)
            mom1 = mom1 + new_f * GPQ_Lmax_weights(l) * GPQ_Lmax_roots(l)
          ENDDO
          
          ! factor of two comes from angular terms
          q_M1(k, i, j, 1) = mom0 * M1_moment_to_distro_inverse(j) / 2.0d0
          ! add back extra momentum causing flux factor to be greater than one
          q_M1(k, i, j, 2) = mom1 * M1_moment_to_distro_inverse(j) / 2.0d0 + extra_p(i,j)
        ENDDO
      ENDDO

      DO i = 1, number_species
        DO j = 1, number_groups
          if (debug_print) WRITE(*,*) 'After FFC', k, i, j, q_M1(k, i, j, 1),q_M1(k, i, j, 2)
        ENDDO
      ENDDO

    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ApplyFFC_Box3D

  SUBROUTINE ComputeDistribution(k, j, i, mu, closure, factor, alpha, eta, extra_momentum)
    INTEGER, INTENT(IN) :: k, j, i
    REAL*8, INTENT(IN) :: mu, factor
    REAL*8, INTENT(OUT) :: alpha, eta, extra_momentum
    REAL*8, PARAMETER :: tol = 1.0d-8
    REAL*8 :: ThisConversion
    CHARACTER(LEN=64), INTENT(IN) :: closure
    
    REAL*8 :: e, f
    REAL*8 :: q_frame(2)
    INTEGER :: l
    LOGICAL :: UseFluidFrame = .false.

    ThisConversion = M1_moment_to_distro(j) / factor
    extra_momentum = 0.0d0

    IF (UseFluidFrame) THEN
      q_frame = q_M1_fluid(k, i, j, 1:2)
    ELSE
      q_frame = q_M1(k, i, j, 1:2)
    ENDIF

    e = ThisConversion * q_frame(1)
    IF (q_frame(2) > q_frame(1)) THEN
      extra_momentum = (q_frame(2) - q_frame(1))*1.01d0
    ENDIF
    f = (q_frame(2) - extra_momentum) / q_frame(1)

    IF (extra_momentum > 0.1d0*q_frame(2)) THEN
      WRITE(*,*) 'Too much extra momentum', k, i, j, q_frame(2), q_frame(1), extra_momentum
      STOP
    ENDIF

    CALL FindAlpha(f, tol, alpha)
    CALL FindEta(alpha, e, eta)

  END SUBROUTINE ComputeDistribution

  SUBROUTINE SurvivalProbability(I_minus, I_plus, n_in_Gamma_minus, P)

    REAL*8, INTENT(IN)  :: I_minus, I_plus
    LOGICAL, INTENT(IN) :: n_in_Gamma_minus(Lmax)
    REAL*8, INTENT(OUT) :: P(Lmax)
    INTEGER :: l
    REAL*8 :: ratio
  
    ! Precompute the ratio to avoid redundant calculations
    IF (I_minus > 1.0d-30) THEN
      ratio = twothirds * I_plus / I_minus
    ELSE
      ratio = 0.0d0
    ENDIF
  
    ! Vectorized loop over Lmax
    !$OMP SIMD
    DO l = 1, Lmax
      IF (I_minus < I_plus) THEN
        IF (n_in_Gamma_minus(l)) THEN
          P(l) = 1.0d0 / 3.0d0
        ELSE
          P(l) = 1.0d0 - ratio
        ENDIF
      ELSE
        IF (.not. n_in_Gamma_minus(l)) THEN
          P(l) = 1.0d0 / 3.0d0
        ELSE
          P(l) = 1.0d0 - ratio
        ENDIF
      ENDIF
    ENDDO
    !$OMP END SIMD
  
  END SUBROUTINE SurvivalProbability

END MODULE FFC
