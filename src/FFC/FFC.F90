MODULE FFC

  USE GR1D_module, ONLY: &
    M1_moment_to_distro, q_M1, q_M1_fluid, &
    ghosts1, M1_imaxradii, number_groups, &
    number_species, Lmax, &
    GQ_weights, GQ_roots, &
    clite, pi4, pi, M1_moment_to_distro_inverse, &
    twothirds, fermi0, energy_gf, length_gf, &
    nulib_energy_gf, hbarc_mevcm, mev_to_erg, &
    Keep_FFC_P_fixed, FFC_Pn_value, X, alp, W, &
    hbarc_cgs, time_gf

  USE nulibtable, ONLY: &
    nulibtable_ewidths, nulibtable_energies

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyFastFlavorConversions
  LOGICAL, parameter :: UseFluidFrame = .true.

CONTAINS

  SUBROUTINE ApplyFastFlavorConversions(closure, method, dts)
    CHARACTER(LEN=64), INTENT(IN) :: closure, method
    REAL*8, INTENT(IN) :: dts

    IF (TRIM(ADJUSTL(method)) == 'Box3D') THEN
      CALL ApplyFFC_Box3D(closure, dts)
    ELSE IF (TRIM(ADJUSTL(method)) == 'FixedProbability') THEN
      CALL ApplyFFC_FixedProb()
    ELSE
      WRITE(*,*) 'FFC Method not supported'
      STOP
    ENDIF

  END SUBROUTINE ApplyFastFlavorConversions

  SUBROUTINE ApplyFFC_Box3D(closure, dts)

    CHARACTER(LEN=64), INTENT(IN) :: closure
    REAL*8, INTENT(IN) :: dts

    INTEGER :: i, j, k, l
    REAL*8 :: alpha(number_species), eta(number_species)
    REAL*8 :: mu, mom0, mom1, new_f, growth_rate, modulation_term
    REAL*8 :: ratio, N_prime, I_plus, I_minus
    REAL*8 :: IntegrationFactor(number_groups)
    REAL*8 :: G_of_n(Lmax), P_of_n(Lmax)
    LOGICAL :: Gamma_minus(Lmax)
    REAL*8 :: ang_distr(number_species, Lmax, number_groups)
    REAL*8 :: ang_distr_asy(number_species, Lmax, number_groups)
    REAL*8 :: q_before(number_species, number_groups, 2)

    LOGICAL :: debug_print = .false.
    
    IntegrationFactor(:) = 2.0d0*pi * SQRT(2.0d0) * fermi0 / hbarc_mevcm * clite &
         / (8.0d0 * pi**3) * nulibtable_ewidths(:) * &
         (nulibtable_energies(:)/nulib_energy_gf)**2 / time_gf


    !$OMP PARALLEL DO PRIVATE(k, ang_distr, ang_distr_asy, &
    !$OMP I_plus, I_minus, Gamma_minus, G_of_n, P_of_n, &
    !$OMP growth_rate, modulation_term, new_f, mom0, mom1, &
    !$OMP i, l, j, mu, alpha, eta, N_prime, ratio, q_before)
    DO k = ghosts1+1, M1_imaxradii-1

      ang_distr = 0.0d0
      ang_distr_asy = 0.0d0
      G_of_n = 0.0d0
      P_of_n = 1.0d0
      Gamma_minus = .false.
      I_plus = 0.0d0
      I_minus = 0.0d0
      
      IF (UseFluidFrame) THEN
        q_before(:,:,1:2) = q_M1_fluid(k,:,:,1:2)
      ELSE
        q_before(:,:,1:2) = q_M1(k,:,:,1:2)
      ENDIF
      DO j = 1, number_groups
        DO i = 1, number_species
          CALL ComputeDistribution(k, j, i, 0.0d0, closure, alpha(i), eta(i))
        ENDDO
        DO l = 1, Lmax
          mu = GQ_roots(l)
          ang_distr(1,l,j) = 1.0d0 / EXP(eta(1) - alpha(1) * mu)
          ang_distr(2,l,j) = 1.0d0 / EXP(eta(2) - alpha(2) * mu)
          ang_distr(3,l,j) = 1.0d0 / EXP(eta(3) - alpha(3) * mu)
          G_of_n(l) = G_of_n(l) + IntegrationFactor(j) * (ang_distr(1,l,j) - ang_distr(2,l,j))
        ENDDO
      ENDDO

      DO l = 1, Lmax
        IF (G_of_n(l) >= 0.0d0) THEN
          I_plus = I_plus + G_of_n(l) * GQ_weights(l)
        ELSE
          I_minus = I_minus - G_of_n(l) * GQ_weights(l)
          Gamma_minus(l) = .true.
        ENDIF
      ENDDO

      IF ((I_minus /= 0.0d0) .OR. (I_plus /= 0.0d0)) THEN
        CALL SurvivalProbability(I_minus, I_plus, Gamma_minus, P_of_n)
      ENDIF

      IF (Keep_FFC_P_fixed) P_of_n(:) = FFC_Pn_value

      DO j = 1, number_groups
        DO l = 1, Lmax
          ang_distr_asy(1,l,j) = P_of_n(l)*ang_distr(1,l,j) + (1.0d0 - P_of_n(l))*ang_distr(3,l,j)
          ang_distr_asy(2,l,j) = P_of_n(l)*ang_distr(2,l,j) + (1.0d0 - P_of_n(l))*ang_distr(3,l,j)
          ang_distr_asy(3,l,j) = 0.25d0*((1.0d0 - P_of_n(l))*ang_distr(1,l,j) + &
                                         (1.0d0 + P_of_n(l))*ang_distr(3,l,j) + &
                                         (1.0d0 - P_of_n(l))*ang_distr(2,l,j) + &
                                         (1.0d0 + P_of_n(l))*ang_distr(3,l,j))
        ENDDO
      ENDDO

      growth_rate = SQRT(I_plus * I_minus)
      modulation_term = MERGE(1.0d0, (1.0d0 - EXP(-growth_rate * dts)), Keep_FFC_P_fixed)

      DO j = 1, number_groups
        DO i = 1, number_species
          mom0 = 0.0d0
          mom1 = 0.0d0
          N_prime = 0.0d0
          DO l = 1, Lmax
            new_f = ang_distr(i,l,j) - modulation_term * (ang_distr(i,l,j) - ang_distr_asy(i,l,j))
            N_prime = N_prime + ang_distr(i,l,j) * GQ_weights(l)
            mom0 = mom0 + new_f * GQ_weights(l)
            mom1 = mom1 + new_f * GQ_weights(l) * GQ_roots(l)
          ENDDO

          IF ( i == 3 ) THEN
            mom0 = mom0 * 4.0d0
            mom1 = mom1 * 4.0d0
          ENDIF

          IF (UseFluidFrame) THEN
            ratio = q_M1_fluid(k,i,j,1) / (N_prime * M1_moment_to_distro_inverse(j) / 2.0d0)
            q_M1(k,i,j,1) = ratio * mom0 * M1_moment_to_distro_inverse(j) / 2.0d0
            q_M1(k,i,j,2) = ratio * mom1 * M1_moment_to_distro_inverse(j) / 2.0d0
          ELSE
            ratio = q_M1(k,i,j,1) / (N_prime * M1_moment_to_distro_inverse(j) / 2.0d0)
            q_M1(k,i,j,1) = ratio * mom0 * M1_moment_to_distro_inverse(j) / 2.0d0
            q_M1(k,i,j,2) = ratio * mom1 * M1_moment_to_distro_inverse(j) * X(k) / 2.0d0
          ENDIF

          IF (ratio < 0.98d0) WRITE(*,*) 'RATIO TOO SMALL', i, j, k, ratio
        ENDDO
      ENDDO
  
      IF (debug_print) THEN
        DO j = 1, number_groups
          DO i = 1, number_species
            IF ( (ABS(q_M1(k,i,j,1) - q_before(i,j,1))/q_before(i,j,1) ) > 1.0d-5 )THEN
              WRITE(*,*) '1', k, i, q_M1(k,i,j,1), q_before(i,j,1)
            ENDIF
            IF ( (ABS(q_M1(k,i,j,2) - q_before(i,j,2))/ABS(q_before(i,j,2)) ) > 1.0d-5 )THEN
              WRITE(*,*) '2', k, i, q_M1(k,i,j,2), q_before(i,j,2), ratio, X(k), &
                ABS(q_M1(k,i,j,2) - q_before(i,j,2))/ABS(q_before(i,j,2))
            ENDIF
          enddo
        enddo
      ENDIF

    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ApplyFFC_Box3D

  SUBROUTINE ApplyFFC_FixedProb()

    INTEGER :: i, j, k, l
    REAL*8  :: P_of_n(Lmax)
    LOGICAL :: Gamma_minus(Lmax)
    
    P_of_n(:) = FFC_Pn_value

    DO k = ghosts1+1, M1_imaxradii-1

      DO j = 1, number_groups
          q_M1(k,1,j,1:2) = P_of_n(l)*q_M1(k,1,j,1:2) + (1.0d0 - P_of_n(l))*q_M1(k,3,j,1:2)/2.0d0
          q_M1(k,2,j,1:2) = P_of_n(l)*q_M1(k,2,j,1:2) + (1.0d0 - P_of_n(l))*q_M1(k,3,j,1:2)/2.0d0
          q_M1(k,3,j,1:2) = 0.5d0*((1.0d0 - P_of_n(l))*q_M1(k,1,j,1:2) + &
                                   (1.0d0 - P_of_n(l))*q_M1(k,2,j,1:2) + &
                                   (1.0d0 + P_of_n(l))*q_M1(k,3,j,1:2))
      ENDDO
    ENDDO

  END SUBROUTINE ApplyFFC_FixedProb

  SUBROUTINE ComputeDistribution(k, j, i, mu, closure, alpha, eta)
    INTEGER, INTENT(IN) :: k, j, i
    REAL*8, INTENT(IN) :: mu
    REAL*8, INTENT(OUT) :: alpha, eta
    REAL*8, PARAMETER :: tol = 1.0d-8
    REAL*8 :: ThisConversion
    CHARACTER(LEN=64), INTENT(IN) :: closure
    REAL*8 :: e, f
    REAL*8 :: q_frame(2)

    IF ( i == 3) THEN
      ThisConversion = M1_moment_to_distro(j) / 4.0d0
    ELSE
      ThisConversion = M1_moment_to_distro(j)
    ENDIF

    IF (UseFluidFrame) THEN
      q_frame = q_M1_fluid(k, i, j, 1:2)
    ELSE
      q_frame = q_M1(k, i, j, 1:2)
      q_frame(2) = q_frame(2) / X(k)
    ENDIF

    e = ThisConversion * q_frame(1)
    IF (q_frame(2) > q_frame(1)) THEN
      WRITE(*,*) q_frame(2), q_frame(1), 'Flux too large in ComputeDistribution'
      STOP
    ENDIF
    f = q_frame(2) / q_frame(1)

    CALL FindAlpha(f, tol, alpha)
    CALL FindEta(alpha, e, eta)

  END SUBROUTINE ComputeDistribution

  SUBROUTINE SurvivalProbability(I_minus, I_plus, n_in_Gamma_minus, P)
    REAL*8, INTENT(IN)  :: I_minus, I_plus
    LOGICAL, INTENT(IN) :: n_in_Gamma_minus(Lmax)
    REAL*8, INTENT(OUT) :: P(Lmax)
    INTEGER :: l

    DO l = 1, Lmax
      IF (I_minus < I_plus) THEN
        P(l) = 1.0d0 / 3.0d0
        IF (.NOT. n_in_Gamma_minus(l)) P(l) = 1.0d0 - twothirds * I_minus / I_plus
      ELSE
        P(l) = 1.0d0 / 3.0d0
        IF (n_in_Gamma_minus(l)) P(l) = 1.0d0 - twothirds * I_plus / I_minus
      ENDIF
    ENDDO

  END SUBROUTINE SurvivalProbability

END MODULE FFC

