MODULe FFC

  USE GR1D_module, ONLY: &
    M1_moment_to_distro, q_M1, q_M1, &
    ghosts1, M1_imaxradii, number_groups, &
    number_species, Lmax, &
    GPQ_Lmax_weights, GPQ_Lmax_roots, &
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
    LOGICAL :: Gamma_minus(Lmax) 
    LOGICAL :: AlwaysSurvive = .false. ! This basically shuts off FFC (use for testing)
    LOGICAL :: debug_print = .false.
    REAL*8 :: ang_distr_asy(number_species,Lmax,number_groups)
    REAL*8 :: ang_distr(number_species,Lmax,number_groups)
    REAL*8 :: q_before(number_species,number_groups,2)

    ! Precompute conversion factors in cgs units!
    IntegrationFactor(:) = SQRT(2.0d0) * fermi0 / (hbarc_mevcm)**3 &
         / (8.0d0 * pi**3) * nulibtable_ewidths(:) * &
         nulibtable_energies(:)**2 / nulib_energy_gf**2 * mev_to_erg
    
    ! Precompute conversion factors in natural units
    ! Frist factor of two pi is for the integration of G_of_n in dphi
    IntegrationFactor(:) = 2.0d0*pi * SQRT(2.0d0) * fermi0 / hbarc_mevcm * clite &
        / (8.0d0 * pi**3) * nulibtable_ewidths(:) * &
        (nulibtable_energies(:) / nulib_energy_gf)**2 &
        / time_gf ! because you need to convert dts

    !$OMP PARALLEL DO PRIVATE(k, ang_distr, ang_distr_asy, &
    !$OMP I_plus, I_minus, Gamma_minus, q_before, &
    !$OMP G_of_n, P_of_n, growth_rate, modulation_term, &
    !$OMP new_f, mom0, mom1, i, l, j, mu, alpha, eta)
    !DO k = 140, 141
    DO k = ghosts1+1, M1_imaxradii
      I_plus = 0.0d0
      I_minus = 0.0d0
      Gamma_minus(:) = .false.
      G_of_n(:) = 0.0d0
      P_of_n(:) = 0.0d0
      ang_distr(:,:,:) = 0.0d0
      ang_distr_asy(:,:,:) = 0.0d0

      !q_before(:,:,:) = q_M1(k,:,:,1:2)
      ! Compute G_of_n values
      DO l = 1, Lmax
        mu = GPQ_Lmax_roots(l)
        DO j = 1, number_groups

          ! Nue 
          CALL ComputeDistribution(k, j, 1, mu, closure, alpha, eta)
          ang_distr(1, l, j) = 1.0d0/EXP(eta - alpha * mu)
          !if (debug_print) WRITE(*,*) 'angdis', l, j, eta, alpha, mu, ang_distr(1, l, j)

          ! Anue
          CALL ComputeDistribution(k, j, 2, mu, closure, alpha, eta)
          ang_distr(2, l, j) = 1.0d0/EXP(eta - alpha * mu)

          ! Nux
          CALL ComputeDistribution(k, j, 3, mu, closure, alpha, eta)
          ang_distr(3, l, j) = 1.0d0/EXP(eta - alpha * mu)

          ! Compute G_of_n. Notice that G_x = 0 since we have 3 species
          G_of_n(l) = G_of_n(l) + IntegrationFactor(j) * &
              (ang_distr(1, l, j) - ang_distr(2, l, j))

        ENDDO
      ENDDO

      ! Compute I_plus and I_minus
      DO l = 1, Lmax
        IF (G_of_n(l) >= 0.0d0) THEN
          I_plus = I_plus + G_of_n(l) * GPQ_Lmax_weights(l)
          Gamma_minus(l) = .false.
        ELSE
          I_minus = I_minus - G_of_n(l) * GPQ_Lmax_weights(l)
          Gamma_minus(l) = .true.
        ENDIF
      ENDDO
     
      IF ( (I_minus == 0.0d0) .and. (I_plus == 0.0d0) ) THEN
        P_of_n(:) = 1.0d0
      ELSE
        CALL SurvivalProbability(I_minus, I_plus, Gamma_minus, P_of_n)
      ENDIF

      IF (Keep_FFC_P_fixed) P_of_n(:) = FFC_Pn_value
      
      ! Compute asymptotic angular distributions
      DO j = 1, number_groups
        DO l = 1, Lmax
          ang_distr_asy(1, l, j) = P_of_n(l) * ang_distr(1, l, j) + (1.0d0 - P_of_n(l)) * ang_distr(3, l, j)
          ang_distr_asy(2, l, j) = P_of_n(l) * ang_distr(2, l, j) + (1.0d0 - P_of_n(l)) * ang_distr(3, l, j)
          ang_distr_asy(3, l, j) = 0.25d0*((1.0d0 - P_of_n(l)) * ang_distr(1, l, j) + &
                                          (1.0d0 + P_of_n(l)) * ang_distr(3, l, j) + &
                                          (1.0d0 - P_of_n(l)) * ang_distr(2, l, j) + &
                                          (1.0d0 + P_of_n(l)) * ang_distr(3, l, j))
        ENDDO
      ENDDO

      growth_rate = SQRT(I_plus * I_minus)
      modulation_term = (1.0d0 - EXP(-growth_rate * dts))
        
      if (Keep_FFC_P_fixed) modulation_term = 1.0d0 

      DO i = 1, number_species
        DO j = 1, number_groups
          mom0 = 0.0d0
          mom1 = 0.0d0
          DO l = 1, Lmax
            new_f = ang_distr(i, l, j) - modulation_term * &
                    (ang_distr(i, l, j) - ang_distr_asy(i, l, j))

            mom0 = mom0 + new_f * GPQ_Lmax_weights(l)
            mom1 = mom1 + new_f * GPQ_Lmax_weights(l) * GPQ_Lmax_roots(l)
          ENDDO
          
          ! factor of two comes from angular terms
          q_M1(k, i, j, 1) = mom0 * M1_moment_to_distro_inverse(j) / 2.0d0
          ! add back extra momentum causing flux factor to be greater than one
          q_M1(k, i, j, 2) = mom1 * M1_moment_to_distro_inverse(j) * X(k) / 2.0d0
          
          !IF ( (ABS((q_M1(k,i,j,1) - q_before(i,j,1)) / q_before(i,j,1)) > 1.0d-2) .OR. &
          !    (ABS((q_M1(k,i,j,2) - q_before(i,j,2)) / q_before(i,j,2)) > 1.0d-2) .AND.&
          !    ABS(q_M1(k,i,j,2)/q_M1(k,i,j,1)) > 1.0d-6 ) THEN
          !    WRITE(*,*) 'E', k,i,j, q_M1(k,i,j,1), q_before(i,j,1)
          !    WRITE(*,*) 'F', k,i,j, q_M1(k,i,j,2), q_before(i,j,2), X(k), alp(k), W(k)
          !    STOP
          !ENDIF

        ENDDO
      ENDDO

    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE ApplyFFC_Box3D

  SUBROUTINE ComputeDistribution(k, j, i, mu, closure, alpha, eta)
    INTEGER, INTENT(IN) :: k, j, i
    REAL*8, INTENT(IN) :: mu
    REAL*8, INTENT(OUT) :: alpha, eta
    REAL*8, PARAMETER :: tol = 1.0d-8
    REAL*8 :: ThisConversion
    CHARACTER(LEN=64), INTENT(IN) :: closure
    
    REAL*8 :: e, f
    REAL*8 :: q_frame(2)
    INTEGER :: l
    LOGICAL, parameter :: UseFluidFrame = .false.

    ThisConversion = M1_moment_to_distro(j)

    IF (UseFluidFrame) THEN
      q_frame = q_M1(k, i, j, 1:2)
    ELSE
      q_frame = q_M1(k, i, j, 1:2)
      q_frame(2) = q_frame(2)/X(k)
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
