subroutine M1_updateeas
  use GR1D_module
  use nulibtable
  use omp_lib
  implicit none

  ! ===== Variables =====
  real*8 :: xrho, xtemp, xye, xeta, energy_x
  integer :: j, k, i, j_prime, keytemp, keyerr
  real*8 :: eosdummy(17)

  ! Spectra arrays
  real*8 :: tempspectrum(nM1,number_species, number_groups, number_eas)
  real*8 :: singlespecies_tempspectrum(nM1,number_groups, number_eas)

  real*8 :: inelastic_tempspectrum(nM1,number_species, number_groups, number_groups, 2)
  real*8 :: singlespecies_inelastic_tempspectrum(nM1,number_groups, number_groups, 2)

  real*8 :: epannihil_tempspectrum(nM1,number_species, number_groups, number_groups, 4)
  real*8 :: singlespecies_epannihil_tempspectrum(nM1,number_groups, number_groups, 4)

  real*8 :: blackbody_spectra(number_species, number_groups)

  real*8 :: t1, t2, t_eas, t_inel, t_ass, t_bb

  t_eas = 0.0d0
  t_inel = 0.0d0
  t_ass  = 0.0d0
  t_bb   = 0.0d0

  t1 = omp_get_wtime()

  CALL CatchOutsideTable()
  CALL CalculateEAS()
  CALL CalculateEta()
  CALL CalculateEmissivity()


     !$OMP PARALLEL DO PRIVATE(xrho,xtemp,xye,tempspectrum,singlespecies_tempspectrum, &
     !$OMP keytemp,keyerr,eosdummy,xeta,inelastic_tempspectrum,singlespecies_inelastic_tempspectrum, &
     !$OMP epannihil_tempspectrum,singlespecies_epannihil_tempspectrum,blackbody_spectra,energy_x,i,j,k,j_prime,t1,t2)
     do k=2, M1_imaxradii+ghosts1-1

        ! --------------------------------------------------------
        ! (1) Get density, temperature, electron fraction
        ! --------------------------------------------------------
        xrho  = rho(k)/rho_gf
        xtemp = temp(k)
        xye   = ye(k)

        ! --------------------------------------------------------
        ! (2) Interpolate EAS spectrum from nulibtable
        ! --------------------------------------------------------
        t1 = omp_get_wtime()
        tempspectrum = 0.0d0
        if (log10(xrho).lt.nulibtable_logrho_min) then
           tempspectrum = 0.0d0
        else
           if (log10(xtemp).lt.nulibtable_logtemp_min) stop "M1_update_eas: temp too low"
           if (xye.lt.nulibtable_ye_min) stop "M1_update_eas: ye too low"
           if (log10(xtemp).gt.nulibtable_logtemp_max) stop "M1_update_eas: temp too high"
           if (xye.gt.nulibtable_ye_max) stop "M1_update_eas: ye too high"

           if (number_species_to_evolve.eq.1) then
              call nulibtable_single_species_range_energy(xrho, xtemp, xye, 1, &
                   singlespecies_tempspectrum, number_groups, number_eas)
              tempspectrum(k,1,:,:) = singlespecies_tempspectrum(k,:,:)
           else if (number_species_to_evolve.eq.3) then
              call nulibtable_range_species_range_energy(xrho, xtemp, xye, tempspectrum, &
                   number_species, number_groups, number_eas)
           else
              stop "set up eas interpolation for this number of species"
           endif
        endif
        t2 = omp_get_wtime()
        t_eas = t_eas + t2 - t1

        ! Reset annihilation contributions if needed
        if (include_epannihil_kernels) then
           tempspectrum(k,3,:,2) = 0.0d0
           tempspectrum(k,3,:,1) = 0.0d0
        endif

        ! --------------------------------------------------------
        ! (3) Update EAS variables
        ! --------------------------------------------------------
        t1 = omp_get_wtime()
        do i=1,number_species
           do j=1,number_groups
              eas(j,k,i,2) = tempspectrum(k,i,j,2)
              eas(j,k,i,3) = tempspectrum(k,i,j,3)
           enddo
        enddo
        t2 = omp_get_wtime()
        t_ass = t_ass + t2 - t1

        ! --------------------------------------------------------
        ! (4) Compute blackbody emissivity (via EOS)
        ! --------------------------------------------------------
        t1 = omp_get_wtime()
        keytemp = 1
        keyerr  = 0
#if HAVE_NUC_EOS
        call nuc_eos_full(xrho, xtemp, xye, eosdummy(1), eosdummy(2), eosdummy(3), &
             eosdummy(4), eosdummy(5), eosdummy(6), eosdummy(7), eosdummy(8), &
             eosdummy(9), eosdummy(10), eosdummy(11), eosdummy(12), &
             eosdummy(13), elechem(k), eosdummy(15), eosdummy(16), eosdummy(17), &
             keytemp, keyerr, eos_rf_prec)
        if(keyerr.ne.0) then
           write(6,*) "EOS PROBLEM in inelastic scatter :"
           write(6,"(i5,1P3E18.9)") k,xrho,xtemp,xye
           stop "This is bad!"
        endif
#else
        stop "Need nuclear EOS for M1 transport"
#endif

        ! Fermi integrals
        xeta = (elechem(k)-eosdummy(17))/xtemp
        do j=1,number_groups
           energy_x = nulibtable_energies(j)/(nulib_energy_gf*xtemp)
           blackbody_spectra(1,j) = clite*(nulibtable_energies(j)/nulib_energy_gf)**3* &
                (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)/(1.0d0+exp(energy_x-xeta))
           blackbody_spectra(2,j) = clite*(nulibtable_energies(j)/nulib_energy_gf)**3* &
                (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)/(1.0d0+exp(energy_x+xeta))
           blackbody_spectra(3,j) = clite*(nulibtable_energies(j)/nulib_energy_gf)**3* &
                (mev_to_erg/(2.0d0*pi*hbarc_mevcm)**3)/(1.0d0+exp(energy_x))
        enddo

        ! Scale emissivity
        eas(:,k,1,1) = (blackbody_spectra(1,:)*tempspectrum(k,1,:,2)* &
             nulibtable_ewidths(:)/nulib_opacity_gf)*nulib_emissivity_gf
        eas(:,k,2,1) = (blackbody_spectra(2,:)*tempspectrum(k,2,:,2)* &
             nulibtable_ewidths(:)/nulib_opacity_gf)*nulib_emissivity_gf
        eas(:,k,3,1) = 4.0d0*(blackbody_spectra(3,:)*tempspectrum(k,3,:,2)* &
             nulibtable_ewidths(:)/nulib_opacity_gf)*nulib_emissivity_gf
        t2 = omp_get_wtime()
        t_bb = t_bb + t2 - t1

        ! --------------------------------------------------------
        ! (5) Inelastic scattering kernels
        ! --------------------------------------------------------
        if (include_Ielectron_imp.or.include_Ielectron_exp) then
           t1 = omp_get_wtime()
           ! Get electron chemical potential again
           keytemp = 1
           keyerr  = 0
#if HAVE_NUC_EOS
           call nuc_eos_full(xrho,xtemp,xye,eosdummy(1),eosdummy(2),eosdummy(3), &
                eosdummy(4),eosdummy(5),eosdummy(6),eosdummy(7),eosdummy(8), &
                eosdummy(9),eosdummy(10),eosdummy(11),eosdummy(12), &
                eosdummy(13),elechem(k),eosdummy(15),eosdummy(16),eosdummy(17), &
                keytemp,keyerr,eos_rf_prec)
           if(keyerr.ne.0) then
              write(6,*) "EOS PROBLEM in M1_updateeas.F90:"
              write(6,"(i5,1P3E18.9)") k,xrho,xtemp,xye
              stop "This is bad!"
           endif
#else
           stop "Need nuclear EOS for M1 transport"
#endif
           xtemp = temp(k)
           xeta  = elechem(k)/temp(k)

           ! Interpolation from nulib
           inelastic_tempspectrum = 0.0d0
           if (log10(xrho).ge.nulibtable_logrho_min) then
              if (number_species_to_evolve.eq.1) then
                 call nulibtable_inelastic_single_species_range_energy2(xtemp,xeta,1, &
                      singlespecies_inelastic_tempspectrum,number_groups,number_groups,2)
                 inelastic_tempspectrum(k,1,:,:,:) = singlespecies_inelastic_tempspectrum(k,:,:,:)
              else if (number_species_to_evolve.eq.3) then
                 call nulibtable_inelastic_range_species_range_energy2(xtemp,xeta, &
                      inelastic_tempspectrum,number_species,number_groups,number_groups,2)
              else
                 stop "set up eas interpolation for this number of species"
              endif
           endif
          t2 = omp_get_wtime()
          t_inel = t_inel + t2 - t1


           ! Store results
           t1 = omp_get_wtime()
           do i=1,number_species
              do j=1,number_groups
                 do j_prime=1,number_groups
                    ies(j_prime,j,k,i,:) = inelastic_tempspectrum(k,i,j,j_prime,:)
                 enddo
              enddo
           enddo
          t2 = omp_get_wtime()
          t_ass = t_ass + t2 - t1

        endif

        ! --------------------------------------------------------
        ! (6) Electron-positron annihilation kernels
        ! --------------------------------------------------------
        if (include_epannihil_kernels) then
           epannihil_tempspectrum = 0.0d0
           if (log10(xrho).ge.nulibtable_logrho_min) then
              i = 3   ! only for nux
              call nulibtable_epannihil_single_species_range_energy2(xtemp,xeta,i, &
                   singlespecies_epannihil_tempspectrum,number_groups,number_groups,4)
              epannihil_tempspectrum(k,i,:,:,:) = singlespecies_epannihil_tempspectrum(k,:,:,:)
           endif

           ! Store results
           do i=1,number_species
              do j=1,number_groups
                 do j_prime=1,number_groups
                    epannihil(j_prime,j,k,i,:) = epannihil_tempspectrum(k,i,j,j_prime,:)
                 enddo
              enddo
           enddo
        endif

     enddo
     !$OMP END PARALLEL DO


  ! ============================================================
  ! OTHER TEST CASE HANDLING
  ! ============================================================
  else if (M1_testcase_number.ge.2.and.M1_testcase_number.le.9) then
     ! already handled elsewhere
  else
     stop "add in eas updating code for this test case"
  endif

 write(*,*) 'TIMES', t_bb, t_ass, t_eas, t_inel

end subroutine M1_updateeas


subroutine CalculateEAS()

     do k=2, M1_imaxradii+ghosts1-1

        ! --------------------------------------------------------
        ! (1) Get density, temperature, electron fraction
        ! --------------------------------------------------------
        xrho  = rho(k)/rho_gf
        xtemp = temp(k)
        xye   = ye(k)

        ! --------------------------------------------------------
        ! (2) Interpolate EAS spectrum from nulibtable
        ! --------------------------------------------------------
        t1 = omp_get_wtime()
        tempspectrum = 0.0d0
        if (log10(xrho).lt.nulibtable_logrho_min) then
           tempspectrum = 0.0d0
        else
           if (log10(xtemp).lt.nulibtable_logtemp_min) stop "M1_update_eas: temp too low"
           if (xye.lt.nulibtable_ye_min) stop "M1_update_eas: ye too low"
           if (log10(xtemp).gt.nulibtable_logtemp_max) stop "M1_update_eas: temp too high"
           if (xye.gt.nulibtable_ye_max) stop "M1_update_eas: ye too high"

           if (number_species_to_evolve.eq.1) then
              call nulibtable_single_species_range_energy(xrho, xtemp, xye, 1, &
                   singlespecies_tempspectrum, number_groups, number_eas)
              tempspectrum(k,1,:,:) = singlespecies_tempspectrum(k,:,:)
           else if (number_species_to_evolve.eq.3) then
              call nulibtable_range_species_range_energy(xrho, xtemp, xye, tempspectrum, &
                   number_species, number_groups, number_eas)
           else
              stop "set up eas interpolation for this number of species"
           endif
        endif
        t2 = omp_get_wtime()
        t_eas = t_eas + t2 - t1

        ! Reset annihilation contributions if needed
        if (include_epannihil_kernels) then
           tempspectrum(k,3,:,2) = 0.0d0
           tempspectrum(k,3,:,1) = 0.0d0
        endif

        ! --------------------------------------------------------
        ! (3) Update EAS variables
        ! --------------------------------------------------------
        t1 = omp_get_wtime()
        do i=1,number_species
           do j=1,number_groups
              eas(j,k,i,2) = tempspectrum(k,i,j,2)
              eas(j,k,i,3) = tempspectrum(k,i,j,3)
           enddo
        enddo
        t2 = omp_get_wtime()
        t_ass = t_ass + t2 - t1
  

end subroutine CalculateEAS
