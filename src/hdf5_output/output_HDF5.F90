!-*-f90-*-
subroutine output_all_HDF5(modeflag)

  use GR1D_module
  use hdf5_output_utils
  use nulibtable
  implicit none

  character*1024 filename
  character*256 basename
  integer modeflag

  integer i,tempghosts
  real*8, allocatable :: cs(:)
  real*8 alp_ana(n1),rho_ana(n1),vel_ana(n1),X_ana(n1)
  real*8 maxr
  
  integer, parameter :: nscalars0 = 256
  real*8 scalars(nscalars0)
  integer nscalars
  integer ishock_radius(1)

  real*8 luminosity(number_species)
  real*8 luminosity_fluid(number_species)
  real*8 luminosity_rad(n1,number_species,2)
  
  real*8 enden_rad(n1,number_species,2)
  real*8 fluxden_rad(n1,number_species,2)
  
  real*8 num_luminosity(number_species)
  real*8 num_luminosity_fluid(number_species)
  real*8 num_luminosity_rad(n1,number_species,2)

  real*8 average_energy(number_species)
  real*8 average_energy_fluid(number_species)
  real*8 average_energy_rad(n1,number_species,2)

  real*8 rms_energy(number_species)
  real*8 rms_energy_fluid(number_species)
  real*8 rms_energy_rad(n1,number_species,2)

  real*8 spectrum(number_groups)

  real*8 fluxfactor_enweighted(n1,number_species)
  real*8 fluxfactor_fluxweighted(n1,number_species)
  real*8 eddingtonfactor_enweighted(n1,number_species)
  real*8 eddingtonfactor_fluxweighted(n1,number_species)
  real*8 total_nu_energy,total_matter_energy,total_energy
  
  integer k,j

  real*8 tau(n1,2) !scattering+absorption/absorption
  real*8 ave_tau(n1,number_species)
  logical nusphere_not_found(number_species,number_groups,2)
  logical ave_nusphere_not_found(number_species)
  real*8 nusphere(number_species,number_groups,2)
  real*8 ave_nusphere(number_species,5)

  integer keyerr,keytemp
  real*8 eosdummy(14)

  call hdf5_initialize()

  if(modeflag.eq.0) then
     
  else if(modeflag.eq.1) then
     
     ! Open xg.h5 file ONCE at the beginning
     call hdf5_open_xg_file()

     do k=ghosts1+1,n1-ghosts1
        keyerr = 0
        keytemp = 0

        ! This is to calculate things like mass fractions and chemical
        ! potentials
        call ApplyEOS_limits
        call eos_full(k,rho(k),temp(k),ye(k),eps(k),press(k),pressth(k), &
          entropy(k),cs2(k),eosdummy(2),&
          eosdummy(3),eosdummy(4),massfrac_a(k),massfrac_h(k), &
          massfrac_n(k),massfrac_p(k),massfrac_abar(k),massfrac_zbar(k), &
          elechem(k),eosdummy(12),eosdummy(13),eosdummy(14), &
          keytemp,keyerr,eoskey,eos_rf_prec)
     enddo

     ! Write time and grid coordinates
     call hdf5_write_grid_data_1d("time", [time], 1)
     
     ! Write velocity
     if (.not.small_output) call hdf5_write_grid_data_1d("v1", v1(1:n1)*clite, n1)
     
     if(do_rotation) then
        call hdf5_write_grid_data_1d("omega", omega(1:n1)*time_gf, n1)
        call hdf5_write_grid_data_1d("ToverW", ToverW(1:n1), n1)
        
        if(GR) then
           call hdf5_write_grid_data_1d("vphi", vphi(1:n1), n1)
        else
           call hdf5_write_grid_data_1d("vphi1", vphi1(1:n1), n1)
        endif
     endif

     if (do_turbulence) then
        call hdf5_write_grid_data_1d("omega2_BV", omega2_BV(1:n1)*time_gf**2, n1)
        call hdf5_write_grid_data_1d("v_turb", v_turb(1:n1)*time_gf/length_gf, n1)
        if (.not. small_output) then
          call hdf5_write_grid_data_1d("dissipated_turb_eps", diss(1:n1)*time_gf/eps_gf, n1)
          call hdf5_write_grid_data_1d("buoyancy_turb_eps", buoy(1:n1)*time_gf/eps_gf, n1)
          call hdf5_write_grid_data_1d("shear_turb_eps", shear(1:n1)*time_gf/eps_gf, n1)
          call hdf5_write_grid_data_1d("Lambda_MLT", lambda_mlt(1:n1)/length_gf, n1)
        endif
     endif   

     call hdf5_write_grid_data_1d("rho", rho(1:n1)/rho_gf, n1)
     
     if (do_nupress.or.do_M1) then
        if (.not.small_output) call hdf5_write_grid_data_1d("nuchem", nuchem(1:n1), n1)
        if (.not.small_output) call hdf5_write_grid_data_1d("press_nu", press_nu(1:n1)/press_gf, n1)
        if (.not.small_output) call hdf5_write_grid_data_1d("energy_nu", energy_nu(1:n1)/press_gf, n1)
        if (.not.small_output) call hdf5_write_grid_data_1d("dnupdr", dnupdr(1:n1)/press_gf*length_gf, n1)
     endif

     call hdf5_write_grid_data_1d("ye", ye(1:n1), n1)
     if (do_M1 .and. (.not.small_output)) then
        call hdf5_write_grid_data_1d("dyedt_hydro", dyedt_hydro(1:n1)*time_gf, n1)
        call hdf5_write_grid_data_1d("depsdt", depsdt(1:n1), n1)
        call hdf5_write_grid_data_1d("ynu", ynu(1:n1), n1)
     endif
     
     call hdf5_write_grid_data_1d("press", press(1:n1)/press_gf, n1)
     
     if (.not.small_output) then
       call hdf5_write_grid_data_1d("eps", eps(1:n1)/eps_gf, n1)
     endif

     if (GR) then
        call hdf5_write_grid_data_1d("mass_grav", mgrav(1:n1)/mass_gf, n1)
        call hdf5_write_grid_data_1d("mass_bary", mass(1:n1)/mass_gf, n1)
     else
        call hdf5_write_grid_data_1d("mass_bary", mass(1:n1)/mass_gf, n1)
     endif
     
     if (eoskey.eq.3) then
        if (.not.small_output) then
           call hdf5_write_grid_data_1d("xn", massfrac_n(1:n1), n1)
           call hdf5_write_grid_data_1d("xp", massfrac_p(1:n1), n1)
           call hdf5_write_grid_data_1d("xa", massfrac_a(1:n1), n1)
           call hdf5_write_grid_data_1d("xh", massfrac_h(1:n1), n1)
           call hdf5_write_grid_data_1d("xabar", massfrac_abar(1:n1), n1)
           call hdf5_write_grid_data_1d("xzbar", massfrac_zbar(1:n1), n1)
        endif
     endif

     if (.not.small_output) then
       if (eoskey.eq.1) then
          call hdf5_write_grid_data_1d("pressth", pressth(1:n1)/press_gf, n1)
       endif

       call hdf5_write_grid_data_1d("eps_kin", eps_kin(1:n1)/eps_gf, n1)

       allocate(cs(n1))
       cs(:) = sqrt(cs2(:))*clite
       call hdf5_write_grid_data_1d("cs", cs, n1)
       deallocate(cs)
     endif

     if(GR) then
        if(initial_data.eq."OSC") then
           call analytic_OSC_alpha(time*time_gf,10.d0,1.0d0, &
                alp_ana,rho_ana,vel_ana,X_ana,maxr)
           call hdf5_write_grid_data_1d("alpha_analytic", alp_ana(1:n1), n1)
           call hdf5_write_grid_data_1d("rho_analytic", rho_ana(1:n1)/rho_gf, n1)
           call hdf5_write_grid_data_1d("vel_analytic", vel_ana(1:n1)*clite, n1)
           call hdf5_write_grid_data_1d("alphamod", alp(1:n1), n1)
           call hdf5_write_grid_data_1d("X_analytic", X_ana(1:n1), n1)
        endif
        if (.not.small_output) call hdf5_write_grid_data_1d("alpha", alp(1:n1), n1)
        
        if (.not.small_output) then
          call hdf5_write_grid_data_1d("X", X(1:n1), n1)
          call hdf5_write_grid_data_1d("W", W(1:n1), n1)
        endif
        call hdf5_write_grid_data_1d("v", v(1:n1)*clite, n1)
      endif

     if (do_effectivepotential) then
        if (.not.small_output) call hdf5_write_grid_data_1d("alpha", alp(1:n1), n1)
     endif

     if (do_M1) then
        luminosity_rad = 0.0d0
        enden_rad = 0.0d0
        fluxden_rad = 0.0d0
        average_energy_rad = 0.0d0
        rms_energy_rad = 0.0d0

        do k=1,number_species
           do i=ghosts1+1,M1_imaxradii
              luminosity_rad(i,k,1) = sum(q_M1_fluid(i,k,:,2))*(4.0d0*pi)**2*x1(i)**2*(clite**5/ggrav)
              enden_rad(i,k,1) = sum(q_M1_fluid(i,k,:,1))*4.0d0*pi*length_gf**3/energy_gf
              fluxden_rad(i,k,1) = sum(q_M1_fluid(i,k,:,2))*4.0d0*pi*length_gf**3/energy_gf
              num_luminosity_rad(i,k,1) = sum(q_M1_fluid(i,k,:,2)/nulibtable_energies(:))* &
                   (4.0d0*pi)**2*x1(i)**2*(clite**3/ggrav*mass_gf)
              average_energy_rad(i,k,1) = sum(q_M1_fluid(i,k,:,1))/ &
                   sum(q_M1_fluid(i,k,:,1)/nulibtable_energies(:)) / (mev_to_erg*energy_gf)
              rms_energy_rad(i,k,1) = sqrt(sum(q_M1_fluid(i,k,:,1)*nulibtable_energies(:))/ &
                   sum(q_M1_fluid(i,k,:,1)/nulibtable_energies(:))) / (mev_to_erg*energy_gf)

              luminosity_rad(i,k,2) = sum(q_M1(i,k,:,2))*(4.0d0*pi)**2*x1(i)**2*(clite**5/ggrav)
              enden_rad(i,k,2) = sum(q_M1(i,k,:,1))*4.0d0*pi*length_gf**3/energy_gf
              fluxden_rad(i,k,2) = sum(q_M1(i,k,:,2))*4.0d0*pi*length_gf**3/energy_gf
              num_luminosity_rad(i,k,2) = sum(q_M1(i,k,:,2)/nulibtable_energies(:))* &
                   (4.0d0*pi)**2*x1(i)**2*(clite**3/ggrav*mass_gf)
              average_energy_rad(i,k,2) = sum(q_M1(i,k,:,1))/ &
                   sum(q_M1(i,k,:,1)/nulibtable_energies(:)) / (mev_to_erg*energy_gf)
              rms_energy_rad(i,k,2) = sqrt(sum(q_M1(i,k,:,1)*nulibtable_energies(:))/ &
                   sum(q_M1(i,k,:,1)/nulibtable_energies(:))) / (mev_to_erg*energy_gf)
           enddo
        enddo

        tau = 0.0d0
        ave_tau = 0.0d0
        do k=1,number_species
           do i=M1_imaxradii,ghosts1+1,-1
              ave_tau(i,k) = ave_tau(i+1,k) + (x1i(i+1)-x1i(i))/ &
                   (sum(q_M1(i,k,:,1)/(eas(i,k,:,2)+eas(i,k,:,3)))/sum(q_M1(i,k,:,1)))
           enddo
        enddo

        fluxfactor_enweighted = 0.0d0
        fluxfactor_fluxweighted = 0.0d0
        eddingtonfactor_enweighted = 0.0d0
        eddingtonfactor_fluxweighted = 0.0d0
        do k=1,number_species
           do i=ghosts1+1,M1_imaxradii
              do j=1,number_groups
                 fluxfactor_enweighted(i,k) = fluxfactor_enweighted(i,k) &
                      + q_M1(i,k,j,1)*q_M1(i,k,j,2)/q_M1(i,k,j,1)
                 fluxfactor_fluxweighted(i,k) = fluxfactor_fluxweighted(i,k) &
                      + q_M1(i,k,j,2)*q_M1(i,k,j,2)/q_M1(i,k,j,1)
                 eddingtonfactor_enweighted(i,k) = eddingtonfactor_enweighted(i,k) &
                      + q_M1(i,k,j,1)*q_M1(i,k,j,3)
                 eddingtonfactor_fluxweighted(i,k) = eddingtonfactor_fluxweighted(i,k) &
                      + q_M1(i,k,j,2)*q_M1(i,k,j,3)
              enddo
              fluxfactor_enweighted(i,k) = fluxfactor_enweighted(i,k) & 
                   / sum(q_M1(i,k,:,1))
              fluxfactor_fluxweighted(i,k) = fluxfactor_fluxweighted(i,k) &
                   / sum(q_M1(i,k,:,2))
              eddingtonfactor_enweighted(i,k) = eddingtonfactor_enweighted(i,k) & 
                   / sum(q_M1(i,k,:,1)) / X(k)**2
              eddingtonfactor_fluxweighted(i,k) = eddingtonfactor_fluxweighted(i,k) &
                   / sum(q_M1(i,k,:,2)) / X(k)**2
           enddo
        enddo

        if (.not.small_output) then
          call hdf5_write_grid_data_1d("M1_fluxfactor_enweighted_nue", fluxfactor_enweighted(:,1), n1)
          call hdf5_write_grid_data_1d("M1_fluxfactor_enweighted_anue", fluxfactor_enweighted(:,2), n1)
          call hdf5_write_grid_data_1d("M1_fluxfactor_enweighted_nux", fluxfactor_enweighted(:,3), n1)
          
          call hdf5_write_grid_data_1d("M1_fluxfactor_fluxweighted_nue", fluxfactor_fluxweighted(:,1), n1)
          call hdf5_write_grid_data_1d("M1_fluxfactor_fluxweighted_anue", fluxfactor_fluxweighted(:,2), n1)
          call hdf5_write_grid_data_1d("M1_fluxfactor_fluxweighted_nux", fluxfactor_fluxweighted(:,3), n1)

          call hdf5_write_grid_data_1d("M1_eddingtonfactor_enweighted_nue", eddingtonfactor_enweighted(:,1), n1)
          call hdf5_write_grid_data_1d("M1_eddingtonfactor_enweighted_anue", eddingtonfactor_enweighted(:,2), n1)
          call hdf5_write_grid_data_1d("M1_eddingtonfactor_enweighted_nux", eddingtonfactor_enweighted(:,3), n1)

          call hdf5_write_grid_data_1d("M1_eddingtonfactor_fluxweighted_nue", eddingtonfactor_fluxweighted(:,1), n1)
          call hdf5_write_grid_data_1d("M1_eddingtonfactor_fluxweighted_anue", eddingtonfactor_fluxweighted(:,2), n1)
          call hdf5_write_grid_data_1d("M1_eddingtonfactor_fluxweighted_nux", eddingtonfactor_fluxweighted(:,3), n1)

          call hdf5_write_grid_data_1d("dyedt_neutrino", dyedt_neutrino(1:n1)*time_gf, n1)
        endif

        call hdf5_write_grid_data_1d("M1_nue_luminosity_fluid_rad", luminosity_rad(:,1,1), n1)
        call hdf5_write_grid_data_1d("M1_nue_luminosity_lab_rad", luminosity_rad(:,1,2), n1)

        if (.not.small_output) then
          call hdf5_write_grid_data_1d("M1_nue_enden_lab_rad", enden_rad(:,1,2), n1)
          call hdf5_write_grid_data_1d("M1_nue_fluxden_lab_rad", fluxden_rad(:,1,2), n1)
          call hdf5_write_grid_data_1d("M1_nue_numluminosity_fluid_rad", num_luminosity_rad(:,1,1), n1)
        endif

        call hdf5_write_grid_data_1d("M1_nue_aveenergy_fluid_rad", average_energy_rad(:,1,1), n1)
        call hdf5_write_grid_data_1d("M1_nue_rmsenergy_fluid_rad", rms_energy_rad(:,1,1), n1)

        call hdf5_write_grid_data_1d("M1_anue_luminosity_fluid_rad", luminosity_rad(:,2,1), n1)
        call hdf5_write_grid_data_1d("M1_anue_luminosity_lab_rad", luminosity_rad(:,2,2), n1)

        if (.not.small_output) then
          call hdf5_write_grid_data_1d("M1_anue_enden_lab_rad", enden_rad(:,2,2), n1)
          call hdf5_write_grid_data_1d("M1_anue_fluxden_lab_rad", fluxden_rad(:,2,2), n1)
          call hdf5_write_grid_data_1d("M1_anue_numluminosity_fluid_rad", num_luminosity_rad(:,2,1), n1)
        endif

        call hdf5_write_grid_data_1d("M1_anue_aveenergy_fluid_rad", average_energy_rad(:,2,1), n1)
        call hdf5_write_grid_data_1d("M1_anue_rmsenergy_fluid_rad", rms_energy_rad(:,2,1), n1)

        call hdf5_write_grid_data_1d("M1_nux_luminosity_fluid_rad", luminosity_rad(:,3,1), n1)
        call hdf5_write_grid_data_1d("M1_nux_luminosity_lab_rad", luminosity_rad(:,3,2), n1)

        if (.not.small_output) then
          call hdf5_write_grid_data_1d("M1_nux_enden_lab_rad", enden_rad(:,3,2), n1)
          call hdf5_write_grid_data_1d("M1_nux_fluxden_lab_rad", fluxden_rad(:,3,2), n1)
          call hdf5_write_grid_data_1d("M1_nux_numluminosity_fluid_rad", num_luminosity_rad(:,3,1), n1)
        endif

        call hdf5_write_grid_data_1d("M1_nux_aveenergy_fluid_rad", average_energy_rad(:,3,1), n1)
        call hdf5_write_grid_data_1d("M1_nux_rmsenergy_fluid_rad", rms_energy_rad(:,3,1), n1)

        if (.not.small_output) then
          call hdf5_write_grid_data_1d("M1_nue_ng1_rad", q_M1(:,1,1,1), n1)
        endif

     endif

     ! Output M1 spectra as HDF5 datasets during modeflag==1
     if (do_M1) then
        ! Extraction radius spectra
        spectrum = q_M1(M1_iextractradii,1,:,2)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_nue_fluxspectra_out", spectrum, number_groups)
        spectrum = q_M1(M1_iextractradii,2,:,2)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_anue_fluxspectra_out", spectrum, number_groups)
        spectrum = q_M1(M1_iextractradii,3,:,2)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_nux_fluxspectra_out", spectrum, number_groups)

        spectrum = q_M1_fluid(M1_iextractradii,1,:,1)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_nue_enspectra_out", spectrum, number_groups)
        spectrum = q_M1_fluid(M1_iextractradii,2,:,1)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_anue_enspectra_out", spectrum, number_groups)
        spectrum = q_M1_fluid(M1_iextractradii,3,:,1)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_nux_enspectra_out", spectrum, number_groups)

        ! Central cell spectra
        spectrum = q_M1(ghosts1+1,1,:,2)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_nue_fluxspectra_cen", spectrum, number_groups)
        spectrum = q_M1(ghosts1+1,2,:,2)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_anue_fluxspectra_cen", spectrum, number_groups)
        spectrum = q_M1(ghosts1+1,3,:,2)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_nux_fluxspectra_cen", spectrum, number_groups)

        spectrum = q_M1_fluid(ghosts1+1,1,:,1)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_nue_enspectra_cen", spectrum, number_groups)
        spectrum = q_M1_fluid(ghosts1+1,2,:,1)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_anue_enspectra_cen", spectrum, number_groups)
        spectrum = q_M1_fluid(ghosts1+1,3,:,1)*M1_moment_to_distro(:)
        call hdf5_write_grid_data_1d("M1_nux_enspectra_cen", spectrum, number_groups)

        ! Capturing factors
        spectrum = alp(ghosts1+1)*(1.0d0-eas(ghosts1+1,1,:,2)* &
             q_M1_fluid(ghosts1+1,1,:,1)/eas(ghosts1+1,1,:,1))
        call hdf5_write_grid_data_1d("capturing_factors", spectrum, number_groups)
     endif

     if(eoskey.eq.3) then
        if (.not.small_output) then
          call hdf5_write_grid_data_1d("entropy", entropy(1:n1), n1)
        endif
        
        call hdf5_write_grid_data_1d("temperature", temp(1:n1), n1)
     endif
     
     ! Close xg.h5 file ONCE at the end of modeflag==1
     call hdf5_close_xg_file()
     
  else if(modeflag.eq.2) then
     
     ! Open dat.h5 file ONCE at the beginning
     call hdf5_open_dat_file()
      
     call hdf5_append_scalar("time", time)
     if (initial_data.eq.'Collapse') then
        !Shock radius
        if (bounce) then
           call hdf5_append_scalar("shock_radius_t", shock_radius/length_gf)
           call hdf5_append_scalar("binding_energy_total", binding_energy_total/energy_gf)
        endif
        
        ! Mass values
        call hdf5_append_scalar("mgrav_Xmax", mgravX)
        call hdf5_append_scalar("mbary_shock", mass(ishock(1)))
        call hdf5_append_scalar("mgrav_shock", mgrav(ishock(1)))
        call hdf5_append_scalar("mgrav_rho1e12", mgrav12)
        call hdf5_append_scalar("mbary_Xmax", mbaryX)
        call hdf5_append_scalar("mbary_rho1e12", mbary12)
        call hdf5_append_scalar("r_Xmax", rXmax/length_gf)
        call hdf5_append_scalar("r_rho1e12", r12max/length_gf)
        call hdf5_append_scalar("r_rho1e11", r11max/length_gf)
        call hdf5_append_scalar("M_innercore", mass_inner_core)
        
        !rotation scalars
        if (do_rotation) then
           call hdf5_append_scalar("total_angular_momentum", angular_momentum/mass_gf*time_gf/length_gf**2)
           call hdf5_append_scalar("ToverW_edge", ToverW(n1-ghosts1-1))
        endif
     endif

     if (do_M1) then

        luminosity = 0.0d0 
        num_luminosity = 0.0d0 
        average_energy = 0.0d0
        rms_energy = 0.0d0
        luminosity_fluid = 0.0d0 
        num_luminosity_fluid = 0.0d0 
        average_energy_fluid = 0.0d0
        rms_energy_fluid = 0.0d0

         do k=1,number_species
           luminosity(k)     = sum(q_M1(M1_iextractradii,k,:,2))&
                * (4.0d0*pi)**2*x1(M1_iextractradii)**2 * (clite**5/ggrav)
           luminosity_fluid(k) = sum(q_M1_fluid(M1_iextractradii,k,:,2))&
                * (4.0d0*pi)**2*x1(M1_iextractradii)**2 * (clite**5/ggrav)
           num_luminosity_fluid(k) = sum(q_M1_fluid(M1_iextractradii,k,:,2)/nulibtable_energies(:))&
                * (4.0d0*pi)**2*x1(M1_iextractradii)**2 * (clite**3/ggrav*mass_gf)
           average_energy_fluid(k) = sum(q_M1_fluid(M1_iextractradii,k,:,2))&
                / sum(q_M1_fluid(M1_iextractradii,k,:,2)/nulibtable_energies(:)) / (mev_to_erg*energy_gf)
           rms_energy_fluid(k)     = sqrt( sum(q_M1_fluid(M1_iextractradii,k,:,2)*nulibtable_energies(:))&
                / sum(q_M1_fluid(M1_iextractradii,k,:,2)/nulibtable_energies(:)) ) / (mev_to_erg*energy_gf)
           average_energy(k) = average_energy_fluid(k)/sqrt(1.0d0-v(M1_iextractradii)**2)*(1.0d0+v(M1_iextractradii))
           rms_energy(k)     = rms_energy_fluid(k)/sqrt(1.0d0-v(M1_iextractradii)**2)*(1.0d0+v(M1_iextractradii))
           num_luminosity(k) = luminosity(k)/(average_energy(k)*mev_to_erg)
        enddo

        total_nu_energy = 0.0d0
        do i=ghosts1+1,M1_imaxradii
           total_nu_energy = total_nu_energy + sum(q_M1(i,1,:,1)*4.0d0*pi*volume(i)/energy_gf) + &
                sum(q_M1(i,2,:,1)*4.0d0*pi*volume(i)/energy_gf) + &
                sum(q_M1(i,3,:,1)*4.0d0*pi*volume(i)/energy_gf)
        enddo

        tau = 0.0d0
        ave_tau = 0.0d0
        nusphere_not_found = .true.
        ave_nusphere_not_found = .true.
        do k=1,number_species
           do i=M1_imaxradii,ghosts1+1,-1
              ave_tau(i,k) = ave_tau(i+1,k) + (x1i(i+1)-x1i(i))/ &
                   (sum(q_M1(i,k,:,1)/(eas(i,k,:,2)+eas(i,k,:,3)))/sum(q_M1(i,k,:,1)))
              if (ave_nusphere_not_found(k).and.ave_tau(i,k).gt.0.66666666d0) then
                 ave_nusphere_not_found(k) = .false.
                 ave_nusphere(k,1) = x1i(i)/length_gf
                 ave_nusphere(k,2) = rho(i)/rho_gf
                 ave_nusphere(k,3) = temp(i)
                 ave_nusphere(k,4) = ye(i)
                 ave_nusphere(k,5) = entropy(i)
              endif
           enddo

           do j=1,number_groups
              tau = 0.0d0
              do i=M1_imaxradii,ghosts1+1,-1
                 tau(i,1) = tau(i+1,1) + (x1i(i+1)-x1i(i))*(eas(i,k,j,2)+eas(i,k,j,3))
                 tau(i,2) = tau(i+1,2) + (x1i(i+1)-x1i(i))*eas(i,k,j,2)
                 if (nusphere_not_found(k,j,1).and.tau(i,1).gt.0.66666666d0) then
                    nusphere_not_found(k,j,1) = .false.
                    nusphere(k,j,1) = x1i(i)
                 endif
                 if (nusphere_not_found(k,j,2).and.tau(i,2).gt.0.66666666d0) then
                    nusphere_not_found(k,j,2) = .false.
                    nusphere(k,j,2) = x1i(i)
                 endif
              enddo
           enddo
        enddo

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = total_nu_energy
        scalars(2) = total_energy_radiated
        scalars(3) = total_energy_absorped
        call hdf5_append_scalar_array("M1_energies", scalars(1:nscalars), nscalars)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = luminosity(1)
        scalars(2) = luminosity(2)
        scalars(3) = luminosity(3)
        call hdf5_append_scalar_array("M1_flux_lum", scalars(1:nscalars), nscalars)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = luminosity_fluid(1)
        scalars(2) = luminosity_fluid(2)
        scalars(3) = luminosity_fluid(3)
        call hdf5_append_scalar_array("M1_flux_lum_fluid", scalars(1:nscalars), nscalars)
        
        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = num_luminosity(1)
        scalars(2) = num_luminosity(2)
        scalars(3) = num_luminosity(3)
        call hdf5_append_scalar_array("M1_flux_numlum", scalars(1:nscalars), nscalars)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = num_luminosity_fluid(1)
        scalars(2) = num_luminosity_fluid(2)
        scalars(3) = num_luminosity_fluid(3)
        call hdf5_append_scalar_array("M1_flux_numlum_fluid", scalars(1:nscalars), nscalars)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 5
        scalars(1) = total_net_heating/(luminosity(1)+luminosity(2))            
        scalars(2) = luminosity(1)+luminosity(2)
        scalars(3) = total_net_heating
        scalars(4) = total_mass_gain
        scalars(5) = total_net_deintdt
        call hdf5_append_scalar_array("M1_net_heating", scalars(1:nscalars), nscalars)
        
        scalars(1:nscalars0) = 0.0d0
        nscalars = number_groups*number_species*2
        if (nscalars.gt.256) stop "increase scalar count"
        do k=1,number_species
           do j=1,number_groups
              scalars((k-1)*number_groups*2+(j-1)*2+1) = nusphere(k,j,1)/length_gf
              scalars((k-1)*number_groups*2+(j-1)*2+2) = nusphere(k,j,2)/length_gf
           enddo
        enddo
        call hdf5_append_scalar_array("M1_nusphere_allE", scalars(1:nscalars), nscalars)

        scalars(1:nscalars0) = 0.0d0
        nscalars = number_species*5
        if (nscalars.gt.256) stop "increase scalar count"
        do k=1,number_species
           scalars((k-1)*5+1) = ave_nusphere(k,1)
           scalars((k-1)*5+2) = ave_nusphere(k,2)
           scalars((k-1)*5+3) = ave_nusphere(k,3)
           scalars((k-1)*5+4) = ave_nusphere(k,4)
           scalars((k-1)*5+5) = ave_nusphere(k,5)
        enddo
        call hdf5_append_scalar_array("M1_nusphere_ave", scalars(1:nscalars), nscalars)
        
        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = average_energy(1)
        scalars(2) = average_energy(2)
        scalars(3) = average_energy(3)
        call hdf5_append_scalar_array("M1_flux_aveenergy_lab", scalars(1:nscalars), nscalars)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = rms_energy_fluid(1)
        scalars(2) = rms_energy_fluid(2)
        scalars(3) = rms_energy_fluid(3)
        call hdf5_append_scalar_array("M1_flux_rmsenergy_fluid", scalars(1:nscalars), nscalars)

        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = average_energy_fluid(1)
        scalars(2) = average_energy_fluid(2)
        scalars(3) = average_energy_fluid(3)
        call hdf5_append_scalar_array("M1_flux_aveenergy_fluid", scalars(1:nscalars), nscalars)
        
        scalars(1:nscalars0) = 0.0d0
        nscalars = 3
        scalars(1) = rms_energy(1)
        scalars(2) = rms_energy(2)
        scalars(3) = rms_energy(3)
        call hdf5_append_scalar_array("M1_flux_rmsenergy_lab", scalars(1:nscalars), nscalars)

        if (.not.small_output) then
          call hdf5_append_scalar("dyedt_neutrino_c_t", dyedt_neutrino(ghosts1+1)*time_gf)
        endif

     endif

     ! central values
     call hdf5_append_scalar("rho_c_t", rho(ghosts1+1)/rho_gf)
     call hdf5_append_scalar("ye_c_t", ye(ghosts1+1))

     if (.not.small_output) then
        call hdf5_append_scalar("dyedt_hydro_c_t", dyedt_hydro(ghosts1+1)*time_gf)
     endif

     call hdf5_append_scalar("ynu_c_t", ynu(ghosts1+1))
     call hdf5_append_scalar("csound_c_t", sqrt(cs2(ghosts1+1)))
     
     if(eoskey.eq.3) then
        call hdf5_append_scalar("entropy_c_t", entropy(ghosts1+1))
        call hdf5_append_scalar("temperature_c_t", temp(ghosts1+1))
     endif
     
     call hdf5_append_scalar("totalmass", totalmass/mass_gf)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           call analytic_OSC_alpha(time*time_gf,10.d0,1.0d0, &
                alp_ana,rho_ana,vel_ana,X_ana,maxr)
           call hdf5_append_scalar("alpha_analytic_c_t", alp_ana(ghosts1+1))
           call hdf5_append_scalar("rho_analytic_c_t", rho_ana(ghosts1+1)/rho_gf)
           call hdf5_append_scalar("vel_analytic_c_t", vel_ana(ghosts1+1)*clite)
        endif
        call hdf5_append_scalar("alpha_c_t", alp(ghosts1+1))
        call hdf5_append_scalar("time_c", time_c)
     endif
     
     if(initial_data.eq."Sedov") then
        ishock_radius = maxloc(abs(v1))
        call hdf5_append_scalar("Sedov_radius", x1(ishock_radius(1)))
        call hdf5_append_scalar("Sedov_velocity", v1(ishock_radius(1)))
        call hdf5_append_scalar("Sedov_press", press(ishock_radius(1)))
        call hdf5_append_scalar("Sedov_density", rho(ishock_radius(1)))
     endif
     
     if (initial_data.eq.'Collapse') then
        call hdf5_append_scalar_array("accretion_rates", accretion_rates(1:11), 11)
        call hdf5_append_scalar_array("accreted_mass", accreted_mass(1:11), 11)
     endif
     
     ! Close dat.h5 file ONCE at the end of modeflag==2
     call hdf5_close_dat_file()
     
  endif
  
  ! Close HDF5 library after all writes complete
  call hdf5_finalize()

end subroutine output_all_HDF5
