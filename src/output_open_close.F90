!-*-f90-*-
subroutine Open_output_files(modeflag)

  use GR1D_module
  implicit none

  character*1024 filename
  integer modeflag, iFile

  if(modeflag.eq.0) then
     
  else if(modeflag.eq.1) then

     iFile = 100

     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/v1.xg"
     call open_output_file(filename,iFile)
     
     if(do_rotation) then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/omega.xg"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/ToverW.xg"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        if(GR) then
           filename = trim(adjustl(outdir))//"/vphi.xg"
           call open_output_file(filename,iFile)
        else
           filename = trim(adjustl(outdir))//"/vphi1.xg"
           call open_output_file(filename,iFile)
        endif
     endif

     if (do_turbulence) then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/omega2_BV.xg"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/v_turb.xg"
        call open_output_file(filename,iFile)
        if (.not. small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/dissipated_turb_eps.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/buoyancy_turb_eps.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/shear_turb_eps.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/Lambda_MLT.xg"
          call open_output_file(filename,iFile)
        endif
     endif   

     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/rho.xg"
     call open_output_file(filename,iFile)
     
     if ( (do_nupress.or.do_M1) .and. (.not.small_output) ) then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/nuchem.xg"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/press_nu.xg"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/energy_nu.xg"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/dnupdr.xg"
        call open_output_file(filename,iFile)
     endif

     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/ye.xg"
     call open_output_file(filename,iFile)
     if (do_M1 .and. (.not.small_output)) then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/dyedt_hydro.xg"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/depsdt.xg"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/ynu.xg"
        call open_output_file(filename,iFile)
     endif
 
     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/press.xg"
     call open_output_file(filename,iFile)
     
     if (.not.small_output) then
       iFile = iFile + 1
       filename = trim(adjustl(outdir))//"/eps.xg"
       call open_output_file(filename,iFile)
     endif

     if (GR) then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mass_grav.xg"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call open_output_file(filename,iFile)
     else
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call open_output_file(filename,iFile)
     endif
     
     if (eoskey.eq.3 .and. (.not.small_output)) then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/xn.xg"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/xp.xg"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1 
        filename = trim(adjustl(outdir))//"/xa.xg"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/xh.xg"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/xabar.xg"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/xzbar.xg"
        call open_output_file(filename,iFile)
     endif

     if (.not.small_output) then
       if (eoskey.eq.1) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/pressth.xg"
          call open_output_file(filename,iFile)
       endif

       iFile = iFile + 1
       filename = trim(adjustl(outdir))//"/eps_kin.xg"
       call open_output_file(filename,iFile)

       iFile = iFile + 1
       filename = trim(adjustl(outdir))//"/cs.xg"
       call open_output_file(filename,iFile)
     endif

     if(GR) then
        if(initial_data.eq."OSC") then
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/alpha_analytic.xg"
           call open_output_file(filename,iFile)
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/rho_analytic.xg"
           call open_output_file(filename,iFile)
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/vel_analytic.xg"
           call open_output_file(filename,iFile)
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/alphamod.xg"
           call open_output_file(filename,iFile)
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/X_analytic.xg"
           call open_output_file(filename,iFile)
        endif
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/alpha.xg"
        call open_output_file(filename,iFile)
        
        if (.not.small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/X.xg"
          call open_output_file(filename,iFile)
          
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/W.xg"
          call open_output_file(filename,iFile)
        
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/v.xg"
          call open_output_file(filename,iFile)
        endif
      endif

     if (do_effectivepotential .and. (.not.small_output)) then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/alpha.xg"
        call open_output_file(filename,iFile)
     endif

     if (do_M1) then

        if (.not.small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_nue.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_anue.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_nux.xg"
          call open_output_file(filename,iFile)
          
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_nue.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_anue.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_nux.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_nue.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_anue.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_nux.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_nue.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_anue.xg"
          call open_output_file(filename,iFile)
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_nux.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/dyedt_neutrino.xg"
          call open_output_file(filename,iFile)
        endif

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nue_luminosity_fluid_rad.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nue_luminosity_lab_rad.xg"
        call open_output_file(filename,iFile)

        if (.not.small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_nue_enden_lab_rad.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_nue_fluxden_lab_rad.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_nue_numluminosity_fluid_rad.xg"
          call open_output_file(filename,iFile)
        endif

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nue_aveenergy_fluid_rad.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nue_rmsenergy_fluid_rad.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_anue_luminosity_fluid_rad.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_anue_luminosity_lab_rad.xg"
        call open_output_file(filename,iFile)

        if (.not.small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_anue_enden_lab_rad.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_anue_fluxden_lab_rad.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_anue_numluminosity_fluid_rad.xg"
          call open_output_file(filename,iFile)
        endif

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_anue_aveenergy_fluid_rad.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_anue_rmsenergy_fluid_rad.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nux_luminosity_fluid_rad.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nux_luminosity_lab_rad.xg"
        call open_output_file(filename,iFile)

        if (.not.small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_nux_enden_lab_rad.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_nux_fluxden_lab_rad.xg"
          call open_output_file(filename,iFile)

          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_nux_numluminosity_fluid_rad.xg"
          call open_output_file(filename,iFile)
        endif

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nux_aveenergy_fluid_rad.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nux_rmsenergy_fluid_rad.xg"
        call open_output_file(filename,iFile)

        if (.not.small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/M1_nue_ng1_rad.xg"
          call open_output_file(filename,iFile)
        endif

     endif

     if(eoskey.eq.3) then
        if (.not.small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/entropy.xg"
          call open_output_file(filename,iFile)
        endif

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/temperature.xg"
        call open_output_file(filename,iFile)
     endif
     
     nOutFiles1 = iFile - 100

  else if(modeflag.eq.2) then
     
     iFile = 1000

     if (initial_data.eq.'Collapse') then
        !Shock radius
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/shock_radius_t.dat"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/binding_energy_total.dat"
        call open_output_file(filename,iFile)
           
        ! Mass values
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mgrav_Xmax.dat"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mbary_shock.dat"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mgrav_shock.dat"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mgrav_rho1e12.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mbary_Xmax.dat"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/mbary_rho1e12.dat"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/r_Xmax.dat"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/r_rho1e12.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/r_rho1e11.dat"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M_innercore.dat"
        call open_output_file(filename,iFile)
        
        !rotation scalars
        if (do_rotation) then
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/total_angular_momentum.dat"
           call open_output_file(filename,iFile)
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/ToverW_edge.dat"
           call open_output_file(filename,iFile)
        endif
     endif


     if (do_M1) then

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_energies.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_flux_lum.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_flux_lum_fluid.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_flux_numlum.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_flux_numlum_fluid.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_net_heating.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nusphere_allE.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nusphere_ave.dat"
        call open_output_file(filename,iFile)
 
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_flux_aveenergy_lab.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_flux_rmsenergy_fluid.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_flux_aveenergy_fluid.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_flux_rmsenergy_lab.dat"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nue_fluxspectra_out.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_anue_fluxspectra_out.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nux_fluxspectra_out.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nue_enspectra_out.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_anue_enspectra_out.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nux_enspectra_out.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nue_fluxspectra_cen.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_anue_fluxspectra_cen.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nux_fluxspectra_cen.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nue_enspectra_cen.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_anue_enspectra_cen.xg"
        call open_output_file(filename,iFile)
        
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/M1_nux_enspectra_cen.xg"
        call open_output_file(filename,iFile)

        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/capturing_factors.xg"
        call open_output_file(filename,iFile)

        if (.not.small_output) then
          iFile = iFile + 1
          filename = trim(adjustl(outdir))//"/dyedt_neutrino_c_t.dat"
          call open_output_file(filename,iFile)
        endif

     endif

     ! central values
     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/rho_c_t.dat"
     call open_output_file(filename,iFile)
        
     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/ye_c_t.dat"
     call open_output_file(filename,iFile)

     if (.not.small_output) then
       iFile = iFile + 1
       filename = trim(adjustl(outdir))//"/dyedt_hydro_c_t.dat"
       call open_output_file(filename,iFile)
     endif

     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/ynu_c_t.dat"
     call open_output_file(filename,iFile)

     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/csound_c_t.dat"
     call open_output_file(filename,iFile)
     
     if(eoskey.eq.3) then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/entropy_c_t.dat"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/temperature_c_t.dat"
        call open_output_file(filename,iFile)
     endif
     
     iFile = iFile + 1
     filename = trim(adjustl(outdir))//"/totalmass.dat"
     call open_output_file(filename,iFile)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/alpha_analytic_c_t.dat"
           call open_output_file(filename,iFile)
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/rho_analytic_c_t.dat"
           call open_output_file(filename,iFile)
           iFile = iFile + 1
           filename = trim(adjustl(outdir))//"/vel_analytic_c_t.dat"
           call open_output_file(filename,iFile)
        endif
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/alpha_c_t.dat"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/time_c.dat"
        call open_output_file(filename,iFile)
     endif
     
     if(initial_data.eq."Sedov") then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/Sedov_radius.dat"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/Sedov_velocity.dat"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/Sedov_press.dat"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/Sedov_density.dat"
        call open_output_file(filename,iFile)
        
     endif
     
     if (initial_data.eq.'Collapse') then
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/accretion_rates.dat"
        call open_output_file(filename,iFile)
        iFile = iFile + 1
        filename = trim(adjustl(outdir))//"/accreted_mass.dat"
        call open_output_file(filename,iFile)
     endif
    
     nOutFiles2 = iFile - 1000

  endif
 
end subroutine Open_output_files
 
subroutine Close_output_files(modeflag)

  use GR1D_module
  implicit none

  integer, intent(in) :: modeflag
  integer :: iFile

  if(modeflag.eq.0) then
     
  else if(modeflag.eq.1) then

    do iFile = 1, nOutFiles1
      close(unit=iFile + 100)
    enddo        

  else if(modeflag.eq.2) then
     
    do iFile = 1, nOutFiles2
      close(unit=iFile + 1000)
    enddo        
     
  endif
 
end subroutine Close_output_files
 
 subroutine open_output_file(filename,nfile)
    
  implicit none
  character(*) filename
  integer nfile
  
  open(unit=nfile,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  
end subroutine open_output_file
