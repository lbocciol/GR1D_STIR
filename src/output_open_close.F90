!-*-f90-*-
subroutine Open_output_files(modeflag)

  use GR1D_module
  implicit none

  character*1024 filename
  integer modeflag

  if(modeflag.eq.0) then
     
  else if(modeflag.eq.1) then

     filename = trim(adjustl(outdir))//"/v1.xg"
     call open_output_file(filename,101)
     
     if(do_rotation) then
        filename = trim(adjustl(outdir))//"/omega.xg"
        call open_output_file(filename,102)
        filename = trim(adjustl(outdir))//"/ToverW.xg"
        call open_output_file(filename,103)
        
        if(GR) then
           filename = trim(adjustl(outdir))//"/vphi.xg"
           call open_output_file(filename,104)
        else
           filename = trim(adjustl(outdir))//"/vphi1.xg"
           call open_output_file(filename,105)
        endif
     endif

     if (do_turbulence) then
        filename = trim(adjustl(outdir))//"/omega2_BV.xg"
        call open_output_file(filename,106)
        filename = trim(adjustl(outdir))//"/v_turb.xg"
        call open_output_file(filename,107)
        filename = trim(adjustl(outdir))//"/dissipated_turb_eps.xg"
        if (.not. small_output) call open_output_file(filename,108)
        filename = trim(adjustl(outdir))//"/buoyancy_turb_eps.xg"
        if (.not. small_output) call open_output_file(filename,109)
        filename = trim(adjustl(outdir))//"/shear_turb_eps.xg"
        if (.not. small_output) call open_output_file(filename,110)
        filename = trim(adjustl(outdir))//"/Lambda_MLT.xg"
        if (.not. small_output) call open_output_file(filename,111)
     endif   

     filename = trim(adjustl(outdir))//"/rho.xg"
     call open_output_file(filename,112)
     
     if (do_nupress.or.do_M1) then
        filename = trim(adjustl(outdir))//"/nuchem.xg"
        if (.not.small_output) call open_output_file(filename,113)
        filename = trim(adjustl(outdir))//"/press_nu.xg"
        if (.not.small_output) call open_output_file(filename,114)
        filename = trim(adjustl(outdir))//"/energy_nu.xg"
        if (.not.small_output) call open_output_file(filename,115) !use press_gf b/c neutrino energy is not specific
        filename = trim(adjustl(outdir))//"/dnupdr.xg"
        if (.not.small_output) call open_output_file(filename,116)
     endif

     filename = trim(adjustl(outdir))//"/ye.xg"
     call open_output_file(filename,117)
     if (do_M1) then
        filename = trim(adjustl(outdir))//"/dyedt_hydro.xg"
        if (.not.small_output) call open_output_file(filename,118)
        filename = trim(adjustl(outdir))//"/depsdt.xg"
        if (.not.small_output) call open_output_file(filename,119)
        filename = trim(adjustl(outdir))//"/ynu.xg"
        if (.not.small_output) call open_output_file(filename,120)
     endif
     
     
     filename = trim(adjustl(outdir))//"/press.xg"
     call open_output_file(filename,121)
     
     filename = trim(adjustl(outdir))//"/eps.xg"
     if (.not.small_output) call open_output_file(filename,122)
     
     if (GR) then
        filename = trim(adjustl(outdir))//"/mass_grav.xg"
        call open_output_file(filename,123)
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call open_output_file(filename,124)
     else
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call open_output_file(filename,125)
     endif
     
     if (eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/xn.xg"
        if (.not.small_output) call open_output_file(filename,126)
        
        filename = trim(adjustl(outdir))//"/xp.xg"
        if (.not.small_output) call open_output_file(filename,127)
        
        filename = trim(adjustl(outdir))//"/xa.xg"
        if (.not.small_output) call open_output_file(filename,128)
        
        filename = trim(adjustl(outdir))//"/xh.xg"
        if (.not.small_output) call open_output_file(filename,129)
        
        filename = trim(adjustl(outdir))//"/xabar.xg"
        if (.not.small_output) call open_output_file(filename,130)
        
        filename = trim(adjustl(outdir))//"/xzbar.xg"
        if (.not.small_output) call open_output_file(filename,131)
     endif

     if (eoskey.eq.1) then
        filename = trim(adjustl(outdir))//"/pressth.xg"
        if (.not.small_output) call open_output_file(filename,132)
     endif

     filename = trim(adjustl(outdir))//"/eps_kin.xg"
     if (.not.small_output) call open_output_file(filename,133)

     filename = trim(adjustl(outdir))//"/cs.xg"
     if (.not.small_output) call open_output_file(filename,134)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           filename = trim(adjustl(outdir))//"/alpha_analytic.xg"
           call open_output_file(filename,135)
           filename = trim(adjustl(outdir))//"/rho_analytic.xg"
           call open_output_file(filename,136)
           filename = trim(adjustl(outdir))//"/vel_analytic.xg"
           call open_output_file(filename,137)
           filename = trim(adjustl(outdir))//"/alphamod.xg"
           call open_output_file(filename,138)
           filename = trim(adjustl(outdir))//"/X_analytic.xg"
           call open_output_file(filename,139)
        endif
        filename = trim(adjustl(outdir))//"/alpha.xg"
        call open_output_file(filename,140)
        filename = trim(adjustl(outdir))//"/X.xg"
        if (.not.small_output) call open_output_file(filename,141)
        filename = trim(adjustl(outdir))//"/W.xg"
        if (.not.small_output) call open_output_file(filename,142)
        filename = trim(adjustl(outdir))//"/v.xg"
        call open_output_file(filename,143)
     endif

     if (do_effectivepotential) then
        filename = trim(adjustl(outdir))//"/alpha.xg"
        if (.not.small_output) call open_output_file(filename,144)
     endif

     if (do_M1) then

        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_nue.xg"
        if (.not.small_output) call open_output_file(filename,145)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_anue.xg"
        if (.not.small_output) call open_output_file(filename,146)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_nux.xg"
        if (.not.small_output) call open_output_file(filename,147)
        
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_nue.xg"
        if (.not.small_output) call open_output_file(filename,148)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_anue.xg"
        if (.not.small_output) call open_output_file(filename,149)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_nux.xg"
        if (.not.small_output) call open_output_file(filename,150)

        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_nue.xg"
        if (.not.small_output) call open_output_file(filename,151)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_anue.xg"
        if (.not.small_output) call open_output_file(filename,152)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_nux.xg"
        if (.not.small_output) call open_output_file(filename,153)

        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_nue.xg"
        if (.not.small_output) call open_output_file(filename,154)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_anue.xg"
        if (.not.small_output) call open_output_file(filename,155)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_nux.xg"
        if (.not.small_output) call open_output_file(filename,156)

        filename = trim(adjustl(outdir))//"/dyedt_neutrino.xg"
        if (.not.small_output) call open_output_file(filename,157)

        filename = trim(adjustl(outdir))//"/M1_nue_luminosity_fluid_rad.xg"
        call open_output_file(filename,158)

        filename = trim(adjustl(outdir))//"/M1_nue_luminosity_lab_rad.xg"
        call open_output_file(filename,159)

        filename = trim(adjustl(outdir))//"/M1_nue_enden_lab_rad.xg"
        if (.not.small_output) call open_output_file(filename,160)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxden_lab_rad.xg"
        if (.not.small_output) call open_output_file(filename,161)

        filename = trim(adjustl(outdir))//"/M1_nue_numluminosity_fluid_rad.xg"
        if (.not.small_output) call open_output_file(filename,162)
        
        filename = trim(adjustl(outdir))//"/M1_nue_aveenergy_fluid_rad.xg"
        call open_output_file(filename,163)

        filename = trim(adjustl(outdir))//"/M1_nue_rmsenergy_fluid_rad.xg"
        call open_output_file(filename,164)


        filename = trim(adjustl(outdir))//"/M1_anue_luminosity_fluid_rad.xg"
        call open_output_file(filename,165)

        filename = trim(adjustl(outdir))//"/M1_anue_luminosity_lab_rad.xg"
        call open_output_file(filename,166)

        filename = trim(adjustl(outdir))//"/M1_anue_enden_lab_rad.xg"
        if (.not.small_output) call open_output_file(filename,167)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxden_lab_rad.xg"
        if (.not.small_output) call open_output_file(filename,168)

        filename = trim(adjustl(outdir))//"/M1_anue_numluminosity_fluid_rad.xg"
        if (.not.small_output) call open_output_file(filename,169)
        
        filename = trim(adjustl(outdir))//"/M1_anue_aveenergy_fluid_rad.xg"
        call open_output_file(filename,170)

        filename = trim(adjustl(outdir))//"/M1_anue_rmsenergy_fluid_rad.xg"
        call open_output_file(filename,171)


        filename = trim(adjustl(outdir))//"/M1_nux_luminosity_fluid_rad.xg"
        call open_output_file(filename,172)

        filename = trim(adjustl(outdir))//"/M1_nux_luminosity_lab_rad.xg"
        call open_output_file(filename,173)

        filename = trim(adjustl(outdir))//"/M1_nux_enden_lab_rad.xg"
        if (.not.small_output) call open_output_file(filename,174)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxden_lab_rad.xg"
        if (.not.small_output) call open_output_file(filename,175)

        filename = trim(adjustl(outdir))//"/M1_nux_numluminosity_fluid_rad.xg"
        if (.not.small_output) call open_output_file(filename,176)

        filename = trim(adjustl(outdir))//"/M1_nux_aveenergy_fluid_rad.xg"
        call open_output_file(filename,177)

        filename = trim(adjustl(outdir))//"/M1_nux_rmsenergy_fluid_rad.xg"
        call open_output_file(filename,178)

            
        filename = trim(adjustl(outdir))//"/M1_nue_ng1_rad.xg"
        if (.not.small_output) call open_output_file(filename,179)

     endif

     if(eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/entropy.xg"
        if (.not.small_output) call open_output_file(filename,180)
        filename = trim(adjustl(outdir))//"/temperature.xg"
        call open_output_file(filename,181)
     endif
     
  else if(modeflag.eq.2) then
     
     if (initial_data.eq.'Collapse') then
        !Shock radius
        filename = trim(adjustl(outdir))//"/shock_radius_t.dat"
        call open_output_file(filename,182)
        filename = trim(adjustl(outdir))//"/binding_energy_total.dat"
        call open_output_file(filename,183)
           
        ! Mass values
        filename = trim(adjustl(outdir))//"/mgrav_Xmax.dat"
        call open_output_file(filename,185)
        
        filename = trim(adjustl(outdir))//"/mbary_shock.dat"
        call open_output_file(filename,186)
        
        filename = trim(adjustl(outdir))//"/mgrav_shock.dat"
        call open_output_file(filename,187)
        
        filename = trim(adjustl(outdir))//"/mgrav_rho1e12.dat"
        call open_output_file(filename,188)

        filename = trim(adjustl(outdir))//"/mbary_Xmax.dat"
        call open_output_file(filename,189)
        
        filename = trim(adjustl(outdir))//"/mbary_rho1e12.dat"
        call open_output_file(filename,190)
        
        filename = trim(adjustl(outdir))//"/r_Xmax.dat"
        call open_output_file(filename,191)
        
        filename = trim(adjustl(outdir))//"/r_rho1e12.dat"
        call open_output_file(filename,192)

        filename = trim(adjustl(outdir))//"/r_rho1e11.dat"
        call open_output_file(filename,193)
        
        filename = trim(adjustl(outdir))//"/M_innercore.dat"
        call open_output_file(filename,194)
        
        !rotation scalars
        if (do_rotation) then
           filename = trim(adjustl(outdir))//"/total_angular_momentum.dat"
           call open_output_file(filename,195)
           filename = trim(adjustl(outdir))//"/ToverW_edge.dat"
           call open_output_file(filename,196)
        endif
     endif


     if (do_M1) then

        filename = trim(adjustl(outdir))//"/M1_energies.dat"
        call open_output_file(filename,197)

        filename = trim(adjustl(outdir))//"/M1_flux_lum.dat"
        call open_output_file(filename,198)

        filename = trim(adjustl(outdir))//"/M1_flux_lum_fluid.dat"
        call open_output_file(filename,199)

        filename = trim(adjustl(outdir))//"/M1_flux_numlum.dat"
        call open_output_file(filename,200)

        filename = trim(adjustl(outdir))//"/M1_flux_numlum_fluid.dat"
        call open_output_file(filename,201)

        filename = trim(adjustl(outdir))//"/M1_net_heating.dat"
        call open_output_file(filename,202)

        filename = trim(adjustl(outdir))//"/M1_nusphere_allE.dat"
        call open_output_file(filename,203)

        filename = trim(adjustl(outdir))//"/M1_nusphere_ave.dat"
        call open_output_file(filename,204)
 
        filename = trim(adjustl(outdir))//"/M1_flux_aveenergy_lab.dat"
        call open_output_file(filename,205)

        filename = trim(adjustl(outdir))//"/M1_flux_rmsenergy_fluid.dat"
        call open_output_file(filename,206)

        filename = trim(adjustl(outdir))//"/M1_flux_aveenergy_fluid.dat"
        call open_output_file(filename,207)

        filename = trim(adjustl(outdir))//"/M1_flux_rmsenergy_lab.dat"
        call open_output_file(filename,208)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxspectra_out.xg"
        call open_output_file(filename,209)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxspectra_out.xg"
        call open_output_file(filename,210)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxspectra_out.xg"
        call open_output_file(filename,211)

        filename = trim(adjustl(outdir))//"/M1_nue_enspectra_out.xg"
        call open_output_file(filename,212)

        filename = trim(adjustl(outdir))//"/M1_anue_enspectra_out.xg"
        call open_output_file(filename,213)

        filename = trim(adjustl(outdir))//"/M1_nux_enspectra_out.xg"
        call open_output_file(filename,214)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxspectra_cen.xg"
        call open_output_file(filename,215)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxspectra_cen.xg"
        call open_output_file(filename,216)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxspectra_cen.xg"
        call open_output_file(filename,217)

        filename = trim(adjustl(outdir))//"/M1_nue_enspectra_cen.xg"
        call open_output_file(filename,218)

        filename = trim(adjustl(outdir))//"/M1_anue_enspectra_cen.xg"
        call open_output_file(filename,219)
        
        filename = trim(adjustl(outdir))//"/M1_nux_enspectra_cen.xg"
        call open_output_file(filename,220)

        filename = trim(adjustl(outdir))//"/capturing_factors.xg"
        call open_output_file(filename,221)

        filename = trim(adjustl(outdir))//"/dyedt_neutrino_c_t.dat"
        if (.not.small_output)  call open_output_file(filename,222)

     endif

     ! central values
     filename = trim(adjustl(outdir))//"/rho_c_t.dat"
     call open_output_file(filename,223)
        
     filename = trim(adjustl(outdir))//"/ye_c_t.dat"
     call open_output_file(filename,224)

     filename = trim(adjustl(outdir))//"/dyedt_hydro_c_t.dat"
     if (.not.small_output) call open_output_file(filename,225)

     filename = trim(adjustl(outdir))//"/ynu_c_t.dat"
     call open_output_file(filename,226)

     filename = trim(adjustl(outdir))//"/csound_c_t.dat"
     call open_output_file(filename,227)
     
     if(eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/entropy_c_t.dat"
        call open_output_file(filename,228)
        filename = trim(adjustl(outdir))//"/temperature_c_t.dat"
        call open_output_file(filename,229)
     endif
     
     filename = trim(adjustl(outdir))//"/totalmass.dat"
     call open_output_file(filename,230)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           filename = trim(adjustl(outdir))//"/alpha_analytic_c_t.dat"
           call open_output_file(filename,231)
           filename = trim(adjustl(outdir))//"/rho_analytic_c_t.dat"
           call open_output_file(filename,232)
           filename = trim(adjustl(outdir))//"/vel_analytic_c_t.dat"
           call open_output_file(filename,233)
        endif
        filename = trim(adjustl(outdir))//"/alpha_c_t.dat"
        call open_output_file(filename,234)
        filename = trim(adjustl(outdir))//"/time_c.dat"
        call open_output_file(filename,235)
     endif
     
     if(initial_data.eq."Sedov") then
        filename = trim(adjustl(outdir))//"/Sedov_radius.dat"
        call open_output_file(filename,236)
        filename = trim(adjustl(outdir))//"/Sedov_velocity.dat"
        call open_output_file(filename,237)
        filename = trim(adjustl(outdir))//"/Sedov_press.dat"
        call open_output_file(filename,238)
        filename = trim(adjustl(outdir))//"/Sedov_density.dat"
        call open_output_file(filename,239)
        
     endif
     
     if (initial_data.eq.'Collapse') then
        filename = trim(adjustl(outdir))//"/accretion_rates.dat"
        call open_output_file(filename,240)
        filename = trim(adjustl(outdir))//"/accreted_mass.dat"
        call open_output_file(filename,241)
     endif
     
  endif
 
end subroutine Open_output_files
 
subroutine Close_output_files(modeflag)

  use GR1D_module
  implicit none

  character*1024 filename
  integer modeflag

  if(modeflag.eq.0) then
     
  else if(modeflag.eq.1) then

     filename = trim(adjustl(outdir))//"/v1.xg"
     call close_output_file(filename,101)
     
     if(do_rotation) then
        filename = trim(adjustl(outdir))//"/omega.xg"
        call close_output_file(filename,102)
        filename = trim(adjustl(outdir))//"/ToverW.xg"
        call close_output_file(filename,103)
        
        if(GR) then
           filename = trim(adjustl(outdir))//"/vphi.xg"
           call close_output_file(filename,104)
        else
           filename = trim(adjustl(outdir))//"/vphi1.xg"
           call close_output_file(filename,105)
        endif
     endif

     if (do_turbulence) then
        filename = trim(adjustl(outdir))//"/omega2_BV.xg"
        call close_output_file(filename,106)
        filename = trim(adjustl(outdir))//"/v_turb.xg"
        call close_output_file(filename,107)
        filename = trim(adjustl(outdir))//"/dissipated_turb_eps.xg"
        if (.not. small_output) call close_output_file(filename,108)
        filename = trim(adjustl(outdir))//"/buoyancy_turb_eps.xg"
        if (.not. small_output) call close_output_file(filename,109)
        filename = trim(adjustl(outdir))//"/shear_turb_eps.xg"
        if (.not. small_output) call close_output_file(filename,110)
        filename = trim(adjustl(outdir))//"/Lambda_MLT.xg"
        if (.not. small_output) call close_output_file(filename,111)
     endif   

     filename = trim(adjustl(outdir))//"/rho.xg"
     call close_output_file(filename,112)
     
     if (do_nupress.or.do_M1) then
        filename = trim(adjustl(outdir))//"/nuchem.xg"
        if (.not.small_output) call close_output_file(filename,113)
        filename = trim(adjustl(outdir))//"/press_nu.xg"
        if (.not.small_output) call close_output_file(filename,114)
        filename = trim(adjustl(outdir))//"/energy_nu.xg"
        if (.not.small_output) call close_output_file(filename,115) !use press_gf b/c neutrino energy is not specific
        filename = trim(adjustl(outdir))//"/dnupdr.xg"
        if (.not.small_output) call close_output_file(filename,116)
     endif

     filename = trim(adjustl(outdir))//"/ye.xg"
     call close_output_file(filename,117)
     if (do_M1) then
        filename = trim(adjustl(outdir))//"/dyedt_hydro.xg"
        if (.not.small_output) call close_output_file(filename,118)
        filename = trim(adjustl(outdir))//"/depsdt.xg"
        if (.not.small_output) call close_output_file(filename,119)
        filename = trim(adjustl(outdir))//"/ynu.xg"
        if (.not.small_output) call close_output_file(filename,120)
     endif
     
     
     filename = trim(adjustl(outdir))//"/press.xg"
     call close_output_file(filename,121)
     
     filename = trim(adjustl(outdir))//"/eps.xg"
     if (.not.small_output) call close_output_file(filename,122)
     
     if (GR) then
        filename = trim(adjustl(outdir))//"/mass_grav.xg"
        call close_output_file(filename,123)
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call close_output_file(filename,124)
     else
        filename = trim(adjustl(outdir))//"/mass_bary.xg"
        call close_output_file(filename,125)
     endif
     
     if (eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/xn.xg"
        if (.not.small_output) call close_output_file(filename,126)
        
        filename = trim(adjustl(outdir))//"/xp.xg"
        if (.not.small_output) call close_output_file(filename,127)
        
        filename = trim(adjustl(outdir))//"/xa.xg"
        if (.not.small_output) call close_output_file(filename,128)
        
        filename = trim(adjustl(outdir))//"/xh.xg"
        if (.not.small_output) call close_output_file(filename,129)
        
        filename = trim(adjustl(outdir))//"/xabar.xg"
        if (.not.small_output) call close_output_file(filename,130)
        
        filename = trim(adjustl(outdir))//"/xzbar.xg"
        if (.not.small_output) call close_output_file(filename,131)
     endif

     if (eoskey.eq.1) then
        filename = trim(adjustl(outdir))//"/pressth.xg"
        if (.not.small_output) call close_output_file(filename,132)
     endif

     filename = trim(adjustl(outdir))//"/eps_kin.xg"
     if (.not.small_output) call close_output_file(filename,133)

     filename = trim(adjustl(outdir))//"/cs.xg"
     if (.not.small_output) call close_output_file(filename,134)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           filename = trim(adjustl(outdir))//"/alpha_analytic.xg"
           call close_output_file(filename,135)
           filename = trim(adjustl(outdir))//"/rho_analytic.xg"
           call close_output_file(filename,136)
           filename = trim(adjustl(outdir))//"/vel_analytic.xg"
           call close_output_file(filename,137)
           filename = trim(adjustl(outdir))//"/alphamod.xg"
           call close_output_file(filename,138)
           filename = trim(adjustl(outdir))//"/X_analytic.xg"
           call close_output_file(filename,139)
        endif
        filename = trim(adjustl(outdir))//"/alpha.xg"
        call close_output_file(filename,140)
        filename = trim(adjustl(outdir))//"/X.xg"
        if (.not.small_output) call close_output_file(filename,141)
        filename = trim(adjustl(outdir))//"/W.xg"
        if (.not.small_output) call close_output_file(filename,142)
        filename = trim(adjustl(outdir))//"/v.xg"
        call close_output_file(filename,143)
     endif

     if (do_effectivepotential) then
        filename = trim(adjustl(outdir))//"/alpha.xg"
        if (.not.small_output) call close_output_file(filename,144)
     endif

     if (do_M1) then

        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_nue.xg"
        if (.not.small_output) call close_output_file(filename,145)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_anue.xg"
        if (.not.small_output) call close_output_file(filename,146)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_enweighted_nux.xg"
        if (.not.small_output) call close_output_file(filename,147)
        
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_nue.xg"
        if (.not.small_output) call close_output_file(filename,148)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_anue.xg"
        if (.not.small_output) call close_output_file(filename,149)
        filename = trim(adjustl(outdir))//"/M1_fluxfactor_fluxweighted_nux.xg"
        if (.not.small_output) call close_output_file(filename,150)

        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_nue.xg"
        if (.not.small_output) call close_output_file(filename,151)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_anue.xg"
        if (.not.small_output) call close_output_file(filename,152)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_enweighted_nux.xg"
        if (.not.small_output) call close_output_file(filename,153)

        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_nue.xg"
        if (.not.small_output) call close_output_file(filename,154)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_anue.xg"
        if (.not.small_output) call close_output_file(filename,155)
        filename = trim(adjustl(outdir))//"/M1_eddingtonfactor_fluxweighted_nux.xg"
        if (.not.small_output) call close_output_file(filename,156)

        filename = trim(adjustl(outdir))//"/dyedt_neutrino.xg"
        if (.not.small_output) call close_output_file(filename,157)

        filename = trim(adjustl(outdir))//"/M1_nue_luminosity_fluid_rad.xg"
        call close_output_file(filename,158)

        filename = trim(adjustl(outdir))//"/M1_nue_luminosity_lab_rad.xg"
        call close_output_file(filename,159)

        filename = trim(adjustl(outdir))//"/M1_nue_enden_lab_rad.xg"
        if (.not.small_output) call close_output_file(filename,160)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxden_lab_rad.xg"
        if (.not.small_output) call close_output_file(filename,161)

        filename = trim(adjustl(outdir))//"/M1_nue_numluminosity_fluid_rad.xg"
        if (.not.small_output) call close_output_file(filename,162)
        
        filename = trim(adjustl(outdir))//"/M1_nue_aveenergy_fluid_rad.xg"
        call close_output_file(filename,163)

        filename = trim(adjustl(outdir))//"/M1_nue_rmsenergy_fluid_rad.xg"
        call close_output_file(filename,164)


        filename = trim(adjustl(outdir))//"/M1_anue_luminosity_fluid_rad.xg"
        call close_output_file(filename,165)

        filename = trim(adjustl(outdir))//"/M1_anue_luminosity_lab_rad.xg"
        call close_output_file(filename,166)

        filename = trim(adjustl(outdir))//"/M1_anue_enden_lab_rad.xg"
        if (.not.small_output) call close_output_file(filename,167)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxden_lab_rad.xg"
        if (.not.small_output) call close_output_file(filename,168)

        filename = trim(adjustl(outdir))//"/M1_anue_numluminosity_fluid_rad.xg"
        if (.not.small_output) call close_output_file(filename,169)
        
        filename = trim(adjustl(outdir))//"/M1_anue_aveenergy_fluid_rad.xg"
        call close_output_file(filename,170)

        filename = trim(adjustl(outdir))//"/M1_anue_rmsenergy_fluid_rad.xg"
        call close_output_file(filename,171)


        filename = trim(adjustl(outdir))//"/M1_nux_luminosity_fluid_rad.xg"
        call close_output_file(filename,172)

        filename = trim(adjustl(outdir))//"/M1_nux_luminosity_lab_rad.xg"
        call close_output_file(filename,173)

        filename = trim(adjustl(outdir))//"/M1_nux_enden_lab_rad.xg"
        if (.not.small_output) call close_output_file(filename,174)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxden_lab_rad.xg"
        if (.not.small_output) call close_output_file(filename,175)

        filename = trim(adjustl(outdir))//"/M1_nux_numluminosity_fluid_rad.xg"
        if (.not.small_output) call close_output_file(filename,176)

        filename = trim(adjustl(outdir))//"/M1_nux_aveenergy_fluid_rad.xg"
        call close_output_file(filename,177)

        filename = trim(adjustl(outdir))//"/M1_nux_rmsenergy_fluid_rad.xg"
        call close_output_file(filename,178)

            
        filename = trim(adjustl(outdir))//"/M1_nue_ng1_rad.xg"
        if (.not.small_output) call close_output_file(filename,179)

     endif

     if(eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/entropy.xg"
        if (.not.small_output) call close_output_file(filename,180)
        filename = trim(adjustl(outdir))//"/temperature.xg"
        call close_output_file(filename,181)
     endif
     
  else if(modeflag.eq.2) then
     
     if (initial_data.eq.'Collapse') then
        !Shock radius
        filename = trim(adjustl(outdir))//"/shock_radius_t.dat"
        call close_output_file(filename,182)
        filename = trim(adjustl(outdir))//"/binding_energy_total.dat"
        call close_output_file(filename,183)
           
        ! Mass values
        filename = trim(adjustl(outdir))//"/mgrav_Xmax.dat"
        call close_output_file(filename,185)
        
        filename = trim(adjustl(outdir))//"/mbary_shock.dat"
        call close_output_file(filename,186)
        
        filename = trim(adjustl(outdir))//"/mgrav_shock.dat"
        call close_output_file(filename,187)
        
        filename = trim(adjustl(outdir))//"/mgrav_rho1e12.dat"
        call close_output_file(filename,188)

        filename = trim(adjustl(outdir))//"/mbary_Xmax.dat"
        call close_output_file(filename,189)
        
        filename = trim(adjustl(outdir))//"/mbary_rho1e12.dat"
        call close_output_file(filename,190)
        
        filename = trim(adjustl(outdir))//"/r_Xmax.dat"
        call close_output_file(filename,191)
        
        filename = trim(adjustl(outdir))//"/r_rho1e12.dat"
        call close_output_file(filename,192)

        filename = trim(adjustl(outdir))//"/r_rho1e11.dat"
        call close_output_file(filename,193)
        
        filename = trim(adjustl(outdir))//"/M_innercore.dat"
        call close_output_file(filename,194)
        
        !rotation scalars
        if (do_rotation) then
           filename = trim(adjustl(outdir))//"/total_angular_momentum.dat"
           call close_output_file(filename,195)
           filename = trim(adjustl(outdir))//"/ToverW_edge.dat"
           call close_output_file(filename,196)
        endif
     endif


     if (do_M1) then

        filename = trim(adjustl(outdir))//"/M1_energies.dat"
        call close_output_file(filename,197)

        filename = trim(adjustl(outdir))//"/M1_flux_lum.dat"
        call close_output_file(filename,198)

        filename = trim(adjustl(outdir))//"/M1_flux_lum_fluid.dat"
        call close_output_file(filename,199)

        filename = trim(adjustl(outdir))//"/M1_flux_numlum.dat"
        call close_output_file(filename,200)

        filename = trim(adjustl(outdir))//"/M1_flux_numlum_fluid.dat"
        call close_output_file(filename,201)

        filename = trim(adjustl(outdir))//"/M1_net_heating.dat"
        call close_output_file(filename,202)

        filename = trim(adjustl(outdir))//"/M1_nusphere_allE.dat"
        call close_output_file(filename,203)

        filename = trim(adjustl(outdir))//"/M1_nusphere_ave.dat"
        call close_output_file(filename,204)
 
        filename = trim(adjustl(outdir))//"/M1_flux_aveenergy_lab.dat"
        call close_output_file(filename,205)

        filename = trim(adjustl(outdir))//"/M1_flux_rmsenergy_fluid.dat"
        call close_output_file(filename,206)

        filename = trim(adjustl(outdir))//"/M1_flux_aveenergy_fluid.dat"
        call close_output_file(filename,207)

        filename = trim(adjustl(outdir))//"/M1_flux_rmsenergy_lab.dat"
        call close_output_file(filename,208)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxspectra_out.xg"
        call close_output_file(filename,209)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxspectra_out.xg"
        call close_output_file(filename,210)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxspectra_out.xg"
        call close_output_file(filename,211)

        filename = trim(adjustl(outdir))//"/M1_nue_enspectra_out.xg"
        call close_output_file(filename,212)

        filename = trim(adjustl(outdir))//"/M1_anue_enspectra_out.xg"
        call close_output_file(filename,213)

        filename = trim(adjustl(outdir))//"/M1_nux_enspectra_out.xg"
        call close_output_file(filename,214)

        filename = trim(adjustl(outdir))//"/M1_nue_fluxspectra_cen.xg"
        call close_output_file(filename,215)

        filename = trim(adjustl(outdir))//"/M1_anue_fluxspectra_cen.xg"
        call close_output_file(filename,216)

        filename = trim(adjustl(outdir))//"/M1_nux_fluxspectra_cen.xg"
        call close_output_file(filename,217)

        filename = trim(adjustl(outdir))//"/M1_nue_enspectra_cen.xg"
        call close_output_file(filename,218)

        filename = trim(adjustl(outdir))//"/M1_anue_enspectra_cen.xg"
        call close_output_file(filename,219)
        
        filename = trim(adjustl(outdir))//"/M1_nux_enspectra_cen.xg"
        call close_output_file(filename,220)

        filename = trim(adjustl(outdir))//"/capturing_factors.xg"
        call close_output_file(filename,221)

        filename = trim(adjustl(outdir))//"/dyedt_neutrino_c_t.dat"
        if (.not.small_output)  call close_output_file(filename,222)

     endif

     ! central values
     filename = trim(adjustl(outdir))//"/rho_c_t.dat"
     call close_output_file(filename,223)
        
     filename = trim(adjustl(outdir))//"/ye_c_t.dat"
     call close_output_file(filename,224)

     filename = trim(adjustl(outdir))//"/dyedt_hydro_c_t.dat"
     if (.not.small_output) call close_output_file(filename,225)

     filename = trim(adjustl(outdir))//"/ynu_c_t.dat"
     call close_output_file(filename,226)

     filename = trim(adjustl(outdir))//"/csound_c_t.dat"
     call close_output_file(filename,227)
     
     if(eoskey.eq.3) then
        filename = trim(adjustl(outdir))//"/entropy_c_t.dat"
        call close_output_file(filename,228)
        filename = trim(adjustl(outdir))//"/temperature_c_t.dat"
        call close_output_file(filename,229)
     endif
     
     filename = trim(adjustl(outdir))//"/totalmass.dat"
     call close_output_file(filename,230)
     
     if(GR) then
        if(initial_data.eq."OSC") then
           filename = trim(adjustl(outdir))//"/alpha_analytic_c_t.dat"
           call close_output_file(filename,231)
           filename = trim(adjustl(outdir))//"/rho_analytic_c_t.dat"
           call close_output_file(filename,232)
           filename = trim(adjustl(outdir))//"/vel_analytic_c_t.dat"
           call close_output_file(filename,233)
        endif
        filename = trim(adjustl(outdir))//"/alpha_c_t.dat"
        call close_output_file(filename,234)
        filename = trim(adjustl(outdir))//"/time_c.dat"
        call close_output_file(filename,235)
     endif
     
     if(initial_data.eq."Sedov") then
        filename = trim(adjustl(outdir))//"/Sedov_radius.dat"
        call close_output_file(filename,236)
        filename = trim(adjustl(outdir))//"/Sedov_velocity.dat"
        call close_output_file(filename,237)
        filename = trim(adjustl(outdir))//"/Sedov_press.dat"
        call close_output_file(filename,238)
        filename = trim(adjustl(outdir))//"/Sedov_density.dat"
        call close_output_file(filename,239)
        
     endif
     
     if (initial_data.eq.'Collapse') then
        filename = trim(adjustl(outdir))//"/accretion_rates.dat"
        call close_output_file(filename,240)
        filename = trim(adjustl(outdir))//"/accreted_mass.dat"
        call close_output_file(filename,241)
     endif
     
  endif
 
end subroutine Close_output_files
 
 subroutine open_output_file(filename,nfile)
    
  implicit none
  character(*) filename
  integer nfile
  
  open(unit=nfile,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  
end subroutine open_output_file
  
 subroutine close_output_file(filename,nfile)
    
  implicit none
  character(*) filename
  integer nfile
  
  open(unit=nfile,file=trim(adjustl(filename)),status="unknown",&
       form='formatted',position="append")
  
end subroutine close_output_file
