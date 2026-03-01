!-*-f90-*-
subroutine turb_diff_terms

   use GR1D_module
   use Grad_module
   implicit none

   integer:: i
   real*8 :: Lambda_diff, Lambda_diff_m, Lambda_diff_p
   real*8 :: D_turb_eps, D_turb_ye, D_turb_K, D_turb_mu_p, D_turb_mu_m
   real*8 :: eps_grad_p(n1), invrho_grad_p(n1)
   real*8 :: ye_grad_p(n1), v2_turb_grad_p(n1)
   real*8 :: v_grad_p(n1), v_grad_m(n1)

   eps_grad_p = Gradient_int_p(eps, x1)
   invrho_grad_p = - 1.0d0/rhop**2*Gradient_int_p(rho, x1)
   ye_grad_p = Gradient_int_p(ye,x1)
   v2_turb_grad_p = Gradient_int_p(v_turb**2,x1)
   v_grad_p = Gradient_int_p(v1, x1)
   v_grad_m = Gradient_int_m(v1, x1)

   !Implement the calculation of the diffusion terms
   do i=ghosts1+1,n1-ghosts1
      call Get_Diffusion_Length(x1(i), pressp(i), rhop(i), dphidr(i), Lambda_diff)

      ! Calculate diffusion coefficients as given in Couch et al. 2019 eqs. 29-31
      D_turb_eps = alpha_turb_e  * v_turbp(i) * Lambda_diff
      D_turb_ye  = alpha_turb_ye * v_turbp(i) * Lambda_diff
      D_turb_K   = alpha_turb_K  * v_turbp(i) * Lambda_diff
   
      diff_term_ye(i)  = qp(i,1) * D_turb_ye * ye_grad_p(i)
      diff_term_K(i)   = qp(i,1) * D_turb_K * v2_turb_grad_p(i)
      diff_term_eps(i) = qp(i,1) * D_turb_eps * (eps_grad_p(i) + &
            press(i)*invrho_grad_p(i))
   enddo 

end subroutine turb_diff_terms

subroutine turbulence_sources

   use GR1D_module
   use Grad_module
   implicit none

   integer i
   real*8 Lambda_mix, Lambda_diss, P_turb, D_turb_mu
   real*8 v_grad(n1), visc_term(n1), visc_grad(n1)

   visc_term = 0.0d0
   visc_grad = 0.0d0
   v_grad = Gradient_3pts(v1,x1)

   !Implement the calculation of the diffusion terms
   do i=ghosts1+1,n1-ghosts1
      call Get_Mixing_Length(x1(i), press(i), rho(i), dphidr(i), Lambda_mix)
      call Get_Dissipation_Length(x1(i), press(i), rho(i), dphidr(i), v_turb(i), omega2_BV(i), Lambda_diss)

      shear(i) = - v_turb(i)**2 * v_grad(i)
      diss(i)  = v_turb(i)**3 / Lambda_diss
      buoy(i)  = v_turb(i) * omega2_BV(i) * Lambda_mix
      if (do_turbulent_viscosity) then
        visc(i)  = 4.0d0/3.0d0 * alpha_turb_mu * Lambda_mix * v_turb(i) * (v_grad(i) - v(i)/x1(i))**2
      endif

      P_turb = 2.0d0*c_turb*rho(i)*v_turb(i)**2
      turb_source(i,2) = 2.0d0*alp(i)*(1 - c_turb)/c_turb*P_turb/X(i)/x1(i)
      turb_source(i,3) = rho(i)* (shear(i) + buoy(i) + visc(i))
      turb_source(i,6) = rho(i)* (shear(i) + buoy(i) - diss(i) )

      if (GR) then
          turb_source(i,3) = alp(i)*X(i)*turb_source(i,3)
          turb_source(i,6) = alp(i)*X(i)*turb_source(i,6)
      else
          turb_source(i,3) = sqrt_gamma(i)*turb_source(i,3)
          turb_source(i,6) = sqrt_gamma(i)*turb_source(i,6)
      endif
      
      if (do_turbulent_viscosity) then
        D_turb_mu  = alpha_turb_mu * v_turb(i) * Lambda_mix
        visc_term(i) = x1(i)**3 * D_turb_mu * q(i,1) * (v_grad(i) - v(i)/x1(i))
      endif

   enddo

   if (do_turbulent_viscosity) then
      visc_grad = Gradient_3pts(visc_term,x1)
      do i=ghosts1+1,n1-ghosts1
        turb_source(i,2) = turb_source(i,2) + alp(i)*X(i)*4.0d0/3.0d0/x1(i)**3*visc_grad(i)
      enddo
   endif

end subroutine turbulence_sources

subroutine Brunt_Vaisala(dts)

   use GR1D_module
   use Grad_module
   implicit none
   
   !local
   integer :: i 
   real*8  :: dts
   real*8  :: v_turb_seed, h, CL
   real*8  :: rho_grad(n1), press_grad(n1), v_grad(n1)
   real*8  :: entropy_grad(n1), yL_grad(n1), lns_grad(n1), lnye_grad(n1)
   real*8  :: this_dpds_rho_yL, this_dpdyL_s_rho, this_dpdrho_s_yL
   real*8  :: this_dlnpdlns_rho_ye, this_dlnpdlnye_s_rho, this_dlnpdlnrho_s_ye
   
   if (GR) then
       rho_grad = Gradient_3pts(rho(:)*(1.0d0 + eps(:)),x1)
       press_grad = Gradient_3pts(press,x1)
       v_grad = Gradient_3pts(v,x1)
   else
       rho_grad = Gradient_3pts(rho,x1)
       press_grad = Gradient_3pts(press,x1)   
       v_grad = Gradient_3pts(v1,x1)
   endif
   
   lns_grad  = Gradient_3pts(log(entropy), x1)
   lnye_grad = Gradient_3pts(log(ye)     , x1)

   do i=ghosts1+1,n1-ghosts1

      call Interpolate_additional_lnderivs( rho(i)/rho_gf, temp(i), ye(i), &
         this_dlnpdlns_rho_ye, this_dlnpdlnye_s_rho, this_dlnpdlnrho_s_ye )
    
      ! This is a more general form
      CL = - rho(i) / this_dlnpdlnrho_s_ye * &
          (this_dlnpdlns_rho_ye * lns_grad(i)  + this_dlnpdlnye_s_rho * lnye_grad(i) )

      ! This is standard STIR
      CL = (rho_grad(i) - press_grad(i)/cs2(i))
      if (GR) then
          h = 1.0d0 + eps(i) + press(i)/rho(i) + v_turb(i)**2

          omega2_BV(i) = alp(i)**2/(rho(i)*h*X(i)**2)* &
                (dphidr(i) - v1(i)* v_grad(i)) * CL
                
      else
          omega2_BV(i) = (dphidr(i) - v(i)* v_grad(i))/rho(i) * CL
      endif

      call Get_Mixing_Length(x1(i), press(i), rho(i), dphidr(i), Lambda_MLT(i))
 
      if (omega2_BV(i) .gt. 0.0d0) then
          v_turb_seed = dts*omega2_BV(i)*Lambda_MLT(i)
          v_turb(i) = max(v_turb(i),v_turb_seed)
      endif

      if ( x1(i) / length_gf .lt. 1.0d5 ) v_turb(i) = 0.0d0

      ! force covective velocity to be zero in the 5 points stencil around the shock,
      ! i.e. do not convect through the shock
      if ( ishock(1)-1 <= i-ghosts1 .and. i-ghosts1 <= ishock(1)+1 ) v_turb(i) = 0.0d0

   enddo

end subroutine Brunt_Vaisala

subroutine Get_Dissipation_Length(r, p, rho, dphidr, vturb, BV2, Lambda)

  use GR1D_module, only: alpha_turb, length_gf, &
    pns_radius, shock_radius, gain_radius

  implicit none

  real(8), intent(in)  :: r, p, rho, dphidr, vturb, BV2
  real(8), parameter   :: delta = 1.0d-10
  real(8), intent(out) :: Lambda
  real(8) :: max_Lambda

  Lambda = alpha_turb * p / (rho * dphidr)

  if (BV2 < 0.0d0) then
    max_Lambda = sqrt(-vturb**2/BV2)
    Lambda = min(Lambda, max_Lambda)
    Lambda = max(Lambda, delta*length_gf)
  endif

end subroutine Get_Dissipation_Length

subroutine Get_Diffusion_Length(r, p, rho, dphidr, Lambda)

  use GR1D_module, only: alpha_turb, length_gf, &
    pns_radius, shock_radius

  implicit none

  real(8), intent(in)  :: r, p, rho, dphidr
  real(8), intent(out) :: Lambda

  CALL Get_Mixing_Length(r, p, rho, dphidr, Lambda)
  RETURN

  ! Lambda = alpha_turb * p / (rho * dphidr)

  ! if (r < pns_radius) then
  !   Lambda = min(Lambda, pns_radius - r, r)
  ! else if (r < shock_radius) then
  !   Lambda = min(Lambda, shock_radius - r, r)
  ! endif

end subroutine Get_Diffusion_Length

subroutine Get_Mixing_Length(r, p, rho, dphidr, Lambda)

  use GR1D_module, only: alpha_turb, length_gf, &
    pns_radius, shock_radius, gain_radius

  implicit none

  real(8), intent(in)  :: r, p, rho, dphidr
  real(8), intent(out) :: Lambda

  Lambda = alpha_turb * p / (rho * dphidr)
  ! Lambda = min(Lambda, r)

end subroutine Get_Mixing_Length
