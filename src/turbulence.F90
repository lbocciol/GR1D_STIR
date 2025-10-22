subroutine turb_diff_terms

   use GR1D_module
   use Grad_module
   implicit none

   integer:: i
   real*8 :: Lambda_mix
   real*8 :: D_turb_eps, D_turb_ye, D_turb_K
   real*8 :: eps_grad_p(n1), invr_grad_p(n1), ye_grad_p(n1), v2_turb_grad_p(n1)
   real*8 :: eps_grad_m(n1), invr_grad_m(n1), ye_grad_m(n1), v2_turb_grad_m(n1)

   eps_grad_p = Gradient_int_p(eps, x1)
   invr_grad_p = Gradient_int_p(1.0d0/rho, x1)
   ye_grad_p = Gradient_int_p(ye,x1)
   v2_turb_grad_p = Gradient_int_p(v_turb**2,x1)

   eps_grad_m = Gradient_int_m(eps, x1)
   invr_grad_m = Gradient_int_m(1.0d0/rho, x1)
   ye_grad_m = Gradient_int_m(ye,x1)
   v2_turb_grad_m = Gradient_int_m(v_turb**2,x1)

   !Implement the calculation of the diffusion terms
   do i=ghosts1+1,n1-ghosts1
      call Get_Mixing_Length(x1i(i), pressp(i), rhop(i), dphidr(i), Lambda_mix)

      ! Calculate diffusion coefficients as given in Couch et al. 2019 eqs. 29-31
      D_turb_eps = alpha_turb_e * v_turbp(i) * Lambda_mix
      D_turb_ye = alpha_turb_ye * v_turbp(i) * Lambda_mix
      D_turb_K = alpha_turb_K * v_turbp(i) * Lambda_mix

      diff_term_K(i) = qp(i,1) * D_turb_K * v2_turb_grad_p(i)
      diff_term_eps(i) = qp(i,1) * D_turb_eps * (eps_grad_p(i) + &
                          press(i)*invr_grad_p(i))
      diff_term_ye(i) = qp(i,1) * D_turb_ye * ye_grad_p(i)
      ! call Get_Mixing_Length(x1i(i-1), pressm(i), rhom(i), dphidr(i), Lambda_mix)

      ! ! Calculate diffusion coefficients as given in Couch et al. 2019 eqs. 29-31
      ! D_turb_eps = alpha_turb_e * v_turbm(i) * Lambda_mix
      ! D_turb_ye = alpha_turb_ye * v_turbm(i) * Lambda_mix
      ! D_turb_K = alpha_turb_K * v_turbm(i) * Lambda_mix

      ! diff_term_K_m(i) = qm(i,1) * D_turb_K * v2_turb_grad_m(i)
      ! diff_term_eps_m(i) = qm(i,1) * D_turb_eps * (eps_grad_m(i) + &
      !                     press(i)*invr_grad_m(i))
      ! diff_term_ye_m(i) = qm(i,1) * D_turb_ye * ye_grad_m(i)
   enddo

end subroutine turb_diff_terms

subroutine turbulence_sources

   use GR1D_module
   use Grad_module
   implicit none

   integer i
   real*8 Lambda_mix, Lambda_diss
   real*8 v_grad(n1)

   v_grad = Gradient_5pts(v1,x1)

   !Implement the calculation of the diffusion terms
   do i=ghosts1+1,n1-ghosts1
      call Get_Mixing_Length(x1(i), press(i), rho(i), dphidr(i), Lambda_mix)
      call Get_Dissipation_Length(x1(i), press(i), rho(i), dphidr(i), v_turb(i), omega2_BV(i), Lambda_diss)

      shear(i) = - v_turb(i)**2 * v_grad(i)
      shear(i) = 0.0d0
      diss(i) = v_turb(i)**3 / Lambda_diss
      buoy(i) = v_turb(i) * omega2_BV(i) * Lambda_mix

      !turb_source(i,3) = rho(i)* (shear(i) + buoy(i))
      turb_source(i,3) = rho(i)* diss(i)
      turb_source(i,6) = rho(i)* (shear(i) + buoy(i) - diss(i))

      if (GR) then
          turb_source(i,3) = alp(i)*X(i)*turb_source(i,3)
          turb_source(i,6) = alp(i)*X(i)*turb_source(i,6)
      else
          turb_source(i,3) = sqrt_gamma(i)*turb_source(i,3)
          turb_source(i,6) = sqrt_gamma(i)*turb_source(i,6)
      endif
   enddo

end subroutine turbulence_sources

subroutine Brunt_Vaisala(dts)

   use GR1D_module
   use Grad_module
   implicit none

   !local
   integer i
   real*8 dts
   real*8 v_turb_seed, h
   real*8 rho_grad(n1), press_grad(n1), v_grad(n1)
   real*8 lnrho_grad(n1), lnP_grad(n1)

   if (GR) then
       rho_grad = Gradient_3pts(rho(:)*(1.0d0 + eps(:)),x1)
       press_grad = Gradient_3pts(press,x1)
       v_grad = Gradient_5pts(v,x1)
   else
       rho_grad = Gradient_3pts(rho,x1)
       press_grad = Gradient_3pts(press,x1)
       v_grad = Gradient_5pts(v1,x1)

       !lnrho_grad = Gradient_3pts(log(rho),x1)
       !lnP_grad = Gradient_3pts(log(press),x1)
   endif

   do i=ghosts1+1,n1-ghosts1

      if (GR) then
          h = 1.0d0 + eps(i) + press(i)/rho(i) + v_turb(i)**2

          omega2_BV(i) = alp(i)**2/(rho(i)*h*X(i)**2)* &
                (dphidr(i) - v1(i)* v_grad(i)) * &
                (rho_grad(i) - press_grad(i)/cs2(i))
      else
          omega2_BV(i) = (dphidr(i) - v(i)* v_grad(i))/rho(i)* &
                (rho_grad(i) - press_grad(i)/cs2(i))
      endif
      call Get_Mixing_Length(x1(i), press(i), rho(i), dphidr(i), Lambda_MLT(i))

      if (omega2_BV(i) .gt. 0.0d0) then
          v_turb_seed = dts*omega2_BV(i)*min(Lambda_MLT(i), x1(i))
          v_turb(i) = max(v_turb(i),v_turb_seed)
      endif
      v_turb(i) = sqrt(omega2_BV(i))*Lambda_MLT(i)

      if ( x1(i) / length_gf .lt. 5.0d5 ) v_turb(i) = 0.0d0

      ! force covective velocity to be zero in the 5 points stencil around the shock,
      ! i.e. do not convect through the shock
      if ( ishock(1)-2 <= i-ghosts1 .and. i-ghosts1 <= ishock(1)+2 ) v_turb(i) = 0.0d0

   enddo

end subroutine Brunt_Vaisala

subroutine Get_Dissipation_Length(r, p, rho, dphidr, vturb, BV2, Lambda)

  use GR1D_module, only: alpha_turb, explosion_reached, length_gf

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


subroutine Get_Mixing_Length(r, p, rho, dphidr, Lambda)

  use GR1D_module, only: alpha_turb, explosion_reached, length_gf

  implicit none

  real(8), intent(in)  :: r, p, rho, dphidr
  real(8), intent(out) :: Lambda
  real(8) :: max_Lambda

  Lambda = alpha_turb * p / (rho * dphidr)

end subroutine Get_Mixing_Length

