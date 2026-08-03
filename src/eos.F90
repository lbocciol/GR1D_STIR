!-*-f90-*-
subroutine eos_full(i,xrhoi,    &
     xtemp,                   &
     xye,                     &
     xenri,                   &
     xprs,                    &
     xprs_th,                 & 
     xent,                    &
     xcs2,                    &
     xdedt,                   &
     xdpderho,                &
     xdpdrhoe,                & 
     xxa,xxh,xxn,xxp,         & 
     xabar,xzbar,             &
     xmu_e,xmu_n,xmu_p,xmunu,& 
     keytemp,keyerr,eoskey,rfeps)

  use GR1D_module, only: clite,rho_gf,press_gf,eps_gf,GR,n1,atmo
  use atmos
  implicit none

  integer,intent(in)    :: i
  real*8, intent(inout) :: xrhoi  ! inout because of atmosphere
  real*8, intent(in)    :: xye
  real*8, intent(inout) :: xtemp,xenri,xent
  real*8, intent(out)   :: xprs,xprs_th,xcs2,xdedt
  real*8, intent(out)   :: xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp
  real*8, intent(out)   :: xabar,xzbar,xmu_e,xmu_n,xmu_p,xmunu
  real*8, intent(in)    :: rfeps
  integer, intent(in)   :: keytemp
  integer, intent(in)   :: eoskey
  integer, intent(out)  :: keyerr

  real*8 :: xrho, xenr
  

  ! convert to CGS:
  xrho = xrhoi/rho_gf
  xenr = xenri/eps_gf

!  if(keyerr.eq.0) then
     if(xrhoi.le.atmo_rho_thr) then
        call atmos_eos(i,xrhoi,xprs,xenri,xcs2)
        return
     endif
!  endif

  if(eoskey.eq.1) then
     ! hybrid EOS
     call hybrid_eos(xrho,xenr,xprs,xprs_th,&
          xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs_th*press_gf
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.4) then
     ! Gamma-Law EOS
     call ideal_eos(xrho,xenr,xprs,xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.2) then
     ! Poly EOS
     call poly_eos(xrho,xenr,xprs,xdpdrhoe,xdpderho,xcs2,keytemp)
     xprs = xprs*press_gf
     xprs_th = xprs
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     xent = 0.0d0
     xdedt = 0.0d0
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
     xmunu = 0.0d0
     if(keytemp.eq.1) then
        xenri = xenr*eps_gf
     endif
     keyerr = 0
  else if(eoskey.eq.3) then
#if HAVE_NUC_EOS
#ifdef HAVE_BURN
     ! composite finite-T nuc_eos <-> Helmholtz EOS (blended by temperature);
     call nuc_helm_eos_short(i,xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)

     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
#else
     call nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt,&
          xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)
     xxa = 0.0d0
     xxh = 0.0d0
     xxp = 0.0d0
     xxn = 0.0d0
     xabar = 0.0d0
     xzbar = 0.0d0
     xmu_e = 0.0d0
     xmu_n = 0.0d0
     xmu_p = 0.0d0
#endif
     xprs = xprs*press_gf
     xprs_th = 0.0d0
     if (GR) then
        xcs2 = xcs2/(clite**2*(1.0d0+xprs/xrhoi+xenri))
     else
        xcs2 = xcs2/clite**2
     endif
     xdpderho = xdpderho*press_gf/eps_gf
     xdpdrhoe = xdpdrhoe*press_gf/rho_gf
     if(keytemp.eq.1.or.keytemp.eq.2) then
        xenri = xenr*eps_gf
     endif
     if(keyerr.eq.667) then
        if(i.ne.n1) then
           if(atmo(i+1).ne.0) then
              atmo(i) = 1
              call atmos_eos(i,xrhoi,xprs,xenri,xcs2)
              keyerr = 0
              stop "eh"
           endif
        endif
     endif
#else
     stop "eoskey=3 impossible, since NUC_EOS not present."
#endif
  else
     write(6,*) "eoskey ",eoskey," not implemented with eos_full"
     stop "This is bad! Fix me please!"
  endif


end subroutine eos_full

subroutine eos(i,ri,tio,y,eio,xx,keytemp,keyerr,eosflag,eoskey,rfeps)

  use GR1D_module, only: clite,rho_gf,press_gf,eps_gf,GR,n1,atmo
  use atmos
  implicit none

  integer, intent(in) :: i
  real*8, intent(in) :: y,rfeps
  real*8, intent(inout) :: ri,eio, tio
  real*8, intent(out)    :: xx

  real*8 r,e,p_th,tp
  real*8 px,sx,cs2x,gammax,dedt,ex
  integer eosflag,eoskey,keytemp,keyerr,keyerrt
  !internal
  real*8 prs,soundsqr
  real*8 dpdrho,dpde

#if HAVE_NUC_EOS
  real*8 xmunu,xent,xdedt
  integer atmo_next
#endif

  xx = 0.0d0

  ! convert to CGS:
  r = ri/rho_gf
  e = eio/eps_gf

! eosflags:                                                                     
! 1 --> pressure                                                                
! 2 --> dpdrhoe                                                                 
! 3 --> dpderho                                                                 
! 4 --> eps                                                                     
! 5 --> temp                                                                    
! 6 --> cs2                                                                     
! 7 --> gamma                                                                   
! 8 --> entropy                                                                 
! 9 --> munu             

! eoskey:
! 1 --> Hybrid EOS
! 2 --> Poly EOS
! 3 --> finite-T EOS
! 4 --> Gamma-Law EOS

     if(ri.le.atmo_rho_thr) then
        call atmos_eos(i,ri,prs,eio,soundsqr)
        
        select case (eosflag)
        case(1)
           !pressure
           xx=prs
        case(2)
           stop "eh 2"

        case(3)
           stop "eh 3"

        case(4)
           xx=eio

        case(5)
           stop "ehbah"

        case(6)
           xx=soundsqr
 
        case(7)
           xx=soundsqr*ri/prs

        case(8)
           xx=xx

        case(9)
           stop "eosflag=9 not implemented for hybrid eos"

        case default
           stop "eosflag not implemented for hybrid eos"

        end select
        
        return
     endif

  if(eoskey.eq.1) then

     select case (eosflag)

        case(1)
           !pressure
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf
           tio=p_th
!           write(*,"(1P10E15.6)") t
        case(2)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf

        case(3)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for hybrid eos"

        case(6)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           call hybrid_eos(r,e,prs,p_th,dpdrho,dpde,soundsqr,keytemp)
           xx=soundsqr*r/prs

        case(8)
           stop "eosflag=8 not implemented for hybrid eos"

        case(9)
           stop "eosflag=9 not implemented for hybrid eos"

        case default
           stop "eosflag not implemented for hybrid eos"

     end select

  else if(eoskey.eq.2) then

     select case (eosflag)

        case(1)
           !pressure
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf

        case(2)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf

        case(3)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for polytropic eos"

        case(6)
           call poly_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           stop "eosflag=7 not implemented for polytropic eos"

        case(8)
           stop "eosflag=8 not implemented for polytropic eos"

        case(9)
           stop "eosflag=9 not implemented for polytropic eos"

        case default
           stop "eosflag not implemented for polytropic eos"

     end select

  else if(eoskey.eq.4) then

     select case (eosflag)

        case(1)
           !pressure
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=prs*press_gf

        case(2)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpdrho*press_gf/rho_gf


        case(3)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=dpde*press_gf/eps_gf

        case(4)
           keytemp=1
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           xx=e*eps_gf
           eio = xx

        case(5)
           stop "eosflag=5 not implemented for ideal fluid"

        case(6)
           call ideal_eos(r,e,prs,dpdrho,dpde,soundsqr,keytemp)
           if (GR) then
              xx=soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx=soundsqr/clite**2
           endif

        case(7)
           stop "eosflag=7 not implemented for ideal fluid"

        case(8)
           stop "eosflag=8 not implemented for ideal fluid"

        case(9)
           stop "eosflag=9 not implemented for ideal fluid"

        case default
           stop "eosflag not implemented for ideal fluid"

     end select

  else if(eoskey.eq.3) then
#if HAVE_NUC_EOS
     ! one backend evaluation, then pick the requested quantity (was 9 separate
     ! calls; consolidated so the composite EOS logic lives in one place)
     if(eosflag.eq.4) keytemp = 1   ! eps requested -> evaluate with T known

#ifdef HAVE_BURN
     ! composite finite-T nuc_eos <-> Helmholtz EOS (blended by temperature)
     call nuc_helm_eos_short(i,r,tio,y,e,prs,xent,soundsqr,xdedt,&
          dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
#else
     call nuc_eos_short(r,tio,y,e,prs,xent,soundsqr,xdedt,&
          dpde,dpdrho,xmunu,keytemp,keyerr,rfeps)
#endif
     if(keytemp.eq.1.or.keytemp.eq.2) eio = e*eps_gf

     select case (eosflag)
        case(1)   !pressure
           xx = prs*press_gf
        case(2)   !dpdrhoe
           xx = dpdrho*press_gf/rho_gf
        case(3)   !dpderho
           xx = dpde*press_gf/eps_gf
        case(4)   !eps
           xx = e*eps_gf
           eio = xx
        case(6)   !cs2
           if (GR) then
              xx = soundsqr/(clite**2*(1.0d0+(prs*press_gf)/ri+eio))
           else
              xx = soundsqr/clite**2
           endif
        case(7)   !gamma
           xx = soundsqr*r/prs
        case(8)   !entropy
           xx = xent
        case(9)   !munu
           xx = xmunu
        case(5)
           stop "eosflag=5 not implemented for eoskey = 3"
        case default
           stop "eosflag not implemented for eoskey = 3"
        end select

        if(keyerr.eq.667) then
           if(i.ne.n1) then
              ! atmo(i+1) belongs to another zone: atomic, since con2prim may
              ! run the zone loop threaded and atmos_eos writes atmo there
              !$omp atomic read
              atmo_next = atmo(i+1)
              if(atmo_next.ne.0) then
                 !$omp atomic write
                 atmo(i) = 1
                 keyerr = 0
                 call atmos_eos(i,ri,prs,eio,soundsqr)
                 select case (eosflag)
                 case(1)
                    !pressure
                    xx=prs
                 case(2)
                    stop "eh 4"
                    
                 case(3)
                    stop "eh 5"

                 case(4)
                    xx=eio
                    
                 case(5)
                    stop "ehbah"
                    
                 case(6)
                    xx=soundsqr
                    
                 case(7)
                    xx=soundsqr*ri/prs
                    
                 case(8)
                    xx=xx
                    
                 case(9)
                    stop "eosflag=9 not implemented for hybrid eos"
                    
                 case default
                    stop "eosflag not implemented for hybrid eos"
                    
                 end select
              endif
           endif
        endif

#endif

  else
     write(6,*) "eoskey: ", eoskey
     call flush(6)
     stop "eos choice not implemented, sorry..."

  endif

end subroutine eos

module hybrid_eos_module

  real*8 gamma1
  real*8 gamma2
  real*8 gammath
  real*8 K1,K2,E1,E2,E3
  real*8 rhonuc

end module hybrid_eos_module

subroutine hybrid_eos(rho,enr,prs,pth,dpdrho_o,dpde_o,soundsqr,keytemp)

  use hybrid_eos_module
  implicit none

! Input
  real*8 rho,enr,prs,soundsqr
  integer keytemp
! Output
  real*8 dpdrho_o,dpde_o


! Local 
!
  real*8 pco,pth,eth,dpth_drho,dpth_denr
  real*8 Gx,Ex,Kx,Ex3,dpco_drho,dpco_denr
  real*8 up,dp_drho,dpde

  if(keytemp.eq.1) then
!	energy wanted
     prs = K1*rho**gamma1 
     enr=prs/rho/(gamma1-1.d0)
     soundsqr=gamma1*prs/rho
  endif

  if(rho .lt. rhonuc) then
     Kx=K1
     Ex=E1
     Gx=gamma1
     Ex3=0.d0
  else
     Kx=K2
     Ex=E2
     Gx=gamma2
     Ex3=E3
  endif

!                Thermal
    up=Ex*rho**Gx+Ex3*rho
    eth = enr*rho - up
    pth = (gammath - 1)*eth
    dpth_drho=(gammath - 1)*(enr-Ex*Gx*rho**(Gx-1.d0)-Ex3)
    dpth_denr=(gammath - 1)*rho

    if (pth.lt.0.0d0) then
       pth = 0.0d0
       dpth_drho = 0.0d0
       dpth_denr = 0.0d0
    endif


!                 Cold 
    pco = Kx*rho**Gx
    dpco_drho=Gx*pco/rho
    dpco_denr=0.d0
!
    prs=pco+pth
    dpde=dpco_denr+dpth_denr
    dp_drho=dpco_drho+dpth_drho
    soundsqr= dp_drho+dpde*prs/rho**2

    dpdrho_o = dp_drho
    dpde_o = dpde

end subroutine hybrid_eos


subroutine init_hybrid_eos

  use hybrid_eos_module

  implicit none

  E1 = K1/(gamma1-1.d0)
  E2 = (gamma1 - 1.d0)/(gamma2-1.d0)*E1*rhonuc**(gamma1-gamma2)
  K2 = (gamma2 - 1.d0)*E2
  E3 = (gamma2 - gamma1)/(gamma2-1.d0)*E1*rhonuc**(gamma1-1.d0)

end subroutine init_hybrid_eos

module poly_eos_module

  real*8 polygamma
  real*8 polyK

end module poly_eos_module

subroutine poly_eos(rho,enr,prs,dpdrho,dpde,soundsqr,keytemp)

  use poly_eos_module
  implicit none

! Input/Output
  real*8 rho,enr,prs,soundsqr
  integer keytemp
  real*8 dpdrho,dpde

! Local 
!
  prs = polyK * rho**polygamma
  soundsqr = polygamma*polyK*rho**(polygamma-1.d0)

  dpdrho = polyK * rho**(polygamma - 1) * polygamma
  dpde = 0.0d0

  if(keytemp.eq.1) then
!	energy wanted
     enr=prs/rho/(polygamma-1.d0)
     soundsqr=polygamma*prs/rho
  endif


end subroutine poly_eos

module ideal_eos_module

  real*8 idealgamma
  real*8 idealK1

end module ideal_eos_module

subroutine ideal_eos(rho,enr,prs,dpdrho,dpde,soundsqr,keytemp)

  use ideal_eos_module
  implicit none

! Input/Output
  real*8 rho,enr,prs,soundsqr
  real*8 dpdrho,dpde
  integer keytemp


  if(keytemp.eq.1) then
!	energy wanted
     prs=idealK1*rho**(idealgamma)
     enr=prs/rho/(idealgamma-1.d0)
     soundsqr=idealgamma*prs/rho

  endif
     
  prs = (idealgamma - 1.0d0) *rho*enr

  dpde = (idealgamma - 1.0d0 ) * rho
  dpdrho = (idealgamma - 1.0d0 ) * enr

  soundsqr= dpdrho+dpde*prs/rho**2


end subroutine ideal_eos

subroutine atmos_eos(i,xrho,xprs,xenr,xcs2)

  use atmos
  use GR1D_module, only: rho_gf,press_gf, &
       v1,v,W,atmo
  implicit none

  integer, intent(in) :: i
  real*8 :: xrho,xprs, xenr, xcs2
  real*8,parameter :: idealK1 =  1.2435d15 * (0.5d0**(4.d0/3.d0))
  real*8,parameter :: idealgamma = 1.66666666666d0

  xrho = atmo_rho
  xprs=idealK1*((atmo_rho/rho_gf)**(idealgamma))*press_gf
  xenr=xprs/atmo_rho/(idealgamma-1.d0)
  xcs2=idealgamma*xprs/atmo_rho
  v(i) = 0.0d0
  v1(i) = 0.0d0
  W(i) = 1.0d0
  ! atomic: atmo(i) may be read as atmo(i+1) by the neighbouring zone's
  ! thread in eos()'s keyerr=667 recovery path
  !$omp atomic write
  atmo(i) = 1

end subroutine atmos_eos

#if HAVE_NUC_EOS
#ifdef HAVE_BURN
! ===========================================================================
! Composite finite-T nuclear EOS <-> Helmholtz EOS (Perego et al. 2015;
! Navo et al. 2023).
!
! Drop-in replacement for nuc_eos_short with ONE extra leading argument: the
! radial zone index i (used only to fetch the local network composition
! Yion(:,i)).  Every other argument and unit (CGS in/out) matches
! nuc_eos_short exactly, so the caller is unaware of the blend.
!
! Dispatch is purely by TEMPERATURE (no density or energy rule; the regimes
! occur at the appropriate densities by construction):
!   T >= T_eos_high                 -> nuc_eos only    (w = 1)
!   T <= T_eos_low                  -> Helmholtz only  (w = 0)
!   T_eos_low < T < T_eos_high      -> blend BOTH with ONE linear weight w(T),
!       applied identically to every returned quantity so the state stays a
!       single consistent thermodynamic point.
! For keytemp=0 the single T whose blended energy e_blend(T) matches the
! target is root-found, starting from the incoming xtemp; at every iterate
! only the backend(s) selected by w(T) are evaluated, so the solver crosses
! the seams naturally as T moves.
!
! Any failure STOPs with a diagnostic: it signals a physics/modeling problem
! (e.g. the network composition drifting far from the table NSE) that must
! not be hidden by a fallback.
! ===========================================================================
subroutine nuc_helm_eos_short(i,xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
     xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)

  use GR1D_module, only: Yion, T_eos_high, T_eos_low, temp_mev_to_kelvin
  use composition, only: nspec, zion
  implicit none

  integer, intent(in)    :: i, keytemp
  real*8,  intent(in)    :: xrho, xye, rfeps
  real*8,  intent(inout) :: xtemp, xenr
  real*8,  intent(out)   :: xprs, xent, xcs2, xdedt
  real*8,  intent(out)   :: xdpderho, xdpdrhoe, xmunu
  integer, intent(out)   :: keyerr

  real*8 :: Y(nspec), abar, zbar, e_offset, w, tk
  real*8 :: t_lo_mev, t_hi_mev, t_mev, f, dfdt
  ! nuc_eos branch outputs (local copies so xtemp/xenr aren't clobbered)
  real*8 :: t_n, e_n, p_n, ent_n, cs2_n, dedt_n, dpde_n, dpdr_n, munu_n
  ! Helmholtz branch outputs
  real*8 :: tk_h, e_h, p_h, ent_h, cs2_h, dedt_h, dpde_h, dpdr_h
  real*8 :: t_pos, t_neg, t_new
  logical :: have_pos, have_neg
  integer :: it
  integer, parameter :: maxit = 200
  real*8 :: t_hist(maxit), f_hist(maxit), w_hist(maxit), dfdt_hist(maxit)

  keyerr = 0
  xmunu  = 0.0d0

  ! entropy-mode inversion (keytemp=2) is wired up only for nuc_eos
  if (keytemp .eq. 2) then
     STOP 'entropy inversion not implemented for Helmholtz; use nuc_eos only'
  end if

  Y    = Yion(:,i)
  abar = 1.0d0 / sum(Y)
  zbar = abar * sum(zion * Y)
  call get_energy_offset(Y, e_offset)

  ! GR1D carries temperature in MeV (nuc_eos convention); the thresholds
  ! T_eos_high/T_eos_low and the Helmholtz table are in Kelvin -> convert for
  ! the weight and for the Helmholtz calls; nuc_eos is always called in MeV.
  tk = xtemp * temp_mev_to_kelvin

  ! =========================================================================
  ! keytemp = 1 (T known): regime selected from the known T
  ! =========================================================================
  if (keytemp .eq. 1) then
     call blend_weight(tk, w)

     if (w .le. 0.0d0) then
        tk_h = tk
        e_h  = xenr
        call eval_helm(xrho,tk_h,abar,zbar,e_offset,e_h, &
             1,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe)
        xenr = e_h
        return
     end if

     if (w .ge. 1.0d0) then
        call nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
             xdpderho,xdpdrhoe,xmunu,1,keyerr,rfeps)
        xmunu = 0.0d0
        return
     end if

     ! transition: evaluate BOTH backends at the known T, blend linearly
     t_n = xtemp
     call nuc_eos_short(xrho,t_n,xye,e_n,p_n,ent_n,cs2_n,dedt_n, &
          dpde_n,dpdr_n,munu_n,1,keyerr,rfeps)
     tk_h = tk
     call eval_helm(xrho,tk_h,abar,zbar,e_offset,e_h, &
          1,p_h,ent_h,cs2_h,dedt_h,dpde_h,dpdr_h)

     xenr     = w*e_n    + (1.0d0-w)*e_h
     xprs     = w*p_n    + (1.0d0-w)*p_h
     xent     = w*ent_n  + (1.0d0-w)*ent_h
     xcs2     = w*cs2_n  + (1.0d0-w)*cs2_h
     xdedt    = w*dedt_n + (1.0d0-w)*dedt_h
     xdpderho = w*dpde_n + (1.0d0-w)*dpde_h
     xdpdrhoe = w*dpdr_n + (1.0d0-w)*dpdr_h
     return
  end if

  ! =========================================================================
  ! keytemp = 0 (eps known): Newton on the piecewise blended energy
  !     f(T) = w(T)*e_nuc(rho,T) + (1-w(T))*e_helm(rho,T) - xenr
  ! over the FULL temperature range, from the incoming xtemp (in evolution
  ! the previous T, which sits next to the root).  The regime is decided by
  ! the temperature the solver is at -- only the backend(s) selected by w(T)
  ! are evaluated, so cold zones never touch nuc_eos and hot zones never
  ! touch Helmholtz -- and the iteration crosses the seams naturally.
  ! For the linear weight dw/dT is constant inside the window and both
  ! backends return dedt, so the derivative is analytic
  !     df/dT = w*dedt_n + (1-w)*dedt_h*K/MeV + (e_n - e_h)/(T_hi - T_lo)
  ! (T in MeV; dedt_n is per MeV, dedt_h = cv is per Kelvin).
  ! =========================================================================
  t_lo_mev = T_eos_low  / temp_mev_to_kelvin
  t_hi_mev = T_eos_high / temp_mev_to_kelvin
  t_mev    = xtemp

  have_pos = .false.
  have_neg = .false.

  do it = 1, maxit
     tk = t_mev * temp_mev_to_kelvin
     call blend_weight(tk, w)

     if (w .gt. 0.0d0) then
        t_n = t_mev
        call nuc_eos_short(xrho,t_n,xye,e_n,p_n,ent_n,cs2_n,dedt_n, &
             dpde_n,dpdr_n,munu_n,1,keyerr,rfeps)
     end if
     if (w .lt. 1.0d0) then
        tk_h = tk
        call eval_helm(xrho,tk_h,abar,zbar,e_offset,e_h, &
             1,p_h,ent_h,cs2_h,dedt_h,dpde_h,dpdr_h)
     end if

     if (w .le. 0.0d0) then
        f    = e_h - xenr
        dfdt = dedt_h * temp_mev_to_kelvin
     else if (w .ge. 1.0d0) then
        f    = e_n - xenr
        dfdt = dedt_n
     else
        f    = w*e_n + (1.0d0-w)*e_h - xenr
        dfdt = w*dedt_n + (1.0d0-w)*dedt_h*temp_mev_to_kelvin &
             + (e_n - e_h)/(t_hi_mev - t_lo_mev)
     end if

     ! Converged when the residual meets the tolerance, or when the sign-
     ! change bracket has collapsed to machine precision in T: there the
     ! residual is pure floating-point cancellation noise (the energies
     ! carry ~1e17 erg/g composition offsets, so f cannot be resolved below
     ! ~1e3 erg/g while xenr itself can pass near zero).  A target with NO
     ! root never brackets and still STOPs below.
     if (abs(f) .le. rfeps*abs(xenr) .or. &
         (have_pos .and. have_neg .and. &
          abs(t_pos - t_neg) .le. 4.0d0*spacing(t_mev))) then
        xtemp = t_mev
        if (w .le. 0.0d0) then
           xprs     = p_h
           xent     = ent_h
           xcs2     = cs2_h
           xdedt    = dedt_h
           xdpderho = dpde_h
           xdpdrhoe = dpdr_h
        else if (w .ge. 1.0d0) then
           xprs     = p_n
           xent     = ent_n
           xcs2     = cs2_n
           xdedt    = dedt_n
           xdpderho = dpde_n
           xdpdrhoe = dpdr_n
        else
           xprs     = w*p_n    + (1.0d0-w)*p_h
           xent     = w*ent_n  + (1.0d0-w)*ent_h
           xcs2     = w*cs2_n  + (1.0d0-w)*cs2_h
           xdedt    = w*dedt_n + (1.0d0-w)*dedt_h
           xdpderho = w*dpde_n + (1.0d0-w)*dpde_h
           xdpdrhoe = w*dpdr_n + (1.0d0-w)*dpdr_h
        end if
        return
     end if

     t_hist(it) = t_mev
     f_hist(it) = f
     w_hist(it) = w
     dfdt_hist(it) = dfdt

     ! Track the sign-change bracket: f(T) is continuous but has derivative
     ! kinks at the window edges (the w' term), where plain Newton can cycle.
     if (f .gt. 0.0d0) then
        t_pos = t_mev
        have_pos = .true.
     else
        t_neg = t_mev
        have_neg = .true.
     end if

     if (have_pos .and. have_neg) then
        ! bracketed: bisect.  Near the window seams the analytic df/dT
        ! (separately tabulated dedt) can underestimate the true slope of
        ! the interpolated e(T), making Newton overshoot the root and cycle;
        ! bisection is immune (root-finder hygiene, not a physics fallback;
        ! a target with no root still has no sign change and STOPs below).
        t_new = 0.5d0*(t_pos + t_neg)
     else
        ! no bracket yet: Newton limited to a factor of 2 per iteration
        t_new = max(0.5d0*t_mev, min(t_mev - f/dfdt, 2.0d0*t_mev))
     end if
     t_mev = t_new
  end do

  write(*,*) 'nuc_helm_eos_short: blended-energy Newton did not converge'
  write(*,*) '  zone i     = ', i
  write(*,*) '  rho        = ', xrho
  write(*,*) '  ye         = ', xye
  write(*,*) '  target eps = ', xenr
  write(*,*) '  guess T    = ', xtemp, ' MeV'
  write(*,*) '  last T     = ', t_mev, ' MeV'
  write(*,*) '  last f     = ', f
  write(*,*) '  abar, zbar = ', abar, zbar
  write(*,*) '  e_offset   = ', e_offset
  write(*,*) '  Yion       = ', Y
  write(*,*) '  iteration history (it, T[MeV], w, f, dfdt):'
  do it = 1, maxit
     write(*,'(i5,1p4e24.15)') it, t_hist(it), w_hist(it), f_hist(it), dfdt_hist(it)
  end do
  STOP 'nuc_helm_eos_short: blended keytemp=0 inversion failed'

end subroutine nuc_helm_eos_short

! ---------------------------------------------------------------------------
! Single Helmholtz evaluation in the SAME CGS in/out convention as
! nuc_eos_short.  The energy offset [erg/g] (from get_energy_offset: the
! nuclear mass-excess per gram of the local composition) is supplied as the
! state's e_offset so the Helmholtz energy zero-point matches the finite-T
! nuclear EOS (see eps-energy convention).
!   keytemp = 1 : temp known -> enr returned
!   keytemp = 0 : enr known  -> temp solved (FullHelmEOS STOPs on any
!                 non-convergence or unreachable target, so the result is
!                 trustworthy)
! cs2 is returned as the NEWTONIAN gam1*p/rho [cm^2/s^2]; the relativistic
! correction is applied uniformly by the caller, exactly as for nuc_eos.
! ---------------------------------------------------------------------------
subroutine eval_helm(rho,temp,abar,zbar,e_offset,enr,keytemp, &
     prs,ent,cs2,dedt,dpderho,dpdrhoe)

  use wlHelmholtzEOS, only: HelmholtzStateType, HelmEOS, &
                            eos_input_rt, eos_input_re
  implicit none

  real*8,  intent(in)    :: rho, abar, zbar, e_offset
  real*8,  intent(inout) :: temp, enr
  integer, intent(in)    :: keytemp
  real*8,  intent(out)   :: prs, ent, cs2, dedt, dpderho, dpdrhoe

  type(HelmholtzStateType) :: st

  st % rho      = rho
  st % T        = temp
  st % abar     = abar
  st % zbar     = zbar
  st % ye       = zbar/abar
  st % e_offset = e_offset

  if (keytemp .eq. 1) then
     call HelmEOS(eos_input_rt, st)
  else
     st % e = enr
     call HelmEOS(eos_input_re, st)
     temp = st % T
  end if
  enr = st % e

  prs     = st % p
  ent     = st % s
  cs2     = st % gam1 * st % p / st % rho   ! Newtonian cs^2 [cm^2/s^2]
  dedt    = st % cv
  dpderho = st % dpde
  dpdrhoe = st % dpdr_e

end subroutine eval_helm

! ---------------------------------------------------------------------------
! Linear temperature blend weight w(T) in [0,1], T in Kelvin
! (Perego et al. 2015; Navo et al. 2023).
!   T >= T_eos_high -> 1 (pure nuc_eos)
!   T <= T_eos_low  -> 0 (pure Helmholtz)
!   between         -> linear
! ---------------------------------------------------------------------------
subroutine blend_weight(tk, w)

  use GR1D_module, only: T_eos_high, T_eos_low
  implicit none

  real*8, intent(in)  :: tk
  real*8, intent(out) :: w

  if (tk .ge. T_eos_high) then
     w = 1.0d0
  else if (tk .le. T_eos_low) then
     w = 0.0d0
  else
     w = (tk - T_eos_low) / (T_eos_high - T_eos_low)
  end if

end subroutine blend_weight

! ---------------------------------------------------------------------------
! Nuclear rest-mass contribution of the composition [erg/g], added to the
! Helmholtz thermal energy so its zero-point matches the tabulated nuc_eos:
!     e_offset = (1/m_u) sum_i Dm_i Y_i ,   Dm_i = m_i - A_i m_u  (mass excess)
! nuclei_mass_excess is in MeV and includes the free nucleons.
subroutine get_energy_offset(Y, e_offset)

  use GR1D_module, only: mev_to_erg, avo
  use composition, only: nspec, nuclei_mass_excess
  implicit none

  real*8, intent(in) :: Y(nspec)
  real*8, intent(out) :: e_offset

  e_offset = mev_to_erg * avo * sum(nuclei_mass_excess * Y)

end subroutine get_energy_offset
#endif
#endif
