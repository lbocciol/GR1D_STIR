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
! ---------------------------------------------------------------------------
! Shared state for the composite nuc_eos <-> Helmholtz EOS.
! ---------------------------------------------------------------------------
module eos_blend_module
  implicit none
  ! regime of the composite EOS, decided ONCE per call (see nuc_helm_eos_short)
  integer, parameter :: REG_COLD = 0   ! pure Helmholtz
  integer, parameter :: REG_BLEND = 1  ! linear blend of both backends
  integer, parameter :: REG_HOT = 2    ! pure tabulated nuc_eos
  ! rate limit for the inverted-seam warning (see warn_inverted_window)
  integer, parameter :: max_inverted_warn = 20
  integer :: n_inverted_warn = 0
end module eos_blend_module

! ===========================================================================
! Composite finite-T nuclear EOS <-> Helmholtz EOS (Perego et al. 2015;
! Navo et al. 2023).
!
! Drop-in replacement for nuc_eos_short with ONE extra leading argument: the
! radial zone index i (used only to fetch the local network composition
! Yion(:,i)).  Every other argument and unit (CGS in/out) matches
! nuc_eos_short exactly, so the caller is unaware of the blend.
!
! The regimes are set by TEMPERATURE:
!   T >= T_eos_high                 -> nuc_eos only    (w = 1)
!   T <= T_eos_low                  -> Helmholtz only  (w = 0)
!   T_eos_low < T < T_eos_high      -> blend BOTH with ONE linear weight w(T),
!       applied identically to every returned quantity so the state stays a
!       single consistent thermodynamic point.
!
! The regime is decided FIRST, then the state is obtained WITHIN that regime:
!   keytemp=1: T is known, so the regime follows directly from it.
!   keytemp=0: the regime follows from the two SEAM energies
!              e_lo = e_helm(rho,T_eos_low) and e_hi = e_nuc(rho,T_eos_high)
!              (see classify_regime), and then
!                REG_COLD  -> Helmholtz inverts itself (eos_input_re); there is
!                             NO outer Newton here at all.  This matters: the
!                             Helmholtz energy carries the composition offset
!                             (~1e17-1e18 erg/g) while eps itself is ~1e11, so
!                             a residual formed as e_helm + e_offset - eps is
!                             quantised far above any sensible tolerance.
!                             eos_input_re strips the offset ONCE and converges
!                             on the step size in T, which is offset-free.
!                REG_HOT   -> nuc_eos_short's own findtemp, i.e. exactly what a
!                             HAVE_BURN=0 build does.
!                REG_BLEND -> bracketed root find on [T_eos_low, T_eos_high],
!                             where classification has already guaranteed a
!                             sign change (see blend_invert).
!
! Two deliberate deviations from a HAVE_BURN=0 build:
!   * xmunu is hard-zeroed in every branch (nuchem stays 0 under HAVE_BURN).
!   * below eos_rhomin*1.2 the cold branch is forced regardless of T, because
!     nuc_eos_short silently substitutes nuc_low_eos there (a K*rho^1.41
!     polytrope with an unrelated energy zero point and xent hardwired to 4),
!     which would poison the seam classification.  Helmholtz with the real
!     composition is strictly better than an arbitrary polytrope.
!
! Any failure STOPs with a diagnostic: it signals a physics/modeling problem
! (e.g. the network composition drifting far from the table NSE) that must
! not be hidden by a fallback.
! ===========================================================================
subroutine nuc_helm_eos_short(i,xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
     xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)

  use GR1D_module, only: Yion, T_eos_high, T_eos_low, temp_mev_to_kelvin
  use composition, only: nspec, zion
  use eosmodule, only: eos_rhomin
  use eos_blend_module, only: REG_COLD, REG_BLEND, REG_HOT
  implicit none

  integer, intent(in)    :: i, keytemp
  real*8,  intent(in)    :: xrho, xye, rfeps
  real*8,  intent(inout) :: xtemp, xenr
  real*8,  intent(out)   :: xprs, xent, xcs2, xdedt
  real*8,  intent(out)   :: xdpderho, xdpdrhoe, xmunu
  integer, intent(out)   :: keyerr

  real*8 :: Y(nspec), abar, zbar, e_offset
  real*8 :: t_lo_mev, t_hi_mev, e_lo, e_hi, e_target
  real*8 :: tk_h, e_h, dedt_h, dfdt
  integer :: regime

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
  ! T_eos_high/T_eos_low and the Helmholtz table are in Kelvin.
  t_lo_mev = T_eos_low  / temp_mev_to_kelvin
  t_hi_mev = T_eos_high / temp_mev_to_kelvin

  e_lo = 0.0d0
  e_hi = 0.0d0

  ! ---------------- regime decision (once, for both keytemps) --------------
  if (xrho .lt. eos_rhomin*1.2d0) then
     ! nuclear table undefined here; nuc_eos_short would return a polytrope
     regime = REG_COLD
  else if (keytemp .eq. 1) then
     if (xtemp .le. t_lo_mev) then
        regime = REG_COLD
     else if (xtemp .ge. t_hi_mev) then
        regime = REG_HOT
     else
        regime = REG_BLEND
     end if
  else
     call classify_regime(i,xrho,xtemp,xye,xenr,abar,zbar,e_offset,Y, &
          t_lo_mev,t_hi_mev,rfeps, regime,e_lo,e_hi)
  end if

  ! ---------------- state, obtained within that regime ---------------------
  select case (regime)

  case (REG_COLD)                       ! pure Helmholtz, its OWN inversion

     if (keytemp .eq. 1) then
        tk_h = xtemp * temp_mev_to_kelvin
        e_h  = xenr
        call eval_helm(xrho,tk_h,abar,zbar,e_offset,e_h, &
             1,xprs,xent,xcs2,dedt_h,xdpderho,xdpdrhoe)
        xenr = e_h
     else
        ! Seed for FullHelmEOS' Newton, in Kelvin.  The upper clamp is valid
        ! because classification has already proved the root is at or below
        ! T_eos_low.  The 1.0d4 floor is required: FullHelmEOS measures its
        ! error as ABS((xnew-x)/x), so a zero seed divides by zero; its
        ! factor-of-2 limiter climbs 1e4 -> 1e9 in 17 of its 100 iterations.
        tk_h = min(max(xtemp*temp_mev_to_kelvin, 1.0d4), T_eos_low)
        e_h  = xenr
        call eval_helm(xrho,tk_h,abar,zbar,e_offset,e_h, &
             0,xprs,xent,xcs2,dedt_h,xdpderho,xdpdrhoe)   ! eos_input_re
        xtemp = tk_h / temp_mev_to_kelvin
     end if
     xdedt = dedt_h * temp_mev_to_kelvin    ! cv is per Kelvin -> per MeV

  case (REG_HOT)                        ! identical to a HAVE_BURN=0 build

     call nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
          xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)
     xmunu = 0.0d0                      ! unchanged behaviour, by decision

  case (REG_BLEND)

     e_target = xenr
     if (keytemp .eq. 0) then
        call blend_invert(i,xrho,xye,xenr,abar,zbar,e_offset,Y, &
             t_lo_mev,t_hi_mev,e_lo,e_hi,rfeps, xtemp)
     end if
     call blend_state(xrho,xtemp,xye,abar,zbar,e_offset,rfeps, &
          xenr,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe,dfdt,keyerr)
     ! keytemp=0: hand back the caller's target energy untouched (the
     ! converged blended energy equals it to within the solver tolerance)
     if (keytemp .eq. 0) xenr = e_target

  end select

end subroutine nuc_helm_eos_short

! ---------------------------------------------------------------------------
! keytemp=0 regime decision, from the two SEAM energies
!     e_lo = e_helm(rho, T_eos_low)  (offset included)
!     e_hi = e_nuc (rho, T_eos_high)
!
! NORMAL window (e_lo <= e_hi).  e_blend(T) is monotone, so the seams partition
! the target space exactly and each regime's reachable energies are
!     cold [.., e_lo]    window [e_lo, e_hi]    hot [e_hi, ..]
! Truth table (cold_ok = xenr<=e_lo, hot_ok = xenr>=e_hi):
!   T F -> cold
!   F T -> hot
!   F F -> e_lo < xenr < e_hi: window, sign change GUARANTEED
!   T T -> impossible unless e_lo = e_hi = xenr
!
! INVERTED window (e_lo > e_hi).  The two backends disagree about the energy
! zero point by more than the width of the window -- a 13-species network at
! ye~0.44 puts the Helmholtz zero ~1.2e18 erg/g above the nuc table -- so
! e_blend(T) FALLS across the window and is no longer monotone.  The seam
! energies then no longer partition anything: a target can have up to three
! roots, and, worse, the energy test can route a window root to a branch that
! cannot reach it at all.  In a degenerate zone e_helm(T) is flat (at
! rho ~ 1e8 it varies by <1e-3 between 1e3 K and 1e8 K), so the cold branch's
! reachable range is a thin sliver just below e_lo; a target below that sliver
! but still <= e_lo has its only root in the window, yet reads as cold_ok.
! There the regime is taken from the incoming T instead, which keeps a zone on
! the branch it was on last step and is exactly right for the map_profile
! round trip (keytemp=1 at temp(i), then keytemp=0 on the eps it just made).
! ---------------------------------------------------------------------------
subroutine classify_regime(i,xrho,xtemp,xye,xenr,abar,zbar,e_offset,Y, &
     t_lo_mev,t_hi_mev,rfeps, regime,e_lo,e_hi)

  use composition, only: nspec
  use eos_blend_module, only: REG_COLD, REG_BLEND, REG_HOT
  implicit none

  integer, intent(in)  :: i
  real*8,  intent(in)  :: xrho, xtemp, xye, xenr, abar, zbar, e_offset
  real*8,  intent(in)  :: Y(nspec), t_lo_mev, t_hi_mev, rfeps
  integer, intent(out) :: regime
  real*8,  intent(out) :: e_lo, e_hi

  logical :: cold_ok, hot_ok

  e_lo = 0.0d0
  e_hi = 0.0d0

  ! BOTH seams, always.  Returning REG_COLD as soon as xenr <= e_lo, without
  ! ever looking at e_hi, cannot tell a normal window from an inverted one --
  ! and in an inverted window that test is not a classification at all.
  call seam_energy_cold(xrho,abar,zbar,e_offset, e_lo)
  call seam_energy_hot (xrho,xye,t_hi_mev,rfeps, e_hi)

  cold_ok = (xenr .le. e_lo)
  hot_ok  = (xenr .ge. e_hi)

  if (e_lo .gt. e_hi) then
     ! inverted: e_blend(T) is non-monotone, so take the regime from the
     ! incoming T (temporal continuity) rather than from the energy, and say
     ! so out loud -- a zone that flips branch between steps jumps in energy.
     call warn_inverted_window(i,xrho,xye,xenr,abar,zbar,e_lo,e_hi,e_offset,Y)
     if (xtemp .le. t_lo_mev) then
        regime = REG_COLD
     else if (xtemp .ge. t_hi_mev) then
        regime = REG_HOT
     else
        regime = REG_BLEND
     end if
  else if (cold_ok) then
     regime = REG_COLD
  else if (hot_ok) then
     regime = REG_HOT
  else
     regime = REG_BLEND      ! e_lo < xenr < e_hi: bracketed by construction
  end if

end subroutine classify_regime

! ---------------------------------------------------------------------------
! Helmholtz energy [erg/g] at the cold seam T_eos_low (offset included).
! ---------------------------------------------------------------------------
subroutine seam_energy_cold(xrho,abar,zbar,e_offset, e_lo)

  use GR1D_module, only: T_eos_low
  implicit none

  real*8, intent(in)  :: xrho, abar, zbar, e_offset
  real*8, intent(out) :: e_lo

  real*8 :: tk, p, ent, cs2, dedt, dpde, dpdr

  tk   = T_eos_low
  e_lo = 0.0d0
  call eval_helm(xrho,tk,abar,zbar,e_offset,e_lo,1,p,ent,cs2,dedt,dpde,dpdr)

end subroutine seam_energy_cold

! ---------------------------------------------------------------------------
! Tabulated nuc_eos energy [erg/g] at the hot seam T_eos_high.
! ---------------------------------------------------------------------------
subroutine seam_energy_hot(xrho,xye,t_hi_mev,rfeps, e_hi)

  implicit none

  real*8, intent(in)  :: xrho, xye, t_hi_mev, rfeps
  real*8, intent(out) :: e_hi

  real*8  :: t, p, ent, cs2, dedt, dpde, dpdr, munu
  integer :: keyerr

  t = t_hi_mev
  call nuc_eos_short(xrho,t,xye,e_hi,p,ent,cs2,dedt,dpde,dpdr,munu, &
       1,keyerr,rfeps)

end subroutine seam_energy_hot

! ---------------------------------------------------------------------------
! Rate-limited warning for the inverted-seam pathology (e_lo > e_hi).  A zone
! that flips branch between steps produces a real energy jump, and this is the
! only signal of it.  Not fatal: the physics fix belongs in the network/NSE
! energy zero points, not here.
! ---------------------------------------------------------------------------
subroutine warn_inverted_window(i,xrho,xye,xenr,abar,zbar,e_lo,e_hi,e_offset,Y)

  use composition, only: nspec
  use eos_blend_module, only: n_inverted_warn, max_inverted_warn
  implicit none

  integer, intent(in) :: i
  real*8,  intent(in) :: xrho, xye, xenr, abar, zbar, e_lo, e_hi, e_offset
  real*8,  intent(in) :: Y(nspec)

  !$OMP CRITICAL (eos_inverted_window_warning)
  if (n_inverted_warn .lt. max_inverted_warn) then
     n_inverted_warn = n_inverted_warn + 1
     write(*,*) 'nuc_helm_eos_short: WARNING inverted blend window (e_lo > e_hi)'
     write(*,*) '  zone i     = ', i
     write(*,*) '  rho        = ', xrho
     write(*,*) '  ye (hydro) = ', xye
     write(*,*) '  ye (zbar/abar) = ', zbar/abar
     write(*,*) '  target eps = ', xenr
     write(*,*) '  e_lo, e_hi = ', e_lo, e_hi
     write(*,*) '  abar, zbar = ', abar, zbar
     write(*,*) '  e_offset   = ', e_offset
     write(*,*) '  Yion       = ', Y
     if (n_inverted_warn .eq. max_inverted_warn) then
        write(*,*) '  (further inverted-window warnings suppressed)'
     end if
  end if
  !$OMP END CRITICAL (eos_inverted_window_warning)

end subroutine warn_inverted_window

! ---------------------------------------------------------------------------
! keytemp=0 solve INSIDE the blend window, confined to [T_eos_low, T_eos_high].
!
! Two exit criteria, neither optional:
!   * the relative-T exit.  The residual test ALONE can be unreachable: the
!     blended energy carries the composition offset (~1e17 erg/g, spacing
!     ~1e2) while the tolerance is rfeps*|eps| ~ 1e2, so |f| can be pure
!     cancellation noise that never falls below it.  This is what hung the old
!     solver.
!   * the ULP floor on the residual, for the same reason: 8*spacing(|e|) is
!     ~1e3 erg/g there, i.e. dT ~ 1e-5 K -- thermodynamically irrelevant.
!
! The seam values are used as a bracket WHEN THEY ARE ONE.  They are not in
! general: e_blend(T) = w*e_nuc + (1-w)*e_helm is only monotone if the two
! backends agree on the energy zero point to better than the window's own
! energy span.  When they do not (see classify_regime), e_blend can dip below
! e_hi just inside the hot seam -- w -> 1 there while e_nuc(T) < e_nuc(T_hi) --
! so a target with a perfectly good interior root reads as "outside the
! bracket".  With a bracket this is rtsafe (Newton, bisect on overshoot); with
! no bracket it is a damped Newton from the incoming guess, which in both the
! map_profile round trip and in evolution starts essentially AT the root, and
! which promotes itself to rtsafe the moment a sign change appears.
! ---------------------------------------------------------------------------
subroutine blend_invert(i,xrho,xye,xenr,abar,zbar,e_offset,Y, &
     t_lo_mev,t_hi_mev,e_lo,e_hi,rfeps, xtemp)

  use composition, only: nspec
  implicit none

  integer, intent(in)    :: i
  real*8,  intent(in)    :: xrho, xye, xenr, abar, zbar, e_offset
  real*8,  intent(in)    :: Y(nspec), t_lo_mev, t_hi_mev, e_lo, e_hi, rfeps
  real*8,  intent(inout) :: xtemp

  real*8  :: a, b, fa, fb, f, t, t_new, e, dfdt, t_prev, f_prev
  real*8  :: p, ent, cs2, dedt, dpde, dpdr
  logical :: bracketed, have_prev
  integer :: it, keyerr
  integer, parameter :: maxit = 100

  a = t_lo_mev ; fa = e_lo - xenr
  b = t_hi_mev ; fb = e_hi - xenr

  bracketed = (fa*fb .le. 0.0d0)
  have_prev = .false.
  t_prev    = 0.0d0
  f_prev    = 0.0d0

  t = min(max(xtemp, a), b)

  do it = 1, maxit

     call blend_state(xrho,t,xye,abar,zbar,e_offset,rfeps, &
          e,p,ent,cs2,dedt,dpde,dpdr,dfdt,keyerr)
     f = e - xenr

     if (abs(f) .le. max(rfeps*abs(xenr), 8.0d0*spacing(abs(e)))) then
        xtemp = t
        return
     end if

     ! a sign change between successive iterates is a bracket too
     if (.not. bracketed .and. have_prev) then
        if (sign(1.0d0,f) .ne. sign(1.0d0,f_prev)) then
           bracketed = .true.
           if (t .lt. t_prev) then
              a = t      ; fa = f
              b = t_prev ; fb = f_prev
           else
              a = t_prev ; fa = f_prev
              b = t      ; fb = f
           end if
        end if
     end if

     if (bracketed) then

        if (sign(1.0d0,f) .eq. sign(1.0d0,fa)) then
           a = t ; fa = f
        else
           b = t ; fb = f
        end if

        if (b - a .le. 1.0d-10*t) then
           xtemp = t
           return
        end if

        ! Newton from the analytic derivative (which includes the dw/dT term);
        ! the bracket makes its accuracy near the seam kinks non-critical.
        if (dfdt .eq. 0.0d0) then
           t_new = 0.5d0*(a + b)
        else
           t_new = t - f/dfdt
           if (t_new .le. a .or. t_new .ge. b) t_new = 0.5d0*(a + b)
        end if

     else

        t_prev = t ; f_prev = f ; have_prev = .true.

        if (dfdt .eq. 0.0d0) then
           t_new = 0.5d0*(a + b)
        else
           t_new = t - f/dfdt
           t_new = max(0.5d0*t, min(t_new, 2.0d0*t))   ! factor-of-2 damping
        end if
        t_new = min(max(t_new, a), b)

        ! converged in T -- but only if the iterate is free to move.  Pinned
        ! against a window edge the step is zero for the wrong reason, and the
        ! root is not in this regime at all; let that run out and STOP.
        if (abs(t_new - t) .le. 1.0d-10*t .and. &
            t_new .gt. a .and. t_new .lt. b) then
           xtemp = t_new
           return
        end if

     end if

     t = t_new

  end do

  write(*,*) 'blend_invert: window root find did not converge'
  call backtrace
  write(*,*) '  guess T    = ', xtemp, ' MeV'
  write(*,*) '  last T     = ', t, ' MeV'
  write(*,*) '  last f     = ', f
  write(*,*) '  bracketed  = ', bracketed
  write(*,*) '  bracket    = ', a, b
  call blend_diagnostics(i,xrho,xye,xenr,abar,zbar,e_offset,e_lo,e_hi,Y)
  STOP 'nuc_helm_eos_short: blended keytemp=0 inversion failed'

end subroutine blend_invert

! ---------------------------------------------------------------------------
! Zone diagnostics shared by the blend-window failure paths.
! ---------------------------------------------------------------------------
subroutine blend_diagnostics(i,xrho,xye,xenr,abar,zbar,e_offset,e_lo,e_hi,Y)

  use composition, only: nspec
  implicit none

  integer, intent(in) :: i
  real*8,  intent(in) :: xrho, xye, xenr, abar, zbar, e_offset, e_lo, e_hi
  real*8,  intent(in) :: Y(nspec)

  write(*,*) '  zone i     = ', i
  write(*,*) '  rho        = ', xrho
  write(*,*) '  ye         = ', xye
  write(*,*) '  target eps = ', xenr
  write(*,*) '  e_lo, e_hi = ', e_lo, e_hi
  write(*,*) '  abar, zbar = ', abar, zbar
  write(*,*) '  e_offset   = ', e_offset
  write(*,*) '  Yion       = ', Y

end subroutine blend_diagnostics

! ---------------------------------------------------------------------------
! The blended thermodynamic state at a KNOWN temperature t_mev.  One weight
! w(T) is applied to every returned quantity, so the state stays a single
! consistent thermodynamic point.  Only the backend(s) selected by w are
! evaluated, so cold zones never touch nuc_eos and hot zones never touch
! Helmholtz.
!
! This is the single evaluator used both by the keytemp=0 root find
! (blend_invert) and to fill the outputs, so the returned state is exactly the
! converged iterate.
!
! dedt is the blended heat capacity in the nuc_eos convention (per MeV; the
! Helmholtz cv is per Kelvin and is converted here).  dfdt additionally carries
! the dw/dT term, i.e. it is the full derivative of the blended energy that the
! root find needs, and differs from dedt inside the window.
! ---------------------------------------------------------------------------
subroutine blend_state(xrho,t_mev,xye,abar,zbar,e_offset,rfeps, &
     e,p,ent,cs2,dedt,dpderho,dpdrhoe,dfdt,keyerr)

  use GR1D_module, only: T_eos_high, T_eos_low, temp_mev_to_kelvin
  implicit none

  real*8,  intent(in)  :: xrho, t_mev, xye, abar, zbar, e_offset, rfeps
  real*8,  intent(out) :: e, p, ent, cs2, dedt, dpderho, dpdrhoe, dfdt
  integer, intent(out) :: keyerr

  real*8 :: w, tk, t_lo_mev, t_hi_mev
  real*8 :: t_n, e_n, p_n, ent_n, cs2_n, dedt_n, dpde_n, dpdr_n, munu_n
  real*8 :: tk_h, e_h, p_h, ent_h, cs2_h, dedt_h, dpde_h, dpdr_h

  keyerr = 0
  tk = t_mev * temp_mev_to_kelvin
  call blend_weight(tk, w)

  if (w .gt. 0.0d0) then
     t_n = t_mev
     call nuc_eos_short(xrho,t_n,xye,e_n,p_n,ent_n,cs2_n,dedt_n, &
          dpde_n,dpdr_n,munu_n,1,keyerr,rfeps)
  end if
  if (w .lt. 1.0d0) then
     tk_h = tk
     e_h  = 0.0d0
     call eval_helm(xrho,tk_h,abar,zbar,e_offset,e_h, &
          1,p_h,ent_h,cs2_h,dedt_h,dpde_h,dpdr_h)
     dedt_h = dedt_h * temp_mev_to_kelvin   ! cv per Kelvin -> per MeV
  end if

  if (w .le. 0.0d0) then
     e       = e_h
     p       = p_h
     ent     = ent_h
     cs2     = cs2_h
     dedt    = dedt_h
     dpderho = dpde_h
     dpdrhoe = dpdr_h
     dfdt    = dedt
  else if (w .ge. 1.0d0) then
     e       = e_n
     p       = p_n
     ent     = ent_n
     cs2     = cs2_n
     dedt    = dedt_n
     dpderho = dpde_n
     dpdrhoe = dpdr_n
     dfdt    = dedt
  else
     e       = w*e_n    + (1.0d0-w)*e_h
     p       = w*p_n    + (1.0d0-w)*p_h
     ent     = w*ent_n  + (1.0d0-w)*ent_h
     cs2     = w*cs2_n  + (1.0d0-w)*cs2_h
     dedt    = w*dedt_n + (1.0d0-w)*dedt_h
     dpderho = w*dpde_n + (1.0d0-w)*dpde_h
     dpdrhoe = w*dpdr_n + (1.0d0-w)*dpdr_h
     ! d/dT of w*e_n + (1-w)*e_h with the linear weight: the extra term is
     ! (e_n - e_h)*dw/dT and dw/dT = 1/(T_hi - T_lo) inside the window
     t_lo_mev = T_eos_low  / temp_mev_to_kelvin
     t_hi_mev = T_eos_high / temp_mev_to_kelvin
     dfdt    = dedt + (e_n - e_h)/(t_hi_mev - t_lo_mev)
  end if

end subroutine blend_state

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
