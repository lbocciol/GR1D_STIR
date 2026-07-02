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
! Composite finite-T nuclear EOS  <->  Helmholtz EOS, blended by temperature.
!
! Drop-in replacement for nuc_eos_short with ONE extra leading argument: the
! radial zone index i (used only to fetch the local network composition
! Yion(:,i)).  Every other argument and unit (CGS in/out) matches
! nuc_eos_short exactly, so the caller is unaware of the blend.
!
! Regimes (T = known temperature for keytemp=1, or the running estimate for
! keytemp=0):
!   T >= T_NSE      -> nuc_eos only              (w = 1)
!   T <= T_interp   -> Helmholtz only            (w = 0)
!   T_interp<T<T_NSE-> blend BOTH with ONE smoothstep weight w(T) in [0,1],
!                      applied identically to every returned quantity so the
!                      state stays a single consistent thermodynamic point.
!
! Robust fallbacks (so a burning build is well-defined in EVERY zone):
!   * no composition loaded here (sum Y ~ 0) -> pure nuc_eos
!   * entropy-mode root find (keytemp=2)     -> pure nuc_eos
!     (the Helmholtz path has no entropy inversion wired up)
! ===========================================================================
subroutine nuc_helm_eos_short(i,xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
     xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)

  use GR1D_module, only: Yion, T_NSE, T_interp, temp_mev_to_kelvin
  use composition, only: nspec, zion
  implicit none

  integer, intent(in)    :: i, keytemp
  real*8,  intent(in)    :: xrho, xye, rfeps
  real*8,  intent(inout) :: xtemp, xenr
  real*8,  intent(out)   :: xprs, xent, xcs2, xdedt
  real*8,  intent(out)   :: xdpderho, xdpdrhoe, xmunu
  integer, intent(out)   :: keyerr

  real*8 :: Y(nspec), ytot, abar, zbar, e_offset, w, s, tk, tk_h
  ! nuc_eos branch outputs (local copies so xtemp/xenr aren't clobbered)
  real*8 :: t_n, e_n, p_n, ent_n, cs2_n, dedt_n, dpde_n, dpdr_n, munu_n
  ! Helmholtz branch outputs
  real*8 :: e_h, p_h, ent_h, cs2_h, dedt_h, dpde_h, dpdr_h
  integer :: mode

  keyerr = 0
  xmunu  = 0.0d0

  ! entropy-mode inversion (keytemp=2) is wired up only for nuc_eos
  if (keytemp .eq. 2) then
     call nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
          xdpderho,xdpdrhoe,xmunu,2,keyerr,rfeps)
     xmunu = 0.0d0
     return
  end if

  Y    = Yion(:,i)
  ytot = sum(Y)

  abar = 1.0d0 / ytot
  zbar = abar * sum(zion * Y)
  ! call get_energy_offset(xtemp, xye, e_offset)
  call get_energy_offset(Y, xye, e_offset)

  ! GR1D carries temperature in MeV (nuc_eos convention); the NSE thresholds
  ! T_NSE/T_interp and the Helmholtz table are in Kelvin -> convert for the
  ! weight and for the Helmholtz calls; nuc_eos is always called in MeV.
  tk = xtemp * temp_mev_to_kelvin

  ! =========================================================================
  ! keytemp = 0 (eps known): invert the BLENDED energy
  !     e_blend(T) = w(T)*e_nuc(rho,T) + (1-w(T))*e_helm(rho,T)
  ! for ONE temperature, so the returned state is a single consistent point
  ! (rather than inverting each backend to its own T and blending across two
  ! different thermodynamic points).  The regime is decided by energy in the
  ! transition window; outside it, only the locally-valid backend is touched
  ! (nuc_eos is not valid in the cold Helmholtz regime, and vice versa).
  ! =========================================================================
  if (keytemp .eq. 0) then
     call blend_weight(tk, w)

     ! pure Helmholtz (cold): invert helm only
     if (w .le. 0.0d0) then
        tk_h = tk
        call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,xenr, &
             0,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe)
        if (tk_h .le. T_interp) then
           xtemp = tk_h / temp_mev_to_kelvin
           xmunu = 0.0d0
           return
        end if
        ! material heated across the seam -> solve consistently below
     end if

     ! pure nuc_eos (hot): invert nuc only
     if (w .ge. 1.0d0) then
        t_n = xtemp
        call nuc_eos_short(xrho,t_n,xye,xenr,xprs,xent,xcs2,xdedt, &
             xdpderho,xdpdrhoe,xmunu,0,keyerr,rfeps)
        if (t_n*temp_mev_to_kelvin .ge. T_NSE) then
           xtemp = t_n
           xmunu = 0.0d0
           return
        end if
        ! material cooled across the seam -> solve consistently below
     end if

     ! transition (both backends valid): root-find the single blended T
     call blend_invert_temp(xrho,xye,abar,zbar,e_offset,xenr,rfeps, &
          xtemp,w,mode,keyerr)
     tk = xtemp * temp_mev_to_kelvin

     if (mode .eq. 0) then
        tk_h = tk
        call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,xenr, &
             1,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe)
        xmunu = 0.0d0
        return
     else if (mode .eq. 1) then
        t_n = xtemp
        call nuc_eos_short(xrho,t_n,xye,xenr,xprs,xent,xcs2,xdedt, &
             xdpderho,xdpdrhoe,xmunu,1,keyerr,rfeps)
        xmunu = 0.0d0
        return
     end if

     ! blended regime: evaluate BOTH backends at the single solved T
     t_n = xtemp
     call nuc_eos_short(xrho,t_n,xye,e_n,p_n,ent_n,cs2_n,dedt_n, &
          dpde_n,dpdr_n,munu_n,1,keyerr,rfeps)
     tk_h = tk
     call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,e_h, &
          1,p_h,ent_h,cs2_h,dedt_h,dpde_h,dpdr_h)

     xprs     = w*p_n    + (1.0d0-w)*p_h
     xent     = w*ent_n  + (1.0d0-w)*ent_h
     xcs2     = w*cs2_n  + (1.0d0-w)*cs2_h
     xdedt    = w*dedt_n + (1.0d0-w)*dedt_h
     xdpderho = w*dpde_n + (1.0d0-w)*dpde_h
     xdpdrhoe = w*dpdr_n + (1.0d0-w)*dpdr_h
     xmunu    = 0.0d0
     return
  end if

  ! ---- temperature blend weight (smoothstep in T, Kelvin) -----------------
  ! keytemp = 1 (T known) or 2 (entropy): regime selected from the known T.
  if (tk .ge. T_NSE) then
     w = 1.0d0
  else if (tk .le. T_interp) then
     w = 0.0d0
  else
     s = (tk - T_interp) / (T_NSE - T_interp)
     w = s*s*(3.0d0 - 2.0d0*s)
  end if

  ! ---- pure Helmholtz -----------------------------------------------------
  if (w .le. 0.0d0) then
     tk_h = tk
     call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,xenr, &
          keytemp,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe)
     if (keytemp .eq. 0) xtemp = tk_h / temp_mev_to_kelvin

     return
  end if

  ! ---- pure nuc_eos -------------------------------------------------------
  if (w .ge. 1.0d0) then
     call nuc_eos_short(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
          xdpderho,xdpdrhoe,xmunu,keytemp,keyerr,rfeps)
     xmunu = 0.0d0
     return
  end if

  ! ---- transition: evaluate BOTH backends at the known T, blend with w ----
  ! Reached only for keytemp = 1 (T known) now; keytemp = 0 is fully handled
  ! by the blended-energy inversion above.
  t_n = xtemp
  call nuc_eos_short(xrho,t_n,xye,e_n,p_n,ent_n,cs2_n,dedt_n, &
       dpde_n,dpdr_n,munu_n,1,keyerr,rfeps)
  tk_h = tk
  call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,e_h, &
       1,p_h,ent_h,cs2_h,dedt_h,dpde_h,dpdr_h)
  xenr = w*e_n + (1.0d0-w)*e_h

  xprs     = w*p_n    + (1.0d0-w)*p_h
  xent     = w*ent_n  + (1.0d0-w)*ent_h
  xcs2     = w*cs2_n  + (1.0d0-w)*cs2_h
  xdedt    = w*dedt_n + (1.0d0-w)*dedt_h
  xdpderho = w*dpde_n + (1.0d0-w)*dpde_h
  xdpdrhoe = w*dpdr_n + (1.0d0-w)*dpdr_h
  xmunu    = 0.0d0

end subroutine nuc_helm_eos_short

! ===========================================================================
! Composite full-vector EOS: mirror of nuc_eos_full (adds the NSE composition
! mass fractions xa/xh/xn/xp, mean nuclei abar/zbar, and chemical potentials
! mu_e/mu_n/mu_p/muhat) on top of the same temperature blend as
! nuc_helm_eos_short.  Drop-in replacement for nuc_eos_full with the extra
! leading zone index i.  This is what eos_full uses for eoskey==3.
!
! In the Helmholtz (network) regime the NSE composition does not apply, so:
!   xa = xh = xn = xp = 0,  abar/zbar from the network composition,
!   mu_e = mu_n = mu_p = muhat = 0  (matching the non-burning short-EOS path).
! Every quantity is blended with the SAME weight w as the thermodynamics.
! ===========================================================================
subroutine nuc_helm_eos_full(i,xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
     xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar, &
     xmu_e,xmu_n,xmu_p,xmuhat,keytemp,keyerr,rfeps)

  use GR1D_module, only: Yion, T_NSE, T_interp, temp_mev_to_kelvin
  use composition, only: nspec, zion
  implicit none

  integer, intent(in)    :: i, keytemp
  real*8,  intent(in)    :: xrho, xye, rfeps
  real*8,  intent(inout) :: xtemp, xenr
  real*8,  intent(out)   :: xprs, xent, xcs2, xdedt
  real*8,  intent(out)   :: xdpderho, xdpdrhoe
  real*8,  intent(out)   :: xxa, xxh, xxn, xxp, xabar, xzbar
  real*8,  intent(out)   :: xmu_e, xmu_n, xmu_p, xmuhat
  integer, intent(out)   :: keyerr

  real*8 :: Y(nspec), ytot, abar, zbar, e_offset, w, s, tk, tk_h
  ! nuc_eos branch outputs (local copies so xtemp/xenr aren't clobbered)
  real*8 :: t_n, e_n, p_n, ent_n, cs2_n, dedt_n, dpde_n, dpdr_n
  real*8 :: a_n, h_n, xn_n, xp_n, abar_n, zbar_n, mue_n, mun_n, mup_n, muhat_n
  ! Helmholtz branch outputs
  real*8 :: e_h, p_h, ent_h, cs2_h, dedt_h, dpde_h, dpdr_h
  integer :: mode

  keyerr = 0

  ! entropy-mode inversion (keytemp=2) is wired up only for nuc_eos
  if (keytemp .eq. 2) then
     call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar, &
          xmu_e,xmu_n,xmu_p,xmuhat,2,keyerr,rfeps)
     return
  end if

  Y    = Yion(:,i)
  ytot = sum(Y)

  abar = 1.0d0 / ytot
  zbar = abar * sum(zion * Y)
  ! call get_energy_offset(xtemp, xye, e_offset)
  call get_energy_offset(Y, xye, e_offset)

  tk = xtemp * temp_mev_to_kelvin

  ! =========================================================================
  ! keytemp = 0 (eps known): invert the BLENDED energy for ONE consistent T
  ! (see nuc_helm_eos_short for the rationale).  Regime is decided by energy
  ! in the transition window; outside it only the locally-valid backend runs.
  ! =========================================================================
  if (keytemp .eq. 0) then
     call blend_weight(tk, w)

     ! pure Helmholtz (cold): invert helm only, no NSE composition
     if (w .le. 0.0d0) then
        tk_h = tk
        call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,xenr, &
             0,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe)
        if (tk_h .le. T_interp) then
           xtemp = tk_h / temp_mev_to_kelvin
           xxa = 0.0d0; xxh = 0.0d0; xxn = 0.0d0; xxp = 0.0d0
           xabar = abar; xzbar = zbar
           xmu_e = 0.0d0; xmu_n = 0.0d0; xmu_p = 0.0d0; xmuhat = 0.0d0
           return
        end if
     end if

     ! pure nuc_eos (hot): invert nuc only
     if (w .ge. 1.0d0) then
        t_n = xtemp
        call nuc_eos_full(xrho,t_n,xye,xenr,xprs,xent,xcs2,xdedt, &
             xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar, &
             xmu_e,xmu_n,xmu_p,xmuhat,0,keyerr,rfeps)
        if (t_n*temp_mev_to_kelvin .ge. T_NSE) then
           xtemp = t_n
           return
        end if
     end if

     ! transition (both backends valid): root-find the single blended T
     call blend_invert_temp(xrho,xye,abar,zbar,e_offset,xenr,rfeps, &
          xtemp,w,mode,keyerr)
     tk = xtemp * temp_mev_to_kelvin

     if (mode .eq. 0) then
        tk_h = tk
        call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,xenr, &
             1,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe)
        xxa = 0.0d0; xxh = 0.0d0; xxn = 0.0d0; xxp = 0.0d0
        xabar = abar; xzbar = zbar
        xmu_e = 0.0d0; xmu_n = 0.0d0; xmu_p = 0.0d0; xmuhat = 0.0d0
        return
     else if (mode .eq. 1) then
        call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
             xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar, &
             xmu_e,xmu_n,xmu_p,xmuhat,1,keyerr,rfeps)
        return
     end if

     ! blended regime: evaluate BOTH backends at the single solved T
     t_n = xtemp
     call nuc_eos_full(xrho,t_n,xye,e_n,p_n,ent_n,cs2_n,dedt_n, &
          dpde_n,dpdr_n,a_n,h_n,xn_n,xp_n,abar_n,zbar_n, &
          mue_n,mun_n,mup_n,muhat_n,1,keyerr,rfeps)
     tk_h = tk
     call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,e_h, &
          1,p_h,ent_h,cs2_h,dedt_h,dpde_h,dpdr_h)

     xprs     = w*p_n    + (1.0d0-w)*p_h
     xent     = w*ent_n  + (1.0d0-w)*ent_h
     xcs2     = w*cs2_n  + (1.0d0-w)*cs2_h
     xdedt    = w*dedt_n + (1.0d0-w)*dedt_h
     xdpderho = w*dpde_n + (1.0d0-w)*dpde_h
     xdpdrhoe = w*dpdr_n + (1.0d0-w)*dpdr_h
     xxa   = w*a_n
     xxh   = w*h_n
     xxn   = w*xn_n
     xxp   = w*xp_n
     xabar = w*abar_n + (1.0d0-w)*abar
     xzbar = w*zbar_n + (1.0d0-w)*zbar
     xmu_e  = w*mue_n
     xmu_n  = w*mun_n
     xmu_p  = w*mup_n
     xmuhat = w*muhat_n
     return
  end if

  ! ---- temperature blend weight (smoothstep in T, Kelvin) -----------------
  ! keytemp = 1 (T known): regime selected from the known temperature.
  if (tk .ge. T_NSE) then
     w = 1.0d0
  else if (tk .le. T_interp) then
     w = 0.0d0
  else
     s = (tk - T_interp) / (T_NSE - T_interp)
     w = s*s*(3.0d0 - 2.0d0*s)
  end if

  ! ---- pure Helmholtz -----------------------------------------------------
  if (w .le. 0.0d0) then
     tk_h = tk
     call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,xenr, &
          keytemp,xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe)
     if (keytemp .eq. 0) xtemp = tk_h / temp_mev_to_kelvin
     xxa = 0.0d0; xxh = 0.0d0; xxn = 0.0d0; xxp = 0.0d0
     xabar = abar; xzbar = zbar
     xmu_e = 0.0d0; xmu_n = 0.0d0; xmu_p = 0.0d0; xmuhat = 0.0d0
     return
  end if

  ! ---- pure nuc_eos -------------------------------------------------------
  if (w .ge. 1.0d0) then
     call nuc_eos_full(xrho,xtemp,xye,xenr,xprs,xent,xcs2,xdedt, &
          xdpderho,xdpdrhoe,xxa,xxh,xxn,xxp,xabar,xzbar, &
          xmu_e,xmu_n,xmu_p,xmuhat,keytemp,keyerr,rfeps)
     return
  end if

  ! ---- transition: evaluate BOTH backends at the known T, blend with w ----
  ! Reached only for keytemp = 1 now; keytemp = 0 is handled by the blended-
  ! energy inversion above.
  t_n = xtemp
  call nuc_eos_full(xrho,t_n,xye,e_n,p_n,ent_n,cs2_n,dedt_n, &
       dpde_n,dpdr_n,a_n,h_n,xn_n,xp_n,abar_n,zbar_n, &
       mue_n,mun_n,mup_n,muhat_n,1,keyerr,rfeps)
  tk_h = tk
  call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,e_h, &
       1,p_h,ent_h,cs2_h,dedt_h,dpde_h,dpdr_h)
  xenr = w*e_n + (1.0d0-w)*e_h

  ! thermodynamics
  xprs     = w*p_n    + (1.0d0-w)*p_h
  xent     = w*ent_n  + (1.0d0-w)*ent_h
  xcs2     = w*cs2_n  + (1.0d0-w)*cs2_h
  xdedt    = w*dedt_n + (1.0d0-w)*dedt_h
  xdpderho = w*dpde_n + (1.0d0-w)*dpde_h
  xdpdrhoe = w*dpdr_n + (1.0d0-w)*dpdr_h
  ! composition: Helmholtz side carries no NSE mass fractions (=0); abar/zbar
  ! cross over from nuc_eos to the network values; mu's cross over to 0.
  xxa   = w*a_n
  xxh   = w*h_n
  xxn   = w*xn_n
  xxp   = w*xp_n
  xabar = w*abar_n + (1.0d0-w)*abar
  xzbar = w*zbar_n + (1.0d0-w)*zbar
  xmu_e  = w*mue_n
  xmu_n  = w*mun_n
  xmu_p  = w*mup_n
  xmuhat = w*muhat_n

end subroutine nuc_helm_eos_full

! ---------------------------------------------------------------------------
! Single Helmholtz evaluation in the SAME CGS in/out convention as
! nuc_eos_short.  The energy offset [erg/g] (from get_energy_offset: the matched
! e_nuc - e_helm table, build_energy_offset_OttEOS) is supplied as the state's e_offset so
! the Helmholtz energy zero-point matches the finite-T nuclear EOS (see
! eps-energy convention / src/CLAUDE.md).
!   keytemp = 1 : temp known -> enr returned
!   keytemp = 0 : enr known  -> temp solved
! cs2 is returned as the NEWTONIAN gam1*p/rho [cm^2/s^2]; the relativistic
! correction is applied uniformly by the caller, exactly as for nuc_eos.
! ---------------------------------------------------------------------------
subroutine eval_helm(rho,temp,ye,abar,zbar,e_offset,enr,keytemp, &
     prs,ent,cs2,dedt,dpderho,dpdrhoe)

  use wlHelmholtzEOS, only: HelmholtzStateType, HelmEOS, &
                            eos_input_rt, eos_input_re
  implicit none

  real*8,  intent(in)    :: rho, ye, abar, zbar, e_offset
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
     enr = st % e
  else
     st % e = enr                 ! target energy; e_offset handled in HelmEOS
     call HelmEOS(eos_input_re, st)
     temp = st % T
  end if

  prs     = st % p
  ent     = st % s
  cs2     = st % gam1 * st % p / st % rho   ! Newtonian cs^2 [cm^2/s^2]
  dedt    = st % cv
  dpderho = st % dpde
  dpdrhoe = st % dpdr_e

end subroutine eval_helm

! ---------------------------------------------------------------------------
! Smoothstep temperature blend weight w(T) in [0,1], T in Kelvin.
!   T >= T_NSE    -> 1 (pure nuc_eos)
!   T <= T_interp -> 0 (pure Helmholtz)
!   between       -> smoothstep s^2(3-2s)
! ---------------------------------------------------------------------------
subroutine blend_weight(tk, w)

  use GR1D_module, only: T_NSE, T_interp
  implicit none

  real*8, intent(in)  :: tk
  real*8, intent(out) :: w
  real*8 :: s

  if (tk .ge. T_NSE) then
     w = 1.0d0
  else if (tk .le. T_interp) then
     w = 0.0d0
  else
     s = (tk - T_interp) / (T_NSE - T_interp)
     w = s*s*(3.0d0 - 2.0d0*s)
  end if

end subroutine blend_weight

! ---------------------------------------------------------------------------
! Invert the BLENDED specific energy
!     e_blend(T) = w(T)*e_nuc(rho,T) + (1-w(T))*e_helm(rho,T)
! for the single temperature T [MeV] at fixed (rho,Ye) matching target_eps
! [erg/g].  Both backends are evaluated in T-input mode, so the result is one
! consistent thermodynamic point (unlike inverting each backend separately and
! blending across two different temperatures).
!
! The transition window [T_interp, T_NSE] gives a guaranteed bracket:
!     e_blend(T_interp) = e_helm(T_interp)  (w=0)
!     e_blend(T_NSE)    = e_nuc(T_NSE)      (w=1)
! If target_eps lies below/above this bracket the solution is in a pure regime
! and the single relevant backend is inverted directly (mode 0/1); otherwise a
! bisection inside the window finds the blended T (mode 2).  Bisection is used
! for robustness: f(T)=e_blend(T)-target_eps has a guaranteed sign change on the
! bracket regardless of any non-monotonicity from the T-dependent weight.
!
! Returns t_mev [MeV], the weight w_out at t_mev, and mode (0/1/2).
! Callers in the transition window only -> both backends are valid here.
! ---------------------------------------------------------------------------
subroutine blend_invert_temp(xrho,xye,abar,zbar,e_offset,target_eps,rfeps, &
     t_mev,w_out,mode,keyerr)

  use GR1D_module, only: T_NSE, T_interp, temp_mev_to_kelvin
  implicit none

  real*8,  intent(in)  :: xrho, xye, abar, zbar, e_offset, target_eps, rfeps
  real*8,  intent(out) :: t_mev, w_out
  integer, intent(out) :: mode, keyerr

  real*8 :: t_lo, t_hi, e_lo, e_hi
  real*8 :: t_a, t_b, t_m, f_m, e_blend, tk, w
  real*8 :: e_n, e_h, t_nuc, tk_h
  ! discarded backend outputs
  real*8 :: p_d, ent_d, cs2_d, dedt_d, dpde_d, dpdr_d, munu_d
  integer :: it
  integer, parameter :: maxit = 100
  real*8,  parameter :: tol_floor = 1.0d-12

  keyerr = 0

  t_lo = T_interp / temp_mev_to_kelvin   ! MeV
  t_hi = T_NSE    / temp_mev_to_kelvin   ! MeV

  ! lower bracket energy: pure Helmholtz at T_interp (w = 0)
  tk_h = T_interp
  call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,e_lo, &
       1,p_d,ent_d,cs2_d,dedt_d,dpde_d,dpdr_d)

  ! upper bracket energy: pure nuc_eos at T_NSE (w = 1)
  t_nuc = t_hi
  call nuc_eos_short(xrho,t_nuc,xye,e_hi,p_d,ent_d,cs2_d,dedt_d, &
       dpde_d,dpdr_d,munu_d,1,keyerr,rfeps)

  if (target_eps .le. e_lo) then
     ! pure Helmholtz regime: invert helm for T
     tk_h = T_interp                 ! initial guess (Kelvin)
     e_h  = target_eps
     call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,e_h, &
          0,p_d,ent_d,cs2_d,dedt_d,dpde_d,dpdr_d)
     t_mev = tk_h / temp_mev_to_kelvin
     w_out = 0.0d0
     mode  = 0
     return
  else if (target_eps .ge. e_hi) then
     ! pure nuc_eos regime: invert nuc for T
     t_nuc = t_hi                    ! initial guess (MeV)
     e_n   = target_eps
     call nuc_eos_short(xrho,t_nuc,xye,e_n,p_d,ent_d,cs2_d,dedt_d, &
          dpde_d,dpdr_d,munu_d,0,keyerr,rfeps)
     t_mev = t_nuc
     w_out = 1.0d0
     mode  = 1
     return
  end if

  ! blended regime: bisection on [t_lo,t_hi].
  ! f(t_lo) = e_lo - target < 0,  f(t_hi) = e_hi - target > 0.
  t_a = t_lo
  t_b = t_hi
  do it = 1, maxit
     t_m = 0.5d0*(t_a + t_b)
     tk  = t_m * temp_mev_to_kelvin
     call blend_weight(tk, w)
     t_nuc = t_m
     call nuc_eos_short(xrho,t_nuc,xye,e_n,p_d,ent_d,cs2_d,dedt_d, &
          dpde_d,dpdr_d,munu_d,1,keyerr,rfeps)
     tk_h = tk
     call eval_helm(xrho,tk_h,xye,abar,zbar,e_offset,e_h, &
          1,p_d,ent_d,cs2_d,dedt_d,dpde_d,dpdr_d)
     e_blend = w*e_n + (1.0d0-w)*e_h
     f_m = e_blend - target_eps
     if (f_m .gt. 0.0d0) then
        t_b = t_m
     else
        t_a = t_m
     end if
     if (abs(t_b - t_a) .le. max(rfeps,tol_floor)*t_m) exit
  end do

  t_mev = 0.5d0*(t_a + t_b)
  call blend_weight(t_mev*temp_mev_to_kelvin, w_out)
  mode = 2

end subroutine blend_invert_temp

! ---------------------------------------------------------------------------
! Build the composite-EOS energy-offset table by MATCHING the two backends at a
! single transition density eos_offset_rho, for every (T,Ye) on the nuc_eos grid
! (only T <= T_NSE, where Helmholtz contributes).  This is the GR1D analog of the
! stellarcollapse `loweos` extension:
!
!     energy_offset(T,Ye) = e_nuc(rho_tr,T,Ye) - e_helm(rho_tr,T,Ye, NSE comp) ,
!
! where the Helmholtz call uses nuc_eos's OWN (NSE) abar/zbar so both sides
! describe the same matter -- the offset is then exactly the nuclear zero-point
! gap, and e_helm + offset meets e_nuc at the seam by construction.  Called once
! at startup (start.F90), after both EOS tables are loaded.
subroutine build_energy_offset_OttEOS()

  use GR1D_module, only: energy_offset_tab, eos_offset_rho, T_NSE, temp_mev_to_kelvin
  use eosmodule,   only: ntemp, nye, logtemp, eos_ye => ye
  implicit none

  integer :: it, iy, keyerr, nfail
  real*8  :: T_mev, tk, ye_v, rfeps
  real*8  :: e_nuc, e_helm
  real*8  :: p, s, cs2, dedt, dpde, dpdr                 ! discarded outputs
  real*8  :: xa, xh, xn, xp, abar, zbar, mue, mun, mup, muhat

  rfeps = 1.0d-9
  if (.not. allocated(energy_offset_tab)) allocate(energy_offset_tab(ntemp,nye))
  energy_offset_tab = 0.0d0
  nfail = 0

  do iy = 1, nye
     ye_v = eos_ye(iy)
     do it = 1, ntemp
        T_mev = 10.0d0**logtemp(it)
        tk    = T_mev * temp_mev_to_kelvin
        if (tk .gt. T_NSE) cycle           ! offset only used where Helmholtz contributes

        keyerr = 0
        call nuc_eos_full(eos_offset_rho, T_mev, ye_v, e_nuc, p,s,cs2,dedt, &
             dpde,dpdr, xa,xh,xn,xp, abar,zbar, mue,mun,mup,muhat, 1, keyerr, rfeps)
        if (keyerr .ne. 0) then
           nfail = nfail + 1
           cycle
        end if

        ! pure Helmholtz (e_offset = 0) with nuc_eos's NSE composition
        call eval_helm(eos_offset_rho, tk, ye_v, 28.0d0, 14.0d0, 0.0d0, e_helm, 1, &
             p, s, cs2, dedt, dpde, dpdr)

        energy_offset_tab(it,iy) = e_nuc - e_helm
     end do
  end do

  write(*,*) "build_energy_offset_OttEOS: matched at rho_tr =", eos_offset_rho, " g/cc"
  write(*,*) "  nuc_eos failures =", nfail, " of", ntemp*nye, " grid points"
  write(*,*) "  energy_offset range [erg/g]:", minval(energy_offset_tab), &
                                               maxval(energy_offset_tab)

end subroutine build_energy_offset_OttEOS

subroutine get_energy_offset(Y, ye, e_offset)

  use GR1D_module, only: mev_to_erg, avo
  use composition, only: nspec, nuclei_binding_energy, nuclei_mass_excess
  implicit none

  real*8, intent(in) :: Y(nspec), ye
  real*8, intent(out) :: e_offset
  
  real*8, parameter  :: Qnp = 1.293333d0

  ! Correctly multiply the slices element-by-element before summing
  ! Notice that nuclei_mass_excess contains nucleons
  e_offset = - mev_to_erg * avo * sum(nuclei_mass_excess * Y)
    
end subroutine get_energy_offset

! ---------------------------------------------------------------------------
! Bilinear lookup of the energy-offset table at (T [MeV], Ye), on the nuc_eos
! log10(T)/Ye grid (clamped to range).  Returns e_offset [erg/g].
! subroutine get_energy_offset(T_mev, ye_in, e_offset)

!   use GR1D_module, only: energy_offset_tab
!   use eosmodule,   only: ntemp, nye, logtemp, eos_ye => ye
!   implicit none

!   real*8, intent(in)  :: T_mev, ye_in
!   real*8, intent(out) :: e_offset
!   real*8  :: lt, yy, dt, dy
!   integer :: it, iy

!   lt = log10(max(T_mev, 1.0d-30))
!   lt = min(max(lt, logtemp(1)), logtemp(ntemp))
!   yy = min(max(ye_in, eos_ye(1)), eos_ye(nye))

!   ! bracketing indices (grids monotonic increasing)
!   it = 1
!   do while (it .lt. ntemp-1 .and. logtemp(it+1) .lt. lt)
!      it = it + 1
!   end do
!   iy = 1
!   do while (iy .lt. nye-1 .and. eos_ye(iy+1) .lt. yy)
!      iy = iy + 1
!   end do

!   dt = (lt - logtemp(it)) / (logtemp(it+1) - logtemp(it))
!   dy = (yy - eos_ye(iy))  / (eos_ye(iy+1)  - eos_ye(iy))

!   e_offset = (1.0d0-dt)*(1.0d0-dy)*energy_offset_tab(it,  iy  ) &
!            +        dt *(1.0d0-dy)*energy_offset_tab(it+1,iy  ) &
!            + (1.0d0-dt)*       dy *energy_offset_tab(it,  iy+1) &
!            +        dt *       dy *energy_offset_tab(it+1,iy+1)

! end subroutine get_energy_offset
#endif
#endif
