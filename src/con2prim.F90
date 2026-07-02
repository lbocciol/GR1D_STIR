!-*-f90-*-
subroutine con2prim
  
  use GR1D_module
  use timers
  implicit none

  real*8 :: t1, t2

  CALL GetThisTime(t1)
  ! call con2prim_1
  call con2prim_grmhd
  CALL GetThisTime(t2)
  timer_c2p = timer_c2p + (t2 - t1)

end subroutine con2prim

!*************************************************************************

subroutine con2prim_1

  use GR1D_module
  use atmos
  use omp_lib
#ifdef HAVE_BURN
  use composition, only: nspec
#endif
  implicit none

#ifdef HAVE_BURN
  integer k
#endif
  real*8 tol, err, h, discrim
  real*8 low_tol
  real*8 pplus,pminus
  real*8 dpdrh, dpde,dedpress,drhodpress,temp1
  real*8 vv,rrho,eeps,pp,ww,op(n1),fp,dfdp
  ! *** for rotation
  real*8 vpv
  ! ***
  real*8 oeps(n1), old_press(n1)
  real*8 :: t1 = 0.0d0
  real*8 :: t2 = 0.0d0
  integer iminb,imaxb
  integer i,j,it
  integer :: max_iterations = 1000
  integer success,pt_counter

  ! dummies for EOS call
  real*8 eosdummy
  integer keyerr,keytemp
  ! end dummies for EOS call

  iminb = ghosts1+1
  imaxb = n1-ghosts1
  tol = 1.0d-10
  low_tol = tol
  err = 1.0d0
  success=0
  pt_counter = 0

  oeps = eps

  !test in shocktube or similar
  if (GR.and.(gravity_active.eqv..false.)) then
     do i=1,n1
        if (X(i).ne.1.0d0) then
           write(*,*) "X is not 1"
           stop
        endif
        if (alp(i).ne.1.0d0) then
           write(*,*) "alp is not 1"
           stop
        endif  
     enddo
  endif

  if (GR.and. .not. do_rotation) then
     do i=iminb,imaxb 

        err = 1.0d0
        low_tol = tol
        ye(i) = q(i,4)/q(i,1)
#ifdef HAVE_BURN
        do k=1,nspec
           Yion(k,i) = max(q(i,6+k)/q(i,1), 1.0d-50)
        enddo
#endif

        if (activate_turbulence) then
            v_turb(i) = sqrt(q(i,6)/q(i,1))
        endif

        if (q(i,1).eq.0.0d0) then
           v1(i) = 0.0d0
           v(i) = 0.0d0
           rho(i) = 0.0d0
           eps(i) = 0.0d0
	   press(i) = 0.0d0
        else
           ! atmosphere handling:
           if(rho(i).eq.atmo_rho) then
              q(i,1) = rho(i)
              q(i,2) = 0.0d0
              q(i,3) = rho(i)*eps(i)
              W(i) = 1.0d0
           endif
 	   if (q(i,2).eq.0.0d0) then
              v1(i) = 0.0d0
              v(i) = 0.0d0
              W(i) = 1.0d0
              rho(i) = q(i,1)/X(i)
              eps(i) = (q(i,3)+q(i,1)-q(i,1)/X(i))/rho(i)
              if (eps(i).lt.1.d-10) then
                 write(*,*) 'Help 1!!',i,eps(i),q(i,3),q(i,1),rho(i)/rho_gf,X(i)
                 eps(i)=1.d-10
                 stop
              endif
              keytemp = 0
              
              ! Limit eos variables, very drastic approach
              call ApplyEOS_limits_zone(i)
              call eos(i,rho(i),temp(i),ye(i),eps(i),pp, keytemp,keyerr,1, eoskey,eos_rf_prec)
              press(i) = pp
           else
              it = 0
              old_press(i) = press(i)
              do while (err.gt.tol.and.it.lt.max_iterations)
                 it = it + 1
                 op(i) = press(i)
                 !here vv is romero's v not v^r
                 vv = q(i,2)/(q(i,3)+q(i,1)+op(i))
                 if (vv.gt.1.0d0.or.vv.lt.-1.0d0) then
                    write(6,*) "We have a problem finding v:" 
                    write(6,*) "Timestep: ", nt
                    write(6,"(2i6,1P10E15.6)") it,i,rrho/rho_gf,rho(i)/rho_gf,eeps/eps_gf,ww
                    write(6,"(2i6,1P10E15.6)") it,i,q(i,2), (q(i,3)+q(i,1)+op(i)),vv,q(i,3), q(i,1), op(i)
                    call flush(6)
                    stop "con2prim problem"
                 endif
                 discrim = 1.0d0-vv**2
                 if (discrim.lt.0.0d0) then
                    write(*,*) "We have a problem with W", discrim, vv
                    stop
                 endif
                 ww = 1.0d0/sqrt(discrim)
                 rrho = q(i,1)/X(i)/ww
                 eeps = (q(i,3)+q(i,1)+op(i)*(1.0d0-ww**2))/(rrho*ww**2)-1.0d0

                 keytemp = 0
                 ! No need to limit EOS variables here
                 call eos_full(i,rrho,temp(i),ye(i),eeps,pp, & 
                      eosdummy,eosdummy,eosdummy,eosdummy,&
                      dpde,dpdrh, &
                      eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,&
                      eosdummy,eosdummy,eosdummy,eosdummy, &
                      keytemp,keyerr,eoskey,eos_rf_prec)
                 if(keyerr.ne.0) then
                    write(6,*) "#############################################"
                    write(6,"(i4,i6,1P10E15.6)") it,i,x1(i)/length_gf,&
                         rrho/rho_gf,rho(i)/rho_gf,temp(i),ye(i),q(i,1)
                    write(6,"(i4,i6,1P10E15.6)") it,i,temp(i),&
                         eps(i)/eps_gf,eeps/eps_gf, &
                         (eeps-eps(i))/eps(i)
                    stop "con2prim: Problem with EOS 1"
                 endif
                 ! atmosphere handling:
                 if(rho(i).eq.atmo_rho) then
                    q(i,1) = rho(i)
                    q(i,2) = 0.0d0
                    q(i,3) = rho(i)*eps(i)
                 endif
                 
                 fp = pp - op(i)
                 temp1 = (q(i,3)+q(i,1)+op(i))**2 - q(i,2)**2
                 if (temp1.lt.0.0d0) then
                    write(*,*) "temp less then zero"
                    stop
                 endif
                 drhodpress = q(i,1)*q(i,2)**2/(sqrt(temp1)*(q(i,3)+q(i,1)+op(i))**2)
                 dedpress = op(i)*q(i,2)**2/(rrho*(q(i,1)+q(i,3)+op(i))*temp1)

                 dfdp = dpdrh*drhodpress+dpde*dedpress-1.0d0

                 if (dfdp.ne.0.0d0) then
                    press(i) = op(i)-fp/dfdp
                 else
                    stop "Problems in dfdp"
                 endif
                 err = abs(1.0d0-press(i)/op(i))
              enddo
              if(it.ge.max_iterations) then
                 do while(success.eq.0.and.pt_counter.lt.7)
                    pt_counter = pt_counter + 1
                    low_tol = low_tol*10.0d0
                    press(i) = old_press(i)
                    call con2prim_pt(low_tol,i,success)
                    if (success.eq.0) then 
                       write(*,*) "con2prim failed with tol: ", low_tol, i, temp(i), eeps/eps_gf, eps(i)/eps_gf
                    endif
                    !if successful, press(i) is set
                 enddo
                 if(pt_counter.eq.7.and.success.eq.0) then
                    write(*,*) "con2prim failed with tol > 1.0d-4", i
                    stop "con2prim problem: iteration on tol"
                 else
                    write(*,*) "Success with smaller tol: ", low_tol, i
                    pt_counter = 0
                    success = 0
                    low_tol = tol
                 endif
              endif
              err = 1.0d0  
              v(i) = q(i,2)/(q(i,3)+q(i,1)+press(i))
              v1(i) = v(i)/X(i)
              if (v1(i).ne.v1(i)) then
                 write(*,*) nt,i,it, q(i,2), q(i,3), q(i,1), press(i), X(i), rho(i), eps(i)
                 stop "NaN in V1"
              endif
              discrim = 1.0d0-v(i)**2
              if (discrim.lt.0.0d0) then
                 stop "We have a problem in con2prim"
              endif
              W(i) = 1.0d0/sqrt(discrim)
              rho(i) = q(i,1)/X(i)/W(i)
              eps(i) = (q(i,3)+q(i,1)+press(i)*(1.0d0-W(i)**2))/(rho(i)*W(i)**2)-1.0d0
	   endif
        endif
     enddo

  else if (GR.and.do_rotation) then
     do i=iminb,imaxb 

        err = 1.0d0
        low_tol = tol
        ye(i) = q(i,4)/q(i,1)
#ifdef HAVE_BURN
        do k=1,nspec
           Yion(k,i) = max(q(i,6+k)/q(i,1), 1.0d-50)
        enddo
#endif

        if (activate_turbulence) then
            v_turb(i) = sqrt(q(i,6)/q(i,1))
        endif

        if (q(i,1).eq.0.0d0) then
           ! D = 0, set everything to zero in this case (can't be good)
           v1(i) = 0.0d0
           v(i) = 0.0d0
           vphi(i) = 0.0d0
           rho(i) = 0.0d0
           eps(i) = 0.0d0
	   press(i) = 0.0d0
           W(i) = 1.0d0
        else
           ! atmosphere handling:
           if(rho(i).eq.atmo_rho) then
              q(i,1) = rho(i)
              q(i,2) = 0.0d0
              q(i,3) = rho(i)*eps(i)
              q(i,5) = 0.0d0
              W(i) = 1.0d0
           endif
           ! special case in which both radial and
           ! angular momenta are zero
 	   if (q(i,2).eq.0.0d0.and.q(i,5).eq.0.0d0) then
              v1(i) = 0.0d0
              v(i) = 0.0d0
              vphi(i) = 0.0d0
              rho(i) = q(i,1)/X(i)
              eps(i) = (q(i,3)+q(i,1)-q(i,1)/X(i))/rho(i)
              W(i) = 1.0d0
              if (eps(i).lt.1.d-10) then
                 eps(i)=1.d-10
                 write(*,*) 'Help 2!!',i
                 stop
              endif
              keytemp = 0

              ! Limit eos variables, very drastic approach
              call ApplyEOS_limits_zone(i)
              call eos(i,rho(i),temp(i),ye(i),eps(i),pp, keytemp,keyerr,1, eoskey,eos_rf_prec)
              press(i) = pp
           else
              ! general case
              it = 0
              old_press(i) = press(i)
              do while (err.gt.tol.and.it.lt.max_iterations)
                 it = it + 1
                 op(i) = press(i)
                 !here vv is romero's v not v^r
                 vv = q(i,2)/(q(i,3)+q(i,1)+op(i))
                 !here vpv is the physical phi velocity
                 vpv = q(i,5)/(q(i,3)+q(i,1)+op(i))/x1(i)
                 if (vv.gt.1.0d0.or.vv.lt.-1.0d0) then
                    write(6,*) "We have a problem finding v:" 
                    write(6,*) "Timestep: ", nt
                    write(6,"(2i6,1P10E15.6)") it,i,rrho/rho_gf,rho(i)/rho_gf,eeps/eps_gf,ww
                    write(6,"(2i6,1P10E15.6)") it,i,q(i,2), (q(i,3)+q(i,1)+op(i)),vv,q(i,3), q(i,1), op(i)
                    call flush(6)
                    stop "con2prim problem"
                 endif
                 if (vpv.gt.1.0d0.or.vpv.lt.-1.0d0) then
                    write(6,*) "We have a problem finding vphi:" 
                    write(6,*) "Timestep: ", nt
                    write(6,"(2i6,1P10E15.6)") it,i,rrho/rho_gf,rho(i)/rho_gf,eeps/eps_gf,ww
                    write(6,"(2i6,1P10E15.6)") it,i,q(i,5), (q(i,3)+q(i,1)+op(i)),vpv,& 
                         q(i,3),q(i,1), op(i)
                    call flush(6)
                    stop "con2prim problem with vphi"
                 endif
                 discrim = 1.0d0-vv**2-twothirds*vpv**2
                 if (discrim.lt.0.0d0) then
                    write(*,*) "We have a problem with W", discrim, vv
                    stop
                 endif
                 ww = 1.0d0/sqrt(discrim)
                 rrho = q(i,1)/X(i)/ww
                 eeps = (q(i,3)+q(i,1)+op(i)*(1.0d0-ww**2))/(rrho*ww**2)-1.0d0
                 keytemp = 0
                 
                 ! No need to limit EOS variables here
                 call eos_full(i,rrho,temp(i),ye(i),eeps,pp, & 
                      eosdummy,eosdummy,eosdummy,eosdummy,&
                      dpde,dpdrh, &
                      eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,&
                      eosdummy,eosdummy,eosdummy,eosdummy, &
                      keytemp,keyerr,eoskey,eos_rf_prec)
                 if(keyerr.ne.0) then
                    write(6,*) "#############################################"
                    write(6,"(i4,i6,1P10E15.6)") it,i,x1(i)/length_gf,&
                         rrho/rho_gf,rho(i)/rho_gf,temp(i),ye(i),q(i,1)
                    write(6,"(i4,i6,1P10E15.6)") it,i,temp(i),&
                         eps(i)/eps_gf,eeps/eps_gf, &
                         (eeps-eps(i))/eps(i)
                    stop "con2prim: Problem with EOS 1"
                 endif
                 ! atmosphere handling:
                 if(rho(i).eq.atmo_rho) then
                    q(i,1) = rho(i)
                    q(i,2) = 0.0d0
                    q(i,3) = rho(i)*eps(i)
                    q(i,5) = 0.0d0
                 endif
                 
                 fp = pp - op(i)
                 temp1 = (q(i,3)+q(i,1)+op(i))**2 - ( q(i,2)**2 + twothirds*(q(i,5)/x1(i))**2 )
                 if (temp1.lt.0.0d0) then
                    write(*,*) "temp less then zero"
                    stop
                 endif
                 drhodpress = q(i,1)*(q(i,2)**2 + twothirds*(q(i,5)/x1(i))**2) / &
                      (sqrt(temp1)*(q(i,3)+q(i,1)+op(i))**2)
                 dedpress = op(i)*(q(i,2)**2 + twothirds*(q(i,5)/x1(i))**2) / &
                      (rrho*(q(i,1)+q(i,3)+op(i))*temp1)

                 dfdp = dpdrh*drhodpress+dpde*dedpress-1.0d0

                 if (dfdp.ne.0.0d0) then
                    press(i) = op(i)-fp/dfdp
                 else
                    stop "Problems in dfdp"
                 endif
                 err = abs(1.0d0-press(i)/op(i))	
              enddo
              if(it.ge.max_iterations) then
                 do while(success.eq.0.and.pt_counter.lt.7)
                    pt_counter = pt_counter + 1
                    low_tol = low_tol*10.0d0
                    press(i) = old_press(i)
                    call con2prim_pt_rot(low_tol,i,success)
                    if (success.eq.0) then 
                       write(*,*) "con2prim failed with tol: ", low_tol, i
                    endif
                    !if successful, press(i) is set
                 enddo
                 if(pt_counter.eq.7.and.success.eq.0) then
                    write(*,*) "con2prim failed with tol > 1.0d-4", i
                    stop "con2prim problem: iteration on tol"
                 else
                    write(*,*) "Success with smaller tol: ", low_tol, i
                    pt_counter = 0
                    success = 0
                    low_tol = tol
                 endif
              endif
              err = 1.0d0  
              v(i) = q(i,2)/(q(i,3)+q(i,1)+press(i))
              v1(i) = v(i)/X(i)
              vphi(i) = q(i,5)/(q(i,3)+q(i,1)+press(i))/x1(i)
              if (v1(i).ne.v1(i)) then
                 write(*,*) nt,i,it, q(i,2), q(i,3), q(i,1), press(i), X(i), rho(i), eps(i)
                 stop "con2prim: NaN in V1"
              endif
              if (vphi(i).ne.vphi(i)) then
                 write(*,*) nt,i,it, q(i,2), q(i,3), q(i,1), press(i), X(i)
                 stop "con2prim: NaN in vphi"
              endif
              discrim = 1.0d0 - (v(i)**2 + twothirds*vphi(i)**2)
              if (discrim.lt.0.0d0) then
                 stop "We have a problem in con2prim"
              endif
              W(i) = 1.0d0/sqrt(discrim)
              rho(i) = q(i,1)/X(i)/W(i)
              eps(i) = (q(i,3)+q(i,1)+press(i)*(1.0d0-W(i)**2))/(rho(i)*W(i)**2)-1.0d0

	   endif
         endif
      enddo
   else   
      rho(iminb:imaxb) = q(iminb:imaxb,1) 
      v1(iminb:imaxb)  =  q(iminb:imaxb,2) / q(iminb:imaxb,1)
      eps(iminb:imaxb) = q(iminb:imaxb,3)/rho(iminb:imaxb) & 
           - 0.5d0*(v1(iminb:imaxb)**2)
      eps_kin(:) = 0.5d0 * v1(:)**2
      if(do_rotation) then
         vphi1(iminb:imaxb) = q(iminb:imaxb,5)/q(iminb:imaxb,1)/x1(iminb:imaxb)
         eps(iminb:imaxb) = eps(iminb:imaxb) &
              - 0.5d0*twothirds*vphi1(iminb:imaxb)**2
         eps_kin(iminb:imaxb) = eps_kin(iminb:imaxb) + 0.5d0 * twothirds*vphi1(iminb:imaxb)**2
      endif
      
      ye(iminb:imaxb) = q(iminb:imaxb,4)/q(iminb:imaxb,1)
      if (activate_turbulence) then
         v_turb(iminb:imaxb) = sqrt(q(iminb:imaxb,6)/q(iminb:imaxb,1))
      endif
   endif
   
   ! a few checks
   if(GR) then
      do i=1,n1
         ! rho better be larger than 0
         if(rho(i).le.0.0d0) then
            write(6,*) "Density <= 0!!!"
            write(6,"(i8,1P10E15.6)") i,x1(i),rho(i)
            stop "Fix me please!"
         endif
      enddo
   endif

end subroutine con2prim_1

!*************************************************************************

subroutine con2prim_pt(tol,i,success)

  use GR1D_module
  use atmos
#ifdef HAVE_BURN
  use composition, only: nspec
#endif
  implicit none

#ifdef HAVE_BURN
  integer k
#endif
  real*8 tol, err, h, discrim
  real*8 pplus,pminus
  real*8 dpdrh, dpde,dedpress,drhodpress,temp1
  real*8 vv,rrho,eeps,pp,ww,op(n1),fp,dfdp
  real*8 :: t1 = 0.0d0
  real*8 :: t2 = 0.0d0
  integer i,j,it
  integer :: max_iterations = 1000
  integer success

  ! dummies for EOS call
  real*8 eosdummy
  integer keyerr,keytemp
  ! end dummies for EOS call

  err = 1.0d0

  !test in shocktube or similar
  if (GR.and.(gravity_active.eqv..false.)) then
     do j=1,n1
        if (X(j).ne.1.0d0) then
           write(*,*) "X is not 1"
        endif
        if (alp(j).ne.1.0d0) then
           write(*,*) "alp is not 1"
        endif  
     enddo
  endif

  if (GR) then
     err = 1.0d0
     ye(i) = q(i,4)/q(i,1)
#ifdef HAVE_BURN
     do k=1,nspec
        Yion(k,i) = max(q(i,6+k)/q(i,1), 1.0d-50)
     enddo
#endif

     if (activate_turbulence) then
        v_turb(i) = sqrt(q(i,6)/q(i,1))
     endif     

     it = 0
     do while (err.gt.tol.and.it.lt.max_iterations)
        it = it + 1
        op(i) = press(i)
        !here vv is romero's v not v^r
        vv = q(i,2)/(q(i,3)+q(i,1)+op(i))
        if (vv.gt.1.0d0.or.vv.lt.-1.0d0) then
           write(6,*) "We have a problem finding v:" 
           write(6,*) "Timestep: ", nt
           write(6,"(2i6,1P10E15.6)") it,i,rrho/rho_gf, &
                rho(i)/rho_gf,eeps/eps_gf,ww
           write(6,"(2i6,1P10E15.6)") it,i,q(i,2), (q(i,3) + &
                q(i,1)+op(i)),vv,q(i,3), q(i,1), op(i)
           call flush(6)
           stop "con2prim problem"
        endif
        discrim = 1.0d0-vv**2
        if (discrim.lt.0.0d0) then
           write(*,*) "We have a problem with W", discrim, vv
           stop
        endif
        ww = 1.0d0/sqrt(discrim)
        rrho = q(i,1)/X(i)/ww
        eeps = (q(i,3)+q(i,1)+op(i)*(1.0d0-ww**2))/(rrho*ww**2)-1.0d0
        
        keytemp = 0
        
        ! No need to limit EOS variables here
        call eos_full(i,rrho,temp(i),ye(i),eeps,pp, & 
             eosdummy,eosdummy,eosdummy,eosdummy,&
             dpde,dpdrh, &
             eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,&
             eosdummy,eosdummy,eosdummy,eosdummy, &
             keytemp,keyerr,eoskey,eos_rf_prec)

        WRITE(*,*) i, rrho,temp(i),ye(i),eeps/eps_gf
        if(keyerr.ne.0) then
           write(6,*) "#############################################"
           write(6,"(i4,i6,1P10E15.6)") it,i,x1(i)/length_gf,&
                rrho/rho_gf,rho(i)/rho_gf,temp(i),ye(i)
           write(6,"(i4,i6,1P10E15.6)") it,i,temp(i),&
                eps(i)/eps_gf,eeps/eps_gf, &
                (eeps-eps(i))/eps(i)
           stop "con2prim: Problem with EOS 1"
        endif
        ! atmosphere handling:
        if(rho(i).eq.atmo_rho) then
           q(i,1) = rho(i)
           q(i,2) = 0.0d0
           q(i,3) = rho(i)*eps(i)
        endif
        
        fp = pp - op(i)
        temp1 = (q(i,3)+q(i,1)+op(i))**2 - q(i,2)**2
        if (temp1.lt.0.0d0) then
           write(*,*) "temp less then zero"
           stop
        endif
        drhodpress = q(i,1)*q(i,2)**2/(sqrt(temp1)*(q(i,3)+q(i,1)+op(i))**2)
        dedpress = op(i)*q(i,2)**2/(rrho*(q(i,1)+q(i,3)+op(i))*temp1)
        
        dfdp = dpdrh*drhodpress+dpde*dedpress-1.0d0
        if (dfdp.ne.0.0d0) then
           press(i) = op(i)-fp/dfdp
        else
           stop "Problems in dfdp"
        endif
        err = abs(1.0d0-press(i)/op(i))	
     enddo
     if(it.ge.max_iterations) then
        success = 0
        return
     endif
     success = 1
     err = 1.0d0  
     v(i) = q(i,2)/(q(i,3)+q(i,1)+press(i))
     v1(i) = v(i)/X(i)
     if (v1(i).ne.v1(i)) then
        write(*,*) i,it, q(i,2), q(i,3), q(i,1), press(i), X(i),eps(i),rho(i)
        stop "Nan in V1"
     endif
     discrim = 1.0d0-v(i)**2
     if (discrim.lt.0.0d0) then
        stop "We have a problem in con2prim"
     endif
     W(i) = 1.0d0/sqrt(discrim)
     rho(i) = q(i,1)/X(i)/W(i)
     eps(i) = (q(i,3)+q(i,1)+press(i)*(1.0d0-W(i)**2))/(rho(i)*W(i)**2)-1.0d0

  else
     stop "Shouldn't be here in con2prim_pt"
  endif

end subroutine con2prim_pt

!*************************************************************************

subroutine con2prim_pt_rot(tol,i,success)

  use GR1D_module
  use atmos
#ifdef HAVE_BURN
  use composition, only: nspec
#endif
  implicit none

#ifdef HAVE_BURN
  integer k
#endif
  real*8 tol, err, h, discrim
  real*8 pplus,pminus
  real*8 dpdrh, dpde,dedpress,drhodpress,temp1
  real*8 vv,vpv,rrho,eeps,pp,ww,op(n1),fp,dfdp
  real*8 :: t1 = 0.0d0
  real*8 :: t2 = 0.0d0
  integer i,j,it
  integer :: max_iterations = 1000
  integer success

  ! dummies for EOS call
  real*8 eosdummy
  integer keyerr,keytemp
  ! end dummies for EOS call

  err = 1.0d0

  if (GR) then
     ye(i) = q(i,4)/q(i,1)
#ifdef HAVE_BURN
     do k=1,nspec
        Yion(k,i) = max(q(i,6+k)/q(i,1), 1.0d-50)
     enddo
#endif

     if (activate_turbulence) then
        v_turb(i) = sqrt(q(i,6)/q(i,1))
     endif
   
     ! general case
     it = 0
     op(i) = press(i)
     do while (err.gt.tol.and.it.lt.max_iterations)
        it = it + 1
        op(i) = press(i)
        !here vv is romero's v not v^r
        vv = q(i,2)/(q(i,3)+q(i,1)+op(i))
        !here vpv is the physical phi velocity
        vpv = q(i,5)/(q(i,3)+q(i,1)+op(i))/x1(i)
        if (vv.gt.1.0d0.or.vv.lt.-1.0d0) then
           write(6,*) "We have a problem finding v:" 
           write(6,*) "Timestep: ", nt
           write(6,"(2i6,1P10E15.6)") it,i,rrho/rho_gf,rho(i)/rho_gf,eeps/eps_gf,ww
           write(6,"(2i6,1P10E15.6)") it,i,q(i,2), (q(i,3)+q(i,1)+op(i)),vv,q(i,3), q(i,1), op(i)
           call flush(6)
           stop "con2prim problem"
        endif
        if (vpv.gt.1.0d0.or.vpv.lt.-1.0d0) then
           write(6,*) "We have a problem finding vphi:" 
           write(6,*) "Timestep: ", nt
           write(6,"(2i6,1P10E15.6)") it,i,rrho/rho_gf,rho(i)/rho_gf,eeps/eps_gf,ww
           write(6,"(2i6,1P10E15.6)") it,i,q(i,5), (q(i,3)+q(i,1)+op(i)),vpv,& 
                q(i,3),q(i,1), op(i)
           call flush(6)
           stop "con2prim problem with vphi"
        endif
        discrim = 1.0d0-vv**2-twothirds*vpv**2
        if (discrim.lt.0.0d0) then
           write(*,*) "We have a problem with W", discrim, vv
           stop
        endif
        ww = 1.0d0/sqrt(discrim)
        rrho = q(i,1)/X(i)/ww
        eeps = (q(i,3)+q(i,1)+op(i)*(1.0d0-ww**2))/(rrho*ww**2)-1.0d0
        keytemp = 0
        
        ! No need to limit EOS variables here
        call eos_full(i,rrho,temp(i),ye(i),eeps,pp, & 
             eosdummy,eosdummy,eosdummy,eosdummy,&
             dpde,dpdrh, &
             eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,eosdummy,&
             eosdummy,eosdummy,eosdummy,eosdummy, &
             keytemp,keyerr,eoskey,eos_rf_prec)
        if(keyerr.ne.0) then
           write(6,*) "#############################################"
           write(6,"(i4,i6,1P10E15.6)") it,i,x1(i)/length_gf,&
                rrho/rho_gf,rho(i)/rho_gf,temp(i),ye(i),q(i,1)
           write(6,"(i4,i6,1P10E15.6)") it,i,temp(i),&
                eps(i)/eps_gf,eeps/eps_gf, &
                (eeps-eps(i))/eps(i)
           stop "con2prim: Problem with EOS 1"
        endif
        
        fp = pp - op(i)
        temp1 = (q(i,3)+q(i,1)+op(i))**2 - ( q(i,2)**2 + twothirds*(q(i,5)/x1(i))**2 )
        if (temp1.lt.0.0d0) then
           write(*,*) "temp less then zero"
           stop
        endif
        drhodpress = q(i,1)*(q(i,2)**2 + twothirds*(q(i,5)/x1(i))**2) / &
             (sqrt(temp1)*(q(i,3)+q(i,1)+op(i))**2)
        dedpress = op(i)*(q(i,2)**2 + twothirds*(q(i,5)/x1(i))**2) / &
             (rrho*(q(i,1)+q(i,3)+op(i))*temp1)
        
        dfdp = dpdrh*drhodpress+dpde*dedpress-1.0d0
        
        if (dfdp.ne.0.0d0) then
           press(i) = op(i)-fp/dfdp
        else
           stop "Problems in dfdp"
        endif
        err = abs(1.0d0-press(i)/op(i))	
     enddo
     if(it.ge.max_iterations) then
        success = 0
        return
     endif
     err = 1.0d0  
     v(i) = q(i,2)/(q(i,3)+q(i,1)+press(i))
     v1(i) = v(i)/X(i)
     vphi(i) = q(i,5)/(q(i,3)+q(i,1)+press(i))/x1(i)
     if (v1(i).ne.v1(i)) then
        write(*,*) nt,i,it, q(i,2), q(i,3), q(i,1), press(i), X(i), rho(i), eps(i)
        stop "con2prim: NaN in V1"
     endif
     if (vphi(i).ne.vphi(i)) then
        write(*,*) nt,i,it, q(i,2), q(i,3), q(i,1), press(i), X(i)
        stop "con2prim: NaN in vphi"
     endif
     discrim = 1.0d0 - (v(i)**2 + twothirds*vphi(i)**2)
     if (discrim.lt.0.0d0) then
        stop "We have a problem in con2prim"
     endif
     W(i) = 1.0d0/sqrt(discrim)
     rho(i) = q(i,1)/X(i)/W(i)
     eps(i) = (q(i,3)+q(i,1)+press(i)*(1.0d0-W(i)**2))/(rho(i)*W(i)**2)-1.0d0

  else
     stop "Shouldn't be here in con2prim_pt_rot"
  endif


end subroutine con2prim_pt_rot

subroutine con2prim_grmhd

  use GR1D_module
  use atmos
  use omp_lib
  implicit none

  ! ---------------------------------------------------------------------------
  ! Primitive recovery via the Kastaun, Kalinani & Ciolfi (2021) master-variable
  ! scheme (Phys. Rev. D 103, 023018; arXiv:2005.01821).  Recovery is reduced to
  ! a single scalar root-find in the master variable  mu = 1/(h*W),  using ONLY
  ! p(rho,eps,Ye) from the EOS -- no thermodynamic derivatives (dpde, dpdrh).
  ! This reproduces the physics of con2prim_1 (Romero pressure iteration) but is
  ! GRMHD-style and derivative-free, which is what we want for a hybrid
  ! tabulated/analytic nuclear EOS where derivatives are noisy or discontinuous.
  !
  ! This routine is self-contained and is NOT wired into the con2prim dispatcher;
  ! swapping it in for con2prim_1 is a separate task (see CLAUDE.md, this dir).
  !
  ! GR1D is pure hydro.  The magnetic field is carried as subroutine-local
  ! variables that are identically zero, and the full Kastaun magnetic
  ! expressions are written out in the hat-functions below.  With B = 0 they
  ! reduce to the hydro limit (Chi = 1, all magnetic terms vanish) and the
  ! compiler constant-folds them away -- so the source reads as a genuine GRMHD
  ! routine at zero runtime cost (see CLAUDE.md sec.6).
  !
  ! GR1D uses the Romero orthonormal formulation, so in the orthonormal (tetrad)
  ! frame the 3-metric is the identity delta_ij: every dot product is Euclidean
  ! and the only non-zero momentum component is the radial one.  The metric
  ! factor X enters in exactly two places: building the undensitized
  ! conservatives, and converting the recovered orthonormal v back to v1 = v^r.
  ! ---------------------------------------------------------------------------

  integer iminb, imaxb
  integer i

  iminb = ghosts1+1
  imaxb = n1-ghosts1

  ! shocktube / no-gravity sanity check (mirrors con2prim_1): in such runs the
  ! metric is trivial, X = 1 and alp = 1, which the orthonormal mapping assumes.
  if (GR.and.(gravity_active.eqv..false.)) then
     do i=1,n1
        if (X(i).ne.1.0d0)   stop "con2prim_grmhd: X is not 1"
        if (alp(i).ne.1.0d0) stop "con2prim_grmhd: alp is not 1"
     enddo
  endif

  ! This GRMHD path is GR, non-rotating only (CLAUDE.md sec.11).
  if (.not.GR .or. do_rotation) then
     stop "con2prim_grmhd: only the GR non-rotating branch is implemented"
  endif

  ! Every zone solve is independent: all shared arrays are touched only at
  ! index i.  The solve itself (with its contained hat-functions) lives in
  ! con2prim_grmhd_zone so that each call -- hence each thread -- owns a
  ! private host frame for the per-zone constants: OpenMP private clauses do
  ! NOT apply to host-associated references inside contained procedures, so
  ! the hat-functions must not live in the loop's own host scope.
  ! schedule(dynamic): per-zone cost varies wildly (closed-form rest branch
  ! vs. many tabulated-EOS root-finds in the Newton solve).
  !$omp parallel do schedule(dynamic)
  do i=iminb,imaxb
     call con2prim_grmhd_zone(i)
  enddo
  !$omp end parallel do

  ! a few checks (mirror con2prim_1)
  if (GR) then
     do i=1,n1
        if (rho(i).le.0.0d0) then
           write(6,*) "Density <= 0!!!"
           write(6,"(i8,1P10E15.6)") i,x1(i),rho(i)
           stop "Fix me please!"
        endif
     enddo
  endif

end subroutine con2prim_grmhd

!*************************************************************************

! Kastaun solve for a single zone i (see the header of con2prim_grmhd).
! Deliberately a separate subroutine rather than the body of the driver's
! loop: the per-zone constants below are locals of THIS call, and the
! contained hat-functions reach them through host association -- which always
! binds to the host's own (per-call, hence per-thread) frame.
subroutine con2prim_grmhd_zone(i)

  use GR1D_module
  use atmos
#ifdef HAVE_BURN
  use composition, only: nspec
#endif
  implicit none

  integer, intent(in) :: i

#ifdef HAVE_BURN
  integer k
#endif

  integer it
  integer keytemp, keyerr

  ! per-zone constants (fixed during the mu solve; seen by the contained
  ! hat-functions through host association -- mirrors CLAUDE.md sec.8)
  real*8 :: D, q_K, r, r2, Ye_fix
  real*8 :: B_i, B_cons_r, B_cons2, rdotB, B_cons2r2_perp
  real*8 :: h_0, v0_2

  ! solver locals
  real*8 :: mu_plus, mu_root, pp, discrim
  real*8 :: lo, hi, mid
  real*8 :: fa, fb, fmu, fmu2, dmu, df, munew
  real*8 :: a_il, b_il
  logical :: converged

  integer, parameter :: nrmax = 200      ! Newton-Raphson iteration cap
  integer, parameter :: ilmax = 200      ! Illinois iteration cap
  real*8,  parameter :: ftol  = 1.0d-11  ! residual tolerance on f(mu)

     ! D == 0: nothing to recover, zero everything and skip the solve.
     if (q(i,1).eq.0.0d0) then
        v1(i)  = 0.0d0
        v(i)   = 0.0d0
        rho(i) = 0.0d0
        eps(i) = 0.0d0
        press(i) = 0.0d0
        W(i)   = 1.0d0
        return
     endif

     ! Ye is fixed for the whole solve; composition / turbulence as con2prim_1.
     ye(i) = q(i,4)/q(i,1)
#ifdef HAVE_BURN
     do k=1,nspec
        Yion(k,i) = max(q(i,6+k)/q(i,1), 1.0d-50)
     enddo
#endif
     if (activate_turbulence) then
        v_turb(i) = sqrt(q(i,6)/q(i,1))
     endif

     ! atmosphere reset (cf. con2prim_1): forces a static, zero-momentum state,
     ! which then routes to the closed-form rest branch just below.
     if (rho(i).eq.atmo_rho) then
        q(i,1) = rho(i)
        q(i,2) = 0.0d0
        q(i,3) = rho(i)*eps(i)
        W(i)   = 1.0d0
     endif

     ! Static, no radial momentum: closed-form rest branch (cheaper than the
     ! root-find and matches con2prim_1).  The Kastaun solve also handles this.
     if (q(i,2).eq.0.0d0) then
        v1(i)  = 0.0d0
        v(i)   = 0.0d0
        W(i)   = 1.0d0
        rho(i) = q(i,1)/X(i)
        eps(i) = (q(i,3)+q(i,1)-q(i,1)/X(i))/rho(i)
        keytemp = 0
        ! Note: no hard eps floor here -- with the composite EOS eps may carry a
        ! (possibly negative) nuclear zero-point offset; ApplyEOS_limits + the
        ! EOS keyerr path handle out-of-range states.
        call ApplyEOS_limits_zone(i)
        call eos(i,rho(i),temp(i),ye(i),eps(i),pp,keytemp,keyerr,1,eoskey,eos_rf_prec)
        if (keyerr.ne.0) then
           write(6,*) "con2prim_grmhd: EOS error (rest branch)", i, &
                rho(i)/rho_gf, temp(i), ye(i), eps(i)/eps_gf
           stop "con2prim_grmhd: Problem with EOS (rest)"
        endif
        press(i) = pp
        return
     endif

     ! ----------------------------------------------------------------------
     ! (i) Rescale: build the undensitized, orthonormal conservatives and the
     !     per-zone constants used throughout the solve (CLAUDE.md sec.4).
     !     Romero identity: q(i,3)+q(i,1) = rho*h*W**2 - press.
     ! ----------------------------------------------------------------------
     D   = q(i,1)/X(i)                              ! = W*rho  (undensitized density)
     q_K = (q(i,3)+q(i,1)-q(i,1)/X(i))/D            ! = tau/D
     r   = X(i)*q(i,2)/q(i,1)                       ! = S_r/D  (only non-zero component)
     r2  = r*r
     Ye_fix = ye(i)

     ! Magnetic field: identically zero in GR1D, but computed for clarity so the
     ! hat-functions below are the genuine GRMHD expressions (CLAUDE.md sec.6).
     B_i            = 0.0d0
     B_cons_r       = B_i/sqrt(D)                   ! orthonormal B^r / sqrt(D)
     B_cons2        = B_cons_r*B_cons_r             ! = 0
     rdotB          = r*B_cons_r                    ! Euclidean dot, = 0
     B_cons2r2_perp = B_cons2*r2 - rdotB**2         ! = 0

     ! Enthalpy floor h_0 > 0 (CLAUDE.md sec.5).  It must be a true lower bound
     ! on the specific enthalpy h = 1 + eps + p/rho for the bracket to contain
     ! the root; a smaller value is conservative (wider mu interval).  For the
     ! analytic EOS (ideal/poly/hybrid) h >= 1, so h_0 = 1 is exact.  For the
     ! tabulated/composite EOS eps may carry a small (possibly negative) nuclear
     ! zero-point offset, so a small strictly-positive floor is used.
     if (eoskey.eq.1 .or. eoskey.eq.2 .or. eoskey.eq.4) then
        h_0 = 1.0d0
     else
        h_0 = 1.0d-2
     endif
     v0_2 = r2/(h_0**2 + r2)                        ! velocity cap v0^2 (precomputed)

     ! ----------------------------------------------------------------------
     ! (ii) Bracket: upper bound mu_plus = root of f_a(mu) in (0, 1/h_0],
     !      f_a(mu) = mu*sqrt(h_0^2 + r_bar2(mu)) - 1   (Kastaun eq. 59).
     !      f_a is monotone increasing there (f_a(0) = -1 < 0, f_a(1/h_0) >= 0)
     !      and needs no EOS call, so a short bisection is robust and cheap.
     ! ----------------------------------------------------------------------
     lo = 0.0d0
     hi = 1.0d0/h_0
     do it=1,100
        mid = 0.5d0*(lo+hi)
        if (fbracket(mid) .gt. 0.0d0) then
           hi = mid
        else
           lo = mid
        endif
        if (hi-lo .le. 1.0d-15*hi) exit
     enddo
     mu_plus = hi

     ! ----------------------------------------------------------------------
     ! (iii) Master solve: root of f(mu) in (0, mu_plus]  (Kastaun eq. 62).
     !       f(mu) = mu - 1/(nu + mu*r_bar2(mu)),  nu = max(vA(mu), vB(mu)).
     !       f changes sign on (0, mu_plus], so the root is bracketed.  Use
     !       Newton-Raphson (numerical derivative -> stays EOS-derivative-free)
     !       and fall back to the Illinois method if NR leaves the bracket or
     !       fails to converge.
     ! ----------------------------------------------------------------------
     fb = fmaster(mu_plus)
     if (fb .le. 0.0d0) then
        ! root is at (or beyond) the velocity cap -> take the cap.
        mu_root = mu_plus
     else
        ! --- Newton-Raphson ---
        mu_root   = 0.5d0*mu_plus
        converged = .false.
        do it=1,nrmax
           fmu = fmaster(mu_root)
           if (abs(fmu) .lt. ftol) then
              converged = .true.
              exit
           endif
           dmu  = max(1.0d-8*mu_root, 1.0d-14)
           fmu2 = fmaster(mu_root+dmu)
           df   = (fmu2-fmu)/dmu
           if (df .eq. 0.0d0) exit
           munew = mu_root - fmu/df
           if (munew .le. 0.0d0 .or. munew .gt. mu_plus) exit   ! left the bracket
           mu_root = munew
        enddo

        ! --- Illinois (modified regula-falsi) fallback on [0, mu_plus] ---
        if (.not. converged) then
           a_il = 0.0d0
           b_il = mu_plus
           fa   = fmaster(a_il)
           fb   = fmaster(b_il)
           munew = b_il
           do it=1,ilmax
              munew = b_il - fb*(b_il-a_il)/(fb-fa)
              fmu   = fmaster(munew)
              if (fmu*fb .lt. 0.0d0) then
                 a_il = b_il
                 fa   = fb
              else
                 fa = 0.5d0*fa            ! Illinois weighting
              endif
              b_il = munew
              fb   = fmu
              if (abs(fmu) .lt. ftol .or. &
                   abs(b_il-a_il) .le. 1.0d-13*abs(b_il)) exit
           enddo
           mu_root = munew
           if (abs(fmaster(mu_root)) .gt. 1.0d-6) then
              write(6,*) "con2prim_grmhd: master solve did not converge", &
                   i, mu_root, fmaster(mu_root)
              stop "con2prim_grmhd: non-convergence"
           endif
        endif
     endif

     ! ----------------------------------------------------------------------
     ! (iv) Recover & store primitives at the root mu_root.  press via a final
     !      EOS call (which also lands temp(i) on its converged value).
     ! ----------------------------------------------------------------------
     W(i)     = W_hat(mu_root)
     rho(i)   = D/W(i)
     eps(i)   = eps_hat(mu_root)
     press(i) = p_hat(mu_root)
     v(i)     = mu_root*r                ! orthonormal velocity (mu*r = v identity)
     v1(i)    = v(i)/X(i)                ! coordinate velocity v^r

     ! Final guards (mirror con2prim_1).
     if (v1(i).ne.v1(i)) then
        write(6,*) nt,i, q(i,2), q(i,3), q(i,1), press(i), X(i), rho(i), eps(i)
        stop "con2prim_grmhd: NaN in V1"
     endif
     discrim = 1.0d0 - v(i)**2
     if (discrim.lt.0.0d0) then
        stop "con2prim_grmhd: v^2 >= 1"
     endif
     if (rho(i).le.0.0d0) then
        write(6,"(i8,1P10E15.6)") i, x1(i), rho(i)
        stop "con2prim_grmhd: density <= 0"
     endif

contains

  ! All "hat" quantities are functions of the single master variable mu; the
  ! per-zone constants (D, q_K, r, r2, h_0, v0_2, B_cons2, B_cons2r2_perp,
  ! rdotB, Ye, i) are fixed during the solve and reach these procedures through
  ! host association.  With B = 0 the magnetic terms vanish (Chi = 1).

  ! Kastaun eq. 61
  real*8 function Chi(mu)
    real*8, intent(in) :: mu
    Chi = 1.0d0/(1.0d0 + mu*B_cons2)
  end function Chi

  ! Kastaun eq. 60   (= r2 when B = 0)
  real*8 function r_bar2(mu)
    real*8, intent(in) :: mu
    real*8 :: c
    c = Chi(mu)
    r_bar2 = r2*c**2 + mu*c*(1.0d0 + c)*rdotB**2
  end function r_bar2

  ! Kastaun eq. 69   (= q_K when B = 0)
  real*8 function q_bar(mu)
    real*8, intent(in) :: mu
    real*8 :: c
    c = Chi(mu)
    q_bar = q_K - 0.5d0*B_cons2 - 0.5d0*mu**2*c**2*B_cons2r2_perp
  end function q_bar

  ! Kastaun eq. 68a  (velocity-capped)
  real*8 function v2_hat(mu)
    real*8, intent(in) :: mu
    v2_hat = min(mu**2*r_bar2(mu), v0_2)
  end function v2_hat

  ! Kastaun eq. 68b
  real*8 function W_hat(mu)
    real*8, intent(in) :: mu
    W_hat = 1.0d0/sqrt(1.0d0 - v2_hat(mu))
  end function W_hat

  ! Kastaun eq. 66c
  real*8 function rho_hat(mu)
    real*8, intent(in) :: mu
    rho_hat = D/W_hat(mu)
  end function rho_hat

  ! Kastaun eq. 67
  real*8 function eps_hat(mu)
    real*8, intent(in) :: mu
    real*8 :: wh
    wh = W_hat(mu)
    eps_hat = wh*(q_bar(mu) - mu*r_bar2(mu)) + v2_hat(mu)*wh**2/(1.0d0 + wh)
  end function eps_hat

  ! Kastaun eq. 66a -- the only EOS evaluation: p(rho_hat, eps_hat, Ye).
  ! Uses the lightweight `eos` wrapper with keytemp = 0 (eps known, temperature
  ! solved internally) and eosflag = 1 (return pressure).  No EOS derivatives.
  ! temp(i) is carried as the solver's seed and updated to the converged value.
  real*8 function p_hat(mu)
    real*8, intent(in) :: mu
    real*8  :: rr, ee, px
    integer :: kt, ke
    rr = rho_hat(mu)
    ee = eps_hat(mu)
    kt = 0
    ! No need to limit EOS variables here
    call eos(i, rr, temp(i), ye(i), ee, px, kt, ke, 1, eoskey, eos_rf_prec)
    if (ke.ne.0) then
       write(6,*) "con2prim_grmhd: EOS error in p_hat", i, &
            rr/rho_gf, temp(i), ye(i), ee/eps_gf, mu
    endif
    p_hat = px
  end function p_hat

  ! Kastaun eq. 66b
  real*8 function a_hat(mu)
    real*8, intent(in) :: mu
    a_hat = p_hat(mu)/(rho_hat(mu)*(1.0d0 + eps_hat(mu)))
  end function a_hat

  ! Kastaun eq. 59 -- bracketing function for mu_plus (no EOS call).
  real*8 function fbracket(mu)
    real*8, intent(in) :: mu
    fbracket = mu*sqrt(h_0**2 + r_bar2(mu)) - 1.0d0
  end function fbracket

  ! Kastaun eq. 62 -- master function.  One EOS call per evaluation (via a_hat).
  ! vA (eq. 64) and vB (eq. 65) are computed inline from a single p_hat so the
  ! root-find does not pay extra EOS evaluations.
  real*8 function fmaster(mu)
    real*8, intent(in) :: mu
    real*8 :: rb2, wh, eh, ah, qm, vA, vB, nu
    rb2 = r_bar2(mu)
    wh  = W_hat(mu)
    eh  = eps_hat(mu)
    ah  = a_hat(mu)
    qm  = q_bar(mu) - mu*rb2
    vA  = (1.0d0 + ah)*(1.0d0 + eh)/wh
    vB  = (1.0d0 + ah)*(1.0d0 + qm)
    nu  = max(vA, vB)
    fmaster = mu - 1.0d0/(nu + mu*rb2)
  end function fmaster

end subroutine con2prim_grmhd_zone
