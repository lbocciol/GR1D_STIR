!-*-f90-*-
subroutine con2prim
  
  use GR1D_module
  use timers
  implicit none

  real*8 :: t1, t2

  CALL GetThisTime(t1)
  call con2prim_1
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
              call ApplyEOS_limits
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
                 ! Limit eos variables, very drastic approach
                 call ApplyEOS_limits
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
              call ApplyEOS_limits
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
                 
                 ! Limit eos variables, very drastic approach
                 call ApplyEOS_limits
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
        
        ! Limit eos variables, very drastic approach
        call ApplyEOS_limits
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
        
        ! Limit eos variables, very drastic approach
        call ApplyEOS_limits
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
#ifdef HAVE_BURN
  use composition, only: nspec
#endif
  implicit none

  ! Notice that I am using Einstein notation here where _i in covariant 
  ! and ^i is contravariant

  ! (i) Calculate the rescaled variables and the following useful
  ! relations which are fixed during the iterations.
  D = X* rho* W ! i.e. q(i,1) in the code as you can see from above
  S_i = rho*h*W**2*v ! i.e. q(i,2) in the code as you can see from above
  tau = rho*h*W**2 - p - D ! i.e. q(i,3) in the code as you can see from above
  ! B^i will be manually set to zero

  q = tau/D
  r_i = S_i/D
  B_cons^i = B^i/SQRT(D) ! This is used for GRMHD, not used in GR1D where it is zero but let's evaluate it for consistency

  ! and then calculate 
  r2 = r^i * r_i
  B_cons2 = B_cons^i * B_cons_i
  B_cons2r2_perp = B_cons2 * r2 - (r^l * B_cons_l)**2

  ! (ii) In the interval (0, h_0**(-1) ], solve the following eq. with Newton–Raphson method,
  ! which is usually more efficient than bracketing methods. In case
  ! the Newton–Raphson method fails to converge, we use the Illinois
  ! algorithm to solve this equation
  f_a = mu * SQRT(h_0**2 + r_bar2(mu)) - 1.0d0 ! eq. 59

  ! (iii) In the interval (0, mu_plus ], solve
  f = mu - 1.0d0/(v + mu*r_bar2(mu)) ! eq. 62
  ! where v = ma(vA(mu, vB(mu)))
  ! we numerically solve equation (62) with the Illinois
  ! algorithm. Note that during the iterations, we enforce the density rho
  ! and the specific energy eps fall within the validity region of the EOS, i.e.
  ! we evaluate the updated rho and eps with rho_hat = max(min(rho_max, rho_hat), rho_min)
  ! and eps_hat = max(min(eps_max(rho_hat), eps_hat), eps_min(rho_hat)).

  ! iv) With the root mu of equation (62), we can then work out the
  ! primitive variables [rho, eps, p] respectively with the equations used in
  ! step 3. The velocity vi can be obtained with mu by

  v_hat^i(mu) = mu*Chi(mu)*(r^i + mu*(r^l B_cons_l)*B_cons^i)
end subroutine con2prim_grmhd

function vA(mu)
  ! eq 64
  (1.0d0 + a_hat(mu)) * ( 1.0d0 + eps_hat(mu))/W_hat(mu)
end function vA

function vB(mu)
  ! eq 65
  (1.0d0 + a_hat(mu)) * ( 1.0d0 + q_bar(mu) - mu*r_bar2(mu))
end function vB

function p_hat(mu)
  ! eq 66a
  p(rho_hat(mu), eps_hat(mu)) ! I guess this is the EOS part? But why is there no Ye?
end function p_hat

function a_hat(mu)
  ! eq 66b
  p_hat(mu) /(rho_hat(mu)*(1.0d0 + eps_hat(mu)))
end function a_hat

function rho_hat(mu)
  ! eq 66c
  D/W_hat(mu)
end function rho_hat

function eps_hat(mu)
  ! eq 67
  W_hat(mu) * (q_bar(mu) - mu*r_bar2(mu)) + v2_hat(mu) * W_hat(mu)**2/(1.0d0 + W_hat(mu))
end function eps_hat

function v2_hat(mu)
  ! eq 68a
  v02 = r2/(h_0**2 + r2)
  MIN(mu**2*r_bar2(mU), v02)
end function v2_hat

function W_hat(mu)
  ! eq 68b
  1.0d0/SQRT(1.0d0 - v2_hat(mu))
end function W_hat

function q_bar(mu)
  ! eq 69
  q - 0.5d0 * B_cons2 - 0.5d0*mu**2*Chi(mu)**2*B_cons2r2_perp
end function q_bar

function r_bar2(mu)
  ! eq 60
  r**2 * Chi(mu)**2 + mu*Chi(mu) * (1.0d0 + Chi(mu)) * (r^l*Beta_l)**2
end function r_bar2

function Chi(mu)
  ! eq 61
  1.0d0/(1.0d0 * mu*B_cons2)
end function Chi