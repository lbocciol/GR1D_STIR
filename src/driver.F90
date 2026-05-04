! -*-f90-*-
subroutine SetTimeStep
  
  use GR1D_module
  implicit none

  real*8 sound,dtnew
  integer keytemp,eosflag,keyerr
  integer i
  logical nan,inf
  
  if (GR) then
     call findnaninf(v,n1,nan,inf)
  else
     call findnaninf(v1,n1,nan,inf)
  endif
  
  if(nan) stop "NaNs found :-/"
    
  if (initial_data.eq."M1test".and.M1_testcase_number.gt.9) then
     stop "make sure dt works with grid"
     return
  endif

  ! sets time step
  ! Note: so far only cfl-stability conditions is implemented
  ! Other conditions may be added
  ! get the speed of sound from the eos
  do i=1,n1
     keytemp = 0 ! not coming in with temperature (that would reset the energy), needs to be zero for the hybrid/poly/ideal EOS
     eosflag = 6 ! we want cs2 to be reset
     keyerr = 0
     call eos(i,rho(i),temp(i),ye(i),eps(i),cs2(i), &
          keytemp,keyerr,eosflag,eoskey,eos_rf_prec)
     if(keyerr.ne.0) then
        stop "problem in eos: cs2 at SetTimeStep"
     endif
  enddo
  dtp = dt

  dtnew = max(1.0d0*time_gf,dt_max*time_gf)
  do i=ghosts1+1, n1-ghosts1
     if (do_M1) then
        sound = 1.0d0
     else
        sound = sqrt(cs2(i))
     endif
     if (GR) then 
        if(.not.do_rotation) then
           dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                max(abs(v(i) + sound),abs(v(i) - sound)))
        else 
           dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                max(max(abs(v(i) + sound),abs(v(i) - sound)),&
                max(abs(vphi(i)+sound),abs(vphi(i)-sound))))
        endif
     else
        if(.not.do_rotation) then
           dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                max(abs(v1(i) + sound),abs(v1(i) - sound)))
        else
           dtnew = min(dtnew,  (x1(i+1)-x1(i)) /  & 
                max(max(abs(v1(i) + sound),abs(v1(i) - sound)),&
                max(abs(vphi1(i)+sound),abs(vphi1(i)-sound))))
        endif
        
     endif
  enddo

  ! if not adaptive, set time step dt
  dt = min(dt_reduction_factor*cffac*dtnew,1.05d0*dtp)
  
  if (dtnew.eq.dt_max*time_gf) then
     dt_max = 1.01*dt_max
  endif

end subroutine SetTimeStep

subroutine handle_output
  
  use timers
  use GR1D_module
#ifdef HAVE_HDF5_OUTPUT
  use hdf5_output_utils
#else
#endif
  implicit none

  real*8 t_now

  call GetThisTime(t_now)
  
  if(dynamic_output_control) then
     call output_control
  endif
  if(mod(nt,ntinfo).eq. 0) then
     call outinfo
  endif
  
  if (nt.ge.ntmax) then
     write(*,*) "Done! :-) ntmax reached"
#ifdef HAVE_HDF5_OUTPUT
        call output_all_HDF5(1)  ! Grid output - opens/closes files for each variable
        call hdf5_increment_output_counter()
        call output_all_HDF5(2)
#else
        call output_all(1)  ! Grid output - opens/closes files for each variable
        call output_all(2)
#endif
     call restart_output_h5
     call PrintTimers()
     stop
  endif
  
  if (time.ge.tend) then
     write(*,*) "Done! :-) tend reached"
#ifdef HAVE_HDF5_OUTPUT
        call output_all_HDF5(1)  ! Grid output - opens/closes files for each variable
        call hdf5_increment_output_counter()
        call output_all_HDF5(2)
#else
        call output_all(1)  ! Grid output - opens/closes files for each variable
        call output_all(2)
#endif
     call restart_output_h5
     call PrintTimers()
     open(unit=666,file=trim(adjustl(outdir))//"/done",status="unknown")
     write(666,*) 1
     close(666)
     call PrintTimers()

     stop
  endif

  if ( (t_now - t_start) .ge. tend_wallclock ) then
    write(*,*) t_now, t_start, tend_wallclock
    write(*,*) "Done! :-) wallclock limit reached"
#ifdef HAVE_HDF5_OUTPUT
        call output_all_HDF5(1)  ! Grid output - opens/closes files for each variable
        call hdf5_increment_output_counter()
        call output_all_HDF5(2)
#else
        call output_all(1)  ! Grid output - opens/closes files for each variable
        call output_all(2)
#endif
    call PrintTimers()
    call restart_output_h5
    stop
  endif

  !!   Output/Checking
  if ( mod(nt,ntout) .eq. 0 .and. ntout.ne.-1) OutputFlag = .true.
  
  if ( mod(nt,ntout_restart) .eq. 0 .and. &
       ntout_restart.ne.-1) OutputFlagRestart = .true.
  
  if ( mod(nt,ntout_scalar).eq.0 .and. ntout_scalar.ne.-1) &
       OutputFlagScalar = .true.
  
  if ( time.ge.tdump) then
     tdump=tdump+dtout
     OutputFlag = .true.
  endif
  
  if ( time.ge.tdump_restart) then
     OutputFlagRestart = .true.
  endif
  
  if ( time.ge.tdump_scalar) then
     tdump_scalar=tdump_scalar+dtout_scalar
     OutputFlagScalar = .true.
  endif
  
  if(nt.eq.0) then
     OutputFlag = .true.
     OutputFlagScalar = .true.
  endif
  
  if (OutputFlag) then
#ifdef HAVE_HDF5_OUTPUT
        call output_all_HDF5(1)  ! Grid output - opens/closes files for each variable
        call hdf5_increment_output_counter()
#else
        call output_all(1)  ! Grid output - opens/closes files for each variable
#endif
     call output_timers
     OutputFlag = .false.
  endif
  
  if (OutputFlagRestart) then
     call restart_output_h5
     tdump_restart=tdump_restart+dtout_restart
     OutputFlagRestart = .false.
  endif
  
  if (OutputFlagScalar) then
#ifdef HAVE_HDF5_OUTPUT
        call output_all_HDF5(2)
#else
        call output_all(2)
#endif
    OutputFlagScalar = .false.
  endif

  if ((time+dt/time_gf).gt.tend) then
    dt = (tend-time)*time_gf
  endif
    
end subroutine handle_output

subroutine postStep_analysis
  
  use GR1D_module
  implicit none

  !for local EOS calls
  real*8 lrho,ltemp,lye,eosdummy(14)
  integer i,keytemp,keyerr
  
  if (initial_data.eq."Collapse".or.initial_data.eq."Collapse_inflow") then
     call get_shock_radius
     call mass_analysis    
     call get_binding_energy
     if (do_rotation) then
        call rotation_analysis
     endif
     if (M1_control_flag) then
        call M1_control
     endif
  endif
  
  nt=nt+1
  time = time+dt/time_gf
  time_c = time_c + dt/time_gf*alp(ghosts1+1)
  if (GR) then
     v_prev = v
  else
     v_prev = v1
  endif
  ye_prev = ye

  !fill mass fractions
#if HAVE_NUC_EOS
  if (eoskey.eq.3) then !if using nuclear EOS
     keytemp = 1 !keep eps constant
     keyerr = 0
     do i=ghosts1+1,n1-ghosts1
        lrho = rho(i)/rho_gf
        ltemp = temp(i)
        lye = ye(i)
        call ApplyEOS_limits
        call nuc_eos_full(lrho,ltemp,lye,eosdummy(1),eosdummy(2),eosdummy(3), &
             eosdummy(4),eosdummy(5),eosdummy(6),eosdummy(7),massfrac_a(i),massfrac_h(i), &
             massfrac_n(i),massfrac_p(i),massfrac_abar(i),massfrac_zbar(i),eosdummy(8), &
             eosdummy(9),eosdummy(10),eosdummy(11),keytemp,keyerr,eos_rf_prec)
        if(keyerr.ne.0) then
           write(6,*) "############################################"
           write(6,*) "EOS PROBLEM in poststep analysis:"
           write(6,*) "timestep number: ",nt
           write(6,"(i4,1P10E15.6)") i,x1(i)/length_gf,lrho,ltemp,lye
           stop "Shouldn't fair here...."
        endif
     enddo
  endif
#endif

end subroutine postStep_analysis

subroutine PrintTimers

  use GR1D_module
  use timers

  real*8 :: tfinal, total_M1

  CALL GetThisTime(tfinal)
  timer_code = tfinal - timer_code

  total_M1   = timer_M1_exp + timer_M1_imp + timer_M1_clo + timer_M1_rec + timer_M1_eas

  print *, '----------------- Timer Summary -----------------'
  print '(A,F12.4)', 'Total code time        = ', timer_code
  print '(A,F12.4)', 'Total step time        = ', timer_step
  print '(A,F12.4)', '  M1 (explicit)        = ', timer_M1_exp
  print '(A,F12.4)', '  M1 (implicit)        = ', timer_M1_imp
  print '(A,F12.4)', '  M1 (closure)         = ', timer_M1_clo
  print '(A,F12.4)', '  M1 (reconstruction)  = ', timer_M1_rec
  print '(A,F12.4)', '  M1 (updateeas)       = ', timer_M1_eas
  print '(A,F12.4)', '  M1 (cons. update)    = ', timer_M1cons
  print '(A,F12.4)', '  M1 (neutrinos only)  = ', total_M1
  print '(A,F12.4)', '  M1                   = ', timer_M1
  print '(A,F12.4)', '  EOS full             = ', timer_eosf
  print '(A,F12.4)', '  EOS short            = ', timer_eoss
  print '(A,F12.4)', '  Hydro                = ', timer_hydro
  print '(A,F12.4)', '  con2prim             = ', timer_c2p
  print '(A,F12.4)', '  con2GR               = ', timer_c2GR
  print '(A,F12.4)', '  reconstruction       = ', timer_rec
  print *
  print '(A,F6.2)',  'Fraction Neutrinos     = ', total_M1 / timer_step
  print '(A,F6.2)',  'Fraction Hydro         = ', timer_hydro / timer_step
  print *, '-------------------------------------------------'

end subroutine PrintTimers
